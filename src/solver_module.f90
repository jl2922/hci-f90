module solver_module

  use constants_module
  use det_module
  use lru_cache_module
  use linked_list_module__int
  use linked_list_module__double
  use spin_det_module
  use types_module
  use utilities_module
  use wavefunction_module

  implicit none

  private

  integer, parameter :: MAX_VAR_ITERATION = 20
  integer, parameter :: MIN_PT_ST_ITERATION = 10

  public :: solver_type

  type solver_type
    private
    logical, public :: dump_var_wf
    integer, public :: n_orb
    integer, public :: n_elec
    integer, public :: n_up
    integer, public :: n_dn
    integer, public :: n_states
    integer, public :: max_n_rs_pairs
    integer, public :: max_connected_dets
    real(DOUBLE), public :: eps_var
    real(DOUBLE), public :: eps_pt
    real(DOUBLE), public :: tol
    real(DOUBLE), public :: hf_energy
    real(DOUBLE), public :: var_energy
    real(DOUBLE), public :: pt_det_energy
    real(DOUBLE), public :: pt_st_energy
    real(DOUBLE), public :: pt_st_uncert
    real(DOUBLE), public :: max_abs_H
    type(optional_int), public :: n_samples
    type(optional_double), public :: eps_pt_det
    type(wavefunction_type), pointer, public :: wf => null()
    contains
      procedure, public :: solve
      procedure :: read_configs
      procedure :: read_system_configs
      procedure :: setup
      procedure :: get_hamiltonian_elem
      procedure :: get_next_dets
      procedure :: find_connected_dets
      procedure :: diagonalize
      procedure :: pt
      procedure :: pt_det
      procedure :: summary
  end type solver_type

  contains

  subroutine read_configs(this, config_file_unit)
    class(solver_type), intent(inout) :: this
    integer, intent(in) :: config_file_unit

    ! Default values.
    integer :: n_up = -1
    integer :: n_dn = -1
    integer :: n_elec
    integer :: n_states = 1
    real(DOUBLE) :: eps_var = 1.0e-3_DOUBLE
    real(DOUBLE) :: eps_pt = 1.0e-5_DOUBLE
    real(DOUBLE) :: tol = 1.0e-3_DOUBLE ! 1 mHa.
    logical :: dump_var_wf = .true.
    real(DOUBLE) :: eps_pt_det = -1 ! eps for the deterministic part of PT.
    real(DOUBLE) :: n_samples = -1 ! For stochastic PT. Real to allow sci format.
    integer :: io_err
    namelist /hci/ n_up, n_dn, n_states, &
        & eps_var, eps_pt, eps_pt_det, tol, dump_var_wf, n_samples

    rewind(config_file_unit)
    read(unit=config_file_unit, nml=hci, iostat=io_err)
    if (io_err > 0) stop 'Cannot read HCI configs.'
    if (n_up < 0 .or. n_dn < 0) then
      stop 'n_up and n_dn must be provided.'
    end if
    this%eps_var = eps_var
    this%eps_pt = eps_pt
    this%tol = tol
    this%dump_var_wf = dump_var_wf
    this%n_up = n_up
    this%n_dn = n_dn
    n_elec = n_up + n_dn
    this%n_elec = n_elec
    this%n_states = n_states
    if (eps_pt_det > 0) then
      this%eps_pt_det%instance = eps_pt_det
      this%eps_pt_det%is_present = .true.
    else
      this%eps_pt_det%is_present = .false.
    end if
    if (n_samples > 0) then
      this%n_samples%instance = nint(n_samples)
      this%n_samples%is_present = .true.
    else
      this%n_samples%is_present = .false.
    end if

    write (6, '(A)') '[CONFIGURATIONS]'
    write (6, '(A, I0, A, I0)') 'n_up/n_dn: ', n_up, '/', n_dn
    write (6, '(A, I0)') 'n_states: ', n_states
    write (6, '(A, F0.10)') 'eps_var: ', eps_var
    write (6, '(A, F0.10)') 'eps_pt: ', eps_pt
    write (6, '(A, F0.10)') 'tol: ', tol
    write (6, '(A, L1)') 'dump_var_wf: ', dump_var_wf
    if (this%eps_pt_det%is_present) then
      write (6, '(A, F0.10)') 'eps_pt_det: ', this%eps_pt_det%instance
    else
      write (6, '(A)') 'eps_pt_det: based on resources available.'
    end if
    if (this%n_samples%is_present) then
      write (6, '(A, I0)') 'n_samples: ', this%n_samples%instance
    else
      write (6, '(A)') 'n_samples: based on resources available.'
    end if

    call this%read_system_configs(config_file_unit)
  end subroutine read_configs

  subroutine read_system_configs(this, config_file_unit)
    class(solver_type), intent(inout) :: this
    integer, intent(in) :: config_file_unit

    write (6, '(A)') 'No system specific configurations.'
  end subroutine read_system_configs

  subroutine solve(this, config_file_unit)
    class(solver_type), intent(inout) :: this
    integer, intent(in) :: config_file_unit
    integer :: i 
    real(DOUBLE) :: energy_prev, energy_cur
    type(wavefunction_type), pointer :: wf_prev

    call this%read_configs(config_file_unit)

    write (6, '(A)') '[SETUP]'
    call this%setup()
    write (6, '()')

    write (6, '(A)') '[VARIATION]'
    write (6, '(A, F0.10)') 'eps_var: ', this%eps_var
    energy_prev = 0
    call build(wf_prev, this%wf)
    do i = 1, MAX_VAR_ITERATION
      write (6, '(A, I0)') 'Variation Iteration #', i
      call this%get_next_dets()

      write (6, '(A, G0.10)') 'Number of dets: ', this%wf%n
      if (this%wf%n < 1.01 * wf_prev%n) then
        this%wf = wf_prev
        write (6, '(A)') 'Number of dets change within 1%. Variation finished.'
        exit
      end if

      call this%diagonalize(energy_cur)
      this%var_energy = energy_cur
      write (6, '(A, F0.10)') 'Energy: ', energy_cur
      write (6, '()')
      call flush(6)

      if (abs(energy_prev - energy_cur) < 1.0e-6) then
        write (6, '(A)') 'Energy change within 1.0e-6. Variation finished.'
        exit
      end if

      energy_prev = energy_cur
      wf_prev = this%wf
    end do
    write (6, '()')

    write (6, '(A)') '[PERTURBATION]'
    write (6, '(A, F0.10)') 'eps_pt: ', this%eps_pt
    call this%pt()
    write (6, '()')

    write (6, '(A)') '[SUMMARY]'
    call this%summary()
  end subroutine solve

  subroutine setup(this)
    class(solver_type), intent(inout) :: this

    write (6, '(A)') 'No system specific setup provided.'
  end subroutine setup

  function get_hamiltonian_elem(this, det_pq, det_rs) result(H)
    class(solver_type), intent(inout) :: this
    type(det_type), pointer, intent(inout) :: det_pq, det_rs
    real(DOUBLE) :: H

    H = 0.0_DOUBLE
    stop 'get_hamiltonian_elem function has not been overloaded.'
  end function

  subroutine get_next_dets(this)
    class(solver_type), intent(inout) :: this
    integer :: i
    integer :: n_dets_old
    type(det_type), pointer :: tmp_det
    type(wavefunction_type), pointer :: connected_dets
    type(wavefunction_type), pointer :: new_dets

    n_dets_old = this%wf%n
    call build(new_dets)
    
    do i = 1, n_dets_old
      tmp_det => this%wf%get_det(i)
      call this%find_connected_dets( &
          & tmp_det, this%eps_var / abs(this%wf%get_coef(i)), connected_dets)
      call connected_dets%sort_dets()
      call new_dets%merge_sorted_dets(connected_dets)
    end do
    call this%wf%merge_sorted_dets(new_dets)
  end subroutine get_next_dets

  subroutine find_connected_dets(this, det, eps_min, connected_dets)
    class(solver_type), intent(inout) :: this
    type(det_type), pointer, intent(inout) :: det
    real(DOUBLE), intent(in) :: eps_min
    type(wavefunction_type), pointer, intent(out) :: connected_dets

    stop 'Default find_connected_dets has not been overloaded.'
  end subroutine find_connected_dets

  subroutine diagonalize(this, lowest_eigenvalue)
    ! Only for the case where dets are represented by integers in binary format.
    class(solver_type), intent(inout) :: this
    real(DOUBLE), intent(out) :: lowest_eigenvalue
    integer :: n
    integer :: n_up
    integer :: n_connections
    integer :: i, j, j_idx
    integer, allocatable :: n_nonzero_elems(:)
    integer, allocatable :: potential_connections(:)
    real(DOUBLE) :: H_ij
    real(DOUBLE), allocatable :: lowest_eigenvalues(:)
    real(DOUBLE), allocatable :: initial_vectors(:, :)
    real(DOUBLE), allocatable :: final_vectors(:, :)
    type(linked_list_type__int), pointer :: H_indices
    type(linked_list_type__double), pointer :: H_values
    type(det_type), pointer :: det_i, det_j

    n = this%wf%n
    n_up = this%n_up
    if (n == 1) then
      lowest_eigenvalue = this%hf_energy
      return
    end if
    allocate(n_nonzero_elems(n))
    call build(H_indices)
    call build(H_values)
    n_nonzero_elems = 0
    write (6, '(A)') 'Generating sparse hamiltonian...'
    ! call generate_sparse_hamiltonian()
    call this%wf%find_potential_connections_setup()
    do i = 1, n
      det_i => this%wf%get_det(i)
      call this%wf%find_potential_connections( &
          & det_i, potential_connections, n_connections)
      do j_idx = 1, n_connections
        j = potential_connections(j_idx)
        if (j < i) cycle
        det_j => this%wf%get_det(j)
        H_ij = this%get_hamiltonian_elem(det_i, det_j)
        if (abs(H_ij) > C%EPS) then
          call H_indices%append(j)
          call H_values%append(H_ij)
          n_nonzero_elems(i) = n_nonzero_elems(i) + 1
        endif
      enddo
    enddo

    allocate(lowest_eigenvalues(1))
    allocate(initial_vectors(n, 1))
    allocate(final_vectors(n, 1))
    do i = 1, n
      initial_vectors(i, 1) = this%wf%get_coef(i)
    enddo
    write (6, '(A)') 'Performing davidson diagonalization...'
    call davidson_sparse(n, 1, final_vectors, lowest_eigenvalues, &
        & H_indices, n_nonzero_elems, H_values, initial_vectors)
    lowest_eigenvalue = lowest_eigenvalues(1)
    do i = 1, n
      call this%wf%set_coef(i, final_vectors(i, 1))
    enddo
    call delete(H_indices)
    call delete(H_values)

    contains

!=============================================================================== 
    subroutine davidson_sparse(n,n_states,final_vector,lowest_eigenvalues,matrix_indices,nelem_nonzero,matrix_values,initial_vector)
    ! Diagonally pre-conditioned Davidson
    ! A Holmes, 17 Nov 2016

    implicit none

    ! Dummy
    integer,intent(in)        :: n,n_states
    integer,intent(in)   :: nelem_nonzero(:)
    type(linked_list_type__int), pointer, intent(inout) :: matrix_indices
    type(linked_list_type__double), pointer, intent(inout) :: matrix_values
    real(DOUBLE),intent(out)      :: final_vector(:,:)
    real(DOUBLE),intent(out)      :: lowest_eigenvalues(:)
    real(DOUBLE),intent(in),optional :: initial_vector(:,:)

    ! Local
    integer                  :: i,it
    !real(DOUBLE)                 :: rannyu
    real(DOUBLE)                 :: energy_shift
    real(DOUBLE)                 :: norm,norm_inv
    real(DOUBLE),allocatable :: residual_norm(:)
    real(DOUBLE),allocatable     :: w(:,:),Hw(:,:),v(:,:),Hv(:,:)
    integer                  :: iterations
    real(DOUBLE),allocatable :: lowest_eigenvalues_prev(:)
    logical                  :: converged=.false.
    integer                  :: len_work,info
    real(DOUBLE),allocatable     :: work(:),eigenvalues(:),h_krylov(:,:),h_overwrite(:,:)
    real(DOUBLE),allocatable     :: diag_elems(:)
    integer :: j
    integer :: niter
    integer :: n_diagonalize

    ! v(:,1:n_states) are input vectors, v(:,i) for i>n_states are residuals
    ! w(:,1:n_states) are the best vectors so far

    iterations=500        ! User option
    iterations=min(n,iterations)
    allocate(lowest_eigenvalues_prev(n_states))
    allocate(v(n,n_states*iterations))
    allocate(Hv(n,n_states*iterations))
    allocate(w(n,n_states))
    allocate(Hw(n,n_states))
    allocate(residual_norm(n_states))

    if (present(initial_vector)) then
      do i=1,n_states
        norm = 1._DOUBLE/sqrt(dot_product(initial_vector(:,i),initial_vector(:,i)))
        v(:,i) = norm*initial_vector(:,i)
        if (i>1) then
          ! Orthogonalize
          do j=1,i-1
            norm=dot_product(v(:,i),v(:,j))
            v(:,i)=v(:,i)-norm*v(:,j)
          enddo
          ! Normalize
          norm=dot_product(v(:,i),v(:,i))
          norm_inv=1._DOUBLE/sqrt(norm)
          v(:,i)=v(:,i)*norm_inv
        endif
      enddo
    else
      ! Start with HF and random dets
      v(:,:)=0
      do i=1,n_states
        v(i,i)=1
      enddo
    endif

    energy_shift=0._DOUBLE
    allocate (h_krylov(n_states*iterations,n_states*iterations))
    allocate (h_overwrite(n_states*iterations,n_states*iterations))

    allocate(eigenvalues(n_states*iterations))
    len_work = 3*n_states*iterations-1
    allocate(work(len_work))

    converged=.false.

    ! w is the lowest energy vector so far
    if (n>1) then
      ! Get diagonal elements
      allocate(diag_elems(n))
      call matrix_values%begin()
      do i = 1, n
        diag_elems(i) = matrix_values%get()
        do j = 1, n_nonzero_elems(i)
          call matrix_values%next()
        enddo
      enddo

      ! First iteration:
      do i=1,n_states
        call sparse_matrix_mul(v(:, i), Hv(:, i))
      enddo

      do i=1,n_states
        lowest_eigenvalues(i) = dot_product(v(:,i),Hv(:,i))
        h_krylov(i,i) = lowest_eigenvalues(i)
        do j=i+1,n_states
          h_krylov(i,j) = dot_product(v(:,i),Hv(:,j))
          h_krylov(j,i) = h_krylov(i,j)
        enddo
      enddo

      write(6,'(''Iteration, Eigenvalues='',i3,10f16.9)') 1, lowest_eigenvalues(:)

      w(:,:) = v(:,1:n_states)
      Hw(:,:) = Hv(:,1:n_states)

      residual_norm(:) = 1._DOUBLE ! so at least one iteration is done

      niter = min(n,n_states*(iterations+1))

      n_diagonalize = 1 ! just for printout

      do it=n_states+1,niter
         ! Compute residual for state i
         i = mod(it-1,n_states)+1
         v(:,it) = (Hw(:,i) - lowest_eigenvalues(i)*w(:,i))/(lowest_eigenvalues(i) - diag_elems(:))

         do j=1,n
           if (abs(lowest_eigenvalues(i)-diag_elems(j))<1e-8_DOUBLE)  v(j,it) = -1._DOUBLE  ! Since denominator could be 0
         enddo

         ! If residual small, converged
         residual_norm(i)=dot_product(v(:,it),v(:,it))
         if (sum(residual_norm)<(1.e-12_DOUBLE))  converged=.true.

         ! Orthogonalize
         do i=1,it-1
           norm=dot_product(v(:,it),v(:,i))
           v(:,it)=v(:,it)-norm*v(:,i)
         enddo

         ! Normalize
         norm=dot_product(v(:,it),v(:,it))
         norm_inv=1._DOUBLE/sqrt(norm)
         v(:,it)=v(:,it)*norm_inv

         ! Apply H once
         call sparse_matrix_mul(v(:, it), Hv(:, it))

         ! Construct Krylov matrix and diagonalize
         do i=1,it
           h_krylov(i,it) = dot_product(v(:,i),Hv(:,it))
           h_krylov(it,i) = h_krylov(i,it)
         enddo

         ! Diagonalize with lapack routine after all new states added
         if (mod(it,n_states)==0) then

           len_work = 3*it-1
           h_overwrite(1:it,1:it) = h_krylov(1:it,1:it)
           call dsyev('V', 'U', it, h_overwrite(1:it,1:it), it, eigenvalues, work, len_work, info)

           lowest_eigenvalues(:)=eigenvalues(1:n_states)

           do i=1,n_states
             w(:,i)=matmul(v(:,1:it),h_overwrite(1:it,i))
             Hw(:,i)=matmul(Hv(:,1:it),h_overwrite(1:it,i))
           enddo

           if (it.gt.1 .and. maxval(abs(lowest_eigenvalues(:)-lowest_eigenvalues_prev(:))) < 1.0e-10_DOUBLE) then
               converged=.true.
               exit
           else
               lowest_eigenvalues_prev(:)=lowest_eigenvalues(:)
               n_diagonalize = n_diagonalize + 1
               write(6,'(''Iteration, Eigenvalues='',i3,10f16.9)') n_diagonalize, lowest_eigenvalues(:)
               call flush(6)
           endif
           if (converged)  exit

         endif ! Diagonalize

      enddo ! it

      write(6,'(''davidson_sparse: n, Lowest eigenvalue ='',i10, 10f17.10)') n, lowest_eigenvalues(:)
      call flush(6)

      if (allocated(eigenvalues)) deallocate(eigenvalues)
      if (allocated(h_krylov))     deallocate(h_krylov)
      if (allocated(work))        deallocate(work)

    else
      stop 'Shall not reach here.'
    endif

    final_vector(:,:)=w(:,:)

    end subroutine davidson_sparse
!=============================================================================== 

    subroutine sparse_matrix_mul(vec, res)
      real(DOUBLE), intent(in) :: vec(:)
      real(DOUBLE), intent(out) :: res(:)
      integer :: i, j, m
      real(DOUBLE) :: x

      call H_indices%begin()
      call H_values%begin()
      do i = 1, this%wf%n
        x = res(i)
        do j = 1, n_nonzero_elems(i)
          m = H_indices%get()
          x = x + H_values%get() * vec(m)
          if (m /= i) then
            res(m) = res(m) + H_values%get() * vec(i)
          end if
          call H_indices%next()
          call H_values%next()
        end do
        res(i) = x
      end do
    end subroutine sparse_matrix_mul

  end subroutine diagonalize

  subroutine pt(this)
    class(solver_type), intent(inout) :: this

    this%pt_det_energy = 0.0_DOUBLE
    this%pt_st_energy = 0.0_DOUBLE
    this%pt_st_uncert = 0.0_DOUBLE

    call this%pt_det()
  end subroutine pt

  subroutine pt_det(this)
    class(solver_type), intent(inout) :: this
    logical :: is_added
    integer :: i, j, a
    integer :: j_idx
    integer :: n_connections
    real(DOUBLE) :: E_a
    real(DOUBLE) :: H_aj
    real(DOUBLE) :: eps_pt
    real(DOUBLE) :: pt_energy
    real(DOUBLE) :: var_energy
    real(DOUBLE) :: sum_a
    real(DOUBLE) :: term
    real(DOUBLE) :: coef_max
    type(det_type), pointer :: det_i, det_j, det_a
    type(wavefunction_type), pointer :: connected_dets_i
    type(wavefunction_type), pointer :: connected_dets_a
    type(lru_cache_type), pointer :: var_dets
    
    integer, allocatable :: potential_connections(:)
    integer :: cnt
    integer :: cnt_tot

    pt_energy = 0.0_DOUBLE
    var_energy = this%var_energy
    eps_pt = this%eps_pt
    call this%wf%sort_dets_by_coefs()
    call this%wf%find_potential_connections_setup()
    cnt = 0
    cnt_tot = 0
    call build(var_dets, this%wf%n * 5)
    do i = 1, this%wf%n
      det_i => this%wf%get_det(i)
      call var_dets%cache(det_i, i)
    enddo
    coef_max = 0.0_DOUBLE
    do i = 1, this%wf%n
      coef_max = max(coef_max, abs(this%wf%get_coef(i)))
    enddo
    do i = 1, this%wf%n
      det_i => this%wf%get_det(i)
      call this%find_connected_dets( &
          & det_i, eps_pt / abs(this%wf%get_coef(i)), connected_dets_i)
      do a = 1, connected_dets_i%n
        det_a => connected_dets_i%get_det(a)
        if (var_dets%has(det_a)) cycle
        if (this%wf%ab_find%lru%has(det_a)) then
          cycle
        endif
        call this%wf%ab_find%lru%cache(det_a)
        sum_a = 0.0_DOUBLE
        is_added = .false.
        cnt_tot = cnt_tot + 1
        call this%find_connected_dets(det_a, eps_pt / coef_max, connected_dets_a)
        do j_idx = 1, connected_dets_a%n
          det_j => connected_dets_a%get_det(j)
          
        enddo
        ! call this%wf%find_potential_connections( &
        !     & det_a, potential_connections, n_connections)
        ! do j_idx = 1, n_connections
        !   j = potential_connections(j_idx)
        !   det_j => this%wf%get_det(j)
        !   if (det_a == det_j) then
        !     cnt = cnt + 1
        !     is_added = .true.
        !     print *, cnt, cnt_tot, 'det_a = det_j'
        !     exit
        !   endif
        !   H_aj = this%get_hamiltonian_elem(det_a, det_j)
        !   if (abs(H_aj) < C%EPS) then
        !     cycle
        !   endif
        !   term = H_aj * this%wf%get_coef(j)
        !   if (abs(term) < eps_pt) then
        !     cycle
        !   else
        !     if (j < i) then
        !       cnt = cnt + 1
        !       is_added = .true.
        !       print *, cnt, cnt_tot, 'det_a = det_j'
        !       exit
        !     else
        !       sum_a = sum_a + term
        !     endif
        !   endif
        ! enddo
        if (is_added) then
          cycle
        endif
        E_a = this%get_hamiltonian_elem(det_a, det_a)
        pt_energy = pt_energy + sum_a**2 / (var_energy - E_a)
        
      enddo
      call delete(connected_dets_i)
    enddo
    this%pt_det_energy = pt_energy
  end subroutine pt_det

  subroutine summary(this)
    class(solver_type), intent(inout) :: this
    real(DOUBLE) :: total_energy
    real(DOUBLE) :: correlation_energy

    total_energy = this%var_energy + this%pt_det_energy + this%pt_st_energy
    correlation_energy = total_energy - this%hf_energy

    write (6, '(A, F0.10)') 'Variational Energy: ', this%var_energy
    write (6, '(A, G0.10)') 'Number of variational dets: ', this%wf%n
    write (6, '(A, F0.10)') 'Deterministic PT Correction: ', this%pt_det_energy
    write (6, '(A, F0.10, A, F0.10)') &
        & 'Stochastic PT Correction: ', this%pt_st_energy, &
        & ' +- ', this%pt_st_uncert
    write (6, '(A, F0.10, A, F0.10)') &
        & 'Final Energy: ', total_energy, ' +- ', this%pt_st_uncert
    write (6, '(A, F0.10, A, F0.10)') &
        & 'Correlation Energy: ', correlation_energy , ' +- ', this%pt_st_uncert
  end subroutine summary

end module solver_module
