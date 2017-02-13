module solver_module

  use det_module
  use dets_module
  use types_module

  implicit none

  private

  integer, parameter :: MAX_VAR_ITERATION = 20

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
    type(optional_int), public :: n_samples
    real(DOUBLE), public :: eps_var
    real(DOUBLE), public :: eps_pt
    real(DOUBLE), public :: tol
    real(DOUBLE), public :: HF_energy
    real(DOUBLE), public :: var_energy
    real(DOUBLE), public :: pt_det_energy
    real(DOUBLE), public :: pt_st_energy
    real(DOUBLE), public :: pt_st_uncert
    real(DOUBLE), public :: max_abs_H
    type(optional_double), public :: eps_pt_det
    real(DOUBLE), allocatable, public :: coefs(:)
    type(dets_type), public :: dets
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

    ! Initialized to -1 for checking availablity in namelist.
    real(DOUBLE) :: eps_pt_det = -1 ! eps for the deterministic part of PT.
    real(DOUBLE) :: n_samples = -1 ! For stochastic PT. Real to allow exponential format.

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

    integer :: iter
    integer :: n_dets_prev, n_dets_cur
    real(DOUBLE) :: energy_prev, energy_cur

    call this%read_configs(config_file_unit)

    write (6, '(A)') '[SETUP]'
    call this%setup()
    write (6, '()')

    write (6, '(A)') '[VARIATION]'
    write (6, '(A, G0.10)') 'eps_var: ', this%eps_var

    energy_prev = 0
    n_dets_prev = 1

    do iter = 1, MAX_VAR_ITERATION
      write (6, '(A, I0)') 'Iteration #', iter
      call this%get_next_dets()
      n_dets_cur = this%dets%n
      if (n_dets_cur < 1.01 * n_dets_prev) then
        this%dets%n = n_dets_prev
        write (6, '(A)') 'Number of dets change within 1%. Variation finished.'
        exit
      end if
      n_dets_prev = n_dets_cur
      energy_cur = this%diagonalize()

      write (6, '(A, F0.10)') 'Energy: ', energy_cur
      write (6, '(A, G0.10)') 'Number of dets: ', this%dets%n
      write (6, '()')
      call flush(6)
      if (abs(energy_prev - energy_cur) < max(1.0e-6, this%tol * 0.1)) then
        this%var_energy = energy_cur
        write (6, '(A)') 'Energy change within tolerance. Variation finished.'
        exit
      end if
      energy_prev = energy_cur
    end do
    write (6, '()')

    write (6, '(A)') '[PERTURBATION]'
    write (6, '(A, G0.10)') 'eps_pt: ', this%eps_pt
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
    type(det_type), intent(in) :: det_pq, det_rs
    real(DOUBLE) :: H
    stop 'get_hamiltonian_elem function has not been overloaded.'
  end function

  subroutine get_next_dets(this)
    class(solver_type), intent(inout) :: this
    integer :: i, j
    type(dets_type) :: connected_dets
    type(dets_type) :: new_dets
    type(det_type) :: tmp_det

    do i = 1, this%dets%n
      call this%find_connected_dets(this%dets%get_det(i), connected_dets, &
          & this%eps_var / abs(this%coefs(i)))
      call connected_dets%sort()
      call connected_dets%filter_by_sorted_dets(this%dets)
      call new_dets%merge_sorted_dets(connected_dets)
    end do
    call this%dets%merge_sorted_dets(new_dets)
  end subroutine get_next_dets

  subroutine find_connected_dets(this, det, connected_dets, eps_min)
    class(solver_type), intent(inout) :: this
    type(det_type), intent(in) :: det
    type(dets_type), intent(out) :: connected_dets
    real(DOUBLE), intent(in) :: eps_min

    stop 'Default find_connected_dets has not been overloaded.'
  end subroutine find_connected_dets

  function diagonalize(this) result(lowest_eigenvalue)
    class(solver_type), intent(inout) :: this
    real(DOUBLE) :: lowest_eigenvalue
    
    call generate_sparse_hamiltonian()
    call davidson_sparse()

    lowest_eigenvalue = 0

    contains

    subroutine generate_sparse_hamiltonian()
      ! Only for the case where dets are represented by integers in binary format.  
       
    end subroutine generate_sparse_hamiltonian

    subroutine davidson_sparse()
    end subroutine davidson_sparse
  end function diagonalize

  subroutine pt(this)
    class(solver_type), intent(inout) :: this

    this%pt_det_energy = 0.0_DOUBLE
    this%pt_st_energy = 0.0_DOUBLE
    this%pt_st_uncert = 0.0_DOUBLE
  end subroutine pt

  subroutine summary(this)
    class(solver_type), intent(inout) :: this
    real(DOUBLE) :: total_energy

    total_energy = this%var_energy + this%pt_det_energy + this%pt_st_energy

    write (6, '(A, F0.10)') 'Variational Energy: ', this%var_energy
    write (6, '(A, G0.10)') 'Number of variational dets: ', this%dets%n
    write (6, '(A, F0.10)') 'Deterministic PT Energy: ', this%pt_det_energy
    write (6, '(A, F0.10, A, F0.10)') &
        & 'Stochastic PT Energy: ', this%pt_st_energy, ' +- ', this%pt_st_uncert
    write (6, '(A, F0.10, A, F0.10)') &
        & 'Total Energy: ', total_energy, ' +- ', this%pt_st_uncert
  end subroutine summary

end module solver_module
