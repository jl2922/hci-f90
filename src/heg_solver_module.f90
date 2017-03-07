module heg_solver_module

  use constants_module
  use det_module
  use solver_module
  use spin_det_module
  use types_module
  use utilities_module
  use wavefunction_module

  implicit none
  
  private

  public :: heg_solver_instance

  type heg_data_type
    integer :: n_diff
    integer :: n_diff_offset
    integer :: n_max ! Max integer n < r_cutoff.
    integer, allocatable :: k_vectors(:, :)
    integer, allocatable :: orbital_lut(:, :, :)
    real(DOUBLE) :: cell_length
    real(DOUBLE) :: k_unit
    real(DOUBLE) :: r_cutoff
    real(DOUBLE) :: r_s
    type(INT_3_DOUBLE), allocatable :: same_spin(:, :, :, :)
    type(INT_3_DOUBLE), allocatable :: opposite_spin(:)
  end type heg_data_type

  type, extends(solver_type) :: heg_solver_type
    type(heg_data_type) :: heg
    contains
      procedure :: find_orb_id
      procedure :: find_connected_dets
      procedure :: generate_k_vectors
      procedure :: get_hamiltonian_elem
      procedure :: generate_orbital_lut
      procedure :: generate_hci_queue
      procedure :: hamiltonian_abs_by_pqrs
      procedure :: read_system_configs
      procedure :: setup
  end type heg_solver_type

  type(heg_solver_type) :: heg_solver_instance

  type(det_type), target :: tmp_det_for_get_hamiltonian_elem
  type(spin_det_type), target :: tmp_up_for_get_hamiltonian_elem
  type(spin_det_type), target :: tmp_dn_for_get_hamiltonian_elem

  contains

  subroutine read_system_configs(this, config_file_unit)
    class(heg_solver_type), intent(inout) :: this
    integer, intent(in) :: config_file_unit
    real(DOUBLE) :: r_cutoff = -1
    real(DOUBLE) :: r_s = -1
    integer :: io_err
    namelist /heg/ r_cutoff, r_s

    rewind(config_file_unit)
    read(unit=config_file_unit, nml=heg, iostat=io_err)
    if (io_err > 0) stop 'Cannot read HEG configs.'
    if (r_cutoff < 0 .or. r_s < 0) then
      stop 'r_cutoff, r_s must be provided for heg system.'
    end if
    this%heg%r_cutoff = r_cutoff
    this%heg%r_s = r_s
    write (6, '(A, F0.10)') 'heg%r_cutoff: ', r_cutoff
    write (6, '(A, F0.10)') 'heg%r_s: ', r_s
    write (6, '()')
  end subroutine read_system_configs

  subroutine setup(this)
    class(heg_solver_type), intent(inout) :: this
    integer :: det_size
    integer :: i
    integer :: n_orb
    integer :: n_up, n_dn
    real(DOUBLE) :: cell_length
    real(DOUBLE) :: density
    real(DOUBLE) :: HF_energy
    real(DOUBLE) :: k_unit
    type(det_type), pointer :: tmp_det

    ! Basic quantities.
    density = 3.0_DOUBLE / (4.0_DOUBLE * C%PI * this%heg%r_s**3)
    cell_length = (this%n_elec / density)**(1.0_DOUBLE / 3)
    k_unit = 2 * C%PI / cell_length
    write (6, '(A, G0.10)') 'Cell length: ', cell_length
    write (6, '(A, G0.10)') 'K unit: ', k_unit
    this%heg%cell_length = cell_length
    this%heg%k_unit = k_unit
    call this%generate_k_vectors()
    n_orb = size(this%heg%k_vectors, 2)
    write (6, '(A, I0)') 'Number of spatial orbitals: ', n_orb
    write (6, '(A, I0)') 'Number of spin orbitals: ', n_orb * 2
    this%n_orb = n_orb
    n_up = this%n_up
    n_dn = this%n_dn
    this%max_connected_dets = &
        & 1 + n_up * n_dn * (n_orb - n_dn) * (n_orb - n_up) &
        & + util%combinatorics%C(n_up, 2) &
        & + util%combinatorics%C(n_orb - n_up, 2) &
        & + util%combinatorics%C(n_dn, 2) &
        & + util%combinatorics%C(n_orb - n_dn, 2)

    ! Setup HF determinant.
    this%wf => new_wavefunction(1)
    tmp_det => new_det(n_orb)
    do i = 1, this%n_up
      call tmp_det%up%set_orbital(i, .true.)
    end do
    do i = 1, this%n_dn
      call tmp_det%dn%set_orbital(i, .true.)
    end do
    call this%wf%append_det(tmp_det)
    call this%wf%set_coef(1, 1.0_DOUBLE)
    HF_energy = this%get_hamiltonian_elem(tmp_det, tmp_det)
    call delete(tmp_det)
    write (6, '(A, F0.10)') 'HF energy: ', HF_energy
    this%HF_energy = HF_energy

    call this%generate_orbital_lut()
    call this%generate_hci_queue()
  end subroutine setup

  subroutine generate_k_vectors(this)
    class(heg_solver_type), intent(inout) :: this

    integer :: n_max
    integer :: n_orb
    integer :: i, j, k
    integer :: length2
    integer, allocatable :: temp_k_vectors(:, :)
    integer, allocatable :: order(:)
    real(DOUBLE), allocatable :: temp_k_length2(:)

    n_max = int(this%heg%r_cutoff + C%EPS)
    this%heg%n_max = n_max
    allocate(temp_k_vectors(3, (2 * n_max + 1)**3))
    allocate(temp_k_length2((2 * n_max + 1)**3))
    n_orb = 0
    do i = -n_max, n_max
      do j = -n_max, n_max
        do k = -n_max, n_max
          length2 = i**2 + j**2 + k**2
          if (length2 > this%heg%r_cutoff**2) then
            cycle
          end if
          n_orb = n_orb + 1
          temp_k_vectors(:, n_orb) = (/i, j, k/)
          temp_k_length2(n_orb) = length2
        end do
      end do
    end do
    allocate(order(n_orb))
    
    call util%sort%arg_sort(temp_k_length2, order, n_orb)
    allocate(this%heg%k_vectors(3, n_orb))
    do i = 1, n_orb
      this%heg%k_vectors(:, i) = temp_k_vectors(:, order(i))
      write (6, '(A, I0, 3F15.10)') &
          & 'K-point #', i, this%heg%k_vectors(:, i) * this%heg%k_unit
    end do
  end subroutine generate_k_vectors

  function get_hamiltonian_elem(this, det_pq, det_rs) result(H)
    class(heg_solver_type), intent(inout) :: this
    type(det_type), pointer, intent(inout) :: det_pq, det_rs
    real(DOUBLE) :: H
    logical :: k_p_set, k_q_set, k_s_set
    integer :: gamma_exp
    integer :: i, j
    integer :: p, q
    integer :: k_change(3)
    integer :: n_occ_up, n_occ_dn
    integer :: n_eor_up, n_eor_dn
    integer :: orb_p, orb_q, orb_s
    integer :: orb_id
    integer, allocatable :: occ_pq_up(:), occ_pq_dn(:), occ_rs_up(:), occ_rs_dn(:)
    integer, allocatable :: eor_up_set_bits(:), eor_dn_set_bits(:)
    type(det_type), pointer :: det_eor
    real(DOUBLE) :: k_unit
    real(DOUBLE) :: one_over_pi_l

    k_unit = this%heg%k_unit
    one_over_pi_l = 1.0_DOUBLE / (C%PI * this%heg%cell_length)

    if (det_pq == det_rs) then
      ! Diagonal elements.
      H = 0.0_DOUBLE
      n_occ_up = det_pq%up%get_n_elec()
      n_occ_dn = det_pq%dn%get_n_elec()
      call det_pq%up%get_elec_orbitals(occ_pq_up, n_occ_up)
      call det_pq%dn%get_elec_orbitals(occ_pq_dn, n_occ_dn)

      ! One electron operator.
      do i = 1, n_occ_up
        p = occ_pq_up(i)
        H = H + sum((this%heg%k_vectors(:, p) * k_unit)**2) * 0.5_DOUBLE
      end do
      do i = 1, n_occ_dn
        p = occ_pq_dn(i)
        H = H + sum((this%heg%k_vectors(:, p) * k_unit)**2) * 0.5_DOUBLE
      end do

      ! Two electron operator.
      do i = 1, n_occ_up
        p = occ_pq_up(i)
        do j = i + 1, n_occ_up
          q = occ_pq_up(j)
          H = H - one_over_pi_l / sum((this%heg%k_vectors(:, p) - this%heg%k_vectors(:, q))**2)
        end do
      end do
      do i = 1, n_occ_dn
        p = occ_pq_dn(i)
        do j = i + 1, n_occ_dn
          q = occ_pq_dn(j)
          H = H - one_over_pi_l / sum((this%heg%k_vectors(:, p) - this%heg%k_vectors(:, q))**2)
        end do
      end do
    else
      ! Off-diagonal elements.
      det_eor => tmp_det_for_get_hamiltonian_elem
      det_eor%up => tmp_up_for_get_hamiltonian_elem
      det_eor%dn => tmp_dn_for_get_hamiltonian_elem
      call det_eor%from_eor(det_pq, det_rs)
      n_eor_up = det_eor%up%get_n_elec()
      n_eor_dn = det_eor%dn%get_n_elec()
      if (n_eor_up + n_eor_dn /= 4) then
        H = 0.0_DOUBLE ! Not double excitation.
        return
      end if
      call det_eor%up%get_elec_orbitals(eor_up_set_bits, n_eor_up)
      call det_eor%dn%get_elec_orbitals(eor_dn_set_bits, n_eor_dn)
      k_change = 0.0_DOUBLE
      k_p_set = .false.
      k_q_set = .false.
      k_s_set = .false.
      do i = 1, n_eor_up
        orb_id = eor_up_set_bits(i)
        if (det_pq%up%get_orbital(orb_id)) then
          k_change = k_change - this%heg%k_vectors(:, orb_id)
          if (.not. k_p_set) then
            orb_p = orb_id
            k_p_set = .true.
          end if
        else
          k_change = k_change + this%heg%k_vectors(:, orb_id)
          if (.not. k_q_set) then
            orb_q = orb_id
            k_q_set = .true.
          else
            orb_s = orb_id
            k_s_set = .true.
          end if
        end if
      end do
      do i = 1, n_eor_dn
        orb_id = eor_dn_set_bits(i)
        if (det_pq%dn%get_orbital(orb_id)) then
          k_change = k_change - this%heg%k_vectors(:, orb_id)
          if (.not. k_p_set) then
            orb_p = orb_id
            k_p_set = .true.
          end if
        else
          k_change = k_change + this%heg%k_vectors(:, orb_id)
          if (.not. k_q_set) then
            orb_q = orb_id
            k_q_set = .true.
          else
            orb_s = orb_id
            k_s_set = .true.
          end if
        end if
      end do
      if (sum(k_change**2) > 0) then
        H = 0.0_DOUBLE ! Momentum not conserved.
        return
      end if

      H = one_over_pi_l / sum( &
          & (this%heg%k_vectors(:, orb_p) - this%heg%k_vectors(:, orb_q))**2)
      if (n_eor_up /= 2) then
        ! Either two up or two dn excitations.
        H = H - one_over_pi_l / sum( &
            & (this%heg%k_vectors(:, orb_p) - this%heg%k_vectors(:, orb_s))**2)
      end if
      call det_pq%up%get_elec_orbitals(occ_pq_up)
      call det_pq%dn%get_elec_orbitals(occ_pq_dn)
      call det_rs%up%get_elec_orbitals(occ_rs_up)
      call det_rs%dn%get_elec_orbitals(occ_rs_dn)
      gamma_exp = &
          & get_gamma_exp(det_pq%up, occ_pq_up, eor_up_set_bits, n_eor_up) &
          & + get_gamma_exp(det_pq%dn, occ_pq_dn, eor_dn_set_bits, n_eor_dn) &
          & + get_gamma_exp(det_rs%up, occ_rs_up, eor_up_set_bits, n_eor_up) &
          & + get_gamma_exp(det_rs%dn, occ_rs_dn, eor_dn_set_bits, n_eor_dn)
      if (iand(gamma_exp, 1) == 1) then
        H = -H
      end if
    end if
  end function get_hamiltonian_elem

  subroutine generate_orbital_lut(this)
    class(heg_solver_type), intent(inout) :: this
    integer :: n_max
    integer :: i
    integer :: k_idx(3)

    n_max = this%heg%n_max
    allocate(this%heg%orbital_lut(2 * n_max + 1, 2 * n_max + 1, 2 * n_max + 1))
    this%heg%orbital_lut = -1
    do i = 1, this%n_orb
      k_idx = this%heg%k_vectors(:, i) + n_max + 1
      this%heg%orbital_lut(k_idx(1), k_idx(2), k_idx(3)) = i
    end do
  end subroutine generate_orbital_lut

  subroutine generate_hci_queue(this)
    class(heg_solver_type), intent(inout) :: this
    logical, allocatable :: same_spin_processed(:, :, :, :, :, :)
    logical, allocatable :: opposite_spin_processed(:, :, :)
    integer :: n_diff
    integer :: n_diff_offset
    integer :: n_max
    integer :: p, q, r, s
    integer :: pq_cnt
    integer :: opposite_spin_cnt
    integer :: max_n_rs_pairs
    integer :: i
    integer, dimension(3) :: diff_pq, diff_pr
    integer, dimension(3) :: diff_pq_idx, diff_pr_idx
    integer, allocatable :: order(:)
    integer, allocatable :: same_spin_pq_cnt(:, :, :)
    real(DOUBLE) :: H_pqrs
    real(DOUBLE) :: max_abs_H
    real(DOUBLE), allocatable :: abs_H_copy(:)
    type(INT_3_DOUBLE), allocatable :: hci_queue_copy(:)

    n_max = this%heg%n_max
    n_diff = 4 * n_max + 1
    n_diff_offset = 2 * n_max + 1
    max_n_rs_pairs = 0
    max_abs_H = 0.0_DOUBLE
    this%heg%n_diff = n_diff
    this%heg%n_diff_offset = n_diff_offset
    allocate(this%heg%same_spin(n_diff, n_diff, n_diff, n_diff**3))
    allocate(this%heg%opposite_spin(n_diff**3))
    allocate(same_spin_pq_cnt(n_diff, n_diff, n_diff))
    allocate(same_spin_processed( &
        & n_diff, n_diff, n_diff, n_diff, n_diff, n_diff))
    allocate(opposite_spin_processed(n_diff, n_diff, n_diff))
    allocate(order(n_diff**3))
    allocate(abs_H_copy(n_diff**3))
    allocate(hci_queue_copy(n_diff**3))
    same_spin_pq_cnt = 0
    same_spin_processed = .false.
    opposite_spin_processed = .false.

    ! Same spin.
    do p = 1, this%n_orb
      do q = p + 1, this%n_orb
        diff_pq = this%heg%k_vectors(:, q) - this%heg%k_vectors(:, p)
        diff_pq_idx = diff_pq + n_diff_offset
        pq_cnt = 0
        do r = 1, this%n_orb
          s = this%find_orb_id( &
              & this%heg%k_vectors(:, p) + this%heg%k_vectors(:, q) &
              & - this%heg%k_vectors(:, r))
          if (s < r) cycle
          diff_pr = this%heg%k_vectors(:, r) - this%heg%k_vectors(:, p)
          diff_pr_idx = diff_pr + n_diff_offset
          if (same_spin_processed( &
              & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3), &
              & diff_pr_idx(1), diff_pr_idx(2), diff_pr_idx(3))) then
            cycle
          end if
          same_spin_processed( &
              & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3), &
              & diff_pr_idx(1), diff_pr_idx(2), diff_pr_idx(3)) = .true.
          H_pqrs = this%hamiltonian_abs_by_pqrs(p, q, r, s)
          if (H_pqrs > C%EPS) then
            pq_cnt = same_spin_pq_cnt( &
                & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3)) + 1
            this%heg%same_spin( &
                & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3), pq_cnt)%int_3 &
                & = diff_pr
            this%heg%same_spin( &
                & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3), pq_cnt)%double &
                & = H_pqrs
            same_spin_pq_cnt( &
                & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3)) = pq_cnt
          end if
        end do
        if (pq_cnt > 0) then
          abs_H_copy(1: pq_cnt) = this%heg%same_spin( &
              & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3), 1: pq_cnt)%double
          hci_queue_copy(1: pq_cnt) = this%heg%same_spin( &
              & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3), 1: pq_cnt)
          call util%sort%arg_sort(-abs_H_copy, order, pq_cnt)
          do i = 1, pq_cnt
            this%heg%same_spin( &
                & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3), i) &
                & = hci_queue_copy(order(i))
          end do
          max_abs_H = max(max_abs_H, this%heg%same_spin( &
              & diff_pq_idx(2), diff_pq_idx(2), diff_pq_idx(3), 1)%double)
          max_n_rs_pairs = max(max_n_rs_pairs, pq_cnt)
        end if
      end do
    end do
    write (6, '(A)') 'Same spin hci queue generated.'

    ! Opposite spin.
    opposite_spin_cnt = 0
    do p = 1, this%n_orb
      do q = p + this%n_orb, this%n_orb * 2
        do r = 1, this%n_orb
          s = this%find_orb_id(this%heg%k_vectors(:, p) &
              & + this%heg%k_vectors(:, q - this%n_orb) &
              & - this%heg%k_vectors(:, r))
          if (s < 0) cycle
          s = s + this%n_orb
          diff_pr = this%heg%k_vectors(:, r) - this%heg%k_vectors(:, p)
          diff_pr_idx = diff_pr + n_diff_offset
          if (opposite_spin_processed( &
              & diff_pr_idx(1), diff_pr_idx(2), diff_pr_idx(3))) then
            cycle
          end if
          opposite_spin_processed( &
              & diff_pr_idx(1), diff_pr_idx(2), diff_pr_idx(3)) = .true.
          H_pqrs = this%hamiltonian_abs_by_pqrs(p, q, r, s)
          if (H_pqrs > C%EPS) then
            opposite_spin_cnt = opposite_spin_cnt + 1
            this%heg%opposite_spin(opposite_spin_cnt)%int_3 = diff_pr
            this%heg%opposite_spin(opposite_spin_cnt)%double = H_pqrs
            abs_H_copy(opposite_spin_cnt) = H_pqrs
          end if
        end do
      end do
    end do
    if (opposite_spin_cnt > 0) then
      call util%sort%arg_sort(-abs_H_copy, order, opposite_spin_cnt)
      hci_queue_copy(1: opposite_spin_cnt) = &
          & this%heg%opposite_spin(1: opposite_spin_cnt)
      do i = 1, opposite_spin_cnt
        this%heg%opposite_spin(i) = hci_queue_copy(order(i))
      end do
      max_abs_H = max(max_abs_H, this%heg%opposite_spin(1)%double)
    else
      this%heg%opposite_spin(1)%double = 0.0_DOUBLE
    end if
    max_n_rs_pairs = max(max_n_rs_pairs, opposite_spin_cnt)
    write (6, '(A)') 'Opposite spin hci queue generated.'
    
    this%max_abs_H = max_abs_H
    this%max_n_rs_pairs = max_n_rs_pairs
  end subroutine generate_hci_queue

  function hamiltonian_abs_by_pqrs(this, p, q, r, s) result(abs_H)
    class(heg_solver_type), intent(inout) :: this
    integer, intent(in) :: p, q, r, s
    real(DOUBLE) :: abs_H
    type(det_type), pointer :: zero_det
    type(det_type), pointer :: det_pq, det_rs
    integer :: n_orb

    if (p == q .or. r == s .or. p == r .or. q == s .or. p == s .or. q == r) then
      abs_H = 0.0_DOUBLE
      return
    end if
    n_orb = this%n_orb
    zero_det => new_det(n_orb)
    det_pq => new_det(zero_det)
    det_rs => new_det(zero_det)
    call delete(zero_det)
    if (p <= n_orb .and. q <= n_orb) then
      call det_pq%up%set_orbital(p, .true.)
      call det_pq%up%set_orbital(q, .true.)
      call det_rs%up%set_orbital(r, .true.)
      call det_rs%up%set_orbital(s, .true.)
    else if (p > n_orb .and. q > n_orb) then
      call det_pq%dn%set_orbital(p - n_orb, .true.)
      call det_pq%dn%set_orbital(q - n_orb, .true.)
      call det_rs%dn%set_orbital(r - n_orb, .true.)
      call det_rs%dn%set_orbital(s - n_orb, .true.)
    else
      call det_pq%up%set_orbital(p, .true.)
      call det_pq%dn%set_orbital(q - n_orb, .true.)
      call det_rs%up%set_orbital(r, .true.)
      call det_rs%dn%set_orbital(s - n_orb, .true.)
    end if
    abs_H = abs(this%get_hamiltonian_elem(det_pq, det_rs))
    call delete(det_pq)
    call delete(det_rs)
  end function hamiltonian_abs_by_pqrs

  function find_orb_id(this, k_vector) result(orb_id)
    class(heg_solver_type), intent(inout) :: this
    integer, intent(in) :: k_vector(3)
    integer :: orb_id
    integer :: k_idx(3)
    integer :: i
    integer :: n_max

    n_max = this%heg%n_max
    do i = 1, 3
      if (k_vector(i) < -n_max .or. k_vector(i) > n_max) then
        orb_id = -1
        return
      end if
    end do
    k_idx = k_vector + n_max + 1
    orb_id = this%heg%orbital_lut(k_idx(1), k_idx(2), k_idx(3))
  end function find_orb_id

  function get_gamma_exp(spin_det, occ, eor_set_bits, n_eor_set_bits) result(res)
    type(spin_det_type), pointer, intent(inout) :: spin_det
    integer, intent(in) :: occ(:)
    integer, intent(in) :: eor_set_bits(:)
    integer, intent(in) :: n_eor_set_bits
    integer :: res
    logical :: is_occupied
    integer :: i
    integer :: ptr
    integer :: orb_id

    res = 0
    ptr = 0
    do i = 1, n_eor_set_bits
      orb_id = eor_set_bits(i)
      is_occupied = spin_det%get_orbital(orb_id)
      if (.not. is_occupied) cycle
      do while (occ(ptr + 1) < orb_id)
        ptr = ptr + 1
      end do
      res = res + ptr
    end do
  end function get_gamma_exp

  subroutine find_connected_dets(this, det, eps_min, connected_dets)
    class(heg_solver_type), intent(inout) :: this
    type(det_type), pointer, intent(inout) :: det
    real(DOUBLE), intent(in) :: eps_min
    integer :: n_elec
    integer :: n_up, n_dn
    integer :: n_orb
    integer :: n_pq_pairs
    integer :: n_rs_pairs
    integer :: p, q, r, s
    integer :: i, j
    type(INT_PAIR), allocatable :: pq_pairs(:)
    type(INT_PAIR), allocatable :: rs_pairs(:)
    type(det_type), pointer :: tmp_det
    type(wavefunction_type), pointer, intent(out) :: connected_dets

    if (this%max_abs_H < eps_min) then
      connected_dets => new_wavefunction(1)
      return
    end if
    n_elec = this%n_elec
    n_up = this%n_up
    n_dn = this%n_dn
    n_orb = this%n_orb
    tmp_det => new_det(n_orb)
    connected_dets => new_wavefunction(this%max_connected_dets)
    call get_pq_pairs(pq_pairs, n_pq_pairs)
    allocate(rs_pairs(this%max_n_rs_pairs + 1))
    do i = 1, n_pq_pairs
      p = pq_pairs(i)%i1
      q = pq_pairs(i)%i2
      call get_rs_pairs(p, q, rs_pairs, n_rs_pairs)
      do j = 1, n_rs_pairs
        r = rs_pairs(j)%i1
        s = rs_pairs(j)%i2
        call process_pqrs(p, q, r, s, connected_dets)
      end do
    end do

    contains

      subroutine get_pq_pairs(pq_pairs, n_pq_pairs)
        type(INT_PAIR), allocatable, intent(out) :: pq_pairs(:)
        integer, intent(out) :: n_pq_pairs
        integer :: i
        integer :: p, q
        integer :: n_occ_up, n_occ_dn
        integer, allocatable :: occ_up(:), occ_dn(:)

        n_pq_pairs = n_up * (n_up - 1) / 2 + n_dn * (n_dn - 1) / 2 + n_up * n_dn
        allocate(pq_pairs(n_pq_pairs))
        call det%up%get_elec_orbitals(occ_up)
        call det%dn%get_elec_orbitals(occ_dn)
        i = 0
        do p = 1, n_up
          do q = p + 1, n_up
            i = i + 1
            pq_pairs(i)%i1 = occ_up(p)
            pq_pairs(i)%i2 = occ_up(q)
          end do
        end do
        do p = 1, n_dn
          do q = p + 1, n_dn
            i = i + 1
            pq_pairs(i)%i1 = occ_dn(p) + n_orb
            pq_pairs(i)%i2 = occ_dn(q) + n_orb
          end do
        end do
        do p = 1, n_up
          do q = 1, n_dn
            i = i + 1
            pq_pairs(i)%i1 = occ_up(p)
            pq_pairs(i)%i2 = occ_dn(q) + n_orb
          end do
        end do
      end subroutine get_pq_pairs

      subroutine get_rs_pairs(p, q, rs_pairs, n_rs_pairs)
        integer, intent(in) :: p, q
        type(INT_PAIR), intent(out) :: rs_pairs(:)
        integer, intent(out) :: n_rs_pairs
        integer :: j
        integer :: r, s
        integer :: tmp
        integer, dimension(3) :: diff_pq, diff_pr
        integer, dimension(3) :: diff_pq_idx, diff_pr_idx
        integer :: p2, q2

        p2 = p
        q2 = q
        if (p > n_orb .and. q > n_orb) then
          p2 = p - n_orb
          q2 = q - n_orb
        else if (p <= n_orb .and. q > n_orb .and. p > q - n_orb) then
          p2 = q - n_orb
          q2 = p + n_orb
        end if

        n_rs_pairs = 0
        do j = 1, this%heg%n_diff**3
          if (p2 <= n_orb .and. q2 <= n_orb) then
            diff_pq = this%heg%k_vectors(:, q2) - this%heg%k_vectors(:, p2)
            diff_pq_idx = diff_pq + this%heg%n_diff_offset
            if (this%heg%same_spin( &
                & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3), j)%double &
                & < eps_min) then
              exit
            end if
            diff_pr = this%heg%same_spin( &
                & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3), j)%int_3
            r = this%find_orb_id(diff_pr + this%heg%k_vectors(:, p2))
            if (r < 0) cycle
            s = this%find_orb_id(this%heg%k_vectors(:, p2) &
                & + this%heg%k_vectors(:, q2) &
                & - this%heg%k_vectors(:, r))
            if (s < r) cycle
          else
            if (this%heg%opposite_spin(j)%double < eps_min) then
              exit
            end if
            diff_pr = this%heg%opposite_spin(j)%int_3
            r = this%find_orb_id(diff_pr + this%heg%k_vectors(:, p2))
            if (r < 0) cycle
            s = this%find_orb_id(this%heg%k_vectors(:, p2) &
                & + this%heg%k_vectors(:, q2 - n_orb) &
                & - this%heg%k_vectors(:, r))
            if (s < 0) cycle
            s = s + n_orb
          end if
          if (p > n_orb .and. q > n_orb) then
            r = r + n_orb
            s = s + n_orb
          else if (p <= n_orb .and. q > n_orb .and. p > q - n_orb) then
            tmp = s - n_orb
            s = r + n_orb
            r = tmp
          end if
          n_rs_pairs = n_rs_pairs + 1
          rs_pairs(n_rs_pairs)%i1 = r
          rs_pairs(n_rs_pairs)%i2 = s
        end do
      end subroutine get_rs_pairs

      subroutine process_pqrs(p, q, r, s, connected_dets)
        integer, intent(in) :: p, q, r, s
        type(wavefunction_type), pointer, intent(inout) :: connected_dets
        integer :: n_orb

        n_orb = this%n_orb
        if (r > n_orb) then
          if (det%dn%get_orbital(r - n_orb)) return
        else
          if (det%up%get_orbital(r)) return
        end if
        if (s > n_orb) then
          if (det%dn%get_orbital(s - n_orb)) return
        else
          if (det%up%get_orbital(s)) return
        end if

        tmp_det = det
        if (p <= n_orb) then
          call tmp_det%up%set_orbital(p, .false.)
        else
          call tmp_det%dn%set_orbital(p - n_orb, .false.)
        end if
        if (q <= n_orb) then
          call tmp_det%up%set_orbital(q, .false.)
        else
          call tmp_det%dn%set_orbital(q - n_orb, .false.)
        end if
        if (r <= n_orb) then
          call tmp_det%up%set_orbital(r, .true.)
        else
          call tmp_det%dn%set_orbital(r - n_orb, .true.)
        end if
        if (s <= n_orb) then
          call tmp_det%up%set_orbital(s, .true.)
        else
          call tmp_det%dn%set_orbital(s - n_orb, .true.)
        end if
        call connected_dets%append_det(tmp_det)
      end subroutine

  end subroutine find_connected_dets

end module heg_solver_module
