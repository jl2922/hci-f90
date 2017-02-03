module heg_solver_class
  
  use types_class
  use constants_class
  use config_class
  use solver_class
  use sort_class

  implicit none

  private

  public :: heg_solver

  type heg_data
    integer :: n_max ! Max integer n < r_cutoff
    integer, allocatable :: k_vectors(:, :)
    integer, allocatable :: orbital_lut(:, :, :)
    real(DOUBLE) :: cell_length
    real(DOUBLE) :: k_unit
    type(INT_3_DOUBLE), allocatable :: same_spin(:, :, :, :)
    type(INT_3_DOUBLE), allocatable :: opposite_spin(:)
  end type heg_data

  type, extends(solver) :: heg_solver
    private
    type(heg_data), public :: heg
    contains
      procedure, public :: setup
      procedure :: find_orb_id
      procedure :: generate_k_vectors
      procedure :: generate_orbital_lut
      procedure :: generate_hci_queue
      procedure :: hamiltonian
      procedure :: hamiltonian_abs_by_orbitals
  end type heg_solver

  contains

    subroutine setup(this, config_instance)
      class(heg_solver), intent(inout) :: this
      type(config), intent(in) :: config_instance

      integer :: n_orb
      integer :: det_size
      integer :: i
      real(DOUBLE) :: density
      real(DOUBLE) :: cell_length
      real(DOUBLE) :: k_unit
      real(DOUBLE) :: HF_energy

      ! Obtain cell size.
      density = 3.0_DOUBLE / (4.0_DOUBLE * PI * config_instance%heg%r_s**3)
      cell_length = (config_instance%heg%n_elec / density)**(1.0_DOUBLE / 3)
      k_unit = 2 * PI / cell_length
      write (6, '(A, G0.10)') 'Cell length: ', cell_length
      write (6, '(A, G0.10)') 'K unit: ', k_unit
      this%heg%cell_length = cell_length
      this%heg%k_unit = k_unit

      call this%generate_k_vectors(config_instance%heg%r_cutoff)

      ! Obtain determinant representation format.
      n_orb = size(this%heg%k_vectors, 2)
      det_size = ceiling(n_orb / 64.0_DOUBLE)
      write (6, '(A, I0)') 'Number of spatial orbitals: ', n_orb
      write (6, '(A, I0)') 'Number of spin orbitals: ', n_orb * 2
      write (6, '(A, I0)') 'Number of INT64 per det: ', det_size
      this%n_orb = n_orb
      this%det_size = det_size

      ! Setup HF determinant.
      this%dets%n = 1
      this%dets%det_size = det_size
      allocate(this%coef(1))
      allocate(this%dets%up(det_size, 1))
      allocate(this%dets%dn(det_size, 1))
      this%coef(1) = 1.0_DOUBLE
      this%dets%up(:, 1) = 0
      this%dets%dn(:, 1) = 0
      do i = 1, config_instance%heg%n_up
        call abset(this%dets%up(:, 1), i)
      end do
      do i = 1, config_instance%heg%n_dn
        call abset(this%dets%dn(:, 1), i)
      end do
      HF_energy = this%hamiltonian(this%dets%up(:, 1), this%dets%dn(:, 1), &
          & this%dets%up(:, 1), this%dets%dn(:, 1))
      write (6, '(A, F0.10)') 'HF energy: ', HF_energy
      this%HF_energy = HF_energy

      call this%generate_orbital_lut()
      
      call this%generate_hci_queue()
    end subroutine setup

    type(integer) function find_orb_id(this, k_vector) result(orb_id)
      class(heg_solver), intent(inout) :: this
      integer, intent(in) :: k_vector(3)

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
    
    subroutine generate_k_vectors(this, r_cutoff)
      class(heg_solver), intent(inout) :: this
      real(DOUBLE), intent(in) :: r_cutoff
      
      integer :: n_max
      integer :: n_orb
      integer :: i, j, k
      integer :: length2
      integer, allocatable :: temp_k_vectors(:, :)
      integer, allocatable :: temp_k_length2(:)
      integer, allocatable :: order(:)

      n_max = int(r_cutoff + EPS)
      this%heg%n_max = n_max

      allocate(temp_k_vectors(3, (2 * n_max + 1)**3))
      allocate(temp_k_length2((2 * n_max + 1)**3))

      n_orb = 0
      do i = -n_max, n_max
        do j = -n_max, n_max
          do k = -n_max, n_max
            length2 = i**2 + j**2 + k**2
            if (length2 > r_cutoff**2) then
              cycle
            end if
            n_orb = n_orb + 1
            temp_k_vectors(:, n_orb) = (/i, j, k/)
            temp_k_length2(n_orb) = length2 
          end do
        end do
      end do

      allocate(order(n_orb))
      call merge_sort_arg(temp_k_length2, order, n_orb) ! Sort order by magnitude.

      allocate(this%heg%k_vectors(3, n_orb))
      do i = 1, n_orb
        this%heg%k_vectors(:, i) = temp_k_vectors(:, order(i))
      end do
    end subroutine generate_k_vectors

    subroutine generate_orbital_lut(this)
      class(heg_solver), intent(inout) :: this

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
      class(heg_solver), intent(inout) :: this
      
      logical, allocatable :: same_spin_processed(:, :, :, :, :, :)
      logical, allocatable :: opposite_spin_processed(:, :, :)
      integer :: n_diff
      integer :: n_diff_offset
      integer :: n_max
      integer :: p, q, r, s
      integer :: pq_cnt
      integer :: opposite_spin_cnt
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
                & this%heg%k_vectors(:, p) + this%heg%k_vectors(:, q) - this%heg%k_vectors(:, r))
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

            H_pqrs = this%hamiltonian_abs_by_orbitals(p, q, r, s)
            
            if (H_pqrs > EPS) then
              pq_cnt = same_spin_pq_cnt( &
                  & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3)) + 1
              this%heg%same_spin( &
                  & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3), pq_cnt)%int_3 = diff_pr
              this%heg%same_spin( &
                  & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3), pq_cnt)%double = H_pqrs
              same_spin_pq_cnt( &
                  & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3)) = pq_cnt
            end if
          end do

          if (pq_cnt > 0) then
            abs_H_copy(1: pq_cnt) = &
                & this%heg%same_spin( &
                    & diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3), 1: pq_cnt)%double
            hci_queue_copy(1: pq_cnt) = &
                & this%heg%same_spin(diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3), 1: pq_cnt)
            call merge_sort_arg(-abs_H_copy, order, pq_cnt)
            do i = 1, pq_cnt
              this%heg%same_spin(diff_pq_idx(1), diff_pq_idx(2), diff_pq_idx(3), i) = &
                  & hci_queue_copy(order(i))
            end do
            max_abs_H = max(max_abs_H, this%heg%same_spin( &
                & diff_pq_idx(2), diff_pq_idx(2), diff_pq_idx(3), 1)%double)
          end if
        end do
      end do

      ! Opposite spin.
      opposite_spin_cnt = 0
      do p = 1, this%n_orb
        do q = p + this%n_orb, this%n_orb * 2
          do r = 1, this%n_orb
            s = this%find_orb_id(this%heg%k_vectors(:, p) + this%heg%k_vectors(:, q - this%n_orb) &
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

            H_pqrs = this%hamiltonian_abs_by_orbitals(p, q, r, s)

            if (H_pqrs > EPS) then
              opposite_spin_cnt = opposite_spin_cnt + 1
              this%heg%opposite_spin(opposite_spin_cnt)%int_3 = diff_pr
              this%heg%opposite_spin(opposite_spin_cnt)%double = H_pqrs
              abs_H_copy(opposite_spin_cnt) = H_pqrs
            end if
          end do
        end do
      end do

      if (opposite_spin_cnt > 0) then
        call merge_sort_arg(-abs_H_copy, order, opposite_spin_cnt)
        hci_queue_copy(1: opposite_spin_cnt) = this%heg%opposite_spin(1: opposite_spin_cnt)
        do i = 1, opposite_spin_cnt
          this%heg%opposite_spin(i) = hci_queue_copy(order(i))
        end do
        max_abs_H = max(max_abs_H, this%heg%opposite_spin(1)%double)
      else
        this%heg%opposite_spin(1)%double = 0.0_DOUBLE
      end if

      this%max_abs_H = max_abs_H
    end subroutine generate_hci_queue

    type(real(DOUBLE)) function hamiltonian( &
        & this, det_pq_up, det_pq_dn, det_rs_up, det_rs_dn) result(H)
      class(heg_solver), intent(inout) :: this
      integer(INT64), intent(in) :: det_pq_up(:), det_pq_dn(:), det_rs_up(:), det_rs_dn(:)

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
      integer(INT64), allocatable :: eor_up(:), eor_dn(:)
      real(DOUBLE) :: k_unit
      real(DOUBLE) :: one_over_pi_l


      k_unit = this%heg%k_unit
      one_over_pi_l = 1.0_DOUBLE / (PI * this%heg%cell_length)

      if (abeq(det_pq_up, det_rs_up) .and. abeq(det_pq_dn, det_rs_dn)) then
        ! Diagonal elements.
        H = 0.0_DOUBLE
        n_occ_up = abcnt(det_pq_up)
        n_occ_dn = abcnt(det_pq_dn)
        call find_set_bits(det_pq_up, occ_pq_up, n_occ_up)
        call find_set_bits(det_pq_dn, occ_pq_dn, n_occ_dn)

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
        call abeor(det_pq_up, det_rs_up, eor_up)
        call abeor(det_pq_dn, det_rs_dn, eor_dn)
        n_eor_up = abcnt(eor_up)
        n_eor_dn = abcnt(eor_dn)

        if (n_eor_up + n_eor_dn /= 4) then
          H = 0.0_DOUBLE
          return
        end if

        call find_set_bits(eor_up, eor_up_set_bits, n_eor_up)
        call find_set_bits(eor_dn, eor_dn_set_bits, n_eor_dn)

        k_change = 0.0_DOUBLE
        k_p_set = .false.
        k_q_set = .false.
        k_s_set = .false.

        do i = 1, n_eor_up
          orb_id = eor_up_set_bits(i)
          if (abtest(det_pq_up, orb_id)) then
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
          if (abtest(det_pq_dn, orb_id)) then
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

        ! Check momentum conservation.
        if (sum(k_change**2) > 0) then
          H = 0.0_DOUBLE
          return
        end if

        H = one_over_pi_l / sum((this%heg%k_vectors(:, orb_p) - this%heg%k_vectors(:, orb_q))**2)

        if (n_eor_up /= 2) then
          ! Either two up or two dn excitations.
          H = H - one_over_pi_l / sum((this%heg%k_vectors(:, orb_p) - this%heg%k_vectors(:, orb_s))**2)
        end if

        call find_set_bits(det_pq_up, occ_pq_up)
        call find_set_bits(det_pq_dn, occ_pq_dn)
        call find_set_bits(det_rs_up, occ_rs_up)
        call find_set_bits(det_rs_dn, occ_rs_dn)
        
        gamma_exp = get_gamma_exp(det_pq_up, occ_pq_up, eor_up_set_bits, n_eor_up) &
            & + get_gamma_exp(det_pq_dn, occ_pq_dn, eor_dn_set_bits, n_eor_dn) &
            & + get_gamma_exp(det_rs_up, occ_rs_up, eor_up_set_bits, n_eor_up) &
            & + get_gamma_exp(det_rs_dn, occ_rs_dn, eor_dn_set_bits, n_eor_dn)

        if (iand(gamma_exp, 1) == 1) then
          H = -H
        end if
      end if
    end function hamiltonian

    type(real(DOUBLE)) function hamiltonian_abs_by_orbitals(this, p, q, r, s) result(abs_H)
      ! TODO: Optimization.
      class(heg_solver), intent(inout) :: this
      integer, intent(in) :: p, q, r, s
      integer(INT64), allocatable :: zero_det(:)
      integer(INT64), allocatable :: det_pq_up(:), det_pq_dn(:), det_rs_up(:), det_rs_dn(:)
      integer :: n_orb

      n_orb = this%n_orb

      if (p == q .or. r == s .or. p == r .or. q == s .or. p == s .or. q == r) then
        abs_H = 0.0_DOUBLE
        return
      end if

      allocate(zero_det(this%det_size))
      zero_det = 0_INT64
      det_pq_up = zero_det
      det_pq_dn = zero_det
      det_rs_up = zero_det
      det_rs_dn = zero_det

      if (p <= n_orb .and. q <= n_orb) then
        call abset(det_pq_up, p)
        call abset(det_pq_up, q)
        call abset(det_rs_up, r)
        call abset(det_rs_up, s)
      else if (p > n_orb .and. q > n_orb) then
        call abset(det_pq_dn, p - n_orb)
        call abset(det_pq_dn, q - n_orb)
        call abset(det_rs_dn, r - n_orb)
        call abset(det_rs_dn, s - n_orb)
      else
        call abset(det_pq_up, p)
        call abset(det_pq_dn, q - n_orb)
        call abset(det_rs_up, r)
        call abset(det_rs_dn, s - n_orb)
      end if
      
      abs_H = abs(this%hamiltonian(det_pq_up, det_pq_dn, det_rs_up, det_rs_dn))
    end function hamiltonian_abs_by_orbitals

    subroutine find_set_bits(det, set_bits_pos, n_set_bits)
      integer(INT64), intent(in) :: det(:)
      integer, allocatable, intent(out) :: set_bits_pos(:)
      integer, optional, intent(in) :: n_set_bits

      integer :: n_set_bits_ ! For subroutine internal use.
      integer :: pos
      integer :: i
      integer(INT64), allocatable :: det_tmp(:)

      if (present(n_set_bits)) then
        n_set_bits_ = n_set_bits
      else
        n_set_bits_ = abcnt(det)
      end if
      
      allocate(set_bits_pos(n_set_bits_))

      det_tmp = det
      do i = 1, n_set_bits_
        pos = abtrailz(det_tmp) + 1
        set_bits_pos(i) = pos 
        call abclr(det_tmp, pos)
      end do
    end subroutine find_set_bits

    type(integer) function get_gamma_exp(det, occ, eor_set_bits, n_eor_set_bits) result(gamma_exp)
      integer(INT64), intent(in) :: det(:)
      integer, intent(in) :: occ(:)
      integer, intent(in) :: eor_set_bits(:)
      integer, intent(in) :: n_eor_set_bits

      integer :: i
      integer :: ptr
      integer :: orb_id

      gamma_exp = 0
      ptr = 0
      do i = 1, n_eor_set_bits
        orb_id = eor_set_bits(i)
        if (.not. abtest(det, orb_id)) then
          cycle
        end if
        do while (occ(ptr + 1) < orb_id)
          ptr = ptr + 1
        end do
        gamma_exp = gamma_exp + ptr
      end do
    end function get_gamma_exp

end module heg_solver_class
