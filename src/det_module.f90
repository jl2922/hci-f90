module det_module

  use constants_module
  use types_module
  use utilities_module

  implicit none

  private

  public :: det_type
  public :: assignment(=), operator(==), operator(.eor.), operator(<)

  type det_type
    private
    integer(LONG), allocatable, public :: up(:)
    integer(LONG), allocatable, public :: dn(:)
    contains
      procedure, public :: initialize
      procedure, public :: get_orbital
      procedure, public :: set_orbital
      procedure, public :: get_n_elec
      procedure, public :: get_elec_orbitals
      procedure, public :: print_det
  end type det_type

  interface assignment(=)
    module procedure assign_det
  end interface

  interface operator(==)
    module procedure equal_det
  end interface

  interface operator(.eor.)
    module procedure eor_det
  end interface

  interface operator(<)
    module procedure lt_det
  end interface

  contains
  
  subroutine initialize(this, det_size)
    class(det_type), intent(inout) :: this
    integer, intent(in) :: det_size
    if (allocated(this%up)) deallocate(this%up)
    if (allocated(this%dn)) deallocate(this%dn)
    allocate(this%up(det_size))
    allocate(this%dn(det_size))
    this%up = 0
    this%dn = 0
  end subroutine initialize

  function get_orbital(this, spin, orbital_idx) result(is_occupied)
    class(det_type), intent(in) :: this
    integer, intent(in) :: spin
    integer, intent(in) :: orbital_idx
    logical :: is_occupied

    if (spin == C%UP_SPIN) then
      is_occupied = util%ab_get(this%up, orbital_idx)
    else
      is_occupied = util%ab_get(this%dn, orbital_idx)
    end if
  end function get_orbital

  subroutine set_orbital(this, spin, orbital_idx, occupation)
    class(det_type), intent(inout) :: this
    integer, intent(in) :: spin
    integer, intent(in) :: orbital_idx
    logical, intent(in) :: occupation

    if (spin == C%UP_SPIN) then
      call util%ab_set(this%up, orbital_idx, occupation)
    else
      call util%ab_set(this%dn, orbital_idx, occupation)
    end if
  end subroutine set_orbital

  function get_n_elec(this, spin) result(n)
    class(det_type), intent(in) :: this
    integer, intent(in) :: spin
    integer :: n
    integer :: i

    n = 0
    if (spin == C%UP_SPIN) then
      do i = 1, size(this%up)
        n = n + popcnt(this%up(i))
      end do
    else
      do i = 1, size(this%dn)
        n = n + popcnt(this%dn(i))
      end do
    end if
  end function get_n_elec

  subroutine get_elec_orbitals(this, spin, orbitals, n_elec_in)
    class(det_type), intent(in) :: this
    integer, intent(in) :: spin
    integer, allocatable, intent(out) :: orbitals(:)
    integer, optional, intent(in) :: n_elec_in
    integer :: n_elec
    integer :: i
    integer :: pos
    integer(LONG), allocatable :: spin_det(:)


    if (present(n_elec_in)) then
      n_elec = n_elec_in
    else
      n_elec = this%get_n_elec(spin)
    end if

    allocate(orbitals(n_elec))

    if (spin == C%UP_SPIN) then
      allocate(spin_det(size(this%up)))
      spin_det(:) = this%up(:)
    else
      allocate(spin_det(size(this%dn)))
      spin_det(:) = this%dn(:)
    end if

    do i = 1, n_elec
      pos = util%ab_trailz(spin_det) + 1
      orbitals(i) = pos
      call util%ab_set(spin_det, pos, .false.)
    end do
  end subroutine get_elec_orbitals

  subroutine print_det(this)
    class(det_type), intent(in) :: this
    integer :: i
    do i = 1, size(this%up)
      write (6, '(A, I0, A, B0.64)') 'UP#', i, ': ', this%up(i)
    end do
    do i = 1, size(this%dn)
      write (6, '(A, I0, A, B0.64)') 'DN#', i, ': ', this%dn(i)
    end do
  end subroutine print_det

  subroutine assign_det(det_dest, det_src)
    class(det_type), intent(in) :: det_src
    class(det_type), intent(out) :: det_dest
    integer :: size_up, size_dn

    size_up = size(det_src%up)
    size_dn = size(det_src%dn)

    ! If not allocated or size mismatch, reallocate. Otherwise, copy over.
    if (allocated(det_dest%up) .and. size(det_dest%up) /= size_up) then
      deallocate(det_dest%up)
    end if
    if (.not. allocated(det_dest%up)) then
      allocate(det_dest%up(size_up))
    end if
    det_dest%up(:) = det_src%up(:)

    if (allocated(det_dest%dn) .and. size(det_dest%dn) /= size_up) then
      deallocate(det_dest%dn)
    end if
    if (.not. allocated(det_dest%dn)) then
      allocate(det_dest%dn(size_up))
    end if
    det_dest%dn(:) = det_src%dn(:)
  end subroutine assign_det

  function equal_det(det1, det2) result(is_equal)
    class(det_type), intent(in) :: det1, det2
    logical :: is_equal
    if (util%ab_eq(det1%up, det2%up) .and. util%ab_eq(det1%dn, det2%dn)) then
      is_equal = .true.
    else
      is_equal = .false.
    end if
  end function equal_det

  function lt_det(det1, det2) result(is_lt)
    class(det_type), intent(in) :: det1, det2
    logical :: is_lt
    if (util%ab_lt(det1%up, det2%up) .or. &
        & (util%ab_eq(det1%up, det2%up) .and. util%ab_lt(det1%dn, det2%dn))) then
      is_lt = .true.
    else
      is_lt = .false.
    end if
  end function lt_det
  
  function eor_det(det1, det2) result(det_res)
    class(det_type), intent(in) :: det1, det2
    type(det_type) :: det_res
    integer :: size1_up, size1_dn, size2_up, size2_dn

    size1_up = size(det1%up)
    size1_dn = size(det1%dn)
    size2_up = size(det2%up)
    size2_dn = size(det2%dn)
    if (size1_up /= size2_up) stop 'eor_det up spin size mismatch.'
    if (size1_dn /= size2_dn) stop 'eor_det dn spin size mismatch.'
    allocate(det_res%up(size1_up))
    allocate(det_res%dn(size1_dn))
    call util%ab_eor(det1%up, det2%up, det_res%up)
    call util%ab_eor(det1%dn, det2%dn, det_res%dn)
  end function eor_det

end module det_module
