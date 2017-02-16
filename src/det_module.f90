module det_module

  use constants_module
  use spin_det_module
  use types_module

  implicit none

  private

  public :: det_type
  public :: assignment(=)
  public :: operator(==)
  public :: operator(<)
  public :: operator(>)

  type det_type
    private
    type(spin_det_type), public :: up, dn
    ! integer, pointer, public :: elec_orbitals(:)
    contains
      procedure, public :: from_eor
      procedure, public :: initialize
      procedure, public :: resize
      procedure, public :: print
  end type det_type

  interface assignment(=)
    module procedure assign_det
  end interface

  interface operator(==)
    module procedure equal_det
  end interface

  interface operator(<)
    module procedure lt_det
  end interface

  interface operator(>)
    module procedure lt_det
  end interface

  contains
  
  subroutine initialize(this, n_orb)
    class(det_type), intent(inout) :: this
    integer, intent(in) :: n_orb
    integer :: det_size
    det_size = ceiling(n_orb * 1.0 / C%TRUNK_SIZE)
    call this%up%initialize_by_det_size(det_size)
    call this%dn%initialize_by_det_size(det_size)
  end subroutine initialize

  subroutine resize(this, det_size)
    class(det_type), intent(inout) :: this
    integer, intent(in) :: det_size

    call this%up%resize(det_size)
    call this%dn%resize(det_size)
  end subroutine resize

  subroutine print(this)
    class(det_type), intent(in) :: this

    call this%up%print()
    call this%dn%print()
  end subroutine print

  subroutine assign_det(dest, src)
    class(det_type), intent(in) :: src
    class(det_type), intent(out) :: dest

    dest%up = src%up
    dest%dn = src%dn
  end subroutine assign_det

  function equal_det(left, right) result(is_equal)
    class(det_type), intent(in) :: left, right
    logical :: is_equal

    is_equal = (left%up == right%up) .and. (left%dn == right%dn)
  end function equal_det

  function lt_det(left, right) result(is_lt)
    class(det_type), intent(in) :: left, right
    logical :: is_lt

    is_lt = (left%up < right%up) .or. &
        & ((left%up == right%up) .and. (left%dn < right%dn))
  end function lt_det
  
  function gt_det(left, right) result(is_gt)
    class(det_type), intent(in) :: left, right
    logical :: is_gt

    is_gt = (left%up < right%up) .or. &
        & ((left%up == right%up) .and. (left%dn < right%dn))
  end function gt_det

  subroutine from_eor(this, left, right)
    class(det_type), intent(inout) :: this
    type(det_type), intent(in) :: left, right

    call this%up%from_eor(left%up, right%up)
    call this%dn%from_eor(left%dn, right%dn)
  end subroutine from_eor

end module det_module
