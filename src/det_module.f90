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
  public :: operator(.eor.)

  type det_type
    private
    type(spin_det_type), public :: up, dn
    contains
      procedure, public :: initialize
      procedure, public :: print
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

  interface operator(>)
    module procedure lt_det
  end interface

  contains
  
  subroutine initialize(this, det_size)
    class(det_type), intent(inout) :: this
    integer, intent(in) :: det_size
    call this%up%initialize(det_size)
    call this%dn%initialize(det_size)
  end subroutine initialize

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

  function eor_det(left, right) result(res)
    class(det_type), intent(in) :: left, right
    type(det_type) :: res
    res%up = left%up .eor. right%up
    res%dn = left%dn .eor. right%dn
  end function eor_det

end module det_module
