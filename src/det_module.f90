module det_module

  use constants_module
  use spin_det_module
  use types_module

  implicit none

  private

  public :: det_type
  public :: new_det
  public :: assignment(=)
  public :: operator(==)
  public :: operator(<)
  public :: operator(>)
  public :: operator(.eor.)

  type det_type
    private
    type(spin_det_type), pointer, public :: up, dn
    ! integer, pointer, public :: elec_orbitals(:)
    contains
      procedure, public :: clean
      procedure, public :: resize
      procedure, public :: print
  end type det_type

  interface new_det
    module procedure new_det_clone
    module procedure new_det_by_n_orb
  end interface new_det

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

  interface operator(.eor.)
    module procedure eor_det
  end interface

  contains

  function new_det_clone(src) result(det)
    type(det_type), pointer, intent(in) :: src
    type(det_type), pointer :: det
    integer :: det_size
    
    allocate(det)
    det%up => new_spin_det(src%up)
    det%dn => new_spin_det(src%dn)
  end function new_det_clone

  function new_det_by_n_orb(n_orb) result(det)
    integer, intent(in) :: n_orb 
    type(det_type), pointer :: det
    integer :: det_size

    allocate(det)
    det_size = ceiling(n_orb * 1.0 / C%TRUNK_SIZE)
    det%up => new_spin_det(det_size)
    det%dn => new_spin_det(det_size)
  end function new_det_by_n_orb

  subroutine clean(this)
    class(det_type), intent(inout) :: this

    call this%up%clean()
    call this%dn%clean()
    deallocate(this%up)
    deallocate(this%dn)
  end subroutine clean

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
    class(det_type), pointer, intent(in) :: src
    class(det_type), pointer, intent(out) :: dest

    dest%up = src%up
    dest%dn = src%dn
  end subroutine assign_det

  function equal_det(left, right) result(is_equal)
    class(det_type), pointer, intent(in) :: left, right
    logical :: is_equal

    is_equal = (left%up == right%up) .and. (left%dn == right%dn)
  end function equal_det

  function lt_det(left, right) result(is_lt)
    class(det_type), pointer, intent(in) :: left, right
    logical :: is_lt

    is_lt = (left%up < right%up) .or. &
        & ((left%up == right%up) .and. (left%dn < right%dn))
  end function lt_det
  
  function gt_det(left, right) result(is_gt)
    class(det_type), pointer, intent(in) :: left, right
    logical :: is_gt

    is_gt = (left%up < right%up) .or. &
        & ((left%up == right%up) .and. (left%dn < right%dn))
  end function gt_det

  function eor_det(left, right) result(res)
    type(det_type), pointer, intent(in) :: left, right
    type(det_type), pointer :: res
    
    allocate(res)
    res%up => left%up .eor. right%up
    res%dn => left%dn .eor. right%dn
  end function eor_det

end module det_module
