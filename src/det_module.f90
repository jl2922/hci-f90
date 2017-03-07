module det_module

  use constants_module
  use spin_det_module
  use types_module

  implicit none

  private

  public :: det_type
  public :: new_det
  public :: delete
  public :: assignment(=)
  public :: operator(==)
  public :: operator(<)
  public :: operator(>)

  type det_type
    private
    type(spin_det_type), pointer, public :: up => null()
    type(spin_det_type), pointer, public :: dn => null()
    contains
      procedure, public :: print
      procedure, public :: from_eor
      procedure, public :: get_n_diff_orbitals
  end type det_type

  interface new_det
    module procedure new_det_clone
    module procedure new_det_by_n_orbs
  end interface new_det

  interface delete
    module procedure delete_det
    module procedure delete_det_arr
  end interface delete

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
    module procedure gt_det
  end interface

  contains

  function new_det_clone(src) result(det)
    type(det_type), pointer, intent(in) :: src
    type(det_type), pointer :: det
    
    allocate(det)
    det%up => new_spin_det(src%up)
    det%dn => new_spin_det(src%dn)
  end function new_det_clone

  function new_det_by_n_orbs(n_orbs) result(det)
    integer, intent(in) :: n_orbs
    type(det_type), pointer :: det

    allocate(det)
    det%up => new_spin_det(n_orbs)
    det%dn => new_spin_det(n_orbs)
  end function new_det_by_n_orbs

  subroutine delete_det(det)
    type(det_type), pointer :: det

    if (.not. associated(det)) return
    call delete(det%up)
    call delete(det%dn)
    deallocate(det)
    nullify(det)
  end subroutine delete_det

  subroutine delete_det_arr(det_arr)
    type(det_type), pointer :: det_arr(:)
    integer :: i

    if (.not. associated(det_arr)) return
    do i = 1, size(det_arr)
      call delete(det_arr(i)%up)
      call delete(det_arr(i)%dn)
    enddo
    deallocate(det_arr)
    nullify(det_arr)
  end subroutine delete_det_arr

  subroutine print(this)
    class(det_type), intent(in) :: this
    integer :: i

    call this%up%print()
    call this%dn%print()
  end subroutine print

  subroutine from_eor(this, det1, det2)
    class(det_type), intent(inout) :: this
    type(det_type), pointer, intent(in) :: det1, det2

    call this%up%from_eor(det1%up, det2%up)
    call this%dn%from_eor(det1%dn, det2%dn)
  end subroutine from_eor

  function get_n_diff_orbitals(this, op) result(n_diff)
    class(det_type), intent(inout) :: this
    type(det_type), pointer, intent(in) :: op
    integer :: n_diff

    n_diff = this%up%get_n_diff_orbitals(op%up) + &
        & this%dn%get_n_diff_orbitals(op%dn)
  end function get_n_diff_orbitals

  subroutine assign_det(dest, src)
    type(det_type), pointer, intent(in) :: src
    type(det_type), pointer, intent(out) :: dest

    dest%up = src%up
    dest%dn = src%dn
  end subroutine assign_det

  function equal_det(left, right) result(is_equal)
    type(det_type), pointer, intent(in) :: left, right
    logical :: is_equal

    is_equal = (left%up == right%up) .and. (left%dn == right%dn)
  end function equal_det

  function lt_det(left, right) result(is_lt)
    type(det_type), pointer, intent(in) :: left, right
    logical :: is_lt

    is_lt = (left%up < right%up) .or. &
        & ((left%up == right%up) .and. (left%dn < right%dn))
  end function lt_det
  
  function gt_det(left, right) result(is_gt)
    type(det_type), pointer, intent(in) :: left, right
    logical :: is_gt

    is_gt = (left%up > right%up) .or. &
        & ((left%up == right%up) .and. (left%dn > right%dn))
  end function gt_det

end module det_module
