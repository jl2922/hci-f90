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
      procedure, public :: resize
      procedure, public :: print
      procedure, public :: from_eor
  end type det_type

  interface new_det
    module procedure new_det_clone
    module procedure new_det_by_n_orb
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
    module procedure lt_det
  end interface

  contains

  function new_det_clone(src) result(det)
    type(det_type), pointer, intent(in) :: src
    type(det_type), pointer :: det
    
    allocate(det)
    det%up => new_spin_det(src%up)
    det%dn => new_spin_det(src%dn)
  end function new_det_clone

  function new_det_by_n_orb(n_orb) result(det)
    integer, intent(in) :: n_orb 
    type(det_type), pointer :: det
    integer :: n_trunks

    allocate(det)
    n_trunks = ceiling(n_orb * 1.0 / C%TRUNK_SIZE)
    det%up => new_spin_det(n_trunks)
    det%dn => new_spin_det(n_trunks)
  end function new_det_by_n_orb

  subroutine delete_det(det)
    type(det_type), pointer :: det

    call delete(det%up)
    call delete(det%dn)
    deallocate(det)
    nullify(det)
  end subroutine delete_det

  subroutine delete_det_arr(arr)
    type(det_type), pointer :: arr(:)
    integer :: i

    do i = 1, size(arr)
      if (associated(arr(i)%up)) then
        call delete(arr(i)%up)
      endif
      if (associated(arr(i)%dn)) then
        call delete(arr(i)%dn)
      endif
    enddo
    deallocate(arr)
    nullify(arr)
  end subroutine delete_det_arr

  subroutine resize(this, n_trunks)
    class(det_type), intent(inout) :: this
    integer, intent(in) :: n_trunks

    call this%up%resize(n_trunks)
    call this%dn%resize(n_trunks)
  end subroutine resize

  subroutine print(this)
    class(det_type), intent(in) :: this
    integer :: i

    do i = 1, this%up%n_trunks
      write (6, '(A, I0, A, B0.63)') 'up#', i, ': ', this%up%trunks(i)
    enddo
    do i = 1, this%dn%n_trunks
      write (6, '(A, I0, A, B0.63)') 'dn#', i, ': ', this%dn%trunks(i)
    enddo
  end subroutine print

  subroutine from_eor(this, det1, det2)
    class(det_type), intent(inout) :: this
    type(det_type), pointer, intent(in) :: det1, det2

    call this%up%from_eor(det1%up, det2%up)
    call this%dn%from_eor(det1%dn, det2%dn)
  end subroutine from_eor

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
