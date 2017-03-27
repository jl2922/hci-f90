module det_module

  use constants_module
  use spin_det_module
  use types_module

  implicit none

  private

  public :: det_type
  public :: build
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
      procedure, public :: is_empty
      procedure, public :: get_hash
      procedure, public :: get_n_diff_orbitals
  end type det_type

  interface build
    module procedure build_det_clone
    module procedure build_det_by_n_orbs
    module procedure build_det_arr
  end interface build

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

  subroutine build_det_clone(this, src)
    type(det_type), pointer, intent(inout) :: this
    type(det_type), pointer, intent(in) :: src
    
    allocate(this)
    call build(this%up, src%up)
    call build(this%dn, src%dn)
  end subroutine build_det_clone

  subroutine build_det_by_n_orbs(this, n_orbs)
    type(det_type), pointer, intent(inout) :: this
    integer, intent(in) :: n_orbs

    allocate(this)
    call build(this%up, n_orbs)
    call build(this%dn, n_orbs)
  end subroutine build_det_by_n_orbs

  subroutine build_det_arr(this, n)
    type(det_type), pointer, intent(inout) :: this(:)
    integer, intent(in) :: n

    allocate(this(n))
  end subroutine build_det_arr

  subroutine delete_det(this)
    type(det_type), pointer :: this

    if (.not. associated(this)) return
    call delete(this%up)
    call delete(this%dn)
    deallocate(this)
    nullify(this)
  end subroutine delete_det

  subroutine delete_det_arr(this)
    type(det_type), pointer :: this(:)
    integer :: i

    if (.not. associated(this)) return
    do i = 1, size(this)
      call delete(this(i)%up)
      call delete(this(i)%dn)
    enddo
    deallocate(this)
    nullify(this)
  end subroutine delete_det_arr

  subroutine print(this)
    class(det_type), intent(in) :: this

    call this%up%print()
    call this%dn%print()
  end subroutine print

  subroutine from_eor(this, det1, det2)
    class(det_type), intent(inout) :: this
    type(det_type), pointer, intent(in) :: det1, det2

    call this%up%from_eor(det1%up, det2%up)
    call this%dn%from_eor(det1%dn, det2%dn)
  end subroutine from_eor

  function is_empty(this)
    class(det_type), intent(inout) :: this
    logical :: is_empty

    is_empty = (.not. associated(this%up)) .and. (.not. associated(this%dn))
  end function is_empty

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

  function get_hash(this, table_size) result(hash)
    class(det_type), intent(inout) :: this
    integer, intent(in) :: table_size
    integer :: hash

    hash = this%up%get_hash(table_size) + this%dn%get_hash(table_size)
    if (hash > table_size) then
      hash = hash - table_size
    endif
  end function get_hash

end module det_module
