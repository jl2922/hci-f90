module linked_list_module__int


  implicit none

  private

  public :: linked_list_type__int
  public :: new_linked_list__int

  integer, parameter :: TRUNK_SIZE = 65536 ! Fits into the CPU L2 cache we have.
  ! integer, parameter :: TRUNK_SIZE = 1024 ! For testing.

  type :: node_type
    private
    integer, dimension(TRUNK_SIZE) :: trunk
    type(node_type), pointer :: next => null()
  end type node_type

  type :: linked_list_type__int
    private
    integer :: n = 0
    integer :: idx = 0
    integer :: trunk_idx = 0 ! Item index inside a trunk.
    type(node_type), pointer :: head => null()
    type(node_type), pointer :: node => null()
    contains     
      procedure, public :: append ! Append at the back.
      procedure, public :: begin ! Returns the iterator to the beginning of the list.
      procedure, public :: get ! Get the current item pointed by the iterator.
      procedure, public :: next ! Advance the iterator.
      procedure, public :: clear
  end type linked_list_type__int

  interface new_linked_list__int
    module procedure new_linked_list__int_default
  end interface

  contains

  function new_linked_list__int_default() result(list)
    type(linked_list_type__int), pointer :: list

    allocate(list)
  end function new_linked_list__int_default

  subroutine append(this, item)
    class(linked_list_type__int), intent(inout) :: this
    integer, intent(in) :: item
    if (.not. associated(this%head)) then
      allocate(this%node)
      this%head => this%node
      this%idx = 0
      this%trunk_idx = 0
      this%n = 0
    else if (this%trunk_idx == TRUNK_SIZE) then
      allocate(this%node%next)
      this%node => this%node%next
      this%trunk_idx = 0
    endif
    this%idx = this%idx + 1
    this%trunk_idx = this%trunk_idx + 1
    this%n = this%n + 1
    this%node%trunk(this%trunk_idx) = item
  end subroutine append

  subroutine begin(this)
    class(linked_list_type__int), intent(inout) :: this
    this%node => this%head
    this%idx = 1
    this%trunk_idx = 1
  end subroutine begin
  
  function get(this) result(item)
    class(linked_list_type__int), intent(in) :: this
    integer :: item
    if (this%idx > this%n) stop 'Index out of bound.'
    item = this%node%trunk(this%trunk_idx)
  end function get

  subroutine next(this, stat)
    class(linked_list_type__int), intent(inout) :: this
    logical, optional, intent(out) :: stat
    if (this%idx >= this%n) then
      if (present(stat)) stat = .false.
      return
    endif
    this%idx = this%idx + 1
    this%trunk_idx = this%trunk_idx + 1
    if (this%trunk_idx > TRUNK_SIZE) then
      this%trunk_idx = 1
      this%node => this%node%next 
    endif
    if(present(stat)) stat = .false.
  end subroutine next

  subroutine clear(this)
    class(linked_list_type__int), intent(inout) :: this
    type(node_type), pointer :: cur_trunk
    type(node_type), pointer :: next_trunk
    if (.not. associated(this%head)) return
    next_trunk => this%head
    do while (associated(next_trunk))
     cur_trunk => next_trunk
     next_trunk => next_trunk%next
     deallocate(cur_trunk)
    enddo
    nullify(this%head)
    nullify(this%node)
    this%n = 0
  end subroutine clear

end module linked_list_module__int
module linked_list_module__double

  use types_module

  implicit none

  private

  public :: linked_list_type__double
  public :: new_linked_list__double

  integer, parameter :: TRUNK_SIZE = 65536 ! Fits into the CPU L2 cache we have.
  ! integer, parameter :: TRUNK_SIZE = 1024 ! For testing.

  type :: node_type
    private
    real(DOUBLE), dimension(TRUNK_SIZE) :: trunk
    type(node_type), pointer :: next => null()
  end type node_type

  type :: linked_list_type__double
    private
    integer :: n = 0
    integer :: idx = 0
    integer :: trunk_idx = 0 ! Item index inside a trunk.
    type(node_type), pointer :: head => null()
    type(node_type), pointer :: node => null()
    contains     
      procedure, public :: append ! Append at the back.
      procedure, public :: begin ! Returns the iterator to the beginning of the list.
      procedure, public :: get ! Get the current item pointed by the iterator.
      procedure, public :: next ! Advance the iterator.
      procedure, public :: clear
  end type linked_list_type__double

  interface new_linked_list__double
    module procedure new_linked_list__double_default
  end interface

  contains

  function new_linked_list__double_default() result(list)
    type(linked_list_type__double), pointer :: list

    allocate(list)
  end function new_linked_list__double_default

  subroutine append(this, item)
    class(linked_list_type__double), intent(inout) :: this
    real(DOUBLE), intent(in) :: item
    if (.not. associated(this%head)) then
      allocate(this%node)
      this%head => this%node
      this%idx = 0
      this%trunk_idx = 0
      this%n = 0
    else if (this%trunk_idx == TRUNK_SIZE) then
      allocate(this%node%next)
      this%node => this%node%next
      this%trunk_idx = 0
    endif
    this%idx = this%idx + 1
    this%trunk_idx = this%trunk_idx + 1
    this%n = this%n + 1
    this%node%trunk(this%trunk_idx) = item
  end subroutine append

  subroutine begin(this)
    class(linked_list_type__double), intent(inout) :: this
    this%node => this%head
    this%idx = 1
    this%trunk_idx = 1
  end subroutine begin
  
  function get(this) result(item)
    class(linked_list_type__double), intent(in) :: this
    real(DOUBLE) :: item
    if (this%idx > this%n) stop 'Index out of bound.'
    item = this%node%trunk(this%trunk_idx)
  end function get

  subroutine next(this, stat)
    class(linked_list_type__double), intent(inout) :: this
    logical, optional, intent(out) :: stat
    if (this%idx >= this%n) then
      if (present(stat)) stat = .false.
      return
    endif
    this%idx = this%idx + 1
    this%trunk_idx = this%trunk_idx + 1
    if (this%trunk_idx > TRUNK_SIZE) then
      this%trunk_idx = 1
      this%node => this%node%next 
    endif
    if(present(stat)) stat = .false.
  end subroutine next

  subroutine clear(this)
    class(linked_list_type__double), intent(inout) :: this
    type(node_type), pointer :: cur_trunk
    type(node_type), pointer :: next_trunk
    if (.not. associated(this%head)) return
    next_trunk => this%head
    do while (associated(next_trunk))
     cur_trunk => next_trunk
     next_trunk => next_trunk%next
     deallocate(cur_trunk)
    enddo
    nullify(this%head)
    nullify(this%node)
    this%n = 0
  end subroutine clear

end module linked_list_module__double