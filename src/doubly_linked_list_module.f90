module doubly_linked_list_module

  use det_module

  implicit none

  private

  public :: doubly_linked_list_node_type
  public :: doubly_linked_list_type
  public :: delete
  public :: new_doubly_linked_list
  public :: new_doubly_linked_list_node_arr

  type doubly_linked_list_node_type
    type(doubly_linked_list_node_type), pointer :: prev => null()
    type(doubly_linked_list_node_type), pointer :: next => null()
    type(det_type), pointer :: item => null()
  end type doubly_linked_list_node_type

  type doubly_linked_list_type
    type(doubly_linked_List_node_type), pointer :: head => null()
    type(doubly_linked_List_node_type), pointer :: tail => null()
    integer :: n
    contains
      generic, public :: push_front => push_front_item, push_front_node
      procedure, public :: pop_back
      procedure, public :: front
      procedure, public :: back
      procedure, public :: remove
      procedure, public :: get_size
      procedure :: push_front_item
      procedure :: push_front_node
  end type doubly_linked_list_type

  interface delete
    module procedure delete_list
    module procedure delete_node
    module procedure delete_node_arr
  end interface

  contains

  function new_node(item) result(node)
    type(det_type), pointer, intent(in) :: item
    type(doubly_linked_list_node_type), pointer :: node

    allocate(node)
    node%item => item
  end function new_node

  function new_doubly_linked_list_node_arr(n) result(arr)
    integer, intent(in) :: n
    type(doubly_linked_list_node_type), pointer :: arr(:)

    allocate(arr(n))
  end function new_doubly_linked_list_node_arr

  function new_doubly_linked_list() result(list)
    type(doubly_linked_list_type), pointer :: list

    allocate(list)
    allocate(list%head)
    allocate(list%tail)
    list%head%next => list%tail
    list%tail%prev => list%head
  end function new_doubly_linked_list

  subroutine delete_list(list)
    type(doubly_linked_list_type), pointer, intent(inout) :: list
    type(doubly_linked_list_node_type), pointer :: cur, next

    if (.not. associated(list)) return
    cur => list%head
    do while (associated(cur))
      next => cur%next
      call delete(cur)
      cur => next
    enddo
    deallocate(list)
    nullify(list)
  end subroutine delete_list

  subroutine delete_node(node)
    type(doubly_linked_list_node_type), pointer, intent(inout) :: node

    if (.not. associated(node)) return
    deallocate(node)
    nullify(node)
  end subroutine delete_node

  subroutine delete_node_arr(arr)
    type(doubly_linked_list_node_type), pointer, intent(inout) :: arr(:)

    if (.not. associated(arr)) return
    deallocate(arr)
    nullify(arr)
  end subroutine delete_node_arr

  subroutine push_front_item(this, item)
    class(doubly_linked_list_type), intent(inout) :: this
    type(det_type), pointer, intent(in) :: item
    type(doubly_linked_list_node_type), pointer :: node

    node => new_node(item)
    call this%push_front(node)
  end subroutine push_front_item

  subroutine push_front_node(this, node)
    class(doubly_linked_list_type), intent(inout) :: this
    type(doubly_linked_list_node_type), pointer :: node

    node%next => this%head%next
    node%prev => this%head
    this%head%next => node
    node%next%prev => node
    this%n = this%n + 1
  end subroutine push_front_node    

  subroutine pop_back(this)
    class(doubly_linked_list_type), intent(inout) :: this
    type(doubly_linked_list_node_type), pointer :: node

    if (this%n <= 0) return
    node => this%tail%prev
    call this%remove(node)
  end subroutine pop_back

  subroutine remove(this, node)
    class(doubly_linked_list_type), intent(inout) :: this 
    type(doubly_linked_list_node_type), pointer, intent(in) :: node

    node%prev%next => node%next
    node%next%prev => node%prev
    this%n = this%n - 1
  end subroutine remove

  function front(this)
    class(doubly_linked_list_type), intent(inout) :: this
    type(doubly_linked_list_node_type), pointer :: front

    front => this%head%next
  end function front

  function back(this)
    class(doubly_linked_list_type), intent(inout) :: this
    type(doubly_linked_list_node_type), pointer :: back

    back => this%tail%prev
  end function back

  function get_size(this) result(res)
    class(doubly_linked_list_type), intent(inout) :: this
    integer :: res

    res = this%n
  end function get_size

end module doubly_linked_list_module
