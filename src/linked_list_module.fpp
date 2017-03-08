#:set shortname {'integer': 'int', 'real': 'real', 'real(DOUBLE)': 'double'}
#:for dtype in ['integer', 'real(DOUBLE)']
module linked_list_module__${shortname[dtype]}$
  ! This is a linked list implementation that optimized for speed but does not 
  ! allow removing items.

#:if dtype == 'real(DOUBLE)'
  use types_module
#:endif

  implicit none

  private

  public :: linked_list_type__${shortname[dtype]}$
  public :: new_linked_list__${shortname[dtype]}$
  public :: new_linked_list_arr__${shortname[dtype]}$
  public :: delete
  public :: assignment(=)

  integer, parameter :: DEFAULT_BLOCK_SIZE = 1024

  type :: node_type
    private
    ${dtype}$, allocatable :: block(:)
    type(node_type), pointer :: next => null()
  end type node_type

  type :: linked_list_type__${shortname[dtype]}$
    private
    integer :: n = 0
    integer :: cur_idx = 0
    integer :: block_idx_cache = 0 ! Item index inside a block.
    integer :: block_size = DEFAULT_BLOCK_SIZE
    type(node_type), pointer :: head => null()
    type(node_type), pointer :: cur_node => null()
    contains     
      procedure, public :: append ! Append at the back.
      procedure, public :: begin ! Returns the iterator to the beginning of the list.
      procedure, public :: get ! Get the current item pointed by the iterator.
      procedure, public :: next ! Advance the iterator.
      procedure, public :: is_empty
      procedure, public :: get_length
  end type linked_list_type__${shortname[dtype]}$

  interface delete
    module procedure delete_linked_list__${shortname[dtype]}$
  end interface delete
  
  interface assignment(=)
    module procedure assign_linked_list
  end interface

  contains

  function new_linked_list__${shortname[dtype]}$( &
      & block_size_opt) result(list)
    type(linked_list_type__${shortname[dtype]}$), pointer :: list
    integer, optional, intent(in) :: block_size_opt

    allocate(list)
    if (present(block_size_opt)) then
      list%block_size = block_size_opt
    endif
  end function new_linked_list__${shortname[dtype]}$

  function new_linked_list_arr__${shortname[dtype]}$( &
      & n, block_size_opt) result(list)
    type(linked_list_type__${shortname[dtype]}$), pointer :: list(:)
    integer, intent(in) :: n
    integer, optional, intent(in) :: block_size_opt
    integer :: i

    allocate(list(n))
    if (present(block_size_opt)) then
      do i = 1, n
        list(i)%block_size = block_size_opt
      enddo
    endif
  end function new_linked_list_arr__${shortname[dtype]}$

  subroutine delete_linked_list__${shortname[dtype]}$(list)
    type(linked_list_type__${shortname[dtype]}$), pointer, intent(inout) :: list
    type(node_type), pointer :: cur_block, next_block
    
    next_block => list%head
    do while(associated(next_block))
      cur_block => next_block
      next_block => next_block%next
      deallocate(cur_block)
    enddo
    deallocate(list)
    nullify(list)
  end subroutine delete_linked_list__${shortname[dtype]}$

  subroutine delete_linked_list_arr__${shortname[dtype]}$(list)
    type(linked_list_type__${shortname[dtype]}$), pointer, intent(inout) :: list(:)
    type(node_type), pointer :: cur_block, next_block
    integer :: i
    
    do i = 1, size(list)
      if (.not. associated(list(i)%head)) cycle
      next_block => list(i)%head
      do while(associated(next_block))
        cur_block => next_block
        next_block => next_block%next
        deallocate(cur_block)
      enddo
    enddo
    deallocate(list)
    nullify(list)
  end subroutine delete_linked_list_arr__${shortname[dtype]}$

  subroutine append(this, item)
    class(linked_list_type__${shortname[dtype]}$), intent(inout) :: this
    ${dtype}$, intent(in) :: item
    if (this%cur_idx /= this%n) then
      call backtrace
      stop 'cur_idx is not at the end of the list.'
    endif
    if (.not. associated(this%head)) then
      allocate(this%cur_node)
      allocate(this%cur_node%block(this%block_size))
      this%head => this%cur_node
    else if (this%block_idx_cache == this%block_size) then
      allocate(this%cur_node%next)
      this%cur_node => this%cur_node%next
      allocate(this%cur_node%block(this%block_size))
      this%block_idx_cache = 0
    endif
    this%cur_idx = this%cur_idx + 1
    this%block_idx_cache = this%block_idx_cache + 1
    this%n = this%n + 1
    this%cur_node%block(this%block_idx_cache) = item
  end subroutine append

  subroutine begin(this)
    class(linked_list_type__${shortname[dtype]}$), intent(inout) :: this
    this%cur_node => this%head
    this%cur_idx = 1
    this%block_idx_cache = 1
  end subroutine begin
  
  function get(this) result(item)
    class(linked_list_type__${shortname[dtype]}$), intent(in) :: this
    ${dtype}$ :: item
    if (this%cur_idx > this%n) then
      call backtrace
      stop 'Index out of bound.'
    endif
    item = this%cur_node%block(this%block_idx_cache)
  end function get

  subroutine next(this, stat)
    class(linked_list_type__${shortname[dtype]}$), intent(inout) :: this
    logical, optional, intent(out) :: stat
    if (this%cur_idx >= this%n) then
      if (present(stat)) stat = .false.
      return
    endif
    this%cur_idx = this%cur_idx + 1
    this%block_idx_cache = this%block_idx_cache + 1
    if (this%block_idx_cache > this%block_size) then
      this%block_idx_cache = 1
      this%cur_node => this%cur_node%next 
    endif
    if(present(stat)) stat = .true.
  end subroutine next

  subroutine assign_linked_list(dest, src)
    type(linked_list_type__${shortname[dtype]}$), pointer, intent(in) :: src
    type(linked_list_type__${shortname[dtype]}$), pointer, intent(out) :: dest

    if (.not. src%is_empty()) then
      stop 'Not implemented yet.'
    endif
    dest%block_size = src%block_size
  end subroutine assign_linked_list

  function is_empty(this)
    class(linked_list_type__${shortname[dtype]}$), intent(inout) :: this
    logical :: is_empty

    is_empty = (this%n == 0)
  end function is_empty

  function get_length(this)
    class(linked_list_type__${shortname[dtype]}$), intent(inout) :: this
    integer :: get_length

    get_length = this%n
  end function get_length

end module linked_list_module__${shortname[dtype]}$
#:endfor
