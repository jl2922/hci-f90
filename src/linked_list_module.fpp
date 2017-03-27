#:set dtypes ['int', 'double']
#:set type_defs {'int': 'integer', 'double': 'real(DOUBLE)'}
#:for dtype in dtypes
#:set type_def type_defs[dtype]
module linked_list_module__${dtype}$
  ! This is a linked list implementation that optimized for speed but does not 
  ! allow removing items.

#:if dtype == 'double'
  use types_module
#:endif

  implicit none

  private

  public :: linked_list_type__${dtype}$
  public :: build
  public :: delete
  public :: assignment(=)

  integer, parameter :: DEFAULT_BLOCK_SIZE = 1024

  type :: node_type
    private
    ${type_def}$, allocatable :: block(:)
    type(node_type), pointer :: next => null()
  end type node_type

  type :: linked_list_type__${dtype}$
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
  end type linked_list_type__${dtype}$

  interface build
    module procedure build_list
    module procedure build_list_arr
  end interface build

  interface delete
    module procedure delete_list
    module procedure delete_list_arr
  end interface delete
  
  interface assignment(=)
    module procedure assign_linked_list
  end interface

  contains

  subroutine build_list(this, block_size_opt)
    type(linked_list_type__${dtype}$), pointer, intent(inout) :: this
    integer, optional, intent(in) :: block_size_opt

    allocate(this)
    if (present(block_size_opt)) then
      this%block_size = block_size_opt
    endif
  end subroutine build_list 

  subroutine build_list_arr(this, n, block_size_opt)
    type(linked_list_type__${dtype}$), pointer, intent(inout) :: this(:)
    integer, intent(in) :: n
    integer, optional, intent(in) :: block_size_opt
    integer :: i

    allocate(this(n))
    if (present(block_size_opt)) then
      do i = 1, n
        this(i)%block_size = block_size_opt
      enddo
    endif
  end subroutine build_list_arr 

  subroutine delete_list(this)
    type(linked_list_type__${dtype}$), pointer, intent(inout) :: this
    type(node_type), pointer :: cur_block, next_block
    
    next_block => this%head
    do while(associated(next_block))
      cur_block => next_block
      next_block => next_block%next
      deallocate(cur_block)
    enddo
    deallocate(this)
    nullify(this)
  end subroutine delete_list

  subroutine delete_list_arr(this)
    type(linked_list_type__${dtype}$), pointer, intent(inout) :: this(:)
    type(node_type), pointer :: cur_block, next_block
    integer :: i
    
    do i = 1, size(this)
      if (.not. associated(this(i)%head)) cycle
      next_block => this(i)%head
      do while(associated(next_block))
        cur_block => next_block
        next_block => next_block%next
        deallocate(cur_block)
      enddo
    enddo
    deallocate(this)
    nullify(this)
  end subroutine delete_list_arr

  subroutine append(this, item)
    class(linked_list_type__${dtype}$), intent(inout) :: this
    ${type_def}$, intent(in) :: item
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
    class(linked_list_type__${dtype}$), intent(inout) :: this
    this%cur_node => this%head
    this%cur_idx = 1
    this%block_idx_cache = 1
  end subroutine begin
  
  function get(this) result(item)
    class(linked_list_type__${dtype}$), intent(in) :: this
    ${type_def}$ :: item
    if (this%cur_idx > this%n) then
      call backtrace
      stop 'Index out of bound.'
    endif
    item = this%cur_node%block(this%block_idx_cache)
  end function get

  subroutine next(this, stat)
    class(linked_list_type__${dtype}$), intent(inout) :: this
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
    type(linked_list_type__${dtype}$), pointer, intent(in) :: src
    type(linked_list_type__${dtype}$), pointer, intent(out) :: dest

    if (.not. src%is_empty()) then
      stop 'Not implemented yet.'
    endif
    dest%block_size = src%block_size
  end subroutine assign_linked_list

  function is_empty(this)
    class(linked_list_type__${dtype}$), intent(inout) :: this
    logical :: is_empty

    is_empty = (this%n == 0)
  end function is_empty

  function get_length(this) result(length)
    class(linked_list_type__${dtype}$), intent(inout) :: this
    integer :: length

    length = this%n
  end function get_length

end module linked_list_module__${dtype}$
#:endfor
