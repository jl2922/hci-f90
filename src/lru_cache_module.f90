module lru_cache_module

  use det_module
  use doubly_linked_list_module
  use hash_table_module__det__det_list_node

  implicit none

  private

  public :: lru_cache_type
  public :: build
  public :: delete

  type lru_cache_type
    integer :: capacity = 0
    type(doubly_linked_list_type), pointer :: list => null()
    type(hash_table_type__det__det_list_node), pointer :: map => null()
    contains
      procedure, public :: has
      procedure, public :: cache
  end type lru_cache_type

  interface build
    module procedure build_lru_cache
  end interface build

  interface delete
    module procedure delete_lru_cache
  end interface delete

  contains

  subroutine build_lru_cache(this, capacity)
    type(lru_cache_type), pointer, intent(inout) :: this
    integer, intent(in) :: capacity

    allocate(this)
    this%capacity = capacity
    call build(this%list)
    call build(this%map, capacity * 5)
  end subroutine build_lru_cache

  function has(this, key)
    class(lru_cache_type), intent(inout) :: this
    type(det_type), pointer, intent(inout) :: key
    logical :: has

    has = this%map%has(key)
  end function has

  subroutine cache(this, key)
    class(lru_cache_type), intent(inout) :: this
    type(det_type), pointer, intent(inout) :: key
    type(doubly_linked_list_node_type), pointer :: node
    type(det_type), pointer :: tmp_det

    if (this%capacity <= 0) then
      call backtrace
      stop 'Invalid LRU capacity'
    endif

    if (this%has(key)) then
      node => this%map%get(key)
      call this%list%remove(node)
      call this%list%push_front(node)
    else
      if (this%list%get_size() >= this%capacity) then
        node => this%list%back()
        call this%map%remove(node%item)
        call this%list%pop_back()
        call delete(node%item)
        call delete(node)
      endif
      call build(tmp_det, key)
      call this%list%push_front(tmp_det)
      node => this%list%front()
      call this%map%set(tmp_det, node)
    endif
  end subroutine cache

  subroutine delete_lru_cache(cache)
    type(lru_cache_type), pointer, intent(inout) :: cache

    call delete(cache%list)
    call delete(cache%map)
    deallocate(cache)
    nullify(cache)
  end subroutine delete_lru_cache

end module lru_cache_module
