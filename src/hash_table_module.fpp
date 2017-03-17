#:set dtypes ['spin_det__int_list', 'det__det_list_node']
#:set key_types { &
    & 'spin_det__int_list': 'type(spin_det_type), pointer', &
    & 'det__det_list_node': 'type(det_type), pointer'}
#:set val_types { &
    & 'spin_det__int_list': 'type(linked_list_type__int), pointer', &
    & 'det__det_list_node': 'type(doubly_linked_list_node_type), pointer'}
#:for dtype in dtypes
module hash_table_module__${dtype}$
  
#:if dtype == 'spin_det__int_list'
  use spin_det_module
  use linked_list_module__int
#:endif
#:if dtype == 'det__det_list_node'
  use det_module
  use doubly_linked_list_module
#:endif
#:set key_type key_types[dtype]
#:set val_type val_types[dtype]

  implicit none

  private

  public :: hash_table_type__${dtype}$
  public :: new_hash_table__${dtype}$
  public :: delete

  type hash_table_type__${dtype}$
    integer :: table_size
    integer :: n = 0
    ${key_type}$ :: keys(:) => null()
    ${val_type}$ :: values(:) => null()
    contains
      procedure :: has
      procedure :: get
      procedure :: set
      procedure :: remove
  end type hash_table_type__${dtype}$

  interface delete
    module procedure delete_hash_table
  end interface delete

  contains

  function new_hash_table__${dtype}$(table_size) result(table)
    integer, intent(in) :: table_size
    type(hash_table_type__${dtype}$), pointer :: table

    allocate(table)
#:if dtype == 'spin_det__int_list'
    table%keys => new_spin_det_arr(table_size)
    table%values => new_linked_list_arr__int(table_size, 8)
#:endif
#:if dtype == 'det__det_list_node'
    table%keys => new_det_arr(table_size)
    table%values => new_doubly_linked_list_node_arr(table_size)
#:endif
    table%table_size = table_size
  end function new_hash_table__${dtype}$

  function has(this, key)
    class(hash_table_type__${dtype}$), intent(inout) :: this
    ${key_type}$, intent(inout) :: key
    logical :: has
    ${val_type}$ :: val

    val => this%get(key)
    has = associated(val)
  end function has

  function get(this, key) result(val)
    class(hash_table_type__${dtype}$), intent(inout) :: this
    ${key_type}$, intent(inout) :: key
    ${val_type}$ :: val
    ${key_type}$ :: tmp_key
    integer :: cnt
    integer :: hash_value

    hash_value = key%get_hash(this%table_size) 
    cnt = 0
    do while (cnt < this%table_size)
      tmp_key => this%keys(hash_value)
      if (tmp_key%is_empty()) then
        exit
      endif
      if (tmp_key == key) then
        val => this%values(hash_value)
        return
      endif
      hash_value = hash_value + 1
      cnt = cnt + 1
      if (hash_value > this%table_size) then
        hash_value = 1
      endif
    enddo
    val => null()
    if (.not. associated(this%keys)) then
      call backtrace
      stop 'illegal state.'
    endif
  end function get

  subroutine set(this, key, val)
    class(hash_table_type__${dtype}$), intent(inout) :: this
    ${key_type}$, intent(inout) :: key
    ${val_type}$, intent(in) :: val
    ${key_type}$ :: tmp_key
    ${val_type}$ :: tmp_val
    integer :: cnt
    integer :: hash_value

    hash_value = key%get_hash(this%table_size)
    cnt = 0
    do while (cnt < this%table_size)
      tmp_key => this%keys(hash_value)
      if (tmp_key%is_empty()) then
        tmp_key = key
        tmp_val => this%values(hash_value)
        tmp_val = val
        return 
      endif
      if (tmp_key == key) then
        tmp_val => this%values(hash_value)
        tmp_val = val
        return
      endif
      hash_value = hash_value + 1
      cnt = cnt + 1
      if (hash_value > this%table_size) then
        hash_value = 1
      endif
    enddo
    stop 'hash table is full.'
  end subroutine set

  subroutine remove(this, key)
    class(hash_table_type__${dtype}$), intent(inout) :: this
    ${key_type}$, intent(in) :: key
    ${key_type}$ :: tmp_key_prev, tmp_key
    ${val_type}$ :: tmp_val_prev, tmp_val
    integer :: hash_value
    integer :: cnt

    hash_value = key%get_hash(this%table_size)
    cnt = 0
    do while (cnt < this%table_size)
      tmp_key => this%keys(hash_value)
      if (tmp_key%is_empty()) then
        exit
      endif
      if (tmp_key == key) then
        do
          tmp_val_prev => this%values(hash_value)
          tmp_key_prev => this%keys(hash_value)
          hash_value = hash_value + 1
          if (hash_value > this%table_size) then
            hash_value = 1
          endif
          tmp_key => this%keys(hash_value)
          if (tmp_key%is_empty()) then
            return
          endif
          if (tmp_key%get_hash(this%table_size) /= hash_value) then
            tmp_val => this%values(hash_value)
            tmp_val_prev = tmp_val
            tmp_key_prev = tmp_key
          endif
        enddo
      endif
      hash_value = hash_value + 1
      cnt = cnt + 1
      if (hash_value > this%table_size) then
        hash_value = 1
      endif
    enddo
    call backtrace
    stop 'Key not found'
  end subroutine remove

  subroutine delete_hash_table(table)
    type(hash_table_type__${dtype}$), pointer, intent(inout) :: table

    if (.not. associated(table)) return
    call delete(table%keys)
    call delete(table%values)
    deallocate(table)
    nullify(table)
  end subroutine delete_hash_table

end module
#:endfor
