module hash_table_module__spin_det__int_list
  
  use spin_det_module
  use linked_list_module__int

  implicit none

  private

  public :: hash_table_type__spin_det__int_list
  public :: build
  public :: delete

  type hash_table_type__spin_det__int_list
    integer :: table_size
    integer :: n = 0
    type(spin_det_type), pointer :: keys(:) => null()
    type(linked_list_type__int), pointer :: values(:) => null()
    contains
      procedure :: has
      procedure :: get
      procedure :: set
      procedure :: remove
  end type hash_table_type__spin_det__int_list

  interface build
    module procedure build_hash_table
  end interface build

  interface delete
    module procedure delete_hash_table
  end interface delete

  contains

  subroutine build_hash_table(this, table_size)
    integer, intent(in) :: table_size
    type(hash_table_type__spin_det__int_list), pointer :: this

    allocate(this)
    call build(this%keys, table_size)
    call build(this%values, table_size, 8)
    this%table_size = table_size
  end subroutine build_hash_table

  function has(this, key)
    class(hash_table_type__spin_det__int_list), intent(inout) :: this
    type(spin_det_type), pointer, intent(inout) :: key
    logical :: has
    type(linked_list_type__int), pointer :: val

    val => this%get(key)
    has = associated(val)
  end function has

  function get(this, key) result(val)
    class(hash_table_type__spin_det__int_list), intent(inout) :: this
    type(spin_det_type), pointer, intent(inout) :: key
    type(linked_list_type__int), pointer :: val
    type(spin_det_type), pointer :: tmp_key
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
    class(hash_table_type__spin_det__int_list), intent(inout) :: this
    type(spin_det_type), pointer, intent(inout) :: key
    type(linked_list_type__int), pointer, intent(in) :: val
    type(spin_det_type), pointer :: tmp_key
    type(linked_list_type__int), pointer :: tmp_val
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
    class(hash_table_type__spin_det__int_list), intent(inout) :: this
    type(spin_det_type), pointer, intent(in) :: key
    type(spin_det_type), pointer :: tmp_key_prev, tmp_key
    type(linked_list_type__int), pointer :: tmp_val_prev, tmp_val
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
    type(hash_table_type__spin_det__int_list), pointer, intent(inout) :: table

    if (.not. associated(table)) return
    call delete(table%keys)
    call delete(table%values)
    deallocate(table)
    nullify(table)
  end subroutine delete_hash_table

end module
module hash_table_module__det__det_list_node
  
  use det_module
  use doubly_linked_list_module

  implicit none

  private

  public :: hash_table_type__det__det_list_node
  public :: build
  public :: delete

  type hash_table_type__det__det_list_node
    integer :: table_size
    integer :: n = 0
    type(det_type), pointer :: keys(:) => null()
    type(doubly_linked_list_node_type), pointer :: values(:) => null()
    contains
      procedure :: has
      procedure :: get
      procedure :: set
      procedure :: remove
  end type hash_table_type__det__det_list_node

  interface build
    module procedure build_hash_table
  end interface build

  interface delete
    module procedure delete_hash_table
  end interface delete

  contains

  subroutine build_hash_table(this, table_size)
    integer, intent(in) :: table_size
    type(hash_table_type__det__det_list_node), pointer :: this

    allocate(this)
    call build(this%keys, table_size)
    call build(this%values, table_size)
    this%table_size = table_size
  end subroutine build_hash_table

  function has(this, key)
    class(hash_table_type__det__det_list_node), intent(inout) :: this
    type(det_type), pointer, intent(inout) :: key
    logical :: has
    type(doubly_linked_list_node_type), pointer :: val

    val => this%get(key)
    has = associated(val)
  end function has

  function get(this, key) result(val)
    class(hash_table_type__det__det_list_node), intent(inout) :: this
    type(det_type), pointer, intent(inout) :: key
    type(doubly_linked_list_node_type), pointer :: val
    type(det_type), pointer :: tmp_key
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
    class(hash_table_type__det__det_list_node), intent(inout) :: this
    type(det_type), pointer, intent(inout) :: key
    type(doubly_linked_list_node_type), pointer, intent(in) :: val
    type(det_type), pointer :: tmp_key
    type(doubly_linked_list_node_type), pointer :: tmp_val
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
    class(hash_table_type__det__det_list_node), intent(inout) :: this
    type(det_type), pointer, intent(in) :: key
    type(det_type), pointer :: tmp_key_prev, tmp_key
    type(doubly_linked_list_node_type), pointer :: tmp_val_prev, tmp_val
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
    type(hash_table_type__det__det_list_node), pointer, intent(inout) :: table

    if (.not. associated(table)) return
    call delete(table%keys)
    call delete(table%values)
    deallocate(table)
    nullify(table)
  end subroutine delete_hash_table

end module
