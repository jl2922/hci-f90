module hash_table_module
  
  use spin_det_module
  use linked_list_module__int
  implicit none

  private

  public :: hash_table_type
  public :: new_hash_table
  public :: delete

  type hash_table_type
    integer :: table_size
    integer :: n = 0
    type(spin_det_type), pointer :: keys(:) => null()
    type(linked_list_type__int), pointer :: values(:) => null()
    contains
      procedure :: get
      procedure :: set
  end type hash_table_type

  type(linked_list_type__int), target :: EMPTY_LIST

  interface delete
    module procedure delete_hash_table
  end interface delete

  contains

  function new_hash_table(table_size) result(table)
    integer, intent(in) :: table_size
    type(hash_table_type), pointer :: table

    allocate(table)
    table%keys => new_spin_det_arr(table_size)
    table%values => new_linked_list_arr__int(table_size, 8)
    table%table_size = table_size
  end function new_hash_table

  function get(this, spin_det) result(int_list)
    class(hash_table_type), intent(inout) :: this
    type(spin_det_type), pointer, intent(inout) :: spin_det
    type(linked_list_type__int), pointer :: int_list
    type(spin_det_type), pointer :: tmp_spin_det
    integer :: cnt
    integer :: hash_value

    hash_value = spin_det%get_hash(this%table_size) 
    cnt = 0
    do while (cnt < this%table_size)
      tmp_spin_det => this%keys(hash_value)
      if (tmp_spin_det%is_empty()) then
        exit
      endif
      if (tmp_spin_det == spin_det) then
        int_list => this%values(hash_value)
        return
      endif
      hash_value = hash_value + 1
      cnt = cnt + 1
      if (hash_value > this%table_size) then
        hash_value = 1
      endif
    enddo
    int_list => EMPTY_LIST
    if (.not. associated(this%keys)) then
      call backtrace
      stop 'illegal state.'
    endif
  end function get

  subroutine set(this, spin_det, int_list)
    class(hash_table_type), intent(inout) :: this
    type(spin_det_type), pointer, intent(inout) :: spin_det
    type(linked_list_type__int), pointer, intent(in) :: int_list
    type(spin_det_type), pointer :: tmp_spin_det
    type(linked_list_type__int), pointer :: tmp_int_list
    integer :: cnt
    integer :: hash_value

    hash_value = spin_det%get_hash(this%table_size)
    cnt = 0
    do while (cnt < this%table_size)
      tmp_spin_det => this%keys(hash_value)
      if (tmp_spin_det%is_empty()) then
        tmp_spin_det = spin_det
        tmp_int_list => this%values(hash_value)
        tmp_int_list = int_list
        return 
      endif
      if (tmp_spin_det == spin_det) then
        tmp_int_list => this%values(hash_value)
        tmp_int_list = int_list
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

  subroutine delete_hash_table(table)
    type(hash_table_type), pointer, intent(inout) :: table

    if (.not. associated(table)) return
    call delete(table%keys)
    call delete(table%values)
    deallocate(table)
    nullify(table)
  end subroutine delete_hash_table

end module
