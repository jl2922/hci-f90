module utilities_module

  use combinatorics_module
  use search_module
  use spin_det_module
  use sort_module
  use types_module

  implicit none

  private

  ! Namespace.
  public :: util

  type utilities_type
    type(sort_type) :: sort
    type(search_type) :: search
    type(combinatorics_type) :: combinatorics
    contains
      procedure, public :: get_free_unit
  end type utilities_type

  type(utilities_type) :: util

  contains

  type(integer) function get_free_unit(this) result(free_unit)
    class(utilities_type), intent(in) :: this
    logical :: is_open

    do free_unit = 100, 999
      inquire(unit=free_unit, opened=is_open)
      if (.not. is_open) then
        exit
      end if
    end do
  end function get_free_unit
  
end module utilities_module
