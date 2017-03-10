module search_module

  use spin_det_module

  implicit none

  private

  ! Namespaces.
  public :: search_type

  type search_type
    contains
      procedure, public :: binary_search_lbound
      procedure, public :: binary_search_rbound
  end type search_type

  contains

  function binary_search_lbound(this, val, arr) result(res)
    class(search_type), intent(in) :: this
    type(spin_det_type), pointer, intent(in) :: val
    type(spin_det_type), pointer, intent(in) :: arr(:)
    integer :: res
    integer :: left, right, mid

    left = 0
    right = size(arr)
    do while (right - left > 1)
      mid = (left + right) / 2
      if (arr(mid) < val) then
        left = mid
      else
        right = mid
      endif
    end do
    res = right
  end function binary_search_lbound

  function binary_search_rbound(this, val, arr) result(res)
    class(search_type), intent(in) :: this
    type(spin_det_type), pointer, intent(in) :: val
    type(spin_det_type), pointer, intent(in) :: arr(:)
    integer :: res
    integer :: left, right, mid

    left = 1
    right = size(arr) + 1
    do while (right - left > 1)
      mid = (left + right) / 2
      if (.not. arr(mid) > val) then
        left = mid
      else
        right = mid
      endif
    end do
    res = left
  end function binary_search_rbound

end module search_module
