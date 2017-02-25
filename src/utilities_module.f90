module utilities_module

  use spin_det_module
  use types_module

  implicit none

  private

  public :: util ! Namespace for all utilility functions.

  type utilities_type
    contains
      procedure, public :: get_free_unit
      generic, public :: arg_sort => arg_sort_int, arg_sort_double
      procedure, public :: n_combinations
      procedure, public :: binary_search_lbound, binary_search_rbound
      procedure :: arg_sort_int, arg_sort_double
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
  
  subroutine arg_sort_double(this, arr, order, n)
    class(utilities_type), intent(in) :: this
    real(DOUBLE), intent(in) :: arr(:)
    integer, allocatable, intent(out) :: order(:)
    integer, intent(in) :: n
    integer :: i
    integer, allocatable :: tmp_order(:)

    if (.not. allocated(order)) then
      allocate(order(n))
    end if
    do i = 1, n
      order(i) = i
    end do
    allocate(tmp_order(n))

    call recur(1, n)

    contains

    recursive subroutine recur(left, right)
      integer, intent(in) :: left, right
      integer :: mid
      integer :: tmp
      integer :: pos
      integer :: ptr1, ptr2

      if (left == right) return
      if (left == right - 1) then
        if (arr(left) > arr(right)) then
          tmp = order(right)
          order(right) = order(left)
          order(left) = tmp
        end if
        return
      end if

      mid = (left + right) / 2
      call recur(left, mid)
      call recur(mid + 1, right)

      ! Merge
      tmp_order(left: mid) = order(left: mid)
      ptr1 = left
      ptr2 = mid + 1
      do pos = left, right
        if (arr(tmp_order(ptr1)) <= arr(order(ptr2))) then
          order(pos) = tmp_order(ptr1)
          ptr1 = ptr1 + 1
        else
          order(pos) = order(ptr2)
          ptr2 = ptr2 + 1
        end if
        if (ptr1 > mid .or. ptr2 > right) then
          exit
        end if
      end do
      if (ptr1 <= mid) then
        order(pos + 1: right) = tmp_order(ptr1: mid)
      end if
    end subroutine recur
  end subroutine arg_sort_double
  subroutine arg_sort_int(this, arr, order, n)
    class(utilities_type), intent(in) :: this
    integer, intent(in) :: arr(:)
    integer, allocatable, intent(out) :: order(:)
    integer, intent(in) :: n
    integer :: i
    integer, allocatable :: tmp_order(:)

    if (.not. allocated(order)) then
      allocate(order(n))
    end if
    do i = 1, n
      order(i) = i
    end do
    allocate(tmp_order(n))

    call recur(1, n)

    contains

    recursive subroutine recur(left, right)
      integer, intent(in) :: left, right
      integer :: mid
      integer :: tmp
      integer :: pos
      integer :: ptr1, ptr2

      if (left == right) return
      if (left == right - 1) then
        if (arr(left) > arr(right)) then
          tmp = order(right)
          order(right) = order(left)
          order(left) = tmp
        end if
        return
      end if

      mid = (left + right) / 2
      call recur(left, mid)
      call recur(mid + 1, right)

      ! Merge
      tmp_order(left: mid) = order(left: mid)
      ptr1 = left
      ptr2 = mid + 1
      do pos = left, right
        if (arr(tmp_order(ptr1)) <= arr(order(ptr2))) then
          order(pos) = tmp_order(ptr1)
          ptr1 = ptr1 + 1
        else
          order(pos) = order(ptr2)
          ptr2 = ptr2 + 1
        end if
        if (ptr1 > mid .or. ptr2 > right) then
          exit
        end if
      end do
      if (ptr1 <= mid) then
        order(pos + 1: right) = tmp_order(ptr1: mid)
      end if
    end subroutine recur
  end subroutine arg_sort_int

  function binary_search_lbound(this, val, arr) result(res)
    class(utilities_type), intent(in) :: this
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
    class(utilities_type), intent(in) :: this
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

  function n_combinations(this, n, k) result(res)
    class(utilities_type), intent(in) :: this
    integer, intent(in) :: n, k
    integer :: res
    integer :: i

    res = 1
    do i = n, n - k + 1
      res = res * i
    end do
    do i = k, 1, -1
      res = res / i
    end do
  end function n_combinations

end module utilities_module
