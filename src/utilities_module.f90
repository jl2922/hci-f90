module utilities_module

  use types_module

  implicit none

  public :: util

  ! Create a namespace for all the utility functions.
  type utilities_type
    contains
      procedure :: get_free_unit
      generic :: arg_sort => arg_sort_int, arg_sort_double
      procedure, private :: arg_sort_int, arg_sort_double
      procedure :: binary_search
      procedure :: ab_get
      procedure :: ab_set
      procedure :: ab_trailz
      procedure :: ab_eq
      procedure :: ab_lt
      procedure :: ab_gt
      procedure :: ab_eor
  end type utilities_type
  type(utilities_type) :: util

  contains

  type(integer) function get_free_unit(this) result(free_unit)
    class(utilities_type), intent(inout) :: this
    logical :: is_open
    do free_unit = 100, 999
      inquire(unit=free_unit, opened=is_open)
      if (.not. is_open) then
        exit
      end if
    end do
  end function get_free_unit
  
  subroutine arg_sort_double(this, arr, order, n)
    class(utilities_type), intent(inout) :: this
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
    class(utilities_type), intent(inout) :: this
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

  function ab_get(this, arr, pos) result(val)
    ! Test whether the bit at pos is set in array representation of bits.
    class(utilities_type), intent(inout) :: this
    integer(LONG), intent(in) :: arr(:)
    integer, intent(in) :: pos
    logical :: val
    integer :: i
    integer :: trunk

    if (pos > size(arr) * 64) then
      call backtrace
      stop 'ab_test with index out of range.'
    end if

    i = pos
    trunk = 1
    do while (i > 64)
      i = i - 64
      trunk = trunk + 1
    end do
    val = btest(arr(trunk), i - 1)
  end function ab_get

  subroutine ab_set(this, arr, pos, val)
    ! Set the bit at pos in array representation of bits to val.
    class(utilities_type), intent(inout) :: this
    integer(LONG), intent(inout) :: arr(:)
    integer, intent(in) :: pos
    logical, optional, intent(in) :: val
    integer :: i
    integer :: trunk

    if (pos > size(arr) * 64) then
      call backtrace
      stop 'ab_test with index out of range.'
    end if

    i = pos
    trunk = 1
    do while (i > 64)
      i = i - 64
      trunk = trunk + 1
    end do
    if ((.not. present(val)) .or. val) then
      arr(trunk) = ibset(arr(trunk), i - 1)
    else
      arr(trunk) = ibclr(arr(trunk), i - 1)
    end if
  end subroutine ab_set

  function ab_trailz(this, arr) result(cnt)
    class(utilities_type), intent(inout) :: this
    integer(LONG), intent(in) :: arr(:)
    integer :: cnt
    integer :: i

    cnt = 0
    do i = 1, size(arr)
      if (arr(i) == 0) then
        cnt = cnt + 64
      else
        cnt = cnt + trailz(arr(i))
        exit
      end if
    end do
  end function ab_trailz

  function ab_eq(this, arr1, arr2) result(is_equal)
    class(utilities_type), intent(inout) :: this
    integer(LONG), intent(in) :: arr1(:), arr2(:)
    logical :: is_equal
    integer :: size1, size2
    integer :: i

    size1 = size(arr1)
    size2 = size(arr2)
    if (size1 /= size2) then
      is_equal = .false.
      return
    end if

    do i = 1, size1
      if (arr1(i) /= arr2(i)) then
        is_equal = .false.
        return
      end if
    end do
    is_equal = .true.
  end function ab_eq

  function ab_lt(this, arr1, arr2) result(is_lt)
    class(utilities_type), intent(inout) :: this
    integer(LONG), intent(in) :: arr1(:), arr2(:)
    logical :: is_lt
    integer :: size1, size2
    integer :: i

    size1 = size(arr1)
    size2 = size(arr2)
    if (size1 /= size2) then
      stop 'ab_lt size mistmatch.'
    end if

    do i = size1, 1, -1
      if (arr1(i) < arr2(i)) then
        is_lt = .true.
        return
      else if (arr1(i) > arr2(i)) then
        is_lt = .false.
        return
      end if
    end do
    is_lt = .false. ! Equal.
  end function ab_lt

  function ab_gt(this, arr1, arr2) result(is_gt)
    class(utilities_type), intent(inout) :: this
    integer(LONG), intent(in) :: arr1(:), arr2(:)
    logical :: is_gt
    integer :: size1, size2
    integer :: i

    size1 = size(arr1)
    size2 = size(arr2)
    if (size1 /= size2) then
      stop 'ab_gt size mistmatch.'
    end if

    do i = size1, 1, -1
      if (arr1(i) > arr2(i)) then
        is_gt = .true.
        return
      else if (arr1(i) < arr2(i)) then
        is_gt = .false.
        return
      end if
    end do
    is_gt = .false. ! Equal.
  end function ab_gt

  subroutine ab_eor(this, arr1, arr2, arr_res)
    class(utilities_type), intent(inout) :: this
    integer(LONG), intent(in) :: arr1(:), arr2(:)
    integer(LONG), intent(out) :: arr_res(:)
    integer :: size1, size2, size_res
    integer :: i

    size1 = size(arr1)
    size2 = size(arr2)
    size_res = size(arr_res)
    if (size1 /= size2 .or. size1 /= size_res) then
      call backtrace
      stop 'ab_eor input size mismatch.'
    end if
    do i = 1, size1
      arr_res(i) = ieor(arr1(i), arr2(i))
    end do
  end subroutine ab_eor

end module utilities_module
