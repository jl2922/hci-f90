module sort_module

  use spin_det_module
  use det_module
  use types_module

  implicit none

  private

  ! Namespaces.
  public :: sort_type
  public :: sort


  type sort_type
    contains
      generic, public :: arg_sort => &
          & arg_sort__int, &
          & arg_sort__double, &
          & arg_sort__spin_det, &
          & arg_sort__det
      procedure :: arg_sort__int 
      procedure :: arg_sort__double 
      procedure :: arg_sort__spin_det 
      procedure :: arg_sort__det 
  end type sort_type

  type(sort_type) :: sort

  contains

  subroutine arg_sort__int(this, arr, order, n)
    class(sort_type), intent(in) :: this
    integer, intent(in) :: arr(:)
    integer, allocatable, intent(out) :: order(:)
    integer, intent(in) :: n
    integer :: i
    integer, allocatable :: tmp_order(:)

    if (n <= 0) then
      return
    endif
    if (allocated(order) .and. size(order) < n) then
      deallocate(order)
    endif
    if (.not. allocated(order)) then
      allocate(order(n))
    endif
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
        if (.not. arr(tmp_order(ptr1)) > arr(order(ptr2))) then
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

  end subroutine arg_sort__int
  subroutine arg_sort__double(this, arr, order, n)
    class(sort_type), intent(in) :: this
    real(DOUBLE), intent(in) :: arr(:)
    integer, allocatable, intent(out) :: order(:)
    integer, intent(in) :: n
    integer :: i
    integer, allocatable :: tmp_order(:)

    if (n <= 0) then
      return
    endif
    if (allocated(order) .and. size(order) < n) then
      deallocate(order)
    endif
    if (.not. allocated(order)) then
      allocate(order(n))
    endif
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
        if (.not. arr(tmp_order(ptr1)) > arr(order(ptr2))) then
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

  end subroutine arg_sort__double
  subroutine arg_sort__spin_det(this, arr, order, n)
    class(sort_type), intent(in) :: this
    type(spin_det_type), pointer, intent(in) :: arr(:)
    integer, allocatable, intent(out) :: order(:)
    integer, intent(in) :: n
    integer :: i
    integer, allocatable :: tmp_order(:)

    if (n <= 0) then
      return
    endif
    if (allocated(order) .and. size(order) < n) then
      deallocate(order)
    endif
    if (.not. allocated(order)) then
      allocate(order(n))
    endif
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
        if (.not. arr(tmp_order(ptr1)) > arr(order(ptr2))) then
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

  end subroutine arg_sort__spin_det
  subroutine arg_sort__det(this, arr, order, n)
    class(sort_type), intent(in) :: this
    type(det_type), pointer, intent(in) :: arr(:)
    integer, allocatable, intent(out) :: order(:)
    integer, intent(in) :: n
    integer :: i
    integer, allocatable :: tmp_order(:)

    if (n <= 0) then
      return
    endif
    if (allocated(order) .and. size(order) < n) then
      deallocate(order)
    endif
    if (.not. allocated(order)) then
      allocate(order(n))
    endif
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
        if (.not. arr(tmp_order(ptr1)) > arr(order(ptr2))) then
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

  end subroutine arg_sort__det

end module sort_module