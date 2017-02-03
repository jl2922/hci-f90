module sort_class

  use types_class, only : DOUBLE

  implicit none

  private

  public :: merge_sort_arg

  interface merge_sort_arg
    module procedure merge_sort_arg_int
    module procedure merge_sort_arg_double
  end interface merge_sort_arg

  contains

    subroutine merge_sort_arg_int(arr, order, n)
      integer, intent(in) :: arr(:)
      integer, intent(out) :: order(:)
      integer, intent(in) :: n

      integer :: i
      integer, allocatable :: tmp_order(:)

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
          
          tmp_order(left: mid) = order(left: mid) 

          ! Merge
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
      
    end subroutine merge_sort_arg_int


    subroutine merge_sort_arg_double(arr, order, n)
      real(DOUBLE), intent(in) :: arr(:)
      integer, intent(out) :: order(:)
      integer, intent(in) :: n

      integer :: i
      integer, allocatable :: tmp_order(:)

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
            if (arr(order(left)) > arr(order(right))) then
              tmp = order(right)
              order(right) = order(left)
              order(left) = tmp
            end if
            return
          end if

          mid = (left + right) / 2

          call recur(left, mid)
          call recur(mid + 1, right)
          
          tmp_order(left: mid) = order(left: mid) 

          ! Merge
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
      
    end subroutine merge_sort_arg_double
end module sort_class
