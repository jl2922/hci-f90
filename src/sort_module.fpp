module sort_module

  use spin_det_module
  use det_module
  use types_module

  implicit none

  private

  ! Namespaces.
  public :: sort_type
  public :: sort

  #:set type_defs { &
      & 'int': 'integer', &
      & 'double': 'real(DOUBLE)', &
      & 'spin_det': 'type(spin_det_type), pointer', &
      & 'det': 'type(det_type), pointer'}
  #:set arg_sort_types ['int', 'double', 'spin_det', 'det']

  type sort_type
    contains
      generic, public :: arg_sort => &
          #:for arg_sort_type in arg_sort_types[0:-1]
          & arg_sort__${arg_sort_type}$, &
          #:endfor
          & arg_sort__${arg_sort_types[-1]}$
      #:for arg_sort_type in arg_sort_types
      procedure :: arg_sort__${arg_sort_type}$ 
      #:endfor
  end type sort_type

  type(sort_type) :: sort

  contains

  #:for arg_sort_type in arg_sort_types
  subroutine arg_sort__${arg_sort_type}$(this, arr, order, n)
    class(sort_type), intent(in) :: this
    ${type_defs[arg_sort_type]}$, intent(in) :: arr(:)
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

  end subroutine arg_sort__${arg_sort_type}$
  #:endfor

end module sort_module
