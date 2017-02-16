module dets_module

  use constants_module
  use det_module
  use types_module

  implicit none

  private

  public :: dets_type
  public :: assignment(=)

  type dets_type
    private
    integer, public :: n = 0 ! n <= 2nd rank of up/dn.
    logical, public :: is_sorted = .true. ! True if sorted in increasing order.
    type(det_type), allocatable, public :: dets(:)
    contains
      procedure, public :: reserve
      procedure, public :: append
      procedure, public :: initialize
      procedure, public :: get_det
      procedure, public :: set_det
      procedure, public :: sort ! Into ascending order.
      procedure, public :: merge_sorted_dets
  end type dets_type

  interface assignment(=)
    module procedure assign_dets
  end interface

  contains

  subroutine reserve(this, capacity)
    class(dets_type), intent(inout) :: this
    integer, intent(in) :: capacity
    if (allocated(this%dets) .and. size(this%dets) /= capacity) then
      deallocate(this%dets)
    end if
    if (.not. allocated(this%dets)) then
      allocate(this%dets(capacity))
    end if
  end subroutine reserve

  subroutine append(this, det)
    class(dets_type), intent(inout) :: this
    class(det_type), intent(in) :: det

    if (size(this%dets) == this%n) then
      stop 'Not enough capacity.'
    endif
    this%n = this%n + 1
    call this%set_det(this%n, det)
  end subroutine append

  subroutine initialize(this, n_orb, capacity)
    class(dets_type), intent(inout) :: this
    integer, intent(in) :: n_orb
    integer, intent(in) :: capacity
    call this%reserve(capacity)
    this%n = 0
  end subroutine initialize

  function get_det(this, idx) result(det)
    class(dets_type), intent(inout) :: this
    integer, intent(in) :: idx
    type(det_type) :: det
    
    if (idx > this%n) then
      call backtrace
      stop 'get_det with idx out of bound.'
    end if
    det = this%dets(idx)
  end function get_det

  subroutine set_det(this, idx, det)
    class(dets_type), intent(inout) :: this
    integer, intent(in) :: idx
    type(det_type), intent(in) :: det

    if (idx > this%n) stop 'set_det with idx out of bound.'
    this%dets(idx) = det
    if (this%n > 1) then
      this%is_sorted = .false.
    end if
  end subroutine set_det

  subroutine sort(this)
    class(dets_type), intent(inout) :: this
    integer :: i
    integer :: n
    integer, allocatable :: order(:)
    integer, allocatable :: tmp_order(:)
    type(det_type), allocatable :: sorted_dets(:)

    n = this%n
    if (n == 0) return
    allocate(order(n))
    allocate(tmp_order(n))
    do i = 1, n
      order(i) = i
    enddo

    call recur(1, n)

    allocate(sorted_dets(n))
    do i = 1, n
      sorted_dets(i) = this%dets(order(i))
    enddo
    deallocate(this%dets)
    call move_alloc(sorted_dets, this%dets)
    this%is_sorted = .true.

    contains

    recursive subroutine recur(left, right)
      integer, intent(in) :: left, right
      integer :: tmp
      integer :: mid
      integer :: ptr1, ptr2
      integer :: pos

      if (left == right) return
      if (left == right - 1) then
        if (this%dets(right) < this%dets(left)) then
          tmp = order(left)
          order(left) = order(right)
          order(right) = tmp
        endif
        return
      endif

      mid = (left + right) / 2
      call recur(left, mid)
      call recur(mid + 1, right)

      tmp_order(left: mid) = order(left: mid)
      ptr1 = left
      ptr2 = mid + 1
      pos = left
      do while (pos < right .and. ptr1 <= mid .and. ptr2 <= right)
        if (this%dets(tmp_order(ptr1)) < this%dets(order(ptr2))) then
          order(pos) = tmp_order(ptr1)
          ptr1 = ptr1 + 1
        else
          order(pos) = order(ptr2)
          ptr2 = ptr2 + 1
        endif
        pos = pos + 1
      enddo
      if (ptr1 <= mid) then
        order(pos: right) = tmp_order(ptr1: mid)
      endif
    end subroutine recur

  end subroutine sort

  subroutine merge_sorted_dets(this, dets, indices_old)
    ! Merge two sorted dets array into one sorted dets array.
    ! Assuming no duplication.
    class(dets_type), intent(inout) :: this
    type(dets_type), intent(inout) :: dets
    integer, allocatable, optional, intent(out) :: indices_old(:)
    integer :: n
    integer :: dets_idx
    integer :: n_merged_dets
    integer :: i
    integer :: idx
    integer, allocatable :: merged_dets_indices(:)
    integer :: ptr_this, ptr_dets
    type(det_type), allocatable :: new_dets(:)

    if (.not. dets%is_sorted) stop 'dets are not sorted.'
    if (.not. this%is_sorted) stop 'target dets are not sorted.'

    n = this%n
    if (dets%n == 0) return
    if (n == 0) then
      this = dets
      return
    endif
    if (present(indices_old)) allocate(indices_old(n))
    allocate(merged_dets_indices(n + dets%n)) ! Indices > n for input dets.
    n_merged_dets = 0
    ptr_this = 1
    ptr_dets = 1
  
    do while (ptr_this <= n .and. ptr_dets <= dets%n)
      n_merged_dets = n_merged_dets + 1
      if (this%dets(ptr_this) < dets%dets(ptr_dets)) then
        merged_dets_indices(n_merged_dets) = ptr_this
        if (present(indices_old)) indices_old(ptr_this) = n_merged_dets
        ptr_this = ptr_this + 1
      else if (this%dets(ptr_this) == dets%dets(ptr_dets)) then
        merged_dets_indices(n_merged_dets) = ptr_this
        if (present(indices_old)) indices_old(ptr_this) = n_merged_dets
        ptr_this = ptr_this + 1
        ptr_dets = ptr_dets + 1
      else
        merged_dets_indices(n_merged_dets) = ptr_dets + n
        ptr_dets = ptr_dets + 1
      endif
    enddo
    do while (ptr_this <= n)
      n_merged_dets = n_merged_dets + 1
      merged_dets_indices(n_merged_dets) = ptr_this
      if (present(indices_old)) indices_old(ptr_this) = n_merged_dets
      ptr_this = ptr_this + 1
    enddo
    do while (ptr_dets <= dets%n)
      n_merged_dets = n_merged_dets + 1
      merged_dets_indices(n_merged_dets) = ptr_dets + n
      ptr_dets = ptr_dets + 1
    enddo

    allocate(new_dets(n_merged_dets))
    do i = 1, n_merged_dets
      idx = merged_dets_indices(i)
      if (idx <= n) then
        new_dets(i) = this%dets(idx)
      else
        new_dets(i) = dets%dets(idx - n)
      endif
    enddo
    deallocate(this%dets)
    call move_alloc(new_dets, this%dets)
    this%n = n_merged_dets
  end subroutine merge_sorted_dets

  subroutine assign_dets(dest, src)
    class(dets_type), intent(in) :: src
    class(dets_type), intent(out) :: dest
    integer :: i

    call dest%reserve(src%n)
    dest%n = src%n
    do i = 1, src%n
      dest%dets(i) = src%dets(i)
    enddo
  end subroutine assign_dets

end module dets_module
