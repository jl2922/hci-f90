module dets_module

  use det_module
  use types_module

  implicit none

  private

  public :: dets_type

  type dets_type
    private
    integer, public :: n = 0
    integer, public :: det_size = -1
    integer(LONG), allocatable, public :: up(:, :)
    integer(LONG), allocatable, public :: dn(:, :)
    logical, public :: is_sorted = .true. ! True if sorted in increasing order.
    contains
      procedure, public :: reserve
      procedure, public :: append
      procedure, public :: initialize
      procedure, public :: get_det
      procedure, public :: set_det
      procedure, public :: sort ! Into ascending order.
      procedure, public :: filter_by_sorted_dets
      procedure, public :: merge_sorted_dets
  end type dets_type

  contains

  subroutine reserve(this, n)
    class(dets_type), intent(inout) :: this
    integer, intent(in) :: n

    if (allocated(this%up)) deallocate(this%up)
    if (allocated(this%dn)) deallocate(this%dn)
    allocate(this%up(this%det_size, n))
    allocate(this%dn(this%det_size, n))
  end subroutine reserve

  subroutine append(this, det)
    class(dets_type), intent(inout) :: this
    class(det_type), intent(in) :: det

    if (size(this%up, 2) == this%n .or. size(this%dn, 2) == this%n) then
      stop 'Not enough capacity.'
    end if
    
    this%n = this%n + 1
    this%up(:, this%n) = det%up(:)
    this%dn(:, this%n) = det%dn(:)
  end subroutine append

  subroutine initialize(this, det_size, n)
    class(dets_type), intent(inout) :: this
    integer, intent(in) :: det_size
    integer, intent(in) :: n

    if (n /= this%n .or. det_size /= this%det_size) then
      if (allocated(this%up)) deallocate(this%up)
      if (allocated(this%dn)) deallocate(this%dn)
      allocate(this%up(det_size, n))
      allocate(this%dn(det_size, n))
      this%n = n
      this%det_size = det_size
    end if

    if (n > 0) then
      this%up = 0
      this%dn = 0
    end if
    if (n == 1) then
      this%is_sorted = .true.
    end if
  end subroutine initialize

  function get_det(this, idx) result(det)
    class(dets_type), intent(in) :: this
    integer, intent(in) :: idx
    type(det_type) :: det
    allocate(det%up(this%det_size))
    allocate(det%dn(this%det_size))
    det%up = this%up(:, idx)
    det%dn = this%dn(:, idx)
  end function get_det

  subroutine set_det(this, idx, det)
    class(dets_type), intent(inout) :: this
    integer, intent(in) :: idx
    type(det_type), intent(in) :: det
    this%up(:, idx) = det%up(:)
    this%dn(:, idx) = det%dn(:)
  end subroutine set_det

  subroutine sort(this)
    class(dets_type), intent(inout) :: this
    integer :: i
    integer :: n
    integer(LONG), allocatable :: sorted_up(:, :)
    integer(LONG), allocatable :: sorted_dn(:, :)
    integer, allocatable :: order(:)
    integer, allocatable :: tmp_order(:)

    n = this%n
    if (n == 0) return
    allocate(order(n))
    allocate(tmp_order(n))
    do i = 1, n
      order(i) = i
    end do

    call recur(1, n)

    allocate(sorted_up(this%det_size, n))
    allocate(sorted_dn(this%det_size, n))
    do i = 1, n
      sorted_up(:, i) = this%up(:, order(i))
      sorted_dn(:, i) = this%dn(:, order(i))
    end do
    deallocate(this%up)
    deallocate(this%dn)
    call move_alloc(sorted_up, this%up)
    call move_alloc(sorted_dn, this%dn)
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
        if (this%get_det(right) < this%get_det(left)) then
          tmp = order(left)
          order(left) = order(right)
          order(right) = tmp
        end if
        return
      end if

      mid = (left + right) / 2
      call recur(left, mid)
      call recur(mid + 1, right)

      tmp_order(left: mid) = order(left: mid)
      ptr1 = left
      ptr2 = mid + 1
      pos = left
      do while (pos < right .and. ptr1 <= mid .and. ptr2 <= right)
        if (this%get_det(tmp_order(ptr1)) < this%get_det(order(ptr2))) then
          order(pos) = tmp_order(ptr1)
          ptr1 = ptr1 + 1
        else
          order(pos) = order(ptr2)
          ptr2 = ptr2 + 1
        end if
        pos = pos + 1
      end do
      if (ptr1 <= mid) then
        order(pos: right) = tmp_order(ptr1: mid)
      end if
    end subroutine recur

  end subroutine sort

  subroutine filter_by_sorted_dets(this, dets)
    ! Filter out the dets in this sorted array if they exist in the sorted dets passed in.
    ! Assuming no duplication.
    class(dets_type), intent(inout) :: this
    type(dets_type), intent(in) :: dets
    integer :: det_size
    integer :: n_left_dets
    integer :: dets_idx
    integer :: i
    integer :: n
    integer, allocatable :: left_dets_indices(:)
    integer(LONG), allocatable :: new_up(:, :), new_dn(:, :)

    if (.not. dets%is_sorted) stop 'dets are not sorted.'
    if (.not. this%is_sorted) stop 'target dets are not sorted.'
    
    n = this%n
    allocate(left_dets_indices(n))
    n_left_dets = 0
    dets_idx = 1
    do i = 1, n
      do while (dets%get_det(dets_idx) < this%get_det(i) .and. dets_idx < dets%n)
        dets_idx = dets_idx + 1
      end do
      if (dets%get_det(dets_idx) == this%get_det(i)) then
        dets_idx = dets_idx + 1
        cycle
      end if
      n_left_dets = n_left_dets + 1
      left_dets_indices(n_left_dets) = i
    end do

    det_size = this%det_size
    allocate(new_up(det_size, n_left_dets))
    allocate(new_dn(det_size, n_left_dets))
    do i = 1, n_left_dets
      new_up(:, i) = this%up(:, left_dets_indices(i))
      new_dn(:, i) = this%dn(:, left_dets_indices(i))
    end do
    deallocate(this%up)
    deallocate(this%dn)
    call move_alloc(new_up, this%up)
    call move_alloc(new_dn, this%dn)
    this%n = n_left_dets
  end subroutine filter_by_sorted_dets

  subroutine merge_sorted_dets(this, dets)
    ! Merge two sorted dets array into one sorted dets array.
    ! Assuming no duplication.
    class(dets_type), intent(inout) :: this
    type(dets_type), intent(in) :: dets
    integer :: n
    integer :: dets_idx
    integer :: det_size
    integer :: n_merged_dets
    integer :: i
    integer :: idx
    integer, allocatable :: merged_dets_indices(:)
    integer :: ptr_this, ptr_dets
    integer(LONG), allocatable :: new_up(:, :), new_dn(:, :)

    if (.not. dets%is_sorted) stop 'dets are not sorted.'
    if (.not. this%is_sorted) stop 'target dets are not sorted.'

    n = this%n
    if (dets%n == 0) return
    if (n == 0) this%det_size = dets%det_size
    allocate(merged_dets_indices(n + dets%n)) ! Indices > n comes from dets passed in.
    n_merged_dets = 0
    ptr_this = 1
    ptr_dets = 1
  
    do while (ptr_this <= n .and. ptr_dets <= dets%n)
      n_merged_dets = n_merged_dets + 1
      if (this%get_det(ptr_this) < dets%get_det(ptr_dets)) then
        merged_dets_indices(n_merged_dets) = ptr_this
        ptr_this = ptr_this + 1
      else if (this%get_det(ptr_this) == dets%get_det(ptr_dets)) then
        merged_dets_indices(n_merged_dets) = ptr_this
        ptr_this = ptr_this + 1
        ptr_dets = ptr_dets + 1
      else
        merged_dets_indices(n_merged_dets) = ptr_dets + n
        ptr_dets = ptr_dets + 1
      end if
    end do
    do while (ptr_this <= n)
      n_merged_dets = n_merged_dets + 1
      merged_dets_indices(n_merged_dets) = ptr_this
      ptr_this = ptr_this + 1
    end do
    do while (ptr_dets <= dets%n)
      n_merged_dets = n_merged_dets + 1
      merged_dets_indices(n_merged_dets) = ptr_dets + n
      ptr_dets = ptr_dets + 1
    end do

    det_size = this%det_size
    allocate(new_up(det_size, n_merged_dets))
    allocate(new_dn(det_size, n_merged_dets))
    do i = 1, n_merged_dets
      idx = merged_dets_indices(i)
      if (idx <= n) then
        new_up(:, i) = this%up(:, idx)
        new_dn(:, i) = this%dn(:, idx)
      else
        new_up(:, i) = dets%up(:, idx - n)
        new_dn(:, i) = dets%dn(:, idx - n)
      end if
    end do
    if(allocated(this%up)) deallocate(this%up)
    if(allocated(this%dn)) deallocate(this%dn)
    call move_alloc(new_up, this%up)
    call move_alloc(new_dn, this%dn)
    this%n = n_merged_dets
  end subroutine merge_sorted_dets

end module dets_module
