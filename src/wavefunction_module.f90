module wavefunction_module

  use constants_module
  use spin_det_module
  use det_module
  use types_module

  implicit none

  private

  public :: wavefunction_type
  public :: new_wavefunction
  public :: assignment(=)

  type wavefunction_type
    private
    integer, public :: n = 0 ! n <= 2nd rank of up/dn.
    logical, public :: is_sorted = .true. ! True if sorted in increasing order.
    type(det_type), pointer, public :: dets(:) => null()
    real(DOUBLE), allocatable, public :: coefs(:)
    contains
      procedure, public :: reserve_dets
      procedure, public :: append_det
      procedure, public :: initialize
      procedure, public :: get_det
      procedure, public :: set_det
      procedure, public :: sort_dets ! Into ascending order.
      procedure, public :: merge_sorted_dets
      procedure, public :: print
  end type wavefunction_type

  interface new_wavefunction
    module procedure new_wavefunction_default
    module procedure new_wavefunction_reserve
    module procedure new_wavefunction_clone
  end interface new_wavefunction

  interface assignment(=)
    module procedure assign_wavefunction
  end interface

  contains

  function new_wavefunction_default() result(wf)
    type(wavefunction_type), pointer :: wf

    allocate(wf)
  end function new_wavefunction_default

  function new_wavefunction_clone(src) result(wf)
    type(wavefunction_type), pointer, intent(in) :: src
    type(wavefunction_type), pointer :: wf

    allocate(wf)
    call assign_wavefunction(wf, src)
  end function new_wavefunction_clone

  function new_wavefunction_reserve(capacity) result(wf)
    integer, intent(in) :: capacity
    type(wavefunction_type), pointer :: wf

    allocate(wf)
    allocate(wf%dets(capacity))
    allocate(wf%coefs(capacity))
    wf%coefs(:) = 0
  end function new_wavefunction_reserve

  subroutine assign_wavefunction(dest, src)
    type(wavefunction_type), pointer, intent(out) :: dest
    type(wavefunction_type), pointer, intent(in) :: src
    integer :: i

    dest%n = src%n
    dest%is_sorted = src%is_sorted
    call dest%reserve_dets(src%n)
    do i = 1, src%n
      call dest%set_det(i, src%get_det(i))
      dest%coefs(i) = src%coefs(i)
    enddo
  end subroutine assign_wavefunction

  subroutine reserve_dets(this, capacity)
    class(wavefunction_type), intent(inout) :: this
    integer, intent(in) :: capacity

    if (associated(this%dets) .and. size(this%dets) /= capacity) then
      deallocate(this%dets)
      nullify(this%dets)
      deallocate(this%coefs)
    end if
    if (.not. associated(this%dets)) then
      allocate(this%dets(capacity))
      allocate(this%coefs(capacity))
    end if
  end subroutine reserve_dets

  subroutine append_det(this, det)
    class(wavefunction_type), intent(inout) :: this
    class(det_type), pointer, intent(in) :: det

    if (size(this%dets) == this%n) then
      stop 'Not enough capacity.'
    endif
    this%n = this%n + 1
    call this%set_det(this%n, det)
  end subroutine append_det

  subroutine initialize(this, n_orb, capacity)
    class(wavefunction_type), intent(inout) :: this
    integer, intent(in) :: n_orb
    integer, intent(in) :: capacity
    call this%reserve_dets(capacity)
    this%n = 0
  end subroutine initialize

  function get_det(this, idx) result(det)
    class(wavefunction_type), intent(inout) :: this
    integer, intent(in) :: idx
    type(det_type), pointer:: det
    
    if (idx > this%n) then
      call backtrace
      stop 'get_det with idx out of bound.'
    end if
    det => this%dets(idx)
  end function get_det

  subroutine set_det(this, idx, det)
    class(wavefunction_type), intent(inout) :: this
    integer, intent(in) :: idx
    type(det_type), pointer, intent(in) :: det
    type(det_type), pointer :: det_idx

    if (idx > this%n) stop 'set_det with idx out of bound.'
    det_idx => this%get_det(idx)
    det_idx%up => new_spin_det(det%up)
    det_idx%dn => new_spin_det(det%dn)
    if (this%n > 1) then
      this%is_sorted = .false.
    end if
  end subroutine set_det

  subroutine sort_dets(this)
    class(wavefunction_type), intent(inout) :: this
    integer :: i
    integer :: n
    integer, allocatable :: order(:)
    integer, allocatable :: tmp_order(:)
    type(det_type), pointer :: sorted_dets(:)

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
    this%dets => sorted_dets
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

  end subroutine sort_dets

  subroutine merge_sorted_dets(this, dets)
    ! Merge two sorted dets array into one sorted dets array.
    ! Coefs of the old dets keep unchanged, new dets 0.
    ! Assuming no duplication.
    class(wavefunction_type), intent(inout) :: this
    type(wavefunction_type), pointer, intent(inout) :: dets
    integer :: n
    integer :: dets_idx
    integer :: n_merged_dets
    integer :: i
    integer :: idx
    integer, allocatable :: merged_dets_indices(:)
    integer :: ptr_this, ptr_dets
    type(det_type), pointer :: new_dets(:)
    real(DOUBLE), allocatable :: new_coefs(:)

    if (.not. dets%is_sorted) stop 'dets are not sorted.'
    if (.not. this%is_sorted) stop 'target dets are not sorted.'

    n = this%n
    if (dets%n == 0) return
    allocate(merged_dets_indices(n + dets%n)) ! Indices > n for input dets.
    n_merged_dets = 0
    ptr_this = 1
    ptr_dets = 1
  
    do while (ptr_this <= n .and. ptr_dets <= dets%n)
      n_merged_dets = n_merged_dets + 1
      if (this%dets(ptr_this) < dets%dets(ptr_dets)) then
        merged_dets_indices(n_merged_dets) = ptr_this
        ptr_this = ptr_this + 1
      else if (this%dets(ptr_this) == dets%dets(ptr_dets)) then
        merged_dets_indices(n_merged_dets) = ptr_this
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
      ptr_this = ptr_this + 1
    enddo
    do while (ptr_dets <= dets%n)
      n_merged_dets = n_merged_dets + 1
      merged_dets_indices(n_merged_dets) = ptr_dets + n
      ptr_dets = ptr_dets + 1
    enddo

    allocate(new_dets(n_merged_dets))
    allocate(new_coefs(n_merged_dets))
    do i = 1, n_merged_dets
      idx = merged_dets_indices(i)
      if (idx <= n) then
        new_dets(i) = this%dets(idx)
        new_coefs(i) = this%coefs(idx)
      else
        new_dets(i) = dets%dets(idx - n)
        new_coefs(i) = 0.0_DOUBLE
      endif
    enddo
    if (associated(this%dets)) deallocate(this%dets)
    this%dets => new_dets
    call move_alloc(new_coefs, this%coefs)
    this%n = n_merged_dets
  end subroutine merge_sorted_dets

  subroutine assign_dets(dest, src)
    type(wavefunction_type), intent(in) :: src
    type(wavefunction_type), intent(out) :: dest
    integer :: i

    call dest%reserve_dets(src%n)
    dest%n = src%n
    do i = 1, src%n
      dest%dets(i) = src%dets(i)
    enddo
  end subroutine assign_dets

  subroutine print(this)
    class(wavefunction_type), intent(inout) :: this
    type(det_type), pointer :: tmp_det
    integer :: i

    do i = 1, this%n
      tmp_det => this%get_det(i)
      print *, '#', i, ' coefs = ', this%coefs(i)
      call tmp_det%print()
    enddo
  end subroutine print

end module wavefunction_module
