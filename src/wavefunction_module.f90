module wavefunction_module

  use constants_module
  use spin_det_module
  use det_module
  use types_module
  use utilities_module

  implicit none

  private

  public :: wavefunction_type
  public :: new_wavefunction
  public :: assignment(=)
  public :: delete

  type wavefunction_type
    private
    integer, public :: n = 0 ! n <= 2nd rank of up/dn.
    integer, public :: n_up = 0
    logical, public :: is_sorted = .true. ! True if sorted in increasing order.
    logical, allocatable :: is_included(:)
    type(det_type), pointer, public :: dets(:) => null()
    real(DOUBLE), allocatable, public :: coefs(:)
    type(spin_det_type), pointer, public :: alpha_m1(:) => null()
    type(spin_det_type), pointer, public :: beta(:) => null()
    integer, allocatable, public :: alpha_m1_idx(:), beta_idx(:)
    contains
      procedure, public :: reserve_dets
      procedure, public :: append_det
      procedure, public :: find_potential_connections
      procedure, public :: find_potential_connections_setup
      procedure, public :: find_potential_connections_post_process
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

  interface delete
    module procedure delete_wavefunction
  end interface delete

  type(spin_det_type), target :: tmp_spin_det_for_find_potential_connections

  contains

  function new_wavefunction_default() result(wf)
    type(wavefunction_type), pointer :: wf

    allocate(wf)
  end function new_wavefunction_default

  function new_wavefunction_clone(src) result(wf)
    type(wavefunction_type), pointer, intent(in) :: src
    type(wavefunction_type), pointer :: wf

    allocate(wf)
    wf = src
  end function new_wavefunction_clone

  function new_wavefunction_reserve(capacity) result(wf)
    integer, intent(in) :: capacity
    type(wavefunction_type), pointer :: wf

    allocate(wf)
    allocate(wf%dets(capacity))
    allocate(wf%coefs(capacity))
    wf%coefs(:) = 0
  end function new_wavefunction_reserve

  subroutine delete_wavefunction(wf)
    type(wavefunction_type), pointer :: wf
    
    call delete(wf%dets)
    if (allocated(wf%coefs)) deallocate(wf%coefs)
    if (associated(wf%alpha_m1)) then
      call delete(wf%alpha_m1)
      call delete(wf%beta)
    endif
    deallocate(wf)
    nullify(wf)
  end subroutine delete_wavefunction

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
      call delete(this%dets)
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

  subroutine set_det(this, idx, src)
    class(wavefunction_type), intent(inout) :: this
    integer, intent(in) :: idx
    type(det_type), pointer, intent(in) :: src
    type(det_type), pointer :: det

    if (idx > this%n) stop 'set_det with idx out of bound.'
    det => this%get_det(idx)
    det%up => new_spin_det(src%up)
    det%dn => new_spin_det(src%dn)
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

  subroutine find_potential_connections_setup(this)
    class(wavefunction_type), intent(inout) :: this
    integer :: n
    integer :: n_up
    integer :: order_i, order_ii, order_ij
    integer :: i, j
    integer, allocatable :: alpha_m1_order(:), beta_order(:)
    integer, allocatable :: up_elec_orbitals(:)
    type(det_type), pointer :: det_i

    n = this%n
    n_up = this%n_up
    if (allocated(this%is_included)) deallocate(this%is_included)
    if (allocated(this%alpha_m1_idx)) deallocate(this%alpha_m1_idx)
    if (allocated(this%beta_idx)) deallocate(this%beta_idx)
    allocate(this%is_included(n))
    this%is_included = .false.
    this%beta => new_spin_det_arr(n)
    allocate(this%beta_idx(n))
    do i = 1, n
      det_i => this%get_det(i)
      this%beta(i) = det_i%dn
    enddo
    call util%arg_sort(this%beta, beta_order, n)
    do i = 1, n
      order_i = beta_order(i)
      det_i => this%get_det(order_i)
      this%beta(i) = det_i%dn
      this%beta_idx(i) = order_i 
    enddo

    det_i => this%get_det(1)
    n_up = det_i%up%get_n_elec()
    this%n_up = n_up
    this%alpha_m1 => new_spin_det_arr(n * n_up)
    allocate(this%alpha_m1_idx(n * n_up))
    allocate(up_elec_orbitals(n_up))
    do i = 1, n
      det_i => this%get_det(i)
      call det_i%up%get_elec_orbitals(up_elec_orbitals, n_up)
      do j = 1, n_up
        this%alpha_m1(n_up * (i - 1) + j) = det_i%up
        call this%alpha_m1(n_up * (i - 1) + j)%set_orbital( &
            & up_elec_orbitals(j), .false.)
      enddo
    enddo
    call util%arg_sort(this%alpha_m1, alpha_m1_order, n * n_up)
    do i = 1, n * n_up
      order_i = alpha_m1_order(i)
      order_ii = (order_i - 1) / n_up + 1
      order_ij = order_i - (order_ii - 1) * n_up
      det_i => this%get_det(order_ii)
      call det_i%up%get_elec_orbitals(up_elec_orbitals, n_up)
      this%alpha_m1(i) = det_i%up
      call this%alpha_m1(i)%set_orbital( &
          & up_elec_orbitals(order_ij), .false.)
      this%alpha_m1_idx(i) = order_ii 
    enddo

    write (6, '(A)') 'alpha_m1 and beta created.'
  end subroutine find_potential_connections_setup

  subroutine find_potential_connections( &
      & this, det, potential_connections, n_connections)
    class(wavefunction_type), intent(inout) :: this
    type(det_type), pointer, intent(in) :: det
    integer, allocatable, intent(out) :: potential_connections(:)
    integer, intent(out) :: n_connections
    integer :: i, j, idx 
    integer :: left, right
    integer :: n, n_up
    integer :: n_nonzero_connections
    integer :: n_diff
    integer, allocatable :: nonzero_connections(:)
    integer, allocatable :: up_elec_orbitals(:)
    integer, allocatable :: order(:)
    type(spin_det_type), pointer :: tmp_spin_det

    n_connections = 0
    n = this%n
    n_up = this%n_up
    allocate(nonzero_connections(n))
    allocate(order(n))
    if (.not. allocated(potential_connections)) then
      allocate(potential_connections(n))
    endif
    ! this%is_included = .false.

    ! beta.
    left = util%binary_search_lbound(det%dn, this%beta)
    right = util%binary_search_rbound(det%dn, this%beta)
    do i = left, right
      idx = this%beta_idx(i)
      if (.not. this%is_included(idx)) then
        n_connections = n_connections + 1
        potential_connections(n_connections) = idx
        this%is_included(idx) = .true.
      endif
    enddo

    ! alpha_m1.
    tmp_spin_det => tmp_spin_det_for_find_potential_connections
    tmp_spin_det = det%up
    call tmp_spin_det%get_elec_orbitals(up_elec_orbitals, n_up)
    do j = 1, n_up
      call tmp_spin_det%set_orbital(up_elec_orbitals(j), .false.)
      left = util%binary_search_lbound(tmp_spin_det, this%alpha_m1)
      right = util%binary_search_rbound(tmp_spin_det, this%alpha_m1)
      do i = left, right
        idx = this%alpha_m1_idx(i)
        if (.not. this%is_included(idx)) then
          n_connections = n_connections + 1
          potential_connections(n_connections) = idx
          this%is_included(idx) = .true.
        endif
      enddo
      call tmp_spin_det%set_orbital(up_elec_orbitals(j), .true.)
    enddo

    n_nonzero_connections = 0
    do i = 1, n_connections
      idx = potential_connections(i)
      this%is_included(idx) = .false.
      n_diff = det%get_n_diff_orbitals(this%get_det(idx))
      ! print *, idx, n_diff
      if (n_diff <= 4) then
        n_nonzero_connections = n_nonzero_connections + 1
        nonzero_connections(n_nonzero_connections) = idx
      endif
    enddo
    ! Sort in order of increasing indices (for deterministic pt).
    call util%arg_sort(nonzero_connections, order, n_nonzero_connections)
    do i = 1, n_nonzero_connections
      ! potential_connections(i) = nonzero_connections(i)
      potential_connections(i) = nonzero_connections(order(i))
      ! this%is_included(order(i)) = .false.
      ! potential_connections(i) = nonzero_connections(i)
      ! this%is_included(i) = .false.
    enddo
    ! call det%print()
    ! print *, potential_connections(1:n_nonzero_connections)
    ! stop
    n_connections = n_nonzero_connections
  end subroutine find_potential_connections

  subroutine find_potential_connections_post_process(this, j)
    class(wavefunction_type), intent(inout) :: this
    integer, intent(in) :: j

    this%is_included(j) = .false.
  end subroutine find_potential_connections_post_process

  subroutine print(this)
    class(wavefunction_type), intent(inout) :: this
    type(det_type), pointer :: tmp_det
    integer :: i

    do i = 1, this%n
      tmp_det => this%get_det(i)
      write (6, '(A, I0, A, F0.10)'), 'det #', i, ' coef = ', this%coefs(i)
      call tmp_det%print()
    enddo
  end subroutine print

end module wavefunction_module
