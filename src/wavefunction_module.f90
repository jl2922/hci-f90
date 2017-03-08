module wavefunction_module

  use constants_module
  use spin_det_module
  use hash_table_module
  use det_module
  use linked_list_module__int
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
    integer, public :: n = 0
    logical, public :: is_sorted = .true. ! True if sorted in increasing order.
    real(DOUBLE), allocatable :: coefs(:)
    type(det_type), pointer :: dets(:) => null()
    type(spin_det_type), pointer :: alpha(:) => null()
    type(spin_det_type), pointer :: alpha_m1(:) => null()
    type(spin_det_type), pointer :: beta(:) => null()
    type(spin_det_type), pointer :: beta_m1(:) => null()
    integer, allocatable :: alpha_idx(:)
    integer, allocatable :: alpha_m1_idx(:)
    integer, allocatable :: beta_idx(:)
    integer, allocatable :: beta_m1_idx(:)
    logical, allocatable :: is_included(:)
    integer :: n_up = 0
    integer :: n_dn = 0
    type(hash_table_type), pointer :: alpha_m1_lut
    type(hash_table_type), pointer :: beta_m1_lut
    logical, allocatable :: is_alpha_m1_connected(:)
    contains
      procedure, public :: reserve_dets
      procedure, public :: append_det
      procedure, public :: find_potential_connections
      procedure, public :: find_potential_connections_setup
      procedure, public :: get_coef
      procedure, public :: get_det
      procedure, public :: set_coef
      procedure, public :: set_det
      procedure, public :: sort_dets ! Into ascending order.
      procedure, public :: merge_sorted_dets
      procedure, public :: print
  end type wavefunction_type

  interface new_wavefunction
    module procedure new_wavefunction_empty
    module procedure new_wavefunction_reserve
    module procedure new_wavefunction_clone
  end interface new_wavefunction

  interface assignment(=)
    module procedure assign_wavefunction
  end interface

  interface delete
    module procedure delete_wavefunction
  end interface delete

  contains

  function new_wavefunction_empty() result(wf)
    type(wavefunction_type), pointer :: wf

    allocate(wf)
  end function new_wavefunction_empty

  function new_wavefunction_clone(src) result(wf)
    type(wavefunction_type), pointer, intent(in) :: src
    type(wavefunction_type), pointer :: wf

    wf => new_wavefunction_empty()
    wf = src
  end function new_wavefunction_clone

  function new_wavefunction_reserve(capacity) result(wf)
    integer, intent(in) :: capacity
    type(wavefunction_type), pointer :: wf

    wf => new_wavefunction_empty()
    allocate(wf%dets(capacity))
    allocate(wf%coefs(capacity))
    wf%coefs(:) = 0
  end function new_wavefunction_reserve

  subroutine delete_wavefunction(wf)
    type(wavefunction_type), pointer :: wf
    
    if (.not. associated(wf)) return
    call delete(wf%dets)
    call delete(wf%alpha)
    call delete(wf%alpha_m1)
    call delete(wf%beta)
    call delete(wf%beta_m1)
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

  function get_coef(this, idx) result(coef)
    class(wavefunction_type), intent(inout) :: this
    integer, intent(in) :: idx
    real(DOUBLE) :: coef

    if (idx > this%n) then
      call backtrace
      stop 'get_coef with idx out of bound.'
    endif
    coef = this%coefs(idx)
  end function get_coef

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

  subroutine set_coef(this, idx, coef)
    class(wavefunction_type), intent(inout) :: this
    integer, intent(in) :: idx
    real(DOUBLE), intent(in) :: coef

    if (idx > this%n) then
      call backtrace
      stop 'set_coef with idx out of bound.'
    endif
    this%coefs(idx) = coef
  end subroutine set_coef

  subroutine sort_dets(this)
    class(wavefunction_type), intent(inout) :: this
    integer :: i
    integer :: n
    integer, allocatable :: order(:)
    type(det_type), pointer :: sorted_dets(:)
    real(DOUBLE), allocatable :: sorted_coefs(:)
    type(det_type), pointer :: det_ptr

    n = this%n
    if (n == 0) return
    call util%sort%arg_sort(this%dets, order, n)
    allocate(sorted_dets(n))
    allocate(sorted_coefs(n))
    do i = 1, n
      det_ptr => sorted_dets(i)
      det_ptr = this%get_det(order(i))
      sorted_coefs(i) = this%get_coef(i)
    enddo
    call delete(this%dets)
    nullify(det_ptr)
    this%dets => sorted_dets
    call move_alloc(sorted_coefs, this%coefs)
    this%is_sorted = .true.
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
    type(det_type), pointer :: det_ptr
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
      det_ptr => new_dets(i)
      if (idx <= n) then
        det_ptr = this%dets(idx)
        new_coefs(i) = this%coefs(idx)
      else
        det_ptr = dets%dets(idx - n)
        new_coefs(i) = 0.0_DOUBLE
      endif
    enddo
    call delete(this%dets)
    nullify(det_ptr)
    this%dets => new_dets
    call move_alloc(new_coefs, this%coefs)
    this%n = n_merged_dets
  end subroutine merge_sorted_dets

  subroutine find_potential_connections_setup(this)
    class(wavefunction_type), intent(inout) :: this
    integer :: n
    integer :: n_up
    integer :: n_dn
    integer :: order_i, order_ii, order_ij
    integer :: i, j
    integer, allocatable :: alpha_order(:), alpha_m1_order(:)
    integer, allocatable :: beta_order(:), beta_m1_order(:)
    integer, allocatable :: up_elec_orbitals(:)
    integer, allocatable :: dn_elec_orbitals(:)
    type(det_type), pointer :: det_i
    type(spin_det_type), pointer :: spin_det_ptr
    type(linked_list_type__int), pointer :: idx_list
    type(linked_list_type__int), target, save :: EMPTY_INT_LIST

    n = this%n
    if (allocated(this%is_included)) deallocate(this%is_included)
    if (allocated(this%alpha_idx)) deallocate(this%alpha_idx)
    if (allocated(this%alpha_m1_idx)) deallocate(this%alpha_m1_idx)
    if (allocated(this%beta_idx)) deallocate(this%beta_idx)
    if (allocated(this%beta_m1_idx)) deallocate(this%beta_m1_idx)
    if (allocated(this%is_alpha_m1_connected)) then
      deallocate(this%is_alpha_m1_connected)
    endif
    allocate(this%is_included(n))
    this%is_included = .false.

    ! Setup alpha and beta.
    this%alpha => new_spin_det_arr(n)
    this%beta => new_spin_det_arr(n)
    allocate(this%alpha_idx(n))
    allocate(this%beta_idx(n))
    do i = 1, n
      det_i => this%get_det(i)
      spin_det_ptr => this%alpha(i)
      spin_det_ptr = det_i%up
      spin_det_ptr => this%beta(i)
      spin_det_ptr = det_i%dn
    enddo
    call util%sort%arg_sort(this%beta, beta_order, n)
    call util%sort%arg_sort(this%alpha, alpha_order, n)
    do i = 1, n
      order_i = alpha_order(i)
      det_i => this%get_det(order_i)
      spin_det_ptr => this%alpha(i)
      spin_det_ptr = det_i%up
      this%alpha_idx(i) = order_i 

      order_i = beta_order(i)
      det_i => this%get_det(order_i)
      spin_det_ptr => this%beta(i)
      spin_det_ptr = det_i%dn
      this%beta_idx(i) = order_i 
    enddo

    ! Setup alpha-m1.
    det_i => this%get_det(1)
    n_up = det_i%up%get_n_elec()
    this%n_up = n_up
    this%alpha_m1 => new_spin_det_arr(n * n_up)
    allocate(this%alpha_m1_idx(n * n_up))
    allocate(up_elec_orbitals(n_up))
    this%alpha_m1_lut => new_hash_table(n * n_up)
    do i = 1, n
      det_i => this%get_det(i)
      call det_i%up%get_elec_orbitals(up_elec_orbitals, n_up)
      do j = 1, n_up
        spin_det_ptr => this%alpha_m1(n_up * (i - 1) + j)
        spin_det_ptr = det_i%up
        call spin_det_ptr%set_orbital(up_elec_orbitals(j), .false.)
        idx_list => this%alpha_m1_lut%get(spin_det_ptr)
        if (idx_list%is_empty()) then
          call this%alpha_m1_lut%set(spin_det_ptr, EMPTY_INT_LIST)
          idx_list => this%alpha_m1_lut%get(spin_det_ptr)
        endif
        call idx_list%append(i)
        ! print *, 'appending ', i, ' to:', spin_det_ptr%get_hash(n * n_up)
        ! call spin_det_ptr%print()
        ! print *, 'current size:', idx_list%get_length()
      enddo
    enddo
    call util%sort%arg_sort(this%alpha_m1, alpha_m1_order, n * n_up)
    do i = 1, n * n_up
      order_i = alpha_m1_order(i)
      order_ii = (order_i - 1) / n_up + 1
      order_ij = order_i - (order_ii - 1) * n_up
      det_i => this%get_det(order_ii)
      call det_i%up%get_elec_orbitals(up_elec_orbitals, n_up)
      spin_det_ptr => this%alpha_m1(i)
      spin_det_ptr = det_i%up
      call spin_det_ptr%set_orbital(up_elec_orbitals(order_ij), .false.)
      this%alpha_m1_idx(i) = order_ii 
    enddo
    allocate(this%is_alpha_m1_connected(n))
    this%is_alpha_m1_connected = .false.

    ! Setup beta-m1.
    det_i => this%get_det(1)
    n_dn = det_i%dn%get_n_elec()
    this%n_dn = n_dn
    this%beta_m1 => new_spin_det_arr(n * n_dn)
    allocate(this%beta_m1_idx(n * n_dn))
    allocate(dn_elec_orbitals(n_dn))
    this%beta_m1_lut => new_hash_table(n * n_up)
    do i = 1, n
      det_i => this%get_det(i)
      call det_i%dn%get_elec_orbitals(dn_elec_orbitals, n_dn)
      do j = 1, n_dn
        spin_det_ptr => this%beta_m1(n_dn * (i - 1) + j)
        spin_det_ptr = det_i%dn
        call spin_det_ptr%set_orbital(dn_elec_orbitals(j), .false.)
        idx_list => this%beta_m1_lut%get(spin_det_ptr)
        if (idx_list%is_empty()) then
          call this%beta_m1_lut%set(spin_det_ptr, EMPTY_INT_LIST)
          idx_list => this%beta_m1_lut%get(spin_det_ptr)
        endif
        call idx_list%append(i)
      enddo
    enddo
    call util%sort%arg_sort(this%beta_m1, beta_m1_order, n * n_dn)
    do i = 1, n * n_dn
      order_i = beta_m1_order(i)
      order_ii = (order_i - 1) / n_dn + 1
      order_ij = order_i - (order_ii - 1) * n_dn
      det_i => this%get_det(order_ii)
      call det_i%dn%get_elec_orbitals(dn_elec_orbitals, n_dn)
      spin_det_ptr => this%beta_m1(i)
      spin_det_ptr = det_i%dn
      call spin_det_ptr%set_orbital(dn_elec_orbitals(order_ij), .false.)
      this%beta_m1_idx(i) = order_ii 
    enddo

    write (6, '(A)') 'alpha and beta created.'
  end subroutine find_potential_connections_setup

  subroutine find_potential_connections( &
      & this, det, potential_connections, n_connections)
    class(wavefunction_type), intent(inout) :: this
    type(det_type), pointer, intent(in) :: det
    integer, allocatable, intent(out) :: potential_connections(:)
    integer, intent(out) :: n_connections
    integer :: i, j, idx 
    integer :: left, right
    integer :: n
    integer :: n_up
    integer :: n_dn
    integer :: n_tmp_connections
    integer :: n_diff
    integer, allocatable :: tmp_connections(:)
    integer, allocatable :: up_elec_orbitals(:)
    integer, allocatable :: dn_elec_orbitals(:)
    integer, allocatable :: order(:)
    integer :: n_alpha_m1_connections
    integer, allocatable :: alpha_m1_connections(:)
    type(spin_det_type), target, save :: tmp_spin_det
    type(spin_det_type), pointer :: spin_det_ptr
    type(linked_list_type__int), pointer :: idx_list
    integer :: list_length

    n_connections = 0
    n = this%n
    n_up = this%n_up
    n_dn = this%n_dn
    allocate(tmp_connections(n))
    allocate(order(n))
    if (.not. allocated(potential_connections)) then
      allocate(potential_connections(n))
    endif

    ! alpha.
    left = util%search%binary_search_lbound(det%up, this%alpha)
    right = util%search%binary_search_rbound(det%up, this%alpha)
    do i = left, right
      idx = this%alpha_idx(i)
      if (.not. this%is_included(idx)) then
        n_connections = n_connections + 1
        potential_connections(n_connections) = idx
        this%is_included(idx) = .true.
      endif
    enddo

    ! beta.
    left = util%search%binary_search_lbound(det%dn, this%beta)
    right = util%search%binary_search_rbound(det%dn, this%beta)
    do i = left, right
      idx = this%beta_idx(i)
      if (.not. this%is_included(idx)) then
        n_connections = n_connections + 1
        potential_connections(n_connections) = idx
        this%is_included(idx) = .true.
      endif
    enddo

    ! alpha_m1.
    n_alpha_m1_connections = 0
    allocate(alpha_m1_connections(n))
    tmp_spin_det = det%up
    call tmp_spin_det%get_elec_orbitals(up_elec_orbitals, n_up)
    spin_det_ptr => tmp_spin_det
    do j = 1, n_up
      call tmp_spin_det%set_orbital(up_elec_orbitals(j), .false.)
      idx_list => this%alpha_m1_lut%get(spin_det_ptr)
      call idx_list%begin()
      list_length = idx_list%get_length()
      do i = 1, list_length
        idx = idx_list%get()
        call idx_list%next()
        if (.not. this%is_included(idx) .and. &
            & .not. this%is_alpha_m1_connected(idx)) then
          n_alpha_m1_connections = n_alpha_m1_connections + 1
          alpha_m1_connections(n_alpha_m1_connections) = idx
          this%is_alpha_m1_connected(idx) = .true.
        endif
      enddo
      ! left = util%search%binary_search_lbound(tmp_spin_det, this%alpha_m1)
      ! right = util%search%binary_search_rbound(tmp_spin_det, this%alpha_m1)
      ! do i = left, right
      !   idx = this%alpha_m1_idx(i)
      !   if (.not. this%is_included(idx)) then
      !     n_alpha_m1_connections = n_alpha_m1_connections + 1
      !     alpha_m1_connections(n_alpha_m1_connections) = idx
      !     is_alpha_m1_connected(idx) = .true.
      !   endif
      ! enddo
      call tmp_spin_det%set_orbital(up_elec_orbitals(j), .true.)
    enddo

    ! beta_m1.
    tmp_spin_det = det%dn
    spin_det_ptr => tmp_spin_det
    call tmp_spin_det%get_elec_orbitals(dn_elec_orbitals, n_dn)
    do j = 1, n_dn
      call tmp_spin_det%set_orbital(dn_elec_orbitals(j), .false.)
      idx_list => this%beta_m1_lut%get(spin_det_ptr)
      call idx_list%begin()
      list_length = idx_list%get_length()
      do i = 1, list_length
        idx = idx_list%get()
        call idx_list%next()
        if (.not. this%is_included(idx) .and. &
            & this%is_alpha_m1_connected(idx)) then
          n_connections = n_connections + 1
          potential_connections(n_connections) = idx
          this%is_included(idx) = .true.
        endif
      enddo
      ! left = util%search%binary_search_lbound(tmp_spin_det, this%beta_m1)
      ! right = util%search%binary_search_rbound(tmp_spin_det, this%beta_m1)
      ! do i = left, right
      !   idx = this%beta_m1_idx(i)
      !   if (.not. this%is_included(idx) .and. this%is_alpha_m1_connected(idx)) then
      !     n_connections = n_connections + 1
      !     potential_connections(n_connections) = idx
      !     this%is_included(idx) = .true.
      !   endif
      ! enddo
      call tmp_spin_det%set_orbital(dn_elec_orbitals(j), .true.)
    enddo
    do i = 1, n_alpha_m1_connections
      idx = alpha_m1_connections(i)
      this%is_alpha_m1_connected(idx) = .false.
    enddo
    
    n_tmp_connections = 0
    do i = 1, n_connections
      idx = potential_connections(i)
      this%is_included(idx) = .false.
      n_tmp_connections = n_tmp_connections + 1
      tmp_connections(n_tmp_connections) = idx
    enddo

    ! Sort in order of increasing indices (for deterministic pt).
    call util%sort%arg_sort(tmp_connections, order, n_tmp_connections)
    do i = 1, n_tmp_connections
      potential_connections(i) = tmp_connections(order(i))
    enddo
    n_connections = n_tmp_connections
  end subroutine find_potential_connections

  subroutine print(this)
    class(wavefunction_type), intent(inout) :: this
    type(det_type), pointer :: det_ptr
    integer :: i

    do i = 1, this%n
      det_ptr => this%get_det(i)
      write (6, '(A, I0, A, F0.10)'), 'det #', i, ', coef = ', this%coefs(i)
      call det_ptr%print()
    enddo
  end subroutine print

end module wavefunction_module
