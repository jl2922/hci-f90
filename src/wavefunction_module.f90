module wavefunction_module

  use constants_module
  use spin_det_module
  use hash_table_module__spin_det__int_list
  use det_module
  use linked_list_module__int
  use lru_cache_module
  use types_module
  use utilities_module

  implicit none

  private

  public :: wavefunction_type
  public :: build
  public :: delete
  public :: assignment(=)

  type ab_find_type
    type(spin_det_type), pointer :: alpha(:) => null()
    type(spin_det_type), pointer :: beta(:) => null()
    type(hash_table_type__spin_det__int_list), pointer :: alpha_m1_lut => null()
    type(hash_table_type__spin_det__int_list), pointer :: beta_m1_lut => null()
    integer, allocatable :: alpha_idx(:)
    integer, allocatable :: beta_idx(:)
    logical, allocatable :: is_alpha_m1_connected(:)
    logical, allocatable :: is_included(:)
    type(lru_cache_type), pointer :: lru => null()
    contains
      procedure, public :: clear => clear_ab_find
  end type ab_find_type

  type wavefunction_type
    integer, public :: n = 0
    logical, public :: is_sorted = .true. ! True if sorted in increasing order.
    real(DOUBLE), allocatable :: coefs(:)
    type(det_type), pointer :: dets(:) => null()
    integer :: n_up = 0
    integer :: n_dn = 0
    type(ab_find_type) :: ab_find
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

  interface build
    module procedure build_wavefunction_empty
    module procedure build_wavefunction_reserve
    module procedure build_wavefunction_clone
  end interface 

  interface delete
    module procedure delete_wavefunction
  end interface delete

  interface assignment(=)
    module procedure assign_wavefunction
  end interface

  contains

  subroutine build_wavefunction_empty(this)
    type(wavefunction_type), pointer, intent(inout) :: this

    allocate(this)
  end subroutine build_wavefunction_empty

  subroutine build_wavefunction_clone(this, src)
    type(wavefunction_type), pointer, intent(inout) :: this
    type(wavefunction_type), pointer, intent(in) :: src

    call build(this)
    this = src
  end subroutine build_wavefunction_clone

  subroutine build_wavefunction_reserve(this, capacity)
    type(wavefunction_type), pointer, intent(inout) :: this
    integer, intent(in) :: capacity

    call build(this)
    allocate(this%dets(capacity))
    allocate(this%coefs(capacity))
    this%coefs(:) = 0
  end subroutine build_wavefunction_reserve

  subroutine delete_wavefunction(this)
    type(wavefunction_type), pointer, intent(inout) :: this
    
    if (.not. associated(this)) return
    call this%ab_find%clear()
    call delete(this%dets)
    deallocate(this)
    nullify(this)
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
    call build(det%up, src%up)
    call build(det%dn, src%dn)
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
    integer :: order_i
    integer :: i, j
    integer, allocatable :: alpha_order(:)
    integer, allocatable :: beta_order(:)
    integer, allocatable :: up_elec_orbitals(:)
    integer, allocatable :: dn_elec_orbitals(:)
    type(det_type), pointer :: det_i
    type(det_type), pointer :: tmp_det
    type(spin_det_type), pointer :: tmp_spin_det
    type(linked_list_type__int), pointer :: idx_list
    type(linked_list_type__int), target, save :: EMPTY_INT_LIST

    n = this%n
    call this%ab_find%clear()
    allocate(this%ab_find%is_included(n))
    this%ab_find%is_included = .false.
    allocate(this%ab_find%is_alpha_m1_connected(n))
    this%ab_find%is_alpha_m1_connected = .false.
    call build(this%ab_find%lru, int(1e7))

    ! Setup alpha and beta.
    call build(this%ab_find%alpha, n)
    call build(this%ab_find%beta, n)
    allocate(this%ab_find%alpha_idx(n))
    allocate(this%ab_find%beta_idx(n))
    do i = 1, n
      det_i => this%get_det(i)
      tmp_spin_det => this%ab_find%alpha(i)
      tmp_spin_det = det_i%up
      tmp_spin_det => this%ab_find%beta(i)
      tmp_spin_det = det_i%dn
    enddo
    call util%sort%arg_sort(this%ab_find%beta, beta_order, n)
    call util%sort%arg_sort(this%ab_find%alpha, alpha_order, n)
    do i = 1, n
      order_i = alpha_order(i)
      det_i => this%get_det(order_i)
      tmp_spin_det => this%ab_find%alpha(i)
      tmp_spin_det = det_i%up
      this%ab_find%alpha_idx(i) = order_i 

      order_i = beta_order(i)
      det_i => this%get_det(order_i)
      tmp_spin_det => this%ab_find%beta(i)
      tmp_spin_det = det_i%dn
      this%ab_find%beta_idx(i) = order_i 
    enddo

    ! Setup alpha-m1.
    tmp_det => this%get_det(1)
    n_up = tmp_det%up%get_n_elec()
    this%n_up = n_up
    allocate(up_elec_orbitals(n_up))
    call build(this%ab_find%alpha_m1_lut, n * n_up)
    do i = 1, n
      det_i => this%get_det(i)
      tmp_spin_det => det_i%up
      call tmp_spin_det%get_elec_orbitals(up_elec_orbitals, n_up)
      do j = 1, n_up
        call tmp_spin_det%set_orbital(up_elec_orbitals(j), .false.)
        idx_list => this%ab_find%alpha_m1_lut%get(tmp_spin_det)
        if (.not. associated(idx_list)) then
          call this%ab_find%alpha_m1_lut%set(tmp_spin_det, EMPTY_INT_LIST)
          idx_list => this%ab_find%alpha_m1_lut%get(tmp_spin_det)
        endif
        call idx_list%append(i)
        call tmp_spin_det%set_orbital(up_elec_orbitals(j), .true.)
      enddo
    enddo

    ! Setup beta-m1.
    tmp_det => this%get_det(1)
    n_dn = tmp_det%dn%get_n_elec()
    this%n_dn = n_dn
    allocate(dn_elec_orbitals(n_dn))
    call build(this%ab_find%beta_m1_lut, n * n_up)
    do i = 1, n
      det_i => this%get_det(i)
      tmp_spin_det => det_i%dn
      call tmp_spin_det%get_elec_orbitals(dn_elec_orbitals, n_dn)
      do j = 1, n_dn
        call tmp_spin_det%set_orbital(dn_elec_orbitals(j), .false.)
        idx_list => this%ab_find%beta_m1_lut%get(tmp_spin_det)
        if (.not. associated(idx_list)) then
          call this%ab_find%beta_m1_lut%set(tmp_spin_det, EMPTY_INT_LIST)
          idx_list => this%ab_find%beta_m1_lut%get(tmp_spin_det)
        endif
        call idx_list%append(i)
        call tmp_spin_det%set_orbital(dn_elec_orbitals(j), .true.)
      enddo
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
    integer, allocatable :: tmp_connections(:)
    integer, allocatable :: up_elec_orbitals(:)
    integer, allocatable :: dn_elec_orbitals(:)
    integer, allocatable :: order(:)
    integer :: n_alpha_m1_connections
    integer, allocatable :: alpha_m1_connections(:)
    type(spin_det_type), pointer :: tmp_spin_det
    type(linked_list_type__int), pointer :: idx_list
    integer :: list_length

    n = this%n
    n_up = this%n_up
    n_dn = this%n_dn
    n_connections = 0
    allocate(tmp_connections(n))
    allocate(order(n))
    if (.not. allocated(potential_connections)) then
      allocate(potential_connections(n))
    endif

    ! alpha.
    left = util%search%binary_search_lbound(det%up, this%ab_find%alpha)
    right = util%search%binary_search_rbound(det%up, this%ab_find%alpha)
    do i = left, right
      idx = this%ab_find%alpha_idx(i)
      if (.not. this%ab_find%is_included(idx)) then
        n_connections = n_connections + 1
        tmp_connections(n_connections) = idx
        this%ab_find%is_included(idx) = .true.
      endif
    enddo

    ! beta.
    left = util%search%binary_search_lbound(det%dn, this%ab_find%beta)
    right = util%search%binary_search_rbound(det%dn, this%ab_find%beta)
    do i = left, right
      idx = this%ab_find%beta_idx(i)
      if (.not. this%ab_find%is_included(idx)) then
        n_connections = n_connections + 1
        tmp_connections(n_connections) = idx
        this%ab_find%is_included(idx) = .true.
      endif
    enddo

    ! alpha_m1.
    n_alpha_m1_connections = 0
    allocate(alpha_m1_connections(n))
    tmp_spin_det => det%up
    call tmp_spin_det%get_elec_orbitals(up_elec_orbitals, n_up)
    do j = 1, n_up
      call tmp_spin_det%set_orbital(up_elec_orbitals(j), .false.)
      idx_list => this%ab_find%alpha_m1_lut%get(tmp_spin_det)
      call tmp_spin_det%set_orbital(up_elec_orbitals(j), .true.)
      if (.not. associated(idx_list)) cycle
      call idx_list%begin()
      list_length = idx_list%get_length()
      do i = 1, list_length
        idx = idx_list%get()
        call idx_list%next()
        if (.not. this%ab_find%is_included(idx) .and. &
            & .not. this%ab_find%is_alpha_m1_connected(idx)) then
          n_alpha_m1_connections = n_alpha_m1_connections + 1
          alpha_m1_connections(n_alpha_m1_connections) = idx
          this%ab_find%is_alpha_m1_connected(idx) = .true.
        endif
      enddo
    enddo

    ! beta_m1.
    tmp_spin_det => det%dn
    call tmp_spin_det%get_elec_orbitals(dn_elec_orbitals, n_dn)
    do j = 1, n_dn
      call tmp_spin_det%set_orbital(dn_elec_orbitals(j), .false.)
      idx_list => this%ab_find%beta_m1_lut%get(tmp_spin_det)
      call tmp_spin_det%set_orbital(dn_elec_orbitals(j), .true.)
      if (.not. associated(idx_list)) cycle
      call idx_list%begin()
      list_length = idx_list%get_length()
      do i = 1, list_length
        idx = idx_list%get()
        call idx_list%next()
        if (this%ab_find%is_alpha_m1_connected(idx)) then
          n_connections = n_connections + 1
          tmp_connections(n_connections) = idx
          this%ab_find%is_included(idx) = .true.
        endif
      enddo
    enddo

    do i = 1, n_alpha_m1_connections
      idx = alpha_m1_connections(i)
      this%ab_find%is_alpha_m1_connected(idx) = .false.
    enddo
    do i = 1, n_connections
      idx = tmp_connections(i)
      this%ab_find%is_included(idx) = .false.
    enddo

    ! Sort in order of increasing indices (for deterministic pt).
    call util%sort%arg_sort(tmp_connections, order, n_connections)
    do i = 1, n_connections
      potential_connections(i) = tmp_connections(order(i))
    enddo
  end subroutine find_potential_connections

  subroutine print(this)
    class(wavefunction_type), intent(inout) :: this
    type(det_type), pointer :: det_ptr
    integer :: i

    do i = 1, this%n
      det_ptr => this%get_det(i)
      write (6, '(A, I0, A, F0.10)') 'det #', i, ', coef = ', this%coefs(i)
      call det_ptr%print()
    enddo
  end subroutine print

  subroutine clear_ab_find(this)
    class(ab_find_type), intent(inout) :: this

    call delete(this%alpha)
    call delete(this%beta)
    if (allocated(this%alpha_idx)) then
      deallocate(this%alpha_idx)
    endif
    if (allocated(this%beta_idx)) then
      deallocate(this%beta_idx)
    endif
    call delete(this%alpha_m1_lut)
    call delete(this%beta_m1_lut)
    if (allocated(this%is_alpha_m1_connected)) then
      deallocate(this%is_alpha_m1_connected)
    endif
    if (allocated(this%is_included)) then
      deallocate(this%is_included)
    endif
  end subroutine clear_ab_find

end module wavefunction_module
