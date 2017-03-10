module spin_det_module
  
  use constants_module
  use types_module

  implicit none

  private

  public :: spin_det_type
  public :: new_spin_det
  public :: new_spin_det_arr
  public :: delete
  public :: assignment(=)
  public :: operator(==)
  public :: operator(<)
  public :: operator(>)

  integer :: spin_det_cnt = 0

  type spin_det_type
    private
    integer(LONG), allocatable :: trunks(:)
    integer :: n_trunks = 0
    integer :: n_elec_cache = -1
    integer, allocatable :: orbitals_cache(:)
    contains
      procedure, public :: get_n_elec
      procedure, public :: get_elec_orbitals
      procedure, public :: set_orbital
      procedure, public :: get_orbital
      procedure, public :: print
      procedure, public :: from_eor ! eor without allocation for performance.
      procedure, public :: get_n_diff_orbitals
      procedure, public :: is_empty
      procedure, public :: get_hash
      procedure :: resize
      procedure :: destroy_orbitals_cache
      procedure :: trail_zero
  end type spin_det_type

  interface assignment(=)
    module procedure assign_spin_det
  end interface

  interface operator(==)
    module procedure equal_spin_det
  end interface

  interface operator(<)
    module procedure lt_spin_det
  end interface

  interface operator(>)
    module procedure gt_spin_det
  end interface

  interface new_spin_det
    module procedure new_spin_det_clone
    module procedure new_spin_det_by_n_orbs
    module procedure new_spin_det_by_trunks
  end interface new_spin_det

  interface delete
    module procedure delete_spin_det
    module procedure delete_spin_det_arr
  end interface delete

  contains

  function new_spin_det_by_n_orbs(n_orbs) result(res)
    integer, intent(in) :: n_orbs
    type(spin_det_type), pointer :: res
    integer :: n_trunks

    n_trunks = ceiling(n_orbs * 1.0 / C%TRUNK_SIZE)
    allocate(res)
    allocate(res%trunks(n_trunks))
    res%trunks(:) = 0
    res%n_trunks = n_trunks
    spin_det_cnt = spin_det_cnt + 1
  end function new_spin_det_by_n_orbs
 
  function new_spin_det_clone(src) result(res)
    type(spin_det_type), pointer, intent(in) :: src
    type(spin_det_type), pointer :: res
    
    allocate(res)
    if (allocated(src%trunks)) then
      allocate(res%trunks(src%n_trunks))
      res%trunks(:) = src%trunks(:)
      res%n_trunks = src%n_trunks
      res%n_elec_cache = src%n_elec_cache
    endif
    spin_det_cnt = spin_det_cnt + 1
  end function new_spin_det_clone

  function new_spin_det_by_trunks(trunks) result(res)
    integer(LONG), intent(in) :: trunks(:)
    type(spin_det_type), pointer :: res
    integer :: n_trunks

    n_trunks = size(trunks)
    allocate(res)
    allocate(res%trunks(n_trunks))
    res%trunks(:) = trunks(:)
    res%n_trunks = n_trunks
    spin_det_cnt = spin_det_cnt + 1
  end function new_spin_det_by_trunks

  function new_spin_det_arr(n) result(res)
    integer, intent(in) :: n
    type(spin_det_type), pointer :: res(:)

    allocate(res(n))
    spin_det_cnt = spin_det_cnt + n
  end function new_spin_det_arr

  subroutine delete_spin_det(spin_det)
    type(spin_det_type), pointer, intent(inout) :: spin_det

    if (.not. associated(spin_det)) return
    spin_det_cnt = spin_det_cnt - 1
    deallocate(spin_det)
    nullify(spin_det)
  end subroutine delete_spin_det

  subroutine delete_spin_det_arr(spin_det_arr)
    type(spin_det_type), pointer, intent(inout) :: spin_det_arr(:)

    if (.not. associated(spin_det_arr)) return
    spin_det_cnt = spin_det_cnt - size(spin_det_arr)
    deallocate(spin_det_arr)
    nullify(spin_det_arr)
  end subroutine delete_spin_det_arr

  subroutine destroy_orbitals_cache(this)
    class(spin_det_type), intent(inout) :: this

    if (allocated(this%orbitals_cache)) then
      deallocate(this%orbitals_cache)
    endif
  end subroutine destroy_orbitals_cache

  function get_orbital(this, orbital_idx) result(is_occupied)
    class(spin_det_type), intent(in) :: this
    integer, intent(in) :: orbital_idx
    logical :: is_occupied
    integer :: i
    integer :: trunk

    if (orbital_idx > this%n_trunks * C%TRUNK_SIZE) then
      call backtrace
      stop 'get_orbital with index out of range.'
    endif

    i = orbital_idx
    trunk = 1
    do while (i > C%TRUNK_SIZE)
      i = i - C%TRUNK_SIZE
      trunk = trunk + 1
    enddo
    is_occupied = btest(this%trunks(trunk), i - 1)
  end function get_orbital
  
  subroutine set_orbital(this, orbital_idx, is_occupied)
    class(spin_det_type), intent(inout) :: this
    integer, intent(in) :: orbital_idx
    logical, optional, intent(in) :: is_occupied ! Default to 1.
    integer :: i
    integer :: trunk

    if (orbital_idx > this%n_trunks * C%TRUNK_SIZE) then
      call backtrace
      stop 'set_orbital with index out of range.'
    endif

    i = orbital_idx
    trunk = 1
    do while (i > C%TRUNK_SIZE)
      i = i - C%TRUNK_SIZE
      trunk = trunk + 1
    enddo
    if ((.not. present(is_occupied)) .or. is_occupied) then
      this%trunks(trunk) = ibset(this%trunks(trunk), i - 1)
    else
      this%trunks(trunk) = ibclr(this%trunks(trunk), i - 1)
    endif
    this%n_elec_cache = -1
    call this%destroy_orbitals_cache()
  end subroutine set_orbital

  function get_n_elec(this) result(n_elec)
    class(spin_det_type), intent(inout) :: this
    integer :: n_elec
    integer :: i

    if (this%n_elec_cache >= 0) then
      n_elec = this%n_elec_cache
      return
    endif
    if (allocated(this%orbitals_cache)) then
      n_elec = size(this%orbitals_cache)
      this%n_elec_cache = n_elec
      return
    endif

    n_elec = 0
    do i = 1, this%n_trunks
      n_elec = n_elec + popcnt(this%trunks(i))
    enddo
    this%n_elec_cache = n_elec
  end function get_n_elec

  subroutine get_elec_orbitals(this, orbitals, n_elec_opt)
    class(spin_det_type), intent(inout) :: this
    integer, allocatable, intent(out) :: orbitals(:)
    integer, optional, intent(in) :: n_elec_opt
    integer :: n_elec
    integer :: i
    integer :: pos
    type(spin_det_type), target, save :: tmp_spin_det

    if (allocated(this%orbitals_cache)) then
      allocate(orbitals(size(this%orbitals_cache)))
      orbitals(:) = this%orbitals_cache(:)
      return
    endif

    if (present(n_elec_opt)) then
      n_elec = n_elec_opt
    else
      n_elec = this%get_n_elec()
    endif
    this%n_elec_cache = n_elec
    allocate(orbitals(n_elec))
    call tmp_spin_det%resize(this%n_trunks)
    tmp_spin_det%trunks(:) = this%trunks(:)
    do i = 1, n_elec
      pos = tmp_spin_det%trail_zero() + 1
      orbitals(i) = pos
      call tmp_spin_det%set_orbital(pos, .false.) 
    enddo

    allocate(this%orbitals_cache(n_elec))
    this%orbitals_cache(:) = orbitals(:)
  end subroutine get_elec_orbitals

  subroutine print(this)
    class(spin_det_type), intent(in) :: this
    integer :: i

    do i = 1, this%n_trunks
      write (6, '(I0, A, B0.63)') i, ': ', this%trunks(i) 
    enddo
  end subroutine print

  function trail_zero(this) result(cnt)
    class(spin_det_type), intent(in) :: this
    integer :: cnt
    integer :: i 

    cnt = 0
    do i = 1, this%n_trunks
      if (this%trunks(i) == 0) then
        cnt = cnt + C%TRUNK_SIZE
      else
        cnt = cnt + trailz(this%trunks(i))
        return
      endif
    enddo
  end function trail_zero

  subroutine resize(this, n_trunks)
    class(spin_det_type), intent(inout) :: this
    integer, intent(in) :: n_trunks

    if (this%n_trunks > 0 .and. this%n_trunks /= n_trunks) then
      deallocate(this%trunks)
    endif
    if (.not. allocated(this%trunks)) then
      allocate(this%trunks(n_trunks))
      this%n_trunks = n_trunks
    endif
    call this%destroy_orbitals_cache()
    this%n_elec_cache = -1
  end subroutine resize

  subroutine from_eor(this, op1, op2)
    class(spin_det_type), intent(inout) :: this
    type(spin_det_type), pointer, intent(in) :: op1, op2
    integer :: i

    if (op1%n_trunks /= op2%n_trunks) then
      call backtrace
      stop 'from eor input size mismatch.'
    endif

    call this%resize(op1%n_trunks)
    do i = 1, this%n_trunks
      this%trunks(i) = ieor(op1%trunks(i), op2%trunks(i))
    enddo
  end subroutine from_eor

  function get_n_diff_orbitals(this, op) result(n_diff)
    class(spin_det_type), intent(inout) :: this
    type(spin_det_type), pointer, intent(inout) :: op
    integer :: n_diff
    integer :: i
    integer(LONG) :: eor_trunk
    
    if (this%n_trunks /= op%n_trunks) then
      call backtrace
      stop 'get_n_diff_orbitals size mismatch.'
    endif

    n_diff = 0
    do i = 1, this%n_trunks
      eor_trunk = ieor(this%trunks(i), op%trunks(i))
      n_diff = n_diff + popcnt(eor_trunk)
    enddo
  end function get_n_diff_orbitals

  subroutine assign_spin_det(dest, src)
    type(spin_det_type), pointer, intent(in) :: src
    type(spin_det_type), pointer, intent(out) :: dest
   
    if (.not. associated(dest)) allocate(dest)
    call dest%resize(src%n_trunks)
    dest%trunks(:) = src%trunks(:)
    dest%n_elec_cache = src%n_elec_cache
    call dest%destroy_orbitals_cache()
  end subroutine assign_spin_det

  function equal_spin_det(left, right) result(is_equal)
    type(spin_det_type), pointer, intent(in) :: left, right
    logical :: is_equal
    integer :: i

    if (left%n_trunks /= right%n_trunks) then
      is_equal = .false.
      return
    endif
    do i = 1, left%n_trunks
      if (left%trunks(i) == right%trunks(i)) then
        cycle
      else
        is_equal = .false.
        return
      endif
    enddo
    is_equal = .true.
  end function equal_spin_det

  function lt_spin_det(left, right) result(is_lt)
    type(spin_det_type), pointer, intent(in) :: left, right
    logical :: is_lt
    integer :: i

    if (left%n_trunks /= right%n_trunks) then
      call backtrace
      stop 'spin det comparison size mismatch.'
    endif
    do i = left%n_trunks, 1, -1
      if (left%trunks(i) == right%trunks(i)) then
        cycle
      else
        is_lt = (left%trunks(i) < right%trunks(i))
        return
      endif
    enddo
    is_lt = .false.
  end function lt_spin_det

  function gt_spin_det(left, right) result(is_gt)
    type(spin_det_type), pointer, intent(in) :: left, right
    logical :: is_gt
    integer :: i

    if (left%n_trunks /= right%n_trunks) then
      call backtrace
      stop 'spin det comparison size mismatch.'
    endif
    do i = left%n_trunks, 1, -1
      if (left%trunks(i) == right%trunks(i)) then
        cycle
      else
        is_gt = (left%trunks(i) > right%trunks(i))
        return
      endif
    enddo
    is_gt = .false.
  end function gt_spin_det

  function is_empty(this) result(res)
    class(spin_det_type), intent(inout) :: this
    logical :: res

    res = (this%n_trunks == 0)
  end function is_empty

  function get_hash(this, table_size) result(hash_value)
    class(spin_det_type), intent(inout) :: this
    integer, intent(in) :: table_size
    integer :: hash_value
    integer :: i
    integer(LONG) :: sum_trunks

    sum_trunks = 0
    do i = 1, this%n_trunks
      sum_trunks = sum_trunks + this%trunks(i)
    enddo
    hash_value = mod(int(sum_trunks), table_size)
    if (hash_value <= 0) then
      hash_value = hash_value + table_size
    endif
  end function get_hash

end module spin_det_module
