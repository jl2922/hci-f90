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
  public :: operator(.eor.)
  public :: tmp_spin_det_instances

  integer :: spin_det_cnt = 0

  type spin_det_type
    private
    integer(LONG), allocatable, public :: trunks(:)
    integer :: n_trunks = 0
    integer, allocatable :: orbitals_cache(:)
    contains
      procedure, public :: get_n_elec
      procedure, public :: get_elec_orbitals
      procedure, public :: set_orbital
      procedure, public :: get_orbital
      procedure, public :: resize
      procedure, public :: print
      procedure, public :: from_eor ! eor without allocation for performance.
      procedure :: clean
      procedure :: destroy_orbitals_cache
      procedure :: trail_zero
  end type spin_det_type

  type(spin_det_type), target :: tmp_spin_det_instances(64)

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

  interface operator(.eor.)
    module procedure eor_spin_det
  end interface

  interface new_spin_det
    module procedure new_spin_det_clone
    module procedure new_spin_det_by_n_trunks
    module procedure new_spin_det_by_trunks
  end interface new_spin_det

  interface delete
    module procedure delete_spin_det
    module procedure delete_spin_det_arr
  end interface delete

  contains

  function new_spin_det_by_n_trunks(n_trunks) result(res)
    integer, intent(in) :: n_trunks
    type(spin_det_type), pointer :: res

    allocate(res)
    allocate(res%trunks(n_trunks))
    res%trunks(:) = 0
    res%n_trunks = n_trunks
    spin_det_cnt = spin_det_cnt + 1
  end function new_spin_det_by_n_trunks
 
  function new_spin_det_clone(src) result(res)
    type(spin_det_type), pointer, intent(in) :: src
    type(spin_det_type), pointer :: res

    allocate(res)
    if (allocated(src%trunks)) then
      allocate(res%trunks(src%n_trunks))
      res%trunks(:) = src%trunks(:)
      res%n_trunks = src%n_trunks
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

    call spin_det%clean()
    deallocate(spin_det)
    nullify(spin_det)
  end subroutine delete_spin_det

  subroutine delete_spin_det_arr(spin_det_arr)
    type(spin_det_type), pointer, intent(inout) :: spin_det_arr(:)
    integer :: i
    do i = 1, size(spin_det_arr)
      call spin_det_arr(i)%clean()
    enddo
    deallocate(spin_det_arr)
    nullify(spin_det_arr)
  end subroutine delete_spin_det_arr

  subroutine clean(this)
    class(spin_det_type), intent(inout) :: this

    deallocate(this%trunks)
    call this%destroy_orbitals_cache()
    spin_det_cnt = spin_det_cnt - 1
    ! print *, spin_det_cnt
    call flush(6)
  end subroutine clean
  
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
  
  subroutine set_orbital(this, orbital_idx, occupancy)
    class(spin_det_type), intent(inout) :: this
    integer, intent(in) :: orbital_idx
    logical, optional, intent(in) :: occupancy ! Default to 1.
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
    if ((.not. present(occupancy)) .or. occupancy) then
      this%trunks(trunk) = ibset(this%trunks(trunk), i - 1)
    else
      this%trunks(trunk) = ibclr(this%trunks(trunk), i - 1)
    endif
    call this%destroy_orbitals_cache()
  end subroutine set_orbital

  function get_n_elec(this) result(n_elec)
    class(spin_det_type), intent(in) :: this
    integer :: n_elec
    integer :: i

    n_elec = 0
    if (allocated(this%orbitals_cache)) then
      n_elec = size(this%orbitals_cache)
      return
    endif
    do i = 1, this%n_trunks
      n_elec = n_elec + popcnt(this%trunks(i))
    enddo
  end function get_n_elec

  subroutine get_elec_orbitals(this, orbitals, n_elec_opt)
    class(spin_det_type), intent(inout) :: this
    integer, allocatable, intent(out) :: orbitals(:)
    integer, optional, intent(in) :: n_elec_opt
    integer :: n_elec
    integer :: i
    integer :: pos
    type(spin_det_type), pointer :: tmp_spin_det

    if (allocated(this%orbitals_cache)) then
      allocate(orbitals(size(this%orbitals_cache)))
      orbitals(:) = this%orbitals_cache
      return
    endif

    if (present(n_elec_opt)) then
      n_elec = n_elec_opt
    else
      n_elec = this%get_n_elec()
    endif
    allocate(orbitals(n_elec))
    tmp_spin_det => tmp_spin_det_instances(1)
    call tmp_spin_det%resize(this%n_trunks)
    do i = 1, this%n_trunks
      tmp_spin_det%trunks(i) = this%trunks(i)
    enddo
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

  subroutine assign_spin_det(dest, src)
    type(spin_det_type), pointer, intent(in) :: src
    type(spin_det_type), pointer, intent(out) :: dest
   
    call dest%resize(src%n_trunks)
    dest%trunks(:) = src%trunks(:)
    call dest%destroy_orbitals_cache()
  end subroutine assign_spin_det

  function equal_spin_det(left, right) result(is_equal)
    type(spin_det_type), pointer, intent(in) :: left, right
    logical :: is_equal
    integer :: i

    do i = left%n_trunks, 1, -1
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

  function eor_spin_det(left, right) result(res)
    type(spin_det_type), pointer, intent(in) :: left, right
    type(spin_det_type), pointer :: res
    integer :: i
    integer :: n_trunks

    n_trunks = left%n_trunks
    if (n_trunks /= right%n_trunks) then
      call backtrace
      stop 'eor_spin_det input size mismatch.'
    endif

    res => new_spin_det(n_trunks)
    do i = 1, n_trunks
      res%trunks(i) = ieor(left%trunks(i), right%trunks(i))
    enddo
  end function eor_spin_det

end module spin_det_module
