module spin_det_module
  
  use constants_module
  use types_module

  implicit none

  private

  public :: spin_det_type
  public :: assignment(=)
  public :: operator(==)
  public :: operator(<)
  public :: operator(>)

  type spin_det_type
    private
    integer(LONG), allocatable, public :: trunks(:)
    integer, pointer :: elec_orbitals(:) => null()
    contains
      procedure, public :: from_eor
      procedure, public :: initialize_by_det_size
      procedure, public :: set_orbital
      procedure, public :: get_orbital
      procedure, public :: get_n_elec
      procedure, public :: get_elec_orbitals
      procedure, public :: destroy_orbitals_cache
      procedure, public :: resize
      procedure, public :: print
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

  contains
 
  subroutine initialize_by_det_size(this, det_size)
    class(spin_det_type), intent(inout) :: this
    integer, intent(in) :: det_size
    call this%resize(det_size)
    this%trunks = 0
    call this%destroy_orbitals_cache()
  end subroutine initialize_by_det_size

  subroutine destroy_orbitals_cache(this)
    class(spin_det_type), intent(inout) :: this
    if (associated(this%elec_orbitals)) then
      deallocate(this%elec_orbitals)
      nullify(this%elec_orbitals)
    endif
  end subroutine destroy_orbitals_cache

  function get_orbital(this, orbital_idx) result(is_occupied)
    class(spin_det_type), intent(in) :: this
    integer, intent(in) :: orbital_idx
    logical :: is_occupied
    integer :: i
    integer :: trunk

    if (orbital_idx > size(this%trunks) * C%TRUNK_SIZE) then
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
    logical, optional, intent(in) :: occupancy
    integer :: i
    integer :: trunk

    if (orbital_idx > size(this%trunks) * C%TRUNK_SIZE) then
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
    do i = 1, size(this%trunks)
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
    type(spin_det_type) :: tmp

    if (associated(this%elec_orbitals)) then
      allocate(orbitals(size(this%elec_orbitals)))
      orbitals(:) = this%elec_orbitals
      return
    endif
    if (present(n_elec_opt)) then
      n_elec = n_elec_opt
    else
      n_elec = this%get_n_elec()
    endif
    if (.not. allocated(orbitals)) then
      allocate(orbitals(n_elec))
    endif
    tmp = this
    do i = 1, n_elec
      pos = tmp%trail_zero() + 1
      orbitals(i) = pos
      call tmp%set_orbital(pos, .false.) 
    enddo
    allocate(this%elec_orbitals(n_elec))
    this%elec_orbitals(:) = orbitals(:)
  end subroutine get_elec_orbitals

  subroutine print(this)
    class(spin_det_type), intent(in) :: this
    integer :: i

    do i = 1, size(this%trunks)
      write (6, '(I0, A, B0.63)') i, ': ', this%trunks(i) 
    enddo
  end subroutine print

  function trail_zero(this) result(cnt)
    class(spin_det_type), intent(in) :: this
    integer :: cnt
    integer :: i 
    cnt = 0
    do i = 1, size(this%trunks)
      if (this%trunks(i) == 0) then
        cnt = cnt + C%TRUNK_SIZE
      else
        cnt = cnt + trailz(this%trunks(i))
        return
      endif
    enddo
  end function trail_zero

  subroutine resize(this, det_size)
    class(spin_det_type), intent(inout) :: this
    integer, intent(in) :: det_size

    if (allocated(this%trunks)) then
      if (size(this%trunks) /= det_size) then
        deallocate(this%trunks)
      end if
    endif
    if (.not. allocated(this%trunks)) then
      allocate(this%trunks(det_size))
    endif
    call this%destroy_orbitals_cache()
  end subroutine resize

  subroutine assign_spin_det(dest, src)
    type(spin_det_type), intent(in) :: src
    type(spin_det_type), intent(out) :: dest
    integer :: src_size
    src_size = size(src%trunks)
    call dest%resize(src_size)
    dest%trunks(:) = src%trunks(:)
    call dest%destroy_orbitals_cache()
  end subroutine assign_spin_det

  function compare_spin_det(left, right) result(res)
    class(spin_det_type), intent(in) :: left, right
    integer :: res
    integer :: i
    if (size(left%trunks) /= size(right%trunks)) then
      call backtrace
      stop 'compare_spin_det size mismatch.'
      return
    endif
    do i = size(left%trunks), 1, -1
      if (left%trunks(i) < right%trunks(i)) then
        res = -1
        return
      else if (left%trunks(i) > right%trunks(i)) then
        res = 1
        return
      endif
    enddo
    res = 0
  end function

  function equal_spin_det(left, right) result(is_equal)
    class(spin_det_type), intent(in) :: left, right
    logical :: is_equal

    is_equal = (compare_spin_det(left, right) == 0)
  end function equal_spin_det

  function lt_spin_det(left, right) result(is_lt)
    class(spin_det_type), intent(in) :: left, right
    logical :: is_lt

    is_lt = (compare_spin_det(left, right) < 0)
  end function lt_spin_det

  function gt_spin_det(left, right) result(is_gt)
    class(spin_det_type), intent(in) :: left, right
    logical :: is_gt

    is_gt = (compare_spin_det(left, right) > 0)
  end function gt_spin_det

  subroutine from_eor(this, left, right)
    class(spin_det_type), intent(inout) :: this
    type(spin_det_type), intent(in) :: left, right
    integer :: i
    integer :: left_size
   
    left_size = size(left%trunks)
    if (left_size /= size(right%trunks)) then
      call backtrace
      stop 'eor_spin_det input size mismatch.'
    endif
    call this%resize(left_size)
    do i = 1, left_size
      this%trunks(i) = ieor(left%trunks(i), right%trunks(i))
    enddo
    call this%destroy_orbitals_cache()
  end subroutine from_eor

end module spin_det_module
