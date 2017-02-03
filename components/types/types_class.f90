module types_class
  
  implicit none

  private

  public :: INT8, INT16, INT32, INT64
  public :: SINGLE, DOUBLE
  public :: LONG

  public :: INT_3_DOUBLE

  public :: optional_int
  public :: optional_double

  public :: abclr, abcnt, abeor, abeq, abprint, abset, abtest, abtrailz

  integer, parameter :: INT8 = selected_int_kind(2)
  integer, parameter :: INT16 = selected_int_kind(4)
  integer, parameter :: INT32 = selected_int_kind(9)
  integer, parameter :: INT64 = selected_int_kind(18)

  integer, parameter :: SINGLE = selected_real_kind(6, 37)
  integer, parameter :: DOUBLE = selected_real_kind(13, 300)

  integer, parameter :: LONG = INT64

  type INT_3_DOUBLE
    integer :: int_3(3)
    real(DOUBLE) :: double
  end type INT_3_DOUBLE

  type optional_instance
    private
    logical, public :: is_present = .false.
  end type optional_instance

  type, extends(optional_instance) :: optional_double
    private
    real(DOUBLE), public :: instance = 0
  end type optional_double

  type, extends(optional_instance) :: optional_int
    private
    integer, public :: instance = 0
  end type optional_int

  contains

    subroutine abclr(arr, pos)
      ! Clear the bit at pos in array representation of bits.
      integer(INT64), intent(inout) :: arr(:)
      integer, intent(in) :: pos

      integer :: i
      integer :: trunk

      if (pos > size(arr) * 64) then
        call backtrace
        stop 'abclr with index out of range.'
      end if

      i = pos
      trunk = 1

      do while (i > 64)
        i = i - 64
        trunk = trunk + 1
      end do
      
      arr(trunk) = ibclr(arr(trunk), i - 1)
    end subroutine abclr

    type(integer) function abcnt(arr) result(cnt)
      ! Count the number of bits set to 1 in array representation of bits.
      integer(INT64), intent(in) :: arr(:)

      integer :: i

      cnt = 0
      do i = 1, size(arr)
        cnt = cnt + popcnt(arr(i))
      end do
    end function abcnt
    
    subroutine abeor(arr1, arr2, arr_res)
      ! Perform eor on array representation of bits.
      integer(INT64), intent(in) :: arr1(:), arr2(:)
      integer(INT64), allocatable, intent(out) :: arr_res(:)

      integer :: size1, size2
      integer :: i

      size1 = size(arr1)
      size2 = size(arr2)

      if (size1 /= size2) then
        call backtrace
        stop 'Size mismatch for abeor.'
      end if

      allocate(arr_res(size1))

      do i = 1, size1
        arr_res(i) = ieor(arr1(i), arr2(i))
      end do
    end subroutine abeor

    type(logical) function abeq(arr1, arr2) result(is_equal)
      ! Check equality for array representation of bits.
      integer(INT64), intent(in) :: arr1(:), arr2(:)

      integer :: size1, size2
      integer :: i

      size1 = size(arr1)
      size2 = size(arr2)

      if (size1 /= size2) then
        is_equal = .false.  
        return
      end if

      do i = 1, size1
        if (.not. (arr1(i) == arr2(i))) then
          is_equal = .false.
          return
        end if
      end do

      is_equal = .true.
    end function abeq

    subroutine abprint(arr)
      ! Print in binary format the array representation of bits.
      integer(INT64), intent(inout) :: arr(:)

      integer :: i

      do i = 1, size(arr)
        write (6, '(I0, A, B0.64)') i, ':', arr(i)
      end do
    end subroutine abprint

    subroutine abset(arr, pos)
      ! Set the bit at pos for array representation of bits.
      integer(INT64), intent(inout) :: arr(:)
      integer, intent(in) :: pos

      integer :: i
      integer :: trunk

      if (pos > size(arr) * 64) then
        call backtrace
        stop 'abset with index out of range.'
      end if

      i = pos
      trunk = 1

      do while (i > 64)
        i = i - 64
        trunk = trunk + 1
      end do
      
      arr(trunk) = ibset(arr(trunk), i - 1)
    end subroutine abset

    type(logical) function abtest(arr, pos) result(is_set)
      ! Test whether the bit at pos is set in array representation of bits.
      integer(INT64), intent(in) :: arr(:)
      integer, intent(in) :: pos

      integer :: i
      integer :: trunk

      if (pos > size(arr) * 64) then
        call backtrace
        stop 'abtest with index out of range.'
      end if

      i = pos
      trunk = 1

      do while (i > 64)
        i = i - 64
        trunk = trunk + 1
      end do
      
      is_set = btest(arr(trunk), i - 1)
    end function abtest

    type(integer) function abtrailz(arr) result(cnt)
      ! Count the number of trailing zeros in array representation of bits.
      integer(INT64), intent(in) :: arr(:)

      integer :: i

      cnt = 0
      do i = 1, size(arr)
        if (arr(i) == 0) then
          cnt = cnt + 64
        else
          cnt = cnt + trailz(arr(i))
          exit
        end if
      end do
    end function abtrailz

end module types_class
