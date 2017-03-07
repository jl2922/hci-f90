module combinatorics_module

  implicit none

  private

  ! Namespaces.
  public :: combinatorics_type
  public :: combinatorics

  type combinatorics_type
    contains
      procedure, public :: C
      procedure, public :: P
      procedure, public :: factorial
  end type combinatorics_type

  type(combinatorics_type) :: combinatorics

  contains

  function C(this, n, k) result(res)
    class(combinatorics_type), intent(in) :: this
    integer, intent(in) :: n, k
    integer :: res
    integer :: i

    res = this%P(n, k) / this%factorial(k)
  end function C

  function P(this, n, k) result(res)
    class(combinatorics_type), intent(in) :: this
    integer, intent(in) :: n, k
    integer :: res
    integer :: i

    res = 1
    do i = n, n - k + 1, -1
      res = res * i
    end do
  end function P

  function factorial(this, k) result(res)
    class(combinatorics_type), intent(in) :: this
    integer, intent(in) :: k
    integer :: res
    integer :: i

    res = 1
    do i = k, 1, -1
      res = res * i
    end do
  end function factorial

end module combinatorics_module
