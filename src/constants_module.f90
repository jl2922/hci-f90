module constants_module

  use types_module

  implicit none

  private

  ! Namespace.
  public :: C
  protected :: C ! Keep invariant.

  type constants_type
    real(DOUBLE) :: EPS = 1.0e-15_DOUBLE
    real(DOUBLE) :: PI = 3.141592653589793238468_DOUBLE
    integer :: TRUNK_SIZE = 63 ! Standard fortran doesn't support unsigned.
  end type constants_type

  type(constants_type) :: C

end module
