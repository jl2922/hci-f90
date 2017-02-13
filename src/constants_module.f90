module constants_module

  use types_module

  implicit none

  private

  public :: C
  protected :: C ! Protected makes it invariant.

  ! Create a namespace for all constants.
  type constants_type
    real(DOUBLE) :: EPS = 1.0e-15
    real(DOUBLE) :: PI = 3.141592653589793238468_DOUBLE
    integer :: UP_SPIN = 1
    integer :: DOWN_SPIN = 0
  end type constants_type
  type(constants_type) :: C

end module
