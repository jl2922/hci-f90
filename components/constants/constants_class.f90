module constants_class

  use types_class

  implicit none

  private

  public :: EPS, PI

  real(DOUBLE), parameter :: EPS = 1.0e-15_DOUBLE
  real(DOUBLE), parameter :: PI = 3.141592653589793238468_DOUBLE
end module constants_class
