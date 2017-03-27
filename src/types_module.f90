module types_module

  use, intrinsic :: iso_fortran_env

  implicit none

  private

  public :: DOUBLE
  public :: LONG

  public :: INT_PAIR
  public :: INT_3_DOUBLE

  public :: optional_int
  public :: optional_double

  integer, parameter :: LONG = INT64
  integer, parameter :: DOUBLE = REAL64

  type INT_PAIR
    integer :: i1 = 0
    integer :: i2 = 0
  end type INT_PAIR

  type INT_3_DOUBLE
    integer :: int_3(3) = 0
    real(DOUBLE) :: double = 0.0_DOUBLE
  end type INT_3_DOUBLE

  type optional_instance
    private
    logical, public :: is_present = .false.
  end type optional_instance

  type, extends(optional_instance) :: optional_double
    private
    real(DOUBLE), public :: instance = 0.0_DOUBLE
  end type optional_double

  type, extends(optional_instance) :: optional_int
    private
    integer, public :: instance = 0
  end type optional_int

end module types_module
