module types_module

  implicit none

  private

  public :: INT8, INT16, INT32, INT64
  public :: SINGLE, DOUBLE
  public :: LONG

  public :: INT_PAIR
  public :: INT_3_DOUBLE

  public :: optional_int
  public :: optional_double

  integer, parameter :: INT8 = selected_int_kind(2)
  integer, parameter :: INT16 = selected_int_kind(4)
  integer, parameter :: INT32 = selected_int_kind(9)
  integer, parameter :: INT64 = selected_int_kind(18)

  integer, parameter :: SINGLE = selected_real_kind(6, 37)
  integer, parameter :: DOUBLE = selected_real_kind(13, 300)

  integer, parameter :: LONG = INT64

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
