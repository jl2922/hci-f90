module types_class
  
  implicit none

  private

  public :: INT8, INT16, INT32, INT64
  public :: SINGLE, DOUBLE
  public :: LONG

  public :: optional_int
  public :: optional_double

  integer, parameter :: INT8 = selected_int_kind(2)
  integer, parameter :: INT16 = selected_int_kind(4)
  integer, parameter :: INT32 = selected_int_kind(9)
  integer, parameter :: INT64 = selected_int_kind(18)

  integer, parameter :: SINGLE = selected_real_kind(6, 37)
  integer, parameter :: DOUBLE = selected_real_kind(13, 300)

  integer, parameter :: LONG = INT64

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

end module types_class
