module assert_class

  use report_class

  implicit none

  private

  public :: assert_true, assert_false

  contains

    subroutine assert_true(x, test)
      logical, intent(in) :: x
      character(len=*), intent(in) :: test
      if (x .neqv. .true.) then
        call report_fail(test)
      end if
      call report_success(test)
    end subroutine assert_true

    subroutine assert_false(x, test)
      logical, intent(in) :: x
      character(len=*), intent(in) :: test
      if (x .neqv. .false.) then
        call report_fail(test)
      end if
      call report_success(test)
    end subroutine assert_false
end module assert_class
