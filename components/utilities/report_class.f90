module report_class
  
  implicit none

  private

  public :: report_success, report_fail

  contains

    subroutine report_success(msg)
      character(len=*), intent(in) :: msg
      write (6, '(A, A)') achar(27)//'[32m[SUCCESS]'//achar(27)//'[0m ', msg
    end subroutine report_success


    subroutine report_fail(msg)
      character(len=*), intent(in) :: msg
      write (6, '(A, A)') achar(27)//'[31m[FAIL]'//achar(27)//'[0m ', msg
      call backtrace
    end subroutine report_fail
end module report_class
