module solver_controller_class

  use config_class
  use heg_solver_class
  use test_system_solver_class

  implicit none

  private

  public :: solver_controller

  type solver_controller
    private
    type(heg_solver) :: heg_solver_instance
    type(test_system_solver), public :: test_system_solver_instance
    ! Not using allocatable for memory locality at the cost of O(1) storage.
    contains
      procedure, public :: solve
  end type solver_controller

  contains

    subroutine solve(this, config_instance)
      class(solver_controller), intent(inout) :: this
      type(config), intent(in) :: config_instance

      select case (config_instance%system)
        case ('heg')
          call this%heg_solver_instance%solve(config_instance)
        case ('test_system')
          call this%test_system_solver_instance%solve(config_instance)
        case default
          stop 'No solver available for this system.'
      end select
    end subroutine solve

end module solver_controller_class
