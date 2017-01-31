module solver_controller_class

  use config_class, only : config
  use heg_solver_class, only : heg_solver

  implicit none

  private

  public :: solver_controller

  type solver_controller
    private
    contains
      procedure, public :: solve
  end type solver_controller

  contains

    subroutine solve(this, config_instance)
      class(solver_controller), intent(inout) :: this
      type(config), intent(in) :: config_instance
      type(heg_solver) :: heg_solver_instance

      select case (config_instance%system)
        case ('heg')
          call heg_solver_instance%solve(config_instance)
        case default
          stop 'No solver available for this system.'
      end select
    end subroutine solve

end module solver_controller_class
