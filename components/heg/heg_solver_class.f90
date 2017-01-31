module heg_solver_class
  
  use config_class, only : config
  use solver_class, only : solver

  implicit none

  private

  public :: heg_solver

  type, extends(solver) :: heg_solver
    integer :: id
    contains
      procedure, public :: setup
  end type heg_solver

  contains

    subroutine setup(this, config_instance)
      class(heg_solver), intent(inout) :: this
      type(config), intent(in) :: config_instance
      print *, 'heg setup procedure...'
    end subroutine setup

end module heg_solver_class
