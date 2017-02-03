module test_system_solver_class

  use config_class, only : config
  use solver_class, only : solver

  implicit none

  private

  public :: test_system_solver

  type test_system_data
    logical :: setup_called = .false.
  end type test_system_data

  type, extends(solver) :: test_system_solver
    private
    type(test_system_data), public :: test_system
    contains
      procedure, public :: setup
  end type test_system_solver

  contains

    subroutine setup(this, config_instance)
      class(test_system_solver), intent(inout) :: this
      type(config), intent(in) :: config_instance

      this%test_system%setup_called = .true.
    end subroutine setup
end module test_system_solver_class
