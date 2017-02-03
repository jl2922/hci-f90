program solver_controller_class_test

  use solver_controller_class, only : solver_controller
  use config_class, only : config
  use assert_class

  implicit none

  call test_solve()

  contains

    subroutine test_solve()
      type(config) :: config_instance
      type(solver_controller) :: solver_controller_instance

      config_instance%system = 'test_system'

      call solver_controller_instance%solve(config_instance)

      call assert_true( &
          & solver_controller_instance%test_system_solver_instance%test_system%setup_called, & 
          & 'solver_controller_solve')
    end subroutine test_solve
end program solver_controller_class_test
