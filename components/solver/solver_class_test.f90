program solver_class_test

  use test_system_solver_class, only : test_system_solver
  use config_class, only : config
  use assert_class
  use report_class

  implicit none

  call test_setup()

  contains

    subroutine test_setup()
      type(config) :: config_instance
      type(test_system_solver) :: test_system_solver_instance

      config_instance%system = 'test_system'

      call test_system_solver_instance%setup(config_instance)

      call assert_true(test_system_solver_instance%test_system%setup_called, 'solver_setup')
    end subroutine test_setup
end program solver_class_test
