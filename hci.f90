program hci

  use config_class
  use solver_controller_class

  implicit none
  
  type(config) :: config_instance
  type(solver_controller) :: solver_controller_instance

  call config_instance%read_configs()
  call solver_controller_instance%solve(config_instance)

end program hci
