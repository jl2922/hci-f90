program hci

  use solver_controller_module

  implicit none

  type(solver_controller_type), pointer :: solver_controller

  solver_controller => new_solver_controller()
  call solver_controller%start()

end program hci
