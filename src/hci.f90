program hci

  use solver_controller_module

  implicit none

  type(solver_controller_type) :: solver_controller

  call solver_controller%start()

end program hci
