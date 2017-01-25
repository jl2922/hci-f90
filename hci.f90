program hci

  use config_class, only : config
  use solver_class, only : solver

  implicit none
  
  type(config) :: hci_config
  type(solver) :: hci_solver

  call hci_config%read_configs()
  call hci_solver%solve(hci_config)
  call hci_solver%print_results()

end program hci
