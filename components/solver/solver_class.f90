module solver_class
  
  use config_class, only : config
  
  implicit none

  private

  public :: solver

  type solver
    private
    contains
      procedure, public :: solve
      procedure, public :: print_results
  end type solver

  contains

    subroutine solve(this, hci_config)
      class(solver), intent(inout) :: this
      type(config), intent(in) :: hci_config
      write (6, *) 'Solving...'
    end subroutine solve

    subroutine print_results(this)
      class(solver), intent(inout) :: this
      write (6, *) 'Printing results...'
    end subroutine print_results

end module solver_class
