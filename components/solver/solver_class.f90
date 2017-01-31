module solver_class
  
  use config_class, only : config
  use types_class
  
  implicit none

  private

  public :: solver

  integer, parameter :: MAX_VAR_ITERATION = 20

  type solver
    private
    integer :: n_dets
    integer :: det_size ! Number of integers used to represent each det.
    real(DOUBLE) :: var_energy
    real(DOUBLE) :: pt_det_energy
    real(DOUBLE) :: pt_st_energy
    real(DOUBLE) :: pt_st_uncert
    integer, allocatable :: dets_up(:, :), dets_dn(:, :)
    contains
      procedure, public :: solve
      procedure :: setup
      procedure :: get_next_dets
      procedure :: diagonalize
      procedure :: pt ! 2nd-order perturabation theory correction.
      procedure :: summary
  end type solver

  contains

    subroutine solve(this, config_instance)
      class(solver), intent(inout) :: this
      type(config), intent(in) :: config_instance
      
      integer :: iteration
      real(DOUBLE) :: energy_prev, energy_cur

      write (6, '(A)') '[SETUP]'
      call this%setup(config_instance)
      write (6, '()')
      
      write (6, '(A)') '[VARIATION]'
      write (6, '(A, G0.10)') 'eps_var: ', config_instance%eps_var
      energy_prev = 0

      do iteration = 1, MAX_VAR_ITERATION
        write (6, '(A, I0)') 'Iteration #', iteration

        call this%get_next_dets(config_instance)
        call this%diagonalize(energy_cur)

        write (6, '(A, F0.10)') 'Energy: ', energy_cur
        write (6, '(A, G0.10)') 'Number of dets: ', this%n_dets
        write (6, '()')
        call flush(6)

        if (abs(energy_prev - energy_cur) < 0.01 * abs(energy_cur)) then
          this%var_energy = energy_cur
          exit
        end if

        energy_prev = energy_cur
      end do

      write (6, '(A)') '[PERTURBATION]' 
      write (6, '(A, G0.10)') 'eps_pt: ', config_instance%eps_pt
      call this%pt(config_instance)
      write (6, '()')

      write (6, '(A)') '[SUMMARY]'
      call this%summary()
    end subroutine solve

    subroutine setup(this, config_instance)
      class(solver), intent(inout) :: this
      type(config), intent(in) :: config_instance

      stop 'Default setup has not been overloaded.'
    end subroutine setup

    subroutine get_next_dets(this, config_instance)
      class(solver), intent(inout) :: this
      type(config), intent(in) :: config_instance
      
      print *, 'getting next dets...'
    end subroutine get_next_dets

    subroutine diagonalize(this, energy)
      class(solver), intent(inout) :: this
      real(DOUBLE), intent(out) :: energy

      energy = -1.0
    end subroutine diagonalize

    subroutine pt(this, config_instance)
      class(solver), intent(inout) :: this
      type(config), intent(in) :: config_instance
    end subroutine pt

    subroutine summary(this)
      class(solver), intent(inout) :: this
      real(DOUBLE) :: total_energy

      write (6, '(A, F0.10)') 'Variational Energy: ', this%var_energy
      write (6, '(A, G0.10)') 'Number of variational dets: ', this%n_dets
      write (6, '(A, F0.10)') 'Deterministic PT Energy: ', this%pt_det_energy
      write (6, '(A, F0.10, A, F0.10)') &
          & 'Stochastic PT Energy: ', this%pt_st_energy, ' +- ', this%pt_st_uncert
      total_energy = this%var_energy + this%pt_det_energy + this%pt_st_energy
      write (6, '(A, F0.10, A, F0.10)') &
          & 'Total Energy: ', total_energy, ' +- ', this%pt_st_uncert
    end subroutine summary

end module solver_class
