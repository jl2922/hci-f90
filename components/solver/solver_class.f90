module solver_class
  
  use config_class
  use types_class
  use dets_class
  
  implicit none

  private

  public :: solver

  integer, parameter :: MAX_VAR_ITERATION = 20

  type solver
    private
    integer, public :: n_orb
    integer, public :: det_size ! Number of integers used to represent each det.
    real(DOUBLE), public :: HF_energy
    real(DOUBLE), public :: var_energy
    real(DOUBLE), public :: pt_det_energy
    real(DOUBLE), public :: pt_st_energy
    real(DOUBLE), public :: pt_st_uncert
    real(DOUBLE), public :: max_abs_H
    real(DOUBLE), allocatable, public :: coef(:)
    type(dets), public :: dets
    contains
      procedure, public :: solve
      procedure :: diagonalize
      procedure :: find_connected_dets
      procedure :: get_next_dets
      procedure :: setup
      procedure :: summary
      procedure :: pt ! 2nd-order perturabation theory correction.
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

        call this%get_next_dets(config_instance%eps_var)
        call this%diagonalize(energy_cur)

        write (6, '(A, F0.10)') 'Energy: ', energy_cur
        write (6, '(A, G0.10)') 'Number of dets: ', this%dets%n
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

    subroutine diagonalize(this, energy)
      class(solver), intent(inout) :: this
      real(DOUBLE), intent(out) :: energy

      energy = -1.0
    end subroutine diagonalize

    subroutine find_connected_dets(this, det_up, det_dn, connected_dets, eps_min)
      class(solver), intent(inout) :: this
      integer(INT64), intent(in) :: det_up(:), det_dn(:)
      type(dets), intent(out) :: connected_dets
      real(DOUBLE), intent(in) :: eps_min
    end subroutine find_connected_dets

    subroutine get_next_dets(this, eps_var)
      class(solver), intent(inout) :: this
      real(DOUBLE), intent(in) :: eps_var

      integer :: i
      type(dets) :: connected_dets
      type(dets) :: new_dets
      
      do i = 1, this%dets%n
        call this%find_connected_dets(this%dets%up(:, i), this%dets%dn(:, i), &
            & connected_dets, eps_var / abs(this%coef(i))) 
        call connected_dets%filter_by_sorted_dets(this%dets)
        call new_dets%merge_sorted_dets(connected_dets)
      end do
      
      call this%dets%merge_sorted_dets(new_dets)
    end subroutine get_next_dets

    subroutine setup(this, config_instance)
      class(solver), intent(inout) :: this
      type(config), intent(in) :: config_instance

      stop 'Default setup has not been overloaded.'
    end subroutine setup

    subroutine summary(this)
      class(solver), intent(inout) :: this

      real(DOUBLE) :: total_energy

      write (6, '(A, F0.10)') 'Variational Energy: ', this%var_energy
      write (6, '(A, G0.10)') 'Number of variational dets: ', this%dets%n
      write (6, '(A, F0.10)') 'Deterministic PT Energy: ', this%pt_det_energy
      write (6, '(A, F0.10, A, F0.10)') &
          & 'Stochastic PT Energy: ', this%pt_st_energy, ' +- ', this%pt_st_uncert
      total_energy = this%var_energy + this%pt_det_energy + this%pt_st_energy
      write (6, '(A, F0.10, A, F0.10)') &
          & 'Total Energy: ', total_energy, ' +- ', this%pt_st_uncert
    end subroutine summary

    subroutine pt(this, config_instance)
      class(solver), intent(inout) :: this
      type(config), intent(in) :: config_instance
    end subroutine pt

end module solver_class
