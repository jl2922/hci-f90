module config_class
  
  use types_class

  implicit none

  private

  public :: config

  integer, parameter :: SYSTEM_NAME_LENGTH = 80
  integer, parameter :: CONFIG_FILE_UNIT = 100

  type heg_config
    private
    integer, public :: n_elec
    integer, public :: n_up
    integer, public :: n_dn
    real(DOUBLE), public :: r_cutoff
    real(DOUBLE), public :: r_s
  end type heg_config
  
  type config
    private
    character(len=SYSTEM_NAME_LENGTH), public :: system
    real(DOUBLE), public :: eps_var
    real(DOUBLE), public :: eps_pt
    real(DOUBLE), public :: tol
    logical, public :: dump_var_wf
    type(optional_double), public :: eps_pt_det
    type(optional_int), public :: n_samples
    type(heg_config), public :: heg
    contains
      procedure, public :: read_configs
      procedure :: read_heg_configs
  end type config

  contains

    subroutine read_configs(this)
      class(config), intent(inout) :: this

      character(len=SYSTEM_NAME_LENGTH) :: system

      ! Default values.
      real(DOUBLE) :: eps_var = 1.0e-3_DOUBLE
      real(DOUBLE) :: eps_pt = 1.0e-5_DOUBLE
      real(DOUBLE) :: tol = 1.0e-3_DOUBLE ! 1 mHa.
      logical :: dump_var_wf = .true.

      ! Initialized to -1 for checking availablity in namelist.
      real(DOUBLE) :: eps_pt_det = -1 ! eps for the deterministic part of PT.
      real(DOUBLE) :: n_samples = -1 ! For stochastic PT. Real to allow exponential format.

      integer :: io_err

      namelist /hci/ system, eps_var, eps_pt, eps_pt_det, tol, dump_var_wf, n_samples

      call create_config_file(CONFIG_FILE_UNIT)

      rewind(config_file_unit)
      read(unit=CONFIG_FILE_UNIT, nml=hci, iostat=io_err)
      if (io_err > 0) stop 'Cannot read HCI configs.'

      this%system = system
      this%eps_var = eps_var
      this%eps_pt = eps_pt
      this%tol = tol
      this%dump_var_wf = dump_var_wf
      
      if (eps_pt_det > 0) then
        this%eps_pt_det%instance = eps_pt_det
        this%eps_pt_det%is_present = .true.
      else
        this%eps_pt_det%is_present = .false.
      end if

      if (n_samples > 0) then
        this%n_samples%instance = nint(n_samples)
        this%n_samples%is_present = .true.
      else
        this%n_samples%is_present = .false.
      end if

      write (6, '(A)') '[CONFIGURATIONS]'
      write (6, '(A, A)') 'system: ', system
      write (6, '(A, G0.10)') 'eps_var: ', eps_var 
      write (6, '(A, G0.10)') 'eps_pt: ', eps_pt
      write (6, '(A, G0.10)') 'tol: ', tol
      write (6, '(A, L1)') 'dump_var_wf: ', dump_var_wf
      if (this%eps_pt_det%is_present) then
        write (6, '(A, G0.10)') 'eps_pt_det: ', this%eps_pt_det%instance
      else
        write (6, '(A)') 'eps_pt_det: based on environment.'
      end if
      if (this%n_samples%is_present) then
        write (6, '(A, I0)') 'n_samples: ', this%n_samples%instance
      else
        write (6, '(A)') 'n_samples: based on environment.'
      end if
      write (6, '()')

      select case (system)
        case ('heg')
          call this%read_heg_configs()
        case default
          stop 'Unrecognized system.'
      end select

      close(config_file_unit)

    end subroutine read_configs

    subroutine read_heg_configs(this)
      class(config), intent(inout) :: this

      real(DOUBLE) :: r_cutoff = -1
      real(DOUBLE) :: r_s = -1
      integer :: n_up = -1
      integer :: n_dn = -1
      integer :: n_elec
      namelist /heg/ r_cutoff, r_s, n_up, n_dn

      integer :: io_err

      rewind(CONFIG_FILE_UNIT)
      read(unit=CONFIG_FILE_UNIT, nml=heg, iostat=io_err)
      if (io_err > 0) stop 'Cannot read HEG configs.'

      if (r_cutoff < 0 .or. r_s < 0 .or. n_up < 0 .or. n_dn < 0) then
        stop 'r_cutoff, r_s, n_up, n_dn must be provided for heg system.'
      end if

      this%heg%r_cutoff = r_cutoff
      this%heg%r_s = r_s
      this%heg%n_up = n_up
      this%heg%n_dn = n_dn
      n_elec = n_up + n_dn
      this%heg%n_elec = n_elec

      write (6, '(A, G0.10)') 'r_cutoff: ', r_cutoff
      write (6, '(A, G0.10)') 'r_s: ', r_s
      write (6, '(A, I0)') 'n_up: ', n_up
      write (6, '(A, I0)') 'n_dn: ', n_dn
      write (6, '(A, I0)') 'n_elec: ', n_elec
      write (6, '()')

    end subroutine read_heg_configs

    subroutine create_config_file(config_file_unit)
      ! Make a snapshot of the config input.
      integer, intent(in) :: config_file_unit
      integer, parameter :: MAX_LINE_LENGTH = 256
      character(MAX_LINE_LENGTH) :: line
      integer :: io_err

      open(unit=config_file_unit, status='scratch', recl=MAX_LINE_LENGTH, delim='APOSTROPHE')

      do
        read(5, '(A)', iostat = io_err) line
        if (io_err < 0) then
          ! End of file reached.
          return
        else if (io_err == 0) then
          write(config_file_unit, '(A)') line
        else 
          stop 'Cannot read config file.'
        end if
      end do
    end subroutine create_config_file

end module config_class
