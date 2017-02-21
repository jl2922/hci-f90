module solver_controller_module

  use heg_solver_module
  use utilities_module

  implicit none

  private

  public :: solver_controller

  integer, parameter :: SYSTEM_NAME_LENGTH = 128

  type solver_controller_type
    contains
      procedure, public :: start
      procedure :: create_config_file
      procedure :: get_system_name
  end type solver_controller_type

  type(solver_controller_type) :: solver_controller

  contains

  subroutine start(this)
    class(solver_controller_type), intent(inout) :: this
    integer :: config_file_unit
    character(len=SYSTEM_NAME_LENGTH) :: system_name

    config_file_unit = this%create_config_file()
    system_name = this%get_system_name(config_file_unit)
    write (6, '(A, A)') 'System: ', system_name
    write (6, '()')

    select case (system_name)
    case ('heg')
      call heg_solver%solve(config_file_unit)
    case default
      stop 'Unrecognized system'
    end select
  end subroutine start

  function create_config_file(this) result(config_file_unit)
    class(solver_controller_type), intent(inout) :: this
    integer :: config_file_unit
    integer, parameter :: MAX_LINE_LENGTH = 1024
    character(MAX_LINE_LENGTH) :: line
    integer :: io_err

    config_file_unit = util%get_free_unit()
    open(unit=config_file_unit, status='scratch', &
        & recl=MAX_LINE_LENGTH, delim='APOSTROPHE')
    do
      read(5, '(A)', iostat = io_err) line
      if (io_err < 0) then
        ! End of file reached.
        exit
      else if (io_err == 0) then
        write(config_file_unit, '(A)') line
      else
        stop 'Cannot read config file.'
      end if
    end do
  end function create_config_file

  function get_system_name(this, config_file_unit) result(system_name)
    class(solver_controller_type), intent(inout) :: this
    integer, intent(in) :: config_file_unit
    integer :: io_err
    character(len=SYSTEM_NAME_LENGTH) :: system_name
    character(len=SYSTEM_NAME_LENGTH) :: name
    namelist /system/ name

    rewind(config_file_unit)
    read(unit=config_file_unit, nml=system, iostat=io_err)
    if (io_err > 0) stop 'Cannot obtain system name.'
    system_name = name
  end function get_system_name

end module solver_controller_module
