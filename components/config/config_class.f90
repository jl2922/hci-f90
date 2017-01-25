module config_class
  
  implicit none

  private
  public :: config

  type config
    private
    contains
      procedure, public :: read_configs
  end type config

  contains

    subroutine read_configs(this)
      class(config), intent(inout) :: this
      write (6, *) 'Reading configs...'
    end subroutine read_configs

end module config_class
