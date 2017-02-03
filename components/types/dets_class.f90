module dets_class

  use types_class
  
  implicit none

  private

  public :: dets

  type dets
    private
    integer, public :: n = 0
    integer, public :: det_size
    integer(INT64), allocatable, public :: up(:, :)
    integer(INT64), allocatable, public :: dn(:, :)
    logical :: is_sorted = .false.
    contains
      procedure, public :: sort
      procedure, public :: filter_by_sorted_dets
      procedure, public :: merge_sorted_dets
  end type dets

  ! No individual det type for keeping consistence and data locality.

  contains

    subroutine sort(this)
      class(dets), intent(inout) :: this
      
      this%is_sorted = .true.
    end subroutine sort

    subroutine filter_by_sorted_dets(this, ref_dets)
      class(dets), intent(inout) :: this
      type(dets), intent(in) :: ref_dets
    end subroutine filter_by_sorted_dets

    subroutine merge_sorted_dets(this, new_dets)
      class(dets), intent(inout) :: this
      type(dets), intent(in) :: new_dets
    end subroutine merge_sorted_dets

end module dets_class
