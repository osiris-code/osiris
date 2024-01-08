#include "os-config.h"
#include "os-preprocess.fpp"

module m_zpulse

use m_system

use m_zpulse_std
use m_zpulse_wall
use m_zpulse_point
use m_zpulse_speckle

use m_input_file

use m_space
use m_emf_define, only: t_emf_bound, t_emf
use m_grid_define
use m_node_conf


implicit none

private

!-----------------------------------------------------------------------------------------
type :: t_zpulse_list

  class( t_zpulse ), pointer :: list => null()

contains

  procedure :: read_input => read_input_zpulse_list
  procedure :: init       => init_zpulse_list
  procedure :: cleanup    => cleanup_zpulse_list

  procedure :: launch => launch_zpulse_list

end type t_zpulse_list
!-----------------------------------------------------------------------------------------


public :: t_zpulse_list

contains


!-----------------------------------------------------------------------------------------
subroutine read_input_zpulse_list( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

  implicit none


  class( t_zpulse_list ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  type( t_space ), intent(in) :: g_space
  class (t_emf_bound), intent(in) :: bnd_con
  logical, dimension(:), intent(in) :: periodic
  class( t_grid ), intent(in) :: grid
  type( t_options ), intent(in) :: sim_options

  class( t_zpulse ), pointer :: tail, item

  tail => null()

  do
    select case (trim(input_file % get_section_name()))
    case('nl_zpulse')
        allocate( t_zpulse          :: item )
    case('nl_zpulse_wall')
        allocate( t_zpulse_wall     :: item )
    case('nl_zpulse_point')
        allocate( t_zpulse_point    :: item )
    case('nl_zpulse_speckle')
        allocate( t_zpulse_speckle  :: item )
    case default
        ! Not a zpulse section, exit
        exit
    end select

    ! Read input for zpulse
    call item % read_input( input_file, g_space, bnd_con, periodic, grid, sim_options )
    call item % check_dimensionality() 
    item % next => null()

    if (associated(tail)) then
      tail % next => item
    else
      this % list => item
    endif

    tail => item

  enddo

end subroutine read_input_zpulse_list
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine init_zpulse_list( this, restart, t, dt, emf, g_space, nx_p_min, no_co, &
                             interpolation )

  implicit none

  class( t_zpulse_list ), intent(inout) :: this
  logical, intent(in) :: restart
  real( p_double ), intent(in) :: t, dt
  class( t_emf ),  intent(in) :: emf
  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: interpolation

  class( t_zpulse ), pointer :: zpulse

  zpulse => this%list

  do
    if ( .not. associated(zpulse)) exit
    call zpulse % init( restart, t, dt, emf%b, g_space, nx_p_min, no_co, interpolation )

    zpulse => zpulse % next
  enddo

end subroutine init_zpulse_list
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine cleanup_zpulse_list( this )

  implicit none

  class( t_zpulse_list ), intent(inout) :: this
  class( t_zpulse ), pointer :: zpulse, next

  zpulse => this % list

  do
    if ( .not. associated(zpulse)) exit
    next => zpulse % next
    call zpulse % cleanup()
    deallocate( zpulse )
    zpulse => next
  enddo

end subroutine cleanup_zpulse_list
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_list( this, emf, g_space, &
                                nx_p_min, g_nx, t, dt, no_co )

  class( t_zpulse_list ), intent(inout) :: this
  class( t_emf ) ,  intent(inout) :: emf

  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t
  real(p_double), intent(in) :: dt
  class( t_node_conf ), intent(in) :: no_co

  class( t_zpulse ), pointer :: zpulse


  zpulse => this % list

  do
    if (.not. associated(zpulse)) exit
    call zpulse % launch( emf, g_space, nx_p_min, g_nx, t, dt, no_co )
    zpulse => zpulse % next
  enddo

end subroutine launch_zpulse_list
!-----------------------------------------------------------------------------------------


end module m_zpulse
