!
! Creating a new simluation type involves:
!
!   1. a 'use' statement to include the new code
!   2. a new entry in the 'case' statement within the 'create_core_simulation' function.
!
module m_simulation_factory

  use m_system
  use m_simulation

  implicit none

  private

  interface create
    module procedure create_simulation
  end interface

  public ::  create, get_algorithm_from_string

  contains

!---------------------------------------------------------------------------------------------------
! Externally called routine that ultimatly passes back a t_simulation pointer.
!---------------------------------------------------------------------------------------------------
  function create_simulation( opts )
!---------------------------------------------------------------------------------------------------
    use m_input_file
    use m_parameters
    implicit none

    type( t_options ), intent(in) :: opts
    class( t_simulation ), pointer :: create_simulation

    class( t_simulation ), pointer :: sim

    sim => create_core_simulation( opts )

    create_simulation => sim
  end function
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Utility function (called from os-main.f03) that maps simulation type strings
!   to the internally used integer parameters
!---------------------------------------------------------------------------------------------------
function get_algorithm_from_string(simulation_label)
  use m_parameters
  implicit none
  character(len=20) :: simulation_label
  integer           :: get_algorithm_from_string

  select case( trim( simulation_label ) )
    case( "standard" )
      get_algorithm_from_string = p_standard
    case default
      call error_and_die( "Error reading global simulation parameters"//NEW_LINE('A')//&
                         &"algorithm: \'"//trim( simulation_label )//"\' is invalid or the code is not compiled to support it. Aborting.")
  end select
end function

!---------------------------------------------------------------------------------------------------
! Core that creates the basic simulation types.
!   This is the code that used to be in 'os-main%initialize'
!---------------------------------------------------------------------------------------------------
function create_core_simulation( opts )
!---------------------------------------------------------------------------------------------------
  use m_input_file
  use m_parameters
  implicit none

  type( t_options ), intent(in) :: opts
  class( t_simulation ), pointer :: create_core_simulation

  class( t_simulation ), pointer :: sim

    ! Allocate simulation object based on algorithm
  select case ( opts%algorithm )
    case( p_standard )
      ! Standard EM-PIC simulation
      allocate( t_simulation :: sim )
    case default
      call error_and_die( '(*error*) Invalid simulation mode, only standard is currently implemented.')
  end select
  create_core_simulation => sim
end function
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Utility function to print out fatal error conditions.
!---------------------------------------------------------------------------------------------------
subroutine error_and_die( message )
!---------------------------------------------------------------------------------------------------
    use m_parameters
    implicit none
    character(len=*) :: message

    if ( mpi_node() == 0 ) then
      write(0,*) message
      call abort_program( p_err_invalid )
    endif

  end subroutine
!---------------------------------------------------------------------------------------------------
end module


!---------------------------------------------------------------------------------------------------
! Version of the simulation factory that rests outside a module
!   to avoid circular depencies in cases where a simulation wants to,
!   inside it's own code, use the factory to create other simulation objects
!   ( an example of this is the 'tandem' simulation class )
!---------------------------------------------------------------------------------------------------
! function global_simulation_factory( opts )

!   use m_system, only: t_options
!   use m_simulation, only: t_simulation
!   use m_simulation_factory

!   implicit none

!   !integer, intent(in) :: algorithm
!   type( t_options ), intent(inout) :: opts

!   class( t_simulation ), pointer :: global_simulation_factory


!   global_simulation_factory => create( opts  )

! end function
