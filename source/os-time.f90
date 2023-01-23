!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     time class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"


module m_time

use m_restart
use m_utilities
use m_parameters

implicit none

private

!       string to id restart data
character(len=*), parameter :: p_time_rst_id = "time rst data - 0x0002"


type :: t_time

  private

  real(p_double) :: tmin, tmax, t

end type t_time


interface read_nml
  module procedure read_nml_time
end interface

interface setup
  module procedure setup_time
end interface

interface restart_write
  module procedure restart_write_time
end interface

interface restart_read
  module procedure restart_read_time
end interface

interface update
  module procedure update_time
end interface

interface test
  module procedure test_time
end interface

interface tmin
  module procedure tmin
end interface

interface tmax
  module procedure tmax
end interface

interface t
  module procedure t
end interface

!       declare things that should be public
public :: t_time, read_nml, setup
public :: restart_write
public :: update, test
public :: tmin, tmax, t


contains 


!---------------------------------------------------
!---------------------------------------------------
subroutine read_nml_time( this, input_file )

  use m_input_file

  implicit none

  type( t_time ), intent(out) :: this
  class( t_input_file ), intent(inout) :: input_file
  
  real(p_double)  ::  tmin, tmax

  namelist /nl_time/ tmin, tmax
  
  integer :: ierr

!       executable statements

  tmin = 0.0_p_double
  tmax = 0.0_p_double

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_time", ierr )
  
  if (ierr /= 0) then
	if (ierr < 0) then
	  print *, "Error reading time parameters"
	else 
	  print *, "Error: time parameters missing"
	endif
	print *, "aborting..."
	stop
  endif
  

  read (input_file%nml_text, nml = nl_time, iostat = ierr)
  if (ierr /= 0) then
	print *, "Error reading time parameters"
	print *, "aborting..."
	stop
  endif

  this%tmin = tmin
  this%tmax = tmax

end subroutine read_nml_time
!---------------------------------------------------

!---------------------------------------------------
subroutine setup_time( this, initial, dt, restart, restart_handle )
!---------------------------------------------------
!       sets up this data structure from the given information
!---------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  type( t_time ), intent( inout )  ::  this
  real(p_double), intent( in )  ::  initial, dt
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

!       local variables

  ! set initial time
  if ( restart ) then
	call restart_read( this, restart_handle )
  else
	this%t = initial
  endif

! Check to assure that this calculation will run for at least one cycle

  if ( this%tmax < ( this%t + 2 * dt ) ) then
	 this%tmax = this%t + 2 * dt
  endif

end subroutine setup_time
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_write_time( this, restart_handle )
!---------------------------------------------------
!       write object information into a restart file
!---------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  type( t_time ),  intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

!       local variables

  integer :: ierr

!       executable statements

  restart_io_wr( p_time_rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then 
	ERROR("error writing restart data for time object.")
	call abort_program(p_err_rstwrt)
  endif

  restart_io_wr( this%t, restart_handle, ierr )
  if ( ierr/=0 ) then
	ERROR('error writing restart data for time object.')
	call abort_program(p_err_rstwrt)
  endif
  restart_io_wr( this%tmin, restart_handle, ierr )
  if ( ierr/=0 ) then
	ERROR('error writing restart data for time object.')
	call abort_program(p_err_rstwrt)
  endif


end subroutine restart_write_time
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_read_time( this, restart_handle )
!---------------------------------------------------
!       read object information from a restart file
!---------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  type( t_time ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

!       local variables

  integer :: ierr
  character(len=len(p_time_rst_id)) :: rst_id
   

!       executable statements

! check if restart file is compatible
  restart_io_rd( rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then 
	 ERROR('error reading restart data for time object.')
	 call abort_program(p_err_rstrd)
  endif

  if ( rst_id /= p_time_rst_id) then
    ERROR('Corrupted restart file, or restart file')
    ERROR('from incompatible binary (time)')
	call abort_program(p_err_rstrd)
  endif

  restart_io_rd( this%t, restart_handle, ierr )
  if ( ierr/=0 ) then
	ERROR('error reading restart data for time object.')
	call abort_program(p_err_rstrd)
  endif
  restart_io_rd( this%tmin, restart_handle, ierr )
  if ( ierr/=0 ) then
	ERROR('error reading restart data for time object.')
	call abort_program(p_err_rstrd)
  endif

end subroutine restart_read_time
!---------------------------------------------------


!---------------------------------------------------
subroutine update_time( this, dt, n )
!---------------------------------------------------
!  update time using new iteration value
!---------------------------------------------------

implicit none

! dummy variables

type( t_time ), intent( inout )  ::  this

real(p_double), intent( in )  ::  dt
integer, intent(in) :: n


this%t = this%tmin + n*dt


end subroutine update_time
!---------------------------------------------------



!---------------------------------------------------
function test_time( this )
!---------------------------------------------------
!       use "this" to determine new value of "test"
!---------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  logical :: test_time

  type( t_time ), intent( in )  ::  this

!       local variables - none

!       executable statements

!         get value of result
  test_time = ( this%t < this%tmax )

end function test_time
!---------------------------------------------------


!---------------------------------------------------
function tmin( this )
!---------------------------------------------------
!       gives the minimum (starting) time
!---------------------------------------------------

  implicit none

!       dummy variables

  real(p_double)  :: tmin

  type( t_time ), intent(in) :: this

!       local variables - none

!       executable statements

  tmin = this%tmin

end function tmin
!---------------------------------------------------


!---------------------------------------------------
function tmax( this )
!---------------------------------------------------
!       gives the maximum (final) time
!---------------------------------------------------

  implicit none

!       dummy variables

  real(p_double)  :: tmax

  type( t_time ), intent(in) :: this

!       local variables - none

!       executable statements

  tmax = this%tmax

end function tmax
!---------------------------------------------------


!---------------------------------------------------
function t( this )
!---------------------------------------------------
!       gives the current time
!---------------------------------------------------

  implicit none

!       dummy variables

  real(p_double)  :: t

  type( t_time ), intent(in) :: this

!       local variables - none

!       executable statements

  t = this%t

end function t
!---------------------------------------------------


end module m_time

