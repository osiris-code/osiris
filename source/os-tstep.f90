!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     time step class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_time_step

use m_parameters
use m_restart
use m_utilities

implicit none

!       restrict access to things explicitly declared public
private


type :: t_time_step

  ! allow access to type components only to module procedures
  private

  real(p_double) :: dt    ! timestep length
  integer :: n            ! timestep counter
  integer :: ndump        ! "master" timestep for reports/data dumps

  integer :: dump_start   ! Initial iteration for diagnostics

end type t_time_step


interface read_nml
  module procedure read_nml_time_step
end interface

interface setup
  module procedure setup_time_step
end interface

interface restart_write
  module procedure restart_write_time_step
end interface

interface restart_read
  module procedure restart_read_time_step
end interface

interface advance
  module procedure advance_time_step
end interface

interface dt
  module procedure dt_time_step
end interface

interface n
  module procedure n_time_step
end interface

interface ndump
  module procedure ndump_time_step
end interface

interface test_if_report
  module procedure test_if_report
end interface

!       declare things that should be public
public :: t_time_step, read_nml, setup
public :: restart_write, advance
public :: dt, n, ndump, test_if_report


contains 


!---------------------------------------------------
!---------------------------------------------------
subroutine read_nml_time_step( this, input_file )

  use m_input_file

  implicit none

  type( t_time_step ), intent(out) :: this
  class( t_input_file ), intent(inout) :: input_file

  real(p_double)  ::  dt
  integer  ::  ndump, dump_start

  namelist /nl_time_step/ dt, ndump, dump_start
  
  integer :: ierr

!       executable statements

  dt = 0.0_p_double
  ndump = 0
  dump_start = 0

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_time_step", ierr )
  
  if ( ierr /= 0) then
	if (ierr < 0) then
	  print *, "Error reading time_step parameters"
	else 
	  print *, "Error time_step parameters missing"
	endif
	print *, "aborting..."
	stop
  endif

  
  read (input_file%nml_text, nml = nl_time_step, iostat = ierr)
  if (ierr /= 0) then
	print *, "Error reading time_step parameters"
	print *, "aborting..."
	stop
  endif

  if (dt == 0.0) then
	print *, "Error: dt = 0.0"
	print *, "aborting..."
	stop
  endif

  this%dt         =  dt
  this%ndump      =  ndump
  this%dump_start =  dump_start

end subroutine read_nml_time_step
!---------------------------------------------------

!---------------------------------------------------
subroutine setup_time_step(this, initial, restart, restart_handle)
!---------------------------------------------------
!       sets up this data structure from the given information
!---------------------------------------------------

  implicit none

!       dummy variables

  type( t_time_step ), intent(inout) :: this

  integer, intent(in) :: initial
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

!       local variables

  if ( restart ) then
	call restart_read( this, restart_handle )
  else
	this%n = initial
  endif

end subroutine setup_time_step
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_write_time_step( this, restart_handle )
!---------------------------------------------------
!       write object information into a restart file
!---------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  type( t_time_step ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

!       local variables
  
  integer :: ierr

!       executable statements

  restart_io_wr( this%dt, restart_handle, ierr )
  if ( ierr/=0 ) then 
	ERROR('error writing restart data for time_step object.')
	call abort_program(p_err_rstwrt)
  endif
  restart_io_wr( this%n, restart_handle, ierr )
  if ( ierr/=0 ) then 
	ERROR('error writing restart data for time_step object.')
	call abort_program(p_err_rstwrt)
  endif


end subroutine restart_write_time_step
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_read_time_step( this, restart_handle )
!---------------------------------------------------
!       read object information from a restart file
!---------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  type( t_time_step ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

!       local variables

  integer :: ierr

!       executable statements
  
  restart_io_rd( this%dt, restart_handle, ierr )
  if ( ierr/=0 ) then 
    ERROR('Error reading time_step restart information, iostat = ', ierr)
	call abort_program(p_err_rstwrt)
  endif
  restart_io_rd( this%n, restart_handle, ierr )
  if ( ierr/=0 ) then 
    ERROR('Error reading time_step restart information, iostat = ', ierr)
	call abort_program(p_err_rstwrt)
  endif

end subroutine restart_read_time_step
!---------------------------------------------------



!---------------------------------------------------
subroutine advance_time_step( this, dn )
!---------------------------------------------------
!       advance time step counter by dn
!---------------------------------------------------

  implicit none

!       dummy variables

  type( t_time_step ), intent( inout )  :: this

  integer, intent(in) :: dn

!       local variables - none

!       executable statements


!         increment time step counter 
!         n      master time step counter
!         dn     step increment
  this%n = this%n + dn


end subroutine advance_time_step
!---------------------------------------------------


!---------------------------------------------------
function dt_time_step( this )
!---------------------------------------------------
!       gives size of timestep
!---------------------------------------------------

  implicit none

!       dummy variables

  real(p_double)  :: dt_time_step

  type( t_time_step ), intent(in) :: this

!       local variables - none

!       executable statements

  dt_time_step = this%dt

end function dt_time_step
!---------------------------------------------------


!---------------------------------------------------
pure function n_time_step( this )
!---------------------------------------------------
!       gives time step counter
!---------------------------------------------------

  implicit none

!       dummy variables

  integer  :: n_time_step

  type( t_time_step ), intent(in) :: this

!       local variables - none

!       executable statements

  n_time_step = this%n

end function n_time_step
!---------------------------------------------------


!---------------------------------------------------
function ndump_time_step( this )
!---------------------------------------------------
!       gives "master" time step number for data dumps/reports
!---------------------------------------------------

  implicit none

!       dummy variables

  integer  :: ndump_time_step

  type( t_time_step ), intent(in) :: this

!       local variables - none

!       executable statements

  ndump_time_step = this%ndump

end function ndump_time_step
!---------------------------------------------------

!---------------------------------------------------
function test_if_report( this, ndump_fac, iter_ )
!---------------------------------------------------
! test if a report should be done at this timestep
!---------------------------------------------------

  implicit none
 
  logical                     :: test_if_report
  
  type( t_time_step ), intent(in) :: this
  integer, intent(in) :: ndump_fac
  integer, intent(in), optional :: iter_
 
  integer :: n, n_del, n_aux
  
  n = this%n
  if ( present( iter_ ) ) n = n + iter_
  
  n_del = this%ndump * ndump_fac
  
  if ( ( n >= this%dump_start ) .and. ( n_del > 0 ) ) then
	n_aux = n / n_del
	test_if_report = (( n - n_aux*n_del ) == 0) 
  else
    test_if_report = .false.
  endif

end function test_if_report
!---------------------------------------------------

end module m_time_step

