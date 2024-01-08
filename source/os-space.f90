!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     space class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_space

use m_system
use m_restart
use m_utilities
use m_parameters

implicit none

private

!       string to id restart data
character(len=*), parameter :: p_space_rst_id = "space rst data - 0x0004"

type :: t_space
 
   ! number of dimensions of space object
   integer :: x_dim = 0
 
 !         first general information on the simulation box that
 !         is independent on how the boundaries of the box move
 !         over longer periods of time
 
 !         space in simulation box
   real(p_double), dimension(2,p_max_dim) :: x_bnd = 0.0
 
 !         space in simulation box before last call of
 !         move_boundary
   real(p_double), dimension(2,p_max_dim) :: x_bnd_old = 0.0
 
 !         number of dx the boundaries of the space move at
 !         current timestep
 !         the sign indicates the direction of the motion of
 !         each boundary
   integer, dimension(2,p_max_dim) :: move_num = 0
 
 !         distances moved since last move of space boundaries
   real(p_double), dimension(2,p_max_dim) :: xmoved = 0.0
 
   ! number of cells moved, for calculating sim. boundaries with less roundoff 
   ! on moving window runs
   integer, dimension( p_max_dim ) :: total_i_moved = 0
 
 !         the following variable(s) contain information about
 !         the specific motion of the simulation box
 !         at this time only motion with c is implemented
 
 !         information on motion of the simulation box
   logical, dimension(p_max_dim) :: if_move = .false.
   real(p_double) :: move_vel = 0.0
   
 
   ! original positions of simulation box at first time step
   real(p_double), dimension(2,p_max_dim) :: x_bnd_initial = 0.0
  
end type t_space


interface read_nml
  module procedure read_nml_space
end interface

interface setup
  module procedure setup_space_global
end interface

interface restart_write
  module procedure restart_write_space
end interface

interface restart_read
  module procedure restart_read_space
end interface

interface move_boundary
  module procedure move_boundary_space_global
end interface

interface xmin
  module procedure xmin_space
  module procedure xmin_dim_space
end interface

interface xmax
  module procedure xmax_space
  module procedure xmax_dim_space
end interface

interface xmin_initial
  module procedure xmin_initial
end interface

interface xmax_initial
  module procedure xmax_initial
end interface

interface x_bnd
  module procedure x_bnd
  module procedure x_bnd_dim          
end interface

interface get_x_bnd
  module procedure get_x_bnd_single
  module procedure get_x_bnd_double
  module procedure get_x_bnd_dim
end interface

interface extend
  module procedure extend_of_space
end interface

interface xmin_old
  module procedure xmin_old
end interface

interface xmax_old
  module procedure xmax_old
end interface

interface if_move
  module procedure if_move_space
  module procedure if_move_dim_space
end interface

interface if_move_any
  module procedure if_move_any
end interface

interface nx_move
  module procedure nx_move_all
  module procedure nx_move_dim
end interface

interface total_xmoved
  module procedure total_xmoved_all
  module procedure total_xmoved_dim_bnd
end interface

interface total_imoved
  module procedure total_imoved_dim
end interface

interface if_center
  module procedure if_center
end interface

interface reshape_obj
  module procedure reshape_space
end interface

interface x_dim
  module procedure x_dim_space
end interface 

interface transpose_obj
  module procedure transpose_space
end interface

interface set_x_bnd
  module procedure set_x_bnd_single
  module procedure set_x_bnd_double
end interface


!       declare things that should be public
public :: t_space, read_nml, setup
public :: restart_write, restart_read, move_boundary
public :: xmin, xmax, x_bnd, get_x_bnd, extend, set_x_bnd
public :: xmin_old, xmax_old, if_move 
public :: if_move_any, if_center, reshape_obj
public :: xmin_initial, xmax_initial, total_xmoved, total_imoved
public :: x_dim, nx_move
public :: transpose_obj

contains 


!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine read_nml_space( this, input_file, periodic, coordinates )

  use m_input_file

  implicit none

!       dummy variables

  type( t_space ), intent(out) :: this
  class( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic
  integer, intent(in) :: coordinates

  real(p_double), dimension(p_max_dim) :: xmin, xmax
  logical,         dimension(p_max_dim) :: if_move
  real(p_double) :: move_u

  namelist /nl_space/ xmin, xmax, if_move, move_u
  
  integer :: i, ierr


  ! this will change
  this%x_dim = p_x_dim          

  xmin = 0.0_p_double
  xmax = 0.0_p_double
  
  if_move = .false.
  move_u = -1.0_p_double

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_space", ierr )

  if (ierr /= 0) then
	if (ierr < 0) then
	  print *, "Error reading space parameters"
	else 
	  print *, "Error: space parameters missing"
	endif
	print *, "aborting..."
	stop
  endif

  read (input_file%nml_text, nml=nl_space, iostat = ierr)
  if (ierr /= 0) then
	print *, "Error reading space parameters"
	print *, "aborting..."
	stop
  endif

  ! check dimensions
  do i = 1, p_x_dim
	if (xmin(i) >= xmax(i)) then
	   print *, "Invalid space parameters, xmin must be < xmax ",& 
				  "for all dimensions."
	   print *, " xmin = ", xmin
	   print *, " xmax = ", xmax
	   print *, "aborting..."
	   stop
	endif
  enddo

  ! check moving window
  do i = 1, p_x_dim
	if ( periodic(i) .and. if_move(i) ) then
	   print *, "Invalid space parameters, you cannot have window ",& 
				  "moving along a periodic direction."
	   print *, "( direction ", i, " )."
	   print *, "aborting..."
	   stop
	endif
  enddo
  
  if ( coordinates == p_cylindrical_b ) then 
	
	! check axial boundary limit
	if ( xmin( p_r_dim ) /= 0.0_p_double ) then
	   print *, "Invalid space parameters, when using cylindrical geometry xmin(2) ", &
				  "(radial coordinate minimum) must be set to 0.0."
	   print *, " xmin(2) = ", xmin(2)
	   print *, "aborting..."
	   stop
	endif
	
	! check moving window
	if ( if_move( p_r_dim ) ) then
	   print *, "Invalid space parameters, when using cylindrical geometry you cannot ", &
				  "have window moving along radial coordinate (x2)"
	   print *, "aborting..."
	   stop
	endif
  endif
  
  this%x_bnd(p_lower,1:this%x_dim)   = xmin(1:this%x_dim)
  this%x_bnd(p_upper,1:this%x_dim)   = xmax(1:this%x_dim)
			
  this%if_move(1:this%x_dim)   = if_move(1:this%x_dim)
  this%x_bnd_initial(p_lower,1:this%x_dim)  = xmin(1:this%x_dim)
  this%x_bnd_initial(p_upper,1:this%x_dim)  = xmax(1:this%x_dim)
  
  if ( move_u > 0.0_p_double ) then
	this%move_vel = move_u/sqrt( 1.0_p_double + move_u**2 )
  else
	this%move_vel = 0.0_p_double
  endif
  
end subroutine read_nml_space
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
subroutine setup_space_global( this, dx, coordinates, restart, restart_handle )
!---------------------------------------------------------------------------------------------------
!       sets up of global space of simulation
!---------------------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  type( t_space ), intent( inout ) :: this

  real(p_double), dimension(:), intent(in) :: dx
  integer ,               intent(in) :: coordinates
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

!       local variables

!       executable statements

  if (restart) then

	 call restart_read( this, restart_handle )

  else

	 ! For cylindrical coordinates the space should start at -0.5*dr
	 if ( coordinates == p_cylindrical_b ) then
		this%x_bnd(p_lower,p_r_dim) = this%x_bnd(p_lower,p_r_dim)-dx(p_r_dim)*0.5_p_double
		this%x_bnd(p_upper,p_r_dim) = this%x_bnd(p_upper,p_r_dim)-dx(p_r_dim)*0.5_p_double
	 endif

	 this%move_num   = 0
	 this%xmoved     = 0.0
	 this%x_bnd_old  = this%x_bnd
	 
	 this%total_i_moved = 0
  endif
			
end subroutine setup_space_global
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
        subroutine restart_write_space( this, restart_handle )
!---------------------------------------------------------------------------------------------------
!       write object information into a restart file
!---------------------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

          implicit none

!       dummy variables

          type( t_space ), intent(in) :: this
          type( t_restart_handle ), intent(inout) :: restart_handle

!       local variables

          integer :: ierr

!       executable statements

          restart_io_wr( p_space_rst_id, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error writing restart data for space object.")
            call abort_program(p_err_rstwrt)
          endif
          
          
          restart_io_wr( this%x_dim, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error writing restart data for space object.")
            call abort_program(p_err_rstwrt)
          endif
          restart_io_wr( this%x_bnd, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error writing restart data for space object.")
            call abort_program(p_err_rstwrt)
          endif
          restart_io_wr( this%x_bnd_old, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error writing restart data for space object.")
            call abort_program(p_err_rstwrt)
          endif
          restart_io_wr( this%move_num, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error writing restart data for space object.")
            call abort_program(p_err_rstwrt)
          endif
          restart_io_wr( this%xmoved, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error writing restart data for space object.")
            call abort_program(p_err_rstwrt)
          endif
          restart_io_wr( this%if_move, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error writing restart data for space object.")
            call abort_program(p_err_rstwrt)
          endif
          restart_io_wr( this%move_vel, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error writing restart data for space object.")
            call abort_program(p_err_rstwrt)
          endif
          restart_io_wr( this%x_bnd_initial, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error writing restart data for space object.")
            call abort_program(p_err_rstwrt)
          endif
          restart_io_wr( this%total_i_moved, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error writing restart data for space object.")
            call abort_program(p_err_rstwrt)
          endif

        end subroutine restart_write_space
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
        subroutine restart_read_space( this, restart_handle )
!---------------------------------------------------------------------------------------------------
!       read object information from a restart file
!---------------------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

          implicit none

!       dummy variables

          type( t_space ), intent(inout) :: this
          type( t_restart_handle ), intent(in) :: restart_handle

!       local variables
           
          integer :: ierr 
          character(len=len(p_space_rst_id)) :: rst_id
           
!       executable statements

		  ! check if restart file is compatible
		  restart_io_rd( rst_id, restart_handle, ierr )
		  if ( ierr/=0 ) then 
			 ERROR('error reading restart data for space object.')
			 call abort_program(p_err_rstrd)
		  endif
	  
		  if ( rst_id /= p_space_rst_id) then
ERROR('Corrupted restart file, or restart file')
ERROR('from incompatible binary (space)')
			call abort_program(p_err_rstrd)
		  endif

          restart_io_rd( this%x_dim, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error reading restart data for space object.")
            call abort_program(p_err_rstrd)
          endif
          restart_io_rd( this%x_bnd, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error reading restart data for space object.")
            call abort_program(p_err_rstrd)
          endif
          restart_io_rd( this%x_bnd_old, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error reading restart data for space object.")
            call abort_program(p_err_rstrd)
          endif
          restart_io_rd( this%move_num, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error reading restart data for space object.")
            call abort_program(p_err_rstrd)
          endif
          restart_io_rd( this%xmoved, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error reading restart data for space object.")
            call abort_program(p_err_rstrd)
          endif
          restart_io_rd( this%if_move, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error reading restart data for space object.")
            call abort_program(p_err_rstrd)
          endif
          restart_io_rd( this%move_vel, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error reading restart data for space object.")
            call abort_program(p_err_rstrd)
          endif
          restart_io_rd( this%x_bnd_initial, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error reading restart data for space object.")
            call abort_program(p_err_rstrd)
          endif
          restart_io_rd( this%total_i_moved, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR("error reading restart data for space object.")
            call abort_program(p_err_rstrd)
          endif

        end subroutine restart_read_space
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine move_boundary_space_global( this, dx, dt )
!---------------------------------------------------------------------------------------------------
!       tests whether space has to move by dx 
!       and move it accordingly
!---------------------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  type( t_space ), intent( inout )  ::  this

  real(p_double), dimension(:), intent(in) :: dx
  real(p_double),               intent(in) :: dt

!       local variables

  integer :: i_dim

!       executable statements

!         move boundaries of space

  this%x_bnd_old = this%x_bnd

  do i_dim=1, this%x_dim
	 if ( this%if_move(i_dim) ) then
	     
		if ( this%move_vel > 0.0_p_double ) then
		  ! move boundary at specified velocity
		  this%xmoved(:,i_dim)   = this%xmoved(:,i_dim) + this%move_vel*dt
		else 
		  ! move boundary at c
		  this%xmoved(:,i_dim)   = this%xmoved(:,i_dim) + dt
		endif
		
		if ( this%xmoved(1,i_dim) >= dx(i_dim) ) then
		   this%move_num(:,i_dim) = 1
		   this%xmoved(:,i_dim)   = this%xmoved(:,i_dim) - dx(i_dim)
		   
		   ! old version, roundoff problems for very long runs
		   !this%x_bnd(:,i_dim)    = this%x_bnd(:,i_dim) + dx(i_dim)
		   
		   ! new version, no roundoff
		   this%total_i_moved(i_dim) = this%total_i_moved(i_dim) + 1
		   this%x_bnd(:,i_dim) = this%x_bnd_initial(:,i_dim) + this%total_i_moved(i_dim) * dx(i_dim) 
		else
		   this%move_num(:,i_dim) = 0
		endif
	 else
		this%move_num(:,i_dim) = 0
	 endif
  enddo

end subroutine move_boundary_space_global
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function xmin_space( this )
!---------------------------------------------------------------------------------------------------
! gives lower boundaries of space
!---------------------------------------------------------------------------------------------------
implicit none

  type( t_space ), intent( in )  ::  this
  
  ! this breaks the intel ifort 9.1 compiler
  !real(p_double), dimension(this%x_dim) :: xmin_space
  
  real(p_double), dimension(p_x_dim) :: xmin_space

  xmin_space(1:this%x_dim) = this%x_bnd(p_lower,1:this%x_dim)

end function xmin_space
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
pure function xmin_dim_space( this, dim )
!---------------------------------------------------------------------------------------------------
! gives lower boundaries of space
!---------------------------------------------------------------------------------------------------
implicit none

  type( t_space ), intent( in )  ::  this
  integer, intent(in) :: dim
  real(p_double) :: xmin_dim_space

  xmin_dim_space = this%x_bnd(p_lower,dim)

end function xmin_dim_space
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function xmax_space( this )
!---------------------------------------------------------------------------------------------------
! gives upper boundaries of space
!---------------------------------------------------------------------------------------------------

  implicit none

  type( t_space ), intent( in )  ::  this
  
  ! this breaks the intel ifort 9.1 compiler
  !real(p_double), dimension(this%x_dim) :: xmax_space
  
  real(p_double), dimension(p_x_dim) :: xmax_space

  xmax_space(1:this%x_dim) = this%x_bnd(p_upper,1:this%x_dim)

end function xmax_space
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
pure function xmax_dim_space( this, dim )
!---------------------------------------------------------------------------------------------------
! gives upper boundaries of space
!---------------------------------------------------------------------------------------------------

  implicit none

  type( t_space ), intent( in )  ::  this
  integer, intent(in) :: dim
  real(p_double) :: xmax_dim_space

  xmax_dim_space = this%x_bnd(p_upper,dim)

end function xmax_dim_space
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function xmin_initial( this )
!---------------------------------------------------------------------------------------------------
!       gives lower boundaries of space
!---------------------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  type( t_space ), intent( in )  ::  this
  real(p_double), dimension(this%x_dim) :: xmin_initial

!       local variables

!       executable statements

  xmin_initial(1:this%x_dim) = this%x_bnd_initial(p_lower,1:this%x_dim)

end function xmin_initial
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function xmax_initial( this )
!---------------------------------------------------------------------------------------------------
!       gives upper boundaries of space
!---------------------------------------------------------------------------------------------------

  implicit none

!       dummy variables


  type( t_space ), intent( in )  ::  this
  real(p_double), dimension(this%x_dim) :: xmax_initial

!       local variables

!       executable statements

  xmax_initial(1:this%x_dim) = this%x_bnd_initial(p_upper,1:this%x_dim)

end function xmax_initial
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function x_bnd( this )
!---------------------------------------------------------------------------------------------------
!       gives boundaries of space
!---------------------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  type( t_space ), intent( in )  ::  this
 
 ! this breaks the intel ifort 9.1 compiler
 ! real(p_double), dimension(2, this%x_dim) :: x_bnd

 real(p_double), dimension(2, p_x_dim) :: x_bnd

!       local variables

!       executable statements

  x_bnd(:,1:this%x_dim) = this%x_bnd(:,1:this%x_dim)

end function x_bnd
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function x_bnd_dim( this , dim)
!---------------------------------------------------------------------------------------------------
!       gives boundaries of space along given direction
!---------------------------------------------------------------------------------------------------

  implicit none

  type( t_space ), intent( in )  ::  this
  integer, intent(in) :: dim
 
  real(p_double), dimension(2) :: x_bnd_dim

  x_bnd_dim = this%x_bnd(:,dim)

end function x_bnd_dim
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine set_x_bnd_single( this, x_bnd_single )
!---------------------------------------------------------------------------------------------------

   implicit none

   type( t_space ), intent( inout )  ::  this
   real(p_single), dimension(:,:), intent(in) :: x_bnd_single

   this%x_bnd(1:2,1:this%x_dim) = x_bnd_single(1:2,1:this%x_dim) 

end subroutine set_x_bnd_single
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine set_x_bnd_double( this, x_bnd_double )
!---------------------------------------------------------------------------------------------------
!       gives boundaries of space
!---------------------------------------------------------------------------------------------------

   implicit none

   type( t_space ), intent( inout )  ::  this
   real(p_double), dimension(:,:), intent(in) :: x_bnd_double

   this%x_bnd(1:2,1:this%x_dim) = x_bnd_double(1:2,1:this%x_dim)

end subroutine set_x_bnd_double
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine get_x_bnd_single( this, x_bnd_single )
!---------------------------------------------------------------------------------------------------
!       gives boundaries of space
!---------------------------------------------------------------------------------------------------

   implicit none

   type( t_space ), intent( in )  ::  this
   real(p_single), dimension(:,:), intent(out) :: x_bnd_single

   x_bnd_single(1:2,1:this%x_dim) = real( this%x_bnd(1:2,1:this%x_dim), p_single )

end subroutine get_x_bnd_single
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine get_x_bnd_double( this, x_bnd_double )
!---------------------------------------------------------------------------------------------------
!       gives boundaries of space
!---------------------------------------------------------------------------------------------------

   implicit none

   type( t_space ), intent( in )  ::  this
   real(p_double), dimension(:,:), intent(out) :: x_bnd_double

   x_bnd_double(1:2,1:this%x_dim) = this%x_bnd(1:2,1:this%x_dim)

end subroutine get_x_bnd_double
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine get_x_bnd_dim( this, dim, x_bnd )
!---------------------------------------------------------------------------------------------------
!       gives boundaries of space for given direction
!---------------------------------------------------------------------------------------------------

   implicit none

!       dummy variables

   type( t_space ), intent( in )  ::  this
   integer, intent( in ) :: dim
   real(p_double), dimension(2), intent(out) :: x_bnd

!       local variables

!       executable statements

   x_bnd(1:2) = this%x_bnd(1:2,dim)

end subroutine get_x_bnd_dim
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function extend_of_space( this )
!---------------------------------------------------------------------------------------------------
!       gives the extend of the space
!---------------------------------------------------------------------------------------------------

  implicit none

  type( t_space ), intent( in )  ::  this
  
  ! this breaks the intel ifort 9.1 compiler
  !real(p_double), dimension(1:this%x_dim) :: extend_of_space
  
  real(p_double), dimension(p_x_dim) :: extend_of_space


  extend_of_space(1:this%x_dim) = this%x_bnd(p_upper,1:this%x_dim) - &
								  this%x_bnd(p_lower,1:this%x_dim)

end function extend_of_space
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function xmin_old( this )
!---------------------------------------------------------------------------------------------------
!       gives lower boundaries of space before last move
!---------------------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  type( t_space ), intent( in )  ::  this
  real(p_double), dimension(this%x_dim) :: xmin_old

!       local variables

!       executable statements

  xmin_old(1:this%x_dim) = this%x_bnd_old(p_lower,1:this%x_dim)

end function xmin_old
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function xmax_old( this )
!---------------------------------------------------------------------------------------------------
!       gives upper boundaries of space before last timestep
!---------------------------------------------------------------------------------------------------

  implicit none

!       dummy variables
  type( t_space ), intent( in )  ::  this
  real(p_double), dimension(this%x_dim) :: xmax_old

!       local variables

!       executable statements

  xmax_old(1:this%x_dim) = this%x_bnd_old(p_upper,1:this%x_dim)

end function xmax_old
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function nx_move_all( this )
!---------------------------------------------------------------------------------------------------
! number of dx space boundaries move now - in this timestep
!---------------------------------------------------------------------------------------------------

  implicit none

  type( t_space ), intent( in )  ::  this
  integer, dimension( this%x_dim ) :: nx_move_all
  
  nx_move_all(1:this%x_dim) = this%move_num( 1, 1:this%x_dim )
  
end function nx_move_all
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function nx_move_dim( this, dim )
!---------------------------------------------------------------------------------------------------
! number of dx space boundaries move now - in this timestep
!---------------------------------------------------------------------------------------------------

  implicit none

  type( t_space ), intent( in )  ::  this
  integer, intent(in) :: dim
  integer :: nx_move_dim
  
  nx_move_dim = this%move_num( 1, dim )
  
end function nx_move_dim
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
pure function x_dim_space( this )
!---------------------------------------------------------------------------------------------------
! number of spatial dimensions of this object
!---------------------------------------------------------------------------------------------------

  implicit none

  integer :: x_dim_space
  type( t_space ), intent( in )  ::  this

  x_dim_space = this%x_dim
  
end function x_dim_space
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function total_xmoved_all( this )
!---------------------------------------------------------------------------------------------------
!       gives total motion of space since initial timestep
!---------------------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  type( t_space ), intent( in )  ::  this
  real(p_double), dimension(2,this%x_dim) :: total_xmoved_all

!       local variables

!       executable statements

  total_xmoved_all(:,1:this%x_dim) = this%x_bnd(:,1:this%x_dim) - &
								 this%x_bnd_initial(:,1:this%x_dim)

end function total_xmoved_all
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function total_xmoved_dim_bnd( this, dim , bnd)
!---------------------------------------------------------------------------------------------------
!       gives total motion of upper boundary since initial timestep
!---------------------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  type( t_space ), intent( in )  ::  this
  integer, intent(in) :: dim, bnd
  real(p_double) :: total_xmoved_dim_bnd

!       local variables

!       executable statements

  total_xmoved_dim_bnd = this%x_bnd(bnd,dim) - &
								 this%x_bnd_initial(bnd, dim)

end function total_xmoved_dim_bnd
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function total_imoved_dim( this, dim )

  implicit none

  type( t_space ), intent( in )  ::  this
  integer, intent(in) :: dim
  integer :: total_imoved_dim

  total_imoved_dim = this%total_i_moved( dim ) 

end function total_imoved_dim
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function if_move_space( this )
!---------------------------------------------------------------------------------------------------
! gives information about general motion of whole space
!---------------------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  type( t_space ), intent( in )  ::  this

! this breaks the intel ifort 9.1 compiler
!  logical, dimension(this%x_dim) :: if_move_space
  
  logical, dimension(p_x_dim) :: if_move_space

  if_move_space(1:this%x_dim) = this%if_move(1:this%x_dim)

end function if_move_space
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
pure function if_move_dim_space( this, dim )
!---------------------------------------------------------------------------------------------------
! gives information about general motion of whole space
!---------------------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  type( t_space ), intent( in )  ::  this
  integer, intent(in) :: dim
  
  logical :: if_move_dim_space

  if_move_dim_space = this%if_move(dim)

end function if_move_dim_space
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
        function if_move_any( this )
!---------------------------------------------------------------------------------------------------
!       tests if space is moving along any coordinate
!---------------------------------------------------------------------------------------------------

          implicit none

!       dummy variables

          logical :: if_move_any

          type( t_space ), intent( in )  ::  this

!       local variables
          
           integer :: i

!       executable statements
          
          if_move_any = this%if_move(1)
          do i = 2, this%x_dim
            if_move_any = if_move_any .or. this%if_move(i)
          enddo
          
        end function if_move_any
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
        function if_center( this, coordinates )
!---------------------------------------------------------------------------------------------------
!       gives information whether r = 0 is part of this space object
!       if cylindrical coordinates are used
!---------------------------------------------------------------------------------------------------

          implicit none

!       dummy variables

          logical :: if_center

         type( t_space ), intent(in) :: this
         integer, intent(in) :: coordinates

!       local variables


!       executable statements

          if ( coordinates == p_cylindrical_b ) then
             if ( this%x_bnd(1,p_r_dim) <= 0.0_p_double ) then
                if_center = .true.
             else
                if_center = .false.
             endif
          else
             if_center = .false.
          endif

        end function if_center
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine reshape_space( this, deltax, dx )
!---------------------------------------------------------------------------------------------------
! Change the shape of the space object by adding/removing dx cells on each boundary
!---------------------------------------------------------------------------------------------------
  implicit none
  
  type( t_space ), intent(inout) :: this
  integer, dimension(:,:), intent(in) :: deltax
  real(p_double), dimension(:), intent(in) :: dx
  
  
  this%x_bnd(1,1:this%x_dim) = this%x_bnd(1,1:this%x_dim) + &
                                        deltax(1,:)*dx(1:this%x_dim)
  this%x_bnd(2,1:this%x_dim) = this%x_bnd(2,1:this%x_dim) + &
                                        deltax(2,:)*dx(1:this%x_dim)
  
end subroutine reshape_space
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine transpose_space( self, dir1, dir2 )
!---------------------------------------------------------------------------------------------------
! Transpose spatial coordinates
!---------------------------------------------------------------------------------------------------
  implicit none
  
  type( t_space ), intent(inout) :: self
  integer, intent(in) :: dir1, dir2

  call swap_value( self%x_bnd, 2, dir1, dir2 )
  call swap_value( self%x_bnd_old, 2, dir1, dir2 )
  
  call swap_value( self%move_num, 2, dir1, dir2 )
  call swap_value( self%xmoved, 2, dir1, dir2 )

  call swap_value( self%if_move, dir1, dir2 )
  call swap_value( self%x_bnd_initial, 2, dir1, dir2 )
  
end subroutine transpose_space
!---------------------------------------------------------------------------------------------------


end module m_space

