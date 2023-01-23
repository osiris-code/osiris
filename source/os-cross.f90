#include "os-config.h"
#include "os-preprocess.fpp"

#ifdef __HAS_IONIZATION__

module m_cross

#include "memory/memory.h"  
  
  use m_parameters
  use m_restart

  use m_utilities

  use stringutil

  implicit none

  private

  real(p_k_fld), parameter, dimension(17) :: &
  Li_ene= (/    .9e-5, 1.0e-5, 1.3125e-5, 1.75e-5, 2e-5, 2.5e-5, &
			 6.125e-5,   1e-4,      1e-3,    1e-2, 1e-1,   5e-1, &
				 8e-1,    1e0,       1e1,     3e2, 1e10/)
				 
  real(p_k_fld), parameter, dimension(17) :: &
  Li_cross=(/  tiny(1.0e0), 1.5433e-16, 1.975e-16, 2.2091e-16, 2.2439e-16, &
		   2.2149e-16, 1.5648e-16,1.1973e-16, 2.0805e-17, 2.7214e-18, &
		   4.0665e-19, 1.8192e-19,1.6436e-19, 1.5948e-19, 1.5465e-19, &
		   1.5546e-19, 1.5546e-19 /)

  type :: t_cross

!         allow access to type components only to module procedures
	private

!         number of points in profile
	integer :: num_ene


!         allocate later: ene,cross( num_ene)
	real(p_k_fld), dimension(:), pointer :: ene => null()
	real(p_k_fld), dimension(:), pointer :: cross => null()


  end type t_cross

  interface read_nml
	module procedure read_nml_cross
  end interface

  interface setup
	module procedure setup_cross
  end interface
  
  interface cleanup
    module procedure cleanup_cross
  end interface

  interface restart_write
	module procedure restart_write_cross
  end interface

  interface restart_read
	module procedure restart_read_cross
  end interface

  interface value
	module procedure value_of_cross
  end interface

!       declare things that should be public
  public :: t_cross, cleanup, read_nml, setup
  public :: restart_write, restart_read
  public :: value


 contains 

!! add gas_type as input so that for known gases, no need to input cross-section
!---------------------------------------------------
subroutine read_nml_cross( this, input_file, gas_type )

  use m_input_file

  implicit none

  type( t_cross ), intent(out) :: this
  class( t_input_file ), intent(inout) :: input_file
  character(len=*), intent(in) :: gas_type

! maximum number of points for the function profiles used here
! necessary because of the limitations in the use of namelists
! namelist can not include variable-size data structures
  integer, parameter :: p_num_ene_max = 100 

! variables that sets the dimensionality of the profile
! and the maximum number of points 
! in case of piecewise-linear profiles
  integer  ::  num_ene
 ! character(20) :: gas_type

! ene is in the unit of Mev, cross is in unit of cm-2
! components needed for piecewise linear profiles
! Here I changed the ene and cross section to be log format to avoid repeating log computation later
  real(p_k_fld), dimension(p_num_ene_max) :: ene, cross    

   integer  ::  i

   namelist /nl_num_ene/ num_ene
   namelist /nl_cross/ cross, ene

  integer::ierr

!       executable statements

   if((lowercase(gas_type)) == "custom") then

	  num_ene = 0  

      call get_namelist( input_file, "nl_num_ene", ierr )
	  read (input_file%nml_text, nml = nl_num_ene)
				
	  if ( num_ene > 0 ) then
		call alloc( this%ene, (/ num_ene /) )
		this%ene = 0.0_p_k_fld
   
		call alloc( this%cross, (/ num_ene /) )
		this%cross = 0.0_p_k_fld
   
		ene  = - huge( 1.0_p_k_fld )
		cross = 0.0_p_k_fld
	  endif

      call get_namelist( input_file, "nl_cross", ierr )
	  read (input_file%nml_text, nml = nl_cross)
   
   else
	 num_ene=17
	 
	 call alloc( this%ene, (/ num_ene /) )
	 this%ene = 0.0_p_k_fld

	 call alloc( this%cross, (/ num_ene /) )
	 this%cross = 0.0_p_k_fld

	 ene  = - huge( 1.0_p_k_fld )
	 cross = 0.0_p_k_fld
   endif
   
   select case ( gas_type )
 
	 case ("h")

	   print *, "   Using hydrogen as neutral gas.finish later"
	   
   
	 case ("he")

	   print *, "   Using helium as neutral gas not finished yet."
		

	 case ("li")

	   print *, "   Using lithium as neutral gas."
	   ene(1:17)= Li_ene
	   cross(1:17)=Li_cross
   
   
	 case ("ar")

	   print *, "   Using argon as neutral gas. not finished yet"
	   
   
	 case ("custom")

	   print *, "   Using a customized neutral gas."
 
	 case default

	   print *, "(*error*) Not a valid neutral gas -> ", gas_type
	   print *, "Available gases:"
	   print *, "    - H"
	   print *, "    - He"                
	   print *, "    - Li"
	   print *, "    - Ar"                
	   print *, "    - Custom"
	   call abort_program()
   end select
   
   this%num_ene=num_ene
   do i=1, this%num_ene
	 this%ene(i)  =  log10( ene(i) )
	 this%cross(i) = log10( cross(i) )
   enddo

   do i=1, num_ene
	 print *, "ene",this%ene(i)  
	 print *, "cross",this%cross(i) 
   enddo
  
end subroutine read_nml_cross
!---------------------------------------------------


!---------------------------------------------------
subroutine setup_cross(this)
!---------------------------------------------------
!       setup profile for later use
!---------------------------------------------------

  implicit none

  type( t_cross ), intent( inout )  ::  this

  integer :: i

  do i =1, this%num_ene-1
     if ( this%ene(i+1) <= this%ene(i) ) then
        this%ene( i+1) = this%ene( i)
        this%cross(i+1) = this%cross(i)
     endif
  enddo

end subroutine setup_cross
!---------------------------------------------------


subroutine cleanup_cross( this )

  implicit none
  type( t_cross ), intent( inout )  ::  this
  
  call freemem( this%ene )
  call freemem( this%cross )

end subroutine cleanup_cross

!---------------------------------------------------
  real(p_k_fld) function value_of_cross( this, ex )
!---------------------------------------------------
!       Computes value at a point from a piecewise linear profile.
!       Here it is the log value of ene and cross can be treated linearly
!		the input ex is in the scale of mc^2
!---------------------------------------------------
	implicit none

!       dummy variables

	type( t_cross ), intent(in) :: this
	real(p_k_fld), intent(in) :: ex

!       local variables

	integer  ::  li
	real(p_k_fld)  :: e_log, cross_log  ! log10 of ex
	
!       executable statements

!         Locate exx between two values of this%ene, and obtain
!         value of density by linear interpolation.
!         If exx lies outside range, density=limit-density.

!         exx     position at which cross desired
!         this%ene  array of energy at which functions are specified
!         this%cross array of corresponding values

	   value_of_cross =10 ** this%cross(1)
	   e_log = log10(ex*.511)
	   if ( e_log <= this%ene(1) )  return
	   do li=2,this%num_ene
		  if ( e_log <= this%ene(li) ) then
			 cross_log     =  this%cross(li-1) +  ( this%cross(li) -  this%cross(li-1) ) / &
			                  (this%ene( li)  - this%ene( li-1) ) * ( e_log - this%ene( li-1)  )
	  
		  value_of_cross = 10**cross_log
	  
		  return
		  endif
	   enddo
	   
	   value_of_cross = 10**this%cross(this%num_ene)
	   
 end function value_of_cross
!---------------------------------------------------


!---------------------------------------------------
  subroutine restart_write_cross( this, restart_handle )
!---------------------------------------------------
!       write object information into a restart file
!---------------------------------------------------

!         <use-statements for this subroutine>

	implicit none

!       dummy variables

	type( t_cross ), intent(in) :: this
    type( t_restart_handle ), intent(inout) :: restart_handle

!       local variables

	logical :: if_x, if_cross
	integer :: ierr

!       executable statements

	if_x  = associated( this%ene  )
	if_cross = associated( this%cross )

    restart_io_wr( this%num_ene, restart_handle, ierr )
    if ( ierr/=0 ) then 
      ERROR('error writing restart data for cross-section object.')
      call abort_program( p_err_rstwrt)
    endif
    restart_io_wr( if_x, restart_handle, ierr )
    if ( ierr/=0 ) then 
      ERROR('error writing restart data for cross-section object.')
      call abort_program( p_err_rstwrt)
    endif
    restart_io_wr( if_cross, restart_handle, ierr )
    if ( ierr/=0 ) then 
      ERROR('error writing restart data for cross-section object.')
      call abort_program( p_err_rstwrt)
    endif
	
	if ( if_x  ) then
      restart_io_wr( size(this%ene), restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error writing restart data for cross-section object.')
        call abort_program( p_err_rstwrt)
      endif
      restart_io_wr( this%ene, restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error writing restart data for cross-section object.')
        call abort_program( p_err_rstwrt)
      endif
	endif
	if ( if_cross ) then 
      restart_io_wr( size(this%cross), restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error writing restart data for cross-section object.')
        call abort_program( p_err_rstwrt)
      endif
      restart_io_wr( this%cross, restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error writing restart data for cross-section object.')
        call abort_program( p_err_rstwrt)
      endif
	endif

  end subroutine restart_write_cross
!---------------------------------------------------

!---------------------------------------------------
  subroutine restart_read_cross( this, restart_handle )
!---------------------------------------------------
!       read object information from a restart file
!---------------------------------------------------

!         <use-statements for this subroutine>

	implicit none

!       dummy variables

	type( t_cross ), intent(inout) :: this
    type( t_restart_handle ), intent(in) :: restart_handle

!       local variables

	integer :: s1, ierr
	logical ::  if_x, if_cross

!       executable statements

	call cleanup(this)

	restart_io_rd( this%num_ene, restart_handle, ierr )
	if ( ierr/=0 ) then 
      ERROR('error reading restart data for cross-section object.')
      call abort_program( p_err_rstwrt)
    endif
	restart_io_rd( if_x, restart_handle, ierr )
	if ( ierr/=0 ) then 
      ERROR('error reading restart data for cross-section object.')
      call abort_program( p_err_rstwrt)
    endif
	restart_io_rd( if_cross, restart_handle, ierr )
	if ( ierr/=0 ) then 
      ERROR('error reading restart data for cross-section object.')
      call abort_program( p_err_rstwrt)
    endif

	if ( if_x  ) then
       restart_io_rd( s1, restart_handle, ierr )
       if ( ierr/=0 ) then 
         ERROR('error reading restart data for cross-section object.')
         call abort_program( p_err_rstwrt)
       endif
	   
	   call alloc( this%ene, (/s1/) )
       restart_io_rd( this%ene, restart_handle, ierr )
       if ( ierr/=0 ) then 
         ERROR('error reading restart data for cross-section object.')
         call abort_program( p_err_rstwrt)
       endif
	endif
	if ( if_cross ) then
       restart_io_rd( s1, restart_handle, ierr )
       if ( ierr/=0 ) then 
         ERROR('error reading restart data for cross-section object.')
         call abort_program( p_err_rstwrt)
       endif
	   
	   call alloc( this%cross, (/s1/) )     
       restart_io_rd( this%cross, restart_handle, ierr )
       if ( ierr/=0 ) then 
         ERROR('error reading restart data for cross-section object.')
         call abort_program( p_err_rstwrt)
       endif
	endif

  end subroutine restart_read_cross
!---------------------------------------------------


end module m_cross

#endif


