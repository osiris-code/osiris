! spline class

module m_spline

#include "memory/memory.h"

use m_system
use m_parameters

implicit none

private

type t_spline
    private
    integer,pointer,dimension(:) :: n_points => null()

    real (p_double), pointer,dimension(:) :: value_x => null(), value_fldx => null()
    real (p_double), pointer,dimension(:) :: value_y => null(), value_fldy => null()
    real (p_double), pointer,dimension(:) :: value_z => null(), value_fldz => null()

end type t_spline

interface delete
    module procedure delete_spline
end interface

interface new
  module procedure new_1d_spline
  module procedure new_2d_spline
  module procedure new_3d_spline
end interface

interface restart_write
  module procedure restart_write_spline
end interface

interface restart_read
  module procedure restart_read_spline
end interface

interface interpolate
  module procedure interpolate_1d
  module procedure interpolate_2d
  module procedure interpolate_3d
end interface

interface write_dbg
  module procedure write_dbg_spline
end interface

public :: new,delete,restart_read,restart_write,interpolate

contains


!---------------------------------------------------------------------------------------------------
subroutine restart_write_spline( this, restart_handle )
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

  implicit none

  type (t_spline), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle
  
  integer :: dim, ierr
  
  if (.not. associated(this%n_points)) then
      restart_io_wr( .FALSE., restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error writing restart data for spline object.')
        call abort_program( p_err_rstwrt)
      endif
  else
      restart_io_wr( .TRUE., restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error writing restart data for spline object.')
        call abort_program( p_err_rstwrt)
      endif
      restart_io_wr( size(this%n_points);, restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error writing restart data for spline object.')
        call abort_program( p_err_rstwrt)
      endif
      restart_io_wr( this%n_points;, restart_handle, ierr )
      if ( ierr/=0 ) then 
        ERROR('error writing restart data for spline object.')
        call abort_program( p_err_rstwrt)
      endif
      dim=size(this%n_points)  ! dim = dimension of the simulation
      select case (dim)
      case (1)
          restart_io_wr( this%value_x, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error writing restart data for spline object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_wr( this%value_fldx, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error writing restart data for spline object.')
            call abort_program( p_err_rstwrt)
          endif
      case (2)
          restart_io_wr( this%value_x, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error writing restart data for spline object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_wr( this%value_fldx, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error writing restart data for spline object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_wr( this%value_y, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error writing restart data for spline object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_wr( this%value_fldy, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error writing restart data for spline object.')
            call abort_program( p_err_rstwrt)
          endif
      case (3)
          restart_io_wr( this%value_x, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error writing restart data for spline object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_wr( this%value_fldx, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error writing restart data for spline object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_wr( this%value_y, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error writing restart data for spline object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_wr( this%value_fldy, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error writing restart data for spline object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_wr( this%value_z, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error writing restart data for spline object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_wr( this%value_fldz, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error writing restart data for spline object.')
            call abort_program( p_err_rstwrt)
          endif
      case default
      write(*,*) "(OS-SPLINE) What kind of crazy system are you running?"
      stop;
      end select
  end if
end subroutine restart_write_spline
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine restart_read_spline(this, restart_handle)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
   implicit none

   type (t_spline), intent(inout) :: this
   type( t_restart_handle ), intent(in) :: restart_handle

   integer :: dim, ierr
   logical :: flag

   restart_io_rd( flag, restart_handle, ierr )
   if ( ierr/=0 ) then 
     ERROR('error reading restart data for XXXXX object.')
     call abort_program( p_err_rstwrt)
   endif
   
   if(flag) then
       restart_io_rd( dim, restart_handle, ierr )
       if ( ierr/=0 ) then 
         ERROR('error reading restart data for XXXXX object.')
         call abort_program( p_err_rstwrt)
       endif
       call alloc(this%n_points, (/ dim /))
       restart_io_rd( this%n_points, restart_handle, ierr )
       if ( ierr/=0 ) then 
         ERROR('error reading restart data for XXXXX object.')
         call abort_program( p_err_rstwrt)
       endif

       select case (dim)
       case (1)
          call alloc(this%value_x, (/ this%n_points(1) /))
          call alloc(this%value_fldx, (/this%n_points(1)/))
          restart_io_rd( this%value_x, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error reading restart data for XXXXX object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_rd( this%value_fldx, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error reading restart data for XXXXX object.')
            call abort_program( p_err_rstwrt)
          endif
      
       case(2)
          call alloc(this%value_x, (/this%n_points(1)/))
          call alloc(this%value_fldx, (/this%n_points(1)/))
          call alloc(this%value_y, (/this%n_points(2)/))
          call alloc(this%value_fldy, (/this%n_points(2)/))
          restart_io_rd( this%value_x, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error reading restart data for XXXXX object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_rd( this%value_fldx, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error reading restart data for XXXXX object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_rd( this%value_y, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error reading restart data for XXXXX object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_rd( this%value_fldy, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error reading restart data for XXXXX object.')
            call abort_program( p_err_rstwrt)
          endif
       
       case(3)
       
          call alloc(this%value_x, (/this%n_points(1)/))
          call alloc(this%value_fldx, (/this%n_points(1)/))
          call alloc(this%value_y, (/this%n_points(2)/))
          call alloc(this%value_fldy, (/this%n_points(2)/))
          call alloc(this%value_z, (/this%n_points(3)/))
          call alloc(this%value_fldz, (/this%n_points(3)/))
          restart_io_rd( this%value_x, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error reading restart data for XXXXX object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_rd( this%value_fldx, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error reading restart data for XXXXX object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_rd( this%value_y, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error reading restart data for XXXXX object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_rd( this%value_fldy, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error reading restart data for XXXXX object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_rd( this%value_z, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error reading restart data for XXXXX object.')
            call abort_program( p_err_rstwrt)
          endif
          restart_io_rd( this%value_fldz, restart_handle, ierr )
          if ( ierr/=0 ) then 
            ERROR('error reading restart data for XXXXX object.')
            call abort_program( p_err_rstwrt)
          endif
       
       case default
         write(*,*) "(OS-SPLINE) What kind of crazy system are you running?"
         stop
       end select

   end if

end subroutine restart_read_spline	  
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine delete_spline(this)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

implicit none;

type (t_spline), intent(inout) :: this

    call freemem(this%n_points) 
    call freemem(this%value_x) 
    call freemem(this%value_fldx) 
    call freemem(this%value_y) 
    call freemem(this%value_fldy) 
    call freemem(this%value_z) 
    call freemem(this%value_fldz) 

end subroutine delete_spline
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine sort(n, position_1d, fld_1d)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
  implicit none

  integer, intent(in) :: n
  real (p_double), dimension(:), intent(inout) :: position_1d, fld_1d
      
  ! local variables
  integer :: i,j
  real (p_double):: current_min,temp_x,temp_fld
  
  do i=1,n-1
    current_min=position_1d(i)
    do j=i+1,n
       if(position_1d(j).lt.position_1d(i)) then
          temp_x=position_1d(i)
          temp_fld=fld_1d(i)
          position_1d(i)=position_1d(j)
          fld_1d(i)=fld_1d(j)
          position_1d(j)=temp_x
          fld_1d(j)=temp_fld
       end if
    end do
  end do
end subroutine sort
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine new_1d_spline(this,n,x,fldx)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

  implicit none
  type (t_spline), intent(inout) :: this
  integer,intent(in)::n
  real (p_double),dimension(:), intent(inout)::x,fldx
  
  !
  if(associated(this%n_points)) then
      call destroy(this)
  end if
  
  call alloc(this%n_points, (/1/))
  this%n_points(1)=n
  
  call alloc(this%value_x, (/n/))
  call alloc(this%value_fldx,(/n/))
  
  call sort(n, x, fldx)
  
  this%value_x=x(1:n)
  this%value_fldx=fldx(1:n)

end subroutine new_1d_spline
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine new_2d_spline(this,nx,ny,x,fldx,y,fldy)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
  implicit none
  type (t_spline), intent(inout) :: this
  integer,intent(in)::nx, ny
  real(p_double),dimension(:)::x,fldx,y,fldy

  if(associated(this%n_points)) then
    call destroy(this)
  end if

  call alloc(this%n_points, (/2/) )
  this%n_points(1)=nx
  this%n_points(2)=ny

  call alloc(this%value_x, (/nx/))
  call alloc(this%value_fldx, (/nx/))

  call alloc(this%value_y, (/ny/))
  call alloc(this%value_fldy, (/ny/))	

  call sort(nx, x, fldx)

  this%value_x=x(1:nx)
  this%value_fldx=fldx(1:nx)

  call sort(ny, y, fldy)

  this%value_y=y(1:ny)
  this%value_fldy=fldy(1:ny)

end subroutine new_2d_spline
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine new_3d_spline(this,nx,ny,nz,x,fldx,y,fldy,z,fldz)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
implicit none
  type (t_spline) :: this
  integer,intent(in)::nx, ny, nz
  real(p_double),dimension(:), intent(inout)::x,fldx,y,fldy, z, fldz

  if(associated(this%n_points)) then
    call destroy(this)
  end if

  call alloc(this%n_points, (/3/))
  this%n_points(1)=nx
  this%n_points(2)=ny
  this%n_points(3)=nz

  call alloc(this%value_x, (/nx/))
  call alloc(this%value_fldx, (/nx/))

  call alloc(this%value_y, (/ny/))
  call alloc(this%value_fldy, (/ny/))

  call alloc(this%value_z, (/nz/))
  call alloc(this%value_fldz, (/nz/))

  call sort(nx, x, fldx)

  this%value_x=x(1:nx)
  this%value_fldx=fldx(1:nx)

  call sort(ny, y, fldy)

  this%value_y=y(1:ny)
  this%value_fldy=fldy(1:ny) 

  call sort(nz, z, fldz)

  this%value_z=z(1:nz)
  this%value_fldz=fldz(1:nz) 

end subroutine new_3d_spline
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function interpolate_basic(value_x, fldx ,x) result (c)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

  implicit none
   
  real(p_double),intent(in)::x
  real(p_double), dimension(:), pointer :: value_x, fldx
  real(p_double) :: c
  
  ! local variables     
  integer i,n
  real (p_double) :: norm  
  
  n=size(value_x)
  if (n.eq.0) then
    c=0
    return
  end if
  
  if(x.le.value_x(1)) then
    c=fldx(1)
    return
  end if
  
  if(x.ge.value_x(n)) then
    c=fldx(n)
    return
  end if
  
  do i=1,n
      if(x.lt.value_x(i)) then
        norm = 1.0/(value_x(i-1)-value_x(i));
        c =norm*((fldx(i-1)-fldx(i))*x + value_x(i-1)*fldx(i) - value_x(i)*fldx(i-1))
        return
      end if
  end do

end function interpolate_basic
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function interpolate_1d(this,x) result (c)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

implicit none

  type (t_spline),intent(in)::this;
  real (p_double),intent(in)::x;
  real (p_double) :: c;
  
  !    if this%n_points is not associated, result=0     
  
  if(.not. associated(this%n_points)) then
    c=0.0
    return
  end if
  
  c=interpolate_basic(this%value_x, this%value_fldx, x)

end function interpolate_1d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function interpolate_2d(this,x,y) result (c)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

  implicit none
  
  type (t_spline),intent(in):: this
  real(p_double),intent(in) :: x,y
  real (p_double) :: c
  
  !    local variables     
  real (p_double) :: cx,cy
  !    if this%n_points is not associated, result=0     
  
  if(.not.associated(this%n_points)) then
    c=0.0
    return
  end if
  
  cx=interpolate_basic(this%value_x, this%value_fldx, x)
  cy=interpolate_basic(this%value_y, this%value_fldy, y)
  c=cx*cy

end function interpolate_2d 
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function interpolate_3d(this,x,y,z) result (c)
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

  implicit none
  
  type (t_spline),intent(in)::this
  real (p_double),intent(in)::x,y,z
  real (p_double) :: c
  
  !    local variables     
  real (p_double) :: cx,cy,cz
  !    if this%n_points is not associated, result=0     
  
  if(.not.associated(this%n_points)) then
    c=0.0
    return
  end if
  
  cx=interpolate_basic(this%value_x, this%value_fldx, x)
  cy=interpolate_basic(this%value_y, this%value_fldy, y)
  cz=interpolate_basic(this%value_z, this%value_fldz, z)
  c=cx*cy*cz

end function interpolate_3d    
!---------------------------------------------------------------------------------------------------


end module m_spline  
