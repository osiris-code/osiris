!! os-lindman.f --> lindman module for OSIRIS
!!  by Frank S. Tsung
!!  version 1.0 -- 9/20/99
!!  Lindman absorbing B.C. can be applied in only one direction.
!!  I guess I am assuming that the perpendicular directions are
!!  periodic for simplicity (I am not familiar with all of OSIRIS
!!  B.C.).
!!
!!
!!  1/10/2000 --> added dtdxi = dt/dx to advance the magnetic fields
!!  on the boundary via your good old Maxwell equation. (FST)
!!
!!
!!  3/10/2000 --> re-write the Lindman B.C. Using the Wall object
!!                implemented by Roy.
!!
!!
!! 12/5/01 --> Summary of Changes
!!             No major changes
!!            (b0 and e0 are no longer needed)

#include "os-config.h"
#include "os-preprocess.fpp"

module m_lindman

use m_emf_define
use m_restart
use m_space
use m_node_conf
use m_grid_define
use m_vdf_define
use m_wall_comm
use m_wall_define
use m_wall
use m_parameters

implicit none

private

! string to id restart data
character(len=*), parameter :: p_lindman_rst_id = "lindman rst data - 0x0004"

interface setup
    module procedure new_lindman
end interface

interface update_e_1d
    module procedure update_lindman_1d_e
end interface

interface update_e_2d
    module procedure update_lindman_2d_e
end interface

interface update_e_3d
    module procedure update_lindman_3d_e
end interface

interface update_b_1d
    module procedure update_lindman_1d_b
end interface

interface update_b_2d
    module procedure update_lindman_2d_b
end interface

interface update_b_3d
    module procedure update_lindman_3d_b
end interface

interface move_window
     module procedure move_window_lindman
end interface

interface update_boundary
     module procedure update_boundary_lindman
end interface

interface cleanup
    module procedure cleanup_lindman
end interface

interface restart_read
    module procedure restart_read_lindman
end interface

interface restart_write
    module procedure restart_write_lindman
end interface


interface is_active
    module procedure is_active_lindman
end interface

interface reshape_obj
  module procedure reshape_lind
end interface

interface update_emfbound
  module procedure update_emfbound_lindman
end interface

public :: t_lindman, setup, cleanup
public :: is_active, restart_write, restart_read
public :: update_boundary, update_emfbound

public :: update_e_1d, update_e_2d, update_e_3d
public :: update_b_1d, update_b_2d, update_b_3d


! load balancing and new moving window routines
public :: reshape_obj, move_window

contains


!---------------------------------------------------
subroutine cleanup_lindman(this)
!---------------------------------------------------
   implicit none

   type (t_lindman), intent(inout) :: this

   this%idir = -1
   call cleanup( this%e_buffer )
   call cleanup( this%b_buffer )

end subroutine cleanup_lindman
!---------------------------------------------------

!---------------------------------------------------
subroutine restart_write_lindman(this, restart_handle)
!---------------------------------------------------

   implicit none

   type (t_lindman), intent(in) :: this
   type( t_restart_handle ), intent(inout) :: restart_handle

   integer :: ierr

   restart_io_wr( p_lindman_rst_id, restart_handle, ierr )
   if ( ierr/=0 ) then
     ERROR('error writing restart data for lindman object.')
     call abort_program(p_err_rstwrt)
   endif

   ! write wall data
   call restart_write(this%e_buffer, restart_handle)
   call restart_write(this%b_buffer, restart_handle)


end subroutine restart_write_lindman
!---------------------------------------------------

!---------------------------------------------------
subroutine restart_read_lindman(this, restart_handle)
!---------------------------------------------------

   implicit none
   type (t_lindman), intent(inout) :: this
   type( t_restart_handle ), intent(in) :: restart_handle

   character(len=len(p_lindman_rst_id)) :: rst_id
   integer :: ierr


   ! clear lindman data
   call cleanup(this)

   ! check if restart file is compatible
   restart_io_rd( rst_id, restart_handle, ierr )
   if ( ierr/=0 ) then
      ERROR('error reading restart data for lindman object.')
      call abort_program(p_err_rstrd)
   endif
   if ( rst_id /= p_lindman_rst_id) then
      ERROR('Corrupted restart file, or restart file ')
      ERROR('from incompatible binary (lindman)')
      call abort_program(p_err_rstrd)
   endif

end subroutine restart_read_lindman
!---------------------------------------------------

!---------------------------------------------------
subroutine new_lindman(this, orientation, i_dim, b, e, dt, restart, restart_handle)
!---------------------------------------------------
!    setup new lindman b.c.
!---------------------------------------------------

  implicit none

  ! dummy variables

  type (t_lindman), intent(inout) :: this
  integer, intent(in) :: orientation, i_dim
  type ( t_vdf ), intent(in) :: b, e
  real (p_double), intent(in) :: dt

  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle


  ! local variables

  integer, dimension(2, p_max_dim) :: lgc_num
  integer, dimension(p_max_dim) :: lnx
  integer, dimension(2) :: range

  integer, dimension(2,p_max_dim) :: i_bounds
  integer :: x_dim

 ! executable statements
 ASSERT(p_f_dim == 3)

 ! Read restart data.
 ! (this is currently empty and just checks the restart id)
 if ( restart ) then
   call restart_read_lindman( this, restart_handle )
 endif

 call cleanup(this)

 this%idir = i_dim
 this%orientation = orientation

 x_dim   = b % x_dim_
 lgc_num(:,1:x_dim) = b%gc_num_(:,1:x_dim)
 lnx(1:x_dim) = b%nx_(1:x_dim)

 ! set simulation parameters
 this%dtdxi(1:x_dim) = real( dt /b%dx_(1:x_dim), p_k_fld )
 this%disp_corr = real( (1.0_p_double - dt/b%dx_(i_dim))/ &
                        (1.0_p_double + dt/b%dx_(i_dim)), p_k_fld )

 ! kludge
 this%dtdxi(1:x_dim) = this%dtdxi(1:x_dim) * 0.5_p_k_fld

 ! set dimensions for buffers
 i_bounds(1,1:x_dim)=1
 i_bounds(2,1:x_dim)=b%nx_(1:x_dim)

 select case(x_dim)

  case (1)
    if (i_dim == 1 ) then
        this%perp(1)=2
        this%perp(2)=3
    else
        ERROR('LINDMAN(NEW): Error in i_dir')
        call abort_program( p_err_invalid )
    endif


  case (2)
    select case(i_dim)
      case (1)
        this%perp(1)=2
        this%perp(2)=3
        this%imin(1)=i_bounds(1,2)
        this%imax(1)=i_bounds(2,2)
      case (2)
        this%perp(1)=3
        this%perp(2)=1
        this%imin(1)=i_bounds(1,1)
        this%imax(1)=i_bounds(2,1)
      case default
        ERROR('LINDMAN(NEW): Error in i_dir')
         call abort_program( p_err_invalid )
   end select

 case (3)

    select case(i_dim)
      case (1)
        this%perp(1)=2
        this%perp(2)=3
        this%imin(1)=i_bounds(1,2)
        this%imin(2)=i_bounds(1,2)
        this%imax(1)=i_bounds(2,2)
        this%imax(2)=i_bounds(2,3)
      case (2)
        this%perp(1)=3
        this%perp(2)=1
        this%imin(1)=i_bounds(1,3)
        this%imin(2)=i_bounds(1,1)
        this%imax(1)=i_bounds(2,3)
        this%imax(2)=i_bounds(2,1)
      case (3)
        this%perp(1)=1
        this%perp(2)=2
        this%imin(1)=i_bounds(1,1)
        this%imin(2)=i_bounds(1,2)
        this%imax(1)=i_bounds(2,1)
        this%imax(2)=i_bounds(2,2)
      case default
        ERROR('LINDMAN(NEW): Error in i_dir')
        call abort_program( p_err_invalid )
    end select

 end select


 ! set range for buffers
 select case (orientation)
   case ( p_lower )
     ! range(1) = 1-lgc_num(p_lower,i_dim)
     range(1) = 0
     range(2) = 1

   case ( p_upper )
     range(1) = lnx(i_dim)
     ! range(2) = lnx(i_dim)  + lgc_num(p_upper,i_dim)
     range(2) = lnx(i_dim)  + 2

   case default
     ERROR('LINDMAN(NEW): Error in orientation')
     call abort_program( p_err_invalid )
 end select
 this%range = range

 ! This will read wall restart data if needed
 call new(this%e_buffer, e, i_dim, orientation, range, restart, restart_handle )
 call new(this%b_buffer, b, i_dim, orientation, range, restart, restart_handle )

end subroutine new_lindman
!---------------------------------------------------


!---------------------------------------------------
subroutine move_window_lindman( this, space )
!---------------------------------------------------
   implicit none

   type (t_lindman), intent(inout) :: this
   type (t_space), intent(in) :: space ! local or global, only nx_move is needed


   if (this%idir > 0) then
      call move_window( this%b_buffer, space )
      call move_window( this%e_buffer, space )
   end if

end subroutine move_window_lindman
!---------------------------------------------------

!---------------------------------------------------
subroutine update_boundary_lindman(this, no_co, nx_move, send_msg, recv_msg)
!---------------------------------------------------
   use m_vdf_comm, only : t_vdf_msg
   implicit none

   type (t_lindman), intent(inout) :: this
   class( t_node_conf ),intent(in) :: no_co
   integer, dimension(:), intent(in):: nx_move
   type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

   ! only update boundaries if boundary is active
   if (this%idir > 0) then
      call update_boundary(this%b_buffer, p_vdf_replace, no_co, nx_move, send_msg, recv_msg)
      call update_boundary(this%e_buffer, p_vdf_replace, no_co, nx_move, send_msg, recv_msg)
   end if

end subroutine update_boundary_lindman
!---------------------------------------------------

!---------------------------------------------------
subroutine update_lindman_1d_e(this, e)
!---------------------------------------------------

   implicit none

  type (t_lindman) , intent(inout) :: this
  type (t_vdf), intent(inout) :: e

  integer:: iperp
  integer, dimension(2) :: i_range

  i_range=range(this%e_buffer)

  select case(this%orientation)

  case(p_lower)
    do iperp=1,p_f_dim - 1
        e%f1(this%perp(iperp),i_range(1)) = &
             this%e_buffer%f1(this%perp(iperp),i_range(1)+1) + &
             this%disp_corr*( this%e_buffer%f1(this%perp(iperp),i_range(1)  )- &
                                          e%f1(this%perp(iperp),i_range(1)+1))

        this%e_buffer%f1(this%perp(iperp), i_range(1)) = &
                    e%f1(this%perp(iperp), i_range(1))

        this%e_buffer%f1(this%perp(iperp), i_range(1)+1) = &
                    e%f1(this%perp(iperp), i_range(1)+1)
    end do
  case(p_upper)
    do iperp=1,p_f_dim - 1
        e%f1(this%perp(iperp),i_range(2))= &
             this%e_buffer%f1(this%perp(iperp),i_range(2)-1) + &
             this%disp_corr*( this%e_buffer%f1(this%perp(iperp),i_range(2)  )-&
                                          e%f1(this%perp(iperp),i_range(2)-1))

        this%e_buffer%f1(this%perp(iperp), i_range(2)-1) = &
                    e%f1(this%perp(iperp), i_range(2)-1)

        this%e_buffer%f1(this%perp(iperp), i_range(2)) = &
                    e%f1(this%perp(iperp), i_range(2))
    end do
  end select

end subroutine update_lindman_1d_e
!---------------------------------------------------

!---------------------------------------------------
subroutine update_lindman_1d_b(this, e)
!---------------------------------------------------

implicit none

  type (t_lindman), intent(inout) :: this
  type ( t_vdf ), intent(in) :: e

  integer :: i_range(2)

  i_range =range(this%b_buffer)

  select case(this%orientation)
  case(p_lower)
    !wall_ptr(1,i_range(1)) = &
    !   wall_ptr(1,i_range(1))

    this%b_buffer%f1(2,i_range(1)) = &
       this%b_buffer%f1(2,i_range(1)) + &
       this%dtdxi(1)*( e%f1(3,i_range(1)+1) - &
                       e%f1(3,i_range(1)))

    this%b_buffer%f1(3,i_range(1)) = &
       this%b_buffer%f1(3,i_range(1)) - &
       this%dtdxi(1)*( e%f1(2,i_range(1)+1) - &
                       e%f1(2,i_range(1)))

  case(p_upper)
   ! do nothing
  end select

end subroutine update_lindman_1d_b
!---------------------------------------------------

!---------------------------------------------------
subroutine update_lindman_2d_e(this, e)
!---------------------------------------------------
!---------------------------------------------------

implicit none

type (t_lindman), intent(inout) :: this
type (t_vdf ), intent(inout) :: e

integer :: iperp
integer :: ind1,i_range(2)

i_range=range(this%e_buffer)

select case (this%idir)
case (1)
   select case (this%orientation)
      case (p_lower)
         do ind1=this%imin(1),this%imax(1)
            do iperp=1,p_f_dim - 1
                e%f2(this%perp(iperp),i_range(1),ind1) = &
                   this%e_buffer%f2(this%perp(iperp),i_range(1)+1,ind1) + &
                   this%disp_corr*(this%e_buffer%f2(this%perp(iperp),i_range(1),ind1) - &
                                          e%f2(this%perp(iperp),i_range(1)+1,ind1))

                this%e_buffer%f2(this%perp(iperp),i_range(1),ind1) = &
                   e%f2(this%perp(iperp), i_range(1),ind1)
                this%e_buffer%f2(this%perp(iperp),i_range(1)+1,ind1) = &
                   e%f2(this%perp(iperp),i_range(1)+1,ind1)
            end do
         end do

      case (p_upper)
         do ind1=this%imin(1),this%imax(1)
           do iperp=1,p_f_dim - 1
             e%f2(this%perp(iperp),i_range(2),ind1) = &
              this%e_buffer%f2(this%perp(iperp),i_range(2)-1,ind1) + &
              this%disp_corr*(this%e_buffer%f2(this%perp(iperp),i_range(2),ind1) - &
                                     e%f2(this%perp(iperp),i_range(2)-1,ind1))

             this%e_buffer%f2(this%perp(iperp), i_range(2)-1,ind1) = &
                    e%f2(this%perp(iperp), i_range(2)-1,ind1)
             this%e_buffer%f2(this%perp(iperp), i_range(2),ind1) = &
                    e%f2(this%perp(iperp), i_range(2),ind1)
           end do
         end do
   end select

case (2)
   select case (this%orientation)
     case (p_lower)
       do ind1=this%imin(1),this%imax(1)
         do iperp=1,p_f_dim - 1
             e%f2(this%perp(iperp),ind1,i_range(1)) = &
               this%e_buffer%f2(this%perp(iperp),ind1,i_range(1)+1) + &
               this%disp_corr*(this%e_buffer%f2(this%perp(iperp),ind1,i_range(1)  ) - &
                                           e%f2(this%perp(iperp),ind1,i_range(1)+1))

             this%e_buffer%f2(this%perp(iperp),ind1,i_range(1)) = &
                         e%f2(this%perp(iperp),ind1,i_range(1))
             this%e_buffer%f2(this%perp(iperp),ind1,i_range(1)+1) = &
                         e%f2(this%perp(iperp),ind1,i_range(1)+1)
         end do
       end do
     case (p_upper)
       do ind1=this%imin(1),this%imax(1)
         do iperp=1,p_f_dim - 1
           e%f2(this%perp(iperp),ind1,i_range(2)) = &
              this%e_buffer%f2(this%perp(iperp),ind1,i_range(2)-1) + &
              this%disp_corr*(this%e_buffer%f2(this%perp(iperp),ind1,i_range(2)  )-&
                                          e%f2(this%perp(iperp),ind1,i_range(2)-1))

           this%e_buffer%f2(this%perp(iperp),ind1,i_range(2)-1)= &
                       e%f2(this%perp(iperp),ind1,i_range(2)-1)
           this%e_buffer%f2(this%perp(iperp),ind1,i_range(2))= &
                       e%f2(this%perp(iperp),ind1,i_range(2))
         end do
       end do
   end select
end select

end subroutine update_lindman_2d_e
!---------------------------------------------------

!---------------------------------------------------
subroutine update_lindman_2d_b(this, e)
!---------------------------------------------------
! Dispite the name, we actually require the values of the electric
! field for this routine
! --------------------------------------------------
implicit none
type (t_lindman), intent(inout) :: this
type (t_vdf), intent(in) :: e

integer :: ind1,i_range(2)

i_range =range(this%b_buffer)


select case (this%idir)

case (1)
select case (this%orientation)

case (p_lower)
  do ind1=this%imin(1),this%imax(1)
      this%b_buffer%f2(1,i_range(1),ind1) = this%b_buffer%f2(1,i_range(1),ind1) - &
          this%dtdxi(2)*(e%f2(3,i_range(1),ind1+1)-e%f2(3,i_range(1),ind1))

      this%b_buffer%f2(2,i_range(1),ind1) = this%b_buffer%f2(2,i_range(1),ind1) + &
          this%dtdxi(1)*(e%f2(3,i_range(1)+1,ind1)-e%f2(3,i_range(1),ind1))

      this%b_buffer%f2(3,i_range(1),ind1) = this%b_buffer%f2(3,i_range(1),ind1) - &
          this%dtdxi(1)*(e%f2(2,i_range(1)+1,ind1)-e%f2(2,i_range(1),ind1)) + &
          this%dtdxi(2)*(e%f2(1,i_range(1),ind1+1)-e%f2(1,i_range(1),ind1))
  end do

case (p_upper)
  ! no update required for upper lindman
  ! do nothing

end select

case (2)
select case (this%orientation)

case (p_lower)
  do ind1=this%imin(1),this%imax(1)
      this%b_buffer%f2(1,ind1,i_range(1)) = this%b_buffer%f2(1,ind1,i_range(1)) - &
          this%dtdxi(2)*(e%f2(3,ind1,i_range(1)+1)-e%f2(3,ind1,i_range(1)))

      this%b_buffer%f2(2,ind1,i_range(1)) = this%b_buffer%f2(2,ind1,i_range(1)) + &
          this%dtdxi(1)*(e%f2(3,ind1+1,i_range(1))-e%f2(3,ind1,i_range(1)))

      this%b_buffer%f2(3,ind1,i_range(1)) = this%b_buffer%f2(3,ind1,i_range(1)) - &
          this%dtdxi(1)*(e%f2(2,ind1+1,i_range(1))-e%f2(2,ind1,i_range(1))) + &
          this%dtdxi(2)*(e%f2(1,ind1,i_range(1)+1)-e%f2(1,ind1,i_range(1)))
  end do

case (p_upper)
  ! no update required for upper lindman
  ! do nothing
end select

end select

end subroutine update_lindman_2d_b
!---------------------------------------------------

!---------------------------------------------------
subroutine update_lindman_3d_e(this, e)
!---------------------------------------------------
!---------------------------------------------------

implicit none
type (t_lindman), intent(inout) ::  this
type (t_vdf), intent(inout) :: e

!! local variables
integer:: iperp, ind1,ind2
integer, dimension(2) :: i_range



if(this%idir < 1) then
ERROR('LINDMAN(UPDATE): The Lindman Object is not active')
call abort_program(p_err_invalid)
end if

i_range = range(this%e_buffer)


select case (this%idir)

case (1)

 select case (this%orientation)

   case (p_lower)
     do ind2=this%imin(2),this%imax(2)
       do ind1=this%imin(1),this%imax(1)
         do iperp=1,p_f_dim - 1
           e%f3(this%perp(iperp),i_range(1),ind1,ind2) = &
              this%e_buffer%f3(this%perp(iperp),i_range(1)+1,ind1,ind2) + &
              this%disp_corr*(this%e_buffer%f3(this%perp(iperp),i_range(1),ind1,ind2) - &
              e%f3(this%perp(iperp),i_range(1)+1,ind1,ind2))

           this%e_buffer%f3(this%perp(iperp),i_range(1),ind1,ind2) = &
              e%f3(this%perp(iperp),i_range(p_lower),ind1,ind2)

           this%e_buffer%f3(this%perp(iperp),i_range(1)+1,ind1,ind2) = &
              e%f3(this%perp(iperp),i_range(p_lower)+1,ind1,ind2)
         end do
       end do
     end do

   case (p_upper)
     do ind2=this%imin(2),this%imax(2)
       do ind1=this%imin(1),this%imax(1)
         do iperp=1,p_f_dim - 1
           e%f3(this%perp(iperp),i_range(2),ind1,ind2) = &
             this%e_buffer%f3(this%perp(iperp),i_range(2)-1,ind1,ind2) + &
             this%disp_corr*(this%e_buffer%f3(this%perp(iperp),i_range(2),ind1,ind2) - &
             e%f3(this%perp(iperp),i_range(2)-1,ind1,ind2))

           this%e_buffer%f3(this%perp(iperp),i_range(p_upper),ind1,ind2) = &
              e%f3(this%perp(iperp),i_range(p_upper),ind1,ind2)
           this%e_buffer%f3(this%perp(iperp),i_range(2)-1,ind1,ind2) = &
              e%f3(this%perp(iperp),i_range(2)-1,ind1,ind2)
         end do
       end do
     end do

 end select

case (2)
 select case (this%orientation)
   case (p_lower)
     do ind1=this%imin(1),this%imax(1)
       do ind2=this%imin(2),this%imax(2)
         do iperp=1,p_f_dim - 1
             e%f3(this%perp(iperp),ind2,i_range(1),ind1) = &
               this%e_buffer%f3(this%perp(iperp),ind2,i_range(1)+1,ind1) + &
               this%disp_corr*(this%e_buffer%f3(this%perp(iperp),ind2,i_range(1),ind1) - &
               e%f3(this%perp(iperp),i_range(1)+1,ind1,ind2))

             this%e_buffer%f3(this%perp(iperp),ind2,i_range(1),ind1) = &
               e%f3(this%perp(iperp),i_range(p_lower),ind1,ind2)

             this%e_buffer%f3(this%perp(iperp),ind2,i_range(1)+1,ind1) = &
               e%f3(this%perp(iperp),ind2,i_range(p_lower)+1,ind1)
         end do
       end do
     end do

   case (p_upper)
     do ind1=this%imin(1),this%imax(1)
       do ind2=this%imin(2),this%imax(2)
         do iperp=1,p_f_dim - 1
           e%f3(this%perp(iperp),ind2,i_range(2),ind1) = &
             this%e_buffer%f3(this%perp(iperp),ind2,i_range(2)-1,ind1) + &
             this%disp_corr*(this%e_buffer%f3(this%perp(iperp),ind2,i_range(2),ind1) - &
             e%f3(this%perp(iperp),ind2,i_range(2)-1,ind1))

           this%e_buffer%f3(this%perp(iperp),ind2,i_range(p_upper),ind1) = &
             e%f3(this%perp(iperp),ind2,i_range(p_upper),ind1)

           this%e_buffer%f3(this%perp(iperp),ind2,i_range(2)-1,ind1) = &
             e%f3(this%perp(iperp),ind2,i_range(2)-1,ind1)
         end do
       end do
     end do
 end select

case (3)
 select case (this%orientation)

   case (p_lower)
     do ind2=this%imin(2),this%imax(2)
       do ind1=this%imin(1),this%imax(1)
         do iperp=1,p_f_dim - 1
           e%f3(this%perp(iperp),ind1,ind2,i_range(1)) = &
             this%e_buffer%f3(this%perp(iperp),ind1,ind2,i_range(1)+1) + &
             this%disp_corr*(this%e_buffer%f3(this%perp(iperp),ind1,ind2,i_range(1)) - &
             e%f3(this%perp(iperp),ind1,ind2,i_range(1)+1))

           this%e_buffer%f3(this%perp(iperp),ind1,ind2,i_range(1)) = &
             e%f3(this%perp(iperp),ind1,ind2,i_range(p_lower))

           this%e_buffer%f3(this%perp(iperp),ind1,ind2,i_range(p_lower)+1) = &
             e%f3(this%perp(iperp),ind1,ind2,i_range(p_lower)+1)
         end do
       end do
     end do

   case (p_upper)
     do ind2=this%imin(2),this%imax(2)
       do ind1=this%imin(1),this%imax(1)
         do iperp=1,p_f_dim - 1
           e%f3(this%perp(iperp),ind1,ind2,i_range(2)) = &
             this%e_buffer%f3(this%perp(iperp),ind1,ind2,i_range(2)-1) + &
             this%disp_corr*(this%e_buffer%f3(this%perp(iperp),ind1,ind2,i_range(2)) - &
             e%f3(this%perp(iperp),ind1,ind2,i_range(2)-1))

           this%e_buffer%f3(this%perp(iperp),ind1,ind2,i_range(p_upper)) = &
             e%f3(this%perp(iperp),ind1,ind2,i_range(p_upper))

           this%e_buffer%f3(this%perp(iperp),ind1,ind2,i_range(p_upper-1)) = &
             e%f3(this%perp(iperp),ind1,ind2,i_range(p_upper)-1)
         end do
       end do
     end do
 end select

end select

end subroutine update_lindman_3d_e
!---------------------------------------------------


!---------------------------------------------------
subroutine update_lindman_3d_b(this, e )
!---------------------------------------------------
!---------------------------------------------------

implicit none

type (t_lindman), intent(inout) :: this
type (t_vdf), intent(in) :: e

! local variables
integer :: ind1,ind2
integer, dimension(2) :: i_range


if( this%idir < 1) then
ERROR('LINDMAN(UPDATE_B): The Lindman Object is not active')
call abort_program(p_err_invalid)
end if

i_range = range(this%b_buffer)

select case (this%idir)
case (1)
  do ind1=this%imin(1),this%imax(1)
    do ind2=this%imin(2),this%imax(2)
      this%b_buffer%f3(1,i_range(1),ind1,ind2) = &
          this%b_buffer%f3(1,i_range(1),ind1,ind2) - &
             this%dtdxi(2)*(e%f3(3,i_range(1),ind1+1,ind2) - &
                            e%f3(3,i_range(1),ind1,ind2)) + &
             this%dtdxi(3)*(e%f3(2,i_range(1),ind1,ind2+1) - &
                            e%f3(2,i_range(1),ind1,ind2))

      this%b_buffer%f3(2,i_range(1),ind1,ind2) = &
          this%b_buffer%f3(2,i_range(1),ind1,ind2) - &
             this%dtdxi(3)*(e%f3(1,i_range(1),ind1,ind2+1) - &
                            e%f3(1,i_range(1),ind1,ind2)) + &
             this%dtdxi(1)*(e%f3(3,i_range(1)+1,ind1,ind2) - &
                            e%f3(3,i_range(1),ind1,ind2))

      this%b_buffer%f3(3,i_range(1),ind1,ind2) = &
          this%b_buffer%f3(3,i_range(1),ind1,ind2) - &
             this%dtdxi(1)*(e%f3(2,i_range(1)+1,ind1,ind2) - &
                            e%f3(2,i_range(1),ind1,ind2)) + &
             this%dtdxi(2)*(e%f3(1,i_range(1),ind1+1,ind2) - &
                            e%f3(1,i_range(1),ind1,ind2))
    end do
  end do

case (2)
  do ind1=this%imin(1),this%imax(1)
      do ind2=this%imin(2),this%imax(2)
        this%b_buffer%f3(1,ind2,i_range(1),ind1) = &
            this%b_buffer%f3(1,ind2,i_range(1),ind1) - &
               this%dtdxi(2)*(e%f3(3,ind2,i_range(1)+1,ind1) - &
                              e%f3(3,ind2,i_range(1),ind1)) + &
               this%dtdxi(3)*(e%f3(2,ind2,i_range(1),ind1+1) - &
                              e%f3(3,ind2,i_range(1),ind1))

        this%b_buffer%f3(2,ind2,i_range(1),ind1) = &
            this%b_buffer%f3(2,ind2,i_range(1),ind1) - &
               this%dtdxi(3)*(e%f3(1,ind2,i_range(1),ind1+1) - &
                              e%f3(1,ind2,i_range(1),ind1)) + &
               this%dtdxi(1)*(e%f3(3,ind2+1,i_range(1),ind1) - &
                              e%f3(3,ind2,i_range(1),ind1))

        this%b_buffer%f3(3,ind2,i_range(1),ind1) = &
            this%b_buffer%f3(3,ind2,i_range(1),ind1) - &
                this%dtdxi(1)*(e%f3(2,ind2+1,i_range(1),ind1) - &
                               e%f3(2,ind2,i_range(1),ind1)) + &
                this%dtdxi(2)*(e%f3(1,ind2,i_range(1)+1,ind1) - &
                               e%f3(1,ind2,i_range(1),ind1))
      end do
  end do

case (3)
  do ind1=this%imin(1),this%imax(1)
      do ind2=this%imin(2),this%imax(2)
        this%b_buffer%f3(1,ind1,ind2,i_range(1)) = &
           this%b_buffer%f3(1,ind1,ind2,i_range(1)) - &
              this%dtdxi(2)*(e%f3(3,ind1,ind2+1,i_range(1)) - &
                             e%f3(3,ind1,ind2,i_range(1))) + &
              this%dtdxi(3)*(e%f3(2,ind1,ind2,i_range(1)+1) - &
                             e%f3(2,ind1,ind2,i_range(1)))

           this%b_buffer%f3(2,ind1,ind2,i_range(1)) = &
              this%b_buffer%f3(2,ind1,ind2,i_range(1)) - &
                 this%dtdxi(3)*(e%f3(1,ind1,ind2,i_range(1)+1) - &
                                e%f3(1,ind1,ind2,i_range(1))) + &
                 this%dtdxi(1)*(e%f3(3,ind1+1,ind2,i_range(1)) - &
                                e%f3(3,ind1,ind2,i_range(1)))

           this%b_buffer%f3(3,ind1,ind2,i_range(1)) = &
              this%b_buffer%f3(3,ind1,ind2,i_range(1)) - &
                 this%dtdxi(1)*(e%f3(2,ind1+1,ind2,i_range(1)) - &
                                e%f3(2,ind1,ind2,i_range(1))) + &
                 this%dtdxi(2)*(e%f3(1,ind1,ind2+1,i_range(1)) - &
                                e%f3(1,ind1,ind2,i_range(1)))
        end do
  end do

end select

end subroutine update_lindman_3d_b
!---------------------------------------------------

!---------------------------------------------------
subroutine update_emfbound_lindman( this, b, e, dim )
!---------------------------------------------------
!
! the actual boundary values for the fields are
! calculated in update_boundary_e and
! update_boundary_b
! here the boundary values calculated for the
! fields in the previous loop in update_boundary_e/b
! (called during the field solve) only need to be
! restored (in case the values in the fields were
! changed since then) in order to be available for
! the particle advance as well as the next field
! solve
!
!---------------------------------------------------

implicit none

type( t_lindman ), intent(in) :: this
type( t_vdf ), intent(inout)  :: b, e
integer, intent(in) :: dim

! local variables
integer :: i1, i2, i3
integer :: lb1, lb2, lb3
integer :: ub1, ub2, ub3


if ( dim /= this%idir ) then
  ERROR("invalid dimension specified")
  call abort_program( p_err_invalid )
endif

select case (p_x_dim)

case(1)
select case ( this%orientation )
 case (p_lower)
      b%f1(:,0) = this%b_buffer%f1(:,0)
      e%f1(:,0) = this%e_buffer%f1(:,0)

 case (p_upper)

      i1 = b % nx_(dim)
      e%f1(:,i1+1) = this%e_buffer%f1(:,i1+1)

end select

case(2)

select case ( this%orientation )
  case (p_lower)

      select case ( dim )
       case (1)
         lb2 = lbound(e%f2,3)
         ub2 = ubound(e%f2,3)
         do i2= lb2, ub2
            b%f2(:,0,i2) = this%b_buffer%f2(:,0,i2)
            e%f2(:,0,i2) = this%e_buffer%f2(:,0,i2)
         enddo
       case (2)
         lb1 = lbound(e%f2,2)
         ub1 = ubound(e%f2,2)
         do i1=lb1, ub1
            e%f2(:,i1,0) = this%e_buffer%f2(:,i1,0)
            b%f2(:,i1,0) = this%b_buffer%f2(:,i1,0)
         enddo
      end select

  case (p_upper)

      select case ( dim )
       case (1)
         i1 = b % nx_(dim)
         lb2 = lbound(e%f2,3)
         ub2 = ubound(e%f2,3)
         do i2=lb2, ub2
            e%f2(:,i1+1,i2) = this%e_buffer%f2(:,i1+1,i2)
         enddo

       case (2)
         i2 = b % nx_(dim)
         lb1 = lbound(e%f2,2)
         ub1 = ubound(e%f2,2)
         do i1=lb1, ub1
            e%f2(:,i1,i2+1) = this%e_buffer%f2(:,i1,i2+1)
         enddo
      end select

end select



case(3)

select case( this%orientation )
 case (p_lower)
      select case( dim )
       case (1)
         lb2 = lbound(e%f3,3)
         ub2 = ubound(e%f3,3)
         lb3 = lbound(e%f3,4)
         ub3 = ubound(e%f3,4)

         do i3 = lb3, ub3
           do i2 = lb2, ub2
              b%f3(:,0,i2,i3) = this%b_buffer%f3(:,0,i2,i3)
              e%f3(:,0,i2,i3) = this%e_buffer%f3(:,0,i2,i3)
           enddo
         enddo
       case (2)
         lb1 = lbound(e%f3,2)
         ub1 = ubound(e%f3,2)
         lb3 = lbound(e%f3,4)
         ub3 = ubound(e%f3,4)

         do i3 = lb3, ub3
           do i1 = lb1, ub1
              b%f3(:,i1,0,i3) = this%b_buffer%f3(:,i1,0,i3)
              e%f3(:,i1,0,i3) = this%e_buffer%f3(:,i1,0,i3)
           enddo
         enddo
       case (3)
         lb1 = lbound(e%f3,2)
         ub1 = ubound(e%f3,2)
         lb2 = lbound(e%f3,3)
         ub2 = ubound(e%f3,3)
         do i2 = lb2,ub2
           do i1 = lb1,ub1
              b%f3(:,i1,i2,0) = this%b_buffer%f3(:,i1,i2,0)
              e%f3(:,i1,i2,0) = this%e_buffer%f3(:,i1,i2,0)
           enddo
         enddo
      end select

 case (p_upper)
      select case( dim )
       case (1)
         i1  = b % nx_(dim)+1
         lb2 = lbound(e%f3,3)
         ub2 = ubound(e%f3,3)
         lb3 = lbound(e%f3,4)
         ub3 = ubound(e%f3,4)

         do i3 = lb3, ub3
           do i2 = lb2, ub2
              e%f3(:,i1,i2,i3) = this%e_buffer%f3(:,i1,i2,i3)
           enddo
         enddo
       case (2)
         lb1 = lbound(e%f3,2)
         ub1 = ubound(e%f3,2)
         i2  = b % nx_(dim)+1
         lb3 = lbound(e%f3,4)
         ub3 = ubound(e%f3,4)

         do i3 = lb3, ub3
           do i1 = lb1, ub1
             e%f3(:,i1,i2,i3) = this%e_buffer%f3(:,i1,i2,i3)
           enddo
         enddo
       case (3)
         lb1 = lbound(e%f3,2)
         ub1 = ubound(e%f3,2)
         lb2 = lbound(e%f3,3)
         ub2 = ubound(e%f3,3)
         i3  = b % nx_(dim)+1

         do i2=lb2, ub2
           do i1=lb1, ub1
              e%f3(:,i1,i2,i3) = this%e_buffer%f3(:,i1,i2,i3)
           enddo
         enddo
      end select
end select


end select


end subroutine update_emfbound_lindman
!---------------------------------------------------


!---------------------------------------------------
function is_active_lindman(this)
!---------------------------------------------------
  implicit none

  logical :: is_active_lindman
  type (t_lindman), intent(in) :: this

  is_active_lindman = ( this%idir > 0 )

end function is_active_lindman
!---------------------------------------------------


!---------------------------------------------------
subroutine reshape_lind( this, old_lb, new_lb, no_co, send_msg, recv_msg )
!---------------------------------------------------
    use m_vdf_comm

  implicit none

  type(t_lindman), intent(inout) :: this
  class(t_grid), intent(in) :: old_lb, new_lb
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer, parameter, dimension(3) :: od1 = (/2,1,1/)
  integer, parameter, dimension(3) :: od2 = (/3,3,2/)
  integer, dimension(p_max_dim) :: nx
  ! exec

  ! no lindman object data needs changing for 1d
  if ( new_lb%x_dim > 1) then

    ! get new values for imin and imax
    nx(1:new_lb%x_dim) = new_lb%my_nx(3,1:new_lb%x_dim)

    this%imin = 1
    this%imax(1)=nx(od1(this%idir))

    if ( new_lb%x_dim == 3 ) then
      this%imax(2)=nx(od2(this%idir))
    endif
  endif

  ! reshape wall object data
  call reshape_copy( this%b_buffer, old_lb, new_lb, no_co, send_msg, recv_msg )
  call reshape_copy( this%e_buffer, old_lb, new_lb, no_co, send_msg, recv_msg )


end subroutine reshape_lind
!---------------------------------------------------

      end module m_lindman
