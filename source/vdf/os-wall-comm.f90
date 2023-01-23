!#define __DEBUG__ 1

#include "os-config.h"
#include "os-preprocess.fpp"

module m_wall_comm

#include "memory/memory.h"

  use m_parameters

  use m_wall_define
  use m_vdf_comm

  use m_space
  use m_node_conf
  use m_grid_define

  private

  interface move_window
    module procedure move_window_wall
  end interface

  interface update_boundary
    module procedure update_boundary_wall
  end interface

  interface reshape_copy
    module procedure reshape_wall_copy
  end interface

  interface reshape_nocopy
    module procedure reshape_wall_nocopy
  end interface

  interface pack_wall_msg
    module procedure pack_wall_msg
  end interface

  interface unpack_wall_msg
    module procedure unpack_wall_msg
  end interface

  ! module variables for normal vdf comm
  ! type( t_vdf_msg ), dimension(2), save, private :: send_msg, recv_msg

  public :: move_window, update_boundary, reshape_copy, reshape_nocopy
  public :: wait_irecv_wall, wait_isend_wall, update_periodic_wall
  public :: pack_wall_msg, unpack_wall_msg

contains

!---------------------------------------------------------------------------------------------------
! Shifts the data in the wall object according to the values supplied in nx_move.
!---------------------------------------------------------------------------------------------------
subroutine move_window_wall( this, space )

  implicit none

  ! dummy variables

  type( t_wall ), intent( inout ) :: this
  type( t_space ), intent( in ) :: space

  ! local variables
  integer, dimension(this%x_dim) :: lnx_move
  integer :: i_dim, i1, i2, i3, lb, ub

  ! executable statements


  ! check if move is legal
  ! this is for debug purposes and should be removed from production code
  if ( x_dim(space) /= this%x_dim) then
    ERROR("The dimensions of the space object, ", x_dim(space) ," do not")
    ERROR("match the dimensions of the wall object, ", this%x_dim)
    call abort_program( p_err_invalid )
  endif

  lnx_move( 1:this%x_dim ) = nx_move( space )

  ! no move perpendicular to wall
  lnx_move( this%idir ) = 0

  do i_dim = 1, this%x_dim
    if(lnx_move(i_dim) < 0) then
      ERROR("Illegal move, only positive moves allowed")
      call abort_program( p_err_invalid )
    endif
  enddo

  do i_dim = 1, this%x_dim
    if ((lnx_move(i_dim) > this%gc_num(1,i_dim)) .or. &
        (lnx_move(i_dim) > this%gc_num(2,i_dim))) then
       ERROR("Illegal move, too many cells along direction ",i_dim)
       ERROR("gc_num(:,",i_dim,") = ",this%gc_num(:,i_dim))
       ERROR("lnx_move(",i_dim,")  = ",lnx_move(i_dim))
       call abort_program(p_err_invalid)
    endif
  enddo

  ! shift the data according to the supplied parameters

  select case(this%x_dim)

    case(1) ! 1D vdf

       if (lnx_move(1) > 0) then ! shift left
         lb = lbound(this%f1,2)
         ub = ubound(this%f1,2) - lnx_move(1)
         do i1 = lb, ub
           this%f1(:,i1) = this%f1(:,i1+lnx_move(1))
         enddo

         ! fill empty cells with 0.0
         this%f1(:,ub+1:) = 0.0_p_k_fld
       endif
    case(2) ! 2D vdf

       if (lnx_move(1) > 0) then ! shift left

         lb = lbound(this%f2,2)
         ub = ubound(this%f2,2) - lnx_move(1)
         do i1 =  lb, ub
           this%f2(:,i1,:) = this%f2(:,i1+lnx_move(1),:)
         enddo
         this%f2(:,ub+1:,:) = 0.0_p_k_fld

       endif

       if (lnx_move(2) > 0) then ! shift left
         lb = lbound(this%f2,3)
         ub = ubound(this%f2,3) - lnx_move(2)
         do i2 = lb, ub
           this%f2(:,:,i2) = this%f2(:,:,i2+lnx_move(2))
         enddo

         this%f2(:,:,ub+1: ) = 0.0_p_k_fld
       endif

    case(3) ! 3D vdf

       if (lnx_move(1) > 0) then ! shift left
         lb = lbound(this%f3,2)
         ub = ubound(this%f3,2) - lnx_move(1)
         do i1 = lb, ub
           this%f3(:,i1,:, :) = this%f3(:,i1+lnx_move(1),:, :)
         enddo
         this%f3(:,ub+1:,:,:) = 0.0_p_k_fld
       endif

       if (lnx_move(2) > 0) then ! shift left
         lb = lbound(this%f3,3)
         ub = ubound(this%f3,3) - lnx_move(2)
         do i2 = lb, ub
           this%f3(:,:,i2,:) = this%f3(:,:,i2+lnx_move(2),:)
         enddo

         this%f3(:,:,ub+1:,:) = 0.0_p_k_fld
       endif

       if (lnx_move(3) > 0) then ! shift left
         lb = lbound(this%f3,4)
         ub = ubound(this%f3,4) - lnx_move(3)
         do i3 = lb, ub
           this%f3(:,:,:,i3) = this%f3(:,:,:,i3+lnx_move(3))
         enddo

         this%f3(:,:,:,ub+1:) = 0.0_p_k_fld
       endif

  end select


end subroutine move_window_wall
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Update single node periodic gc / grid values for a vdf along a given direction
!---------------------------------------------------------------------------------------------------
subroutine update_periodic_wall( wall, dim, update_type )

  implicit none

  type( t_wall ), intent(inout) :: wall
  integer, intent(in) :: dim, update_type

  integer, dimension(2, p_max_dim) :: lrange, urange
  integer, dimension( p_max_dim ) :: shift

  integer :: i, i1, i2, i3
  real( p_k_fld ) :: val

  ! Number of cells to shift
  do i = 1, wall%x_dim
    if ( i == dim ) then
      shift(i) = wall%nx(i)
    else
      shift(i) = 0
    endif
  enddo

  do i = 1, wall%x_dim
    if ( i == dim ) then
      ! Select only guard cells along update direction
      lrange( p_lower, i ) = 1          - wall%gc_num( p_lower, i )
      lrange( p_upper, i ) = 0
      urange( p_lower, i ) = wall%nx(i) + 1
      urange( p_upper, i ) = wall%nx(i) + wall%gc_num( p_upper, i )
    else
      ! full cell range along other directions
      if ( i == wall%idir ) then
        lrange( p_lower, i ) = wall%range( p_lower )
        lrange( p_upper, i ) = wall%range( p_upper )
      else
        lrange( p_lower, i ) = 1          - wall%gc_num( p_lower, i )
        lrange( p_upper, i ) = wall%nx(i) + wall%gc_num( p_upper, i )
      endif
      urange( p_lower, i ) = lrange( p_lower, i )
      urange( p_upper, i ) = lrange( p_upper, i )
    endif
  enddo

  select case ( update_type )

    case ( p_vdf_add )

      ! Add grid values to corresponding gc values and store in both grid and gc
      select case ( wall%x_dim )
         ! case (1)
         ! No 1D code required

         case (2)
           ! update lower gc
           do i2 = lrange( p_lower, 2 ), lrange( p_upper, 2 )
             do i1 = lrange( p_lower, 1 ), lrange( p_upper, 1 )
               do i = 1, wall%f_dim
                 val =  wall%f2( i, i1, i2 ) + wall%f2( i, i1 + shift(1), i2 + shift(2) )
                 wall%f2( i, i1, i2 )                       = val
                 wall%f2( i, i1 + shift(1), i2 + shift(2) ) = val
               enddo
             enddo
           enddo

           ! update upper gc
           do i2 = urange( p_lower, 2 ), urange( p_upper, 2 )
             do i1 = urange( p_lower, 1 ), urange( p_upper, 1 )
               do i = 1, wall%f_dim
                 val = wall%f2( i, i1, i2 ) + wall%f2( i, i1 - shift(1), i2 - shift(2) )
                 wall%f2( i, i1 - shift(1), i2 - shift(2) ) = val
                 wall%f2( i, i1, i2 )                       = val
               enddo
             enddo
           enddo

         case (3)
           ! update lower gc
           do i3 = lrange( p_lower, 3 ), lrange( p_upper, 3 )
             do i2 = lrange( p_lower, 2 ), lrange( p_upper, 2 )
               do i1 = lrange( p_lower, 1 ), lrange( p_upper, 1 )
                 do i = 1, wall%f_dim
                   val =  wall%f3( i,i1,i2,i3 ) + wall%f3( i, i1+shift(1), i2+shift(2), i3+shift(3) )
                   wall%f3( i, i1, i2, i3 )                                     = val
                   wall%f3( i, i1 + shift(1), i2 + shift(2), i3 + shift(3) ) = val
                 enddo
               enddo
             enddo
           enddo

           ! update upper gc
           do i3 = urange( p_lower, 3 ), urange( p_upper, 3 )
             do i2 = urange( p_lower, 2 ), urange( p_upper, 2 )
               do i1 = urange( p_lower, 1 ), urange( p_upper, 1 )
                 do i = 1, wall%f_dim
                   val =  wall%f3( i,i1,i2,i3 ) + wall%f3( i, i1-shift(1),i-shift(2),i3-shift(3) )
                   wall%f3( i, i1 - shift(1), i2 - shift(2), i3 - shift(3) ) = val
                   wall%f3( i, i1, i2, i3 )                                  = val
                 enddo
               enddo
             enddo
           enddo

      end select

    case ( p_vdf_replace )

      ! Copy gc values from opposing grid values
      select case ( wall%x_dim )

         ! case (1)
         ! No 1D code required

         case (2)
           ! update lower gc
           do i2 = lrange( p_lower, 2 ), lrange( p_upper, 2 )
             do i1 = lrange( p_lower, 1 ), lrange( p_upper, 1 )
               do i = 1, wall%f_dim
                 wall%f2( i, i1, i2 ) = wall%f2( i, i1 + shift(1), i2 + shift(2) )
               enddo
             enddo
           enddo

           ! update upper gc
           do i2 = urange( p_lower, 2 ), urange( p_upper, 2 )
             do i1 = urange( p_lower, 1 ), urange( p_upper, 1 )
               do i = 1, wall%f_dim
                 wall%f2( i, i1, i2 ) = wall%f2( i, i1 - shift(1), i2 - shift(2) )
               enddo
             enddo
           enddo

         case (3)
           ! update lower gc
           do i3 = lrange( p_lower, 3 ), lrange( p_upper, 3 )
             do i2 = lrange( p_lower, 2 ), lrange( p_upper, 2 )
               do i1 = lrange( p_lower, 1 ), lrange( p_upper, 1 )
                 do i = 1, wall%f_dim
                   wall%f3( i, i1, i2, i3 ) = wall%f3( i, i1 + shift(1), &
                                                        i2 + shift(2), &
                                                        i3 + shift(3) )
                 enddo
               enddo
             enddo
           enddo

           ! update upper gc
           do i3 = urange( p_lower, 3 ), urange( p_upper, 3 )
             do i2 = urange( p_lower, 2 ), urange( p_upper, 2 )
               do i1 = urange( p_lower, 1 ), urange( p_upper, 1 )
                 do i = 1, wall%f_dim
                   wall%f3( i, i1, i2, i3 ) = wall%f3( i, i1 - shift(1), &
                                                          i2 - shift(2), &
                                                          i3 - shift(3) )
                 enddo
               enddo
             enddo
           enddo

      end select

  end select

end subroutine update_periodic_wall
!---------------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Updates the VPML wall objects boundary
!! @param   this            Wall object
!! @param   update_type     Type of boundary update, can be p_vdf_add or 
!!                          p_vdf_replace
!! @param   no_co           Parallel universe configuration
!! @param   move_num        Number of cells to move in this iteration (used
!!                          when moving window is active)
!! @param   send_msg        Send message buffers to use
!! @param   recv_msg        Receive message buffers to use
subroutine update_boundary_wall( wall, update_type, no_co, move_num, send_msg, recv_msg )

  implicit none

  type( t_wall ), intent(inout) :: wall
  integer,                 intent(in) :: update_type
  class( t_node_conf ),             intent(in) :: no_co
    type(t_vdf_msg), dimension(:), intent(inout) :: send_msg, recv_msg
  integer, dimension(:), intent(in) :: move_num

  integer :: i, my_node, lneighbor, uneighbor

  my_node = no_co % my_aid()

  do i = 1, wall%x_dim

    ! skip the direction perpendicular to the wall
    if (i /= wall%idir) then

      if ( my_node == neighbor(no_co,p_lower,i) ) then

         !single node periodic
         call update_periodic_wall( wall, i, update_type )

      else

         ! communication with other nodes
         lneighbor = neighbor( no_co, p_lower, i )
         uneighbor = neighbor( no_co, p_upper, i )

         ! post receives
         if ( lneighbor > 0 ) call irecv_wall( wall, update_type, i, p_lower, no_co, move_num(i), recv_msg )
         if ( uneighbor > 0 ) call irecv_wall( wall, update_type, i, p_upper, no_co, move_num(i), recv_msg )

         ! post sends
         ! Note: sends MUST be posted in the opposite order of receives to account for a 2 node
         ! periodic partition, where 1 node sends 2 messages to the same node
         if ( uneighbor > 0 ) call isend_wall( wall, update_type, i, p_upper, no_co, move_num(i), send_msg )
         if ( lneighbor > 0 ) call isend_wall( wall, update_type, i, p_lower, no_co, move_num(i), send_msg )

         ! wait for receives and unpack data
         if ( lneighbor > 0 ) call wait_irecv_wall( wall, p_lower, recv_msg )
         if ( uneighbor > 0 ) call wait_irecv_wall( wall, p_upper, recv_msg )

         ! wait for sends to complete
         if ( uneighbor > 0 ) call wait_isend_wall( p_upper, send_msg )
         if ( lneighbor > 0 ) call wait_isend_wall( p_lower, send_msg )

      endif
    endif
 enddo

end subroutine update_boundary_wall
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! post send message for boundary communication
!---------------------------------------------------------------------------------------------------
subroutine isend_wall( wall, update_type, dim, bnd, no_co, nx_move, send_msg )

  implicit none

  type( t_wall ), intent(in) :: wall
  integer, intent(in) :: update_type, dim, bnd
  class( t_node_conf ), intent(in) :: no_co

  integer, intent(in) :: nx_move
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg

  integer, dimension(2,p_max_dim) :: ns
  integer :: i, bsize, msg_size

  ! sanity check
  if ( send_msg(bnd)%handle /= MPI_REQUEST_NULL ) then
    ERROR( 'isend_msg_wall called when message still waiting to complete' )
    call abort_program( p_err_invalid )
  endif

    if ( neighbor( no_co, bnd, dim ) < 0 ) then
        write(0,*) 'isend_msg_wall called without a valid neighbor, dim, bnd = ',dim, bnd, &
        trim(__FILE__), ':', __LINE__
        call abort_program( p_err_invalid )
    endif

  ! Get range of cells to send
  msg_size = wall%f_dim
  do i = 1, wall%x_dim
    if ( i == dim ) then
       if ( update_type == p_vdf_add ) then
         select case ( bnd )
           case ( p_upper )
                    ns( p_lower, i ) = wall%nx(i) + 1 - wall%gc_num_comm( p_lower, i )
                    ns( p_upper, i ) = wall%nx(i) + wall%gc_num_comm( p_upper, i )
           case ( p_lower )
                    ns( p_lower, i ) = 1 - wall%gc_num_comm( p_lower, i )
                    ns( p_upper, i ) = wall%gc_num_comm( p_upper, i )
         end select
       else
         select case ( bnd )
           case ( p_upper )
                    ns( p_lower, i ) = wall%nx(i) + 1 - wall%gc_num_comm( p_lower, i )
             ns( p_upper, i ) = wall%nx(i) - nx_move
           case ( p_lower )
             ns( p_lower, i ) = 1 - nx_move
                    ns( p_upper, i ) = wall%gc_num_comm( p_upper, i )
         end select
       endif
    else
      ! full cell range along other directions
      ns( p_lower, i ) = lbound( wall, i+1 )
      ns( p_upper, i ) = ubound( wall, i+1 )
    endif
    msg_size = msg_size * (ns( p_upper, i ) - ns( p_lower,i ) + 1 )
  enddo

  ! Use vdf-comm send buffers
  send_msg(bnd)%type                       = update_type
  send_msg(bnd)%range( :, 1 : wall%x_dim ) = ns ( :, 1 : wall%x_dim )
  send_msg(bnd)%msg_size                   = msg_size
  send_msg(bnd)%node                       = neighbor( no_co, bnd, dim )
  send_msg(bnd)%tag                        = bnd

  ! grow message buffer if necessary
  if ( associated( send_msg(bnd)%buffer ) ) then
    bsize = size( send_msg(bnd)%buffer )
  else
    bsize = 0
  endif

  if ( msg_size > bsize ) then
        if ( bsize > 0 ) call freemem( send_msg(bnd)%buffer )
    call alloc( send_msg(bnd)%buffer, (/msg_size/) )
  endif

  ! pack message data
  call pack_wall_msg( send_msg(bnd), wall )

  ! post send using vdf-comm isend
  call isend( send_msg(bnd), no_co )

end subroutine isend_wall
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Wait for a wall BC send message to complete
!---------------------------------------------------------------------------------------------------
subroutine wait_isend_wall( bnd, send_msg )

  implicit none

  integer, intent(in) :: bnd
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg

  call wait( send_msg( bnd ) )

end subroutine wait_isend_wall
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! post recv message for boundary communication
!---------------------------------------------------------------------------------------------------
subroutine irecv_wall( wall, update_type, dim, bnd, no_co, nx_move, recv_msg )

  implicit none

  type(t_wall), intent(inout) :: wall
  integer, intent(in) :: update_type, dim, bnd
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: nx_move
  type(t_vdf_msg), dimension(2), intent(inout) :: recv_msg

  integer, dimension(2,p_max_dim) :: nr
  integer :: i, bsize, msg_size

  ! Check if message is still waiting to complete
  if ( recv_msg( bnd )%handle /= MPI_REQUEST_NULL ) then
    ERROR( 'irecv_wall called when message still waiting to complete' )
    call abort_program( p_err_invalid )
  endif

    ! Sanity check
    if ( neighbor( no_co, bnd, dim ) < 0 ) then
        write(0,*) 'irecv_wall called without a valid neighbor, dim, bnd = ',dim, bnd, &
        trim(__FILE__), ':', __LINE__
        call abort_program( p_err_invalid )
    endif


  ! Get range of cells to receive
    ! Note that there are no messages in the wall normal direction
  msg_size = wall%f_dim
  do i = 1, wall%x_dim
    if ( i == dim ) then
       if ( update_type == p_vdf_add ) then
         select case ( bnd )
           case ( p_upper )
                    nr( p_lower, i ) = wall%nx(i) + 1 - wall%gc_num_comm( p_lower, i )
                    nr( p_upper, i ) = wall%nx(i) + wall%gc_num_comm( p_upper, i )
           case ( p_lower )
                    nr( p_lower, i ) = 1 - wall%gc_num_comm( p_lower, i )
                    nr( p_upper, i ) = wall%gc_num_comm( p_upper, i )
         end select
       else
         select case ( bnd )
           case ( p_upper )
             nr( p_lower, i ) = wall%nx(i) + 1 - nx_move
                    nr( p_upper, i ) = wall%nx(i) + wall%gc_num_comm( p_upper, i )
           case ( p_lower )
                    nr( p_lower, i ) = 1 - wall%gc_num_comm( p_lower, i )
             nr( p_upper, i ) = 0 - nx_move
         end select
       endif
    else
      ! full cell range along other directions
      nr( p_lower, i ) = lbound( wall, i+1 )
      nr( p_upper, i ) = ubound( wall, i+1 )
    endif
    msg_size = msg_size * (nr( p_upper, i ) - nr( p_lower,i ) + 1 )
  enddo

  ! (*debug*)
  if ( msg_size == 0 ) then
    ERROR('Invalid message size, aborting.')
    call abort_program( p_err_invalid )
  endif

  ! Use vdf-comm recv buffers
  recv_msg(bnd)%type                       = update_type
  recv_msg(bnd)%range( :, 1 : wall%x_dim ) = nr ( :, 1 : wall%x_dim )
  recv_msg(bnd)%msg_size                   = msg_size
  recv_msg(bnd)%node                       = neighbor( no_co, bnd, dim )

  ! receive tag is the opposing boundary
  select case (bnd)
    case ( p_upper )
      recv_msg(bnd)%tag = p_lower
    case ( p_lower )
      recv_msg(bnd)%tag = p_upper
  end select

  ! grow message buffer if necessary
  if ( associated( recv_msg(bnd)%buffer ) ) then
    bsize = size( recv_msg(bnd)%buffer )
  else
    bsize = 0
  endif

  if ( msg_size > bsize ) then
        if ( bsize > 0 ) call freemem( recv_msg(bnd)%buffer )
    call alloc( recv_msg(bnd)%buffer, (/msg_size/) )
  endif

  ! post receive
  call irecv( recv_msg(bnd), no_co )

end subroutine irecv_wall
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Wait for a BC vdf recv message to complete
!---------------------------------------------------------------------------------------------------
subroutine wait_irecv_wall( wall, bnd, recv_msg )

  implicit none

  type(t_wall), intent(inout) :: wall
  integer, intent(in) :: bnd
  type(t_vdf_msg), dimension(2), intent(inout) :: recv_msg

  ! wait for message to complete
  call wait( recv_msg( bnd ) )

  ! unpack message data
  call unpack_wall_msg( recv_msg( bnd ), wall )

end subroutine wait_irecv_wall
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Pack wall data into vdf message buffer
!---------------------------------------------------------------------------------------------------
subroutine pack_wall_msg( msg, wall )

  implicit none

  type( t_vdf_msg ), intent(inout) :: msg
  type( t_wall ), intent(in) :: wall

  integer :: k
  integer :: i, i1, i2, i3

  k = 1

  select case ( wall%x_dim )

    ! case (1)
    !   There is no need to communicate in 1D because there is only 1 node in the wall plane

    case (2)
       do i2 = msg%range( p_lower, 2 ), msg%range( p_upper, 2 )
         do i1 = msg%range( p_lower, 1 ), msg%range( p_upper, 1 )
           do i = 1, wall%f_dim
             msg%buffer( k ) = wall%f2( i, i1, i2 )
             k = k + 1
           enddo
         enddo
       enddo

    case (3)
       do i3 = msg%range( p_lower, 3 ), msg%range( p_upper, 3 )
         do i2 = msg%range( p_lower, 2 ), msg%range( p_upper, 2 )
           do i1 = msg%range( p_lower, 1 ), msg%range( p_upper, 1 )
             do i = 1, wall%f_dim
               msg%buffer( k ) = wall%f3( i, i1, i2, i3 )
               k = k + 1
             enddo
           enddo
         enddo
       enddo
  end select

end subroutine pack_wall_msg
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Unapack message buffer data into vdf
!---------------------------------------------------------------------------------------------------
subroutine unpack_wall_msg( msg, wall )

  implicit none

  type( t_vdf_msg ), intent(inout) :: msg
  type( t_wall ), intent(inout) :: wall

  integer :: k
  integer :: i, i1, i2, i3

  k = 1

  if ( msg%type == p_vdf_add ) then

    ! Add message values to local values

    select case ( wall%x_dim )
    ! case (1)
    !   There is no need to communicate in 1D because there is only 1 node in the wall plane

      case (2)
         do i2 = msg%range( p_lower, 2 ), msg%range( p_upper, 2 )
           do i1 = msg%range( p_lower, 1 ), msg%range( p_upper, 1 )
             do i = 1, wall%f_dim
               wall%f2( i, i1, i2 ) = wall%f2( i, i1, i2 ) + msg%buffer( k )
               k = k + 1
             enddo
           enddo
         enddo

      case (3)
         do i3 = msg%range( p_lower, 3 ), msg%range( p_upper, 3 )
           do i2 = msg%range( p_lower, 2 ), msg%range( p_upper, 2 )
             do i1 = msg%range( p_lower, 1 ), msg%range( p_upper, 1 )
               do i = 1, wall%f_dim
                 wall%f2( i, i1, i2 ) = wall%f2( i, i1, i2 ) + msg%buffer( k )
                 k = k + 1
               enddo
             enddo
           enddo
         enddo
    end select

  else

    ! replace local values with message values

    select case ( wall%x_dim )
    ! case (1)
    !   There is no need to communicate in 1D because there is only 1 node in the wall plane

      case (2)
         do i2 = msg%range( p_lower, 2 ), msg%range( p_upper, 2 )
           do i1 = msg%range( p_lower, 1 ), msg%range( p_upper, 1 )
             do i = 1, wall%f_dim
               wall%f2( i, i1, i2 ) = msg%buffer( k )
               k = k + 1
             enddo
           enddo
         enddo

      case (3)
         do i3 = msg%range( p_lower, 3 ), msg%range( p_upper, 3 )
           do i2 = msg%range( p_lower, 2 ), msg%range( p_upper, 2 )
             do i1 = msg%range( p_lower, 1 ), msg%range( p_upper, 1 )
               do i = 1, wall%f_dim
                 wall%f3( i, i1, i2, i3 ) = msg%buffer( k )
                 k = k + 1
               enddo
             enddo
           enddo
         enddo
    end select
  endif

end subroutine unpack_wall_msg
!---------------------------------------------------------------------------------------------------


!***************************************************************************************************
!***************************************************************************************************
!**                                                                                               **
!**                                Dynamic load balancing code                                    **
!**                                                                                               **
!***************************************************************************************************
!***************************************************************************************************


!---------------------------------------------------------------------------------------------------
!> Reshapes wall object when node grids change. 
!! No data is copied from previous wall object.
!! If dbg parameter is set then then the new vdf is initialized with NaNs, otherwise the wall
!! data is left undefined
!!
!! @param   this            t_wall object
!! @param   new_grid          t_grid object specifying new boundaries
!! @param   dbg             (optional) set this variable to initialize the new vdf with NaN values
subroutine reshape_wall_nocopy( this, new_grid, dbg )

  implicit none

  type (t_wall), intent(inout) :: this
    class( t_grid ), intent(in) :: new_grid

    ! Initialize the new array with NaN values
    logical,  optional, intent(in) :: dbg

  ! local variables
    integer :: i
  integer, dimension(this%x_dim + 1) :: lower, upper

  integer :: x_dim

  ! executable statements

  ! reshape wall if any data in it
  if (this%x_dim > 0) then

     x_dim = this%x_dim ! copy the value for simplicity

     lower(1)          = 1
        upper(1) = this%f_dim

        do i = 1, x_dim
            if ( i == this%idir ) then
     ! along wall normal use wall range
                this%nx(i)  = this%range(2) - this%range(1) + 1
                lower(i+1)  = this%range(1)
                upper(i+1)  = this%range(2)
            else
                this%nx(i) = new_grid%my_nx(3,i)
                lower(i+1) = 1          - this%gc_num(p_lower,i)
                upper(i+1) = this%nx(i) + this%gc_num(p_upper,i)
            endif
        enddo

     ! deallocate the old data and allocate new one
     select case (x_dim)
        case(1)
          call freemem(this%f1)
          call alloc( this%f1, lower, upper )

            if ( present(dbg) ) then
                call setnan( this%f1 )
          endif

        case(2)
          call freemem( this%f2 )
          call alloc( this%f2, lower, upper )

            if ( present(dbg) ) then
                call setnan( this%f2 )
          endif

        case(3)
          call freemem( this%f3 )
          call alloc( this%f3, lower, upper )

            if ( present(dbg) ) then
                call setnan( this%f3 )
          endif
     end select
  endif

end subroutine reshape_wall_nocopy
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! allocate message buffers for all the messages being sent/received in msg_patt
! and copy send data to msg buffers
!---------------------------------------------------------------------------------------------------
subroutine init_msg_buffers_wall( this, msg_patt, send, recv, my_node )

  implicit none

  type (t_wall), intent(inout) :: this
  type (t_msg_patt), intent(in) :: msg_patt
  type (t_vdf_msg), dimension(:), pointer :: send, recv
  integer, intent(in) :: my_node

    integer, parameter :: p_dyn_lb_tag = 1001
  integer :: i, j, bsize, self_send_id

  ! initiate send buffers if message is sent to other nodes
  self_send_id = -1
  do i = 1, msg_patt%n_send

        send(i)%type = p_vdf_replace
    send(i)%node = msg_patt%send(i)%node
        send(i)%tag  = p_dyn_lb_tag

    bsize = this%f_dim
    do j = 1, this%x_dim
      if ( j == this%idir ) then
        send(i)%range(p_lower,j) = this%range( p_lower )
        send(i)%range(p_upper,j) = this%range( p_upper )
      else
        send(i)%range(p_lower,j) = msg_patt%send(i)%cells(p_lower,j)
        send(i)%range(p_upper,j) = msg_patt%send(i)%cells(p_upper,j)
      endif
      bsize = bsize * ( send(i)%range(p_upper,j) -  send(i)%range(p_lower,j) + 1 )
    enddo

    call alloc( send(i)%buffer, (/bsize/) )
        send(i)%msg_size = bsize

    ! Pack data into message
    call pack_wall_msg( send(i), this )

    if (send(i)%node == my_node) then
      ! sanity check
      if ( self_send_id >= 1 ) then
        ERROR( 'More than 1 message to be sent to self' )
        call abort_program( p_err_invalid )
      endif

      ! Store id of message being sent to self
      self_send_id = i
    endif
  enddo

  ! initiate receive buffers
  do i = 1, msg_patt%n_recv

        recv(i)%type = p_vdf_replace
    recv(i)%node = msg_patt%recv(i)%node
        recv(i)%tag  = p_dyn_lb_tag

    bsize = this%f_dim
    do j = 1, this%x_dim
        if ( j == this%idir ) then
        recv(i)%range(p_lower,j) = this%range(p_lower)
        recv(i)%range(p_upper,j) = this%range(p_upper)
      else
        recv(i)%range(p_lower,j) = msg_patt%recv(i)%cells(p_lower,j)
        recv(i)%range(p_upper,j) = msg_patt%recv(i)%cells(p_upper,j)
      endif
      bsize = bsize * ( recv(i)%range(p_upper,j) -  recv(i)%range(p_lower,j) + 1 )
    enddo
        recv(i)%msg_size = bsize

    ! if message comes from local node just point the recv buffer to the send buffer that already
    ! has the data
    if ( recv(i)%node == my_node ) then

      ! sanity check
      if ( self_send_id < 1 ) then
        ERROR( 'No message is being sent to self' )
        call abort_program( p_err_invalid )
      endif

      recv(i)%buffer => send( self_send_id )%buffer
      send( self_send_id )%buffer => null()
    else
      ! otherwise allocate space for incoming message
      call alloc( recv(i)%buffer, (/bsize/) )
    endif

  enddo

end subroutine init_msg_buffers_wall
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!> Reshapes wall object when node grids change and redistribute the data through all nodes
!! @param   this        t_wall object
!! @param   old_lb      t_grid object describing the old node grid distribution
!! @param   new_grid    t_grid object describing the new node grid distribution
!! @param   no_co       t_node_conf object with parallel informantion
subroutine reshape_wall_copy( this, old_lb, new_grid, no_co, send_msg, recv_msg )

  use m_grid_parallel
    use m_vdf_comm

  implicit none

  type (t_wall), intent(inout) :: this
  class(t_grid), intent(in) :: old_lb, new_grid
  class ( t_node_conf ), intent(in) :: no_co
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  type (t_msg_patt) :: msg_patt
    integer :: i

  type ( t_vdf_msg ), dimension(:), pointer :: send, recv
    integer, dimension( 2, p_max_dim ) :: deltax, gc
    integer, dimension( p_max_dim ) :: move
    integer, dimension( 2 ) :: slice

#ifdef __DEBUG__
    SCR_ROOT('reshape_wall_copy started')
#endif

  ! inactive walls and 1D walls do not need to be reshaped
  if ( this%x_dim > 1 ) then

    ! get message pattern for nodes along the wall only

    ! slice position in global parallel partition
    if ( this%ibnd == p_lower ) then
      i = 1
    else
      i = nx(no_co,this%idir)
    endif

        ! Get message pattern, only consider nodes in the same slice of the
        ! parallel grid partition
        slice(1) = this % idir
        slice(2) = i
        call new_msg_patt( msg_patt, old_lb, new_grid, no_co, this%gc_num, slice )

    ! allocate send/recv buffers
    if ( msg_patt%n_send > 0 ) call alloc( send, (/msg_patt%n_send/) )
    if ( msg_patt%n_recv > 0 ) call alloc( recv, (/msg_patt%n_recv/) )

    ! Initiatlize message buffers and copy data to send/recv buffers
    call init_msg_buffers_wall( this, msg_patt, send, recv, no_co%my_aid() )

    ! start communication
    do i = 1, msg_patt%n_recv
       call irecv( recv(i), no_co )
    enddo

    do i = 1, msg_patt%n_send
       call isend( send(i), no_co )
    enddo

        ! reshape wall

#ifdef __DEBUG__
        ! (* debug *) - This sets the initial array values to NaN so the code can check if any value
        !               was left unset
        call reshape_wall_nocopy( this, new_grid, dbg = .true. )
#else
        ! (* production *)
        call reshape_wall_nocopy( this, new_grid )
#endif

    ! wait for receive messages to complete
    if ( msg_patt%n_recv > 0) then
       call wait( recv )
    endif

        ! copy data to reshaped wall
     do i = 1, msg_patt%n_recv
       call unpack_wall_msg( recv(i), this )
     enddo

    ! clear message pattern data
    call msg_patt%cleanup()

    ! clear recv buffers
    if ( msg_patt%n_recv > 0) then
      call cleanup( recv )
      call freemem( recv )
    endif

    ! wait for send messages to complete
    if ( msg_patt%n_send > 0) then
       call wait( send )

       ! clear send buffers
       call cleanup(send)
       call freemem( send )
    endif

  endif

    if ( this%x_dim > 0 ) then
    ! if grid changed in the same direction as the wall and
    ! the wall is upper wall then we need to shift the range
    if ( this%ibnd == p_upper ) then
            deltax( :, 1:this%x_dim ) = my_deltax( old_lb, new_grid )
      if ( deltax( p_upper, this%idir ) /= deltax( p_lower, this%idir ) ) then
        call shift_range_wall(this, deltax(p_upper,this%idir) - deltax(p_lower,this%idir))
      endif
    endif

        do i = 1, this % x_dim
            gc(p_lower,i) = this % gc_num(p_lower,i)
            gc(p_upper,i) = this % gc_num(p_upper,i)
            this % gc_num(p_lower,i) = 1
            this % gc_num(p_upper,i) = 2
        enddo

        move = 0
        call update_boundary_wall( this, p_vdf_replace, no_co, move, send_msg, recv_msg )

        do i = 1, this % x_dim
            this % gc_num(p_lower,i) = gc(p_lower,i)
            this % gc_num(p_upper,i) = gc(p_upper,i)
        enddo

  endif

#ifdef __DEBUG__
    !(*debug*)  
    call check_nan( this, "Checking Wall grid after reshape_wall_copy" )
  SCR_ROOT('reshape_wall_copy finished ok')
#endif

end subroutine reshape_wall_copy


!---------------------------------------------------------------------------------------------------
! Shift the wall range
!---------------------------------------------------------------------------------------------------
subroutine shift_range_wall( this, shift )

  implicit none

  type( t_wall ), intent(inout) :: this
  integer, intent(in) :: shift

  real(p_k_fld), dimension(:,:), pointer :: old_f1
  real(p_k_fld), dimension(:,:,:), pointer :: old_f2
  real(p_k_fld), dimension(:,:,:,:), pointer :: old_f3

  integer, dimension( p_max_dim + 1 ) :: lb, ub

  ! set new range
  this%range = this%range + shift

  ! reallocate the buffers with new boundaries and copy data
  select case ( this%x_dim )
    case(1)
       ! allocate temp buffer
       lb(1:2) = lbound(this%f1)
       ub(1:2) = ubound(this%f1)

       call alloc( old_f1, lb, ub )

       ! copy values to temp buffer
       old_f1 = this%f1

       ! reallocate wall buffer
       call freemem( this%f1 )

       lb(1 + this%idir) = this%range(1)
       ub(1 + this%idir) = this%range(2)

       call alloc( this%f1, lb, ub )

       ! copy values back to wall buffer
       this%f1 = old_f1

       ! cleanup temp buffer
       call freemem( old_f1 )

    case(2)

       ! allocate temp buffer
       lb(1:3) = lbound(this%f2)
       ub(1:3) = ubound(this%f2)
       call alloc( old_f2, lb, ub )

       ! copy values to temp buffer
       old_f2 = this%f2

       ! reallocate wall buffer
       call freemem( this%f2 )

       lb(1 + this%idir) = this%range(1)
       ub(1 + this%idir) = this%range(2)

       call alloc( this%f2, lb, ub )

       ! copy values back to wall buffer
       this%f2 = old_f2

       ! cleanup temp buffer
       call freemem( old_f2 )

    case(3)
       ! allocate temp buffer
       lb(1:4) = lbound(this%f3)
       ub(1:4) = ubound(this%f3)

       call alloc( old_f3, lb, ub )

       ! copy values to temp buffer
       old_f3 = this%f3

       ! reallocate wall buffer
       call freemem( this%f3 )

       lb(1 + this%idir) = this%range(1)
       ub(1 + this%idir) = this%range(2)

       call alloc( this%f3, lb, ub )

       ! copy values back to wall buffer
       this%f3 = old_f3

       ! cleanup temp buffer
       call freemem( old_f3 )

  end select


end subroutine shift_range_wall
!---------------------------------------------------------------------------------------------------



end module m_wall_comm

