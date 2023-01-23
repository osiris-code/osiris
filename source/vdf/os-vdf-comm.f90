!#define __DEBUG__ 1

#include "os-config.h"
#include "os-preprocess.fpp"

module m_vdf_comm
    
#include "memory/memory.h"

use m_vdf_define, only : t_vdf
use m_node_conf
use m_grid_define

use m_parameters

implicit none

private

type t_vdf_msg

! comm type (add or replace)
integer :: type

! cell range involved in communication
integer, dimension( 2, p_max_dim ) :: range

! node communicating to/from
integer :: node = -1

! message tag (optional)
integer :: tag = 0

! message request id
integer :: handle = MPI_REQUEST_NULL

! pointer to field data
real(p_k_fld), dimension(:), pointer :: buffer => null()

! message size ( buffer may be larger than message size )
integer :: msg_size

! parameters for use with the tiles module
! step, dimension in tiled boundary update
integer :: step, dim

! add potential to have a linked list of these
type(t_vdf_msg), pointer :: next => null()

end type t_vdf_msg

! module variables for normal vdf comm
! type( t_vdf_msg ), dimension(2), save, private :: send_msg, recv_msg


interface update_boundary
    module procedure update_boundary_vdf
end interface

interface reshape_copy
    module procedure reshape_vdf_copy
end interface

interface reshape_nocopy
    module procedure reshape_vdf_nocopy
end interface

interface alloc
    module procedure alloc_vdf_msg
    module procedure alloc_1d_vdf_msg
    module procedure alloc_bound_1d_vdf_msg
end interface

interface freemem
    module procedure freemem_vdf_msg
    module procedure freemem_1d_vdf_msg
end interface

interface wait
    module procedure wait_vdf_msg
    module procedure wait_vdf_msg_array
end interface

interface isend
    module procedure isend_vdf
    module procedure isend_vdf_msg
end interface

interface isend_wait
    module procedure wait_isend_vdf
end interface

interface irecv
    module procedure irecv_vdf
    module procedure irecv_vdf_msg
end interface

interface irecv_wait
    module procedure wait_irecv_vdf
end interface

interface cleanup
    module procedure cleanup_vdf_msg_1
    module procedure cleanup_vdf_msg_array
end interface

interface update_periodic_vdf
    module procedure update_periodic_vdf
end interface

interface pack_vdf_msg
    module procedure pack_vdf_msg
end interface

interface unpack_vdf_msg
    module procedure unpack_vdf_msg
end interface

public :: t_vdf_msg
public :: update_boundary, reshape_copy, reshape_nocopy
public :: alloc, freemem, wait, isend, isend_wait, irecv, irecv_wait
public :: cleanup, update_periodic_vdf
public :: pack_vdf_msg, unpack_vdf_msg

contains

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines for t_vdf_msg

#define __TYPE__ type( t_vdf_msg )
#define __TYPE_STR__ "t_vdf_msg"
#define FNAME( a )  a ## _vdf_msg
#define __MAX_DIM__ 1
#include "memory/mem-template.h"

!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Update single node periodic gc / grid values for a vdf along a given direction
!---------------------------------------------------------------------------------------------------
subroutine update_periodic_vdf( vdf, dim, update_type )
    
    implicit none
    
    class( t_vdf ), intent(inout) :: vdf
    integer, intent(in) :: dim, update_type
    
    integer, dimension(2, p_max_dim) :: lrange, urange
    integer, dimension( p_max_dim ) :: shift
    
    integer :: i, i1, i2, i3
    real( p_k_fld ) :: val
    
    ! Number of cells to shift
    do i = 1, vdf%x_dim_
        if ( i == dim ) then
            shift(i) = vdf%nx_(i)
        else
            shift(i) = 0
        endif
    enddo
    
    do i = 1, vdf%x_dim_
        if ( i == dim ) then
            ! Select only guard cells along update direction
            lrange( p_lower, i ) = 1         - vdf%gc_num_( p_lower, i )
            lrange( p_upper, i ) = 0
            urange( p_lower, i ) = vdf%nx_(i) + 1
            urange( p_upper, i ) = vdf%nx_(i) + vdf%gc_num_( p_upper, i )
        else
            ! full cell range along other directions
            lrange( p_lower, i ) = 1         - vdf%gc_num_( p_lower, i )
            lrange( p_upper, i ) = vdf%nx_(i) + vdf%gc_num_( p_upper, i )
            urange( p_lower, i ) = lrange( p_lower, i )
            urange( p_upper, i ) = lrange( p_upper, i )
        endif
    enddo
    
    select case ( update_type )
        
    case ( p_vdf_add )
        
        ! Add grid values to corresponding gc values and store in both grid and gc
        select case ( vdf%x_dim_ )
        case (1)
            ! update lower gc
            do i1 = lrange( p_lower, 1 ), lrange( p_upper, 1 )
                do i = 1, vdf%f_dim_
                    val =  vdf%f1( i, i1 ) + vdf%f1( i, i1 + shift(1) )
                    vdf%f1( i, i1 )             = val
                    vdf%f1( i, i1 + shift(1) ) = val
                enddo
            enddo
            
            ! update upper gc
            do i1 = urange( p_lower, 1 ), urange( p_upper, 1 )
                do i = 1, vdf%f_dim_
                    val = vdf%f1( i, i1 ) + vdf%f1( i, i1 - shift(1) )
                    vdf%f1( i, i1 - shift(1) ) = val
                    vdf%f1( i, i1 )             = val
                enddo
            enddo
            
        case (2)
            ! update lower gc
            do i2 = lrange( p_lower, 2 ), lrange( p_upper, 2 )
                do i1 = lrange( p_lower, 1 ), lrange( p_upper, 1 )
                    do i = 1, vdf%f_dim_
                        val =  vdf%f2( i, i1, i2 ) + vdf%f2( i, i1 + shift(1), i2 + shift(2) )
                        vdf%f2( i, i1, i2 )                         = val
                        vdf%f2( i, i1 + shift(1), i2 + shift(2) ) = val
                    enddo
                enddo
            enddo
            
            ! update upper gc
            do i2 = urange( p_lower, 2 ), urange( p_upper, 2 )
                do i1 = urange( p_lower, 1 ), urange( p_upper, 1 )
                    do i = 1, vdf%f_dim_
                        val = vdf%f2( i, i1, i2 ) + vdf%f2( i, i1 - shift(1), i2 - shift(2) )
                        vdf%f2( i, i1 - shift(1), i2 - shift(2) ) = val
                        vdf%f2( i, i1, i2 )                         = val
                    enddo
                enddo
            enddo
            
        case (3)
            ! update lower gc
            do i3 = lrange( p_lower, 3 ), lrange( p_upper, 3 )
                do i2 = lrange( p_lower, 2 ), lrange( p_upper, 2 )
                    do i1 = lrange( p_lower, 1 ), lrange( p_upper, 1 )
                        do i = 1, vdf%f_dim_
                            val =  vdf%f3( i,i1,i2,i3 ) + vdf%f3( i, i1+shift(1), i2+shift(2), i3+shift(3) )
                            vdf%f3( i, i1, i2, i3 )                                     = val
                            vdf%f3( i, i1 + shift(1), i2 + shift(2), i3 + shift(3) ) = val
                        enddo
                    enddo
                enddo
            enddo
            
            ! update upper gc
            do i3 = urange( p_lower, 3 ), urange( p_upper, 3 )
                do i2 = urange( p_lower, 2 ), urange( p_upper, 2 )
                    do i1 = urange( p_lower, 1 ), urange( p_upper, 1 )
                        do i = 1, vdf%f_dim_
                            val =  vdf%f3( i,i1,i2,i3 ) + vdf%f3( i, i1-shift(1),i2-shift(2),i3-shift(3) )
                            vdf%f3( i, i1 - shift(1), i2 - shift(2), i3 - shift(3) ) = val
                            vdf%f3( i, i1, i2, i3 )                                     = val
                        enddo
                    enddo
                enddo
            enddo
            
        end select
        
    case ( p_vdf_replace )
        
        ! Copy gc values from opposing grid values
        select case ( vdf%x_dim_ )
        case (1)
            ! update lower gc
            do i1 = lrange( p_lower, 1 ), lrange( p_upper, 1 )
                do i = 1, vdf%f_dim_
                    vdf%f1( i, i1 ) = vdf%f1( i, i1 + shift(1) )
                enddo
            enddo
            
            ! update upper gc
            do i1 = urange( p_lower, 1 ), urange( p_upper, 1 )
                do i = 1, vdf%f_dim_
                    vdf%f1( i, i1 ) = vdf%f1( i, i1 - shift(1) )
                enddo
            enddo
            
        case (2)
            ! update lower gc
            do i2 = lrange( p_lower, 2 ), lrange( p_upper, 2 )
                do i1 = lrange( p_lower, 1 ), lrange( p_upper, 1 )
                    do i = 1, vdf%f_dim_
                        vdf%f2( i, i1, i2 ) = vdf%f2( i, i1 + shift(1), i2 + shift(2) )
                    enddo
                enddo
            enddo
            
            ! update upper gc
            do i2 = urange( p_lower, 2 ), urange( p_upper, 2 )
                do i1 = urange( p_lower, 1 ), urange( p_upper, 1 )
                    do i = 1, vdf%f_dim_
                        vdf%f2( i, i1, i2 ) = vdf%f2( i, i1 - shift(1), i2 - shift(2) )
                    enddo
                enddo
            enddo
            
        case (3)
            ! update lower gc
            do i3 = lrange( p_lower, 3 ), lrange( p_upper, 3 )
                do i2 = lrange( p_lower, 2 ), lrange( p_upper, 2 )
                    do i1 = lrange( p_lower, 1 ), lrange( p_upper, 1 )
                        do i = 1, vdf%f_dim_
                            vdf%f3( i, i1, i2, i3 ) = vdf%f3( i, i1 + shift(1), &
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
                        do i = 1, vdf%f_dim_
                            vdf%f3( i, i1, i2, i3 ) = vdf%f3( i, i1 - shift(1), &
                            i2 - shift(2), &
                            i3 - shift(3) )
                        enddo
                    enddo
                enddo
            enddo
            
        end select
        
    end select
    
end subroutine update_periodic_vdf
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Pack vdf data into message buffer
!---------------------------------------------------------------------------------------------------
subroutine pack_vdf_msg( msg, vdf, k )
    
    implicit none
    
    type( t_vdf_msg ), intent(inout) :: msg
    class( t_vdf ), intent(in) :: vdf
    
    integer, intent(inout) :: k
    
    integer :: i, i1, i2, i3
    
    select case ( vdf%x_dim_ )
    case (1)
        do i1 = msg%range( p_lower, 1 ), msg%range( p_upper, 1 )
            do i = 1, vdf%f_dim_
                msg%buffer( k ) = vdf%f1( i, i1 )
                k = k + 1
            enddo
        enddo
        
    case (2)
        do i2 = msg%range( p_lower, 2 ), msg%range( p_upper, 2 )
            do i1 = msg%range( p_lower, 1 ), msg%range( p_upper, 1 )
                do i = 1, vdf%f_dim_
                    msg%buffer( k ) = vdf%f2( i, i1, i2 )
                    k = k + 1
                enddo
            enddo
        enddo
        
    case (3)
        do i3 = msg%range( p_lower, 3 ), msg%range( p_upper, 3 )
            do i2 = msg%range( p_lower, 2 ), msg%range( p_upper, 2 )
                do i1 = msg%range( p_lower, 1 ), msg%range( p_upper, 1 )
                    do i = 1, vdf%f_dim_
                        msg%buffer( k ) = vdf%f3( i, i1, i2, i3 )
                        k = k + 1
                    enddo
                enddo
            enddo
        enddo
    end select
    
    
end subroutine pack_vdf_msg
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Unapack message buffer data into vdf
!---------------------------------------------------------------------------------------------------
subroutine unpack_vdf_msg( msg, vdf, k )
    
    implicit none

    type( t_vdf_msg ), intent(inout) :: msg
    class( t_vdf ), intent(inout) :: vdf
    
    integer, intent(inout) :: k
    
    integer :: i, i1, i2, i3
    
    select case ( msg%type )
        
    case ( p_vdf_add )
        
        ! Add message values to local values
        
        select case ( vdf%x_dim_ )
        case (1)
            do i1 = msg%range( p_lower, 1 ), msg%range( p_upper, 1 )
                do i = 1, vdf%f_dim_
                    vdf%f1( i, i1 ) = vdf%f1( i, i1 ) + msg%buffer( k )
                    k = k + 1
                enddo
            enddo
            
        case (2)
            do i2 = msg%range( p_lower, 2 ), msg%range( p_upper, 2 )
                do i1 = msg%range( p_lower, 1 ), msg%range( p_upper, 1 )
                    do i = 1, vdf%f_dim_
                        vdf%f2( i, i1, i2 ) = vdf%f2( i, i1, i2 ) + msg%buffer( k )
                        k = k + 1
                    enddo
                enddo
            enddo
            
        case (3)
            do i3 = msg%range( p_lower, 3 ), msg%range( p_upper, 3 )
                do i2 = msg%range( p_lower, 2 ), msg%range( p_upper, 2 )
                    do i1 = msg%range( p_lower, 1 ), msg%range( p_upper, 1 )
                        do i = 1, vdf%f_dim_
                            vdf%f3( i, i1, i2, i3 ) = vdf%f3( i, i1, i2, i3 ) + msg%buffer( k )
                            k = k + 1
                        enddo
                    enddo
                enddo
            enddo
        end select
        
    case ( p_vdf_replace )
        
        ! replace local values with message values
        
        select case ( vdf%x_dim_ )
        case (1)
            do i1 = msg%range( p_lower, 1 ), msg%range( p_upper, 1 )
                do i = 1, vdf%f_dim_
                    vdf%f1( i, i1 ) = msg%buffer( k )
                    k = k + 1
                enddo
            enddo
            
        case (2)
            do i2 = msg%range( p_lower, 2 ), msg%range( p_upper, 2 )
                do i1 = msg%range( p_lower, 1 ), msg%range( p_upper, 1 )
                    do i = 1, vdf%f_dim_
                        vdf%f2( i, i1, i2 ) = msg%buffer( k )
                        k = k + 1
                    enddo
                enddo
            enddo
            
        case (3)
            do i3 = msg%range( p_lower, 3 ), msg%range( p_upper, 3 )
                do i2 = msg%range( p_lower, 2 ), msg%range( p_upper, 2 )
                    do i1 = msg%range( p_lower, 1 ), msg%range( p_upper, 1 )
                        do i = 1, vdf%f_dim_
                            vdf%f3( i, i1, i2, i3 ) = msg%buffer( k )
                            k = k + 1
                        enddo
                    enddo
                enddo
            enddo
        end select
        
    case default
        ERROR('unpack_vdf_msg called with an invalid msg%type')
        call abort_program()
        
    end select
    
end subroutine unpack_vdf_msg
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! cleanup msg buffers for a linked list
!---------------------------------------------------------------------------------------------------
recursive subroutine cleanup_vdf_msg_1( msg )

    implicit none

    type( t_vdf_msg ), intent(inout) :: msg

    if (associated(msg%next)) then
        call cleanup_vdf_msg_1(msg%next)
        deallocate(msg%next)
    endif

    call freemem( msg%buffer )

    msg%handle = MPI_REQUEST_NULL
    msg%msg_size = 0

end subroutine cleanup_vdf_msg_1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! clear msg buffers
!---------------------------------------------------------------------------------------------------
subroutine cleanup_vdf_msg_array( vdf_msg )

    implicit none

    type (t_vdf_msg), dimension(:), intent(inout) :: vdf_msg

    integer :: i

    do i=1, size(vdf_msg)
        
        vdf_msg(i)%handle = MPI_REQUEST_NULL
        vdf_msg(i)%msg_size = 0
        call freemem( vdf_msg(i)%buffer )
        
    enddo

end subroutine cleanup_vdf_msg_array
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Post receive for message
!---------------------------------------------------------------------------------------------------
subroutine irecv_vdf_msg( recv, no_co )

    implicit none

    type (t_vdf_msg), intent(inout)  :: recv
    class (t_node_conf ), intent(in) :: no_co

    integer :: ierr
    integer :: my_node
    integer :: mpi_type

    mpi_type = mpi_real_type( p_k_fld )
    my_node = no_co % my_aid()

    if ( recv%node /= my_node ) then
        call mpi_irecv( recv%buffer, recv%msg_size, mpi_type, recv%node - 1, recv%tag, comm( no_co ), &
        recv%handle, ierr )
        
        if (ierr/=0) then
            ERROR("MPI error")
            call abort_program( p_err_mpi )
        endif
        !    print *, '[', my_node, '] <- ', recv%node , ' ', recv%msg_size, ' reals, handle = ', recv%handle, ' tag = ', recv%tag
    else
        recv%handle = MPI_REQUEST_NULL
    endif

end subroutine irecv_vdf_msg
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Post send for message
!---------------------------------------------------------------------------------------------------
subroutine isend_vdf_msg( send, no_co )

    implicit none

    type (t_vdf_msg), intent(inout)  :: send
    class (t_node_conf ), intent(in) :: no_co

    integer :: ierr
    integer :: my_node
    integer :: mpi_type

    mpi_type = mpi_real_type( p_k_fld )
    my_node = no_co % my_aid()

    if ( send%node /= my_node ) then
        call mpi_isend( send%buffer, send%msg_size, mpi_type, send%node-1, send%tag, comm(no_co), &
        send%handle, ierr )
        if (ierr/=0) then
            ERROR("MPI error")
            call abort_program( p_err_mpi )
        endif
        
        !     print *, '[', my_node, '] -> ', send%node , ' ', send%msg_size, ' reals, handle = ', send%handle, ' tag = ', send%tag
        
    else
        send%handle = MPI_REQUEST_NULL
    endif

end subroutine isend_vdf_msg
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Wait for message set to finish
!---------------------------------------------------------------------------------------------------
subroutine wait_vdf_msg_array( vdf_msg )

    implicit none

    type (t_vdf_msg), dimension(:), intent(inout) :: vdf_msg

    integer :: i, j, ierr
    integer, dimension(:,:), pointer :: status
    integer, dimension(:), pointer :: request
    integer :: count

    count = size( vdf_msg )

    call alloc(status,  (/ MPI_STATUS_SIZE, count /) )
    call alloc(request, (/ count /) )

    j = 0
    do i = 1, count
        if (vdf_msg(i)%handle /= MPI_REQUEST_NULL) then
            j = j+1
            request(j) = vdf_msg(i)%handle
        endif
    enddo

    ! Check if there are any handles waiting
    if (j > 0) then
        call mpi_waitall(j, request, status, ierr)
        if (ierr/=0) then
            ERROR("MPI error")
            call abort_program( p_err_mpi )
        endif
    endif

    ! clear message handles
    do i = 1, count
        vdf_msg(i)%handle = MPI_REQUEST_NULL
    enddo

    call freemem(request)
    call freemem(status)

end subroutine wait_vdf_msg_array
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Wait for single message to finish
!---------------------------------------------------------------------------------------------------
subroutine wait_vdf_msg( vdf_msg )

    implicit none

    type (t_vdf_msg), intent(inout) :: vdf_msg

    integer :: ierr
    integer, dimension(mpi_status_size):: stat

    ! if message has no valid handle return quietly
    if ( vdf_msg%handle /= MPI_REQUEST_NULL ) then
        ! wait for message to complete
        call mpi_wait( vdf_msg%handle, stat, ierr )
        if (ierr/=0) then
            ERROR("MPI error")
            call abort_program( p_err_mpi )
        endif
        
        ! clear handle
        vdf_msg%handle = MPI_REQUEST_NULL
    endif

end subroutine wait_vdf_msg
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Update boundary values of a vdf on non-physical boundaries ( communication and parallel )
!---------------------------------------------------------------------------------------------------
subroutine update_boundary_vdf( vdf, update_type, no_co, send_msg, recv_msg, move_num_ )

    implicit none

    class( t_vdf ), intent(inout) :: vdf
    integer,                 intent(in) :: update_type
    class( t_node_conf ),             intent(in) :: no_co
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
    integer, dimension(:), intent(in), optional :: move_num_

    integer, dimension(p_max_dim) :: move_num
    integer :: i, my_node, lneighbor, uneighbor

    if ( present( move_num_ ) ) then
        move_num( 1 : vdf%x_dim_ ) = move_num_( 1 : vdf%x_dim_ )
    else
        move_num = 0
    endif

    my_node = no_co % my_aid()

    do i = 1, vdf%x_dim_
        
        if ( my_node == neighbor(no_co,p_lower,i) ) then
            
            !single node periodic
            call update_periodic_vdf( vdf, i, update_type )
            
        else
            
            ! communication with other nodes
            lneighbor = neighbor( no_co, p_lower, i )
            uneighbor = neighbor( no_co, p_upper, i )
            
            ! post receives
            if ( lneighbor > 0 ) call irecv_vdf( vdf, update_type, i, p_lower, no_co, move_num(i), recv_msg )
            if ( uneighbor > 0 ) call irecv_vdf( vdf, update_type, i, p_upper, no_co, move_num(i), recv_msg )
            
            ! post sends
            ! Note: sends MUST be posted in the opposite order of receives to account for a 2 node
            ! periodic partition, where 1 node sends 2 messages to the same node
            if ( uneighbor > 0 ) call isend_vdf( vdf, update_type, i, p_upper, no_co, move_num(i), send_msg )
            if ( lneighbor > 0 ) call isend_vdf( vdf, update_type, i, p_lower, no_co, move_num(i), send_msg )
            
            ! wait for receives and unpack data
            if ( lneighbor > 0 ) call wait_irecv_vdf( vdf, p_lower, recv_msg )
            if ( uneighbor > 0 ) call wait_irecv_vdf( vdf, p_upper, recv_msg )
            
            ! wait for sends to complete
            if ( uneighbor > 0 ) call wait_isend_vdf( vdf, p_upper, send_msg )
            if ( lneighbor > 0 ) call wait_isend_vdf( vdf, p_lower, send_msg )
            
        endif
        
    enddo


end subroutine update_boundary_vdf
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! post send message for boundary communication
!---------------------------------------------------------------------------------------------------
subroutine isend_vdf( vdf, update_type, dim, bnd, no_co, nx_move, send_msg, vdf_b )

    implicit none

    class( t_vdf ), intent(in) :: vdf
    class( t_vdf ), intent(in), optional :: vdf_b

    integer, intent(in) :: update_type, dim, bnd
    class( t_node_conf ), intent(in) :: no_co

    integer, intent(in) :: nx_move
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg

    integer, dimension(2,p_max_dim) :: ns
    integer :: i, bsize, msg_size

    ! Sanity checks
    if ( send_msg(bnd)%handle /= MPI_REQUEST_NULL ) then
        ERROR( 'isend_msg_vdf called when message still waiting to complete' )
        call abort_program( p_err_invalid )
    endif

    if ( neighbor( no_co, bnd, dim ) < 0 ) then
        write(0,*) 'isend_msg_vdf called without a valid neighbor, dim, bnd = ',dim, bnd, &
        trim(__FILE__), ':', __LINE__
        call abort_program( p_err_invalid )
    endif

    ! Get range of cells to send
    msg_size = vdf%f_dim_
    do i = 1, vdf%x_dim_
        if ( i == dim ) then
            if ( update_type == p_vdf_add ) then
                select case ( bnd )
                case ( p_upper )
                    ns( p_lower, i ) = vdf%nx_(i) + 1 - vdf%gc_num_( p_lower, i )
                    ns( p_upper, i ) = vdf%nx_(i) + vdf%gc_num_( p_upper, i )
                case ( p_lower )
                    ns( p_lower, i ) = 1 - vdf%gc_num_( p_lower, i )
                    ns( p_upper, i ) = vdf%gc_num_( p_upper, i )
                end select
            else
                select case ( bnd )
                case ( p_upper )
                    ns( p_lower, i ) = vdf%nx_(i) + 1 - vdf%gc_num_( p_lower, i )
                    ns( p_upper, i ) = vdf%nx_(i) - nx_move
                case ( p_lower )
                    ns( p_lower, i ) = 1 - nx_move
                    ns( p_upper, i ) = vdf%gc_num_( p_upper, i )
                end select
            endif
        else
            ns( p_lower, i ) = vdf%lbound( i+1 )
            ns( p_upper, i ) = vdf%ubound( i+1 )
        endif
        msg_size = msg_size * (ns( p_upper, i ) - ns( p_lower,i ) + 1 )
    enddo

    ! If second vdf is present double message size
    if ( present( vdf_b ) ) msg_size = msg_size * 2

    send_msg(bnd)%type                      = update_type ! unnecessary
    send_msg(bnd)%range( :, 1 : vdf%x_dim_ ) = ns ( :, 1 : vdf%x_dim_ )
    send_msg(bnd)%msg_size                  = msg_size
    send_msg(bnd)%node                      = neighbor( no_co, bnd, dim )
    send_msg(bnd)%tag                       = bnd

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
    i = 1
    call pack_vdf_msg( send_msg(bnd), vdf, i )
    if ( present( vdf_b ) )  call pack_vdf_msg( send_msg(bnd), vdf_b, i )

    ! post send
    call isend_vdf_msg( send_msg(bnd), no_co )

end subroutine isend_vdf
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Wait for a BC vdf send message to complete
!---------------------------------------------------------------------------------------------------
subroutine wait_isend_vdf( this, bnd, send_msg )

    implicit none

    class( t_vdf ), intent(in) :: this
    integer, intent(in) :: bnd
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg

    call wait_vdf_msg( send_msg( bnd ) )

end subroutine wait_isend_vdf
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! post recv message for boundary communication
!---------------------------------------------------------------------------------------------------
subroutine irecv_vdf( vdf, update_type, dim, bnd, no_co, nx_move, recv_msg, vdf_b )

    implicit none

    class( t_vdf ), intent(in) :: vdf
    integer, intent(in) :: update_type, dim, bnd
    class( t_node_conf ), intent(in) :: no_co
    integer, intent(in) :: nx_move
    type(t_vdf_msg), dimension(2), intent(inout) :: recv_msg
    class( t_vdf ), intent(in), optional :: vdf_b


    integer, dimension(2,p_max_dim) :: nr
    integer :: i, bsize, msg_size

    ! Check if message is still waiting to complete
    if ( recv_msg( bnd )%handle /= MPI_REQUEST_NULL ) then
        ERROR( 'irecv_vdf called when message still waiting to complete' )
        call abort_program( p_err_invalid )
    endif

    ! Sanity check
    if ( neighbor( no_co, bnd, dim ) < 0 ) then
        write(0,*) 'irecv_vdf called without a valid neighbor, dim, bnd = ',dim, bnd, &
        trim(__FILE__), ':', __LINE__
        call abort_program( p_err_invalid )
    endif

    ! Get range of cells to receive
    msg_size = vdf%f_dim_
    do i = 1, vdf%x_dim_
        if ( i == dim ) then
            if ( update_type == p_vdf_add ) then
                select case ( bnd )
                case ( p_upper )
                    nr( p_lower, i ) = vdf%nx_(i) + 1 - vdf%gc_num_( p_lower, i )
                    nr( p_upper, i ) = vdf%nx_(i) + vdf%gc_num_( p_upper, i )
                case ( p_lower )
                    nr( p_lower, i ) = 1 - vdf%gc_num_( p_lower, i )
                    nr( p_upper, i ) = vdf%gc_num_( p_upper, i )
                end select
            else
                select case ( bnd )
                case ( p_upper )
                    nr( p_lower, i ) = vdf%nx_(i) + 1 - nx_move
                    nr( p_upper, i ) = vdf%nx_(i) + vdf%gc_num_( p_upper, i )
                case ( p_lower )
                    nr( p_lower, i ) = 1 - vdf%gc_num_( p_lower, i )
                    nr( p_upper, i ) = 0 - nx_move
                end select
            endif
        else
            nr( p_lower, i ) = vdf%lbound( i+1 )
            nr( p_upper, i ) = vdf%ubound( i+1 )
        endif
        msg_size = msg_size * (nr( p_upper, i ) - nr( p_lower,i ) + 1 )
    enddo

    ! (* debug *)
    if ( msg_size == 0 ) then
        ERROR('Invalid message size, aborting.')
        call abort_program( p_err_invalid )
    endif

    ! If second vdf is present double message size
    if ( present( vdf_b ) ) msg_size = msg_size * 2

    ! Get send buffer and grow it if necessary
    recv_msg(bnd)%type                       = update_type
    recv_msg(bnd)%range( :, 1 : vdf%x_dim_ ) = nr ( :, 1 : vdf%x_dim_ )
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
    call irecv_vdf_msg( recv_msg(bnd), no_co )

end subroutine irecv_vdf
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Wait for a BC vdf recv message to complete
!---------------------------------------------------------------------------------------------------
subroutine wait_irecv_vdf( vdf, bnd, recv_msg, vdf_b )

    implicit none

    class( t_vdf ), intent(inout) :: vdf
    integer, intent(in) :: bnd
    type(t_vdf_msg), dimension(2), intent(inout) :: recv_msg
    class( t_vdf ), intent(inout), optional :: vdf_b

    integer :: i

    ! wait for message to complete
    call wait_vdf_msg( recv_msg( bnd ) )

    ! unpack message data
    i = 1
    call unpack_vdf_msg( recv_msg( bnd ), vdf, i )
    if ( present( vdf_b ) ) call unpack_vdf_msg( recv_msg( bnd ), vdf_b, i )

end subroutine wait_irecv_vdf
!---------------------------------------------------------------------------------------------------


!***************************************************************************************************
!***************************************************************************************************
!**                                                                                               **
!**                                Dynamic load balancing code                                    **
!**                                                                                               **
!***************************************************************************************************
!***************************************************************************************************

!---------------------------------------------------------------------------------------------------
!> Reshape vdf object when node grids change. No data is copied from previous vdf object.
!! If dbg parameter is set then then the new vdf is initialized with NaNs, otherwise the vdf
!! data is left undefined
!!
!! @param   this        vdf object
!! @param   new_grid    grid object with the new parallel grid parameters
!! @param   dbg         (optional) set this variable to initialize the new vdf with NaN values
subroutine reshape_vdf_nocopy( this, new_grid, dbg )

    implicit none

    type (t_vdf), intent(inout) :: this
    class( t_grid ), intent(in) :: new_grid

    ! Initialize the new array with NaN values
    logical,  optional, intent(in) :: dbg

    ! local variables
    integer :: i
    integer, dimension(p_max_dim+1) :: lower, upper

    type(c_ptr) :: tmp
    integer, dimension(1) :: shape
  
    integer :: x_dim

    ! executable statements

    ! reshape vdf if any data in it
    if (this%x_dim_ > 0) then
        
        x_dim = this%x_dim_ ! copy the value for simplicity
        
        lower(1) = 1
        upper(1) = this%f_dim_
        
        do i = 1, x_dim
            this%nx_(i) = new_grid%my_nx(3,i)
            lower(i+1)  = 1           - this%gc_num_(p_lower,i)
            upper(i+1)  = this%nx_(i) + this%gc_num_(p_upper,i)
        enddo
        
        ! deallocate the old data and allocate new one
        select case (x_dim)
        case(1)
            call freemem(this%f1)
            call alloc( this%f1, lower, upper )
            tmp = c_loc( this % f1( lower(1), lower(2) ) )
            
            if ( present(dbg) ) then
                call setnan( this%f1 )
            endif
            
        case(2)
            
            call freemem( this%f2 )
            call alloc( this%f2, lower, upper )
            tmp = c_loc( this % f2( lower(1), lower(2), lower(3) ) )

            if ( present(dbg) ) then
                call setnan( this%f2 )
            endif
            
        case(3)
            
            call freemem( this%f3 )
            call alloc( this%f3, lower, upper )
            tmp = c_loc( this % f3( lower(1), lower(2), lower(3), lower(4) ) )

            if ( present(dbg) ) then
                call setnan( this%f3 )
            endif
            
        end select
        
    endif

    ! Associate buffer pointer with allocated data
    shape(1) = this%f_dim_
    do i = 1, x_dim
        shape(1) = shape(1) * (upper(i+1)-lower(i+1)+1)
    enddo
    call c_f_pointer( tmp, this % buffer, shape )

end subroutine reshape_vdf_nocopy
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! allocate message buffers for all the messages being sent/received in msg_patt
! and copy send data to msg buffers
!---------------------------------------------------------------------------------------------------
subroutine init_msg_buffers( this, msg_patt, send, recv, my_node )
    implicit none

    ! dummy variables
    type (t_vdf), intent(inout) :: this
    type (t_msg_patt), intent(in) :: msg_patt

    type (t_vdf_msg), dimension(:), pointer :: send, recv

    integer, intent(in) :: my_node

    ! local variables
    integer :: i, j, bsize, self_send_id
    integer, parameter :: p_dyn_lb_tag = 1001


    ! initiate send buffers if message is sent to other nodes
    self_send_id = -1
    do i = 1, msg_patt%n_send
        
        send(i)%type = p_vdf_replace
        send(i)%node = msg_patt%send(i)%node
        send(i)%tag  = p_dyn_lb_tag
        
        bsize = this%f_dim_
        do j = 1, this%x_dim_
            send(i)%range(p_lower,j) = msg_patt%send(i)%cells(p_lower,j)
            send(i)%range(p_upper,j) = msg_patt%send(i)%cells(p_upper,j)
            bsize = bsize * ( send(i)%range(p_upper,j) -  send(i)%range(p_lower,j) + 1 )
        enddo
        
        call alloc( send(i)%buffer, (/bsize/) )
        
        ! In this case the msg_size is the same as the buffer size
        send(i)%msg_size = bsize
        
        ! Pack data into message
        j = 1
        call pack_vdf_msg( send(i), this, j )
        
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
        
        bsize = this%f_dim_
        do j = 1, this%x_dim_
            recv(i)%range(p_lower,j) = msg_patt%recv(i)%cells(p_lower,j)
            recv(i)%range(p_upper,j) = msg_patt%recv(i)%cells(p_upper,j)
            bsize = bsize * ( recv(i)%range(p_upper,j) -  recv(i)%range(p_lower,j) + 1 )
        enddo
        
        ! In this case the msg_size is the same as the buffer size
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

end subroutine init_msg_buffers
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
! reshape vdf object when node grids change and redistribute the data through all nodes
! Waits for all messages to arive and then unpacks them.
!---------------------------------------------------------------------------------------------------
subroutine reshape_vdf_copy( this, old_lb, new_grid, no_co, send_msg, recv_msg )

    use m_grid_parallel

    implicit none

    ! dummy variables
    type (t_vdf), intent(inout) :: this
    class (t_grid ), intent(in) :: old_lb, new_grid
    class (t_node_conf ), intent(in) :: no_co
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

    ! local variables
    type (t_msg_patt) :: msg_patt
    integer :: i, j

    type (t_vdf_msg), dimension(:), pointer :: send, recv

    ! get message pattern
    call new_msg_patt( msg_patt, old_lb, new_grid, no_co, this%gc_num_ )

    ! allocate send/recv buffers
    if ( msg_patt%n_send > 0 ) call alloc( send, (/msg_patt%n_send/) )
    if ( msg_patt%n_recv > 0 ) call alloc( recv, (/msg_patt%n_recv/) )

    ! Initiatlize message buffers and copy data to send/recv buffers
    call init_msg_buffers( this, msg_patt, send, recv, no_co%my_aid() )

    ! start communication
    do i = 1, msg_patt%n_recv
        call irecv_vdf_msg( recv(i), no_co )
    enddo

    do i = 1, msg_patt%n_send
        call isend_vdf_msg( send(i), no_co )
    enddo

    ! reshape vdf

#ifdef __DEBUG__
    ! (* debug *) - This sets the initial array values to NaN so the code can check if any value
    !               was left unset
    call reshape_nocopy( this, new_grid, dbg = .true. )
#else
    ! (* production *)
    call reshape_nocopy( this, new_grid )
#endif
    
    ! wait for receive messages to complete
    if ( msg_patt%n_recv > 0) then
        call wait_vdf_msg_array( recv )
    endif
    
    ! copy data to reshaped vdf
    do i = 1, msg_patt%n_recv
        j = 1
        call unpack_vdf_msg( recv(i), this, j )
    enddo
    
    ! clear message pattern data
    call msg_patt%cleanup()
    
    ! clear recv buffers
    if ( msg_patt%n_recv > 0) then
        call cleanup_vdf_msg_array( recv )
        
        call freemem( recv )
    endif
    
    ! wait for send messages to complete
    if ( msg_patt%n_send > 0) then
        call wait_vdf_msg_array( send )
        
        ! clear send buffers
        call cleanup_vdf_msg_array(send)
        
        call freemem( send )
    endif
    
    call update_boundary_vdf( this, p_vdf_replace, no_co, send_msg, recv_msg )

#ifdef __DEBUG__
    !(*debug*)  
    call this % check_nan( "Checking VDF grid after reshape_vdf_copy" )
#endif
    
end subroutine reshape_vdf_copy
!---------------------------------------------------------------------------------------------------



end module m_vdf_comm
                                                                                    