#include "os-config.h"
#include "os-preprocess.fpp"


module m_emf_psi

#include "memory/memory.h"

use m_parameters

use m_emf_define

use m_vdf
use m_vdf_define
use m_vdf_comm
use m_node_conf

use m_logprof

implicit none

private

interface get_psi
  module procedure get_psi
end interface

interface cleanup
  module procedure cleanup_psi
end interface

public :: get_psi, cleanup

contains

!---------------------------------------------------------------------------------------------------
! Initialize psi object
!---------------------------------------------------------------------------------------------------
subroutine init_psi( emf, no_co )

  implicit none

  class( t_emf ), intent(inout) :: emf
  class( t_node_conf ), intent(in) :: no_co

  integer :: nodes_x, ngp_x
  integer :: color, key, ierr

  ! Initialize psi vdf
  call emf%psi%data%new( emf%e, f_dim = 1 )

  ! Initialize communicator if necessary
  nodes_x = nx( no_co, 1 )

  if ( nodes_x > 1 ) then
    ngp_x = my_ngp( no_co, 1 )

    ! Communicator for MPI_SCAN use

    ! Reverse the order of the individual ranks so that the scan function adds the values
    ! right to left
    key = nodes_x - ngp_x

    ! split communicators along x2 and x3 directions so that nodes in the same communicator
    ! share the same x2 and x3 coordinates
    select case( emf%psi%data%x_dim_ )
      case(1)
        color = 1
      case(2)
        color = my_ngp( no_co, 2 )
      case(3)
        color = ( my_ngp( no_co, 2 ) - 1 ) + nx( no_co, 2 ) * ( my_ngp( no_co, 3 ) - 1 )
    end select

    call MPI_COMM_SPLIT( comm( no_co ), color, key, emf%psi%comm, ierr )

  else

    emf%psi%comm = mpi_comm_null

  endif

end subroutine init_psi
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Cleanup psi object
!---------------------------------------------------------------------------------------------------
subroutine cleanup_psi( psi )

  implicit none

  type( t_emf_psi ), intent(inout) :: psi

  integer :: ierr

  call psi%data%cleanup()
  psi%n = -1

  if ( psi%comm /= mpi_comm_null ) then
    call MPI_COMM_FREE( psi%comm, ierr )
    psi%comm = mpi_comm_null
  endif

end subroutine cleanup_psi
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine get_psi( emf, n, no_co, psi, send_msg, recv_msg )

  implicit none

  class( t_emf ), intent(inout), target :: emf
  integer, intent(in) :: n
  class( t_node_conf ), intent(in) :: no_co

  type( t_vdf ), pointer :: psi
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  ! Initialize the Psi object if required
  if ( emf%psi%n < 0 ) then
    call init_psi( emf, no_co )
  endif

  ! If data for this iteration is not available calculate it
  if ( n /= emf%psi%n ) then
    call update_psi( emf, no_co, send_msg, recv_msg )
    emf%psi%n = n
  endif

  ! return pointer to psi data
  psi => emf%psi%data

end subroutine get_psi
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Calculate the Psi diagnostic in parallel, MPI_SCAN version
!---------------------------------------------------------------------------------------------------
subroutine update_psi( emf, no_co, send_msg, recv_msg )

  implicit none

  class( t_emf ), intent(inout) :: emf
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg


  integer :: i1, i2, i3, nx1
  real( p_k_fld ) :: dx

  real( p_k_fld ) :: acc1D, acc_out1D
  real( p_k_fld ), dimension(:), pointer   :: acc2D, acc_out2D
  real( p_k_fld ), dimension(:,:), pointer :: acc3D, acc_out3D

  integer :: nodes_x, ngp_x
  integer :: mpi_type, count, ierr

  !---

  ! sychronize nodes (this is mostly for benchmarking purposes)
  call no_co % barrier( )

  call begin_event(getpsi_ev)

  ! calculate local psi
  dx = real( emf%e%dx_(1), p_k_fld )
  nx1 = emf%psi%data%nx_( 1 )
  select case ( emf%psi%data%x_dim_ )

    case (1)
      emf%psi%data%f1( 1, nx1 ) = emf%e%f1( 1, nx1 ) * dx
      do i1 = nx1-1, 1, -1
        emf%psi%data%f1( 1, i1 ) = emf%psi%data%f1( 1, i1+1 ) + emf%e%f1( 1, i1 ) * dx
      enddo

    case (2)
      do i2 = 1, emf%psi%data%nx_( 2 )
        emf%psi%data%f2( 1, nx1, i2 ) = emf%e%f2( 1, nx1, i2 ) * dx
        do i1 = nx1-1, 1, -1
          emf%psi%data%f2( 1, i1, i2 ) = emf%psi%data%f2( 1, i1+1, i2 ) + emf%e%f2( 1, i1, i2 ) * dx
        enddo
      enddo

    case (3)
      do i3 = 1, emf%psi%data%nx_( 3 )
        do i2 = 1, emf%psi%data%nx_( 2 )
          emf%psi%data%f3( 1, nx1, i2, i3 ) = emf%e%f3( 1, nx1, i2, i3 ) * dx
          do i1 = nx1-1, 1, -1
            emf%psi%data%f3( 1, i1, i2, i3 ) = emf%psi%data%f3( 1, i1+1, i2, i3 ) + emf%e%f3( 1, i1, i2, i3 ) * dx
          enddo
        enddo
      enddo

  end select

  ! Add contribution from other nodes to the right
  if ( no_num( no_co ) > 1 ) then

    nodes_x = nx( no_co, 1 )

    if ( nodes_x > 1 ) then

       ngp_x = my_ngp( no_co, 1 )

       select case ( p_k_fld )
         case( p_single )
           mpi_type = MPI_REAL
         case( p_double )
           mpi_type = MPI_DOUBLE_PRECISION
       end select

       ! Version using MPI_SCAN
       ! SCR_ROOT('Calculating PSI using MPI_SCAN')

       select case ( emf%psi%data%x_dim_ )
         case(1)

           acc1D = emf%psi%data%f1( 1, 1 )

           ! add contributions from all nodes in communicator
           ! (the openmpi-1.2.9 implementation of MPI_EXSCAN + MPI_IN_PLACE is broken )
           count = 1
           call MPI_EXSCAN( acc1D, acc_out1D, count, mpi_type, MPI_SUM, emf%psi%comm, ierr )

           ! add contribution to local data
           if ( ngp_x < nodes_x ) then
             do i1 = 1, emf%psi%data%nx_( 1 )
               emf%psi%data%f1( 1, i1 ) = emf%psi%data%f1( 1, i1 ) + acc_out1D
             enddo
           endif

         case(2)
           ! create temp buffer with leftmost values
           call alloc( acc2D, (/ emf%psi%data%nx_( 2 ) /))
           call alloc( acc_out2D, (/ emf%psi%data%nx_( 2 ) /) )

           do i2 = 1, emf%psi%data%nx_( 2 )
             acc2D( i2 ) = emf%psi%data%f2( 1, 1, i2 )
           enddo

           ! add contributions from all nodes in communicator
           ! (the openmpi-1.2.9 implementation of MPI_EXSCAN + MPI_IN_PLACE is broken )
           count = emf%psi%data%nx_( 2 )
           call MPI_EXSCAN( acc2D, acc_out2D, count, mpi_type, MPI_SUM, emf%psi%comm, ierr )

           ! add contribution to local data
           if ( ngp_x < nodes_x ) then
             do i2 = 1, emf%psi%data%nx_( 2 )
               do i1 = 1, emf%psi%data%nx_( 1 )
                 emf%psi%data%f2( 1, i1, i2 ) = emf%psi%data%f2( 1, i1, i2 ) + acc_out2D( i2 )
               enddo
             enddo
           endif

           ! delete temp buffer with leftmost values
           call freemem( acc2D )
           call freemem( acc_out2D )

         case(3)

           ! create temp buffer with leftmost values
           call alloc( acc3D, (/ emf%psi%data%nx_( 2 ), emf%psi%data%nx_( 3 ) /) )
           call alloc( acc_out3D, (/ emf%psi%data%nx_( 2 ), emf%psi%data%nx_( 3 ) /) )

           do i3 = 1, emf%psi%data%nx_( 3 )
             do i2 = 1, emf%psi%data%nx_( 2 )
               acc3D( i2, i3 ) = emf%psi%data%f3( 1, 1, i2, i3 )
             enddo
           enddo

           ! add contributions from all nodes in communicator
           ! (the openmpi-1.2.9 implementation of MPI_EXSCAN + MPI_IN_PLACE is broken )
           count = emf%psi%data%nx_( 2 ) * emf%psi%data%nx_( 3 )
           call MPI_EXSCAN( acc3D, acc_out3D, count, mpi_type, MPI_SUM, emf%psi%comm, ierr )

           ! add contribution to local data
           if ( ngp_x < nodes_x ) then
             do i3 = 1, emf%psi%data%nx_( 3 )
               do i2 = 1, emf%psi%data%nx_( 2 )
                 do i1 = 1, emf%psi%data%nx_( 1 )
                   emf%psi%data%f3( 1, i1, i2, i3 ) = emf%psi%data%f3( 1, i1, i2, i3 ) + acc_out3D( i2, i3 )
                 enddo
               enddo
             enddo
           endif

           ! delete temp buffer with leftmost values
           call freemem( acc3D )
           call freemem( acc_out3D )

       end select

    endif

    ! fill guard cells with values from another node
    ! (this should be made optional)
    call update_boundary( emf%psi%data, p_vdf_replace, no_co, send_msg, recv_msg )
  endif

  call end_event(getpsi_ev)


end subroutine update_psi
!---------------------------------------------------------------------------------------------------



end module m_emf_psi
