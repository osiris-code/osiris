! if __TEMPLATE__ is defined then read the template definition at the end of the file
#ifndef __TEMPLATE__

#include "os-config.h"
#include "os-preprocess.fpp"

module m_vdf_reportfile

#include "memory/memory.h"

use m_vdf_define
use m_space
use m_system

private

interface init_diag_file_vdf
  module procedure init_diag_file_vdf
end interface

interface write_vdf_1d_r4
  module procedure write_vdf_1d_r4
end interface

interface write_vdf_2d_r4
  module procedure write_vdf_2d_r4
end interface

interface write_vdf_3d_r4
  module procedure write_vdf_3d_r4
end interface

interface write_vdf_1d_r8
  module procedure write_vdf_1d_r8
end interface

interface write_vdf_2d_r8
  module procedure write_vdf_2d_r8
end interface

interface write_vdf_3d_r8
  module procedure write_vdf_3d_r8
end interface

public :: init_diag_file_vdf
public :: write_vdf_1d_r4, write_vdf_2d_r4, write_vdf_3d_r4
public :: write_vdf_1d_r8, write_vdf_2d_r8, write_vdf_3d_r8

contains

!---------------------------------------------------------------------------------------------------
! Initializes t_diag_file file object with supplied information
!---------------------------------------------------------------------------------------------------
subroutine init_diag_file_vdf( this, report, g_space, grid, diagFile )

  use m_space
  use m_grid_define
  use m_diagnostic_utilities

  implicit none

  class( t_vdf ),        intent(in) :: this
  type( t_vdf_report ), intent(in) :: report
  type( t_space ),      intent(in) :: g_space
  class( t_grid ),      intent(in) :: grid
  class( t_diag_file ),  allocatable :: diagFile

  integer :: i

  ! Initialize grid output file with default parameters
  call create_diag_file( diagFile )

  diagFile % ftype   = p_diag_grid

  ! Fill in parameters from specific report options
  diagFile%name           = report%name

  ! Iteration info
  diagFile%iter%n         = report%n
  diagFile%iter%t         = report%t
  diagFile%iter%time_units = report%time_units
  diagFile%grid%offset_t  = report%offset_t

  ! Grid info
  diagFile%grid%ndims = this % x_dim_
  diagFile%grid%name  = report % name
  diagFile%grid%label = report % label
  diagFile%grid%units = report % units

  do i = 1, this % x_dim_
    diagFile%grid%count(i) = grid % g_nx(i)
    diagFile%grid%axis(i)%type  = diag_axis_linear
    diagFile%grid%axis(i)%min  = xmin( g_space, i )
    diagFile%grid%axis(i)%max  = xmax( g_space, i )
    diagFile%grid%axis(i)%name  = report%xname(i)
    diagFile%grid%axis(i)%label = report%xlabel(i)
    diagFile%grid%axis(i)%units = report%xunits(i)
    diagFile%grid%offset_x(i) = report%offset_x(i)
  enddo

  ! Set file name and path
  diagFile%filename = trim(report%filename)
  diagFile%filepath = trim(report%path)

end subroutine init_diag_file_vdf
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!  Generate specific template functions for single and double precision datatypes.
!  Note that the module interfaces are not generated automatically and must be explicity written
!  in the module header above.
!
!  - write_vdf_1d_r4
!  - write_vdf_2d_r4
!  - write_vdf_3d_r4
!  - write_vdf_1d_r8
!  - write_vdf_2d_r8
!  - write_vdf_3d_r8
!---------------------------------------------------------------------------------------------------


#define __TEMPLATE__

! single precision real
#define __RKIND__ p_single
#define __MPI_TYPE__ MPI_REAL
#define FNAME(a) a ## _r4
#include __FILE__

! double precision real
#define __RKIND__ p_double
#define __MPI_TYPE__ MPI_DOUBLE_PRECISION
#define FNAME(a) a ## _r8
#include __FILE__

end module m_vdf_reportfile


!---------------------------------------------------------------------------------------------------
! Save vdf grid to file
!---------------------------------------------------------------------------------------------------
subroutine write_vdf( vdf, report, fc, g_space, grid  )

  use m_vdf_define,     only : t_vdf, t_vdf_report
  use m_space
  use m_grid_define
  use m_node_conf
  use m_system
  use m_parameters
  use m_vdf_reportfile

  implicit none

  class( t_vdf ),       intent(in) :: vdf
  type( t_vdf_report ), intent(in) :: report
  integer,              intent(in) :: fc
  type( t_space ),      intent(in) :: g_space
  class( t_grid ),      intent(in) :: grid

  ! Precision to be used for diagnostics
  integer :: prec

  ! Don't do diagnostics with better precision than the one used for storing the values
  if ( report%prec > p_k_fld ) then
    prec = p_k_fld
  else
    prec = report%prec
  endif

  ! Call appropriate routine according to precision / dimensions
  select case( prec )
    case ( p_single )
      select case (vdf % x_dim_)
      case(1)
        call write_vdf_1d_r4( vdf, report, fc, g_space, grid  )
      case(2)
        call write_vdf_2d_r4( vdf, report, fc, g_space, grid  )
      case(3)
        call write_vdf_3d_r4( vdf, report, fc, g_space, grid  )
      case default
        ERROR('Invalid value for vdf % x_dim_')
        call abort_program( p_err_invalid )
      end select

    case ( p_double )
      select case (vdf % x_dim_)
      case(1)
        call write_vdf_1d_r8( vdf, report, fc, g_space, grid  )
      case(2)
        call write_vdf_2d_r8( vdf, report, fc, g_space, grid  )
      case(3)
        call write_vdf_3d_r8( vdf, report, fc, g_space, grid  )
      case default
        ERROR('Invalid value for vdf % x_dim_')
        call abort_program( p_err_invalid )
      end select

    case default
      ERROR('Invalid precision selected for VDF diagnostics')
      call abort_program( p_err_invalid )
  end select

end subroutine write_vdf


#else

!---------------------------------------------------------------------------------------------------
!  Template Function definitions
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Merge / write a 1D dataset
!---------------------------------------------------------------------------------------------------
subroutine FNAME(write_vdf_1d)( vdf, report, fc, g_space, grid  )

  use m_space
  use m_grid_define
  use m_node_conf
  use m_diagnostic_utilities
  use m_parameters
  use zdf
  use m_diagfile

  implicit none

  class( t_vdf ),       intent(in) :: vdf
  type( t_vdf_report ), intent(in) :: report
  integer,              intent(in) :: fc
  type( t_space ),      intent(in) :: g_space
  class( t_grid ),      intent(in) :: grid

  ! Buffer for merging data - only group leaders will allocate this
  real( __RKIND__ ), dimension(:), pointer :: write_buffer

  ! Buffer for MPI communications - all nodes allocate this
  real( __RKIND__ ), dimension(:), pointer :: comm_buffer

  real( __RKIND__ ), dimension(:), pointer :: gather_buffer

  ! MPI datatype for receiving data on group root nodes
  integer :: recv_type

  ! Sizes and offsets for MPI datatype
  integer, dimension(1) :: sizes, starts

  ! Rank of node sending data
  integer :: rank

  ! Status variable for messages
  integer, dimension( MPI_STATUS_SIZE ) :: stat

  ! Ping variables
  integer :: ping, ping_handle

  ! Total data size being received
  integer :: data_size

  ! Sizes and displacements of data from each node
  integer, dimension(:), pointer :: recvcounts, displs

  ! Diagnostic file object
  class(t_diag_file), allocatable :: diagFile

  integer :: i, i1, k, ierr

  type( t_diag_dataset )   :: parDset
  type( t_diag_chunk ) :: parChunk

  ! Dummy receive buffer for MPI_GATHERV
  integer, dimension(1) :: tmp

  ! Merge data if necessary
  select case ( grid%io%merge_type )
  case (p_none)

    ! No merging, just copy local data to output buffer
    call alloc( write_buffer, [ vdf % nx_(1) ] )

    do i1 = 1, vdf % nx_(1)
      write_buffer( i1 ) = real( vdf % f1 ( fc, i1 ), __RKIND__ )
    enddo

  case (p_point2point)

    ! Group leaders receive data from other group members and store it in merge_buffer
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc( write_buffer, grid%io%group_tile( 1, 1:1 ) )

      ! Copy local data to merge buffer
      starts(1) = grid%io%tiles( 2, 1, 1 )
      do i1 = 1, vdf % nx_(1)
        write_buffer( starts(1)+i1 ) = real( vdf % f1 ( fc, i1 ), __RKIND__ )
      enddo

      ! Get data from other node
      do rank = 1, grid%io%group_size-1

        ! get size of grid tile on source node
        sizes(1) = grid%io%tiles( 1, 1, rank + 1)

        ! get start position of grid tile on merged data grid
        starts(1) = grid%io%tiles( 2, 1, rank + 1)

        ! Create dataype describing the tile to be received
        call mpi_type_create_subarray( 1, grid%io%group_tile( 1, : ), sizes, starts, &
              MPI_ORDER_FORTRAN, __MPI_TYPE__, recv_type, ierr )
        call mpi_type_commit( recv_type, ierr )

        ! notify node that we are ready
        call mpi_isend( ping, 1, MPI_INTEGER, rank, 0, &
                        grid%io%group_comm, ping_handle, ierr )

        ! Receive tile data
        call mpi_recv( write_buffer, 1, recv_type, rank, 1, &
                       grid%io%group_comm, stat, ierr )

        ! Free the datatype
        call mpi_type_free( recv_type, ierr )

        ! Wait for ping message to complete (this is just for cleanup, the ping must
        ! have completed by now)
        call mpi_wait( ping_handle, stat, ierr )

      enddo

    else

      ! Allocate buffer for messaging
      call alloc( comm_buffer, [ vdf % nx_(1) ] )

      ! Pack data for sending
      do i1 = 1, vdf % nx_(1)
        comm_buffer( i1 ) = real( vdf % f1 ( fc, i1 ), __RKIND__ )
      enddo

      ! receive ping from group leader
      call mpi_recv( ping, 1, MPI_INTEGER, 0, 0, &
                     grid%io%group_comm, stat, ierr )

      ! send data to group leader
      call mpi_send( comm_buffer, size(comm_buffer), __MPI_TYPE__, 0, &
                     1, grid%io%group_comm, ierr )

      ! free message buffer
      call freemem( comm_buffer )
    endif

  case (p_gather)

    ! Find counts and displacements for gathering the data
    if ( grid%io%group_rank == 0 ) then
      call alloc( recvcounts, [ grid % io % group_size ] )
      call alloc( displs, [ grid % io % group_size ] )

      data_size = 0

      do i = 1, grid%io%group_size
        recvcounts(i) = grid%io%tiles( 1, 1, i )
        data_size = data_size + recvcounts(i)
      enddo

      displs(1) = 0
      do i = 2, grid%io%group_size
        displs(i) = displs(i-1) + recvcounts(i-1)
      enddo

      call alloc( gather_buffer, [ data_size ] )
    else
      recvcounts   => null()
      displs       => null()
      gather_buffer => null()
    endif

    ! Copy the data to communication buffer
    call alloc( comm_buffer, [ vdf % nx_(1) ] )
    do i1 = 1, vdf % nx_(1)
      comm_buffer( i1 ) = real( vdf % f1 ( fc, i1 ), __RKIND__ )
    enddo

    ! Gather the data on the group root
    if ( grid%io%group_rank == 0 ) then
      call mpi_gatherv( comm_buffer, size( comm_buffer ), __MPI_TYPE__, &
                        gather_buffer, recvcounts, displs, __MPI_TYPE__, &
                        0, grid%io%group_comm, ierr )
    else
      call mpi_gatherv( comm_buffer, size( comm_buffer ), __MPI_TYPE__, &
                        tmp, tmp, tmp, __MPI_TYPE__, &
                        0, grid%io%group_comm, ierr )
    endif

    ! free the communication buffer
    call freemem( comm_buffer )

    ! Rearrange data into properly formed 2D arrays
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc( write_buffer, grid%io%group_tile( 1, : ) )

      ! Transpose data from merge_buffer into data_buffer
      k = 0
      do i = 1, grid%io%group_size
        do i1 = grid%io%tiles( 2, 1, i ) + 1, &
                grid%io%tiles( 2, 1, i ) + grid%io%tiles( 1, 1, i )
          k = k+1
          write_buffer( i1 ) = gather_buffer( k )
        enddo
      enddo

      ! Free allocated memory on group root
      call freemem( gather_buffer )
      call freemem( displs )
      call freemem( recvcounts )

    endif

  end select

  ! Write merged data to disk

  ! Only group leaders participate in writing data to disk
  ! Note than when not using merging all nodes have group_rank = 0
  if ( grid%io%group_rank == 0 ) then

    ! create file
    call init_diag_file_vdf( vdf, report, g_space, grid, diagFile )

    ! Open the file for (possible) parallel I/O
    call diagFile % open( p_diag_create, grid%io%comm )

    call diagFile % start_cdset( report%name, 1, grid % g_nx, freal_to_diagtype( __RKIND__ ), &
                    parDset, grid % io % chunk_size )

    parChunk % count(1) = grid % io % group_tile(  1 , 1 )
    parChunk % start(1) = grid % io % group_tile(  2 , 1 )
    parChunk % stride(1) = 1
    parChunk % data = c_loc( write_buffer )

    ! Write the data
    call diagFile % write_par_cdset( parDset, parChunk )

    ! Close dataset
    call diagFile % end_cdset( parDset )

    ! Close file
    call diagFile % close( )

    deallocate( diagFile )

    ! free merge data buffer
    call freemem( write_buffer )

  endif

end subroutine FNAME(write_vdf_1d)


!---------------------------------------------------------------------------------------------------
! Merge / write a 2D dataset
!---------------------------------------------------------------------------------------------------
subroutine FNAME(write_vdf_2d)( vdf, report, fc, g_space, grid  )

  use m_space
  use m_grid_define
  use m_node_conf
  use m_diagnostic_utilities
  use m_parameters
  use m_diagfile

  implicit none

  class( t_vdf ),       intent(in) :: vdf     ! vdf object
  type( t_vdf_report ), intent(in) :: report
  integer,              intent(in) :: fc       ! field component to report
  type( t_space ),      intent(in) :: g_space    ! spatial information
  class( t_grid ),      intent(in) :: grid


  ! Buffer for merging data - only group leaders will allocate this
  real( __RKIND__ ), dimension(:,:), pointer :: write_buffer => null()

  ! Buffer for MPI communications - all nodes allocate this
  real( __RKIND__ ), dimension(:,:), pointer :: comm_buffer

  real( __RKIND__ ), dimension(:), pointer :: gather_buffer

  ! MPI datatype for receiving data on group root nodes
  integer :: recv_type

  ! Sizes and offsets for MPI datatype
  integer, dimension(2) :: sizes, starts

  ! Rank of node sending data
  integer :: rank

  ! Status variable for messages
  integer, dimension( MPI_STATUS_SIZE ) :: stat

  ! Ping variables
  integer :: ping, ping_handle

  ! Total data size being received
  integer :: data_size

  ! Sizes and displacements of data from each node
  integer, dimension(:), pointer :: recvcounts, displs

  ! Diagnostic file object
  class(t_diag_file), allocatable :: diagFile

  integer :: i, i1, i2, k, ierr

  type( t_diag_dataset )   :: parDset
  type( t_diag_chunk ) :: parChunk

  ! Dummy receive buffer for MPI_GATHERV
  integer, dimension(1) :: tmp

  ! Merge data if necessary
  select case ( grid%io%merge_type )
  case (p_none)

    ! No merging, just copy local data to output buffer
    call alloc( write_buffer, vdf % nx_(1:2))

    do i2 = 1, vdf % nx_(2)
      do i1 = 1, vdf % nx_(1)
        write_buffer( i1, i2 ) = real( vdf % f2 ( fc, i1, i2 ), __RKIND__ )
      enddo
    enddo

  case (p_point2point)

    ! Group leaders receive data from other group members and store it in merge_buffer
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc( write_buffer, grid%io%group_tile( 1, 1:2 ) )

      ! Copy local data to merge buffer
      starts(1) = grid%io%tiles( 2, 1, 1 )
      starts(2) = grid%io%tiles( 2, 2, 1 )
      do i2 = 1, vdf % nx_(2)
        do i1 = 1, vdf % nx_(1)
          write_buffer( starts(1)+i1, starts(2)+i2 ) = real( vdf % f2 ( fc, i1, i2 ), __RKIND__ )
        enddo
      enddo

      ! Get data from other node
      do rank = 1, grid%io%group_size-1

        ! get size of grid tile on source node
        sizes(1) = grid%io%tiles( 1, 1, rank + 1)
        sizes(2) = grid%io%tiles( 1, 2, rank + 1)

        ! get start position of grid tile on merged data grid
        starts(1) = grid%io%tiles( 2, 1, rank + 1)
        starts(2) = grid%io%tiles( 2, 2, rank + 1)

        ! Create dataype describing the tile to be received
        call mpi_type_create_subarray( 2, grid%io%group_tile( 1, 1:2 ), sizes, starts, &
              MPI_ORDER_FORTRAN, __MPI_TYPE__, recv_type, ierr )
        call mpi_type_commit( recv_type, ierr )

        ! notify node that we are ready
        call mpi_isend( ping, 1, MPI_INTEGER, rank, 0, &
                        grid%io%group_comm, ping_handle, ierr )

        ! Receive tile data
        call mpi_recv( write_buffer, 1, recv_type, rank, 1, &
                       grid%io%group_comm, stat, ierr )

        ! Free the datatype
        call mpi_type_free( recv_type, ierr )

        ! Wait for ping message to complete (this is just for cleanup, the ping must
        ! have completed by now)
        call mpi_wait( ping_handle, stat, ierr )

      enddo

    else

      ! Allocate buffer for messaging
      call alloc( comm_buffer, vdf % nx_(1:2))

      ! Pack data for sending
      do i2 = 1, vdf % nx_(2)
        do i1 = 1, vdf % nx_(1)
          comm_buffer( i1, i2 ) = real( vdf % f2 ( fc, i1, i2 ), __RKIND__ )
        enddo
      enddo

      ! receive ping from group leader
      call mpi_recv( ping, 1, MPI_INTEGER, 0, 0, &
                     grid%io%group_comm, stat, ierr )

      ! send data to group leader
      call mpi_send( comm_buffer, size(comm_buffer), __MPI_TYPE__, 0, &
                     1, grid%io%group_comm, ierr )

      ! free message buffer
      call freemem( comm_buffer )
    endif

  case (p_gather)

    ! Find counts and displacements for gathering the data
    if ( grid%io%group_rank == 0 ) then
      call alloc( recvcounts, [ grid % io % group_size ] )
      call alloc( displs, [ grid % io % group_size ] )

      data_size = 0

      do i = 1, grid%io%group_size
        recvcounts(i) = grid%io%tiles( 1, 1, i ) * &
                        grid%io%tiles( 1, 2, i )
        data_size = data_size + recvcounts(i)
      enddo

      displs(1) = 0
      do i = 2, grid%io%group_size
        displs(i) = displs(i-1) + recvcounts(i-1)
      enddo

      call alloc( gather_buffer, [ data_size ] )
    else
      recvcounts   => null()
      displs       => null()
      gather_buffer => null()
    endif

    ! Copy the data to communication buffer
    call alloc( comm_buffer, vdf % nx_(1:2) )
    do i2 = 1, vdf % nx_(2)
      do i1 = 1, vdf % nx_(1)
        comm_buffer( i1, i2 ) = real( vdf % f2 ( fc, i1, i2 ), __RKIND__ )
      enddo
    enddo

    ! Gather the data on the group root
    if ( grid%io%group_rank == 0 ) then
      call mpi_gatherv( comm_buffer, size( comm_buffer ), __MPI_TYPE__, &
                        gather_buffer, recvcounts, displs, __MPI_TYPE__, &
                        0, grid%io%group_comm, ierr )
    else
      call mpi_gatherv( comm_buffer, size( comm_buffer ), __MPI_TYPE__, &
                        tmp, tmp, tmp, __MPI_TYPE__, &
                        0, grid%io%group_comm, ierr )
    endif

    ! free the communication buffer
    call freemem( comm_buffer )

    ! Rearrange data into properly formed 2D arrays
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc( write_buffer, grid%io%group_tile( 1, : ) )

      ! Transpose data from merge_buffer into data_buffer
      k = 0
      do i = 1, grid%io%group_size
        do i2 = grid%io%tiles( 2, 2, i ) + 1, &
                grid%io%tiles( 2, 2, i ) + grid%io%tiles( 1, 2, i )
          do i1 = grid%io%tiles( 2, 1, i ) + 1, &
                  grid%io%tiles( 2, 1, i ) + grid%io%tiles( 1, 1, i )
            k = k+1
            write_buffer( i1, i2 ) = gather_buffer( k )
          enddo
        enddo
      enddo

      ! Free allocated memory on group root
      call freemem( gather_buffer )
      call freemem( displs )
      call freemem( recvcounts )

    endif

  end select

  ! Write merged data to disk

  ! Only group leaders participate in writing data to disk
  ! Note than when not using merging all nodes have group_rank = 0
  if ( grid%io%group_rank == 0 ) then

    ! create file
    call init_diag_file_vdf( vdf, report, g_space, grid, diagFile )

    ! Open the file for (possible) parallel I/O
    call diagFile % open( p_diag_create, grid%io%comm )

    call diagFile % start_cdset( report%name, 2, grid % g_nx, freal_to_diagtype( __RKIND__ ), &
                    parDset, grid % io % chunk_size )

    do i = 1, 2
      parChunk % count(i) = grid % io % group_tile(  1 , i )
      parChunk % start(i) = grid % io % group_tile(  2 , i )
      parChunk % stride(i) = 1
    enddo
    parChunk % data = c_loc( write_buffer )

    ! Write the data
    call diagFile % write_par_cdset( parDset, parChunk )

    ! Close dataset
    call diagFile % end_cdset( parDset )

    ! Close file
    call diagFile % close( )

    deallocate( diagFile )

    ! free merge data buffer
    call freemem( write_buffer )

  endif

end subroutine FNAME(write_vdf_2d)


!---------------------------------------------------------------------------------------------------
! Merge / write a 3D dataset
!---------------------------------------------------------------------------------------------------
subroutine FNAME(write_vdf_3d)( vdf, report, fc, g_space, grid  )

  use m_space
  use m_grid_define
  use m_node_conf
  use m_diagnostic_utilities
  use m_parameters
  use m_diagfile

  implicit none

  class( t_vdf ),       intent(in) :: vdf     ! vdf object
  type( t_vdf_report ), intent(in) :: report
  integer,              intent(in) :: fc       ! field component to report
  type( t_space ),      intent(in) :: g_space    ! spatial information
  class( t_grid ),      intent(in) :: grid


  ! Buffer for merging data - only group leaders will allocate this
  real( __RKIND__ ), dimension(:,:,:), pointer :: write_buffer

  ! Buffer for MPI communications - all nodes allocate this
  real( __RKIND__ ), dimension(:,:,:), pointer :: comm_buffer

  real( __RKIND__ ), dimension(:), pointer :: gather_buffer

  ! MPI datatype for receiving data on group root nodes
  integer :: recv_type

  ! Sizes and offsets for MPI datatype
  integer, dimension(3) :: sizes, starts

  ! Rank of node sending data
  integer :: rank

  ! Status variable for messages
  integer, dimension( MPI_STATUS_SIZE ) :: stat

  ! Ping variables
  integer :: ping, ping_handle

  ! Total data size being received
  integer :: data_size

  ! Sizes and displacements of data from each node
  integer, dimension(:), pointer :: recvcounts, displs

  ! Diagnostic file object
  class(t_diag_file), allocatable :: diagFile

  integer :: i, i1, i2, i3, k, ierr

  type( t_diag_dataset )   :: parDset
  type( t_diag_chunk ) :: parChunk

  ! Dummy receive buffer for MPI_GATHERV
  integer, dimension(1) :: tmp

  ! Merge data if necessary
  select case ( grid%io%merge_type )
  case (p_none)

    ! No merging, just copy local data to output buffer
    call alloc( write_buffer, vdf % nx_(1:3))

    do i3 = 1, vdf % nx_(3)
      do i2 = 1, vdf % nx_(2)
        do i1 = 1, vdf % nx_(1)
          write_buffer( i1, i2, i3 ) = real( vdf % f3 ( fc, i1, i2, i3 ), __RKIND__ )
        enddo
      enddo
    enddo

  case (p_point2point)

    ! Group leaders receive data from other group members and store it in merge_buffer
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc( write_buffer, grid%io%group_tile( 1, 1:3 ) )

      ! Copy local data to merge buffer
      starts(1) = grid%io%tiles( 2, 1, 1 )
      starts(2) = grid%io%tiles( 2, 2, 1 )
      starts(3) = grid%io%tiles( 2, 3, 1 )

      do i3 = 1, vdf % nx_(3)
        do i2 = 1, vdf % nx_(2)
          do i1 = 1, vdf % nx_(1)
            write_buffer( starts(1)+i1, starts(2)+i2, starts(3)+i3 ) = &
              real( vdf % f3 ( fc, i1, i2, i3 ), __RKIND__ )
          enddo
        enddo
      enddo

      ! Get data from other node
      do rank = 1, grid%io%group_size-1

        ! get size of grid tile on source node
        sizes(1) = grid%io%tiles( 1, 1, rank + 1)
        sizes(2) = grid%io%tiles( 1, 2, rank + 1)
        sizes(3) = grid%io%tiles( 1, 3, rank + 1)

        ! get start position of grid tile on merged data grid
        starts(1) = grid%io%tiles( 2, 1, rank + 1)
        starts(2) = grid%io%tiles( 2, 2, rank + 1)
        starts(3) = grid%io%tiles( 2, 3, rank + 1)

        ! Create dataype describing the tile to be received
        call mpi_type_create_subarray( 3, grid%io%group_tile( 1, 1:3 ), sizes, starts, &
              MPI_ORDER_FORTRAN, __MPI_TYPE__, recv_type, ierr )
        call mpi_type_commit( recv_type, ierr )

        ! notify node that we are ready
        call mpi_isend( ping, 1, MPI_INTEGER, rank, 0, &
                        grid%io%group_comm, ping_handle, ierr )

        ! Receive tile data
        call mpi_recv( write_buffer, 1, recv_type, rank, 1, &
                       grid%io%group_comm, stat, ierr )

        ! Free the datatype
        call mpi_type_free( recv_type, ierr )

        ! Wait for ping message to complete (this is just for cleanup, the ping must
        ! have completed by now)
        call mpi_wait( ping_handle, stat, ierr )

      enddo

    else

      ! Allocate buffer for messaging
      call alloc( comm_buffer, vdf % nx_(:))

      ! Pack data for sending
      do i3 = 1, vdf % nx_(3)
        do i2 = 1, vdf % nx_(2)
          do i1 = 1, vdf % nx_(1)
            comm_buffer( i1, i2, i3 ) = real( vdf % f3 ( fc, i1, i2, i3 ), __RKIND__ )
          enddo
        enddo
      enddo

      ! receive ping from group leader
      call mpi_recv( ping, 1, MPI_INTEGER, 0, 0, &
                     grid%io%group_comm, stat, ierr )

      ! send data to group leader
      call mpi_send( comm_buffer, size(comm_buffer), __MPI_TYPE__, 0, &
                     1, grid%io%group_comm, ierr )

      ! free message buffer
      call freemem( comm_buffer )
    endif

  case (p_gather)

    ! Find counts and displacements for gathering the data
    if ( grid%io%group_rank == 0 ) then
      call alloc( recvcounts, [ grid % io % group_size ] )
      call alloc( displs,  [ grid % io % group_size ] )

      data_size = 0

      do i = 1, grid%io%group_size
        recvcounts(i) = grid%io%tiles( 1, 1, i ) * &
                        grid%io%tiles( 1, 2, i ) * &
                        grid%io%tiles( 1, 3, i )
        data_size = data_size + recvcounts(i)
      enddo

      displs(1) = 0
      do i = 2, grid%io%group_size
        displs(i) = displs(i-1) + recvcounts(i-1)
      enddo

      call alloc( gather_buffer, [ data_size ] )
    else
      recvcounts   => null()
      displs       => null()
      gather_buffer => null()
    endif

    ! Copy the data to communication buffer
    call alloc( comm_buffer, vdf % nx_(1:3) )
    do i3 = 1, vdf % nx_(3)
      do i2 = 1, vdf % nx_(2)
        do i1 = 1, vdf % nx_(1)
          comm_buffer( i1, i2, i3 ) = real( vdf % f3 ( fc, i1, i2, i3 ), __RKIND__ )
        enddo
      enddo
    enddo

    ! Gather the data on the group root
    if ( grid%io%group_rank == 0 ) then
      call mpi_gatherv( comm_buffer, size( comm_buffer ), __MPI_TYPE__, &
                        gather_buffer, recvcounts, displs, __MPI_TYPE__, &
                        0, grid%io%group_comm, ierr )
    else
      call mpi_gatherv( comm_buffer, size( comm_buffer ), __MPI_TYPE__, &
                        tmp, tmp, tmp, __MPI_TYPE__, &
                        0, grid%io%group_comm, ierr )
    endif

    ! free the communication buffer
    call freemem( comm_buffer )

    ! Rearrange data into properly formed 2D arrays
    if ( grid%io%group_rank == 0 ) then

      ! Allocate buffer for merging data
      call alloc( write_buffer, grid%io%group_tile( 1, : ) )

      ! Transpose data from merge_buffer into data_buffer
      k = 0
      do i = 1, grid%io%group_size
        do i3 = grid%io%tiles( 2, 3, i ) + 1, &
                grid%io%tiles( 2, 3, i ) + grid%io%tiles( 1, 3, i )
          do i2 = grid%io%tiles( 2, 2, i ) + 1, &
                  grid%io%tiles( 2, 2, i ) + grid%io%tiles( 1, 2, i )
            do i1 = grid%io%tiles( 2, 1, i ) + 1, &
                    grid%io%tiles( 2, 1, i ) + grid%io%tiles( 1, 1, i )
              k = k+1
              write_buffer( i1, i2, i3 ) = gather_buffer( k )
            enddo
          enddo
        enddo
      enddo

      ! Free allocated memory on group root
      call freemem( gather_buffer )
      call freemem( displs )
      call freemem( recvcounts )

    endif

  end select

!  SCR_MPINODE('Writing merged data...')

  ! Write merged data to disk

  ! Only group leaders participate in writing data to disk
  ! Note than when not using merging all nodes have group_rank = 0
  if ( grid%io%group_rank == 0 ) then

    ! create file
    call init_diag_file_vdf( vdf, report, g_space, grid, diagFile )

!  SCR_MPINODE('Opening file...')

    ! Open the file for (possible) parallel I/O
    call diagFile % open( p_diag_create, grid % io % comm )

!  SCR_MPINODE('Starting CDSET...')

    call diagFile % start_cdset( report%name, 3, grid % g_nx, freal_to_diagtype( __RKIND__ ), &
                    parDset, grid % io % chunk_size )

    do i = 1, 3
      parChunk % count(i) = grid % io % group_tile(  1 , i )
      parChunk % start(i) = grid % io % group_tile(  2 , i )
      parChunk % stride(i) = 1
    enddo
    parChunk % data = c_loc( write_buffer )

!  SCR_MPINODE('Writing CDSET...')

    ! Write the data
    call diagFile % write_par_cdset( parDset, parChunk )

    ! Close dataset
    call diagFile % end_cdset( parDset )

    ! Close file
    call diagFile % close( )

    deallocate( diagFile )

    ! free merge data buffer
    call freemem( write_buffer )

  endif

end subroutine FNAME(write_vdf_3d)

!---------------------------------------------------------------------------------------------------
!  End of template Functions, do not change below
!---------------------------------------------------------------------------------------------------

#undef __RKIND__
#undef __MPI_TYPE__
#undef FNAME

#endif
