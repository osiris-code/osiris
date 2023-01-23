!-----------------------------------------------------------------------------------------
! Species tracks module
!
! This file contains the routines for handling the t_tracks objects
!
!-----------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

#define TRACK_FILE_FORMAT2 1
! #define DEBUG_IO 1

#ifdef DEBUG_IO
#define DEBUGMSG(...) print *, __VA_ARGS__
#define CHECK_IOERROR( ierr ) if (ierr<0) write(0,*) "(*io error*) ", __FILE__, __LINE__-1
#else
#define DEBUGMSG(...)
#define CHECK_IOERROR( ierr )
#endif

module m_species_tracks

#include "memory/memory.h"

use m_system
use m_parameters

use stringutil
use m_utilities
use m_node_conf
use m_species_define

use m_emf_define
use m_emf_psi
use m_vdf_interpolate

use m_vdf_define

implicit none

private

integer, parameter :: p_max_tagstring_len = 32

! string to id restart data
character(len=*), parameter :: p_track_rst_id = "track rst data - 0x0003"

interface missing_particles
  module procedure missing_particles_list
end interface

interface inbound_particle
  module procedure inbound_particle
end interface

interface new_particles
  module procedure new_particles_bounds
  module procedure new_particles_single
end interface

interface add_track_data
  module procedure add_track_set_data
end interface

interface add_single_track
  module procedure add_single_track_data
end interface

interface write_tracks
  module procedure write_tracks
end interface

interface create_file
  module procedure create_track_set_file
end interface

interface setup
  module procedure setup_single_track
  module procedure setup_track_set
end interface

interface cleanup
  module procedure cleanup_single_track
  module procedure cleanup_track_set
end interface

interface update_indexes
  module procedure update_indexes_single_track
  module procedure update_indexes_track_set
end interface

interface change_index
  module procedure change_index_track_set
end interface

interface restart_write
  module procedure restart_write_tracks
  module procedure restart_write_single_track
end interface

interface restart_read
  module procedure restart_read_tracks
  module procedure restart_read_single_track
end interface

interface cleanup_buffers_tracks
  module procedure cleanup_buffers_tracks
end interface

public :: setup, restart_write, cleanup, cleanup_buffers_tracks
public :: create_file, write_tracks, add_track_data, add_single_track
public :: update_indexes, change_index
public :: missing_particles, new_particles, inbound_particle

interface alloc
  module procedure alloc_track
  module procedure alloc_1d_track
  module procedure alloc_bound_1d_track
end interface

interface freemem
  module procedure freemem_track
  module procedure freemem_1d_track
end interface

contains

!-----------------------------------------------------------------------------------------
! rearrange indexes after main species buffer has been sorted
!-----------------------------------------------------------------------------------------
subroutine update_indexes_single_track( track, new_idx )

   implicit none

   type( t_track ), intent(inout) :: track
   integer, dimension(:), intent(in) :: new_idx

   track%part_idx = new_idx(track%part_idx)

end subroutine update_indexes_single_track

!-----------------------------------------------------------------------------------------
! sort track according to iteration (needed after gathering from all nodes)
!-----------------------------------------------------------------------------------------
subroutine sort_track( track )

   implicit none

   type( t_track ), intent(inout) :: track

   integer, dimension(:), pointer :: idx => null(), temp_int => null()
   real(p_k_part), dimension(:), pointer :: temp_double => null()
   integer :: i, j, nquants

   ! if just 1 point there's no need to sort
   if ( track%npoints > 1 ) then
      call alloc( idx, (/track%npoints/))

      ! get sorted index
      call heapsort( track%n(1:track%npoints), idx )

      ! reorder iterations
      call alloc( temp_int, (/track%npoints/))

      do i = 1, track%npoints
        temp_int(i) = track%n(idx(i))
      enddo
      do i = 1, track%npoints
        track%n(i) = temp_int(i)
      enddo
      call freemem( temp_int )

      ! reorder remaining data
      call alloc( temp_double, (/track%npoints/))

      nquants = size( track % data, 1 )

      do j = 1, nquants
         do i = 1, track%npoints
           temp_double(i) = track%data(j,idx(i))
         enddo
         do i = 1, track%npoints
           track%data(j,i) = temp_double(i)
         enddo
      enddo
      call freemem( temp_double )

      call freemem( idx )
   endif

end subroutine sort_track

!-----------------------------------------------------------------------------------------
! Pack single track data into comm buffer
!-----------------------------------------------------------------------------------------
subroutine pack_data_single_track( track, no_co, buffer, bufsize, position )

   !use mpi

   implicit none

   type( t_track ), intent(inout) :: track
   class( t_node_conf ), intent(in) :: no_co

   integer( p_byte ), dimension(:), intent(inout) :: buffer
   integer, intent(in) :: bufsize
   integer, intent(inout) :: position

   integer :: mpi_type
   integer :: i, nquants, ierr

   ! pack number of points
   call mpi_pack( track%npoints, 1, MPI_INTEGER, &
                  buffer, bufsize, position, &
                  comm(no_co), ierr )

   if ( track%npoints > 0 ) then
      ! add iterations
      call mpi_pack( track%n(1:track%npoints), track%npoints, MPI_INTEGER, &
                     buffer, bufsize, position, &
                     comm(no_co), ierr )

      mpi_type = mpi_real_type( p_k_part )

      ! add datasets
      nquants = size( track % data, 1 )
      do i = 1, nquants
         call mpi_pack( track%data(i, 1:track%npoints), track%npoints, mpi_type, &
                        buffer, bufsize, position, &
                        comm(no_co), ierr )
      enddo

      ! free data points
      track%npoints = 0

   endif

end subroutine pack_data_single_track

!-----------------------------------------------------------------------------------------
! Unpack single track data from comm buffer
!-----------------------------------------------------------------------------------------
subroutine unpack_data_single_track( track, buffer, bufsize, position )

   !use mpi

   implicit none

   type( t_track ), intent(inout) :: track

   integer( p_byte ), intent(inout), dimension(:) :: buffer
   integer, intent(in) :: bufsize
   integer, intent(inout) :: position

   integer :: nquants, npoints, i, ierr
   integer :: mpi_type
   real( p_k_part ), dimension(:), pointer :: quant_buffer => null()
   integer, dimension(:), pointer :: n_buffer => null()

   ! unpack number of points
   call mpi_unpack( buffer, bufsize, position, &
                    npoints, 1, MPI_INTEGER, &

                    mpi_comm_world, ierr )

   if ( npoints > 0 ) then
      call alloc( n_buffer, (/ npoints /) )
      ! unpack iterations
      call mpi_unpack( buffer, bufsize, position, &
                       n_buffer, npoints, MPI_INTEGER, &
                       mpi_comm_world, ierr )
      track%n(track%npoints+1: track%npoints+npoints) = n_buffer
      call freemem( n_buffer )

      ! add datasets
      mpi_type = mpi_real_type( p_k_part )

      call alloc( quant_buffer, (/ npoints /) )
      nquants = size( track % data, 1 )
      do i = 1, nquants
         call mpi_unpack( buffer, bufsize, position, &
                          quant_buffer, npoints, mpi_type, &
                          mpi_comm_world, ierr )
         track%data(i, track%npoints+1: track%npoints+npoints) = quant_buffer
      enddo
      call freemem( quant_buffer )

      track%npoints = track%npoints + npoints
   endif

end subroutine unpack_data_single_track

!-----------------------------------------------------------------------------------------
! Add current iteration point to track
!-----------------------------------------------------------------------------------------
subroutine add_single_track_data( track, species, n, t, fld )

!-----------------------------------------------------------------------------------------

   implicit none

   type( t_track ), intent(inout) :: track
   class( t_species ), intent(in) :: species
   integer, intent(in) :: n
   real(p_double), intent(in) :: t

   real(p_k_part), dimension(:), intent(in), optional :: fld

   ! local variables
   real(p_k_part), dimension(:), pointer :: pos
   real(p_k_part) :: u2

  integer :: i, nquants, pp, n_x_dim

   n_x_dim = species%get_n_x_dims()
   call alloc( pos ,(/n_x_dim/))
   track%npoints = track%npoints + 1

   track%n( track%npoints ) = n

   track%data(1, track%npoints ) = real( t, p_k_part )
   track%data(2, track%npoints ) = species%q( track%part_idx )

   ! energy must be calculated here
   u2 = species%p( 1, track%part_idx )**2 + &
        species%p( 2, track%part_idx )**2 + &
        species%p( 3, track%part_idx )**2

  track%data(3, track%npoints ) = u2 / ( sqrt(1 + u2) + 1 )

   pp = 3
   call species % get_position( track%part_idx, pos )
   track%data( pp+1:pp+n_x_dim , track%npoints ) = pos
   pp = pp + n_x_dim

   track%data( pp+1:pp+p_p_dim, track%npoints ) = species%p( 1:p_p_dim, track%part_idx )
   pp = pp + p_p_dim

#ifdef __HAS_SPIN__

   track%data( pp+1:pp+p_s_dim, track%npoints ) = species%s( 1:p_s_dim, track%part_idx )
   pp = pp + p_s_dim

#endif

   ! Add extra field data, if any
   if (present(fld)) then
    nquants = size( track % data, 1 )
    do i = 1, nquants - pp
       track%data( pp+i, track%npoints ) = fld(i)
     enddo
   endif
   call freemem(pos)

end subroutine add_single_track_data
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleanup single track object
!-----------------------------------------------------------------------------------------
subroutine cleanup_single_track( track )

!-----------------------------------------------------------------------------------------
   implicit none

   type( t_track ), intent(inout) :: track

   track%savedpoints = 0
   track%npoints = 0
   track%tag = 0
   call freemem( track%data )
   call freemem( track%n )

end subroutine cleanup_single_track
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Setup single track object
!-----------------------------------------------------------------------------------------
subroutine setup_single_track( track, max_points, nquants, tag )

!-----------------------------------------------------------------------------------------
   implicit none

   type( t_track ), intent(inout) :: track
  integer, intent(in) :: max_points, nquants
   integer, dimension(:) :: tag


   call alloc( track%n, (/max_points/))
  call alloc( track%data, (/nquants, max_points/))

   track%savedpoints = 0
   track%npoints = 0
   track%tag(1) = tag(1)
   track%tag(2) = tag(2)
   track%part_idx = -1

end subroutine setup_single_track
!-----------------------------------------------------------------------------------------

!*****************************************************************************************
!************************************ Track set routines *********************************
!*****************************************************************************************

!-----------------------------------------------------------------------------------------
! Change particle pointer if particle has moved inside particle buffer (usually
! because another one was deleted, sort is handled below)
!-----------------------------------------------------------------------------------------
subroutine change_index_track_set( track_set, old_idx, new_idx )

   implicit none

   type( t_track_set ), intent(inout) :: track_set

   integer, intent(in) :: old_idx, new_idx

   integer :: j

   do j = 1, track_set%npresent
     if ( track_set%tracks(track_set%present(j))%part_idx == old_idx ) then
       track_set%tracks(track_set%present(j))%part_idx = new_idx
       exit
     endif
   enddo

end subroutine change_index_track_set

!-----------------------------------------------------------------------------------------
! rearrange indexes after main species buffer has been sorted
!-----------------------------------------------------------------------------------------
subroutine update_indexes_track_set( track_set, new_idx )

   implicit none

   type( t_track_set ), intent(inout) :: track_set
   integer, dimension(:), intent(in) :: new_idx

   integer :: i

   do i = 1, track_set%npresent
      call update_indexes( track_set%tracks( track_set%present(i) ), new_idx )
   enddo

end subroutine update_indexes_track_set

!-----------------------------------------------------------------------------------------
! Particles have left the local node, update present/missing lists
!-----------------------------------------------------------------------------------------
subroutine missing_particles_list( track_set, idx, n_idx, tags )

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  integer, dimension(:) :: idx
  integer, dimension(:,:), intent(in) :: tags

  integer, intent(in) :: n_idx

  integer :: i, j

  ! loop through all the missing particles
  do i = 1, n_idx

    ! if no tags are present we can end the search
    if (track_set%npresent == 0) exit

    ! if particle was being tracked its tag(1) will be < 0
    if ( tags(1,idx(i)) < 0 ) then

     ! loop through all missing tags

     j = 1
     do
        ! The particle must be on the missing list, if not there is an error
        ASSERT( j <= track_set%npresent )

        if (idx(i) == track_set%tracks(track_set%present(j))%part_idx) then

           ! delete particle index
           ! (not necessary for production but useful for debugging)
           track_set%tracks(track_set%present(j))%part_idx = -1

           ! add track to missing list
           track_set%nmissing = track_set%nmissing + 1
           track_set%missing( track_set%nmissing ) = track_set%present(j)

           ! remove track from present list
           track_set%present(j) = track_set%present(track_set%npresent)

           track_set%npresent = track_set%npresent - 1

           exit
        else

           j = j + 1
        endif
     enddo

    endif

  enddo

end subroutine missing_particles_list

!-----------------------------------------------------------------------------------------
! New particle is present in local nodes, check if its tag is in the
! missing list and if so store their index and move it to the present list
!-----------------------------------------------------------------------------------------
subroutine new_particles_single( track_set, tag, idx )

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  integer, dimension(2), intent(inout) :: tag
  integer, intent(in) :: idx

  integer :: j

  if ( track_set%nmissing > 0 ) then
    j = 1
    do
      if (j > track_set%nmissing) exit
      if ( track_set%tracks(track_set%missing(j))%tag(1) == tag(1) .and. &
         track_set%tracks(track_set%missing(j))%tag(2) == tag(2) ) then

         ! store particle index
         track_set%tracks( track_set%missing(j))%part_idx = idx

         ! add tag to present list
         track_set%npresent = track_set%npresent + 1
         track_set%present( track_set%npresent ) = track_set%missing(j)

         ! remove tag from missing list
         track_set%missing(j) = track_set%missing(track_set%nmissing)

         track_set%nmissing = track_set%nmissing - 1

         ! mark the particle as being tracked by flipping the sign of tag(1)
         tag(1) = - tag(1)

         exit
      else
         j = j + 1
      endif
    enddo
  endif

end subroutine new_particles_single

!-----------------------------------------------------------------------------------------
! New particles are creates in local nodes, check if their tags are in the
! missing list and if so store their index and move them to the present list
!-----------------------------------------------------------------------------------------
subroutine new_particles_bounds( track_set, spec, idx0, idx1 )

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  class( t_species ), intent(inout) :: spec
  integer, intent(in) :: idx0, idx1

  integer, dimension(2) :: tag
  integer :: i, j

  ! loop through all the new particles
  do i = idx0, idx1

     ! loop through all missing tags

     j = 1
     do
        if (j > track_set%nmissing) exit

        tag(1) = track_set%tracks(track_set%missing(j))%tag(1)
        tag(2) = track_set%tracks(track_set%missing(j))%tag(2)

        if (( spec%tag(1, i) == tag(1)) .and. ( spec%tag(2, i) == tag(2))) then

           ! store particle index
           track_set%tracks( track_set%missing(j))%part_idx = i

           ! add tag to present list
           track_set%npresent = track_set%npresent + 1
           track_set%present( track_set%npresent ) = track_set%missing(j)

           ! remove tag from missing list
           track_set%missing(j) = track_set%missing(track_set%nmissing)

           track_set%nmissing = track_set%nmissing - 1

           ! mark the particle as being tracked by flipping the sign of tag(1)
           spec%tag(1, i) = -spec%tag(1, i)

           exit
        else

           j = j + 1
        endif
     enddo

     ! if no tags are missing we can end the search
     if (track_set%nmissing == 0) exit
  enddo

end subroutine new_particles_bounds

!-----------------------------------------------------------------------------------------
! Check if particles coming from another node are in the missing list. Note that particles
! being tracked will have a negative tag(1), and that this function will only be called
! for those particles.
!-----------------------------------------------------------------------------------------

subroutine inbound_particle( track_set, tag, idx )

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  integer, dimension(2), intent(in) :: tag
  integer, intent(in) :: idx

  integer :: j

  ! The tag(1) of the inbound particle must be < 0

  ASSERT( tag(1) < 0 )

  ! At least 1 particle must be missing from this node
  ASSERT( track_set%nmissing > 0 )

  j = 1
  do

    ! The particle must be on the missing list, if not there is an error
    ASSERT( j <= track_set%nmissing )

    if ( track_set%tracks(track_set%missing(j))%tag(1) == - tag(1) .and. &
         track_set%tracks(track_set%missing(j))%tag(2) ==   tag(2) ) then

       ! debug
       !write(0,*) mpi_node(), ' particle has entered node node, tag = ', &
       !                       tag

       ! store particle index
       track_set%tracks( track_set%missing(j))%part_idx = idx

       ! add tag to present list
       track_set%npresent = track_set%npresent + 1
       track_set%present( track_set%npresent ) = track_set%missing(j)

       ! remove tag from missing list
       track_set%missing(j) = track_set%missing(track_set%nmissing)

       track_set%nmissing = track_set%nmissing - 1
       exit
    else
       j = j + 1
    endif
  enddo

end subroutine inbound_particle

!-----------------------------------------------------------------------------------------
! Add current values to track data
!-----------------------------------------------------------------------------------------
subroutine add_track_set_data( track_set, spec, no_co, emf, n, t, send_msg, recv_msg )

  use m_vdf_interpolate
  use m_vdf_comm, only : t_vdf_msg

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  class( t_species ), intent(in) :: spec
  class( t_node_conf ), intent(in) :: no_co
  class( t_emf ),  intent(inout) :: emf
  integer, intent(in) :: n
  real( p_double ), intent(in) :: t
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  ! local variables
  integer :: i, k, idx, n_x_dim

  type( t_vdf ), pointer :: psi

  psi => null()

  n_x_dim = spec % get_n_x_dims()
  ! This needs to be done even if no particles are present in this node
  if ( track_set%ifdmp_tracks_psi ) then
    call get_psi( emf, n, no_co, psi, send_msg, recv_msg )
  endif

  if ( track_set%npresent > 0 ) then

    if ( track_set%nfields > 0 ) then

      if ( track_set%npresent > track_set%size_buffer ) then
        if ( track_set%size_buffer > 0 ) then
          call freemem( track_set%field_buffer )
          call freemem( track_set%pos_buffer )
          call freemem( track_set%ipos_buffer )
        endif

        ! grow buffer in sizes of 1024 tracks
        track_set%size_buffer = ceiling( track_set%npresent / 1024.0 ) * 1024

        call alloc( track_set%pos_buffer, (/ n_x_dim, track_set%size_buffer /) )
        call alloc( track_set%ipos_buffer, (/ p_x_dim, track_set%size_buffer /) )

        ! If interpolating any fields also grow field_buffer
        if ( track_set%nfields > 0 ) &
          call alloc( track_set%field_buffer, (/ track_set%size_buffer, track_set%nfields /) )

      endif

      ! Get positions of particles being tracked
      do i = 1, track_set%npresent

        idx = track_set%tracks( track_set%present(i) )%part_idx
        track_set%ipos_buffer(1:p_x_dim, i) = spec%ix( 1:p_x_dim, idx )
        track_set%pos_buffer(1:n_x_dim, i)  = spec%x( 1:n_x_dim, idx )
      enddo

      k = 1
      do i = 1, 3
        if ( track_set%ifdmp_tracks_efl(i) ) then
          call interpolate( emf%e, i, p_emf_hc_e(i), track_set%ipos_buffer, track_set%pos_buffer, track_set%npresent, &
                            spec%interpolation, track_set%field_buffer( :, k ) )
          k = k + 1
        endif
      enddo

      do i = 1, 3
        if ( track_set%ifdmp_tracks_bfl(i) ) then
          call interpolate( emf%b, i, p_emf_hc_b(i), track_set%ipos_buffer, track_set%pos_buffer, track_set%npresent, &
                            spec%interpolation, track_set%field_buffer( :, k ) )
          k = k + 1
        endif
      enddo

      if ( track_set%ifdmp_tracks_psi ) then
        call interpolate( psi, 1, 0, track_set%ipos_buffer, track_set%pos_buffer, track_set%npresent, &
                          spec%interpolation, track_set%field_buffer( :, k ) )
        k = k + 1
      endif

      ! Add track data including interpolated E and B fields, and Psi diagnostic
      do i = 1, track_set%npresent
        call add_single_track_data( track_set%tracks( track_set%present(i) ), &
                                    spec, n, t, track_set%field_buffer( i, : ) )
      enddo

    else

      ! Add track data
      do i = 1, track_set%npresent
        call add_single_track_data( track_set%tracks( track_set%present(i) ), &
                                    spec, n, t )
      enddo

    endif

  endif

end subroutine add_track_set_data

!-----------------------------------------------------------------------------------------
! Cleanup temporary buffers
!-----------------------------------------------------------------------------------------
subroutine cleanup_buffers_tracks()

  implicit none


  ! not needed now
  ! if ( size_buffer > 0 ) then
  !   call freemem( field_buffer )
  !   call freemem( pos_buffer )
  !   call freemem( ipos_buffer )

  !   size_buffer = 0
  ! endif

end subroutine cleanup_buffers_tracks
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Gather single track data from all nodes and sort it
!-----------------------------------------------------------------------------------------
subroutine gather_data( track_set, no_co, spec )

   !use mpi

   implicit none

   type( t_track_set ), intent(inout) :: track_set
   class( t_node_conf ), intent(in) :: no_co
   class( t_species ), intent(in) :: spec

   integer( p_byte ), dimension(:), pointer :: buffer => null()
   integer :: bufsize, position, totalpoints, nquants
   integer :: pack_size_int
   integer :: pack_size_double

   integer :: mpi_type
   integer, dimension(mpi_status_size):: stat
   integer :: handle

   integer :: i, j, ierr

   nquants = 3 + spec%get_n_x_dims() + p_p_dim + track_set%nfields
#ifdef __HAS_SPIN__
   nquants = nquants + p_s_dim
#endif

   if ( no_num( no_co ) > 1 ) then

      ! get pack sizes
      call mpi_pack_size( 1, MPI_INTEGER, comm(no_co), pack_size_int, ierr)
      mpi_type = mpi_real_type( p_k_part )
      call mpi_pack_size( 1, mpi_type, comm(no_co), pack_size_double, ierr)

      if ( root( no_co ) ) then
         ! allocate comm bufer to the max size
         ! npoints, iteration data, remaining data
         bufsize = track_set%ntracks * ( pack_size_int * (1 + track_set%maxpoints ) + &
                        pack_size_double * track_set%maxpoints * nquants )

         call alloc( buffer, (/ bufsize /))

         do j = 2, no_num( no_co )
            ! print *, ' 1 - sending ping to ', j
            call send_ping( no_co, j-1, 0 )

            ! print *, ' 1 - posted receive from ', j
            call MPI_RECV( buffer, bufsize, MPI_PACKED, j-1, 0, &
                           comm( no_co ), stat, ierr )

            position = 0
            ! print *, ' 1 - unpacking data from ', j
            do i = 1, track_set%ntracks
              call unpack_data_single_track( track_set%tracks(i),buffer, bufsize, position )
            enddo
            ! print *, ' 1 - finished with data from node ',j
         enddo

         ! sort the data
         do i = 1, track_set%ntracks
            call sort_track( track_set%tracks(i) )
         enddo
      else

         totalpoints = track_set%tracks(1)%npoints
         do i = 2, track_set%ntracks
           totalpoints = totalpoints + track_set%tracks(i)%npoints
         enddo

         ! allocate comm bufer to the size of points in this node
         bufsize = track_set%ntracks * pack_size_int + &
                   totalpoints*(pack_size_int + nquants * pack_size_double )

         call alloc( buffer, (/ bufsize /))

         position = 0
         do i = 1, track_set%ntracks
           call pack_data_single_track( track_set%tracks(i), no_co, &
                                        buffer, bufsize, position )
         enddo

         ! wait for ping
         ! print *, my_aid( no_co ) , ' - waiting for ping from node 1'
         call recv_ping( no_co, 0, 0 )

         ! send data
         !print *, no_co%my_aid() , ' - sending data to node 1'
         call MPI_ISEND( buffer, bufsize, MPI_PACKED, 0, 0, &
                         comm( no_co ), handle, ierr )
         call MPI_WAIT( handle, stat, ierr )

      endif

      call freemem( buffer )

   endif

end subroutine gather_data

!-----------------------------------------------------------------------------------------
! setup track set object
!-----------------------------------------------------------------------------------------
subroutine setup_track_set( track_set, spec_name, n_x_dim, ndump, restart, restart_handle )

   use m_restart

   implicit none

   type( t_track_set ), intent(inout) :: track_set
   character(len = *), intent(in) :: spec_name
   integer, intent(in) :: n_x_dim
   integer, intent(in) :: ndump ! number of iterations between each write
   logical, intent(in) :: restart
   type( t_restart_handle ), intent(in) :: restart_handle

   character(len = 256) :: linebuffer
   integer :: i, nquants, tag1,tag2, ierr


   if ( restart ) then

      call restart_read( track_set, n_x_dim, restart_handle )

      if ( track_set%ntracks > 0 ) then
        ! Check if the user did not change ndump_fac_tracks in the input file
        if ( track_set%maxpoints /= ndump/track_set%niter + 1 ) then
          write(0,*) '(*error*) ndump_fac and/or ndump_fac_tracks set in the input file are not compatible with restart data, aborting.'
          call abort_program()
        endif
      endif

   else if ( ndump > 0 ) then

      ! read tags
      open( file_id_tem, file = track_set%file_tags, status = 'old', &
            access = 'sequential', form = 'formatted', iostat = ierr )
      if ( ierr /= 0 ) then
         write(0,*) 'Unable to open particle tags file "', trim(track_set%file_tags),'"'
         call abort_program()
      endif

      ! skip comments until first valid value
      do
         read( unit = file_id_tem, fmt= '(A)') linebuffer
         if ( linebuffer(1:1) /= '!' ) then
            backspace file_id_tem
            exit
         endif
      enddo

      read( unit = file_id_tem, fmt= '(I12)') track_set%ntracks

      if ( track_set%ntracks < 1 ) then
         write(0,*) 'Error reading particle tags file "',trim(track_set%file_tags),'"'
         write(0,*) 'Invalid number of tracks: ',track_set%ntracks
         call abort_program()
      endif

      ! allocate track objects
      call alloc( track_set%tracks, (/track_set%ntracks/) )
      call alloc( track_set%present, (/track_set%ntracks/) )
      call alloc( track_set%missing, (/ track_set%ntracks /) )

      ! set numeber of present/missing tracks
      track_set%npresent = 0
      track_set%nmissing = track_set%ntracks

      ! skip comments until first valid value
      do
         read( unit = file_id_tem, fmt= '(A)') linebuffer
         if ( linebuffer(1:1) /= '!' ) then
            backspace file_id_tem
            exit
         endif
      enddo

      track_set%maxpoints = ndump/track_set%niter + 1

      ! read tags and setup individual track objects
      nquants = 3 + n_x_dim + p_p_dim + track_set%nfields
#ifdef __HAS_SPIN__
      nquants = nquants + p_s_dim
#endif

      do i = 1, track_set%ntracks
         read( unit = file_id_tem, fmt= '(I12,I12)') tag1, tag2
         call setup( track_set%tracks(i), track_set%maxpoints, nquants, (/tag1,tag2/) )

         ! store all tracks in the missing list
         track_set%missing(i) = i
      enddo

      ! close tags file
      close( file_id_tem )

      ! set output file name and path
      track_set % file_path = trim(path_mass) // 'TRACKS' // p_dir_sep
      track_set % file_name = replace_blanks(trim(spec_name)) // '-tracks'

   endif

end subroutine setup_track_set

!-----------------------------------------------------------------------------------------
! destroy track set object
!-----------------------------------------------------------------------------------------
subroutine cleanup_track_set( track_set )

   implicit none

   type( t_track_set ), intent(inout) :: track_set

   integer :: i

   do i = 1, track_set%ntracks
     call cleanup( track_set%tracks(i) )
   enddo

   if ( track_set%ntracks > 0 ) then

      call freemem( track_set%tracks )
      call freemem( track_set%present )
      call freemem( track_set%missing )
   endif

   if ( track_set%size_buffer > 0 ) then
    call freemem( track_set%field_buffer )
    call freemem( track_set%pos_buffer )
    call freemem( track_set%ipos_buffer )

    track_set%size_buffer = 0
  endif

end subroutine

!-----------------------------------------------------------------------------------------
! write track restart information
!-----------------------------------------------------------------------------------------
subroutine restart_write_tracks( track_set, restart_handle )

   use m_restart

   implicit none

   type( t_track_set ), intent(in) :: track_set
   type( t_restart_handle ), intent(inout) :: restart_handle

   character(len=*), parameter :: err_msg = 'error writing restart data for track object.'
   integer :: i, ierr

   restart_io_wr( p_track_rst_id, restart_handle, ierr )
   CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

   restart_io_wr( track_set%ntracks, restart_handle, ierr )
   CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

   if ( track_set%ntracks > 0 ) then

      restart_io_wr( track_set%niter, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      restart_io_wr( track_set%ifdmp_tracks_efl, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      restart_io_wr( track_set%ifdmp_tracks_bfl, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      restart_io_wr( track_set%ifdmp_tracks_psi, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      restart_io_wr( track_set%nfields, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      restart_io_wr( track_set%maxpoints, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      restart_io_wr( track_set%npresent, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      restart_io_wr( track_set%present, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      restart_io_wr( track_set%nmissing, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      restart_io_wr( track_set%missing, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      restart_io_wr( track_set%file_path, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      restart_io_wr( track_set%file_name, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      do i = 1, track_set%ntracks
        call restart_write( track_set%tracks(i), restart_handle )
      enddo

   endif

end subroutine restart_write_tracks
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! read track restart information
!-----------------------------------------------------------------------------------------
subroutine restart_read_tracks( track_set, n_x_dim, restart_handle )

   use m_restart

   implicit none

   type( t_track_set ), intent(inout) :: track_set
   integer, intent(in) :: n_x_dim
   type( t_restart_handle ), intent(in) :: restart_handle

   character(len=*), parameter :: err_msg = 'error reading restart data for track object.'
   character(len=len(p_track_rst_id)) :: rst_id
   integer :: i, nquants, ierr

   restart_io_rd( rst_id, restart_handle, ierr )
   CHECK_ERROR( ierr, err_msg, p_err_rstrd )

   ! check if restart file is compatible
   if ( rst_id /= p_track_rst_id) then
     ERROR('Corrupted restart file, or restart file ')
     ERROR('from incompatible binary (tracks)')
     call abort_program(p_err_rstrd)
   endif

   restart_io_rd( track_set%ntracks, restart_handle, ierr )
   CHECK_ERROR( ierr, err_msg, p_err_rstrd )

   if ( track_set%ntracks > 0 ) then
      ! allocate track objects
      call alloc( track_set%tracks, (/track_set%ntracks/) )
      call alloc( track_set%present, (/ track_set%ntracks /) )
      call alloc( track_set%missing, (/ track_set%ntracks /) )

      ! read tracking lists
      restart_io_rd( track_set%niter, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstrd )

      restart_io_rd( track_set%ifdmp_tracks_efl, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstrd )

      restart_io_rd( track_set%ifdmp_tracks_bfl, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstrd )

      restart_io_rd( track_set%ifdmp_tracks_psi, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstrd )

      restart_io_rd( track_set%nfields, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstrd )

      restart_io_rd( track_set%maxpoints, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstrd )

      restart_io_rd( track_set%npresent, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstrd )

      restart_io_rd( track_set%present, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstrd )

      restart_io_rd( track_set%nmissing, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstrd )

      restart_io_rd( track_set%missing, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstrd )

      restart_io_rd( track_set%file_path, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstrd )

      restart_io_rd( track_set%file_name, restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstrd )

      ! read individual tracks
      nquants = 3 + n_x_dim + p_p_dim + track_set%nfields
#ifdef __HAS_SPIN__
      nquants = nquants + p_s_dim
#endif

      do i = 1, track_set%ntracks
        call restart_read( track_set%tracks(i), track_set%maxpoints, &
                           nquants, restart_handle )
      enddo

   endif

end subroutine restart_read_tracks
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! write track restart information
!-----------------------------------------------------------------------------------------
subroutine restart_write_single_track( track, restart_handle )

   use m_restart

   implicit none

   type( t_track ), intent(in) :: track
   type( t_restart_handle ), intent(inout) :: restart_handle

   character(len=*), parameter :: err_msg = 'error writing restart data for single_track object.'
   integer :: ierr

   restart_io_wr( track%npoints, restart_handle, ierr )
   CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

   restart_io_wr( track%savedpoints, restart_handle, ierr )
   CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

   restart_io_wr( track%tag, restart_handle, ierr )
   CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

   restart_io_wr( track%part_idx, restart_handle, ierr )
   CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

   if ( track%npoints > 0 ) then
      restart_io_wr( track%data( :, 1:track%npoints ), restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      restart_io_wr( track%n(1:track%npoints ), restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
   endif

end subroutine restart_write_single_track
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! read track restart information
!-----------------------------------------------------------------------------------------
subroutine restart_read_single_track( track, max_points, nquants, restart_handle )

   use m_restart

   implicit none

   type( t_track ), intent(inout) :: track
   integer, intent(in) :: max_points, nquants
   type( t_restart_handle ), intent(in) :: restart_handle

   character(len=*), parameter :: err_msg = 'error reading restart data for single_track object.'
   integer :: ierr

   ! allocate required memory
   call setup_single_track( track, max_points, nquants, (/0,0/) )

   ! read information from file
   restart_io_rd( track%npoints, restart_handle, ierr )
   CHECK_ERROR( ierr, err_msg, p_err_rstrd )

   restart_io_rd( track%savedpoints, restart_handle, ierr )
   CHECK_ERROR( ierr, err_msg, p_err_rstrd )

   restart_io_rd( track%tag, restart_handle, ierr )
   CHECK_ERROR( ierr, err_msg, p_err_rstrd )

   restart_io_rd( track%part_idx, restart_handle, ierr )
   CHECK_ERROR( ierr, err_msg, p_err_rstrd )

   if ( track%npoints > 0 ) then
      restart_io_rd( track%data( :, 1:track%npoints ), restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstrd )

      restart_io_rd( track%n( 1:track%npoints ), restart_handle, ierr )
      CHECK_ERROR( ierr, err_msg, p_err_rstrd )
   endif

end subroutine restart_read_single_track
!-----------------------------------------------------------------------------------------


!*****************************************************************************************
!*********************************** File output routines ********************************
!*****************************************************************************************


#ifdef TRACK_FILE_FORMAT2
!-----------------------------------------------------------------------------------------
! create tracks file - format 2
!-----------------------------------------------------------------------------------------
subroutine create_track_set_file( track_set, species, ndump )
  use m_diagnostic_utilities

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  class( t_species ), intent(in) :: species
  integer, intent(in) :: ndump

  class( t_diag_file ), allocatable   :: diagFile
  type( t_diag_dataset ) :: dset

  integer, parameter :: p_max_txt_len = 16

  if (track_set%ntracks > 0 ) then

    if ( mpi_node() == 0 ) then

      call create_diag_file( diagFile )
      diagFile % ftype   = p_diag_tracks

      diagFile % filepath = track_set % file_path
      diagFile % filename = track_set % file_name

      diagFile % name = species%name

      diagFile % tracks % name = species%name
      ! Species object doesn't have a label (yet)
      diagFile % tracks % label = species%name
      diagFile % tracks % ntracks = track_set%ntracks
      diagFile % tracks % ndump   = ndump
      diagFile % tracks % niter   = track_set%niter

      call species % enumerate_quants( diagFile, track_set)

      ! Open file
      call diagFile % open( p_diag_create )

      ! create iteration map
      call diagFile % start_cdset( "itermap", 2, [3,0], p_diag_int32, dset, &
                        chunk_size = [3, track_set % ntracks] )

      ! create main dataset
      ! The iteration is not stored explicitly in the file so the number of quantities is nquants-1
      call diagFile % start_cdset( "data", 2, [diagFile % tracks %nquants-1,0], diag_type_real( p_k_part ), dset, &
                        chunk_size = [diagFile % tracks %nquants-1, track_set % ntracks] )

      ! Close file
      call diagFile % close()

      deallocate( diagFile )


    endif

  endif

end subroutine create_track_set_file
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Write tracks data to file - format 2
!-----------------------------------------------------------------------------------------
subroutine write_tracks( track_set, no_co, spec )

  use m_diagnostic_utilities

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  class( t_node_conf ), intent(in) :: no_co
  class( t_species ), intent(in) :: spec
  class( t_diag_file ), allocatable :: diagFile
  type( t_diag_dataset ) :: dset
  type( t_diag_chunk ) :: chunk
  integer :: i, k, np, np_total, nquants, ntracks

  real(p_k_part), dimension(:,:), pointer :: track_data
  integer, dimension(:,:), pointer :: iter_map
  integer :: idx

  if (track_set%ntracks > 0 ) then

    ! gather data from all nodes
    call gather_data( track_set, no_co, spec )

    if ( root( no_co ) ) then

      ! open file
      call create_diag_file( diagFile )
      diagFile % filepath = track_set % file_path
      diagFile % filename = track_set % file_name

      call diagFile % open( p_diag_update )
      ! build iteration map
      call alloc( iter_map, (/ 3, track_set%ntracks /) )

      k = 0
      np_total = 0
      do i = 1, track_set%ntracks
        np = track_set%tracks(i)%npoints
        if ( np > 0 ) then

          k = k + 1
          np_total = np_total + np
          iter_map( 1, k ) = i                          ! track id
          iter_map( 2, k ) = np                         ! number of points
          iter_map( 3, k ) = track_set%tracks(i)%n(1)   ! start iter
        endif

      enddo

      ntracks = k

      if ( ntracks > 0 ) then

        ! append iteration map
        dset % name = "itermap"
        call diagFile % open_cdset( dset )

        chunk % count(1:2) = [3,ntracks]
        chunk % start(1:2) = [0_p_int64, dset % count(2)]
        chunk % stride(1:2) = [1,1]
        chunk % data = c_loc( iter_map )
        call diagFile % extend_cdset( dset, [ dset % count(1), dset % count(2) + ntracks ] )
        call diagFile % write_cdset( dset, chunk )
        call diagFile % close_cdset( dset )

        ! build track_data
        nquants = size( track_set % tracks(1) % data, 1 )
        call alloc( track_data, (/ nquants, np_total /) )

        k = 1
        do i = 1, ntracks
          idx = iter_map( 1, i )
          track_data( 1 : nquants, k : k + iter_map( 2, i ) - 1 ) = &
                    track_set % tracks( idx ) % data( 1 : nquants, 1 : iter_map( 2, i ))

          k = k + iter_map( 2, i )

          ! update saved points data
          track_set % tracks( idx ) % savedpoints = &
                    track_set % tracks( idx ) % savedpoints + &
                    track_set % tracks( idx ) % npoints

          track_set % tracks( idx ) % npoints = 0

        enddo

        ! append track data
        dset % name = "data"
        call diagFile % open_cdset( dset )
        chunk % count(1:2) = [ nquants, np_total ]
        chunk % start(1:2) = [ 0_p_int64, dset % count(2) ]
        chunk % stride(1:2) = [ 1, 1 ]
        chunk % data = c_loc( track_data )
        call diagFile % extend_cdset( dset, [ dset % count(1), dset % count(2) + np_total ] )
        call diagFile % write_cdset( dset, chunk )
        call diagFile % close_cdset( dset )

        ! free temp memory
        call freemem( track_data )

      endif

      ! free temp memory
      call freemem( iter_map )

      ! close file
      call diagFile % close()

      deallocate( diagFile )
    endif

  endif

end subroutine write_tracks
!-----------------------------------------------------------------------------------------
#else
!-----------------------------------------------------------------------------------------
! create tracks file - format 1
!-----------------------------------------------------------------------------------------
subroutine create_track_set_file( track_set, species, ndump )
  use m_diagnostic_utilities

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  class( t_species ), intent(in) :: species
  integer, intent(in) :: ndump

  integer :: i
  class( t_diag_file ), allocatable   :: diagFile
  type( t_diag_dataset ) :: dset
  integer, parameter :: p_max_txt_len = 16
  character( len = p_max_tagstring_len ) :: tagstring



  if (track_set%ntracks > 0 ) then
    if ( mpi_node() == 0 ) then

      call create_diag_file( diagFile )
      diagFile % ftype   = p_diag_tracks

      diagFile % filepath = track_set % file_path
      diagFile % filename = track_set % file_name

      diagFile % name = species%name

      diagFile % tracks % name = species%name
      ! Species object doesn't have a label (yet)
      diagFile % tracks % label = species%name
      diagFile % tracks % ntracks = track_set%ntracks
      diagFile % tracks % ndump   = ndump
      diagFile % tracks % niter   = track_set%niter

      call species % enumerate_quants( diagFile, track_set)

      ! Open file
      call diagFile % open( p_diag_create )

      do i = 1, track_set%ntracks

        tagstring = trim(tostring(track_set%tracks(i)%tag(1)))// &
                    '-'//trim(tostring(track_set%tracks(i)%tag(2)))
        ! create particle group
        call diagFile % start_cdset(tagstring, 2, [diagFile % tracks %nquants,0], diag_type_real( p_k_part ), dset, &
                        chunk_size = [diagFile % tracks %nquants, track_set% maxpoints] )
        call diagFile % close_cdset( dset )
      enddo

      ! Close file
      call diagFile % close()

      deallocate( diagFile )


    endif

  endif

end subroutine create_track_set_file
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Write tracks data to file - format 1
!-----------------------------------------------------------------------------------------
subroutine write_tracks( track_set, no_co, spec )

  use m_diagnostic_utilities

  implicit none

  type( t_track_set ), intent(inout) :: track_set
  class( t_node_conf ), intent(in) :: no_co
  class( t_species ), intent(in) :: spec

  class( t_diag_file ), allocatable :: diagFile
  type( t_diag_dataset ) :: dset
  type( t_diag_chunk ) :: chunk
  integer :: i, npoints, nquants
  integer, parameter :: p_max_txt_len = 16

  real(p_k_part), dimension(:,:), pointer :: track_data

  if (track_set%ntracks > 0 ) then

    ! gather data from all nodes
    call gather_data( track_set, no_co, spec )

    if ( root( no_co ) ) then

      ! open file
      call create_diag_file( diagFile )
      diagFile % filepath = track_set % file_path
      diagFile % filename = track_set % file_name

      call diagFile % open( p_diag_update )
      nquants = size( track_set % tracks(1) % data, 1 ) + 1
      !iterate over tracks
      do i = 1, track_set%ntracks

        npoints = track_set%tracks(i)%npoints
        if(npoints > 0) then

          call alloc( track_data, (/ nquants, npoints /) )

          ! open track dset
          dset % name = trim(tostring(track_set%tracks(i)%tag(1)))// &
                        '-'//trim(tostring(track_set%tracks(i)%tag(2)))

          call diagFile % open_cdset( dset )

          ! store data
          track_data(1,:) = track_set%tracks(i)%n(1:npoints)
          track_data(2:nquants,:) = track_set%tracks(i)%data(:,1:npoints)

          chunk % count(1:2) = [ nquants, npoints ]
          chunk % start(1:2) = [ 0_p_int64, dset % count(2) ]
          chunk % stride(1:2) = [ 1, 1 ]
          chunk % data = c_loc( track_data )

          call diagFile % extend_cdset( dset, [ dset % count(1), dset % count(2) + npoints ] )
          call diagFile % write_cdset( dset, chunk )
          call diagFile % close_cdset( dset )

          track_set % tracks(i) % savedpoints  = & 
                                                track_set % tracks(i) % savedpoints + &
                                                track_set % tracks(i) % npoints
          track_set % tracks(i) % npoints = 0

          call freemem(track_data)

        endif

      enddo

      ! close file
      call diagFile % close()

      deallocate( diagFile )
    endif

  endif

end subroutine write_tracks
!-----------------------------------------------------------------------------------------


#endif


!-----------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines

#define __TYPE__ type( t_track )
#define __TYPE_STR__ "t_track"
#define FNAME( a )  a ## _track
#define __MAX_DIM__ 1
#include "memory/mem-template.h"

!-----------------------------------------------------------------------------------------

end module m_species_tracks

!-----------------------------------------------------------------------------------------
! Write metadata for track files
!-----------------------------------------------------------------------------------------
subroutine enumerate_quants_spec( this, diagFile, track_set )

  use m_parameters
  use m_diagnostic_utilities, only : t_diag_file
  use m_species_define, only : t_species, t_track_set
#ifdef __HAS_SPIN__
  use m_parameters, only : p_s_dim
#endif

  implicit none

  class( t_species ), intent(in) :: this
  class( t_diag_file ),intent(inout) :: diagFile
  type( t_track_set ), intent(inout) :: track_set
  integer :: i, j

  ! enumerate quants
  j = 0

  j = j+1
  diagFile % tracks % quants(j) = 'n'
  diagFile % tracks % qlabels(j) = 'n'
  diagFile % tracks % qunits(j)  = ''
  diagFile % tracks % offset_t(j)= 0.0_p_double

  j = j+1
  diagFile % tracks % quants(j) = 't'
  diagFile % tracks % qlabels(j) = 't'
  diagFile % tracks % qunits(j)  = '1/\omega_p'
  diagFile % tracks % offset_t(j)= 0.0_p_double

  j = j+1
  diagFile % tracks % quants(j) = 'q'
  diagFile % tracks % qlabels(j) = 'q'
  diagFile % tracks % qunits(j)  = 'e'
  diagFile % tracks % offset_t(j)= 0.0_p_double

  j = j+1
  diagFile % tracks % quants(j) = 'ene'
  diagFile % tracks % qlabels(j) = 'Ene'
  diagFile % tracks % qunits(j)  = 'm c^2'
  diagFile % tracks % offset_t(j)= -0.5_p_double

  do i = 1, this % get_n_x_dims()
    j = j+1
    diagFile % tracks % quants(j) = 'x'//(char(iachar('0')+i))
    diagFile % tracks % qlabels(j) = 'x_'//(char(iachar('0')+i))
    diagFile % tracks % qunits(j) = 'c/\omega_p'
    diagFile % tracks % offset_t(j)= 0.0_p_double
  enddo

  do i = 1, p_p_dim
    j = j+1
    diagFile % tracks % quants(j) = 'p'//(char(iachar('0')+i))
    diagFile % tracks % qlabels(j) = 'p_'//(char(iachar('0')+i))
    diagFile % tracks % qunits(j) = 'm c'
    diagFile % tracks % offset_t(j)= -0.5_p_double
  enddo

#ifdef __HAS_SPIN__

  do i = 1, p_s_dim
    j = j+1
    diagFile % tracks % quants(j) = 's'//(char(iachar('0')+i))
    diagFile % tracks % qlabels(j) = 's_'//(char(iachar('0')+i))
    diagFile % tracks % qunits(j) = 'h/4\pi'
    diagFile % tracks % offset_t(j)= -0.5_p_double
  enddo

#endif

  do i = 1, p_f_dim
    if (track_set%ifdmp_tracks_efl(i)) then
      j = j+1
      diagFile % tracks % quants(j) = 'E'//(char(iachar('0')+i))
      diagFile % tracks % qlabels(j) = 'E_'//(char(iachar('0')+i))
      diagFile % tracks % qunits(j) = 'm_e c \omega_p/e'
      diagFile % tracks % offset_t(j)= 0.0_p_double
    endif
  enddo

  do i = 1, p_f_dim
    if (track_set%ifdmp_tracks_bfl(i)) then
      j = j+1
      diagFile % tracks % quants(j) = 'B'//(char(iachar('0')+i))
      diagFile % tracks % qlabels(j) = 'B_'//(char(iachar('0')+i))
      diagFile % tracks % qunits(j) = 'm_e c \omega_p/e'
      diagFile % tracks % offset_t(j)= 0.0_p_double
    endif
  enddo

  if (track_set%ifdmp_tracks_psi) then
    j = j+1
    diagFile % tracks % quants(j) = 'psi'
    diagFile % tracks % qlabels(j) = 'Psi'
    diagFile % tracks % qunits(j) = ''
    diagFile % tracks % offset_t(j)= 0.0_p_double
  endif

  diagFile % tracks % nquants = j

end subroutine enumerate_quants_spec
!-----------------------------------------------------------------------------------------

