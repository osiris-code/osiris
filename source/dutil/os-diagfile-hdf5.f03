#ifndef __TEMPLATE__

#define DEBUG_IO 1

#ifdef DEBUG_IO
#define DEBUGMSG(...) print *, __VA_ARGS__
#define CHECK_ERROR( ierr ) if (ierr<0) write(0,"('(*io error*)',A,':',I0)") __FILE__,__LINE__-1
#else
#define DEBUGMSG(...)
#define CHECK_ERROR( ierr )
#endif

#include "os-config.h"

module m_diagfile_hdf5

#include "memory/memory.h"

use m_diagfile
use m_parameters
use hdf5
use mpi

implicit none

private

! TODO: use the ones that are in os-diagfile.f03 rather then cut-n-paste here
!integer, private, parameter :: p_single = kind(1.0e0)
!integer, private, parameter :: p_double = kind(1.0d0)
!integer, private, parameter :: p_byte   = selected_int_kind(2)
!integer, private, parameter :: p_int64  = selected_int_kind(10)

type, extends( t_diag_file ) :: t_diag_file_hdf5

  ! File id
  integer(hid_t)   :: id

  ! Parallel communicator
  integer :: comm

  ! Parallel I/O mode
  integer :: iomode

  ! flag to determine if local parallel node writes to the file
  logical :: write

contains

  procedure :: open_diag_file     => open_file_hdf5
  procedure :: open_par_diag_file => open_par_file_hdf5

  procedure :: close => close_file_hdf5

  procedure :: add_dataset_1D_r4 => add_dataset_1D_r4_hdf5
  procedure :: add_dataset_2D_r4 => add_dataset_2D_r4_hdf5
  procedure :: add_dataset_3D_r4 => add_dataset_3D_r4_hdf5

  procedure :: add_dataset_1D_r8 => add_dataset_1D_r8_hdf5
  procedure :: add_dataset_2D_r8 => add_dataset_2D_r8_hdf5
  procedure :: add_dataset_3D_r8 => add_dataset_3D_r8_hdf5

  ! Gfortran throws the strangest error in 'os-paritcles-define.f90' if a enable this method.
  !   for disable for now since it is unused in the rest of the code
  !procedure :: add_dataset_loc_   => add_dataset_loc_hdf5_

  procedure :: start_cdset => start_cdset_hdf5
  procedure :: write_cdset => write_simple_cdset_hdf5
  procedure :: write_par_cdset  => write_par_cdset_hdf5
  procedure :: end_cdset   => end_cdset_hdf5

  procedure :: open_cdset   => open_cdset_hdf5
  procedure :: extend_cdset => extend_cdset_hdf5
  procedure :: close_cdset  => close_cdset_hdf5

  procedure :: file_extension => file_extension_hdf5

end type t_diag_file_hdf5

public :: t_diag_file_hdf5

contains

pure function file_extension_hdf5( this )

  implicit none
  class( t_diag_file_hdf5 ), intent(in) :: this
  character(len=3) :: file_extension_hdf5
  file_extension_hdf5 = 'h5'

end function

! --------------------------------------------------------------------------------------------------
! Close diagnostics file
! --------------------------------------------------------------------------------------------------
subroutine close_file_hdf5( this )

  implicit none
  class( t_diag_file_hdf5 ), intent(inout) :: this
  integer :: ierr

  if ( this % write ) then
    call h5fclose_f( this%id, ierr)
    CHECK_ERROR( ierr )
  endif

end subroutine close_file_hdf5
! --------------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------------------
! Adds general file attributes
! --------------------------------------------------------------------------------------------------
subroutine add_file_attr_hdf5( this )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this

  integer(hid_t) :: rootID
  integer :: ierr

  if ( this % write ) then

    ! Open root group to add metadata
    call h5gopen_f( this%id, '/', rootID, ierr )
    CHECK_ERROR(ierr)

    ! Add name property
    call hdf5_add( rootID, 'NAME', this%name )

    ! Add general metadata
    select case ( this%ftype )
    case (p_diag_grid)
      call hdf5_add( rootID, "TYPE", "grid" )
      call hdf5_add( rootID, this % grid )
      call hdf5_add( rootID, this % iter )

    case (p_diag_particles)
      call hdf5_add( rootID, "TYPE", "particles");
      call hdf5_add( rootID, this % particles );
      call hdf5_add( rootID, this % iter );

    case (p_diag_tracks)
      call hdf5_add( rootID, "TYPE", "tracks-2");
      call hdf5_add( rootID, this % tracks );

    case default
      write(0,*) '(*error*) Unsupported file type, aborting.'
      call abort_program( p_err_invalid )
    end select

    ! Add simulation info
    call hdf5_add( rootID, "SIMULATION", sim_info )

    ! Close root group
    call h5gclose_f( rootID, ierr )
    CHECK_ERROR( ierr )

  endif
end subroutine add_file_attr_hdf5

! --------------------------------------------------------------------------------------------------
! Open Diagnostic file
! --------------------------------------------------------------------------------------------------
subroutine open_file_hdf5( this, amode )

  use m_grid_define
  use m_system
  use m_parameters
  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  integer, intent(in) :: amode

  character(len=1024) :: path
  integer :: ierr
  integer(hid_t) :: plistID
  integer(hsize_t) :: threshold, alignment

  ! integer(size_t) :: rdcc_nslots, rdcc_nbytes
  ! real :: rdcc_w0

  ! Create directory
  call mkdir( this%filepath, ierr )

  this % write = .true.

  ! open the file
  path = trim(this%filepath)//trim(this%filename)//'.h5'

  ! set to true for all amode
  this % write = .true.

  if ( amode == p_diag_create ) then

    call h5pcreate_f(H5P_FILE_ACCESS_F, plistID, ierr)
    CHECK_ERROR( ierr )

    ! Sets alignment properties
    threshold = 0
    alignment = 262144
    !call h5pset_alignment_f( plistID, threshold, alignment, ierr )
    CHECK_ERROR( ierr )

    call h5fcreate_f( path, H5F_ACC_TRUNC_F, this % id, ierr, access_prp = plistID )
    CHECK_ERROR( ierr )

    call add_file_attr_hdf5( this )

    ! Close the property list
    call h5pclose_f(plistID, ierr)
    CHECK_ERROR( ierr )

  else
    call h5fopen_f( path, H5F_ACC_RDWR_F, this % id, ierr)
    CHECK_ERROR( ierr )

  endif

  this % comm=MPI_COMM_NULL

end subroutine open_file_hdf5

subroutine open_par_file_hdf5( this, amode, comm, iomode )

  use m_grid_define
  use m_system
  use m_parameters
  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  integer, intent(in) :: amode
  integer, intent(in) :: comm
  integer, intent(in), optional :: iomode

  integer :: ierr, rank

#ifdef H5_HAVE_PARALLEL
  character(len=1024) :: path
  integer(hid_t) :: plistID
  integer :: fileInfo
  integer(hsize_t) :: threshold, alignment
#endif

  ! Variables for setting the data chunk cache
  ! integer(hsize_t) :: rdcc_nslots, rdcc_nbytes
  ! real :: rdcc_w0

  if ( present(iomode)) then
    this % iomode = iomode
  else
    this % iomode = diag_file_default_iomode
  endif

  ! If invalid comm just use the serial file open
  ! Ricardo: Look here... does this fit with how you were thinking of things
  !          (I speak of the the ' comm == MPI_COMM_SELF' portion)
  if ( comm == MPI_COMM_NULL .or. comm == MPI_COMM_SELF) then
    call open_file_hdf5( this, amode )
    return
  endif

  ! Only root node creates the directory
  call mpi_comm_rank( comm, rank, ierr )
  if ( rank == 0 ) then
    call mkdir( this%filepath, ierr )
    CHECK_ERROR( ierr )
  endif

  select case ( this % iomode )
  case ( DIAG_MPI )
    ! For MPI access mode only root node writes to the file
    if ( rank == 0 ) then
      this % write = .true.
      call open_file_hdf5( this, amode )
    else
      this % write = .false.
    endif

  case ( DIAG_MPIIO_INDEPENDENT, DIAG_MPIIO_COLLECTIVE )

#ifdef H5_HAVE_PARALLEL

    ! For other parallel I/O methods all nodes are involved in writing the file
    this % write = .true.

    call h5pcreate_f(H5P_FILE_ACCESS_F, plistID, ierr)
    CHECK_ERROR( ierr )

    ! Sets alignment properties
    threshold = 0
    alignment = 524288
    !call h5pset_alignment_f( plistID, threshold, alignment, ierr )
    CHECK_ERROR( ierr )

    fileInfo = MPI_INFO_NULL
    call MPI_INFO_CREATE( fileInfo, ierr )
    CHECK_ERROR( ierr )

    call MPI_INFO_SET( fileInfo, "romio_cb_write", "enable", ierr)
    CHECK_ERROR( ierr )

    call h5pset_fapl_mpio_f(plistID, comm, fileInfo, ierr)
    CHECK_ERROR( ierr )

    path = trim(this%filepath)//trim(this%filename)//'.h5'
    if ( amode == p_diag_create ) then
      call h5fcreate_f( path, H5F_ACC_TRUNC_F, this % id, ierr, access_prp = plistID )
      CHECK_ERROR( ierr )
      call add_file_attr_hdf5( this )
    else
      call h5fopen_f( path, H5F_ACC_RDWR_F, this % id, ierr, access_prp = plistID )
      CHECK_ERROR( ierr )
    endif

    call MPI_INFO_FREE( fileInfo, ierr )
    CHECK_ERROR( ierr )

    ! Close the property list
    call h5pclose_f(plistID, ierr)
    CHECK_ERROR( ierr )
#else
  if ( rank == 0 ) then
    write(0,'(A)') '(*error*) Unable to create HDF5 file, parallel i/o mode not supported.'
    write(0,'(A)') '(*error*) HDF5 library was not compiled with parallel i/o support please'
    write(0,'(A)') '(*error*) use "mpi" mode instead.'
  endif
#endif

  case default

    if ( rank == 0 ) then
      write(0,'(A)') '(*error*) Unable to create HDF5 file, "posix" parallel i/o mode not supported.'
    endif

  end select

  ! Store the communicator
  this % comm = comm

end subroutine open_par_file_hdf5

subroutine add_dataset_loc_hdf5_( this, name, buffer, ndims, count, data_type )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  character(len=*), intent(in) :: name
  type(c_ptr), intent(in) :: buffer
  integer, intent(in) :: data_type, ndims
  integer, dimension(:), intent(in) :: count

  integer(p_byte), dimension(:), pointer :: tmp
  integer(hsize_t), dimension(ndims) :: dims
  integer(hid_t) :: dataspaceID, datasetID, datatypeID
  integer :: bsize, i, ierr

  bsize = diag_sizeof( data_type )
  do i = 1, ndims
    dims(i) = count(i)
    bsize = bsize * count(i)
  enddo
  call c_f_pointer( buffer, tmp, [bsize] )

  ! Convert data type to HDF5 equivalent
  datatypeID = diag_hdf5_type( data_type )

  ! create the dataset
  call h5screate_simple_f( ndims, dims, dataspaceID, ierr )
  CHECK_ERROR( ierr )

  call h5dcreate_f( this % id, name, datatypeID, dataspaceID, datasetID, ierr )
  CHECK_ERROR( ierr )

  ! write the data
  call h5dwrite_f( datasetID, datatypeID, tmp, dims, ierr )
  CHECK_ERROR( ierr )

  ! close the dataset
  call h5sclose_f( dataspaceID, ierr )
  CHECK_ERROR( ierr )

  call h5dclose_f( datasetID, ierr )
  CHECK_ERROR( ierr )


end subroutine add_dataset_loc_hdf5_

! --------------------------------------------------------------------------------------------------
! Start chunked dataset
! --------------------------------------------------------------------------------------------------
subroutine start_cdset_hdf5( this, name, ndims, count, data_type, dset, chunk_size )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  character(len = *), intent(in) :: name
  integer, intent(in) :: ndims
  integer, dimension(:), intent(in) :: count
  integer, intent(in) :: data_type
  type( t_diag_dataset ), intent(inout) :: dset

  integer, dimension(:), intent(in), optional :: chunk_size

  integer(hsize_t), dimension(diag_max_dims) :: dims, maxDims, chunkDims
  integer(hid_t) :: filespaceID, dsetID, propertyID
  integer :: i, ierr

  dset % ndims = ndims
  dset % count(1:ndims) = count(1:ndims)
  dset % data_type = data_type

  if ( this % write ) then

    if ( present(chunk_size) ) then
      ! Create chunked dataset
      do i = 1, dset%ndims
        dims(i) = dset % count(i)
        maxDims(i) = dset % count(i)
        chunkDims(i) = chunk_size(i)
        if ( chunkDims(i) > maxDims(i) .and. maxDims(i) > 0 ) chunkDims(i) = maxDims(i)
      enddo
      maxDims( dset%ndims ) = H5S_UNLIMITED_F

      call h5screate_simple_f( dset % ndims, dims, filespaceID, ierr, &
                               maxdims = maxDims )
      CHECK_ERROR( ierr )

      ! Using unlimited dataspaces requires chunking
      call h5pcreate_f( H5P_DATASET_CREATE_F, propertyID, ierr )
      CHECK_ERROR( ierr )

      call h5pset_chunk_f( propertyID, dset%ndims, chunkDims, ierr )
      CHECK_ERROR( ierr )

      call h5dcreate_f( this%id, trim(name), diag_hdf5_type( dset%data_type), filespaceID, dsetID, ierr, &
                        dcpl_id = propertyID )
      CHECK_ERROR( ierr )

      call h5pclose_f( propertyID, ierr )
      CHECK_ERROR( ierr )
    else
      ! Create contiguos dataset
      do i = 1, dset%ndims
        dims(i) = dset % count(i)
      enddo

      if ( dims(1) > 0 ) then
        call h5screate_simple_f( dset % ndims, dims, filespaceID, ierr )
        CHECK_ERROR( ierr )
      else
        call h5screate_f( H5S_SCALAR_F, filespaceID, ierr )
        CHECK_ERROR( ierr )
      endif

      call h5dcreate_f( this%id, trim(name), diag_hdf5_type( dset%data_type), filespaceID, dsetID, ierr )
      CHECK_ERROR( ierr )
    endif

    call h5sclose_f(filespaceID, ierr)
    CHECK_ERROR( ierr )

    ! Store dataset ID
    dset % id = dsetID


  endif

end subroutine start_cdset_hdf5

! --------------------------------------------------------------------------------------------------
! Write part of a chunked dataset (not parallel)
! --------------------------------------------------------------------------------------------------
subroutine write_simple_cdset_hdf5( this, dset, chunk )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset
  type( t_diag_chunk ), intent(inout) :: chunk

  integer(hid_t) :: filespaceID, memspaceID, datasetID
  integer(hsize_t), dimension(diag_max_dims) :: start, count
  integer(p_byte), dimension(:), pointer :: buffer
  integer :: bsize, i, ierr

  datasetID = dset % id

  bsize = diag_sizeof( dset % data_type )
  do i = 1, dset % ndims
    start(i) = chunk % start(i)
    count(i) = chunk % count(i)
    bsize = bsize * chunk % count(i)
  enddo
  call c_f_pointer( chunk % data, buffer, [bsize] )

  ! create memory dataspace
  call h5screate_simple_f( dset%ndims, count, memspaceID, ierr)
  CHECK_ERROR( ierr )

  ! select hyperslab in the file
  call h5dget_space_f( datasetID, filespaceID, ierr )
  CHECK_ERROR( ierr )
  call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr )
  CHECK_ERROR( ierr )

  ! write data
  call h5dwrite_f(datasetID, diag_hdf5_type( dset % data_type ), buffer, count, ierr, &
      file_space_id = filespaceID, mem_space_id = memspaceID)

  ! close resources
  call h5sclose_f(filespaceID, ierr)
  CHECK_ERROR( ierr )
  call h5sclose_f(memspaceID, ierr)
  CHECK_ERROR( ierr )

end subroutine write_simple_cdset_hdf5

! --------------------------------------------------------------------------------------------------
! Write parallel chunked dataset using mpi/io
! --------------------------------------------------------------------------------------------------
#ifdef H5_HAVE_PARALLEL

subroutine write_cdset_mpio( this, dset, chunk, mode )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset
  type( t_diag_chunk ), intent(inout) :: chunk
  integer, intent(in) :: mode

  integer(hid_t) :: dcplid, xferid, datasetID, filespaceID, memspaceID
  integer(hsize_t), dimension(diag_max_dims) :: dimsChunk, dimsFile, start, count
  integer(p_byte), dimension(:), pointer :: buffer
  integer :: bsize, i, ierr

  datasetID = dset % id

  ! crete property list for datset creation
  call h5pcreate_f(H5P_DATASET_CREATE_F, dcplID, ierr)
  CHECK_ERROR( ierr )

  ! create property list for dataset transfer
  call h5pcreate_f(H5P_DATASET_XFER_F, xferID, ierr)
  CHECK_ERROR( ierr )

  bsize = diag_sizeof( dset % data_type )
  do i = 1, dset % ndims
    dimsChunk(i) = chunk % count(i)
    dimsFile(i)  = dset % count(i)
    bsize = bsize * chunk % count(i)
  enddo
  call c_f_pointer( chunk % data, buffer, [bsize] )

  ! create file dataspace
  call h5screate_simple_f(dset % ndims, dimsFile, filespaceID, ierr )
  CHECK_ERROR( ierr )

  ! create memory dataspace
  call h5screate_simple_f(dset % ndims, dimsChunk, memspaceID, ierr )
  CHECK_ERROR( ierr )

  ! close resources
  call h5sclose_f(filespaceID, ierr)
  CHECK_ERROR( ierr )

  do i = 1, dset % ndims
    ! hyperslab coordinates are 0 indexed
    start(i) = chunk % start(i)
    count(i) = chunk % count(i)
  enddo

  ! select hyperslab in the file
  call h5dget_space_f( datasetID, filespaceID, ierr )
  CHECK_ERROR( ierr )
  call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr )
  CHECK_ERROR( ierr )

  ! write the dataset collectively or independently.
  call h5pset_dxpl_mpio_f(xferID, mode, ierr)
  CHECK_ERROR( ierr )

  call h5dwrite_f(datasetID, diag_hdf5_type( dset%data_type ), buffer, dimsFile, ierr, &
     file_space_id = filespaceID, mem_space_id = memspaceID, xfer_prp = xferID)

  ! close resources
  call h5sclose_f(filespaceID, ierr)
  CHECK_ERROR( ierr )
  call h5sclose_f(memspaceID, ierr)
  CHECK_ERROR( ierr )

end subroutine write_cdset_mpio

#endif

subroutine write_cdset_mpi( this, dset, chunk )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset
  type( t_diag_chunk ), intent(inout) :: chunk

  integer, parameter :: ping_tag = 1001, tile_tag = 1002, data_tag = 1003
  integer :: ping, data_size, comm_rank, comm_size, ping_req, data_req
  integer :: step, src, i, ierr, bsize
  integer, dimension( 2 * diag_max_dims ) :: cinfo0, cinfo1
  integer(p_byte), dimension(:), pointer :: buffer0, buffer1, bufferw
  integer :: bsize0, bsize1, bsizew
  integer(hsize_t), dimension(diag_max_dims) :: start, count
  integer(hid_t) :: datasetID, memspaceID, filespaceID

  datasetID = dset % id

  ping = 1234

  ! Get size of data element
  data_size = diag_sizeof( dset % data_type )

  ! Point buffer0 to local data
  bsize0 = 1
  do i = 1, dset % ndims
    cinfo0(                  i ) = chunk % count(i)
    cinfo0(   dset % ndims + i ) = chunk % start(i)
    bsize0 = bsize0 * chunk % count(i)
  enddo
  call c_f_pointer( chunk % data, buffer0, [bsize0*data_size] )

  call MPI_COMM_RANK( this%comm, comm_rank, ierr )

  if ( comm_rank == 0 ) then

    ! Get filespace
    call h5dget_space_f(datasetID, filespaceID, ierr)
    CHECK_ERROR( ierr )

    ! Loop over all nodes
    call MPI_COMM_SIZE( this%comm, comm_size, ierr )

    bsize0 = 0
    bsize1 = 0
    step = 0
    do src = 0, comm_size-1
      ! Start receiving data from next node
      if ( src < comm_size - 1 ) then
        ! send ping
        call MPI_ISEND( ping, 1, MPI_INTEGER, src+1, ping_tag, this % comm, ping_req, ierr )

        if ( step == 0 ) then
          ! Receive on buffer1
          call MPI_RECV( cinfo1, 2 * dset%ndims, MPI_INTEGER, src+1, tile_tag, this%comm, MPI_STATUS_IGNORE, ierr)
          bsize = cinfo1(1)
          do i = 2, dset % ndims
            bsize = bsize * cinfo1(i)
          enddo
          if ( bsize > bsize1 ) then
            if ( bsize1 > 0 ) call freemem( buffer1 )
            call alloc( buffer1, [bsize * data_size])
            bsize1 = bsize
          endif

          if( bsize .ne. 0 ) call MPI_IRECV( buffer1, bsize, diag_mpi_type(dset%data_type), src+1, data_tag, this%comm, data_req, ierr)

        else
          ! Receive on buffer0
          call MPI_RECV( cinfo0, 2 * dset%ndims, MPI_INTEGER, src+1, tile_tag, this%comm, MPI_STATUS_IGNORE, ierr)
          bsize = cinfo0(1)
          do i = 2, dset % ndims
            bsize = bsize * cinfo0(i)
          enddo
          if ( bsize > bsize0 ) then
            if ( bsize0 > 0 ) call freemem( buffer0 )
            call alloc( buffer0, [bsize * data_size])
            bsize0 = bsize
          endif

          if( bsize .ne. 0 ) call MPI_IRECV( buffer0, bsize, diag_mpi_type(dset%data_type), src+1, data_tag, this%comm, data_req, ierr)
        endif
      endif

      ! write data from current node
      if ( step == 0 ) then
        bufferw => buffer0
        do i = 1, dset % ndims
          count(i)  = cinfo0(                i )
          start(i)  = cinfo0( dset % ndims + i )
        enddo
      else
        bufferw => buffer1
        do i = 1, dset % ndims
          count(i)  = cinfo1(                i )
          start(i)  = cinfo1( dset % ndims + i )
        enddo
      endif

      bsizew = 1
      do i = 1, dset % ndims
        bsizew = bsizew * count(i)
      enddo

      if ( bsizew .ne. 0 ) then

        call h5screate_simple_f( dset % ndims, count, memspaceID, ierr)
        CHECK_ERROR( ierr )

        call h5sselect_hyperslab_f( filespaceID, H5S_SELECT_SET_F, start, count, ierr )
        CHECK_ERROR( ierr )

        call h5dwrite_f( datasetID,  diag_hdf5_type( dset % data_type ), bufferw, count, ierr, &
          mem_space_id = memspaceID, file_space_id = filespaceID)
        CHECK_ERROR( ierr )

        ! close resources
        call h5sclose_f(memspaceID, ierr)
        CHECK_ERROR( ierr )

      endif

      ! Wait for messages to complete
      if ( src < comm_size - 1 ) then
        call MPI_WAIT( ping_req, MPI_STATUS_IGNORE, ierr )
        if ( bsize .ne. 0 ) call MPI_WAIT( data_req, MPI_STATUS_IGNORE, ierr )
      endif

      step = 1 - step
    enddo

    ! Free communication buffers
    if ( bsize0 > 0 ) call freemem( buffer0 )
    if ( bsize1 > 0 ) call freemem( buffer1 )

    ! close resources
    call h5sclose_f( filespaceID, ierr )

  else
    ! Wait for ping from root node
    call MPI_RECV( ping, 1, MPI_INTEGER, 0, ping_tag, this%comm, MPI_STATUS_IGNORE, ierr)

    ! send local chunk information to root node
    call MPI_SEND( cinfo0, 2 * dset%ndims, MPI_INTEGER, 0, tile_tag, this%comm, ierr)

    if( bsize0 > 0 ) then
      ! Send chunk data to root node
      call MPI_SEND( buffer0, bsize0, diag_mpi_type( dset % data_type ), 0, data_tag, this%comm, ierr )
    endif

  endif

end subroutine write_cdset_mpi

subroutine write_par_cdset_hdf5( this, dset, chunk, offset_ )

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset
  type( t_diag_chunk ), intent(inout) :: chunk

  ! Offset is not used by HDF5 routines
  integer, intent(in), optional :: offset_

  integer :: rank, ierr

  ! RICARDO LOOK HERE
  ! If comm is MPI_COMM_NULL switch to the serial write code,
  !   but the code should not get here so this is a bit of a hack.
  !   But not so much of a hack since the parallel open function
  !   switches back to serial... so it's kinda ok but maybe not..
  !
  ! Also, Ricardo: does the MPI_COMM_SELF fit with your thinking here?
  if ( this%comm == MPI_COMM_NULL .or. this%comm == MPI_COMM_SELF ) then
    call this%write_cdset( dset, chunk )
    return
  endif

  select case ( this % iomode )
  case(DIAG_MPI)
    call write_cdset_mpi( this, dset, chunk )
#ifdef H5_HAVE_PARALLEL
  case(DIAG_MPIIO_INDEPENDENT)
    call write_cdset_mpio( this, dset, chunk, H5FD_MPIO_INDEPENDENT_F )
  case(DIAG_MPIIO_COLLECTIVE)
    call write_cdset_mpio( this, dset, chunk, H5FD_MPIO_COLLECTIVE_F )
#endif
  case default
    call mpi_comm_rank( this % comm, rank, ierr )
    if ( rank == 0 ) then
      write(0,'(A)') '(*error*) Invalid parallel IO mode, unable to write dataset.'
    endif
  end select

end subroutine write_par_cdset_hdf5

! --------------------------------------------------------------------------------------------------
! Write end marker for a chunked dataset
! --------------------------------------------------------------------------------------------------
subroutine end_cdset_hdf5( this, dset )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset

  integer :: ierr
  integer(hid_t) :: dsetID

  if ( this % write ) then
    dsetID = dset % id
    call h5dclose_f( dsetID, ierr )
    CHECK_ERROR( ierr )
  endif

end subroutine end_cdset_hdf5

! --------------------------------------------------------------------------------------------------
! Open a chunked dataset for updating
! --------------------------------------------------------------------------------------------------
subroutine open_cdset_hdf5( this, dset )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset

  integer :: i, ierr
  integer(hsize_t), dimension(diag_max_dims) :: dims, maxdims
  integer(hid_t) :: datasetID, datatypeID, filespaceID

  if ( .not. this % write ) return

  call h5dopen_f( this%id, trim(dset%name), datasetID, ierr )
  CHECK_ERROR( ierr )

  ! Get datatype
  call h5dget_type_f(datasetID, datatypeID, ierr )
  CHECK_ERROR( ierr )
  dset % data_type = hdf5_diag_type( datatypeID )
  call h5tclose_f( datatypeID, ierr )
  CHECK_ERROR( ierr )

  ! Get dimensions
  call h5dget_space_f( datasetID, filespaceID, ierr )
  CHECK_ERROR( ierr )

  call h5sget_simple_extent_ndims_f(filespaceID, dset % ndims, ierr )
  CHECK_ERROR( ierr )

  call h5sget_simple_extent_dims_f( filespaceID, dims, maxdims, ierr )
  CHECK_ERROR( ierr )

  do i = 1, dset % ndims
    dset % count(i) = dims(i)
  enddo

  ! Store dataset ID
  dset % id = datasetID

  ! The offset is not used for HDF5 files
  dset % offset = -1

end subroutine open_cdset_hdf5

! --------------------------------------------------------------------------------------------------
! Extend dimensions of chunked dataset
! --------------------------------------------------------------------------------------------------
subroutine extend_cdset_hdf5( this, dset, dims )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset
  integer(p_int64), dimension(:), intent(in) :: dims

  integer :: i, ierr
  integer(hid_t) :: datasetID

  datasetID = dset % id
  call h5dset_extent_f( datasetID, dims, ierr )
  CHECK_ERROR( ierr )

  do i = 1, dset % ndims
    dset % count(i) = dims(i)
  enddo

end subroutine extend_cdset_hdf5

! --------------------------------------------------------------------------------------------------
! Close previously open chunked dataset
! --------------------------------------------------------------------------------------------------
subroutine close_cdset_hdf5( this, dset )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset

  integer :: ierr
  integer(hid_t) :: datasetID

  datasetID = dset % id
  call h5dclose_f( datasetID, ierr )
  CHECK_ERROR( ierr )

  dset % ndims = 0
  dset % id = -1

end subroutine close_cdset_hdf5

! --------------------------------------------------------------------------------------------------
! add_dataset routines
! --------------------------------------------------------------------------------------------------

#define __TEMPLATE__

#define __TYPE__ real(p_single)
#define __DIMS__ :
#define __NDIMS__ 1
#define __H5_TYPE__ H5T_NATIVE_REAL
#define FNAME(a) a ## _1d_r4_hdf5
#include __FILE__

#define __TYPE__ real(p_single)
#define __DIMS__ :,:
#define __NDIMS__ 2
#define __H5_TYPE__ H5T_NATIVE_REAL
#define FNAME(a) a ## _2d_r4_hdf5
#include __FILE__

#define __TYPE__ real(p_single)
#define __DIMS__ :,:,:
#define __NDIMS__ 3
#define __H5_TYPE__ H5T_NATIVE_REAL
#define FNAME(a) a ## _3d_r4_hdf5
#include __FILE__

#define __TYPE__ real(p_double)
#define __DIMS__ :
#define __NDIMS__ 1
#define __H5_TYPE__ H5T_NATIVE_DOUBLE
#define FNAME(a) a ## _1d_r8_hdf5
#include __FILE__

#define __TYPE__ real(p_double)
#define __DIMS__ :,:
#define __NDIMS__ 2
#define __H5_TYPE__ H5T_NATIVE_DOUBLE
#define FNAME(a) a ## _2d_r8_hdf5
#include __FILE__

#define __TYPE__ real(p_double)
#define __DIMS__ :,:,:
#define __NDIMS__ 3
#define __H5_TYPE__ H5T_NATIVE_DOUBLE
#define FNAME(a) a ## _3d_r8_hdf5
#include __FILE__


end module m_diagfile_hdf5

#else

!---------------------------------------------------------------------------------------------------
!  Template Function definitions
!---------------------------------------------------------------------------------------------------

subroutine FNAME(add_dataset)( this, name, buffer )

  use hdf5_util

  implicit none

  class( t_diag_file_hdf5 ), intent(inout) :: this
  character(len=*), intent(in) :: name
  __TYPE__, dimension(__DIMS__), intent(inout), target :: buffer

  integer(hsize_t), dimension(__NDIMS__) :: dims
  integer(hid_t) :: dataspaceID, datasetID
  integer :: i, ierr

  do i = 1, __NDIMS__
    dims(i) = size(buffer, i)
  enddo

  ! create the dataset
  call h5screate_simple_f( __NDIMS__, dims, dataspaceID, ierr )
  CHECK_ERROR( ierr )

  call h5dcreate_f( this % id, name, __H5_TYPE__, dataspaceID, datasetID, ierr )
  CHECK_ERROR( ierr )

  ! write the data
  call h5dwrite_f( datasetID, __H5_TYPE__, buffer, dims, ierr )
  CHECK_ERROR( ierr )

  ! close the dataset
  call h5sclose_f( dataspaceID, ierr )
  CHECK_ERROR( ierr )

  call h5dclose_f( datasetID, ierr )
  CHECK_ERROR( ierr )


end subroutine FNAME(add_dataset)


#undef __TYPE__
#undef __DIMS__
#undef __NDIMS__
#undef __H5_TYPE__
#undef FNAME

#endif
