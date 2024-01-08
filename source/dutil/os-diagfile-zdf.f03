#ifndef __TEMPLATE__

#include "os-config.h"
#include "os-preprocess.fpp"

module m_diagfile_zdf

use m_diagfile
use zdf
use zdf_parallel
use m_system
use m_parameters

private

! Diagnostic file in ZDF format
type, extends( t_diag_file ) :: t_diag_file_zdf

    type( t_zdf_file ) :: file
    type( t_zdf_par_file ) :: parFile

contains

  procedure :: open_diag_file     => open_file_zdf
  procedure :: open_par_diag_file => open_par_file_zdf
  procedure :: close              => close_file_zdf

  procedure :: add_dataset_1D_r4 => add_dataset_1D_r4_zdf
  procedure :: add_dataset_2D_r4 => add_dataset_2D_r4_zdf
  procedure :: add_dataset_3D_r4 => add_dataset_3D_r4_zdf

  procedure :: add_dataset_1D_r8 => add_dataset_1D_r8_zdf
  procedure :: add_dataset_2D_r8 => add_dataset_2D_r8_zdf
  procedure :: add_dataset_3D_r8 => add_dataset_3D_r8_zdf

  ! Gfortran throws the strangest error in 'os-paritcles-define.f90' if a enable this method.
  !   for disable for now since it is unused in the rest of the code
  !procedure :: add_dataset_loc_   => add_dataset_loc_zdf_

  procedure :: start_cdset       => start_cdset_zdf
  procedure :: write_cdset       => write_cdset_zdf
  procedure :: write_par_cdset   => write_par_cdset_zdf
  procedure :: end_cdset         => end_cdset_zdf

  procedure :: open_cdset        => open_cdset_zdf
  procedure :: extend_cdset      => extend_cdset_zdf
  procedure :: close_cdset       => close_cdset_zdf

  procedure :: file_extension    => file_extension_zdf

end type t_diag_file_zdf


public :: t_diag_file_zdf

interface diag_to_zdf
    module procedure diag_to_zdf_dataset
    module procedure diag_to_zdf_grid
    module procedure diag_to_zdf_axis
    module procedure diag_to_zdf_part
    module procedure diag_to_zdf_tracks
    module procedure diag_to_zdf_iteration
end interface

interface zdf_to_diag
    module procedure zdf_to_diag_dataset
end interface

contains

! -----------------------------------------------------------------------------
!> Gets ZDF data type from t_diag_file
!@ param    t_diag_file data type
!@ return   ZDF data type
function diag_zdf_type( diag_type ) result( zdf_type )

    implicit none

    integer, intent(in) :: diag_type
    integer :: zdf_type

    select case ( diag_type )
    case( diag_null )
        zdf_type = zdf_null
    case( diag_int8 )
        zdf_type = zdf_int8
    case( diag_uint8 )
        zdf_type = zdf_uint8
    case( diag_int16 )
        zdf_type = zdf_int16
    case( zdf_uint16 )
        zdf_type = zdf_uint16
    case( diag_int32 )
        zdf_type = zdf_int32
    case( diag_uint32 )
        zdf_type = zdf_uint32
    case( diag_int64 )
        zdf_type = zdf_int64
    case( diag_uint64 )
        zdf_type = zdf_uint64
    case( diag_float32 )
        zdf_type = zdf_float32
    case( diag_float64 )
        zdf_type = zdf_float64
    case default
        zdf_type = zdf_null
    end select

end function

! -----------------------------------------------------------------------------
!> Converts t_diag_dataset object to a t_zdf_dataset object
!@ param    diag    t_diag_object
!@ param    zdf     t_zdf_object
subroutine diag_to_zdf_dataset( diag_dset, zdf_dset )

    implicit none

    type(t_diag_dataset), intent(in) :: diag_dset
    type(t_zdf_dataset), intent(inout) :: zdf_dset

    integer :: i

    zdf_dset % name      = diag_dset % name
    zdf_dset % data_type = diag_zdf_type( diag_dset%data_type) 
    zdf_dset % ndims     = diag_dset % ndims
    do i = 1, diag_dset % ndims
        zdf_dset % count(i) = diag_dset % count(i)
    enddo
    zdf_dset % data = diag_dset % data
    zdf_dset % id   = diag_dset % id
    zdf_dset % offset = diag_dset % offset

end subroutine diag_to_zdf_dataset

! -----------------------------------------------------------------------------
!> Converts t_diag_axis_info object to a t_zdf_axis_info object
!@ param    diag    t_diag_grid_axis
!@ param    zdf     t_zdf_grid_axis
subroutine diag_to_zdf_axis( diag_axis, zdf_axis )

    implicit none

    type(t_diag_grid_axis), intent(in) :: diag_axis
    type(t_zdf_grid_axis), intent(inout) :: zdf_axis

    zdf_axis % type  = diag_axis % type
    zdf_axis % min   = diag_axis % min
    zdf_axis % max   = diag_axis % max
    zdf_axis % name  = diag_axis % name
    zdf_axis % label = diag_axis % label
    zdf_axis % units = diag_axis % units

end subroutine diag_to_zdf_axis

! -----------------------------------------------------------------------------
!> Converts t_diag_grid_info object to a t_zdf_grid_info object
!@ param    diag    t_diag_grid_info
!@ param    zdf     t_zdf_grid_info
subroutine diag_to_zdf_grid( diag_grid, zdf_grid )

    implicit none

    type(t_diag_grid_info), intent(in) :: diag_grid
    type(t_zdf_grid_info), intent(inout) :: zdf_grid

    integer :: i

    zdf_grid % ndims = diag_grid % ndims
    zdf_grid % name = diag_grid % name
    zdf_grid % label = diag_grid % label
    zdf_grid % units = diag_grid % units

    do i = 1, zdf_grid % ndims
        zdf_grid % count(i) = diag_grid % count(i)
        call diag_to_zdf( diag_grid % axis(i), zdf_grid % axis(i) )
    enddo

end subroutine diag_to_zdf_grid


! -----------------------------------------------------------------------------
!> Converts t_diag_part_info object to a t_zdf_part_info object
!@ param    diag    t_diag_part_info
!@ param    zdf     t_zdf_part_info
subroutine diag_to_zdf_part( diag_part, zdf_part )

    implicit none

    type(t_diag_part_info), intent(in) :: diag_part
    type(t_zdf_part_info), intent(inout) :: zdf_part

    integer :: i

    zdf_part % name = diag_part % name
    zdf_part % label = diag_part % label
    zdf_part % np = diag_part % np
    zdf_part % nquants = diag_part % nquants
    zdf_part % has_tags = diag_part % has_tags

    do i = 1, zdf_part % nquants
        zdf_part % quants(i) = diag_part % quants(i)
        zdf_part % qlabels(i) = diag_part % qlabels(i)
        zdf_part % qunits(i) = diag_part % qunits(i)
    enddo
    
end subroutine diag_to_zdf_part

! -----------------------------------------------------------------------------
!> Converts t_diag_tracks_info object to a t_zdf_tracks_info object
!@ param    diag    t_diag_tracks_info
!@ param    zdf     t_zdf_tracks_info
subroutine diag_to_zdf_tracks( diag_tracks, zdf_tracks )

    implicit none

    type(t_diag_track_info), intent(in) :: diag_tracks
    type(t_zdf_track_info), intent(inout) :: zdf_tracks

    integer :: i

    zdf_tracks % name      = diag_tracks % name
    zdf_tracks % label     = diag_tracks % label
    zdf_tracks % ntracks   = diag_tracks % ntracks
    zdf_tracks % ndump     = diag_tracks % ndump
    zdf_tracks % niter     = diag_tracks % niter
    zdf_tracks % nquants   = diag_tracks % nquants

    do i = 1, zdf_tracks % nquants
        zdf_tracks % quants(i)  = diag_tracks % quants(i)
        zdf_tracks % qlabels(i) = diag_tracks % qlabels(i)
        zdf_tracks % qunits(i)  = diag_tracks % qunits(i)
    enddo

end subroutine diag_to_zdf_tracks

! -----------------------------------------------------------------------------
!> Converts t_diag_iteration_info object to a t_zdf_iteration_info object
!@ param    diag    t_diag_iteration_info
!@ param    zdf     t_zdf_iteration_info
subroutine diag_to_zdf_iteration( diag_iteration, zdf_iteration )

    implicit none

    type(t_diag_iteration), intent(in) :: diag_iteration
    type(t_zdf_iteration), intent(inout) :: zdf_iteration

    
    zdf_iteration % name = diag_iteration % name
    zdf_iteration % n    = diag_iteration % n
    zdf_iteration % t    = diag_iteration % t
    zdf_iteration % time_units = diag_iteration % time_units

end subroutine diag_to_zdf_iteration

! -----------------------------------------------------------------------------
!> Converts t_zdf_dataset object to a t_diag_dataset object
!@ param    zdf     t_zdf_object
!@ param    diag    t_diag_object
subroutine zdf_to_diag_dataset( zdf_dset, diag_dset )

    implicit none

    type(t_zdf_dataset), intent(in) :: zdf_dset
    type(t_diag_dataset), intent(inout) :: diag_dset

    integer :: i

    diag_dset % name      = zdf_dset % name
    diag_dset % data_type = diag_zdf_type( zdf_dset%data_type ) 
    diag_dset % ndims     = zdf_dset % ndims
    do i = 1, zdf_dset % ndims
        diag_dset % count(i) = zdf_dset % count(i)
    enddo
    diag_dset % data = zdf_dset % data
    diag_dset % id   = zdf_dset % id
    diag_dset % offset = zdf_dset % offset

end subroutine zdf_to_diag_dataset


! -----------------------------------------------------------------------------
pure function file_extension_zdf( this ) result(ext)

  implicit none

  class( t_diag_file_zdf ), intent(in) :: this

  character(len=3) :: ext
  ext = 'zdf'

end function

! --------------------------------------------------------------------------------------------------
! Open Diagnostic file
! --------------------------------------------------------------------------------------------------
subroutine open_file_zdf( this, amode )

  use m_grid_define
  use m_system
  use m_parameters

  implicit none

  class( t_diag_file_zdf ), intent(inout) :: this
  integer, intent(in) :: amode

  character(len=1024) :: path
  integer :: ierr
  type( t_zdf_grid_info ) :: grid
  type( t_zdf_part_info ) :: particles
  type( t_zdf_iteration ) :: iter
  type( t_zdf_track_info ) :: tracks


  ! Serial I/O
  this % parallel = .false.

  call mkdir( this%filepath, ierr )

  path = trim(this%filepath)//trim(this%filename)//'.zdf'
  call zdf_open_file( this % file, path, amode )

  if ( amode == p_diag_create ) then
    select case ( this%ftype )
    case (p_diag_grid)
      call zdf_add( this%file, "TYPE", "grid" )
      call diag_to_zdf( this % grid, grid ); call zdf_add( this%file, grid )
      call diag_to_zdf( this % iter, iter ); call zdf_add( this%file, iter )

    case (p_diag_particles)
      call zdf_add( this%file, "TYPE", "particles");
      call diag_to_zdf( this % particles, particles ); call zdf_add( this%file, particles );
      call diag_to_zdf( this % iter, iter ); call zdf_add( this%file, iter );

    case (p_diag_tracks)
      call zdf_add( this%file, "TYPE", "tracks-2");
      call diag_to_zdf( this % tracks, tracks ); call zdf_add( this%file, tracks );

    case default
      write(0,'(A,I0)') "(*error*) Unsupported file type: ",  this%ftype
      write(0,'(A,A)') "(*error*) path: ", trim(path)
      call abort_program( p_err_invalid )
    end select
  endif

end subroutine open_file_zdf

subroutine open_par_file_zdf( this, amode, comm, iomode )

  use m_grid_define
  use m_system
  use m_parameters

  implicit none

  class( t_diag_file_zdf ), intent(inout) :: this
  integer, intent(in) :: amode
  integer, intent(in) :: comm
  integer, intent(in), optional :: iomode

  character(len=1024) :: path
  integer :: ierr, rank, iomode_

  type( t_zdf_grid_info ) :: grid
  type( t_zdf_part_info ) :: particles
  type( t_zdf_iteration ) :: iter

  if ( comm == MPI_COMM_NULL ) then
    ! Use serial IO
    this % parallel = .false.
    call open_file_zdf( this, amode )
  else
    ! Use parallel IO
    this % parallel = .true.

    if ( present(iomode)) then
      iomode_ = iomode
    else
      iomode_ = diag_file_default_iomode
    endif

    ! Create directory
    call mpi_comm_rank( comm, rank, ierr )
    if ( rank == 0 ) then
      call mkdir( this%filepath, ierr )
    endif
    ! Create file
    path = trim(this%filepath)//trim(this%filename)//'.zdf'
    call zdf_open_file( this % parFile, path, amode, comm, iomode_ )

    ! Add metadata
    if ( amode == p_diag_create ) then
      select case ( this%ftype )
      case (p_diag_grid)
        call zdf_add( this%parFile, "TYPE", "grid" )
        call diag_to_zdf( this % grid, grid ); call zdf_add( this%parFile, grid )
        call diag_to_zdf( this % iter, iter ); call zdf_add( this%parFile, iter )
  
      case (p_diag_particles)
        call zdf_add( this%parFile, "TYPE", "particles");
        call diag_to_zdf( this % particles, particles ); call zdf_add( this%parFile, particles );
        call diag_to_zdf( this % iter, iter ); call zdf_add( this%parFile, iter );
  
      case default
        write(0,'(A,I0)') "(*error*) Unsupported file type: ",  this%ftype
        write(0,'(A,A)') "(*error*) path: ", trim(path)
        call abort_program( p_err_invalid )
      end select
    endif
  endif

end subroutine open_par_file_zdf

! --------------------------------------------------------------------------------------------------
! Close diagnostics file
! --------------------------------------------------------------------------------------------------
subroutine close_file_zdf( this )

  use m_grid_define

  implicit none
  class( t_diag_file_zdf ), intent(inout) :: this

  if ( this % parallel ) then
    call zdf_close_file( this % parFile )
  else
    call zdf_close_file( this % file )
  endif

end subroutine close_file_zdf
! --------------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------------------
! Add dataset to file 
! --------------------------------------------------------------------------------------------------
subroutine add_dataset_loc_zdf_( this, name, buffer, ndims, count, data_type )

  use zdf_parallel

  implicit none

  class( t_diag_file_zdf ), intent(inout) :: this
  character(len=*), intent(in) :: name
  type(c_ptr), intent(in) :: buffer
  integer, intent(in) :: data_type, ndims
  integer, dimension(:), intent(in) :: count

  type(t_zdf_dataset) :: dset

  dset % name = name
  dset % data_type =   data_type
  dset % ndims     = ndims
  dset % count(1:ndims) = count(1:ndims)
  dset % data = buffer
  dset % id = -1
  dset % offset = 0

  call zdf_add( this%file, dset )

end subroutine add_dataset_loc_zdf_

! --------------------------------------------------------------------------------------------------
! Start chunked dataset
! --------------------------------------------------------------------------------------------------
subroutine start_cdset_zdf( this, name, ndims, count, data_type, dset, chunk_size )

  implicit none

  class( t_diag_file_zdf ), intent(inout) :: this
  character(len = *), intent(in) :: name
  integer, intent(in) :: ndims
  integer, dimension(:), intent(in) :: count
  integer, intent(in) :: data_type
  type( t_diag_dataset ), intent(inout) :: dset
  integer, dimension(:), intent(in), optional :: chunk_size

  type( t_zdf_dataset ):: zdf_dset

  zdf_dset % name = name
  zdf_dset % ndims = ndims
  zdf_dset % count(1:ndims) = count(1:ndims)
  zdf_dset % data_type = data_type
  zdf_dset % data = c_null_ptr
  zdf_dset % id = -1
  zdf_dset % offset = 0

  if ( present(chunk_size) ) then
    ! chunk size is ignored for ZDF files
    continue
  endif

  if ( this % parallel ) then
    call zdf_start_cdset( this % parFile, zdf_dset )
  else
    call zdf_start_cdset( this % file, zdf_dset )
  endif

  call zdf_to_diag( zdf_dset, dset )

end subroutine start_cdset_zdf

! --------------------------------------------------------------------------------------------------
! Write part of a chunked dataset
! --------------------------------------------------------------------------------------------------
subroutine write_cdset_zdf( this, dset, chunk )

  implicit none

  class( t_diag_file_zdf ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset
  type( t_diag_chunk ), intent(inout) :: chunk

  type(t_zdf_chunk) :: zdf_chunk
  type( t_zdf_dataset ) :: zdf_dset

    ! Check if the file is a serial file
  if ( this % parallel ) then
    write(0,'(A)') "(*error*) Cannot write simple chunk on parallel file, aborting..."
    call abort_program( p_err_invalid )
  endif

  call diag_to_zdf( dset, zdf_dset )

  zdf_chunk % count  = chunk % count 
  zdf_chunk % start  = chunk % start 
  zdf_chunk % stride = chunk % stride
  zdf_chunk % data   = chunk % data  

  call zdf_write_cdset( this % file, zdf_dset, zdf_chunk )

end subroutine write_cdset_zdf


! --------------------------------------------------------------------------------------------------
! Write all parts of a chunked dataset distributed over a parallel partition
! --------------------------------------------------------------------------------------------------
subroutine write_par_cdset_zdf( this, dset, chunk, offset_ )

  implicit none

  class( t_diag_file_zdf ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset
  type( t_diag_chunk ), intent(inout) :: chunk
  integer, intent(in), optional :: offset_

  integer(c_int64_t) :: offset
  type( t_zdf_dataset ) :: zdf_dset
  type(t_zdf_chunk) :: zdf_chunk

  call diag_to_zdf( dset, zdf_dset )
  zdf_chunk % count  = chunk % count 
  zdf_chunk % start  = chunk % start 
  zdf_chunk % stride = chunk % stride
  zdf_chunk % data   = chunk % data  

  ! Check if the file is a parallel file
  if ( this % parallel ) then

    if ( present(offset_) ) then
        ! Use the supplied offset
        offset = offset_
    else
        ! Calculate offset if required
        offset = ZDF_GET_OFFSET
    endif

    call zdf_par_write_cdset( this % parFile, zdf_dset, zdf_chunk, offset )

  else
    call zdf_write_cdset( this % file, zdf_dset, zdf_chunk )
  endif

end subroutine write_par_cdset_zdf

! --------------------------------------------------------------------------------------------------
! Write end marker for a chunked dataset
! --------------------------------------------------------------------------------------------------
subroutine end_cdset_zdf( this, dset )

  implicit none

  class( t_diag_file_zdf ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset

  type( t_zdf_dataset ) :: zdf_dset

  call diag_to_zdf( dset, zdf_dset )

  if ( this % parallel ) then
    call zdf_end_cdset( this % parFile, zdf_dset )
  else
    call zdf_end_cdset( this % file, zdf_dset )
  endif

  call zdf_to_diag( zdf_dset, dset )

end subroutine end_cdset_zdf

! --------------------------------------------------------------------------------------------------
! Open a chunked dataset for updating
! --------------------------------------------------------------------------------------------------
subroutine open_cdset_zdf( this, dset )

  implicit none

  class( t_diag_file_zdf ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset

  type( t_zdf_dataset ) :: zdf_dset

  call diag_to_zdf( dset, zdf_dset )
  call zdf_open_cdset( this % file, zdf_dset )
  call zdf_to_diag( zdf_dset, dset )

end subroutine open_cdset_zdf

! --------------------------------------------------------------------------------------------------
! Extend dimensions of chunked dataset
! --------------------------------------------------------------------------------------------------
subroutine extend_cdset_zdf( this, dset, dims )

  implicit none

  class( t_diag_file_zdf ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset
  integer(p_int64), dimension(:), intent(in) :: dims

  type( t_zdf_dataset ) :: zdf_dset

  call diag_to_zdf( dset, zdf_dset )
  if ( this % parallel ) then
    call zdf_extend_cdset( this % parFile, zdf_dset, dims )
  else
    call zdf_extend_cdset( this % file, zdf_dset, dims )
  endif

  call zdf_to_diag( zdf_dset, dset )

end subroutine extend_cdset_zdf


! --------------------------------------------------------------------------------------------------
! Close previously open chunked dataset
! --------------------------------------------------------------------------------------------------
subroutine close_cdset_zdf( this, dset )

  implicit none

  class( t_diag_file_zdf ), intent(inout) :: this
  type( t_diag_dataset ), intent(inout) :: dset

  ! For ZDF files we don't need to do anything
  dset % ndims = 0
  dset % id = 0
  dset % offset = -1

end subroutine close_cdset_zdf

! --------------------------------------------------------------------------------------------------
! add_dataset routines
! --------------------------------------------------------------------------------------------------

#define __TEMPLATE__

#define __TYPE__ real(p_single)
#define __DIMS__ :
#define __NDIMS__ 1
#define __ZDF_TYPE__ zdf_float32
#define FNAME(a) a ## _1d_r4_zdf
#include __FILE__

#define __TYPE__ real(p_single)
#define __DIMS__ :,:
#define __NDIMS__ 2
#define __ZDF_TYPE__ zdf_float32
#define FNAME(a) a ## _2d_r4_zdf
#include __FILE__

#define __TYPE__ real(p_single)
#define __DIMS__ :,:,:
#define __NDIMS__ 3
#define __ZDF_TYPE__ zdf_float32
#define FNAME(a) a ## _3d_r4_zdf
#include __FILE__

#define __TYPE__ real(p_double)
#define __DIMS__ :
#define __NDIMS__ 1
#define __ZDF_TYPE__ zdf_float64
#define FNAME(a) a ## _1d_r8_zdf
#include __FILE__

#define __TYPE__ real(p_double)
#define __DIMS__ :,:
#define __NDIMS__ 2
#define __ZDF_TYPE__ zdf_float64
#define FNAME(a) a ## _2d_r8_zdf
#include __FILE__

#define __TYPE__ real(p_double)
#define __DIMS__ :,:,:
#define __NDIMS__ 3
#define __ZDF_TYPE__ zdf_float64
#define FNAME(a) a ## _3d_r8_zdf
#include __FILE__

end module m_diagfile_zdf

#else

!---------------------------------------------------------------------------------------------------
!  Template Function definitions
!---------------------------------------------------------------------------------------------------

subroutine FNAME(add_dataset)( this, name, buffer )

  use zdf_parallel

  implicit none

  class( t_diag_file_zdf ), intent(inout) :: this
  character(len=*), intent(in) :: name
  __TYPE__, dimension(__DIMS__), intent(inout), target :: buffer

  type(t_zdf_dataset) :: dset
  integer :: i

  dset % name = name
  dset % data_type = diag_zdf_type(__ZDF_TYPE__)
  dset % ndims     = __NDIMS__
  do i = 1, __NDIMS__
    dset % count(i) = size(buffer, i)
  enddo
  dset % data = c_loc( buffer )

  call zdf_add( this%file, dset )

end subroutine FNAME(add_dataset)


#undef __TYPE__
#undef __DIMS__
#undef __NDIMS__
#undef __ZDF_TYPE__
#undef FNAME

#endif
