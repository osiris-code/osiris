
#include "os-config.h"
#include "os-preprocess.fpp"

module m_diagnostic_utilities

#include "memory/memory.h"

use m_parameters

use m_diagfile
use m_diagfile_zdf

#ifdef HDF5
use m_diagfile_hdf5
#endif

implicit none

private

! Precision for diagnostics
! some diagnostics (e.g. field grids) can be controlled in the input file using the prec
! parameter

integer, parameter :: p_diag_prec = p_single
!integer, parameter :: p_diag_prec = p_double

! Size of buffer for diagnostics (to be used by phasespaces & averaged grids)
#if defined( __MIC__ )

#warning Setting smaller diagnostics buffer for MIC platform

! on Intel MIC use a smaller buffer
! (this sets it to 1 MB in single precision)
integer, parameter :: p_max_diag_buffer_size = 1024*1024/4

#elif defined( __KNL__ )

#warning Setting smaller diagnostics buffer for KNL platform

! on Intel MIC use a smaller buffer
! (this sets it to 2 MB in single precision)
integer, parameter :: p_max_diag_buffer_size = 1024*1024*2/4

#else

! (this sets it to 64 MB in single precision)
integer, parameter :: p_max_diag_buffer_size = 1024*1024*(64/4)

#endif

integer :: diag_buffer_size = -1
real(p_diag_prec), dimension(:), pointer :: diag_buffer


integer, parameter :: p_zdf_format  = 0
integer, parameter :: p_hdf5_format = 1

! Default file format to use
integer, private :: file_format_default

interface get_filename
  module procedure get_filename
end interface

interface init_dutil
  module procedure init_diagutil
end interface

interface cleanup_dutil
  module procedure cleanup_diagutil
end interface

interface get_diag_buffer
  module procedure get_diag_buffer
end interface

interface diag_type_real
  module procedure diag_type_real
end interface

interface create_diag_file
  module procedure create_diag_file
end interface

public :: t_diag_file, create_diag_file
public :: p_diag_grid, p_diag_particles, p_diag_tracks
public :: get_filename

public :: init_dutil, cleanup_dutil, get_diag_buffer

public :: p_diag_prec, diag_type_real, p_max_diag_buffer_size

! pass though some parameters/symbols from os-diagfile for convenience 
public :: diag_axis_linear, diag_axis_log10, diag_axis_log2
public :: t_diag_dataset, t_diag_chunk
public :: p_diag_create, p_diag_read, p_diag_update
public :: p_diag_float32, p_diag_float64, p_diag_int32, p_diag_int64, p_diag_null


public :: p_hdf5_format, p_zdf_format
public :: DIAG_MPI, DIAG_INDEPENDENT, DIAG_MPIIO_INDEPENDENT, DIAG_MPIIO_COLLECTIVE, DIAG_POSIX
public :: sim_info

contains

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
pure function diag_type_real( kind )

  implicit none

  integer :: diag_type_real
  integer, intent(in) :: kind

  select case (kind)
  case (p_double)
    diag_type_real = p_diag_float64
  case (p_single)
    diag_type_real = p_diag_float32
  case default
    diag_type_real = -1
  end select

end function diag_type_real

!--------------------------------------------------------------------------------------------------
! Initialize module variables
!--------------------------------------------------------------------------------------------------
subroutine init_diagutil( sim_options, buffer_size )

#ifdef HDF5
  use hdf5_util
#endif

  implicit none

  type( t_options ), intent(in) :: sim_options
  integer, intent(in) :: buffer_size

#ifdef HDF5
  ! Initialize hdf5 subsystem
  call open_hdf5( )
#endif

  ! Set default file format
  file_format_default = sim_options % file_format

  ! Set default parallel io algorithm
  call diag_file_set_def_iomode( sim_options % parallel_io )

  ! allocate diagnostics buffer
  if (buffer_size > p_max_diag_buffer_size) then
    diag_buffer_size = p_max_diag_buffer_size
  else
    diag_buffer_size = buffer_size
  endif

  if ( diag_buffer_size > 0 ) then
    call alloc( diag_buffer, (/diag_buffer_size/) )
  endif

end subroutine init_diagutil
!--------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
subroutine cleanup_diagutil( )
!--------------------------------------------------------------------------------------------------
! Cleanup module
!--------------------------------------------------------------------------------------------------
  implicit none

  if ( diag_buffer_size > 0 ) then
    call freemem( diag_buffer )
  endif

  call freemem( sim_info % nx_x1_node )
  call freemem( sim_info % nx_x2_node )
  call freemem( sim_info % nx_x3_node )

end subroutine cleanup_diagutil
!--------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
subroutine get_diag_buffer( buffer, bsize )
!--------------------------------------------------------------------------------------------------
  implicit none

  real(p_diag_prec), dimension(:), pointer :: buffer
  integer, intent(out), optional :: bsize

  buffer => diag_buffer
  if ( present(bsize) ) then
    bsize = diag_buffer_size
  endif

end subroutine get_diag_buffer
!--------------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------------
! Generate a file name from the information provided
!--------------------------------------------------------------------------------------------------
function get_filename( n, prefix, suffix, node )

  use stringutil

  implicit none

  ! arguments and return value

  character(len = 80) :: get_filename
  integer, intent(in) :: n
  character(len = *), intent(in) :: prefix

  character(len = *), intent(in), optional :: suffix
  integer,    intent(in), optional :: node

  get_filename = trim(prefix)//'-'//trim(idx_string(n,p_time_length))

  if (present(node)) get_filename = trim(get_filename)// '.'//trim(idx_string(node,3))

  if (present(suffix)) get_filename = trim(get_filename)//trim(suffix)

end function get_filename
!--------------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------------
! Return new diagfile
!--------------------------------------------------------------------------------------------------
subroutine create_diag_file(diagFile, type )

  implicit none

  class( t_diag_file ), allocatable :: diagFile
  integer, intent(in), optional :: type
  
  integer :: ltype

  if ( present(type) ) then
    ltype = type
  else
    ltype = file_format_default
  endif
  
  if ( allocated(diagFile) ) deallocate(diagFile)
  
  select case( ltype )
  case( p_zdf_format )
    allocate( t_diag_file_zdf  :: diagFile )

#ifdef HDF5
  case( p_hdf5_format )
    allocate( t_diag_file_hdf5 :: diagFile )
#endif

  case default
    ERROR("Invalid diagnostics file format requested, aborting.")
    call abort_program( p_err_invalid )
  end select

end subroutine
!--------------------------------------------------------------------------------------------------

end module m_diagnostic_utilities
