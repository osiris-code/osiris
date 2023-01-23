#include "os-config.h"
#include "os-preprocess.fpp"

module m_diagfile
    
use iso_c_binding
use m_system, only: t_sim_info

implicit none

private

! abstracted datatypes.
enum, bind(c)
    enumerator :: diag_null, diag_int8, diag_uint8, diag_int16, diag_uint16, diag_int32, &
    diag_uint32, diag_int64, diag_uint64, diag_float32, diag_float64
endenum

public :: diag_null, diag_int8, diag_uint8, diag_int16, diag_uint16, diag_int32, &
diag_uint32, diag_int64, diag_uint64, diag_float32, diag_float64

! constants to specifiy I/O types
enum, bind(c)
    enumerator :: DIAG_MPI, DIAG_INDEPENDENT, DIAG_MPIIO_INDEPENDENT, DIAG_MPIIO_COLLECTIVE, DIAG_POSIX
endenum

! some general Fortran type info
integer, private, parameter :: p_single = kind(1.0e0)
integer, private, parameter :: p_double = kind(1.0d0)
integer, private, parameter :: p_byte   = selected_int_kind(2)
integer, private, parameter :: p_int64  = selected_int_kind(10)

! General simulation info to be added to diagnostic files
type( t_sim_info ), save, public :: sim_info

integer, parameter :: p_diag_grid      = 1
integer, parameter :: p_diag_particles = 2
integer, parameter :: p_diag_tracks    = 3

! Needs to be this large for tracking all field quantities with 15 cyl_modes
integer, parameter, private :: p_max_quants = 114

integer, parameter :: p_diag_float32 = diag_float32
integer, parameter :: p_diag_float64 = diag_float64
integer, parameter :: p_diag_int32   = diag_int32
integer, parameter :: p_diag_int64   = diag_int64
integer, parameter :: p_diag_null    = diag_null


integer, parameter :: p_diag_create = 0
integer, parameter :: p_diag_read   = 1
integer, parameter :: p_diag_update = 2


integer, parameter :: diag_axis_linear = 0
integer, parameter :: diag_axis_log10  = 1
integer, parameter :: diag_axis_log2   = 2

integer, parameter :: diag_max_dims = 3

integer, parameter :: p_str_maxlen = 256

! Number of digits for time step in diag files
integer, parameter :: p_time_length = 6

real(p_double), parameter :: p_no_offset = - huge(1.0_p_double)

! Default parallel io mode
integer :: diag_file_default_iomode = DIAG_MPI


type :: t_diag_dataset
    character(len = p_str_maxlen ) :: name = ''
    integer(c_int32_t) :: data_type = 0
    integer(c_int32_t) :: ndims = 0
    integer(c_int64_t), dimension(diag_max_dims) :: count
    type(c_ptr) :: data
    integer(c_int64_t) :: id = -1
    integer(c_int64_t) :: offset = 0
end type


type :: t_diag_grid_axis
    integer :: type = diag_axis_linear
    real(p_double) :: min = 0.0
    real(p_double) :: max = 1.0
    character(len = p_str_maxlen ) :: name = ''
    character(len = p_str_maxlen ) :: label = ''
    character(len = p_str_maxlen ) :: units = ''
end type

type :: t_diag_iteration
    character(len = p_str_maxlen ) :: name = 'ITERATION'
    integer :: n
    real(p_double) :: t
    character(len = p_str_maxlen ) :: time_units
end type

type :: t_diag_grid_info
    integer :: ndims = -1
    integer, dimension(diag_max_dims) :: count = 0
    character(len = p_str_maxlen ) :: name = ''
    character(len = p_str_maxlen ) :: label = ''
    character(len = p_str_maxlen ) :: units = ''
    type(t_diag_grid_axis), dimension(diag_max_dims) :: axis
    real(p_double), dimension(diag_max_dims) :: offset_x = p_no_offset
    real(p_double) :: offset_t = p_no_offset
    ! This offset is for phasespace axes
    real(p_double), dimension(diag_max_dims) :: offset_t_ax = p_no_offset
end type

type :: t_diag_part_info
    character(len = p_str_maxlen ) :: name = ''
    character(len = p_str_maxlen ) :: label = ''
    integer(p_int64) :: np
    integer :: nquants
    character(len = p_str_maxlen ), dimension(p_max_quants) :: quants = ''
    character(len = p_str_maxlen ), dimension(p_max_quants) :: qlabels = ''
    character(len = p_str_maxlen ), dimension(p_max_quants) :: qunits = ''
    real(p_double), dimension(p_max_quants) :: offset_t = p_no_offset
    logical :: has_tags
end type

type, bind(C) :: t_diag_chunk
    integer(c_int64_t), dimension(diag_max_dims) :: count
    integer(c_int64_t), dimension(diag_max_dims) :: start
    integer(c_int64_t), dimension(diag_max_dims) :: stride
    type(c_ptr) :: data
end type

type :: t_diag_track_info
    character(len = p_str_maxlen ) :: name = ''
    character(len = p_str_maxlen ) :: label = ''
    integer :: ntracks
    integer :: ndump
    integer :: niter
    integer :: nquants
    character(len = p_str_maxlen ), dimension(p_max_quants) :: quants = ''
    character(len = p_str_maxlen ), dimension(p_max_quants) :: qlabels = ''
    character(len = p_str_maxlen ), dimension(p_max_quants) :: qunits = ''
    real(p_double), dimension(p_max_quants) :: offset_t = p_no_offset
end type


!----------------------------------------------------------------------------
! Diagnostic File abstract defintion
!----------------------------------------------------------------------------
type, abstract :: t_diag_file

    ! type of file (grid / particles / tracks)
    integer :: ftype

    ! file handles for serial and parallel formats
    ! type(t_diag_file_c)     :: file
    ! type(t_diag_par_file_c) :: parFile
    logical                 :: parallel

    ! name
    character(len=256) :: name

    ! file name and path
    character( len = 1024 ) :: filename = '', filepath = ''

    ! time information (not used for tracks)
    type( t_diag_iteration ) :: iter

    ! Grid information
    type( t_diag_grid_info ) :: grid

    ! Particle information
    type( t_diag_part_info ) :: particles

    ! Tracks information
    type( t_diag_track_info ) :: tracks

    ! selection (particles)
    real(p_double)     :: raw_gamma_limit = 0.0_p_double
    real(p_double)     :: raw_fraction = 1.0_p_double
    character(len=1024):: raw_math_expr = ''

contains

    procedure(open_diag_file), deferred     :: open_diag_file
    procedure(open_par_diag_file), deferred :: open_par_diag_file
    generic   :: open  => open_diag_file, open_par_diag_file

    procedure(close_file), deferred        :: close

    procedure(add_dataset_1D_r4), deferred :: add_dataset_1D_r4
    procedure(add_dataset_2D_r4), deferred :: add_dataset_2D_r4
    procedure(add_dataset_3D_r4), deferred :: add_dataset_3D_r4

    procedure(add_dataset_1D_r8), deferred :: add_dataset_1D_r8
    procedure(add_dataset_2D_r8), deferred :: add_dataset_2D_r8
    procedure(add_dataset_3D_r8), deferred :: add_dataset_3D_r8

    ! Gfortran throws the strangest error in 'os-paritcles-define.f90' if a enable this method.
    !   for disable for now since it is unused in the rest of the code
    !procedure(add_dataset_loc_), deferred :: add_dataset_loc_

    generic :: add_dataset => add_dataset_1D_r4, add_dataset_2D_r4, add_dataset_3D_r4, &
        add_dataset_1D_r8, add_dataset_2D_r8, add_dataset_3D_r8 !, add_dataset_loc_

    procedure(start_cdset), deferred :: start_cdset

    procedure(write_cdset), deferred :: write_cdset
    procedure(write_par_cdset), deferred :: write_par_cdset

    procedure(end_cdset), deferred :: end_cdset

    procedure(open_cdset), deferred :: open_cdset
    procedure(extend_cdset), deferred :: extend_cdset
    procedure(close_cdset), deferred :: close_cdset

    procedure(file_extension), deferred :: file_extension

end type t_diag_file


abstract interface
    subroutine open_diag_file( this, amode )
        import t_diag_file
        class( t_diag_file ), intent(inout) :: this
        integer, intent(in) :: amode
    end subroutine
end interface

abstract interface
    subroutine open_par_diag_file( this, amode, comm, iomode )
        import t_diag_file
        class( t_diag_file ), intent(inout) :: this
        integer, intent(in) :: amode
        integer, intent(in) :: comm
        integer, intent(in), optional :: iomode
    end subroutine
end interface

abstract interface
subroutine close_file( this )
    import t_diag_file 
        class( t_diag_file ), intent(inout) :: this
    end subroutine
end interface

abstract interface
    subroutine add_dataset_1D_r4(this, name, buffer )
        import t_diag_file, p_single
        class( t_diag_file ), intent(inout) :: this
        character(len=*), intent(in) :: name
        real(p_single), dimension(:), intent(inout), target :: buffer
    end subroutine
end interface

abstract interface
    subroutine add_dataset_2D_r4(this, name, buffer )
        import t_diag_file, p_single
        class( t_diag_file ), intent(inout) :: this
        character(len=*), intent(in) :: name
        real(p_single), dimension(:,:), intent(inout), target :: buffer
    end subroutine
end interface

abstract interface
    subroutine add_dataset_3D_r4(this, name, buffer )
        import t_diag_file, p_single
        class( t_diag_file ), intent(inout) :: this
        character(len=*), intent(in) :: name
        real(p_single), dimension(:,:,:), intent(inout), target :: buffer
    end subroutine
end interface


abstract interface
    subroutine add_dataset_1D_r8(this, name, buffer )
        import t_diag_file, p_double
        class( t_diag_file ), intent(inout) :: this
        character(len=*), intent(in) :: name
        real(p_double), dimension(:), intent(inout), target :: buffer
    end subroutine
end interface

abstract interface
    subroutine add_dataset_2D_r8(this, name, buffer )
        import t_diag_file, p_double
        class( t_diag_file ), intent(inout) :: this
        character(len=*), intent(in) :: name
        real(p_double), dimension(:,:), intent(inout), target :: buffer
    end subroutine
end interface

abstract interface
    subroutine add_dataset_3D_r8(this, name, buffer )
        import t_diag_file, p_double
        class( t_diag_file ), intent(inout) :: this
        character(len=*), intent(in) :: name
        real(p_double), dimension(:,:,:), intent(inout), target :: buffer
    end subroutine
end interface


abstract interface
    subroutine add_dataset_loc_( this, name, buffer, ndims, count, data_type )
        import t_diag_file, c_ptr
        class( t_diag_file ), intent(inout) :: this
        character(len=*), intent(in) :: name
        type(c_ptr), intent(in) :: buffer
        integer, intent(in) :: data_type, ndims
        integer, dimension(:), intent(in) :: count
    end subroutine
end interface

abstract interface
    subroutine start_cdset( this, name, ndims, count, data_type, dset, chunk_size )
        import t_diag_file, t_diag_dataset
        class( t_diag_file ), intent(inout) :: this
        character(len = *), intent(in) :: name
        integer, intent(in) :: ndims
        integer, dimension(:), intent(in) :: count
        integer, intent(in) :: data_type
        type( t_diag_dataset ), intent(inout) :: dset

        integer, dimension(:), intent(in), optional :: chunk_size
    end subroutine
end interface

abstract interface
    subroutine write_cdset( this, dset, chunk )
        import t_diag_file, t_diag_dataset, t_diag_chunk
        class( t_diag_file ), intent(inout) :: this
        type( t_diag_dataset ), intent(inout) :: dset
        type( t_diag_chunk ), intent(inout) :: chunk
    end subroutine
end interface

abstract interface
    subroutine write_par_cdset( this, dset, chunk, offset_ )
        import t_diag_file, t_diag_dataset, t_diag_chunk
        class( t_diag_file ), intent(inout) :: this
        type( t_diag_dataset ), intent(inout) :: dset
        type( t_diag_chunk ), intent(inout) :: chunk
        ! Offset is not used by HDF5 routines
        integer, intent(in), optional :: offset_
    end subroutine
end interface

abstract interface
    subroutine end_cdset( this, dset )
        import t_diag_file, t_diag_dataset
        class( t_diag_file ), intent(inout) :: this
        type( t_diag_dataset ), intent(inout) :: dset
    end subroutine
end interface

abstract interface
    subroutine open_cdset( this, dset )
        import t_diag_file, t_diag_dataset
        class( t_diag_file ), intent(inout) :: this
        type( t_diag_dataset ), intent(inout) :: dset
    end subroutine
end interface

abstract interface
    subroutine extend_cdset( this, dset, dims )
        import t_diag_file, t_diag_dataset, p_int64
        class( t_diag_file ), intent(inout) :: this
        type( t_diag_dataset ), intent(inout) :: dset
        integer(p_int64), dimension(:), intent(in) :: dims
    end subroutine
end interface

abstract interface
    subroutine close_cdset( this, dset )
        import t_diag_file, t_diag_dataset
        class( t_diag_file ),   intent(inout) :: this
        type( t_diag_dataset ), intent(inout) :: dset
    end subroutine
end interface

abstract interface
    pure function file_extension( this )
        import t_diag_file
        class( t_diag_file ),   intent(in) :: this
        character(len=3) :: file_extension
    end function
end interface

! Fortran interfaces
interface diag_sizeof
    module procedure diag_sizeof
end interface

interface freal_to_diagtype
    module procedure freal_to_diagtype
end interface

interface diag_mpi_type
    module procedure diag_mpi_type
end interface

interface diag_iomode_name
    module procedure diag_iomode_name
end interface

interface diag_file_set_def_iomode
    module procedure diag_file_set_def_iomode
end interface


public :: DIAG_MPI, DIAG_INDEPENDENT, DIAG_MPIIO_INDEPENDENT, DIAG_MPIIO_COLLECTIVE, DIAG_POSIX
public :: diag_sizeof, freal_to_diagtype, diag_mpi_type, diag_iomode_name
public :: diag_file_set_def_iomode
public :: t_diag_file
public :: p_diag_grid, p_diag_particles, p_diag_tracks
public :: p_diag_float32, p_diag_float64, p_diag_int32, p_diag_int64, p_diag_null
public :: p_diag_create, p_diag_read, p_diag_update
public :: diag_axis_linear, diag_axis_log10, diag_axis_log2
public :: diag_max_dims
public :: p_str_maxlen, p_no_offset, p_time_length
public :: diag_file_default_iomode

! public :: t_diag_file_c, t_diag_par_file_c
public :: t_diag_iteration, t_diag_grid_info, t_diag_part_info, t_diag_grid_axis
public :: t_diag_dataset, t_diag_chunk, t_diag_track_info

contains


subroutine diag_file_set_def_iomode( iomode )
    
    implicit none
    
    integer, intent(in) :: iomode
    
    diag_file_default_iomode = iomode
    
end subroutine


function diag_sizeof( data_type )
    
    implicit none
    
    integer, intent(in) :: data_type
    integer :: diag_sizeof
    
    select case ( data_type )
    case( diag_null )
        diag_sizeof = 0
    case( diag_int8 )
        diag_sizeof = 1
    case( diag_uint8 )
        diag_sizeof = 1
    case( diag_int16 )
        diag_sizeof = 2
    case( diag_uint16 )
        diag_sizeof = 2
    case( diag_int32 )
        diag_sizeof = 4
    case( diag_uint32 )
        diag_sizeof = 4
    case( diag_float32 )
        diag_sizeof = 4
    case( diag_int64 )
        diag_sizeof = 8
    case( diag_uint64 )
        diag_sizeof = 8
    case( diag_float64 )
        diag_sizeof = 8
    case default
        diag_sizeof = 0
    end select
    
end function

function freal_to_diagtype( rkind )
    
    implicit none
    
    integer, intent(in) :: rkind
    integer :: freal_to_diagtype
    
    select case(rkind)
    case(p_single)
        freal_to_diagtype = diag_float32
    case(p_double)
        freal_to_diagtype = diag_float64
    case default
        freal_to_diagtype = diag_null
    end select
    
end function freal_to_diagtype


function diag_iomode_name( iomode )
    
    implicit none
    
    integer, intent(in) :: iomode
    character(len=18) :: diag_iomode_name
    
    select case(iomode)
    case (DIAG_MPI)
        diag_iomode_name = "mpi"
    case (DIAG_INDEPENDENT)
        diag_iomode_name = "independent"
    case (DIAG_MPIIO_INDEPENDENT)
        diag_iomode_name = "mpi-io independent"
    case (DIAG_MPIIO_COLLECTIVE)
        diag_iomode_name = "mpi-io collective"
    case default
        diag_iomode_name = "invalid"
    end select
    
end function diag_iomode_name

function diag_mpi_type( data_type )
    
    use mpi
    
    implicit none
    
    integer, intent(in) :: data_type
    integer :: diag_mpi_type
    
    select case ( data_type )
    case( diag_null )
        diag_mpi_type = MPI_DATATYPE_NULL
    case( diag_int8 )
        diag_mpi_type = MPI_BYTE
    case( diag_uint8 )
        diag_mpi_type = MPI_BYTE
    case( diag_int16 )
        diag_mpi_type = MPI_INTEGER2
    case( diag_uint16 )
        diag_mpi_type = MPI_INTEGER2
    case( diag_int32 )
        diag_mpi_type = MPI_INTEGER4
    case( diag_uint32 )
        diag_mpi_type = MPI_INTEGER4
    case( diag_float32 )
        diag_mpi_type = MPI_REAL
    case( diag_int64 )
        diag_mpi_type = MPI_INTEGER8
    case( diag_uint64 )
        diag_mpi_type = MPI_INTEGER8
    case( diag_float64 )
        diag_mpi_type = MPI_DOUBLE_PRECISION
    case default
        diag_mpi_type = MPI_DATATYPE_NULL
    end select
    
end function

end module m_diagfile