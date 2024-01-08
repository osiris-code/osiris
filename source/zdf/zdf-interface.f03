! Fortran90 interface routines for ZDF
module zdf

use iso_c_binding

integer, private, parameter :: p_single = kind(1.0e0)
integer, private, parameter :: p_double = kind(1.0d0)
integer, private, parameter :: p_byte   = selected_int_kind(2)
integer, private, parameter :: p_int64  = selected_int_kind(10)

integer, parameter :: zdf_str_maxlen = 256

integer, parameter :: zdf_max_dims = 3

integer, parameter :: zdf_file_create = 0
integer, parameter :: zdf_file_read   = 1
integer, parameter :: zdf_file_update = 2

integer, parameter :: zdf_axis_linear = 0
integer, parameter :: zdf_axis_log10  = 1
integer, parameter :: zdf_axis_log2   = 2

enum, bind(c)
    enumerator :: zdf_null, zdf_int8, zdf_uint8, zdf_int16, zdf_uint16, zdf_int32, &
                   zdf_uint32, zdf_int64, zdf_uint64, zdf_float32, zdf_float64
endenum

! Needs to be this large for tracking all field quantities with 15 cyl_modes
integer, private, parameter :: p_max_quants = 114

type, bind(C) :: t_zdf_file
    type(c_ptr) :: fp = C_NULL_PTR
    integer(c_int) :: mode
    integer(c_int32_t) :: ndatasets
end type

type :: t_zdf_dataset
    character(len = zdf_str_maxlen ) :: name = 'DATASET'
    integer(c_int32_t) :: data_type = 0
    integer(c_int32_t) :: ndims = 0
    integer(c_int64_t), dimension(zdf_max_dims) :: count
    type(c_ptr) :: data = C_NULL_PTR
    integer(c_int64_t) :: id = -1
    integer(c_int64_t) :: offset = 0
end type

type :: t_zdf_grid_axis
    integer :: type = zdf_axis_linear
    real(p_double) :: min = 0.0
    real(p_double) :: max = 1.0
    character(len = zdf_str_maxlen ) :: name = ''
    character(len = zdf_str_maxlen ) :: label = ''
    character(len = zdf_str_maxlen ) :: units = ''
end type

type :: t_zdf_iteration
    character(len = zdf_str_maxlen ) :: name = 'ITERATION'
    integer :: n
    real(p_double) :: t
    character(len = zdf_str_maxlen ) :: time_units = ''
end type

type :: t_zdf_grid_info
    integer :: ndims = -1
    integer, dimension(zdf_max_dims) :: count = 0
    character(len = zdf_str_maxlen ) :: name = 'GRID'
    character(len = zdf_str_maxlen ) :: label = ''
    character(len = zdf_str_maxlen ) :: units = ''
    type(t_zdf_grid_axis), dimension(zdf_max_dims) :: axis
end type

type :: t_zdf_part_info
    character(len = zdf_str_maxlen ) :: name = 'PARTICLES'
    character(len = zdf_str_maxlen ) :: label = ''
    integer(p_int64) :: np
    integer :: nquants
    character(len = zdf_str_maxlen ), dimension(p_max_quants) :: quants  = ''
    character(len = zdf_str_maxlen ), dimension(p_max_quants) :: qlabels = ''
    character(len = zdf_str_maxlen ), dimension(p_max_quants) :: qunits  = ''
    logical :: has_tags
end type

type, bind(C) :: t_zdf_chunk
    integer(c_int64_t), dimension(zdf_max_dims) :: count
    integer(c_int64_t), dimension(zdf_max_dims) :: start
    integer(c_int64_t), dimension(zdf_max_dims) :: stride
    type(c_ptr) :: data
end type

type :: t_zdf_track_info
    character(len = zdf_str_maxlen ) :: name = "TRACKS"
    character(len = zdf_str_maxlen ) :: label = ""
    integer :: ntracks
    integer :: ndump
    integer :: niter
    integer :: nquants
    character(len = zdf_str_maxlen ), dimension(p_max_quants) :: quants  = ''
    character(len = zdf_str_maxlen ), dimension(p_max_quants) :: qlabels = ''
    character(len = zdf_str_maxlen ), dimension(p_max_quants) :: qunits  = ''
end type

! C types to be used in C function calls
type, bind(C) :: t_zdf_dataset_c
    type(c_ptr) :: name = C_NULL_PTR
    integer(c_int32_t) :: data_type = 0
    integer(c_int32_t) :: ndims = 0
    integer(c_int64_t), dimension(zdf_max_dims) :: count
    type(c_ptr) :: data = C_NULL_PTR
    integer(c_int64_t) :: id = -1
    integer(c_int64_t) :: offset = 0
end type

type, bind(C) :: t_zdf_iteration_c
    type(c_ptr) :: name = C_NULL_PTR
    integer(c_int32_t) :: n
    real(c_double) :: t
    type(c_ptr) :: time_units = C_NULL_PTR
end type

type, bind(C) :: t_zdf_grid_axis_c
    type(c_ptr) :: name = C_NULL_PTR
    integer(c_int) :: type
    real(c_double) :: min, max
    type(c_ptr) :: label = C_NULL_PTR
    type(c_ptr) :: units = C_NULL_PTR
end type

type, bind(C) :: t_zdf_grid_info_c
    type(c_ptr) :: name = C_NULL_PTR
    integer(c_int32_t) :: ndims
    integer(c_int64_t), dimension(zdf_max_dims) :: count
    type(c_ptr) :: label = C_NULL_PTR
    type(c_ptr) :: units = C_NULL_PTR
    type(c_ptr) :: axis = C_NULL_PTR
end type

type, bind(C) :: t_zdf_part_info_c
    type(c_ptr) :: name = C_NULL_PTR
    type(c_ptr) :: label = C_NULL_PTR
    integer(c_int64_t) :: np
    integer(c_int32_t) :: nquants
    type(c_ptr) :: quants = C_NULL_PTR
    type(c_ptr) :: qlabels = C_NULL_PTR
    type(c_ptr) :: qunits = C_NULL_PTR
end type

type, bind(C) :: t_zdf_track_info_c
    type(c_ptr) :: name = C_NULL_PTR
    type(c_ptr) :: label = C_NULL_PTR
    integer(c_int32_t) :: ntracks
    integer(c_int32_t) :: ndump
    integer(c_int32_t) :: niter
    integer(c_int32_t) :: nquants
    type(c_ptr) :: quants = C_NULL_PTR
    type(c_ptr) :: qlabels = C_NULL_PTR
    type(c_ptr) :: qunits = C_NULL_PTR
end type

! Fortran interfaces
interface zdf_sizeof
    module procedure zdf_sizeof
end interface

interface zdf_open_file
    module procedure zdf_open_file
end interface

interface zdf_close_file
    module procedure zdf_close_file
end interface

interface zdf_open_grid_file
    module procedure zdf_open_grid_file
end interface

interface zdf_open_part_file
    module procedure zdf_open_part_file
end interface

interface zdf_add
    module procedure zdf_add_string
    module procedure zdf_add_int32
    module procedure zdf_add_r8
    module procedure zdf_add_iteration
    module procedure zdf_add_grid_info
    module procedure zdf_add_part_info
    module procedure zdf_add_track_info
    module procedure zdf_add_dataset
end interface

! Chunked dataset Interface
interface zdf_start_cdset
    module procedure zdf_start_cdset
end interface

interface zdf_write_cdset
    module procedure zdf_write_cdset
end interface

interface zdf_end_cdset
    module procedure zdf_end_cdset
end interface

interface zdf_open_cdset
    module procedure zdf_open_cdset
end interface

interface zdf_extend_cdset
    module procedure zdf_extend_cdset
end interface

! High level interfaces
interface zdf_save_grid
    module procedure zdf_save_grid_1D_r4
    module procedure zdf_save_grid_2D_r4
    module procedure zdf_save_grid_3D_r4
end interface


interface zdf_add_quant_part_file
    module procedure zdf_add_quant_part_file
end interface

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_open_file( this, filename, access_mode )

    implicit none

    interface
        subroutine zdf_open_file_c( handle, filename, zdf_file_access ) bind(C, name="zdf_open_file")
            import t_zdf_file, c_char, c_int
            type(t_zdf_file)       :: handle
            character(kind=c_char) :: filename(*)
            integer(c_int), value  :: zdf_file_access
        end subroutine zdf_open_file_c
    end interface

    type( t_zdf_file ), intent(inout) :: this
    character(len = *), intent(in)  :: filename
    integer, intent(in)             :: access_mode

    call zdf_open_file_c( this, trim(filename) // C_NULL_CHAR, access_mode )

end subroutine zdf_open_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_close_file( this )

    implicit none

    interface
        subroutine zdf_close_file_c( handle ) bind(C, name="zdf_close_file")
            import t_zdf_file
            type(t_zdf_file)     :: handle
        end subroutine zdf_close_file_c
    end interface

    type( t_zdf_file ), intent(inout) :: this

    call zdf_close_file_c( this )

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_add_string( file, name, str )

    implicit none

    type(t_zdf_file), intent(in) :: file
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: str

    interface
        subroutine zdf_add_string_c( file, name, str ) bind(C, name="zdf_add_string")
            import t_zdf_file, c_char
            type(t_zdf_file)       :: file
            character(kind=c_char) :: name(*)
            character(kind=c_char) :: str(*)
        end subroutine zdf_add_string_c
    end interface

    call zdf_add_string_c( file, trim(name) // C_NULL_CHAR, str // C_NULL_CHAR )

end subroutine zdf_add_string

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_add_int32( file, name, val )

    implicit none

    type(t_zdf_file), intent(in) :: file
    character(len=*), intent(in) :: name
    integer, intent(in) :: val

    interface
        subroutine zdf_add_int32_c( file, name, val ) bind(C, name="zdf_add_int32")
            import t_zdf_file, c_char, c_int32_t
            type(t_zdf_file)       :: file
            character(kind=c_char) :: name(*)
            integer(c_int32_t), value  :: val
        end subroutine zdf_add_int32_c
    end interface

    integer(c_int32_t) :: val_c

    val_c = val
    call zdf_add_int32_c( file, trim(name) // C_NULL_CHAR, val_c )

end subroutine zdf_add_int32

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_add_r8( file, name, val )

    implicit none

    type(t_zdf_file), intent(in) :: file
    character(len=*), intent(in) :: name
    real(p_double), intent(in) :: val

    interface
        subroutine zdf_add_double_c( file, name, val ) bind(C, name="zdf_add_double")
            import t_zdf_file, c_char, c_double
            type(t_zdf_file)       :: file
            character(kind=c_char) :: name(*)
            real(c_double), value  :: val
        end subroutine zdf_add_double_c
    end interface

    real(c_double) :: val_c

    val_c = val
    call zdf_add_double_c( file, trim(name) // C_NULL_CHAR, val_c )

end subroutine zdf_add_r8

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_add_iteration( file, iter )

    implicit none

    type(t_zdf_file), intent(in) :: file
    type(t_zdf_iteration), intent(in) :: iter

    interface
        subroutine zdf_add_iteration_c( file, iter ) bind(C, name="zdf_add_iteration")
            import t_zdf_file, c_char, t_zdf_iteration_c
            type(t_zdf_file)        :: file
            type(t_zdf_iteration_c) :: iter
        end subroutine zdf_add_iteration_c
    end interface

    type( t_zdf_iteration_c ) :: iter_c
    character(len = len( iter % time_units ) + 1 ), target :: time_units
    character(len = len(trim(iter % name))+1), target :: lname

    ! Prepare iteration info
    lname = trim(iter % name) // C_NULL_CHAR
    iter_c % name = c_loc(lname)
  
    time_units = trim( iter % time_units ) // C_NULL_CHAR
    iter_c % n = iter % n
    iter_c % t = iter % t
    iter_c % time_units = c_loc( time_units )

    call zdf_add_iteration_c( file, iter_c )

end subroutine zdf_add_iteration

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_add_grid_info( file, info )

    implicit none

    type(t_zdf_file), intent(in) :: file
    type(t_zdf_grid_info), intent(in) :: info

    interface
        subroutine zdf_add_grid_info_c( file, info ) bind(C, name="zdf_add_grid_info")
            import t_zdf_file, c_char, t_zdf_grid_info_c
            type(t_zdf_file)        :: file
            type(t_zdf_grid_info_c) :: info
        end subroutine zdf_add_grid_info_c
    end interface

    type(t_zdf_grid_info_c) :: info_c
    character(len = zdf_str_maxlen ), target :: name,label, units

    type(t_zdf_grid_axis_c), dimension(info % ndims), target :: axis
    character(len = zdf_str_maxlen ), dimension(info % ndims), target :: axis_name
    character(len = zdf_str_maxlen ), dimension(info % ndims), target :: axis_label
    character(len = zdf_str_maxlen ), dimension(info % ndims), target :: axis_units

    integer :: i

    ! Prepare grid info
    info_c % ndims = info % ndims
    do i = 1, info % ndims
        info_c % count(i) = info % count(i)
    enddo

    ! Label and units
    name  = trim( info % name ) // C_NULL_CHAR
    label = trim( info % label ) // C_NULL_CHAR
    units = trim( info % units ) // C_NULL_CHAR
    info_c % name  = c_loc( name )
    info_c % label = c_loc( label )
    info_c % units = c_loc( units )

    ! Axis
    do i = 1, info % ndims
        axis(i) % type = info % axis(i) % type
        axis(i) % min  = info % axis(i) % min
        axis(i) % max  = info % axis(i) % max
        axis_name(i) = trim(info % axis(i) % name) // C_NULL_CHAR
        axis_label(i) = trim(info % axis(i) % label) // C_NULL_CHAR
        axis_units(i) = trim(info % axis(i) % units) // C_NULL_CHAR
        axis(i) % name  = c_loc( axis_name(i) )
        axis(i) % label = c_loc( axis_label(i) )
        axis(i) % units = c_loc( axis_units(i) )
    enddo

    info_c % axis = c_loc( axis )

    call zdf_add_grid_info_c( file, info_c )

end subroutine zdf_add_grid_info

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_add_dataset( file, dataset )

    implicit none

    type(t_zdf_file), intent(in) :: file
    type(t_zdf_dataset), intent(inout) :: dataset

    interface
        subroutine zdf_add_dataset_c( file, dataset ) bind(C, name="zdf_add_dataset")
            import t_zdf_file, c_char, t_zdf_dataset_c
            type(t_zdf_file)       :: file
            type(t_zdf_dataset_c)    :: dataset
        end subroutine
    end interface

    type(t_zdf_dataset_c) :: dataset_c
    character(len = zdf_str_maxlen), target :: name
    integer :: i

    name = trim(dataset % name) // c_null_char
    dataset_c % name = c_loc( name )
    dataset_c % data_type = dataset % data_type
    dataset_c % ndims = dataset % ndims
    do i = 1, dataset % ndims
        dataset_c % count(i) = dataset % count(i)
    enddo
    dataset_c % data = dataset % data
    dataset_c % id = dataset % id
    dataset_c % offset = dataset % offset

    call zdf_add_dataset_c( file, dataset_c )

    ! Copy updated id to fortran variable
    dataset % id = dataset_c % id

end subroutine zdf_add_dataset


!------------------------------------------------------------------------------
!> Start ZDF chunked dataset
!! @param   file        ZDF file
!! @param   dataset     ZDF chunked dataset
subroutine zdf_start_cdset( file, dataset )

    implicit none

    type(t_zdf_file), intent(in) :: file
    type(t_zdf_dataset), intent(inout) :: dataset

    interface
        subroutine zdf_start_cdset_c( file, dataset ) bind(C, name="zdf_start_cdset")
            import t_zdf_file, c_char, t_zdf_dataset_c
            type(t_zdf_file)       :: file
            type(t_zdf_dataset_c)    :: dataset
        end subroutine zdf_start_cdset_c
    end interface

    type(t_zdf_dataset_c) :: dataset_c
    character(len = zdf_str_maxlen), target :: name
    integer :: i

    name = trim(dataset % name) // c_null_char
    dataset_c % name = c_loc( name )
    dataset_c % data_type = dataset % data_type
    dataset_c % ndims = dataset % ndims
    do i = 1, dataset % ndims
        dataset_c % count(i) = dataset % count(i)
    enddo
    dataset_c % data = dataset % data
    dataset_c % id = dataset % id
    dataset_c % offset = dataset % offset

    call zdf_start_cdset_c( file, dataset_c )

    ! Copy updated id to fortran variable
    dataset % id = dataset_c % id

end subroutine zdf_start_cdset

!------------------------------------------------------------------------------
!> Write ZDF chunked dataset
!> @param file      ZDF file
!> @param cdset     Chunked dataset
!> @param chunk     Data chunk
subroutine zdf_write_cdset( file, cdset, chunk )

    implicit none

    type(t_zdf_file), intent(in)  :: file
    type(t_zdf_dataset), intent(in) :: cdset
    type(t_zdf_chunk), intent(in)   :: chunk

    interface
        subroutine zdf_write_cdset_c( file, cdset, chunk ) bind(C, name="zdf_write_cdset")
            import t_zdf_file, t_zdf_dataset_c, t_zdf_chunk
            type(t_zdf_file)    :: file
            type(t_zdf_dataset_c) :: cdset
            type(t_zdf_chunk)   :: chunk
        end subroutine zdf_write_cdset_c
    end interface

    type(t_zdf_dataset_c) :: cdset_c
    character(len = zdf_str_maxlen), target :: name
    integer :: i

    name = trim(cdset % name) // c_null_char
    cdset_c % name = c_loc( name )
    cdset_c % data_type = cdset % data_type
    cdset_c % ndims = cdset % ndims
    do i = 1, cdset % ndims
        cdset_c % count(i) = cdset % count(i)
    enddo
    cdset_c % data = cdset % data
    cdset_c % id = cdset % id
    cdset_c % offset = cdset % offset

    call zdf_write_cdset_c( file, cdset_c, chunk )

end subroutine zdf_write_cdset

!------------------------------------------------------------------------------
!> Ends access to chunked dataset
!!
!! @param   file    ZDF file
!! @param   cdset   ZDF chunked dataset
subroutine zdf_end_cdset( file, cdset )

    implicit none

    type(t_zdf_file), intent(in)  :: file
    type(t_zdf_dataset), intent(in) :: cdset

    interface
        subroutine zdf_end_cdset_c( file, cdset ) bind(C, name="zdf_end_cdset")
            import t_zdf_file, t_zdf_dataset_c
            type(t_zdf_file)    :: file
            type(t_zdf_dataset_c) :: cdset
        end subroutine zdf_end_cdset_c
    end interface

    type(t_zdf_dataset_c) :: cdset_c
    character(len = zdf_str_maxlen), target :: name
    integer :: i

    name = trim(cdset % name) // c_null_char
    cdset_c % name = c_loc( name )
    cdset_c % data_type = cdset % data_type
    cdset_c % ndims = cdset % ndims
    do i = 1, cdset % ndims
        cdset_c % count(i) = cdset % count(i)
    enddo
    cdset_c % data = cdset % data
    cdset_c % id = cdset % id
    cdset_c % offset = cdset % offset

    call zdf_end_cdset_c( file, cdset_c )

end subroutine 

!------------------------------------------------------------------------------
!> Extends dataset dimensions
!!
!! @param   file        ZDF file
!! @param   cdset       ZDF chunked dataset
!! @param   new_count   New dataset dimensions
subroutine zdf_extend_cdset( file, cdset, new_count )

    implicit none

    type(t_zdf_file), intent(in) :: file
    type(t_zdf_dataset), intent(inout) :: cdset
    integer(c_int64_t), dimension(:), intent(in) :: new_count

    interface
        subroutine zdf_extend_cdset_c( file, cdset, new_count ) bind(C, name="zdf_extend_dataset")
            import t_zdf_file, t_zdf_dataset_c, c_int64_t
            type(t_zdf_file)    :: file
            type(t_zdf_dataset_c) :: cdset
            integer(c_int64_t)  :: new_count(*)
        end subroutine
    end interface

    type(t_zdf_dataset_c) :: cdset_c
    character(len = zdf_str_maxlen), target :: name
    integer :: i

    name = trim(cdset % name) // c_null_char
    cdset_c % name = c_loc( name )
    cdset_c % data_type = cdset % data_type
    cdset_c % ndims = cdset % ndims
    do i = 1, cdset % ndims
        cdset_c % count(i) = cdset % count(i)
    enddo
    cdset_c % data = cdset % data
    cdset_c % id = cdset % id
    cdset_c % offset = cdset % offset

    call zdf_extend_cdset_c( file, cdset_c, new_count )

    do i = 1, cdset % ndims
        cdset % count(i) = cdset_c % count(i)
    enddo

end subroutine


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_open_cdset( file, dataset )

    implicit none

    type(t_zdf_file), intent(in) :: file
    type(t_zdf_dataset), intent(inout) :: dataset

    interface
        subroutine zdf_open_cdset_c( file, dataset ) bind(C, name="zdf_open_dataset")
            import t_zdf_file, c_char, t_zdf_dataset_c
            type(t_zdf_file)       :: file
            type(t_zdf_dataset_c)    :: dataset
        end subroutine zdf_open_cdset_c
    end interface

    type(t_zdf_dataset_c) :: dataset_c
    character(len = zdf_str_maxlen), target :: name
    integer :: i

    name = trim(dataset % name) // c_null_char
    dataset_c % name = c_loc( name )
    dataset_c % data_type = dataset % data_type
    dataset_c % ndims = dataset % ndims
    do i = 1, dataset % ndims
        dataset_c % count(i) = dataset % count(i)
    enddo
    dataset_c % data = dataset % data
    dataset_c % id = dataset % id
    dataset_c % offset = dataset % offset

    call zdf_open_cdset_c( file, dataset_c )

    ! Copy updated id to fortran variable
    dataset % id = dataset_c % id

end subroutine zdf_open_cdset


!------------------------------------------------------------------------------
! *************************** High Level Interfaces ***************************
!------------------------------------------------------------------------------

subroutine zdf_save_grid_generic( data, type, info, iteration, path )

    implicit none

    type(c_ptr), intent(in) :: data
    integer, intent(in)     :: type
    type(t_zdf_grid_info), intent(inout) :: info
    type(t_zdf_iteration), intent(in) :: iteration
    character(len=*), intent(in) :: path

    type( t_zdf_iteration_c ) :: iter_c
    character(len = zdf_str_maxlen ), target :: time_units

    type(t_zdf_grid_info_c) :: info_c
    character(len = zdf_str_maxlen ), target :: label, units

    type(t_zdf_grid_axis_c), dimension(info % ndims), target :: axis
    character(len = zdf_str_maxlen ), dimension(info % ndims), target :: axis_label
    character(len = zdf_str_maxlen ), dimension(info % ndims), target :: axis_units

    integer :: i

    interface
        subroutine zdf_save_grid_c( data, type, info, iteration, path ) bind(C, name="zdf_save_grid")
            import t_zdf_grid_info_c, t_zdf_iteration_c, c_ptr, c_int32_t, c_char
            type(c_ptr), value      :: data
            integer(c_int32_t), value :: type
            type(t_zdf_grid_info_c) :: info
            type(t_zdf_iteration_c) :: iteration
            character(kind=c_char)  :: path(*)
        end subroutine zdf_save_grid_c
    end interface


    ! Prepare iteration info
    time_units = trim( iteration % time_units ) // C_NULL_CHAR
    iter_c % n = iteration % n
    iter_c % t = iteration % t
    iter_c % time_units = c_loc( time_units )

    ! Prepare grid info
    info_c % ndims = info % ndims
    do i = 1, info % ndims
        info_c % count(i) = info % count(i)
    enddo

    ! Label and units
    label = trim( info % label ) // C_NULL_CHAR
    units = trim( info % units ) // C_NULL_CHAR
    info_c % label = c_loc( label )
    info_c % units = c_loc( units )

    ! Axis
    do i = 1, info % ndims
        axis(i) % type = info % axis(i) % type
        axis(i) % min  = info % axis(i) % min
        axis(i) % max  = info % axis(i) % max
        axis_label(i) = trim(info % axis(i) % label) // C_NULL_CHAR
        axis_units(i) = trim(info % axis(i) % units) // C_NULL_CHAR
        axis(i) % label = c_loc( axis_label(i) )
        axis(i) % units = c_loc( axis_units(i) )
    enddo

    info_c % axis = c_loc( axis )

    call zdf_save_grid_c( data, type, info_c, iter_c, trim(path) // C_NULL_CHAR )

end subroutine zdf_save_grid_generic
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_save_grid_1D_r4( data, info, iteration, path )

    implicit none

    real(p_single), dimension( : ), target, intent(inout) :: data
    type(t_zdf_grid_info), intent(inout) :: info
    type(t_zdf_iteration), intent(in) :: iteration
    character(len=*), intent(in) :: path

    info % ndims = 1
    info % count(1) = size( data )

    call zdf_save_grid_generic( c_loc(data), zdf_float32, info, iteration, path )

end subroutine zdf_save_grid_1D_r4

subroutine zdf_save_grid_2D_r4( data, info, iteration, path )

    implicit none

    real(p_single), dimension( :, : ), target, intent(inout) :: data
    type(t_zdf_grid_info), intent(inout) :: info
    type(t_zdf_iteration), intent(in) :: iteration
    character(len=*), intent(in) :: path

    info % ndims = 2
    info % count(1) = size( data, 1 )
    info % count(2) = size( data, 2 )

    call zdf_save_grid_generic( c_loc(data), zdf_float32, info, iteration, path )

end subroutine zdf_save_grid_2D_r4

subroutine zdf_save_grid_3D_r4( data, info, iteration, path )

    implicit none

    real(p_single), dimension( :, :, : ), target, intent(inout) :: data
    type(t_zdf_grid_info), intent(inout) :: info
    type(t_zdf_iteration), intent(in) :: iteration
    character(len=*), intent(in) :: path

    info % ndims = 3
    info % count(1) = size( data, 1 )
    info % count(2) = size( data, 2 )
    info % count(3) = size( data, 3 )

    call zdf_save_grid_generic( c_loc(data), zdf_float32, info, iteration, path )

end subroutine zdf_save_grid_3D_r4
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_open_grid_file( file, info, iteration, path )

    implicit none

    interface
        subroutine zdf_open_grid_file_c( file, info, iteration, path ) bind(C, name="zdf_open_grid_file")
            import t_zdf_file, t_zdf_grid_info_c, t_zdf_iteration_c, c_char
            type(t_zdf_file)        :: file
            type(t_zdf_grid_info_c) :: info
            type(t_zdf_iteration_c) :: iteration
            character(kind=c_char)  :: path(*)
        end subroutine zdf_open_grid_file_c
    end interface

    type(t_zdf_file), intent(in) :: file
    type(t_zdf_grid_info), intent(in) :: info
    type(t_zdf_iteration), intent(in) :: iteration
    character(len=*), intent(in) :: path

    type( t_zdf_iteration_c ) :: iter_c
    character(len = zdf_str_maxlen ), target :: time_units

    type(t_zdf_grid_info_c) :: info_c
    character(len = zdf_str_maxlen ), target :: label, units

    type(t_zdf_grid_axis_c), dimension(info % ndims), target :: axis
    character(len = zdf_str_maxlen ), dimension(info % ndims), target :: axis_label
    character(len = zdf_str_maxlen ), dimension(info % ndims), target :: axis_units

    integer :: i

    ! Prepare iteration info
    time_units = trim( iteration % time_units ) // C_NULL_CHAR
    iter_c % n = iteration % n
    iter_c % t = iteration % t
    iter_c % time_units = c_loc( time_units )

    ! Prepare grid info
    info_c % ndims = info % ndims
    do i = 1, info % ndims
        info_c % count(i) = info % count(i)
    enddo

    ! Label and units
    label = trim( info % label ) // C_NULL_CHAR
    units = trim( info % units ) // C_NULL_CHAR

    info_c % label = c_loc( label )
    info_c % units = c_loc( units )

    ! Axis
    do i = 1, info % ndims
        axis(i) % type = info % axis(i) % type
        axis(i) % min  = info % axis(i) % min
        axis(i) % max  = info % axis(i) % max
        axis_label(i) = trim(info % axis(i) % label) // C_NULL_CHAR
        axis_units(i) = trim(info % axis(i) % units) // C_NULL_CHAR
        axis(i) % label = c_loc( axis_label(i) )
        axis(i) % units = c_loc( axis_units(i) )
    enddo

    info_c % axis = c_loc( axis )

    call zdf_open_grid_file_c( file, info_c, iter_c, trim(path) // C_NULL_CHAR )

end subroutine zdf_open_grid_file
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_add_part_info( file, info )

    implicit none

    type(t_zdf_file), intent(in) :: file
    type(t_zdf_part_info), intent(in) :: info

    interface
        subroutine zdf_add_part_info_c( file, info ) bind(C, name="zdf_add_part_info")
            import t_zdf_file, c_char, t_zdf_part_info_c
            type(t_zdf_file)        :: file
            type(t_zdf_part_info_c) :: info
        end subroutine zdf_add_part_info_c
    end interface

    type( t_zdf_part_info_c ) :: info_c
    character(len = zdf_str_maxlen ), target :: name, label
    character(len = zdf_str_maxlen ), dimension(:), pointer :: quants, qlabels, qunits
    type(c_ptr), dimension(64), target :: pquants, plabels, punits

    integer :: i

    ! Prepare particles info
    name = trim(info % name) // C_NULL_CHAR
    label = trim(info % label) // C_NULL_CHAR
    info_c % name = c_loc(name)
    info_c % label = c_loc(label)
    info_c % np = info % np
    info_c % nquants = info % nquants

    allocate( quants( info % nquants ))
    allocate( qlabels( info % nquants ))
    allocate( qunits( info % nquants ))

    do i = 1, info % nquants
        quants(i) = trim(info % quants(i)) // C_NULL_CHAR
        qlabels(i) = trim(info % qlabels(i)) // C_NULL_CHAR
        qunits(i)  = trim(info % qunits(i)) // C_NULL_CHAR

        pquants(i) = c_loc( quants(i) )
        plabels(i) = c_loc( qlabels(i) )
        punits(i)  = c_loc( qunits(i) )
    enddo

    info_c % quants = c_loc( pquants )
    info_c % qlabels = c_loc( plabels )
    info_c % qunits  = c_loc( punits )

    call zdf_add_part_info_c( file, info_c )

    deallocate( quants )
    deallocate( qlabels )
    deallocate( qunits )

end subroutine zdf_add_part_info


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_add_track_info( file, info )

    implicit none

    type(t_zdf_file), intent(in) :: file
    type(t_zdf_track_info), intent(in) :: info

    interface
        subroutine zdf_add_track_info_c( file, info ) bind(C, name="zdf_add_track_info")
            import t_zdf_file, c_char, t_zdf_track_info_c
            type(t_zdf_file)        :: file
            type(t_zdf_track_info_c) :: info
        end subroutine zdf_add_track_info_c
    end interface

    type( t_zdf_track_info_c ) :: info_c
    character(len = zdf_str_maxlen ), target :: lname, llabel
    character(len = zdf_str_maxlen ), dimension(:), pointer :: quants, labels, units
    type(c_ptr), dimension(64), target :: pquants, plabels, punits

    integer :: i

    ! Prepare particles info
    lname = trim(info % name) // C_NULL_CHAR
    llabel = trim(info % label) // C_NULL_CHAR
    info_c % name = c_loc(lname)
    info_c % label = c_loc(llabel)
    info_c % ntracks = info % ntracks
    info_c % ndump   = info % ndump
    info_c % niter   = info % niter
    info_c % nquants = info % nquants

    allocate( quants( info % nquants ))
    allocate( labels( info % nquants ))
    allocate( units( info % nquants ))

    do i = 1, info % nquants
        quants(i) = trim(info % quants(i)) // C_NULL_CHAR
        labels(i) = trim(info % qlabels(i)) // C_NULL_CHAR
        units(i)  = trim(info % qunits(i)) // C_NULL_CHAR

        pquants(i) = c_loc( quants(i) )
        plabels(i) = c_loc( labels(i) )
        punits(i)  = c_loc( units(i) )
    enddo

    info_c % quants = c_loc( pquants )
    info_c % qlabels = c_loc( plabels )
    info_c % qunits  = c_loc( punits )

    call zdf_add_track_info_c( file, info_c )

    deallocate( units )
    deallocate( labels )
    deallocate( quants )

end subroutine zdf_add_track_info

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_open_part_file( file, info, iteration, path )

    implicit none

    interface
        subroutine zdf_open_part_file_c( file, info, iteration, path ) bind(C, name="zdf_open_part_file")
            import t_zdf_file, t_zdf_part_info_c, t_zdf_iteration_c, c_char
            type(t_zdf_file)        :: file
            type(t_zdf_part_info_c) :: info
            type(t_zdf_iteration_c) :: iteration
            character(kind=c_char)  :: path(*)
        end subroutine zdf_open_part_file_c
    end interface

    type(t_zdf_file), intent(in) :: file
    type(t_zdf_part_info), intent(in) :: info
    type(t_zdf_iteration), intent(in) :: iteration
    character(len=*), intent(in) :: path

    type( t_zdf_iteration_c ) :: iter_c
    character(len = zdf_str_maxlen ), target :: time_units

    type( t_zdf_part_info_c ) :: info_c
    character(len = zdf_str_maxlen ), target :: name, label
    character(len = zdf_str_maxlen ), dimension(:), pointer :: quants, qlabels, qunits
    type(c_ptr), dimension(64), target :: pquants, plabels, punits

    integer :: i

    ! Prepare iteration info
    time_units = trim( iteration % time_units ) // C_NULL_CHAR
    iter_c % n = iteration % n
    iter_c % t = iteration % t
    iter_c % time_units = c_loc( time_units )

    ! Prepare particles info
    name = trim(info % name) // C_NULL_CHAR
    label = trim(info % label) // C_NULL_CHAR
    info_c % name = c_loc(name)
    info_c % label = c_loc(label)
    info_c % nquants = info % nquants

    allocate( quants( info % nquants ))
    allocate( qlabels( info % nquants ))
    allocate( qunits( info % nquants ))

    do i = 1, info % nquants
        quants(i) = trim(info % quants(i)) // C_NULL_CHAR
        qlabels(i) = trim(info % qlabels(i)) // C_NULL_CHAR
        qunits(i)  = trim(info % qunits(i)) // C_NULL_CHAR

        pquants(i) = c_loc( quants(i) )
        plabels(i) = c_loc( qlabels(i) )
        punits(i)  = c_loc( qunits(i) )
    enddo

    info_c % quants  = c_loc( pquants )
    info_c % qlabels = c_loc( plabels )
    info_c % qunits  = c_loc( punits )
    info_c % np = info % np

    call zdf_open_part_file_c( file, info_c, iter_c, trim(path) // C_NULL_CHAR )

    deallocate( qunits )
    deallocate( qlabels )
    deallocate( quants )

end subroutine zdf_open_part_file
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine zdf_add_quant_part_file( file, name, data, np )

    implicit none

    interface
        subroutine zdf_add_quant_part_file_c( file, name, data, np ) bind(C, name="zdf_add_quant_part_file")
            import t_zdf_file, c_char, c_int, c_float
            type(t_zdf_file)       :: file
            character(kind=c_char) :: name(*)
            real(c_float)          :: data(*)
            integer(c_int), value  :: np
        end subroutine zdf_add_quant_part_file_c
    end interface

    type(t_zdf_file), intent(in) :: file
    character(len=*), intent(in) :: name
    real(p_single),dimension(:), intent(in) :: data
    integer, intent(in) :: np

    call zdf_add_quant_part_file_c( file, trim(name) // C_NULL_CHAR, data, np )

end subroutine zdf_add_quant_part_file
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Returns storage size in bytes for each datatype
!! @param   zdf_type
!! @return  Storage size in bytes
function zdf_sizeof( data_type ) result(sizeof)
    
    implicit none
    
    integer, intent(in) :: data_type
    integer :: sizeof
    
    select case ( data_type )
    case( zdf_null )
        sizeof = 0
    case( zdf_int8 )
        sizeof = 1
    case( zdf_uint8 )
        sizeof = 1
    case( zdf_int16 )
        sizeof = 2
    case( zdf_uint16 )
        sizeof = 2
    case( zdf_int32 )
        sizeof = 4
    case( zdf_uint32 )
        sizeof = 4
    case( zdf_float32 )
        sizeof = 4
    case( zdf_int64 )
        sizeof = 8
    case( zdf_uint64 )
        sizeof = 8
    case( zdf_float64 )
        sizeof = 8
    case default
        sizeof = 0
    end select
    
end function

!------------------------------------------------------------------------------
!> Returns storage size in bytes for each datatype
!! @param   rkind   Real kind type parameter
!! @return          Equivalent ZDF datatype
function freal_to_zdftype( rkind )
    
    implicit none
    
    integer, intent(in) :: rkind
    integer :: freal_to_zdftype
    
    select case(rkind)
    case(p_single)
        freal_to_zdftype = zdf_float32
    case(p_double)
        freal_to_zdftype = zdf_float64
    case default
        freal_to_zdftype = zdf_null
    end select
    
end function


end module zdf

!------------------------------------------------------------------------------


#ifdef __ZDF_TEST__

program zdf_test

    use zdf
    use m_diagfile
    implicit none

    !type(t_zdf_file) :: file

    type(t_zdf_grid_info) :: grid
    type(t_zdf_part_info) :: part
    type(t_zdf_iteration) :: iter
    type(t_zdf_file) :: file
    type(t_zdf_dataset) :: cdset
    type(t_zdf_chunk) :: chunk

    real(p_single), dimension(16,8) :: data
    real(p_single), dimension(8,4), target :: tile
    character(len=zdf_str_maxlen), dimension(2), target :: quants, units

    integer :: i, j

    iter % n = 1
    iter % t = 0.0
    iter % time_units = 'time units'

    do j = 1, size(data,2)
        do i = 1, size(data,1)
            data(i,j) = i + (j-1) * size(data,1)
        enddo
    enddo

    ! Grid
    grid % ndims = 2
    grid % count(1:2) = [16,8]

    grid % label = 'label'
    grid % units = 'units'

    grid % axis(1) % type  = zdf_axis_linear
    grid % axis(1) % min   = -1.0
    grid % axis(1) % max   = +1.0
    grid % axis(1) % label = 'x1'
    grid % axis(1) % units = 'x1 units'

    grid % axis(2) % type  = zdf_axis_linear
    grid % axis(2) % min   = -1.0
    grid % axis(2) % max   = +1.0
    grid % axis(2) % label = 'x2'
    grid % axis(2) % units = 'x2 units'

    ! Save grid all at once
    !call zdf_save_grid( data, grid, iter, '.')

    call zdf_open_grid_file( file, grid, iter, '.' )

    cdset%data_type = zdf_float32
    cdset%ndims = 2
    cdset%count(1:2) = [16,8]
    call zdf_start_cdset( file, "DATA", cdset )

    chunk%count(1:2)  = [8,4]
    chunk%stride(1:2) = [1,1]
    chunk%data = c_loc( tile )

    chunk%start(1:2) = [0,0]
    tile = data( 1:8, 1:4 )
    call zdf_write_cdset( file, cdset, chunk )

    chunk%start(1:2) = [8,0]
    tile = data( 9:16, 1:4 )
    call zdf_write_cdset( file, cdset, chunk )

    chunk%start(1:2) = [0,4]
    tile = data( 1:8, 5:8 )
    call zdf_write_cdset( file, cdset, chunk )

!    chunk%start(1:2) = [8,4]
!    tile = data( 9:16, 5:8 )
!    call zdf_write_cdset( file, cdset, chunk )

    call zdf_end_cdset( file, cdset )

    ! Particles

    ! part % name = 'electrons'
    ! part % nquants = 2
    ! part % np = 16
    ! quants = ['x1','x2']
    ! part % quants => quants
    ! units = ['x1 units','x2 units']
    ! part % units  => units

    ! call zdf_open_part_file( file, part, iter, '.')
    ! call zdf_add_quant_part_file( file, "x1", data(:,1), 16 )
    ! call zdf_add_quant_part_file( file, "x2", data(:,2), 16 )
    ! call zdf_close_file( file )

end program zdf_test

#endif
