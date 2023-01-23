!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     system module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_system

use iso_c_binding
use mpi

implicit none

! give general access to all information in this module
public

! flag that knows if mpi has started
logical, private :: os_mpi_started = .false.

! maximum number of spatial dimensions
integer, parameter :: p_max_dim = 3

!==============================================================
! parameters essential at compile time that depend on the system
! the program is running on
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Constants defining datatypes
integer, public, parameter :: p_single = kind(1.0e0)
integer, public, parameter :: p_double = kind(1.0d0)
integer, public, parameter :: p_byte   = selected_int_kind(2)
integer, public, parameter :: p_int64  = selected_int_kind(10)

#ifdef HAS_QUAD_PREC

! 128 bit floating point. This definition ensures portability between all
! compilers tested
integer, parameter :: p_quad   = max( selected_real_kind( 18 ), &
                                      selected_real_kind( 19 ) )

! MPI type for quad precision communication
integer :: mpi_quad = -1

! Custom quad precision operations for reduce operations
! These are defined at the end of the file after the end of the module
integer, external  :: sum_quad, max_quad

! operation codes
integer  :: mpi_quad_sum, mpi_quad_max

#endif

! Vector width for SIMD code. This affects:
! - particle buffer sizes (must be a multiple of this)
! - OpenMP code: the number of particles in each thread must be a multiple of this
! This should be done on a system by system basis, but for now just set to the max size
! on any current system (16 on single precision MIC code)
integer, public, parameter :: p_vecwidth = 16


! directory seperator for this system
character, parameter :: p_dir_sep  = '/'      ! unix like systems

integer, parameter :: p_stderr = 0
integer, parameter :: p_stdin  = 5
integer, parameter :: p_stdout = 6

! max size of filename
integer, parameter :: p_max_filename_len = 256
integer, parameter :: file_id_tem        =  10
integer, parameter :: p_err_assert       = -21 ! assertion

!=================================================================================================
! High resolution timer interfaces
! these routines are defined in os-sys-multi-c.c
!-------------------------------------------------------------------------------------------------

! number of seconds elapsed since a fixed time in the past

!real(p_double), external :: timer_cpu_seconds
interface
  real( c_double ) function timer_cpu_seconds() bind(c, name="timer_cpu_seconds")
    use iso_c_binding
  end function timer_cpu_seconds
end interface

! minimal difference between successive high resolution timer
! calls (note: this is usual much less than the actual resolution
! due to the conversion to seconds)

!real(p_double), external :: timer_resolution
interface
  real( c_double ) function timer_resolution() bind(c, name="timer_resolution")
    use iso_c_binding
  end function timer_resolution
end interface


! number of ticks since a fixed time in the past
!integer(p_int64), external :: timer_ticks
interface
  integer( C_INT64_T ) function timer_ticks() bind(c, name="timer_ticks")
    use iso_c_binding
  end function timer_ticks
end interface

! Convert tick interval to seconds
!real(p_double), external :: timer_interval_seconds
interface
  real( c_double ) function timer_interval_seconds(a,b) bind(c, name="timer_interval_seconds")
    use iso_c_binding
    integer(C_INT64_T), intent(in) :: a,b
  end function timer_interval_seconds
end interface

interface
  subroutine timer_init() bind(c, name="timer_init")
    use iso_c_binding
  end subroutine timer_init
end interface

interface
  subroutine get_vmrss( a ) bind(c, name="get_vmrss")
    use iso_c_binding
    integer(C_INT64_T), intent(inout) :: a
  end subroutine get_vmrss
end interface

interface
  subroutine get_vmsize( a ) bind(c, name="get_vmsize")
    use iso_c_binding
    integer(C_INT64_T), intent(inout) :: a
  end subroutine get_vmsize
end interface

! buffer variables for debug/error routines
! used by the DEBUG, LOG, ERROR and WARNING macros
character(len = 1024) :: dbg_buf__, err_buf__, wrn_buf__

! Structure to hold command line / environment options and general simulation options
type :: t_options

  logical :: test    = .false.
  logical :: restart = .false.
  character( len = p_max_filename_len ) :: input_file = ''
  character( len = p_max_filename_len ) :: work_dir  = ''

  integer :: random_seed = 0      ! Random seed for random number generator
  integer :: random_class = 0
  logical :: enforce_rng_constancy = .false.

  integer :: algorithm   = 0      ! Algorithm for the simulation (normal / hybrid / pgc)
  real( p_double ) :: omega_p0 = -1.0      ! Reference frequency
  real( p_double ) ::  n0 = -1.0        ! Reference density

  real( p_double ) ::  gamma = -1.0       ! Lorentz boosted frame gamma

  integer :: ndump_prof = 0        ! Frequency for profiling dumps

  real( p_double ) :: wall_clock_start = 0       ! Wall clock time at the beggining of the run (seconds)
  real( p_double ) :: wall_clock_limit = -1     ! Wall clock limit time (seconds)
  integer          :: wall_clock_check = -1     ! Frequency at which to check if wall time limit
                                              ! has been reached

  logical :: wall_clock_checkpoint = .true.  ! Dump checkpointing information when stopping due to
                                        ! wall clock limit

  ! integer at which time steps timings will be started and ended to be measured
  integer :: prof_start = 0
  integer :: prof_end   = 0

  ! File format to use for diagnostics output
  integer :: file_format

  ! Parallel I/O algorithm to use
  integer :: parallel_io

end type t_options

! General simulation information
! Used to tag the diagnostics file
type :: t_sim_info

  character(len = 1024) :: git_version
  character(len = 1024) :: compile_time

  character(len = 1024) :: input_file ! Name of the input file used
  character(len = 1024) :: timestamp  ! output of date_and_time command at the beggining of the sim.

  real(p_double)        :: input_file_crc32 ! a simple checksum so that the input file can be matched.
  real(p_double)        :: dt               ! timestep size

  integer                                 :: ndims
  integer, dimension(p_max_dim)           :: box_nx   ! Simulation box size in grid cells
  real( p_double ), dimension(p_max_dim)  :: box_min  ! coordinates of the lower box corner
  real( p_double ), dimension(p_max_dim)  :: box_max  ! coordinates of the upper box corner
  logical, dimension(p_max_dim)           :: periodic ! Periodic simulation
  logical, dimension(p_max_dim)           :: move_c   ! Moving window simulation

  integer, dimension(p_max_dim)  :: node_conf  ! parallel node configuration
  integer, dimension(:), pointer :: nx_x1_node => null() ! nodes along x1
  integer, dimension(:), pointer :: nx_x2_node => null() ! nodes along x2
  integer, dimension(:), pointer :: nx_x3_node => null() ! nodes along x3

end type t_sim_info


interface isnan
  module procedure isnan_single
  module procedure isnan_double
end interface

interface isinf
  module procedure isinf_single
  module procedure isinf_double
end interface

interface ffloor
  module procedure fastfloor_single
  module procedure fastfloor_double
end interface

interface file_link
  module procedure file_link
end interface

interface abort_program
  module procedure abort_program
end interface

interface mpi_real_type
  module procedure mpi_real_type
end interface

interface file_exists
   module procedure file_exists
end interface

interface mpi_node
   module procedure mpi_node
end interface

interface system_init
  module procedure system_init
end interface system_init

interface system_cleanup
  module procedure system_cleanup
end interface system_cleanup

#ifdef __bgp__
interface bgp_core_memory
  module procedure bgp_core_memory
end interface

interface bgp_used_memory
  module procedure bgp_used_memory
end interface

#endif

interface getopt
  module procedure getopt
end interface getopt

interface memcpy
  module procedure memcpy_1d_int
  module procedure memcpy_2d_int

  module procedure memcpy_1d_r4
  module procedure memcpy_2d_r4
  module procedure memcpy_3d_r4
  module procedure memcpy_4d_r4
  module procedure memcpy_3d_1d_r4

  module procedure memcpy_1d_r8
  module procedure memcpy_2d_r8
  module procedure memcpy_3d_r8
  module procedure memcpy_4d_r8
  module procedure memcpy_3d_1d_r8
end interface

interface setnan
  module procedure setnan_1d_r4
  module procedure setnan_2d_r4
  module procedure setnan_3d_r4
  module procedure setnan_4d_r4

  module procedure setnan_1d_r8
  module procedure setnan_2d_r8
  module procedure setnan_3d_r8
  module procedure setnan_4d_r8
end interface

! interfaces to system routines
interface
  subroutine memset(b,c,len) bind(c, name="memset")
    use iso_c_binding
    type(c_ptr), intent(in),       value :: b
    integer(c_int), intent(in),    value :: c
    integer(c_size_t), intent(in), value :: len
  end subroutine
end interface

interface
  subroutine zero(a,len) bind(c, name="zero")
    use iso_c_binding
    type(c_ptr), intent(in),       value :: a
    integer(c_size_t), intent(in), value :: len
  end subroutine
end interface

interface
  subroutine memcpy_c(dst,src,n) bind(c, name="memcopy")
    use iso_c_binding
    type(c_ptr), intent(in),       value :: dst
    type(c_ptr), intent(in),       value :: src
    integer(c_size_t), intent(in), value :: n
  end subroutine
end interface

! interfaces to wrappers in os-sys-multi-c
interface
  subroutine setnan_r4(s1, n) bind(c, name="setnan_r4")
    use iso_c_binding
    type(c_ptr), intent(in), value :: s1
    integer(c_int), intent(in), value  :: n
  end subroutine setnan_r4
end interface

interface
  subroutine setnan_r8(s1, n) bind(c, name="setnan_r8")
    use iso_c_binding
    type(c_ptr), intent(in), value :: s1
    integer(c_int) :: n
  end subroutine setnan_r8
end interface

! profiling flag
!   - if true, `begin_event` and `end_event` will get the timings
!   - if false, `begin_event` and `end_event` will not do anything
!   - state can only be changed by subroutine `change_profiling`
logical, private :: p_profiling = .false.

interface profiling
  module procedure set_profiling
end interface

interface if_profiling
  module procedure if_profiling
end interface

contains

!---------------------------------------------------------------------------------------------------
! Initialize system settings
!---------------------------------------------------------------------------------------------------
subroutine system_init( options )
!---------------------------------------------------------------------------------------------------

  implicit none

  type( t_options ), intent(inout) :: options

  integer :: ierr
  character(len = 1024) :: osiris_wdir

  ! if only testing the input file the following is skipped
  if ( .not. options%test ) then

    call mpi_init( ierr )
    if ( ierr /= MPI_SUCCESS ) then
      write(0,'(A)') "(*error*) Unable to initialize MPI"
      stop
    endif

    ! set the os_mpi_started flag to true
    os_mpi_started = .true.

    ! printout host running osiris
    SCR_ROOT("Osiris running on host ",trim(hostname()) )

    ! set umask to O'022' => file permissions will be O0755
    call sys_umask(18) ! 18 (int) is 022 (octal)

    ! Change working directory if needed
    if ( trim(options%work_dir) /= '' ) then

      if ( mpi_node()==0 ) then
        print *, 'Changing working directory to "',trim(options%work_dir),'"'
      endif

      call sys_chdir(options%work_dir, ierr)

      if (ierr/=0) then
        SCR_MPINODE("(*error*) ", trim(strerror(ierr)) )
        SCR_MPINODE("(*error*) Unable to change directory to ",trim(options%work_dir))
        call abort_program()
      endif
    endif

    call sys_getcwd( osiris_wdir, ierr )
    if (ierr /= 0) then
      SCR_MPINODE("(*error*) Unable to get current working directory:")
      SCR_MPINODE("(*error*) ", trim(strerror(ierr)))
      SCR_MPINODE("(*error*) Aborting...")
      call abort_program()
    endif

    SCR_ROOT("Current working directory = ", trim(osiris_wdir))

    SCR_ROOT('Timer resolution is ', timer_resolution(), ' seconds')

#ifdef HAS_QUAD_PREC

    ! Create quad mpi datatype
    call mpi_type_contiguous( 16, MPI_BYTE, mpi_quad, ierr )
    call mpi_type_commit( mpi_quad, ierr )

    ! Create sum operation for quads
    call mpi_op_create( sum_quad, .true., mpi_quad_sum, ierr )
    call mpi_op_create( max_quad, .true., mpi_quad_max, ierr )

#endif

  endif

end subroutine system_init
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine system_cleanup()

  implicit none

  integer :: ierr

#ifdef HAS_QUAD_PREC

  ! Release the custom operation handles
  call mpi_op_free( mpi_quad_sum, ierr )
  call mpi_op_free( mpi_quad_max, ierr )

  ! Release the type handle
  call mpi_type_free( mpi_quad, ierr )

#endif

  ! close MPI
  call mpi_finalize( ierr )

end subroutine system_cleanup
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function file_exists( filename )
!---------------------------------------------------------------------------------------------------
! Check is supplied filename exists and is readable
! - This is implemented in pure fortran, a posix version using fstat() or access() is also
!   possible
!---------------------------------------------------------------------------------------------------

  implicit none

  logical :: file_exists
  character(len = *), intent(in) :: filename

  integer :: ierr

  open( file_id_tem, file = trim(filename), status = 'old', iostat = ierr, action = 'read' )
  if ( ierr == 0 ) then
    close( file_id_tem )
    file_exists = .true.
  else
    file_exists = .false.
  endif

end function file_exists
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! returns the correct mpi type for the type kind k
!---------------------------------------------------------------------------------------------------
function mpi_real_type( k )

  implicit none

  integer, intent(in) :: k
  integer :: mpi_real_type

  mpi_real_type = -1
  select case (k)
    case (p_double)
      mpi_real_type = MPI_DOUBLE_PRECISION
    case (p_single)
      mpi_real_type = MPI_REAL

#ifdef HAS_QUAD_PREC
    case (p_quad)
      mpi_real_type = mpi_quad
#endif

  end select


end function mpi_real_type
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine file_link( target, name, ierr_out )
!---------------------------------------------------------------------------------------------------
! create a symbolic link "name" to the file "target"
!---------------------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  character(len=*), intent(in) :: target, name
  integer, intent(out), optional :: ierr_out

  ! local variables
  integer :: ierr

  ! remove old link
  call remove( name, ierr )

  ! create new one
  call symlink( target, name, ierr )

  ! return result if necessary
  if (present(ierr_out)) then
    ierr_out = ierr
  endif

end subroutine file_link
!---------------------------------------------------------------------------------------------------


!************************************************************************
! Fortran wrappers for os-sys-multi-c.c
!************************************************************************

!---------------------------------------------------------------------------------------------------
function hostname()
!---------------------------------------------------------------------------------------------------
!        gets the hostname
!---------------------------------------------------------------------------------------------------

  implicit none

  interface
    subroutine gethostname_f(hostname, ierr)  bind(c, name="gethostname_f")
      use iso_c_binding
      character(kind=c_char), intent(out) :: hostname(*)
      integer(c_int), intent(out) :: ierr
    end subroutine gethostname_f
  end interface

  character(len=256) :: hostname
  integer :: ierr

  character(len=256) :: lhostname

  lhostname = " "
  call gethostname_f(lhostname, ierr)
  call cleanup_cstring(lhostname)

  hostname = trim(lhostname)

end function hostname
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! gets the mpi node number (for MPI_COMM_WORLD communicator)
!---------------------------------------------------------------------------------------------------
function mpi_node()

  implicit none

  integer :: mpi_node

  integer:: ierr

  if (os_mpi_started) then
    call MPI_COMM_RANK( mpi_comm_world, mpi_node, ierr )
  else
    mpi_node = 0
  endif

end function mpi_node
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! gets the mpi communicator
!---------------------------------------------------------------------------------------------------
function mpi_world_comm()

  implicit none

  integer :: mpi_world_comm

  if (os_mpi_started) then
    mpi_world_comm = MPI_COMM_WORLD
  else
    mpi_world_comm = MPI_COMM_NULL
  endif

end function mpi_world_comm
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! set the current umask
!---------------------------------------------------------------------------------------------------
subroutine sys_umask( cmask )
!---------------------------------------------------------------------------------------------------
  implicit none

  interface
    subroutine umask_f(cmask)  bind(c, name="umask_f")
      use iso_c_binding
      integer(c_int), intent(in), value :: cmask
    end subroutine umask_f
  end interface

  integer, intent(in) :: cmask

  call umask_f( cmask )

end subroutine sys_umask
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!        create a directory
!---------------------------------------------------------------------------------------------------
subroutine mkdir(path, ierr)

  implicit none

  interface
    subroutine mkdir_f(path, ierr)  bind(c, name="mkdir_f")
      use iso_c_binding
      character(kind=c_char) :: path(*)
      integer(c_int) :: ierr
    end subroutine mkdir_f
  end interface

  character(len=*), intent(in) :: path
  integer, intent(out) :: ierr

  if ( trim(path) /= '' ) then
    call mkdir_f( trim(path) // C_NULL_CHAR, ierr)
  else
    ierr = 0
  endif

end subroutine mkdir
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!        check to see if directory exists.
!---------------------------------------------------------------------------------------------------
function dir_exists(path)

  implicit none

  interface
    subroutine does_dir_exist_f(path, ierr)  bind(c, name="does_dir_exist_f")
      use iso_c_binding
      character(kind=c_char) :: path(*)
      integer(c_int) :: ierr
    end subroutine does_dir_exist_f
  end interface
  
  logical :: dir_exists
  character(len=*), intent(in) :: path
  integer :: ierr
  
  ierr = 0
  if ( trim(path) /= '' ) then
    call does_dir_exist_f(  trim(path) // C_NULL_CHAR, ierr )
    if(ierr == 0) then
      dir_exists = .false.
    else
      dir_exists = .true.
    endif
  else
    dir_exists = .false.
  endif

end function dir_exists
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
!        change current working directory
!---------------------------------------------------------------------------------------------------
subroutine sys_chdir(path, ierr)
!---------------------------------------------------------------------------------------------------

  implicit none

  interface
    subroutine chdir_f(path, ierr)  bind(c, name="chdir_f")
      use iso_c_binding
      character(kind=c_char) :: path(*)
      integer(c_int) :: ierr
    end subroutine chdir_f
  end interface

  character(len=*), intent(in) :: path
  integer, intent(out) :: ierr

  call chdir_f(trim(path) // C_NULL_CHAR, ierr)

end subroutine sys_chdir
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!        get current working directory
!---------------------------------------------------------------------------------------------------
subroutine sys_getcwd(path, ierr)

  implicit none

  interface
    subroutine getcwd_f(path, len, ierr)  bind(c, name="getcwd_f")
      use iso_c_binding
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), intent(in), value :: len
      integer(c_int), intent(out) ::  ierr
    end subroutine getcwd_f
  end interface

  character(len=*), intent(out) :: path
  integer, intent(out) :: ierr


  path = " "
  call getcwd_f(path, len(path), ierr)
  call cleanup_cstring(path)


end subroutine sys_getcwd
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
!        get the name of a system error
!---------------------------------------------------------------------------------------------------
function strerror(ierr)
!---------------------------------------------------------------------------------------------------

  implicit none

  interface
    subroutine strerror_f(errnum, err)  bind(c, name="strerror_f")
      use iso_c_binding
      integer(c_int), intent(in), value :: errnum
      character(kind=c_char), intent(out) :: err(*)
    end subroutine strerror_f
  end interface

  integer, intent(in) :: ierr
  character(len = 256) :: strerror

  strerror = ""
  call strerror_f(ierr, strerror)
  call cleanup_cstring(strerror)

end function strerror
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
!        create a symbolic link
!---------------------------------------------------------------------------------------------------
subroutine symlink(targetf, linkf, ierr)
!---------------------------------------------------------------------------------------------------

  implicit none

  interface
    subroutine symlink_f(name1, name2, ierr)  bind(c, name="symlink_f")
      use iso_c_binding
      character(kind=c_char) :: name1(*)
      character(kind=c_char) :: name2(*)
      integer(c_int) :: ierr
    end subroutine symlink_f
  end interface

  character(len=*), intent(in) :: targetf, linkf
  integer, intent(out) :: ierr

  call symlink_f( trim(targetf) // C_NULL_CHAR, trim(linkf) // C_NULL_CHAR, ierr)

 end subroutine symlink
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!        remove a file
!---------------------------------------------------------------------------------------------------
subroutine remove(path, ierr)
!---------------------------------------------------------------------------------------------------

  implicit none

  interface
    subroutine remove_f(path, ierr)  bind(c, name="remove_f")
      use iso_c_binding
      character(kind=c_char) :: path(*)
      integer(c_int) :: ierr
    end subroutine remove_f
  end interface

  character(len=*), intent(in) :: path
  integer, intent(out) :: ierr

  call remove_f(trim(path)//C_NULL_CHAR, ierr)

 end subroutine remove
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! get option character from command line argument list
!
! optstring  - string specifying options
! optidx    - command line argument index to process (start with 1)
! opt    - next known option character in optstring
! optarg  - option argument, if it is anticipated
! opt is set to '\0' when finished processing arguments
!---------------------------------------------------------------------------------------------------
subroutine getopt( optstring, optidx, opt, optarg  )
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(in) :: optstring
  integer, intent(inout) :: optidx
  character, intent(out) :: opt
  character(len = *), intent(out) :: optarg

  integer :: argc, idx
  character(len = p_max_filename_len) :: argv
  logical :: valid

  argc = command_argument_count()

  ! If over the last command line argument return '\0'
  if ( optidx > argc ) then
    opt = achar(0)
    return
  endif

  ! get argument value
  call get_command_argument( optidx, argv )

  if ( argv(1:1) == '-' ) then

    ! termination with '--'
    if ( len_trim(argv) == 2 .and. argv(2:2) == '-' ) then
    opt = achar(0)
    return
    endif

    ! look for options in optstring
    idx = 1
    valid = .false.
    do while ( idx <= len_trim(optstring) )
      ! found option
      if ( optstring(idx:idx) == argv(2:2) ) then
        ! check for option argument
        if ( idx + 1 <= len_trim(optstring) ) then
      if ( optstring(idx+1:idx+1) == ':' ) then
      optidx = optidx + 1
      if ( optidx > argc ) then
        exit
      else
        call get_command_argument( optidx, optarg )
      endif
      endif
        endif
        ! set opt to option and finish search
        opt = argv(2:2)
        valid = .true.
        exit
      endif
      idx = idx + 1
    enddo

    if ( .not. valid ) then
      opt = '?'
      return
    endif

    optidx = optidx + 1

  else
    opt = achar(0)
  endif


end subroutine getopt
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
function isnan_single( x )
!---------------------------------------------------------------------------------------------------
! Check if supplied value is a NAN
!---------------------------------------------------------------------------------------------------

  use, intrinsic :: ieee_arithmetic

  implicit none

  real( p_single ), intent(in) :: x
  logical :: isnan_single

  isnan_single = ieee_is_nan( x )

end function isnan_single
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function isnan_double( x )
!---------------------------------------------------------------------------------------------------
! Check if supplied value is a NAN
!---------------------------------------------------------------------------------------------------

  use, intrinsic :: ieee_arithmetic

  implicit none

  real( p_double ), intent(in) :: x
  logical :: isnan_double

  isnan_double = ieee_is_nan( x )

end function isnan_double
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function isinf_single( x )
!---------------------------------------------------------------------------------------------------
!  Check if number is infinity
!---------------------------------------------------------------------------------------------------

  use, intrinsic :: IEEE_ARITHMETIC

  implicit none

  real( p_single ), intent(in) :: x
  logical :: isinf_single

  isinf_single = ( .not. ieee_is_finite(x) )

end function isinf_single
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function isinf_double( x )
!---------------------------------------------------------------------------------------------------
!  Check if number is infinity
!---------------------------------------------------------------------------------------------------

  use, intrinsic :: IEEE_ARITHMETIC

  implicit none

  real( p_double ), intent(in) :: x
  logical :: isinf_double

  isinf_double = ( .not. ieee_is_finite(x) )

end function isinf_double
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! For their usability range these are actually (MUCH) faster than the native floor implementation
! (tested in the range [-2,2])
!---------------------------------------------------------------------------------------------------
function fastfloor_single(x)
!---------------------------------------------------------------------------------------------------
  implicit none
  real(p_single), intent(in) :: x

  integer :: fastfloor_single

  fastfloor_single = int(x)
  if ( x < 0 ) fastfloor_single = fastfloor_single - 1

end function fastfloor_single
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function fastfloor_double(x)
!---------------------------------------------------------------------------------------------------
  implicit none
  real(p_double), intent(in) :: x

  integer :: fastfloor_double

  fastfloor_double = int(x)
  if ( x < 0 ) fastfloor_double = fastfloor_double - 1

end function fastfloor_double
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine cleanup_cstring( str )
!---------------------------------------------------------------------------------------------------
! Cleans up a C string (null terminated) so that it is a proper fortran string
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(inout) :: str
  integer :: pos
  ! a C string will be terminated by a null character

  pos = index( str, char(0) )

  if (pos > 0) then
    str = str( 1:pos-1)
  endif


end subroutine cleanup_cstring
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine linebuf_stdio()
!---------------------------------------------------------------------------------------------------
! Force stdout and stderr to be line buffered
!---------------------------------------------------------------------------------------------------

  implicit none

  interface
    subroutine linebuf_stdio_f()  bind(c, name="linebuf_stdio_f")
      use iso_c_binding
    end subroutine linebuf_stdio_f
  end interface

  call linebuf_stdio_f()

 end subroutine linebuf_stdio
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!        compute the CRC32 checksum of a file.
!---------------------------------------------------------------------------------------------------
function compute_crc32( path )

  implicit none

  interface
    subroutine compute_crc32_f(path, crc32)  bind(c, name="compute_crc32_f")
      use iso_c_binding
      character(kind=c_char) :: path(*)
      real( c_double ) :: crc32
    end subroutine compute_crc32_f
  end interface

  character(len=*), intent(in) :: path
  real(p_double) :: compute_crc32, temp
  
  if ( trim(path) /= '' ) then
    call compute_crc32_f( trim(path) // C_NULL_CHAR, temp)
    compute_crc32 = temp
  else
    compute_crc32 = 0
  endif

end function compute_crc32
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!        Debug/Error routines
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine wrn__( file, line )
!---------------------------------------------------------------------------------------------------
!  Write warning message to screen and log file
!  Used by the WARNING macro
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(in) :: file
  integer, intent(in)          :: line

  write(0,*) '[',mpi_node(),'] (*warning*) ', trim(wrn_buf__), &
                                    " [", trim(adjustl(file)), &
                                    ":", line, "]"

end subroutine wrn__
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine err__( file, line )
!---------------------------------------------------------------------------------------------------
!  Write error message to screen
!  Used by the ERROR macro
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(in) :: file
  integer, intent(in)          :: line

  write(0,'(A,I0,A,A,A,A,A,I0,A)') ' (*error*) [', mpi_node(), ']', trim(err_buf__), &
                                    " [", trim(adjustl(file)), &
                                    ":", line, "]"

end subroutine err__
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine check_error( ierr, msg, code, file, line )
!---------------------------------------------------------------------------------------------------
!  Check for error and write error message to screen and abort when needed
!  Used by the CHECK_ERROR macro
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, intent(in) :: ierr
  character(len=*), intent(in) :: msg
  integer, intent(in) :: code
  character(len=*), intent(in) :: file
  integer, intent(in)          :: line

  if ( ierr /= 0 ) then

  write(0,'(A,I0,A,A,A,A,A,I0,A)') '[',mpi_node(),'] (*error*) ', trim(msg), &
                    " [", trim(adjustl(file)), &
                    ":", line, "]"

    call abort_program( code )

  endif

end subroutine check_error
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine dbg__( file, line )
!---------------------------------------------------------------------------------------------------
!  Write debug message to screen
!  Used by the DEBUG macro
!---------------------------------------------------------------------------------------------------

  implicit none

  character(len=*), intent(in) :: file
  integer, intent(in)          :: line

  write(0,'(A,I0,A,A,A,A,A,I0,A)') '[',mpi_node(),'] (*debug*) ', trim(dbg_buf__), &
                  " [", trim(adjustl(file)), &
                  ":", line, "]"

end subroutine dbg__
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine asrt__( test, message, file, line )
!---------------------------------------------------------------------------------------------------
!  Assertion routine. If test is false issue error message to screen and and stop program
!  Used by the ASSERT macro
!---------------------------------------------------------------------------------------------------

  implicit none

  logical, intent(in)          :: test
  character(len=*), intent(in) :: message
  character(len=*), intent(in) :: file
  integer, intent(in)          :: line

  if (.not. test) then
    call abort_program( p_err_assert, 'Assertion failed: '//message, file, line )
  endif

end subroutine asrt__
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!       abort the paralell program
!---------------------------------------------------------------------------------------------------
subroutine abort_program( errorcode, message, file, line )

  implicit none

  integer, optional, intent(in) :: errorcode
  character(len=*), intent(in), optional :: message
  character(len=*), intent(in), optional :: file
  integer, intent(in), optional :: line

  !       local variables
  integer :: ierr, lerrorcode, node
  character( len = 256 ) :: lmessage

  if(present(errorcode)) then
  lerrorcode = errorcode
  else
    lerrorcode = -1
  endif

  lmessage = ""

  if (present(line)) write(lmessage,*) " line: ", line
  if (present(file)) lmessage = " file: "// trim(adjustl(file)) // trim(lmessage)
  if (present(message)) lmessage = trim(message) // trim(lmessage)

  node = mpi_node()

  if (lmessage /= '') write(0,'(A)') " (*error*) ", lmessage
  write(0,'(A,I0,A,I0)') " (*error*) abort command received on node ", node, &
         ", error code: ", lerrorcode
  write(0,'(A)') ' (*error*) shutting down.'

  if (os_mpi_started) then
  call MPI_ABORT( mpi_comm_world, lerrorcode, ierr )
  endif
  call exit(lerrorcode)

end subroutine abort_program
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! memcpy routines
!
! Syntax is:
!
! call memcpy_c( to, from, size )
!
! where:
!    to   - Array to copy to
!    from - Array to copy from
!    size - Number of elements (not bytes) to copy
!---------------------------------------------------------------------------------------------------
subroutine memcpy_1d_int( s1, s2, n )

  implicit none

  integer, dimension(:), intent(inout), target :: s1
  integer, dimension(:), intent(in), target :: s2
  integer, intent(in) :: n

  integer :: tmp
  call memcpy_c( c_loc(s1), c_loc(s2), n * c_sizeof(tmp) )

end subroutine memcpy_1d_int

subroutine memcpy_2d_int( s1, s2, n )

  implicit none

  integer, dimension(:,:), intent(inout), target :: s1
  integer, dimension(:,:), intent(in), target :: s2
  integer, intent(in) :: n

  integer :: tmp
  call memcpy_c( c_loc(s1), c_loc(s2), n * c_sizeof(tmp) )

end subroutine memcpy_2d_int

subroutine memcpy_3d_int( s1, s2, n )

  implicit none

  integer, dimension(:,:,:), intent(inout), target :: s1
  integer, dimension(:,:,:), intent(in), target :: s2
  integer, intent(in) :: n

  integer :: tmp
  call memcpy_c( c_loc(s1), c_loc(s2), n * c_sizeof(tmp) )

end subroutine memcpy_3d_int

subroutine memcpy_4d_int( s1, s2, n )

  implicit none

  integer, dimension(:,:,:), intent(inout), target :: s1
  integer, dimension(:,:,:), intent(in), target :: s2
  integer, intent(in) :: n

  integer :: tmp
  call memcpy_c( c_loc(s1), c_loc(s2), n * c_sizeof(tmp) )

end subroutine memcpy_4d_int

subroutine memcpy_1d_r4( s1, s2, n )

  implicit none

  real(p_single), dimension(:), intent(inout), target :: s1
  real(p_single), dimension(:), intent(in), target :: s2
  integer, intent(in) :: n

  real(p_single) :: tmp
  call memcpy_c( c_loc(s1), c_loc(s2), n * c_sizeof(tmp) )

end subroutine memcpy_1d_r4

subroutine memcpy_2d_r4( s1, s2, n )

  implicit none

  real(p_single), dimension(:,:), intent(inout), target :: s1
  real(p_single), dimension(:,:), intent(in), target :: s2
  integer, intent(in) :: n

  real(p_single) :: tmp
  call memcpy_c( c_loc(s1), c_loc(s2), n * c_sizeof(tmp) )

end subroutine memcpy_2d_r4

subroutine memcpy_3d_r4( s1, s2, n )

  implicit none

  real(p_single), dimension(:,:,:), intent(inout), target :: s1
  real(p_single), dimension(:,:,:), intent(in), target :: s2
  integer, intent(in) :: n

  real(p_single) :: tmp
  call memcpy_c( c_loc(s1), c_loc(s2), n * c_sizeof(tmp) )

end subroutine memcpy_3d_r4

subroutine memcpy_4d_r4( s1, s2, n )

  implicit none

  real(p_single), dimension(:,:,:,:), intent(inout), target :: s1
  real(p_single), dimension(:,:,:,:), intent(in), target :: s2
  integer, intent(in) :: n

  real(p_single) :: tmp
  call memcpy_c( c_loc(s1), c_loc(s2), n * c_sizeof(tmp) )

end subroutine memcpy_4d_r4

subroutine memcpy_3d_1d_r4( s1, s2, n )

  implicit none

  real(p_single), dimension(:,:,:), intent(inout), target :: s1
  real(p_single), dimension(:), intent(in), target :: s2
  integer, intent(in) :: n

  real(p_single) :: tmp
  call memcpy_c( c_loc(s1), c_loc(s2), n * c_sizeof(tmp) )

end subroutine memcpy_3d_1d_r4

subroutine memcpy_1d_r8( s1, s2, n )

  implicit none

  real(p_double), dimension(:), intent(inout), target :: s1
  real(p_double), dimension(:), intent(in), target :: s2
  integer, intent(in) :: n

  real(p_double) :: tmp
  call memcpy_c( c_loc(s1), c_loc(s2), n * c_sizeof(tmp) )

end subroutine memcpy_1d_r8

subroutine memcpy_2d_r8( s1, s2, n )

  implicit none

  real(p_double), dimension(:,:), intent(inout), target :: s1
  real(p_double), dimension(:,:), intent(in), target :: s2
  integer, intent(in) :: n

  real(p_double) :: tmp
  call memcpy_c( c_loc(s1), c_loc(s2), n * c_sizeof(tmp) )

end subroutine memcpy_2d_r8

subroutine memcpy_3d_r8( s1, s2, n )

  implicit none

  real(p_double), dimension(:,:,:), intent(inout), target :: s1
  real(p_double), dimension(:,:,:), intent(in), target :: s2
  integer, intent(in) :: n

  real(p_double) :: tmp
  call memcpy_c( c_loc(s1), c_loc(s2), n * c_sizeof(tmp) )

end subroutine memcpy_3d_r8

subroutine memcpy_4d_r8( s1, s2, n )

  implicit none

  real(p_double), dimension(:,:,:,:), intent(inout), target :: s1
  real(p_double), dimension(:,:,:,:), intent(in), target :: s2
  integer, intent(in) :: n

  real(p_double) :: tmp
  call memcpy_c( c_loc(s1), c_loc(s2), n * c_sizeof(tmp) )

end subroutine memcpy_4d_r8


subroutine memcpy_3d_1d_r8( s1, s2, n )

  implicit none

  real(p_double), dimension(:,:,:), intent(inout), target :: s1
  real(p_double), dimension(:), intent(in), target :: s2
  integer, intent(in) :: n

  real(p_double) :: tmp
  call memcpy_c( c_loc(s1), c_loc(s2), n * c_sizeof(tmp) )

end subroutine memcpy_3d_1d_r8

!---------------------------------------------------------------------------------------------------
! setnan routines
!---------------------------------------------------------------------------------------------------

subroutine setnan_1d_r4( s1 )

  implicit none

  real(p_single), dimension(:), intent(inout), target :: s1

  call setnan_r4( C_LOC(s1), size(s1) )

end subroutine setnan_1d_r4

subroutine setnan_2d_r4( s1 )

  implicit none

  real(p_single), dimension(:,:), intent(inout), target :: s1

  call setnan_r4( C_LOC(s1), size(s1) )

end subroutine setnan_2d_r4

subroutine setnan_3d_r4( s1 )

  implicit none

  real(p_single), dimension(:,:,:), intent(inout), target :: s1

  call setnan_r4( C_LOC(s1), size(s1) )

end subroutine setnan_3d_r4

subroutine setnan_4d_r4( s1 )

  implicit none

  real(p_single), dimension(:,:,:,:), intent(inout), target :: s1

  call setnan_r4( C_LOC(s1), size(s1) )

end subroutine setnan_4d_r4

subroutine setnan_1d_r8( s1 )

  implicit none

  real(p_double), dimension(:), intent(inout), target :: s1

  call setnan_r8( C_LOC(s1), size(s1) )

end subroutine setnan_1d_r8

subroutine setnan_2d_r8( s1 )

  implicit none

  real(p_double), dimension(:,:), intent(inout), target :: s1

  call setnan_r8( C_LOC(s1), size(s1) )

end subroutine setnan_2d_r8

subroutine setnan_3d_r8( s1 )

  implicit none

  real(p_double), dimension(:,:,:), intent(inout), target :: s1

  call setnan_r8( C_LOC(s1), size(s1) )

end subroutine setnan_3d_r8

subroutine setnan_4d_r8( s1 )

  implicit none

  real(p_double), dimension(:,:,:,:), intent(inout), target :: s1

  call setnan_r8( C_LOC(s1), size(s1) )

end subroutine setnan_4d_r8

!-------------------------------------------------------------------------------
! Profiling checks and state changes
!-------------------------------------------------------------------------------
subroutine set_profiling(flag)

  implicit none

  logical :: flag

  p_profiling = flag

end subroutine set_profiling
!-------------------------------------------------------------------------------
function if_profiling()

  implicit none

  logical :: if_profiling

  if_profiling = p_profiling

end function if_profiling
!-------------------------------------------------------------------------------

end module m_system


#ifdef HAS_QUAD_PREC

!---------------------------------------------------------------------------------------------------
! Reduce operations for quad dataytypes. These must be defined after the end of the
! of the module so that they are defined as "external"
!---------------------------------------------------------------------------------------------------

! sum operation for quad datatypes
function sum_quad( invec, inoutvec, len, datatype )

  implicit none

  integer :: sum_quad
  real(p_single), intent(in) :: len
  real( p_quad ), dimension(:), intent(in) :: invec
  real( p_quad ), dimension(:), intent(inout) :: inoutvec
  real(p_single), intent(in) :: datatype

  integer :: i

  do i = 1, len
    inoutvec(i) = invec(i) + inoutvec(i)
  enddo

  sum_quad = 0

end function sum_quad

! max operation for quad datatypes
function max_quad( invec, inoutvec, len, datatype )

  implicit none

  integer :: max_quad
  real(p_single), intent(in) :: len
  real( p_quad ), dimension(:), intent(in) :: invec
  real( p_quad ), dimension(:), intent(inout) :: inoutvec
  real(p_single), intent(in) :: datatype

  integer :: i

  do i = 1, len
    if ( inoutvec(i) < invec(i) ) inoutvec(i) = invec(i)
  enddo

  max_quad = 0

end function max_quad

# endif
