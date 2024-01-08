!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     restart operation class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#include "os-config.h"
#include "os-preprocess.fpp"

!-----------------------------------------------------------------------------------------
! Set the default restart type. Possible values are:
!   __RST_BINARY__  : Use the standard fortran binary output
!   __RST_SION__  : Use the SION library (not available on all systems)
!   __RST_HDF5__  : Use hdf5 output (not implemented yet)
!-----------------------------------------------------------------------------------------

#if (! defined(__RST_BINARY__)) &&  (! defined(__RST_SION__))

! default to binary output
#define __RST_BINARY__

#endif

!-----------------------------------------------------------------------------------------


module m_restart

  use m_parameters
  use stringutil
  use m_diagfile, only: p_time_length

  implicit none

  ! restrict access to things explicitly declared public
  private

  ! string to id restart data
  character(len=*), parameter :: p_rest_rst_id = "restart rst data - 0x0002"
  ! string length for time/node id
  integer, parameter :: p_node_length = 5

  type :: t_restart

    ! frequency of restart dumps (set to -1 to use timed restart)
    integer :: ndump_fac

    ! switch to restart from a restart file
    logical :: if_restart

    ! switch to remove the previous restart file when finished writing
    ! the new one
    logical :: if_remold

    ! frequency of restart dumps in wall clock time (seconds)
    real(p_double) :: ndump_time

    ! time at wich the last restart was done (seconds)
    real(p_double) :: last_restart_time

    ! dump number for timed dumps
    integer :: time_dump_number

    ! force writing a restart at a given iteration
    integer :: debug_iter
  
    character(len = 80) :: old_file_name = ''

    ! if writing restart from hitting a wall clock limit
    logical :: wall_clock_limit = .false.
 
  end type t_restart

  type  :: t_restart_handle

    ! should restart write only accumulate the data size that would be written?
    logical   ::   if_acc_size

    ! Restart file size calculation (needs to be a 64bit integer)
    integer(p_int64)  ::  data_size

    ! file handle for restart file
#if defined( __RST_BINARY__ )
    integer  ::  file_id
#elif defined(  __RST_SION__ )
    integer  ::  sid
#endif

    ! filename currently open
    character(len = 80) :: current_file_name = ''

  end type t_restart_handle


  interface read_nml
    module procedure read_nml_restart
  end interface

  interface setup
    module procedure setup_restart
  end interface

  interface restart_write
    module procedure restart_write_restart
  end interface

  interface restart_read
    module procedure restart_read_restart
  end interface

  interface if_restart_write
    module procedure if_restart_write
  end interface

  interface if_restart_read
    module procedure if_restart_read
  end interface

  interface restart_write_open
    module procedure restart_write_open
  end interface

  interface restart_read_open
    module procedure restart_read_open
  end interface

  interface restart_write_close
    module procedure restart_write_close
  end interface

  interface restart_read_close
    module procedure restart_read_close
  end interface

  interface restart_io_write
    module  procedure restart_io_write_string
    module  procedure restart_io_write_string1d
    module  procedure restart_io_write_string2d
    module  procedure restart_io_write_logical
    module  procedure restart_io_write_logical1d
    module  procedure restart_io_write_integer
    module  procedure restart_io_write_integer1d
    module  procedure restart_io_write_integer2d
    module  procedure restart_io_write_integer3d
    module  procedure restart_io_write_integer4d
    module  procedure restart_io_write_int64
    module  procedure restart_io_write_int641d
    module  procedure restart_io_write_single
    module  procedure restart_io_write_single1d
    module  procedure restart_io_write_single2d
    module  procedure restart_io_write_single3d
    module  procedure restart_io_write_single4d
    module  procedure restart_io_write_double
    module  procedure restart_io_write_double1d
    module  procedure restart_io_write_double2d
    module  procedure restart_io_write_double3d
    module  procedure restart_io_write_double4d
  end interface

  interface restart_io_read
    module  procedure restart_io_read_string
    module  procedure restart_io_read_string1d
    module  procedure restart_io_read_string2d
    module  procedure restart_io_read_logical
    module  procedure restart_io_read_logical1d
    module  procedure restart_io_read_integer
    module  procedure restart_io_read_integer1d
    module  procedure restart_io_read_integer2d
    module  procedure restart_io_read_integer3d
    module  procedure restart_io_read_integer4d
    module  procedure restart_io_read_int64
    module  procedure restart_io_read_int641d
    module  procedure restart_io_read_single
    module  procedure restart_io_read_single1d
    module  procedure restart_io_read_single2d
    module  procedure restart_io_read_single3d
    module  procedure restart_io_read_single4d
    module  procedure restart_io_read_double
    module  procedure restart_io_read_double1d
    module  procedure restart_io_read_double2d
    module  procedure restart_io_read_double3d
    module  procedure restart_io_read_double4d
  end interface


  ! declare things that should be public
  public :: t_restart, t_restart_handle
  public :: read_nml, setup
  public :: if_restart_write, if_restart_read
  public :: restart_write_open,  restart_read_open
  public :: restart_write_close, restart_read_close
  public :: restart_write
  public :: restart_io_write, restart_io_read

contains


!---------------------------------------------------
!---------------------------------------------------
subroutine read_nml_restart( this, input_file, restart )

  use m_input_file

  implicit none

  type( t_restart ), intent(out) :: this
  class( t_input_file ), intent(inout) :: input_file
  logical, intent(in) :: restart  ! override if_restart value

  integer :: ndump_fac
  logical         :: if_restart
  logical         :: if_remold
  real(p_double) :: ndump_time
  integer :: debug_iter

  namelist /nl_restart/ ndump_fac, if_restart, &
            if_remold, ndump_time, debug_iter

  integer::ierr

  ndump_fac  = 0          ! restart information is off
  if_restart = .false.
  if_remold  = .false.
  ndump_time = 0.0
  debug_iter = -1

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_restart", ierr )

  if (ierr == 0) then
    read ( input_file%nml_text, nml = nl_restart, iostat = ierr)
    if (ierr /= 0) then
      print *, "Error reading restart parameters"
      print *, "aborting..."
      stop
    endif
  else
    if (ierr < 0) then
      write(*,*) "Error reading restart parameters"
      write(*,*) "aborting..."
      stop
    else
      if (disp_out(input_file)) then
        SCR_ROOT(" - restart parameters missing, using default")
      endif
    endif
  endif

  this%ndump_fac  = ndump_fac
  this%if_restart = if_restart
  this%if_remold = if_remold

  this%ndump_time = ndump_time
  this%debug_iter = debug_iter

  ! override input file if restart is set to true
  if ( restart ) this%if_restart = .true.

end subroutine read_nml_restart
!---------------------------------------------------


!---------------------------------------------------
subroutine setup_restart( this, restart, restart_handle)
!---------------------------------------------------
!       sets up this data structure from the given information
!---------------------------------------------------


  implicit none

  ! dummy variables

  type( t_restart ), intent( inout )  ::  this
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle


  ! executable statements
  if ((this%ndump_fac == -1) .and. &
    (this%ndump_time <= 0.0)) then
    write(*,*) "(*error*) Writing of restart information is set for timed restart"
    write(*,*) "(*error*) and the frequency is set to ",this%ndump_time," which is"
    write(*,*) "(*error*) less than 0.0"
    call abort_program()
  endif

  ! set the last restart wall time to the initial simulation wall time
  this%last_restart_time = timer_cpu_seconds()

  ! read restart values if necessary
  if ( restart ) then
    call restart_read( this, restart_handle )
  else
    this%time_dump_number = 1
  endif

end subroutine setup_restart
!---------------------------------------------------

!---------------------------------------------------
subroutine restart_write_restart( this, restart_handle )
!---------------------------------------------------
!       write object information into a restart file
!---------------------------------------------------

  implicit none

  type( t_restart ), intent(in)  ::  this
  type( t_restart_handle ), intent(inout) :: restart_handle

  integer :: ierr

  restart_io_wr( p_rest_rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then
    ERROR('error writing restart data for restart object.')
    call abort_program(p_err_rstwrt)
  endif

  ! we need to write time_dump_number+1 because
  ! this variable only gets incremented after the
  ! write_restart is done
  restart_io_wr( this%time_dump_number+1, restart_handle, ierr )
  if ( ierr/=0 ) then
    ERROR('error writing restart data for restart object.')
    call abort_program(p_err_rstwrt)
  endif

  restart_io_wr( restart_handle%current_file_name, restart_handle, ierr )
  if ( ierr/=0 ) then
    ERROR('error writing restart data for restart object.')
    call abort_program(p_err_rstwrt)
  endif

end subroutine restart_write_restart
!---------------------------------------------------

!---------------------------------------------------
subroutine restart_read_restart( this, restart_handle )
!---------------------------------------------------
!       read object information from a restart file
!---------------------------------------------------

  implicit none

  type( t_restart ), intent(inout)  ::  this
  type( t_restart_handle ), intent(in) :: restart_handle

  integer :: ierr
  character(len=len(p_rest_rst_id)) :: rst_id

  restart_io_rd( rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then
    ERROR('error reading restart data for restart object.')
    call abort_program(p_err_rstrd)
  endif

  ! check if restart file is compatible
  if ( rst_id /= p_rest_rst_id) then
    ERROR('Corrupted restart file, or restart file from incompatible binary (restart)')
    call abort_program(p_err_rstrd)
  endif

  restart_io_rd( this%time_dump_number, restart_handle, ierr )
  if ( ierr/=0 ) then
    ERROR('error writing restart data for restart object.')
    call abort_program(p_err_rstwrt)
  endif

  restart_io_rd( this%old_file_name, restart_handle, ierr )
  if ( ierr/=0 ) then
    ERROR('error writing restart data for restart object.')
    call abort_program(p_err_rstwrt)
  endif

end subroutine restart_read_restart
!---------------------------------------------------


!---------------------------------------------------
function if_restart_write( this, n, ndump, gcom, ngp_id )
!---------------------------------------------------
!       check whether to write a restart file
!---------------------------------------------------

  implicit none

  logical :: if_restart_write

  type( t_restart ), intent(in) :: this
  integer,   intent(in) :: n, ndump, gcom, ngp_id

  ! local variables

  integer :: n_del, n_aux, ierr

  integer :: overtime, sendcount, root

  ! restart dumps are turned off if ndump_fac is 0
  if ( ndump > 0) then
    if ( this%ndump_fac > 0 )  then  ! check for number of iterations restart dumps

      n_del = ndump * this%ndump_fac
      n_aux = n / n_del

      !            dump only when n is multiples of ndump*this%ndump_fac
      if_restart_write = ( ( n - n_aux*n_del ) .eq. 0 )

      !            no restart file needed for n = 0
      if ( n == 0 ) if_restart_write = .false.

    elseif ( this%ndump_fac == -1 ) then ! check for timed restart dumps

      ! only node 0 wall clock time is used
      if (ngp_id == 1) then
        if (timer_cpu_seconds()-this%last_restart_time >= this%ndump_time) then

          write(*,*) "More than ", this%ndump_time, "seconds have passed, "
          write(*,*) "wrinting restart information..."

          overtime = 1
        else
          overtime = 0
        endif
        if ( gcom /= MPI_COMM_NULL) then
          sendcount = 1
          root = 0
          call MPI_BCAST( overtime, sendcount, MPI_INTEGER, root, gcom, ierr)
        endif
      else
        if ( gcom /= MPI_COMM_NULL) then
          sendcount = 1
          root = 0
          call MPI_BCAST( overtime, sendcount, MPI_INTEGER, root, gcom, ierr)
        endif
      endif

      if_restart_write = (overtime == 1)

    else

      if_restart_write = .false.

    endif
  else

    if_restart_write = .false.

  endif

  ! allow for debug restarts i.e. restarts at a specific time step
  if_restart_write = ( if_restart_write .or. ( n == this%debug_iter ))

end function if_restart_write
!---------------------------------------------------


!---------------------------------------------------
function if_restart_read( this )
!---------------------------------------------------
!       check whether to to read the restart file
!---------------------------------------------------

  implicit none

  logical :: if_restart_read

  type( t_restart ), intent(in) :: this

  if_restart_read = this%if_restart

end function if_restart_read
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_write_open( this, comm, ngp_id, n, ndump, file_id, restart_handle )
!---------------------------------------------------
!       open new restart file for writing
!---------------------------------------------------

  implicit none

  type( t_restart ),   intent(inout) :: this
  integer, intent(in) :: comm, ngp_id
  integer, intent(in) :: n, ndump
  integer, intent(in) :: file_id
  type( t_restart_handle ), intent(inout) :: restart_handle

  integer :: n_aux
  character(len = 128) :: current_file_name
  integer :: ierr

! SIONlib variables
#ifdef __RST_SION__
  character(len = 128) :: newfname
  integer :: nfiles, fsblksize, gcomm, lcomm
  integer*8 :: chunksize
#endif

  ! prepare path and file names

  if ( n /= this%debug_iter ) then

    if ( this%wall_clock_limit ) then
      current_file_name = 'wall-'
    else
      current_file_name = 'rst-'
    endif

    if (ndump <= 0) then
      ! executing a timed restart, number restart files sequentially
      n_aux = this%time_dump_number
    else
      n_aux =  n / ndump
    endif

  else

    current_file_name = 'debug-'
    n_aux =  n

  endif

#if defined( __RST_BINARY__ )

  call mkdir( trim(path_rest), ierr )

#elif defined( __RST_SION__ )

  if ( ngp_id == 1 ) then
    call mkdir( trim(path_rest), ierr )
  else
    ierr=0
  endif

#else
#error no valid RST_IO defined!
#endif

  if (ierr /= 0) then
    ERROR('Failed to create directory ', trim(path_rest))
    ERROR('Error code ', ierr, " ", strerror(ierr))
    call abort_program(p_err_rstwrt)
  endif

#if defined( __RST_BINARY__)

  current_file_name = trim(current_file_name)//trim(idx_string( n_aux, p_time_length )) &
                      //'.n'//idx_string( ngp_id, p_node_length )
  open (unit=file_id, file=trim(path_rest)//current_file_name, &
        form='unformatted', iostat = ierr)
  restart_handle%file_id = file_id

#elif defined( __RST_SION__)

  if ( comm == MPI_COMM_NULL ) then
    write(0,*) '(*error*) Using SIONlib for checkpointing is not supported for serial runs'
    write(0,*) '(*error*) For checkpointing serial runs please recompile the code without SIONlib support'
    call abort_program( p_err_invalid)
  endif

  current_file_name = trim(current_file_name)//trim(idx_string( n_aux, p_time_length ))

  nfiles = 0
  gcomm  = comm
  lcomm  = comm
  fsblksize = -1
  chunksize = restart_handle%data_size

  call fsion_paropen_mpi( trim(path_rest)//trim(current_file_name), 'bw', nfiles, &
                          gcomm, lcomm, chunksize, fsblksize, &
                          ngp_id, newfname, restart_handle%sid )

  if ( restart_handle%sid == -1 ) then
    ierr = 1
  else
    ierr = 0
  endif

#else
#error no valid RST_IO defined!
#endif

  if (ierr /= 0) then
    ERROR('Unable to create restart file ', trim(current_file_name) )
    call abort_program(p_err_rstwrt)
  endif

  restart_handle%current_file_name = current_file_name
  restart_handle%if_acc_size = .false.

end subroutine restart_write_open
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_read_open( this, comm, ngp_id, file_id, restart_handle )
!---------------------------------------------------
!       open restart file for reading
!---------------------------------------------------

implicit none

  type( t_restart ),   intent(in) :: this
  integer, intent(in) :: comm, ngp_id
  integer,     intent(in) :: file_id
  type( t_restart_handle ), intent(out) :: restart_handle

  character(80)   :: path
  character(80)   :: full_name
  character(20)   :: ch_me

#ifdef __RST_SION__
  character(len = 128) :: newfname
  integer :: nfiles, fsblksize, gcomm, lcomm
  integer*8 :: chunksize
#endif

  integer :: ierr

  path  =  trim(path_rest)
  ch_me  = '.n'//idx_string( ngp_id, p_node_length )

#if defined( __RST_BINARY__)

  full_name = trim(path) // 'rstrtdat' // trim(ch_me)
  open (unit=file_id, file=full_name, status='old', action='read',&
      form='unformatted', iostat = ierr)
  restart_handle%file_id = file_id

#elif defined( __RST_SION__)

  full_name = trim(path) // 'rstrtdat'

  nfiles = 0
  gcomm  = comm
  lcomm  = comm
  fsblksize = -1
  chunksize = 0

  call fsion_paropen_mpi( trim(full_name), 'rb', nfiles, &
                          gcomm, lcomm, chunksize, fsblksize, &
                          ngp_id, newfname, restart_handle%sid )

  if ( restart_handle%sid == -1 ) then
    ierr = 1
  else
    ierr = 0
  endif

#else
#error no valid RST_IO defined!
#endif

  if (ierr /= 0) then
    ERROR(" Error opening restart file : ", trim(full_name) )
    call abort_program(p_err_rstrd)
  endif

  restart_handle%if_acc_size = .false.
  restart_handle%data_size = 0

end subroutine restart_read_open
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_write_close( this, gcom, ngp_id, n, restart_handle )
!---------------------------------------------------
!  close newly written restart files and make them accessible
!---------------------------------------------------

  use m_system

  implicit none

  type( t_restart ),   intent(inout) :: this
  integer, intent(in) :: gcom, ngp_id
  integer, intent(in) :: n
  type( t_restart_handle ), intent(inout) :: restart_handle

  logical   :: remold
  character(len = 80)   :: target_name, name, oldfile
  integer         :: ierr

  remold = ( (n /= this%debug_iter) .and. this%if_remold .and. &
           (this%old_file_name /= '' ))

#if defined( __RST_BINARY__)
  close( restart_handle%file_id, iostat = ierr )
#elif defined( __RST_SION__)
  call fsion_parclose_mpi(restart_handle%sid,ierr)
  ierr=ierr-1
#endif

  if ( ierr/= 0) then
    ERROR("error closing restart file")
    call abort_program( p_err_rstwrt )
  endif

  ! link rstrtdat.xxx -> rst-yyyy.xxx
  target_name    = restart_handle%current_file_name

#if defined( __RST_BINARY__)
  name      = trim(path_rest) // 'rstrtdat' // '.n' // idx_string( ngp_id, p_node_length )
  call file_link( trim(target_name), trim(name), ierr )
#elif defined( __RST_SION__)
  if ( ngp_id == 1 ) then
    name      = trim(path_rest) // 'rstrtdat'
    call file_link( trim(target_name), trim(name), ierr )
  else
    ierr = 0
  endif
#endif

  if (ierr /= 0) then
    ERROR("error creating restart file link")
    call abort_program( p_err_rstwrt )
  endif

  ! Synchronize nodes (in case some of them fail to write/link)
  if( remold ) then
    if ( gcom /= MPI_COMM_NULL ) then
      call MPI_BARRIER( gcom, ierr )
    endif
  endif

  if( remold ) then
    ! remove old file
    oldfile = trim(path_rest) // trim(this%old_file_name)
    call remove(trim(oldfile),ierr)
  endif
  if (n /= this%debug_iter) then
    this%old_file_name = restart_handle%current_file_name
  endif
  restart_handle%current_file_name = ''

  ! Synchronize nodes (in case some of them fail to write/link)
  if ( gcom /= MPI_COMM_NULL ) then
    call MPI_BARRIER( gcom, ierr )
  endif
  this%last_restart_time = timer_cpu_seconds()
  this%time_dump_number = this%time_dump_number + 1

end subroutine restart_write_close
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_read_close( this, restart_handle )
!---------------------------------------------------
!       close restart file after reading it
!---------------------------------------------------

  implicit none

  type( t_restart ), intent(in) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

#if defined( __RST_BINARY__)

  close( restart_handle%file_id )

#elif defined( __RST_SION__)

  integer :: ierr

  ! if it fails here, we already have all the data, so we won't bail out.
  call fsion_parclose_mpi(restart_handle%sid,ierr)

#else
#error no valid RST_IO defined!
#endif


end subroutine restart_read_close
!---------------------------------------------------



!---------------------------------------------------
! include i/o procedures for the appropriate file backend
!---------------------------------------------------

#if defined( __RST_BINARY__ )

#include "os-restart-io-binary.f90"

#elif defined( __RST_SION__ )

#include "os-restart-io-sion.f90"

#else

#error no valid RST_IO defined!

#endif
!---------------------------------------------------

end module m_restart





