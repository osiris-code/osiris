!#define DEBUG_FILE 1

program osiris

#include "os-config.h"
#include "os-preprocess.fpp"
#include "memory/memory.h"

  use m_system
  use m_simulation
  use m_node_conf
  use m_time_step
  use m_logprof

  implicit none

  class ( t_simulation ), pointer :: sim => null()
  type( t_options ) :: opts

  ! Initialize timers
  call timer_init()

  ! Initialize runtime options and store wall clock time
  call init_options( opts )

  ! Initialize system code (this also initializes MPI)
  call system_init( opts )

#ifdef __USE_PAPI__
  ! Initialize PAPI library and counters
  call init_papi()
#endif

  ! Initialize dynamic memory system
  call init_mem()

  ! Print Initialization banner
  call print_init_banner()

  ! Read input file and setup simulation structure
  call initialize( opts, sim )

  ! Print main algorithm parameters
  if ( root(sim%no_co) ) then
    call sim%list_algorithm()
  endif

  ! Print dynamic memory status
  call status_mem( comm( sim%no_co ), vmrss=.false., vmsize=.false. )

  ! write information about this run's configuration to directory of the executable
  call print_run_info_to_file( sim, '.' )

  ! Run simulation
  call run_sim( sim )

  ! Write detailed timing information to disk
  call list_total_event_times( n(sim%tstep), comm(sim%no_co), label = 'final' )

  ! write information about this run's configuration to directory holding diagnostics
  !   if any were created during simulation.
  call print_run_info_to_file( sim, 'MS' )

  ! Cleanup and shut down
  call sim%cleanup( )
  deallocate( sim )

  ! Print final message and timings overview
  if ( mpi_node() == 0 ) then
    print *, ''
    print *, 'Osiris run completed normally'
  endif

  ! Finalize dynamic memory system
  call finalize_mem()

  ! this calls mpi_finalize
  call system_cleanup()

!-----------------------------------------------------------------------------------------

contains

subroutine set_diag_sim_info( opts, sim )

  use m_diagnostic_utilities
  use m_space
  use m_grid_define
  use m_parameters

  implicit none

  type( t_options ), intent(in) :: opts
  class ( t_simulation ), intent(in) :: sim
  integer :: i

  ! The preprocessor gets confused if you try to concatenate these macros
  ! so define proper Fortran parameters
  character(len=*), parameter :: date = __DATE__
  character(len=*), parameter :: time = __TIME__


#ifdef OS_REV
  sim_info % git_version  = OS_REV
#else
  sim_info % git_version  = "n/a"
#endif
  sim_info % input_file_crc32 = 0
  sim_info % compile_time = date // " " // time
  sim_info % input_file   = trim( opts%input_file )

  sim_info % dt = dt( sim % tstep )

  sim_info % ndims = sim % grid % x_dim

  do i = 1, p_x_dim
    sim_info % box_nx(i)    = sim % grid % g_nx(i)
    sim_info % box_min(i)   = xmin( sim % g_space, i )
    sim_info % box_max(i)   = xmax( sim % g_space, i )
    sim_info % periodic(i)  = sim % no_co % ifpr_(i)
    sim_info % move_c(i)    = if_move( sim % g_space, i )
    sim_info % node_conf(i) = sim % grid % nnodes(i)
  enddo

  call alloc( sim_info%nx_x1_node, [sim % grid % nnodes(1)] )
  call sim % grid % get_part_width( sim_info%nx_x1_node, 1 )

  if (sim_info % ndims >= 2) then
    call alloc(sim_info%nx_x2_node, [ sim % grid % nnodes(2)])
    call sim % grid % get_part_width( sim_info%nx_x2_node, 2 )

    if (sim_info % ndims == 3) then
      call alloc(sim_info%nx_x3_node, [sim % grid % nnodes(3)])
      call sim % grid % get_part_width( sim_info%nx_x3_node, 3 )
    endif
  endif

  ! the root node will compute the crc32 checksum of the input deck,
  if ( root(sim%no_co) ) then
    sim_info % input_file_crc32 = compute_crc32( opts%input_file )
  endif
end subroutine set_diag_sim_info

!-----------------------------------------------------------------------------------------
!
subroutine initialize( opts, sim )

  use m_input_file
  use m_parameters
  use m_simulation_factory

  implicit none

  type( t_options ), intent(inout) :: opts
  class( t_simulation ), pointer :: sim

  class( t_input_file ), pointer :: input_file

  allocate( input_file )

  if ( opts%test ) then
    print '(A,A)', 'Testing input file: ', trim(opts%input_file)
  endif

  ! Open input file
  SCR_ROOT('')
  SCR_ROOT('Reading input file, full name: ', trim( opts%input_file ) )
  call open_input( input_file, trim( opts%input_file ), mpi_world_comm() )
  SCR_ROOT('')

  ! Read global options
  SCR_ROOT('Reading global simulation options... ')
  call read_sim_options( opts, input_file )

  ! Allocate simulation object based on algorithm
  sim => create( opts )

  ! Store input file in simulation object
  allocate( sim%input_file, source=input_file )
  deallocate( input_file )

  ! Allocate polymorphic object
  call sim % allocate_objs( )

  ! Copy global options to sim object
  sim%options = opts

  ! Read in the rest of the input file
  call sim%read_input( sim%input_file )

  ! Test Courant condition and parallel configuration and other input parameters
  call sim%test_input()

  if ( opts%test ) then
     print '(A)', '35711 - Input file reads ok!'
     stop
  endif

  ! finished reading information from input file, setup simulation data
  call sim%init(sim%min_gc_sim())

  ! store simulation info in system module
  call set_diag_sim_info( opts, sim )

end subroutine initialize
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!
subroutine run_sim( sim )

  use m_time
  use m_restart
  use m_dynamic_loadbalance

  implicit none

  class( t_simulation ), intent(inout) :: sim

  ! --

  if ( root(sim%no_co) ) then
    print *, ''
    print *, 'Starting main simulation loop:'
    print *, ''
  endif

  ! synchronize all nodes
  DEBUG("Sychronizing nodes before main simulation loop")
  call sim%no_co % barrier()

  do
    ! check if profiling should be added
    if( (n(sim%tstep) >= sim%options%prof_start) .and. &
        (n(sim%tstep) <= sim%options%prof_end) ) then
      call profiling(.true.)
    else
      call profiling(.false.)
    endif

     ! check if wall clock time limit has been exceeded
     if ( wall_clock_over_limit( sim%options, n(sim%tstep), sim%no_co ) ) then
       SCR_ROOT( '')
       SCR_ROOT( ' Wall clock limit reached, stopping gracefully')
       SCR_ROOT( '')
       if ( sim%options%wall_clock_checkpoint ) then
         ! write checkpoint data
         sim%restart%wall_clock_limit = .true.
         call write_restart( sim )
       endif
       exit
     endif

     ! If simulation finished exit the loop
     if ( .not. test( sim%time ) ) then
       SCR_ROOT( '')
       SCR_ROOT( ' Simulation completed.')
       SCR_ROOT( '')
       exit
     endif

     ! Do 1 iteration
     call sim%iter()

     ! do any per-iteration maintenance
     call sim%iter_finished()

     ! write restart files if necessary
     if ( if_restart_write( sim%restart, n(sim%tstep), ndump(sim%tstep), &
                                    comm(sim%no_co), sim%no_co%ngp_id() ) ) then
        call write_restart( sim )
     endif

     ! Simulation timings
     if ( ( n(sim%tstep) > 0 ) .and. ( sim%options%ndump_prof > 0 ) ) then
       if ( mod( n(sim%tstep), sim%options%ndump_prof ) == 0) then
          call list_total_event_times( n(sim%tstep), comm(sim%no_co) )
       endif
     endif

     ! Report on computational load
     DEBUG('Before report load')
     call sim % report_load()

     ! Repartition simulation to get optimal load balance
     DEBUG('Before dynamic load balance')
     call sim % dlb()


  enddo

end subroutine run_sim
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Write checkpoint information
!-----------------------------------------------------------------------------------------
subroutine write_restart( sim )

  use m_parameters
  use m_restart
  use m_logprof

  implicit none

  class ( t_simulation ), intent(inout) :: sim

  type( t_restart_handle )    ::  restart_handle

  call begin_event(restart_write_ev)

  if ( mpi_node() == 0 ) then
    print *, ''
    print *, ' Writing checkpoint information for timestep n =', n(sim%tstep)
    print *, ''
  endif

#ifdef __RST_SION__

  ! Measure expected file size, only required by the sion library
  restart_handle%if_acc_size = .true.
  restart_handle%data_size = 0
  call sim%write_checkpoint( restart_handle )
  restart_handle%if_acc_size = .false.

#endif

  ! Open checkpoint files
  call restart_write_open( sim%restart, comm(sim%no_co), sim%no_co%ngp_id(), &
                           n(sim%tstep), ndump(sim%tstep), &
                           file_id_rst, restart_handle )

  ! Write checkpoint data
  call sim%write_checkpoint( restart_handle )

  ! Close checkpoint files
  call restart_write_close( sim%restart, comm(sim%no_co), sim%no_co%ngp_id(), &
                            n(sim%tstep), restart_handle )

  call end_event(restart_write_ev)

end subroutine write_restart

!-----------------------------------------------------------------------------------------
! Set runtime options using command line arguments and/or environment variables
!-----------------------------------------------------------------------------------------
subroutine init_options( options )

  use m_parameters

  implicit none

  type(t_options), intent(inout) :: options

  integer :: idx, ierr
  character :: opt
  character(len=p_max_filename_len) :: arg

  !
  options%wall_clock_start = timer_cpu_seconds()

  ! Parse command line options
  options%work_dir   = ""
  options%input_file = p_default_input_file

  idx = 1
  do
    call getopt( 'vrtw:', idx, opt, arg )
    if ( opt == achar(0) ) exit

    select case( opt )
      case('r')
        options%restart = .true.
      case('t')
        options%test = .true.
      case('w')
        print *, 'Setting work_dir to "',trim(arg), '"'
        options%work_dir = trim(arg)
      case('v')
        call print_version_info(6)
        stop
      case default
        print *, '(*error*) Invalid command line parameter'
        call print_usage()
        stop
    end select

  enddo

  ! Check if more than 1 argument was given
  if ( command_argument_count() > idx ) then
     print *, '(*error*) Invalid command line parameters'
     print *, '(*error*) Only one input file may be specified'
     call print_usage()
     stop
  endif

  ! Get optional input file name
  if ( command_argument_count() == idx ) then
    call get_command_argument( idx, options%input_file )
  endif

  ! Parse environment variables
  call get_environment_variable( p_env_osiris_restart, status = ierr)
  if ( ierr == 0 ) options%restart = .true.

  call get_environment_variable( p_env_osiris_test, status = ierr)
  if ( ierr == 0 ) options%test = .true.

  call get_environment_variable( p_env_osiris_wdir, value = arg, status = ierr)
  if ( ierr == 0 ) options%work_dir = trim(arg)

end subroutine init_options
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Read global options from the input file
!-----------------------------------------------------------------------------------------
subroutine read_sim_options( options, input_file )

   use stringutil
   use m_parameters
   use m_input_file
   use m_random
   use m_simulation_factory
   use m_diagnostic_utilities

   implicit none

   type( t_options ), intent(inout) :: options
   class( t_input_file ), intent(inout) :: input_file

   integer :: random_seed                ! seed for random number generator
   character(len=20) :: random_algorithm ! algorithm to use for global random number generation
   logical :: enforce_rng_constancy      ! If false, an identical iput deck will give diffrent result due to rnadom number generation 
                                         !   depending on a given run's spatial dcomposition (e.g. different number/load of MPI NODES)
                                         ! If true, an identical input deck identical results (within roundoff) no matter the run's MPI setup
                                         !  (i.e. running a simulation on 1 node will give same results no matter how many MPI nodes used. )
                                         !  This is comes at a small perforance cost, but is usefull for testing, debugging the and subtraction-technique

   character(len=20) :: algorithm        ! standard, hd-hybrid
   real(p_double) :: omega_p0, n0        ! reference plasma frequency and/or density
   real(p_double) :: gamma              ! relativistic factor for boosted frame

   integer :: ndump_prof                ! frequency for profiling information dumps

   character(len=16):: wall_clock_limit ! Wall clock limit time (h:m:s)
   integer :: wall_clock_check          ! Frequency at which to check if wall time limit
                                        ! has been reached

   logical :: wall_clock_checkpoint     ! Dump checkpointing information when stopping due to
                                        ! wall clock limit
  integer :: prof_start      ! time step to begin profiling
  integer :: prof_end        ! time step to end profiling

  character(len=20) :: file_format
  character(len=20) :: parallel_io

   integer :: ierr

   ! namelist for reading
   namelist /nl_simulation/ random_seed, random_algorithm, algorithm, omega_p0, n0, gamma, ndump_prof, &
                            wall_clock_limit, wall_clock_check, wall_clock_checkpoint, &
    prof_start, prof_end, file_format, parallel_io, enforce_rng_constancy

   algorithm   = "standard"

   random_seed = 0
   random_algorithm = 'default'
   enforce_rng_constancy = .false.
   
   omega_p0    = 0
   n0          = 0
   gamma       = 1
   ndump_prof  = 0

   wall_clock_limit        = ""
   wall_clock_check      = -1
   wall_clock_checkpoint = .true.

  prof_start = -1
  prof_end   = huge(prof_start)

  ! If the code was compiled with hdf5 default to that format
#ifdef HDF5
  file_format = 'hdf5'
#else
  file_format = 'zdf'
#endif

  parallel_io = 'mpi'

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_simulation", ierr )

   if ( ierr == 0 ) then
     read (input_file%nml_text, nml = nl_simulation, iostat = ierr)
     if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error reading global simulation parameters"
        write(0,*) "aborting..."
      endif
       stop
     endif
   else
     if (ierr < 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error reading global simulation parameters"
        write(0,*) "aborting..."
      endif
       stop
     else
       SCR_ROOT(" - global simulation parameters missing, using defaults")
     endif
   endif

   ! random number generator options
   options%random_seed  = random_seed
   options%random_class = random_type( random_algorithm )
   options%enforce_rng_constancy = enforce_rng_constancy

   if ( options%random_class < 0 ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading global simulation parameters"
      write(0,*) "invalid random_algorithm, aborting..."
    endif
    stop
   endif

   ! 'get_algorithm_from_string' is in os-simulation-factory 
   options%algorithm = get_algorithm_from_string( algorithm )

   ! set reference plasma frequency

   ! default undefined values
   options%omega_p0 = -1
   options%n0 = -1

   if ( omega_p0 < 0. ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading global simulation parameters"
      write(0,*) "Invalid value for omega_p0, must be > 0, aborting..."
    endif
      stop
   else
      options%omega_p0 = omega_p0

       ! cgs units
      options%n0 = 3.142077870426918d-10 * omega_p0**2   ! n0 in [cm^-3]
   endif

   ! set reference plasma frequency using plasma density ( this overrides setting
   ! the plasma frequency directly )
   if ( n0 < 0. ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading global simulation parameters"
      write(0,*) "Invalid value for n0, must be > 0, aborting..."
    endif
      stop
   else if ( n0 > 0. ) then

      options%n0 = n0

      ! cgs units
      options%omega_p0 = 56414.60192463558d0*sqrt(n0)     ! n0 in [cm^-3]

   endif

   if ( gamma < 1. ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading global simulation parameters"
      write(0,*) "Invalid value for gamma, must be > 1, aborting..."
    endif
      stop
   else
      options%gamma = gamma
   endif

   if ( ndump_prof > 0 ) then
     options%ndump_prof = ndump_prof
   else
     options%ndump_prof = 0
   endif

   if ( wall_clock_limit /= "" ) then
     options%wall_clock_limit = time2seconds( wall_clock_limit )
     if ( options%wall_clock_limit < 0 ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error reading global simulation parameters"
        write(0,*) "Invalid value for wall_clock_limit, must be in the form h:m:s aborting..."
      endif
        stop
     endif
     options%wall_clock_check      = wall_clock_check
     options%wall_clock_checkpoint = wall_clock_checkpoint
   endif

  if( prof_end >= prof_start ) then
    options%prof_start = prof_start
    options%prof_end   = prof_end
  else
    if ( mpi_node() == 0 ) then
      write(0,*) "Error in profiling parameters, requires prof_end >= prof_start"
      write(0,*) "aborting..."
    endif
    stop
  endif

  select case( trim( file_format ) )
  case( "zdf" )
    options%file_format = p_zdf_format
  case( "hdf5" )
    options%file_format = p_hdf5_format
  case default
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading global simulation parameters"
      write(0,*) "file_format must be one of 'zdf' or 'hdf5', aborting..."
    endif
    stop
  end select

  select case( trim( file_format ) )
  case( "zdf" )
    options%file_format = p_zdf_format
  case( "hdf5" )
    options%file_format = p_hdf5_format
  case default
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading global simulation parameters"
      write(0,*) "file_format must be one of 'zdf' or 'hdf5', aborting..."
    endif
    stop
  end select

  select case( trim( file_format ) )
  case( "zdf" )
    options%file_format = p_zdf_format
  case( "hdf5" )
#ifdef HDF5
    options%file_format = p_hdf5_format
#else
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading global simulation parameters"
      write(0,*) "Invalid file_format, the code was not compiled with HDF5 support, aborting..."
    endif
    stop
#endif
  case default
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading global simulation parameters"
      write(0,*) "file_format must be one of 'zdf' or 'hdf5', aborting..."
    endif
    stop
  end select

  select case( trim( parallel_io ) )
  case( "mpi" )
    options%parallel_io = DIAG_MPI
  case( "indep" )
    if ( options%file_format ==  p_hdf5_format ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error in global simulation parameters"
        write(0,*) "parallel_io type 'indep' not supported by 'hdf5' file format, aborting..."
      endif
      stop
  endif
    options%parallel_io = DIAG_INDEPENDENT
  case( "mpiio-indep" )
    options%parallel_io = DIAG_MPIIO_INDEPENDENT
  case( "mpiio-coll" )
    options%parallel_io = DIAG_MPIIO_COLLECTIVE
  case default
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading global simulation parameters"
      write(0,*) "parallel_io must be one of 'mpi', 'indep', 'mpiio-indep' or 'mpiio-coll', aborting..."
    endif
    stop
  end select

end subroutine read_sim_options
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Check wall time
!-----------------------------------------------------------------------------------------
function wall_clock_over_limit( options, n, no_co )

  implicit none

  type( t_options ), intent(in) :: options
  integer, intent(in) :: n
  class( t_node_conf ), intent(in) :: no_co

  logical :: wall_clock_over_limit
  integer :: t

  if ( options%wall_clock_check > 0 ) then

    if ( (n > 0) .and. (mod( n, options%wall_clock_check) == 0) ) then

      ! Check if wall time was exceeded on root node
      if ( root( no_co ) ) then

        ! Compare elapsed time with limit
        if (  timer_cpu_seconds() >= options%wall_clock_start + options%wall_clock_limit ) then
          t = 1
        else
          t = 0
        endif
      endif

      ! Broadcast result to all nodes
      call broadcast( no_co, t )

      ! set the function value on all nodes
      wall_clock_over_limit = ( t == 1 )
    else
      wall_clock_over_limit = .false.
    endif
  else
    wall_clock_over_limit = .false.
  endif

end function wall_clock_over_limit
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Print Initial banner and code revision
!-----------------------------------------------------------------------------------------
subroutine print_init_banner()

  implicit none

  character(len = 6) :: prec

    if ( mpi_node() == 0 ) then
       print *, ''
       print *, '***************************************************************'

      ! SIMD Code
#ifdef SIMD

#if defined( PRECISION_SINGLE )
      prec = 'single'
#elif defined( PRECISION_DOUBLE )
      prec = 'double'
#else
      ! This should never happen, it should be caught by the configure script
      prec = 'unknow'
#endif

#if defined(SIMD_SSE)
      print *, '*          Using ', prec, ' precision SSE optimized Code          *'
#elif defined(SIMD_AVX)
      print *, '*          Using ', prec, ' precision AVX optimized Code          *'
#elif defined(SIMD_AVX2)
      print *, '*          Using ', prec, ' precision AVX2 optimized Code          *'
#elif defined(SIMD_BGQ)
      print *, '*        Using ', prec, ' precision BlueGene/Q optimized Code     *'
#elif defined(SIMD_MIC)
      print *, '*         Using ', prec, ' precision Intel MIC optimized Code     *'
#endif


#else
       print *, '*                     Using Fortran version                   *'
#endif
       print *, '***************************************************************'
       print *, ''


#ifdef OS_REV
       ! Print revision version is available
       print *, "Software revision: ", OS_REV
       print *, ''
#endif

    endif



end subroutine print_init_banner
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Print help on command line parameters
!-----------------------------------------------------------------------------------------
subroutine print_usage()

  implicit none

  print *, 'Usage: osiris [-t] [-r] [-w work_dir] [input_file]'
  print *, ''
  print *, 'Arguments:'
  print *, 'input_file    -   Input file to use. If not specified defaults to "os-stdin"'
  print *, ''
  print *, 'Options:'
  print *, ''
  print *, '-t            -   Test only. If set osiris will only test the validity of the'
  print *, '                  specified file, and all other options are ignored.'
  print *, '-r            -   Restart. If set it will force osiris to restart from'
  print *, '                  checkpoint data.'
  print *, '-w work_dir   -   Work directory. If set osiris will change to this directory'
  print *, '                  before starting the simulation.'
  print *, '-v            -   Print version information and exit (no simulation run).'

end subroutine print_usage
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Print out the build parameters used to create this particular Osiris executable
!-----------------------------------------------------------------------------------------
subroutine print_version_info(output_unit)

  implicit none

  integer :: output_unit

  ! '__generated_config_info.f90' is created by the Makefile and contains a series of 'write' statments that print out
  !   various build parameters. See 'Makefile.config-info'
  #include "__generated_config_info.f90"

end subroutine

subroutine print_run_info(output_unit, input_deck_str)
  use m_diagnostic_utilities, only: sim_info
  implicit none
  integer :: output_unit
  character(len=*), intent(in) :: input_deck_str

  write(output_unit,'(A)')         'run_info = {'
  write(output_unit,'(A,F20.0,A)') '  "input_deck_crc":', (sim_info%input_file_crc32), ','
  write(output_unit,'(A)')         '  "input_deck": """'
  write(output_unit,'(A)')         input_deck_str
  write(output_unit,'(A)')         '"""'
  write(output_unit,'(A)')         '}'
end subroutine
!-----------------------------------------------------------------------------------------
! Print out the build parameters used to create this particular Osiris executable 
!    to a file named 'build-info'
!-----------------------------------------------------------------------------------------
subroutine print_run_info_to_file( sim, subdir )

  ! use m_system,               only: dir_exists
  ! use m_input_file,           only: p_max_input_file_size
  use m_diagnostic_utilities, only: sim_info

  implicit none
  class(t_simulation) :: sim
  character(len=*)    :: subdir

  integer :: ierr
  character(len=:), allocatable :: input_deck_str
  integer :: istat, filesize

  if ( root(sim%no_co) ) then 

    ! -----
    ! read in the input deck using 'stream access'
    !   so that we can add it the 'run-info' output file.
    ! -----
    open(unit=67,file=trim(sim_info%input_file),status='OLD',&
        form='UNFORMATTED',access='STREAM',iostat=istat)
    if( istat .eq. 0 ) then
      inquire(file=trim(sim_info%input_file), size=filesize)
      allocate( character(len=filesize) :: input_deck_str )
      read(67,pos=1,iostat=istat) input_deck_str
      close(67, iostat=istat)
    endif
    ! if something happened and we couldn't read the input deck
    !  then allocated a 1 character string and set to ' ' as a placeholder.
    if (.not. ALLOCATED(input_deck_str)) then
      allocate( character(len=1) :: input_deck_str )
      input_deck_str = ' '
    endif
    close(67, iostat=istat) 

    ! check to see the passed-in directory exists
    !  If so, place a copy of 'run-info' in there.
    if( dir_exists(subdir) ) then
      open(unit = 66, iostat = ierr, file = trim(subdir)// p_dir_sep //"run-info", action = "write")
      if ( ierr /= 0 ) then
        ! do nothing.. not major error if we can't write this.
      endif

      call print_version_info( 66 )
      call print_run_info( 66, input_deck_str )
      close(66)    
    endif
    deallocate( input_deck_str )
  endif

end subroutine

end program osiris
! -----------------------------------------------------------------
