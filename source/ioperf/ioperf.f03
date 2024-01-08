! Test program to measure IO performance!

!---------------------------------------------------------------------------------------------------

program ioperf

#include "os-config.h"
#include "os-preprocess.fpp"
#include "memory/memory.h"

use m_system
use m_parameters
use m_simulation
use m_node_conf
use m_time_step
use m_logprof

use hdf5_util
use hdf5

use m_grid_define
use m_grid
use m_grid_parallel
use m_space
use m_vdf_define
use m_input_file
use m_vdf_define

use m_restart

use stringutil

use m_diagnostic_utilities

implicit none

type :: t_test
  integer, dimension(p_x_dim) :: grid_size
  integer :: file_format
  integer :: parallel_io
  integer, dimension(p_x_dim) :: n_merge
  integer :: merge_type
end type t_test


type :: t_test_set
  type(t_test), dimension(:), pointer :: tests
  integer :: n_tests
  integer :: repeat
  logical :: keepFiles
  character(len=1024) :: path
end type t_test_set

!---------------------------------------------------------------------------------------------------
! IO Perf
!---------------------------------------------------------------------------------------------------

class( t_node_conf ) ::  no_co
type( t_options ) :: opts
type( t_test_set ) :: test_set
integer :: i
integer, parameter :: fid = 10

! Initialize timers
call timer_init()

! Initialize system code (this also initializes MPI)
call system_init( opts )
call set_diag_sim_info()

! Initialize hdf5 subsystem
call open_hdf5( )

call ReadInput( no_co, test_set )

! Initialize communications module
call setup( no_co )

! Open file for test summary
open( fid, file = "ioperf.out", status = 'replace', form = 'formatted')

! write header
select case (p_x_dim)
case(1)
  write(fid,"(A)") "   Size     |  Dim. | File |    IO    | merge  | nmrg |       time        |                Bandwidth             "
  write(fid,"(A)") "------------+-------+------+----------+--------+------+-------------------+--------------------------------------"
case(2)
  write(fid,"(A)") "   Size     | Dimensions | File |    IO    | merge  | nmerge  |       time        |                Bandwidth             "
  write(fid,"(A)") "------------+------------+------+----------+--------+---------+-------------------+--------------------------------------"
case(3)
  write(fid,"(A)") "   Size     |    Dimensions     | File |    IO    | merge  |    nmerge    |       time        |                Bandwidth             "
  write(fid,"(A)") "------------+-------------------+------+----------+--------+--------------+-------------------+--------------------------------------"
end select

! run tests
do i = 1, test_set % n_tests
  call RunTest( i, fid, test_set % tests(i), test_set % path, no_co, test_set % repeat, test_set % keepFiles )
enddo

close( fid )

! cleanup data structures
call cleanup( no_co )

call cleanup_dutil()

! Finalize dynamic memory system
call finalize_mem()

! this calls mpi_finalize
call system_cleanup()

!---------------------------------------------------------------------------------------------------

contains

!---------------------------------------------------------------------------------------------------
! Runs individual test
!---------------------------------------------------------------------------------------------------
subroutine RunTest( id, fid, test, path, no_co, repeat, keepFiles )

  use m_vdf_report
  use m_diagnostic_utilities

  implicit none

  integer, intent(in) :: id, fid
  type(t_test), intent(in) :: test
  character(len=*), intent(in) :: path
  class( t_node_conf), intent(in) :: no_co
  integer, intent(in) :: repeat
  logical, intent(in) :: keepFiles

  class(t_grid) :: grid
  type(t_space) :: space
  integer, dimension(2, p_x_dim) :: gc_num
  real(p_double), dimension(2,p_x_dim) :: g_box
  real(p_double), dimension(p_x_dim) :: dx
  integer, dimension(p_x_dim) :: nx
  type( t_restart_handle ) :: restart_handle
  type( t_options ) :: sim_options
  type(t_vdf) :: charge
  type(t_vdf_report) :: charge_report

  character(len=*), dimension(2), parameter :: file_format = ['zdf ','hdf5']
  character(len=*), dimension(4), parameter :: parallel_io = ['mpi     ', 'indep   ', 'mpiio(i)', 'mpiio(c)']
  character(len=*), dimension(3), parameter :: merge_type  = ['none  ','p2p   ','gather']

  real(p_double) :: delta, bw, bw_s1, bw_s2, mindelta, maxdelta, t_s1, t_s2
  real(p_double) :: s0, t_stdev, bw_stdev
  integer(p_int64) :: t0, t1
  real(p_k_fld) :: val
  integer :: i, ierr
  integer(p_int64) :: datasize

  ! create grid object from test parameters
  gc_num = 1
  g_box( 1, : ) = 0
  g_box( 2, : ) = 1

  grid % x_dim = p_x_dim
  grid % load_balance = .false.
  grid % io % merge_type = test % merge_type
  grid % io % n_merge(1:p_x_dim) = test % n_merge(1:p_x_dim)

  do i = 1, p_x_dim
    grid % g_nx(i) = test % grid_size(i)
    grid % io % n_merge(i) = test % n_merge(i)
  enddo
  call grid % init( no_co, gc_num, g_box, .false., restart_handle )

  ! initialize parallel partition and io object
  call grid % parallel_partition( no_co, 0 )

  ! Create space object
  dx = 1.0
  do i = 1, p_x_dim
    space % x_bnd( 1, i ) = g_box( 1, i )
    space % x_bnd( 2, i ) = g_box( 2, i )
  enddo
  call setup( space, dx, 0, .false., restart_handle )

  ! Initialize diagnostic file format
  sim_options % file_format = test % file_format
  sim_options % parallel_io = test % parallel_io
  call init_dutil( sim_options, 0 )

  ! create vdf and fill it with local node value
  nx  = grid % my_nx( p_size, 1 : p_x_dim )
  val = no_co % local_rank()

  call charge % new( p_x_dim, 1, nx, gc_num, dx, .false. )
  charge = val

  ! set diagnostic metadata
  charge_report%xname  = (/'x1', 'x2', 'x3'/)

  charge_report%xlabel = (/'x_1', 'x_2', 'x_3'/)
  charge_report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

  charge_report%t = 1.0
  charge_report%n = 1

  charge_report%label = 'Charge Density'
  charge_report%name  = 'Charge Density '

  charge_report%units  = 'n.a.'

  charge_report%prec = p_double

  ! full grid diagnostic
  charge_report%path = path

  dataSize = test % grid_size(1)
  do i = 2, p_x_dim
    dataSize = dataSize * test % grid_size(i)
  end do

  if ( charge_report%prec == p_single ) then
    dataSize = dataSize * 4
  else
    dataSize = dataSize * 8
  endif

  if ( no_co % local_rank() == 0 ) then
    print *, ""
    print "('Test id: ',I0)", id
    print *, ""
    select case (p_x_dim)
    case(1)
      print "(A)", " Size       |  Grid | file | io       | merge  |  merge#     "
      if ( test % merge_type /= p_none ) then
        print '(A12,"|",I6," | ",A," | ",A," | ",A, " | ",I2)', &
          strmem(dataSize), test % grid_size(1), file_format( test % file_format + 1 ), parallel_io( test % parallel_io + 1 ), &
          merge_type( test % merge_type + 1 ), test % n_merge(1)
      else
        print '(A12,I6," | ",A," | ",A," | ",A)', &
          strmem(dataSize), test % grid_size(1), file_format( test % file_format + 1 ), parallel_io( test % parallel_io + 1 ), &
          merge_type( test % merge_type + 1 )
      endif
    case(2)
      print "(A)", " Size       | Grid dims  | file | io       | merge  | merge#"
      if ( test % merge_type /= p_none ) then
        print '(A12,"|",I4," x ",I4," | ",A," | ",A," | ",A, " | ",I2," x ",I2)', &
          strmem(dataSize), test % grid_size(1:2), file_format( test % file_format + 1 ), parallel_io( test % parallel_io + 1 ), &
          merge_type( test % merge_type + 1 ), test % n_merge(1:2)
      else
        print '(A12,"|",I4," x ",I4," | ",A," | ",A," | ",A, " | ---")', &
          strmem(dataSize), test % grid_size(1:2), file_format( test % file_format + 1 ), parallel_io( test % parallel_io + 1 ), &
          merge_type( test % merge_type + 1 )
      endif
    case(3)
      print "(A)", " Size       | Grid dims         | file | io       | merge  |  merge#     "
      if ( test % merge_type /= p_none ) then
      print '(A12,"|",I4," x ",I4," x ",I4," | ",A," | ",A," | ",A, " | ",I2," x ",I2," x ",I2)', &
        strmem(dataSize), test % grid_size(1:3), file_format( test % file_format + 1 ), parallel_io( test % parallel_io + 1 ), &
        merge_type( test % merge_type + 1 ), test % n_merge(1:3)
      else
      print '(A12,"|",I4," x ",I4," x ",I4," | ",A," | ",A," | ",A, " |   ---")', &
        strmem(dataSize),test % grid_size(1:3), file_format( test % file_format + 1 ), parallel_io( test % parallel_io + 1 ), &
        merge_type( test % merge_type + 1 )
      endif
    end select
    print *, ""
    print "(A)", "           Time      Bandwidth"
  endif

  ! Initialize bandwidth measurements
  t_s1 = 0
  t_s2 = 0
  bw_s1 = 0
  bw_s2 = 0
  maxdelta = 0
  mindelta = 0
  s0 = 0

  ! synchronize nodes
  call no_co  % barrier()

  do i = 1, repeat

    ! set file name
    charge_report%filename = 'ioperf_test_'//trim(tostring_int(id))//'_'// trim(tostring_int(i))

    ! save file
    t0 = timer_ticks()

    call charge % write( charge_report, 1, space, grid, no_co )

    call no_co  % barrier()
    t1 = timer_ticks()

    delta = timer_interval_seconds(t0,t1)
    bw = (dataSize / delta) / (1024*1024)

    s0 = i

    t_s1 = t_s1 + delta
    t_s2 = t_s2 + delta**2

    bw_s1 = bw_s1 + bw
    bw_s2 = bw_s2 + bw**2

    if ( i == 1 ) then
      maxdelta = delta
      mindelta = delta
    else
      if ( delta > maxdelta ) maxdelta = delta
      if ( delta < mindelta ) mindelta = delta
    endif

    if ( root(no_co) ) then
      print "('  [',I2,'/',I2,'] ',F6.3,' s, ',F8.2, ' MB/s')", i, repeat, delta, bw

      ! remove the file
      if ( .not. keepFiles ) then
        if ( trim(charge_report%path) /= '' ) then
          if ( test % file_format == p_hdf5_format ) then
            call remove( trim(charge_report%path) // '/' // trim(charge_report%filename) // '.h5', ierr )
          else
            call remove( trim(charge_report%path) // '/' // trim(charge_report%filename) // '.zdf', ierr )
          endif
        else
          if ( test % file_format == p_hdf5_format ) then
            call remove( trim(charge_report%filename) // '.h5', ierr )
          else
            call remove( trim(charge_report%filename) // '.zdf', ierr )
          endif
        endif
      endif
    endif

  enddo

  if ( root(no_co) ) then
    if ( repeat > 1 ) then
      t_stdev = sqrt( (s0*t_s2 - t_s1**2) / (s0 * (s0-1)) )
      delta = t_s1 / s0
      bw_stdev = sqrt( (s0*bw_s2 - bw_s1**2) / (s0 * (s0-1)) )
      bw = bw_s1 / s0
    else
      t_stdev = 0
      delta = t_s1
      bw_stdev = 0
      bw = bw_s1
    endif
    print *, ""
    print "('[final] ',F6.3,' ± ',F6.3,' s, ', F8.2,' ± ',F8.2,' MB/s, ∆ = ',F5.1,' %')", &
      delta, t_stdev, bw, bw_stdev, 100 * bw_stdev / bw
    select case (p_x_dim)
    case(1)
      if ( test % merge_type /= p_none ) then
        write (fid,'(A12,"|",I6," | ",A," | ",A," | ",A, " | ",I4," | ",F6.3," ± ",F6.3," s | ", F8.2," ± ",F8.2," MB/s, ∆ = ",F5.1," %")') &
          strmem(dataSize), test % grid_size(1), file_format( test % file_format + 1 ), parallel_io( test % parallel_io + 1 ), &
          merge_type( test % merge_type + 1 ), test % n_merge(1), delta, t_stdev, bw, bw_stdev, 100 * bw_stdev / bw
      else
        write (fid,'(A12,"|",I6," | ",A," | ",A," | ",A, " |  --- | ",F6.3," ± ",F6.3," s | ", F8.2," ± ",F8.2," MB/s, ∆ = ",F5.1," %")') &
          strmem(dataSize),test % grid_size(1), file_format( test % file_format + 1 ), parallel_io( test % parallel_io + 1 ), &
          merge_type( test % merge_type + 1 ), delta, t_stdev, bw, bw_stdev, 100 * bw_stdev / bw
      endif
    case(2)
      if ( test % merge_type /= p_none ) then
        write (fid,'(A12,"|",I4," x ",I4," | ",A," | ",A," | ",A, " | ",I2," x ",I2," | ",F6.3," ± ",F6.3," s | ", F8.2," ± ",F8.2," MB/s, ∆ = ",F5.1," %")') &
          strmem(dataSize), test % grid_size(1:2), file_format( test % file_format + 1 ), parallel_io( test % parallel_io + 1 ), &
          merge_type( test % merge_type + 1 ), test % n_merge(1:2), delta, t_stdev, bw, bw_stdev, 100 * bw_stdev / bw
      else
        write (fid,'(A12,"|",I4," x ",I4," | ",A," | ",A," | ",A, " |   ---   | ",F6.3," ± ",F6.3," s | ", F8.2," ± ",F8.2," MB/s, ∆ = ",F5.1," %")') &
          strmem(dataSize),test % grid_size(1:2), file_format( test % file_format + 1 ), parallel_io( test % parallel_io + 1 ), &
          merge_type( test % merge_type + 1 ), delta, t_stdev, bw, bw_stdev, 100 * bw_stdev / bw
      endif
    case(3)
      if ( test % merge_type /= p_none ) then
        write (fid,'(A12,"|",I4," x ",I4," x ",I4," | ",A," | ",A," | ",A, " | ",I2," x ",I2," x ",I2, " | ",F6.3," ± ",F6.3," s | ", F8.2," ± ",F8.2," MB/s, ∆ = ",F5.1," %")') &
          strmem(dataSize), test % grid_size(1:3), file_format( test % file_format + 1 ), parallel_io( test % parallel_io + 1 ), &
          merge_type( test % merge_type + 1 ), test % n_merge(1:3), delta, t_stdev, bw, bw_stdev, 100 * bw_stdev / bw
      else
        write (fid,'(A12,"|",I4," x ",I4," x ",I4," | ",A," | ",A," | ",A, " |      ---     | ",F6.3," ± ",F6.3," s | ", F8.2," ± ",F8.2," MB/s, ∆ = ",F5.1," %")') &
          strmem(dataSize),test % grid_size(1:3), file_format( test % file_format + 1 ), parallel_io( test % parallel_io + 1 ), &
          merge_type( test % merge_type + 1 ), delta, t_stdev, bw, bw_stdev, 100 * bw_stdev / bw
      endif
    end select

  endif

  call charge % cleanup( )

  call grid % cleanup( )

end subroutine RunTest

!---------------------------------------------------------------------------------------------------
! Read input file
!---------------------------------------------------------------------------------------------------
subroutine ReadInput( no_co, test_set )

  use zdf_parallel
  use m_diagfile
  use m_diagnostic_utilities

  implicit none
  class( t_node_conf), intent(inout) :: no_co
  type(t_test_set), intent(inout) :: test_set

  character(len=*), parameter :: input_fname = "ioperf.params"
  class( t_input_file ) :: input_file

  character(len = 16), dimension(16) :: grid
  integer :: n_grid
  integer, dimension(p_x_dim) :: igrid
  character(len = 16), dimension(16) :: format
  integer :: n_format
  integer :: iformat
  character(len = 16), dimension(16) :: parallel_io
  integer :: n_parallel_io
  integer :: iparallel_io
  character(len = 16), dimension(16) :: n_merge
  integer :: n_n_merge
  integer, dimension(p_x_dim) :: in_merge
  character(len = 16), dimension(16) :: merge_type
  integer :: n_merge_type
  integer :: imerge_type

  logical :: keepFiles

  integer :: repeat, max_tests
  integer :: i5, i4, i3, i2, i1, i, j, k, ierr
  character(len=1024) :: path

  namelist /nl_test_set/ grid, format, parallel_io, merge_type, n_merge, repeat, path, keepFiles

  ! open input file
  call open_input( input_file, trim( input_fname ), mpi_world_comm() )

  SCR_ROOT('Reading parallel node configuration... ')
  call no_co % read_input( input_file, p_x_dim )

  SCR_ROOT('Reading test set parameters... ')

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_test_set", ierr )

  if ( ierr /= 0 ) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        write(0,*) "Error reading test_set parameters"
      else
        write(0,*) "Error: test_set parameters missing"
      endif
      write(0,*) "aborting..."
    endif
    stop
  endif

  grid = "-"
  format = "-"
  parallel_io = "-"
  n_merge = "-"
  merge_type = "-"
  repeat = 3
  path = ''
  keepFiles = .false.

  read (input_file%nml_text, nml = nl_test_set, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading test_set parameters"
      write(0,*) "aborting..."
    endif
    stop
  endif

  SCR_ROOT('Finished reading input file. ')

  test_set % repeat = repeat
  test_set % path = path
  test_set % keepFiles = keepFiles

  ! Create test_set object
  n_grid = 0
  do i = 1, 16
    if ( grid(i) == "-" ) exit
    n_grid = n_grid + 1
  enddo

  if ( n_grid == 0 ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading test_set parameters, no grid dimensions supplied"
      write(0,*) "aborting..."
    endif
    stop
  endif

  n_format = 0
  do i = 1, 16
    if ( format(i) == "-" ) exit
    n_format = n_format + 1
  enddo
  if (n_format == 0) then
    format(1) = 'zdf'
    n_format = 1
  endif

  n_parallel_io = 0
  do i = 1, 16
    if ( parallel_io(i) == "-" ) exit
    n_parallel_io = n_parallel_io + 1
  enddo
  if (n_parallel_io == 0) then
    parallel_io(1) = 'mpi'
    n_parallel_io = 1
  endif

  n_n_merge = 0
  do i = 1, 16
    if ( n_merge(i) == "-" ) exit
    n_n_merge = n_n_merge + 1
  enddo
  if (n_n_merge == 0) then
    select case (p_x_dim)
    case(1)
      n_merge(1) = "1"
    case(2)
      n_merge(1) = "1,1"
    case(3)
      n_merge(1) = "1,1,1"
    end select
    n_n_merge = 1
  endif

  n_merge_type = 0
  do i = 1, 16
    if ( merge_type(i) == "-" ) exit
    n_merge_type = n_merge_type + 1
  enddo
  if (n_merge_type == 0) then
    merge_type(1) = 'point2point'
    n_merge_type = 1
  endif

  max_tests = n_grid * n_format * n_parallel_io * n_n_merge * n_merge_type

  allocate( test_set%tests( max_tests ) )

  k = 1
  do i5 = 1, n_grid
    ! get grid size
    read( grid(i5), * ) igrid(1:p_x_dim)

    do i4 = 1, n_format
      ! get file format
      select case (lowercase(format(i4)))
      case('zdf')
        iformat = p_zdf_format
      case('hdf5')
        iformat = p_hdf5_format
      case default
        if ( mpi_node() == 0 ) then
          write(0,*) "(*error*) Invalid format value, must be one of 'zdf' or 'hdf5'."
        endif
        stop
        iformat = -1
      end select

      do i3 = 1, n_parallel_io
        ! get parallel io
        select case (lowercase(parallel_io(i3)))
        case('mpi')
          iparallel_io = DIAG_MPI
        case('posix')
          if ( iformat ==  p_hdf5_format ) then
            if ( mpi_node() == 0 ) then
              write(0,*) "Error in global simulation parameters"
              write(0,*) "parallel_io type 'posix' not supported by 'hdf5' file format, aborting..."
            endif
            stop
          endif

          iparallel_io = DIAG_POSIX
        case('independent')
          iparallel_io = DIAG_MPIIO_INDEPENDENT
        case('collective')
          iparallel_io = DIAG_MPIIO_COLLECTIVE
        case default
          if ( mpi_node() == 0 ) then
            write(0,*) "(*error*) Invalid parallel_io value, must be one of 'mpi', 'posix', 'independent' or 'collective'."
          endif
          stop
          iparallel_io = DIAG_MPI
        end select

        do i2 = 1, n_merge_type

          ! Get merge_type
          select case (lowercase(merge_type(i2)))
          case('none')
            imerge_type = p_none
          case('point2point')
            imerge_type = p_point2point
          case('gather')
            imerge_type = p_gather
          case default
            if ( mpi_node() == 0 ) then
              write(0,*) "(*error*) Invalid merge_type value, must be one of 'none', 'point2point' or 'gather'."
            endif
            stop
            imerge_type = p_none
          end select
          
          if ( imerge_type /= p_none ) then
            do i1 = 1, n_n_merge
              ! get n_merge
              read( n_merge(i1), * ) in_merge(1:p_x_dim)

              ! Validate n_merge
              j = 1
              do i = 1, p_x_dim
                if (mod( no_co % nx(i), in_merge(i) ) /= 0 ) then
                  if ( mpi_node() == 0 ) then
                    write(0,*) '(*error*) Invalid n_merge value. n_merge must divide partition sizes evenly.'
                  endif
                  stop
                endif
                j = j * in_merge(i)
              enddo
              if ( j == 1 ) then
                if ( mpi_node() == 0 ) then
                  write(0,*) '(*error*) Invalid n_merge value. n_merge must be >1 in at least one direction.'
                endif
                stop
              endif

              ! Store test
              test_set % tests(k) % grid_size = igrid
              test_set % tests(k) % file_format = iformat
              test_set % tests(k) % parallel_io = iparallel_io
              test_set % tests(k) % merge_type = imerge_type
              test_set % tests(k) % n_merge = in_merge
  
              k = k+1
  
            enddo

          else

            ! Store test
            test_set % tests(k) % grid_size = igrid
            test_set % tests(k) % file_format = iformat
            test_set % tests(k) % parallel_io = iparallel_io
            test_set % tests(k) % merge_type = p_none
            test_set % tests(k) % n_merge = 1

            k = k+1
          endif
          
        enddo
      enddo
    enddo
  enddo

  test_set % n_tests = k - 1

  ! Print info regarding selected tests
  if ( mpi_node() == 0 ) then
    print "(A)",     "Number of parameter combinations selected"
    print "(A,I0)",  "Grid              : ", n_grid
    print "(A,I0)",  "File format       : ", n_format
    print "(A,I0)",  "Parallel I/O      : ", n_parallel_io
    print "(A,I0)",  "Merge nodes       : ", n_n_merge
    print "(A,I0)",  "Merge type        : ", n_merge_type
    print "(A,I0)",  "            Total : ", test_set % n_tests
  endif


end subroutine ReadInput
!---------------------------------------------------------------------------------------------------

subroutine set_diag_sim_info( )

  use m_diagnostic_utilities

  implicit none

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

  sim_info % compile_time = date // " " // time
  sim_info % input_file   = "ioperf.params"

  sim_info % dt = 0

  sim_info % ndims = p_x_dim

  do i = 1, p_x_dim
    sim_info % box_nx(i)    = 0
    sim_info % box_min(i)   = 0
    sim_info % box_max(i)   = 1
    sim_info % periodic(i)  = .false.
    sim_info % move_c(i)    = .false.
    sim_info % node_conf(i) = 1
  enddo

  call alloc( sim_info%nx_x1_node, [1] )
  sim_info%nx_x1_node = 0

  if (sim_info % ndims >= 2) then
    call alloc(sim_info%nx_x2_node, [1])
    sim_info%nx_x2_node = 0

    if (sim_info % ndims == 3) then
      call alloc(sim_info%nx_x3_node, [1])
      sim_info%nx_x3_node = 0
    endif
  endif

end subroutine set_diag_sim_info

end program
