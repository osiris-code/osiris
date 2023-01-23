!-------------------------------------------------------------------------------
! Grid definition module
!
! This file contains the class definition for the following classes:
!
!  t_grid
!  t_msg
!  t_msg_patt
!-------------------------------------------------------------------------------

#include "os-preprocess.fpp"
#include "os-config.h"

module m_grid_define

#include "memory/memory.h"

  use m_system
  use m_parameters
  use m_fparser
  use m_input_file
  use m_node_conf
  use zdf_parallel
  use m_diagfile

  implicit none

  private

  integer, parameter :: p_no_load_balance       = 0  ! evenly distribute nodes
  integer, parameter :: p_static_load_balance   = 1  ! static load balance
  integer, parameter :: p_dynamic_load_balance  = 2  ! dynamic load balance
  integer, parameter :: p_test_load_balance     = 3  ! test non-even load balance
  integer, parameter :: p_expr_load_balance     = 4  ! specified expression load balance

  ! parameters for particle load type
  integer, parameter :: p_density = 0
  integer, parameter :: p_particles = 1

  integer, parameter :: p_global = 1
  integer, parameter :: p_node   = 2
  integer, parameter :: p_grid   = 3

  ! string to id restart data
  character(len=*), parameter :: p_grid_rst_id = "grid rst data - 0x0009"


  ! timing variables
  integer, save :: dynlb_sum_up_load_ev = -1

  ! *******************************************************************************
  ! * Grid I/O object
  ! *******************************************************************************

  ! Type of merging to use
  integer, parameter, public :: p_none        = 0
  integer, parameter, public :: p_point2point = 1
  integer, parameter, public :: p_gather      = 2

  type :: t_grid_io

    ! Chunk size to use for the file
    integer, dimension( p_max_dim ) :: max_chunk_size

    ! -- Grid merging data

    ! Type of merging to use
    integer :: merge_type

    ! number of nodes to merge for data output
    integer, dimension( p_max_dim ) :: n_merge

    ! Communicator involving nodes in the same group
    integer :: group_comm = MPI_COMM_NULL
    ! Size of group
    integer :: group_size
    ! rank in group communicator
    integer :: group_rank

    ! information on all node tiles in group grid ( count|start, dir, rank in group_comm )
    integer, dimension(:,:,:), pointer :: tiles => null()

    ! information on local group tile relative to global grid (count|start, dir)
    integer, dimension( 2, p_max_dim ) :: group_tile

    ! Communicator involving group leaders
    integer :: comm = MPI_COMM_NULL

    ! Chunk size to use for some file formats
    integer, dimension(p_max_dim) :: chunk_size
  contains

    procedure :: init    => init_grid_io
    procedure :: cleanup => cleanup_grid_io

  end type

  ! *******************************************************************************
  ! * Grid object
  ! *******************************************************************************

  type :: t_grid

    integer :: x_dim = 0                      ! number of spatial dimensions

    integer :: ndump_global_load = 0          ! frequency to dump global load info (min, max, avg parts/node)
    integer :: ndump_node_load   = 0          ! frequency to dump node load info (parts/node)
    integer :: ndump_grid_load   = 0          ! frequency to dump grid load info (parts/grid)

    integer, dimension(p_max_dim) :: g_nx     ! global number of cells

    integer, dimension(p_max_dim) :: min_nx   ! minimum number of cells in each partition
                                              ! (used by load balancing algorithm)

    integer, dimension(p_max_dim) :: max_nx   ! maximum number of cells in each partition
                                              ! (used for diagnostics)

    integer :: coordinates ! variable to select type of coordinates for the space

    ! number of nodes in each direction
    integer, dimension(p_max_dim) :: nnodes

    ! number of cells for each partition for all directions
    ! partition( bound, part_idx )
    ! - bound is 1:3 (lower/upper/#cells)
    ! - part_idx is organized as [ 1 .. np1, 1 .. np2, 1 .. np3 ] i.e.
    integer, pointer, dimension(:,:) :: nx_node => null()

    ! global position and number of cells for local node
    integer, dimension(3, p_max_dim) :: my_nx



    integer :: lb_type = p_no_load_balance    ! load balance type

    logical :: lb_gather_max = .false.        ! use maximum value when gathering load
                                              ! from multiple nodes

    ! load balance along given direction
    logical, dimension(p_max_dim) :: load_balance = .false.

    integer :: n_dynamic = 0                  ! number of timesteps between
                                              ! each dynamic load balance

    integer :: start_load_balance = -1        ! number of iterations at which to start load balance

    logical :: balance_on_start = .false.     ! when n == start_load_balance, do a dynamic load balance

    real(p_single) :: max_imbalance = 0.0     ! threshold to trigger reshaping
                                              ! simulation.

    real(p_single) :: cell_weight = 0.0       ! weight for grid cells in computacionial load

    ! array holding the integral of the load
    real(p_double), pointer, dimension(:) :: int_load => null()

    ! simulation box dimensions
    real(p_double), dimension(2, p_max_dim) :: g_box
    
    ! This is parser which will hold the equation specified in the input file
    type(t_fparser) :: spatial_loaddensity

    ! I/O parameters
    type(t_grid_io) :: io

  contains

    procedure :: init          => init_grid
    procedure :: copy          => copy_grid
    procedure :: cleanup       => cleanup_grid
    procedure :: nx_p_part     => nx_p_part_grid
    procedure :: read_input    => read_input_grid
    procedure :: write_checkpoint => write_checkpoint_grid
    procedure :: read_checkpoint => read_checkpoint_grid

    procedure :: num_partitions
    procedure :: get_part_width

    procedure :: parallel_partition => parallel_partition_grid
    procedure :: needs_int_load => needs_int_load_grid

    procedure :: gather_load        => gather_load_grid
    procedure :: if_report          => if_report_grid
    procedure :: ndump              => ndump_report_grid
    procedure :: test_partition     => test_partition_grid
    procedure :: get_node_limits    => get_node_limits_grid

  end type t_grid

  interface
    subroutine parallel_partition_grid( this, no_co, n )
      import t_grid, t_node_conf
      class( t_grid ), intent(inout) :: this
      class( t_node_conf ), intent(in)  :: no_co
      integer, intent(in) :: n
    end subroutine
  end interface

  interface
    subroutine gather_load_grid( this, no_co )
      import t_grid, t_node_conf
      class( t_grid ), intent(inout) :: this
      class( t_node_conf ), intent(in)  :: no_co
    end subroutine
  end interface

  ! aux. types for dynamic load balancing
  type :: t_msg

    ! node to send/receive messages from
    integer :: node

    ! cells to send/receive
    ! cells( lb:ub, x_dim )
    integer, dimension(2,p_max_dim) :: cells

    ! extra parameters used for tiles module
    ! message request id
    integer :: request = MPI_REQUEST_NULL

    ! list of tiles to send/receive in terms of their g_til_aid's
    integer, dimension(:), pointer :: tils => null()

    ! list of integers containing MPI devired types used in tile diag communication
    integer, dimension(:), pointer :: sub_type => null()

    ! buffer for messages
    integer(p_byte), dimension(:), pointer :: buffer => null()

    ! pointer to field data
    ! TODO: make this local variable in send_til_msg_diag() freed immediately after
    ! blocking send
    real(p_k_fld), dimension(:), pointer :: vdf_diag_buffer => null()

  end type t_msg

  type :: t_msg_patt

    ! number of nodes to send messages to
    integer :: n_send = 0

    ! description of messages to send
    type(t_msg), dimension(:), pointer :: send => null()

    ! number of nodes to receive messages from
    integer :: n_recv = 0

    ! description of messages to receive
    type(t_msg), dimension(:), pointer :: recv => null()

  contains

    procedure :: cleanup => cleanup_msg_patt

  end type t_msg_patt

  ! exported symbols
  public :: t_grid, t_grid_io, t_msg, t_msg_patt
  public :: p_no_load_balance, p_static_load_balance, p_dynamic_load_balance, &
            p_test_load_balance, p_expr_load_balance
  public :: p_density, p_particles
  public :: p_global, p_node, p_grid

  public :: p_grid_rst_id, dynlb_sum_up_load_ev

  interface alloc
    module procedure alloc_msg
    module procedure alloc_1d_msg
    module procedure alloc_bound_1d_msg
    end interface

    interface freemem
    module procedure freemem_msg
    module procedure freemem_1d_msg
    end interface

    public :: alloc, freemem

contains

!-----------------------------------------------------------------------------------------
! clear message pattern object
!-----------------------------------------------------------------------------------------
subroutine cleanup_msg_patt( patt )

  ! use m_grid_memory

  implicit none

  class( t_msg_patt ), intent(inout) :: patt

  integer :: i

  if ( patt%n_send>0 ) then
    do i = 1, patt%n_send
      call freemem( patt%send(i)%tils )
      call freemem( patt%send(i)%sub_type )
      call freemem( patt%send(i)%buffer )
    enddo

    call freemem( patt%send )
  endif

  if ( patt%n_recv>0 ) then
    do i = 1, patt%n_recv
      call freemem( patt%recv(i)%tils )
      call freemem( patt%recv(i)%sub_type )
      call freemem( patt%recv(i)%buffer )
    enddo

    call freemem( patt%recv )
  endif

end subroutine cleanup_msg_patt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get partition sizes for given direction
!-----------------------------------------------------------------------------------------
subroutine get_part_width( this, part_width, dir )

  implicit none

  class( t_grid ), intent(in) :: this
  integer, dimension(:), intent(inout) :: part_width
  integer, intent(in) :: dir

  integer :: i, idx

  idx = 1
  do i = 2, dir
    idx = idx + this%nnodes(i-1)
  enddo

  ! part_width = this%nx_node( 3, idx : idx + this%nnodes(dir) - 1 )

  do i = 1, this%nnodes(dir)
    part_width(i) = this%nx_node( 3, idx + i - 1 )
  enddo


end subroutine get_part_width
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!  gets the number of partitions in each direction
!-----------------------------------------------------------------------------------------
function num_partitions( this )

  implicit none

  class( t_grid ), intent(in)     :: this
  integer, dimension(this%x_dim) :: num_partitions

  num_partitions = this%nnodes(1:this%x_dim)

end function num_partitions
!-----------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Gets position of lower/upper boundary of a given partition / direction on the global simulation
! grid. For convenience it also returns the size of the specified partition.
!---------------------------------------------------------------------------------------------------
function nx_p_part_grid( this, bound, dir, part_idx  )

  implicit none

  class( t_grid ), intent(in) :: this
  integer, intent(in) :: bound        ! 1 - lower boundary, 2 - upper boundary, 3 - size
  integer, intent(in) :: dir          ! direction
  integer, intent(in) :: part_idx     ! partition index, 1 to nnodes(dir)

  integer :: nx_p_part_grid
  integer :: shift

  select case (dir)
    case(2)
      shift = this % nnodes(1)
    case(3)
      shift = this % nnodes(1) + this % nnodes(2)
    case default
      shift = 0
  end select

  nx_p_part_grid = this % nx_node( bound, shift + part_idx )

end function nx_p_part_grid
!---------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine init_grid( this, no_co, gc_num, g_box, restart, restart_handle )

  use m_restart
  use m_logprof

  implicit none

  ! dummy variables
  class( t_grid ), intent(inout) :: this
  class( t_node_conf ), intent(in)  :: no_co    ! node configuration

  integer, dimension(:, :), intent(in) :: gc_num
  real(p_double), dimension(:,:), intent(in) :: g_box

  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle


  ! local variables
  integer :: i, x_dim_, nx_node_size

  if ( .not. restart ) then

    x_dim_ = this%x_dim
    this%nnodes(1:x_dim_) = nx( no_co )

    nx_node_size = 0
    do i = 1, x_dim_
      this%g_box( :, i ) = g_box( :, i )

      ! Increase nx_node buffer size
      nx_node_size   = nx_node_size + this%nnodes(i)

      ! Minimum number of cells per node
      this%min_nx(i) = max( 1, gc_num(p_lower,i)+gc_num(p_upper,i))
    enddo

    call alloc( this%nx_node, (/ 3, nx_node_size /) )

  else

    call this % read_checkpoint( restart_handle )

  endif

  ! create event for timing
  if (dynlb_sum_up_load_ev < 0) dynlb_sum_up_load_ev = create_event('dynlb sum up load')

end subroutine init_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine read_input_grid( this, input_file, partition )

  implicit none

  class( t_grid ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in), dimension(:) :: partition

  ! local variables

  integer, dimension(p_x_dim) :: nx_p
  character(len=16)           :: coordinates

  logical, dimension(p_x_dim) :: load_balance

  logical :: any_lb

  character(len=16) :: lb_type, lb_gather
  integer :: n_dynamic, start_load_balance
  logical :: balance_on_start
  real( p_single ) :: max_imbalance, cell_weight

  integer :: ndump_global_load, ndump_node_load, ndump_grid_load
  character (len=p_max_expr_len) :: spatial_loaddensity

  integer, dimension(p_x_dim) :: io_nmerge
  character(len=16) :: io_merge_type


  namelist /nl_grid/ nx_p, coordinates, io_nmerge, io_merge_type, &
                     load_balance, lb_type, lb_gather, n_dynamic, start_load_balance, &
                     balance_on_start, max_imbalance, cell_weight, &
                     ndump_global_load, ndump_node_load, ndump_grid_load, spatial_loaddensity


  integer :: i, ierr

  ! executable statements

  this%x_dim = p_x_dim


  nx_p = 0
  coordinates = "cartesian"

  io_merge_type = "point2point"
  io_nmerge = 1

  load_balance = .false.
  lb_type   = "none"
  lb_gather = "sum"
  n_dynamic = -1
  spatial_loaddensity = "NO FUNCTION DEFINED"

  ! Default is to start load balance at timestep 0 (initialization)
  start_load_balance = -1
  balance_on_start = .false.

  max_imbalance = 0.0
  cell_weight   = 0.0

  ndump_global_load = 0
  ndump_node_load   = 0
  ndump_grid_load   = 0

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_grid", ierr )

  if ( ierr /= 0 ) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        write(0,*) "Error reading grid parameters"
      else
        write(0,*) "Error - grid parameters missing"
      endif
      write(0,*) "aborting..."
    endif
    stop
  endif

  read (input_file%nml_text, nml = nl_grid, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading grid parameters"
      write(0,*) "aborting..."
    endif
    stop
  endif

  ! validate global grid parameters
  do i = 1, this%x_dim
    if (nx_p(i) <= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error in grid parameters"
        write(0,*) "Invalid grid parameters, nx_p must be >= 0 for all dimensions."
        write(0,*) "nx_p(",p_x_dim,") = ", nx_p
        write(0,*) "(also check the number of dimensions of the binary)"
        write(0,*) "aborting..."
      endif
      stop
    endif
  enddo

  ! store global grid size
  this%g_nx(1:this%x_dim) = nx_p(1:this%x_dim)

  select case ( trim( coordinates ) )
  case ( "cylindrical" )
    SCR_ROOT('cylindrical coordinates - B1 on axis')
    this%coordinates = p_cylindrical_b
    if ( p_x_dim /= 2 ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error in grid parameters"
        write(0,*) "Invalid coordinates, cylindrical coordinates "
        write(0,*) "are only available in 2D."
      endif
      stop
    endif

  case ( "cartesian" )
    SCR_ROOT('cartesian coordinates')
    this%coordinates = p_cartesian

  case default
    if ( mpi_node() == 0 ) then
      write(0,*) "Error in grid parameters"
      write(0,*) "Invalid coordinates, coordinates must be:"
      write(0,*) "'cartesian' (1D, 2D, 3D) or 'cylindrical' (2D)"
    endif
    stop

  end select

  this%start_load_balance = start_load_balance
  this%balance_on_start = balance_on_start

  select case (trim(lb_type))
  case("none")
    this%lb_type = p_no_load_balance
    this%load_balance = .false.

  case("test")
    this%lb_type = p_test_load_balance
    this%load_balance = .true.
    this%n_dynamic = n_dynamic

  case("static")
    this%lb_type = p_static_load_balance
    this%cell_weight = cell_weight

    ! The start load balance parameter is ignored, since this occurs a timestep 0
    this%start_load_balance = -1

  case("dynamic")
    this%lb_type = p_dynamic_load_balance
    if (n_dynamic < 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error in grid parameters"
        write(0,*) "n_dynamic must be >= 1 when using dynamic load balancing"
        write(0,*) "aborting..."
      endif
      stop
    endif
    this%n_dynamic = n_dynamic

    this%max_imbalance = max_imbalance
    this%cell_weight = cell_weight

    select case (trim( lb_gather ))
    case ("sum")
      this%lb_gather_max = .false.
    case ("max")
      this%lb_gather_max = .true.
    case default
      if ( mpi_node() == 0 ) then
        write(0,*) "Error in grid parameters"
        write(0,*) "lb_gather must be 'sum' or 'max' when using dynamic load balancing"
        write(0,*) "aborting..."
      endif
      stop
    end select

  case("expression")
    this%lb_type = p_expr_load_balance
    this%cell_weight = cell_weight

    select case (p_x_dim)
      case (1)
        call setup(this%spatial_loaddensity, trim(spatial_loaddensity), (/'x1'/), ierr)
      case (2)
        call setup(this%spatial_loaddensity, trim(spatial_loaddensity), (/'x1','x2'/), ierr)
      case (3)
        call setup(this%spatial_loaddensity, trim(spatial_loaddensity), (/'x1','x2','x3'/), ierr)
    end select

    ! check if function compiled ok
    if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error compiling spatial_loaddensity"
        write(0,*) "aborting..."
      endif
      stop
    endif

  case default
    if ( mpi_node() == 0 ) then
      write(0,*) "Error in grid parameters"
      write(0,*) "Invalid lb_type selected."
      write(0,*) "aborting..."
    endif
    stop
  end select

  if ( this%lb_type /= p_no_load_balance ) then

    ! Store load balance directions
    any_lb = .false.
    do i = 1, this%x_dim
      this%load_balance(i) = load_balance(i)
      if ( load_balance(i) ) any_lb = .true.
    enddo

    if ( .not. any_lb ) then
      if ( mpi_node() == 0 ) then
        write(0,*) ""
        write(0,*) "Error in grid parameters"
        write(0,*) "Load balance type '", trim(lb_type), "' selected, but no load balance"
        write(0,*) "direction specified"
        write(0,*) "aborting..."
      endif
      stop
    endif

    ! store cell weight
    this%cell_weight = cell_weight
    if ( cell_weight < 0.0 ) then
      if ( mpi_node() == 0 ) then
        write(0,*) ""
        write(0,*) "Error in grid parameters"
        write(0,*) "cell_weight must be >= 0."
        write(0,*) "aborting..."
      endif
      stop
    endif
  endif


  ! grid reports
  this%ndump_global_load = ndump_global_load
  this%ndump_node_load   = ndump_node_load
  this%ndump_grid_load   = ndump_grid_load


  ! Data merging for file output
  select case (trim(io_merge_type))
  case("none")
    this % io % merge_type = p_none
  case("point2point")
    this % io % merge_type = p_point2point
  case("gather")
    this % io % merge_type = p_gather
  case default
    if ( mpi_node() == 0 ) then
      write(0,*) "Error in grid parameters"
      write(0,*) "io_merge_type must be one of 'none', 'point2point', or 'gather'"
      write(0,*) "aborting..."
    endif
    stop
  end select

  do i = 1, this%x_dim

    if ( io_nmerge(i) < 1 ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error in grid parameters"
        write(0,*) "io_nmerge(:) must be >= 1 in all directions."
        write(0,*) "aborting..."
      endif
      stop
    endif

    if ( mod(partition(i), io_nmerge(i) ) /= 0 ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error in grid parameters"
        write(0,*) "io_nmerge(:) must divide number of nodes in parallel partition exactly."
        write(0,*) "aborting..."
      endif
      stop
    endif

    this % io % n_merge(i) = io_nmerge(i)
  enddo

  if ( this % io % merge_type == p_none .and. product(io_nmerge(1:this%x_dim)) > 1 ) then
    if (mpi_node() == 0) then
      write(0,'(A,I0,A)') " (*warning*) io_nmerge(*) is > 1 but io_merge_type was set to 'none'"
      write(0,'(A)')      " (*warning*) no data merging for I/O will be performed"
    endif
  endif


end subroutine read_input_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Checks if creating the parallel partition requires a computational load estimate
!-----------------------------------------------------------------------------------------
function needs_int_load_grid( this, n )

  implicit none

  logical :: needs_int_load_grid
  class( t_grid ), intent(in)     :: this
  integer, intent(in) :: n

  needs_int_load_grid = ((this%lb_type == p_static_load_balance) .or. &
                    (this%lb_type == p_dynamic_load_balance) .or. &
                    (this%lb_type == p_expr_load_balance)) .and. &
                   ( n >= this%start_load_balance )

end function needs_int_load_grid


!-----------------------------------------------------------------------------------------
! Write checkpoint information for this object
!-----------------------------------------------------------------------------------------
subroutine write_checkpoint_grid( this, restart_handle )

  use m_restart

  implicit none

  class( t_grid ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for grid object.'
  integer :: ierr

  ! executable statements

  restart_io_wr( p_grid_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%x_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%g_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%coordinates, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%nnodes, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%g_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%min_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%max_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%nx_node, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%my_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )


  ! - Dynamic load balance parameters
  restart_io_wr( this%lb_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%lb_gather_max, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%load_balance, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%n_dynamic, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%start_load_balance, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%balance_on_start, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%max_imbalance, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%cell_weight, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! - Global box (phsysical) dimensions
  restart_io_wr( this%g_box, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

end subroutine write_checkpoint_grid
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Read checkpoint information for grid object
!-----------------------------------------------------------------------------------------
subroutine read_checkpoint_grid( this, restart_handle )

  use m_restart

  implicit none

  class( t_grid ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=len(p_grid_rst_id)) :: rst_id
  character(len=*), parameter :: err_msg = 'error reading restart data for grid object.'
  integer :: i, ierr, nx_node_size


  call this % cleanup()

  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  ! check if restart file is compatible
  if ( rst_id /= p_grid_rst_id) then
  write(0,*) '(*error*) Corrupted restart file, or restart file from incompatible binary (grid)'
  call abort_program(p_err_rstrd)
  endif

  restart_io_rd( this%x_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%g_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%coordinates, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%nnodes, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%g_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%min_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%max_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  nx_node_size = 0
  do i = 1, this%x_dim
    nx_node_size = nx_node_size + this%nnodes(i)
  enddo

  call alloc( this%nx_node, (/ 3, nx_node_size /))
  restart_io_rd( this%nx_node, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%my_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )


  ! - Dynamic load balance parameters
  restart_io_rd( this%lb_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%lb_gather_max, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%load_balance, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%n_dynamic, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%start_load_balance, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%balance_on_start, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%max_imbalance, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_rd( this%cell_weight, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! - Global box (phsysical) dimensions
  restart_io_rd( this%g_box, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )


end subroutine read_checkpoint_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! clear object data
!-----------------------------------------------------------------------------------------
subroutine cleanup_grid(this)

  implicit none

  class( t_grid ), intent(inout) :: this

  call freemem( this%nx_node ); this % nx_node => null()
  call freemem( this%int_load ); this % int_load => null()

  call this % io % cleanup()

end subroutine cleanup_grid
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine cleanup_grid_io( this )

  implicit none

  class( t_grid_io ),   intent(inout) :: this

  integer :: ierr


  call freemem( this % tiles )
  this % tiles => null()

  if ( this % merge_type /= p_none ) then
    if ( this % comm /= MPI_COMM_NULL ) then
      call mpi_comm_free( this % comm, ierr )
      this % comm = MPI_COMM_NULL
    endif
  endif

end subroutine cleanup_grid_io
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine init_grid_io( this, grid, no_co )

  use m_node_conf

  implicit none

  class( t_grid_io ),   intent(inout) :: this
  class( t_grid ),      intent(in)    :: grid
  class( t_node_conf ), intent(in)    :: no_co

  ! local node rank on global MPI partition
  integer :: global_rank

  ! local node coordinates on global MPI partition
  integer, dimension( grid%x_dim ) :: global_coords

  ! local group leader coordinates on global MPI partition
  integer, dimension( grid%x_dim ) :: group_leader_coords

  ! Group leader rank on global MPI partition
  integer :: group_leader_rank

  ! Color and key for creating new MPI communicators
  integer :: color, key

  ! count | start  of local chunk on group
  integer, dimension(2, grid%x_dim) :: tile

  ! Variables for calculating the chunk_size
  integer :: chunk_size, max_chunk_size

  integer ::  i, j, k, ierr

  ! Disable merging if n_merge(:) == 1
  if ( product( this%n_merge(1:grid%x_dim)) == 1 ) then
    this % merge_type = p_none
  endif

  ! If merge_type is none make n_merge(:) = 1
  if ( this % merge_type == p_none ) this%n_merge = 1

  ! Get rank and coordinates of local node on global MPI partition
  ! These values are 0 based
  global_rank   = no_co % local_rank()
  global_coords = no_co % local_coords()

  do i = 1, grid % x_dim
    ! Find global coordinates of group leader
    group_leader_coords(i) = (global_coords(i)/this % n_merge(i)) * this % n_merge(i)

    ! Find information on group tile relative to global grid (count|start, dir)
    this%group_tile( 1, i ) = grid%nx_p_part( 2,       i, group_leader_coords(i) + this % n_merge(i) ) - &
                              grid%nx_p_part( 1,       i, group_leader_coords(i) + 1  ) + 1
    this%group_tile( 2, i ) = grid%nx_p_part( p_lower, i, group_leader_coords(i) + 1  ) - 1
  enddo

  ! Find rank of local merged partition root node
  group_leader_rank = no_co % coords_rank( group_leader_coords )

  if ( this % merge_type /= p_none .and. no_co % comm_size() > 1) then

    ! Create communicator with members of the same group
    color = group_leader_rank

    ! Make sure group leader will be the root of new communicator
    if ( global_rank == group_leader_rank ) then
      key = 0
    else
      key = 1 + global_rank
    endif
    call mpi_comm_split( no_co % comm, color, key, this % group_comm, ierr )

    ! Get group size (must match n_merge(1) * n_merge(2) )
    call mpi_comm_size( this % group_comm, this % group_size, ierr )

    ! Sanity check
    j = this%n_merge(1)
    do i = 2, grid % x_dim
      j = j * this%n_merge(i)
    enddo
    if ( this % group_size /= j ) then
      SCR_MPINODE('(*error*) Invalid group size ')
      call abort_program( p_err_invalid )
    endif

    ! Get tile rank on group
    call mpi_comm_rank( this % group_comm, this % group_rank, ierr )

    ! Allocate on group root only
    if ( this % group_rank == 0 ) then
      ! Sanity check
      if ( global_rank /= group_leader_rank ) then
        SCR_MPINODE('(*error*) Group rank 0 is not group leader ')
        call abort_program( p_err_invalid )
      endif

      call alloc( this % tiles, [ 2, grid % x_dim, this % group_size ])
    else
      this%tiles => null()
    endif

    ! Gather tile limits on group root
    do i = 1, grid % x_dim
      tile(1, i) = grid % my_nx( 3, i )
      tile(2, i) = (grid % my_nx( p_lower, i )-1) - this % group_tile( 2, i )
    enddo

    call mpi_gather( tile,       2 * grid % x_dim, MPI_INTEGER, &
                     this%tiles, 2 * grid % x_dim, MPI_INTEGER, &
                     0, this%group_comm, ierr )

    ! create an MPI group containing only the group leaders
    if ( this%group_rank == 0 ) then
      color = 0
    else
      color = 1
    endif

    ! Make sure global root will be the root of new communicator
    if ( global_rank == 0 ) then
      key = 0
    else
      key = 1 + global_rank
    endif

    ! Create communicator
    call mpi_comm_split( no_co % comm, color, key, this % comm, ierr )

    ! Find chunk size for file formats that support chunking
    do i = 1, grid % x_dim
      max_chunk_size = 0

      do j = 1, grid % nnodes(i), this % n_merge(i)
        chunk_size = 0
        do k = 0, this % n_merge(i)-1
          chunk_size = chunk_size + grid % nx_p_part( 3, i, j + k )
        enddo
        if ( chunk_size > max_chunk_size ) max_chunk_size = chunk_size
      enddo
      this % chunk_size(i) = max_chunk_size
    enddo

  else
    ! No merging being performed
    this % merge_type = p_none

    ! Group holds only one node
    this % group_comm = MPI_COMM_SELF
    this % group_size = 1
    this % group_rank = 0

    ! No need to store tile information
    this%tiles => null()

    ! The group of group leaders is the same as the global group
    this % comm = no_co % comm

    do i = 1, grid % x_dim
      this % chunk_size(i) = grid % max_nx( i )
    enddo
  endif

end subroutine init_grid_io
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get dump frequency of required diagnostic
!-----------------------------------------------------------------------------------------
function ndump_report_grid( this, rep_type )
  
  implicit none
  
  class( t_grid ), intent(in) :: this
  integer, intent(in) :: rep_type
  integer :: ndump_report_grid
  
  select case( rep_type )
    case ( p_global )
      ndump_report_grid = this%ndump_global_load
    case ( p_node )
      ndump_report_grid = this%ndump_node_load
    case ( p_grid )
      ndump_report_grid = this%ndump_grid_load
    case default
      ndump_report_grid = 0
  end select
    
end function ndump_report_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Check if load diagnostics are required at this timestep
!-----------------------------------------------------------------------------------------
function if_report_grid( this, n, rep_type )
  
  implicit none
  
  class( t_grid ), intent(in) :: this
  integer, intent(in) :: n, rep_type
  logical :: if_report_grid
  
  select case( rep_type )
    case ( p_global )
      if_report_grid = test_if_report( n, this%ndump_global_load )
    case ( p_node )
      if_report_grid = test_if_report( n, this%ndump_node_load )
    case ( p_grid )
      if_report_grid = test_if_report( n, this%ndump_grid_load )
    case default
      if_report_grid = .false.
  end select
  
contains

  function test_if_report( n, ndump )
    
    implicit none
    
    integer, intent(in) :: n, ndump
    logical :: test_if_report
    
    if ( ndump > 0 ) then
      test_if_report = ( mod( n, ndump ) == 0 )
    else
      test_if_report = .false.
    endif
  
  end function test_if_report
  
end function if_report_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Test parallel partition
!  - The minimum number of grid cells per node is gc_num( p_lower ) + gc_num( p_upper )
!-----------------------------------------------------------------------------------------
subroutine test_partition_grid( this, no_co, gc_num )

  implicit none

  class( t_grid ), intent(inout) :: this
  class( t_node_conf ), intent(in)  :: no_co
  integer, dimension(:,:), intent(in) :: gc_num

  integer :: i, min_nx

  do i = 1, p_x_dim

    ! Minimum number of cells along given direction
    min_nx = max( 1, gc_num( p_lower, i ) + gc_num( p_upper, i ) )

    if ( min_nx * nx( no_co, i ) > this%g_nx(i) ) then
     write(0,'(A,I0)') '(*error*) Too many partitions along direction ', i
     write(0,       *) '(*error*) nx = ', this%g_nx(i), ' n_nodes = ', nx( no_co, i ), &
             ' mininum grid cells/node = ', min_nx
     stop
    endif

  enddo

end subroutine test_partition_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!  Returns the grid limits for all nodes
!  The nx_p_nodes needs to be already allocated and has the following structure:
!  nx_p_nodes(bound, dim, pe)
!  where bound is 1:2 (lower/upper boundary), dim is 1:p_x_dim (dimension),
!  and pe is 1:no_num(no_co) the node number
!-----------------------------------------------------------------------------------------
subroutine get_node_limits_grid( this, no_co, nx_p_nodes  )

  implicit none

  class( t_grid ), intent(in)            :: this
  class( t_node_conf ), intent(in)       :: no_co
  integer, dimension(:,:,:), intent(out) :: nx_p_nodes

  ! local variables

  integer :: node, i, idx

  do node = 1, no_num(no_co)
    idx = 0
    do i = 1, this%x_dim
      nx_p_nodes(1,i,node) = this%nx_node( 1, idx + no_co%ngp( node, i ) )
      nx_p_nodes(2,i,node) = this%nx_node( 2, idx + no_co%ngp( node, i ) )
      idx = idx + this%nnodes(i)
    enddo
  enddo

end subroutine get_node_limits_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Copy data from another object
!-----------------------------------------------------------------------------------------
subroutine copy_grid(this, from)

    implicit none
  
    class( t_grid ), intent(in) :: from
    class( t_grid ), intent(inout) :: this
  
    ! Cleanup any dynamically allocated data
    call this % cleanup()
  
    this%x_dim             = from%x_dim
  
    this%ndump_global_load = from%ndump_global_load
    this%ndump_node_load   = from%ndump_node_load
    this%ndump_grid_load   = from%ndump_grid_load
  
    this%g_nx      = from%g_nx
    this%min_nx    = from%min_nx
    this%max_nx    = from%max_nx
  
    this%coordinates = from%coordinates
  
    this%nnodes    = from%nnodes
  
    call alloc( this%nx_node, (/ 3, size( from%nx_node, 2 ) /) )
    this%nx_node = from%nx_node
  
    this%my_nx         = from%my_nx
    this%lb_type       = from%lb_type
    this%lb_gather_max = from%lb_gather_max
    this%load_balance  = from%load_balance
  
    this%n_dynamic     = from%n_dynamic
    this%start_load_balance     = from%start_load_balance
    this%balance_on_start = from%balance_on_start
    this%max_imbalance = from%max_imbalance
    this%cell_weight   = from%cell_weight
  
    if ( associated( from%int_load) ) then
      call alloc( this%int_load, [ size(from%int_load) ])
      this%int_load = from%int_load
    endif
  
    this%g_box = from%g_box
  
    ! Not sure if these will copy correctly
    this%spatial_loaddensity = from%spatial_loaddensity
    this%io = from%io
  
  end subroutine copy_grid
  
!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation for t_msg

#define __TYPE__ type( t_msg )
#define __TYPE_STR__ "t_msg"
#define FNAME( a )  a ## _msg
#define __MAX_DIM__ 1
#include "memory/mem-template.h"

!---------------------------------------------------------------------------------------------------


end module m_grid_define
