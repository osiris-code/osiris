!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     node configuration class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#define DEBUG_FILE 1

#include "os-config.h"
#include "os-preprocess.fpp"

module m_node_conf

#include "memory/memory.h"

!use mpi

use m_parameters
use m_system
use m_restart
use m_logprof

implicit none

! restrict access to things explicitly declared public
private

! string to id restart data
character(len=*), parameter :: p_nconf_rst_id = "node conf rst data - 0x0002"


! ping constants
integer, parameter :: ping_size = 4
integer, dimension(ping_size), parameter :: ping_code = (/3,5,7,11/)

! reduce operations
integer, parameter :: p_sum = MPI_SUM
integer, parameter :: p_max = MPI_MAX
integer, parameter :: p_min = MPI_MIN
integer, parameter :: p_prod = MPI_PROD

! topology calculation
integer, parameter :: p_topology_mpi = 0
integer, parameter :: p_topology_old = 1
integer, parameter :: p_topology_bgq = 2

type :: t_node_conf

  ! number of dimensions
  integer :: x_dim

  ! mpi communicator
  integer :: comm = mpi_comm_null
  logical :: free_comm = .true.

  ! total number of nodes
  integer :: no_num

  ! number of nodes in each dimension
  integer, dimension(p_max_dim) :: nx

  ! parallel rank of this node (Fortran ordering)
  integer :: my_aid_

  ! Node position in global partition
  integer, dimension(p_max_dim) :: my_ngp

  ! parallel rank of neighboring nodes
  integer, dimension( 2, p_max_dim ) :: neighbor

  ! flags for periodic node grid (periodic boundary conditions)
  logical, dimension(p_max_dim) :: ifpr_

  ! number of threads per node
  integer :: n_threads

  ! Type of topology to use ( mpi, old, bgq )
  integer :: topology_type

contains

    procedure :: init => init_node_conf  
    procedure :: cleanup => cleanup_node_conf  
    procedure :: read_input => read_input_node_conf
    procedure :: write_checkpoint => write_checkpoint_node_conf
    procedure :: restart_read => restart_read_node_conf

  ! This is depreacated and will be removed in the near future
  ! Please use local_rank instead
  procedure, public :: my_aid
  
  procedure :: on_edge => on_edge_node_conf
  procedure :: physical_boundary_local => physical_boundary_local_nconf
  procedure :: physical_boundary_node => physical_boundary_node_nconf
  generic   :: phys_bound => physical_boundary_local, physical_boundary_node
  procedure :: ngp => ngp_nconf
  procedure :: ngp_id => ngp_id_nconf
  procedure :: ifpr => relevant_ifpr_nconf
  
      ! These values are 0 based
    procedure :: local_rank      ! MPI rank of local node
    procedure :: local_coords    ! Coordinates of local node in parallel partition
    procedure :: rank_coords     ! Coordinates of specified rank in parallel partition
    procedure :: coords_rank     ! MPI rank of the node with the specified coordinates
    procedure :: barrier         ! synchronize nodes

    procedure :: comm_size

end type t_node_conf

interface no_num
  module procedure no_num_node_conf
end interface

interface nx
  module procedure get_nx
  module procedure get_nx_dim
end interface

interface neighbor
  module procedure neighbor
end interface

interface my_ngp
  module procedure my_ngp_all
  module procedure my_ngp_dim
end interface

interface send_ping
  module procedure send_ping
end interface

interface recv_ping
  module procedure recv_ping
end interface

interface reduce_array
  module procedure reduce_array_int_1d

  module procedure reduce_array_rsing_1d
  module procedure reduce_array_rsing_2d
  module procedure reduce_array_rsing_3d

  module procedure reduce_array_rdoub_1d
  module procedure reduce_array_rdoub_2d
  module procedure reduce_array_rdoub_3d
end interface

interface reduce_array_size
  module procedure reduce_array_size_rsing_1d
  module procedure reduce_array_size_rdoub_1d
end interface


interface reduce
  module procedure reduce_int
  module procedure reduce_int64
  module procedure reduce_rsing
  module procedure reduce_rdoub
end interface

interface gather_array
  module procedure gather_array_int_1d
end interface

interface allgather
  module procedure allgather_int
end interface

interface broadcast
  module procedure broadcast_int
end interface

interface periodic
  module procedure periodic_dir
  module procedure periodic_all
end interface

interface comm
  module procedure comm
end interface
interface n_threads
  module procedure n_threads_node_conf
end interface

interface root
  module procedure root
end interface

! declare things that should be public
public :: t_node_conf
public :: no_num, n_threads, nx, neighbor, my_ngp
public :: periodic, root

public :: gather_array, reduce, reduce_array, reduce_array_size
public :: send_ping, recv_ping, comm
public :: allgather, broadcast

! global operation codes
public :: p_sum, p_min, p_max, p_prod
public :: p_lower, p_upper
public :: p_topology_mpi

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine read_input_node_conf( this, input_file, x_dim )

  use m_input_file

  implicit none

  class( t_node_conf ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: x_dim

  integer, dimension(p_x_dim) ::  node_number
  logical, dimension(p_x_dim) ::  if_periodic
  integer :: n_threads

  character(len=16) :: topology

  namelist /nl_node_conf/ node_number, if_periodic, n_threads, topology

  integer :: i, ierr

  this%x_dim = x_dim

  node_number = 1
  if_periodic = .false.
  n_threads = 1

  ! Default topology routines
  topology = "mpi"

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_node_conf", ierr )

  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        write(0,*) "Error reading node_conf parameters"
      else
        write(0,*) "Error: node_conf parameters missing"
      endif
      write(0,*) "aborting..."
    endif
    stop
  endif

  read (input_file%nml_text, nml = nl_node_conf, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) "aborting..."
    endif
    stop
  endif

  this%no_num = 1
  do i=1, x_dim
    this%nx(i)   = node_number(i)
    this%ifpr_(i) = if_periodic(i)
    this%no_num  = this%no_num * this%nx(i)
  enddo

  ! process topology type
  select case ( trim( topology ))
  case ( "mpi" )
    this%topology_type = p_topology_mpi
  case ( "old" )
    this%topology_type = p_topology_old

#ifdef __bgq__
  case ( "bgq" )
    this%topology_type = p_topology_bgq
#endif

  case default
    if ( mpi_node() == 0 ) then
      write(0,*) 'Error reading node_conf parameters'
      write(0,*) 'invalid topology type: "', trim( topology ), '"'
#ifdef __bgq__
      write(0,*) 'Valid values are "bgq", "mpi" and "old".'
#else
      write(0,*) 'Valid values are "mpi" and "old".'
#endif
      write(0,*) 'aborting...'
    endif
    stop
  end select


  if ( this%no_num < 1 ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) 'Number of nodes (node_number) must be > 0 in all directions'
      write(0,*) "aborting..."
    endif
    stop
  endif


  ! Multithread support
  this%n_threads = n_threads

  if ( n_threads < 1 ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) "Invalid number of threads per node (n_threads)."
      write(0,*) "aborting..."
    endif
    stop
  endif

#ifndef _OPENMP
  if ( n_threads /= 1 ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) "Multiple threads per node are only supported when using OpenMP."
      write(0,*) "aborting..."
    endif
    stop
  endif
#endif


end subroutine read_input_node_conf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine cleanup_node_conf( this, barrier )
!-----------------------------------------------------------------------------------------
! Object destructor, cleans up all dynamically
! allocated memory
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_node_conf ), intent( inout )  ::  this
  logical, intent(in), optional :: barrier
  integer :: ierr

  ! Synchronize nodes
  if (present(barrier)) then
    if (barrier) then
      call this % barrier()
    endif
  else
    call this % barrier()
  endif

  ! Free communicator if necessary
  if ((this%free_comm) .and. ( this%comm /= mpi_comm_null )) then
    call mpi_comm_free( this%comm, ierr )
    if (ierr /= 0) then
      ERROR("Unable to free communicator")
      call abort_program( p_err_dealloc )
    endif

  endif

  this%comm = mpi_comm_null

end subroutine cleanup_node_conf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine init_node_conf( this )
!-----------------------------------------------------------------------------------------
! sets up this data structure from the given information
!-----------------------------------------------------------------------------------------

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  class( t_node_conf ), intent( inout )  ::  this

  integer  ::  started, id, i, ierr, n_threads

  ! Initialize MPI if required
  if ( this%no_num > 1 ) then

    call MPI_COMM_SIZE( mpi_comm_world, started, ierr )
    if ( ierr /= MPI_SUCCESS ) then
      ERROR("Unable to get MPI size.")
      call abort_program( p_err_mpi )
    endif

    if ( started /= this%no_num ) then
      write(0,'(A,I0,A)') "(*error*) ", started, " nodes started up"
      write(0,'(A,I0,A)') "(*error*) ", this%no_num, " nodes should have started"
      call abort_program( p_err_invalid )
    endif

    call MPI_COMM_RANK( mpi_comm_world, id, ierr )
    if ( ierr /= MPI_SUCCESS ) then
      ERROR("Unable to get MPI rank.")
      call abort_program( p_err_mpi )
    endif

    ! store using fortran ordering
    this%my_aid_ = id + 1

  else
    this%my_aid_ = 1
  endif


#ifdef _OPENMP

  ! Initilize OpenMP if required
  if ( this%n_threads > 1 ) then

    if ( this%n_threads > omp_get_max_threads() ) then
      write(0,*) '(*error*) The number of threads requested, ', this%n_threads , ' is more than the'
      write(0,*) '(*error*) maximum number of threads available, ', omp_get_max_threads()
      call abort_program()
    endif

  endif

  ! Set the number of threads to use
  call omp_set_num_threads( this%n_threads )

  ! Check the operation was successfull (on some systems you need to set the number of threads
  ! externally e.g. in the job script)

  !$omp parallel shared(n_threads)
  n_threads = omp_get_num_threads()
  !$omp end parallel

  if ( this%n_threads /= n_threads ) then
    write(0,*) '(*error*) Unable to set the number of threads requested'
    write(0,*) '(*error*) You may need to set the number of threads outside OSIRIS e.g. ', &
                          'in the job script'
    call abort_program()
  endif

  ! Disable dynamic variation of number of threads
  call omp_set_dynamic (.false.)

  ! Print parallel information
  if ( this%my_aid_ == 1 ) then
    print '(A,I0,A,I0,A)', "  - Running on ", this%no_num, " processes, ", this%n_threads, &
               " thread(s) per process"
  endif

#else

  ! Print parallel information
  if ( this%my_aid_ == 1 ) then
    print '(A,I0,A)', "  - Running on ", this%no_num, " processes"
  endif

#endif


  ! setup comm, ngp and neighbor data
  SCR_ROOT(" setting up topology")
  if ( this%no_num > 1 ) then

  select case ( this%topology_type )
    case ( p_topology_mpi )
    if ( this%my_aid_ == 1 ) then
      print '(A)', "  - Using MPI_Cart_create based topology"
    endif
      call create_topology_mpi( this )
    case ( p_topology_old )
    if ( this%my_aid_ == 1 ) then
      print '(A)', "  - Using old topology routines (debug only)"
    endif
      call create_topology_old( this )
#ifdef __bgq__
    case ( p_topology_bgq )
    if ( this%my_aid_ == 1 ) then
      print '(A)', "  - Using BG/Q topology routines"
    endif
      call create_topology_bgq( this )
#endif
    end select

  else

    ! Serial run, no need to setup the topology
    this%my_ngp   = 1
    do i = 1, this%x_dim
      if ( this%ifpr_(i) ) then
        this%neighbor(:,i) = 1
      else
        this%neighbor(:,i) = -1
      endif
    enddo

  endif

end subroutine init_node_conf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Create simulation topology using MPI topology routines
!-----------------------------------------------------------------------------------------
subroutine create_topology_mpi( this )
!-----------------------------------------------------------------------------------------
  implicit none

  class( t_node_conf ), intent(inout)  ::  this

  integer :: ndim, dims(p_max_dim), coords(p_max_dim), myid, i, ierr
  integer :: source, dest
  logical :: isperiodic(p_max_dim+1)

  this%free_comm = .true.
  ndim = this%x_dim

  ! Get a cartesian topology through MPI
  dims(1:ndim) = this%nx(1:ndim)
  isperiodic(1:ndim) = this%ifpr_(1:ndim)

  call mpi_cart_create( mpi_comm_world, ndim, dims, isperiodic, &
           .true., this%comm, ierr )
  if ( ierr /= MPI_SUCCESS ) then
    ERROR("mpi_cart_create failed")
    call abort_program( p_err_mpi )
  endif

  ! this will probably reorder the nodes so we need to get our new node id
  call mpi_comm_rank( this%comm, myid, ierr )
  if ( ierr /= MPI_SUCCESS ) then
    ERROR("mpi_comm_rank failed")
    call abort_program( p_err_mpi )
  endif

  this%my_aid_ = myid + 1 ! we use 1 for the 1st node in osiris (i.e. fortran numbering)

  ! get local ngp
  call mpi_cart_coords( this%comm, myid, ndim, coords, ierr )
  if ( ierr /= MPI_SUCCESS ) then
    ERROR("mpi_cart_coords failed")
    call abort_program( p_err_mpi )
  endif
  this%my_ngp(1:this%x_dim) = coords(1:this%x_dim) + 1


  ! get neighbors; direction is 0 -> x1, 1 -> x2, 2 -> x3
  do i = 1, this%x_dim
    call mpi_cart_shift( this%comm, i-1, 1, source, dest, ierr)
    if ( ierr /= MPI_SUCCESS ) then
      ERROR("mpi_cart_shift failed")
      call abort_program( p_err_mpi )
    endif

    this%neighbor( p_lower, i ) = source + 1
    this%neighbor( p_upper, i ) = dest + 1
  enddo

  ! (* debug *) Test topology
  !call test_topology( this )

end subroutine create_topology_mpi
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Create simulation topology using original OSIRIS algorithm (not network optimized)
! Nodes will be distributed by x1 first, then x2, then x3.
!-----------------------------------------------------------------------------------------
subroutine create_topology_old( this )
!-----------------------------------------------------------------------------------------
  implicit none

  class( t_node_conf ), intent(inout)  ::  this

  integer :: n, stride, i

  this%comm = MPI_COMM_WORLD
  this%free_comm = .false.


  ! get local ngp
  n = this%my_aid() - 1
  do i=1, this%x_dim
    this%my_ngp(i) =  mod( n, this%nx(i) ) + 1
    n = n / this%nx(i)
  enddo

  ! Locate neighbors
  stride = 1
  do i = 1, this%x_dim
    ! lower
    if (this%my_ngp(i) == 1) then
      if ( this%ifpr_(i) ) then
        this%neighbor( p_lower, i ) = this%my_aid() + (this%nx(i)-1)*stride
      else
        this%neighbor( p_lower, i ) = -1
      endif
    else
      this%neighbor( p_lower, i ) = this%my_aid() - stride
    endif

    ! upper
    if (this%my_ngp(i) == this%nx(i)) then
      if ( this%ifpr_(i) ) then
        this%neighbor( p_upper, i ) = this%my_aid() - (this%nx(i)-1)*stride
      else
        this%neighbor( p_upper, i ) = -1
      endif
    else
      this%neighbor( p_upper, i ) = this%my_aid() + stride
    endif

    stride = stride*this%nx(i)
  enddo

  ! (* debug *) Test topology
  ! call test_topology( this )

end subroutine create_topology_old

#ifdef __bgq__

!-----------------------------------------------------------------------------------------
! Create simulation topology that matches the BlueGene/Q hardware topology
!-----------------------------------------------------------------------------------------
subroutine create_topology_bgq( this )
!-----------------------------------------------------------------------------------------
  implicit none

  class( t_node_conf ), intent(inout)  ::  this

  integer :: i,j,rank,ierr
  integer, dimension(p_max_dim) :: coords

  ! Initialize BlueGene/Q partition data
  call bgq_initpart( this%x_dim, this%nx, ierr )

  ! Failsafe, if for some reason the previous failed default to mpi topology
  if ( ierr /= 0 ) then
    if ( this%my_aid() == 1 ) then
      write(0,*) "(*warning*) Unable to initialize BlueGene/Q torus mapping, &
                  & defaulting to MPI topology"
    endif
    this%topology_type = p_topology_mpi
    call create_topology_mpi( this )
    return
  endif

  ! We use the global communicator
  this%comm = MPI_COMM_WORLD
  this%free_comm = .false.

  ! get local ngp
  call bgq_rank2sim( this%my_aid - 1, coords, ierr )
  if ( ierr /= 0 ) then
    ERROR("bgq_rank2sim failed")
    call abort_program( p_err_mpi )
  endif
  do i=1, this%x_dim
    this%my_ngp(i) = coords(i)
  enddo

  ! setup neighbours
  do i = 1, this%x_dim
    ! lower
    coords(i) = this%my_ngp(i) - 1
    if ( this%ifpr_(i) .and. coords(i) < 0 ) coords(i) = coords(i) + this%nx(i)
    if ( coords(i) >= 0 ) then
      call bgq_sim2rank( coords, rank, ierr )
      if ( ierr /= 0 ) then
        ERROR("bgq_sim2rank failed")
        call abort_program( p_err_mpi )
      endif
      this%neighbor( p_lower, i ) = rank + 1
    else
      this%neighbor( p_lower, i ) = -1
    endif

    ! upper
    coords(i) = this%my_ngp(i) + 1
    if ( this%ifpr_(i) .and. coords(i) >= this%nx(i) ) coords(i) = coords(i) - this%nx(i)
    if ( coords(i) < this%nx(i) ) then
      call bgq_sim2rank( coords, rank, ierr )
      if ( ierr /= 0 ) then
        ERROR("bgq_sim2rank failed")
        call abort_program( p_err_mpi )
      endif
      this%neighbor( p_upper, i ) = rank + 1
    else
      this%neighbor( p_upper, i ) = -1
    endif

    ! Clear shift
    coords(i) = this%my_ngp(i)
  enddo

  ! Convert to 1 indexed coordinates
  do i=1, this%x_dim
    this%my_ngp(i) = this%my_ngp(i) + 1
  enddo

  ! (* debug *) Test topology
  call test_topology( this )

end subroutine create_topology_bgq


#endif

!-----------------------------------------------------------------------------------------
! Test the topology for debug purposes
!-----------------------------------------------------------------------------------------
subroutine test_topology( this )

  implicit none

  class( t_node_conf ), intent(in)  ::  this

  integer :: i, j, ierr

  ! MPI variables
  integer, dimension(p_max_dim) :: lrecv, urecv, send
  integer :: lrreq, urreq, lsreq, usreq
  integer, dimension( MPI_STATUS_SIZE ) :: status

  if ( this%my_aid() == 1 ) then
    print *, '(*debug*) Testing topology...'
  endif

  ! Prepare data to send
  do i = 1, this%x_dim
    send(i) = this%my_ngp(i)
  enddo

  ! send / receive messages from neighbors
  do i = 1, this%x_dim

    if ( this%my_aid() == 1 ) then
      print *, '(*debug*) Testing dim ', i
    endif

    ! receive from lower neighbor
    if ( this%neighbor( p_lower, i ) > 0 ) then
      if ( this%neighbor( p_lower, i ) < 1 .or. &
        this%neighbor( p_lower, i ) > this%no_num ) then
        print *, '(*error*) Node: ', this%my_ngp(1:this%x_dim), &
        'invalid rank for lower neighbor : ', this%neighbor( p_lower, i )
        call abort_program( p_err_mpi )
      endif
      call mpi_irecv( lrecv, this%x_dim, MPI_INTEGER, this%neighbor( p_lower, i )-1, 0, &
        this%comm, lrreq, ierr )
      CHECK_ERROR( ierr, "Unable to post recv(lower)", p_err_mpi )
    endif
    ! receive from upper neighbor
    if ( this%neighbor( p_upper, i ) > 0 ) then
      if ( this%neighbor( p_upper, i ) < 1 .or. &
        this%neighbor( p_upper, i ) > this%no_num ) then
        print *, '(*error*) Node: ', this%my_ngp(1:this%x_dim), &
        'invalid rank for upper neighbor : ', this%neighbor( p_upper, i )
        call abort_program( p_err_mpi )
      endif
      call mpi_irecv( urecv, this%x_dim, MPI_INTEGER, this%neighbor( p_upper, i )-1, 0, &
        this%comm, urreq, ierr )
      CHECK_ERROR( ierr, "Unable to post recv(upper)", p_err_mpi )
    endif

    ! send to lower neighbor
    if ( this%neighbor( p_lower, i ) > 0 ) then
      call mpi_isend( send, this%x_dim, MPI_INTEGER, this%neighbor( p_lower, i )-1, 0, &
        this%comm, lsreq, ierr )
      CHECK_ERROR( ierr, "Unable to post send(lower)", p_err_mpi )
    endif
    ! send to upper neighbor
    if ( this%neighbor( p_upper, i ) > 0 ) then
      call mpi_isend( send, this%x_dim, MPI_INTEGER, this%neighbor( p_upper, i )-1, 0, &
        this%comm, usreq, ierr )
      CHECK_ERROR( ierr, "Unable to post send(upper)", p_err_mpi )
    endif

    ! Verify lower neighbor
    if ( this%neighbor( p_lower, i ) > 0 ) then
      call mpi_wait( lrreq, status, ierr )
      CHECK_ERROR( ierr, "Wait recv(lower) failed", p_err_mpi )
      lrecv( i ) = lrecv( i ) + 1
      if ( lrecv( i ) > this%nx(i) ) lrecv( i ) = lrecv( i ) - this%nx(i)

      do j = 1, this%x_dim
        if ( lrecv(j) /= this%my_ngp(j) ) then
          print *, '(*error*) Node: ', this%my_ngp(1:this%x_dim),' bad lower neighbor : ', &
                                       lrecv(1:this%x_dim)
          call abort_program( p_err_mpi )
        endif
      enddo
    endif

    ! Verify upper neighbor
    if ( this%neighbor( p_upper, i ) > 0 ) then
      call mpi_wait( urreq, status, ierr )
      CHECK_ERROR( ierr, "Wait recv(upper) failed", p_err_mpi )
      urecv( i ) = urecv( i ) - 1
      if ( urecv( i ) < 1 ) urecv( i ) = urecv( i ) + this%nx(i)

      do j = 1, this%x_dim
        if ( urecv(j) /= this%my_ngp(j) ) then
          print *, '(*error*) Node: ', this%my_ngp(1:this%x_dim),' bad upper neighbor.', &
                                       urecv(1:this%x_dim)
          call abort_program( p_err_mpi )
        endif
      enddo
    endif

    ! clear send messages
    if ( this%neighbor( p_lower, i ) > 0 ) call mpi_wait( lsreq, status, ierr )
    if ( this%neighbor( p_upper, i ) > 0 ) call mpi_wait( usreq, status, ierr )

  enddo

  ! Synchronize all nodes
  call this % barrier()

  if ( this%my_aid() == 1 ) then
    print *, '(*debug*) Topology is correct.'
  endif

end subroutine test_topology

!-----------------------------------------------------------------------------------------
subroutine write_checkpoint_node_conf( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       write object information into a restart file
!-----------------------------------------------------------------------------------------

  implicit none


  class( t_node_conf ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  integer :: ierr

  restart_io_wr( p_nconf_rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then
    ERROR('error writing restart data for node conf object.')
    call abort_program(p_err_rstwrt)
  endif


  ! no data is needed for restart, this is just used to insure
  ! that the input deck has not changed
  restart_io_wr( this%nx, restart_handle, ierr )
  if ( ierr/=0 ) then
    ERROR('error writing restart data for node_conf object.')
    call abort_program( p_err_rstwrt)
  endif

  restart_io_wr( this%ifpr_, restart_handle, ierr )
  if ( ierr/=0 ) then
    ERROR('error writing restart data for node_conf object.')
    call abort_program( p_err_rstwrt)
  endif

  restart_io_wr( this%my_ngp, restart_handle, ierr )
  if ( ierr/=0 ) then
    ERROR('error writing restart data for node_conf object.')
    call abort_program( p_err_rstwrt)
  endif


end subroutine write_checkpoint_node_conf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine restart_read_node_conf( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       read object information from a restart file
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_node_conf ), intent(inout)  ::  this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=len(p_nconf_rst_id)) :: rst_id
  integer, dimension(p_max_dim) :: nx, my_ngp
  logical, dimension(p_max_dim) :: ifpr
  integer :: i, ierr

  restart_io_rd( rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then
    ERROR('error reading restart data for node conf object.')
    ERROR('rst_id = ', rst_id)
    call abort_program(p_err_rstrd)
  endif

  ! check if restart file is compatible
  if ( rst_id /= p_nconf_rst_id) then
    ERROR('Corrupted restart file, or restart file ')
    ERROR('from incompatible binary (node conf)')
    ERROR('rst_id = ', rst_id)
    call abort_program(p_err_rstrd)
  endif

  restart_io_rd( nx, restart_handle, ierr )
  if ( ierr/=0 ) then
    ERROR('error reading restart data for node_conf object.')
    call abort_program(p_err_rstrd)
  endif

  restart_io_rd( ifpr, restart_handle, ierr )
  if ( ierr/=0 ) then
    ERROR('error reading restart data for node_conf object.')
    call abort_program(p_err_rstrd)
  endif

  restart_io_rd( my_ngp, restart_handle, ierr )
  if ( ierr/=0 ) then
    ERROR('error reading restart data for node_conf object.')
    call abort_program(p_err_rstrd)
  endif

  ! check if the user messed up the input file
  do i = 1, p_x_dim
    if ( nx(i) /= this%nx(i) ) then
      ERROR('The node_number specified on the input deck does')
      ERROR('not match the restart data')
      call abort_program(p_err_rstrd)
    endif

    if ( ifpr(i) .neqv. this%ifpr_(i) ) then
      ERROR('The if_periodic specified on the input deck does')
      ERROR('not match the restart data')
      call abort_program(p_err_rstrd)
    endif

    if ( my_ngp(i) /= this%my_ngp(i) ) then
      ERROR('NGP mismatch, reading data for different node')
      call abort_program(p_err_rstrd)
    endif

  enddo

end subroutine restart_read_node_conf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function n_threads_node_conf(this)
!-----------------------------------------------------------------------------------------
! result is number of threads per node
!-----------------------------------------------------------------------------------------

  implicit none

  integer :: n_threads_node_conf

  class( t_node_conf ), intent( in )  ::  this

#ifdef _OPENMP
  n_threads_node_conf = this%n_threads
#else
  n_threads_node_conf = 1
#endif

end function n_threads_node_conf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
pure function no_num_node_conf(this)
!-----------------------------------------------------------------------------------------
! result is total number of nodes
!-----------------------------------------------------------------------------------------

  implicit none

  integer :: no_num_node_conf

  class( t_node_conf ), intent( in )  ::  this

  no_num_node_conf=this%no_num

end function no_num_node_conf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
function get_nx(this)
!-----------------------------------------------------------------------------------------
! result is number of nodes in each direction
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_node_conf ), intent( in )  ::  this
  integer, dimension(this%x_dim) :: get_nx

  get_nx(1:this%x_dim)=this%nx(1:this%x_dim)

end function get_nx
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_nx_dim(this, dim)
!-----------------------------------------------------------------------------------------
! result is number of nodes in direction given by dim
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_node_conf ), intent( in )  ::  this
  integer, intent(in) :: dim

  integer :: get_nx_dim

  get_nx_dim=this%nx(dim)

end function get_nx_dim
!-----------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Returns true if node boundary is on the edge of the simulation box
!---------------------------------------------------------------------------------------------------
function on_edge_node_conf( this, dim, bnd )

  implicit none

  class( t_node_conf ), intent(in) :: this
  integer, intent(in) :: dim, bnd

  logical :: on_edge_node_conf

  if ( bnd == p_lower ) then
    on_edge_node_conf = ( this%my_ngp(dim) == 1 )
  else
    on_edge_node_conf = ( this%my_ngp(dim) == this%nx(dim) )
  endif

end function on_edge_node_conf
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Returns true if local node boundary is a physical boundary (non communication / periodic )
!---------------------------------------------------------------------------------------------------
function physical_boundary_local_nconf( this, dim, bnd )

  implicit none

  class( t_node_conf ), intent(in) :: this
  integer, intent(in) :: dim, bnd

  logical :: physical_boundary_local_nconf

  if ( this%ifpr_( dim )  ) then
    ! Periodic boundaries are not physical boundaries
    physical_boundary_local_nconf = .false.
  else
    ! To be a physical boundary the boundary must be on the edge of the simulation box
    if ( bnd == p_lower ) then
      physical_boundary_local_nconf = ( this%my_ngp(dim) == 1 )
    else
      physical_boundary_local_nconf = ( this%my_ngp(dim) == this%nx(dim) )
    endif
  endif

end function physical_boundary_local_nconf
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Returns true if  node boundary is a physical boundary (non communication / periodic )
!---------------------------------------------------------------------------------------------------
function physical_boundary_node_nconf( this, node, dim, bnd )

  implicit none

  class( t_node_conf ), intent(in) :: this
  integer, intent(in) :: node, dim, bnd

  logical :: physical_boundary_node_nconf

  if ( this%ifpr_( dim )  ) then
    ! Periodic boundaries are not physical boundaries
    physical_boundary_node_nconf = .false.
  else
    ! To be a physical boundary the boundary must be on the edge of the simulation box
    if ( bnd == p_lower ) then
      physical_boundary_node_nconf = ( this%ngp( node, dim ) == 1 )
    else
      physical_boundary_node_nconf = ( this%ngp( node, dim ) == this%nx(dim) )
    endif
  endif

end function physical_boundary_node_nconf
!---------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!  the result is the array id of the node on the
!  boundary described by bnd_idx and i_dim
!-----------------------------------------------------------------------------------------
function neighbor( this, bnd_idx, i_dim )

  implicit none

  integer :: neighbor

  class( t_node_conf ), intent(in) :: this
  integer, intent(in) :: bnd_idx, i_dim

  neighbor = this%neighbor( bnd_idx, i_dim )

end function neighbor
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! parallel rank of this node (Fortran ordering)
! ! This is depreacated and will be removed in the near future
! Please use local_rank instead
!-----------------------------------------------------------------------------------------
pure function my_aid( this )

  implicit none

  integer :: my_aid

  class( t_node_conf ), intent( in )  ::  this

  my_aid=this%my_aid_

end function my_aid
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Returns the position of the requested node on the parallel partition
!-----------------------------------------------------------------------------------------
function ngp_nconf( this, node, dim )

    implicit none

    class( t_node_conf ), intent( in ) :: this
    integer, intent( in ) :: node
    integer, intent(in) :: dim

    integer :: ngp_nconf

    integer, dimension( p_max_dim ) :: coords
    integer :: i, nt, ierr

    if ( this%no_num > 1 ) then

      select case ( this%topology_type )
      case ( p_topology_mpi )
        call mpi_cart_coords( this%comm, node - 1, this%x_dim, coords, ierr )
        if ( ierr /= MPI_SUCCESS ) then
          ERROR("mpi_cart_coords failed")
          call abort_program( p_err_mpi )
        endif
        ngp_nconf = coords(dim) + 1

      case ( p_topology_old )
        nt = node - 1
        do i=1, this%x_dim
          coords(i) = mod( nt, this%nx(i) ) + 1
          nt = nt / this%nx(i)
        enddo
        ngp_nconf = coords(dim)

#ifdef __bgq__
      case ( p_topology_bgq )
        call bgq_rank2sim( node-1, coords, ierr )
        if ( ierr /= 0 ) then
          ERROR("bgq_rank2sim failed")
          call abort_program( p_err_mpi )
        endif
        ngp_nconf = coords(dim) + 1
#endif
      case default
        ngp_nconf = 1
      end select
    else

      ngp_nconf = 1

    endif


end function ngp_nconf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! my N_ode G_rid P_osition
! the result is the position of this processor
! in the grid of processors
!-----------------------------------------------------------------------------------------
function my_ngp_all( this )

  implicit none

  class( t_node_conf ), intent( in )  ::  this
  integer, dimension(this%x_dim) :: my_ngp_all

  my_ngp_all(1:this%x_dim)=this%my_ngp(1:this%x_dim)

end function my_ngp_all
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!  my N_ode G_rid P_osition
!  the result is the position of this processor
!  in the grid of processors
!-----------------------------------------------------------------------------------------
function my_ngp_dim( this, dim )

  implicit none

  class( t_node_conf ), intent( in )  ::  this
  integer, intent(in) :: dim
  integer :: my_ngp_dim

  my_ngp_dim = this % my_ngp(dim)

end function my_ngp_dim
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function relevant_ifpr_nconf( this )
!-----------------------------------------------------------------------------------------
! gives flags for the periodic node configuration WHERE IT IS RELEVANT for the rest of
! the code. That means in directions which are set for periodicity and that have only one
! node if there is more then one node in that directions the rest of the code does not
! notice the periodicity in that direction since it just changes which nodes data are
! shipped to
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_node_conf ), intent(in)  ::  this
  logical, dimension(this%x_dim)    ::  relevant_ifpr_nconf

  integer :: i

  do i=1, this%x_dim
    if ( this%ifpr_(i) .and. ( this%nx(i) .eq. 1 ) ) then
      relevant_ifpr_nconf(i) = .true.
    else
      relevant_ifpr_nconf(i) = .false.
    endif
  enddo

end function relevant_ifpr_nconf
!-----------------------------------------------------------------------------------------



#ifdef HAS_MPI_IN_PLACE

!---------------------------------------------------------------------------------------------------
! MPI 2 reduce operations, using MPI_IN_PLACE, use less memory
!---------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reduce_array_int_1d( this, buffer, operation, all )
!-----------------------------------------------------------------------------------------
!       Execute specified global reduce operation and return result
!       on node 1. If no operation is specified the routine defaults
!       to p_sum i.e. sums up the data from all the nodes
!-----------------------------------------------------------------------------------------

  implicit none


  class( t_node_conf ), intent(in) :: this
  integer, dimension(:), intent(inout) :: buffer
  integer, intent(in), optional :: operation   ! which operation to do
  logical, intent(in), optional :: all         ! do a MPI_ALLREDUCE
  ! i.e. return result to
  ! all nodes


  integer, dimension(1) :: s

  integer :: loperation
  logical :: lall
  integer :: ierr
  integer :: count

  ! if single node run do nothing

  if (this%no_num > 1) then

    if (present(operation)) then
      loperation = operation
    else
      loperation = p_sum
    endif

    if (present(all)) then
      lall = all
    else
      lall = .false.
    endif

    ! allocate temporary buffer

    s(1) = size(buffer,1)
    count = s(1)

    if (lall) then
      ! return result to all nodes
      call MPI_ALLREDUCE(MPI_IN_PLACE, buffer, count, MPI_INTEGER, loperation, &
                          this%comm, ierr)
    else
      if ( this%my_aid() == 1 ) then
        call MPI_REDUCE(MPI_IN_PLACE, buffer, count, MPI_INTEGER, loperation, &
                          0, this%comm, ierr)
      else
        call MPI_REDUCE(buffer, 0, count, MPI_INTEGER, loperation, &
                          0, this%comm, ierr)
      endif
    endif

  endif


end subroutine reduce_array_int_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reduce_array_rsing_1d( this, buffer, operation )
!-----------------------------------------------------------------------------------------
!       Execute specified global reduce operation and return result
!       on node 1. If no operation is specified the routine defaults
!       to p_sum i.e. sums up the data from all the nodes
!-----------------------------------------------------------------------------------------

  implicit none


  class( t_node_conf ), intent(in) :: this
  real(p_single), dimension(:), intent(inout) :: buffer
  integer, intent(in), optional :: operation



  integer :: loperation
  integer :: ierr
  integer :: count

  ! if single node run do nothing

  if (this%no_num > 1) then

    if (present(operation)) then
      loperation = operation
    else
      loperation = p_sum
    endif

    ! allocate temporary buffer

    count = size(buffer,1)

    if ( this%my_aid() == 1 ) then
      call MPI_REDUCE(MPI_IN_PLACE, buffer, count, MPI_REAL, loperation, &
        0, this%comm, ierr)
    else
      call MPI_REDUCE(buffer, 0, count, MPI_REAL, loperation, &
        0, this%comm, ierr)
    endif

  endif


end subroutine reduce_array_rsing_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reduce_array_size_rsing_1d( this, buffer, size, operation )
!-----------------------------------------------------------------------------------------

implicit none

class( t_node_conf ), intent(in) :: this
real(p_single), dimension(:), intent(inout) :: buffer
integer, intent(in) :: size
integer, intent(in), optional :: operation

integer :: loperation
integer :: ierr

if (this%no_num > 1) then

  if (present(operation)) then
    loperation = operation
  else
    loperation = p_sum
  endif

  ! allocate temporary buffer

  if ( this%my_aid() == 1 ) then
    call MPI_REDUCE(MPI_IN_PLACE, buffer, size, MPI_REAL, loperation, &
      0, this%comm, ierr)
  else
    call MPI_REDUCE(buffer, 0, size, MPI_REAL, loperation, &
      0, this%comm, ierr)
  endif

endif

end subroutine reduce_array_size_rsing_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reduce_array_rdoub_1d( this, buffer, operation, all )
!-----------------------------------------------------------------------------------------
!       Execute specified global reduce operation and return result
!       on node 1. If no operation is specified the routine defaults
!       to p_sum i.e. sums up the data from all the nodes
!-----------------------------------------------------------------------------------------

implicit none


class( t_node_conf ), intent(in) :: this
real(p_double), dimension(:), intent(inout) :: buffer
integer, intent(in), optional :: operation
logical, intent(in), optional :: all


logical :: doall

! communication variables

integer :: loperation
integer :: ierr
integer :: count


! if single node run do nothing

if (this%no_num > 1) then

  if ( .not. present(all) ) then
    doall = .false.
  else
    doall = all
  endif

  if (present(operation)) then
    loperation = operation
  else
    loperation = p_sum
  endif

  ! allocate temporary buffer

  count = size(buffer,1)

  if ( doall ) then
    call MPI_ALLREDUCE(MPI_IN_PLACE, buffer, count, MPI_DOUBLE_PRECISION, loperation, &
      0, this%comm, ierr)
  else
    if ( this%my_aid() == 1 ) then
      call MPI_REDUCE(MPI_IN_PLACE, buffer, count, MPI_DOUBLE_PRECISION, loperation, &
        0, this%comm, ierr)
    else
      call MPI_REDUCE(buffer, 0, count, MPI_DOUBLE_PRECISION, loperation, &
        0, this%comm, ierr)
    endif
  endif
endif


end subroutine reduce_array_rdoub_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reduce_array_rsing_2d( this, buffer, operation )
!-----------------------------------------------------------------------------------------
!       Execute specified global reduce operation and return result
!       on node 1. If no operation is specified the routine defaults
!       to p_sum i.e. sums up the data from all the nodes
!-----------------------------------------------------------------------------------------

! this version uses MPI_IN_PLACE to save memory!

  implicit none


  class( t_node_conf ), intent(in) :: this
  real(p_single), dimension(:,:), intent(inout) :: buffer
  integer, intent(in), optional :: operation

  integer, dimension(2) :: s

  !       communication variables

  integer :: loperation
  integer :: ierr
  integer :: count


  ! if single node run do nothing

  if (this%no_num > 1) then

    if (present(operation)) then
      loperation = operation
    else
      loperation = p_sum
    endif

    ! allocate temporary buffer

    s(1) = size(buffer,1)
    s(2) = size(buffer,2)
    count = product(s)

    if ( this%my_aid() == 1 ) then
      call MPI_REDUCE(MPI_IN_PLACE, buffer, count, MPI_REAL, loperation, &
        0, this%comm, ierr)
    else
      call MPI_REDUCE(buffer, 0, count, MPI_REAL, loperation, &
        0, this%comm, ierr)
    endif

  endif

end subroutine reduce_array_rsing_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reduce_array_rdoub_2d( this, buffer, operation, all )
!-----------------------------------------------------------------------------------------
!       Execute specified global reduce operation and return result
!       on node 1. If no operation is specified the routine defaults
!       to p_sum i.e. sums up the data from all the nodes
!-----------------------------------------------------------------------------------------

! this version uses MPI_IN_PLACE to save memory!

  implicit none


  class( t_node_conf ), intent(in) :: this
  real(p_double), dimension(:,:), intent(inout) :: buffer
  integer, intent(in), optional :: operation
  logical, intent(in), optional :: all

  logical :: doall
  integer, dimension(2) :: s

  !       communication variables

  integer :: loperation
  integer :: ierr
  integer :: count


  ! if single node run do nothing

  if (this%no_num > 1) then

    if ( .not. present(all) ) then
      doall = .false.
    else
      doall = all
    endif

    if (present(operation)) then
      loperation = operation
    else
      loperation = p_sum
    endif

    ! allocate temporary buffer

    s(1) = size(buffer,1)
    s(2) = size(buffer,2)
    count = product(s)

    if ( doall ) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, buffer, count, MPI_DOUBLE_PRECISION, loperation, &
        0, this%comm, ierr)
    else
      if ( this%my_aid_== 1 ) then
        call MPI_REDUCE(MPI_IN_PLACE, buffer, count, MPI_DOUBLE_PRECISION, loperation, &
          0, this%comm, ierr)
      else
        call MPI_REDUCE(buffer, 0, count, MPI_DOUBLE_PRECISION, loperation, &
          0, this%comm, ierr)
      endif
    endif
  endif

end subroutine reduce_array_rdoub_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reduce_array_rsing_3d( this, buffer, operation )
!-----------------------------------------------------------------------------------------
!       Execute specified global reduce operation and return result
!       on node 1. If no operation is specified the routine defaults
!       to p_sum i.e. sums up the data from all the nodes
!-----------------------------------------------------------------------------------------

  implicit none


  class( t_node_conf ), intent(in) :: this
  real(p_single), dimension(:,:,:), intent(inout) :: buffer
  integer, intent(in), optional :: operation


  integer, dimension(3) :: s

  integer :: loperation
  integer :: ierr, count

  ! if single node run do nothing

  if (this%no_num > 1) then

    if (present(operation)) then
      loperation = operation
    else
      loperation = p_sum
    endif

    s(1) = size(buffer,1)
    s(2) = size(buffer,2)
    s(3) = size(buffer,3)
    count = product(s)

    if ( this%my_aid() == 1 ) then
      call MPI_REDUCE(MPI_IN_PLACE, buffer, count, MPI_REAL, loperation, &
        0, this%comm, ierr)
    else
      call MPI_REDUCE(buffer, 0, count, MPI_REAL, loperation, &
        0, this%comm, ierr)
    endif
  endif

end subroutine reduce_array_rsing_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reduce_array_rdoub_3d( this, buffer, operation, all )
!-----------------------------------------------------------------------------------------
!       Execute specified global reduce operation and return result
!       on node 1. If no operation is specified the routine defaults
!       to p_sum i.e. sums up the data from all the nodes
!-----------------------------------------------------------------------------------------

  implicit none


  class( t_node_conf ), intent(in) :: this
  real(p_double), dimension(:,:,:), intent(inout) :: buffer
  integer, intent(in), optional :: operation
  logical, intent(in), optional :: all

  logical :: doall
  integer, dimension(3) :: s

  integer :: loperation
  integer :: ierr, count

  ! if single node run do nothing

  if (this%no_num > 1) then

    if ( .not. present(all) ) then
      doall = .false.
    else
      doall = all
    endif

    if (present(operation)) then
      loperation = operation
    else
      loperation = p_sum
    endif

    s(1) = size(buffer,1)
    s(2) = size(buffer,2)
    s(3) = size(buffer,3)
    count = product(s)

    if ( doall ) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, buffer, count, MPI_DOUBLE_PRECISION, loperation, &
        0, this%comm, ierr)
    else
      if ( this%my_aid_== 1 ) then
        call MPI_REDUCE(MPI_IN_PLACE, buffer, count, MPI_DOUBLE_PRECISION, loperation, &
          0, this%comm, ierr)
      else
        call MPI_REDUCE(buffer, 0, count, MPI_DOUBLE_PRECISION, loperation, &
          0, this%comm, ierr)
      endif
    endif
  endif

end subroutine reduce_array_rdoub_3d
!-----------------------------------------------------------------------------------------


#else

!---------------------------------------------------------------------------------------------------
! MPI 1 reduce operations, not using MPI_IN_PLACE => use more memory
!---------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine reduce_array_int_1d( this, buffer_in, operation, all )
!-----------------------------------------------------------------------------------------
!       Execute specified global reduce operation and return result
!       on node 1. If no operation is specified the routine defaults
!       to p_sum i.e. sums up the data from all the nodes
!-----------------------------------------------------------------------------------------

  implicit none


  class( t_node_conf ), intent(in) :: this
  integer, dimension(:), intent(inout) :: buffer_in
  integer, intent(in), optional :: operation   ! which operation to do
  logical, intent(in), optional :: all         ! do a MPI_ALLREDUCE
  ! i.e. return result to
  ! all nodes


  integer, dimension(:), pointer :: buffer_out
  integer, dimension(1) :: s

  integer :: loperation
  logical :: lall
  integer :: ierr
  integer :: count



  ! if single node run do nothing

  if (this%no_num > 1) then

    if (present(operation)) then
      loperation = operation
    else
      loperation = p_sum
    endif

    if (present(all)) then
      lall = all
    else
      lall = .false.
    endif

    ! allocate temporary buffer

    s(1) = size(buffer_in,1)
    call alloc( buffer_out, s)
    count = s(1)

    if (lall) then
      ! return result to all nodes
      call MPI_ALLREDUCE(buffer_in, buffer_out, count, MPI_INTEGER, loperation, &
        this%comm, ierr)
      buffer_in=buffer_out
    else
      ! return result to node 0 only
      call MPI_REDUCE(buffer_in, buffer_out, count, MPI_INTEGER, loperation, &
        0, this%comm, ierr)
      if (this%my_aid() == 1) then
        buffer_in=buffer_out
      endif
    endif

    call freemem( buffer_out )
  endif

end subroutine reduce_array_int_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reduce_array_rsing_1d( this, buffer_in, operation )
!-----------------------------------------------------------------------------------------
!       Execute specified global reduce operation and return result
!       on node 1. If no operation is specified the routine defaults
!       to p_sum i.e. sums up the data from all the nodes
!-----------------------------------------------------------------------------------------

  implicit none


  class( t_node_conf ), intent(in) :: this
  real(p_single), dimension(:), intent(inout) :: buffer_in
  integer, intent(in), optional :: operation



  real(p_single), dimension(:), pointer :: buffer_out
  integer, dimension(1) :: s

  integer :: loperation
  integer :: ierr
  integer :: count



  ! if single node run do nothing

  if (this%no_num > 1) then

    if (present(operation)) then
      loperation = operation
    else
      loperation = p_sum
    endif

    ! allocate temporary buffer

    s(1) = size(buffer_in,1)
    call alloc( buffer_out, s )

    count = s(1)

    call MPI_REDUCE(buffer_in, buffer_out, count, MPI_REAL, loperation, &
      0, this%comm, ierr)
    if (this%my_aid() == 1) then
      buffer_in=buffer_out
    endif

    call freemem( buffer_out )
  endif

end subroutine reduce_array_rsing_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reduce_array_size_rsing_1d( this, buffer, size, operation )

  implicit none

  class( t_node_conf ), intent(in) :: this
  real(p_single), dimension(:), intent(inout) :: buffer
  integer, intent(in) :: size
  integer, intent(in), optional :: operation

  real(p_single), dimension(:), pointer :: buffer_out
  integer :: loperation
  integer :: ierr

  if (this%no_num > 1) then

    if (present(operation)) then
      loperation = operation
    else
      loperation = p_sum
    endif

    ! allocate temporary buffer

    call alloc( buffer_out, (/size/) )

    call MPI_REDUCE(buffer, buffer_out, size, MPI_REAL, loperation, &
      0, this%comm, ierr)
    if (this%my_aid() == 1) then
      buffer(1:size)=buffer_out(1:size)
    endif

    call freemem( buffer_out )
  endif

end subroutine reduce_array_size_rsing_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reduce_array_rdoub_1d( this, buffer_in, operation, all )
!-----------------------------------------------------------------------------------------
!       Execute specified global reduce operation and return result
!       on node 1. If no operation is specified the routine defaults
!       to p_sum i.e. sums up the data from all the nodes
!-----------------------------------------------------------------------------------------

  implicit none


  class( t_node_conf ), intent(in) :: this
  real(p_double), dimension(:), intent(inout) :: buffer_in
  integer, intent(in), optional :: operation
  logical, intent(in), optional :: all


  logical :: doall

  ! communication variables

  integer :: loperation
  real( p_double ), dimension(:), pointer :: buffer_out
  integer :: ierr
  integer :: count


  ! if single node run do nothing

  if (this%no_num > 1) then

    if ( .not. present(all) ) then
      doall = .false.
    else
      doall = all
    endif

    if (present(operation)) then
      loperation = operation
    else
      loperation = p_sum
    endif

    ! allocate temporary buffer

    count = size(buffer_in,1)
    call alloc( buffer_out, (/ count /))

    if ( doall ) then
      ! return result to all nodes
      call MPI_ALLREDUCE(buffer_in, buffer_out, count, MPI_DOUBLE_PRECISION, loperation, &
        this%comm, ierr)
      buffer_in = buffer_out

    else
      ! return result to node 0 only
      call MPI_REDUCE(buffer_in, buffer_out, count, MPI_DOUBLE_PRECISION, loperation, &
        0, this%comm, ierr)

      if (this%my_aid() == 1) then
        buffer_in=buffer_out
      endif

    endif

    call freemem( buffer_out )
  endif


end subroutine reduce_array_rdoub_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reduce_array_size_rdoub_1d( this, buffer, size, operation )

    implicit none

    class( t_node_conf ), intent(in) :: this
    real(p_double), dimension(:), intent(inout) :: buffer
    integer, intent(in) :: size
    integer, intent(in), optional :: operation

    real(p_double), dimension(:), pointer :: buffer_out
    integer :: loperation
    integer :: ierr

    if (this%no_num > 1) then

      if (present(operation)) then
        loperation = operation
      else
        loperation = p_sum
      endif

      ! allocate temporary buffer

      call alloc( buffer_out, (/size/) )

      call MPI_REDUCE(buffer, buffer_out, size, MPI_DOUBLE_PRECISION, loperation, &
        0, this%comm, ierr)
      if (this%my_aid() == 1) then
        buffer(1:size)=buffer_out(1:size)
      endif

      call freemem( buffer_out )
    endif

end subroutine reduce_array_size_rdoub_1d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine reduce_array_rsing_2d( this, buffer_in, operation )
!-----------------------------------------------------------------------------------------
!       Execute specified global reduce operation and return result
!       on node 1. If no operation is specified the routine defaults
!       to p_sum i.e. sums up the data from all the nodes
!-----------------------------------------------------------------------------------------

implicit none


class( t_node_conf ), intent(in) :: this
real(p_single), dimension(:,:), intent(inout) :: buffer_in
integer, intent(in), optional :: operation

integer, dimension(2) :: s

!       communication variables
real(p_single), dimension(:,:), pointer :: buffer_out

integer :: loperation
integer :: ierr
integer :: count


! if single node run do nothing

if (this%no_num > 1) then

  if (present(operation)) then
    loperation = operation
  else
    loperation = p_sum
  endif

  ! allocate temporary buffer

  s(1) = size(buffer_in,1)
  s(2) = size(buffer_in,2)
  call alloc( buffer_out, s )

  count = product(s)

  call MPI_REDUCE(buffer_in, buffer_out, count, MPI_REAL, loperation, &
    0, this%comm, ierr)
  if (this%my_aid() == 1) then
    buffer_in=buffer_out
  endif

  call freemem( buffer_out )

endif

end subroutine reduce_array_rsing_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reduce_array_rdoub_2d( this, buffer_in, operation, all )
!-----------------------------------------------------------------------------------------
!       Execute specified global reduce operation and return result
!       on node 1. If no operation is specified the routine defaults
!       to p_sum i.e. sums up the data from all the nodes
!-----------------------------------------------------------------------------------------

  implicit none


  class( t_node_conf ), intent(in) :: this
  real(p_double), dimension(:,:), intent(inout) :: buffer_in
  integer, intent(in), optional :: operation
  logical, intent(in), optional :: all

  logical :: doall
  integer, dimension(2) :: s

  ! communication variables
  real(p_double), dimension(:,:), pointer :: buffer_out

  integer :: loperation
  integer :: ierr
  integer :: count


  ! if single node run do nothing

  if (this%no_num > 1) then

    if ( .not. present(all) ) then
      doall = .false.
    else
      doall = all
    endif

    if (present(operation)) then
      loperation = operation
    else
      loperation = p_sum
    endif

    ! allocate temporary buffer

    s(1) = size(buffer_in,1)
    s(2) = size(buffer_in,2)
    call alloc( buffer_out, s )

    count = product(s)

    if ( doall ) then
      ! return result to all nodes
      call MPI_ALLREDUCE(buffer_in, buffer_out, count, MPI_DOUBLE_PRECISION, loperation, &
        this%comm, ierr)
      buffer_in = buffer_out

    else
      ! return result to node 0 only
      call MPI_REDUCE(buffer_in, buffer_out, count, MPI_DOUBLE_PRECISION, loperation, &
        0, this%comm, ierr)
      if (this%my_aid_== 1) then
        buffer_in=buffer_out
      endif

    endif

    call freemem( buffer_out )

  endif

end subroutine reduce_array_rdoub_2d
!-----------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Execute specified global reduce operation on a 3D single precision array and return result on node
! 1. If no operation is pecified the routine defaults to p_sum i.e. sums up the data from all the
! nodes.
! This version will split the calculation over sets of x1*x2 slices to conserve memory if the total
! buffer size is larger than 4 MB
!---------------------------------------------------------------------------------------------------
subroutine reduce_array_rsing_3d( this, data, operation )

  implicit none

  real(p_single), dimension(:,:,:), intent(inout) :: data
  class( t_node_conf ), intent(in) :: this
  integer, intent(in), optional :: operation

  ! Target Buffer size to use. This corresponds to 4 MB
  integer, parameter :: p_buf_size = 1048576
  integer :: loperation
  integer :: i, bsize, chunk, ierr
  integer, dimension(3) :: s
  real(p_single), dimension(:), pointer :: buffer


  ! if single node run do nothing
  if (this%no_num > 1) then

    ! Choose reduce operation, default to sum
    if (present(operation)) then
      loperation = operation
    else
      loperation = p_sum
    endif

    ! Get size of required buffer and chunk
    s(1) = size(data,1)
    s(2) = size(data,2)
    s(3) = size(data,3)

    if ( size(data) <= p_buf_size ) then
      chunk = s(3)
      bsize  = size(data)
    else
      chunk = p_buf_size / ( s(1)*s(2) )
      if ( chunk < 1 ) then
        chunk = 1
      endif
    endif
    bsize  = chunk * s(1)*s(2)

    ! allocate temporary buffer
    call alloc( buffer, (/bsize/) )

    ! Loop through all chunks
    i = 1
    do
      ! if all chunks have been processed then exit
      if ( i > s(3) ) exit

      ! last chunk may be smaller
      if ( i + chunk - 1 > s(3) ) then
        bsize = s(3) - i + 1
      else
        bsize = chunk
      endif
      bsize = bsize * s(1) * s(2)

      call mpi_reduce( data(:,:,i:), buffer, bsize, MPI_REAL, loperation, 0, this%comm, ierr )
      if ( ierr /= MPI_SUCCESS ) then
        ERROR('MPI Error')
        call abort_program( p_err_mpi )
      endif

      ! copy data to output array
      if ( this%my_aid() == 1 ) then
        call memcpy( data(:,:,i:), buffer, bsize )
      endif

      i = i + chunk
    enddo

    call freemem( buffer )

  endif

end subroutine reduce_array_rsing_3d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Execute specified global reduce operation on a 3D double precision array and return result on node
! 1. If no operation is pecified the routine defaults to p_sum i.e. sums up the data from all the
! nodes.
! This version will split the calculation over sets of x1*x2 slices to conserve memory if the total
! buffer size is larger than 4 MB
!---------------------------------------------------------------------------------------------------
subroutine reduce_array_rdoub_3d( this, data, operation, all )

  implicit none

  real(p_double), dimension(:,:,:), intent(inout) :: data
  class( t_node_conf ), intent(in) :: this
  integer, intent(in), optional :: operation
  logical, intent(in), optional :: all

  ! Target Buffer size to use. This corresponds to 4 MB
  integer, parameter :: p_buf_size = 1048576
  integer :: loperation
  integer :: i, bsize, chunk, ierr
  logical :: doall
  integer, dimension(3) :: s
  real(p_double), dimension(:), pointer :: buffer


  ! if single node run do nothing
  if (this%no_num > 1) then

    if ( .not. present(all) ) then
      doall = .false.
    else
      doall = all
    endif

    ! Choose reduce operation, default to sum
    if (present(operation)) then
      loperation = operation
    else
      loperation = p_sum
    endif

    ! Get size of required buffer and chunk
    s(1) = size(data,1)
    s(2) = size(data,2)
    s(3) = size(data,3)

    if ( size(data) <= p_buf_size ) then
      chunk = s(3)
      bsize  = size(data)
    else
      chunk = p_buf_size / ( s(1)*s(2) )
      if ( chunk < 1 ) then
        chunk = 1
      endif
    endif
    bsize  = chunk * s(1)*s(2)

    ! allocate temporary buffer
    call alloc( buffer, (/bsize/) )

    ! Loop through all chunks
    i = 1
    do
      ! if all chunks have been processed then exit
      if ( i > s(3) ) exit

      ! last chunk may be smaller
      if ( i + chunk - 1 > s(3) ) then
        bsize = s(3) - i + 1
      else
        bsize = chunk
      endif
      bsize = bsize * s(1) * s(2)

      if ( doall ) then
        ! return result to all nodes
        call mpi_allreduce( data(:,:,i:), buffer, bsize, MPI_DOUBLE_PRECISION, loperation, this%comm, ierr )
        if ( ierr /= MPI_SUCCESS ) then
          ERROR('MPI Error')
          call abort_program( p_err_mpi )
        endif
        call memcpy( data(:,:,i:), buffer, bsize )

      else

        ! return result to node 0 only
        call mpi_reduce( data(:,:,i:), buffer, bsize, MPI_DOUBLE_PRECISION, loperation, 0, this%comm, ierr )
        if ( ierr /= MPI_SUCCESS ) then
          ERROR('MPI Error')
          call abort_program( p_err_mpi )
        endif

        ! copy data to output array
        if ( this%my_aid_== 1 ) then
          call memcpy( data(:,:,i:), buffer, bsize )
        endif

      endif

      i = i + chunk
    enddo

    call freemem( buffer )

  endif

end subroutine reduce_array_rdoub_3d
!---------------------------------------------------------------------------------------------------


#endif



!---------------------------------------------------------------------------------------------------
! Parallel reduction of integer scalar value
!---------------------------------------------------------------------------------------------------
subroutine reduce_int( no_co, buffer, operation, all )

  implicit none

  class( t_node_conf ), intent(in) :: no_co
  integer,  intent(inout) :: buffer
  integer, intent(in) :: operation             ! which operation to do
  logical, intent(in), optional :: all         ! do a MPI_ALLREDUCE
  ! i.e. return result to
  ! all nodes

  integer, dimension(1) :: outval, inval
  integer :: ierr
  logical :: lall

  if ( no_co%no_num > 1 ) then

    if (present(all)) then
      lall = all
    else
      lall = .false.
    endif

    inval = buffer

    if ( lall ) then
      ! return result to all nodes
      call MPI_ALLREDUCE(inval, outval, 1, MPI_INTEGER, operation, &
        no_co%comm, ierr)
      buffer=outval(1)
    else
      ! return result to root node only
      call MPI_REDUCE(inval, outval, 1, MPI_INTEGER, operation, &
        0, no_co%comm, ierr)
      if ( no_co%my_aid() == 1 ) buffer = outval(1)
    endif

  endif

end subroutine reduce_int
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Parallel reduction of 64 bits integer scalar value
!---------------------------------------------------------------------------------------------------
subroutine reduce_int64( no_co, buffer, operation, all )

  implicit none

  class( t_node_conf ), intent(in) :: no_co
  integer( p_int64),  intent(inout) :: buffer
  integer, intent(in) :: operation             ! which operation to do
  logical, intent(in), optional :: all         ! do a MPI_ALLREDUCE
  ! i.e. return result to
  ! all nodes

  integer( p_int64), dimension(1) :: outval, inval
  integer :: ierr
  logical :: lall

  if ( no_co%no_num > 1 ) then

    if (present(all)) then
      lall = all
    else
      lall = .false.
    endif

    inval = buffer

    if ( lall ) then
      ! return result to all nodes
      call MPI_ALLREDUCE(inval, outval, 1, MPI_INTEGER8, operation, &
        no_co%comm, ierr)
      buffer=outval(1)
    else
      ! return result to root node only
      call MPI_REDUCE(inval, outval, 1, MPI_INTEGER8, operation, &
        0, no_co%comm, ierr)
      if ( no_co%my_aid() == 1 ) buffer = outval(1)
    endif

  endif

end subroutine reduce_int64
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Parallel reduction of single precision scalar value
!---------------------------------------------------------------------------------------------------
subroutine reduce_rsing( no_co, buffer, operation, all )

  implicit none

  class( t_node_conf ), intent(in) :: no_co
  real(p_single),  intent(inout) :: buffer
  integer, intent(in) :: operation            ! which operation to do
  logical, intent(in), optional :: all        ! do a MPI_ALLREDUCE
                                              ! i.e. return result to
                                              ! all nodes

  ! local variables
  real( p_single ), dimension(1) :: outval, inval
  integer :: ierr
  logical :: lall

  if ( no_co%no_num > 1 ) then

    if (present(all)) then
      lall = all
    else
      lall = .false.
    endif

    inval = buffer

    if ( lall ) then
      ! return result to all nodes
      call MPI_ALLREDUCE(inval, outval, 1, MPI_REAL, operation, &
        no_co%comm, ierr)
      buffer=outval(1)
    else
      ! return result to root node only
      call MPI_REDUCE(inval, outval, 1, MPI_REAL, operation, &
        0, no_co%comm, ierr)
      if ( no_co%my_aid() == 1 ) buffer = outval(1)
    endif

  endif

end subroutine reduce_rsing
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Parallel reduction of double precision scalar value
!---------------------------------------------------------------------------------------------------
subroutine reduce_rdoub( no_co, buffer, operation, all )

  implicit none

  class( t_node_conf ), intent(in) :: no_co
  real(p_double),  intent(inout) :: buffer
  integer, intent(in) :: operation            ! which operation to do
  logical, intent(in), optional :: all        ! do a MPI_ALLREDUCE
                                              ! i.e. return result to
                                              ! all nodes

  ! local variables
  real( p_double ), dimension(1) :: outval, inval
  integer :: ierr
  logical :: lall

  if ( no_co%no_num > 1 ) then

    if (present(all)) then
      lall = all
    else
      lall = .false.
    endif

    inval = buffer

    if ( lall ) then
      ! return result to all nodes
      call MPI_ALLREDUCE(inval, outval, 1, MPI_DOUBLE_PRECISION, operation, &
        no_co%comm, ierr)
      buffer=outval(1)
    else
      ! return result to root node only
      call MPI_REDUCE(inval, outval, 1, MPI_DOUBLE_PRECISION, operation, &
        0, no_co%comm, ierr)
      if ( no_co%my_aid() == 1 ) buffer = outval(1)
    endif

  endif

end subroutine reduce_rdoub
!---------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Synchronize all nodes
!-----------------------------------------------------------------------------------------
subroutine barrier( this )

  implicit none

  class( t_node_conf ), intent(in)  ::  this
  integer :: ierr

  if ( this%no_num > 1 ) then
    call MPI_BARRIER( this%comm, ierr )
  endif

end subroutine barrier
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! routine that sends a short message (ping) to another node
!-----------------------------------------------------------------------------------------
subroutine send_ping( this, dest_node, tag )

  implicit none

  class( t_node_conf ), intent(in)  ::  this
  integer, intent(in) :: dest_node
  integer, intent(in) :: tag

  integer :: ierr

  call MPI_SEND(ping_code, ping_size, mpi_integer, dest_node, tag, this%comm, ierr)

end subroutine send_ping
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       routine that receives a short message (ping)
!       from another node
!-----------------------------------------------------------------------------------------
subroutine recv_ping( this, source_node, tag )

  implicit none

  class( t_node_conf ), intent(in)  ::  this
  integer, intent(in) :: source_node
  integer, intent(in) :: tag

  integer ::  ierr
  integer, dimension(ping_size) :: ping_buffer
  integer, dimension(mpi_status_size):: stat

  call MPI_RECV(ping_buffer, ping_size, mpi_integer, source_node, tag, this%comm, stat, ierr)

end subroutine recv_ping
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       routine that gathers an integer array from all
!       nodes on node 0
!-----------------------------------------------------------------------------------------
subroutine gather_array_int_1d( this, local, all )

  implicit none


  class( t_node_conf ), intent(in)  ::  this
  integer, dimension(:),    intent(in)  :: local
  integer, dimension(:,:),  intent(out) :: all


  integer :: root, sendcount, recvcount, ierr


  if ( this%no_num > 1 ) then

    sendcount = size(local)
    recvcount = sendcount
    root = 0

    call MPI_GATHER( local, sendcount, MPI_INTEGER, all, recvcount, MPI_INTEGER, &
      root, this%comm, ierr)
  else

    all(:,1) = local(:)

  endif

end subroutine gather_array_int_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine allgather_int( this, send, recv )

  implicit none

  class( t_node_conf ), intent(in)  ::  this
  integer, intent(in)  :: send
  integer, dimension(:),  intent(out) :: recv


  integer :: sendcount, recvcount, ierr


  if ( this%no_num > 1 ) then

    sendcount = 1
    recvcount = 1

    call MPI_ALLGATHER( send, sendcount, MPI_INTEGER, &
      recv, recvcount, MPI_INTEGER, &
      this%comm, ierr)
  else

    recv(1) = send

  endif

end subroutine allgather_int
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!       routine that broadcasts an integer from
!       node 0 to all nodes
!-----------------------------------------------------------------------------------------
subroutine broadcast_int( this, n )

  implicit none

  class( t_node_conf ), intent(in)     ::  this
  integer,             intent(inout)  ::  n

  integer :: root, sendcount, ierr

  if ( this%no_num > 1) then

    sendcount = 1
    root = 0

    call MPI_BCAST( n, sendcount, MPI_INTEGER, root, this%comm, ierr)

  endif

end subroutine broadcast_int
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure function periodic_dir( this, dir )

  implicit none

  class( t_node_conf ), intent(in)     ::  this
  integer, intent(in) :: dir
  logical :: periodic_dir

  periodic_dir = this%ifpr_( dir )

end function periodic_dir
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function periodic_all( this )

  implicit none

  class( t_node_conf ), intent(in)     ::  this

  ! this breaks the intel ifort 9.1 compiler
  ! logical, dimension( this%x_dim ) :: periodic_all

  logical, dimension( p_x_dim ) :: periodic_all

  periodic_all(1:this%x_dim) = this%ifpr_(1:this%x_dim)

end function periodic_all
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! returns the MPI communicator for this object
!-----------------------------------------------------------------------------------------
pure function comm(this)

  implicit none

  class( t_node_conf ), intent(in)     ::  this
  integer :: comm

  comm = this%comm

end function comm
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! returns a unique id integer based on the node grid position
!-----------------------------------------------------------------------------------------
pure function ngp_id_nconf(this)

  implicit none

  class( t_node_conf ), intent(in)     ::  this
  integer :: i,st,ngp_id_nconf

  st = 1
  ngp_id_nconf = 1

  do i = 1, this%x_dim
    ngp_id_nconf = ngp_id_nconf + ( this%my_ngp(i)-1 ) * st
    st = st * this%nx(i)
  enddo

end function ngp_id_nconf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
pure function root(this)

  implicit none

  class( t_node_conf ), intent(in) :: this
  logical :: root

  root = ( this%my_aid() == 1 )

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! New class mehtods
!
! The results are 0 indexed, as is used by mpi functions
! Currently they only support MPI topologies
!
! local_rank      ! MPI rank of local node
! local_coords    ! Coordinates of local node in parallel partition
! rank_coords     ! Coordinates of specified rank in parallel partition
! coords_rank     ! MPI rank of the node with the specified coordinates
!-----------------------------------------------------------------------------------------

function local_rank( this )

  implicit none

  class(t_node_conf), intent(in) :: this
  integer :: local_rank, ierr

  if ( this % comm /= MPI_COMM_NULL ) then
    call mpi_comm_rank( this % comm, local_rank, ierr )
  else
    local_rank = 0
  endif

end function local_rank

function local_coords( this )

  implicit none

  class(t_node_conf), intent(in) :: this
  integer, dimension(p_x_dim) :: local_coords
  integer :: rank, ierr

  if ( this % comm /= MPI_COMM_NULL ) then
    call mpi_comm_rank( this % comm, rank, ierr )
    call mpi_cart_coords( this%comm, rank, p_x_dim, local_coords, ierr )
  else
    local_coords = 0
  endif

end function local_coords

function rank_coords( this, rank )

  implicit none

  class(t_node_conf), intent(in) :: this
  integer, intent(in) :: rank
  integer, dimension(p_x_dim) :: rank_coords
  integer :: ierr

  if ( this%comm /= MPI_COMM_NULL ) then
    call mpi_cart_coords( this%comm, rank, p_x_dim, rank_coords, ierr )
  else
    rank_coords = 0
  endif

end function rank_coords

function coords_rank( this, coords )

  implicit none

  class(t_node_conf), intent(in) :: this
  integer, dimension(:), intent(in) :: coords
  integer :: coords_rank, ierr

  if ( this%comm /= MPI_COMM_NULL ) then
    call mpi_cart_rank( this%comm, coords, coords_rank, ierr )
  else
    coords_rank = 0
  endif

end function coords_rank

function comm_size( this )

  implicit none

  class(t_node_conf), intent(in) :: this
  integer :: comm_size
  integer :: ierr

  if ( this%comm /= MPI_COMM_NULL ) then
    call mpi_comm_size( this%comm, comm_size, ierr )
  else
    comm_size = 1
  endif

end function comm_size

end module m_node_conf
