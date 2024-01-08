#include "os-config.h"
#include "os-preprocess.fpp"

module m_psource_file

#include "memory/memory.h"

use m_species_define
use m_parameters
use m_current_define
use m_psource_std
use m_node_conf

implicit none

private


type, extends(t_psource_std) :: t_psource_file

  ! flag for RAW file import
  logical :: raw_read_done = .false.

  ! parameter to shift x1 offset
  real(p_k_part) :: x1_offset

  ! q scaling factor for different grid resolutions
  real(p_k_part) :: q_mult

  ! path to file
  character(len = p_max_filename_len) :: file_name

  ! node_maps
  integer, dimension(:), pointer :: node_map_1d => null()
  integer, dimension(:,:), pointer :: node_map_2d => null()
  integer, dimension(:,:,:), pointer :: node_map_3d => null()

contains

  procedure :: if_inject => if_inject_file
  procedure :: read_input => read_input_file
  procedure :: cleanup => cleanup_file
  procedure :: inject => inject_file
  procedure :: read_raw_chunk => read_raw_chunk
  procedure :: send_particle_chunk => send_particle_chunk
  procedure :: inject_particle_chunk => inject_particle_chunk
end type t_psource_file

public :: t_psource_file

contains

!-----------------------------------------------------------------------------------------
! Read information from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_file( this, input_file, coordinates )

  use m_input_file

  implicit none

  class( t_psource_file ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: coordinates

  character(len = p_max_filename_len) :: file_name
  real(p_k_part) :: x1_offset, q_mult
  namelist /nl_profile/ file_name, x1_offset, q_mult
  integer :: ierr

  x1_offset = 0
  q_mult = 1.0_p_k_part
  ! Get namelist text from input file
  call get_namelist( input_file, "nl_profile", ierr )

  if ( ierr == 0 ) then
    read (input_file%nml_text, nml = nl_profile, iostat = ierr)
    if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) ""
        write(0,*) "   Error reading profile information"
        write(0,*) "   aborting..."
      endif
      stop
    endif
  else
    SCR_ROOT("   - profile parameters missing, using uniform density")
  endif

  this%x1_offset = x1_offset
  this%q_mult = q_mult
  this%file_name = file_name



end subroutine read_input_file

!-----------------------------------------------------------------------------------------
! Cleanup t_psource object
!-----------------------------------------------------------------------------------------
subroutine cleanup_file( this )

  implicit none

  class( t_psource_file ), intent( inout ) :: this

  ! nothing to clean

end subroutine cleanup_file
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Returns .true. if this psource injects particles
!-----------------------------------------------------------------------------------------
function if_inject_file( this ) result( if_inject )

  implicit none

  logical :: if_inject
  class( t_psource_file ), intent( in ) :: this

  if_inject = .not. this % raw_read_done

end function if_inject_file
!-----------------------------------------------------------------------------------------


#ifdef HDF5

!-----------------------------------------------------------------------------------------
! RAW Initialization HDF5 routines
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Inject particles using information from external file
!-----------------------------------------------------------------------------------------
function inject_file( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, &
                      send_msg, recv_msg ) result(num_inj)

  use hdf5_util
  use m_parameters
  use m_system
  use m_diagfile_hdf5
  use m_diagfile
  use hdf5
  use m_diagnostic_utilities
  use stringutil
  implicit none

  class(t_psource_file), intent(inout) :: this
  class(t_species), intent(inout), target :: species
  integer, dimension(:, :), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout)   :: jay
  class( t_node_conf ), intent(in)     :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: num_inj, chunk_inj, total_part_num, niter, npar
  integer :: part_chunk_size = 1e7

  integer, dimension(:, :), allocatable :: ig_xbnd_inj_temp
  integer :: ierr, nquants
  integer(hid_t) :: file_id      ! File identifier
  integer(hid_t) :: dspace, dset_id
  integer :: ping, handle, rank, i, j, k
  ! Status variable for messages

  integer, dimension( MPI_STATUS_SIZE ) :: stat
  real(p_double), dimension(p_x_dim) :: dx
  real(p_double), dimension(2, p_x_dim) :: g_xbnd
  real(p_k_part), dimension(p_x_dim) :: xpos
  integer, dimension(p_x_dim) :: g_nx, idxpos
  integer :: num_read
  real(p_k_part), dimension(:,:), pointer :: chunk_data, node_data
  integer(hsize_t), dimension(1) :: dims , maxdims
  integer :: broadcast_type, ltype
  real :: t1,t2

  chunk_data => null(); node_data => null()

  if(this%raw_read_done) then

    num_inj = this % t_psource_std % inject( species, ig_xbnd_inj, jay, no_co, &
                                             bnd_cross, node_cross, send_msg, recv_msg )
  else
    ! x, p, q quantities needed
    nquants = species % get_n_x_dims() + p_p_dim + 1
    do i = 1, p_x_dim
      ! get cell size
      dx(i) = species%dx(i)
      g_xbnd( :, i ) = species%g_box( :, i ) - 0.5_p_k_part * dx(i)
      g_nx(i) = species%g_nx(i)


    enddo
    ! first pass create node map & particle count for each node
    allocate(ig_xbnd_inj_temp(2,p_x_dim))
    if(mpi_node() == 0) then

      !initialize node map
      select case (p_x_dim)
      case (1)
        call freemem(this%node_map_1d)
        call alloc(this%node_map_1d, g_nx)
        this%node_map_1d(:) = 0
      case (2)
        call freemem(this%node_map_2d)
        call alloc(this%node_map_2d, g_nx)
        this%node_map_2d(:,:) = 0
      case (3)
        call freemem(this%node_map_3d)
        call alloc(this%node_map_3d, g_nx)
        this%node_map_3d(:,:,:) = 0
      case default
        write(0,*) 'Error initializing node map'
        call abort_program()
      end select

      ! initialize particle per node & particle count arrays
      do rank = 2, no_num( no_co )

        call mpi_recv( ig_xbnd_inj_temp, 2 * p_x_dim, MPI_INTEGER, rank-1, 0, &
                       comm(no_co), stat, ierr )

        select case (p_x_dim)
        case (1)
          do i = ig_xbnd_inj_temp(p_lower,1), ig_xbnd_inj_temp(p_upper,1)
            this%node_map_1d(i) = rank-1
          enddo
        case (2)
          do j = ig_xbnd_inj_temp(p_lower,2), ig_xbnd_inj_temp(p_upper,2)
            do i = ig_xbnd_inj_temp(p_lower,1), ig_xbnd_inj_temp(p_upper,1)
              this%node_map_2d(i,j) = rank-1
            enddo
          enddo
        case (3)
          do k = ig_xbnd_inj_temp(p_lower,3), ig_xbnd_inj_temp(p_upper,3)
            do j = ig_xbnd_inj_temp(p_lower,2), ig_xbnd_inj_temp(p_upper,2)
              do i = ig_xbnd_inj_temp(p_lower,1), ig_xbnd_inj_temp(p_upper,1)
                this%node_map_3d(i,j,k) = rank-1
              enddo
            enddo
          enddo
        case default
          write(0,*) 'Error initializing node map'
          call abort_program()
        end select

      enddo

    else

      call mpi_ssend( ig_xbnd_inj(p_lower:p_upper, 1:p_x_dim), 2 * p_x_dim, MPI_INTEGER, &
                      0, 0, comm(no_co), ierr )

    endif


    if( mpi_node() == 0 ) then
      call cpu_time(t1)
      print *," - Initializing file read of ", trim(this%file_name), " ..."
      ! open hdf5 file
      call open_hdf5( )
      call h5fopen_f (this%file_name, H5F_ACC_RDONLY_F, file_id, ierr)

      ! get dataspace parameters
      call h5dopen_f(file_id, 'x1', dset_id, ierr)

      call h5dget_space_f(dset_id,dspace,ierr)

      ! get dims of RAW file
      call h5sget_simple_extent_dims_f( dspace, dims, maxdims, ierr )

      ! close dataset
      call h5dclose_f(dset_id, ierr)

      ! setup particle chunk reading
      total_part_num = dims(1)
      niter = total_part_num/part_chunk_size

      ! determinte num of mpisends/recvs
      if(mod(total_part_num,part_chunk_size) /= 0) then
        niter = niter + 1
      endif

      do rank = 2, no_num( no_co )
        call mpi_send( niter, 1, MPI_INTEGER, rank-1, 0, &
                       comm(no_co), ierr )
      enddo

      print *, " - Reading in ", trim(tostring_int(total_part_num)), " particles ..."
      ! allocate buffer & read data (all x,p,q quants)
      if( total_part_num > part_chunk_size) then
        allocate(chunk_data(nquants,part_chunk_size))
      else
        allocate(chunk_data(nquants,total_part_num))
      endif

      do j =1, total_part_num , part_chunk_size

        ! check if last chunk of particles
        if(j + part_chunk_size > total_part_num) then
          npar = total_part_num - j + 1
        else
          npar = part_chunk_size
        endif
        ! read part chunk npar
        num_read = read_raw_chunk(this, species, file_id, j-1, j-1 + npar,chunk_data)
        ! post process chunk_data
        do i = 1, npar
          chunk_data(1,i) = chunk_data(1,i) - this%x1_offset
          chunk_data(nquants,i) = chunk_data(nquants,i) * this%q_mult
        enddo
        ! map and send parts to nodes
        ! returns chunk_inj to node 0
        ! node_data contains chunk_inj parts for rank 0
        chunk_inj = send_particle_chunk(this, node_data, chunk_data, no_co, g_xbnd, &
                                        g_nx, dx, npar, nquants)
        ! inject chunk
        if(chunk_inj > 0) then
          call inject_particle_chunk(this, species, node_data, chunk_inj, g_xbnd, dx)
        endif

      enddo


      ! deallocate node map/ chunk data
      deallocate(chunk_data)
      select case (p_x_dim)
      case (1)
        call freemem(this%node_map_1d)
      case (2)
        call freemem(this%node_map_2d)
      case (3)
        call freemem(this%node_map_3d)
      case default
        write(0,*) 'Error initializing node map'
        call abort_program()
      end select
      !close file
      CALL h5fclose_f(file_id, ierr)
      call cpu_time(t2)
      print *," - Raw file initialization completed in",t2-t1, "seconds"

    else
      ! mpi ranks >= 1 receiving data

      ! number of chunks to receive
      call mpi_recv( niter, 1, MPI_INTEGER, 0, 0, &
                     comm(no_co), stat, ierr )

      ! recieve niter chunks
      do i = 1, niter

        ! recv num of particles
        call mpi_recv( ping, 1, MPI_INTEGER, 0, 0, &
                       comm(no_co), stat, ierr )
        chunk_inj = ping
        ! recv and inject chunk if nonzero # of particles
        if(chunk_inj > 0) then
          ! allocate buffer
          allocate(node_data(nquants,chunk_inj))
          call mpi_recv( node_data, chunk_inj*nquants, mpi_real_type(p_k_part), 0, 0, &
                         comm(no_co), stat, ierr )
          call inject_particle_chunk(this, species, node_data, chunk_inj, g_xbnd, dx)

          ! free data
          deallocate(node_data)
        endif

      enddo

    endif
    ! injection number update handled by inject_particle_chunk --  return 0
    num_inj = 0
    this%raw_read_done = .true.
  endif



end function inject_file

!-----------------------------------------------------------------------------------------
! read chunk from RAW hdf5
!-----------------------------------------------------------------------------------------
function read_raw_chunk( this, species, file_id, start, last, &
                              data_out) result(num_read)

  use hdf5_util
  use m_parameters
  use m_system
  use m_diagfile_hdf5
  use m_diagfile
  use hdf5
  use m_diagnostic_utilities
  use stringutil
  implicit none

  class( t_psource_file ), intent( inout ) :: this
  class(t_species), intent(inout), target :: species
  integer(hid_t), intent(inout) :: file_id       ! File identifier
  integer, intent(in) :: start, last
  real(p_k_part), dimension(:,:), pointer, intent(inout) :: data_out
  integer :: num_read
  real :: t1, t2, bytes_read



  integer(hid_t) :: dset_id       ! Dataset identifier
  integer(hid_t) :: datatypeID, dspace, rootID , memspaceID
  real(p_double), dimension(:,:), pointer :: data_dp
  real(p_single), dimension(:,:), pointer :: data_sp

  integer(hsize_t), dimension(1) :: dims,maxdims
  integer(hsize_t), dimension(1) :: count, offset, offset_out, count_out
  character(len = p_str_maxlen ), dimension(:), allocatable :: datasets
  integer :: ierr, nquants, i, j, k
  character(len=24) :: data, bandwidth, num_parts

  data_dp => null(); data_sp => null()

  ! set nquants to read in xdims + pdims + charge q
  nquants = species % get_n_x_dims() + p_p_dim + 1
  ! set offset in file to start reading
  offset(1) = start
  count(1) = last-start
  offset_out(1) = 0
  count_out(1) = last-start
  ! allocate buffers
  allocate(datasets(nquants))
  allocate(data_dp(count(1),nquants))
  allocate(data_sp(count(1),nquants))
  data_dp = 0
  data_sp = 0
  k = 1
  do i = 1, species % get_n_x_dims()
    datasets(k) = 'x'//(char(iachar('0')+i))
    k = k + 1
  enddo
  do i = 1, p_p_dim
    datasets(k) = 'p'//(char(iachar('0')+i))
    k = k + 1
  enddo
  datasets(k) = 'q'
  call cpu_time(t1)

  do k = 1, nquants

    call h5dopen_f(file_id, datasets(k), dset_id, ierr)

    call h5dget_type_f(dset_id, datatypeID, ierr )

    call h5dget_space_f(dset_id,dspace,ierr)

    call h5sget_simple_extent_dims_f( dspace, dims, maxdims, ierr )

    if(start > maxdims(1)) then
      num_read = 0
      call h5dclose_f(dset_id, ierr)
      return
    else
      count(1) = min(maxdims(1) - start, count(1))
      num_read = count(1)
      num_parts = tostring_int(num_read)
    endif

    ! select slab in dataspace
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, &
                               offset, count, ierr)

    ! create memory dspace
    call h5screate_simple_f(1,count, memspaceID, ierr)

    ! select hyperslab in memory
    call h5sselect_hyperslab_f(memspaceID, H5S_SELECT_SET_F, &
                               offset_out, count_out, ierr)

    if(hdf5_diag_type(datatypeID) == p_diag_float64) then
      call h5dread_f(dset_id, datatypeID, data_dp(:,k),count, ierr, &
                     file_space_id=dspace, mem_space_id = memspaceID )
      do j = 1, count(1)
        data_out(k,j) = real(data_dp(j,k), p_k_part)
      enddo

    else
      call h5dread_f(dset_id, datatypeID , data_sp(:,k),count, ierr, &
                     file_space_id=dspace, mem_space_id = memspaceID)
      do j = 1, count(1)
        data_out(k,j) = real(data_sp(j,k), p_k_part)
      enddo


    endif
    call h5sclose_f(dspace, ierr)
    call h5sclose_f(memspaceID, ierr)
    call h5dclose_f(dset_id, ierr)
  enddo

  call cpu_time(t2)
  bytes_read = num_read *nquants
  if(hdf5_diag_type(datatypeID) == p_diag_float64) then
    bytes_read = bytes_read * 8
  else
    bytes_read = bytes_read * 4
  endif
  data = tostring_r4(real(bytes_read/1024.0**2, p_single))
  bandwidth = tostring_r4(real(bytes_read/1024.0**2/(t2-t1), p_single))
  !print *, " - Read in ", trim(data)," MB @ ", trim(bandwidth), " MB/s"
  ! clean up
  deallocate(data_dp)
  deallocate(data_sp)
  deallocate(datasets)



end function read_raw_chunk

#else
!-----------------------------------------------------------------------------------------
! RAW Initialization ZDF routines
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Inject particles using information from external ZDF file
!-----------------------------------------------------------------------------------------
function inject_file( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, &
                      send_msg, recv_msg ) result(num_inj)

  use m_parameters
  use m_system
  use m_diagfile
  use m_diagnostic_utilities
  use stringutil
  implicit none

  class(t_psource_file), intent(inout) :: this
  class(t_species), intent(inout), target :: species
  integer, dimension(:, :), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout)   :: jay
  class( t_node_conf ), intent(in)     :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: num_inj 

  write(0,*) 'RAW ZDF particle initialization not supported yet. Please use HDF5'
  call abort_program()

end function inject_file


!-----------------------------------------------------------------------------------------
! read chunk from RAW ZDF file
!-----------------------------------------------------------------------------------------
function read_raw_chunk( this, species, file_id, start, last, &
                              data_out) result(num_read)

  use m_parameters
  use m_system
  use m_diagfile
  use m_diagnostic_utilities
  use stringutil
  implicit none

  class( t_psource_file ), intent( inout ) :: this
  class(t_species), intent(inout), target :: species
  character(len=1024), intent(inout) :: file_id       ! File identifier
  integer, intent(in) :: start, last
  real(p_k_part), dimension(:,:), pointer, intent(inout) :: data_out

  integer :: num_read 


  write(0,*) 'RAW ZDF particle initialization not supported yet. Please use HDF5'
  call abort_program()


end function read_raw_chunk



#endif

!-----------------------------------------------------------------------------------------
! send particles
!-----------------------------------------------------------------------------------------
function send_particle_chunk( this, node_data, chunk_data, no_co, g_xbnd, g_nx, dx, npar,&
                              nquants) result(num_inj)

  use m_parameters
  use m_system
  use m_diagfile
  use m_diagnostic_utilities
  use stringutil
  implicit none

  class(t_psource_file), intent(inout) :: this
  real(p_k_part), dimension(:,:), pointer, intent(in) :: chunk_data
  real(p_k_part), dimension(:,:), pointer, intent(inout) :: node_data
  class( t_node_conf ), intent(in)     :: no_co
  real(p_double), dimension(p_x_dim), intent(in) :: dx
  real(p_double), dimension(2, p_x_dim), intent(in) :: g_xbnd
  integer, dimension(p_x_dim), intent(in) :: g_nx
  integer, intent(in) :: npar, nquants

  integer, dimension(:), pointer :: displacements, blocks, part_per_node
  integer, dimension(:), pointer :: node_starts, node_counters, particle_tags
  integer :: handle, rank, node_number, tag_index
  integer :: i, j, k, ltype, ierr, num_inj
  integer, dimension(p_x_dim) :: idxpos
  real(p_k_part), dimension(p_x_dim) :: xpos
  integer, dimension(2) :: dims

  displacements => null(); blocks => null(); part_per_node => null()
  node_starts => null(); node_counters => null(); particle_tags => null()

  call freemem( part_per_node)
  call alloc(part_per_node, (/ no_num(no_co) /))
  part_per_node(:) = 0
  ! first pass; calculate # parts to send to each node
  j=0
  do i = 1, npar
    xpos = chunk_data(1:p_x_dim,i)
    ! calculate cell idx
    do k = 1, p_x_dim
      idxpos(k) = int((xpos(k)-g_xbnd(p_lower,k))/dx(k)) + 1
    enddo

    ! check if inside grid bounds
    if(any(idxpos < 1) .or. any(idxpos > g_nx) ) then
      cycle
    endif
    j = j + 1


    select case (p_x_dim)
    case (1)
      node_number = this%node_map_1d(idxpos(1))
      part_per_node(node_number+1) = part_per_node(node_number +1) + 1
    case (2)
      node_number = this%node_map_2d(idxpos(1),idxpos(2))
      part_per_node(node_number+1) = part_per_node(node_number +1) + 1
    case (3)
      node_number = this%node_map_3d(idxpos(1),idxpos(2),idxpos(3))
      part_per_node(node_number+1) = part_per_node(node_number +1) + 1
    case default
      write(0,*) 'Error setting particle counts for first pass'
      call abort_program()
    end select


  enddo
  ! second pass; generate part list for each node in order (node 0, node 1, etc..)
  ! allocate counter/start/tag arrays
  call freemem(node_counters)
  call alloc(node_counters, (/ no_num(no_co) /))
  node_counters = 0

  call freemem(node_starts)
  call alloc(node_starts, (/ no_num(no_co) /))
  node_starts = 0

  call freemem(particle_tags)
  call alloc(particle_tags, (/ j /))
  particle_tags = 0

  do i = 2, no_num(no_co)
    node_starts(i) = node_starts(i-1) + part_per_node(i-1)
  enddo


  do i = 1, npar
    xpos = chunk_data(1:p_x_dim,i)
    do k = 1, p_x_dim
      idxpos(k) = int((xpos(k)-g_xbnd(p_lower,k))/dx(k)) + 1
      ! check if inside grid bounds

    enddo
    if(any(idxpos < 1) .or. any(idxpos > g_nx) ) then
      cycle
    endif


    select case (p_x_dim)
    case (1)
      node_number = this%node_map_1d(idxpos(1))
      tag_index = node_counters(node_number+1) + node_starts(node_number+1)
      particle_tags(tag_index + 1) = i
      node_counters(node_number+1) = node_counters(node_number+1) + 1
    case (2)
      node_number = this%node_map_2d(idxpos(1),idxpos(2))
      tag_index = node_counters(node_number+1) + node_starts(node_number+1)
      particle_tags(tag_index + 1) = i
      node_counters(node_number+1) = node_counters(node_number+1) + 1
    case (3)
      node_number = this%node_map_3d(idxpos(1),idxpos(2),idxpos(3))
      tag_index = node_counters(node_number+1) + node_starts(node_number+1)
      particle_tags(tag_index + 1) = i
      node_counters(node_number+1) = node_counters(node_number+1) + 1
    case default
      write(0,*) 'Error setting particle counts for first pass'
      call abort_program()
    end select


  enddo
  ! disperse particles to each node
  do rank = 2, no_num( no_co )

    call mpi_send( part_per_node(rank), 1, MPI_INTEGER, rank-1, 0, &
                   comm(no_co), ierr )

    if(part_per_node(rank) > 0) then
      allocate(blocks(part_per_node(rank)))
      allocate(displacements(part_per_node(rank)))
      do i = 1, part_per_node(rank)
        tag_index = i + node_starts(rank)
        displacements(i) = nquants * (particle_tags(tag_index)-1)
        blocks(i) = nquants
      enddo
      call mpi_type_indexed( part_per_node(rank), blocks, displacements, &
                             mpi_real_type(p_k_part), ltype, ierr)
      call mpi_type_commit(ltype,ierr)
      call mpi_send( chunk_data, 1, ltype, rank-1, 0, &
                     comm(no_co), ierr )

      deallocate(blocks)
      deallocate(displacements)
    endif
  enddo

  ! set root node data
  num_inj = part_per_node(1)
  if(part_per_node(1) > 0) then
    allocate(node_data(nquants,part_per_node(1)))

    ! copy chunk data to node data
    do i = 1, part_per_node(1)
      tag_index = i + node_starts(1)
      do j = 1, nquants
        node_data(j,i) = chunk_data(j,particle_tags(tag_index))
      enddo
    enddo
  endif

  call freemem(part_per_node)
  call freemem(node_counters)
  call freemem(node_starts)
  call freemem(particle_tags)



end function send_particle_chunk


!-----------------------------------------------------------------------------------------
! inject particles
!-----------------------------------------------------------------------------------------
subroutine inject_particle_chunk(this, species, node_data, num_inj, g_xbnd, dx)

  use m_parameters
  use m_system
  use m_diagfile
  use m_diagnostic_utilities
  use stringutil
  use m_species_tag, only : set_tags
  implicit none

  class(t_psource_file), intent( inout ) :: this
  class(t_species), intent(inout), target :: species
  real(p_k_part), dimension(:,:), pointer, intent(in) :: node_data
  real(p_double), dimension(p_x_dim), intent(in) :: dx
  real(p_double), dimension(2, p_x_dim), intent(in) :: g_xbnd
  integer, intent(in) :: num_inj

  integer :: i, j, k


  if ( num_inj > species%num_par_max - species%num_par ) then
    call species % grow_buffer( species%num_par_max + num_inj + p_spec_buf_block )
  endif
  do i = 1, num_inj

    ! set grid cell locations
    do j = 1, p_x_dim
      species%ix(j, species%num_par + 1) = int((node_data(j,i)-g_xbnd(p_lower,j))/dx(j)) &
                                           + 1 - species%my_nx_p( p_lower, j ) + 1
    enddo

    ! set displacements
    k = 1
    do j = 1, p_x_dim
      species%x(j, species%num_par + 1) = (node_data(j,i)-species%g_box(p_lower,j))/dx(j)&
                                          + 1 - species%ix(j, species%num_par + 1) &
                                          - species%my_nx_p( p_lower, j ) + 1
      k = k + 1
    enddo

    ! copy residual data, i.e. from quasi-3D
    do j = p_x_dim + 1, species% get_n_x_dims()
      species%x(j, species%num_par + 1) = node_data(j,i)
      k = k + 1
    enddo

    ! copy momentum
    do j = 1, p_p_dim
      species%p(j, species%num_par + 1) = node_data(k,i)
      k = k + 1
    enddo
    species%q(species%num_par + 1) = node_data(k, i)
    

    ! if required set tags of particles
    if(species%add_tag) then
      call set_tags( species, species%num_par + 1)
    endif

    ! update species counters
    species%num_par = species%num_par + 1
    species%num_created = species%num_created + 1
  enddo



end subroutine inject_particle_chunk




end module m_psource_file
