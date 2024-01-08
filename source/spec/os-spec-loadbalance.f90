! Uncomment below to enable debug code
!#define __DEBUG__ 1

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_loadbalance

#include "memory/memory.h"

use m_species_define
use m_species_comm

use m_node_conf

use m_grid_define
use m_parameters

implicit none

private

! Amount of space to grow idx lists when necessary
integer, parameter :: p_list_idx_block = 65536

interface deposit_cell_load
  module procedure deposit_cell_load
end interface

interface add_density_load
  module procedure add_density_load
end interface

interface add_particle_load
  module procedure add_particle_load
end interface

interface init_send_idx
  module procedure init_send_idx
end interface

interface pack_particles
  module procedure pack_particles
end interface

public :: add_particle_load, add_density_load, deposit_cell_load, init_send_idx
public :: pack_particles

contains


!---------------------------------------------------------------------------------------------------
! Deposit the number of particles on each cell
!---------------------------------------------------------------------------------------------------
subroutine deposit_cell_load( this, vdf )

  use m_vdf_define

  implicit none

  class(t_species), intent(in) :: this
  type(t_vdf), intent(inout)  :: vdf

  integer :: i

  select case ( p_x_dim )

    case(1)
      do i = 1, this%num_par
        vdf%f1(1, this%ix(1,i) ) = vdf%f1(1, this%ix(1,i) ) + 1
      enddo
    case(2)
      do i = 1, this%num_par
        vdf%f2(1, this%ix(1,i), this%ix(2,i)) = vdf%f2(1, this%ix(1,i), this%ix(2,i)) + 1
      enddo

    case(3)
      do i = 1, this%num_par
        vdf%f3(1, this%ix(1,i), this%ix(2,i), this%ix(3,i)) = &
                                             vdf%f3(1, this%ix(1,i), this%ix(2,i), this%ix(3,i)) + 1
      enddo

  end select


end subroutine deposit_cell_load
!---------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Determines the number of particles per cell along all directions based on
! the density profile for the species
!-------------------------------------------------------------------------------
subroutine add_density_load( this, grid )
!-------------------------------------------------------------------------------

  use m_psource_std

  implicit none

  ! this must be inout because of the den_value calls
  class( t_species ), intent(inout) :: this
  class( t_grid ), intent(inout) :: grid


  real( p_double ), dimension(:), pointer :: int_load

  integer, dimension(p_max_dim) :: g_nx
  integer, dimension(2,p_max_dim) :: my_nx
  integer :: npar_cell

  integer :: i

  real(p_k_part), dimension(p_max_dim) :: dx
  real(p_k_part), dimension(p_max_dim ) :: g_xmin
  real(p_k_part) :: x1, x2

  integer, dimension( p_max_dim ) :: dim_off

  integer :: i1, i2, i3

  ! particle positions (global / inside cell )
  real(p_k_part), dimension(:,:), pointer :: ppos, ppos_cell

  ! particle charges
  real(p_k_part), dimension(:), pointer :: pcharge

  ! increment in position for the particles
  real(p_k_part), dimension(3) :: dxpart

  real( p_double ) :: cshift

  class(t_psource), pointer :: profile

  profile => this % source

  ! This will only work if using a t_profile_std (or subclass) profile
  select type (profile)

  class is (t_psource_std)

  ! get pointer to int_load array
  int_load => grid%int_load

  ! Get dimensional offset for int_load array
  dim_off(1) = 0
  do i = 2, grid%x_dim
    dim_off(i) = dim_off(i-1) + grid%g_nx(i-1)
  enddo

  ! get grid sizes and distance between particles
  do i = 1, p_x_dim
    g_nx(i) = grid%g_nx(i)
    my_nx(p_lower:p_upper, i) = grid%my_nx(p_lower:p_upper, i)
    g_xmin( i ) = grid%g_box( p_lower, i )

    ! the species%dx variable has not been defined yet
    dx(i) = ( grid%g_box( p_upper, i ) - grid%g_box( p_lower, i ) ) / g_nx(i)

    dxpart(i) = 1.0_p_k_part / this%num_par_x(i)
  enddo

  ! get total number of particles per cell
  npar_cell = this%num_par_x(1)
  do i = 2, p_x_dim
    npar_cell = npar_cell*this%num_par_x(i)
  enddo

  ! initialize temp buffers
  call alloc( ppos, (/ p_x_dim, npar_cell /) )
  call alloc( ppos_cell, (/ p_x_dim, npar_cell /) )
  call alloc( pcharge, (/ npar_cell /))


  if ( this%pos_type == p_cell_near ) then

    ! positions inside cell are from -0.5 to +0.5
    cshift = 0.5_p_double

    if ( this%coordinates /= p_cylindrical_b ) then

      ! The particle global position is considered to be calculated from the center of the cell, so
      ! that a particle with gix = 1, x = -0.5 will be at the left edge of the simulation

      ! To simplify absolute position calculations we can just shift the global min by half a cell
      do i = 1, p_x_dim
        g_xmin(i) = g_xmin(i) + 0.5_p_double * this%dx(i)
      enddo

      ! all quantities depending on global positions must (should, it's just half a cell) be
      ! adjusted

    else
      ! In cylindrical coordinates we cannot do this for the radial direction, because the algorithm
      ! puts the radial boundary at position r = 0, which is placed at the center of cell 1.

      ! In this case, for even interpolation, particles with gix = 2, x = -0.5 will be in the axis,
      ! and for odd interpolation, this would be gix = 1, x = +0.5

      g_xmin(1) = g_xmin(1) + 0.5_p_double * this%dx(1)
    endif

  else

    cshift = 0.0_p_double
  endif

  ! count number of particles to inject in each direction
  select case( p_x_dim )

  case (1)
    i = 0
    do i1 = 1, this%num_par_x(1)
      i = i + 1
      ppos_cell(1,i) = real( dxpart(1) * ( i1 - 0.5_p_double ) - cshift, p_k_part )
    enddo

    do i1 = my_nx(p_lower,1), my_nx(p_upper,1)
      ppos(1,:) = (ppos_cell(1,:) + (i1-1))*dx(1) + g_xmin(1)
      call profile % get_den_value( ppos, npar_cell, pcharge )

      do i = 1, npar_cell
        if (pcharge(i) >= profile%den_min) then
          int_load(i1) = int_load(i1) + 1
        endif
      enddo
    enddo

  case (2)
    i = 0
    do i1 = 1, this%num_par_x(1)
      x1 = real( dxpart(1) * ( i1 - 0.5_p_double ) - cshift, p_k_part )
      do i2 = 1, this%num_par_x(2)
        i = i + 1
        ppos_cell(1,i) = x1
        ppos_cell(2,i) = real( dxpart(2) * ( i2 - 0.5_p_double ) - cshift, p_k_part )
      enddo
    enddo

    do i1 = my_nx(p_lower,1), my_nx(p_upper,1)
      ppos(1,:) = (ppos_cell(1,:) + (i1-1))*dx(1) + g_xmin(1)

      do i2 = my_nx(p_lower,2), my_nx(p_upper,2)
        ppos(2,:) = (ppos_cell(2,:) + (i2-1))*dx(2) + g_xmin(2)

        call profile % get_den_value( ppos, npar_cell, pcharge )
        do i = 1, npar_cell
          if ( pcharge(i) >= profile%den_min) then
            int_load(i1+dim_off(1)) = int_load(i1+dim_off(1)) + 1
            int_load(i2+dim_off(2)) = int_load(i2+dim_off(2)) + 1
          endif
        enddo

      enddo

    enddo

  case (3)
    i = 0
    do i1 = 1, this%num_par_x(1)
      x1 = real( dxpart(1) * ( i1 - 0.5_p_double) - cshift, p_k_part )
      do i2 = 1, this%num_par_x(2)
        x2 = real( dxpart(2) * ( i2 - 0.5_p_double ) - cshift, p_k_part )
        do i3 = 1, this%num_par_x(3)
          i = i + 1
          ppos_cell(1,i) = x1
          ppos_cell(2,i) = x2
          ppos_cell(3,i) = real( dxpart(3) * ( i3 - 0.5_p_double ) - cshift, p_k_part )
        enddo
      enddo
    enddo


    do i1 = my_nx(p_lower,1), my_nx(p_upper,1)
      ppos(1,:) = (ppos_cell(1,:) + (i1-1))*dx(1) + g_xmin(1)

      do i2 = my_nx(p_lower,2), my_nx(p_upper,2)
        ppos(2,:) = (ppos_cell(2,:) + (i2-1))*dx(2) + g_xmin(2)

        do i3 = my_nx(p_lower,3), my_nx(p_upper,3)
          ppos(3,:) = (ppos_cell(3,:) + (i3-1))*dx(3) + g_xmin(3)

          call profile % get_den_value( ppos, npar_cell, pcharge )
          do i = 1, npar_cell
            if (pcharge(i) >= profile%den_min) then
              int_load(i1+dim_off(1)) = int_load(i1+dim_off(1)) + 1
              int_load(i2+dim_off(2)) = int_load(i2+dim_off(2)) + 1
              int_load(i3+dim_off(3)) = int_load(i3+dim_off(3)) + 1
            endif
          enddo

        enddo

      enddo

    enddo

  end select

  call freemem( ppos )
  call freemem( ppos_cell )
  call freemem( pcharge )

  class default

  if ( mpi_node() == 0 ) then
    write(0,*) "(*warning*) add_density_load is not supported by the selected &
    &species initialization type"
  endif

  end select

end subroutine add_density_load
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Determines the number of particles per cell along all directions based on
! the particle positions for the species
!-------------------------------------------------------------------------------
subroutine add_particle_load( this, grid )

   implicit none

   class( t_species ), intent(in) :: this
   class( t_grid ), intent(inout) :: grid

   integer, dimension( p_max_dim ) :: dim_off
   real( p_double ), dimension(:), pointer :: int_load
   integer :: i, i1, i2, i3

   ! Get the offset of the int_load array for each direction
   dim_off(1) = 0
   do i = 2, grid%x_dim
     dim_off(i) = dim_off(i-1) + grid%g_nx(i-1)
   enddo

   ! Add the correction for global cell positions to dim_off
   do i = 1, grid%x_dim
     dim_off(i) = dim_off(i) + ( this%my_nx_p( p_lower, i ) - 1 )
   enddo

   ! just for clarity
   int_load => grid%int_load

   ! Doing it like this allows us to unroll the inner loop over the number of dimensions
   select case ( grid%x_dim )
     case(1)
	   do i=1, this%num_par
		 i1 = dim_off(1) + this%ix(1,i)

		 int_load(i1) = int_load(i1) + 1
	   enddo

     case(2)
	   do i=1, this%num_par
		 i1 = dim_off(1) + this%ix(1,i)
		 i2 = dim_off(2) + this%ix(2,i)

		 int_load(i1) = int_load(i1) + 1
		 int_load(i2) = int_load(i2) + 1
	   enddo

     case(3)
	   do i=1, this%num_par
		 i1 = dim_off(1) + this%ix(1,i)
		 i2 = dim_off(2) + this%ix(2,i)
		 i3 = dim_off(3) + this%ix(3,i)

		 int_load(i1) = int_load(i1) + 1
		 int_load(i2) = int_load(i2) + 1
		 int_load(i3) = int_load(i3) + 1
	   enddo

   end select

end subroutine add_particle_load
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine init_send_idx( this, msg_patt, new_grid, send_idx, local_idx   )
!-------------------------------------------------------------------------------
! Find indexes of particles to be sent to other nodes / correct particle positions for
! particles staying in current node.
!-------------------------------------------------------------------------------
  implicit none

  ! dummy variables
  class (t_species), intent(inout) :: this
  type (t_msg_patt), intent(in) :: msg_patt
  class (t_grid), intent(in) :: new_grid

  type(t_part_idx), dimension(:), intent(inout) :: send_idx
  type(t_part_idx), intent(inout) :: local_idx


  ! local variables

  ! New coordinates of local grid relative to current grid
  integer, dimension(2,p_max_dim) :: new_nx_p

  ! integer shift for particles staying in the same node
  integer, dimension(3) :: lshift

  integer :: i, j

  ! position of new local grid coordinates relative to old local grid coordinates, and
  ! cell shift for particles staying in current node
  do i = 1, p_x_dim
    new_nx_p( 1, i ) = new_grid % my_nx( 1, i ) - this%my_nx_p( p_lower, i ) + 1
    new_nx_p( 2, i ) = new_grid % my_nx( 2, i ) - this%my_nx_p( p_lower, i ) + 1

	lshift(i)        = this%my_nx_p( p_lower, i ) - new_grid % my_nx( p_lower, i )
  enddo

  ! Initialize index lists (this should not be necessary since these are the default
  ! values for the object)
  do j = 1, msg_patt%n_send
    send_idx(j)%start   = 1
    send_idx(j)%nidx    = 0
    send_idx(j)%buf_size = 0
  enddo

  local_idx%start = 1
  local_idx%nidx = 0
  local_idx%buf_size = 0

  ! Process particles
  select case (p_x_dim)
	case(1)
	   do i = 1, this%num_par
		 if ((this%ix(1,i) < new_nx_p(1,1)) .or. &  ! exit x1 lower
			 (this%ix(1,i) > new_nx_p(2,1))) then   ! exit x1 uppper

			! loop through all nodes and see where the particle needs to
			! be sent
			do j = 1, msg_patt%n_send
			   if ((this%ix(1,i) >=  msg_patt%send(j)%cells(1,1)) .and. &
				   (this%ix(1,i) <=  msg_patt%send(j)%cells(2,1)) ) then

                  call add_particle_index( send_idx(j),  i )
				  exit

			   endif
			enddo

		 else
			! correct indexes of particle
			this%ix(1,i) = this%ix(1,i) + lshift(1)
            call add_particle_index( local_idx, i )
		 endif
	   enddo

	case(2)

	   do i = 1, this%num_par

		 if ((this%ix(1,i) < new_nx_p(1,1)) .or. &  ! exit x1 lower
			 (this%ix(1,i) > new_nx_p(2,1)) .or. &  ! exit x1 uppper
			 (this%ix(2,i) < new_nx_p(1,2)) .or. &  ! exit x2 lower
			 (this%ix(2,i) > new_nx_p(2,2)) ) then  ! exit x2 uppper


			! loop through all nodes and see where the particle needs to
			! be sent
			do j = 1, msg_patt%n_send
			   if ((this%ix(1,i) >=  msg_patt%send(j)%cells(1,1)) .and. &
				   (this%ix(1,i) <=  msg_patt%send(j)%cells(2,1)) .and. &
				   (this%ix(2,i) >=  msg_patt%send(j)%cells(1,2)) .and. &
				   (this%ix(2,i) <=  msg_patt%send(j)%cells(2,2)) ) then

                  call add_particle_index( send_idx(j),  i )
				  exit

			   endif
			enddo

		 else

		   ! correct indexes of particle
		   this%ix(1,i) = this%ix(1,i) + lshift(1)
		   this%ix(2,i) = this%ix(2,i) + lshift(2)
           call add_particle_index( local_idx, i )
		 endif

	   enddo

	case(3)
	   do i = 1, this%num_par
		 if ((this%ix(1,i) < new_nx_p(1,1)) .or. &  ! exit x1 lower
			 (this%ix(1,i) > new_nx_p(2,1)) .or. &  ! exit x1 uppper
			 (this%ix(2,i) < new_nx_p(1,2)) .or. &  ! exit x2 lower
			 (this%ix(2,i) > new_nx_p(2,2)) .or. &  ! exit x2 uppper
			 (this%ix(3,i) < new_nx_p(1,3)) .or. &  ! exit x3 lower
			 (this%ix(3,i) > new_nx_p(2,3)) ) then  ! exit x3 uppper

			! loop through all nodes and see where the particle needs to
			! be sent

			do j = 1, msg_patt%n_send
			   if ((this%ix(1,i) >=  msg_patt%send(j)%cells(1,1)) .and. &
				   (this%ix(1,i) <=  msg_patt%send(j)%cells(2,1)) .and. &
				   (this%ix(2,i) >=  msg_patt%send(j)%cells(1,2)) .and. &
				   (this%ix(2,i) <=  msg_patt%send(j)%cells(2,2)) .and. &
				   (this%ix(3,i) >=  msg_patt%send(j)%cells(1,3)) .and. &
				   (this%ix(3,i) <=  msg_patt%send(j)%cells(2,3)) ) then

                  call add_particle_index( send_idx(j),  i )

				  exit

			   endif
			enddo

		 else
			! correct indexes of particle
			this%ix(1,i) = this%ix(1,i) + lshift(1)
			this%ix(2,i) = this%ix(2,i) + lshift(2)
			this%ix(3,i) = this%ix(3,i) + lshift(3)
            call add_particle_index( local_idx, i )

		 endif
	   enddo

  end select


  contains

  ! ----

  subroutine add_particle_index( list, i )

     implicit none

     type(t_part_idx), intent(inout) :: list
     integer, intent(in)              :: i

     integer, dimension(:), pointer :: tmp_idx

     ! grow buffer if necessary
     if ( list%nidx + 1 > list%buf_size ) then
        list%buf_size = list%buf_size + p_spec_buf_block
        call alloc(  tmp_idx, (/ list%buf_size /))
        if ( list%nidx > 0 ) call memcpy( tmp_idx, list%idx, list%nidx )
        call freemem( list%idx )
        list%idx => tmp_idx
     endif

     ! add particle index to send index list
     list%nidx = list%nidx + 1
     list%idx( list%nidx ) = i

  end subroutine add_particle_index

end subroutine init_send_idx
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
! Pack all particles in local_idx at the beggining of the particle buffer. Particles not
! in local_idx are destroyed, and the species%num_par is corrected
!-----------------------------------------------------------------------------------------
subroutine pack_particles( species, local_idx )

  implicit none

  class( t_species ), intent(inout) :: species
  type(t_part_idx), intent(inout) :: local_idx

  integer :: i, j, k, l, n_x_dim

  k = 0
  n_x_dim = species%get_n_x_dims()

  do i = local_idx%start, local_idx%nidx
    j = local_idx%idx(i)

    k = k+1

    species%x (1:n_x_dim, k) = species%x (1:n_x_dim, j)
    species%ix(1:p_x_dim, k) = species%ix(1:p_x_dim, j)

    species%q( k  ) = species%q(   j)

    species%p(1:p_p_dim, k) = species%p(1:p_p_dim, j)

#ifdef __HAS_SPIN__
    species%s(1:p_s_dim, k) = species%s(1:p_s_dim, j)
#endif

    ! remove tracking data if necessary
    if ( species%add_tag ) then
      species%tag(1, k) = species%tag(1, j)
      species%tag(2, k) = species%tag(2, j)

#ifdef __HAS_TRACKS__

      ! change tracked particle indexes if needed
      if ( species%diag%ndump_fac_tracks > 0 .and. species%tag(1,k) < 0 ) then
        do l = 1, species%diag%tracks%npresent
          if (species%diag%tracks%tracks( species%diag%tracks%present(l) )%part_idx == j) then
            species%diag%tracks%tracks( species%diag%tracks%present(l) )%part_idx = k
            exit
          endif
        enddo
      endif

#endif

    endif

  enddo

  species%num_par = k

#ifdef __DEBUG__
  ! (* debug *) Invalidate charge for rest of the buffer
  do i = species%num_par + 1, species%num_par_max
    species%q(i) = 0
  enddo
#endif

end subroutine pack_particles
!-----------------------------------------------------------------------------------------


end module

!-----------------------------------------------------------------------------------------
! redistribute the species particle data through all nodes when node grids change and
! store new grid boundaries
!-----------------------------------------------------------------------------------------
subroutine reshape_spec( this, old_grid, new_grid, msg_patt, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

  use m_species_define, only : t_species, t_part_idx, t_spec_msg
  use m_species_comm
  use m_node_conf, only : t_node_conf
  use m_grid_define, only : t_grid, t_msg_patt
  use m_parameters
  use m_species_memory
  use m_species_loadbalance
  use m_vdf_comm, only : t_vdf_msg, reshape_copy
  use m_vdf_define, only : t_vdf_report

  implicit none

  class( t_species ), intent(inout) :: this
  type( t_msg_patt ), intent(in) :: msg_patt
  class( t_grid ), intent(in) :: old_grid, new_grid
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  ! local variables

  ! send and receive buffers
  type (t_spec_msg), dimension(:), pointer :: send, recv

  ! index of particles to send
  type(t_part_idx) :: local_idx
  type(t_part_idx), dimension(:), pointer :: send_idx

  integer :: i,j

  integer :: particle_size, my_node
  integer, parameter, dimension(p_max_dim) :: shift = 0
  integer, parameter :: p_recv = 1, p_send = 2
  integer, parameter :: p_dyn_lb_tag = 1000
  type(t_vdf_report), pointer :: report

#ifdef __DEBUG__
  integer(p_int64) :: npar, totpar0, totpar1
  integer :: ierr

  npar = this%num_par
  call MPI_ALLREDUCE( npar, totpar0, 1, MPI_INTEGER8, MPI_SUM, comm(no_co), ierr)
#endif

  my_node = no_co%my_aid()
  particle_size = particle_pack_size( this )

  ! if any messages to receive
  if ( msg_patt%n_recv > 0 ) then

    ! allocate recv buffers
    call alloc( recv, (/msg_patt%n_recv/))

    ! post reveives for 1st message block

    do j = 1, msg_patt%n_recv

      if ( msg_patt%recv(j)%node /= my_node ) then

        recv(j)%comm = comm( no_co )
        recv(j)%node = msg_patt%recv(j)%node
        recv(j)%tag  = p_dyn_lb_tag
        recv(j)%particle_size  = particle_size

        call irecv_1( recv(j) )

      else

        recv(j)%request = MPI_REQUEST_NULL

      endif

    enddo

  endif

  ! if any messages to send
  if ( msg_patt%n_send > 0 ) then

    ! allocate send buffers and lists
    call alloc(send,     (/msg_patt%n_send/))
    call alloc(send_idx, (/msg_patt%n_send/))

    ! get index of particles to send
    call init_send_idx( this, msg_patt, new_grid, send_idx, local_idx )

    ! post sends for 1st message block
    do j = 1, msg_patt%n_send

      if ( msg_patt%send(j)%node /= my_node ) then
        send(j)%node = msg_patt%send(j)%node
        send(j)%comm = comm( no_co )
        send(j)%tag  = p_dyn_lb_tag
        send(j)%particle_size = particle_size

        call isend_1( send(j), this, send_idx(j)%idx, send_idx(j)%nidx, shift )

      else

        send(j)%request = MPI_REQUEST_NULL

      endif
    enddo

    ! cleanup send_idx lists
    do j = 1, msg_patt%n_send
      call freemem( send_idx(j)%idx )
    enddo
    call freemem(send_idx)

  endif


  ! Move local particles to the beginning of the particle buffer and delete all particles
  ! leaving the node
  call pack_particles( this, local_idx )

  ! Cleanup the local_idx buffer
  if ( local_idx%buf_size > 0 ) call freemem( local_idx%idx )
  local_idx%start = 1
  local_idx%nidx = 0
  local_idx%buf_size = 0

  ! store new local grid boundaries
  ! this must happen before unpacking received particles
  do j = 1, p_max_dim
    do i = 1, 3
      this%my_nx_p(i,j) = new_grid%my_nx(i,j)
    enddo
  enddo


#ifdef __DEBUG__
  ! (* debug *) Validate can only be called after storing the new local grid boundaries
  call this % validate( "reshape_spec : before 2nd communication stage" )
#endif

  ! Complete 2nd communication stage

! Wait for each stage 1 receive message to complete and post 2nd receive
  ! if necessary
  do j = 1, msg_patt%n_recv
    if ( msg_patt%recv(j)%node /= my_node ) call irecv_2( recv( j ) )
  enddo

  ! post sends for 2nd message block
  ! If message has a 2nd part, wait for 1st part to complete and send it
  ! otherwise just skip it
  do j = 1, msg_patt%n_send
    if ( msg_patt%send(j)%node /= my_node ) call isend_2( send( j ) )
  enddo

  ! Unpack received data (this waits for the 2nd stage message to complete if needed)
  do j = 1, msg_patt%n_recv
    if ( msg_patt%recv(j)%node /= my_node ) call unpack( recv(j), this, local_idx )
    call cleanup( recv(j) )
  enddo

  ! Cleanup message array
  if ( msg_patt%n_recv > 0 ) call freemem( recv )

  ! Wait for any remaining send messages and cleanup message buffers
  if ( msg_patt%n_send > 0 ) then
    call wait_spec_msg( send )
    call cleanup( send )
    call freemem( send )
  endif

  ! Reshape tavg_data
  report => this%diag%reports
  do
    if ( .not. associated(report) ) exit
    if ( report%tavg_data%x_dim_ > 0 ) then
      call reshape_copy( report%tavg_data, old_grid, new_grid, no_co, send_msg, recv_msg )
    endif
    report => report%next
  enddo

  report => this%diag%rep_cell_avg
  do
    if ( .not. associated(report) ) exit
    if ( report%tavg_data%x_dim_ > 0 ) then
      call reshape_copy( report%tavg_data, old_grid, new_grid, no_co, send_msg, recv_msg )
    endif
    report => report%next
  enddo

  report => this%diag%rep_udist
  do
    if ( .not. associated(report) ) exit
    if ( report%tavg_data%x_dim_ > 0 ) then
      call reshape_copy( report%tavg_data, old_grid, new_grid, no_co, send_msg, recv_msg )
    endif
    report => report%next
  enddo

#ifdef __DEBUG__
  ! (* debug *) Validate particles after reshape
  call this % validate( "reshape_spec : reshape complete" )

  ! (* debug *) Check total number of particles
  npar = this%num_par
  call MPI_ALLREDUCE( npar, totpar1, 1, MPI_INTEGER8, MPI_SUM, comm(no_co), ierr)
  if ( totpar0 /= totpar1 ) then
    write(0,'(A,I0,A)') '[', mpi_node(), '] (* error *) reshaping species failed:'
    write(0,*) 'Initial number of particles : ', totpar0
    write(0,*) 'Final number of particles   : ', totpar1
    call abort_program()
  endif
#endif

end subroutine reshape_spec
!-----------------------------------------------------------------------------------------
