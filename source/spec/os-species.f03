!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     particle species class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species

#include "memory/memory.h"

  use m_parameters
  use m_math

  use m_species_define
  use m_species_boundary
  use m_species_diagnostics
  use m_species_memory

  use m_species_charge
  use m_species_current

  use m_species_udist
#ifdef __HAS_SPIN__
  use m_species_sdist
#endif

  use m_species_radcool

#ifdef __HAS_TRACKS__
  use m_species_tracks
#endif

  use m_logprof
  use m_vdf_define

  use m_emf
  use m_emf_define
  use m_emf_interpolate

  use m_node_conf
  use m_grid_define
  use m_space

  use m_utilities
  use m_diagnostic_utilities

  use m_piston

  use m_parameters

  implicit none

  private

  ! timing events
  integer :: sortev, sort_genidx_ev, sort_rearrange_ev

  interface test
    module procedure test_species
  end interface

  interface move_window
    module procedure move_window_species
  end interface

  interface sort
    module procedure sort_species
  end interface

  interface name
    module procedure name_species
  end interface

  interface push_start_time
    module procedure push_start_time
  end interface

  interface bctype
    module procedure bctype_spec
  end interface bctype

  interface cleanup_buffers_spec
    module procedure cleanup_buffers_spec
  end interface

  interface init_buffers_spec
    module procedure init_buffers_spec
  end interface

  interface rearrange
    module procedure rearrange_species
  end interface

  public :: bctype, restart_write, restart_read

  public :: sortev, sort_genidx_ev, sort_rearrange_ev

  public :: move_window

  public :: init_buffers_spec, cleanup_buffers_spec

  public :: sort, rearrange




contains

!-----------------------------------------------------------------------------------------
subroutine test_species( this, no_co  )
!-----------------------------------------------------------------------------------------
!       tests the species for invalid values
!-----------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

  implicit none

!       dummy variables

  class( t_species ), intent(in) :: this
  class( t_node_conf ), intent(in) :: no_co


!       local variables

  integer  :: j, i

!       executable statements

     do i = 1, this%num_par
       do j= 1, p_x_dim
         if ( isnan(this%x(j,i)) .or. isinf(this%x(j,i))) then
ERROR('x(',j,',',i,') is invalid ',this%x(j,i), 'in node ', no_co%my_aid())
           call abort_program( -24 )
         endif
       enddo

       do j= 1, p_p_dim
         if ( isnan(this%p(j,i)) .or. isinf(this%p(j,i))) then
ERROR('p(',j,',',i,') is invalid ',this%p(j,i), 'in node ', no_co%my_aid())
           call abort_program( -25 )
         endif
       enddo

       if ( isnan(this%q(i)) .or. isinf(this%q(i))) then
ERROR('q(',i,') is invalid ',this%q(i), 'in node ', no_co%my_aid())
           call abort_program( -26 )
         endif

#ifdef __HAS_SPIN__

        do j= 1, p_s_dim
          if ( isnan(this%s(j,i)) .or. isinf(this%s(j,i))) then
ERROR('s(',j,',',i,') is invalid ',this%s(j,i), 'in node ', no_co%my_aid())
            call abort_program( -27 )
          endif
        enddo

#endif

     enddo

end subroutine test_species
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------
subroutine move_window_species( this, g_space, grid, jay, no_co, bnd_cross, node_cross, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

  use m_current_define

  implicit none

  ! dummy variables
  class( t_species ), intent(inout) :: this

  type( t_space ), intent(in) :: g_space
  class( t_grid ),  intent(in) :: grid
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  ! local variables
  integer, dimension(3, p_x_dim) :: ig_xbnd_inj
  integer :: i, k, nmove

  do i = 1, p_x_dim

    nmove = nx_move( g_space, i )

    this%move_num(i) = nmove

    if ( nmove > 0 ) then

       ! Get the new global boundaries

       ! The internal global boundaries are shifted from the simulation global boundaries
       ! by half a cell to simplify global particle position calculations
       this%g_box( p_lower, i ) = xmin( g_space, i ) + 0.5_p_double * this%dx(i)
       this%g_box( p_upper, i ) = xmax( g_space, i ) + 0.5_p_double * this%dx(i)

       this%total_xmoved( i ) = total_xmoved( g_space, i, p_lower )

       ! Shift particle positions

       !$omp parallel do
       do k = 1, this%num_par
         this%ix( i, k ) = this%ix( i, k ) - nmove
       enddo
       !$omp end parallel do

       if ( type( this%bnd_con, p_upper, i ) == p_bc_move_c .and. this%init_part ) then

          ig_xbnd_inj = grid%my_nx(:,1:p_x_dim)
          ig_xbnd_inj(p_lower, i) = ig_xbnd_inj(p_upper, i) + 1 - nmove
          ig_xbnd_inj(3,i) = nmove

          call this % inject( ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, send_msg, recv_msg )
       endif
    endif
  enddo

end subroutine move_window_species
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Grow a particle index list
!-----------------------------------------------------------------------------------------
subroutine grow_part_idx( buffer, new_size )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_part_idx ), intent(inout) :: buffer
  integer, intent(in) :: new_size

  integer, dimension(:), pointer :: tmp

  if ( buffer%nidx > 0 ) then

    ! Allocate new buffer
    call alloc( tmp, (/ new_size /) )

    ! Copy existing data to new buffer
    call memcpy( tmp, buffer % idx, buffer%nidx )

    ! Free old buffer and point it to the new buffer
    call freemem( buffer % idx )
    buffer % idx => tmp

  else
    ! Just reallocate the buffer with the new size
    call freemem( buffer%idx )
    call alloc( buffer%idx, (/ new_size /) )
  endif


end subroutine grow_part_idx
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine sort_species( this, n, t )
!-----------------------------------------------------------------------------------------
!     Selects the appropriate sorting routine
!-----------------------------------------------------------------------------------------

    implicit none

!       dummy variables

    class( t_species ), intent(inout) :: this

    integer, intent(in) :: n
    real(p_double), intent(in) :: t


!       local variables

    integer, dimension(:), pointer :: idx

!       executable statements

    ! sort particles at certain timesteps if required

    idx => null()

    if ( (this%n_sort > 0) .and. (this%num_par > 0) .and. (t > this%push_start_time)) then

       if (( mod( n, this%n_sort ) == 0 ) .and. (n > 0)) then

          call begin_event(sortev)

          call alloc( idx, (/ this%num_par /) )

          call begin_event( sort_genidx_ev )

          select case (p_x_dim)
           case (1)
             call generate_sort_idx_1d(this, idx)

           case (2)
             call generate_sort_idx_2d(this, idx)

           case (3)
             call generate_sort_idx_3d(this, idx)

          end select

          call end_event( sort_genidx_ev )

          call begin_event( sort_rearrange_ev )

          call rearrange_species( this, idx )

          call end_event( sort_rearrange_ev )

          call freemem( idx )

          call end_event(sortev)

       endif

    endif

  end subroutine sort_species
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine generate_sort_idx_1d( this, ip )
!-----------------------------------------------------------------------------------------
!       generate sort indexes for a 1d run
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species ),             intent(in) :: this
  integer, dimension(:), intent(out) :: ip

  integer, pointer, dimension(:) :: npic
  integer :: i
  integer :: index
  integer :: isum,ist
  integer :: n_grid

  npic => null()

  n_grid = this%my_nx_p( 3, 1 )

  call alloc(npic, (/n_grid/) )
  npic = 0

  ! sort by interpolation cell
  do i = 1,this%num_par
    index = this%ix(1,i)
    npic(index) = npic( index ) + 1
    ip(i)=index
  end do

  isum=0
  do i=1,n_grid
    ist=npic(i)
    npic(i)=isum
    isum=isum+ist
  end do

  do i=1,this%num_par
    index=ip(i)
    npic(index)=npic(index)+1
    ip(i)=npic(index)
  end do

  call freemem( npic )



end subroutine generate_sort_idx_1d

!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine generate_sort_idx_2d( this, ip )
!-----------------------------------------------------------------------------------------
!       generate sort indexes for a 2d run
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species ),             intent(in) :: this
  integer, dimension(:), intent(out) :: ip

  integer, pointer, dimension(:) :: npic
  integer :: i
  integer :: index
  integer :: isum,ist
  integer :: n_grid, n_grid_x, n_grid_y

  !       executable statements
  npic => null()

  n_grid_x = this%my_nx_p(3, 1)
  n_grid_y = this%my_nx_p(3, 2)
  n_grid = n_grid_x * n_grid_y

  call alloc(npic, (/n_grid/) )
  npic = 0

  ! sort by interpolation cell
  do i=1,this%num_par
     index = this%ix(1,i) + n_grid_x * (this%ix(2,i)-1)       ! sort by y then x
     npic(index) = npic(index) + 1
     ip(i)=index
  end do

  isum=0
  do i=1,n_grid
     ist = npic(i)
     npic(i) = isum
     isum = isum + ist
  end do

  ! isum must be the same as the total number of particles
  ! ASSERT(isum == this%num_par)

  do i=1,this%num_par
     index=ip(i)
     npic(index) = npic(index) + 1
     ip(i) = npic(index)
  end do

  call freemem( npic )



end subroutine generate_sort_idx_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine generate_sort_idx_3d( this, ip )
!-----------------------------------------------------------------------------------------
!       generate sort indexes for a 2d run
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  class( t_species ),             intent(in) :: this
  integer, dimension(:), intent(out) :: ip

  ! local variables

  integer, pointer, dimension(:) :: npic
  integer :: i
  integer :: index
  integer :: isum,ist
  integer :: n_grid, n_grid_x, n_grid_y, n_grid_xy, n_grid_z

  npic => null()

  n_grid_x = this%my_nx_p( 3, 1 )
  n_grid_y = this%my_nx_p( 3, 2 )
  n_grid_z = this%my_nx_p( 3, 3 )

  n_grid_xy = n_grid_x * n_grid_y
  n_grid    = n_grid_xy * n_grid_z

  call alloc(npic, (/n_grid/) )
  npic = 0

  ! sort by interpolation cell ( fastest push )
  do i=1,this%num_par
     ! sort by z then y then x
     index = this%ix(1,i) + &
             n_grid_x * (this%ix(2,i)-1) + &
             n_grid_xy * (this%ix(3,i)-1)

     npic(index) = npic(index) + 1
     ip(i)=index
  end do

  isum=0
  do i=1,n_grid
     ist = npic(i)
     npic(i) = isum
     isum = isum + ist
  end do

  do i=1,this%num_par
     index=ip(i)
     npic(index) = npic(index) + 1
     ip(i) = npic(index)
  end do

  call freemem( npic )



end subroutine generate_sort_idx_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Rearrange the order of particles using the array of new indices new_idx
! (OpenMP parallel))
!-----------------------------------------------------------------------------------------
subroutine rearrange_species( spec, new_idx )

  implicit none

  ! dummy variables

  class( t_species ), intent(inout) :: spec
  integer, dimension(:), intent(in) :: new_idx

  integer :: j, n_x_dim
  real(p_k_part), dimension(:), pointer   :: temp1_r
  real(p_k_part), dimension(:,:), pointer :: temp2_r
  integer, dimension(:,:), pointer        :: temp2_i


  ! Check if we are close to filling the particle buffer and grow it if so. Doing it here has zero
  ! overhead and may prevent having to do it during communications.
!  if ( spec%num_par > spec%reshape_grow_threshold * spec%num_par_max ) then
!
!    ! Get new size
!    spec%num_par_max = max( spec%num_par_max / spec%reshape_grow_threshold, &
!                            spec%num_par_max + p_spec_buf_block )
!
!    ! Make sure it is a multiple of 4 because of SSE code
!    spec%num_par_max = (( spec%num_par_max + 3 ) / 4) * 4
!  endif

  temp1_r => null()
  temp2_r => null()
  temp2_i => null()

  ! rearrange positions
  n_x_dim = spec%get_n_x_dims()
  call alloc( temp2_r, (/ n_x_dim, spec%num_par_max /) )

  ! Specific version for each dimensionality to manually unroll dimension loop
  select case ( n_x_dim )
    case (1)
      !$omp parallel do
      do j = 1, spec%num_par
        temp2_r( 1, new_idx(j) ) = spec%x(1,j)
      enddo
      !$omp end parallel do
    case (2)
      !$omp parallel do
      do j = 1, spec%num_par
        temp2_r( 1, new_idx(j) ) = spec%x(1,j)
        temp2_r( 2, new_idx(j) ) = spec%x(2,j)
      enddo
      !$omp end parallel do
    case (3)
      !$omp parallel do
      do j = 1, spec%num_par
        temp2_r( 1, new_idx(j) ) = spec%x(1,j)
        temp2_r( 2, new_idx(j) ) = spec%x(2,j)
        temp2_r( 3, new_idx(j) ) = spec%x(3,j)
      enddo
      !$omp end parallel do
    case (4)
      !$omp parallel do
      do j = 1, spec%num_par
        temp2_r( 1, new_idx(j) ) = spec%x(1,j)
        temp2_r( 2, new_idx(j) ) = spec%x(2,j)
        temp2_r( 3, new_idx(j) ) = spec%x(3,j)
        temp2_r( 4, new_idx(j) ) = spec%x(4,j)
      enddo
      !$omp end parallel do
    case (5)
      !$omp parallel do
      do j = 1, spec%num_par
        temp2_r( 1, new_idx(j) ) = spec%x(1,j)
        temp2_r( 2, new_idx(j) ) = spec%x(2,j)
        temp2_r( 3, new_idx(j) ) = spec%x(3,j)
        temp2_r( 4, new_idx(j) ) = spec%x(4,j)
        temp2_r( 5, new_idx(j) ) = spec%x(5,j)
      enddo
      !$omp end parallel do
  end select

  call freemem( spec%x )
  spec%x => temp2_r
  temp2_r => null()

  call alloc( temp2_i, (/ p_x_dim, spec%num_par_max /) )

  ! Specific version for each dimensionality to manually unroll dimension loop
  select case ( p_x_dim )
    case (1)
      !$omp parallel do
      do j = 1, spec%num_par
        temp2_i( 1, new_idx(j) ) = spec%ix(1,j)
      enddo
      !$omp end parallel do
    case (2)
      !$omp parallel do
      do j = 1, spec%num_par
        temp2_i( 1, new_idx(j) ) = spec%ix(1,j)
        temp2_i( 2, new_idx(j) ) = spec%ix(2,j)
      enddo
      !$omp end parallel do
    case (3)
      !$omp parallel do
      do j = 1, spec%num_par
        temp2_i( 1, new_idx(j) ) = spec%ix(1,j)
        temp2_i( 2, new_idx(j) ) = spec%ix(2,j)
        temp2_i( 3, new_idx(j) ) = spec%ix(3,j)
      enddo
      !$omp end parallel do
  end select

  call freemem( spec%ix )
  spec%ix => temp2_i
  temp2_i => null()

  ! rearrange momenta
  ! We always have 3 momenta components
  call alloc( temp2_r, (/ 3, spec%num_par_max /) )

  !$omp parallel do
  do j = 1, spec%num_par
     temp2_r( 1, new_idx(j) ) = spec%p(1,j)
     temp2_r( 2, new_idx(j) ) = spec%p(2,j)
     temp2_r( 3, new_idx(j) ) = spec%p(3,j)
  enddo
  !$omp end parallel do

  call freemem( spec%p )
  spec%p => temp2_r
  temp2_r => null()

#ifdef __HAS_SPIN__

  ! rearrange spin
  ! We always have 3 spin components
  call alloc( temp2_r, (/ 3, spec%num_par_max /) )

  !$omp parallel do
  do j = 1, spec%num_par
     temp2_r( 1, new_idx(j) ) = spec%s(1,j)
     temp2_r( 2, new_idx(j) ) = spec%s(2,j)
     temp2_r( 3, new_idx(j) ) = spec%s(3,j)
  enddo
  !$omp end parallel do

  call freemem( spec%s )
  spec%s => temp2_r
  temp2_r => null()

#endif

  ! rearrange charge
  call alloc( temp1_r, (/ spec%num_par_max /) )

  !$omp parallel do
  do j = 1, spec%num_par
     temp1_r( new_idx(j) ) = spec%q(j)
  enddo
  !$omp end parallel do

  call freemem( spec%q )
  spec%q => temp1_r
  temp1_r => null()

  ! reorder tags if present
  if ( spec%add_tag ) then
     call alloc( temp2_i, (/ 2, spec%num_par_max /) )

     !$omp parallel do
     do j = 1, spec%num_par
         temp2_i(1, new_idx(j)) = spec%tag(1,j)
       temp2_i(2, new_idx(j)) = spec%tag(2,j)
     enddo
     !$omp end parallel do

     call freemem( spec%tag )
     spec%tag => temp2_i
     temp2_i => null()

  endif

#ifdef __HAS_TRACKS__

  ! change tracked particle indexes if needed
  if ( spec%diag%ndump_fac_tracks > 0 ) then
    call update_indexes( spec%diag%tracks, new_idx )
  endif

#endif


end subroutine rearrange_species
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function bctype_spec( this, bnd, dim )
!-----------------------------------------------------------------------------------------
   implicit none

   class( t_species ), intent(in) :: this
   integer, intent(in) :: bnd, dim
   integer :: bctype_spec

   bctype_spec = type( this%bnd_con, bnd, dim )

end function bctype_spec
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine init_buffers_spec( )
!-----------------------------------------------------------------------------------------
!       initializes cache buffers for advance_deposit
!       routines
!-----------------------------------------------------------------------------------------
  implicit none

  ! init buffers for pistons - currently offline
  ! call init_buffers_piston()

end subroutine init_buffers_spec
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine cleanup_buffers_spec( )
!-----------------------------------------------------------------------------------------

  implicit none

  ! cleanup buffers for pistons  - currently offline
  ! call cleanup_buffers_piston()

#ifdef __HAS_TRACKS__
  ! cleanup buffers for tracks
  call cleanup_buffers_tracks()
#endif

end subroutine cleanup_buffers_spec
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!  gives the name of this species
!-----------------------------------------------------------------------------------------
function name_species( this )

  implicit none

  class(t_species), intent(in) :: this
  character(len=len_trim(this%name)) :: name_species

  name_species = trim(this%name)

end function name_species
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! returns this%push_start_time, the time at
! which you should begin pushing this species
!-----------------------------------------------------------------------------------------
function push_start_time( this )

  implicit none

  real( p_double ) :: push_start_time
  class( t_species ), intent(in) :: this

  push_start_time = this%push_start_time

end function push_start_time
!-----------------------------------------------------------------------------------------


end module m_species

!-----------------------------------------------------------------------------------------
! wrap the get_emf function 
! this function may be reimplemented for t_species subclasses in different
! simulation modes
!-----------------------------------------------------------------------------------------
subroutine get_emf_spec( this, emf, bp, ep, np, ptrcur, t )

  use m_species_define, only : t_species
  use m_emf_define, only : t_emf
  use m_emf_interpolate, only : get_emf
  use m_parameters, only : p_k_part, p_double

  implicit none 

  class( t_species ), intent(in) :: this
  class( t_emf ), intent(in), target :: emf
  real(p_k_part), dimension(:,:), intent(out) :: bp, ep
  integer, intent(in) :: np, ptrcur
  real(p_double), intent(in), optional :: t

  call get_emf( emf, bp, ep, this%ix(:,ptrcur:), this%x(:,ptrcur:), np, &
                this%interpolation )

end subroutine get_emf_spec

!-----------------------------------------------------------------------------------------
! Read information from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_species( this, input_file, def_name, periodic, if_move, grid, &
                             dt, read_prof, sim_options )

  use m_system
  use m_parameters
  use m_input_file

  use m_species_define

  use m_piston

  use m_species_diagnostics
  use m_species_boundary
  use m_species_udist
#ifdef __HAS_SPIN__
  use m_species_sdist
#endif
  use m_grid_define

  use m_profile

  implicit none

  class( t_species ),  intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  character(len = *),    intent(in) :: def_name
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ),               intent(in) :: grid

  type(t_options), intent(in) :: sim_options

  real(p_double), intent(in) :: dt
  logical, intent(in) :: read_prof

  character(len = p_max_spname_len) :: name
  character(len = 10 ) :: push_type

  integer                      :: num_par_max
  integer                      :: n_sort
  real(p_k_part)                      :: rqm

  real(p_k_part)                      :: q_real
#ifdef __HAS_SPIN__
  real(p_k_part)                       :: anom_mag_moment
#endif
  integer, dimension(p_x_dim)  :: num_par_x
  integer(p_int64), dimension(p_x_dim)  :: tot_par_x

  integer                              :: num_pistons

  logical                              :: add_tag
  logical                              :: free_stream

  real(p_double)                      :: push_start_time

  logical             :: if_collide
  logical             :: if_like_collide

  logical :: init_fields

  character(len = 16) :: init_type

  real(p_double) :: iter_tol
  logical :: rad_react

#ifdef __HAS_SPIN__

  namelist /nl_species/ name, num_par_max, n_sort, rqm, q_real, &
                        num_par_x, tot_par_x, push_type, &
                        push_start_time, num_pistons, &
                        add_tag, free_stream, init_fields, &
                        if_collide, if_like_collide, init_type, iter_tol, rad_react, &
                        anom_mag_moment

#else

  namelist /nl_species/ name, num_par_max, n_sort, rqm, q_real, &
                        num_par_x, tot_par_x, push_type, &
                        push_start_time, num_pistons, &
                        add_tag, free_stream, init_fields, &
                        if_collide, if_like_collide, init_type, iter_tol, rad_react
#endif

  integer :: i, ierr, piston_id

  ! Copy in relevant simulation parameters
  this%dt = dt
  this%g_nx( 1:p_x_dim ) = grid%g_nx( 1:p_x_dim )
  this%coordinates = grid%coordinates

#if defined(SIMD)
  push_type = "simd"
#else
  push_type = "standard"
#endif

  ! Initialization type
  init_type = "standard"

  num_par_max = 0
  rqm = 0
  q_real = 0
  num_par_x = -1
  tot_par_x = -1
  n_sort = 25

  push_start_time = -1.0_p_double
  name = trim(adjustl(def_name))

  num_pistons = 0

  add_tag = .false.

  free_stream = .false.
#ifdef __HAS_SPIN__
  anom_mag_moment = 0.00115965218076_p_k_part
#endif

  if_collide = .false.
  if_like_collide = .false.
  init_fields = .false.

  ! convergence threshold of root finding routines in exact pusher
  iter_tol = 1.0d-3

  ! if include radiation reaction (only valid for exact pusher)
  rad_react = .false.

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_species", ierr )
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        write(0,*) "Error reading species parameters"
      else
        write(0,*) "Error: species parameters missing"
      endif
      write(0,*) "aborting..."
    endif
    stop
  endif


  read (input_file%nml_text, nml = nl_species, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*) "   Error reading species parameters"
      write(0,*) "   aborting..."
    endif
    stop
  endif

  if (disp_out(input_file)) then
    SCR_ROOT("   Species name : ", trim(name))
  endif

  ! -------------- Allocate any sub-objects needed --------------------
  call this%allocate_objs()

  this%name = name
  this%free_stream = free_stream

  this%num_par_max      = num_par_max
  this%rqm              = rqm

  ! Validate num_par_x and tot_par_x
  do i = 1, p_x_dim
    if ( (num_par_x(i) <= 0) .and. (tot_par_x(i) <= 0) ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "   Error reading species parameters"
        write(0,*) "   Number of particles per cell must be >= 1 in all directions:"
        write(0,"(A,I0,A,I0,A,I0,A,I0)") "    -> num_par_x(",i,") = ", &
                 num_par_x(i), ", tot_par_x(",i,") = ", tot_par_x(i)
        write(0,*) "   aborting..."
      endif
      stop
    endif

    ! num_par_x overrides tot_par_x
    if ( num_par_x(i) <= 0 ) then
      this%num_par_x(i) = tot_par_x(i) / this % g_nx(i)
      this%tot_par_x(i) = tot_par_x(i)
    else
      this%num_par_x(i) = num_par_x(i)
      this%tot_par_x(i) = num_par_x(i) * this % g_nx(i)
    endif
  enddo

  this%n_sort           = n_sort

  ! Pusher type
  select case ( trim( push_type ) )
    case('standard')
       this%push_type = p_std
    case('vay')
      this%push_type = p_vay
    case('fullrot')
      this%push_type = p_fullrot
    case('euler')
      this%push_type = p_euler
    case('cond_vay')
      this%push_type = p_cond_vay
    case('cary')
      this%push_type = p_cary
    case('exact', 'analytic')
      this%push_type = p_exact
      this%iter_tol = iter_tol
    case('exact-rr', 'analytic-rr')
      this%push_type = p_exact_rr
      this%iter_tol = iter_tol

    case('simd')

#if defined(SIMD)

       this%push_type = p_simd

#else
       if ( mpi_node() == 0) then
        write(0,*) "   The code was not compiled with SIMD code support,"
         write(0,*) "   Please recompile the code with the appropriate flags."
         write(0,*) "   Error reading species parameters"
         write(0,*) "   aborting..."
       endif
       stop
#endif

    case('radcool')
       this%push_type = p_radcool

       if ( sim_options%omega_p0 <= 0. ) then
          if ( mpi_node() == 0) then
            write(0,*) "   Error: omega_p0 missing or invalid"
            write(0,*) "   When choosing push_type = 'radcool' you must also set omega_p0 in"
            write(0,*) "   the simulation section of the input file"
          endif
          stop
       endif

    case default
       if ( mpi_node() == 0) then
         write(0,*) "   Invalid push type:'", trim(push_type), "'"
         write(0,*) "   Allowed values are 'standard', 'vay', 'fullrot', 'euler', 'cond_vay', 'cary',"
         write(0,*) "   'exact'/'analytic', 'exact-rr'/'analytic-rr', 'simd' and 'radcool'"
         write(0,*) "   Error reading species parameters"
         write(0,*) "   aborting..."
       endif
       stop
  end select

  ! parameters for RR
  this%rad_react = rad_react
  if ( this%rad_react ) then
    if ( sim_options%omega_p0 <= 0. ) then
      if ( mpi_node() == 0) then
        write(0,*) "   Error: omega_p0 missing or invalid"
        write(0,*) "   When choosing pusher with radiation reaction you must"
        write(0,*) "   also set omega_p0 in the simulation section of the input file"
      endif
      stop
    else
      ! the leading coefficient is a characteristic time 2*e^2/(3*m_e*c^3)
      ! (e = elementary charge, m_e = electron static mass) in cgs unit.
      ! see Eq. (16.3) in Chapter 16 of Jackson
      this%k_rr = 6.266424752633625d-24 * sim_options%omega_p0
    endif

    if ( this%push_type /= p_std .and. this%push_type /= p_exact ) then
      if ( mpi_node() == 0) then
        write(0,*) "   Error: radiation reaction correction is only implemented for"
        write(0,*) "   'standard' and 'exact'/'analytic' push types"
      endif
      stop
    endif
  endif

  ! Free streaming is not implemented in simd code
  if ( free_stream ) then
    this%push_type = p_std
  endif

  this%push_start_time = push_start_time

  this%add_tag = add_tag

  this%num_pistons = num_pistons

  ! this is also used by some people even without collisions
  ! q_real and rqm having different signs is inconsistent, give rqm precedence
  this%q_real           = sign(q_real, rqm)

  if ( (this%push_type == p_exact .and. this%rad_react) .or. &
    this%push_type == p_exact_rr ) then
    if ( this%q_real == 0. ) then
      if ( mpi_node() == 0) then
        write(0,*) "   Error: q_real missing or invalid"
        write(0,*) "   When choosing push_type = 'exact' with radiation reaction"
        write(0,*) "   or push_type = 'exact-rr', you must also set q_real"
        write(0,*) "   aborting..."
      endif
      stop
    endif
  endif

#ifdef __HAS_COLLISIONS__

  this%if_collide       = if_collide
  this%if_like_collide  = if_like_collide

  if (( if_collide .or. if_like_collide ) .and. q_real == 0.0 ) then
     if ( mpi_node() == 0) then
       write(0,*) "   Error reading species parameters"
       write(0,*) "   When using collisions the user must also set q_real"
       write(0,*) "   aborting..."
     endif
     stop
  endif

#else

  if ( if_collide .or. if_like_collide ) then
     if ( mpi_node() == 0) then
       write(0,*) "   Error reading species parameters"
       write(0,*) "   Collisions are not supported in this version"
       write(0,*) "   aborting..."
     endif
     stop
  endif

#endif

  ! read momentum distribution
  call read_nml( this%udist, input_file )

#ifdef __HAS_SPIN__

  this%anom_mag_moment = anom_mag_moment
  ! read spin distribution
  call read_nml( this%sdist, input_file )

#endif

  ! Initialize fields from initial density / momentum profile
  this%init_fields = init_fields

  ! read profile definition unless read_prof was set to false
  if ( read_prof ) then
    this % source => new_source( init_type, input_file, this%coordinates )
    if ( .not. associated( this%source )) then
      if ( mpi_node() == 0) then
        write(0,*) "   Error reading species parameters"
        write(0,*) "   Invalid 'init_type' value (", trim(init_type), ")"
        write(0,*) "   aborting..."
      endif
      stop
    endif
  else
    this % source => null()
  endif

  ! read boundary condtions
  call read_nml( this%bnd_con, input_file, periodic, if_move, this%coordinates )

  ! read pistons
  if ( num_pistons > 0 ) then
    call alloc( this%pistons, (/ num_pistons /) )
    do piston_id=1, num_pistons
      call read_nml_piston( this%pistons(piston_id), input_file )
    enddo

    call check_piston( this%pistons, dt )

  endif

  ! read diagnostics
  call this%diag%read_input( input_file, sim_options%gamma )

end subroutine read_input_species
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine init_species( this, sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
                          no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
                          recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmax )
!-----------------------------------------------------------------------------------------
!       sets up this data structure from the given information
!-----------------------------------------------------------------------------------------

  use m_system
  use m_species_define
  use m_grid_define, only : t_grid
  use m_space
  use m_emf_define, only : t_emf
  use m_current_define, only : t_current
  use m_vdf_define
  use m_vdf_memory, only : alloc
  use m_vdf_comm, only : t_vdf_msg
  use m_node_conf
  use m_restart

  use m_parameters

  !use m_species, only : restart_read
  use m_species_push
  use m_species_extpush
  use m_species_current
  use m_species_diagnostics
  use m_species_boundary
  use m_species_initfields
  use m_profile
  use m_psource_file
  use m_time_step

  implicit none

  class( t_species ), intent(inout) :: this

  integer,     intent(in) :: sp_id
  integer, intent(in) :: interpolation
  logical, intent(in) :: grid_center
  class( t_grid ), intent(in)     :: grid
  type( t_space ),     intent(in) :: g_space
  class( t_emf ), intent(inout) :: emf
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_vdf, recv_vdf
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_spec, recv_spec
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  class(t_psource), pointer :: src
  type( t_restart_handle ), intent(in) :: restart_handle
  type(t_options), intent(in) :: sim_options
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: tmax 

  integer :: i

  ! Interpolation level is now defined by the global particles object
  this%interpolation = interpolation

  ! Field values are grid centered
  this%grid_center = grid_center

  if ( this%grid_center .and. this%push_type == p_simd ) then
    if (mpi_node()==0) then
      write(0,*) '(*warning*) Grid centered field values are not supported by SIMD code'
      write(0,*) '(*warning*) defaulting to the standard fortran pusher.'
    endif
    this%push_type = p_std
  endif

  ! initialize diagnostics
  ! diagnostics must be setup before inject_area because of particle tracking
  call this%diag%init(this%name,this% get_n_x_dims(), ndump_fac, interpolation, restart, &
                      restart_handle )

  ! allocate aux. array for openMP particle energy calculations
  ! (if OpenMP support was not compiled the function n_threads always returns 1)
  call alloc( this%energy, (/ n_threads(no_co) /) )

  ! initialize boundary conditions
  ! If the particle is a photon, characterized by this%rqm=0, we add it to the list of arguments
  if ( this % rqm /= 0.0_p_k_part ) then
    call setup( this%bnd_con, g_space, no_co, grid%coordinates, restart, restart_handle )
  else
    call setup( this%bnd_con, g_space, no_co, grid%coordinates, restart, restart_handle, this%rqm )
  endif

  ! Get unique id based on node grid position (independent of mpi rank)
  ! This needs to happen before calling inject_area
  this%ngp_id = no_co % ngp_id()


  ! Ensure that the particles buffer size, if set in the input file, is a multiple of 8 because
  ! of vector code
  this%num_par_max = ((this%num_par_max + 7) / 8 ) *8


  if ( restart ) then
     call this % restart_read( restart_handle )

    ! turn off raw read during restart
    src => this%source
    select type(src)
    class is ( t_psource_file )
      src%raw_read_done = .true.
    class default
      ! do nothing
    end select
  else
    this%sp_id = sp_id

    this%coordinates = grid%coordinates

    ! process position type
    select case ( interpolation )
    case ( p_linear, p_cubic )
      this%pos_type = p_cell_low
    case ( p_quadratic, p_quartic )
      this%pos_type = p_cell_near
    end select


    ! Global box boundaries are shifted from the simulation box boundaries by + 0.5*dx to simplify
    ! global particle position calculations, which are only used for injection (density profile)
    ! and diagnostics
    do i = 1, p_x_dim
      this%dx( i ) = emf%e%dx(i)
      this%g_box( p_lower, i ) = xmin( g_space, i ) + 0.5_p_double * this%dx(i)
      this%g_box( p_upper, i ) = xmax( g_space, i ) + 0.5_p_double * this%dx(i)
    enddo

    ! Cylindrical coordinates
    if ( this%coordinates == p_cylindrical_b .and. this%pos_type == p_cell_near ) then

      ! In cylindrical coordinates the algorithm puts the cylindrical symmetry axis (r = 0) at the
      ! center of cell 1, and the global simulation edge at -0.5*dx(2)

      ! This means that, for even interpolation, particles with gix = 2, x = -0.5 will be on
      ! the cylindrical axis, and for odd interpolation, this would be gix = 1, x = +0.0

      ! For even interpolation we therefore need to keep the minimum to be -0.5*dx(2)
      this%g_box( p_lower, p_r_dim ) = xmin( g_space, p_r_dim )
      this%g_box( p_upper, p_r_dim ) = xmax( g_space, p_r_dim )
    endif



    this%g_nx( 1:p_x_dim ) = grid%g_nx( 1:p_x_dim )

    do i = 1, p_x_dim
      this%total_xmoved( i ) = total_xmoved( g_space, i, p_lower )
    enddo

    this%my_nx_p(:,1:p_x_dim) = grid%my_nx(:,1:p_x_dim)

    ! allocate needed buffers
    if ( this%num_par_max > 0 .and. this%init_part ) then
      ! Check if a superclass did not already allocate the memory
      if (.not. associated(this%x)) call this%init_buffer( this%num_par_max )
    endif

    ! actual initialization of particles of this species
    this%num_par = 0
    this%num_created = 0

    call this%inject( this%my_nx_p, jay, no_co, bnd_cross, node_cross, send_spec, recv_spec )

    ! Initialize fields associated with the initial distribution / momentum
    if ( this%init_fields ) then
      call init_fields( this, emf, no_co, grid, send_vdf, recv_vdf )
    endif

  endif

  ! Setup data that is not in the restart file
  call this % set_dudt()

  ! Set pointer to current deposition functions
  select case ( this%interpolation )
  case( p_linear )
    this % dep_current_1d => dep_current_1d_s1
    this % dep_current_2d => dep_current_2d_s1
    this % dep_current_3d => dep_current_3d_s1
  case( p_quadratic )
    this % dep_current_1d => dep_current_1d_s2
    this % dep_current_2d => dep_current_2d_s2
    this % dep_current_3d => dep_current_3d_s2
  case( p_cubic )
    this % dep_current_1d => dep_current_1d_s3
    this % dep_current_2d => dep_current_2d_s3
    this % dep_current_3d => dep_current_3d_s3
  case( p_quartic )
    this % dep_current_1d => dep_current_1d_s4
    this % dep_current_2d => dep_current_2d_s4
    this % dep_current_3d => dep_current_3d_s4
  case default
    ERROR("Invalid interpolation type")
    call abort_program( p_err_invalid )
  end select

  ! calculate m_real (this is currently only used by collisions)
  this%m_real = this%q_real * this%rqm
  if( this%q_real /= 0.0_p_k_part ) this%rq_real = 1.0_p_k_part / this%q_real

end subroutine init_species
!-----------------------------------------------------------------------------------------

subroutine set_dudt_spec( this ) 

  use m_species_define, only : t_species
  use m_species_push
  use m_species_extpush
  use m_parameters

  implicit none

  class( t_species ), intent(inout) :: this

  ! Set pointer to dudt function
  select case ( this % push_type )
  case(p_std)
    this % dudt => dudt_boris
  case(p_vay)
    this % dudt => dudt_vay
  case(p_fullrot)
    this % dudt => dudt_fullrot
  case(p_euler)
    this % dudt => dudt_euler
  case(p_cond_vay)
    this % dudt => dudt_cond_vay
  case(p_cary)
    this % dudt => dudt_cary
  case(p_exact)
    this % dudt => dudt_exact
  case(p_exact_rr)
    this % dudt => dudt_exact_rr
  case(p_simd)
    ! SIMD pusher uses its own dudt (currently boris)
    this % dudt => null()
  case(p_radcool)
    ! Rad. cooling pusher uses its own dudt
    this % dudt => null()
  case default
    ERROR("Invalid pusher type")
    call abort_program( p_err_invalid )
  end select

#ifdef __HAS_SPIN__
  ! Set pointer to dsdt function
  this % dsdt => dsdt_vieira
#endif

end subroutine set_dudt_spec

!-----------------------------------------------------------------------------------------
! This routine dispatches the code to the appropriate push.
! Currently the following of push are available:
!  i)   Standard PIC push (Fortran or Hardware Optimized)
!  ii)  Accelerate pusher (used for initializning beams )
!  iii)  Radiation cooling pusher
!-----------------------------------------------------------------------------------------
subroutine push_species( this, emf, jay, t, tstep, tid, n_threads, options )

  use m_system
  use m_parameters

  use m_species_define, only : t_species
  use m_emf_define, only : t_emf
  use m_vdf_define, only : t_vdf

  use m_time_step
  use m_species_push
  use m_species_vpush
  use m_species_accelerate
  use m_species_udist, only : use_accelerate, use_q_incr
  use m_species_radcool

  implicit none

  class( t_species ), intent(inout) :: this
  class( t_emf ), intent(inout)  ::  emf
  type( t_vdf ), intent(inout) :: jay

  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  type( t_options ), intent(in) :: options

  integer, intent(in) :: tid        ! local thread id
  integer, intent(in) :: n_threads  ! total number of threads

  integer :: push_type

  ! if before push start time return silently
  if ( t >= this%push_start_time ) then

    push_type = this%push_type
    if ( use_accelerate( this%udist, n(tstep) ) ) push_type = p_beam_accel
    if ( use_q_incr( this%udist, n(tstep) ) ) push_type = p_beam_accel

    select case ( push_type )

#ifdef SIMD
      case ( p_simd )
        call vadvance_deposit( this, emf, jay, t, tstep, tid, n_threads )
#endif

      case ( p_beam_accel )
        call accelerate_deposit( this, jay, t, tstep, tid, n_threads )

      case ( p_radcool )
         call advance_deposit_radcool( this, emf, jay, t, tstep, options%omega_p0, tid, n_threads )

      case default
         call advance_deposit( this, emf, jay, t, tstep, tid, n_threads )
  
    end select

  endif

  ! TODO: REMOVE: this%n_current_spec = this%n_current_spec + 1
end subroutine push_species
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! cleanup all dynamic memory allocated by this object
!-----------------------------------------------------------------------------------------
subroutine cleanup_species( this )

#include "memory/memory.h"

  use m_system
  use m_parameters
  use m_species_define
  use m_piston, only : freemem
  use m_species_diagnostics, only : cleanup

  implicit none

  class( t_species ), intent(inout) :: this

  ! cleanup memory used by the particle buffers
  call freemem( this%x )
  call freemem( this%ix )
  call freemem( this%p )
  call freemem( this%q )
  call freemem( this%tag )
  call freemem( this%energy )
#ifdef __HAS_SPIN__
  call freemem( this%s )
#endif

  ! these are just pointers to data structures in the emf object so
  ! there is no need to deallocate
  nullify( this%E )
  nullify( this%B )
  ! cleanup subobjects
  if(associated(this%source )) then
    call this%source%cleanup( )
  endif
  
  call cleanup( this%diag )
  ! cleanup piston data
  if ( this%num_pistons > 0 ) then
    call freemem( this%pistons )
  endif

end subroutine cleanup_species
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Inject particles in the grid area specified by ig_xbnd_inj
!-----------------------------------------------------------------------------------------
subroutine species_inject( this, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, send_msg, recv_msg )

  use m_species_define
  use m_species_tag
  use m_current_define
  use m_node_conf

  implicit none

  class( t_species ), intent(inout) :: this
  integer, dimension(:,:), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout)   :: jay
  class( t_node_conf ), intent(in)     :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: num_inj

  if ( associated( this % source ) .and. this%init_part ) then
    ! Inject particles using the initialization object
    ! This will also set particle momentum
    num_inj = this % source % inject( this, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, send_msg, recv_msg )

    if ( num_inj > 0 ) then

      ! if required set tags of particles
      if ( this%add_tag ) then
        call set_tags( this, this%num_par+1, this%num_par + num_inj )
      endif

      ! increase num_created counter
      this%num_created = this%num_created + num_inj
      this%num_par = this%num_par + num_inj

    endif
  endif

end subroutine species_inject
!-----------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine create_particle_single_cell( this, ix, x, q )
!---------------------------------------------------------------------------------------------------
! creates a single particle in this this
! does not initialize particle momentum (to be done with set_momentum)
! ix and x are expected to match the this pos_type:
!  - cell_lower => 0 <= x < 1
!  - cell_near  => -0.5 <= x < 0.5
!---------------------------------------------------------------------------------------------------
  use m_system

  use m_species_define
  use m_species_udist
  use m_species_memory
  use m_species_tag

  #ifdef __HAS_TRACKS__
  use m_species_tracks
  #endif

  use m_space
  use m_math
  use m_random

  use m_parameters

   implicit none

  ! dummy variables

   class(t_species), intent(inout) :: this
   integer, dimension(:), intent(in) :: ix
   real(p_k_part), dimension(:), intent(in) :: x
   real(p_k_part),               intent(in) :: q


   ! local variables

   integer :: bseg_size
   integer :: n_x_dim

   ! executable statements

   if ( this%num_par + 1 > this%num_par_max ) then
     ! buffer is full, resize it
     bseg_size = max(this%num_par_max / 16, p_spec_buf_block )
     call this%grow_buffer( this%num_par_max + bseg_size )
   endif

   this%num_par = this%num_par + 1

   ! if required set tags of particles
   if (this%add_tag) then
     call set_tags( this, this%num_par )
   endif

   n_x_dim = this%get_n_x_dims()
   this%x(1:n_x_dim,this%num_par)  = x(1:n_x_dim)
   this%ix(1:p_x_dim,this%num_par) = ix(1:p_x_dim)

   this%q(  this%num_par) = q
   this%num_created = this%num_created + 1


 end subroutine create_particle_single_cell
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine create_particle_single_cell_p( this, ix, x, p, q )
!---------------------------------------------------------------------------------------------------
! creates a single particle in this species and initialized momentum to the given value
! ix and x are expected to match the species pos_type:
!  - cell_lower => 0 <= x < 1
!  - cell_near  => -0.5 <= x < 0.5
!---------------------------------------------------------------------------------------------------
   use m_system

   use m_species_define
   use m_species_tag
   use m_species_udist
   use m_species_memory

   #ifdef __HAS_TRACKS__
   use m_species_tracks
   #endif

   use m_space
   use m_math
   use m_random

   use m_parameters

   implicit none

  ! dummy variables

   class(t_species), intent(inout) :: this
   integer, dimension(:), intent(in) :: ix
   real(p_k_part), dimension(:), intent(in) :: x
   real(p_k_part), dimension(:), intent(in) :: p
   real(p_k_part),               intent(in) :: q


   ! local variables

   integer :: bseg_size
   integer :: n_x_dim

   ! executable statements

#if 0

  ! (* debug *)
  integer :: i

  do i = 1, p_x_dim
    if ( x(i) < -0.5_p_k_part .or. x(i) >= +0.5_p_k_part ) then
      write(0,*) '(*error*) Invalid particle being created, position x = ', x(1:p_x_dim)
      call abort_program()
    endif
  enddo

  if ( q == 0 ) then
     write(0,*) '(*error*) Invalid particle being created, charge q = ', q
     call abort_program()
  endif
  ! (* debug *)

#endif

   if ( this%num_par + 1 > this%num_par_max ) then
     ! buffer is full, resize it
     bseg_size = max(this%num_par_max / 16, p_spec_buf_block )
     call this%grow_buffer( this%num_par_max + bseg_size )
   endif

   this%num_par = this%num_par + 1

   ! if required set tags of particles
   if (this%add_tag) then
     call set_tags( this, this%num_par )
   endif

   n_x_dim = this%get_n_x_dims()
   this%x(1:n_x_dim,this%num_par)  = x(1:n_x_dim)
   this%ix(1:p_x_dim,this%num_par) = ix(1:p_x_dim)
   this%q(  this%num_par) = q
   this%p(1,this%num_par) = p(1)
   this%p(2,this%num_par) = p(2)
   this%p(3,this%num_par) = p(3)

   this%num_created = this%num_created + 1


 end subroutine create_particle_single_cell_p
!-------------------------------------------------------------------------------

#ifdef __HAS_SPIN__
!-------------------------------------------------------------------------------
subroutine create_particle_single_cell_p_s( this, ix, x, p, s, q )
!---------------------------------------------------------------------------------------------------
! creates a single particle in this species and initialized momentum to the given value
! ix and x are expected to match the species pos_type:
!  - cell_lower => 0 <= x < 1
!  - cell_near  => -0.5 <= x < 0.5
!---------------------------------------------------------------------------------------------------
   use m_system

   use m_species_define
   use m_species_tag
   use m_species_udist
   use m_species_memory

   #ifdef __HAS_TRACKS__
   use m_species_tracks
   #endif

   use m_space
   use m_math
   use m_random

   use m_parameters

   implicit none

  ! dummy variables

   class(t_species), intent(inout) :: this
   integer, dimension(:), intent(in) :: ix
   real(p_k_part), dimension(:), intent(in) :: x
   real(p_k_part), dimension(:), intent(in) :: p
   real(p_k_part), dimension(:), intent(in) :: s
   real(p_k_part),               intent(in) :: q


   ! local variables

   integer :: bseg_size

   ! executable statements

#if 0

  ! (* debug *)
  integer :: i

  do i = 1, p_x_dim
    if ( x(i) < -0.5_p_k_part .or. x(i) >= +0.5_p_k_part ) then
      write(0,*) '(*error*) Invalid particle being created, position x = ', x(1:p_x_dim)
      call abort_program()
    endif
  enddo

  if ( q == 0 ) then
     write(0,*) '(*error*) Invalid particle being created, charge q = ', q
     call abort_program()
  endif
  ! (* debug *)

#endif

   if ( this%num_par + 1 > this%num_par_max ) then
     ! buffer is full, resize it
     bseg_size = max(this%num_par_max / 16, p_spec_buf_block )
     call this%grow_buffer( this%num_par_max + bseg_size )
   endif

   this%num_par = this%num_par + 1

   ! if required set tags of particles
   if (this%add_tag) then
     call set_tags( this, this%num_par )
   endif

   this%x(1:p_x_dim,this%num_par)  = x(1:p_x_dim)
   this%ix(1:p_x_dim,this%num_par) = ix(1:p_x_dim)
   this%q(  this%num_par) = q
   this%p(1,this%num_par) = p(1)
   this%p(2,this%num_par) = p(2)
   this%p(3,this%num_par) = p(3)
   this%s(1,this%num_par) = s(1)
   this%s(2,this%num_par) = s(2)
   this%s(3,this%num_par) = s(3)

   this%num_created = this%num_created + 1


 end subroutine create_particle_single_cell_p_s
!-------------------------------------------------------------------------------
#endif

!---------------------------------------------------------------------------------------------------
! Used for CUDA (and similar) situations where data may have to be copied from another location.
!---------------------------------------------------------------------------------------------------
subroutine fill_data_species(this, tstep)
  use m_parameters
  use m_time_step
  use m_species_define

  implicit none
  class( t_species ), intent(inout) :: this
  type( t_time_step ),   intent(in) :: tstep

  ! If used, the general format for this function would be along the lines of:
  !
  ! if( n(tstep) /= this%n_current_spec) then
  !   ! get the data needed 
  !   call ....
  !   ! update the status counter so that the data is only copied as needed
  !   this%n_current_spec = n(tstep)
  ! endif

  ! put a fatal error here in case any subclasses forgot to increment 'n_current_spec'
  !   this error is most likely from a subclass not having the line
  !      "this%n_current_spec=this%n_current_spec+1" 
  !   after the a species 'push' (e.g. 'push_species'). Ask Ricardo or Anton or Adam if needed.
  if( n(tstep) /= this%n_current_spec) then
    ERROR( "Error in fill_data_species.. see 'fill_data_species' is os-species.f03")
    call abort_program(p_err_rstrd)
  endif

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_write_species( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       write object information into a restart file
!-----------------------------------------------------------------------------------------

  use m_restart
  use m_parameters
  use m_species_define
  use m_species_boundary
  use m_species_diagnostics

  implicit none

  class( t_species ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for species object.'
  integer :: ierr

  ! these need to happen before the species object data is saved
  ! because of the structure of the setup routine
  call restart_write( this%diag, restart_handle )
  call restart_write( this%bnd_con, restart_handle )

  ! write restart information for species

  restart_io_wr( p_spec_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%name, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%sp_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%num_par_max, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%num_par, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%rqm, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%num_par_x, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%tot_par_x, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%dx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%coordinates, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%push_start_time, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%add_tag, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%free_stream, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%energy, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%pos_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%g_box, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%g_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%dx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%total_xmoved, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%my_nx_p, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%push_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! Collision data
  restart_io_wr( this%if_collide, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%if_like_collide, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%q_real, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%num_created, restart_handle, ierr)
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  if ( this%num_par > 0 ) then
     restart_io_wr( this%x(:,1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

     restart_io_wr( this%ix(:,1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

     restart_io_wr( this%p(:,1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

     restart_io_wr( this%q(1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

#ifdef __HAS_SPIN__
     restart_io_wr( this%s(:,1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
#endif

     if ( this%add_tag ) then
        restart_io_wr( this%tag(:,1:this%num_par), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
     endif
  endif

end subroutine restart_write_species
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine restart_read_species( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       read object information from a restart file
!-----------------------------------------------------------------------------------------

  use m_restart
  use m_parameters
  use m_species_define

  implicit none

  class( t_species ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for species object.'
  character(len=len(p_spec_rst_id)) :: rst_id
  integer :: ierr
  integer :: num_par_max

  ! The restart_read calls for diag (tracks only) and bnd_con happen in
  ! the respective setup functions

  num_par_max = this%num_par_max

  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  ! check if restart file is compatible
  if ( rst_id /= p_spec_rst_id) then
    ERROR('Corrupted restart file, or restart file ')
    ERROR('from incompatible binary (species)')
    call abort_program(p_err_rstrd)
  endif

  ! read restart data
  restart_io_rd( this%name, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%sp_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%num_par_max, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%num_par, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%rqm, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%num_par_x, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%tot_par_x, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%dx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%coordinates, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%push_start_time, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%add_tag, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%free_stream, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%energy, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%pos_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%g_box, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%g_nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%dx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%total_xmoved, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%my_nx_p, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%push_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  ! Collision data
  restart_io_rd( this%if_collide, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%if_like_collide, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%q_real, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%num_created, restart_handle, ierr)
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  ! set the values that can be changed in the input deck between restarts

  ! num_par_max - The user may only request to grow this

  if (num_par_max > this%num_par_max) then
    this%num_par_max = num_par_max
  endif

  ! initialize buffers
  call this % init_buffer( this%num_par_max )


  ! read particle data, if available
  if ( this%num_par > 0 ) then

     restart_io_rd( this%x(:,1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstrd )

     restart_io_rd( this%ix(:,1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstrd )

     restart_io_rd( this%p(:,1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstrd )

     restart_io_rd( this%q(1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstrd )

#ifdef __HAS_SPIN__
     restart_io_rd( this%s(:,1:this%num_par), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstrd )
#endif

     if ( this%add_tag ) then
        restart_io_rd( this%tag(:,1:this%num_par), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstrd )
     endif
  endif

end subroutine restart_read_species
!-----------------------------------------------------------------------------------------
