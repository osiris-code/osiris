!#define __DEBUG__ 1

#include "os-preprocess.fpp"
#include "os-config.h"

module m_species_boundary

#include "memory/memory.h"

use m_species_define
use m_species_current
use m_species_comm

#ifdef __HAS_TRACKS__
use m_species_tracks
#endif

use m_space
use m_node_conf

use m_system
use m_parameters
use m_utilities
use m_random

use m_vdf_define
use m_current_define

private

! string to id restart data
character(len=*), parameter :: p_spbnd_rst_id = "species bound rst data - 0x0000"

integer, parameter :: p_idx_block = 16384

! indexes of particles crossing lower/upper boundary
! type( t_part_idx ), dimension(2), save :: bnd_cross

! indexes of particles leaving the node
! type( t_part_idx ), save :: node_cross

interface read_nml
  module procedure read_nml_spe_bound
end interface

interface setup
  module procedure setup_spe_bound
end interface

interface restart_write
  module procedure restart_write_spe_bound
end interface

interface restart_read
  module procedure restart_read_spe_bound
end interface

interface type
  module procedure type_spe_bound
  module procedure type_bnd_dim_spe_bound
end interface

interface check_boundary
  module procedure check_boundary
end interface

interface process_phys_boundary
  module procedure process_phys_boundary
end interface

interface check_node_cross
  module procedure check_node_cross
end interface

interface thermal_bath
  module procedure thermal_bath
end interface

! declare things that should be public
public :: read_nml, setup
public :: restart_write, restart_read
public :: type

public :: check_boundary, check_node_cross, process_phys_boundary, thermal_bath

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------
! Read necessary information from input file
!-------------------------------------------------------------------------------
subroutine read_nml_spe_bound( this, input_file, periodic, if_move, coordinates, rqm )

  use m_input_file

  implicit none

  class( t_spe_bound ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic, if_move
  integer,             intent(in) :: coordinates
  real(p_k_part), intent(in), optional  :: rqm

  ! integer, dimension(2, p_x_dim) :: type
  character(20), dimension(2, p_x_dim) :: type

  real(p_k_part), dimension(p_p_dim,2, p_x_dim) :: uth_bnd
  real(p_k_part), dimension(p_p_dim,2, p_x_dim) :: ufl_bnd
  character(len=20), dimension(2,p_max_dim) :: thermal_type

  namelist /nl_spe_bound/ type, uth_bnd, ufl_bnd, thermal_type

  integer :: ierr, i, j
  logical :: if_photons

!       executable statements

  ! type = p_bc_invalid
  type = "---"

  thermal_type = "thermal"
  uth_bnd = 0
  ufl_bnd = 0

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_spe_bound", ierr )
  if ( ierr /= 0 ) then
    if (ierr < 0) then
      print *, "Error reading spe_bound parameters"
    else
      print *, "Error: spe_bound parameters missing"
    endif
    print *, "aborting..."
    stop
  endif

  read (input_file%nml_text, nml = nl_spe_bound, iostat = ierr)
  if (ierr /= 0) then
    print *, ""
    print *, "   Error reading spe_bound parameters"
    print *, "   aborting..."
    stop
  endif

  ! Read boundary condition types
  do i = 1, p_x_dim
  do j = 1, 2
    select case ( trim( type(j,i) ) )
      case ( "absorbing","open" )
        this%type(j,i) = p_bc_absorbing
      case ( "axial" )
        this%type(j,i) = p_bc_cyl_axis
      case ( "specular","reflecting" )
        this%type(j,i) = p_bc_specular
      case ( "thermal" )
        this%type(j,i) = p_bc_thermal
      case ( "periodic" )
        this%type(j,i) = p_bc_periodic

      case ( "---" )
        this%type(j,i) = p_bc_invalid
      case default
        print "(A,I0,A,I0)", "   Error reading spe_bound parameters, dir = ", i, ", wall = ", j
      print *, '   Error: unknown boundary condition type, "', trim(type(i,j)), '"', &
               ' valid values are "axial","absorbing"/"open", "specular"/"reflecting", and "thermal"'
      print *, '   aborting...'
      stop
    end select
  enddo
  enddo

  ! Process thermal bath boundaries
  this%uth_bnd = uth_bnd
  this%ufl_bnd = ufl_bnd

  ! process thermal_type
  do i = 1, p_x_dim
    do j = p_lower, p_upper
    select case ( trim( thermal_type(j,i) ) )

     case ( "thermal", "half max" )
       this%thermal_type(j,i) = p_thermal

     case ( "random dir" )
       this%thermal_type(j,i) = p_random_dir

     case ( "waterbag" )
       this%thermal_type(j,i) = p_waterbag

     case default
       print *, "   Error reading spe_bound parameters"
       print *, "   invalid thermal_type, valid values are 'thermal', 'random dir',"
       print *, "   'waterbag', and 'half max',"
       print *, "   aborting..."
       stop
    end select
    enddo
  enddo

  !boolean to know if the species is a photon : this%rqm = 0
  if_photons = .false.
  if( present(rqm) ) then
    if_photons = ( (rqm==0.0_p_k_part) .and. (coordinates==p_cylindrical_b) )
  endif

  do i = 1, p_x_dim
    !periodic boundaries and moving windows override local settings
    if ((.not. periodic(i)) .and. (.not. if_move(i))) then
       do j = 1, 2
         if (this%type(j,i) == p_bc_invalid) then
           print *, ""
           print *, "   Error reading spe_bound parameters, dir = ", i, ", wall = ", j
           print *, "   Invalid or no boundary conditions specified."
           print *, "   aborting..."
           stop
         elseif (this%type(j,i) == p_bc_thermal .and. uth_bnd(i,j,i) <= 0 ) then
           print *, ""
           print *, "   Error reading spe_bound parameters, dir = ", i, ", wall = ", j
           print *, "   When specifying thermal bath boundary conditions the pth_bnd "
           print *, "   parameter must be > 0.0 for the specified boundary"
           print *, "   aborting..."
           stop
         elseif (this%type(j,i) == p_bc_periodic) then
           if ( .not. if_photons ) then
            print *, ""
            print *, "   Error reading spe_bound parameters, dir = ", i
            print *, "   Periodic boundaries cannot be specified unless global periodic"
            print *, "   boundaries were set in the node_conf section"
            print *, "   aborting..."
            stop
           else
            SCR_ROOT("   Warning, the code is running with a specific BC for photons ")
            SCR_ROOT("   It is only meant to be used with the qed_cyl simulation mode")
            SCR_ROOT("   The BC is periodical in the radial direction, only for photons")
           endif
         endif
       enddo
    endif
  enddo

  ! check on cylindrical boundaries
  do i = 1, p_x_dim
   do j = 1, 2
     if ((this%type(j,i) == p_bc_cyl_axis) .and. &
       ((i /= p_r_dim) .or. ( j /= p_lower ) .or. (coordinates /= p_cylindrical_b))) then
       print *, ""
       print *, "   Error reading spe_bound parameters, dir = ", i
       print *, "   Cylindrical axial boundaries can only be specified for"
       print *, "   the lower boundary of dimension ",p_r_dim
       print *, "   aborting..."
       stop
     endif
   enddo
  enddo


end subroutine read_nml_spe_bound
!-------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Cleanup
!---------------------------------------------------------------------------------------------------
!subroutine cleanup_spe_bound( this )
!
!  implicit none
!
!  class ( t_spe_bound ), intent(inout) :: this
!
!end subroutine cleanup_spe_bound
!---------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! sets up this data structure from the given information
!-------------------------------------------------------------------------------
subroutine setup_spe_bound( this, g_space, no_co, coordinates, restart, restart_handle, rqm )

  use m_restart

  implicit none

  class( t_spe_bound ), intent( inout )  ::  this

  type( t_space ),     intent(in) :: g_space
  class( t_node_conf ), intent(in) :: no_co
  integer,             intent(in) :: coordinates
  logical,             intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  real( p_k_part ), intent(in), optional :: rqm

  ! local variables

  logical, dimension(p_x_dim) :: ifpr_l, if_move_l
  logical :: dep_current
  integer :: i_bnd, i_x_dim
  logical :: if_photons

  ! executable statements



  if ( restart ) then

     call restart_read( this, restart_handle )

  else

   ifpr_l = no_co%ifpr()
   if_move_l = if_move(g_space)

   ! overwrite b.c. flags if motion of simulation space
   do i_x_dim = 1, p_x_dim
    if ( if_move_l(i_x_dim) ) then
       this%type( p_lower , i_x_dim ) = p_bc_absorbing
       this%type( p_upper , i_x_dim ) = p_bc_move_c
    endif
   enddo

   ! overwrite b.c. flags for boundaries to other nodes
   do i_x_dim = 1, p_x_dim
    do i_bnd=1, 2
       if ( neighbor( no_co, i_bnd, i_x_dim ) > 0 ) then
        this%type( i_bnd, i_x_dim ) = p_bc_other_node
       endif
    enddo
   enddo

   ! overwrite b.c.-flags if periodicity
   do i_x_dim = 1, p_x_dim
    if ( ifpr_l(i_x_dim) ) then
       this%type( : , i_x_dim ) = p_bc_periodic
    endif
   enddo

   !boolean to know if the species is a photon : this%rqm = 0
   if_photons = .false.
   if( present(rqm) ) then
     if_photons = ( (rqm==0.0_p_k_part) .and. (coordinates==p_cylindrical_b) )
   endif

   ! do additional cylindrical coordinates setup
   if ( coordinates == p_cylindrical_b ) then
     if ( no_co%on_edge( p_r_dim, p_lower )) then
      this%type(p_lower, p_r_dim) = p_bc_cyl_axis
     endif

     ! check if we have periodic boundaries along radial direction
     if ((this%type(p_lower, p_r_dim) == p_bc_periodic) .or. &
       (this%type(p_upper, p_r_dim) == p_bc_periodic)) then
       if (.not. if_photons) then
        ERROR('no periodic b.c. for radial direction')
        call abort_program( p_err_invalid )
       endif
     endif

   endif

  endif

  ! Check if any boundary condition deposits current and initialize tmp buffers if needed

  dep_current = .false.
  do i_x_dim = 1, p_x_dim
    if ( this%type(p_lower, i_x_dim) == p_bc_specular .or. &
         this%type(p_lower, i_x_dim) == p_bc_thermal .or. &
         this%type(p_upper, i_x_dim) == p_bc_specular .or. &
         this%type(p_upper, i_x_dim) == p_bc_thermal .or. &
         this%type(p_upper, i_x_dim) == p_bc_move_c ) then
       dep_current = .true.
       exit
    endif
  enddo

  if ( dep_current ) call this % init_tmp_buf_current()

end subroutine setup_spe_bound
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Write object information into a restart file
!-------------------------------------------------------------------------------
subroutine restart_write_spe_bound( this, restart_handle )

  use m_restart

  implicit none

  ! dummy variables

  class ( t_spe_bound ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  ! local variables
  integer :: ierr

  ! executable statements
  restart_io_wr( p_spbnd_rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then
    ERROR('error writing restart data for species boundary object.')
  call abort_program(p_err_rstwrt)
  endif

  restart_io_wr( this%type, restart_handle, ierr )
  if ( ierr/=0 ) then
     ERROR('error writing restart data for spe_bound object.')
    call abort_program(p_err_rstwrt)
  endif
  restart_io_wr( this%uth_bnd, restart_handle, ierr )
  if ( ierr/=0 ) then
     ERROR('error writing restart data for spe_bound object.')
    call abort_program(p_err_rstwrt)
  endif
  restart_io_wr( this%ufl_bnd, restart_handle, ierr )
  if ( ierr/=0 ) then
     ERROR('error writing restart data for spe_bound object.')
    call abort_program(p_err_rstwrt)
  endif

end subroutine restart_write_spe_bound
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Read object information from a restart file
!-------------------------------------------------------------------------------
subroutine restart_read_spe_bound( this, restart_handle )

  use m_restart

  implicit none

  ! dummy variables

  class ( t_spe_bound ), intent(inout)  ::  this
  type( t_restart_handle ), intent(in) :: restart_handle

  ! local variables
  character(len=len(p_spbnd_rst_id)) :: rst_id
  integer :: ierr

  ! executable statements
  restart_io_rd( rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then
  ERROR('error reading restart data for species boundary object.')
  call abort_program(p_err_rstrd)
  endif

  ! check if restart file is compatible
  if ( rst_id /= p_spbnd_rst_id) then
  ERROR('Corrupted restart file, or restart file ')
  ERROR('from incompatible binary (species boundary)')
  call abort_program(p_err_rstrd)
  endif

  restart_io_rd( this%type, restart_handle, ierr )
  if ( ierr/=0 ) then
     ERROR('error reading restart data for spe_bound object.')
    call abort_program(p_err_rstrd)
  endif
  restart_io_rd( this%uth_bnd, restart_handle, ierr )
  if ( ierr/=0 ) then
     ERROR('error reading restart data for spe_bound object.')
    call abort_program(p_err_rstrd)
  endif
  restart_io_rd( this%ufl_bnd, restart_handle, ierr )
  if ( ierr/=0 ) then
     ERROR('error reading restart data for spe_bound object.')
    call abort_program(p_err_rstrd)
  endif

end subroutine restart_read_spe_bound
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Process physical boundary conditions
!-------------------------------------------------------------------------------
subroutine process_phys_boundary( species, jay, dt, i_dim, bnd_idx, bnd_cross )

  implicit none

!       dummy variables

  class( t_species ), intent(inout)  :: species
  class( t_current ),      intent(inout) :: jay
  real(p_double), intent(in)        :: dt
  integer, intent(in)               :: i_dim    ! dimension being processed
  integer, intent(in)               :: bnd_idx  ! upper or lower boundary


  type( t_part_idx ), intent(inout) :: bnd_cross

  ! local variables
  integer   :: bnd_type

  ! executable statements

  if ( bnd_cross%nidx > 0 ) then

    bnd_type = species%bnd_con%type( bnd_idx, i_dim )

    if ( bnd_type == p_bc_absorbing ) then

      ! particles will be removed automatically by the update_boundary routine

#ifdef __HAS_TRACKS__

      ! if tracking remove particles from tracking list
      if (( species%diag%ndump_fac_tracks > 0 ) .and. (bnd_cross % nidx > 0))  then
        call missing_particles( species%diag%tracks, bnd_cross % idx, bnd_cross % nidx, species%tag )
      endif

#endif

    else

      call species % phys_boundary( jay, dt, i_dim, bnd_idx, &
                                    bnd_cross % idx, bnd_cross % nidx )

      ! keep all the particles
      bnd_cross % nidx = 0

    endif
  endif

end subroutine process_phys_boundary
!-------------------------------------------------------------------------------


#if 1

!---------------------------------------------------------------------------------------------------
! Check particles with index in cross buffer if outside boundary.
! V1 - Check only particles in the node_cross list, and particles placed at the end of the buffer
!
! The node_cross list must be created before calling this function, and min_npar has the minimum
! number of particles the species buffer has held during the update boundary process
!---------------------------------------------------------------------------------------------------
subroutine check_boundary( species, min_npar, i_dim, bnd_cross, node_cross )

  implicit none

  class( t_species ), intent(in) :: species
  integer, intent(in) :: min_npar, i_dim
  type( t_part_idx ), dimension(:), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross

  integer :: ixmax, i, k
  integer :: ncross, num_par_l, num_par_u
  integer, dimension(:), pointer :: cross, par_idx_l, par_idx_u
  integer :: par_idx_size

  cross => null(); par_idx_l => null(); par_idx_u => null();

  ! get maximum number of particles possibly crossing the node
  ncross = node_cross%nidx
  if ( min_npar < species%num_par ) ncross = ncross + ( species%num_par - min_npar )

  ! check if idx buffers are big enough
  if ( associated( bnd_cross( p_lower ) % idx ) ) then
    par_idx_size = size( bnd_cross( p_lower ) % idx )
  else
    par_idx_size = 0
  endif

  if ( par_idx_size < ncross ) then
    ! free existing buffers if needed
    call freemem( bnd_cross( p_lower ) % idx )
    call freemem( bnd_cross( p_upper ) % idx )

    ! give it a little head room
    call alloc( bnd_cross( p_lower ) % idx, (/ ncross + p_idx_block /) )
    call alloc( bnd_cross( p_upper ) % idx, (/ ncross + p_idx_block /) )
  endif

  num_par_l = 0
  num_par_u = 0
  par_idx_l => bnd_cross( p_lower ) % idx
  par_idx_u => bnd_cross( p_upper ) % idx

  ixmax = species%my_nx_p( 3, i_dim )

  ! Check particles in the node_cross index list
  cross  => node_cross%idx
  ncross =  node_cross%nidx
  do i = 1, ncross
  k = cross( i )

  ! Check if the particle buffer is now smaller
  if ( k > min_npar ) then
    node_cross%nidx = i - 1
    exit
  endif

  if ( species%ix(i_dim, k) < 1 ) then ! lower boundary
    num_par_l = num_par_l + 1
    par_idx_l(num_par_l) = k
  else if ( species%ix( i_dim , k ) > ixmax ) then ! upper boundary
    num_par_u = num_par_u + 1
    par_idx_u(num_par_u) = k
  endif
  enddo ! num_par

  ! Check particles that have arrived from a different node and have been placed at the end
  ! of the particle buffer
  do k = min_npar + 1, species%num_par

  if ( species%ix(i_dim, k) < 1 ) then ! lower boundary
    num_par_l = num_par_l + 1
    par_idx_l(num_par_l) = k
  else if ( species%ix( i_dim , k ) > ixmax ) then ! upper boundary
    num_par_u = num_par_u + 1
    par_idx_u(num_par_u) = k
  endif
  enddo

  ! Set values on bnd_cross index lists
  bnd_cross( p_lower ) % start = 1
  bnd_cross( p_upper ) % start = 1

  bnd_cross( p_lower ) % nidx = num_par_l
  bnd_cross( p_upper ) % nidx = num_par_u

end subroutine check_boundary
!---------------------------------------------------------------------------------------------------

#else

!---------------------------------------------------------------------------------------------------
! Check particles with index in cross buffer if outside boundary
! - Old version, checks all particles in particle buffer, keep around for debug purposes
!---------------------------------------------------------------------------------------------------
subroutine check_boundary( species, init_npar, i_dim, bnd_cross, node_cross )

  implicit none

  class( t_species ), intent(in) :: species
  integer, intent(in) :: init_npar, i_dim
  type( t_part_idx ), dimension(:), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross

  integer :: ixmax, i, k, i1
  integer :: ncross, num_par_l, num_par_u
  integer, dimension(:), pointer :: cross => null(), par_idx_l => null(), par_idx_u => null()
  integer :: par_idx_size


  ! get maximum number of particles possibly crossing the node
  ncross = species%num_par

  ! check if idx buffers are big enough
  if ( associated( bnd_cross( p_lower ) % idx ) ) then
    par_idx_size = size( bnd_cross( p_lower ) % idx )
  else
    par_idx_size = 0
  endif

  if ( par_idx_size < ncross ) then
  ! free existing buffers if needed
  call freemem( bnd_cross( p_lower ) % idx )
  call freemem( bnd_cross( p_upper ) % idx )

    ! give it a little head room
  call alloc( bnd_cross( p_lower ) % idx, (/ ncross + p_idx_block /) )
  call alloc( bnd_cross( p_upper ) % idx, (/ ncross + p_idx_block /) )
  endif

  num_par_l = 0
  num_par_u = 0
  par_idx_l => bnd_cross( p_lower ) % idx
  par_idx_u => bnd_cross( p_upper ) % idx

  ixmax = species%my_nx_p( 3, i_dim )

  do k = 1, species%num_par
  if ( species%ix(i_dim, k) < 1 ) then ! lower boundary
    num_par_l = num_par_l + 1
    par_idx_l(num_par_l) = k
  else if ( species%ix( i_dim , k ) > ixmax ) then ! upper boundary
    num_par_u = num_par_u + 1
    par_idx_u(num_par_u) = k
  endif
  enddo

  ! Set values on bnd_cross index lists
  bnd_cross( p_lower ) % start = 1
  bnd_cross( p_upper ) % start = 1

  bnd_cross( p_lower ) % nidx = num_par_l
  bnd_cross( p_upper ) % nidx = num_par_u

end subroutine check_boundary
!---------------------------------------------------------------------------------------------------

#endif

#if 1

!---------------------------------------------------------------------------------------------------
! Checks which particles have crossed the node boundaries. The result is stored in the module
! variable node_cross. This version supports OpenMP parallelism
!---------------------------------------------------------------------------------------------------
subroutine check_node_cross( this, n_threads, node_cross )

  implicit none

  class( t_species ), intent(in) :: this
  integer, intent(in):: n_threads
  type( t_part_idx ), intent(inout) :: node_cross

  integer, dimension( p_max_dim ) :: nx
  integer, dimension(:), pointer :: cross
  integer :: i
  integer :: ncross, cross_bsize

  integer :: chunk, tid, i0, i1
  integer, dimension( n_threads ) :: ncross_th

  ! Get node edges
  do i = 1, p_x_dim
    nx(i) = this%my_nx_p( 3, i )
  enddo

  if ( associated( node_cross % idx ) ) then
    cross_bsize = size( node_cross % idx )
  else
    cross_bsize = 0
  endif

  ! Grow buffer to be able to hold all particles in species
  if ( this%num_par > cross_bsize ) then
    call freemem( node_cross % idx )
    call alloc( node_cross % idx, (/ this%num_par /) )
  endif

  ! (* debug *)
  ! if ( associated( node_cross % idx ) ) node_cross % idx = -1


  ! Reset number of particles leaving the node
  cross => node_cross % idx

  ! Select range of particles per thread
  chunk = ( this%num_par + n_threads - 1 ) / n_threads

  !$omp parallel do private(tid, i0, i1, i, ncross)
  do tid = 1, n_threads

    i0    = (tid-1) * chunk + 1
    i1    = min( tid * chunk, this%num_par )

    ncross = 0

    select case ( p_x_dim )

    case (1)
      do i = i0, i1
        if ( (this%ix(1,i) < 1) .or. (this%ix(1,i) > nx(1)) ) then
          cross( i0 + ncross ) = i
          ncross = ncross + 1
        endif
      enddo

    case (2)
      do i = i0, i1
        if ( (this%ix(1,i) < 1) .or. (this%ix(1,i) > nx(1)) ) then
          cross( i0 + ncross ) = i
          ncross = ncross + 1
        else if ( (this%ix(2,i) < 1) .or. (this%ix(2,i) > nx(2)) ) then
          cross( i0 + ncross ) = i
          ncross = ncross + 1
        endif
      enddo

    case (3)
      do i = i0, i1
        if ( (this%ix(1,i) < 1) .or. (this%ix(1,i) > nx(1)) ) then
          cross( i0 + ncross ) = i
          ncross = ncross + 1
        else if ( (this%ix(2,i) < 1) .or. (this%ix(2,i) > nx(2)) ) then
          cross( i0 + ncross ) = i
          ncross = ncross + 1
        else if ( (this%ix(3,i) < 1) .or. (this%ix(3,i) > nx(3)) ) then
          cross( i0 + ncross ) = i
          ncross = ncross + 1
        endif
      enddo

    end select

    ncross_th( tid ) = ncross
  enddo

  ! compact cross idx
  ncross = ncross_th(1)
  do tid = 2, n_threads
    i0    = (tid-1) * chunk + 1
    do i = 1, ncross_th( tid )
      cross( ncross + i ) = cross( i0 + i - 1)
    enddo
    ncross = ncross + ncross_th( tid )
  enddo

  ! (*debug*)
  !do i = 1, ncross
  !  if ( cross(i) < 1 ) then
  !    write(0,*) '(* error *) Bad cross list'
  !    write(0,*) 'cross = ', cross(1:ncross)
  !    call abort_program()
  !  endif
  !enddo

  ! store the number of crossings
  node_cross % start = 1
  node_cross % nidx  = ncross

end subroutine check_node_cross
!---------------------------------------------------------------------------------------------------

#else

!---------------------------------------------------------------------------------------------------
! Checks which particles have crossed the node boundaries. The result is stored in the module
! variable node_cross. This version does not support OpenMP parallelism
!---------------------------------------------------------------------------------------------------

! This version is kept here for debug purposes only

! TODO: get rid of #warning.. it's not standard and if we really want to warn the user,
!       then move it up to the build/configure script
#if defined(_MSC_VER)
#else
#warning check_node_cross() using non OpenMP version
#endif


subroutine check_node_cross( this, n_threads, node_cross )

  implicit none

  class( t_species ), intent(in) :: this
  integer, intent(in):: n_threads ! This is not used in this version

  integer, dimension( p_x_dim ) :: nx
  integer, dimension(:), pointer :: cross
  type( t_part_idx ), intent(inout) :: node_cross
  integer :: i
  integer :: ncross, cross_bsize

  ! Get node edges
  do i = 1, p_x_dim
    nx(i) = this%my_nx_p( 3, i )
  enddo

  if ( associated( node_cross % idx ) ) then
    cross_bsize = size( node_cross % idx )
  else
    cross_bsize = 0
  endif

  ! Grow buffer to be able to hold all particles in species
  if ( this%num_par > cross_bsize ) then
    call freemem( node_cross % idx )
    call alloc( node_cross % idx, (/ this%num_par /) )
  endif

  ! Reset number of particles leaving the node
  cross => node_cross % idx

  ncross = 0

  select case ( p_x_dim )

  case (1)
    do i = 1, this % num_par
      if ( (this%ix(1,i) < 1) .or. (this%ix(1,i) > nx(1)) ) then
        ncross = ncross + 1
        cross( ncross ) = i
      endif
    enddo

  case (2)
    do i = 1, this % num_par
      if ( (this%ix(1,i) < 1) .or. (this%ix(1,i) > nx(1)) ) then
        ncross = ncross + 1
        cross( ncross ) = i
      else if ( (this%ix(2,i) < 1) .or. (this%ix(2,i) > nx(2)) ) then
        ncross = ncross + 1
        cross( ncross ) = i
      endif
    enddo

  case (3)
    do i = 1, this % num_par
      if ( (this%ix(1,i) < 1) .or. (this%ix(1,i) > nx(1)) ) then
        ncross = ncross + 1
        cross( ncross ) = i
      else if ( (this%ix(2,i) < 1) .or. (this%ix(2,i) > nx(2)) ) then
        ncross = ncross + 1
        cross( ncross ) = i
      else if ( (this%ix(3,i) < 1) .or. (this%ix(3,i) > nx(3)) ) then
        ncross = ncross + 1
        cross( ncross ) = i
      endif
    enddo

  end select

  ! store the number of crossings
  node_cross % start = 1
  node_cross % nidx  = ncross

end subroutine check_node_cross
!---------------------------------------------------------------------------------------------------


#endif


!---------------------------------------------------
function type_spe_bound( this )
!---------------------------------------------------
!       gives the boundary type variable for this boundary condition
!---------------------------------------------------

  implicit none

!       dummy variables

  integer, dimension(2,p_x_dim)  ::  type_spe_bound

  class ( t_spe_bound ),intent(in) :: this

!       local variables - none

!       executable statements

  type_spe_bound(:,1:p_x_dim) = this%type(:,1:p_x_dim)

end function type_spe_bound
!---------------------------------------------------

!---------------------------------------------------
function type_bnd_dim_spe_bound( this, bnd, dim )
!---------------------------------------------------
!       gives the boundary type variable for this boundary condition
!---------------------------------------------------

  implicit none

!       dummy variables

  integer  ::  type_bnd_dim_spe_bound

  class ( t_spe_bound ),intent(in) :: this
  integer, intent(in) :: bnd, dim

!       local variables - none

!       executable statements

  type_bnd_dim_spe_bound = this%type( bnd, dim )

end function type_bnd_dim_spe_bound
!---------------------------------------------------




!-------------------------------------------------------------------------------
! Set the momenta of particles coming from thermal bath boundary.
!-------------------------------------------------------------------------------
subroutine thermal_bath( this, i_dim, bnd_idx, u, npar, tmp_ix )

  implicit none

  class ( t_spe_bound ),intent(in) :: this
  integer, intent(in) :: i_dim, bnd_idx
  real(p_k_part), intent(inout), dimension(:,:) :: u
  integer, intent(in) :: npar
  integer, pointer, dimension(:,:) :: tmp_ix

  integer :: k
  real(p_k_part) :: usign, u1, u2, u3, u_tmp, upar, uperp, theta
  real(p_k_part) :: a, b, gth

  real(p_double), parameter :: pi      = 3.14159265358979323846264_p_double


  ! Make sure the component normal to the boundary points inward
  select case( bnd_idx )
  case( p_lower )
    usign = 1
  case( p_upper )
    usign = -1
  case default
    usign = 0
  end select

  u1 = this%uth_bnd( 1, bnd_idx, i_dim )
  u2 = this%uth_bnd( 2, bnd_idx, i_dim )
  u3 = this%uth_bnd( 3, bnd_idx, i_dim )

  select case ( this%thermal_type(bnd_idx, i_dim) )

    case( p_thermal )

      u_tmp = usign * this%uth_bnd(i_dim,bnd_idx,i_dim)

      do k = 1, npar
        u(1,k) = u1 * real( rng % genrand_gaussian( tmp_ix(:, k)), p_k_part )
        u(2,k) = u2 * real( rng % genrand_gaussian( tmp_ix(:, k)), p_k_part )
        u(3,k) = u3 * real( rng % genrand_gaussian( tmp_ix(:, k)), p_k_part )

        ! This is required to maintain a null flux across the boundary
        ! Higer momentum particles are more likely to exit the boundary, so re-injecting
        ! from a normal gaussian distribution will generally yield lower energy particles
        ! and cool down the plasma
        u(i_dim, k) = u_tmp * real(rng % genrand_half_max(tmp_ix(:, k)), p_k_part) 
      enddo

    case ( p_waterbag )
      do k = 1, npar
        u(1,k) = u1 * real( rng % genrand_real3( tmp_ix(:, k)) - 0.5, p_k_part )
        u(2,k) = u2 * real( rng % genrand_real3( tmp_ix(:, k)) - 0.5, p_k_part )
        u(3,k) = u3 * real( rng % genrand_real3( tmp_ix(:, k)) - 0.5, p_k_part )

        u(i_dim, k) = usign * abs( u(i_dim, k) )
      enddo

    case ( p_random_dir )

      u_tmp = sqrt( u1**2 + u2**2 + u3**2 )

      do k = 1, npar
    ! get p3
    upar = real(2 * rng % genrand_real3( tmp_ix(:, k)) - 1.0_p_double, p_k_part)

    ! get p2 and p1
    uperp = sqrt(1 - upar)*u_tmp
    theta = real(2 * pi * rng % genrand_real3( tmp_ix(:, k)), p_k_part)
    u(1,k) = cos(theta) * uperp
    u(2,k) = sin(theta) * uperp
    u(3,k) = upar       * u_tmp

        u(i_dim, k) = usign * abs( u(i_dim, k) )
      enddo

  end select

  ! Add fluid momenta if required
  u1 = this%ufl_bnd( 1, bnd_idx, i_dim )
  u2 = this%ufl_bnd( 2, bnd_idx, i_dim )
  u3 = this%ufl_bnd( 3, bnd_idx, i_dim )

  if ( u1 /= 0.0 .and. u2 /= 0.0 .and. u3 /= 0.0 ) then
    a = 1.0_p_k_part / (1.0_p_k_part + sqrt( 1 + u1**2 + u2**2 + u3**2 ) )

    do k = 1, npar
      gth = sqrt( 1 + u(1,k)**2 + u(2,k)**2 + u(3,k)**2 )
      b = gth + a * (u1*u(1,k) + u2*u(2,k) + u3*u(3,k))
      u(1,k) = u(1,k) + b * u1
      u(2,k) = u(2,k) + b * u2
      u(3,k) = u(3,k) + b * u3
    enddo
  endif

end subroutine thermal_bath
!-------------------------------------------------------------------------------

function ntrim(x)

  implicit none

  real(p_k_part), intent(in) :: x
  integer :: ntrim, a, b

  if ( x < -.5 ) then
  a = -1
  else
    a = 0
  endif

  if ( x >= .5 ) then
  b = +1
  else
    b = 0
  endif

  ntrim = a+b

end function ntrim

end module m_species_boundary


!---------------------------------------------------------------------------------------------------
! Update species boundary conditions
!---------------------------------------------------------------------------------------------------
subroutine update_boundary_species( spec, jay, no_co, dt, bnd_cross, node_cross, send_msg, recv_msg )

  use m_species_define
  use m_species_boundary
  use m_species_comm
  use m_vdf_define
  use m_node_conf
  use m_parameters
  use m_current_define

  implicit none

  ! species
  class( t_species ), intent(inout) :: spec
  ! current
  class( t_current ),     intent(inout) :: jay
  ! node configuration
  class( t_node_conf ), intent(in) :: no_co
  real(p_double), intent(in) :: dt
  ! boundary information
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i
  integer :: my_node, comms, min_npar


  ! Get list of particles crossing the node boundaries. The result is stored in the module
  ! variable node_cross
  call check_node_cross( spec, n_threads(no_co), node_cross )

  my_node = no_co%my_aid()
  min_npar = spec%num_par

  do i = 1, p_x_dim

    if ( spec%num_par < min_npar ) min_npar = spec%num_par


    ! communication with other nodes
    comms = 0
    if ( my_node /= neighbor(no_co,p_lower,i) ) then
      if ( neighbor( no_co, p_lower, i ) > 0 ) comms = comms + p_lower
      if ( neighbor( no_co, p_upper, i ) > 0 ) comms = comms + p_upper
    endif

    ! Post receives for number of particles and particle data
    select case ( comms )
      case (0)
         ! no comms
      case (p_lower)
         ! communication with lower node only
         call irecv_1( spec, i, p_lower, no_co, recv_msg )
      case (p_upper)
         ! communication with upper node only
         call irecv_1( spec, i, p_upper, no_co, recv_msg )
      case (p_lower + p_upper)
         ! communication with both nodes
         call irecv_1( spec, i, p_lower, no_co, recv_msg )
         call irecv_1( spec, i, p_upper, no_co, recv_msg )
    end select

    ! check which particles crossed the node boundary and store their indexes
    call check_boundary( spec, min_npar, i, bnd_cross, node_cross )

    ! Process particles crossing nodes
    select case ( comms )
      case (0)
         ! no comms
         call process_phys_boundary( spec, jay, spec%dt, i, p_lower, bnd_cross(p_lower) )
         call process_phys_boundary( spec, jay, spec%dt, i, p_upper, bnd_cross(p_upper) )

      case (p_lower)
         ! Pack data and send it to lower node
         call isend_1( spec, i, p_lower, no_co, bnd_cross(p_lower), send_msg )

         ! process physical boundary
         call process_phys_boundary( spec, jay, spec%dt, i, p_upper, bnd_cross(p_upper) )

         ! Process step2 of communications
         call irecv_2( p_lower, recv_msg )
         call isend_2( p_lower, send_msg )

         ! wait recv step2 and unpack data
         call unpack( spec, p_lower, bnd_cross, recv_msg )

      case (p_upper)
         ! communication with upper node only
         call isend_1( spec, i, p_upper, no_co, bnd_cross( p_upper ), send_msg )

         ! process physical boundary
         call process_phys_boundary( spec, jay, spec%dt, i, p_lower, bnd_cross( p_lower )  )

         ! Process step2 of communications
         call irecv_2( p_upper, recv_msg )
         call isend_2( p_upper, send_msg )

         ! wait recv step2 and unpack data
         call unpack( spec, p_upper, bnd_cross, recv_msg )

      case (p_lower+p_upper)
         ! communication with both nodes

         ! Pack data and send it to nodes. Messages are posted in opposite order of receives.
         call isend_1( spec, i, p_upper, no_co, bnd_cross( p_upper ), send_msg )
         call isend_1( spec, i, p_lower, no_co, bnd_cross( p_lower ), send_msg )

         ! Process step2 of communications (this handles the 4 messages)
         call irecv_2( p_lower, recv_msg )
         call irecv_2( p_upper, recv_msg )

         call isend_2( p_upper, send_msg )
         call isend_2( p_lower, send_msg )

         ! wait recv step2 and unpack data
         call unpack( spec, p_lower + p_upper, bnd_cross, recv_msg )

    end select

    ! remove particles whose index is still in par_idx_l or par_idx_u and pack the species buffer
    call remove_particles( spec, bnd_cross )

    ! Wait for any leftover send messages to complete
    if ( comms > 0 ) call wait_send( comms, send_msg )

  enddo

#ifdef __DEBUG__
  call spec % validate( 'after update boundary')
  call no_co % barrier()
#endif

end subroutine update_boundary_species
!---------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Process physical boundary conditions
!-----------------------------------------------------------------------------------------
subroutine phys_boundary_species( this, current, dt, i_dim, bnd_idx, par_idx, npar )

  use m_species_define, only : t_species
  use m_current_define, only : t_current
  use m_parameters
  use m_vdf_define, only : t_vdf
  use m_species_current, only : tmp_xold, tmp_ix, tmp_q, tmp_p, tmp_rg, tmp_xnew, tmp_dxi
  use m_random, only : rng
  use m_species_boundary, only : thermal_bath
  use m_species_current, only : deposit_current_1d, deposit_current_2d, deposit_current_3d

  implicit none

  class( t_species ), intent(inout) :: this
  class( t_current ), intent(inout) :: current
  real(p_double), intent(in) :: dt

  integer, intent(in) :: i_dim, bnd_idx
  integer, dimension(:), intent(in) :: par_idx
  integer, intent(in) :: npar

  ! local variables
  type(t_vdf), pointer :: jay
  real(p_k_part), dimension(p_x_dim) :: dt_dx
  integer :: k, i, l, nbuf, bc_type

  integer :: shift, ipos
  real(p_k_part) :: delta, vperp, xpos

  jay => current%pf(1)

  ! dt / dx
  dt_dx(1:p_x_dim) = real( dt/this%dx(1:p_x_dim), p_k_part )

  ! Boundary condition type
  bc_type = this%bnd_con%type( bnd_idx, i_dim )

  ! On a moving window run, the upper boundary behaves as specular
  if ( bnd_idx == p_upper .and. bc_type == p_bc_move_c ) bc_type = p_bc_specular

  select case ( bc_type )

  case( p_bc_periodic )  ! Single node periodic boundaries

    if ( bnd_idx == p_lower ) then
      shift = + this%my_nx_p( 3, i_dim )
    else
      shift = - this%my_nx_p( 3, i_dim )
    endif

    do  k = 1, npar
      this%ix( i_dim , par_idx(k) ) = this%ix( i_dim , par_idx(k) ) + shift
    enddo


  case ( p_bc_specular ) ! specular reflection

    if ( bnd_idx == p_lower ) then
      ipos = 1
    else
      ipos = this%my_nx_p( 3, i_dim )
    endif

    do k = 1, npar
      ! flip momentum
      this%p(i_dim,par_idx(k) )  = -this%p(i_dim,par_idx(k) )

      ! position after reflection
      this%ix(i_dim,par_idx(k) ) = ipos
      this%x(i_dim,par_idx(k) )  = - this%x(i_dim,par_idx(k))
    enddo

    ! There is no current being deposited from the reflected particles here.
    ! If the user is using some form of conducting boundaries for the fields then
    ! the electric current object will reflect the current deposited out of the
    ! simulation box back into box

  case ( p_bc_thermal ) ! thermal bath

    ! Get position of boundary
    if ( bnd_idx == p_lower ) then
      ipos = 1
      xpos = -0.5_p_k_part
    else
      ipos = this%my_nx_p( 3, i_dim )
      ! to inject inside the current grid cell, xpos should be samller than 0.5_p_k_part
      xpos = nearest(+0.5_p_k_part, -1.0_p_k_part)
    endif

    do i = 1, npar, p_cache_size

      if( i + p_cache_size > npar ) then
        nbuf = npar - i + 1
      else
        nbuf = p_cache_size
      endif

      do  k= 1, nbuf
        this%ix(i_dim,par_idx(i+k-1) ) = ipos
        tmp_ix(:,k) = this%ix(:,par_idx(i+k-1) )
      end do

      call thermal_bath( this%bnd_con, i_dim, bnd_idx, tmp_p, nbuf, tmp_ix )


      do  k= 1, nbuf

        ! Store new momenta
        this%p (1,par_idx(i+k-1) ) = tmp_p(1,k)
        this%p (2,par_idx(i+k-1) ) = tmp_p(2,k)
        this%p (3,par_idx(i+k-1) ) = tmp_p(3,k)

        tmp_rg(k) = 1.0 / sqrt( 1.0 + tmp_p(1,k)**2+ tmp_p(2,k)**2+ tmp_p(3,k)**2 )

        call rng % harvest_real3( delta, this%ix(:,par_idx(i+k-1)))

        vperp = tmp_rg(k) * tmp_p(i_dim,k) * dt_dx(i_dim)

        this%x (i_dim,par_idx(i+k-1) ) = xpos + vperp * delta

        ! store values for current deposition, change the sign of charge and
        ! reverse the in-plane momentum
        tmp_q( k ) = - this%q( par_idx(i+k-1) )
        ! we don't really need to change the sign of tmp_p because
        ! the in-plane momentum is not used in deposit_current_{1,2,3}d directly.
        ! The particle motion is reversed in rev_current below
        !do l = p_p_dim-p_x_dim, p_p_dim
        !  tmp_p(l,k) = -tmp_p(l,k)
        !enddo

        do l = 1, p_x_dim
          ! the following line is commented out since we do this
          !   operation already above (right before calling 'thermal_bath')
          !   since the IX data is sometimes needed for random number generation.
          !tmp_ix( l, k )   = this%ix(l, par_idx(i+k-1) )
          tmp_xold( l, k ) = this%x(l, par_idx(i+k-1) )
        enddo

      enddo

      ! deposit current as if particles were coming from another node
      call rev_current( this, jay, tmp_xold, tmp_ix, tmp_q, tmp_p, tmp_rg, nbuf, dt_dx, dt )
    enddo

  case default

    ERROR("invalid particle bc")
    print *, 'type(', bnd_idx,', ', i_dim, ' ) = ', bc_type

  end select


  contains

  subroutine rev_current( this, jay, xold, ix, q, p, rg, nbuf, dt_dx, dt )

    implicit none

    class( t_species ), intent(inout)  :: this
    type( t_vdf ),      intent(inout) :: jay
    real(p_k_part), dimension(:,:), intent(inout) :: xold
    integer, dimension(:,:), intent(inout) :: ix
    real(p_k_part), dimension( : ), intent(inout) :: q
    real(p_k_part), dimension(:,:), intent(inout) :: p
    real(p_k_part), dimension( : ), intent(inout) :: rg
    integer, intent(in) :: nbuf
    real(p_k_part), dimension(:), intent(in)        :: dt_dx
    real(p_double), intent(in)        :: dt

    integer :: k,l

#if 1

    ! get previous position
    do k = 1, nbuf
      do l = 1, p_x_dim
        tmp_xnew(l,k) = xold(l,k) - p(l,k)*rg(k)*dt_dx(l)
        tmp_dxi(l,k)  = ntrim( tmp_xnew(l,k) )
      enddo
    enddo

    select case (p_x_dim)
    case (1)
      call deposit_current_1d( this, jay, tmp_dxi, tmp_xnew, ix, xold, q, &
        rg, p, nbuf, dt)
    case (2)
      call deposit_current_2d( this, jay, tmp_dxi, tmp_xnew, ix, xold, q, &
        rg, p, nbuf, dt)
    case (3)
      call deposit_current_3d( this, jay, tmp_dxi, tmp_xnew, ix, xold, q, &
        nbuf, dt)
    end select

#else

    SCR_ROOT('(*warning*) current deposition by particle boundary routines')
    SCR_ROOT('(*warning*) has been turned off')
    return

#endif


  end subroutine rev_current

  function ntrim(x)

    implicit none

    real(p_k_part), intent(in) :: x
    integer :: ntrim, a, b

    if ( x < -.5 ) then
    a = -1
    else
      a = 0
    endif

    if ( x >= .5 ) then
    b = +1
    else
      b = 0
    endif

    ntrim = a+b

  end function ntrim

end subroutine phys_boundary_species
!-----------------------------------------------------------------------------------------

