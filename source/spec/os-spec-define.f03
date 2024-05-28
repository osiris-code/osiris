!-----------------------------------------------------------------------------------------

! Species definition module
!
! This file contains the class definition for the following classes:
!
!  t_spe_bound
!  t_piston
!  t_species
!-----------------------------------------------------------------------------------------


#include "os-preprocess.fpp"
#include "os-config.h"

module m_species_define

#include "memory/memory.h"

use m_system
use m_parameters

use m_fparser


use m_time_avg

use m_diagnostic_utilities
use m_diagfile, only : p_no_offset

use m_time_step,     only : t_time_step
use m_space,         only : t_space
use m_grid_define,   only : t_grid, t_msg_patt
use m_node_conf,     only : t_node_conf
use m_restart,       only : t_restart_handle
use m_input_file,    only : t_input_file
use m_emf_define,    only : t_emf
use m_vdf_define,    only : t_vdf, t_vdf_report
use m_vdf_comm,      only : t_vdf_msg
use m_current_define,only : t_current

!-----------------------------------------------------------------------------------------

implicit none

private

! string to id restart data
character(len=*), parameter :: p_spec_rst_id = "species rst data - 0x000F"
public :: p_spec_rst_id

! minimal block size to be used when growing particle buffers
integer, parameter :: p_spec_buf_block = 131072 ! 128k

! species position definition
integer, parameter :: p_cell_low  = 1
integer, parameter :: p_cell_near = 2

!-----------------------------------------------------------------------------------------
! Momentum distribution definitions
!-----------------------------------------------------------------------------------------

! constants for velocity distribution
integer, parameter :: p_none         = -1  ! all particles initialized at rest
integer, parameter :: p_thermal      = 1    ! normal distribution of momentum
integer, parameter :: p_random_dir   = 2    ! fixed velocity, random direction
integer, parameter :: p_waterbag     = 3    ! waterbag
integer, parameter :: p_relmax       = 4    ! Relativistic maxwellian
integer, parameter :: p_waterbag_rel = 5    ! Relativistic waterbag
integer, parameter :: p_relmax_boosted = 6  ! Relativistic juettner with a boost
integer, parameter :: p_uth_expr     = 7    ! Math function input for the momentum

integer, parameter :: p_half_max     = 10   ! this is only used for thermal boundaries


integer, parameter :: p_uniform = 0    ! uniform distribution
integer, parameter :: p_spatial = 1    ! spatially dependent distribution

#ifdef __HAS_SPIN__

! constants for spin distribution
integer, parameter :: p_isotropic   = 0  ! isotropic spin distribution
integer, parameter :: p_fixed_dir   = 1  ! fixed spin direction
integer, parameter :: p_spherical_gauss = 2
public :: p_isotropic, p_fixed_dir, p_spherical_gauss

#endif

! force recalculation of energy after particle push
integer, parameter :: p_ene_recalc = -1

! constant for default relmax distribution cut-off (parameter x vdist_T)
real, parameter :: p_relmax_umax = 20.0

! constants for momentum and charge increase for beam loading
integer, parameter :: incr_linear = 1 ! increase linear within n steps, default
integer, parameter :: incr_regressive = 2 ! increase regressive  (first strong, than weaker)
!integer, parameter :: incr_progressive = 3 ! increase progressive (first weak, than stronger)

! Maximum buffer size for step 1 communication (in bytes)

#if defined( __MIC__ )

#warning Setting smaller particle comm buffer for MIC platform

! integer, parameter :: p_max_buffer1_size = 65536
 integer, parameter :: p_max_buffer1_size = 131072

#elif defined( __KNL__ )

#warning Setting smaller particle comm buffer for KNL platform

 integer, parameter :: p_max_buffer1_size = 262144

#else

integer, parameter :: p_max_buffer1_size = 4194304
!integer, parameter :: p_max_buffer1_size = 1048576
!integer, parameter :: p_max_buffer1_size = 128

#endif

public :: p_none, p_ene_recalc
public :: p_thermal, p_random_dir, p_waterbag, p_relmax
public :: p_relmax_boosted, p_waterbag_rel, p_uth_expr
public :: p_half_max, p_relmax_umax
public :: p_uniform, p_spatial
public :: incr_linear, incr_regressive !,incr_progressive
public :: p_max_buffer1_size

type t_udist

  integer :: uth_type                             ! distribution type
  real(p_k_part), dimension(p_p_dim) :: uth       ! thermo-momentum

  ! Parameters for relativistic maxwellian temperature
  real(p_k_part)                     :: relmax_T         ! relmax temperature
  real(p_k_part)                     :: relmax_umax      ! relmax max temperature

  real(p_k_part)                     :: relmax_boost     ! \beta * \gamma of boost
  integer                            :: relmax_boost_dir ! direction of boost

  logical :: use_spatial_uth
  type(t_fparser), dimension(p_p_dim) :: spatial_uth

  integer :: ufl_type
  real(p_k_part), dimension(p_p_dim) :: ufl       ! fluid-momentum
  type(t_fparser), dimension(p_p_dim) :: spatial_ufl

  type(t_fparser), dimension(p_p_dim) :: math_func_uth

  logical :: use_classical_uadd
  logical :: use_particle_uacc

  ! Gradual acceleration
  integer :: n_accelerate = -1
  integer :: n_accelerate_type
  
  ! initialize with a free streaming pusher with increasing charge
  ! within the first n steps
  integer :: n_q_incr = -1
  integer :: n_q_incr_type

  real(p_k_part), dimension(p_p_dim) :: umin
  real(p_k_part), dimension(p_p_dim) :: umax
  logical, dimension(p_p_dim) :: math_func_use_thermal

end type t_udist

#ifdef __HAS_SPIN__

!-----------------------------------------------------------------------------------------
! Species spin distribution classes definitions
!-----------------------------------------------------------------------------------------

type t_sdist

  integer :: sdist_type                            ! distribution type

  real(p_k_part), dimension(:,:), pointer :: sdir => null() ! spin orientation, (/ n_dirs, p_s_dim /)
  real(p_k_part), dimension(:), pointer :: pdf => null()  ! probability dist. function for each spin direction
  integer :: n_dirs
  real(p_k_part), dimension(p_s_dim) :: gauss_mu
  real(p_k_part) :: gauss_theta

end type t_sdist

public :: t_sdist

#endif

!-----------------------------------------------------------------------------------------
! Species diagnostic classes definitions
!-----------------------------------------------------------------------------------------

#include "diagnostics/os-spec-diag-def.f03"

!-----------------------------------------------------------------------------------------
! Species Boundary Condition Class
!-----------------------------------------------------------------------------------------

type :: t_spe_bound

  ! Boundary condition type
  integer, dimension(2,p_max_dim) :: type

  ! Type of momenta distribution for each wall
  integer, dimension(2,p_max_dim) :: thermal_type

  ! thermal and fluid momenta of the thermal bath BC
  real(p_k_part), dimension( p_p_dim,2,p_x_dim ) :: uth_bnd
  real(p_k_part), dimension( p_p_dim,2,p_x_dim ) :: ufl_bnd

contains

  procedure, nopass :: init_tmp_buf_current => init_tmp_buf_current_spe_bnd
  procedure, nopass :: cleanup_tmp_buf_current => cleanup_tmp_buf_current_spe_bnd

end type t_spe_bound
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! species message object
!-----------------------------------------------------------------------------------------

type t_spec_msg

  ! communicator
  integer :: comm = MPI_COMM_NULL

  ! message request id
  integer :: request = MPI_REQUEST_NULL

  ! node communicating to/from
  integer :: node = -1

  ! message tag ( optional )
  integer :: tag = 0

  ! number of particles to send/receive
  integer :: n_part = 0

  ! size of packed data for single particle
  integer :: particle_size

  ! buffer sizes for communication
  integer :: max_buffer1_size = p_max_buffer1_size
  integer :: n_part1 = 0, n_part2 = 0
  integer :: size1 = 0, size2 = 0

  ! buffers for messages
  integer(p_byte), dimension(:), pointer :: buffer1 => null()
  integer(p_byte), dimension(:), pointer :: buffer2 => null()

  ! attributes for tile module
  ! step, dim, and num_par in tiled boundary update
  integer :: step, dim

  ! number of tils grouped into first message buffer
  integer :: n_tils1 = 0

  ! add potential to have a linked list of these
  type(t_spec_msg), pointer :: next => null()

end type t_spec_msg
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Piston Class
!-----------------------------------------------------------------------------------------

type :: t_piston
  ! parameters from inputdeck (processed)
  integer          :: dim        ! dimension in which piston moves
  integer          :: updown     ! edge from whitch piston starts
  real(p_k_part)  :: u          ! proper velocity ( gamma v) of piston
  real(p_k_part)  :: start_pos  ! position from which piston starts
  real(p_k_part)  :: start_time ! time when piston launched of edge
  real(p_k_part)  :: stop_time  ! time when piston disapears
  real(p_k_part)  :: opacity_factor  ! uniform opacity or factor for profile
  integer          :: profile_type    ! type of profile
  type(t_fparser)  :: profile_parser  ! profile function (A, B: transverse variables

  !auxiliary variabes
  real(p_k_part)  :: v          ! v
  real(p_k_part)  :: squ        ! u^2 = (gamma v)^2
  real(p_k_part)  :: g          ! gamma
  real(p_k_part)  :: sqg        ! gamma^2
  ! real(p_k_part)  :: sqv        ! v^2

  real(p_k_part)  :: pos_after  ! position of piston at t
  real(p_k_part)  :: pos_before ! position of piston at t-dt

  logical          :: inbox      ! true if piston is in this node
end type t_piston


!-----------------------------------------------------------------------------------------
! Species particle index class
!  - Stores a list of particle indexes
!-----------------------------------------------------------------------------------------

type t_part_idx

  ! indexes of particles
  integer, dimension(:), pointer :: idx => null()

  ! size of buffer
  integer :: buf_size = 0

  ! Total number of particles
  integer :: nidx    = 0

  ! starting position
  integer :: start   = 1
end type t_part_idx
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Species Class
!-----------------------------------------------------------------------------------------

! max. length of species name
integer, parameter :: p_max_spname_len = 64

type :: t_species
  
  ! current iteration that is currently held by this
  !   t_current object. Used with CUDA to keep alert
  !   if jay needs to be copied from the GPU to be
  !   output for diagnostics.
  integer :: n_current_spec = 0

  ! species name
  character(len = p_max_spname_len) :: name

  ! species number - internal id
  integer :: sp_id

  ! particles per cell
  integer, dimension(p_max_dim) :: num_par_x

  ! total number of particle in one direction. this can be really big.
  integer(p_int64), dimension(p_max_dim) :: tot_par_x

  ! simulation box dimensions

  ! This is shifted from the simulation box dimensions by + 0.5*dx to simplify global particle
  ! position calculations, which are only used for injection (density profile) and diagnostics

  ! The only exception is for radial cylindrical coordinates with even order interpolation, where the
  ! algorithm also uses it during the push and requires that this is set to the same value value as
  ! the global value

  real(p_double), dimension(2, p_max_dim) :: g_box

  ! simulation parameters
  real(p_double)                          :: dt             ! simulation time step
  real(p_double), dimension(p_max_dim)    :: dx             ! cell size

  integer, dimension(p_max_dim)           :: g_nx           ! simulation box grid size
  real(p_double), dimension(p_x_dim)      :: total_xmoved   ! motion of the simulation box
  integer, dimension(p_x_dim)             :: move_num = 0   ! number of dx moved for current timestep
  integer, dimension(3, p_max_dim)        :: my_nx_p        ! local grid position in global grid
  integer :: coordinates                    ! coordinate system in use
                                  ! cartesian / cylindrical

  ! local node position in global grid as a single int
  ! (used for particle tagging)
  integer :: ngp_id

  ! particle buffer size
  integer :: num_par_max

  ! number of particles in buffer
  integer :: num_par

  ! number of particles saved in update_boundary for tiles module
  integer :: num_par_save

  ! number of particles that have been created in this node
  integer :: num_created = 0

  !  mass to charge ratio
  real(p_k_part)                     :: rqm ! [me/e]
  real(p_k_part)                     :: rq_real, q_real ! charge of a real particle [1/e], [e]
  real(p_k_part)                     :: m_real ! mass of a real particle [me]

  ! essential particle data position, momentum, charge, and spin
  real(p_k_part), dimension(:,:), pointer :: x => null()
  real(p_k_part), dimension(:,:), pointer :: p => null()
  real(p_k_part), dimension(:),   pointer :: q => null()

#ifdef __HAS_SPIN__

  real(p_k_part), dimension(:,:), pointer :: s => null()

  ! spin distribution parameters
  type( t_sdist ) :: sdist

  ! anomalous magnetic moment
  real(p_double) :: anom_mag_moment

  ! Pointer to dsdt function used by some pushers
  procedure(dsdt_species), pointer :: dsdt => null()

#endif

  ! define particle positions related to current cell or simulation box
  integer                          :: pos_type
  integer, dimension(:,:), pointer :: ix => null()

  ! type of current deposition / field interpolation
  integer :: interpolation
  logical :: grid_center

  ! Source
  class( t_psource ), pointer :: source => null()

  ! velocity distribution parameters
  type( t_udist ) :: udist

  ! initialize fields from initial density / momentum distribution
  logical :: init_fields = .false.

  ! initialize particle distribution (only false for tile dlb)
  logical :: init_part = .true.

  ! particle tags
  logical :: add_tag = .false.
  integer, dimension(:,:), pointer :: tag => null()


  ! Free streaming species (i.e. constant velocity). dgam can still be used
  logical :: free_stream = .false.

  ! boundary conditions for this species
  class( t_spe_bound ), pointer :: bnd_con => null()

  ! diagnostic for this species
  class( t_diag_species ), pointer  :: diag  => null()

  ! number of timesteps between sorting of particles
  ! n_sort = 0 turns sorting off
  integer :: n_sort

  ! switch wether to collide this species
  logical :: if_collide
  logical :: if_like_collide

  ! local E and B fields for n_push > 1
  ! averaged for n_push times
  type(t_vdf), pointer :: E => NULL()
  type(t_vdf), pointer :: B => NULL()


  ! Numerical piston data
  integer                                 :: num_pistons
  type( t_piston ), dimension(:), pointer :: pistons  => NULL()

  ! time-centered total energy diagnostic (this is always done in double precision)
  ! this is set to an array to allow calculations in OpenMP parallel runs
  real(p_double), dimension(:), pointer :: energy => null()

  ! Push options
  integer             :: push_type          ! type of push (standard/simd)
  real(p_double)      :: push_start_time    ! delayed push

  ! Iteration tolerance of exact pusher
  real(p_double) :: iter_tol = 1.0d-3

  ! radiation reaction parameters for exact pusher
  real(p_k_part) :: k_rr ! leading coefficient of RR
  logical :: rad_react = .false.

  ! Pointer to dudt function used by some pushers
  procedure(dudt_species), pointer :: dudt => null()

  ! Pointers to current deposition routines (not used by simd code)
  procedure(dep_current_1d), nopass, pointer :: dep_current_1d => null()
  procedure(dep_current_2d), nopass, pointer :: dep_current_2d => null()
  procedure(dep_current_3d), nopass, pointer :: dep_current_3d => null()

  ! Next species (for storing groups of species as linked lists)
  class(t_species), pointer :: next => null()

contains

  procedure :: allocate_objs        => allocate_objs_spec
  procedure :: get_n_x_dims         => get_n_x_dims_spec
  procedure :: init_buffer          => init_buffer_spec
  procedure :: grow_buffer          => grow_buffer_spec

  procedure :: init                 => init_species
  procedure :: set_dudt             => set_dudt_spec
  procedure :: push                 => push_species
  procedure :: get_emf              => get_emf_spec 
  procedure :: read_input           => read_input_species
  procedure :: get_diag_buffer_size => get_diag_buffer_size_spec
  procedure :: report               => report_species
  procedure :: report_energy        => report_energy_species
  procedure :: fill_data            => fill_data_species
  procedure :: get_phasespace_axis  => get_phasespace_axis_spec
  procedure :: reshape              => reshape_spec

  procedure :: get_quant            => get_quant_spec
  procedure :: deposit_density      => deposit_density_spec
  procedure :: enumerate_quants     => enumerate_quants_spec
  procedure :: update_boundary      => update_boundary_species
  procedure :: phys_boundary        => phys_boundary_species
  procedure :: cleanup              => cleanup_species

  procedure :: inject                         => species_inject
  procedure :: create_particle_single_cell_p  => create_particle_single_cell_p
  procedure :: create_particle_single_cell    => create_particle_single_cell

  procedure :: list_algorithm       => list_algorithm_spec
  procedure :: position_single_4
  procedure :: position_single_8
  procedure :: position_range_comp_4
  procedure :: position_range_comp_8
  procedure :: position_idx_comp_4
  procedure :: position_idx_comp_8
  procedure :: position_ref_box_idx_comp_4
  procedure :: position_ref_box_idx_comp_8
  procedure :: position_range_full_4
  procedure :: position_range_full_8
  generic   :: get_position => position_single_4, position_single_8, &
                               position_range_comp_4, position_range_comp_8, &
                               position_idx_comp_4, position_idx_comp_8, &
                               position_ref_box_idx_comp_4, position_ref_box_idx_comp_8, &
                               position_range_full_4, position_range_full_8 
                               
  procedure :: get_energy

#ifdef __HAS_SPIN__
  procedure :: create_particle_single_cell_p_s
  generic   :: create_particle => create_particle_single_cell_p, create_particle_single_cell, &
                                  create_particle_single_cell_p_s
#else
  generic   :: create_particle => create_particle_single_cell_p, create_particle_single_cell
#endif
  procedure :: validate           => validate_spec

  procedure :: restart_write    => restart_write_species
  procedure :: restart_read     => restart_read_species

end type t_species
!-----------------------------------------------------------------------------------------

abstract interface 
  subroutine dudt_species( this, emf, dt, i0, i1, energy, time )
    import t_species, t_emf, p_double

    class( t_species ), intent(inout) :: this
    class( t_emf ), intent(in) ::  emf
    real(p_double), intent(in) :: dt
    integer, intent(in) :: i0, i1
    real(p_double), intent(inout) :: energy
    real(p_double), intent(in) :: time
  end subroutine
end interface

#ifdef __HAS_SPIN__

interface
  subroutine dsdt_species( this, ep, bp, p_old, gamma, dt, ptrcur, chunk_size )
    import t_species, p_double, p_k_part

    class( t_species ), intent(inout) :: this
    real(p_k_part), intent(in), dimension(:,:) :: bp, ep, p_old
    real(p_k_part), intent(in), dimension(:) :: gamma
    real(p_double), intent(in) :: dt
    integer, intent(in) :: ptrcur, chunk_size
  end subroutine
end interface

#endif

interface 
  subroutine dep_current_1d( jay, dxi, xnew, ixold, xold, q, rgamma, u, np, dt )
    import t_vdf, p_k_part, p_double

    type( t_vdf ), intent(inout) :: jay
    integer, dimension(:,:), intent(inout) :: dxi, ixold
    real(p_k_part), dimension(:,:), intent(inout) :: xnew, xold ,u
    real(p_k_part), dimension(:), intent(inout) :: q,rgamma
    integer, intent(in) :: np
    real(p_double), intent(in) :: dt
  end subroutine
end interface

interface
  subroutine dep_current_2d(  jay, dxi, xnew, ixold, xold, q, rgamma, u, np, dt ) 
    import t_vdf, p_k_part, p_double

    type( t_vdf ), intent(inout) :: jay
    integer, dimension(:,:), intent(inout) :: dxi, ixold
    real(p_k_part), dimension(:,:), intent(inout) :: xnew, xold ,u
    real(p_k_part), dimension(:), intent(inout) :: q,rgamma
    integer, intent(in) :: np
    real(p_double), intent(in) :: dt
  end subroutine
end interface

interface
  subroutine dep_current_3d(  jay, dxi, xnew, ixold, xold, q, np, dt )
    import t_vdf, p_k_part, p_double

    type( t_vdf ), intent(inout) :: jay
    integer, dimension(:,:), intent(inout) :: dxi, ixold
    real(p_k_part), dimension(:,:), intent(inout) :: xnew, xold
    real(p_k_part), dimension(:), intent(inout) :: q
    integer, intent(in) :: np
    real(p_double), intent(in) :: dt
  end subroutine
end interface

interface
  subroutine restart_write_species( this, restart_handle )
    import t_species, t_restart_handle
    
    class( t_species ), intent(in) :: this
    type( t_restart_handle ), intent(inout) :: restart_handle
  end subroutine
end interface

interface
  subroutine restart_read_species( this, restart_handle )
    import t_species, t_restart_handle
    
    class( t_species ), intent(inout) :: this
    type( t_restart_handle ), intent(in) :: restart_handle
  end subroutine
end interface

! this extra type necessary to make an array of pointers in fortran
type :: t_spec_arr

  class( t_species ), pointer :: s => null()

end type t_spec_arr

interface
  subroutine init_tmp_buf_current_spe_bnd()
  end subroutine
end interface

interface
  subroutine cleanup_tmp_buf_current_spe_bnd()
  end subroutine
end interface

interface
  subroutine init_species( this, sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
              no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
              recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmax )

  import t_species, t_emf, t_current, t_space, t_grid, t_node_conf, t_restart_handle, &
         t_options, t_vdf_msg, t_part_idx, t_spec_msg, p_double, t_time_step

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
  type( t_restart_handle ), intent(in) :: restart_handle
  type(t_options), intent(in) :: sim_options
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: tmax 

  end subroutine
end interface

interface 
  subroutine set_dudt_spec( this ) 

    import t_species 
    
    class( t_species ), intent(inout) :: this
  
  end subroutine
end interface

interface
  subroutine push_species( this, emf, current, t, tstep, tid, n_threads, options )

  import t_species, t_emf, t_vdf, p_double, t_time_step, t_options

  class( t_species ), intent(inout) :: this
  class( t_emf ), intent( inout )  ::  emf
  type( t_vdf ), intent(inout) :: current

  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  integer, intent(in) :: tid    ! local thread id
  integer, intent(in) :: n_threads  ! total number of threads
  type( t_options ), intent(in) :: options

  end subroutine
end interface

interface 
  subroutine get_emf_spec( this, emf, bp, ep, np, ptrcur, t )

    import t_species, t_emf, p_k_part, p_double

    class( t_species ), intent(in) :: this
    class( t_emf ), intent(in), target :: emf
    real(p_k_part), dimension(:,:), intent(out) :: bp, ep
    integer, intent(in) :: np, ptrcur
    real(p_double), intent(in), optional :: t

  end subroutine
end interface

interface
  subroutine read_input_species( this, input_file, def_name, periodic, if_move, grid, &
               dt, read_prof, sim_options )

  import t_species, t_input_file, t_options, p_double, t_grid

  class( t_species ),  intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  character(len = *),    intent(in) :: def_name
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ),               intent(in) :: grid

  real(p_double), intent(in) :: dt
  logical, intent(in) :: read_prof

  type(t_options), intent(in) :: sim_options

  end subroutine
end interface

interface
  subroutine get_diag_buffer_size_spec( spec, gnx, diag_buffer_size )

  import t_species

  class( t_species ), intent(in) :: spec
  integer, dimension(:), intent(in) :: gnx
  integer, intent(inout) :: diag_buffer_size
  end subroutine
end interface

interface
  subroutine report_species( this, emf, g_space, grid, no_co, tstep, t, send_msg, recv_msg )

  import t_species, t_emf, t_space, t_grid, t_node_conf, t_time_step, p_double, t_vdf_msg

  class( t_species ),  intent(inout) :: this
  class( t_emf ),      intent(inout) :: emf
  type( t_space ),       intent(in) :: g_space
  class( t_grid ),        intent(in) :: grid
  class( t_node_conf ),   intent(in) :: no_co
  type( t_time_step ),   intent(in) :: tstep
  real(p_double),        intent(in) :: t
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  end subroutine
end interface

interface
  subroutine report_energy_species( this, no_co, tstep, t )

  import t_species, t_node_conf, t_time_step, p_double

  class( t_species ),  intent(in) :: this
  class( t_node_conf ),   intent(in) :: no_co
  type( t_time_step ),   intent(in) :: tstep
  real(p_double),        intent(in) :: t

  end subroutine
end interface

interface
  subroutine fill_data_species(this, tstep)
  import t_species, t_time_step 
  class( t_species ), intent(inout)    :: this
  type( t_time_step ),   intent(in) :: tstep
  end subroutine
end interface

interface
  subroutine get_phasespace_axis_spec( spec, xp, l, lp, x_or_p, xp_dim )

  import t_species, p_diag_prec 

  class(t_species), intent(in) :: spec
  real(p_diag_prec), dimension(:), intent(out) :: xp
  integer, intent(in) :: l, lp
  integer, intent(in) :: x_or_p
  integer, intent(in) :: xp_dim

  end subroutine
end interface

interface
  subroutine reshape_spec( this, old_grid, new_grid, msg_patt, no_co, send_msg, recv_msg )

  import t_species, t_msg_patt, t_grid, t_node_conf, t_vdf_msg

  class( t_species ), intent(inout) :: this
  type( t_msg_patt ), intent(in) :: msg_patt
  class( t_grid ), intent(in) :: old_grid, new_grid
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  end subroutine
end interface

interface
  subroutine get_quant_spec( this, i1, i2, quant, q )

  import t_species, p_k_part

  class( t_species ), intent(in) :: this
  integer, intent(in) :: i1, i2
  integer, intent(in) :: quant
  real(p_k_part), dimension(:), intent(out) :: q

  end subroutine
end interface

interface
  subroutine deposit_density_spec( this, charge, i1, i2, q )

  import t_species, t_vdf, p_k_part
  
  class( t_species ), intent(in) :: this
  type( t_vdf ), intent(inout) :: charge
  integer, intent(in) :: i1, i2
  real(p_k_part), dimension(:), intent(in) :: q
  end subroutine
end interface

interface
  subroutine enumerate_quants_spec( this, diagFile, track_set )

  import t_species, t_diag_file, t_track_set

  class( t_species ), intent(in) :: this
  class( t_diag_file ),intent(inout) :: diagFile
  type( t_track_set ), intent(inout) :: track_set

  end subroutine
end interface

interface
  subroutine update_boundary_species( this, jay, no_co, dt, bnd_cross, node_cross, send_msg, recv_msg )

    import t_species, t_vdf, t_node_conf, p_double, t_current, t_part_idx, t_spec_msg

    class( t_species ), intent(inout) :: this
    class( t_current ),      intent(inout) :: jay
    class( t_node_conf ),   intent(in) :: no_co
    real(p_double),        intent(in) :: dt
    type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
    type( t_part_idx ), intent(inout) :: node_cross
    type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  end subroutine
end interface

interface
  subroutine phys_boundary_species( this, current, dt, i_dim, bnd_idx, par_idx, npar )

    import t_species, t_current, p_double

    class( t_species ), intent(inout) :: this
    class( t_current ), intent(inout) :: current
    real(p_double), intent(in) :: dt
    integer, intent(in) :: i_dim, bnd_idx
    integer, dimension(:), intent(in) :: par_idx
    integer, intent(in) :: npar

  end subroutine
end interface

interface
  subroutine cleanup_species( this )
  import t_species
  class( t_species ), intent(inout) :: this
  end subroutine
end interface

!-----------------------------------------------------------------------------------------
! Fortran 2003 detritus required for the species utility proceedures
!-----------------------------------------------------------------------------------------
interface
  subroutine species_inject( this, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, send_msg, recv_msg )

    import t_species, t_current, t_node_conf, t_part_idx, t_spec_msg

    class( t_species ), intent(inout) :: this
    integer, dimension(:, :), intent(in) :: ig_xbnd_inj
    class( t_current ), intent(inout)   :: jay
    class( t_node_conf ), intent(in)     :: no_co
    type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
    type( t_part_idx ), intent(inout) :: node_cross
    type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  end subroutine
end interface

interface
  subroutine create_particle_single_cell( this, ix, x, q)

    import t_species, p_k_part

    class(t_species), intent(inout) :: this
    integer, dimension(:), intent(in) :: ix
    real(p_k_part), dimension(:), intent(in) :: x
    real(p_k_part),               intent(in) :: q

  end subroutine
end interface

interface
  subroutine create_particle_single_cell_p( this, ix, x, p, q )

    import t_species, p_k_part

    class(t_species), intent(inout) :: this
    integer, dimension(:), intent(in) :: ix
    real(p_k_part), dimension(:), intent(in) :: x
    real(p_k_part), dimension(:), intent(in) :: p
    real(p_k_part),               intent(in) :: q

  end subroutine
end interface

#ifdef __HAS_SPIN__
interface
  subroutine create_particle_single_cell_p_s( this, ix, x, p, s, q )

    import t_species, p_k_part

    class(t_species), intent(inout) :: this
    integer, dimension(:), intent(in) :: ix
    real(p_k_part), dimension(:), intent(in) :: x
    real(p_k_part), dimension(:), intent(in) :: p
    real(p_k_part), dimension(:), intent(in) :: s
    real(p_k_part),               intent(in) :: q

  end subroutine
end interface
#endif

!-----------------------------------------------------------------------------------------

public :: t_diag_species, t_track, t_track_set
public :: t_phasespace_diagnostics, t_phasespace, t_phasespace_list, t_phasespace_params
public :: t_spe_bound, t_spec_msg, t_piston, t_species, t_part_idx, t_udist
public :: t_spec_arr

!public :: p_none, p_uniform, p_pw_linear, p_gaussian, p_channel, p_sphere, p_func
public :: p_cell_near, p_cell_low
public :: p_spec_buf_block

public :: p_max_spname_len

!-----------------------------------------------------------------------------------------
! Density profile definitions
!-----------------------------------------------------------------------------------------
#include "psource/os-psource-def.f03"


contains

subroutine allocate_objs_spec( this )
!-----------------------------------------------------------------------------------------
!       Allocate any objects contained within the species (t_species) object
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species ), intent(inout) :: this

  ! Allocate default t_diag_species object class
  if ( .not. associated( this%diag ) ) then
    allocate( this%diag )
  endif

  ! Allocate default t_spe_bound object class
  if ( .not. associated( this%bnd_con ) ) then
    allocate( this%bnd_con )
  endif

  ! If needed, 'allocate_objs' on any of classes we just made so that any subobjects are created
  call this%diag%allocate_objs()
  !call this%bnd_con%allocate_objs()

end subroutine

!-----------------------------------------------------------------------------------------
! Return the number of spatial dimensions on the x buffer
! - This will generally just be the same as p_x_dim
! - Some simulation modes may override this
!-----------------------------------------------------------------------------------------
pure function get_n_x_dims_spec( this )

  implicit none

  class( t_species ), intent(in) :: this
  integer :: get_n_x_dims_spec

  get_n_x_dims_spec = p_x_dim

end function get_n_x_dims_spec

!-----------------------------------------------------------------------------------------
! Return species kinetic energy on local node
!-----------------------------------------------------------------------------------------
function get_energy( this )

  implicit none

  class( t_species ), intent(in) :: this
  real( p_double ) :: get_energy

  integer :: i
  real( p_double ) :: energy, u2, gamma, kin, cell_volume
  
  ! Check if energy was calculated during the push (time-centered), otherwise
  ! calculate it here
  if ( this % energy(1) == p_ene_recalc ) then
     ! Energy not available, recalculate it
     energy = 0

     !$omp parallel do private(u2, gamma, kin) reduction(+ : energy)
     do i = 1, this%num_par
        u2 = this%p(1,i)**2 + this%p(2,i)**2 + this%p(3,i)**2
        gamma = sqrt( u2 + 1 )
        kin = u2 / (gamma + 1)
        energy = energy + this%q(i) * kin
     enddo
     !$omp end parallel do

     ! Store value so it can be reused
     this%energy(1) = energy
   #ifdef _OPENMP
     do i = 2, ubound(this % energy, 1)
        this % energy(i) = 0
     enddo
   #endif

  else

     ! Get energy calculated previously
     energy = this % energy(1)

   #ifdef _OPENMP
     ! If doing a multi-threaded run, add contribution from all threads
     do i = 2, ubound(this % energy, 1)
        energy = energy + this % energy(i)
     enddo
   #endif

  endif

  ! Normalize to cell size and charge over mass ratio
  cell_volume = this % dx(1)
  do i = 2, p_x_dim
    cell_volume = cell_volume * this % dx(i)
  enddo

  get_energy = energy * this%rqm * cell_volume


end function get_energy
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! get_position routines
!
! These routines return particle positions indexed to the global simulation box
!  - For these to work the species%g_box( p_lower, : ) needs to be shifted 0.5 cells
!    with respect to the global simulation box. See the t_species definition for details
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine position_idx_comp_4( this, comp, idx, np, pos )
!-----------------------------------------------------------------------------------------


   implicit none

   class( t_species ), intent(in) :: this
   integer, intent(in) :: comp
   integer, intent(in), dimension(:) :: idx
   integer, intent(in) :: np
   real( p_single ), dimension(:), intent(out) :: pos

   integer :: j, ixmin
   real( p_single ) :: xmin, dx

   xmin  = real( this%g_box( p_lower, comp ), p_single )
   ixmin = this%my_nx_p( 1, comp ) - 2
   dx    = real( this%dx( comp ), p_single )

   do j = 1, np
    pos(j) = xmin + dx * real((this%ix( comp, idx(j) ) + ixmin) + &
                               this%x ( comp, idx(j) ) , p_single )
   enddo

end subroutine position_idx_comp_4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine position_idx_comp_8( this, comp, idx, np, pos )

   implicit none

   class( t_species ), intent(in) :: this
   integer, intent(in) :: comp
   integer, intent(in), dimension(:) :: idx
   integer, intent(in) :: np
   real( p_double ), dimension(:), intent(out) :: pos

   integer :: j, ixmin
   real( p_double ) :: xmin, dx

   xmin  = this%g_box( p_lower, comp )
   ixmin = this%my_nx_p( 1, comp ) - 2
   dx = this%dx( comp )

   do j = 1, np
    pos(j) = xmin + dx * (real(( this%ix( comp, idx(j) ) + ixmin), p_double ) + &
                                 this%x( comp, idx(j) ) )
   enddo

end subroutine position_idx_comp_8
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine position_ref_box_idx_comp_4( this, comp, idx, np, pos, if_pos_ref_box )

   implicit none

   class( t_species ), intent(in) :: this
   integer, intent(in) :: comp
   integer, intent(in), dimension(:) :: idx
   integer, intent(in) :: np
   real( p_single ), dimension(:), intent(out) :: pos
   logical, intent(in) :: if_pos_ref_box

   integer :: j, ixmin
   real( p_single ) :: xmin, dx

   if (if_pos_ref_box) then
      xmin = 0.0_p_single
   else
      xmin  = real( this%g_box( p_lower, comp ), p_single )
   endif
   ixmin = this%my_nx_p( 1, comp ) - 2
   dx    = real( this%dx( comp ), p_single )

   do j = 1, np
    pos(j) = xmin + dx * real((this%ix( comp, idx(j) ) + ixmin) + &
                               this%x ( comp, idx(j) ) , p_single )
   enddo

end subroutine position_ref_box_idx_comp_4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine position_ref_box_idx_comp_8( this, comp, idx, np, pos, if_pos_ref_box )

   implicit none

   class( t_species ), intent(in) :: this
   integer, intent(in) :: comp
   integer, intent(in), dimension(:) :: idx
   integer, intent(in) :: np
   real( p_double ), dimension(:), intent(out) :: pos
   logical, intent(in) :: if_pos_ref_box

   integer :: j, ixmin
   real( p_double ) :: xmin, dx

   if (if_pos_ref_box) then
      xmin = 0.0_p_double
   else
      xmin  = this%g_box( p_lower, comp )
   endif
   ixmin = this%my_nx_p( 1, comp ) - 2
   dx = this%dx( comp )

   do j = 1, np
    pos(j) = xmin + dx * (real(( this%ix( comp, idx(j) ) + ixmin), p_double ) + &
                                 this%x( comp, idx(j) ) )
   enddo

end subroutine position_ref_box_idx_comp_8
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------

subroutine position_single_4( this, idx, pos )

   implicit none

   class( t_species ), intent(in) :: this
   integer, intent(in) :: idx
   real( p_single ), dimension(:), intent(out) :: pos

   integer :: j

   do j = 1, p_x_dim
    pos(j) = real( ( (this%ix( j, idx ) + this%my_nx_p(p_lower, j) - 2) + &
                      this%x( j, idx ) ) * this%dx( j ) + &
                      this%g_box( p_lower, j ), p_single )
   enddo

end subroutine position_single_4
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine position_single_8( this, idx, pos )

   implicit none

   class( t_species ), intent(in) :: this
   integer, intent(in) :: idx
   real( p_double ), dimension(:), intent(out) :: pos

   integer :: j

   do j = 1, p_x_dim
    pos(j) = this%g_box( p_lower, j ) + this%dx( j ) * &
               ( real((this%ix( j, idx ) + this%my_nx_p(p_lower, j) - 2), p_double) + &
                       this%x( j, idx ) )
   enddo

end subroutine position_single_8
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------

subroutine position_range_comp_4( this, comp, idx0, idx1, pos )
!-----------------------------------------------------------------------------------------


   implicit none

   class( t_species ), intent(in) :: this
   integer, intent(in) :: comp, idx0, idx1
   real( p_single ), dimension(:), intent(out) :: pos

   integer :: j, ixmin
   real( p_single ) :: xmin, dx

   xmin  = real( this%g_box( p_lower, comp ), p_single )
   ixmin = this%my_nx_p( 1, comp ) - 2
   dx = real( this%dx( comp ), p_single )

   do j = idx0, idx1
    pos(j-idx0+1) = xmin + dx * (real(( this%ix( comp, j ) + ixmin ) + &
                                        this%x( comp, j ), p_single ) )
   enddo

end subroutine position_range_comp_4
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------

subroutine position_range_comp_8( this, comp, idx0, idx1, pos )
!-----------------------------------------------------------------------------------------


   implicit none

   class( t_species ), intent(in) :: this
   integer, intent(in) :: comp, idx0, idx1
   real( p_double ), dimension(:), intent(out) :: pos

   integer :: j, ixmin
   real( p_double ) :: xmin, dx

   xmin  = this%g_box( p_lower, comp )
   ixmin = this%my_nx_p( 1, comp ) - 2
   dx    = this%dx( comp )

   do j = idx0, idx1
    pos(j-idx0+1) = xmin + dx * (real(( this%ix( comp, j ) + ixmin ), p_double ) + &
                                        this%x( comp, j ) )
   enddo

end subroutine position_range_comp_8
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------

subroutine position_range_full_4( this, idx0, idx1, pos )
!-----------------------------------------------------------------------------------------


   implicit none

   class( t_species ), intent(in) :: this
   integer, intent(in) :: idx0, idx1
   real( p_single ), dimension(:,:), intent(out) :: pos

   integer :: j 
   integer ixmin(p_x_dim)
   
   real( p_single ) :: xmin( p_x_dim ), dx( p_x_dim )

   xmin(1:p_x_dim)  = real( this%g_box( p_lower, 1:p_x_dim ), p_single )
   ixmin(1:p_x_dim) = this%my_nx_p( 1, 1:p_x_dim ) - 2
   dx(1:p_x_dim) = real( this%dx( 1:p_x_dim ), p_single )

   do j = idx0, idx1
    pos(:,j-idx0+1) = xmin(:) + dx(:) * (real(( this%ix( :, j ) + ixmin(:) ) + &
                                        this%x( :, j ), p_single ) )
   enddo

end subroutine position_range_full_4
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------

subroutine position_range_full_8( this, idx0, idx1, pos )
!-----------------------------------------------------------------------------------------


   implicit none

   class( t_species ), intent(in) :: this
   integer, intent(in) :: idx0, idx1
   real( p_double ), dimension(:,:), intent(out) :: pos

   integer :: j, ixmin(p_x_dim)
   real( p_double ) :: xmin(p_x_dim), dx(p_x_dim)

   xmin(1:p_x_dim)  = this%g_box( p_lower, 1:p_x_dim )
   ixmin(1:p_x_dim) = this%my_nx_p( 1, 1:p_x_dim ) - 2
   dx(1:p_x_dim)    = this%dx(1:p_x_dim)

   do j = idx0, idx1
    pos(:,j-idx0+1) = xmin(:) + dx(:) * (real(( this%ix( :, j ) + ixmin( : ) ) + &
                                        this%x( :, j ), p_double ) )
   enddo

end subroutine position_range_full_8
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Printout the algorithm used by the pusher
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_spec( this )

  implicit none
  class( t_species ), intent(in) :: this

  print *, ' '
  print *, trim(this%name),' :'

  ! Push type
  select case ( this%push_type )
  case (p_std)
    print *, '- Standard (Boris) pusher'
  case (p_simd)
    print *, '- SIMD optimized pusher'
  case (p_radcool)
    print *, '- Radiation cooling'
  case (p_vay)
    print *, '- Vay velocity push'
  case (p_fullrot)
    print *, '- Full magnetic rotation pusher'
  case (p_euler)
    print *, '- Euler-Rodriguez magnetic rotation pusher'
  case (p_cond_vay)
    print *, '- Conditional (gamma>5) standard/Vay pusher'
  case (p_cary)
    print *, '- Higuera-Cary pusher'
  case (p_exact)
    print *, '- Exact pusher'
  case (p_exact_rr)
    print *, '- Exact pusher with radiation reaction'
  end select

  if ( this%free_stream ) then
    print *, '- Free streaming particles (no dudt)'
  write(0,*) '- (*warning*) ', trim(this%name), ' are free streaming!'
  endif

end subroutine list_algorithm_spec
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Allocate buffers for particle quantities
!-----------------------------------------------------------------------------------------
subroutine init_buffer_spec( this, num_par_req )

  implicit none

  class(t_species), intent(inout) :: this
  integer, intent(in) :: num_par_req

  integer :: n_x_dim

  ! The buffer size must always be a multiple of vector width in size because of SIMD code
  this%num_par_max = (( num_par_req + p_vecwidth - 1 ) / p_vecwidth) * p_vecwidth

  ! Under some compilers / configurations (e.g. gfortran 4.9.1, OS X, single precision)
  ! if the total q buffer size is below 1Kb it may not be allocated to a 32bit boundary
  ! which is required by AVX code
  if ( this%num_par_max < 256 ) this%num_par_max = 256

  ! setup position buffer
  n_x_dim = this%get_n_x_dims()
  call freemem( this%x )
  call alloc(  this%x, (/ n_x_dim, this%num_par_max /))

  ! initialize particle cell information
  call freemem( this%ix )
  call alloc( this%ix, (/ p_x_dim, this%num_par_max /) )

  ! setup momenta buffer
  call freemem( this%p )
  call alloc( this%p, (/ p_p_dim, this%num_par_max /) )

#ifdef __HAS_SPIN__

  ! setup spin buffer
  call freemem( this%s )
  call alloc( this%s, (/ p_s_dim, this%num_par_max /) )

#endif

  ! setup particle charge buffer
  call freemem( this%q )
  call alloc( this%q, (/ this%num_par_max /))

  ! initialize tracking data if necessary
  if ( this%add_tag ) then
     call freemem( this%tag )
     call alloc( this%tag, (/ 2, this%num_par_max /))
  endif

end subroutine init_buffer_spec
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Grow the particle buffers
!-------------------------------------------------------------------------------
subroutine grow_buffer_spec( this, num_par_req )

  implicit none

  class(t_species), intent(inout) :: this
  integer, intent(in) :: num_par_req

  real(p_k_part), dimension(:), pointer   :: temp1_r
  real(p_k_part), dimension(:,:), pointer :: temp2_r
  integer, dimension(:,:), pointer :: temp2_i

  integer :: num_par_old, num_par_new, n_x_dim

  ! The buffer size must always be a multiple of vector width in size because of SIMD code
  num_par_new = (( num_par_req + p_vecwidth - 1 ) / p_vecwidth) * p_vecwidth

  if ( this%num_par > 0 ) then

     num_par_old = this%num_par

     write(0,'(A,I0,A,A)') '[', mpi_node(), '] (* warning *) resizing particle buffers for this ', &
                            trim(this%name)
     write(0,'(A,I0,A,I0,A,I0)') '[', mpi_node(), '] (* warning *) Buffer size: ', this%num_par_max, ' -> ', num_par_new
     write(0,'(A,I0,A,I0)') '[', mpi_node(), &
             '] (* warning *) Number of particles currently in buffer: ', this%num_par

     if ( num_par_new <= num_par_old ) then
       ERROR('Invalid size for new buffer')
       call abort_program( p_err_invalid )
     endif

     ! particle positions (may not be p_x_dim)
     n_x_dim = this%get_n_x_dims()
     call alloc( temp2_r, (/ n_x_dim, num_par_new /))
     call memcpy( temp2_r, this%x, n_x_dim * num_par_old )
     call freemem( this%x )
     this%x => temp2_r

     ! particle cell index
     call alloc( temp2_i, (/ p_x_dim, num_par_new /))
     call memcpy( temp2_i, this%ix, p_x_dim * num_par_old )
     call freemem( this%ix )
     this%ix => temp2_i

     ! particle momenta
     call alloc( temp2_r, (/ p_p_dim, num_par_new /))
     call memcpy( temp2_r, this%p, p_p_dim * num_par_old )
     call freemem( this%p )
     this%p => temp2_r

     ! particle charge
     call alloc( temp1_r, (/ num_par_new /) )
     call memcpy( temp1_r, this%q, num_par_old )
     call freemem( this%q )
     this%q => temp1_r

#ifdef __HAS_SPIN__

     ! particle spin
    call alloc( temp2_r, (/ p_s_dim, num_par_new /))
    call memcpy( temp2_r, this%s, p_s_dim * num_par_old )
    call freemem( this%s )
    this%s => temp2_r

#endif

     ! Resize tracking data if necessary
     if ( this%add_tag ) then
        call alloc( temp2_i, (/ 2, num_par_new /) )
        call memcpy( temp2_i, this%tag, 2 * num_par_old )
        call freemem( this%tag )
        this%tag => temp2_i
     endif

     ! resize ionization data if necessary
     !(...)

     this%num_par_max = num_par_new
     write(0,'(A,I0,A)') '[', mpi_node(), '] (* warning *) resize successfull!'

  else

    ! no particles in buffer, simply reallocate the buffers
    call this % init_buffer( num_par_new )

  endif

end subroutine grow_buffer_spec
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Checks if all particle values are ok.
!-------------------------------------------------------------------------------
subroutine validate_spec( this, msg, over )

  implicit none

  ! dummy variables

  class(t_species), intent(in) :: this
  character( len = * ) , intent(in) :: msg
  logical, intent(in), optional :: over

  ! local variables
  integer :: i_dim, k
  integer, dimension(p_x_dim) :: ilb, iub
  logical :: over_

  ! Verify buffers
  if ((.not. associated(this%x)) .or. (.not. associated(this%ix)) .or. &
      (.not. associated(this%p)) .or. (.not. associated(this%q))) then

    write(0,'(A,I0,A,A)') '[', mpi_node(), '] ', trim(msg)
    write(0,*)  '[', mpi_node(), '] buffers are not allocated'
    call abort_program()
  endif

  if ((size(this%x,1)/=this%get_n_x_dims()) .or. (size(this%ix,1)/=p_x_dim) .or. &
      (size(this%p,1)/=p_p_dim)) then
        write(0,'(A,I0,A,A)') '[', mpi_node(), '] ', trim(msg)
        write(0,*)  '[', mpi_node(), '] buffer dimensions are invalid'
        call abort_program()
      endif

    if ((size(this%x,2)/=this%num_par_max) .or. (size(this%ix,2)/=this%num_par_max) .or. &
        (size(this%p,2)/=this%num_par_max) .or. (size(this%q)/=this%num_par_max)) then
        write(0,'(A,I0,A,A)') '[', mpi_node(), '] ', trim(msg)
        write(0,*)  '[', mpi_node(), '] buffer size is invalid'
        call abort_program()
    endif

  ! get grid boundaries
  do k = 1, p_x_dim
    ilb(k) = 1 - this%move_num(k)
    iub(k) = this%my_nx_p(3,k)
  enddo

  if ( present( over) ) then
    over_ = over
  else
    over_ = .false.
  endif

  ! allow for 1 cell overflow (ok before update boundary)
  if ( over_ ) then
    ilb = ilb - 1
    iub = iub + 1
  endif

  ! validate positions
  do i_dim = 1, p_x_dim
    do k = 1, this%num_par
       if ( ( this%x(i_dim, k)  < -0.5_p_k_part ) .or. &
            ( this%x(i_dim, k)  >= 0.5_p_k_part ) .or. &
            ( this%ix(i_dim, k) <  ilb( i_dim ) ) .or. &
            ( this%ix(i_dim, k) >  iub( i_dim)  ) ) then

                SCR_MPINODE('over = ', over_)
                call bad_particle( k, this, ilb, iub, msg // " - Invalid position " )

       endif
    enddo
  enddo

  ! validate momenta
  do k = 1, this%num_par
    do i_dim = 1, p_p_dim
       if ( isinf( this%p(i_dim, k) ) .or. isnan( this%p(i_dim, k) ) ) then

          call bad_particle( k, this, ilb, iub, msg // " - Invalid momenta " )

       endif
    enddo
  enddo

#ifdef __HAS_SPIN__

  ! validate spin
  do k = 1, this%num_par
    do i_dim = 1, p_s_dim
       if ( isinf( this%s(i_dim, k) ) .or. isnan( this%s(i_dim, k) ) ) then

          call bad_particle( k, this, ilb, iub, msg // " - Invalid spin " )

       endif
    enddo
  enddo

#endif

  ! validate charge
  do k = 1, this%num_par
     if ( this%q(k) == 0 .or. isinf( this%q(k) ) .or. isnan( this%q(k) ) ) then

        call bad_particle( k, this, ilb, iub, msg // " - Invalid charge " )

     endif
  enddo

contains

subroutine bad_particle( k, this, ilb, iub, msg )

  implicit none

  integer, intent(in) :: k

  class(t_species), intent(in) :: this
  integer, dimension(:), intent(in) :: ilb, iub
  character( len = * ), intent(in) :: msg

  write(0,'(A,I0,A,A)') '[', mpi_node(), '] ', trim(msg)

  write(0,'(A,I0,A,I0,A,I0)') "[", mpi_node(), "] Bad particle ", k, " of ", this%num_par
  write(0,*) "[", mpi_node(), "] p (:)  =", this%p(:, k)
  write(0,*) "[", mpi_node(), "] x (:)  =", this%x(:, k)
  write(0,*) "[", mpi_node(), "] ix (:) =", this%ix(:, k)
  write(0,*) "[", mpi_node(), "] q      =", this%q(k)

#ifdef __HAS_SPIN__
  write(0,*) "[", mpi_node(), "] s (:)  =", this%s(:, k)
#endif

  write(0,*) "[", mpi_node(), "] ilb(:) = ", ilb
  write(0,*) "[", mpi_node(), "] iub(:) = ", iub

  write(0,'(A,I0,A,A)') '[', mpi_node(), '] (* error *) Validate species failed for ', &
                                            trim(this%name), ' aborting...'

  call abort_program( p_err_invalid )


end subroutine bad_particle

end subroutine validate_spec
!-------------------------------------------------------------------------------

end module m_species_define
