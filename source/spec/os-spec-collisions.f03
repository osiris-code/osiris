!-------------------------------------------------------------------------------
! Collision module
!-------------------------------------------------------------------------------
!
!
!
! ------ToDo------
!
! Errata and incomplete work includes but is not limited to
!  - Implement linear and cubic interpolator for Nanbu equation
!    - Verify/refactor Newton solver for Nanbu equation
!  - Communicate upstream to the load balancer to prevent node breaks through collision cells
!    - Verify the grid shape at read_nml() instead of setup()
!
!
!
! ---Bibliography---
!
! On collisions in PIC
!
! - Takizuka, T., & Abe, H. (1977). A binary collision model for plasma simulation with a 
!   particle code. Journal of Computational Physics, 25(3), 205–219. doi:10.1016/0021-9991(77)90099-7
! - Nanbu, K. (1997). Theory of cumulative small-angle collisions in plasmas. Physical Review E,
!    55(4), 4642–4652. doi:10.1103/PhysRevE.55.4642
! - Nanbu, K., & Yonemura, S. (1998). Weighted Particles in Coulomb Collision Simulations
!   Based on the Theory of a Cumulative Scattering Angle. Journal of Computational Physics, 145, 639–654.
! - Sentoku, Y., Mima, K., Kishimoto, Y., & Honda, M. (1998). Effects of Relativistic 
!   Binary Collisions on PIC Simulations of Laser Plasmas. Journal of the Physical Society
!   of Japan, 67(12), 4084–4088.
! - Sentoku, Y., & Kemp, A. J. (2008). Numerical methods for particle simulations at extreme
!   densities and temperatures: Weighted particles, relativistic collisions and reduced currents.
!   Journal of Computational Physics, 227(14), 6846–6861. doi:10.1016/j.jcp.2008.03.043
! - Peano, F., Marti, M., Silva, L., & Coppa, G. (2009). Statistical kinetic treatment of
!   relativistic binary collisions. Physical Review E, 79(2), 025701. doi:10.1103/PhysRevE.79.025701
! - Perez, F., Gremillet, L., Decoster, A., Drouin, M., & Lefebvre, E. (2012) Improved
!   modeling of relativistic collisions and collisional ionization in particle-in-cell codes.
!   Physics of Plasmas, 19, 083104, http://dx.doi.org/10.1063/1.4742167
!
! On the physics of collisions
!
! - Spitzer, Lyman jr., Physics of Fully Ionized Gases, Second Edition (1962), Dover
! - Lifshitz, L.M., and Pitaevskii, L.P., Physical Kinetics (1981), Pergamon Press
! - Rosenbluth, M. N., MacDonald, W. M., & Judd, D. L. (1957). Fokker-Planck Equation
!   for an Inverse-Square Force. Physical Review, Vol. 107, No.1, pp. 1–6
! - Huba, J.D., NRL Plasma Formulary (2004), Office of Naval Research
!
!
!
! ---Back-of-the-envelope Calculations---
!
! As a very rough rule of thumb, for electron-electron and electron-ion collisions,
! the ratio of a collision frequency to the plasma frequency,
! up to a geometric factor of order 10, is given by
! \nu_{coll} / \omega_{pe} = Z ln(g) / g
! where Z is the ion charge state, omega_{pe} is the electron plasma frequency, and g is the plasma parameter
! g = n_e \lambda_d^3
! 
! For electron-ion energy exchange (as opposed to elastic scattering) there is a reduction by a factor of \sqrt{m/M}
!
! Numerically, the plasma parameter (as defined above) is 
! g = u_{th}^3 \sqrt{\frac{2.25 \times 10^{34} cm^{-3}}{n_e}}
! where u_{th} is the normalized thermal momentum ( u = p/mc = \gamma \beta .)

#include "os-config.h"
#include "os-preprocess.fpp"

#ifdef __HAS_COLLISIONS__

module m_species_collisions

#include "memory/memory.h"
  
  ! I can use either m_random or m_random_class to get t_random it seems
  ! Why is m_random able to share t_random if it's defined in m_random_class?
  !use m_random
  use m_random_class
  use m_random_hash
  use m_random_kiss
  use m_random_mt
  use m_restart
  use m_species
  use m_space
  use m_math
  use m_parameters
  use m_species_define
  use m_species
  use m_system
  use m_logprof
  use stringutil ! why is there no m_ in the name? Oh well.
  
  implicit none
  
  private

!-------------------------------------------------------------------------------
!  Collision parameters and types
!-------------------------------------------------------------------------------
  
  character(len=*), parameter :: p_collisions_rst_id = "collisions rst data - 0x0000"
  
  ! Parameters for Nanbu collision model
  integer, parameter :: p_max_lookup_path_length = 128
  ! The values are chosen to have an exact binary representation
  real(p_k_part), parameter  :: p_useA_min        = 1.40988e-3_p_k_part
  real(p_k_part), parameter  :: p_useA_max        = 6.0_p_k_part
  real(p_k_part), parameter  :: p_lookup_min      = 0.009765625e0_p_k_part
  real(p_k_part), parameter  :: p_lookup_max      = 3.0_p_k_part
  real(p_k_part), parameter  :: p_lookup_delta    = 1.220703125e-4_p_k_part
  real(p_k_part), parameter  :: p_r_lookup_delta  = 8192.0_p_k_part
  integer, parameter  :: p_lookup_size     = 24498
  ! These can be input variables - it might be best to put them in the type.
  real(p_k_part), dimension(:), pointer :: lookup_table_A
  real(p_k_part), dimension(:), pointer :: lookup_table_Ap
  
  ! Physical parameters in cgs units
  real(p_k_part), parameter  ::  cgs_c  = 2.9979e10_p_k_part ! [units of cm/sec]
  real(p_k_part), parameter  ::  cgs_e  = 4.8032e-10_p_k_part ! [units of statcoulomb]
  real(p_k_part), parameter  ::  cgs_me = 9.1094e-28_p_k_part ! [units of grams]
  real(p_k_part), parameter  ::  cgs_hbar = 1.0546e-27_p_k_part ! [units of erg-seconds]
  
  ! Parameters for switches in input deck
  
  ! Nanbu Collision Model
  ! cf Nanbu, K. (1997), 'Theory of cumulative small-angle collisions in plasmas',
  ! Physical Review E, 55(4), 4642–4652. doi:10.1103/PhysRevE.55.4642
  integer, parameter :: p_CollisionModel_Nanbu = 1
  
  ! Sentoku Collision Model
  ! cf Sentoku et al., "Effects of Relativistic Binary Collisions on PIC Simulations of Laser Plasmas"
  ! Journal of the Physical Society of Japan, Vol. 67, no. 12, pp. 4084-4088 (1998)
  integer, parameter :: p_CollisionModel_Sentoku = 2
  
  ! Nanbu Collision Model, with correction from theta_CoM to theta_lab
  ! c.f. Anything? Or just optimism?
  integer, parameter :: p_CollisionModel_Nanbu_Rel = 3
  
  ! For testing only! Scatters into 4pi at all timesteps!
  integer, parameter :: p_CollisionModel_Isotropic = 4
  
  ! Takizuka non-relativistic model
  integer, parameter :: p_CollisionModel_Takizuka = 5
  
  ! Perez
  integer, parameter :: p_CollisionModel_Perez = 6
  
  ! Correction for Monte-Carlo collision statistic to came the number of collisions relativistically invariant
  ! c.f. Peano et al., "Statistical kinetic treatment of relativistic binary collisions"
  ! Physical Review E (Rapid Communication), 79(2), 025701. doi:10.1103/PhysRevE.79.025701
  integer, parameter :: p_RelativisticCrossSectionCorrection_None = 0
  integer, parameter :: p_RelativisticCrossSectionCorrection_Rejection = 1
  integer, parameter :: p_RelativisticCrossSectionCorrection_Scattering = 2 ! Only inside Sentoku collisions right now (don't know why)
  
  ! Correction for boosted frame(?)
  ! This should probably be gamma_fluid, not gamma_particle
  ! (And Josh hasn't convinced himself the exponent has the right sign)
  integer, parameter :: p_FrameCorrection_None = 0
  integer, parameter :: p_FrameCorrection_dt = 1
  integer, parameter :: p_FrameCorrection_dt_n = 2
  
  ! Solver used to invert the Nanbu equation (Nanbu 1997 eq. 13)
  integer, parameter :: p_RootFinder_None = 0
  integer, parameter :: p_RootFinder_Newton = 1
  integer, parameter :: p_RootFinder_Linear = 2
  integer, parameter :: p_RootFinder_Cubic = 3
  
  integer, parameter :: p_LongestName = 128 ! Long enough right? :)
  
  ! 4th root of the largest subnormal number, used to avoid dividing by zero
  ! (3rd root would naively be sufficinent, but an extra power for saftey seems wise)
  real(p_k_part), parameter :: mr_small = sqrt(sqrt(tiny(0.0_p_k_part)))
  
  !-------------------------------------------------------------------------------
  ! Collision class
  !-------------------------------------------------------------------------------
  
  type :: t_collisions
    
    ! Number of dt per collision cycle (set to 0 to turn off collisions)
    integer :: n_collide
    
    ! Size of collision cell (in PIC cells) in each direction
    integer, dimension(p_x_dim) :: nx_collision_cells
    
    ! Total number of PIC cells per collision cell
    integer :: collision_cell_size
    
    ! Normalization for charge density n_0
    ! An input deck parameter
    real(p_k_part) :: norm_charge_density ![e cm^-3]
    
    ! The plasma frequency for the normalization charge density
    ! Calculated automatically, not an input deck parameter
    real(p_k_part) :: omega_pe_0 ![rad/sec]
    
    ! Input deck switch to turn on (or off) calculation of the Coulomb Logarithm
    logical        :: coulomb_logarithm_automatic = .true.
    
    ! Flags for collision model, cross section correction methods, etc.
    integer :: CollisionModel = p_CollisionModel_Perez
    integer :: CrossCorrection = p_RelativisticCrossSectionCorrection_None
    logical :: CrossCorrTime = .false.
    integer :: FrameCorrection = p_FrameCorrection_None
    integer :: RootFinder = p_RootFinder_None
    
    ! Value of ln(lambda) if the code isn't calculating it
    real(p_k_part) :: coulomb_logarithm_value = 0.0_p_k_part
    
    ! file names of lookup tables
    character(len=p_max_lookup_path_length) :: lookup_file_A, lookup_file_Ap
    
    integer :: number_of_collision_cells
    class (t_random), pointer :: node_rng => NULL()
    integer :: rng_seed
    class (t_random), dimension(:), pointer :: cell_rngs => NULL()
    
    logical :: PerezLowTempCorrection = .false.
    
  end type t_collisions
  
  !-------------------------------------------------------------------------------
  !  (some) Fluid quantities (and number density) in the collision cell for a given species
  !-------------------------------------------------------------------------------
  type :: t_local_fluid_state
    
    real(p_k_part) :: e_kin_ave   ! Mean kinetic energy (in fluid reference frame) [units of m_e c**2]
    real(p_k_part) :: q_tot       ! Total charge in cell
    real(p_k_part), dimension(p_p_dim) :: u_fluid  ! Fluid velocity (NOT the mean of all velocities) [units of c]
    real(p_k_part) :: rsq_ldebye ! (1/lambda_debeye^2) [units of omega_p0^2/c^2]
    real(p_k_part) :: charge_density ! Charge density [units of e/cm^3]
    real(p_k_part) :: number_density ! Number density
    
    real(p_k_part) :: p_th2 ! the scalar 'thermal momentum' (squared to avoid having to take roots when comparing)
  end type t_local_fluid_state
  
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  End parameters and types
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
  
  integer :: collide_part_ev=0, sort_coll_ev, sort_coll_genidx_ev, sort_coll_rearrange_ev
  
  interface list_algorithm_collisions
    module procedure list_collision_algorithms
  end interface
  
  interface read_nml
    module procedure read_nml_coll
  end interface
  
  interface setup
    module procedure setup_coll
  end interface
  
  interface reshape
    module procedure reshape_collisions
  end interface
  
  interface cleanup
    module procedure cleanup_coll
  end interface
  
  interface collide_particles
    module procedure collide_particles
  end interface
  
  interface if_collide
    module procedure if_collide
  end interface
  
  interface A_lookup_setup
    module procedure A_lookup_setup_linear
    module procedure A_lookup_setup_cubic
  end interface A_lookup_setup
  
  interface lorentz_transform_proper
    module procedure lorentz_transform_proper_double
    module procedure lorentz_transform_proper_single
  end interface lorentz_transform_proper
  
  interface alloc
    module procedure alloc_2d_local_fluid_state
  end interface

  interface freemem
    module procedure freemem_2d_local_fluid_state
  end interface
  
  interface restart_write
    module procedure restart_write_collisions
  end interface restart_write
  
  interface restart_read
    module procedure restart_read_collisions
  end interface restart_read
  
  public :: read_nml, setup, cleanup, collide_particles, if_collide
  public :: reshape
  public :: p_max_lookup_path_length
  public :: list_algorithm_collisions
  public :: restart_write, restart_read
  
  public :: t_collisions
       
contains 

!-------------------------------------------------------------------------------
subroutine read_nml_coll( self, input_file, species_head, num_species, sim_options )
!-------------------------------------------------------------------------------
! Read input file parameters 
!-------------------------------------------------------------------------------
  
  use m_input_file
  
  implicit none
  
  type( t_collisions ), intent(inout) :: self
  class( t_input_file ), intent(inout) :: input_file
  class ( t_species ), pointer :: species_head
  integer, intent(in) :: num_species
  type( t_options ), intent(in) :: sim_options
  
  integer  :: n_collide
  integer, dimension(p_x_dim) :: nx_collision_cells
  logical          :: coulomb_logarithm_automatic
  real(p_k_part)  :: coulomb_logarithm_value
  
  character( len = p_LongestName ) :: collision_model
  character( len = p_LongestName ) :: cross_section_correction
  logical                          :: timestep_cross_corr
  character( len = p_LongestName ) :: frame_correction
  character( len = p_LongestName ) :: root_finder
  
  integer :: rng_seed
  
  logical :: perez_low_temp_correction
  
  namelist /nl_collisions/ n_collide, nx_collision_cells, &
                         coulomb_logarithm_automatic, &
                         coulomb_logarithm_value, &
                         collision_model, cross_section_correction, &
                         timestep_cross_corr, frame_correction, &
                         root_finder, rng_seed, &
                         perez_low_temp_correction
                         
  integer :: ierr
  
  n_collide = 0
  nx_collision_cells = 1
  
  coulomb_logarithm_automatic = .true.
  coulomb_logarithm_value = 0.0_p_k_part
  
  collision_model = ""
  cross_section_correction = ""
  timestep_cross_corr = .false.
  frame_correction = ""
  root_finder = ""
  
  rng_seed = 1337*(mpi_node()+1)
  
  perez_low_temp_correction = .false.
  
  ! Get namelist text from input file
  call get_namelist( input_file, "nl_collisions", ierr )
  
  if ( ierr == 0 ) then
    read (input_file%nml_text, nml = nl_collisions, iostat = ierr)
    if (ierr /= 0) then
      print *, "Error reading collisions parameters"
      print *, "aborting..."
      stop
    endif
    
    if ( num_species == 0 ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "(*warning*) 'collisions' section found, but no species were defined"

        write(0,*) "(*warning*) 'collisions' section will be ignored"

      endif

      self%n_collide = 0
      return
    endif
    
  else if (ierr < 0) then
    print *, "Error reading collisions parameters"
    print *, "aborting..."
    stop
  endif
  
  self%n_collide = n_collide
  self%nx_collision_cells = nx_collision_cells
  self%coulomb_logarithm_automatic = coulomb_logarithm_automatic
  self%coulomb_logarithm_value = coulomb_logarithm_value
  
  ! reference frequency / density is now set globally
  self%omega_pe_0          = sim_options%omega_p0
  self%norm_charge_density = sim_options%n0
  
  self%rng_seed = rng_seed
  
  self%PerezLowTempCorrection = perez_low_temp_correction
  
  select case ( lowercase(trim(collision_model)) )
    case ("nanbu")
      ERROR("Nanbu model is poorly implemented at this point, so much so that I'm not going to let you use it")
      ERROR("Perez model is recommended, otherwise please review the code and update the Nanbu model yourself")
      ERROR("Shutting down")
      call abort_program(p_err_notimplemented)
      self%CollisionModel = p_CollisionModel_Nanbu
    case ("sentoku")
      WARNING("Sentoku model currently seems to have a factor of 2 error")
      WARNING("(Yes we could just fix it, but haven't actually found the error, it's just off by 2")
      WARNING("(so, for now, it persists)")
      self%CollisionModel = p_CollisionModel_Sentoku
    case ("nanbu_rel")
      ERROR("'Nanbu relativistic' model is a bad idea")
      ERROR("Use Perez Model")
      call abort_program(p_err_notimplemented)
      self%CollisionModel = p_CollisionModel_Nanbu_Rel
    case ("isotropic")
      self%CollisionModel = p_CollisionModel_Isotropic
    case ("takizuka", "takizuka and abe", "ta", "non-rel", "non-relativistic")
      self%CollisionModel = p_CollisionModel_Takizuka
    case ("perez")
      self%CollisionModel = p_CollisionModel_Perez
    case ("")
      self%CollisionModel = p_CollisionModel_Perez
    case default
      ERROR("Error: Can't parse Collision Model in Input Deck")
      call abort_program(p_err_invalid)
  end select
  
  select case ( lowercase(trim(cross_section_correction)) )
    case ("none")
      self%CrossCorrection = p_RelativisticCrossSectionCorrection_None
    case ("rejection")
      self%CrossCorrection = p_RelativisticCrossSectionCorrection_Rejection
    case ("scattering")
      self%CrossCorrection = p_RelativisticCrossSectionCorrection_Scattering
    case ("")
      self%CrossCorrection = p_RelativisticCrossSectionCorrection_None
    case default
      ERROR("Error: Can't parse Relativistic Cross Section Correction in Input Deck")
      call abort_program(p_err_invalid)
  end select
   
  self%CrossCorrTime = timestep_cross_corr
  
  select case ( lowercase(trim(frame_correction)) )
    case ("none")
      self%FrameCorrection = p_FrameCorrection_None
    case ("time")
      self%FrameCorrection = p_FrameCorrection_dt
    case ("time_denisty")
      self%FrameCorrection = p_FrameCorrection_dt_n
    case ("")
      self%FrameCorrection = p_FrameCorrection_None
    case default
      ERROR("Error: Can't parse Relativistic Frame Correction in Input Deck")
      call abort_program(p_err_invalid)
  end select
  
  select case ( lowercase(trim(root_finder)) )
    case ("none")
      self%RootFinder = p_RootFinder_None
    case ("newton")
      self%RootFinder = p_RootFinder_Newton 
    case ("linear")
      self%RootFinder = p_RootFinder_Linear
    case ("cubic")
      self%RootFinder = p_RootFinder_Cubic
    case ("")
      self%RootFinder = p_RootFinder_None
    case default
      ERROR("Error: Can't parse Root Finding Alogrithm in Input Deck")
      call abort_program(p_err_invalid)
  end select
  
  if( n_collide > 0 ) then
    call verify_input_deck( self, species_head )
  endif !n_collide
  
end subroutine read_nml_coll
!------------------------------------------------------------------------------

subroutine verify_input_deck( self, species_head )
  
  implicit none
  
  ! Dummy variables
  type( t_collisions ), intent(in) :: self
  class ( t_species ), pointer :: species_head
  
  ! Local variables
  class( t_species ), pointer :: sp1, sp2
  
  real(p_double) :: push_start_time
  logical :: any_species_to_collide
  
  ! Executable statements
    
    ! Make sure species is always sorted at collision time

    sp1 => species_head
    do
      if( .not. associated(sp1)) exit
      if ( MODULO(self%n_collide, sp1%n_sort) /= 0 ) then
        SCR_ROOT('(*warning*) n_collide is not a multiple of n_sort for species ', trim(sp1%name))
        SCR_ROOT('(*warning*) setting n_sort to n_collide for this species.')
        sp1%n_sort = self%n_collide
      endif
      sp1 => sp1%next
    enddo
  
  ! Verify that there actually is a colliding species
  ! If not and n_collide is set, it's probably a mistake so abort
  any_species_to_collide = .false.
  sp1 => species_head
  do
    if (.not. associated(sp1)) exit
    if ( sp1%if_collide ) then
      any_species_to_collide = .true.
      exit
    endif
    sp1 => sp1%next
  enddo
  if ( .not. any_species_to_collide ) then
    ERROR('n_collide set, but no species to be collided')
    ERROR('Abort !')
    call abort_program()
  endif
  
  ! Verify the colliding species has something to collide with
  ! (This could be done much simpler and clearer with a few logical variables,
  ! but, well, I didn't do that)
  sp1 => species_head
  outer_loop: do
    if (.not. associated(sp1) ) exit outer_loop
    if ( sp1%if_collide ) then
      if( sp1%if_like_collide) exit outer_loop ! Redundant but clearer
      sp2 => sp1%next
      inner_loop: do
        if (.not. associated(sp2)) exit inner_loop
        if( sp2%if_collide ) exit outer_loop
        sp2 => sp2%next
      enddo inner_loop
      ! At this point there is only one species to collide and it doesn't collide itself, so...
      ERROR('Only one species to collide, and if_like_collide == .false.')
      call abort_program()
    endif
    sp1 => sp1%next
  enddo outer_loop
  
  ! Find start time of first colliding species
  ! (They all have to be the same so this is sufficient)
  sp1 => species_head
  do
    if (.not. associated(sp1)) exit
    if( sp1%if_collide ) then
      push_start_time = sp1%push_start_time
      exit
    endif
    sp1 => sp1%next
  enddo
  
  ! And now verify that all push_start_time's are consistent
  sp1 => species_head
  do
    if (.not. associated(sp1)) exit
    if( sp1%if_collide ) then
      if( push_start_time /= sp1%push_start_time ) then
        ERROR( "push_start_time of colliding species not compatible:" )
        ERROR( "all colliding species must have same push_start_time!" )
        ERROR( "ABORT!" )
        call abort_program(p_err_invalid)
      endif
    else
      if( push_start_time < sp1%push_start_time ) then
        ERROR( "push_start_time of non colliding species not compatible:" )
        ERROR( "non-colliding species must have push_start_time not greater than colliding species!" )
        ERROR( "ABORT!" )
        call abort_program(p_err_invalid)
      endif
    endif
    sp1 => sp1%next
  enddo
  
  ! If we're not calculating ln(lambda), make sure it was set
  if ( .not.(self%coulomb_logarithm_automatic) .and. &
       (self%coulomb_logarithm_value == 0.0_p_k_part) ) then
    ERROR( "It is not wise to set the Coulomb Logarithm to 0 if it is not calculated automatically." )
    ERROR( "ABORT!" )
    call abort_program(p_err_invalid)
  endif
  
  ! Make sure that q_real (i.e. Z) was set for each species
  ! (We could just assume Z=1, but safer to require the user be explicit)
  sp1 => species_head
  do
    if (.not. associated(sp1)) exit
    if ( sp1%if_collide .and. ( sp1%q_real == 0 ) ) then
      ERROR( 'No q_real was given for species ', sp1%name )
      ERROR( 'Charge (q_real) needed for collisions!' )
      call abort_program(p_err_invalid)
    endif
    sp1 => sp1%next
  enddo
  
  ! If we need it, make sure a root finder was specified
  if( self%CollisionModel == p_CollisionModel_Nanbu &
      .or. self%CollisionModel == p_CollisionModel_Nanbu_Rel) then
    if( self%RootFinder == p_RootFinder_None ) then
      ERROR( "Cumulative scattering model selected, but not root finding method specified" )
      call abort_program(p_err_invalid)
    endif
  endif
  
  ! A normalization density is strictly needed
  if ( self%norm_charge_density <= 0.0_p_k_part ) then
    ERROR('No norm number density given.')
    ERROR('Needed for collisions. ABORT!')
    call abort_program()
  endif
  
  if ( self%CrossCorrection == p_RelativisticCrossSectionCorrection_None &
       .and. self%CrossCorrTime) then
    ERROR('Timestep adjusted for Cross Section Correction, but no correction set.')
    call abort_program()
  endif
  
  if ( self%PerezLowTempCorrection .and. (self%CollisionModel .ne. p_CollisionModel_Perez) ) then
    ERROR('Low-temperature correction only implemented for Perez collision Model')
    call abort_program(p_err_invalid)
  endif
  
end subroutine verify_input_deck

!-------------------------------------------------------------------------------
subroutine setup_coll( self, g_nx, l_nx, restart, restart_handle)
!-------------------------------------------------------------------------------
! Setup collision data
!-------------------------------------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  type( t_collisions ), intent(inout) :: self
  integer, dimension(:), intent(in) :: g_nx, l_nx
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  
  ! Local Variables
  integer :: dim, i
  
  ! Executable statements
  
  if ( self%n_collide > 0 ) then
    
    ! reference frequency and density are now set globally
    
    do dim=1, p_x_dim
      !check compatibility of pic grid vs. collision grid
      if ( modulo( g_nx(dim), self%nx_collision_cells(dim) ) /= 0 ) then
        write(0,'(A,I0)') ' (*error*) Global PIC grid (nx_p) is not a multiple of collision grid for dimension ', dim
        call abort_program( p_err_invalid )
      endif
      
      !check compatibility of local pic grid vs. collision grid
      if ( modulo( l_nx(dim), self%nx_collision_cells(dim) ) /= 0) then
        write(0,'(A,A,I0)') ' (*error*) Local node grid (nx_p/node_num) is not a multiple ', &
                   'of collision grid for dimension ', dim
        call abort_program( p_err_invalid )
      endif
    enddo
    
    ! Setup lookup table for A
    select case ( self%RootFinder )
      
      case (p_RootFinder_Linear)
        call A_lookup_setup()
      case (p_RootFinder_Cubic)
        call A_lookup_setup(self%lookup_file_A, self%lookup_file_Ap)
      case default
        continue
    end select
    
    self%collision_cell_size = 1
    do dim = 1, p_x_dim
      self%collision_cell_size = self%collision_cell_size * self%nx_collision_cells(dim)
    enddo
    
    ! Calculate the number of collision cells
    self%number_of_collision_cells = 1
    do dim=1, p_x_dim
      self%number_of_collision_cells = self%number_of_collision_cells * (l_nx(dim) / self%nx_collision_cells(dim))
    enddo
    
    !call init_genrand( self%node_rng, 1337*(mpi_node()+1) )
    ! IF I REMOVE THIS LINE CAN I NOT USE m_random?
    !self%node_rng => new_random( p_random_mt )
    allocate( t_random_mt :: self%node_rng )
    call self%node_rng % init( self%rng_seed )
    !call alloc(self%rng_array, (/self%number_of_collision_cells/))
    allocate(t_random_kiss :: self%cell_rngs( self%number_of_collision_cells ))
    
    do i = 1, self%number_of_collision_cells
      !call init_genrand( self%rng_array(i), genrand_int32(self%node_rng))
      !self%cell_rngs(i) => new_random( p_random_hash )
      call self%cell_rngs(i) % init( self%node_rng%genrand_int32() )
    enddo
    
    ! Create collision events
    if(collide_part_ev==0) then
      collide_part_ev        = create_event('particle collisions (total)')
      sort_coll_ev           = create_event('particle sort_coll (total)')
      sort_coll_genidx_ev    = create_event('particle sort_coll, gen. idx')
      sort_coll_rearrange_ev = create_event('particle sort_coll, rearrange particles')
    endif
    
    if( restart) then
      call restart_read(self, restart_handle)
    endif
    
  endif !collisions
  
end subroutine setup_coll
!-------------------------------------------------------------------------------

subroutine reshape_collisions( self, l_nx )
  
  implicit none
  
  ! Dummy variables
  type( t_collisions ), intent(inout) :: self
  integer, dimension(:), intent(in) :: l_nx
  
  ! Local Variables
  integer :: dim, i
  
  if ( self%n_collide > 0 ) then
  
    do dim=1, p_x_dim
      !check compatibility of local pic grid vs. collision grid
      if ( modulo( l_nx(dim), self%nx_collision_cells(dim) ) /= 0) then
        ERROR('Local node grid (nx_p/node_num) is not a multiple of collision grid for dimension ', dim)
        call abort_program( p_err_invalid )
      endif
    enddo
  
    self%collision_cell_size = 1
    do dim = 1, p_x_dim
      self%collision_cell_size = self%collision_cell_size * self%nx_collision_cells(dim)
    enddo
    
    ! Calculate the number of collision cells
    self%number_of_collision_cells = 1
    do dim=1, p_x_dim
      self%number_of_collision_cells = self%number_of_collision_cells * (l_nx(dim) / self%nx_collision_cells(dim))
    enddo
    
    ! I kind of want to communicate and save the rng seeds and indexes, so the state wont change at rebalance,
    ! But the system rng definitely changes anyway, so I think we just have to make sure the rebalances
    ! happen at the same time if we're running tracks or anything.
    !call freemem(self%rng_array)
    deallocate(self%cell_rngs)
    !call alloc(self%rng_array, (/self%number_of_collision_cells/))
    allocate(t_random_hash :: self%cell_rngs(self%number_of_collision_cells))
    
    do i = 1, self%number_of_collision_cells
      !call init_genrand( self%rng_array(i), genrand_int32(self%node_rng))
      !self%cell_rngs(i) => new_random( p_random_hash )
      call self%cell_rngs(i) % init( self%node_rng%genrand_int32() )
    enddo
    
  endif
  
end subroutine reshape_collisions

!-------------------------------------------------------------------------------
subroutine list_collision_algorithms( self )
!-------------------------------------------------------------------------------
! Print out the algorithm details for collisions
!-------------------------------------------------------------------------------
  
  implicit none
  
  type( t_collisions), intent(in) :: self
  
  SCR_ROOT( '' )
  SCR_ROOT( 'Collisions:' )
  
  ! Report collision model
  select case ( self%CollisionModel )
     case(p_CollisionModel_Nanbu)
        SCR_ROOT( "- Model: Nanbu, cumulative" )
     case(p_CollisionModel_Sentoku)
        SCR_ROOT( "- Model: Sentoku, relativistic" )
     case(p_CollisionModel_Nanbu_Rel)
        SCR_ROOT( "- Model: Nanbu, relativistic" )
     case(p_CollisionModel_Isotropic)
        SCR_ROOT( "- Model: Uniform scattering angle" )
     case(p_CollisionModel_Takizuka)
        SCR_ROOT( "- Model: Non-Relativistic Takizuka&Abe Model" )
     case(p_CollisionModel_Perez)
        SCR_ROOT( "- Model: Perez Relativistic Extension of Nanbu Model" )
     case default
        SCR_ROOT( "- Model: UNKNOWN - Please update os-spec-collisions%list_collision_algorithms" )
  end select
  
  ! Report correction type for crosssection
  select case ( self%CrossCorrection )
     case(p_RelativisticCrossSectionCorrection_None)
        SCR_ROOT( "- Type of crosssection correction: none" )
     case(p_RelativisticCrossSectionCorrection_Rejection)
        SCR_ROOT( "- Type of crosssection correction: rejection method" )
     case(p_RelativisticCrossSectionCorrection_Scattering)
        SCR_ROOT( "- Type of crosssection correction: scattering angle distribution" )
     case default
        SCR_ROOT( "- Type of crosssection correction: UNKNOWN - Please update os-spec-collisions%list_collision_algorithms" )
  end select
  
  ! Report correction type for crossection time
  select case ( self%CrossCorrTime )
     case (.false.)
        SCR_ROOT( "- Type of crosssection TIME correction: none" )
     case (.true.)
        SCR_ROOT( "- Type of crosssection TIME correction: dt = dt * 2" )
     !case default
     !   SCR_ROOT( "- Type of crosssection TIME correction: huh?" )
     ! Right, no need for case default here. Doesn't feel right but okay.
  end select
  
  ! Report correction type for timestep
  select case( self%FrameCorrection )
     case(p_FrameCorrection_None)  
        SCR_ROOT( "- Type of frame correction: none" )
     case(p_FrameCorrection_dt)  
        SCR_ROOT( "- Type of frame correction: correction for delta T only" )
     case(p_FrameCorrection_dt_n)  
        SCR_ROOT( "- Type of frame correction: correction for delta T and density" )
     case default
        SCR_ROOT( "- Type of frame correction: UNKNOWN - Please update os-spec-collisions%list_collision_algorithms" )
  end select
  
  ! Report which solver is used
  select case( self%RootFinder )
     case(p_RootFinder_None)
       continue
     case(p_RootFinder_Newton)
        SCR_ROOT( "- Type of solver used for A(s): newton" )
     case(p_RootFinder_Linear)
        SCR_ROOT( "- Type of solver used for A(s): linear interpolation" )
     case(p_RootFinder_Cubic) 
        SCR_ROOT( "- Type of solver used for A(s): cubic interpolation" )
     case default
        SCR_ROOT( "- Type of solver used for A(s): UNKNOWN - Please update os-spec-collisions%list_collision_algorithms" )
  end select
  
end subroutine list_collision_algorithms
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine collide_particles( self, species_head, num_species, t, dt, n, n_threads )
!-------------------------------------------------------------------------------
!       sort and collide requesting particles 
!-------------------------------------------------------------------------------
  
#ifdef _OPENMP
  use omp_lib
#endif
  
  type( t_collisions ), intent(inout) :: self
  ! head of species linked list (don't point to anything else, use local pointer to iterate)
  class ( t_species ), pointer :: species_head
  
  integer, intent(in) :: num_species 
  
  real(p_double), intent(in) :: t
  real(p_double), intent(in) :: dt
  integer, intent(in) :: n
  integer, intent(in) :: n_threads
  
  ! local variables
  
  ! Iterators for species objects
  integer :: sp_id, sp_id2
  ! Pointers for species objects I guess now
  class (t_species), pointer :: sp1, sp2
  
  ! Iterator for collision cells
  integer ::  i_cell
  
  ! Time interval between two collision cycles = n_collide*dt [1/w0]
  real(p_k_part) :: dt_n
  
  real(p_k_part) :: lambda_debye
  
  !total number of collision cells (on local node) and their size in dx^3
  integer :: number_of_collision_cells, collision_cell_size
  
  ! Max number of particles one in collision cell
  ! Needed to allocate sufficiently large permutation array
  integer :: num_par_max_coll_cell
  
  ! Array for permutation vectors for each species for this cell
  ! (Simply unused for non-colliding species)
  ! dimension( num_species, num_par_max_coll_cell + 1)
  integer, dimension(:,:,:), pointer :: permutation_vector
  
  ! Array with local fluid quantities for each species
  ! (Even non-colliding species will contribute to the Debye length,
  ! if this is physical or not is, well, a good question)
  type(t_local_fluid_state), dimension(:,:), pointer :: cell_fluid_states
  
  ! Array for ids of last particle for each collision cell for each species
  ! dimension(num_species, number_of_collision_cells)
  integer, dimension(:,:), pointer :: last_particle_list
  
  ! Number of particles for each species for current cell
  integer, dimension(:,:), pointer :: coll_cell_num_par
  
  procedure(like_collide), pointer :: like_collide_p
  procedure(unlike_collide), pointer :: unlike_collide_p
  
  integer :: tid
  
  class( t_random ), pointer :: cell_rng
  
  ! point variables to null after initialization to avoid making them save variables
  sp1                => null()
  sp2                => null()
  permutation_vector => null()
  cell_fluid_states  => null()
  last_particle_list => null()
  coll_cell_num_par  => null()
  like_collide_p     => null()
  unlike_collide_p   => null()
  cell_rng           => null()
  
  ! Executable statements
  
  !num_species = size( species )
  
  ! Loop over all species, check if any exist that should collide and haven't been started
  ! If there are, loop over all species again sorting any that *have* started,
  ! then return without running collisions
!  do sp_id = 1, num_species
!    if( species( sp_id )%if_collide .and. ( species( sp_id )%push_start_time > t ) ) then 
!      do sp_id2 = 1, num_species        
!        if( t > species( sp_id2 )%push_start_time ) &
!          call sort( species( sp_id2 ), n, t )
!      enddo
!      return
!    endif
!  enddo


  ! Loop over all species, check if any exist that should collide but have not been started
  ! If so, then only sorting happens now;
  ! loop over all species again, sorting any that *have* started,
  ! and then return
  sp1 => species_head
  sortlevel1: do
    if (.not. associated(sp1)) exit sortlevel1
    if (sp1%if_collide .and. sp1%push_start_time > t) then
      sp2 => species_head
      sortlevel2: do
        if (.not. associated(sp2)) exit sortlevel2
        if ( t > sp2%push_start_time) call sort(sp2, n, t)
        sp2 => sp2%next
      enddo sortlevel2
      return
    endif
    sp1 => sp1%next
  enddo sortlevel1
  
  
  call begin_event( collide_part_ev )
  
  ! Calculate the number of collision cells and their size
  ! Obvious not though, must have happened here before, now just aliases
  number_of_collision_cells = self%number_of_collision_cells
  collision_cell_size = self%collision_cell_size
  
!  do i_cell = 1, number_of_collision_cells
!    !if ( self%cell_rngs(i_cell)%index .lt. 1 ) then
!    if ( self%cell_rngs(i_cell) % has_cycled() ) then
!      call self%cell_rngs(i_cell) % init(self%node_rng % genrand_int32() )
!    endif
!  enddo
  
  do i_cell = 1, number_of_collision_cells
    select type( cell_rng => self%cell_rngs(i_cell) )
      class is (t_random_hash)
        if ( cell_rng % has_cycled() ) then
          call cell_rng % init(self%node_rng % genrand_int32() )
        endif
      class default
        ! This isn't necessarily an error, just do nothing I think
      continue
    end select
  enddo
  
  cell_rng => null()
  
  call alloc( last_particle_list, (/ number_of_collision_cells +1, num_species /) )
  last_particle_list(:,:) = 0
  
  num_par_max_coll_cell=0
  sp1 => species_head; sp_id = 1
  do
    if ( .not. associated(sp1)) exit
    call sort_coll_species( sp1, n, t, self%nx_collision_cells, &
                            last_particle_list(2:,sp_id) )
    ! Find the greatest number of colliding particles for any cell or species
    if ( sp1%if_collide ) then
      do i_cell = 1, number_of_collision_cells
        num_par_max_coll_cell = max( num_par_max_coll_cell, &
          last_particle_list( i_cell+1, sp_id ) - last_particle_list( i_cell, sp_id ) )
      enddo 
    endif
    sp1 => sp1%next; sp_id = sp_id + 1
  enddo
  
  ! Time interval between 2 collision cycles
  dt_n = real( dt, p_k_part ) * self%n_collide
  
  call alloc( coll_cell_num_par, (/ num_species, n_threads /) )
  call alloc( cell_fluid_states, (/ num_species, n_threads /) )
  call alloc( permutation_vector, (/ num_par_max_coll_cell + 1, num_species, n_threads /) )
  
  ! Loop over all collision cells
  ! Oh crap, can I keep the species pointers private???? Seems like yeah?
  !$omp parallel do schedule(dynamic) private(sp_id, lambda_debye, sp_id2, tid, sp1, sp2)
  do i_cell = 1, number_of_collision_cells
    
    #ifdef _OPENMP
      tid = omp_get_thread_num()+1
    #else
      tid = 1
    #endif
    
    ! Calculate fluid quantities
    sp1 => species_head; sp_id = 1
    fluid_loop: do
      if( .not. associated(sp1)) exit fluid_loop
      
      coll_cell_num_par(sp_id, tid) = last_particle_list(i_cell+1, sp_id) - last_particle_list( i_cell, sp_id )
      
      cell_fluid_states(sp_id, tid) = find_t_local_fluid_state( &
                                      sp1, &
                                      last_particle_list( i_cell, sp_id ), &
                                      coll_cell_num_par(sp_id, tid), &
                                      collision_cell_size, &
                                      calculate_lambda_debye = self%coulomb_logarithm_automatic &
                                                              )
      
      sp1 => sp1%next; sp_id = sp_id + 1
    enddo fluid_loop
    
    ! Generate permutation list
    sp1 => species_head; sp_id = 1
    perm_loop: do
      if (.not. associated(sp1)) exit perm_Loop
      if(sp1%if_collide) then
        call generate_shuffle_list( permutation_vector(:,sp_id,tid), &
                                    last_particle_list(i_cell,sp_id), &
                                    coll_cell_num_par(sp_id,tid), &
                                    num_par_max_coll_cell + 1, &
                                    self%cell_rngs(i_cell) &
                                   )
      endif
      sp1 => sp1%next; sp_id = sp_id + 1
      
    enddo perm_loop
    
    ! Calculate total Debye length
    if( self%coulomb_logarithm_automatic ) then
      lambda_debye = 0.0
      do sp_id = 1, num_species
        ! if number of particles in cell is not enough then rsq_ldebye will be set to -1
        ! and should be ignored
        if ( cell_fluid_states(sp_id, tid)%rsq_ldebye > 0 ) then
          lambda_debye = lambda_debye + cell_fluid_states(sp_id, tid)%rsq_ldebye
        endif
      enddo
      if( lambda_debye > 0.0 ) lambda_debye = 1.0/sqrt(lambda_debye)
    else
      lambda_debye = -1.0_p_k_part
    endif
    
    select case (self%CollisionModel)
      case (p_CollisionModel_Takizuka)
        like_collide_p => like_collide_non_rel
        unlike_collide_p => unlike_collide_non_rel
      case (p_CollisionModel_Sentoku)
        like_collide_p => like_collide_sentoku
        unlike_collide_p => unlike_collide_sentoku
      case (p_CollisionModel_Perez)
        like_collide_p => like_collide_perez
        unlike_collide_p => unlike_collide_perez
      case default
        like_collide_p => like_collide
        unlike_collide_p => unlike_collide
    end select
    
    ! Call collision subroutines for this cell
    sp1 => species_head; sp_id = 1
    call_routines_ol: do
      if (.not. associated(sp1)) exit call_routines_ol
      ! I'm being somewhat stubborn, to avoid more nested ifs
      ! With linked lists this makes a lot less sense though
      if( .not. sp1%if_collide) then
        sp1 => sp1%next; sp_id = sp_id + 1
        cycle call_routines_ol
      endif
      
      ! Unlike collisions with the species after species( sp_id )
      sp2 => sp1%next; sp_id2 = sp_id + 1
      call_routines_il: do
        if( .not. associated(sp2)) exit call_routines_il
        if( .not. sp2%if_collide) then
          sp2 => sp2% next; sp_id2 = sp_id2 + 1
          cycle call_routines_il
        endif
            
        call unlike_collide_p( &
                                 self, &
                                 sp1, &
                                 sp2, &
                                 sp_id, &
                                 sp_id2, &
                                 last_particle_list(i_cell,:), &
                                 coll_cell_num_par(:,tid), &
                                 permutation_vector(:,:,tid), &
                                 cell_fluid_states(:,tid), &
                                 lambda_debye, &
                                 dt_n, &
                                 self%cell_rngs(i_cell) &
                             )
        
        sp2 => sp2%next; sp_id2 = sp_id2 + 1
      enddo call_routines_il ! Unlike collisions ( species with index > sp_id )
      
      ! Like collision
      if( sp1%if_like_collide ) then
        
        call like_collide_p( &
                               self, &
                               sp1, &
                               coll_cell_num_par(sp_id,tid), &
                               permutation_vector(:,sp_id,tid), &
                               cell_fluid_states(sp_id,tid), &
                               lambda_debye, &
                               dt_n, &
                               self%cell_rngs(i_cell) &
                           )
        
      endif
      
      sp1 => sp1%next; sp_id = sp_id + 1
    enddo call_routines_ol ! Species to collide
  
  enddo ! All collision cells
  !$omp end parallel do
  
  call freemem( last_particle_list )
  call freemem( coll_cell_num_par )
  call freemem( cell_fluid_states )
  call freemem( permutation_vector )
    
  call end_event( collide_part_ev )
  
end subroutine collide_particles
!-------------------------------------------------------------------------------

!---------------------------------------------------
subroutine like_collide( &
                         self, &
                         spec, &
                         num_par, &
                         shuffle_vec, &
                         cell_fluid_state, &
                         l_debye, &
                         dt_n, &
                         local_rng &
                         )
!---------------------------------------------------
!       update momenta for collision
!---------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  type( t_collisions ), intent(in) :: self
  class ( t_species ), intent(inout) :: spec
  integer, intent(in) :: num_par
  integer, dimension(:) :: shuffle_vec
  type(t_local_fluid_state), intent(in) :: cell_fluid_state
  ! total debye length [c/w0]
  real(p_k_part), intent(in) :: l_debye
  ! collision cycle: dt * n_collide [1/w0]
  real(p_k_part), intent(in) :: dt_n
  class (t_random), intent(inout) :: local_rng
  
  ! Local variables
  
  ! Masses
  real(p_k_part)   :: m_real, rq_real ! Masses of real particles (in simulation units) (just aliases)
  real(p_k_part)   :: mu_ab ! Reduced mass (=m/2, can probably prune) (unless we want to use the relativistically correct reduced mass)
  
  ! Weight of particle a, particle b, maximum or minimum of those (used for statistical scattering)
  real(p_k_part) :: w_a, w_b, w_m
  ! The sum of all w_m for each pair, used to normalize weighted averages
  !
  real(p_k_part) :: n_alpha_alpha
  
  ! Random number between zero and one ( [0,1) )
  real(p_k_part) :: random_unit_interval
  
  ! Coulomb logarithm, and its argument
  real(p_k_part) :: lambda, ln_lambda
  
  ! Average of u_rel**2 (used in log(lambda))
  real(p_k_part) :: u_sq_avg
  
  ! Collision 'probability'
  ! i.e. (4 pi (e_1 e_2)**2 n_p ln(lambda) / g_rel**2 beta_rel**3 m_e**2 c**3 ) * delta_t
  ! i.e. (w_pe**2 e**2 ln(lambda) / g_rel**2 beta_rel**3 m_c c**3) (just work it out, it should be true)
  real(p_k_part) :: s
  
  ! 'Base' collision 'probability'
  ! i.e. s = s_first / g_rel**2 beta_rel**3
  ! (With a possible relativistic correction factor which Josh doesn't trust)
  real(p_k_part) :: s_first
  
  ! Adjusted timestep (dt * n_collide * {weighting for statistics} * {possible factor of 2} )
  real(p_k_part) :: dt_s
  
  ! The value needed for the Nanbu collision models
  ! coth(A) - 1/A = exp(-s)
  ! 0 < A_Nanbu < inf
  ! A_Nanbu( s -> inf ) = 3 * exp(-s)
  ! A_Nanbu( s -> 0 )   = 1 / s
  ! A_Nanbu( s <= 0 ) invalid
  real(p_k_part)   :: A_Nanbu
  
  ! Reference plasma frequency (just an alias to self%omega_pe_0)
  real(p_k_part) :: omega_pe_0
  
  ! Proper 4-velocity of particles in LAB frame
  real(p_k_part), dimension(p_p_dim) :: u_a_lab, u_b_lab
  real(p_k_part)                     :: g_a_lab, g_b_lab
  
  ! Proper 4-velocity of either particle in the One Particle at Rest (OPR) frame
  ! Also referred to as the velocity of part_a RELative to part_b (or vice versa)
  real(p_k_part), dimension(p_p_dim) :: u_rel
  real(p_k_part)                     :: g_rel
  ! Magnitude of space-like proper-velocity, true velocity of a relative to b (or vice versa)
  real(p_k_part)                     :: au_rel, av_rel
  
  ! Proper 4-velocity of the Center-of-Momentum (CoM) frame (as seen in the LAB Frame)
  real(p_k_part), dimension(p_p_dim) :: ucm
  real(p_k_part)                     :: gcm
  
  ! Proper 4-velocity of either particle as seen in the CoM frame (equal because m_a = m_b)
  ! Note this is therefore also the velocity *OF* the CoM frame as seen in either OPR frame
  ! Note also the way it's calculated below, actually u_b_in_cm = -u_in_cm
  real(p_k_part), dimension(p_p_dim) :: u_in_cm
  real(p_k_part)                     :: g_in_cm
  ! Magnitude of the space-like proper velocity of either particle as seen in the CoM frame
  real(p_k_part)                     :: au_in_cm
  
  ! Scattering angles
  ! 0<phi<2pi, with equal distribution 
  ! -inf < tan(theta_CoM) < inf, with a gaussian distribution and variance given by the collision frequency
  ! (Though if tan(theta_CoM)!~theta you're well outside the realm of legitimacy of these models)
  real(p_k_part)   :: phi, theta_CoM
  
  ! Trigonometric functions of argument the scattering angles
  ! ( cos(phi), sin(phi), cos(theta_CoM), sin(theta_CoM) )
  real(p_k_part) :: cphi, sphi, ctheta, stheta
  
  ! theta in the One Particle at Rest frame
  real(p_k_part)  :: theta_opr
  
  ! Variables used for actual momenta update
  !
  ! Normalized CoM proper velocity
  ! (For like collisions this isn't strictly necessary, but it is for unlike so we use it here for consistency)
  real(p_k_part), dimension(p_p_dim) :: u_norm_a_cm
  
  ! The components of the normalized CoM proper velocity
  ! ( n_t = sqrt(n_x**2 + n_y**2) i.e. "n_transverse"; rn_t = 1/n_t)
  real(p_k_part) :: n_x, n_y, n_z, n_t, rn_t
  
  ! The change in momentum (u_new - u_old), still normalized
  real(p_k_part), dimension(p_p_dim) :: delta
  
  ! Updated space-like proper velocity after collision, in CoM
  real(p_k_part), dimension(p_p_dim) :: u_after_collision
  
  ! Numerical variables ( indexing )
  ! Iterator used for looping over pairs in the shuffle vector
  ! Used in loop to calculate the statistical weights 
  ! and in the actual collision loop
  integer :: particle_pair_loop_iterator
  
  ! Index of last pair to collide (needed in case num_par is odd)
  integer :: last_pair
  
  ! Translation from shuffle vector array to particle array
  integer :: id_a, id_b
  
  ! Self-explanatory
  real(p_k_part) :: devnull
  
  ! Executable statements
  
  ! Return if no pairs are in cell
  if( num_par < 2 ) return
  
  ! If there are an odd number of particles, collide the last with the first
  ! ( shuffle_vector( num_par + 1 ) = shuffle_vector(1) by construction )
  ! c.f. Nanbu 1998 Fig. 3(b)
  last_pair = num_par + mod( num_par, 2 )
  
  ! Aliases
  m_real  = spec%m_real
  rq_real = spec%rq_real
  omega_pe_0 = self%omega_pe_0
  ! It might be better to calculate the relativistic reduced mass (for each pair)
  mu_ab = m_real / 2.0_p_k_part
  
  ! Calculate statistical values in collision cell
  
  u_sq_avg = 0
  n_alpha_alpha = 0
  if( self%coulomb_logarithm_automatic ) then
    
    do particle_pair_loop_iterator = 1, last_pair, 2
      id_a = shuffle_vec(particle_pair_loop_iterator)
      id_b = shuffle_vec(particle_pair_loop_iterator + 1)
      
      ! Accumulate statistical weight
      w_m = min( abs(spec%q(id_a)), abs(spec%q(id_b)) ) 
      n_alpha_alpha = n_alpha_alpha + w_m  
      
      ! Get u_a_lab, u_b_lab, g_a_lab, g_b_lab
      u_a_lab = spec%p(:,id_a)
      u_b_lab = spec%p(:,id_b)
      g_b_lab = sqrt( 1.0_p_k_part + dot_product( u_b_lab, u_b_lab ) )
      
      ! Transform u_b from lab frame to ref frame of a
      call lorentz_transform_proper( u_a_lab, u_b_lab, g_b_lab, u_rel, devnull)
      ! c.f. ???
      ! I can't say right now why we're doing the weighting this way
      u_sq_avg = u_sq_avg + w_m * dot_product(u_rel,u_rel)
    enddo 
    u_sq_avg = u_sq_avg / n_alpha_alpha
    
  else ! Nothing to do but find the weighting
    
    do particle_pair_loop_iterator = 1, last_pair, 2
      id_a = shuffle_vec(particle_pair_loop_iterator)
      id_b = shuffle_vec(particle_pair_loop_iterator + 1)
      n_alpha_alpha = n_alpha_alpha + min( abs(spec%q(id_a)), abs(spec%q(id_b)) )
    enddo 
    
  endif
  
  if( self%coulomb_logarithm_automatic ) then
    
    ! For the lower limit, use either the distance of classical closest approach, or the first Born approximation
    ! of the quantum mechanical result for close scattering (see Landau and Lifshitz 'Quantum Mechanics'
    ! chapter 127) depending on which is greater (and of course the Debye length for the upper limit.)
    ! (Note that this line is equivalent to a simple test of u < fine structure constant, but for physical clarity it's written out explicitly )
    ! (Note also, this should technically be u*v in the classical limit, but gamma < 1.000027 so it doesn't matter)
    lambda = min( ( mu_ab*cgs_me * l_debye*(cgs_c/omega_pe_0) * u_sq_avg * cgs_c**2 ) / &
                  ( 2.0_p_k_part * spec%q_real**2 * cgs_e**2 ), &
                  ( mu_ab*cgs_me * l_debye*(cgs_c/omega_pe_0) * sqrt(u_sq_avg) * cgs_c ) / &
                  ( cgs_hbar ) )
    
    ! c.f. Spitzer 'Physics of Fully Ionized Gases' equation 5-12
    ! $\int_0^p dx \frac{x^3}{\sqrt{1+x^2}} = \frac{1}{2}\left(\frac{1}{1+p^2}-1+ln(1+p^2)\right)$
    ! Note that since this theory is only logarithmically correct, we could also use ln_lambda = log(lambda)
    ! here without loss of validity (as is often done)
    ln_lambda = 1./2 * ( 1./(1 + lambda**2) - 1 + log(1 + lambda**2) )
    
    ! Clamp log(lambda) at 2 (chosen, with some arbitrariness, to be consistent with Lee&More 1984 (c.f. discussion after eq. 17))
    ln_lambda = max( ln_lambda, 2._p_k_part )
    
  else
    ln_lambda = self%coulomb_logarithm_value
  endif
  
  ! Adjust time-step for collision statistics for unequal particle weights
  ! c.f. Nanbu 1998, eq.14 and the one immediately preceding it
  dt_s = dt_n * ( abs(cell_fluid_state%q_tot) / ( 2.0_p_k_part * n_alpha_alpha ) )
  if ( self%CrossCorrTime ) then
    ! Josh believes this is to correct for the Peano rejection method, as in the non-relativistic limit that will skip 50% of the collisions
    dt_s = 2.0_p_k_part * dt_s
  endif
  
  ! prepare s calculations
  !s_first = ln_lambda * ( spec%q_real**2 / mu_ab )**2 * &
  !          cell_fluid_state%number_density * dt_s * &
  !          ((cgs_e**2 * omega_pe_0)/(cgs_me * cgs_c**3))
  ! The line above and the calculation below are equivalent, but below is rearranged for clarity
  
  s_first = ln_lambda * &
            ( spec%q_real**2 / mu_ab )**2 * &
            (dt_s/omega_pe_0) * &
            cell_fluid_state%number_density * omega_pe_0**2 * &
            cgs_e**2 / (cgs_me * cgs_c**3)
  
  ! Loop over pairs of particles and 'collide' them
  do particle_pair_loop_iterator = 1, last_pair, 2
    ! Get particle indexes
    id_a = shuffle_vec(particle_pair_loop_iterator)
    id_b = shuffle_vec(particle_pair_loop_iterator + 1)
    
    ! Read u, g in lab frame
    u_a_lab = spec%p(:,id_a)
    u_b_lab = spec%p(:,id_b)
    g_a_lab = sqrt( 1.0_p_k_part + dot_product( u_a_lab, u_a_lab ) )
    g_b_lab = sqrt( 1.0_p_k_part + dot_product( u_b_lab, u_b_lab ) )
    
    ! Relative velocity
    call lorentz_transform_proper( u_a_lab, u_b_lab, g_b_lab, u_rel, g_rel)
    au_rel = sqrt( dot_product( u_rel, u_rel ) )
    av_rel = au_rel/ g_rel
    
    ! Cycle if au_rel too small (au_rel**3 < tiny, plus one extra power for safety)
    ! (Could be made more efficient by comparing to the 4th root of zero,
    ! (but would need a select case switch for the precision )
    if( au_rel**4 <= TINY(au_rel) ) cycle
    
    if ( self%CrossCorrection == p_RelativisticCrossSectionCorrection_Rejection ) then
      ! Correction for relativistic crosssection
      ! c.f. Peano 2009
      if( 2.0_p_k_part * local_rng % genrand_real2() >=    &
          (1.0_p_k_part - dot_product(u_a_lab,u_b_lab)/(g_a_lab*g_b_lab)) ) cycle
    endif
    
    ! Calculate center-of-momentum 4-velocity
    ucm =  ( u_a_lab + u_b_lab ) / &
          sqrt( 2.0_p_k_part * ( 1.0_p_k_part + ( g_a_lab*g_b_lab - dot_product(u_a_lab,u_b_lab) ) ) )
    gcm = sqrt(1.0_p_k_part + dot_product(ucm,ucm))
    
    ! Transform u_i to the com frame
    call lorentz_transform_proper( ucm, u_a_lab, g_a_lab, u_in_cm, g_in_cm )
    au_in_cm = sqrt( dot_product( u_in_cm, u_in_cm ) )
    
    ! This is code from the original module
    ! At this point there's no reference for it and we recommend p_FrameCorrection_None
    select case ( self%FrameCorrection )
      case ( p_FrameCorrection_None )
        ! c.f. e.g. Sentoku 1998 eq. 2.6
        s = (s_first * g_rel) / (au_rel**3)
      case ( p_FrameCorrection_dt )
        ! c.f. nothing ???
        s = (g_b_lab * s_first * g_rel) / (au_rel**3)
      case ( p_FrameCorrection_dt_n )
        ! c.f. nothing ???
        s = (g_b_lab**2 * s_first * g_rel) / (au_rel**3)
      case default
        ERROR("Somehow frame correction was not set")
        ERROR("This is a bug")
        call abort_program(p_err_invalid)
    end select
    
    ! Calculate theta based on collision theory ( though actually cos(theta) )
    select case (self%CollisionModel)
      
      case( p_CollisionModel_Nanbu, p_CollisionModel_Nanbu_Rel )
        
        call local_rng%harvest_real2( random_unit_interval )
         
        if( s < p_useA_min ) then
          ! c.f. Nanbu, discussion just after eq. 17
          ctheta = 1.0_p_k_part + s * log(random_unit_interval)
          ctheta = max(-1.0_p_k_part, ctheta)
        elseif( s > p_useA_max ) then
          ! Isotropic scattering
          ! c.f. ibid
          ctheta = 1.0_p_k_part - 2.0_p_k_part * random_unit_interval
        else
          
          select case ( self%RootFinder )
          
            case(p_RootFinder_Newton)
              A_Nanbu = solve_nanbu_equation(s)
              
            case(p_RootFinder_Linear)
              if( s < p_lookup_min ) then
                A_Nanbu = 1.0_p_k_part / s
              elseif( s > p_lookup_max ) then
                A_Nanbu = 3.0_p_k_part * exp( -s )
              else
                A_Nanbu = A_lin_inter(s)
              endif
              
            case(p_RootFinder_Cubic)
              if( s < p_lookup_min ) then
                A_Nanbu = 1.0_p_k_part / s
              elseif( s > p_lookup_max ) then
                A_Nanbu = 3.0_p_k_part * exp( -s )
              else
                A_Nanbu = A_cubic_inter(s)
              endif
              
            case default
              ERROR("Somehow no root finder was set")
              ERROR("This is a bug")
              call abort_program(p_err_invalid)
              
          end select
          
          ! c.f. Nanbu 1997 eq. 17
          ctheta = log( exp(-A_Nanbu) + 2.0_p_k_part*random_unit_interval*sinh(A_Nanbu) ) / A_Nanbu
          ctheta = max(-1.0_p_k_part, ctheta)
          ctheta = min( 1.0_p_k_part, ctheta)
        endif
        
        if( self%CollisionModel == p_CollisionModel_Nanbu_Rel ) then
          theta_opr = acos( ctheta )
          ! c.f. Sentoku 2008 eq. 14
          ! ( This is, at best, only roughly the correct relativistic correction.
          ! ( This should be replaced with the calculation of Perez et al. 2012 .)
          theta_CoM = atan2( sin(theta_opr), g_in_cm*cos(theta_opr) - au_in_cm/av_rel )
          ctheta = cos( theta_CoM )
        endif
        
      case (p_CollisionModel_Sentoku)
        
        if ( self%CrossCorrection == p_RelativisticCrossSectionCorrection_Scattering ) then
          ! c.f. Peano 2009 ?
          ! Peano actually recommends using a rejection method (as above), currently this method isn't verified
          theta_opr = atan( real(local_rng%genrand_gaussian(), p_k_part) * sqrt(s) ) * &
                     ( 2.0_p_k_part - ( av_rel * (1.0_p_k_part - dot_product(u_a_lab,u_b_lab)/(g_a_lab*g_b_lab)) ) )
        else
          ! c.f Sentoku 2008 discussion after eq. 11
          ! (though really we should have a warning/error if this inequality is being hit often, as it means we're not getting the physics right either way)
          !s = min( s, 0.02_p_k_part )
          
          ! c.f. Sentoku 2008 eq. 11
          ! <{tan{theta/2}**2> = nu_ab * n_collide * dt = nu_ab * dt_s = s
          ! => tan(theta/2) = genrand_gaussian() * sqrt( s )
          ! => theta/2 = atan( genrand_gaussian() * sqrt( s ) )
          ! All this done in the One Particle at Rest frame
          theta_opr =  2.0_p_k_part * atan( real(local_rng%genrand_gaussian(), p_k_part) * sqrt( s ) )
        endif
        
        ! c.f. Sentoku 2008 eq. 14  (lower case 'cm' in that equation means 'velocity relative to the CoM', not 'velocity of the CoM')
        theta_CoM = atan2( sin(theta_opr), g_in_cm*cos(theta_opr) - au_in_cm/av_rel )
        ctheta = cos( theta_CoM )
        
      case( p_CollisionModel_Isotropic)
        
        ! c.f. http://mathworld.wolfram.com/SpherePointPicking.html
        ctheta = 2.0_p_k_part * real( local_rng%genrand_real2(), p_k_part ) - 1.0_p_k_part
        
      case default
        ERROR("Somehow no collision model was selected")
        ERROR("This is a bug")
        call abort_program( p_err_invalid)
    end select
    
    stheta = sqrt(1.0_p_k_part - ctheta**2)
    
    call local_rng%harvest_real2( random_unit_interval )
    phi = 2.0_p_k_part * random_unit_interval * real( pi, p_k_part )
    sphi = sin( phi )
    cphi = cos( phi )
    
    ! Now that theta, phi are found, calculate delta_p (in CoM frame, but without rotation)
    u_norm_a_cm = u_in_cm / au_in_cm !vector with length 1 and same direction as u_a_cm
    
    n_x = u_norm_a_cm(1)
    n_y = u_norm_a_cm(2)
    n_z = u_norm_a_cm(3)
    
    if (( n_x /= 0) .or. ( n_y /= 0)) then
      ! c.f. Sentoku 1998 eq.s 2.11-2.13
      ! Their derivation involves solving the matrix equation ibid eq. 2.5, then inverting that matrix
      ! (note Sentoku 2008 has a typo in this matrix, be sure to use Sentoku 1998)
      n_t = sqrt( n_x**2 + n_y**2 )
      rn_t = 1.0_p_k_part / n_t
      
      delta(1) = n_x * rn_t * n_z * stheta * cphi - &
                 n_y * rn_t * stheta * sphi + &
                 n_x * ( ctheta - 1.0_p_k_part)
      delta(2) = n_y * rn_t * n_z * stheta * cphi + &
                 n_x * rn_t * stheta * sphi + &
                 n_y * ( ctheta - 1.0_p_k_part )
      delta(3) = n_z * ( ctheta - 1.0_p_k_part ) - n_t * stheta * cphi
    else
      delta(1) = stheta * cphi
      delta(2) = stheta * sphi
      delta(3) = ctheta - 1.0_p_k_part
    endif
    
    ! Find size of macro-particle again, to weight statistical collision probability
    w_a = abs(spec%q(id_a))
    w_b = abs(spec%q(id_b))
    w_m = max( w_a, w_b )
    call local_rng%harvest_real2( random_unit_interval )
     
    ! Add delta to proper velocities (in CoM frame) (based on statistics)
    if( random_unit_interval < w_b/w_m) then
      u_after_collision = u_in_cm + delta * au_in_cm
      call lorentz_transform_proper( -ucm, u_after_collision, g_in_cm, spec%p(:,id_a), devnull)
    endif
    
    if( random_unit_interval < w_a/w_m) then
      u_after_collision = -u_in_cm - delta * au_in_cm
      call lorentz_transform_proper( -ucm, u_after_collision, g_in_cm, spec%p(:,id_b), devnull)
    endif
    
  enddo !all particles in cell
  
end subroutine like_collide
!---------------------------------------------------

!---------------------------------------------------
subroutine unlike_collide( &
                           self, &
                           sp1, &
                           sp2, &
                           sp_id1, &
                           sp_id2, &
                           coll_cell_id0, &
                           num_par, &
                           shuffle_vec, &
                           cell_fluid_state, &
                           l_debye, &
                           dt_n, &
                           local_rng &
                          )
!---------------------------------------------------
!       update momenta for collision
!---------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  type( t_collisions ), intent(in) :: self
  class ( t_species ), pointer :: sp1, sp2
  integer, intent(in) :: sp_id1, sp_id2
  ! Currently the code collides the shuffled b_vector off the a_species in order
  ! If we were to collide the shuffled b_vector off the shuffled a_vector I think we wouldn't need coll_cell_id0
  ! Doesn't make much difference though, I don't think
  integer, dimension(:), intent(in) :: coll_cell_id0
  integer, dimension(:), intent(in) :: num_par
  integer, dimension(:,:), target :: shuffle_vec
  type( t_local_fluid_state ), dimension(:), intent(in) :: cell_fluid_state
  ! total debye length [units of c/w0]
  real(p_k_part), intent(in) :: l_debye
  ! collision cycle: dt * n_collide [units of 1/w0]
  real(p_k_part), intent(in) :: dt_n
  class (t_random), intent(inout) :: local_rng
  
  ! Local variables
  
  ! Aliases for input parameters (depending on which species has more particles to collide)
  class ( t_species ), pointer :: spec_a, spec_b
  integer :: coll_cell_id0_a, coll_cell_id0_b
  integer :: num_par_a, num_par_b
  integer, dimension(:), pointer :: shuffle_vec_a, shuffle_vec_b
  type(t_local_fluid_state) :: cell_fluid_state_a, cell_fluid_state_b
  
  ! Masses [units of m_e]
  real(p_k_part)   :: m_real_a, m_real_b, rq_real_a, rq_real_b
  real(p_k_part)   :: mu_ab ! Reduced mass
  
  ! Weight of particle a, particle b, maximum of those (used for statistical scattering)
  real(p_k_part) :: w_a, w_b, w_m
  ! The sum of all w_m for each pair, used to normalize weighted averages
  ! c.f. Nanbu 1998, first equation after eq. 12
  real(p_k_part) :: n_alpha_beta
  
  ! Random number between zero and one
  real(p_k_part) :: random_unit_interval
  
  ! Coulomb logarithm, and its argument
  real(p_k_part) :: lambda, ln_lambda
  
  ! Average of u_rel**2 (used in log(lambda))
  real(p_k_part) :: u_sq_avg
  
  ! Collision 'probability'
  real(p_k_part) :: s
  
  ! 'Base' collision 'probability'
  real(p_k_part) :: s_first
  
  ! Adjusted timestep (dt * n_collide * {weighting for statistics} * {possible factor of 2} )
  real(p_k_part) :: dt_s
  
  ! The value needed for the Nanbu collision models
  ! coth(A) - 1/A = exp(-s)
  ! 0 < A_Nanbu < inf
  ! A_Nanbu( s -> inf ) = 3 * exp(-s)
  ! A_Nanbu( s -> 0 )   = 1 / s
  ! A_Nanbu( s <= 0 ) invalid
  real(p_k_part)   :: A_Nanbu
  
  ! Reference plasma frequency (just an alias to self%omega_pe_0)
  real(p_k_part) :: omega_pe_0
  
  ! Proper 4-velocity of particles in LAB frame
  real(p_k_part), dimension(p_p_dim) :: u_a_lab, u_b_lab
  real(p_k_part)                     :: g_a_lab, g_b_lab
  
  ! Proper 4-velocity of either particle in the One Particle at Rest (OPR) frame
  ! Also referred to as the velocity of part_a RELative to part_b (or vice versa)
  real(p_k_part), dimension(p_p_dim) :: u_rel
  real(p_k_part)                     :: g_rel
  ! Magnitude of space-like proper-velocity, true velocity of a relative to b (or vice versa)
  real(p_k_part)                     :: au_rel, av_rel
  
  ! Proper 4-velocity of the Center-of-Momentum (CoM) frame (as seen in the LAB Frame)
  real(p_k_part), dimension(p_p_dim) :: ucm
  real(p_k_part)                     :: gcm
  
  ! Proper 4-velocities of the particles as seen in the CoM frame
  ! Note this is therefore also the velocity *OF* the CoM frame as seen in the respective OPR frame
  ! (The frame of the more massive particle is the one we boost out of below)
  real(p_k_part), dimension(p_p_dim) :: u_a_in_cm, u_b_in_cm
  real(p_k_part)                     :: g_a_in_cm, g_b_in_cm
  ! Magnitude of the proper velocities
  real(p_k_part)                     :: au_a_in_cm, au_b_in_cm
  
  ! Scattering angles
  ! 0<phi<2pi, with equal distribution 
  ! -inf < tan(theta_CoM) < inf, with a gaussian distribution and variance given by the collision frequency
  real(p_k_part)   :: phi, theta_CoM
  
  ! Trigonometric functions of argument the scattering angles
  ! ( cos(phi), sin(phi), cos(theta_CoM), sin(theta_CoM) )
  real(p_k_part) :: cphi, sphi, ctheta, stheta
  
  ! theta in the One Particle at Rest frame
  real(p_k_part)                   :: theta_opr
  
  ! gamma factor and magnitude of the space-like proper-velocity of the One Particle at Rest frame
  ! as seen from the CoM Frame, where the 'One Particle' is the more massive one
  real(p_k_part)                   :: au_opr, g_opr
  
  ! Variables used for actual momenta update
  !
  ! Normalized CoM proper velocity
  ! (I think this may be unneeded but I have to check)
  real(p_k_part), dimension(p_p_dim) :: u_norm_a_cm
  
  ! The components of the normalized CoM proper velocity
  ! ( n_t = sqrt(n_x**2 + n_y**2) i.e. "n_transverse"; rn_t = 1/n_t)
  real(p_k_part) :: n_x, n_y, n_z, n_t, rn_t
  
  ! The change in momentum (p_new - p_old), still normalized
  real(p_k_part), dimension(p_p_dim) :: delta
  
  ! Updated space-like proper velocity after collision, in CoM
  real(p_k_part), dimension(p_p_dim) :: u_after_collision
  
  ! Numerical variables ( indexing )
  integer :: id_a, id_b, i
  
  ! Self-explanatory
  real(p_k_part) :: devnull
  
  ! Executable statements
  
  ! Return if no pairs are in cell
  if( (num_par(sp_id1) < 1) .or. (num_par(sp_id2) < 1) ) return
  
  ! Order the two species by the number of sim. particles (num_par_a >= num_par_b)
  if ( num_par(sp_id1) .ge. num_par(sp_id2) ) then
    !
    spec_a => sp1
    coll_cell_id0_a = coll_cell_id0(       sp_id1 )
    num_par_a = num_par(                   sp_id1 )
    shuffle_vec_a => shuffle_vec(    :   , sp_id1 ) ! This is never used
    cell_fluid_state_a = cell_fluid_state( sp_id1 )
    !
    spec_b => sp2
    coll_cell_id0_b = coll_cell_id0(       sp_id2 ) ! This is never used
    num_par_b = num_par(                   sp_id2 )
    shuffle_vec_b => shuffle_vec(    :   , sp_id2 )
    cell_fluid_state_b = cell_fluid_state( sp_id2 )
  else
    spec_a => sp2
    coll_cell_id0_a = coll_cell_id0(       sp_id2 )
    num_par_a = num_par(                   sp_id2 )
    shuffle_vec_a => shuffle_vec(    :   , sp_id2 ) ! This is never used
    cell_fluid_state_a = cell_fluid_state( sp_id2 )
    !
    spec_b => sp1
    coll_cell_id0_b = coll_cell_id0(       sp_id1 ) ! This is never used
    num_par_b = num_par(                   sp_id1 )
    shuffle_vec_b => shuffle_vec(    :   , sp_id1 )
    cell_fluid_state_b = cell_fluid_state( sp_id1 )
  endif
  
  m_real_a  = spec_a%m_real
  m_real_b  = spec_b%m_real
  rq_real_a = spec_a%rq_real
  rq_real_b = spec_b%rq_real
  ! Reduced mass (not relativistically correct...)
  mu_ab = ( m_real_a * m_real_b ) / ( m_real_a + m_real_b )
  omega_pe_0 = self%omega_pe_0
  
  u_sq_avg = 0
  n_alpha_beta = 0.0_p_k_part
  if( self%coulomb_logarithm_automatic ) then
    do i=1, num_par_a
      id_a = coll_cell_id0_a + i
      id_b = shuffle_vec_b(i) ! b has less particles, but shuffle_vec is periodic
      
      ! Collect sn_ab
      w_m = min( spec_a%q(id_a)*rq_real_a, spec_b%q(id_b)*rq_real_b )
      n_alpha_beta = n_alpha_beta + w_m
      
      !get u_a, u_b, g_a, g_b
      u_a_lab = spec_a%p(:,id_a)
      u_b_lab = spec_b%p(:,id_b)
      g_b_lab = sqrt( 1.0_p_k_part + dot_product( u_b_lab, u_b_lab ) )
      
      !get au_rel (transform a_b to ref frame of a)
      call lorentz_transform_proper( u_a_lab, u_b_lab, g_b_lab, u_rel, devnull)
      
      u_sq_avg = u_sq_avg + w_m * dot_product(u_rel,u_rel)
    enddo !all a
    u_sq_avg = u_sq_avg / n_alpha_beta
    
  else ! Nothing to do but find the weighting
    do i=1, num_par_a
      id_a = coll_cell_id0_a + i
      id_b = shuffle_vec_b(i)
      n_alpha_beta = n_alpha_beta + min( spec_a%q(id_a)*rq_real_a, spec_b%q(id_b)*rq_real_b )
    enddo
  endif
  
  if( self%coulomb_logarithm_automatic ) then
    
    ! Except for the q's this is the same as the like-like code
    ! Josh *thinks* that's right but, but if someone wants to extra-think that Josh wouldn't be offended
    lambda = min( ( mu_ab*cgs_me * l_debye*(cgs_c/omega_pe_0) * u_sq_avg * cgs_c**2 ) / &
                  ( 2.0_p_k_part * abs( spec_a%q_real * spec_b%q_real ) * cgs_e**2 ), &
                  ( mu_ab*cgs_me * l_debye*(cgs_c/omega_pe_0) * sqrt(u_sq_avg) * cgs_c ) / &
                  ( cgs_hbar ) )
    
    ln_lambda = 1./2 * ( 1./(1 + lambda**2) - 1 + log(1 + lambda**2) )
    
    ln_lambda = max( ln_lambda, 2._p_k_part)
    
  else
    ln_lambda = self%coulomb_logarithm_value
  endif
  
  ! Adjust time-step for collision statistics for unequal particle weights
  ! c.f. Nanbu 1998, second equation after eq. 12
  dt_s = dt_n * ((rq_real_a * cell_fluid_state_a%q_tot) / n_alpha_beta )
  if (self%CrossCorrTime ) then
    dt_s = 2.0_p_k_part * dt_s
  endif
  
  !s_first = ln_lambda * ( (spec_a%q_real * spec_b%q_real) / mu_ab )**2 * &
  !          cell_fluid_state_b%number_density * &
  !          dt_s * &
  !          ((cgs_e**2 * omega_pe_0)/(cgs_me * cgs_c**3))
  
  s_first = ln_lambda * &
            ( spec_a%q_real * spec_b%q_real / mu_ab )**2 * &
            cell_fluid_state_b%q_tot * rq_real_b / self%collision_cell_size * omega_pe_0**2  * & ! 'number density'
            cell_fluid_state_a%q_tot * rq_real_a / n_alpha_beta * dt_n/omega_pe_0 * & ! 'timestep'
            cgs_e**2 /(cgs_me * cgs_c**3)
  
  if (self%CrossCorrTime ) then
    s_first = 2.0_p_k_part * s_first
  endif
  
  ! Loop over particles 'a' and 'collide' them with particles 'b'
  do i=1, num_par_a
    ! Get particle indexes
    id_a = coll_cell_id0_a + i
    id_b = shuffle_vec_b(i) ! b has less particles, but shuffle_vec is periodic
    
    ! Read u, g in lab frame
    u_a_lab = spec_a%p(:,id_a)
    u_b_lab = spec_b%p(:,id_b)
    g_a_lab = sqrt( 1.0_p_k_part + dot_product( u_a_lab, u_a_lab ) )
    g_b_lab = sqrt( 1.0_p_k_part + dot_product( u_b_lab, u_b_lab ) )
    
    ! Relative velocity:
    call lorentz_transform_proper( u_a_lab, u_b_lab, g_b_lab, u_rel, g_rel)
    au_rel = sqrt( dot_product( u_rel, u_rel ) )
    av_rel = au_rel/ g_rel
    
    ! Cycle if au_rel too small (au_rel**3 < tiny, plus one extra power for safety)
    if( au_rel**4 <= TINY(au_rel) ) cycle
    
    if ( self%CrossCorrection == p_RelativisticCrossSectionCorrection_Rejection ) then
    ! Correction for relativistic crosssection
      if( 2.0_p_k_part * local_rng % genrand_real2() >=      &
          (1.0_p_k_part - dot_product(u_a_lab,u_b_lab)/(g_a_lab*g_b_lab)) ) cycle
    endif              
    
    ! Calculate Center-of-Momentum 4-velocity
    ucm = (m_real_a * u_a_lab + m_real_b * u_b_lab) / &
          sqrt( m_real_a**2 + m_real_b**2 + 2.0_p_k_part * &
                m_real_a * m_real_b * ( g_a_lab*g_b_lab - dot_product(u_a_lab,u_b_lab) ))
    gcm = sqrt(1.0_p_k_part + dot_product(ucm,ucm))
    
    ! Transform u_i to the CoM frame
    ! m_real_a * u_a_in_cm = -m_real_b * u_b_in_cm
    call lorentz_transform_proper( ucm, u_a_lab, g_a_lab, u_a_in_cm, g_a_in_cm )
    call lorentz_transform_proper( ucm, u_b_lab, g_b_lab, u_b_in_cm, g_b_in_cm )
    
    ! Calculate magnitude of the proper velocities in the CoM frame
    au_a_in_cm = sqrt( dot_product( u_a_in_cm, u_a_in_cm ) )
    au_b_in_cm = sqrt( dot_product( u_b_in_cm, u_b_in_cm ) )
    
    select case ( self%FrameCorrection )
    
      case ( p_FrameCorrection_None )
        s = (s_first * g_rel) / (au_rel**3)
        
      case ( p_FrameCorrection_dt )
        if( m_real_b > m_real_a ) then
          s = (g_b_lab * s_first * g_rel) / (au_rel**3)
        else
          s = (g_a_lab * s_first * g_rel) / (au_rel**3)
        endif
        
      case ( p_FrameCorrection_dt_n )
        if( m_real_b > m_real_a ) then
          s = (g_b_lab**2 * s_first * g_rel) / (au_rel**3)
        else
          s = (g_a_lab**2 * s_first * g_rel) / (au_rel**3)
        endif
        
      case default
        ERROR("Somehow frame correction was not set")
        ERROR("This is a bug")
        call abort_program(p_err_invalid)
        
    end select
    
    ! Calculate theta based on collision theory ( though actually cos(theta) )
    select case (self%CollisionModel)
      
      case( p_CollisionModel_Nanbu, p_CollisionModel_Nanbu_Rel )
        
        call local_rng % harvest_real2( random_unit_interval )
        
        if( s < p_useA_min ) then
          ctheta = 1.0_p_k_part + s * log(random_unit_interval)
          ctheta = max(-1.0_p_k_part, ctheta)
        elseif( s > p_useA_max ) then
          ctheta = 1.0_p_k_part - 2.0_p_k_part * random_unit_interval
        else
          
          select case( self%RootFinder )
          
            case(p_RootFinder_Newton)
              A_Nanbu = solve_nanbu_equation(s)
              
            case(p_RootFinder_Linear)
              if( s < p_lookup_min ) then
                A_Nanbu = 1.0_p_k_part / s
              elseif( s > p_lookup_max ) then
                A_Nanbu = 3.0_p_k_part * exp( -s )
              else
                A_Nanbu = A_lin_inter(s)
              endif
              
            case(p_RootFinder_Cubic)
              if( s < p_lookup_min ) then
                A_Nanbu = 1.0_p_k_part / s
              elseif( s > p_lookup_max ) then
                A_Nanbu = 3.0_p_k_part * exp( -s )
              else
                A_Nanbu = A_cubic_inter(s)
              endif
              
            case default
              ERROR("Somehow no root finder was set")
              ERROR("This is a bug")
              call abort_program(p_err_invalid)
              
          end select
          
          ctheta = log( exp(-A_Nanbu) + 2.0_p_k_part*random_unit_interval*sinh(A_Nanbu) ) / A_Nanbu
          ctheta = max(-1.0_p_k_part, ctheta)
          ctheta = min( 1.0_p_k_part, ctheta)
        endif
        
        if( self%CollisionModel == p_CollisionModel_Nanbu_Rel ) then
          if( m_real_b > m_real_a ) then
            au_opr = au_b_in_cm
            g_opr = g_b_in_cm
          else
            au_opr = au_a_in_cm
            g_opr = g_a_in_cm
          endif
          
          theta_opr = acos( ctheta )
          ! c.f. Sentoku 2008 eq. 14
          theta_CoM = atan2( sin(theta_opr), g_opr*cos(theta_opr) - au_opr/av_rel )
          ctheta = cos( theta_CoM )
          
        endif
        
      case (p_CollisionModel_Sentoku)
        
        ! The boost from OPR to CM isn't unique if m_a != m_b
        ! That is to say, it's not clear which OPR frame the scattering angle is calculated in
        ! Here we assume the OPR frame is the rest frame of the more massive particle, 
        ! as that particle's frame will change the least in the collision
        if( m_real_b > m_real_a ) then
          au_opr = au_b_in_cm
          g_opr = g_b_in_cm
        else
          au_opr = au_a_in_cm
          g_opr = g_a_in_cm
        endif
        
        if (self%CrossCorrection == p_RelativisticCrossSectionCorrection_Scattering) then
          theta_opr = atan( real(local_rng%genrand_gaussian(), p_k_part)  * sqrt(s) ) * &
                     ( 2.0_p_k_part - ( av_rel * (1.0_p_k_part - dot_product(u_a_lab,u_b_lab)/(g_a_lab*g_b_lab)) ) )
        else
          !s = min( s, 0.02_p_k_part )
          theta_opr =  2.0_p_k_part * atan( real(local_rng%genrand_gaussian(), p_k_part) * sqrt( s ) )
        endif
        
        ! c.f. Sentoku 2008 eq. 14
        theta_CoM = atan2( sin(theta_opr), g_opr*cos(theta_opr) - au_opr/av_rel )
        ctheta = cos( theta_CoM )
        
      case( p_CollisionModel_Isotropic)
        ctheta = 2.0_p_k_part * real(local_rng%genrand_real2(), p_k_part) - 1.0_p_k_part
        
      case default
        ERROR("Somehow no collision model was selected")
        ERROR("This is a bug")
        call abort_program( p_err_invalid)
        
    end select
    
    stheta = sqrt(1.0_p_k_part - ctheta**2)
    
    call local_rng%harvest_real2( random_unit_interval )
    phi = 2.0_p_k_part * random_unit_interval * real(pi,p_k_part)
    sphi = sin( phi )
    cphi = cos( phi )
    
    ! Now that theta, phi are found, calculate delta_p (in CoM frame, but without rotation)
    u_norm_a_cm = u_a_in_cm / au_a_in_cm ! Vector with length 1 and same direction as u_a_in_cm
    
    n_x = u_norm_a_cm(1)
    n_y = u_norm_a_cm(2)
    n_z = u_norm_a_cm(3)
    
    if (( n_x /= 0) .or. ( n_y /= 0)) then
      n_t = sqrt( n_x**2 + n_y**2 )
      rn_t = 1.0_p_k_part / n_t
      
      delta(1) = n_x * rn_t * n_z * stheta * cphi - &
                 n_y * rn_t * stheta * sphi + &
                 n_x * ( ctheta - 1.0_p_k_part)
      delta(2) = n_y * rn_t * n_z * stheta * cphi + &
                 n_x * rn_t * stheta * sphi + &
                 n_y * ( ctheta - 1.0_p_k_part )
      delta(3) = n_z * ( ctheta - 1.0_p_k_part ) - n_t * stheta * cphi
    else
      delta(1) = stheta * cphi
      delta(2) = stheta * sphi
      delta(3) = ctheta - 1.0_p_k_part
    endif
    
    ! Find size of macro-particle again, to weight statistical collision probability
    w_a = spec_a%q(id_a)*rq_real_a
    w_b = spec_b%q(id_b)*rq_real_b
    w_m = max( w_a, w_b )
    call local_rng%harvest_real2( random_unit_interval )
    
    ! Add delta to proper velocities (in CoM frame) (based on statistics)
    if( random_unit_interval < w_b/w_m) then
      u_after_collision = u_a_in_cm + delta * au_a_in_cm
      call lorentz_transform_proper( -ucm, u_after_collision, g_a_in_cm, spec_a%p(:,id_a), devnull)
    endif
    
    if( random_unit_interval < w_a/w_m) then
      u_after_collision = u_b_in_cm - delta * au_b_in_cm
      call lorentz_transform_proper( -ucm, u_after_collision, g_b_in_cm, spec_b%p(:,id_b), devnull)
    endif
    
  enddo ! All species_a particles in cell
  
end subroutine unlike_collide
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Strictly Non-Relativistic Collision Routines
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine like_collide_non_rel( &
                                 self, &
                                 spec, &
                                 num_par, &
                                 shuffle_vec, &
                                 cell_fluid_state, &
                                 l_debye, &
                                 dt_n, &
                                 local_rng &
                               )
!---------------------------------------------------
!       update momenta for collision
!---------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  type( t_collisions ), intent(in) :: self
  class ( t_species ), intent(inout) :: spec
  integer, intent(in) :: num_par
  integer, dimension(:) :: shuffle_vec
  type(t_local_fluid_state), intent(in) :: cell_fluid_state
  ! total debye length [c/w0]
  real(p_k_part), intent(in) :: l_debye
  ! collision cycle: dt * n_collide [1/w0]
  real(p_k_part), intent(in) :: dt_n
  class (t_random), intent(inout) :: local_rng
  
  ! Local variables
  
  ! Masses
  real(p_k_part)   :: m_real, rq_real ! Masses of real particles (in simulation units) (just aliases)
  real(p_k_part)   :: mu_ab ! Reduced mass (=m/2)
  
  ! Weight of particle a, particle b, maximum or minimum of those (used for statistical scattering)
  real(p_k_part) :: w_a, w_b, w_m
  ! The sum of all w_m for each pair, used to normalize weighted averages
  real(p_k_part) :: n_alpha_alpha
  
  ! Random number between zero and one ( [0,1) )
  real(p_k_part) :: random_unit_interval
  
  ! Coulomb logarithm, and its argument
  real(p_k_part) :: lambda, ln_lambda
  
  ! Average of u_rel**2 (used in log(lambda))
  real(p_k_part) :: u_sq_avg
  
  ! Collision 'probability'
  real(p_k_part) :: s
  
  ! 'Base' collision 'probability'
  ! i.e. s = s_first / u_rel**3
  real(p_k_part) :: s_first
  
  ! Takizuka and Abe equation 8a
  real(p_k_part) :: delta_ta
  
  ! Adjusted timestep (dt * n_collide * {weighting for statistics} )
  real(p_k_part) :: dt_s
  
  ! Reference plasma frequency (just an alias to self%omega_pe_0)
  real(p_k_part) :: omega_pe_0
  
  ! Relative velocity of the particles
  real(p_k_part), dimension(p_p_dim) :: u_rel
  real(p_k_part)                     :: au_rel
  
  ! Trigonometric functions of argument the scattering angles
  ! ( cos(phi), sin(phi), cos(theta), sin(theta) )
  real(p_k_part) :: cphi, sphi, ctheta, stheta
  
  ! Variables used for actual momenta update
  !
  ! Normalized velocity
  real(p_k_part), dimension(p_p_dim) :: u_norm
  
  ! The components of the normalized velocity
  ! ( n_t = sqrt(n_x**2 + n_y**2) i.e. "n_transverse"; rn_t = 1/n_t)
  real(p_k_part) :: n_x, n_y, n_z, n_t, rn_t
  
  ! The change in momentum (u_new - u_old), still normalized
  real(p_k_part), dimension(p_p_dim) :: delta
  
  ! Numerical variables ( indexing )
  ! Iterator used for looping over pairs in the shuffle vector
  ! Used in loop to calculate the statistical weights 
  ! and in the actual collision loop
  integer :: particle_pair_loop_iterator
  
  ! Index of last pair to collide (needed in case num_par is odd)
  integer :: last_pair
  
  ! Translation from shuffle vector array to particle array
  integer :: id_a, id_b
  
  ! Executable statements
  
  ! Return if no pairs are in cell
  if( num_par .lt. 2 ) return
  
  ! If there are an odd number of particles, collide the last with the first
  ! ( shuffle_vector( num_par + 1 ) = shuffle_vector(1) by construction )
  ! c.f. Nanbu 1998 Fig. 3(b)
  last_pair = num_par + mod( num_par, 2 )
  
  ! Aliases
  m_real  = spec%m_real
  rq_real = spec%rq_real
  omega_pe_0 = self%omega_pe_0
  mu_ab = m_real / 2.0_p_k_part
  
  ! Calculate statistical values in collision cell
  
  u_sq_avg = 0
  n_alpha_alpha = 0
  if( self%coulomb_logarithm_automatic ) then
    
    do particle_pair_loop_iterator = 1, last_pair, 2
      id_a = shuffle_vec(particle_pair_loop_iterator)
      id_b = shuffle_vec(particle_pair_loop_iterator + 1)
      
      ! Accumulate statistical weight
      w_m = min( abs(spec%q(id_a)), abs(spec%q(id_b)) ) 
      n_alpha_alpha = n_alpha_alpha + w_m  
      u_rel = spec%p(:,id_a) - spec%p(:,id_b)
      u_sq_avg = u_sq_avg + w_m * dot_product(u_rel,u_rel)
    enddo 
    u_sq_avg = u_sq_avg / n_alpha_alpha
    
  else ! Nothing to do but find the weighting
    
    do particle_pair_loop_iterator = 1, last_pair, 2
      id_a = shuffle_vec(particle_pair_loop_iterator)
      id_b = shuffle_vec(particle_pair_loop_iterator + 1)
      n_alpha_alpha = n_alpha_alpha + min( abs(spec%q(id_a)), abs(spec%q(id_b)) )
    enddo 
    
  endif
  
  if( self%coulomb_logarithm_automatic ) then
    
    lambda = ( mu_ab*cgs_me * l_debye*(cgs_c/omega_pe_0) * u_sq_avg * cgs_c**2 ) / &
             ( 2.0_p_k_part * spec%q_real**2 * cgs_e**2 )
    ln_lambda = log( lambda )
    ln_lambda = max( ln_lambda, 2._p_k_part )
    
  else
    ln_lambda = self%coulomb_logarithm_value
  endif
  
  ! Adjust time-step for collision statistics for unequal particle weights
  ! c.f. Nanbu 1998, eq.14 and the one immediately preceding it
  dt_s = dt_n * ( abs(cell_fluid_state%q_tot) / ( 2.0_p_k_part * n_alpha_alpha ) )
  
  ! prepare s calculations
  
  s_first = 1./2 * &
            omega_pe_0**2 * &                      ! Plasma frequency of reference density
            cell_fluid_state%number_density * &    ! Local density normalized to reference density
            ( spec%q_real**2 / mu_ab )**2 * &      ! Charge and reduced mass in normalized units
            ln_lambda * &
            cgs_e**2 / (cgs_me * cgs_c**3) * &
            dt_s/omega_pe_0                        ! Time step in cgs units
  
  ! Loop over pairs of particles and 'collide' them
  do particle_pair_loop_iterator = 1, last_pair, 2
    ! Get particle indexes
    id_a = shuffle_vec(particle_pair_loop_iterator)
    id_b = shuffle_vec(particle_pair_loop_iterator + 1)
    
    ! Relative velocity
    u_rel = spec%p(:,id_a) - spec%p(:,id_b)
    au_rel = sqrt( dot_product( u_rel, u_rel ) )
    if( au_rel**4 .le. TINY(au_rel) ) cycle ! Safety Cycle
    
    ! Calculate ctheta, stheta
    s = s_first / (au_rel**3)
    delta_TA = sqrt(s)*local_rng%genrand_gaussian()
    ctheta = 1 - 2*delta_TA**2/(1+delta_TA**2)
    stheta = 2*delta_TA/(1+delta_TA**2)
    
    call local_rng%harvest_real2( random_unit_interval )
    ! pi because sphi is always positive anyway (which is fine because theta is symmetric around zero)
    cphi = cos( random_unit_interval * real( pi, p_k_part ) )
    sphi = sqrt( 1. - cphi**2)
    
    ! Now that theta, phi are found, calculate delta_p
    u_norm = u_rel/au_rel
    n_x = u_norm(1)
    n_y = u_norm(2)
    n_z = u_norm(3)
    
    if (( n_x .ne. 0) .or. ( n_y .ne. 0)) then
      ! c.f. TA Equation 4
      
      n_t = sqrt( n_x**2 + n_y**2 )
      rn_t = 1.0_p_k_part / n_t
      
      delta(1) = n_x * rn_t * n_z * stheta * cphi - &
                 n_y * rn_t * stheta * sphi + &
                 n_x * ( ctheta - 1.0_p_k_part)
      delta(2) = n_y * rn_t * n_z * stheta * cphi + &
                 n_x * rn_t * stheta * sphi + &
                 n_y * ( ctheta - 1.0_p_k_part )
      delta(3) = n_z * ( ctheta - 1.0_p_k_part ) - n_t * stheta * cphi
    else
      delta(1) = stheta * cphi
      delta(2) = stheta * sphi
      delta(3) = ctheta - 1.0_p_k_part
    endif
    
    ! Find size of macro-particle again, to weight statistical collision probability
    w_a = abs(spec%q(id_a))
    w_b = abs(spec%q(id_b))
    w_m = max( w_a, w_b )
    call local_rng%harvest_real2( random_unit_interval )
     
    ! Add delta to proper velocities (based on statistics)
    if( random_unit_interval .lt. w_b/w_m) then
      spec%p(:,id_a) = spec%p(:,id_a) + au_rel * delta * mu_ab/m_real
    endif
    
    if( random_unit_interval .lt. w_a/w_m) then
      spec%p(:,id_b) = spec%p(:,id_b) - au_rel * delta * mu_ab/m_real
    endif
    
  enddo !all particles in cell
  
end subroutine like_collide_non_rel
!---------------------------------------------------

!---------------------------------------------------
subroutine unlike_collide_non_rel( &
                                   self, &
                                   sp1, &
                                   sp2, &
                                   sp_id1, &
                                   sp_id2, &
                                   coll_cell_id0, &
                                   num_par, &
                                   shuffle_vec, &
                                   cell_fluid_state, &
                                   l_debye, &
                                   dt_n, &
                                   local_rng &
                                 )
!---------------------------------------------------
!       update momenta for collision
!---------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  type( t_collisions ), intent(in) :: self
  class ( t_species ), pointer :: sp1, sp2
  integer, intent(in) :: sp_id1, sp_id2
  ! Currently the code collides the shuffled b_vector off the a_species in order
  ! If we were to collide the shuffled b_vector off the shuffled a_vector I think we wouldn't need coll_cell_id0
  ! Doesn't make much difference though, I don't think
  integer, dimension(:), intent(in) :: coll_cell_id0
  integer, dimension(:), intent(in) :: num_par
  integer, dimension(:,:), target :: shuffle_vec
  type( t_local_fluid_state ), dimension(:), intent(in) :: cell_fluid_state
  ! total debye length [units of c/w0]
  real(p_k_part), intent(in) :: l_debye
  ! collision cycle: dt * n_collide [units of 1/w0]
  real(p_k_part), intent(in) :: dt_n
  class (t_random), intent(inout) :: local_rng
  
  ! Local variables
  
  ! Aliases for input parameters (depending on which species has more particles to collide)
  class ( t_species ), pointer :: spec_a, spec_b
  integer :: coll_cell_id0_a, coll_cell_id0_b
  integer :: num_par_a, num_par_b
  integer, dimension(:), pointer :: shuffle_vec_a, shuffle_vec_b
  type(t_local_fluid_state) :: cell_fluid_state_a, cell_fluid_state_b
  
  ! Masses [units of m_e]
  real(p_k_part)   :: m_real_a, m_real_b, rq_real_a, rq_real_b
  real(p_k_part)   :: mu_ab ! Reduced mass
  
  ! Weight of particle a, particle b, maximum of those (used for statistical scattering)
  real(p_k_part) :: w_a, w_b, w_m
  ! The sum of all w_m for each pair, used to normalize weighted averages
  ! c.f. Nanbu 1998, first equation after eq. 12
  real(p_k_part) :: n_alpha_beta
  
  ! Random number between zero and one
  real(p_k_part) :: random_unit_interval
  
  ! Coulomb logarithm, and its argument
  real(p_k_part) :: lambda, ln_lambda
  
  ! Average of u_rel**2 (used in log(lambda))
  real(p_k_part) :: u_sq_avg
  
  ! Collision 'probability'
  real(p_k_part) :: s
  
  ! 'Base' collision 'probability'
  real(p_k_part) :: s_first
  
  real(p_k_part) :: delta_TA
  
  ! Adjusted timestep (dt * n_collide * {weighting for statistics} * {possible factor of 2} )
  real(p_k_part) :: dt_s
  
  ! Reference plasma frequency (just an alias to self%omega_pe_0)
  real(p_k_part) :: omega_pe_0
  
  ! Relative velocity
  real(p_k_part), dimension(p_p_dim) :: u_rel
  real(p_k_part)                     :: au_rel
  
  ! Trigonometric functions of argument the scattering angles
  ! ( cos(phi), sin(phi), cos(theta), sin(theta) )
  real(p_k_part) :: cphi, sphi, ctheta, stheta
  
  ! Variables used for actual momenta update
  !
  ! Normalized velocity
  real(p_k_part), dimension(p_p_dim) :: u_norm
  
  ! The components of the normalized velocity
  ! ( n_t = sqrt(n_x**2 + n_y**2) i.e. "n_transverse"; rn_t = 1/n_t)
  real(p_k_part) :: n_x, n_y, n_z, n_t, rn_t
  
  ! The change in momentum (p_new - p_old), still normalized
  real(p_k_part), dimension(p_p_dim) :: delta
  
  ! Numerical variables ( indexing )
  integer :: id_a, id_b, i
  
  ! Executable statements
  
  ! Return if no pairs are in cell
  if( (num_par(sp_id1) .lt. 1) .or. (num_par(sp_id2) .lt. 1) ) return
  
  ! Order the two species by the number of sim. particles (num_par_a >= num_par_b)
  if ( num_par(sp_id1) .ge. num_par(sp_id2) ) then
    !
    spec_a => sp1
    coll_cell_id0_a = coll_cell_id0(       sp_id1 )
    num_par_a = num_par(                   sp_id1 )
    shuffle_vec_a => shuffle_vec(    :   , sp_id1 ) ! This is never used
    cell_fluid_state_a = cell_fluid_state( sp_id1 )
    !
    spec_b => sp2
    coll_cell_id0_b = coll_cell_id0(       sp_id2 ) ! This is never used
    num_par_b = num_par(                   sp_id2 )
    shuffle_vec_b => shuffle_vec(    :   , sp_id2 )
    cell_fluid_state_b = cell_fluid_state( sp_id2 )
  else
    spec_a => sp2
    coll_cell_id0_a = coll_cell_id0(       sp_id2 )
    num_par_a = num_par(                   sp_id2 )
    shuffle_vec_a => shuffle_vec(    :   , sp_id2 ) ! This is never used
    cell_fluid_state_a = cell_fluid_state( sp_id2 )
    !
    spec_b => sp1
    coll_cell_id0_b = coll_cell_id0(       sp_id1 ) ! This is never used
    num_par_b = num_par(                   sp_id1 )
    shuffle_vec_b => shuffle_vec(    :   , sp_id1 )
    cell_fluid_state_b = cell_fluid_state( sp_id1 )
  endif
  
  m_real_a  = spec_a%m_real
  m_real_b  = spec_b%m_real
  rq_real_a = spec_a%rq_real
  rq_real_b = spec_b%rq_real
  ! Reduced mass
  mu_ab = ( m_real_a * m_real_b ) / ( m_real_a + m_real_b )
  omega_pe_0 = self%omega_pe_0
  
  u_sq_avg = 0
  n_alpha_beta = 0.0_p_k_part
  if( self%coulomb_logarithm_automatic ) then
    do i=1, num_par_a
      id_a = coll_cell_id0_a + i
      id_b = shuffle_vec_b(i) ! b has less particles, but shuffle_vec is periodic
      
      ! Collect sn_ab
      w_m = min( spec_a%q(id_a)*rq_real_a, spec_b%q(id_b)*rq_real_b )
      n_alpha_beta = n_alpha_beta + w_m
      
      u_rel = spec_a%p(:,id_a) - spec_b%p(:,id_b)
      u_sq_avg = u_sq_avg + w_m * dot_product(u_rel,u_rel)
    enddo !all a
    u_sq_avg = u_sq_avg / n_alpha_beta
    
  else ! Nothing to do but find the weighting
    do i=1, num_par_a
      id_a = coll_cell_id0_a + i
      id_b = shuffle_vec_b(i)
      n_alpha_beta = n_alpha_beta + min( spec_a%q(id_a)*rq_real_a, spec_b%q(id_b)*rq_real_b )
    enddo
  endif
  
  if( self%coulomb_logarithm_automatic ) then
    
    lambda = ( mu_ab*cgs_me * l_debye*(cgs_c/omega_pe_0) * u_sq_avg * cgs_c**2 ) / &
             ( 2.0_p_k_part * abs( spec_a%q_real * spec_b%q_real ) * cgs_e**2 )
    ln_lambda = log( lambda )
    ln_lambda = max( ln_lambda, 2._p_k_part)
    
  else
    ln_lambda = self%coulomb_logarithm_value
  endif
  
  ! Adjust time-step for collision statistics for unequal particle weights
  ! c.f. Nanbu 1998, second equation after eq. 12
  dt_s = dt_n * ((rq_real_a * cell_fluid_state_a%q_tot) / n_alpha_beta )
  
  s_first = 1./2 * &
            omega_pe_0**2 * &                                ! Plasma frequency of reference density
            cell_fluid_state_b%number_density * &            ! Local density of b species, normalized to reference density
            ( spec_a%q_real * spec_b%q_real / mu_ab )**2 * & ! Charge and reduced mass in normalized units
            ln_lambda * &
            cgs_e**2 / (cgs_me * cgs_c**3) * &
            dt_s/omega_pe_0                                  ! Time step in cgs units
  
  
  ! Loop over particles 'a' and 'collide' them with particles 'b'
  do i=1, num_par_a
    ! Get particle indexes
    id_a = coll_cell_id0_a + i
    id_b = shuffle_vec_b(i) ! b has less particles, but shuffle_vec is periodic
    
    ! Relative velocity:
    u_rel = spec_a%p(:,id_a) - spec_b%p(:,id_b)
    au_rel = sqrt( dot_product( u_rel, u_rel ) )
    if( au_rel**4 .le. TINY(au_rel) ) cycle ! Safety cycle
    
    ! Calculate ctheta, stheta
    s = s_first / (au_rel**3)
    delta_TA = sqrt(s)*local_rng%genrand_gaussian()
    ctheta = 1 - 2*delta_TA**2/(1+delta_TA**2)
    stheta = 2*delta_TA/(1+delta_TA**2)
    
    call local_rng%harvest_real2( random_unit_interval )
    cphi = cos( random_unit_interval * real( pi, p_k_part ) )
    sphi = sqrt( 1. - cphi**2)
    
    ! Now that theta, phi are found, calculate delta_p
    u_norm = u_rel / au_rel
    
    n_x = u_norm(1)
    n_y = u_norm(2)
    n_z = u_norm(3)
    
    if (( n_x .ne. 0) .or. ( n_y .ne. 0)) then
      n_t = sqrt( n_x**2 + n_y**2 )
      rn_t = 1.0_p_k_part / n_t
      
      delta(1) = n_x * rn_t * n_z * stheta * cphi - &
                 n_y * rn_t * stheta * sphi + &
                 n_x * ( ctheta - 1.0_p_k_part)
      delta(2) = n_y * rn_t * n_z * stheta * cphi + &
                 n_x * rn_t * stheta * sphi + &
                 n_y * ( ctheta - 1.0_p_k_part )
      delta(3) = n_z * ( ctheta - 1.0_p_k_part ) - n_t * stheta * cphi
    else
      delta(1) = stheta * cphi
      delta(2) = stheta * sphi
      delta(3) = ctheta - 1.0_p_k_part
    endif
    
    ! Find size of macro-particle again, to weight statistical collision probability
    w_a = spec_a%q(id_a)*rq_real_a
    w_b = spec_b%q(id_b)*rq_real_b
    w_m = max( w_a, w_b )
    call local_rng%harvest_real2( random_unit_interval )
    
    if( random_unit_interval .lt. w_b/w_m) then
      spec_a%p(:,id_a) = spec_a%p(:,id_a) + au_rel * delta * mu_ab/m_real_a
    endif
    
    if( random_unit_interval .lt. w_a/w_m) then
      spec_b%p(:,id_b) = spec_b%p(:,id_b) - au_rel * delta * mu_ab/m_real_b
    endif
    
  enddo ! All species_a particles in cell
  
end subroutine unlike_collide_non_rel
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Sentoku Collision Routines
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine like_collide_sentoku( &
                         self, &
                         spec, &
                         num_par, &
                         shuffle_vec, &
                         cell_fluid_state, &
                         l_debye, &
                         dt_n, &
                         local_rng &
                         )
!---------------------------------------------------
!       update momenta for collision
!---------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  type( t_collisions ), intent(in) :: self
  class ( t_species ), intent(inout) :: spec
  integer, intent(in) :: num_par
  integer, dimension(:) :: shuffle_vec
  type(t_local_fluid_state), intent(in) :: cell_fluid_state
  ! total debye length [c/w0]
  real(p_k_part), intent(in) :: l_debye
  ! collision cycle: dt * n_collide [1/w0]
  real(p_k_part), intent(in) :: dt_n
  class (t_random), intent(inout) :: local_rng
  
  ! Local variables
  
  ! Masses
  real(p_k_part)   :: m_real, rq_real ! Masses of real particles (in simulation units) (just aliases)
  real(p_k_part)   :: mu_ab ! Reduced mass (=m/2)
  
  ! Weight of particle a, particle b, maximum or minimum of those (used for statistical scattering)
  real(p_k_part) :: w_a, w_b, w_m
  ! The sum of all w_m for each pair, used to normalize weighted averages
  !
  real(p_k_part) :: n_alpha_alpha
  
  ! Random number between zero and one ( [0,1) )
  real(p_k_part) :: random_unit_interval
  
  ! Coulomb logarithm, and its argument
  real(p_k_part) :: lambda, ln_lambda
  
  ! Average of u_rel**2 (used in log(lambda))
  real(p_k_part) :: u_sq_avg
  
  ! Collision 'probability'
  real(p_k_part) :: s
  
  ! 'Base' collision 'probability'
  ! i.e. s = s_first / g_rel**2 beta_rel**3
  real(p_k_part) :: s_first
  
  ! Adjusted timestep (dt * n_collide * {weighting for statistics} )
  real(p_k_part) :: dt_s
  
  ! Reference plasma frequency (just an alias to self%omega_pe_0)
  real(p_k_part) :: omega_pe_0
  
  ! Proper 4-velocity of particles in LAB frame
  real(p_k_part), dimension(p_p_dim) :: u_a_lab, u_b_lab
  real(p_k_part)                     :: g_a_lab, g_b_lab
  
  ! Proper 4-velocity of either particle in the One Particle at Rest (OPR) frame
  ! Also referred to as the velocity of part_a RELative to part_b (or vice versa)
  real(p_k_part), dimension(p_p_dim) :: u_rel
  real(p_k_part)                     :: g_rel
  ! Magnitude of space-like proper-velocity, true velocity of a relative to b (or vice versa)
  real(p_k_part)                     :: au_rel, av_rel
  
  ! Proper 4-velocity of the Center-of-Momentum (CoM) frame (as seen in the LAB Frame)
  real(p_k_part), dimension(p_p_dim) :: ucm
  real(p_k_part)                     :: gcm
  
  ! Proper 4-velocity of either particle as seen in the CoM frame (equal because m_a = m_b)
  ! Note this is therefore also the velocity *OF* the CoM frame as seen in either OPR frame
  ! Note also the way it's calculated below, actually u_b_in_cm = -u_in_cm
  real(p_k_part), dimension(p_p_dim) :: u_in_cm
  real(p_k_part)                     :: g_in_cm
  ! Magnitude of the space-like proper velocity of either particle as seen in the CoM frame
  real(p_k_part)                     :: au_in_cm
  
  ! Scattering angles
  ! 0<phi<2pi, with equal distribution 
  ! -inf < tan(theta_CoM) < inf, with a gaussian distribution and variance given by the collision frequency
  real(p_k_part)   :: phi, theta_CoM
  
  ! Trigonometric functions of argument the scattering angles
  ! ( cos(phi), sin(phi), cos(theta_CoM), sin(theta_CoM) )
  real(p_k_part) :: cphi, sphi, ctheta, stheta
  
  ! theta in the One Particle at Rest frame
  real(p_k_part)  :: theta_opr
  
  ! Variables used for actual momenta update
  !
  real(p_k_part), dimension(p_p_dim) :: p_in_cm
  
  ! The components of the normalized CoM proper velocity
  ! ( n_t = sqrt(n_x**2 + n_y**2) i.e. "n_transverse"; rn_t = 1/n_t)
  real(p_k_part) :: n_x, n_y, n_z, n_t, rn_t, p
  
  ! The change in momentum (u_new - u_old)
  real(p_k_part), dimension(p_p_dim) :: delta
  
  ! Updated space-like proper velocity after collision, in CoM
  real(p_k_part), dimension(p_p_dim) :: u_after_collision
  
  ! Numerical variables ( indexing )
  ! Iterator used for looping over pairs in the shuffle vector
  ! Used in loop to calculate the statistical weights 
  ! and in the actual collision loop
  integer :: particle_pair_loop_iterator
  
  ! Index of last pair to collide (needed in case num_par is odd)
  integer :: last_pair
  
  ! Translation from shuffle vector array to particle array
  integer :: id_a, id_b
  
  ! Self-explanatory
  real(p_k_part) :: devnull
  
  ! Executable statements
  
  ! Return if no pairs are in cell
  if( num_par .lt. 2 ) return
  
  ! If there are an odd number of particles, collide the last with the first
  ! ( shuffle_vector( num_par + 1 ) = shuffle_vector(1) by construction )
  ! c.f. Nanbu 1998 Fig. 3(b)
  last_pair = num_par + mod( num_par, 2 )
  
  ! Aliases
  m_real  = spec%m_real
  rq_real = spec%rq_real
  omega_pe_0 = self%omega_pe_0
  ! It might be better to calculate the relativistic reduced mass (for each pair)
  mu_ab = m_real / 2.0_p_k_part
  
  ! Calculate statistical values in collision cell
  
  u_sq_avg = 0
  n_alpha_alpha = 0
  if( self%coulomb_logarithm_automatic ) then
    
    do particle_pair_loop_iterator = 1, last_pair, 2
      id_a = shuffle_vec(particle_pair_loop_iterator)
      id_b = shuffle_vec(particle_pair_loop_iterator + 1)
      
      ! Accumulate statistical weight
      w_m = min( abs(spec%q(id_a)), abs(spec%q(id_b)) ) 
      n_alpha_alpha = n_alpha_alpha + w_m  
      
      ! Get u_a_lab, u_b_lab, g_a_lab, g_b_lab
      u_a_lab = spec%p(:,id_a)
      u_b_lab = spec%p(:,id_b)
      g_b_lab = sqrt( 1.0_p_k_part + dot_product( u_b_lab, u_b_lab ) )
      
      ! Transform u_b from lab frame to ref frame of a
      call lorentz_transform_proper( u_a_lab, u_b_lab, g_b_lab, u_rel, devnull)
      u_sq_avg = u_sq_avg + w_m * dot_product(u_rel,u_rel)
    enddo 
    u_sq_avg = u_sq_avg / n_alpha_alpha
    
  else ! Nothing to do but find the weighting
    
    do particle_pair_loop_iterator = 1, last_pair, 2
      id_a = shuffle_vec(particle_pair_loop_iterator)
      id_b = shuffle_vec(particle_pair_loop_iterator + 1)
      n_alpha_alpha = n_alpha_alpha + min( abs(spec%q(id_a)), abs(spec%q(id_b)) )
    enddo 
    
  endif
  
  if( self%coulomb_logarithm_automatic ) then
    
    ! For the lower limit, use either the distance of classical closest approach, or the first Born approximation
    ! of the quantum mechanical result for close scattering (see Landau and Lifshitz 'Quantum Mechanics'
    ! chapter 127) depending on which is greater (and of course the Debye length
    ! for the upper limit.)
    ! (Note that this line is equivalent to a simple test of u < fine structure constant, but, well, that's not how it's written )
    ! (Note also, this should technically be u*v in the classical limit, but gamma < 1.000027 so it doesn't matter)
    lambda = min( ( mu_ab*cgs_me * l_debye*(cgs_c/omega_pe_0) * u_sq_avg * cgs_c**2 ) / &
                  ( 2.0_p_k_part * spec%q_real**2 * cgs_e**2 ), &
                  ( mu_ab*cgs_me * l_debye*(cgs_c/omega_pe_0) * sqrt(u_sq_avg) * cgs_c ) / &
                  ( cgs_hbar ) )
    
    ! c.f. Spitzer 'Physics of Fully Ionized Gases' equation 5-12
    ! $\int_0^p dx \frac{x^3}{\sqrt{1+x^2}} = \frac{1}{2}\left(\frac{1}{1+p^2}-1+ln(1+p^2)\right)$
    ln_lambda = 1./2 * ( 1./(1 + lambda**2) - 1 + log(1 + lambda**2) )
    
    ! Clamp log(lambda) at 2 (an arbitrary choice)
    ln_lambda = max( ln_lambda, 2._p_k_part )
    
  else
    ln_lambda = self%coulomb_logarithm_value
  endif
  
  ! Adjust time-step for collision statistics for unequal particle weights
  ! c.f. Nanbu 1998, eq.14 and the one immediately preceding it
  dt_s = dt_n * ( abs(cell_fluid_state%q_tot) / ( 2.0_p_k_part * n_alpha_alpha ) )
  
  ! prepare s calculations
  s_first = ( spec%q_real**2 / m_real )**2 * &
            omega_pe_0**2 * &
            ln_lambda * &
            cell_fluid_state%number_density * &
            cgs_e**2 / (cgs_me * cgs_c**3) * &
            (dt_s/omega_pe_0)  
  
  ! Loop over pairs of particles and 'collide' them
  do particle_pair_loop_iterator = 1, last_pair, 2
    ! Get particle indexes
    id_a = shuffle_vec(particle_pair_loop_iterator)
    id_b = shuffle_vec(particle_pair_loop_iterator + 1)
    
    ! Read u, g in lab frame
    u_a_lab = spec%p(:,id_a)
    u_b_lab = spec%p(:,id_b)
    g_a_lab = sqrt( 1.0_p_k_part + dot_product( u_a_lab, u_a_lab ) )
    g_b_lab = sqrt( 1.0_p_k_part + dot_product( u_b_lab, u_b_lab ) )
    
    ! Relative velocity
    call lorentz_transform_proper( u_a_lab, u_b_lab, g_b_lab, u_rel, g_rel)
    au_rel = sqrt( dot_product( u_rel, u_rel ) )
    av_rel = au_rel/ g_rel
    
    ! Cycle if au_rel too small (au_rel**3 < tiny, plus one extra power for safety)
    ! (Could be made more efficient by comparing to the 4th root of zero,
    ! (but would need a select case switch for the precision )
    if( au_rel**4 .le. TINY(au_rel) ) cycle
    
    ! Calculate center-of-momentum 4-velocity
    ucm =  ( u_a_lab + u_b_lab ) / &
          sqrt( 2.0_p_k_part * ( 1.0_p_k_part + ( g_a_lab*g_b_lab - dot_product(u_a_lab,u_b_lab) ) ) )
    gcm = sqrt(1.0_p_k_part + dot_product(ucm,ucm))
    
    ! Transform u_i to the com frame
    call lorentz_transform_proper( ucm, u_a_lab, g_a_lab, u_in_cm, g_in_cm )
    au_in_cm = sqrt( dot_product( u_in_cm, u_in_cm ) )
    
    s = (s_first * g_rel) / (au_rel**3)
    
    ! Calculate theta
          
    ! c.f. Sentoku 2008 eq. 11
    ! <{tan{theta/2}**2> = nu_ab * n_collide * dt = nu_ab * dt_s = s
    ! => tan(theta/2) = genrand_gaussian() * sqrt( s )
    ! => theta/2 = atan( genrand_gaussian() * sqrt( s ) )
    ! All this done in the One Particle at Rest frame
    theta_opr =  2.0_p_k_part * atan( real(local_rng%genrand_gaussian(), p_k_part) * sqrt( s ) )
    
    ! c.f. Sentoku 2008 eq. 14  (lower case 'cm' in that equation means 'velocity relative to the CoM', not 'velocity of the CoM')
    theta_CoM = atan2( sin(theta_opr), g_in_cm*cos(theta_opr) - au_in_cm/av_rel )
    ctheta = cos( theta_CoM )
    
    stheta = sqrt(1.0_p_k_part - ctheta**2)
    
    call local_rng%harvest_real2( random_unit_interval )
    phi = 2.0_p_k_part * random_unit_interval * real( pi, p_k_part )
    sphi = sin( phi )
    cphi = cos( phi )
    
    ! Now that theta, phi are found, calculate delta_p (in CoM frame, but without rotation)
    p_in_cm = u_in_cm * m_real
    
    n_x = p_in_cm(1)
    n_y = p_in_cm(2)
    n_z = p_in_cm(3)
    p = sqrt(n_x**2 + n_y**2 + n_z**2)
    
    if (( n_x .ne. 0) .or. ( n_y .ne. 0)) then
      ! c.f. Sentoku 1998 eq.s 2.11-2.13
      ! Their derivation involves solving the matrix equation ibid eq. 2.5, then inverting that matrix
      ! (note Sentoku 2008 has a typo in this matrix, be sure to use Sentoku 1998)
      n_t = sqrt( n_x**2 + n_y**2 )
      rn_t = 1.0_p_k_part / n_t
      
      delta(1) = n_x * rn_t * n_z * stheta * cphi - &
                 n_y * rn_t * p * stheta * sphi + &
                 n_x * ( ctheta - 1.0_p_k_part)
      delta(2) = n_y * rn_t * n_z * stheta * cphi + &
                 n_x * rn_t * p * stheta * sphi + &
                 n_y * ( ctheta - 1.0_p_k_part )
      delta(3) = n_z * ( ctheta - 1.0_p_k_part ) - n_t * stheta * cphi
    else
      delta(1) = p * stheta * cphi
      delta(2) = p * stheta * sphi
      delta(3) = p * (ctheta - 1.0_p_k_part)
    endif
    
    ! Find size of macro-particle again, to weight statistical collision probability
    w_a = abs(spec%q(id_a))
    w_b = abs(spec%q(id_b))
    w_m = max( w_a, w_b )
    call local_rng%harvest_real2( random_unit_interval )
     
    ! Add delta to proper velocities (in CoM frame) (based on statistics)
    if( random_unit_interval .lt. w_b/w_m) then
      u_after_collision = u_in_cm + delta / m_real
      call lorentz_transform_proper( -ucm, u_after_collision, g_in_cm, spec%p(:,id_a), devnull)
    endif
    
    if( random_unit_interval .lt. w_a/w_m) then
      u_after_collision = -u_in_cm - delta / m_real
      call lorentz_transform_proper( -ucm, u_after_collision, g_in_cm, spec%p(:,id_b), devnull)
    endif
    
  enddo !all particles in cell
  
end subroutine like_collide_sentoku
!---------------------------------------------------

!---------------------------------------------------
subroutine unlike_collide_sentoku( &
                           self, &
                           sp1, &
                           sp2, &
                           sp_id1, &
                           sp_id2, &
                           coll_cell_id0, &
                           num_par, &
                           shuffle_vec, &
                           cell_fluid_state, &
                           l_debye, &
                           dt_n, &
                           local_rng &
                          )
!---------------------------------------------------
!       update momenta for collision
!---------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  type( t_collisions ), intent(in) :: self
  class ( t_species ), pointer :: sp1, sp2
  integer, intent(in) :: sp_id1, sp_id2
  ! Currently the code collides the shuffled b_vector off the a_species in order
  ! If we were to collide the shuffled b_vector off the shuffled a_vector I think we wouldn't need coll_cell_id0
  ! Doesn't make much difference though, I don't think
  integer, dimension(:), intent(in) :: coll_cell_id0
  integer, dimension(:), intent(in) :: num_par
  integer, dimension(:,:), target :: shuffle_vec
  type( t_local_fluid_state ), dimension(:), intent(in) :: cell_fluid_state
  ! total debye length [units of c/w0]
  real(p_k_part), intent(in) :: l_debye
  ! collision cycle: dt * n_collide [units of 1/w0]
  real(p_k_part), intent(in) :: dt_n
  class (t_random), intent(inout) :: local_rng
  
  ! Local variables
  
  ! Aliases for input parameters (depending on which species has more particles to collide)
  class ( t_species ), pointer :: spec_a, spec_b
  integer :: coll_cell_id0_a, coll_cell_id0_b
  integer :: num_par_a, num_par_b
  integer, dimension(:), pointer :: shuffle_vec_a, shuffle_vec_b
  type(t_local_fluid_state) :: cell_fluid_state_a, cell_fluid_state_b
  
  ! Masses [units of m_e]
  real(p_k_part)   :: m_real_a, m_real_b, rq_real_a, rq_real_b
  real(p_k_part)   :: mu_ab ! Reduced mass
  real(p_k_part)   :: m_rel ! The relative mass (the mass associated with the 'relative momentum' {yeah I don't know why you would call it that either})
  
  ! Weight of particle a, particle b, maximum of those (used for statistical scattering)
  real(p_k_part) :: w_a, w_b, w_m
  ! The sum of all w_m for each pair, used to normalize weighted averages
  ! c.f. Nanbu 1998, first equation after eq. 12
  real(p_k_part) :: n_alpha_beta
  
  ! Random number between zero and one
  real(p_k_part) :: random_unit_interval
  
  ! Coulomb logarithm, and its argument
  real(p_k_part) :: lambda, ln_lambda
  
  ! Average of u_rel**2 (used in log(lambda))
  real(p_k_part) :: u_sq_avg
  
  ! Collision 'probability'
  real(p_k_part) :: s
  
  ! 'Base' collision 'probability'
  real(p_k_part) :: s_first

  ! Adjusted timestep (dt * n_collide * {weighting for statistics} )
  real(p_k_part) :: dt_s
  
  ! Reference plasma frequency (just an alias to self%omega_pe_0)
  real(p_k_part) :: omega_pe_0
  
  ! Proper 4-velocity of particles in LAB frame
  real(p_k_part), dimension(p_p_dim) :: u_a_lab, u_b_lab
  real(p_k_part)                     :: g_a_lab, g_b_lab
  
  ! Proper 4-velocity of either particle in the One Particle at Rest (OPR) frame
  ! Also referred to as the velocity of part_a RELative to part_b (or vice versa)
  real(p_k_part), dimension(p_p_dim) :: u_rel
  real(p_k_part)                     :: g_rel
  ! Magnitude of space-like proper-velocity, true velocity of a relative to b (or vice versa)
  real(p_k_part)                     :: au_rel, av_rel
  
  ! Proper 4-velocity of the Center-of-Momentum (CoM) frame (as seen in the LAB Frame)
  real(p_k_part), dimension(p_p_dim) :: ucm
  real(p_k_part)                     :: gcm
  
  ! Proper 4-velocities of the particles as seen in the CoM frame
  ! Note this is therefore also the velocity *OF* the CoM frame as seen in the respective OPR frame
  ! (The frame of the more massive particle is the one we boost out of below)
  real(p_k_part), dimension(p_p_dim) :: u_a_in_cm, u_b_in_cm
  real(p_k_part)                     :: g_a_in_cm, g_b_in_cm
  ! Magnitude of the proper velocities
  real(p_k_part)                     :: au_a_in_cm, au_b_in_cm
  
  ! Scattering angles
  ! 0<phi<2pi, with equal distribution 
  ! -inf < tan(theta_CoM) < inf, with a gaussian distribution and variance given by the collision frequency
  real(p_k_part)   :: phi, theta_CoM
  
  ! Trigonometric functions of argument the scattering angles
  ! ( cos(phi), sin(phi), cos(theta_CoM), sin(theta_CoM) )
  real(p_k_part) :: cphi, sphi, ctheta, stheta
  
  ! theta in the One Particle at Rest frame
  real(p_k_part)                   :: theta_opr
  
  ! gamma factor and magnitude of the space-like proper-velocity of the One Particle at Rest frame
  ! as seen from the CoM Frame, where the 'One Particle' is the more massive one
  real(p_k_part)                   :: au_opr, g_opr
  
  ! Variables used for actual momenta update
  !
  real(p_k_part), dimension(p_p_dim) :: p_in_cm
  
  ! The components of the CoM momentum
  ! ( n_t = sqrt(n_x**2 + n_y**2) i.e. "n_transverse"; rn_t = 1/n_t)
  real(p_k_part) :: n_x, n_y, n_z, n_t, rn_t, p
  
  ! The change in momentum (p_new - p_old)
  real(p_k_part), dimension(p_p_dim) :: delta
  
  ! Updated space-like proper velocity after collision, in CoM
  real(p_k_part), dimension(p_p_dim) :: u_after_collision
  
  ! Numerical variables ( indexing )
  integer :: id_a, id_b, i
  
  ! Self-explanatory
  real(p_k_part) :: devnull
  real(p_k_part) :: avg
  integer :: count
  avg =0.
  count = 0
  
  ! Executable statements
  
  ! Return if no pairs are in cell
  if( (num_par(sp_id1) .lt. 1) .or. (num_par(sp_id2) .lt. 1) ) return
  
  ! Order the two species by the number of sim. particles (num_par_a >= num_par_b)
  if ( num_par(sp_id1) .ge. num_par(sp_id2) ) then
    !
    spec_a => sp1
    coll_cell_id0_a = coll_cell_id0(       sp_id1 )
    num_par_a = num_par(                   sp_id1 )
    shuffle_vec_a => shuffle_vec(    :   , sp_id1 ) ! This is never used
    cell_fluid_state_a = cell_fluid_state( sp_id1 )
    !
    spec_b => sp2
    coll_cell_id0_b = coll_cell_id0(       sp_id2 ) ! This is never used
    num_par_b = num_par(                   sp_id2 )
    shuffle_vec_b => shuffle_vec(    :   , sp_id2 )
    cell_fluid_state_b = cell_fluid_state( sp_id2 )
  else
    spec_a => sp2
    coll_cell_id0_a = coll_cell_id0(       sp_id2 )
    num_par_a = num_par(                   sp_id2 )
    shuffle_vec_a => shuffle_vec(    :   , sp_id2 ) ! This is never used
    cell_fluid_state_a = cell_fluid_state( sp_id2 )
    !
    spec_b => sp1
    coll_cell_id0_b = coll_cell_id0(       sp_id1 ) ! This is never used
    num_par_b = num_par(                   sp_id1 )
    shuffle_vec_b => shuffle_vec(    :   , sp_id1 )
    cell_fluid_state_b = cell_fluid_state( sp_id1 )
  endif
  
  m_real_a  = spec_a%m_real
  m_real_b  = spec_b%m_real
  rq_real_a = spec_a%rq_real
  rq_real_b = spec_b%rq_real
  ! Reduced mass (not relativistically correct...)
  mu_ab = ( m_real_a * m_real_b ) / ( m_real_a + m_real_b )
  m_rel = min(m_real_a, m_real_b)
  omega_pe_0 = self%omega_pe_0
  
  u_sq_avg = 0
  n_alpha_beta = 0.0_p_k_part
  if( self%coulomb_logarithm_automatic ) then
    do i=1, num_par_a
      id_a = coll_cell_id0_a + i
      id_b = shuffle_vec_b(i) ! b has less particles, but shuffle_vec is periodic
      
      ! Collect sn_ab
      w_m = min( spec_a%q(id_a)*rq_real_a, spec_b%q(id_b)*rq_real_b )
      n_alpha_beta = n_alpha_beta + w_m
      
      !get u_a, u_b, g_a, g_b
      u_a_lab = spec_a%p(:,id_a)
      u_b_lab = spec_b%p(:,id_b)
      g_b_lab = sqrt( 1.0_p_k_part + dot_product( u_b_lab, u_b_lab ) )
      
      !get au_rel (transform a_b to ref frame of a)
      call lorentz_transform_proper( u_a_lab, u_b_lab, g_b_lab, u_rel, devnull)
      
      u_sq_avg = u_sq_avg + w_m * dot_product(u_rel,u_rel)
    enddo !all a
    u_sq_avg = u_sq_avg / n_alpha_beta
    
  else ! Nothing to do but find the weighting
    do i=1, num_par_a
      id_a = coll_cell_id0_a + i
      id_b = shuffle_vec_b(i)
      n_alpha_beta = n_alpha_beta + min( spec_a%q(id_a)*rq_real_a, spec_b%q(id_b)*rq_real_b )
    enddo
  endif
  
  if( self%coulomb_logarithm_automatic ) then
    
    ! Except for the q's this is the same as the like-like code
    ! Josh *thinks* that's right but, but if someone wants to extra-think that Josh wouldn't be offended
    lambda = min( ( mu_ab*cgs_me * l_debye*(cgs_c/omega_pe_0) * u_sq_avg * cgs_c**2 ) / &
                  ( 2.0_p_k_part * abs( spec_a%q_real * spec_b%q_real ) * cgs_e**2 ), &
                  ( mu_ab*cgs_me * l_debye*(cgs_c/omega_pe_0) * sqrt(u_sq_avg) * cgs_c ) / &
                  ( cgs_hbar ) )
    
    ln_lambda = 1./2 * ( 1./(1 + lambda**2) - 1 + log(1 + lambda**2) )
    
    ln_lambda = max( ln_lambda, 2._p_k_part)
    
  else
    ln_lambda = self%coulomb_logarithm_value
  endif
  
  ! Adjust time-step for collision statistics for unequal particle weights
  ! c.f. Nanbu 1998, second equation after eq. 12
  dt_s = dt_n * ((rq_real_a * cell_fluid_state_a%q_tot) / n_alpha_beta )
  
  s_first = ( spec_a%q_real * spec_b%q_real / m_rel )**2 * &
            omega_pe_0**2 * &
            ln_lambda * &
            cell_fluid_state_b%number_density * &
            cgs_e**2 /(cgs_me * cgs_c**3) * &
            dt_s / omega_pe_0
  
  ! Loop over particles 'a' and 'collide' them with particles 'b'
  do i=1, num_par_a
    ! Get particle indexes
    id_a = coll_cell_id0_a + i
    id_b = shuffle_vec_b(i) ! b has less particles, but shuffle_vec is periodic
    
    ! Read u, g in lab frame
    u_a_lab = spec_a%p(:,id_a)
    u_b_lab = spec_b%p(:,id_b)
    g_a_lab = sqrt( 1.0_p_k_part + dot_product( u_a_lab, u_a_lab ) )
    g_b_lab = sqrt( 1.0_p_k_part + dot_product( u_b_lab, u_b_lab ) )
    
    ! Relative velocity:
    call lorentz_transform_proper( u_a_lab, u_b_lab, g_b_lab, u_rel, g_rel)
    au_rel = sqrt( dot_product( u_rel, u_rel ) )
    av_rel = au_rel/ g_rel
    
    ! Cycle if au_rel too small (au_rel**3 < tiny, plus one extra power for safety)
    if( au_rel**4 .le. TINY(au_rel) ) cycle         
    
    ! Calculate Center-of-Momentum 4-velocity
    ucm = (m_real_a * u_a_lab + m_real_b * u_b_lab) / &
          sqrt( m_real_a**2 + m_real_b**2 + 2.0_p_k_part * &
                m_real_a * m_real_b * ( g_a_lab*g_b_lab - dot_product(u_a_lab,u_b_lab) ))
    gcm = sqrt(1.0_p_k_part + dot_product(ucm,ucm))
    
    ! Transform u_i to the CoM frame
    ! m_real_a * u_a_in_cm = -m_real_b * u_b_in_cm
    call lorentz_transform_proper( ucm, u_a_lab, g_a_lab, u_a_in_cm, g_a_in_cm )
    call lorentz_transform_proper( ucm, u_b_lab, g_b_lab, u_b_in_cm, g_b_in_cm )
    
    ! Calculate magnitude of the proper velocities in the CoM frame
    au_a_in_cm = sqrt( dot_product( u_a_in_cm, u_a_in_cm ) )
    au_b_in_cm = sqrt( dot_product( u_b_in_cm, u_b_in_cm ) )
    
    s = (s_first * g_rel) / (au_rel**3)
    
    ! Calculate theta based on collision theory ( though actually cos(theta) )
    if( m_real_b .gt. m_real_a ) then
      au_opr = au_b_in_cm
      g_opr = g_b_in_cm
    else
      au_opr = au_a_in_cm
      g_opr = g_a_in_cm
    endif
    
    theta_opr = 2.0_p_k_part * atan( real(local_rng%genrand_gaussian(), p_k_part) * sqrt( s ) )
    
    ! c.f. Sentoku 2008 eq. 14
    theta_CoM = atan2( sin(theta_opr), g_opr*cos(theta_opr) - au_opr/av_rel )
    ctheta = cos( theta_CoM )
    
    stheta = sqrt(1.0_p_k_part - ctheta**2)
    
    call local_rng%harvest_real2( random_unit_interval )
    phi = 2.0_p_k_part * random_unit_interval * real(pi,p_k_part)
    sphi = sin( phi )
    cphi = cos( phi )
    
    ! Now that theta, phi are found, calculate delta_p (in CoM frame, but without rotation)
    !u_norm_a_cm = u_a_in_cm / au_a_in_cm ! Vector with length 1 and same direction as u_a_in_cm
    p_in_cm = u_a_in_cm * m_real_a
    p = sqrt(p_in_cm(1)**2 + p_in_cm(2)**2 + p_in_cm(3)**2)
    
    n_x = p_in_cm(1)
    n_y = p_in_cm(2)
    n_z = p_in_cm(3)
    
    if (( n_x .ne. 0) .or. ( n_y .ne. 0)) then
      n_t = sqrt( n_x**2 + n_y**2 )
      rn_t = 1.0_p_k_part / n_t
      
      delta(1) = n_x * rn_t * n_z * stheta * cphi - &
                 n_y * rn_t * p * stheta * sphi + &
                 n_x * ( ctheta - 1.0_p_k_part)
      delta(2) = n_y * rn_t * n_z * stheta * cphi + &
                 n_x * rn_t * p * stheta * sphi + &
                 n_y * ( ctheta - 1.0_p_k_part )
      delta(3) = n_z * ( ctheta - 1.0_p_k_part ) - n_t * stheta * cphi
    else
      delta(1) = p * stheta * cphi
      delta(2) = p * stheta * sphi
      delta(3) = p * (ctheta - 1.0_p_k_part)
    endif
    
    ! Find size of macro-particle again, to weight statistical collision probability
    w_a = spec_a%q(id_a)*rq_real_a
    w_b = spec_b%q(id_b)*rq_real_b
    w_m = max( w_a, w_b )
    call local_rng%harvest_real2( random_unit_interval )
    
    ! Add delta to proper velocities (in CoM frame) (based on statistics)
    if( random_unit_interval .lt. w_b/w_m) then
      u_after_collision = u_a_in_cm + delta / m_real_a
      call lorentz_transform_proper( -ucm, u_after_collision, g_a_in_cm, spec_a%p(:,id_a), devnull)
    endif
    
    if( random_unit_interval .lt. w_a/w_m) then
      u_after_collision = u_b_in_cm - delta / m_real_b
      call lorentz_transform_proper( -ucm, u_after_collision, g_b_in_cm, spec_b%p(:,id_b), devnull)
    endif
    
  enddo ! All species_a particles in cell
  
end subroutine unlike_collide_sentoku
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Perez Collision Routines
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine like_collide_perez( &
                         self, &
                         spec, &
                         num_par, &
                         shuffle_vec, &
                         cell_fluid_state, &
                         l_debye, &
                         dt_n, &
                         local_rng &
                         )
!---------------------------------------------------
!       update momenta for collision
!---------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  type( t_collisions ), intent(in) :: self
  class ( t_species ), intent(inout) :: spec
  integer, intent(in) :: num_par
  integer, dimension(:) :: shuffle_vec
  type(t_local_fluid_state), intent(in) :: cell_fluid_state
  ! total debye length [c/w0]
  real(p_k_part), intent(in) :: l_debye
  ! collision cycle: dt * n_collide [1/w0]
  real(p_k_part), intent(in) :: dt_n
  class (t_random), intent(inout) :: local_rng
  
  ! Local variables
  
  ! Masses
  real(p_k_part)   :: m_real, rq_real ! Masses of real particles (in simulation units) (just aliases)
  real(p_k_part)   :: mu_ab ! Reduced mass (=m/2)
  
  ! Weight of particle a, particle b, maximum or minimum of those (used for statistical scattering)
  real(p_k_part) :: w_a, w_b, w_m
  ! The sum of all w_m for each pair, used to normalize weighted averages
  !
  real(p_k_part) :: n_alpha_alpha
  
  ! Random number between zero and one ( [0,1) )
  real(p_k_part) :: random_unit_interval
  
  ! Coulomb logarithm, and its argument
  real(p_k_part) :: lambda, ln_lambda
  
  ! Average of u_rel**2 (used in log(lambda))
  real(p_k_part) :: u_sq_avg
  
  ! Collision 'probability'
  real(p_k_part) :: s
  
  ! 'Base' collision 'probability'
  real(p_k_part) :: s_first
  
  ! Collision probability variables, for low temperature correction
  real(p_k_part) :: s_prime, s_prime_first
  ! Normal density for low temp correction (just from omega_pe_0)
  real(p_k_part) :: n_0
  
  ! Adjusted timestep (dt * n_collide * {weighting for statistics} )
  real(p_k_part) :: dt_s
  
  ! Reference plasma frequency (just an alias to self%omega_pe_0)
  real(p_k_part) :: omega_pe_0
  
  ! Proper 4-velocity of particles in LAB frame
  real(p_k_part), dimension(p_p_dim) :: u_a_lab, u_b_lab
  real(p_k_part)                     :: g_a_lab, g_b_lab
  
  ! Proper 4-velocity of either particle in the One Particle at Rest (OPR) frame
  ! Also referred to as the velocity of part_a RELative to part_b (or vice versa)
  real(p_k_part), dimension(p_p_dim) :: u_rel
  real(p_k_part)                     :: g_rel
  ! Magnitude of space-like proper-velocity, true velocity of a relative to b (or vice versa)
  real(p_k_part)                     :: au_rel, av_rel
  
  ! Proper 4-velocity of the Center-of-Momentum (CoM) frame (as seen in the LAB Frame)
  real(p_k_part), dimension(p_p_dim) :: ucm
  real(p_k_part)                     :: gcm
  
  ! Proper 4-velocity of either particle as seen in the CoM frame (equal because m_a = m_b)
  ! Note this is therefore also the velocity *OF* the CoM frame as seen in either OPR frame
  ! Note also the way it's calculated below, actually u_b_in_cm = -u_in_cm
  real(p_k_part), dimension(p_p_dim) :: u_in_cm
  real(p_k_part)                     :: g_in_cm
  ! Magnitude of the space-like proper velocity of either particle as seen in the CoM frame
  real(p_k_part)                     :: au_in_cm
  
  ! Nanbu parameter
  real(p_k_part) :: A
  
  ! Scattering angle phi
  !real(p_k_part)   :: phi
  
  ! Trigonometric functions of argument the scattering angles
  ! ( cos(phi), sin(phi), cos(theta), sin(theta) )
  ! (we never solve for theta)
  real(p_k_part) :: cphi, sphi, ctheta, stheta
  ! try skipping phi as well
  real(p_k_part) :: x_phi, y_phi
  
  ! Variables used for actual momenta update
  !
  real(p_k_part), dimension(p_p_dim) :: p_in_cm
  
  ! The components of the normalized CoM proper velocity
  ! ( n_t = sqrt(n_x**2 + n_y**2) i.e. "n_transverse"; rn_t = 1/n_t)
  real(p_k_part) :: n_x, n_y, n_z, n_t, rn_t, p
  
  ! The change in momentum (u_new - u_old)
  real(p_k_part), dimension(p_p_dim) :: delta
  
  ! Updated space-like proper velocity after collision, in CoM
  real(p_k_part), dimension(p_p_dim) :: u_after_collision
  
  ! Numerical variables ( indexing )
  ! Iterator used for looping over pairs in the shuffle vector
  ! Used in loop to calculate the statistical weights 
  ! and in the actual collision loop
  integer :: particle_pair_loop_iterator
  
  ! Index of last pair to collide (needed in case num_par is odd)
  integer :: last_pair
  
  ! Translation from shuffle vector array to particle array
  integer :: id_a, id_b
  
  ! Self-explanatory
  real(p_k_part) :: devnull
  
  ! Executable statements
  
  ! Return if no pairs are in cell
  if( num_par .lt. 2 ) return
  
  ! If there are an odd number of particles, collide the last with the first
  ! ( shuffle_vector( num_par + 1 ) = shuffle_vector(1) by construction )
  ! c.f. Nanbu 1998 Fig. 3(b)
  last_pair = num_par + mod( num_par, 2 )
  
  ! Aliases
  m_real  = spec%m_real
  rq_real = spec%rq_real
  omega_pe_0 = self%omega_pe_0
  ! Classical reduced mass (only for log(lambda))
  mu_ab = m_real / 2.0_p_k_part
  
  ! Calculate statistical values in collision cell
  
  u_sq_avg = 0
  n_alpha_alpha = 0
  if( self%coulomb_logarithm_automatic ) then
    
    do particle_pair_loop_iterator = 1, last_pair, 2
      id_a = shuffle_vec(particle_pair_loop_iterator)
      id_b = shuffle_vec(particle_pair_loop_iterator + 1)
      
      ! Accumulate statistical weight
      w_m = min( abs(spec%q(id_a)), abs(spec%q(id_b)) ) 
      n_alpha_alpha = n_alpha_alpha + w_m  
      
      ! Get u_a_lab, u_b_lab, g_a_lab, g_b_lab
      u_a_lab = spec%p(:,id_a)
      u_b_lab = spec%p(:,id_b)
      g_b_lab = sqrt( 1.0_p_k_part + dot_product( u_b_lab, u_b_lab ) )
      
      ! Transform u_b from lab frame to ref frame of a
      call lorentz_transform_proper( u_a_lab, u_b_lab, g_b_lab, u_rel, devnull)
      u_sq_avg = u_sq_avg + w_m * dot_product(u_rel,u_rel)
    enddo 
    u_sq_avg = u_sq_avg / n_alpha_alpha
    
  else ! Nothing to do but find the weighting
    
    do particle_pair_loop_iterator = 1, last_pair, 2
      id_a = shuffle_vec(particle_pair_loop_iterator)
      id_b = shuffle_vec(particle_pair_loop_iterator + 1)
      n_alpha_alpha = n_alpha_alpha + min( abs(spec%q(id_a)), abs(spec%q(id_b)) )
    enddo 
    
  endif
  
  if( self%coulomb_logarithm_automatic ) then
    
    ! For the lower limit, use either the distance of classical closest approach, or the first Born approximation
    ! of the quantum mechanical result for close scattering (see Landau and Lifshitz 'Quantum Mechanics'
    ! chapter 127) depending on which is greater (and of course the Debye length
    ! for the upper limit.)
    ! (Note that this line is equivalent to a simple test of u < fine structure constant, but, well, that's not how it's written )
    ! (Note also, this should technically be u*v in the classical limit, but gamma < 1.000027 so it doesn't matter)
    lambda = min( ( mu_ab*cgs_me * l_debye*(cgs_c/omega_pe_0) * u_sq_avg * cgs_c**2 ) / &
                  ( 2.0_p_k_part * spec%q_real**2 * cgs_e**2 ), &
                  ( mu_ab*cgs_me * l_debye*(cgs_c/omega_pe_0) * sqrt(u_sq_avg) * cgs_c ) / &
                  ( cgs_hbar ) )
    
    ! c.f. Spitzer 'Physics of Fully Ionized Gases' equation 5-12
    ! $\int_0^p dx \frac{x^3}{\sqrt{1+x^2}} = \frac{1}{2}\left(\frac{1}{1+p^2}-1+ln(1+p^2)\right)$
    ln_lambda = 1./2 * ( 1./(1 + lambda**2) - 1 + log(1 + lambda**2) )
    
    ! Clamp log(lambda) at 2 (an arbitrary choice)
    ln_lambda = max( ln_lambda, 2._p_k_part )
    
  else
    ln_lambda = self%coulomb_logarithm_value
  endif
  
  ! Adjust time-step for collision statistics for unequal particle weights
  ! c.f. Nanbu 1998, eq.14 and the one immediately preceding it
  dt_s = dt_n * ( abs(cell_fluid_state%q_tot) / ( 2.0_p_k_part * n_alpha_alpha ) )
  
  ! prepare s calculations
  s_first = ( spec%q_real**2 )**2 * &
            omega_pe_0**2 * &
            ln_lambda * &
            cell_fluid_state%number_density * &
            cgs_e**2 / (cgs_me * cgs_c**3) * &
            (dt_s/omega_pe_0)  
  
  if (self%PerezLowTempCorrection) then
    n_0 = cgs_me * omega_pe_0**2/(4._p_k_part * pi * cgs_e**2)
    s_prime_first = &
      (4._p_k_part * pi/3._p_k_part)**(1./3) * &
      (abs(cell_fluid_state%q_tot)) / (2.0_p_k_part * n_alpha_alpha) * &
      cell_fluid_state%number_density * n_0 * &
      dt_n / omega_pe_0 * &
      (2._p_k_part) / &
      (cell_fluid_state%number_density*n_0)**(2./3)
  endif
  
  ! Loop over pairs of particles and 'collide' them
  do particle_pair_loop_iterator = 1, last_pair, 2
    ! Get particle indexes
    id_a = shuffle_vec(particle_pair_loop_iterator)
    id_b = shuffle_vec(particle_pair_loop_iterator + 1)
    
    ! Read u, g in lab frame
    u_a_lab = spec%p(:,id_a)
    u_b_lab = spec%p(:,id_b)
    g_a_lab = sqrt( 1.0_p_k_part + dot_product( u_a_lab, u_a_lab ) )
    g_b_lab = sqrt( 1.0_p_k_part + dot_product( u_b_lab, u_b_lab ) )
    
    ! Relative velocity
    call lorentz_transform_proper( u_a_lab, u_b_lab, g_b_lab, u_rel, g_rel)
    au_rel = sqrt( dot_product( u_rel, u_rel ) )
    av_rel = au_rel/ g_rel
    
    ! Calculate center-of-momentum 4-velocity
    ucm =  ( u_a_lab + u_b_lab ) / &
          sqrt( 2.0_p_k_part * ( 1.0_p_k_part + ( g_a_lab*g_b_lab - dot_product(u_a_lab,u_b_lab) ) ) )
    gcm = sqrt(1.0_p_k_part + dot_product(ucm,ucm))
    
    ! Transform u_i to the com frame
    call lorentz_transform_proper( ucm, u_a_lab, g_a_lab, u_in_cm, g_in_cm )
    au_in_cm = sqrt( dot_product( u_in_cm, u_in_cm ) )
    
    ! Cycle if au_in_cm too small (au_in_cm < 4th root of a subnormal number)
    if( au_in_cm .le. mr_small ) cycle
    
    s = s_first * &
        1/(m_real * g_a_lab * m_real * g_b_lab) * &
        gcm * m_real * au_in_cm / (m_real * g_a_lab + m_real * g_b_lab) * &
        (m_real * g_in_cm * m_real * g_in_cm/(m_real * au_in_cm)**2 + 1 )**2
    
    if (self%PerezLowTempCorrection) then
      s_prime = s_prime_first * au_rel * cgs_c
      s = min(s, s_prime)
    endif
    
    call local_rng%harvest_real3( random_unit_interval )
    
    if (s < 0.1) then
      ctheta = 1 + s * log(random_unit_interval)
      !ctheta = max(ctheta, -1._p_k_part)
    elseif( s < 3 ) then
      A = 0.0056958 + &
          0.9560202*s - &
          0.508139*s**2 + &
          0.47913906*s**3 - &
          0.12788975*s**4 + &
          0.02389567*s**5
      A = 1./A
      ctheta = 1./A * log(exp(-A) + 2 * random_unit_interval * sinh(A))
    elseif( s < 6 ) then
      A = 3 * exp(-s)
      ctheta = 1./A * log(exp(-A) + 2 * random_unit_interval * sinh(A))
    else
      ctheta = 2*random_unit_interval - 1
    endif
    
    ! Although it should be mathematically impossible for this to be out of range,
    ! it appears with rounding errors it can happen sometimes.
    ctheta = min(ctheta, 1.0_p_k_part)
    ctheta = max(ctheta, -1.0_p_k_part)
    
    stheta = sqrt(1.0_p_k_part - ctheta**2)
    
    !call local_rng%harvest_real2( random_unit_interval )
    !phi = 2.0_p_k_part * random_unit_interval * real( pi, p_k_part )
    !sphi = sin( phi )
    !cphi = cos( phi )
    !
    ! c.f. http://mathworld.wolfram.com/CirclePointPicking.html
    ! and https://mcnp.lanl.gov/pdf_files/nbs_vonneumann.pdf
    ! (avoiding the root seems like overkill, but since von Neumann already did it)
    x_phi = 1.0_p_k_part
    y_phi = 1.0_p_k_part
    do while (x_phi**2 + y_phi**2 .gt. 1.0)
      call local_rng%harvest_real3( random_unit_interval )
      x_phi = 2*random_unit_interval - 1.0
      call local_rng%harvest_real3( random_unit_interval )
      y_phi = 2*random_unit_interval - 1.0
    end do
    cphi = 1/(x_phi**2 + y_phi**2)
    sphi = cphi * 2 * x_phi * y_phi
    cphi = cphi * (x_phi**2 - y_phi**2)
    
    ! Now that theta, phi are found, calculate delta_p (in CoM frame, but without rotation)
    p_in_cm = u_in_cm * m_real
    
    n_x = p_in_cm(1)
    n_y = p_in_cm(2)
    n_z = p_in_cm(3)
    p = sqrt(n_x**2 + n_y**2 + n_z**2)
    
    if (( n_x .ne. 0) .or. ( n_y .ne. 0)) then
      n_t = sqrt( n_x**2 + n_y**2 )
      rn_t = 1.0_p_k_part / n_t
      
      delta(1) = n_x * rn_t * n_z * stheta * cphi - &
                 n_y * rn_t * p * stheta * sphi + &
                 n_x * ( ctheta - 1.0_p_k_part)
      delta(2) = n_y * rn_t * n_z * stheta * cphi + &
                 n_x * rn_t * p * stheta * sphi + &
                 n_y * ( ctheta - 1.0_p_k_part )
      delta(3) = n_z * ( ctheta - 1.0_p_k_part ) - n_t * stheta * cphi
    else
      delta(1) = p * stheta * cphi
      delta(2) = p * stheta * sphi
      delta(3) = p * (ctheta - 1.0_p_k_part)
    endif
    
    ! Find size of macro-particle again, to weight statistical collision probability
    w_a = abs(spec%q(id_a))
    w_b = abs(spec%q(id_b))
    w_m = max( w_a, w_b )
    call local_rng%harvest_real2( random_unit_interval )
     
    ! Add delta to proper velocities (in CoM frame) (based on statistics)
    ucm = -ucm
    if( random_unit_interval .lt. w_b/w_m) then
      u_after_collision = u_in_cm + delta / m_real
      !call lorentz_transform_proper( -ucm, u_after_collision, g_in_cm, spec%p(:,id_a), devnull)
      call lorentz_transform_proper( ucm, u_after_collision, g_in_cm, u_a_lab, devnull)
      spec%p(:,id_a)= u_a_lab
    endif
    
    if( random_unit_interval .lt. w_a/w_m) then
      u_after_collision = -u_in_cm - delta / m_real
      !call lorentz_transform_proper( -ucm, u_after_collision, g_in_cm, spec%p(:,id_b), devnull)
      call lorentz_transform_proper( ucm, u_after_collision, g_in_cm, u_b_lab, devnull)
      spec%p(:,id_b) = u_b_lab
    endif
    
  enddo !all particles in cell
  
end subroutine like_collide_perez
!---------------------------------------------------

subroutine unlike_collide_perez( &
                           self, &
                           sp1, &
                           sp2, &
                           sp_id1, &
                           sp_id2, &
                           coll_cell_id0, &
                           num_par, &
                           shuffle_vec, &
                           cell_fluid_state, &
                           l_debye, &
                           dt_n, &
                           local_rng &
                          )
!---------------------------------------------------
!       update momenta for collision
!---------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  type( t_collisions ), intent(in) :: self
  class ( t_species ), pointer :: sp1, sp2
  integer, intent(in) :: sp_id1, sp_id2
  ! Currently the code collides the shuffled b_vector off the a_species in order
  ! If we were to collide the shuffled b_vector off the shuffled a_vector I think we wouldn't need coll_cell_id0
  ! Doesn't make much difference though, I don't think
  integer, dimension(:), intent(in) :: coll_cell_id0
  integer, dimension(:), intent(in) :: num_par
  integer, dimension(:,:), target :: shuffle_vec
  type( t_local_fluid_state ), dimension(:), intent(in) :: cell_fluid_state
  ! total debye length [units of c/w0]
  real(p_k_part), intent(in) :: l_debye
  ! collision cycle: dt * n_collide [units of 1/w0]
  real(p_k_part), intent(in) :: dt_n
  class (t_random), intent(inout) :: local_rng
  
  ! Local variables
  
  ! Aliases for input parameters (depending on which species has more particles to collide)
  class ( t_species ), pointer :: spec_a, spec_b
  integer :: coll_cell_id0_a, coll_cell_id0_b
  integer :: num_par_a, num_par_b
  integer, dimension(:), pointer :: shuffle_vec_a, shuffle_vec_b
  type(t_local_fluid_state) :: cell_fluid_state_a, cell_fluid_state_b
  
  ! Masses [units of m_e]
  real(p_k_part)   :: m_real_a, m_real_b, rq_real_a, rq_real_b
  real(p_k_part)   :: mu_ab ! Reduced mass
  
  ! Weight of particle a, particle b, maximum of those (used for statistical scattering)
  real(p_k_part) :: w_a, w_b, w_m
  ! The sum of all w_m for each pair, used to normalize weighted averages
  ! c.f. Nanbu 1998, first equation after eq. 12
  real(p_k_part) :: n_alpha_beta
  
  ! Random number between zero and one
  real(p_k_part) :: random_unit_interval
  
  ! Coulomb logarithm, and its argument
  real(p_k_part) :: lambda, ln_lambda
  
  ! Average of u_rel**2 (used in log(lambda))
  real(p_k_part) :: u_sq_avg
  
  ! Collision 'probability'
  real(p_k_part) :: s
  
  ! 'Base' collision 'probability'
  real(p_k_part) :: s_first
  
  ! Collision probability variables, for low temperature correction
  real(p_k_part) :: s_prime, s_prime_first
  ! Normal density for low temp correction (just from omega_pe_0)
  real(p_k_part) :: n_0
  
  ! Adjusted timestep (dt * n_collide * {weighting for statistics} )
  real(p_k_part) :: dt_s
  
  ! Reference plasma frequency (just an alias to self%omega_pe_0)
  real(p_k_part) :: omega_pe_0
  
  ! Proper 4-velocity of particles in LAB frame
  real(p_k_part), dimension(p_p_dim) :: u_a_lab, u_b_lab
  real(p_k_part)                     :: g_a_lab, g_b_lab
  
  ! Proper 4-velocity of either particle in the One Particle at Rest (OPR) frame
  ! Also referred to as the velocity of part_a RELative to part_b (or vice versa)
  real(p_k_part), dimension(p_p_dim) :: u_rel
  real(p_k_part)                     :: g_rel
  ! Magnitude of space-like proper-velocity, true velocity of a relative to b (or vice versa)
  real(p_k_part)                     :: au_rel, av_rel
  
  ! Proper 4-velocity of the Center-of-Momentum (CoM) frame (as seen in the LAB Frame)
  real(p_k_part), dimension(p_p_dim) :: ucm
  real(p_k_part)                     :: gcm
  
  ! Proper 4-velocities of the particles as seen in the CoM frame
  ! Note this is therefore also the velocity *OF* the CoM frame as seen in the respective OPR frame
  ! (The frame of the more massive particle is the one we boost out of below)
  real(p_k_part), dimension(p_p_dim) :: u_a_in_cm, u_b_in_cm
  real(p_k_part)                     :: g_a_in_cm, g_b_in_cm
  ! Magnitude of the proper velocities
  real(p_k_part)                     :: au_a_in_cm, au_b_in_cm
  
  ! Nanbu parameter
  real(p_k_part) :: A
  
  ! Scattering angle phi
  !real(p_k_part)   :: phi
  
  ! Trigonometric functions of argument the scattering angles
  ! ( cos(phi), sin(phi), cos(theta), sin(theta) )
  real(p_k_part) :: cphi, sphi, ctheta, stheta
  ! try skipping phi as well
  real(p_k_part) :: x_phi, y_phi
  
  ! Variables used for actual momenta update
  !
  real(p_k_part), dimension(p_p_dim) :: p_in_cm
  
  ! The components of the CoM momentum
  ! ( n_t = sqrt(n_x**2 + n_y**2) i.e. "n_transverse"; rn_t = 1/n_t)
  real(p_k_part) :: n_x, n_y, n_z, n_t, rn_t, p
  
  ! The change in momentum (p_new - p_old)
  real(p_k_part), dimension(p_p_dim) :: delta
  
  ! Updated space-like proper velocity after collision, in CoM
  real(p_k_part), dimension(p_p_dim) :: u_after_collision
  
  ! Numerical variables ( indexing )
  integer :: id_a, id_b, i
  
  ! Self-explanatory
  real(p_k_part) :: devnull
  
  ! Executable statements
  
  ! Return if no pairs are in cell
  if( (num_par(sp_id1) .lt. 1) .or. (num_par(sp_id2) .lt. 1) ) return
  
  ! Order the two species by the number of sim. particles (num_par_a >= num_par_b)
  if ( num_par(sp_id1) .ge. num_par(sp_id2) ) then
    !
    spec_a => sp1
    coll_cell_id0_a = coll_cell_id0(       sp_id1 )
    num_par_a = num_par(                   sp_id1 )
    shuffle_vec_a => shuffle_vec(    :   , sp_id1 ) ! This is never used
    cell_fluid_state_a = cell_fluid_state( sp_id1 )
    !
    spec_b => sp2
    coll_cell_id0_b = coll_cell_id0(       sp_id2 ) ! This is never used
    num_par_b = num_par(                   sp_id2 )
    shuffle_vec_b => shuffle_vec(    :   , sp_id2 )
    cell_fluid_state_b = cell_fluid_state( sp_id2 )
  else
    spec_a => sp2
    coll_cell_id0_a = coll_cell_id0(       sp_id2 )
    num_par_a = num_par(                   sp_id2 )
    shuffle_vec_a => shuffle_vec(    :   , sp_id2 ) ! This is never used
    cell_fluid_state_a = cell_fluid_state( sp_id2 )
    !
    spec_b => sp1
    coll_cell_id0_b = coll_cell_id0(       sp_id1 ) ! This is never used
    num_par_b = num_par(                   sp_id1 )
    shuffle_vec_b => shuffle_vec(    :   , sp_id1 )
    cell_fluid_state_b = cell_fluid_state( sp_id1 )
  endif
  
  m_real_a  = spec_a%m_real
  m_real_b  = spec_b%m_real
  rq_real_a = spec_a%rq_real
  rq_real_b = spec_b%rq_real
  ! Classical reduced mass
  mu_ab = ( m_real_a * m_real_b ) / ( m_real_a + m_real_b )
  
  omega_pe_0 = self%omega_pe_0
  
  u_sq_avg = 0
  n_alpha_beta = 0.0_p_k_part
  if( self%coulomb_logarithm_automatic ) then
    do i=1, num_par_a
      id_a = coll_cell_id0_a + i
      id_b = shuffle_vec_b(i) ! b has less particles, but shuffle_vec is periodic
      
      ! Collect sn_ab
      w_m = min( spec_a%q(id_a)*rq_real_a, spec_b%q(id_b)*rq_real_b )
      n_alpha_beta = n_alpha_beta + w_m
      
      !get u_a, u_b, g_a, g_b
      u_a_lab = spec_a%p(:,id_a)
      u_b_lab = spec_b%p(:,id_b)
      g_b_lab = sqrt( 1.0_p_k_part + dot_product( u_b_lab, u_b_lab ) )
      
      !get au_rel (transform a_b to ref frame of a)
      call lorentz_transform_proper( u_a_lab, u_b_lab, g_b_lab, u_rel, devnull)
      
      u_sq_avg = u_sq_avg + w_m * dot_product(u_rel,u_rel)
    enddo !all a
    u_sq_avg = u_sq_avg / n_alpha_beta
    
  else ! Nothing to do but find the weighting
    do i=1, num_par_a
      id_a = coll_cell_id0_a + i
      id_b = shuffle_vec_b(i)
      n_alpha_beta = n_alpha_beta + min( spec_a%q(id_a)*rq_real_a, spec_b%q(id_b)*rq_real_b )
    enddo
  endif
  
  if( self%coulomb_logarithm_automatic ) then
    
    ! Except for the q's this is the same as the like-like code
    lambda = min( ( mu_ab*cgs_me * l_debye*(cgs_c/omega_pe_0) * u_sq_avg * cgs_c**2 ) / &
                  ( 2.0_p_k_part * abs( spec_a%q_real * spec_b%q_real ) * cgs_e**2 ), &
                  ( mu_ab*cgs_me * l_debye*(cgs_c/omega_pe_0) * sqrt(u_sq_avg) * cgs_c ) / &
                  ( cgs_hbar ) )
    
    ln_lambda = 1./2 * ( 1./(1 + lambda**2) - 1 + log(1 + lambda**2) )
    
    ln_lambda = max( ln_lambda, 2._p_k_part)
    
  else
    ln_lambda = self%coulomb_logarithm_value
  endif
  
  ! Adjust time-step for collision statistics for unequal particle weights
  ! c.f. Nanbu 1998, second equation after eq. 12
  dt_s = dt_n * ((rq_real_a * cell_fluid_state_a%q_tot) / n_alpha_beta )
  
  s_first = ( spec_a%q_real * spec_b%q_real )**2 * &
            omega_pe_0**2 * &
            ln_lambda * &
            cell_fluid_state_b%number_density * &
            cgs_e**2 /(cgs_me * cgs_c**3) * &
            dt_s / omega_pe_0
  
  if (self%PerezLowTempCorrection) then
    n_0 = cgs_me * omega_pe_0**2/(4._p_k_part * pi * cgs_e**2)
    s_prime_first = &
      (4._p_k_part * pi/3._p_k_part)**(1./3) * &
      (rq_real_a * cell_fluid_state_a%q_tot) / n_alpha_beta * &
      cell_fluid_state_b%number_density * n_0 * &
      dt_n / omega_pe_0 * &
      (m_real_a + m_real_b) / &
      (max(m_real_a *cell_fluid_state_a%number_density**(2./3), m_real_b *cell_fluid_state_b%number_density**(2./3)) * &
      (n_0)**(2./3))
  endif
  
  ! Loop over particles 'a' and 'collide' them with particles 'b'
  do i=1, num_par_a
    ! Get particle indexes
    id_a = coll_cell_id0_a + i
    id_b = shuffle_vec_b(i) ! b has less particles, but shuffle_vec is periodic
    
    ! Read u, g in lab frame
    u_a_lab = spec_a%p(:,id_a)
    u_b_lab = spec_b%p(:,id_b)
    g_a_lab = sqrt( 1.0_p_k_part + dot_product( u_a_lab, u_a_lab ) )
    g_b_lab = sqrt( 1.0_p_k_part + dot_product( u_b_lab, u_b_lab ) )
    
    ! Relative velocity:
    call lorentz_transform_proper( u_a_lab, u_b_lab, g_b_lab, u_rel, g_rel)
    au_rel = sqrt( dot_product( u_rel, u_rel ) )
    av_rel = au_rel/ g_rel
            
    
    ! Calculate Center-of-Momentum 4-velocity
    ucm = (m_real_a * u_a_lab + m_real_b * u_b_lab) / &
          sqrt( m_real_a**2 + m_real_b**2 + 2.0_p_k_part * &
                m_real_a * m_real_b * ( g_a_lab*g_b_lab - dot_product(u_a_lab,u_b_lab) ))
    gcm = sqrt(1.0_p_k_part + dot_product(ucm,ucm))
    
    ! Transform u_i to the CoM frame
    ! m_real_a * u_a_in_cm = -m_real_b * u_b_in_cm
    call lorentz_transform_proper( ucm, u_a_lab, g_a_lab, u_a_in_cm, g_a_in_cm )
    call lorentz_transform_proper( ucm, u_b_lab, g_b_lab, u_b_in_cm, g_b_in_cm )
    
    ! Calculate magnitude of the proper velocities in the CoM frame
    au_a_in_cm = sqrt( dot_product( u_a_in_cm, u_a_in_cm ) )
    au_b_in_cm = sqrt( dot_product( u_b_in_cm, u_b_in_cm ) )
    
    ! Cycle if au_a_in_cm too small (au_a_in_cm < 4th root of a subnormal number)
    if( au_a_in_cm .le. mr_small ) cycle
    
    s = s_first * &
        1/(m_real_a * g_a_lab * m_real_b * g_b_lab) * &
        gcm * m_real_a * au_a_in_cm / (m_real_a * g_a_lab + m_real_b * g_b_lab) * &
        (m_real_a * g_a_in_cm * m_real_b * g_b_in_cm/(m_real_a * au_a_in_cm)**2 + 1 )**2
    
    if (self%PerezLowTempCorrection) then
      s_prime = s_prime_first * au_rel * cgs_c
      s = min(s, s_prime)
    endif
    
    call local_rng%harvest_real3( random_unit_interval )
    
    if (s < 0.1) then
      ctheta = 1 + s * log(random_unit_interval)
      !ctheta = max(ctheta, -1._p_k_part)
    elseif( s < 3 ) then
      A = 0.0056958 + &
          0.9560202*s - &
          0.508139*s**2 + &
          0.47913906*s**3 - &
          0.12788975*s**4 + &
          0.02389567*s**5
      A = 1./A
      ctheta = 1./A * log(exp(-A) + 2 * random_unit_interval * sinh(A))
    elseif( s < 6 ) then
      A = 3 * exp(-s)
      ctheta = 1./A * log(exp(-A) + 2 * random_unit_interval * sinh(A))
    else
      ctheta = 2*random_unit_interval - 1
    endif
    
    ctheta = min(ctheta, 1.0_p_k_part)
    ctheta = max(ctheta, -1.0_p_k_part)
    
    stheta = sqrt(1.0_p_k_part - ctheta**2)
    
    !call local_rng%harvest_real2( random_unit_interval )
    !phi = 2.0_p_k_part * random_unit_interval * real(pi,p_k_part)
    !sphi = sin( phi )
    !cphi = cos( phi )
    x_phi = 1.0_p_k_part
    y_phi = 1.0_p_k_part
    do while (x_phi**2 + y_phi**2 .gt. 1.0)
      call local_rng%harvest_real3( random_unit_interval )
      x_phi = 2*random_unit_interval - 1.0
      call local_rng%harvest_real3( random_unit_interval )
      y_phi = 2*random_unit_interval - 1.0
    end do
    cphi = 1/(x_phi**2 + y_phi**2)
    sphi = cphi * 2 * x_phi * y_phi
    cphi = cphi * (x_phi**2 - y_phi**2)
    
    ! Now that theta, phi are found, calculate delta_p (in CoM frame, but without rotation)
    p_in_cm = u_a_in_cm * m_real_a
    p = sqrt(p_in_cm(1)**2 + p_in_cm(2)**2 + p_in_cm(3)**2)
    
    n_x = p_in_cm(1)
    n_y = p_in_cm(2)
    n_z = p_in_cm(3)
    
    if (( n_x .ne. 0) .or. ( n_y .ne. 0)) then
      n_t = sqrt( n_x**2 + n_y**2 )
      rn_t = 1.0_p_k_part / n_t
      
      delta(1) = n_x * rn_t * n_z * stheta * cphi - &
                 n_y * rn_t * p * stheta * sphi + &
                 n_x * ( ctheta - 1.0_p_k_part)
      delta(2) = n_y * rn_t * n_z * stheta * cphi + &
                 n_x * rn_t * p * stheta * sphi + &
                 n_y * ( ctheta - 1.0_p_k_part )
      delta(3) = n_z * ( ctheta - 1.0_p_k_part ) - n_t * stheta * cphi
    else
      delta(1) = p * stheta * cphi
      delta(2) = p * stheta * sphi
      delta(3) = p * (ctheta - 1.0_p_k_part)
    endif
    
    ! Find size of macro-particle again, to weight statistical collision probability
    w_a = spec_a%q(id_a)*rq_real_a
    w_b = spec_b%q(id_b)*rq_real_b
    w_m = max( w_a, w_b )
    call local_rng%harvest_real2( random_unit_interval )
    
    ! Add delta to proper velocities (in CoM frame) (based on statistics)
    ucm = -ucm
    if( random_unit_interval .lt. w_b/w_m) then
      u_after_collision = u_a_in_cm + delta / m_real_a
      !call lorentz_transform_proper( -ucm, u_after_collision, g_a_in_cm, spec_a%p(:,id_a), devnull)
      call lorentz_transform_proper( ucm, u_after_collision, g_a_in_cm, u_a_lab, devnull)
      spec_a%p(:,id_a) = u_a_lab
    endif
    
    if( random_unit_interval .lt. w_a/w_m) then
      u_after_collision = u_b_in_cm - delta / m_real_b
      !call lorentz_transform_proper( -ucm, u_after_collision, g_b_in_cm, spec_b%p(:,id_b), devnull)
      call lorentz_transform_proper( ucm, u_after_collision, g_b_in_cm, u_b_lab, devnull)
      spec_b%p(:,id_b) = u_b_lab
    endif
    
  enddo ! All species_a particles in cell
  
end subroutine unlike_collide_perez
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Collision cell sorting routines
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine sort_coll_species( species, n, t, nx_collision_cells, coll_aid )
!---------------------------------------------------------------------------------------------------
! Sort the species along collision cells
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  
  class ( t_species ), intent(inout) :: species

  integer, intent(in) :: n
  real(p_double), intent(in) :: t
  integer, dimension(p_x_dim), intent(in) :: nx_collision_cells
  integer, dimension(:), intent(out) :: coll_aid    
  
  ! Local variables
  
  integer, dimension(:), pointer :: idx
  
  ! Executable statements
  
  ! sort particles at certain timesteps if required
  
  ! need to sort for collisions even if numpart==0
  if ( (species%n_sort > 0) .and. (t > species%push_start_time) .and. (species%num_par .gt. 0) ) then
    
    if ( mod( n, species%n_sort ) == 0 ) then
      
      ! Currently we need to sort even non-colliding species at every collision step,
      ! as they are still needed for the calculation of the collision cell fluid state
      
      call begin_event(sort_coll_ev)
      
      call alloc( idx, (/ species%num_par /) )
      
      call begin_event( sort_coll_genidx_ev )
      
      select case (p_x_dim)
        case (1)
          call generate_sort_coll_idx_1d(species, idx, nx_collision_cells, coll_aid )
        case (2)
          call generate_sort_coll_idx_2d(species, idx, nx_collision_cells, coll_aid )
        case (3)
          call generate_sort_coll_idx_3d(species, idx, nx_collision_cells, coll_aid )
      end select          
      
      call end_event( sort_coll_genidx_ev )
      
      call begin_event( sort_coll_rearrange_ev )
      
      call rearrange( species, idx )
     
      call end_event( sort_coll_rearrange_ev )
       
      call freemem( idx )
      
      call end_event(sort_coll_ev)
      
    endif
    
  endif
    
end subroutine sort_coll_species
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine generate_sort_coll_idx_1d( species, ip, nx_collision_cells, coll_aid )
!---------------------------------------------------------------------------------------------------
!       generate sort indexes for a 1d run
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  !
  ! intent(in)
  class ( t_species ), intent(in) :: species

  integer, dimension(:), intent(in) :: nx_collision_cells
  ! intent(out)
  integer, dimension(:), intent(out) :: ip
  integer, dimension(:), intent(out) :: coll_aid
  
  ! Local variables
  
  integer, pointer, dimension(:) :: npic
  integer :: i
  integer :: index
  
  ! The total number of collision cells on the node
  integer :: number_of_collision_cells
  
  ! The size of a collision cell, in grid cells
  integer :: collision_cell_size
  
  integer :: isum,ist
  integer :: n_grid
  
  ! Executable statements
  
  n_grid = species%my_nx_p( 3, 1 )
  
  collision_cell_size = nx_collision_cells(1)
  number_of_collision_cells = n_grid / collision_cell_size
  
  call alloc(npic, (/n_grid/) )
  npic = 0
  
  do i=1,species%num_par
    index = species%ix(1,i)
    npic(index) = npic( index ) + 1
    ip(i)=index
  enddo

  ! At this point npic( n ) is the number of particles in the nth cell
  ! And ip( m ) is the index of the cell the mth particle will be sorted into,
  ! that is it is an index referencing the grid
  
  isum=0
  do i=1,n_grid
    ist=npic(i)
    npic(i)=isum
    isum=isum+ist
  end do
  ! At this point npic( n ) is the total number of particles in cells before the nth cell
  
  do i=1,species%num_par
    index=ip(i)
    npic(index)=npic(index)+1
    ip(i)=npic(index)
  end do
  ! At this point npic( n ) is the total number of particles in cells up to and including the nth cell
  ! And ip( m ) is the index within the particle arrays where the mth particle will be put,
  ! that is it is an index referencing the species object
  
  ! Return array with the last indexes per collision cell
  ! i*collision_cell_size is the index of the last grid cell within the ith collision cell
  ! npic(i*collision_cell_size) is the index within the particle array of the last particle
  ! in the i*collision_cell_sizeth grid cell, and so is also the index of the last particle
  ! within the ith collision cell
  do i=1,number_of_collision_cells
     coll_aid(i) = npic(i*collision_cell_size)
  enddo
  
  call freemem( npic )
  
end subroutine generate_sort_coll_idx_1d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine generate_sort_coll_idx_2d( species, ip, nx_collision_cells, coll_aid )
!---------------------------------------------------------------------------------------------------
!       generate sort indexes for a 2d run
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  
  class ( t_species ), intent(in) :: species
  integer, dimension(:), intent(out) :: ip
    
  integer, dimension(:), intent(in) :: nx_collision_cells
  integer, dimension(:) :: coll_aid
  
  ! Local variables
  
  integer, pointer, dimension(:) :: npic
  integer :: i
  integer :: index
  
  integer :: index_xpic, index_ypic
  integer :: index_xcoll, index_ycoll
  integer :: index_dxpic, index_dypic
  integer :: n_coll_cells_x, n_coll_cells_y, n_coll_cells
  integer :: s_coll_cell_x, s_coll_cell_y, s_coll_cell
  
  integer :: isum,ist
  integer :: n_grid, n_grid_x, n_grid_y
  
  ! Executable statements
  
  n_grid_x = species%my_nx_p(3, 1)
  n_grid_y = species%my_nx_p(3, 2)

  n_grid = n_grid_x * n_grid_y
  
  s_coll_cell_x = nx_collision_cells(1)
  s_coll_cell_y = nx_collision_cells(2)
  s_coll_cell = s_coll_cell_x * s_coll_cell_y
  n_coll_cells_x = n_grid_x / s_coll_cell_x
  n_coll_cells_y = n_grid_y / s_coll_cell_y
  n_coll_cells = n_coll_cells_x * n_coll_cells_y
  
  call alloc( npic, (/ n_grid /) )
  npic = 0
  
  do i=1,species%num_par
    ! The simulation cell the particle is in
     ! Less one because of the way index is calculated below
     index_xpic = species%ix(1,i) - 1
     index_ypic = species%ix(2,i) - 1
     ! The collision cell the particle is in
     index_xcoll = index_xpic / s_coll_cell_x
     index_ycoll = index_ypic / s_coll_cell_y
     ! The index of the simulation cell within the collision cell (i.e. modulo)
     index_dxpic = index_xpic - index_xcoll * s_coll_cell_x
     index_dypic = index_ypic - index_ycoll * s_coll_cell_y
     ! The first line is the first index within the collision cell,
     ! and the second line is the offset of the index of the simulation cell within the collision cell
     index = s_coll_cell * ( index_xcoll + n_coll_cells_x * index_ycoll ) + &
             index_dxpic + s_coll_cell_x * index_dypic + 1
             
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
  ! ASSERT(isum == species%num_par)
  
  do i=1,species%num_par
     index=ip(i)
     npic(index) = npic(index) + 1
     ip(i) = npic(index)
  end do
  
  !return array with the last indexes per collision cell
  do i=1,n_coll_cells
     coll_aid(i) = npic(i*s_coll_cell)
  enddo
  
  call freemem( npic )
  
end subroutine generate_sort_coll_idx_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine generate_sort_coll_idx_3d( species, ip, nx_collision_cells, coll_aid )
!---------------------------------------------------------------------------------------------------
!       generate sort indexes for a 3d run
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  
  class ( t_species ), intent(in) :: species
  integer, dimension(:), intent(out) :: ip
    
  integer, dimension(:), intent(in) :: nx_collision_cells
  integer, dimension(:) :: coll_aid
  
  ! Local variables
  
  integer, pointer, dimension(:) :: npic
  integer :: i
  integer :: index
  
  integer :: index_xpic, index_ypic, index_zpic
  integer :: index_xcoll, index_ycoll, index_zcoll
  integer :: index_dxpic, index_dypic, index_dzpic
  integer :: n_coll_cells_x, n_coll_cells_y, n_coll_cells_z
  integer :: n_coll_cells
  integer :: s_coll_cell_x, s_coll_cell_y, s_coll_cell_z
  integer :: s_coll_cell
  
  integer :: isum, ist
  integer :: n_grid_x, n_grid_y, n_grid_z
  integer :: n_grid
  
  ! Executable statements
  
  n_grid_x = species%my_nx_p(3, 1)
  n_grid_y = species%my_nx_p(3, 2)
  n_grid_z = species%my_nx_p(3, 3)

  n_grid = n_grid_x * n_grid_y * n_grid_z
  
  s_coll_cell_x = nx_collision_cells(1)
  s_coll_cell_y = nx_collision_cells(2)
  s_coll_cell_z = nx_collision_cells(3)
  s_coll_cell = s_coll_cell_x * s_coll_cell_y * s_coll_cell_z
  n_coll_cells_x = n_grid_x / s_coll_cell_x
  n_coll_cells_y = n_grid_y / s_coll_cell_y
  n_coll_cells_z = n_grid_z / s_coll_cell_z
  n_coll_cells = n_coll_cells_x * n_coll_cells_y * n_coll_cells_z
  
  call alloc( npic, (/n_grid/) )
  npic = 0
  
  do i=1,species%num_par
     ! The simulation cell the particle is in
     ! Less one because of the way index is calculated below
     index_xpic = species%ix(1,i) - 1
     index_ypic = species%ix(2,i) - 1
     index_zpic = species%ix(3,i) - 1
     ! The collision cell the particle is in
     index_xcoll = index_xpic / s_coll_cell_x
     index_ycoll = index_ypic / s_coll_cell_y
     index_zcoll = index_zpic / s_coll_cell_z
     ! The index of the simulation cell within the collision cell (i.e. modulo)
     index_dxpic = index_xpic - index_xcoll * s_coll_cell_x
     index_dypic = index_ypic - index_ycoll * s_coll_cell_y
     index_dzpic = index_zpic - index_zcoll * s_coll_cell_z
     
     ! If anyone can confirm that Josh did this right that would be great
     index = s_coll_cell * ( index_xcoll + n_coll_cells_x * index_ycoll + &
                             n_coll_cells_x * n_coll_cells_y * index_zcoll) + &
             index_dxpic + s_coll_cell_x * index_dypic + s_coll_cell_x * s_coll_cell_y * index_dzpic + 1
             
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
  ! ASSERT(isum == species%num_par)
  
  do i=1,species%num_par
     index=ip(i)
     npic(index) = npic(index) + 1
     ip(i) = npic(index)
  end do
  
  !return array with the last indexes per collision cell
  do i=1,n_coll_cells
     coll_aid(i) = npic(i*s_coll_cell)
  enddo
  
  call freemem(npic)
  
end subroutine generate_sort_coll_idx_3d
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------
function find_t_local_fluid_state( this, coll_cell_id0, num_par, num_cells, calculate_lambda_debye )
!---------------------------------------------------
  
  ! Dummy variables
  class ( t_species ), intent(in) :: this
  integer, intent(in) :: coll_cell_id0
  integer, intent(in) :: num_par
  integer, intent(in) :: num_cells
  logical, intent(in) :: calculate_lambda_debye
  
  type( t_local_fluid_state) :: find_t_local_fluid_state
  
  ! Local variables
  ! Numerical
  integer :: i
  ! Always equal to i plus an offset of coll_cell_id0
  integer:: part_id
  
  ! Physics
  ! Total four-momentum in the lab frame, used to find Lorentz transform to center-of-momentum frame
  real(p_k_part), dimension(p_p_dim) :: p_tot
  real(p_k_part) :: e_tot
  ! Current particle weight, total weight of all particles in collision cell (note Lorentz invariant)
  real(p_k_part) :: w_i, w_tot
  ! Current particle momentum, energy (here we're calling total energy g instead of e, to keep you on your toes)
  real(p_k_part), dimension(p_p_dim) :: u_i
  real(p_k_part) :: g_i
  ! The fluid velocity of the plasma in the lab-frame, and the Lorentz factor associated with that velocity
  real(p_k_part), dimension(p_p_dim) :: u_fluid
  real(p_k_part) :: g_fluid
  ! Current particle momentum, energy, *in the fluid rest frame*
  real(p_k_part), dimension(p_p_dim) :: up_i
  real(p_k_part) :: gp_i
  ! "Total gamma minus 1", the weighted sum of g-1, used to calculate average energy per particle
  real(p_k_part) :: g_1_tot
  ! Average kinetic energy per particle
  real(p_k_part) :: e_kin_ave
  ! Charge density, number density, and the inverse of the square of the Debye length for this species
  real(p_k_part) :: charge_density, number_density, rsq_ldebye
  
  ! If not calculating lambda debye, just get total charge, don't need to get weights of anything
  real(p_k_part) :: total_q
  
  ! Executable statements
  
  w_tot = 0._p_k_part
  p_tot = 0._p_k_part
  e_tot = 0._p_k_part
  
  if ( num_par > 0 .and. calculate_lambda_debye ) then
    
    do i=1, num_par
      part_id = coll_cell_id0 + i
      
      ! fetch velocity of particles and store it locally
      u_i = this%p(:,part_id)
      g_i = sqrt(1.0_p_k_part + dot_product(u_i,u_i))
      ! fetch weight of particles
      w_i = abs(this%q(part_id))
      
      ! total weight
      w_tot = w_tot + w_i
      
      ! total energy and momentum four vector
      p_tot = p_tot + w_i * u_i
      e_tot = e_tot + w_i * g_i
    enddo
    
    ! Fluid velocity
    ! p = gamma*m*beta*c
    ! u = gamma*beta*c
    ! u = p/m
    ! ('m' is the rest-mass of the system of particles, not the sum of their rest masses)
    ! total four-momentum = ( e_tot/c, p_tot)
    ! four-momentum inner product : e_tot^2/c^2 - p_tot^2 = m^2 c^2
    ! u = p/m = p/sqrt(e^2 - p^2)
    ! (c = 1 sometimes)
    
    u_fluid = p_tot / sqrt( e_tot**2 - dot_product(p_tot,p_tot))
    g_fluid = sqrt( 1.0_p_k_part + dot_product(u_fluid, u_fluid) )
    
    ! Transform all the velocities into the fluid reference frame
    ! Sum up all (gamma-1) = u^2/(1+gamma^2)
    g_1_tot = 0._p_k_part
    
    do i=1, num_par
      part_id = i + coll_cell_id0
      u_i = this%p(:,part_id)
      g_i = sqrt(1.0_p_k_part + dot_product(u_i,u_i))
      
      w_i = abs(this%q(part_id))
      
      ! transform to com frame
      call lorentz_transform_proper(u_fluid, u_i, g_i, up_i, gp_i)
      
      ! sum[wp_i ( gp_i -1)]
      g_1_tot = g_1_tot + ( ( w_i * dot_product(up_i,up_i) ) / ( gp_i + 1.0_p_k_part ) )
    enddo
    
    
    ! Average kinetic energy per real particle
    if (w_tot > 0._p_k_part ) then
      e_kin_ave = g_1_tot * this%m_real / w_tot ![units of m_e c^2]
    else
      if (num_par .gt. 0) then
        ERROR("Any (in fact all) particles in collision cell have zero charge, I don't see how that's right")
        call abort_program(-1)
      else
        ERROR("Honestly I'm just confused")
        call abort_program(-1)
      endif
    endif
    
    ! Densties
    number_density = w_tot / num_cells / abs(this%q_real) ! [units of n0]
    charge_density = number_density * this%q_real ! [units of |e| n0]
    
    ! Debye length
    ! w^2 = w0^2 * number_density * Z^2 / m_real = w0^2 * charge_density * Z / m_real
    ! v_th^2 = 2/3 e_kin_ave / m_real * c^2 (assuming a Maxwellian, which I mean)
    if ( e_kin_ave .gt. 0._p_k_part ) then
      rsq_ldebye = (3.0_p_k_part * charge_density * this%q_real) / &
                   (2.0_p_k_part * e_kin_ave  )  ! [units of w0^2/c^2]
    else
      rsq_ldebye = 0._p_k_part
    endif
                 
    ! Return values
    find_t_local_fluid_state%e_kin_ave = e_kin_ave
    find_t_local_fluid_state%q_tot = sign(w_tot, this%q_real)
    find_t_local_fluid_state%u_fluid = u_fluid
    find_t_local_fluid_state%charge_density = charge_density
    find_t_local_fluid_state%number_density = number_density
    find_t_local_fluid_state%rsq_ldebye = rsq_ldebye
    find_t_local_fluid_state%p_th2 = e_kin_ave*e_kin_ave + 2*e_kin_ave
    
  else
     
    total_q = 0._p_k_part
    
    do i=1, num_par
      part_id = coll_cell_id0 + i
      ! total charge
      total_q = total_q + this%q(part_id)
    enddo
    
    charge_density = total_q / num_cells ![units of |e| n0]
    number_density = charge_density / this%q_real ![units of n0]
    
    find_t_local_fluid_state%q_tot = total_q
    find_t_local_fluid_state%charge_density = charge_density
    find_t_local_fluid_state%number_density = number_density
    
    find_t_local_fluid_state%e_kin_ave = -1.0_p_k_part
    find_t_local_fluid_state%u_fluid = -1.0_p_k_part
    find_t_local_fluid_state%rsq_ldebye = -1.0_p_k_part
    find_t_local_fluid_state%p_th2 = -1.0_p_k_part
    
  endif
  
end function find_t_local_fluid_state
!---------------------------------------------------

!-------------------------------------------------------------------------------
function if_collide( self, n )
!-------------------------------------------------------------------------------
! Check if collisions are to occur at this timestep
!-------------------------------------------------------------------------------
  
  implicit none
  
  type( t_collisions ), intent(in) :: self
  integer, intent(in) :: n
  logical :: if_collide
  
  if_collide = .false.
  if ( self%n_collide > 0 ) then
    if ( mod( n, self%n_collide ) == 0 ) if_collide = .true.
  endif
  
end function if_collide
!-------------------------------------------------------------------------------

!---------------------------------------------------
subroutine generate_shuffle_list( list, base0, num, length, local_rng)
!---------------------------------------------------
!       generate a shuffeld list with indexes
!       [base0+1,base0+num]
!---------------------------------------------------
  
  use m_random
  
  ! Dummy variables
  integer, dimension(:) :: list
  integer, intent(in) :: base0, num, length
  class (t_random), intent(inout) :: local_rng
  
  ! Local variables
  integer :: i, j, temp, f, r, p0
  
  ! Executable statements
  
  ! If no particles, nothing to do
  if( num == 0 ) return
  
  ! Populate vector with the particle indexes
  do i=1, num
    list(i)=i+base0
  enddo
  
  ! Fisher-Yates
  ! http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
  !do i=2, num
  do i = num, 2, -1
    ! Get random between 1 and i
    j = int(local_rng%genrand_real2()*i)+1
    temp = list(i)
    list(i) = list(j)
    list(j) = temp
  enddo
  
  ! Periodically extend the vector to its size
  ! This if statement seems to solve no problems, and perhaps should be removed
  if( num == 1 ) then
    !temp = list(1)
    !list = temp
    list = list(1)
    return
  endif
   
  f = length / num
  r = length - num * f
  p0 = num
  
  do j=2, f
    do i=1, num
      list(p0+i) = list(i)
    enddo
    p0 = p0 + num
  enddo
  do i=1, r
    list(p0+i) = list(i)
  enddo
  
end subroutine generate_shuffle_list
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine lorentz_transform_proper_double( u, a, a0, b, b0 )
!---------------------------------------------------------------------------------------------------
! Transforms the 4-vector (a,a0) into (b,b0) using the proper velocity u
! Will transform any 4-vector correctly, it doesn't also have to be a 4-velocity
!
! Note unlike the transform in os-util, this does not take the gamma factor of the transform
! as an input; doing so (and not sanitizing) serves no purpose and runs the risk that (g,u)
! isn't a 4-velocity. This will always give a proper Lorentz transform (I can't of course
! guarantee it's the one the user wanted though).
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  
  real(p_double), dimension(p_p_dim), intent(in) :: u
  
  real(p_double), dimension(p_p_dim), intent(in) :: a
  real(p_double), intent(in) :: a0
  real(p_double), dimension(p_p_dim), intent(out) :: b
  real(p_double), intent(out) :: b0
  
  real(p_double) :: g
  real(p_double) :: a_u ! dot_product(a,u)
  real(p_double) :: eta ! dot_product(a,u) / (1 + g) - a0
  
  g = sqrt( 1.0_p_double + u(1)**2 + u(2)**2 + u(3)**2 )
  
  a_u = a(1)*u(1) + a(2)*u(2) + a(3)*u(3)
  eta = a_u / ( g + 1.0_p_double ) - a0
  
  b(1) = a(1) + eta*u(1)
  b(2) = a(2) + eta*u(2)
  b(3) = a(3) + eta*u(3)
  
  ! b0 = g * a0 - dot_product(a,u)
  b0 = g * a0 - a_u
  
end subroutine lorentz_transform_proper_double
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine lorentz_transform_proper_single( u, a, a0, b, b0 )
!---------------------------------------------------------------------------------------------------
! Transforms the 4-vector (a,a0) into (b,b0) using the proper velocity u
! Will transform any 4-vector correctly, it doesn't also have to be a 4-velocity
!---------------------------------------------------------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  
  real(p_single), dimension(p_p_dim), intent(in) :: u
  
  real(p_single), dimension(p_p_dim), intent(in) :: a
  real(p_single), intent(in) :: a0
  real(p_single), dimension(p_p_dim), intent(out) :: b
  real(p_single), intent(out) :: b0
  
  real(p_single) :: g
  real(p_single) :: a_u ! dot_product(a,u)
  real(p_single) :: eta ! dot_product(a,u) / (1 + g) - a0
  
  g = sqrt( 1.0_p_single + u(1)**2 + u(2)**2 + u(3)**2 )
  
  a_u = a(1)*u(1) + a(2)*u(2) + a(3)*u(3)
  eta = a_u / ( g + 1.0_p_single ) - a0
  
  b(1) = a(1) + eta*u(1)
  b(2) = a(2) + eta*u(2)
  b(3) = a(3) + eta*u(3)
  
  ! b0 = g * a0 - dot_product(a,u)
  b0 = g * a0 - a_u
  
end subroutine lorentz_transform_proper_single
!---------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine cleanup_coll( this )
!-------------------------------------------------------------------------------
! Cleanup collsions object
!-------------------------------------------------------------------------------
  implicit none
  
  type(t_collisions), intent(inout) :: this
  
  !call freemem(this%rng_array)
  if( associated (this%cell_rngs) ) deallocate( this%cell_rngs )
  if( associated (this%node_rng)  ) deallocate( this%node_rng  )
  
end subroutine cleanup_coll
!-------------------------------------------------------------------------------

subroutine restart_write_collisions(this, restart_handle)
  
  implicit none
  
  type(t_collisions), intent(in) :: this
  type(t_restart_handle), intent(inout) :: restart_handle
  
  integer :: ierr, i
  
  if( this%n_collide .lt. 1) return
  
  restart_io_wr( p_collisions_rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then 
   ERROR('error writing restart data for collisions.')
   call abort_program(p_err_rstwrt)
  endif
  
  restart_io_wr( this%number_of_collision_cells, restart_handle, ierr)
  if ( ierr/=0 ) then 
   ERROR('error writing restart data for collisions.')
   call abort_program(p_err_rstwrt)
  endif
  
  !call restart_write(this%node_rng, restart_handle)
  call this%node_rng%write_checkpoint( restart_handle )
  
  do i = 1, this%number_of_collision_cells
    !call restart_write(this%rng_array(i), restart_handle)
    call this%cell_rngs(i)%write_checkpoint( restart_handle )
  enddo
  
end subroutine restart_write_collisions

subroutine restart_read_collisions(this, restart_handle)
  
  implicit none
  
  type(t_collisions), intent(inout) :: this
  type(t_restart_handle), intent(in) :: restart_handle
  
  character(len=len(p_collisions_rst_id)) :: rst_id
  character(len=*), parameter :: err_msg = 'error reading restart data for collisions object.'
  integer :: ierr, i, num_coll_cells
  
  restart_io_rd( rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then 
    ERROR('error reading restart data for collisions.')
    call abort_program(p_err_rstwrt)
  endif
  
  ! check if restart file is compatible
  if ( rst_id /= p_collisions_rst_id) then
    ERROR('Corrupted restart file, or restart file ')
    ERROR('from incompatible binary (collisions)')
    call abort_program(p_err_rstrd)
  endif
  
  restart_io_rd( num_coll_cells, restart_handle, ierr)
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )
  
  if( this%number_of_collision_cells .ne. num_coll_cells) then
    ERROR("Collision cells in restart data and input deck don't match")
    ERROR("Expected", this%number_of_collision_cells,"found",num_coll_cells)
    ERROR("This is non-recoverable")
    call abort_program(p_err_invalid)
  endif
  
  !call restart_read(this%node_rng, restart_handle)
  call this%node_rng%read_checkpoint( restart_handle )
  
  do i = 1, this%number_of_collision_cells
    !call restart_read(this%rng_array(i), restart_handle)
    call this%cell_rngs(i)%read_checkpoint( restart_handle )
  enddo

end subroutine restart_read_collisions

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Nanbu solver routines
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!---------------------------------------------------
function Nanbu_Func( a, ems )
!---------------------------------------------------
!       nonlinear function to solve
! It's not a non-linear function, it's a non-invertible function
! (Well, it is non-linear, but that's not relevant)
!---------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  real( p_k_part ), intent(in) :: a
  real( p_k_part ), intent(in) :: ems
  real( p_k_part ) :: Nanbu_Func
  
  ! Executable statements
  
  ! a/3 has higher accuracy than the actual functions for very small y
  ! (The result is 0 for a<=1.0e-8)
  if ( a < 1.0e-5 ) then
    Nanbu_Func = a/3.0_p_k_part - exp(-ems)
  else
    Nanbu_Func = 1.0_p_k_part/tanh(a) - 1.0_p_k_part/a - exp(-ems)
  endif
  
end function Nanbu_Func
!---------------------------------------------------

!---------------------------------------------------
function Nanbu_Derivative( a )
!---------------------------------------------------
!       derivate of the nonlinear function
!---------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  real( p_k_part ), intent(in) :: a
  real( p_k_part ) :: Nanbu_Derivative
  
  ! The result of the function is 4e-9 off the real answer for a = 1.0e-4
  ! 1/3 is 7e-10 off.
  if( a < 1.0e-4 ) then
    Nanbu_Derivative = 1.0_p_k_part / 3.0_p_k_part
  else
    Nanbu_Derivative = a**(-2) - sinh(a)**(-2)
  endif
  
end function Nanbu_Derivative
!---------------------------------------------------

!---------------------------------------------------
function solve_nanbu_equation( ems )
!---------------------------------------------------
!
! The idea here is to solve for y(x) where we have the equation
! coth(y) - 1/y = exp(-x)
! (c.f. Nanbu, K. (1997),
! ('Theory of cumulative small-angle collisions in plasmas',
! (Physical Review E, 55(4), 4642–4652.
! (doi:10.1103/PhysRevE.55.4642 )
!
! It attempts to accomplish this by creating the function
! f(y,x) (or f_x(y) if you will) = coth(y) -1/y - exp(-x)
! and then finding the root of that equation in y
! It uses the Newton-Raphson Method (c.f. Numerical Recipes in C section 9.4) to attempt this.
!
! I have edited it for readability, but have not confirmed that it works properly -Josh
!---------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  real( p_k_part ), intent(in) :: ems
  real( p_k_part ) :: solve_nanbu_equation
  
  ! Local variables
  ! for interval prediction
  real( p_k_part ) :: a1, a2
  real( p_k_part ) :: f1, f2
  real( p_k_part ) :: factor
  
  !for solve
  integer, PARAMETER :: MAXIT=10000
  real( p_k_part ), PARAMETER :: aacc=1.0E-32_p_k_part
  
  integer :: j
  real( p_k_part ) :: df, da, daold, f, fh, fl, temp, ah, al
  
  ! Executable statements
  j=0
  ! Guessing inital interval
  if (ems == 0) then
    a1=1E32_p_k_part
    a2=1E33_p_k_part
  elseif (ems > 1.58525_p_k_part) then
    a1 = 0.0_p_k_part
    a2 = 1.00215_p_k_part / ems
  else
    a1 = 1.00215_p_k_part / ems
    a2 = a1 + 0.1_p_k_part
  endif
  
  factor=2
  f1 = Nanbu_Func( a1, ems )
  f2 = Nanbu_Func( a2, ems )
  
  do
    if( (f1 >= 0.0_p_k_part .and. f2 <= 0.0_p_k_part) .or. &
        (f1 <= 0.0_p_k_part .and. f2 >= 0.0_p_k_part) &
      ) exit
    
    if (abs(f1) < abs(f2)) then
      temp = a1
      a1 = a1 + factor * (a1-a2)
      a2 = temp
      f2 = f1
      f1 = Nanbu_Func(a1, ems)
    else
      temp = a2
      a2 = a2 + factor * (a2-a1)
      a1 = temp
      f1 = f2
      f2 = Nanbu_Func(a2, ems)
    end if
    
  enddo
  
  fl = Nanbu_Func(a1, ems)
  fh = Nanbu_Func(a2, ems)
  
  if (fl == 0.0_p_k_part) then
    solve_nanbu_equation=a1
    RETURN
  else if (fh == 0.0_p_k_part) then
    solve_nanbu_equation = a2
    RETURN
  else if (fl < 0.0_p_k_part) then !Orient the search so that f(al) < 0.0_p_k_part
    al=a1
    ah=a2
  else
    ah=a1
    al=a2
  end if
  
  solve_nanbu_equation = a1 !0.5_p_k_part * (a1+a2) !Initialize the guess for root,
  daold=abs(a2-a1)      !the “stepsize before last,”
  da=daold              !and the last step.
  
  f  = Nanbu_Func(solve_nanbu_equation, ems)
  df = Nanbu_Derivative(solve_nanbu_equation )
  
  do j=1,MAXIT !Loop over allowed iterations.
    
    if (((solve_nanbu_equation-ah)*df-f)*((al-solve_nanbu_equation)*df-f) >= 0.0_p_k_part .or. &
    abs(2.0_p_k_part*f) > abs(daold*df) ) then
    
      !Bisect if Newton out of range, or not decreasing fast enough.
      daold=da
      da=0.5_p_k_part*(ah-al)
      temp=solve_nanbu_equation
      solve_nanbu_equation=al+da
      if (temp == solve_nanbu_equation) then !Change in root is negligible.
        RETURN
      endif
    else !Newton step acceptable. Take it.
      daold=da
      da=f/df
      temp=solve_nanbu_equation
      solve_nanbu_equation=solve_nanbu_equation-da
      if (temp == solve_nanbu_equation) then
        RETURN
      endif
    end if
    
    !One new function evaluation per iteration.
    f  = Nanbu_Func(solve_nanbu_equation, ems)
    df = Nanbu_Derivative(solve_nanbu_equation )
    if( f == 0.0_p_k_part ) then
      RETURN
    endif
    
    if (f < 0.0_p_k_part) then !Maintain the bracket on the root.
      al=solve_nanbu_equation
    else
      ah=solve_nanbu_equation
    end if
  
  end do
  
end function solve_nanbu_equation
!---------------------------------------------------

!---------------------------------------------------
subroutine A_lookup_setup_linear()
!---------------------------------------------------
!       sets up then lookup table
!---------------------------------------------------
  implicit none
  
  !       executable statements
  
  ! allocate table
  call alloc( lookup_table_A, (/ p_lookup_size /) )
  
  ! write values
  !     #include "os-lookup.f90"
  !          include "os-lookup.f90"
     continue
  
  ERROR("Linear lookup table not written.")
  ERROR("If you want to use that feature you'll have to write one.")
  call abort_program(p_err_notimplemented)
  
end subroutine A_lookup_setup_linear
!---------------------------------------------------

!---------------------------------------------------
function A_lin_inter( s )
!---------------------------------------------------
!       nonlinear function lookup
!       lookup_table_A is suposed to contain values of A
!       starting with s = p_lookup_min ending with s = p_lookup_max + p_lookup_delta
!       and with a spacing of delta_s = p_lookup_delta = 1/p_r_lookup_delta
!---------------------------------------------------
  implicit none
  
  ! Dummy variables
  real( p_k_part ), intent(in)            :: s
  real( p_k_part )                        :: A_lin_inter
   
  ! Local variables
  real( p_k_part )                        :: normalized
  real( p_k_part )                        :: w
  
  integer                          :: index
  
  ! Executable statements
  
  ERROR("Linear lookup table not written.")
  ERROR("If you want to use that feature you'll have to write one.")
  call abort_program(p_err_notimplemented)
  
  normalized = ( s - p_lookup_min ) * p_r_lookup_delta
  index = int( normalized )
  w = normalized - real( index, p_k_part )
  
  A_lin_inter = (1.0_p_k_part - w) * lookup_table_A(index+1) + w * lookup_table_A(index+2)
  
end function A_lin_inter
!---------------------------------------------------

!---------------------------------------------------
subroutine A_lookup_setup_cubic( lookup_file_A, lookup_file_Ap)
!---------------------------------------------------
!       sets up then lookup tables
!---------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  character(len=*), intent(in)       :: lookup_file_A, lookup_file_Ap
  
  ! Executable statements
  
  ! Allocate table
  call alloc( lookup_table_A, (/p_lookup_size/) )
  call alloc( lookup_table_Ap, (/p_lookup_size/) )
  
  call lookup_read_file(lookup_file_A, lookup_table_A)          
  call lookup_read_file(lookup_file_Ap, lookup_table_Ap)
  
end subroutine A_lookup_setup_cubic
!---------------------------------------------------

!---------------------------------------------------
function A_cubic_inter( s )
!---------------------------------------------------
!       nonlinear function lookup
!       lookup_table_A is suposed to contain values of A
!       starting with s = p_lookup_min ending with s = p_lookup_max + p_lookup_delta
!       and with a spacing of delta_s = p_lookup_delta = 1/p_r_lookup_delta
!
!       lookup_table_Ap is suposed to contain values of dA/ds * p_lookup_delta
!       starting with s = p_lookup_min ending with s = p_lookup_max + p_lookup_delta
!       and with a spacing of delta_s = p_lookup_delta = 1/p_r_lookup_delta
!---------------------------------------------------
  implicit none
  
  ! Dummy variables
  real(p_k_part), intent(in)            :: s
  real(p_k_part)                        :: A_cubic_inter
  
  ! Local variables
  real(p_k_part)                        :: normalized
  real(p_k_part)                        :: w
  integer                        :: index
  
  real(p_k_part)                        :: f0, f1, fd0, fd1
  real(p_k_part)                        :: a, b, c, d
  
  ! Executable statements
  normalized = ( s - p_lookup_min ) * p_r_lookup_delta
  index = int( normalized )
  w = normalized - real( index, p_k_part)
  
  f0  = lookup_table_A(index+1)
  f1  = lookup_table_A(index+2)
  fd0 = lookup_table_Ap(index+1)
  fd1 = lookup_table_Ap(index+2)
  
  a = f0
  b = fd0
  c = -3.0_p_k_part * f0 + 3.0_p_k_part * f1 - 2.0_p_k_part * fd0 - fd1
  d = 2.0_p_k_part * f0 - 2.0_p_k_part * f1 + fd0 + fd1
  
  A_cubic_inter = ((d * w + c) * w + b) * w + a ! == a + b*w + c*w^2 + d*w^3
  
end function A_cubic_inter
!---------------------------------------------------

!---------------------------------------------------
subroutine lookup_read_file( fname, tab )
!---------------------------------------------------
!       sets up then lookup tables
!---------------------------------------------------
  
  implicit none
  
  ! Dummy variables
  character(len=*), intent(in)       :: fname
  real(p_k_part), dimension(:), intent(inout)    :: tab
  
  ! Local variables
  integer                        :: ierr, i
  
  ! Executable statements
  SCR_ROOT( "Reading lookuptable: ", fname )
  open(file_id_tem, file=fname, iostat = ierr, action = "read")
  if( ierr /= 0 ) then
    ERROR("Problem opening lookup table: ", fname)
    ERROR("Abort!")
    call abort_program()
  endif
  
  do i=1, size(tab)
    read(file_id_tem, '(E21.5)', iostat=ierr) tab(i)
    if( ierr /= 0 ) then
      ERROR("Problem reading value (", i, ") from lookup table: ", fname)
      ERROR("Abort!")
      call abort_program()
    endif
  enddo
  
  close(file_id_tem, iostat = ierr)
  if( ierr /= 0 ) then
    ERROR("Problem closing lookup table: ", fname)
    ERROR("Abort!")
    call abort_program()
  endif
  
end subroutine lookup_read_file
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines

#define __TYPE__ type( t_local_fluid_state )
#define __TYPE_STR__ "t_local_fluid_state"
#define FNAME( a )  a ## _local_fluid_state
#include "memory/mem-template.h"


end module m_species_collisions

#endif
      ! def __HAS_COLLISIONS__
