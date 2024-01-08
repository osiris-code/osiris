!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     parameter module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"

module m_parameters

use m_system

implicit none

public

!==============================================================
! simulation compile time parameters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! no. of spatial dimensions. Currently set to the preprocessor macro P_X_DIM
integer, parameter :: p_x_dim  = P_X_DIM

! Data precision
! p_k_fld   : Field quantities
! p_k_part  : Particle quantities

#if defined(PRECISION_SINGLE)

! Single precision code
integer, parameter :: p_k_part = p_single
integer, parameter :: p_k_fld  = p_single

#elif defined(PRECISION_DOUBLE)

! Double precision code
integer, parameter :: p_k_part = p_double
integer, parameter :: p_k_fld  = p_double

#else

! Choose custom numerical precision for each of the quantities

! precision of particle quantities
!integer, parameter :: p_k_part = p_single
integer, parameter :: p_k_part = p_double

! precision of field quantities
!integer, parameter :: p_k_fld = p_single
integer, parameter :: p_k_fld = p_double

#endif

!==============================================================
! parameters describing algorithms
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

integer, parameter :: p_standard          = 0  ! standard EM-PIC algorithm
integer, parameter :: p_sim_pgc           = 3  ! PGC algorithm
integer, parameter :: p_sim_qed           = 4  ! QED algorithm
integer, parameter :: p_sim_cyl_modes     = 5  ! quasi-3D (cylindrical modes) algorithm
integer, parameter :: p_sim_shear         = 6  ! SHEAR algorithm
integer, parameter :: p_sim_rad           = 7  ! radiat
integer, parameter :: p_sim_tiles         = 8  ! Tiling algorithm
integer, parameter :: p_sim_overdense     = 9  ! Overdense algorithm (damping and splitting)
integer, parameter :: p_sim_overdense_cyl = 10 ! Overdense-cyl algorithm
integer, parameter :: p_sim_xxfel         = 11 ! More percise external fields for particle motion (e.g. xxfel)
integer, parameter :: p_sim_qed_cyl       = 12 ! quasi-3D (cylindrical modes) + QED  algorithm
integer, parameter :: p_sim_neut_spin     = 13 ! TDSE ionization (spin-dependent)
integer, parameter :: p_sim_gr            = 14 ! GR algorithm
integer, parameter :: p_sim_rad_cyl       = 15 ! radiat cyl

integer, parameter :: p_std        = 0   ! standard pusher
integer, parameter :: p_simd       = 1   ! Hardware optimized pusher

integer, parameter :: p_vay        = 2
integer, parameter :: p_fullrot    = 3
integer, parameter :: p_euler      = 4

integer, parameter :: p_beam_accel = 5   ! Beam acceleration pusher

integer, parameter :: p_pgc        = 6   ! Ponderomotive guiding center pusher
integer, parameter :: p_radcool    = 7   ! Radiation cooling pusher
integer, parameter :: p_shear      = 8   ! shear frame pusher
integer, parameter :: p_cond_vay   = 9   ! Conditional (gamma>5) standard/Vay pusher
integer, parameter :: p_cary       = 10  ! Cary pusher
integer, parameter :: p_exact      = 11  ! Exact pusher
integer, parameter :: p_exact_rr   = 12  ! Exact pusher with radiation reaction


!==============================================================
! parameters describing coordinates/dimensions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
integer, parameter :: p_p_dim   = 3 !no. of dim. of momentum
integer, parameter :: p_f_dim   = 3 !no. of comp. of a vectorfield
integer, parameter :: p_s_dim   = 3 !no. of dim. of spin

! axial, radial and azimuthal directions
integer, parameter :: p_z_dim = 1
integer, parameter :: p_r_dim = 2
integer, parameter :: p_t_dim = 3

!==============================================================
! parameters to identify coordinate systems
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
integer, parameter :: p_cartesian     = 0
integer, parameter :: p_cylindrical_b = 1 ! B1 on axis

!==============================================================
! parameters to identify type of current deposition / field interpolation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
integer, parameter :: p_ngp           = 0 ! nearest grid point (not implemented)
integer, parameter :: p_linear        = 1 ! linear interpolation (i.e. area weighting)
integer, parameter :: p_quadratic     = 2 ! quadratic splines
integer, parameter :: p_cubic         = 3 ! cubic splines
integer, parameter :: p_quartic       = 4 ! quartic (4th order) splines

!==============================================================
! parameters to identify boundary conditions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Invalid boundary conditions (used to ensure the user specifies bc)
integer, parameter :: p_bc_invalid         = -2
! boundary to other node
integer, parameter :: p_bc_other_node      = -1
! Periodic b.c.
integer, parameter :: p_bc_periodic        =  0
! boundary moving with c
integer, parameter :: p_bc_move_c          =  1
! Particle absorption / open boundary
integer, parameter :: p_bc_absorbing       =  5

! Conducting (perfect electric conductor)
integer, parameter :: p_bc_pec             =  6
! Reflecting (perfect magnetic conductor)
integer, parameter :: p_bc_pmc             =  7

! Internally PEC and PMC boundaries are different for odd or even interpolation levels
integer, parameter :: p_bc_pec_odd         =  6
integer, parameter :: p_bc_pmc_odd         =  7
integer, parameter :: p_bc_pec_even        =  8
integer, parameter :: p_bc_pmc_even        =  9

! Axial boundary for r-z simulations
integer, parameter :: p_bc_cyl_axis        = 20

! Lindman absorbing boundary
integer, parameter :: p_bc_lindman         = 30

! Particle specular reflection
integer, parameter :: p_bc_specular        = 40

! Particle thermal bath
integer, parameter :: p_bc_thermal         = 50

! VPML absorbing boundaries
integer, parameter :: p_bc_vpml            = 60

integer, parameter :: p_bc_none            = 100

!==============================================================
! parameters to identify pair production models
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
integer, parameter :: p_pairprod_qed  = 1
integer, parameter :: p_pairprod_gthr = 2

!==============================================================
!       parameters to define boundary operations at internal boundaries
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  vdf boundary type/operations with regard to periodic or internode b.c.
!  the following possibilities are currently given:
!  p_vdf_replace := 1 replacement of guard cell values by values
!                     from other node or opposite boundary
!                     (example: E and B fields)
!  p_vdf_add     := 2 summing up of guard cell values and values
!                     from corresponding cells on the other node
!                     or opposing boundary
!                     (example: current or charge densities)

integer, parameter :: p_vdf_replace = 1 ! replacement
integer, parameter :: p_vdf_add     = 2 ! summing up

! boundary indexes
integer, parameter :: p_lower = 1
integer, parameter :: p_upper = 2
integer, parameter :: p_size  = 3

!=================================================================================================
! id-numbers for outputfiles
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! temp file is file_id_tem        =  10, defined in os-system

integer, parameter, public :: file_id_rst  =  40

! load balance
integer, parameter, public :: file_id_part_load = 50

! energy reports
integer, parameter, public :: file_id_fldene  =  70
integer, parameter, public :: file_id_parene  =  71
integer, parameter, public :: file_id_partemp =  74

integer, parameter, public :: file_id_prof =  73

!=================================================================================================
!       parameters to define error codes
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


integer, parameter :: p_err_alloc           = -10 ! memory allocation failed
integer, parameter :: p_err_dealloc         = -11 ! memory deallocation failed
integer, parameter :: p_err_range           = -12 ! index out of range
integer, parameter :: p_err_invalid         = -13 ! invalid parameters
integer, parameter :: p_err_msg             = -14 ! message related errors
integer, parameter :: p_err_rstwrt          = -15 ! error writing restart file
integer, parameter :: p_err_rstrd           = -16 ! error reading restart file
integer, parameter :: p_err_nullp           = -17 ! null pointer or var. not properly initialized
integer, parameter :: p_err_diagfile        = -18 ! error writing diagnostic file
integer, parameter :: p_err_notimplemented  = -19 ! feature not implemented
integer, parameter :: p_err_nan             = -20 ! nan or infinity found
integer, parameter :: p_err_mpi             = -21 ! error in mpi call


! select number of particles to be held in sofware cache
! this number should be optimized to be compatiable with the
! hardware cache
integer, parameter :: p_cache_size = 1024

! environment variable holding working directory
character(len=*), parameter :: p_env_osiris_wdir    = "OSIRIS_WDIR"
character(len=*), parameter :: p_env_osiris_test    = "OSIRIS_TEST"
character(len=*), parameter :: p_env_osiris_restart = "OSIRIS_RESTART"

! default input file name
character(len = *), parameter :: p_default_input_file =  'os-stdin'

! Default paths
character(len = *), parameter :: p_default_mass       =  'MS'
character(len = *), parameter :: p_default_rest       =  'RE'
character(len = *), parameter :: p_default_hist       =  'HIST'
character(len = *), parameter :: p_default_time       =  'TIMINGS'

! Actual paths
character(len = p_max_filename_len) :: path_mass = trim( p_default_mass // p_dir_sep)
character(len = p_max_filename_len) :: path_rest = trim( p_default_rest // p_dir_sep)
character(len = p_max_filename_len) :: path_hist = trim( p_default_hist // p_dir_sep)
character(len = p_max_filename_len) :: path_time = trim( p_default_time // p_dir_sep)


interface list_algorithm
  module procedure list_algorithm_options
end interface

contains

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine list_algorithm_options( options )

  implicit none

  type( t_options ), intent(in) :: options

  if ( options%wall_clock_check > 0 ) then
    print *, ' '
    print *, 'Simulation will run for a maximum wall clock time of ', options%wall_clock_limit, ' s,'
    print *, 'checking at every ', options%wall_clock_check, ' iterations.'
  endif

  print *, ' '
  print *, 'Numerical precision :'

  select case ( p_k_fld )
    case( p_single )
      print *, '- Fields    : single precision'
    case( p_double )
      print *, '- Fields    : double precision'
  end select

  select case ( p_k_part)
    case( p_single )
       print *, '- Particles : single precision'
    case( p_double )
       print *, '- Particles : double precision'
  end select

end subroutine list_algorithm_options
!---------------------------------------------------------------------------------------------------



end module m_parameters
