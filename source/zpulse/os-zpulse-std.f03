
#include "os-config.h"
#include "os-preprocess.fpp"

module m_zpulse_std

#include "memory/memory.h"

use m_node_conf
use m_restart
use m_parameters 
use m_math

use m_vdf_define
use m_vdf_math
use m_vdf

use m_emf_define, only: t_emf_bound, t_emf

use m_space
use m_grid_define
use m_fparser
use m_utilities

implicit none

private

! parameters for zpulses

integer, parameter :: p_zpulse_bnormal = 0
integer, parameter :: p_zpulse_bint = 1


integer, parameter :: p_forward    = 1  ! forward propagation
integer, parameter :: p_backward   = 2  ! backward propagation

integer, parameter :: p_plane          = 0              ! plane wave
integer, parameter :: p_hermite_gaussian  = 1           ! hermite-gaussian beam
integer, parameter :: p_hermite_gaussian_astigmatic = 2 ! astigmatic beam
integer, parameter :: p_bessel            = 3              ! bessel
integer, parameter :: p_gaussian_asym     = 4              ! asymetric gaussian
integer, parameter :: p_laguerre_gaussian = 5              ! laguerre-gaussian beam

integer, parameter :: p_gaussian       = 1  ! gaussian
integer, parameter :: p_polynomial     = 2  ! polynomial
integer, parameter :: p_func           = 3  ! math function
integer, parameter :: p_sin2           = 4  ! sin^2
integer, parameter :: p_const          = 5 


integer, parameter :: p_max_chirp_order = 4

! string to id restart data
character(len=*), parameter :: p_zpulse_rst_id = "zpulse rst data - 0x0002"

type :: t_zpulse

  ! Defines how to create B Field values ( see launch_zpulse_1d_bint for details )
  integer :: b_type = p_zpulse_bnormal

  ! EM Wave Parameters
  real(p_double) :: a0
  real(p_double) :: omega0
  real(p_double) :: phase0
  integer         :: pol_type
  real(p_double) :: pol
  integer         :: propagation
  integer         :: direction

  ! k chirp parameters
  integer         :: chirp_order
  real(p_double), dimension(p_max_chirp_order) :: chirp_coefs

  ! Longitudinal Profile

  integer        :: lon_type
  real(p_double) :: lon_start

  real(p_double) :: lon_rise, lon_flat, lon_fall
  real(p_double) :: lon_duration, lon_range, lon_x0

  type(t_fparser) :: lon_math_func

  real(p_double), dimension(2) :: lon_tilt
  logical :: if_lon_tilt


  ! Perpendicular Profile
  integer         :: per_type
  real(p_double), dimension(2) :: per_center

  ! perpendicular chirp parameters
  integer         :: per_chirp_order
  real(p_double), dimension(p_max_chirp_order) :: per_chirp_coefs

  ! Hermite / Laguerre gaussian parameters
  real(p_double), dimension(2) :: per_w0, per_focus
  integer, dimension(2) :: per_tem_mode

  ! assymetric gaussian parameters
  real(p_double), dimension(2,2) :: per_w0_asym, per_focus_asym  ! gaussian parameters
  real(p_double), dimension(2) :: per_asym_trans

  ! bessel parameters
  integer         :: per_n        ! Bessel mode
  real(p_double)  :: per_kt       ! Transverse wavenumber
  integer         :: per_0clip    ! number of Zero to start clipping
  real(p_double)  :: per_clip_pos ! actual clipping position
  ! real(p_double) :: per_kbessel  ! k*sqrt(1-kt^2/k^2)

  ! Extra parameters
  logical         :: if_launch
  real(p_double)  :: launch_time
  logical         :: no_div_corr
  integer         :: interpolation

  ! boost variables
  real(p_double)  :: gamma
  real(p_double)  :: gamma_times_beta
  real(p_double)  :: boost_per_fac
  logical         :: if_boost

  ! Pointer for next item on the list
  class( t_zpulse ), pointer :: next


contains

  procedure :: init => init_zpulse
  procedure :: cleanup => cleanup_zpulse
  procedure :: read_input => read_input_zpulse
  procedure :: check_dimensionality => check_dimensionality_zpulse

  procedure :: launch => launch_zpulse

  procedure :: lon_envelope    => lon_envelope_zpulse
  procedure :: lon_center      => lon_center_zpulse
  procedure :: per_envelope_2d => per_envelope_2d_zpulse
  procedure :: per_envelope_3d => per_envelope_3d_zpulse
  procedure :: apply_lorentz_transformation => apply_lorentz_transformation_zpulse


end type t_zpulse

public :: t_zpulse

public :: p_zpulse_bnormal
public :: p_max_chirp_order, p_polynomial, p_sin2, p_gaussian, p_func, p_const
public :: p_forward, p_backward
public :: p_hermite_gaussian, p_hermite_gaussian_astigmatic, p_laguerre_gaussian, p_plane, p_gaussian_asym

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine read_input_zpulse( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

  use m_input_file

  implicit none

  class( t_zpulse ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  type( t_space ), intent(in) :: g_space
  class (t_emf_bound), intent(in) :: bnd_con
  logical, dimension(:), intent(in) :: periodic
  class( t_grid ), intent(in) :: grid
  type( t_options ), intent(in) :: sim_options

  ! local variables
  character(len=16) :: b_type
  real(p_double) :: a0
  real(p_double) :: omega0
  real(p_double) :: phase
  integer         :: pol_type
  real(p_double) :: pol
  character(len=16) :: propagation
  integer         :: direction

  integer         :: chirp_order
  real(p_double), dimension(p_max_chirp_order) :: chirp_coefs

  character(len=16) :: lon_type
  real(p_double) :: lon_start

  real(p_double) :: lon_fwhm

  real(p_double) :: lon_rise, lon_flat, lon_fall
  real(p_double) :: lon_duration, lon_x0, lon_range
  character(len = p_max_expr_len) :: lon_math_func

  real(p_double), dimension(2) :: lon_tilt

  character(len=16) :: per_type
  real(p_double), dimension(2) :: per_center

  integer            :: per_chirp_order
  real(p_double), dimension(p_max_chirp_order) :: per_chirp_coefs

  real(p_double), dimension(2) :: per_w0, per_fwhm, per_focus
  integer, dimension(2) :: per_tem_mode

  real(p_double), dimension(2,2) :: per_w0_asym, per_fwhm_asym, per_focus_asym
  real(p_double), dimension(2) :: per_asym_trans

  integer         :: per_n, per_0clip
  real(p_double) :: per_kt

  real(p_double) :: launch_time
  logical        :: if_launch, no_div_corr

  namelist /nl_zpulse/ &
    if_launch, b_type, a0, omega0,phase, pol_type, pol, propagation, direction, &
    chirp_order, chirp_coefs, lon_fwhm, lon_type, lon_start, lon_rise, lon_flat, lon_fall, &
    lon_duration, lon_x0, lon_range, lon_math_func, lon_tilt, per_type, per_center, per_w0, per_fwhm, per_focus, &
    per_chirp_order, per_chirp_coefs, per_w0_asym, per_fwhm_asym, per_focus_asym, per_asym_trans, &
    per_n, per_0clip, per_kt, launch_time, no_div_corr, per_tem_mode


  integer :: i, j, ierr

  if (disp_out(input_file)) then
    SCR_ROOT(" - Reading zpulse configuration...")
  endif

  if_launch = .true.
  b_type = "normal"

  a0             = 1.0_p_double
  omega0         = 10.0_p_double
  phase          = 0.0_p_double
  pol_type       = 0
  pol            = 90.0_p_double

  chirp_order    = 0
  chirp_coefs    = 0.0_p_double

  propagation    = "forward"
  direction      = 1

  lon_type       = "polynomial"
  lon_start      = -huge(1.0_p_double)

  lon_fwhm       = -1.0_p_double

  lon_rise       = 0.0_p_double
  lon_flat       = 0.0_p_double
  lon_fall       = 0.0_p_double

  lon_duration   = 0.0_p_double
  lon_x0         = -huge(1.0_p_double)
  lon_range      = -huge(1.0_p_double)
  lon_math_func  = "NO_FUNCTION_SUPPLIED!"

  lon_tilt       = 0.0_p_double

  per_type       = "plane"
  per_center     = -huge(1.0_p_double)

  per_chirp_order    = 0
  per_chirp_coefs    = 0.0_p_double

  per_w0         = 0.0_p_double
  per_fwhm       = 0.0_p_double
  per_focus      = -huge(1.0_p_double)
  per_tem_mode = 0 ! tem mode, integer, default = 0,0

  per_w0_asym     = 0.0_p_double
  per_fwhm_asym   = 0.0_p_double
  per_focus_asym  = -huge(1.0_p_double)
  per_asym_trans  = 0.0_p_double

  per_n          = 0
  per_0clip      = 1
  per_kt         = 0.0_p_double


  if_launch       = .true.
  launch_time           = 0.0_p_double
  no_div_corr    = .false.


  ! Get namelist text from input file
  call get_namelist( input_file, "nl_zpulse", ierr )

  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      if ( ierr < 0 ) then
        print *, "Error reading zpulse parameters"
      else
        print *, "Error: zpulse parameters missing"
      endif
      print *, "aborting..."
    endif
    stop
  endif

  read (input_file%nml_text, nml = nl_zpulse, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error reading zpulse parameters"
      write(0,*)  "   aborting..."
    endif
    stop
  endif


  select case ( trim( b_type ) )
  case ( "normal" )
    this%b_type = p_zpulse_bnormal
  case ( "int" )
    this%b_type = p_zpulse_bint

  case default
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse parameters"
      write(0,*)  "   btype must be either 'normal' or 'int'"
      write(0,*)  "   aborting..."
    endif
    stop
  end select

  this%a0         = a0
  this%omega0     = omega0
  this%phase0     = real( phase * pi_180, p_double )

  if ((pol_type <-1) .or. (pol_type >1)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse parameters"
      write(0,*)  "   pol must be in the range [-1,0,+1]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%pol_type  = pol_type
  if ( pol /= 0.0_p_double ) then
    this%pol = real( pol * pi_180, p_double )
  else
    this%pol = 0.0_p_double
  endif

  if ((chirp_order < 0) .or. (chirp_order > p_max_chirp_order)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse parameters"
      write(0,*)  "   chirp_order must be in the range [0,",p_max_chirp_order,"]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%chirp_order    = chirp_order
  this%chirp_coefs    = chirp_coefs

  if ((per_chirp_order < 0) .or. (per_chirp_order > p_max_chirp_order)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse parameters"
      write(0,*)  "   per_chirp_order must be in the range [0,",p_max_chirp_order,"]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%per_chirp_order    = per_chirp_order
  this%per_chirp_coefs    = per_chirp_coefs


  select case ( trim(propagation) )
  case ("forward")
    this%propagation = p_forward
  case ("backward")
    this%propagation = p_backward
  case default
    if ( mpi_node() == 0 ) then
      write(0,*)  ''
      write(0,*)  '   Error in zpulse parameters'
      write(0,*)  '   propagation must be either "forward" or "backward"'
      write(0,*)  '   aborting...'
    endif
    stop
  end select

  if ((direction < 1) .or. (direction > p_x_dim)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse parameters"
      write(0,*)  "   direction must be in the range [1 .. x_dim]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif

  this%direction = direction

  this%if_lon_tilt = .false.
  if ( p_x_dim > 1 ) then
    do i = 1, p_x_dim-1
      if ( lon_tilt(i) /= 0.0_p_double ) then
        if (( lon_tilt(i) < -45.0_p_double ) .or. (lon_tilt(i) > 45.0_p_double )) then
          if ( mpi_node() == 0 ) then
            write(0,*)  ""
            write(0,*)  "   Error in zpulse parameters"
            write(0,*)  "   lon_tilt must be in the range [ -45 , 45 ] deg"
            write(0,*)  "   aborting..."
          endif
          stop
        endif
        this%lon_tilt(i) = tan( real( lon_tilt(i) * pi_180, p_double ))
        if ( this%propagation == p_forward ) this%lon_tilt(i) = -this%lon_tilt(i)
        this%if_lon_tilt = .true.
      else
        this%lon_tilt(i) = 0.0_p_double
      endif
    enddo
  endif


  ! longitudinal profile

  select case ( trim(lon_type) )
  case ("polynomial")
    this%lon_type = p_polynomial

    if ( lon_fwhm > 0.0_p_double ) then
      this%lon_rise       = lon_fwhm/2.0_p_double
      this%lon_flat       = 0.0_p_double
      this%lon_fall       = lon_fwhm/2.0_p_double
    else
      this%lon_rise       = lon_rise
      this%lon_flat       = lon_flat
      this%lon_fall       = lon_fall
    endif

  case ("sin2")
    this%lon_type = p_sin2

    if ( lon_fwhm > 0.0_p_double ) then
      this%lon_rise       = lon_fwhm/2.0_p_double
      this%lon_flat       = 0.0_p_double
      this%lon_fall       = lon_fwhm/2.0_p_double
    else
      this%lon_rise       = lon_rise
      this%lon_flat       = lon_flat
      this%lon_fall       = lon_fall
    endif

  case ("gaussian")
    this%lon_type = p_gaussian

    if ( lon_fwhm > 0.0_p_double  ) then
      this%lon_duration   = lon_fwhm/sqrt( 2.0_p_double * log(2.0_p_double))
    else
      this%lon_duration   = lon_duration
    endif
    this%lon_range      = lon_range
    this%lon_x0         = lon_x0

    if ( lon_range == -huge(1.0_p_double) ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse parameters'
        write(0,*)  '   When using "gaussian" longitudinal envelopes, '
        write(0,*)  '   lon_range must also be defined'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

    if ( lon_x0 == -huge(1.0_p_double) ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse parameters'
        write(0,*)  '   When using "gaussian" longitudinal envelopes, '
        write(0,*)  '   lon_x0 must also be defined'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

  case ("math")
    this%lon_type = p_func
    this%lon_range      = lon_range
    this%lon_x0         = lon_x0

    if ( lon_range == -huge(1.0_p_double) ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse parameters'
        write(0,*)  '   When using "math" longitudinal envelopes, '
        write(0,*)  '   lon_range must also be defined'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

    if ( lon_x0 == -huge(1.0_p_double) ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse parameters'
        write(0,*)  '   When using "math" longitudinal envelopes, '
        write(0,*)  '   lon_x0 must also be defined'
        write(0,*)  '   aborting...'
      endif
      stop
    endif


    call setup(this%lon_math_func, trim(lon_math_func), (/'x'/), ierr)

    ! check if function compiled ok

    if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse parameters'
        write(0,*)  '   Supplied lon_math_func failed to compile'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

  case default
    if ( mpi_node() == 0 ) then
      write(0,*)  ''
      write(0,*)  '   Error in zpulse parameters'
      write(0,*)  '   lon_type must be either "polynomial", "gaussian",'
      write(0,*)  '   "sin2", "math" or "hermite"'
      write(0,*)  '   aborting...'
    endif
    stop
  end select

  if ( this%lon_type == p_polynomial .or. this%lon_type == p_sin2 ) then
    if ( lon_start == -huge(1.0_p_double) ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse parameters'
        write(0,*)  '   When using "polynomial" or "sin2" longitudinal envelopes, '
        write(0,*)  '   lon_start must also be defined'
        write(0,*)  '   aborting...'
      endif
      stop
    else
      this%lon_start      = lon_start
    endif
  endif


  ! perpendiular profile

  select case ( trim(per_type) )
  case ("plane")
    this%per_type       = p_plane

  case ("gaussian", "hermite")

    this%per_type       = p_hermite_gaussian
    if ( p_x_dim > 1 ) then
      this%per_tem_mode = 0
      do i = 1, p_x_dim - 1
        this%per_tem_mode(i) = per_tem_mode(i)
        if ( per_tem_mode(i) > 7 ) then
          if ( mpi_node() == 0 ) then
            write(0,*)  ''
            write(0,*)  '   Error in zpulse parameters'
            write(0,*)  '   Only TEM modes up to 7 have been implemented. Please contact the'
            write(0,*)  '   development team if you need higher order modes.'
            write(0,*)  '   aborting...'
          endif
          stop
        endif
      enddo
    endif

    ! get main spot size and focal plane
    this%per_w0 = 0.0_p_double
    if (per_fwhm(1) > 0.0_p_double) then
      this%per_w0 = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
    else
      this%per_w0 = per_w0(1)
    endif

    if ( this%per_w0(1) <= 0.0_p_double ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse parameters'
        write(0,*)  '   for gaussian pulses per_w0 or per_fwhm must be > 0'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

    if ( per_focus(1) == -huge(1.0_p_double) ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse parameters'
        write(0,*)  '   for gaussian pulses per_focus must be defined'
        write(0,*)  '   aborting...'
      endif
      stop
    endif
    this%per_focus      = per_focus(1)

    ! check for astigmatic beam
    if (p_x_dim == 3) then
      if (per_fwhm(2) > 0.0_p_double) then
        this%per_w0(2) = per_fwhm(2)/sqrt(4.0_p_double * log(2.0_p_double))
      else if ( per_w0(2) > 0.0_p_double ) then
        this%per_w0(2) = per_w0(2)
      endif

      if ( per_focus(2) /= -huge(1.0_p_double) ) then
        this%per_focus(2)      = per_focus(2)
      else
        this%per_focus(2)      = this%per_focus(1)
      endif

      if ((this%per_w0(2) /= this%per_w0(1)) .or. (this%per_focus(2) /= this%per_focus(1))) then
        this%per_type = p_hermite_gaussian_astigmatic
      endif

    endif

  case ("laguerre")

    this%per_type       = p_laguerre_gaussian
    if ( p_x_dim > 1 ) then
      this%per_tem_mode = 0
      do i = 1, 2
        this%per_tem_mode(i) = per_tem_mode(i)
        if ( per_tem_mode(i) > 7 ) then
          if ( mpi_node() == 0 ) then
            print *, ''
            print *, '   Error in zpulse parameters'
            print *, '   Only Laguerre modes up to 7 have been implemented. Please contact the'
            print *, '   development team if you need higher order modes.'
            print *, '   aborting...'
          endif
          stop
        endif
      enddo
    endif

    ! get main spot size and focal plane
    this%per_w0 = 0.0_p_double
    if (per_fwhm(1) > 0.0_p_double) then
      this%per_w0 = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
    else
      this%per_w0 = per_w0(1)
    endif

    if ( this%per_w0(1) <= 0.0_p_double ) then
      if ( mpi_node() == 0 ) then
        print *, ''
        print *, '   Error in zpulse parameters'
        print *, '   for Laguerre-Gaussian pulses per_w0 or per_fwhm must be > 0'
        print *, '   aborting...'
      endif
      stop
    endif

    if ( per_focus(1) == -huge(1.0_p_double) ) then
      if ( mpi_node() == 0 ) then
        print *, ''
        print *, '   Error in zpulse parameters'
        print *, '   for Laguerre-Gaussian pulses per_focus must be defined'
        print *, '   aborting...'
      endif
      stop
    endif
    this%per_focus      = per_focus(1)

  case ("asymmetric")

    if ( chirp_order /= 0 .or. per_chirp_order /= 0 ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse parameters'
        write(0,*)  '   asymmetric Gaussian beams cannot be used with chirped pulses.'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

    this%per_type       = p_gaussian_asym

    ! default parameters
    if (per_fwhm(1) > 0.0_p_double) then
      this%per_w0(1) = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
    else
      this%per_w0(1) = per_w0(1)
    endif
    this%per_focus(1)      = per_focus(1)

    ! assymetric parameters
    do i = 1, 2
      do j = 1, 2
        if (per_fwhm_asym(i,j) > 0.0_p_double) then
          this%per_w0_asym(i,j) = per_fwhm_asym(i,j)/sqrt(4.0_p_double * log(2.0_p_double))
        else
          this%per_w0_asym(i,j) = per_w0_asym(i,j)
        endif

        if ( per_focus_asym(i,j) /= -huge(1.0_p_double) ) then
          this%per_focus_asym(i,j) = per_focus_asym(i,j)
        else
          this%per_focus_asym(i,j) = this%per_focus(1)
        endif
      enddo
    enddo

    ! width of transition region
    this%per_asym_trans = per_asym_trans
    do i = 1, p_x_dim-1
      if ( this%per_asym_trans(i) <= 0.0_p_double  ) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ''
          write(0,*)  '   Error in zpulse parameters'
          write(0,*)  '   Width of transition for asymmetric Gaussian beams must be > 0.'
          write(0,*)  '   aborting...'
        endif
        stop
      endif
    enddo



    if (chirp_order > 0) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse parameters'
        write(0,*)  '   Chirped Bessel beams are not allowed.'
      endif
      stop
    endif

    this%per_type       = p_bessel
    this%per_n          = per_n
    this%per_0clip      = per_0clip
    this%per_kt         = per_kt

  case default
    if ( mpi_node() == 0 ) then
      write(0,*)  ''
      write(0,*)  '   Error in zpulse parameters'
      write(0,*)  '   per_type must be one of the following:'
      write(0,*)  '   "plane", "gaussian"/"hermite", "bessel", or "asymmetric"'
      write(0,*)  '   aborting...'
    endif
    stop
  end select

  ! Set default per_center values (this needs to happen after reading the direction )
  if ( p_x_dim > 1 ) then
    if ( per_center(1) == -huge(1.0_p_double) ) then
      if ( direction == 1 ) then
        per_center(1) = 0.5_p_double * ( xmin( g_space, 2 ) + xmax( g_space, 2 ) )
      else
        per_center(1) = 0.5_p_double * ( xmin( g_space, 1 ) + xmax( g_space, 1 ) )
      endif
    endif

    if ( p_x_dim == 3 ) then
      if ( per_center(2) == -huge(1.0_p_double) ) then
        if ( direction /= 3 ) then
          per_center(2) = 0.5_p_double * ( xmin( g_space, 3 ) + xmax( g_space, 3 ) )
        else
          per_center(2) = 0.5_p_double * ( xmin( g_space, 2 ) + xmax( g_space, 2 ) )
        endif
      endif
    endif
  endif

  this%per_center     = per_center

  ! extra parameters

  this%if_launch      = if_launch
  this%launch_time    = launch_time
  this%no_div_corr    = no_div_corr

  ! Boost pulse
  this%if_boost = .false.
  
  if (sim_options%gamma > 1.0_p_double) then

    this % gamma          = sim_options%gamma
    this % gamma_times_beta = sqrt(this%gamma**2-1.0_p_double)

    if( this % propagation .eq. p_forward ) then
      this % boost_per_fac = 1.0_p_double / ( this%gamma + this%gamma_times_beta )
    else
      this % boost_per_fac = this%gamma + this%gamma_times_beta
    endif 
    SCR_ROOT('   boosting zpulse with gamma = ', this%gamma )

    this%if_boost = .true.
  else
    this % boost_per_fac = 1.0_p_double
  endif

end subroutine read_input_zpulse
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Everything related to dimensionality is checked here so that
! read_input can still be called in different dimensionality subclasses
! (e.g. cyl_modes) while this subroutine can be overridden to establish
! new dimensionality rules 
!-----------------------------------------------------------------------------------------
subroutine check_dimensionality_zpulse( this ) 

  implicit none 

  class( t_zpulse ), intent(inout) :: this

  select case( this%per_type )

  case (p_bessel)
  
    if (p_x_dim /= 3) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse parameters'
        write(0,*)  '   Bessel beams only exist in 3D geometries.'
      endif
      stop
    endif

  case( p_laguerre_gaussian )

    if (p_x_dim < 3) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse parameters'
        write(0,*)  '   Laguerre-Gaussian modes are available in 3D geometries only'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

  case( p_gaussian_asym )

    if (p_x_dim < 2 ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse parameters'
        write(0,*)  '   asymmetric Gaussian beams can only be used in 2D or 3D geometries.'
        write(0,*)  '   aborting...'
      endif
      stop
    endif


  end select

end subroutine check_dimensionality_zpulse
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine init_zpulse( this, restart, t, dt, b, g_space, nx_p_min, no_co, interpolation )

  implicit none

  class( t_zpulse ), intent(inout) :: this
  logical, intent(in) :: restart
  real( p_double ), intent(in) :: t, dt
  type( t_vdf ) ,  intent(in)  :: b
  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: interpolation

  ! When restarting, check if the pulse was previously launched and disable it if so
  if ( restart .and. this%if_launch ) then
    if ( t >= this%launch_time ) this%if_launch = .false.
  endif

  if ( this%if_launch ) then

    select case ( this%per_type )

    case (p_plane)
      ! plane wave pulses do not require divergence correction
      this%no_div_corr = .true.

    case (p_gaussian)
      ! this needs to be done locally because of chirp

      ! get rayleigh range
      !this%per_z0(1) = this%omega0 * this%per_w0(1)**2 /2

    case (p_gaussian_asym)
      ! no initialization required

    case (p_bessel)
      if (p_x_dim /= 3) then
        ERROR('Bessel beams only exist in 3D geometries.')
        call abort_program( p_err_invalid )
      endif

      if (this%per_kt/this%omega0 > 0.5_p_double ) then
        if ( mpi_node() == 0 ) then
          write(0,*) "(* warning *) kt / k > 0.5 in bessel beam, ", &
          "breaking paraxial approximation"
        endif
      endif

      ! this needs to be done locally because of chirp
      !this%per_kbessel = sqrt(this%omega0**2 + this%per_kt**2)

      if (this%per_0clip > 0) then
        this%per_clip_pos = besselJnZero(this%per_n, this%per_0clip)
      else
        this%per_clip_pos = huge(0.0_p_double)
      endif
    end select

  endif

  this % interpolation = interpolation

end subroutine init_zpulse
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine cleanup_zpulse( this )

  implicit none

  class( t_zpulse ), intent(inout) :: this

  ! no cleanup required
  continue

end subroutine cleanup_zpulse
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_1d( this, b, e, t, nx_p_min, g_space )

  implicit none

  integer, parameter :: rank = 1

  class( t_zpulse ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  real(p_double), intent(in) :: t
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space

  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  integer :: ik, i1

  real(p_double) :: zmin, z, z_2, dz_2
  real(p_double) :: lenv, lenv_2

  real(p_double) :: cos_pol, sin_pol, prop_sign
  real(p_double) :: z_center, k_z, k_z_2, amp

  ! Boost variables
  real(p_double) :: gam_one_beta
  real(p_k_fld)  :: gam_one_beta_fld

  ! executable statements
  !g_x_range = x_bnd( g_space )
  call get_x_bnd( g_space, g_x_range )
  ldx(1) = real( b%dx_(1), p_double )
  zmin = real( g_x_range(p_lower,1), p_double ) + real(nx_p_min(1)-1, p_double )*ldx(1)

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
  dz_2 = ldx(1)*0.5_p_double

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + dz_2
  end select

  if ( this%if_boost ) then
    ! equivalent to gamma * (1 - beta) = 1/[ gamma * (1 + beta) ]
    if (this%propagation==p_forward) then
      gam_one_beta = this%gamma - sqrt(this%gamma**2-1.0_p_double)
    else
      gam_one_beta = this%gamma + sqrt(this%gamma**2-1.0_p_double)
    endif
    this%omega0 = this%omega0*gam_one_beta
  else
    gam_one_beta = 0
  endif

  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  ! get pulse center for chirp and phase calculations
  z_center = this % lon_center( )

  ! loop through all grid cells and add pulse
  do i1 = b%lbound(2), b%ubound(2)
    z = zmin + real(i1-1,p_double) * ldx(1)
    z_2 = z + dz_2

    lenv   = amp* this % lon_envelope( z  , t )
    lenv_2 = amp* this % lon_envelope( z_2, t )

    ! get wavenumber
    k_z  = this%omega0
    k_z_2  = this%omega0

    ! add chirp
    do ik = 1, this%chirp_order
      k_z   = k_z   + this%chirp_coefs(ik) * (prop_sign*(z   - z_center))**ik
      k_z_2 = k_z_2 + this%chirp_coefs(ik) * (prop_sign*(z_2 - z_center))**ik
    enddo


    ! b%f1(1,i1) = b%f1(1,i1) + 0.0_p_double
    b%f1(2,i1) = b%f1(2,i1) - real( lenv_2 * cos(k_z_2*(z_2 - z_center - t)+this%phase0) * &
      sin_pol * prop_sign, p_k_fld )
    b%f1(3,i1) = b%f1(3,i1) + real( lenv_2 * cos(k_z_2*(z_2 - z_center - t)+this%phase0) * &
      cos_pol * prop_sign, p_k_fld )

    ! e%f1(1,i1) = e%f1(1,i1) + 0.0_p_double
    e%f1(2,i1) = e%f1(2,i1) + real( lenv   * cos( k_z *(z - z_center - t) + this%phase0) * &
      cos_pol, p_k_fld )
    e%f1(3,i1) = e%f1(3,i1) + real( lenv   * cos( k_z *(z - z_center - t) + this%phase0) * &
      sin_pol, p_k_fld )

  enddo

  ! Boost fields
  if (this%if_boost) then

    gam_one_beta_fld = real( gam_one_beta, p_k_fld )

    ! also boost 1 guard cell
    do i1=0, b%nx_(1)+1

      e%f1(2,i1) = gam_one_beta_fld * e%f1(2,i1)
      e%f1(3,i1) = gam_one_beta_fld * e%f1(3,i1)
      b%f1(2,i1) = gam_one_beta_fld * b%f1(2,i1)
      b%f1(3,i1) = gam_one_beta_fld * b%f1(3,i1)

    enddo

  endif


end subroutine launch_zpulse_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_2d( this, b, e, t, nx_p_min, g_space )

  implicit none

  integer, parameter :: rank = 2

  class( t_zpulse ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  real(p_double), intent(in) :: t
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space

  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  integer :: i1, i2

  real(p_double) :: zmin, z, z_2, dz_2
  real(p_double) :: rmin, r, r_2, dr_2, r_center
  real(p_double) :: lenv, lenv_2

  real(p_double) :: cos_pol, sin_pol, prop_sign

  real(p_double) :: amp

  ! Boost variables
  real(p_k_fld) :: gam_one_beta

  ! executable statements
  call get_x_bnd( g_space, g_x_range )
  ldx(1) = b%dx_(1)
  ldx(2) = b%dx_(2)

  zmin = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1)
  rmin = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2)

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0

  dz_2 = ldx(1)/2.0_p_double
  dr_2 = ldx(2)/2.0_p_double

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + dz_2
    rmin = rmin + dr_2
  end select

  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  ! loop through all grid cells and initialize
  if ( .not. this%if_lon_tilt ) then
    ! normal envelope
    do i1 = b%lbound(2), b%ubound(2)
      z = zmin + real(i1-1, p_double) * ldx(1)
      z_2 = z + dz_2

      lenv   = amp * this % lon_envelope( z  , t )
      lenv_2 = amp * this % lon_envelope( z_2, t )

      do i2 = b%lbound(3), b%ubound(3)

        r = rmin + real(i2-1, p_double) * ldx(2)
        r_2 = r + dr_2


        b%f2(1,i1,i2) = 0.0_p_k_fld
        b%f2(2,i1,i2) = real( - lenv_2 * this % per_envelope_2d( z_2,r,   t   ) * sin_pol * prop_sign, p_k_fld )
        b%f2(3,i1,i2) = real( + lenv_2 * this % per_envelope_2d( z_2,r_2, t ) * cos_pol * prop_sign, p_k_fld )

        e%f2(1,i1,i2) = 0.0_p_k_fld
        e%f2(2,i1,i2) = real( lenv   * this % per_envelope_2d( z  ,r_2, t ) * cos_pol, p_k_fld )
        e%f2(3,i1,i2) = real( lenv   * this % per_envelope_2d( z  ,r  , t ) * sin_pol, p_k_fld )

      enddo
    enddo
  else
    ! tilted envelope
    r_center = this%per_center(1)

    do i1 = b%lbound(2), b%ubound(2)
      z = zmin + real(i1-1,p_double) * ldx(1)
      z_2 = z + dz_2

      do i2 = b%lbound(3), b%ubound(3)

        r = rmin + real(i2-1,p_double) * ldx(2)
        r_2 = r + dr_2

        b%f2(1,i1,i2) = 0.0_p_k_fld
        b%f2(2,i1,i2) = real( - amp * this % lon_envelope( z_2 + this%lon_tilt(1) * ( r - r_center ), t ) * &
          this % per_envelope_2d( z_2,r, t   ) * sin_pol * prop_sign, p_k_fld )
        b%f2(3,i1,i2) = real( + amp * this % lon_envelope( z_2 + this%lon_tilt(1) * (r_2 - r_center), t ) * &
          this % per_envelope_2d( z_2,r_2, t ) * cos_pol * prop_sign, p_k_fld )

        e%f2(1,i1,i2) = 0.0_p_k_fld
        e%f2(2,i1,i2) = real( + amp * this % lon_envelope( z + this%lon_tilt(1) * ( r_2 - r_center ), t ) * &
          this % per_envelope_2d( z  ,r_2, t ) * cos_pol, p_k_fld )
        e%f2(3,i1,i2) = real( + amp * this % lon_envelope( z + this%lon_tilt(1) * ( r - r_center ), t ) * &
          this % per_envelope_2d( z  ,r  , t ) * sin_pol, p_k_fld )

      enddo
    enddo
  endif



  ! Boost fields
  if (this%if_boost) then

    SCR_ROOT('Boosting zpulse initialization...')

    ! equivalent to gamma * (1 - beta)
    if (this%propagation == p_forward) then
      gam_one_beta = real( this%gamma - sqrt(this%gamma**2-1.0_p_double), p_k_fld )
    else
      gam_one_beta = real( this%gamma + sqrt(this%gamma**2-1.0_p_double), p_k_fld )
    endif
    ! also boost 1 guard cell
    do i1=0, b%nx_(1)+1
      do i2=0, b%nx_(2)+1

        e%f2(2,i1,i2) = gam_one_beta * e%f2(2,i1,i2)
        e%f2(3,i1,i2) = gam_one_beta * e%f2(3,i1,i2)
        b%f2(2,i1,i2) = gam_one_beta * b%f2(2,i1,i2)
        b%f2(3,i1,i2) = gam_one_beta * b%f2(3,i1,i2)

      enddo
    enddo

  endif

end subroutine launch_zpulse_2d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_3d( this, b, e, t, nx_p_min, g_space )

  implicit none

  integer, parameter :: rank = 3

  class( t_zpulse ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  real(p_double), intent(in) :: t
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space

  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  integer :: i1, i2, i3

  real(p_double) :: zmin, z, z_2, dz_2
  real(p_double) :: r1min, r1, r1_2, dr1_2, r1_center
  real(p_double) :: r2min, r2, r2_2, dr2_2, r2_center
  real(p_double) :: lenv, lenv_2

  real(p_double) :: cos_pol, sin_pol, prop_sign

  real(p_double) :: amp

  ! Boost variables
  real(p_k_fld) :: gam_one_beta

  ! executable statements
  call get_x_bnd( g_space, g_x_range )
  ldx(1) = b%dx_(1)
  ldx(2) = b%dx_(2)
  ldx(3) = b%dx_(3)

  zmin  = g_x_range(1,1) + real(nx_p_min(1)-1,p_double)*ldx(1)
  r1min = g_x_range(1,2) + real(nx_p_min(2)-1,p_double)*ldx(2)
  r2min = g_x_range(1,3) + real(nx_p_min(3)-1,p_double)*ldx(3)

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0

  dz_2  = ldx(1)/2.0_p_double
  dr1_2 = ldx(2)/2.0_p_double
  dr2_2 = ldx(3)/2.0_p_double

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + dz_2
    r1min = r1min + dr1_2
    r2min = r2min + dr2_2
  end select

  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  if ( .not. this%if_lon_tilt ) then
    ! loop through all grid cells and initialize
    do i1 = b%lbound(2), b%ubound(2)
      z = zmin + real(i1-1, p_double) * ldx(1)
      z_2 = z + dz_2

      lenv   = amp * this % lon_envelope( z  , t )
      lenv_2 = amp * this % lon_envelope( z_2, t )

      do i2 = b%lbound(3), b%ubound(3)

        r1 = r1min + real(i2-1, p_double) * ldx(2)
        r1_2 = r1 + dr1_2

        do i3 = b%lbound(4), b%ubound(4)

          r2 = r2min + real(i3-1, p_double) * ldx(3)
          r2_2 = r2 + dr2_2

          !b(1,i1,i2,i3) = 0.0_p_double
          b%f3(2,i1,i2,i3) = b%f3(2,i1,i2,i3) - real( lenv_2 * this % per_envelope_3d( z_2,r1  ,r2_2, t ) * &
            sin_pol * prop_sign, p_k_fld )
          b%f3(3,i1,i2,i3) = b%f3(3,i1,i2,i3) + real( lenv_2 * this % per_envelope_3d( z_2,r1_2,r2, t   ) * &
            cos_pol * prop_sign, p_k_fld )

          !e%f3(1,i1,i2,i3) = 0.0_p_double
          e%f3(2,i1,i2,i3) = e%f3(2,i1,i2,i3) + real( lenv   * this % per_envelope_3d( z  ,r1_2, r2  , t   ) * cos_pol, p_k_fld )
          e%f3(3,i1,i2,i3) = e%f3(3,i1,i2,i3) + real( lenv   * this % per_envelope_3d( z  ,r1  , r2_2, t ) * sin_pol, p_k_fld )

        enddo
      enddo
    enddo
  else
    ! tilted envelope
    r1_center = this%per_center(1)
    r2_center = this%per_center(2)

    ! loop through all grid cells and initialize
    do i1 = b%lbound(2), b%ubound(2)
      z = zmin + real(i1-1, p_double) * ldx(1)
      z_2 = z + dz_2

      do i2 = b%lbound(3), b%ubound(3)

        r1 = r1min + real(i2-1, p_double) * ldx(2)
        r1_2 = r1 + dr1_2

        do i3 = b%lbound(4), b%ubound(4)

          r2 = r2min + real(i3-1, p_double) * ldx(3)
          r2_2 = r2 + dr2_2

          !b(1,i1,i2,i3) = 0.0_p_double
          b%f3(2,i1,i2,i3) = real( - amp * this % lon_envelope( z_2 + &
            this%lon_tilt(1) * ( r1 - r1_center ) + &
            this%lon_tilt(2) * ( r2_2 - r2_center ), t ) * &
            this % per_envelope_3d( z_2,r1  ,r2_2, t ) * sin_pol * prop_sign, p_k_fld )
          b%f3(3,i1,i2,i3) = real( +  amp * this % lon_envelope( z_2 + &
            this%lon_tilt(1) * ( r1_2 - r1_center ) + &
            this%lon_tilt(2) * ( r2 - r2_center ), t ) * &
            this % per_envelope_3d( z_2,r1_2,r2, t   ) * cos_pol * prop_sign, p_k_fld )

          !e%f3(1,i1,i2,i3) = 0.0_p_double
          e%f3(2,i1,i2,i3) = real( + amp * this % lon_envelope( z + &
            this%lon_tilt(1) * ( r1_2 - r1_center ) + &
            this%lon_tilt(2) * ( r2 - r2_center ), t ) * &
            this % per_envelope_3d( z  ,r1_2, r2, t   ) * cos_pol, p_k_fld )
          e%f3(3,i1,i2,i3) = real( + amp * this % lon_envelope( z + &
            this%lon_tilt(1) * ( r1 - r1_center ) + &
            this%lon_tilt(2) * ( r2_2 - r2_center ), t ) * &
            this % per_envelope_3d( z  ,r1  , r2_2, t ) * sin_pol, p_k_fld )

        enddo
      enddo
    enddo

  endif

  ! Boost fields
  if (this%if_boost) then

    if (this%propagation == p_forward) then
      ! equivalent to gamma * (1 - beta)
      gam_one_beta = real( this%gamma - sqrt(this%gamma**2-1.0_p_double), p_k_fld )
    else
      ! equivalent to gamma * (1 + beta)
      gam_one_beta = real( this%gamma + sqrt(this%gamma**2-1.0_p_double), p_k_fld )

    endif

    ! also boost 1 guard cell
    do i1=0, b%nx_(1)+1
      do i2=0, b%nx_(2)+1
        do i3=0, b%nx_(3)+1

          e%f3(2,i1,i2,i3) = gam_one_beta * e%f3(2,i1,i2,i3)
          e%f3(3,i1,i2,i3) = gam_one_beta * e%f3(3,i1,i2,i3)
          b%f3(2,i1,i2,i3) = gam_one_beta * b%f3(2,i1,i2,i3)
          b%f3(3,i1,i2,i3) = gam_one_beta * b%f3(3,i1,i2,i3)

        enddo
      enddo
    enddo

  endif

end subroutine launch_zpulse_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_1d_bint( this, b, e, t, nx_p_min, g_space )

  implicit none

  integer, parameter :: rank = 1

  class( t_zpulse ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  real(p_double), intent(in) :: t
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space

  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  integer :: ik, i1

  real(p_double) :: zmin, z
  real(p_double) :: lenv

  real(p_double) :: cos_pol, sin_pol, prop_sign
  real(p_double) :: z_center, k_z, amp

  real(p_k_fld) :: prop_sign_fld

  ! executable statements
  !g_x_range = x_bnd( g_space )
  call get_x_bnd( g_space, g_x_range )
  ldx(1) = real( b%dx_(1), p_double )

  zmin = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1)

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1) * 0.5_p_double
  end select

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0

  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  ! get pulse center for chirp and phase calculations
  z_center = this % lon_center( )

  ! loop through all grid cells and add pulse
  do i1 = b%lbound(2), b%ubound(2)
    z = zmin + real(i1-1, p_double) * ldx(1)

    ! get wavenumber
    k_z  = this%omega0

    ! add chirp
    do ik = 1, this%chirp_order
      k_z   = k_z   + this%chirp_coefs(ik) * (prop_sign*(z   - z_center))**ik
    enddo

    lenv   = amp * this % lon_envelope( z, t  ) * cos( k_z *(z-z_center )+this%phase0)

    ! e%f1(1,i1) = e%f1(1,i1) + 0.0_p_double
    e%f1(2,i1) = real( lenv  * cos_pol, p_k_fld )
    e%f1(3,i1) = real( lenv  * sin_pol, p_k_fld )

  enddo

  ! Set B field values by interpolating E field values
  ! This is better than calculating the B values directly, because it is faster
  ! and eliminates the ghost pulse that propagates in the oppositie direction
  ! when dt is equal to courant condition, even for low resolution

  prop_sign_fld = real( prop_sign, p_k_fld )
  do i1 = b%lbound(2), b%ubound(2)-1
    ! b%f1(1,i1) =  0.0_p_double
    b%f1(2,i1) =  - 0.5_p_k_fld * (e%f1(3,i1) + e%f1(3,i1+1)) * prop_sign_fld
    b%f1(3,i1) =  + 0.5_p_k_fld * (e%f1(2,i1) + e%f1(2,i1+1)) * prop_sign_fld
  enddo

end subroutine launch_zpulse_1d_bint
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_2d_bint( this, b, e, t, nx_p_min, g_space )

  implicit none

  ! dummy variables

  integer, parameter :: rank = 2

  class( t_zpulse ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  real(p_double), intent(in) :: t
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space

  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  integer :: i1, i2

  real(p_double) :: zmin, z
  real(p_double) :: rmin, r, r_2, dr_2
  real(p_double) :: lenv

  real(p_double) :: cos_pol, sin_pol, prop_sign
  real(p_double) :: amp

  real(p_k_fld) :: prop_sign_fld

  ! executable statements
  call get_x_bnd( g_space, g_x_range )
  ldx(1) = real( b%dx_(1), p_double )
  ldx(2) = real( b%dx_(2), p_double )

  zmin = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double )*ldx(1)
  rmin = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double )*ldx(2)

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0

  dr_2 = ldx(2)/2.0_p_double

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1) * 0.5_p_double
    rmin = rmin + dr_2
  end select

  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  ! loop through all grid cells and initialize E field
  do i1 = e%lbound(2), e%ubound(2)
    z = zmin + real(i1-1, p_double) * ldx(1)

    lenv   = amp * this % lon_envelope( z, t )

    do i2 = e%lbound(3), e%ubound(3)

      r = rmin + real(i2-1, p_double) * ldx(2)
      r_2 = r + dr_2

      e%f2(1,i1,i2) = 0.0_p_k_fld
      e%f2(2,i1,i2) = real( lenv   * this % per_envelope_2d( z  ,r_2, t ) * cos_pol, p_k_fld )
      e%f2(3,i1,i2) = real( lenv   * this % per_envelope_2d( z  ,r  , t ) * sin_pol, p_k_fld )

    enddo
  enddo

  ! Initialize B field by interpolation (see note in launch_zpulse_1d)
  prop_sign_fld = real( prop_sign, p_k_fld )

  do i1 = b%lbound(2), b%ubound(2) - 1
    do i2 = b%lbound(3), b%ubound(3)
      b%f2(1,i1,i2) =  0.0_p_k_fld
      b%f2(2,i1,i2) =  - 0.5_p_k_fld * (e%f2(3,i1,i2) + e%f2(3,i1+1,i2)) * prop_sign_fld
      b%f2(3,i1,i2) =  + 0.5_p_k_fld * (e%f2(2,i1,i2) + e%f2(2,i1+1,i2)) * prop_sign_fld
    enddo
  enddo



end subroutine launch_zpulse_2d_bint
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_3d_bint( this, b, e, t, nx_p_min, g_space )

  implicit none

  integer, parameter :: rank = 3

  class( t_zpulse ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  real(p_double), intent(in) :: t
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space

  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  integer :: i1, i2, i3

  real(p_double) :: zmin, z
  real(p_double) :: r1min, r1, r1_2, dr1_2
  real(p_double) :: r2min, r2, r2_2, dr2_2
  real(p_double) :: lenv

  real(p_double) :: cos_pol, sin_pol, prop_sign
  real(p_double) :: amp

  real(p_k_fld) :: prop_sign_fld

  ! executable statements
  call get_x_bnd( g_space, g_x_range )
  ldx(1) = b%dx_(1)
  ldx(2) = b%dx_(2)
  ldx(3) = b%dx_(3)

  zmin  = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1)
  r1min = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2)
  r2min = g_x_range(p_lower,3) + real(nx_p_min(3)-1, p_double)*ldx(3)

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0

  dr1_2 = ldx(2)/2.0_p_double
  dr2_2 = ldx(3)/2.0_p_double

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1) * 0.5_p_double
    r1min = r1min + dr1_2
    r2min = r2min + dr2_2
  end select

  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  ! loop through all grid cells and initialize E field
  do i1 = b%lbound(2), b%ubound(2)
    z = zmin + real(i1-1, p_double) * ldx(1)

    lenv   = amp * this % lon_envelope( z, t )

    do i2 = b%lbound(3), b%ubound(3)

      r1 = r1min + real(i2-1, p_double) * ldx(2)
      r1_2 = r1 + dr1_2

      do i3 = b%lbound(4), b%ubound(4)

        r2 = r2min + real(i3-1, p_double) * ldx(3)
        r2_2 = r2 + dr2_2

        e%f3(1,i1,i2,i3) = 0.0_p_k_fld
        e%f3(2,i1,i2,i3) = e%f3(2,i1,i2,i3) + real( lenv * this % per_envelope_3d( z  ,r1_2, r2,   t   ) * cos_pol, p_k_fld )
        e%f3(3,i1,i2,i3) = e%f3(3,i1,i2,i3) + real( lenv * this % per_envelope_3d( z  ,r1  , r2_2, t ) * sin_pol, p_k_fld )

      enddo
    enddo
  enddo

  ! Initialize B field by interpolation (see note in launch_zpulse_1d)
  prop_sign_fld = real( prop_sign, p_k_fld )
  do i1 = b%lbound(2), b%ubound(2) - 1
    do i2 = b%lbound(3), b%ubound(3)
      do i3 = b%lbound(4), b%ubound(4)
        b%f3(1,i1,i2,i3) =  0.0_p_k_fld
        b%f3(2,i1,i2,i3) =  - 0.5_p_k_fld * (e%f3(3,i1,i2,i3) + e%f3(3,i1+1,i2,i3)) * prop_sign_fld
        b%f3(3,i1,i2,i3) =  + 0.5_p_k_fld * (e%f3(2,i1,i2,i3) + e%f3(2,i1+1,i2,i3)) * prop_sign_fld
      enddo
    enddo
  enddo


end subroutine launch_zpulse_3d_bint
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! correct the divergence of a laser pulse by adding
! longitudinal components
!-----------------------------------------------------------------------------------------
subroutine div_corr_zpulse_2d( this, b, e, t, g_x_range, nx_p_min, g_nx)

  implicit none

  integer, parameter :: rank = 2

  class( t_zpulse ), intent(inout)             :: this
  type( t_vdf ), intent(inout)                :: e, b
  real(p_double), intent(in) :: t
  real(p_double), dimension(:,:), intent(in) :: g_x_range
  integer, dimension(:), intent(in)           :: nx_p_min, g_nx

  ! local variables
  integer, dimension(2, rank) :: lgc_num

  ! Boost variables
  real(p_double) :: gam_one_beta

  real(p_double), dimension(rank) :: ldx

  real(p_double) :: e1, b1    ! e1, b1 field correction
  real(p_double) :: e2p, e2m, b2p, b2m
  real(p_double) :: x1min, dx1
  real(p_double) :: x2min, dx2

  real(p_double) :: lenv, amp, cos_pol, sin_pol, prop_sign

  integer :: i1, i2
  integer :: idx1 = 1, idx2 = 2

  ldx(1) = real( b%dx_(1), p_double )
  ldx(2) = real( b%dx_(2), p_double )

  lgc_num(:,1:rank) = b%gc_num()

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  if (this%if_boost) then
    if (this%propagation == p_forward) then
      ! equivalent to gamma * (1 - beta)
      gam_one_beta = real( this%gamma - sqrt(this%gamma**2-1.0_p_double), p_k_fld )
    else
      ! equivalent to gamma * (1 + beta)
      gam_one_beta = real( this%gamma + sqrt(this%gamma**2-1.0_p_double), p_k_fld )

    endif
  else
    gam_one_beta = 0
  endif

  dx1   = ldx(idx1)
  dx2   = ldx(idx2)

  x1min = g_x_range(1,p_lower) + real(nx_p_min(idx1)-1, p_double)*dx1
  x2min = g_x_range(1,p_upper) + real(nx_p_min(idx2)-1, p_double)*dx2

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    x1min = x1min + dx1 * 0.5_p_double
    x2min = x2min + dx2 * 0.5_p_double
  end select

  ! apply the divergenge corr. line by line
  ! only inside the vertical physical area
  ! (no need to correct guard cells in this direction)
  do i2 = 1, b%nx_(idx2)

    e1 = 0.0_p_double
    b1 = 0.0_p_double

    ! calculate correction from space on other nodes
    ! coordinates are referenced to local grid coordinates

    do i1 = g_nx(1) - nx_p_min(1) + 1, b%nx_(1) + 1, -1
      ! correct e1
      !e2p = e(2,i1+1,i2)
      !e2m = e(2,i1+1,i2-1)

      lenv = amp * this % lon_envelope( x1min + i1*dx1, t )
      if (lenv > 0.0_p_double) then
        e2p  = lenv * this % per_envelope_2d( x1min + i1*dx1, x2min + (i2-0.5)*dx2, t ) * cos_pol
        e2m  = lenv * this % per_envelope_2d( x1min + i1*dx1, x2min + (i2-1.5)*dx2, t ) * cos_pol

        e1 = e1 + dx1*(e2p - e2m)/dx2
      else
        e1 = 0.0_p_double
      endif

      ! correct b1
      !b2p = b(2,i1,i2+1)
      !b2m = b(2,i1,i2)

      lenv = amp * this % lon_envelope( x1min + (i1-0.5)*dx1, t )
      if (lenv > 0.0_p_double) then
        b2p = - lenv * this % per_envelope_2d( x1min + (i1-0.5)*dx1, x2min + (i2)*dx2  , t ) * sin_pol * prop_sign
        b2m = - lenv * this % per_envelope_2d( x1min + (i1-0.5)*dx1, x2min + (i2-1)*dx2, t ) * sin_pol * prop_sign

        b1 = b1 + dx1*(b2p - b2m)/dx2
      else
        b1 = 0.0_p_double
      endif
    enddo

    ! ! Boost fields
    ! if (this%if_boost) then
    !   e1 = e1*gam_one_beta
    !   b1 = b1*gam_one_beta
    ! endif

    e%f2(1, b%nx_(1)+1, i2) = real( e1, p_k_fld )
    b%f2(1, b%nx_(1)+1, i2) = real( b1, p_k_fld )

    ! local node correction

    do i1 = b%nx_(1), 1-lgc_num(1,1), -1

      if ( this % lon_envelope( x1min + i1*dx1, t ) > 0.0_p_double ) then
        ! correct e1
        e2p = e%f2(2,i1+1,i2)
        e2m = e%f2(2,i1+1,i2-1)

        e1 = e1 + dx1*(e2p - e2m)/dx2
        e%f2(1, i1, i2) = real( e1, p_k_fld )
      else
        e1 = 0.0_p_double
      endif

      if ( this % lon_envelope( x1min + (i1-0.5)*dx1, t ) > 0.0_p_double ) then
        ! correct b1
        b2p = b%f2(2,i1,i2+1)
        b2m = b%f2(2,i1,i2)

        b1 = b1 + dx1*(b2p - b2m)/dx2

        b%f2(1, i1, i2) = real(b1, p_k_fld )
      else
        b1 = 0.0_p_double
      endif

    enddo

  enddo



end subroutine div_corr_zpulse_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! correct the divergence of a laser pulse by adding
! longitudinal components
!-----------------------------------------------------------------------------------------
subroutine div_corr_zpulse_2d_bint( this, b, e, t, g_x_range, nx_p_min, g_nx)

  implicit none


  integer, parameter :: rank = 2

  class( t_zpulse ), intent(inout)           :: this
  type( t_vdf ), intent(inout)               :: e, b
  real(p_double), intent(in)                 :: t 
  real(p_double), dimension(:,:), intent(in) :: g_x_range
  integer, dimension(:), intent(in)          :: nx_p_min, g_nx

  ! local variables
  integer, dimension(2, rank) :: lgc_num

  real(p_double), dimension(rank) :: ldx

  real(p_double) :: e1, b1    ! e1, b1 field correction
  real(p_double) :: e2p, e2m, b2p, b2m
  real(p_double) :: x1min, dx1
  real(p_double) :: x2min, dx2

  real(p_double) :: lenv, lenvp1, amp, cos_pol, sin_pol, prop_sign

  integer :: i1, i2
  integer :: idx1 = 1, idx2 = 2

  ldx(1) = real( b%dx_(1), p_double )
  ldx(2) = real( b%dx_(2), p_double )

  lgc_num(:,1:rank) = b%gc_num()

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  dx1   = ldx(idx1)
  dx2   = ldx(idx2)

  x1min = real( g_x_range(1,1), p_double ) + (nx_p_min(idx1)-1)*dx1
  x2min = real( g_x_range(1,2), p_double ) + (nx_p_min(idx2)-1)*dx2

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    x1min = x1min + dx1 * 0.5_p_double
    x2min = x2min + dx2 * 0.5_p_double
  end select

  ! apply the divergenge corr. line by line
  ! only inside the vertical physical area
  ! (no need to correct guard cells in this direction)

  !  do i2 = 1, b%nx_(idx2)+1
  do i2 = e%lbound(3)+1, e%ubound(3)-1

    e1 = 0.0_p_double
    b1 = 0.0_p_double

    ! calculate correction from space on other nodes
    ! coordinates are referenced to local grid coordinates

    do i1 = g_nx(1) - nx_p_min(1) + 1, b%nx_(1) + 2, -1
      ! correct e1
      !e2p = e(2,i1+1,i2)
      !e2m = e(2,i1+1,i2-1)

      lenv = amp * this % lon_envelope( x1min + i1*dx1, t )
      if (lenv > 0.0_p_double) then
        e2p  = lenv * this % per_envelope_2d( x1min + i1*dx1, x2min + (i2-0.5)*dx2, t ) * cos_pol
        e2m  = lenv * this % per_envelope_2d( x1min + i1*dx1, x2min + (i2-1.5)*dx2, t ) * cos_pol

        e1 = e1 + dx1*(e2p - e2m)/dx2
      else
        e1 = 0.0_p_double
      endif

      ! correct b1
      !b2p = b(2,i1,i2+1)
      !b2m = b(2,i1,i2)

      ! faster implementation
      if ( lenv > 0.0_p_double) then
        b2p = lenv * this % per_envelope_2d( x1min + i1*dx1, x2min + i2*dx2,     t )
        b2m = lenv * this % per_envelope_2d( x1min + i1*dx1, x2min + (i2-1)*dx2, t )

        lenvp1 = amp * this % lon_envelope( x1min + (i1-1)*dx1, t )
        if ( lenvp1 > 0.0_p_double ) then
          b2p = -0.5_p_double * prop_sign * sin_pol * ( b2p + &
            lenvp1 * this % per_envelope_2d( x1min + (i1-1)*dx1,x2min + i2*dx2,     t ))
          b2m = -0.5_p_double * prop_sign * sin_pol * ( b2m + &
            lenvp1 * this % per_envelope_2d( x1min + (i1-1)*dx1,x2min + (i2-1)*dx2, t ))
        else
          b2p = -0.5_p_double * prop_sign * sin_pol * b2p
          b2m = -0.5_p_double * prop_sign * sin_pol * b2m
        endif
        b1 = b1 + dx1*(b2p - b2m)/dx2
      else
        b1 = 0.0_p_double
      endif

      ! reference implementation
      !      lenvp1 = amp*lon_envelope( this, x1min + (i1-1)*dx1 )
      !
      !      b2p = -0.5_p_double * prop_sign * sin_pol * &
      !              ( lenv  *per_envelope_2d(this, x1min+   i1  *dx1,x2min + i2*dx2 ) + &
      !                lenvp1*per_envelope_2d(this, x1min+ (i1-1)*dx1,x2min + i2*dx2 ))
      !
      !      b2m = -0.5_p_double * prop_sign * sin_pol * &
      !              ( lenv  *per_envelope_2d(this, x1min+  i1 * dx1,x2min + (i2-1)*dx2 ) + &
      !                lenvp1*per_envelope_2d(this, x1min+(i1-1)*dx1,x2min + (i2-1)*dx2 ))
      !
      !      b1 = b1 + dx1*(b2p - b2m)/dx2

    enddo

    e%f2(1, b%nx_(1)+2, i2) = real( e1, p_k_fld )
    b%f2(1, b%nx_(1)+2, i2) = real( b1, p_k_fld )

    ! local node correction

    do i1 = b%nx_(1)+1, 1-lgc_num(1,1), -1

      if (( this % lon_envelope( x1min + (i1+1)*dx1, t ) > 0.0_p_double ) .or. &
          ( this % lon_envelope( x1min +   i1 * dx1, t ) > 0.0_p_double )) then
        ! correct e1
        e2p = e%f2(2,i1+1,i2)
        e2m = e%f2(2,i1+1,i2-1)

        e1 = e1 + dx1*(e2p - e2m)/dx2
        e%f2(1, i1, i2) = real( e1, p_k_fld )
      else
        e1 = 0.0_p_double
      endif

      if (( this % lon_envelope( x1min + (i1+0.5)*dx1, t ) > 0.0_p_double ) .or. &
          ( this % lon_envelope( x1min + (i1-0.5)*dx1, t ) > 0.0_p_double )) then
        ! correct b1
        b2p = b%f2(2,i1,i2+1)
        b2m = b%f2(2,i1,i2)

        b1 = b1 + dx1*(b2p - b2m)/dx2

        b%f2(1, i1, i2) = real( b1, p_k_fld )
      else
        b1 = 0.0_p_double
      endif

    enddo

  enddo



end subroutine div_corr_zpulse_2d_bint
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! correct the divergence of a laser pulse by adding longitudinal components
! calculations are always done in double precision
!-----------------------------------------------------------------------------------------
subroutine div_corr_zpulse_3d( this, b, e, t, g_x_range, nx_p_min, g_nx )

  implicit none

  ! dummy variables

  integer, parameter :: rank = 3

  class( t_zpulse ), intent(inout)             :: this
  type( t_vdf ), intent(inout)                :: e, b
  real(p_double), intent(in)                  :: t
  real(p_double), dimension(:,:), intent(in)  :: g_x_range
  integer, dimension(:), intent(in)           :: nx_p_min, g_nx

  ! local variables
  integer, dimension(2, rank) :: lgc_num

  ! Boost variables
  real(p_double) :: gam_one_beta

  real(p_double), dimension(rank) :: ldx

  real(p_double) :: e1, b1    ! e1, b1 field correction
  real(p_double) :: e2p, e2m, b2p, b2m
  real(p_double) :: e3p, e3m, b3p, b3m
  real(p_double) :: x1min, dx1
  real(p_double) :: x2min, dx2
  real(p_double) :: x3min, dx3

  real(p_double) :: lenv, amp, cos_pol, sin_pol, prop_sign

  integer :: i1, i2, i3
  integer :: idx1 = 1, idx2 = 2, idx3 = 3

  ldx(1) = real( b%dx_(1), p_double )
  ldx(2) = real( b%dx_(2), p_double )
  ldx(3) = real( b%dx_(3), p_double )

  lgc_num(:, 1:idx3) = b%gc_num()

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0
  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  if (this%if_boost) then
    if (this%propagation == p_forward) then
      ! equivalent to gamma * (1 - beta)
      gam_one_beta = this%gamma - sqrt(this%gamma**2-1.0_p_double)
    else
      ! equivalent to gamma * (1 + beta)
      gam_one_beta = this%gamma + sqrt(this%gamma**2-1.0_p_double)
    endif
  else
    gam_one_beta = 0
  endif

  dx1   = ldx(idx1)
  dx2   = ldx(idx2)
  dx3   = ldx(idx3)

  x1min = g_x_range(1,idx1) + (nx_p_min(idx1)-1)*dx1
  x2min = g_x_range(1,idx2) + (nx_p_min(idx2)-1)*dx2
  x3min = g_x_range(1,idx3) + (nx_p_min(idx3)-1)*dx3

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    x1min = x1min + dx1 * 0.5_p_double
    x2min = x2min + dx2 * 0.5_p_double
    x3min = x3min + dx3 * 0.5_p_double
  end select

  ! apply the divergenge corr. line by line
  ! only inside the transvers physical area
  ! (no need to correct guard cells in these directions)
  do i3 = 1, b%nx_(idx3)
    do i2 = 1, b%nx_(idx2)

      e1 = 0.0_p_double
      b1 = 0.0_p_double

      ! calculate correction from space on other nodes
      ! coordinates are referenced to local grid coordinates

      do i1 = g_nx(1) - nx_p_min(1) + 1, b%nx_(1) + 1, -1

        ! correct e1
        !e2p = e(2,i1+1,i2  ,i3  )
        !e2m = e(2,i1+1,i2-1,i3  )
        !e3p = e(3,i1+1,i2  ,i3  )
        !e3m = e(3,i1+1,i2  ,i3-1)

        lenv = amp * this % lon_envelope( x1min + i1*dx1, t )
        if (lenv > 0.0_p_double) then
          e2p  = lenv * this % per_envelope_3d( x1min + i1      *dx1, &
                                                x2min + (i2-0.5)*dx2, &
                                                x3min + (i3-1  )*dx3, t  ) * cos_pol
          e2m  = lenv * this % per_envelope_3d( x1min + i1      *dx1, &
                                                x2min + (i2-1.5)*dx2, &
                                                x3min + (i3-1  )*dx3, t  ) * cos_pol

          e3p  = lenv * this % per_envelope_3d( x1min + i1      *dx1, &
                                                x2min + (i2-1  )*dx2, &
                                                x3min + (i3-0.5)*dx3, t  ) * sin_pol
          e3m  = lenv * this % per_envelope_3d( x1min + i1      *dx1, &
                                                x2min + (i2-1  )*dx2, &
                                                x3min + (i3-1.5)*dx3, t  ) * sin_pol

          e1 = e1 + dx1*((e2p - e2m)/dx2 + (e3p - e3m)/dx3)
        else
          e1 = 0.0_p_double
        endif

        ! correct b1
        !b2p = b(2,i1,i2+1,i3  )
        !b2m = b(2,i1,i2  ,i3  )
        !b3p = b(3,i1,i2  ,i3+1)
        !b3m = b(3,i1,i2  ,i3  )

        lenv = amp * this % lon_envelope( x1min + (i1-0.5)*dx1, t )
        if (lenv > 0.0_p_double) then
          b2p = - lenv * this % per_envelope_3d( x1min + (i1-0.5)*dx1, &
                                                 x2min + (i2    )*dx2, &
                                                 x3min + (i3-0.5)*dx3, t  ) * sin_pol * prop_sign
          b2m = - lenv * this % per_envelope_3d( x1min + (i1-0.5)*dx1, &
                                                 x2min + (i2-1  )*dx2, &
                                                 x3min + (i3-0.5)*dx3, t  ) * sin_pol * prop_sign

          b3p =   lenv * this % per_envelope_3d( x1min + (i1-0.5)*dx1, &
                                                 x2min + (i2-0.5)*dx2, &
                                                 x3min + (i3    )*dx3, t  ) * cos_pol * prop_sign
          b3m =   lenv * this % per_envelope_3d( x1min + (i1-0.5)*dx1, &
                                                 x2min + (i2-0.5)*dx2, &
                                                 x3min + (i3-1  )*dx3, t  ) * cos_pol * prop_sign

          b1 = b1 + dx1*((b2p - b2m)/dx2 + (b3p - b3m)/dx3)
        else
          b1 = 0.0_p_double
        endif
      enddo

      ! ! Boost fields
      ! if (this%if_boost) then
      !   e1 = e1*gam_one_beta
      !   b1 = b1*gam_one_beta
      ! endif

      e%f3(1, b%nx_(idx1)+1, i2, i3) = real( e1, p_k_fld )
      b%f3(1, b%nx_(idx1)+1, i2, i3) = real( b1, p_k_fld )

      ! local node correction

      do i1 = b%nx_(1), 1-lgc_num(1,1), -1

        if ( this % lon_envelope( x1min + i1*dx1, t ) > 0.0_p_double ) then
          ! correct e1
          e2p = e%f3(2,i1+1,i2  ,i3  )
          e2m = e%f3(2,i1+1,i2-1,i3  )
          e3p = e%f3(3,i1+1,i2  ,i3  )
          e3m = e%f3(3,i1+1,i2  ,i3-1)

          e1 = e1 + dx1*((e2p - e2m)/dx2 + (e3p - e3m)/dx3)

          e%f3(1, i1, i2, i3) = real( e1, p_k_fld )
        else
          e1 = 0.0_p_double
        endif

        ! correct b1
        if ( this % lon_envelope( x1min + (i1-0.5)*dx1, t ) > 0.0_p_double ) then
          b2p = b%f3(2,i1,i2+1,i3  )
          b2m = b%f3(2,i1,i2  ,i3  )
          b3p = b%f3(3,i1,i2  ,i3+1)
          b3m = b%f3(3,i1,i2  ,i3  )


          b1 = b1 + dx1*((b2p - b2m)/dx2 + (b3p - b3m)/dx3)

          b%f3(1, i1, i2, i3) = real( b1, p_k_fld )
        else
          b1 = 0.0_p_double
        endif

      enddo

    enddo

  enddo

end subroutine div_corr_zpulse_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine div_corr_zpulse( this, g_space, nx_p_min, g_nx, &
                                b_pulse, e_pulse, t )
!-----------------------------------------------------------------------------------------
! select the routine for the appropriate dimension
! for correcting divergence
!
!-----------------------------------------------------------------------------------------
  implicit none

  class( t_zpulse ), intent(inout) :: this
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in)  :: nx_p_min, g_nx
  type( t_vdf ), intent(inout) :: b_pulse, e_pulse
  real(p_double), intent(in)                  :: t

  real(p_double), dimension(:,:), pointer :: lg_x_range

  ! divergence correction is not required in plane waves
  if ((this%per_type /= p_plane) .and. (.not. this%no_div_corr)) then


  ! grid parameters are stored in local variables
  ! so they can be rearanged for launching pulses in
  ! different directions

  call alloc( lg_x_range, (/2,b_pulse%x_dim_/) )

  lg_x_range = real( x_bnd( g_space ), p_double )

  ! call the launch zpulse routine with the
  ! appropriate dimensions
  select case ( this%b_type )
    case ( p_zpulse_bnormal )
     select case (b_pulse%x_dim_)
       case(1)
       ! divergence correction is not necessary in 1D

       case(2)

       call div_corr_zpulse_2d(this, b_pulse, e_pulse, t, lg_x_range, nx_p_min, g_nx)

       case(3)

       call div_corr_zpulse_3d(this, b_pulse, e_pulse, t, lg_x_range, nx_p_min, g_nx)
     end select

    case ( p_zpulse_bint )

     select case (b_pulse%x_dim_)
       case(1)
       ! divergence correction is not necessary in 1D

       case(2)

       call div_corr_zpulse_2d_bint(this, b_pulse, e_pulse, t, lg_x_range, nx_p_min, g_nx)

       case(3)

       !call div_corr_zpulse_3d_bint(this, b_pulse, e_pulse, &
       !             lg_x_range, nx_p_min, g_nx)

ERROR('Divergence correction of zpulses using interpolated B-Field')
ERROR('is not yet implemented')
       call abort_program( p_err_notimplemented )
     end select

  end select

  ! clear local values variables
  call freemem( lg_x_range )
  endif

end subroutine div_corr_zpulse
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine create_zpulse_field( this, g_space, nx_p_min,  &
                                b_pulse, e_pulse, t )
!-----------------------------------------------------------------------------------------
! select the routine for the appropriate dimension
! currently only launching pulses along x1 is allowed
!
!-----------------------------------------------------------------------------------------
  implicit none

  class( t_zpulse ), intent(inout) :: this
  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min
  type( t_vdf ), intent(inout) :: b_pulse, e_pulse
  real(p_double), intent(in)                  :: t

  ! call the launch zpulse routine with the
  ! appropriate dimensions

  select case ( this%b_type )

    case ( p_zpulse_bnormal )

     select case (b_pulse%x_dim_)
     case(1)
       call launch_zpulse_1d( this, b_pulse, e_pulse, t, &
                  nx_p_min, g_space )

     case(2)
       call launch_zpulse_2d( this, b_pulse, e_pulse, t, &
                  nx_p_min, g_space )

     case(3)
       call launch_zpulse_3d( this, b_pulse, e_pulse, t, &
                  nx_p_min, g_space )

     end select

    case ( p_zpulse_bint )

     select case (b_pulse%x_dim_)
     case(1)
       call launch_zpulse_1d_bint( this, b_pulse, e_pulse, t, &
                  nx_p_min, g_space )

     case(2)
       call launch_zpulse_2d_bint( this, b_pulse, e_pulse, t, &
                  nx_p_min, g_space )

     case(3)
       call launch_zpulse_3d_bint( this, b_pulse, e_pulse, t, &
                  nx_p_min, g_space )

     end select


  end select



end subroutine create_zpulse_field
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_zpulse( this, emf, g_space, nx_p_min, g_nx, t, dt, no_co )
!-----------------------------------------------------------------------------------------
! applies pulse to e and b fields
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_zpulse ), intent(inout) :: this
  class( t_emf ) ,  intent(inout) :: emf

  type( t_space ), intent(in) :: g_space

  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t
  real(p_double), intent(in) :: dt
  class( t_node_conf ), intent(in) :: no_co

  ! local variables
  type( t_vdf ) :: e_pulse, b_pulse
  type( t_space ) :: lg_space
  integer, dimension(p_x_dim) :: lnx_p_min, lg_nx

  ! executable statements

  if ( this%if_launch .and. ( t >= this%launch_time ) ) then

   ! create new vdfs to hold the pulse field
   call b_pulse % new( emf%b, zero = .true. )
   call e_pulse % new( emf%e, zero = .true. )

   lg_space = g_space
   lnx_p_min = nx_p_min(1:p_x_dim)
   lg_nx = g_nx(1:p_x_dim)

   ! transpose data for launching pulses along x2 and x3
   if ( this%direction > 1 ) then

     call swap_value( lnx_p_min, this%direction, 1 )
     call swap_value( lg_nx, this%direction, 1 )
     call transpose_obj( lg_space, this%direction, 1 )

     ! 3D needs additional swapping, since each transpose
     ! is obtained with 2 rotations (see details in transpose_vdf routine)
     if (b_pulse%x_dim_ == 3) then
     call swap_value( lnx_p_min, 3, 2 )
     call swap_value( lg_nx, 3, 2 )
     call transpose_obj( lg_space, 3, 2 )
       endif

     ! For some reason the ifort compiler breaks if the  optional
     ! tback is not there (only in parallel!)
     call transpose_obj( b_pulse, this%direction, 1, copy = .false., tback = .false. )
     call transpose_obj( e_pulse, this%direction, 1, copy = .false., tback = .false. )

   endif

   ! create pulse field in b_pulse and e_pulse
   call create_zpulse_field( this, lg_space, lnx_p_min, &
                             b_pulse, e_pulse, t )

   ! correct divergence
   call div_corr_zpulse( this, lg_space, lnx_p_min, lg_nx, &
                             b_pulse, e_pulse, t )

   ! transpose data back if necessary
   if ( this%direction > 1 ) then
     call transpose_obj( b_pulse, this%direction, 1, tback = .true. )
     call transpose_obj( e_pulse, this%direction, 1, tback = .true. )
   endif

   ! add laser pulse to emf
   call add( emf%b, b_pulse )
   call add( emf%e, e_pulse )

   ! take care of circular polarization
   if (this%pol_type /= 0) then
     ! clear pulse fields
     call b_pulse % zero()
     call e_pulse % zero()

     ! Get phase and polarization of second pulse
     this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )) , p_double )
     this%pol    = this%pol    + real( pi_2, p_double )

     ! transpose data for launching pulses along x2 and x3
     if ( this%direction > 1 ) then
     call transpose_obj( b_pulse, this%direction, 1, copy = .false., tback = .false. )
     call transpose_obj( e_pulse, this%direction, 1, copy = .false., tback = .false. )
     endif

     ! create pulse field in b_pulse and e_pulse
     call create_zpulse_field( this, lg_space, lnx_p_min, &
                 b_pulse, e_pulse, t )

     ! correct divergence
     call div_corr_zpulse( this, lg_space, lnx_p_min, lg_nx, &
                 b_pulse, e_pulse, t )

     ! transpose data back if necessary
     if ( this%direction > 1 ) then
     call transpose_obj( b_pulse, this%direction, 1, tback = .true. )
     call transpose_obj( e_pulse, this%direction, 1, tback = .true. )
     endif

     ! add laser pulse to emf
     call add( emf%b, b_pulse )
     call add( emf%e, e_pulse )
   endif

   ! free b_pulse and e_pulse vdfs
   call b_pulse % cleanup()
   call e_pulse % cleanup()

   ! pulse has been launched, turn off if_launch
   this%if_launch = .false.

  endif


end subroutine launch_zpulse
!-----------------------------------------------------------------------------------------





!-----------------------------------------------------------------------------------------
!       returns the value of the longitudinal
!       envelope at the requested position
!-----------------------------------------------------------------------------------------
function lon_envelope_zpulse( this, z_sim, tin_sim ) result( lon_envelope )

  implicit none

  class( t_zpulse ), intent(inout) :: this ! intent must be inout because of
                                           ! the eval( f_parser ) function
  real(p_double), intent(in) :: z_sim
  real(p_double), intent(in), optional :: tin_sim

  real(p_double) :: lon_envelope

  ! lab frame coordinates
  real(p_double) :: z, t 

  real(p_double) :: t_sim

  real(p_k_fparse), dimension(1) ::z_dbl

  ! executable statements
  if( present( tin_sim ) ) then 
    t_sim = tin_sim
  else 
    t_sim = 0.0_p_double
  endif 

  ! t is ignored unless we are boosted the simulation
  if( this % if_boost ) then 
    call this % apply_lorentz_transformation( z_sim, t_sim, z, t )
    z = z - t 
  else 
    z = z_sim 
  endif 

  select case ( this%lon_type )

  case (p_polynomial)
    if ( this%propagation == p_forward ) then
      lon_envelope = fenv_poly( this%lon_start - z, &
        this%lon_rise, this%lon_flat, this%lon_fall)
    else
      lon_envelope = fenv_poly( z - this%lon_start, &
        this%lon_rise, this%lon_flat, this%lon_fall)
    endif

  case (p_sin2)
    if ( this%propagation == p_forward ) then
      lon_envelope = fenv_sin2( this%lon_start - z, &
        this%lon_rise, this%lon_flat, this%lon_fall)
    else
      lon_envelope = fenv_sin2( z - this%lon_start, &
        this%lon_rise, this%lon_flat, this%lon_fall)
    endif

  case (p_gaussian)

    if (abs(z-this%lon_x0) > this%lon_range / 2) then
      lon_envelope = 0.0
    else
      lon_envelope = exp( - 2*((z-this%lon_x0)/this%lon_duration)**2 )
    endif

  case (p_func)

    if (abs(z-this%lon_x0) > this%lon_range / 2) then
      lon_envelope = 0.0
    else
      if ( this%propagation == p_forward ) then
        z_dbl(1) = real( z - this%lon_x0, p_k_fparse )
      else
        z_dbl(1) = real( this%lon_x0 - z, p_k_fparse )
      endif

      lon_envelope = eval( this%lon_math_func, z_dbl )
    endif

  case default

    lon_envelope = 0

  end select

end function lon_envelope_zpulse
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Total duration of the laser pulse
!-----------------------------------------------------------------------------------------
function lon_duration( this )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  class( t_zpulse ), intent(in) :: this
  real(p_double) :: lon_duration


  ! executable statements

  select case ( this%lon_type )

  case (p_polynomial, p_sin2)

      lon_duration = this%lon_rise + this%lon_flat + this%lon_fall

  case (p_gaussian, p_func)

    lon_duration = this%lon_range

  end select

end function lon_duration
!----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!       returns the value of the perpendicular
!       envelope at the requested position
!-----------------------------------------------------------------------------------------
function per_envelope_2d_zpulse( this, x1_sim, x2, tin_sim ) result( per_envelope_2d )

  implicit none

  class( t_zpulse ), intent(in) :: this

  integer, parameter :: rank = 2

  ! x1 is the longitudinal coordinate, x2 the perpendicular one
  real(p_double), intent(in) :: x1_sim, x2
  real(p_double), intent(in), optional :: tin_sim

  ! lab frame coordinates 
  real(p_double) :: x1, t

  ! simulation time 
  real(p_double) :: t_sim 

  real(p_double) :: per_envelope_2d

  ! local variables
  real(p_double) :: z_center, r_center, prop_sign
  real(p_double) :: k, z, rWl2, Wl2, z0, rho, rho2
  real(p_double) :: gouy_shift, rWu, curv

  ! asymetric laser pulses variables
  real(p_double) :: w0_asym, focus_asym, shp

  integer :: i

  ! ! Boost variables
  ! real(p_double) :: uvel, beta, rg1pb

  if( present( tin_sim ) ) then 
    t_sim = tin_sim
  else 
    t_sim = 0.0_p_double
  endif 

  if( this % if_boost ) then 
    call this % apply_lorentz_transformation( x1_sim, t_sim, x1, t )
  else 
    x1 = x1_sim 
    t = t_sim 
  endif 

  ! executable statements

  ! get pulse center for chirp and phase calculations
  z_center = this%lon_center()
  r_center = this%per_center(1)

  ! get wavenumber
  k  = this%omega0

  ! add chirp
  if ( this%propagation == p_forward ) then
    do i = 1, this%chirp_order
      k = k + this%chirp_coefs(i) * ((x1 - z_center)**(i))
    enddo
    do i = 1, this%per_chirp_order
      k = k + this%per_chirp_coefs(i) * ((x2 - r_center)**(i))
    enddo
  else
    do i = 1, this%chirp_order
      k = k + this%chirp_coefs(i) * ((z_center - x1)**(i))
    enddo
    do i = 1, this%per_chirp_order
      k = k + this%per_chirp_coefs(i) * ((r_center - x2)**(i))
    enddo
  endif

  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  ! get envelope
  select case ( this%per_type )

  case (p_plane)
    per_envelope_2d = cos( k * ( x1 - z_center - prop_sign * t ) + this%phase0 )

  case (p_hermite_gaussian)

    ! get rayleigh range
    z0 = k * this%per_w0(1)**2 * 0.5_p_double

    ! calculate radius
    rho = x2 - r_center
    rho2 = rho**2

    ! calculate gaussian beam parameters
    z  = x1 - this%per_focus(1)

    ! In some situations (i.e. chirps) z0 may be 0 so z = 0 must be treated
    ! as a special case
    if ( z /= 0 ) then
      rWl2 = z0**2 / (z0**2 + z**2)
      curv = 0.5_p_double*rho2*z/(z**2 + z0**2)
      gouy_shift = ( this%per_tem_mode(1) + 1 ) * atan2( z, z0 )
    else
      rWl2 = 1
      curv = 0
      gouy_shift = 0
    endif
    rWu = sqrt(2 * rWl2) / this%per_w0(1)

    ! note that this is different in 2D and 3D
    per_envelope_2d = sqrt( sqrt(rWl2) ) * &
    hermite( this%per_tem_mode(1), rho * rWu ) * &
    exp( - rho2 * rWl2/this%per_w0(1)**2 ) * &
    cos( k * ( x1 - z_center - prop_sign * t ) + k * curv - gouy_shift + this%phase0 )

  case (p_gaussian_asym)

    ! THIS DOES NOT WORK WITH CHIRPED PULSES
    ! (there is a possibility of a division by 0 in that case)

    ! calculate radius
    rho = x2 - r_center
    rho2 = rho**2

    ! get w0 and focsus as a smooth transition between the 2 regimes
    shp = smooth_heaviside( rho/this%per_asym_trans(1) )
    w0_asym = this%per_w0_asym(1, 1) + shp*( this%per_w0_asym(2, 1) - this%per_w0_asym(1, 1))
    focus_asym = this%per_focus_asym(1,1) + &
    shp*(this%per_focus_asym(2,1) - this%per_focus_asym(1,1))

    ! get rayleigh range
    z0 = k * w0_asym**2 * 0.5_p_double

    ! calculate gaussian beam parameters
    z  = x1 - focus_asym

    Wl2 = 1.0_p_double + (z/z0)**2

    gouy_shift = atan2( z, z0 )

    ! note that this is different in 2D and 3D
    per_envelope_2d = sqrt( 1.0_p_double /  sqrt(Wl2) ) * &
    exp( - rho2/(w0_asym**2 * Wl2) ) * &
    cos( k * (x1 - z_center - prop_sign * t ) + &
      0.5_p_double*k*rho2*z/(z**2 + z0**2) - &
      gouy_shift + this%phase0 )

  case default

    ! This should never happen
    per_envelope_2d = 0

  end select

end function per_envelope_2d_zpulse
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       returns the value of the perpendicular
!       envelope at the requested position
!-----------------------------------------------------------------------------------------
function per_envelope_3d_zpulse( this, x1_sim, x2, x3, tin_sim ) result( per_envelope_3d )

  implicit none

  integer, parameter :: rank = 3

  real(p_double) :: per_envelope_3d

  class( t_zpulse ), intent(in) :: this

  ! x1, x2, x3 are organized as lon, per1, per2
  real(p_double), intent(in) :: x1_sim, x2, x3
  real(p_double), intent(in), optional :: tin_sim

  ! lab-frame x1 and t
  real(p_double) :: x1, t

  ! simulation time 
  real(p_double) :: t_sim 

  ! local variables
  real(p_double) :: z_center, rx_center, ry_center, prop_sign
  real(p_double) :: kbessel, k, z, rWl2, z0, rho2, phi, ox, oy
  real(p_double) :: gouy_shift, rWu, curv, theta
  real(p_single) :: ktr

  ! variables for astigmatic gaussian beams
  real(p_double) :: r_wx, r_wy, wx, wy, zx, zy, z0x, z0y, r_Rx, r_Ry, ox2, oy2
  real(p_double) :: gouy_angx, gouy_angy

  ! variables for asymmetric gaussian beams
  real(p_double), dimension(2) :: w0_asym, focus_asym, shp

  integer :: i

  if( present( tin_sim ) ) then 
    t_sim = tin_sim
  else 
    t_sim = 0.0_p_double
  endif 

  if( this % if_boost ) then 
    call this % apply_lorentz_transformation( x1_sim, t_sim, x1, t )
  else 
    x1 = x1_sim 
    t = t_sim 
  endif 

  ! get pulse center for chirp and phase calculations
  z_center  = this % lon_center()
  rx_center = this%per_center(1)
  ry_center = this%per_center(2)

  !calculate axial distances
  ox = x2 - rx_center
  oy = x3 - ry_center

  ! get wavenumber
  k  = this%omega0

  ! add chirp
  if ( this%propagation == p_forward ) then
    do i = 1, this%chirp_order
      k = k + this%chirp_coefs(i) * ((x1 - z_center)**(i))
    enddo
    do i = 1, this%per_chirp_order
      k = k + this%per_chirp_coefs(i) * (ox**(i))
    enddo
  else
    do i = 1, this%chirp_order
      k = k + this%chirp_coefs(i) * ((z_center - x1)**(i))
    enddo
    do i = 1, this%per_chirp_order
      k = k + this%per_chirp_coefs(i) * ((-ox)**(i))
    enddo
  endif

  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif
  
  select case ( this%per_type )

  case (p_plane)
    per_envelope_3d = cos( k * (x1-z_center - prop_sign * t ) + this%phase0 )

  case (p_hermite_gaussian)

    !calculate radius
    rho2 = ox**2 + oy**2

    z  = x1 - this%per_focus(1)

    z0 = k * this%per_w0(1)**2 /2

    ! In some situations (i.e. chirps) z0 may be 0 so z = 0 must be treated
    ! as a special case
    if ( z /= 0 ) then
      rWl2 = z0**2 / (z0**2 + z**2)
      curv = 0.5_p_double*rho2*z/(z**2 + z0**2)
      gouy_shift = (this%per_tem_mode(1) + this%per_tem_mode(2) + 1) * atan2( z, z0 )
    else
      rWl2 = 1
      curv = 0
      gouy_shift = 0
    endif

    rWu = sqrt(2 * rWl2) / this%per_w0(1)

    ! for mode(1) = 0 and mode(2) = 0 this reproduces the gaussian beam
    per_envelope_3d = sqrt(rWl2) * &
    hermite( this%per_tem_mode(1), ox * rWu ) * &
    hermite( this%per_tem_mode(2), oy * rWu ) * &
    exp( - rho2 * rWl2 /this%per_w0(1)**2 ) * &
    cos( k * ( x1 - z_center - prop_sign * t ) + k * curv - gouy_shift + this%phase0 )

  case (p_laguerre_gaussian)

    !calculate radius
    rho2 = ox**2 + oy**2

    !calculate phi
    theta = atan2(oy,ox)

    z  = x1 - this%per_focus(1)

    z0 = k * this%per_w0(1)**2 /2

    ! In some situations (i.e. chirps) z0 may be 0 so z = 0 must be treated
    ! as a special case
    if ( z /= 0 ) then
      rWl2 = z0**2 / (z0**2 + z**2)
      curv = 0.5_p_double*rho2*z/(z**2 + z0**2)
      gouy_shift = ( 2*this%per_tem_mode(1) + this%per_tem_mode(2) + 1) * atan2( z, z0 )
    else
      rWl2 = 1
      curv = 0
      gouy_shift = 0
    endif

    rWu = sqrt(2 * rWl2) / this%per_w0(1)

    ! for mode(1) = 0 and mode(2) = 0 this reproduces the gaussian beam
    ! if (.not. this%if_boost) then
    per_envelope_3d = sqrt(rWl2) * &
    ( sqrt(rho2)*sqrt(rWl2)/this%per_w0(1) )**this%per_tem_mode(2) * &
    laguerre( this%per_tem_mode(1), this%per_tem_mode(2), &
      2*rho2*rWl2/this%per_w0(1)**2  ) * &
    exp( - rho2 * rWl2/this%per_w0(1)**2  ) * &
    cos( k * ( x1 - z_center - prop_sign * t ) + k * curv - gouy_shift + this%phase0 + &
      this%per_tem_mode(2) * theta)


  case (p_hermite_gaussian_astigmatic )

    ! get rayleigh ranges
    z0x = k * this%per_w0(1)**2 *0.5_p_double ! Rayleigh range for x
    z0y = k * this%per_w0(2)**2 *0.5_p_double ! Rayleigh range for y

    ! calculate transverse coordinates squared
    ox = x2 - rx_center
    oy = x3 - ry_center
    ox2 = ox**2
    oy2 = oy**2

    ! calculate distance to focal planes
    zx  = x1 - this%per_focus(1) ! distance to x focal plane
    zy  = x1 - this%per_focus(2) ! distance to y focal plane

    ! In some situations (i.e. chirps) z0x or z0y may be 0 so zx = 0 and zy = 0 must be treated
    ! as a special case

    if ( zx /= 0 ) then
      r_wx = z0x / sqrt( zx**2 + z0x**2 ) / this%per_w0(1) ! 1/local spot, x direction
      r_Rx = zx/(zx**2 + z0x**2) ! 1 / ( x curvature radius )
      gouy_angx = atan2( zx, z0x )
    else
      r_wx = 1/this%per_w0(1)
      r_Rx = 0
      gouy_angx = 0
    endif

    if ( zy /= 0 ) then
      r_wy = z0y / sqrt( zy**2 + z0y**2 ) / this%per_w0(2) ! 1/local spot, y direction
      r_Ry = zy/(zy**2 + z0y**2) ! 1 / ( y curvature radius )
      gouy_angy = atan2( zy, z0y )
    else
      r_wy = 1/this%per_w0(2)
      r_Ry = 0
      gouy_angy = 0
    endif

    ! gouy phase shift
    ! Half of the Guoy phase comes from each transverse coord. See
    ! Siegman, Lasers, 16.4 Higher Order Gaussian Modes, Astigmatic
    ! mode functions
    gouy_shift = 0.5_p_double*( gouy_angx + gouy_angy ) * &
    (this%per_tem_mode(1) + this%per_tem_mode(2) + 1)

    ! asymmetric gaussian envelope
    per_envelope_3d = sqrt(this%per_w0(1) * this%per_w0(2) * r_wx * r_wy) * &
    hermite( this%per_tem_mode(1), real(sqrt_2,p_double) * ox * r_wx ) * &
    hermite( this%per_tem_mode(2), real(sqrt_2,p_double) * oy * r_wy ) * &
    exp( - ox2 * r_wx**2 - oy2 * r_wy**2 ) * &
    cos( k * ( x1 - z_center - prop_sign * t ) + &
      0.5_p_double*k*(ox2*r_Rx + oy2*r_Ry) - &
      gouy_shift + this%phase0 )


  case (p_bessel)
    ! calculate axial distances
    ox = x2 - rx_center
    oy = x3 - ry_center

    ! calculate Kt * rho
    ! currently bessel functions only work in single precision
    ktr = real(sqrt(ox**2+oy**2) * this%per_kt, p_single)

    ! get kbessel
    kbessel = sqrt(k**2 + this%per_kt**2)

    ! clip to the specified 0
    if (ktr >= this%per_clip_pos) then
      per_envelope_3d = 0.0_p_double
    else
      !kbessel = k * sqrt(1.0 - (kt/k)**2), calculated at setup
      phi = this%per_n * atan2(oy,ox) + kbessel * (x1 - z_center)
      per_envelope_3d = besseljn(this%per_n, ktr)*cos( phi + this%phase0 )
    endif

  case (p_gaussian_asym)

    ! THIS DOES NOT WORK WITH CHIRPED PULSES
    ! (there is a possibility of a division by 0 in that case)

    ! calculate transverse coordinates squared
    ox = x2 - rx_center
    oy = x3 - ry_center
    ox2 = ox**2
    oy2 = oy**2

    ! detect which quadrant we are in:
    ! get w0 and focsus as a smooth transition between the 2 regimes
    shp(1) = smooth_heaviside( ox/this%per_asym_trans(1) )
    shp(2) = smooth_heaviside( oy/this%per_asym_trans(2) )

    w0_asym(1) = this%per_w0_asym(1, 1) + shp(1)*( this%per_w0_asym(2, 1) - this%per_w0_asym(1, 1))
    focus_asym(1) = this%per_focus_asym(1,1) + &
    shp(1)*(this%per_focus_asym(2,1) - this%per_focus_asym(1,1))
    w0_asym(2) = this%per_w0_asym(1, 2) + shp(2)*( this%per_w0_asym(2, 2) - this%per_w0_asym(1, 2))
    focus_asym(2) = this%per_focus_asym(1,2) + &
    shp(1)*(this%per_focus_asym(2,2) - this%per_focus_asym(1,2))

    ! get rayleigh ranges
    z0x = k * w0_asym(1)**2 *0.5_p_double ! Rayleigh range for x
    z0y = k * w0_asym(2)**2 *0.5_p_double ! Rayleigh range for y

    ! calculate distance to focal planes
    zx  = x1 - focus_asym(1) ! distance to x focal plane
    zy  = x1 - focus_asym(2) ! distance to y focal plane

    wx = w0_asym(1) *sqrt( 1.0_p_double + (zx/z0x)**2 ) ! local spot, x direction
    wy = w0_asym(2) *sqrt( 1.0_p_double + (zy/z0y)**2 ) ! local spot, y direction

    r_Rx = zx/(zx**2 + z0x**2) ! 1 / ( x curvature radius )
    r_Ry = zy/(zy**2 + z0y**2) ! 1 / ( y curvature radius )

    ! gouy phase shift
    ! Half of the Guoy phase comes from each transverse coord. See
    ! Siegman, Lasers, 16.4 Higher Order Gaussian Modes, Astigmatic
    ! mode functions
    gouy_shift = 0.5_p_double*(atan2( zx, z0x ) + atan2( zy, z0y ))

    ! asymmetric gaussian envelope
    per_envelope_3d = sqrt( w0_asym(1) * w0_asym(2) / (wx * wy)) * &
    exp( - ox2/(wx**2) - oy2/(wy**2) ) * &
    cos( k * ( x1 - z_center - prop_sign * t ) + 0.5_p_double*k*(ox2*r_Rx + oy2*r_Ry) - &
      gouy_shift + this%phase0 )

  case default

    ! This should never happen
    per_envelope_3d = 0

  end select


end function per_envelope_3d_zpulse
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
function lon_center_zpulse( this ) result( lon_center )
!-----------------------------------------------------------------------------------------
!  find the center of the pulse longitudinally
!-----------------------------------------------------------------------------------------

  implicit none

     real(p_double) :: lon_center
     class( t_zpulse ), intent(in)  :: this

   select case ( this%lon_type )

     case (p_polynomial, p_sin2)

         ! the center of the pulse is define as being
         ! half way in the flat region
     if ( this%propagation == p_forward ) then
      lon_center = this%lon_start - this%lon_rise - this%lon_flat/2.0
     else
      lon_center = this%lon_start + this%lon_rise + this%lon_flat/2.0
     endif

     case (p_gaussian)

     ! for this longitudinal profile the center is
     ! user defined

     lon_center = this%lon_x0

     case (p_func)

     ! for this longitudinal profile the center is
     ! user defined

     lon_center = this%lon_x0

     case default

       lon_center = 0.0_p_double

   end select

end function lon_center_zpulse
!-----------------------------------------------------------------------------------------


! return the lab-frame coordinates given coordinates 
! in a boosted frame moving at speed beta relative to the lab frame
subroutine apply_lorentz_transformation_zpulse( this, xsim, tsim, xlab, tlab )
  class( t_zpulse ), intent(in)  :: this
  real(p_double), intent(in) :: xsim, tsim
  real(p_double), intent(out) :: xlab, tlab

  xlab = this % gamma * xsim + this % gamma_times_beta * tsim 
  tlab = this % gamma * tsim + this % gamma_times_beta * xsim 

end subroutine apply_lorentz_transformation_zpulse


!-----------------------------------------------------------------------------------------
! ** EXPERIMENTAL **
! Electric current based antennae
!-----------------------------------------------------------------------------------------

#if 0

!-----------------------------------------------------------------------------------------
! Get weights and attenuation for moving antenna
!-----------------------------------------------------------------------------------------
subroutine movwall_param( omega0, dx, s, att )

  implicit none

  real(p_double), intent(in) :: omega0, dx
  real(p_double), dimension( -1 - 2*smooth_order : 2 + 2*smooth_order ), intent(inout) :: s
  real(p_double), intent(out) :: att

  integer :: m, i, k
  real(p_double), dimension( -2 - 2*smooth_order : 3 + 2*smooth_order ) :: tmp
  real(p_double) :: c0, c1, f0

  ! Get wire shape. Cubic interpolation gives good results
  s = 0
  s(-1) = -(-1 + dx)**3/6.
  s(0)  = (4 - 6*dx**2 + 3*dx**3)/6.
  s(1)  = (1 + 3*dx + 3*dx**2 - 3*dx**3)/6.
  s(2)  = dx**3/6.

  ! Smooth current
  m = 2 * smooth_order - 1
  tmp = 0

  do k = 0, smooth_order-1
  do i = -2 - 2*k, 3 + 2*k
    tmp(i) = 0.25*s(i-1) + 0.5*s(i) + 0.25*s(i+1)
  enddo
  if ( k < smooth_order-1 ) then
    do i = -3 - 2*k, 4 + 2*k
    s(i) = 0.25*tmp(i-1) + 0.5*tmp(i) + 0.25*tmp(i+1)
    enddo
  else
    ! In the final pass apply a compensator
      c0 = (2.0 + m)/2.0
      c1 = - m / 4.0
    do i = -3 - 2*k, 4 + 2*k
    s(i) = c1*tmp(i-1) + c0*tmp(i) + c1*tmp(i+1)
    enddo
    endif
  enddo

  ! Get attenuation of wire shape + filter at laser frequency
  att = (sin( omega0 ) / ( omega0 ) ) ** 4        ! wire shape
  att = att * ( 0.5*( 1 + cos( omega0 ) ) )**m    ! binomial filters
  att = att * 0.5 * ( 2 + m - m * cos( omega0 ) ) ! compensator

end subroutine movwall_param
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Launch pulse from a moving wall
!-----------------------------------------------------------------------------------------
subroutine launch_current_movwall_1d( this, jay, g_space, nx_p_min, t )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  integer, parameter :: rank = 1

  class( t_zpulse ), pointer :: this
  type( t_vdf ), intent(inout) :: jay
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double) :: t, dt


  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(2), parameter :: t_range = 0
  real(p_double), dimension(rank) :: ldx

  real(p_double) :: zmin

  real(p_double) :: cos_pol, sin_pol
  real(p_double) :: amp, pha

  integer :: ipos
  real(p_double) :: t_0, duration

  integer :: i
  real(p_double) :: x, dx, att
  real(p_double), dimension( -1 - 2*smooth_order : 2 + 2*smooth_order ) :: s

  integer :: gix

  duration = this%lon_rise + this%lon_flat + this%lon_fall

  ! Turn off antenna if pulse has ended
  if ( t > duration ) then
    this%if_launch = .false.
    return
  endif

  call get_x_bnd( g_space, g_x_range )
  ldx(1) = real( jay%dx_(1), p_double )
  zmin = real( g_x_range(p_lower,1), p_double )

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1)*0.5_p_double
  end select

  x   = ( this%wall_pos - this%wall_vel*t - zmin ) / ldx(1)
  gix = int(x)

  i     = gix - nx_p_min( 1 ) + 1
  dx    = x - gix

  ! Get current parameters
  call movwall_param( this%omega0 * ldx(1), dx, s, att )

  ! Calculate current at each position, compensating for attenuation
  t_0 = t * (1 + this%wall_vel)
  amp =  4 * this%omega0 * this%a0 / ldx(1) / att
  pha = amp*lon_envelope( this, -t_0  , t_range ) * cos( this%omega0*t_0 + this%phase0 )

  cos_pol = cos( this%pol ) * pha
  sin_pol = sin( this%pol ) * pha

  do ipos = -1 - 2*smooth_order, 2 + 2*smooth_order
  if ( i+ipos >= -1 .and. i+ipos <= jay%nx_(1) + 2 ) then
    jay%f1(2, i+ipos)  = jay%f1(2, i+ipos)  + s(ipos) * cos_pol
    jay%f1(3, i+ipos)  = jay%f1(3, i+ipos)  + s(ipos) * sin_pol
  endif
  enddo

end subroutine launch_current_movwall_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Launch pulse from a moving wall
!-----------------------------------------------------------------------------------------
subroutine launch_current_movwall_2d( this, jay, g_space, nx_p_min, t )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  integer, parameter :: rank = 2

  class( t_zpulse ), pointer :: this
  type( t_vdf ), intent(inout) :: jay
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double) :: t, dt


  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range

  real(p_double), dimension(rank) :: xmin0

  real(p_double), dimension(2), parameter :: t_range = 0
  real(p_double), dimension(rank) :: ldx

  real(p_double) :: zmin, rmin

  real(p_double) :: cos_pol, sin_pol
  real(p_double) :: amp, pha

  integer :: ipos, i1, i2
  real(p_double) :: t_0, duration

  real(p_double) :: x, dx, att
  real(p_double), dimension( -1 - 2*smooth_order : 2 + 2*smooth_order ) :: s

   real(p_double) :: r, r_2, inj_time, j2, j3

  integer :: gix

  duration = this%lon_rise + this%lon_flat + this%lon_fall

  ! Turn off antenna if pulse has ended
  if ( t > duration ) then
    this%if_launch = .false.
    return
  endif

  call get_x_bnd( g_space, g_x_range )
  xmin0 = xmin_initial( g_space )

  ldx(1) = real( jay%dx_(1), p_double )
  zmin = real( g_x_range(p_lower,1), p_double )
  rmin = real( g_x_range(p_lower,2), p_double )

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1)*0.5_p_double
    rmin = rmin + jay%dx_(2)*0.5_p_double
  end select

  x   = ( this%wall_pos - this%wall_vel*t - zmin ) / ldx(1)
  gix = int(x)
  dx  = x - gix
  i1  = gix - nx_p_min( 1 ) + 1

  ! Get current parameters
  call movwall_param( this%omega0 * ldx(1), dx, s, att )

  ! Calculate current at each position, compensating for attenuation
  t_0 = t * (1 + this%wall_vel)
  amp =  4 * this%omega0 * this%a0 / ldx(1) / att * lon_envelope( this, -t_0  , t_range )

  cos_pol = cos( this%pol ) * amp
  sin_pol = sin( this%pol ) * amp


  g_x_range(p_lower, 1) = 0
  g_x_range(p_upper, 1) = 0

  do i2 = lbound( jay%f2, 3 ), ubound( jay%f2, 3 )

  inj_time = -t + this%wall_pos
  r = rmin + ( i2 - 1 + nx_p_min( 2 ) ) * jay%dx_(2)
  r_2 = r + 0.5 * jay%dx_(2)

  j2 = per_envelope_2d( this, inj_time, r_2 ) * cos_pol
  j3 = per_envelope_2d( this, inj_time, r   ) * sin_pol

  do ipos = -1 - 2*smooth_order, 2 + 2*smooth_order
    if ( i1+ipos >= -1 .and. i1+ipos <= jay%nx_(1) + 2 ) then
      jay%f2(2, i1+ipos, i2)  = jay%f2(2, i1+ipos, i2)  + s(ipos) * j2
      jay%f2(3, i1+ipos, i2)  = jay%f2(3, i1+ipos, i2)  + s(ipos) * j3
    endif
  enddo
  enddo


end subroutine launch_current_movwall_2d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_current_wall_1d( this, jay, t, dt, no_co )
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  class( t_zpulse ), pointer     :: this
  type( t_vdf ) ,  intent(inout)  :: jay      ! electric current
  real( p_double ),   intent(in)  :: t        ! simulation time
  real( p_double ),   intent(in)  :: dt       ! time step
  class( t_node_conf ), intent(in) :: no_co    ! node configuration

  ! local variables
  integer :: inj_node, inpos
  real( p_double ) :: inj_time
  real( p_double ) :: amp
  real(p_double)  :: ldx, rdtdx
  real(p_double) :: cos_pol, sin_pol

  ! executable statements

  ! check if correct node for launching particles
  if (this%propagation == p_forward) then
    inj_node = 1
    inpos = 1
  else
    inj_node = nx( no_co, this%direction )
    inpos = nx( jay, this%direction)
  endif

  if ( my_ngp( no_co, this%direction) == inj_node ) then

  ldx = dx(jay, this%direction)
  rdtdx = real( dt/ldx, p_k_fld )

  inj_time = real(t, p_double) - this%launch_time

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = 2.0_p_k_fld * this%omega0 * this%a0 &
        * fenv_poly( inj_time, this%lon_rise, this%lon_flat, this%lon_fall) &
        * rdtdx * cos( this%omega0*inj_time + this%phase0)

  jay%f1(2, inpos) = jay%f1(2, inpos) + amp * sin_pol
  jay%f1(3, inpos) = jay%f1(3, inpos) + amp * cos_pol

  endif

end subroutine launch_current_wall_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Launch EM pulses by setting an electric current
!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_list_current( list, jay, g_space, &
                                nx_p_min, g_nx, t, dt, no_co )
!-----------------------------------------------------------------------------------------

  use m_current

  class( t_zpulse_list ), intent(inout) :: list
  class( t_current ) ,  intent(inout) :: jay

  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t
  real(p_double), intent(in) :: dt
  class( t_node_conf ), intent(in) :: no_co

!       local variables

  class( t_zpulse ), pointer :: pulse

!       executable statements

   pulse => list%head

   do
   if (.not. associated(pulse)) exit

   select case ( pulse%type )

!    case ( p_zpulse_box )
!     call launch_zpulse( pulse, emf, g_space, nx_p_min, g_nx, t )

     case ( p_zpulse_wall )

      if ( p_x_dim /= 1 ) then
        ERROR('Only implemented in 1D')
        call abort_program( p_err_notimplemented )
      endif
      call launch_current_wall_1d( pulse, jay%pf(1), t, dt, no_co )

     case ( p_zpulse_mov_wall )

      select case ( p_x_dim )
        case (1)
          call launch_movwall_1d( pulse, jay%pf(1), g_space, nx_p_min, t )
        case (2)
          call launch_movwall_2d( pulse, jay%pf(1), g_space, nx_p_min, t )
        case default
          ERROR('Not implemented yet')
          call abort_program( p_err_notimplemented )

          end select

   end select


   pulse => pulse%next
   enddo


end subroutine launch_zpulse_list_current
!-----------------------------------------------------------------------------------------

#endif



end module m_zpulse_std
