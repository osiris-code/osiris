#include "os-config.h"
#include "os-preprocess.fpp"

module m_zpulse_wall

use m_zpulse_std
use m_system
use m_parameters
use m_math

use m_fparser
use m_input_file

use m_node_conf
use m_space
use m_grid_define
use m_vdf_define

use m_emf_define, only: t_emf_bound, t_emf
use m_emf_bound, only: is_open

implicit none

private

type, extends(t_zpulse) :: t_zpulse_wall

  integer :: tenv_type

  real(p_double)  :: inj_time

  real(p_double)  :: tenv_rise, tenv_flat, tenv_fall
  real(p_double)  :: tenv_duration, tenv_range
  type(t_fparser) :: tenv_math_func

  logical         :: if_tenv_tilt
  real(p_double), dimension(2) :: tenv_tilt

contains

  procedure :: init       => init_zpulse_wall
  procedure :: read_input => read_input_zpulse_wall
  procedure :: launch     => launch_zpulse_wall
  procedure :: t_envelope => t_env_zpulse_wall
  procedure :: t_duration => t_duration_zpulse_wall

  ! Overriding this allows for the per_envelope routines from
  ! t_zpulse to be used here
  procedure :: lon_center => lon_center_zpulse_wall

end type t_zpulse_wall

public :: t_zpulse_wall

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine init_zpulse_wall( this, restart, t, dt, b, g_space, nx_p_min, no_co, interpolation )

  implicit none

  class( t_zpulse_wall ), intent(inout) :: this
  logical, intent(in) :: restart
  real( p_double ), intent(in) :: t, dt
  type( t_vdf ) ,  intent(in)  :: b
  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: interpolation

  logical :: if_launch

  ! Need to preserve if_launch value for this zpulse
  if_launch = this%if_launch

  call this % t_zpulse % init( restart, t, dt, b, g_space, nx_p_min, no_co, interpolation )

  this%if_launch = if_launch

  if ( this%if_launch .and. ( t > this%launch_time + this % t_duration() ) ) then
    this%if_launch = .false.
  endif

end subroutine init_zpulse_wall
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine read_input_zpulse_wall( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

  implicit none

  class( t_zpulse_wall ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  type( t_space ), intent(in) :: g_space
  class (t_emf_bound), intent(in) :: bnd_con
  logical, dimension(:), intent(in) :: periodic
  class( t_grid ), intent(in) :: grid
  type( t_options ), intent(in) :: sim_options

  ! local variables
  real(p_double) :: a0
  real(p_double) :: omega0
  real(p_double) :: phase
  integer        :: pol_type
  real(p_double) :: pol
  character(len=16) :: propagation
  integer         :: direction

  integer            :: chirp_order
  real(p_double), dimension(p_max_chirp_order) :: chirp_coefs

  ! Temporal profile
  character(len=16) :: tenv_type

  real(p_double) :: tenv_fwhm

  real(p_double) :: tenv_rise, tenv_flat, tenv_fall
  real(p_double) :: tenv_duration, tenv_range
  character(len = p_max_expr_len) :: tenv_math_func
  real(p_double), dimension(2) :: tenv_tilt

  ! Beam profile
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
  logical        :: if_launch
  integer :: pos

  namelist /nl_zpulse_wall/ &
        a0, omega0,phase, pol_type, pol, propagation, direction, &
        chirp_order, chirp_coefs, tenv_fwhm, &
        tenv_type, tenv_rise, tenv_flat, tenv_fall, &
        tenv_duration, tenv_range, tenv_math_func, tenv_tilt, &
        per_type, per_center, per_w0, per_fwhm, per_focus, &
        per_chirp_order, per_chirp_coefs, &
        per_w0_asym, per_fwhm_asym, per_focus_asym, per_asym_trans, &
        per_n, per_0clip, per_kt, per_tem_mode, &
        launch_time, if_launch

  integer :: i, ierr

  if (disp_out(input_file)) then
    SCR_ROOT(" - Reading zpulse_wall configuration...")
  endif

  if_launch = .true.

  a0             = 1.0_p_double
  omega0         = 10.0_p_double
  phase          = 0.0_p_double
  pol_type       = 0
  pol            = 90.0_p_double

  chirp_order    = 0
  chirp_coefs    = 0.0_p_double

  propagation    = "forward"
  direction      = 1

  tenv_type       = "polynomial"

  tenv_fwhm       = -1.0_p_double

  tenv_rise       = 0.0_p_double
  tenv_flat       = 0.0_p_double
  tenv_fall       = 0.0_p_double

  tenv_duration   = 0.0_p_double
  tenv_range      = -huge(1.0_p_double)
  tenv_math_func  = "NO_FUNCTION_SUPPLIED!"

  tenv_tilt       = 0.0_p_double

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


  launch_time           = 0.0_p_double


  ! read data from file
  call get_namelist( input_file, "nl_zpulse_wall", ierr )

  if ( ierr /= 0 ) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        print *, "Error reading zpulse_wall parameters"
      else
        print *, "Error: zpulse_wall parameters missing"
      endif
      print *, "aborting..."
    endif
    stop
  endif

  read (input_file%nml_text, nml = nl_zpulse_wall, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error reading zpulse_wall parameters"
      write(0,*)  "   aborting..."
    endif
    stop
  endif

  this%a0         = a0
  this%omega0     = omega0
  this%phase0     = real( phase * pi_180, p_double )

  ! Polarization parameters
  if ((pol_type <-1) .or. (pol_type >1)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_wall parameters"
      write(0,*)  "   pol must be in the range [-1,0,+1]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%pol_type  = pol_type
  this%pol = real( pol * pi_180, p_double )

  ! Frequency chirp parameters
  if ((chirp_order < 0) .or. (chirp_order > p_max_chirp_order)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_wall parameters"
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
      write(0,*)  "   Error in zpulse_wall parameters"
      write(0,*)  "   per_chirp_order must be in the range [0,",p_max_chirp_order,"]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%per_chirp_order    = per_chirp_order
  this%per_chirp_coefs    = per_chirp_coefs


  ! Propagation direction
  select case ( trim(propagation) )
  case ("forward")
    this%propagation = p_forward
  case ("backward")
    this%propagation = p_backward
  case default
    if ( mpi_node() == 0 ) then
      write(0,*)  ''
      write(0,*)  '   Error in zpulse_wall parameters'
      write(0,*)  '   propagation must be either "forward" or "backward"'
      write(0,*)  '   aborting...'
    endif
    stop
  end select

  ! Launch direction
  if ((direction < 1) .or. (direction > p_x_dim)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_wall parameters"
      write(0,*)  "   direction must be in the range [1 .. x_dim]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%direction = direction

  ! Temporal envelope tilt
  this%if_tenv_tilt = .false.
  if ( p_x_dim > 1 ) then
    do i = 1, p_x_dim-1
      if ( tenv_tilt(i) /= 0.0_p_double ) then
        if (( tenv_tilt(i) < -45.0_p_double ) .or. (tenv_tilt(i) > 45.0_p_double )) then
          if ( mpi_node() == 0 ) then
            write(0,*)  ""
            write(0,*)  "   Error in zpulse_wall parameters"
            write(0,*)  "   tenv_tilt must be in the range [ -45 , 45 ] deg"
            write(0,*)  "   aborting..."
          endif
          stop
        endif
        this%tenv_tilt(i) = tan( real( tenv_tilt(i) * pi_180, p_double ))
        if ( this%propagation == p_forward ) this%tenv_tilt(i) = -this%tenv_tilt(i)
        this%if_tenv_tilt = .true.
      else
        this%tenv_tilt(i) = 0.0_p_double
      endif
    enddo
  endif


  ! Temporal envelope parameters
  select case ( trim(tenv_type) )
  case ("polynomial")
    this%tenv_type = p_polynomial

    if ( tenv_fwhm > 0.0_p_double ) then
      this%tenv_rise       = tenv_fwhm/2.0_p_double
      this%tenv_flat       = 0.0_p_double
      this%tenv_fall       = tenv_fwhm/2.0_p_double
    else
      this%tenv_rise       = tenv_rise
      this%tenv_flat       = tenv_flat
      this%tenv_fall       = tenv_fall
    endif

  case ("sin2")
    this%tenv_type = p_sin2

    if ( tenv_fwhm > 0.0_p_double ) then
      this%tenv_rise       = tenv_fwhm/2.0_p_double
      this%tenv_flat       = 0.0_p_double
      this%tenv_fall       = tenv_fwhm/2.0_p_double
    else
      this%tenv_rise       = tenv_rise
      this%tenv_flat       = tenv_flat
      this%tenv_fall       = tenv_fall
    endif

  case ("gaussian")
    this%tenv_type = p_gaussian

    if ( tenv_fwhm > 0.0_p_double  ) then
      this%tenv_duration   = tenv_fwhm/sqrt( 2.0_p_double * log(2.0_p_double))
    else
      this%tenv_duration   = tenv_duration
    endif
    this%tenv_range      = tenv_range

    if ( tenv_range == -huge(1.0_p_double) ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_wall parameters'
        write(0,*)  '   When using "gaussian" temporal envelope, '
        write(0,*)  '   tenv_range must also be defined'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

  case ("math")
    this%tenv_type = p_func
    this%tenv_range      = tenv_range

    if ( tenv_range == -huge(1.0_p_double) ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_wall parameters'
        write(0,*)  '   When using "math" temporal envelopes, '
        write(0,*)  '   tenv_range must also be defined'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

    call setup(this%tenv_math_func, trim(tenv_math_func), (/'t'/), ierr)
    if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_wall parameters'
        write(0,*)  '   Supplied tenv_math_func failed to compile'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

  case default
    if ( mpi_node() == 0 ) then
      write(0,*)  ''
      write(0,*)  '   Error in zpulse_wall parameters'
      write(0,*)  '   lon_type must be either "polynomial", "gaussian",'
      write(0,*)  '   "sin2", "math" or "hermite"'
      write(0,*)  '   aborting...'
    endif
    stop
  end select



  ! perpendiular profile
  select case ( trim(per_type) )
  case ("plane")
    this%per_type       = p_plane

  case ("gaussian","hermite")

    this%per_type       = p_hermite_gaussian
    if ( p_x_dim > 1 ) then
      this%per_tem_mode = 0
      do i = 1, p_x_dim - 1
           this%per_tem_mode(i) = per_tem_mode(i)
           if ( per_tem_mode(i) > 7 ) then
            if (mpi_node() == 0) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_wall parameters'
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
      if (mpi_node() == 0) then
     write(0,*)  ''
     write(0,*)  '   Error in zpulse_wall parameters'
     write(0,*)  '   for gaussian pulses per_w0 or per_fwhm must be > 0'
     write(0,*)  '   aborting...'
   endif
     stop
    endif

    if ( per_focus(1) == -huge(1.0_p_double) ) then
      if (mpi_node() == 0) then
     write(0,*)  ''
     write(0,*)  '   Error in zpulse_wall parameters'
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
      do i = 1, p_x_dim - 1
           this%per_tem_mode(i) = per_tem_mode(i)
           if ( per_tem_mode(i) > 7 ) then
            if (mpi_node() == 0) then
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

      if (p_x_dim < 3) then
        if (mpi_node() == 0) then
        print *, ''
        print *, '   Error in zpulse parameters'
        print *, '   Laguerre-Gaussian modes are available in 3D geometries only'
        print *, '   aborting...'
      endif
        stop
    endif

    ! get main spot size and focal plane
    this%per_w0 = 0.0_p_double
    if (per_fwhm(1) > 0.0_p_double) then
      this%per_w0 = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
    else
      this%per_w0 = per_w0(1)
    endif

    if ( this%per_w0(1) <= 0.0_p_double ) then
      if (mpi_node() == 0) then
     print *, ''
     print *, '   Error in zpulse parameters'
     print *, '   for Laguerre-Gaussian pulses per_w0 or per_fwhm must be > 0'
     print *, '   aborting...'
   endif
     stop
    endif

    if ( per_focus(1) == -huge(1.0_p_double) ) then
      if (mpi_node() == 0) then
     print *, ''
     print *, '   Error in zpulse parameters'
     print *, '   for Laguerre-Gaussian pulses per_focus must be defined'
     print *, '   aborting...'
   endif
     stop
    endif
    this%per_focus      = per_focus(1)

    ! default parameters
    if (per_fwhm(1) > 0.0_p_double) then
      this%per_w0(1) = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
    else
      this%per_w0(1) = per_w0(1)
    endif
    this%per_focus(1)      = per_focus(1)

  case default
if (mpi_node() == 0) then
      write(0,*)  ''
      write(0,*)  '   Error in zpulse_wall parameters'
      write(0,*)  '   per_type must be one of the following:'
      write(0,*)  '   "plane", "gaussian"/"hermite"'
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

  ! do not allow wall type without open boundaries
  pos = p_lower
  if ( this%propagation == p_backward) pos = p_upper

  if ( periodic(this%direction) .or. (.not. is_open(bnd_con, this%direction, pos)) ) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_wall parameters"
      write(0,*)  "   Boundary condition must be open ('vpml' or 'lindman') in the wall launching the pulse."
      write(0,*)  "   aborting..."
    endif
    stop
  endif

  ! extra parameters
  this%launch_time = launch_time
  this%if_launch   = if_launch
  
  ! Wall pulses don't need divergence correction, this just insures that
  ! the parameter has the same value in all directions
  this%no_div_corr = .false.

  ! Use in boosted frame simulations is currently unsupported
  this%if_boost    = .false.

  ! To reuse superclass per_envelope methods we must flip the sign of per_focus
  ! when propagating backwards
  if (this%propagation == p_backward) then
    this%per_focus(1) = - this%per_focus(1)
    this%per_focus(2) = - this%per_focus(2)
  endif

end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_wall( this, emf, g_space, &
                               nx_p_min, g_nx, t, dt, no_co )

  implicit none

  class( t_zpulse_wall ), intent(inout) :: this
  class( t_emf ) ,  intent(inout) :: emf

  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t
  real(p_double), intent(in) :: dt
  class( t_node_conf ), intent(in) :: no_co

  real(p_double) :: phase0_temp, pol_temp

  if ( this%if_launch .and. ( t > this%launch_time + this % t_duration() ) ) then
    this%if_launch = .false.
  endif

  if ( this%if_launch .and. ( t >= this%launch_time ) ) then

    ! call the routine with the appropriate dimensions

    select case ( emf % b % x_dim_)
    case(1)
      call launch_em_wall_1d( this, emf%b, t, dt, no_co )

    case(2)
      call launch_em_wall_2d( this, emf%b, g_space, nx_p_min, t, dt, no_co )

    case(3)
      call launch_em_wall_3d( this, emf%b, g_space, nx_p_min, t, dt, no_co )

    end select

    if (this%pol_type /= 0) then

      ! Store phase and polarization of first pulse
      phase0_temp = this%phase0
      pol_temp    = this%pol

      ! Get phase and polarization of second pulse
      this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )) , p_double )
      this%pol    = this%pol    + real( pi_2, p_double )

      ! Launch second pulse
      select case ( emf % b % x_dim_)
      case(1)
        call launch_em_wall_1d( this, emf%b, t, dt, no_co )

      case(2)
        call launch_em_wall_2d( this, emf%b, g_space, nx_p_min, t, dt, no_co )

      case(3)
        call launch_em_wall_3d( this, emf%b, g_space, nx_p_min, t, dt, no_co )

      end select

      ! Restore phase and polarization for first pulse
      this%phase0 = phase0_temp
      this%pol = pol_temp

    endif

  endif

end subroutine launch_zpulse_wall
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
function t_env_zpulse_wall( this, t ) result(t_env)

  implicit none

  class( t_zpulse_wall ), intent(inout) :: this ! intent must be inout because of
                                                ! the eval( f_parser ) function
  real(p_double), intent(in) :: t
  real(p_double) :: t_env

  ! local variables

  real(p_k_fparse), dimension(1) ::t_dbl

  ! executable statements

  select case ( this%tenv_type )

  case (p_polynomial)
    t_env = fenv_poly( t, this%tenv_rise, this%tenv_flat, this%tenv_fall )

  case (p_sin2)
    t_env = fenv_sin2( t, this%tenv_rise, this%tenv_flat, this%tenv_fall )

  case (p_gaussian)
    t_env = exp( - 2*((t - 0.5*this % tenv_range)/this%tenv_duration)**2 )


  case (p_func)

    if ( t > this % tenv_range ) then
      t_env = 0
    else
      t_dbl(1) = t
      t_env = eval( this%lon_math_func, t_dbl )
    endif

  case default
    t_env = 0

  end select

end function t_env_zpulse_wall
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Total duration of the laser pulse
!-----------------------------------------------------------------------------------------
function t_duration_zpulse_wall( this ) result( t_duration )

  implicit none

  class( t_zpulse_wall ), intent(in) :: this
  real(p_double) :: t_duration

  select case ( this%tenv_type )
  case (p_polynomial, p_sin2)
    t_duration = this%tenv_rise + this%tenv_flat + this%tenv_fall

  case (p_gaussian, p_func)
    t_duration = this%tenv_range

  case default
    t_duration = 0

  end select

end function t_duration_zpulse_wall
!----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!  find the center of the pulse longitudinally
!  - This function is used by the per_envelope routines for phase calculations
!  - Adding this % t_inj_time allows for the use of the t_zpulse_std per_envelope methods
!-----------------------------------------------------------------------------------------

function lon_center_zpulse_wall( this ) result( lon_center )

  implicit none

  class( t_zpulse_wall ), intent(in)  :: this
  real(p_double) :: lon_center

  select case ( this%tenv_type )

  case (p_polynomial, p_sin2)
    ! the center of the pulse is define as being
    ! half way in the flat region
    lon_center = this % tenv_rise + this % tenv_flat / 2

  case (p_gaussian, p_func)
    lon_center = 0.5*this % tenv_range

  case default
    lon_center = 0

  end select

  lon_center = lon_center + this % inj_time

end function lon_center_zpulse_wall


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine launch_em_wall_1d( this, b, t, dt, no_co )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_zpulse_wall ), intent(inout) :: this
  type( t_vdf ) ,  intent(inout)  :: b           ! magnetic field
  real( p_double ),   intent(in)  :: t           ! simulation time
  real( p_double ),   intent(in)  :: dt         ! time step
  class( t_node_conf ), intent(in) :: no_co       ! node configuration

  ! local variables
  integer :: inj_node, inpos
  real( p_double ) :: inj_time, prop_sign
  real( p_double ) :: amp
  real(p_double)  :: ldx, rdtdx
  real(p_double) :: cos_pol, sin_pol

  ! executable statements

  ! check if correct node for launching particles
  if (this%propagation == p_forward) then
    inj_node = 1
    inpos = 1
    prop_sign = +1.0_p_double
  else
    inj_node = nx( no_co, this%direction )
    inpos = b %nx_( this%direction)
    prop_sign = -1.0_p_double
  endif

  if ( my_ngp( no_co, this%direction) == inj_node ) then

  ldx = b%dx(this%direction)
  rdtdx = real( dt/ldx, p_k_fld )

  inj_time = t - this%launch_time

  cos_pol = + cos( this%pol ) * prop_sign
  sin_pol = - sin( this%pol ) * prop_sign
  amp = 2.0_p_k_fld * this%omega0 * this%a0 * this % t_envelope( inj_time ) * &
           rdtdx * cos( this%omega0*inj_time + this%phase0)

  b%f1(2, inpos) = b%f1(2, inpos) + real( amp * sin_pol, p_k_fld )
  b%f1(3, inpos) = b%f1(3, inpos) + real( amp * cos_pol, p_k_fld )

  endif

end subroutine launch_em_wall_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! The trick to using the standard per_envelope routines is to use the
! lon_start /
!-----------------------------------------------------------------------------------------
subroutine launch_em_wall_2d( this, b, g_space, nx_p_min, t, dt, no_co )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  class( t_zpulse_wall ), intent(inout) :: this
  type( t_vdf ) ,  intent(inout)  :: b           ! magnetic field
  type( t_space ),     intent(in) :: g_space     ! global space information
  integer, intent(in), dimension(:) :: nx_p_min
  real( p_double ),   intent(in)  :: t           ! simulation time
  real( p_double ),   intent(in)  :: dt         ! time step
  class( t_node_conf ), intent(in) :: no_co       ! node configuration

  integer :: inj_node, inpos
  integer :: ip, iperp
  real( p_double ) :: prop_sign
  real( p_double ) :: amp, lenv, rmin, z, r, r_2, dr_2
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double)  :: ldx, ldxp, rdtdx, b12, b3
  real(p_double) :: cos_pol, sin_pol


  call get_x_bnd( g_space, g_x_range )

  if (this%propagation == p_forward) then
   inj_node = 1
   inpos = 1
   prop_sign = 1.0_p_double

   z =  g_x_range(p_lower,this%direction)

   this % inj_time = t - this%launch_time
  else
   inj_node = nx( no_co, this%direction )
   inpos = b %nx_( this%direction)
   prop_sign = -1.0_p_double

   z =  -g_x_range(p_upper,this%direction)
   this % inj_time = t - this%launch_time
  endif

  if ( my_ngp( no_co, this%direction) == inj_node ) then

  ldx = b%dx(this%direction)
  rdtdx = real( dt/ldx, p_k_fld )


  ! perpendicular direction
  iperp = 3 - this%direction
  ldxp = b%dx(iperp)

  ! cells for wall injection
  rmin =  g_x_range(p_lower,iperp) + (nx_p_min(iperp)-1)*ldxp

  cos_pol = + cos( this%pol ) * prop_sign
  sin_pol = - sin( this%pol ) * prop_sign
  amp = 2.0_p_k_fld * this%omega0 * this%a0

  dr_2 = ldxp/2.0_p_k_fld

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    z = z + ldx*0.5_p_double
    rmin = rmin + dr_2
  end select

  ! inject column of cells
  lenv = amp * this % t_envelope( this % inj_time   )

  do ip = 0, b%nx_(iperp)+2

    r = rmin + (ip-1) * ldxp
    r_2 = r + dr_2

    b12 = lenv   * rdtdx * this % per_envelope_2d( z, r   ) * sin_pol
    b3  = lenv   * rdtdx * this % per_envelope_2d( z, r_2 ) * cos_pol

    if (this%direction == 1) then
      ! launching along x1
      b%f2(2, inpos, ip) = b%f2(2, inpos, ip) + real( b12, p_k_fld )
      b%f2(3, inpos, ip) = b%f2(3, inpos, ip) + real( b3, p_k_fld )
    else
      ! launching along x2
      b%f2(1, ip, inpos) = b%f2(1, ip, inpos) + real( b12, p_k_fld )
      b%f2(3, ip, inpos) = b%f2(3, ip, inpos) + real( b3, p_k_fld )
    endif

  enddo

  endif

end subroutine launch_em_wall_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine launch_em_wall_3d( this, b, g_space, nx_p_min, t, dt, no_co )

  implicit none

  integer, parameter :: rank = 3

  class( t_zpulse_wall ), intent(inout) :: this
  type( t_vdf ) ,  intent(inout)  :: b           ! magnetic field
  type( t_space ),     intent(in) :: g_space     ! global space information
  integer, intent(in), dimension(:) :: nx_p_min
  real( p_double ),   intent(in)  :: t           ! simulation time
  real( p_double ),   intent(in)  :: dt         ! time step
  class( t_node_conf ), intent(in) :: no_co       ! node configuration

  ! local variables
  integer :: inj_node, inpos
  integer :: ip1, ip2, iperp1, iperp2
  real( p_double ) :: z, prop_sign
  real( p_double ) :: amp, lenv, rmin1, rmin2, r1, r2, r1_2, r2_2, dr1_2, dr2_2
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double)  :: ldx, ldxp1, ldxp2, rdtdx, b0, b90
  real(p_double) :: cos_pol, sin_pol

  ! executable statements

  this % inj_time = t - this%launch_time
  call get_x_bnd( g_space, g_x_range )

  if (this%propagation == p_forward) then
    inj_node = 1
    inpos = 1
    prop_sign = 1.0_p_double

    z =  g_x_range(p_lower,this%direction)
  else
    inj_node = nx( no_co, this%direction )
    inpos = b % nx_( this%direction)
    prop_sign = -1.0_p_double

    z =  -g_x_range(p_upper,this%direction)
  endif

  if ( my_ngp( no_co, this%direction) == inj_node ) then

    ldx = b%dx(this%direction)
    rdtdx = real( dt/ldx, p_k_fld )

    cos_pol =  + cos( this%pol ) * prop_sign
    sin_pol =  - sin( this%pol ) * prop_sign
    amp = 2.0_p_k_fld * this%omega0 * this%a0

    ! perpendicular directions
    select case ( this%direction )
    case (1)
      iperp1 = 2
      iperp2 = 3
    case (2)
      iperp1 = 1
      iperp2 = 3
    case (3)
      iperp1 = 1
      iperp2 = 2
    end select

    ldxp1 = b%dx(iperp1)
    ldxp2 = b%dx(iperp2)

    dr1_2 = ldxp1/2.0_p_k_fld
    dr2_2 = ldxp2/2.0_p_k_fld

    rmin1 = real( g_x_range(p_lower,iperp1), p_k_fld ) + (nx_p_min(iperp1)-1)*ldxp1
    rmin2 = real( g_x_range(p_lower,iperp2), p_k_fld ) + (nx_p_min(iperp2)-1)*ldxp2

    select case( this%interpolation )
    case( p_linear, p_cubic )
      ! Nothing to change here
    case( p_quadratic, p_quartic )
      z = z + ldx*0.5_p_double
      rmin1 = rmin1 + dr1_2
      rmin2 = rmin2 + dr2_2
    end select

    ! inject column of cells
    lenv = amp * this % t_envelope( this % inj_time )

    do ip2 = 0, b%nx_(iperp2)+2
      do ip1 = 0, b%nx_(iperp1)+2

        r1 = rmin1 + (ip1-1) * ldxp1
        r2 = rmin2 + (ip2-1) * ldxp2
        r1_2 = r1 + dr1_2
        r2_2 = r2 + dr2_2

        b0  = lenv * rdtdx * this % per_envelope_3d( z, r1,   r2_2 ) * sin_pol
        b90 = lenv * rdtdx * this % per_envelope_3d( z, r1_2, r2   ) * cos_pol

        select case (this%direction)

        case (1)
          b%f3(2, inpos, ip1, ip2) = b%f3(2, inpos, ip1, ip2) + real( b0, p_k_fld )
          b%f3(3, inpos, ip1, ip2) = b%f3(3, inpos, ip1, ip2) + real( b90, p_k_fld )

        case (2)
          b%f3(1, ip1, inpos, ip2) = b%f3(1, ip1, inpos, ip2) + real( b0, p_k_fld )
          b%f3(3, ip1, inpos, ip2) = b%f3(3, ip1, inpos, ip2) + real( b90, p_k_fld )

        case (3)
          b%f3(1, ip1, ip2, inpos) = b%f3(1, ip1, ip2, inpos) + real( b0, p_k_fld )
          b%f3(2, ip1, ip2, inpos) = b%f3(2, ip1, ip2, inpos) + real( b90, p_k_fld )

        end select

      enddo
    enddo

  endif

end subroutine launch_em_wall_3d
!-----------------------------------------------------------------------------------------


end module m_zpulse_wall

