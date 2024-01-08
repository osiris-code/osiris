#include "os-config.h"
#include "os-preprocess.fpp"

module m_zpulse_point

use m_parameters
use m_zpulse_std
use m_zpulse_wall
use m_node_conf

use m_space
use m_emf_define, only: t_emf_bound, t_emf
use m_grid_define

use m_fparser
use m_math

implicit none

private

type, extends( t_zpulse_wall ) :: t_zpulse_point

  real(p_double), dimension(p_x_dim)  :: point_pos
  real(p_double), dimension(3)        :: pol_vector

  contains

  procedure :: launch     => launch_zpulse_point
  procedure :: read_input => read_input_zpulse_point

end type t_zpulse_point

public :: t_zpulse_point

contains


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine read_input_zpulse_point( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

  use m_input_file

  implicit none

  class( t_zpulse_point ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  type( t_space ), intent(in) :: g_space
  class (t_emf_bound), intent(in) :: bnd_con
  logical, dimension(:), intent(in) :: periodic
  class( t_grid ), intent(in) :: grid
  type( t_options ), intent(in) :: sim_options

  real(p_double) :: a0, omega0, phase

  integer        :: chirp_order
  real(p_double), dimension(p_max_chirp_order) :: chirp_coefs

  character(len=16) :: tenv_type

  real(p_double) :: tenv_fwhm

  real(p_double) :: tenv_rise, tenv_flat, tenv_fall
  real(p_double) :: tenv_duration, tenv_range
  character(len = p_max_expr_len) :: tenv_math_func

  logical        :: if_launch
  real(p_double) :: launch_time, pol_norm
  real(p_double), dimension( p_x_dim ) :: point_pos

  real(p_double), dimension(3) :: pol_vector

  namelist /nl_zpulse_point/ &
    if_launch, a0, omega0, phase, pol_vector, chirp_order, chirp_coefs, &
    tenv_type, tenv_fwhm, tenv_rise, tenv_flat, tenv_fall, &
    tenv_duration, tenv_range, tenv_math_func, launch_time, point_pos

  integer :: i, ierr

  if (disp_out(input_file)) then
    SCR_ROOT(" - Reading zpulse_point configuration...")
  endif

  if_launch = .true.

  a0             = 1.0_p_double
  omega0         = 10.0_p_double
  phase          = 0.0_p_double

  chirp_order    = 0
  chirp_coefs    = 0.0_p_double

  tenv_type       = "polynomial"

  tenv_fwhm       = -1.0_p_double

  tenv_rise       = 0.0_p_double
  tenv_flat       = 0.0_p_double
  tenv_fall       = 0.0_p_double

  tenv_duration   = 0.0_p_double
  tenv_range      = -huge(1.0_p_double)

  tenv_math_func  = "NO_FUNCTION_SUPPLIED!"

  launch_time = 0.0
  point_pos      = -huge(1.0_p_double)
  pol_vector(1:3) = (/ 0.0d0, 0.0d0, 1.0d0 /)


  ! Get namelist text from input file
  call get_namelist( input_file, "nl_zpulse_point", ierr )

  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      if ( ierr < 0 ) then
        print *, "Error reading zpulse_point parameters"
      else
        print *, "Error: zpulse_point parameters missing"
      endif
      print *, "aborting..."
    endif
    stop
  endif

  ! read data from file
  read (input_file%nml_text, nml = nl_zpulse_point, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error reading zpulse_point parameters"
      write(0,*)  "   aborting..."
    endif
    stop
  endif

  this%a0      = a0
  this%omega0  = omega0
  this%phase0  = phase * pi_180


  if ((chirp_order < 0) .or. (chirp_order > p_max_chirp_order)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_point parameters"
      write(0,*)  "   chirp_order must be in the range [0,",p_max_chirp_order,"]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%chirp_order    = chirp_order
  this%chirp_coefs    = chirp_coefs

  if (( chirp_order < 0) .or. ( chirp_order > p_max_chirp_order)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_point parameters"
      write(0,*)  "   per_chirp_order must be in the range [0,",p_max_chirp_order,"]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%chirp_order    = chirp_order
  this%chirp_coefs    = chirp_coefs


  ! Temporal profile

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
        write(0,*)  '   Error in zpulse_point parameters'
        write(0,*)  '   When using "gaussian" longitudinal envelopes, '
        write(0,*)  '   lon_range must also be defined'
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
        write(0,*)  '   Error in zpulse_point parameters'
        write(0,*)  '   When using "math" longitudinal envelopes, '
        write(0,*)  '   lon_range must also be defined'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

    call setup(this%lon_math_func, trim(tenv_math_func), (/'t'/), ierr)

    ! check if function compiled ok

    if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_point parameters'
        write(0,*)  '   Supplied tenv_math_func failed to compile'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

  case default
    if ( mpi_node() == 0 ) then
      write(0,*)  ''
      write(0,*)  '   Error in zpulse_point parameters'
      write(0,*)  '   lon_type must be either "polynomial", "gaussian",'
      write(0,*)  '   "sin2", "math" or "hermite"'
      write(0,*)  '   aborting...'
    endif
    stop
  end select

  this % if_launch   = if_launch
  this % launch_time = launch_time

  ! point source parameters
  do i = 1, p_x_dim
    if ( point_pos(i) == -huge(1.0_p_double) ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ""
        write(0,*)  "   Error in zpulse_point parameters"
        write(0,*)  "   The point position (point_pos) must be set"
        write(0,*)  "   aborting..."
      endif
      stop
    endif
  enddo

  this%point_pos   = point_pos

  pol_norm = pol_vector(1)**2 + pol_vector(2)**2 + pol_vector(3)**2
  if ( pol_norm == 0 ) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_point parameters"
      write(0,*)  "   pol_vector must not be a null vector."
      write(0,*)  "   aborting..."
    endif
    stop
  endif

  ! Normalize and store polarization vector
  this % pol_vector = pol_vector / pol_norm



end subroutine read_input_zpulse_point
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_point( this, emf, g_space, nx_p_min, g_nx, t, dt, no_co )

  implicit none

  class( t_zpulse_point ), intent(inout) :: this
  class( t_emf ) ,  intent(inout)  :: emf
  type( t_space ),     intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real( p_double ),   intent(in)  :: t
  real( p_double ),   intent(in)  :: dt
  class( t_node_conf ), intent(in) :: no_co

  ! local variables
  real(p_double), dimension(p_x_dim) :: xmin0
  real(p_double) :: x1, x2, x3
  real(p_double) :: dx1, dx2, dx3
  integer :: gix1, gix2, gix3
  integer :: i1, i2, i3

  real( p_double ) :: amp, inj_time

  if ( this%if_launch .and. ( t > this%launch_time + this % t_duration() ) ) then
    this%if_launch = .false.
  endif

  if ( this%if_launch .and. ( t >= this%launch_time ) ) then

    xmin0 = xmin_initial( g_space )

    ! call the routine with the appropriate dimensions

    select case (p_x_dim)
    case(1)
      x1 = ( this%point_pos(1) - xmin0(1) ) / emf%e%dx_(1)
      gix1 = int(x1)
      dx1 = x1 - gix1

      i1  = gix1 - nx_p_min( 1 ) + 1

      if ( i1 >= 1 .and. i1 <= emf%e%nx_(1) ) then

        inj_time = t - this%launch_time

        amp = 2.0_p_k_fld * this%omega0 * this%a0 * this % t_envelope( inj_time ) * &
        cos( this%omega0*inj_time + this%phase0)

        emf%e%f1(1, i1) = emf%e%f1(1, i1) + real( amp * this % pol_vector(1), p_k_fld )
        emf%e%f1(2, i1) = emf%e%f1(2, i1) + real( amp * this % pol_vector(2), p_k_fld )
        emf%e%f1(3, i1) = emf%e%f1(3, i1) + real( amp * this % pol_vector(3), p_k_fld )

      endif

    case(2)

      x1 = ( this%point_pos(1) - xmin0(1) ) / emf%e%dx_(1)
      x2 = ( this%point_pos(2) - xmin0(2) ) / emf%e%dx_(2)

      gix1 = int(x1)
      gix2 = int(x2)

      dx1 = x1 - gix1
      dx2 = x2 - gix2

      i1  = gix1 - nx_p_min( 1 ) + 1
      i2  = gix2 - nx_p_min( 2 ) + 1

      if ( ( i1 >= 1 .and. i1 <= emf%e%nx_(1) ) .and. &
           ( i2 >= 1 .and. i2 <= emf%e%nx_(2) ) ) then

        inj_time = t - this%launch_time

        amp = 2.0_p_k_fld * this%omega0 * this%a0 * this % t_envelope( inj_time ) * &
        cos( this%omega0*inj_time + this%phase0)


        emf%e%f2(1, i1, i2) = emf%e%f2(1, i1, i2) + real( amp * this % pol_vector(1), p_k_fld )
        emf%e%f2(2, i1, i2) = emf%e%f2(2, i1, i2) + real( amp * this % pol_vector(2), p_k_fld )
        emf%e%f2(3, i1, i2) = emf%e%f2(3, i1, i2) + real( amp * this % pol_vector(3), p_k_fld )

      endif

    case(3)

      x1 = ( this%point_pos(1) - xmin0(1) ) / emf%e%dx_(1)
      x2 = ( this%point_pos(2) - xmin0(2) ) / emf%e%dx_(2)
      x3 = ( this%point_pos(3) - xmin0(3) ) / emf%e%dx_(3)

      gix1 = int(x1)
      gix2 = int(x2)
      gix3 = int(x3)

      dx1 = x1 - gix1
      dx2 = x2 - gix2
      dx3 = x3 - gix3

      i1  = gix1 - nx_p_min( 1 ) + 1
      i2  = gix2 - nx_p_min( 2 ) + 1
      i3  = gix3 - nx_p_min( 3 ) + 1

      if ( ( i1 >= 1 .and. i1 <= emf%e%nx_(1) ) .and. &
           ( i2 >= 1 .and. i2 <= emf%e%nx_(2) ) .and. &
           ( i3 >= 1 .and. i3 <= emf%e%nx_(3) ) ) then

        inj_time = t - this%launch_time

        amp = 2.0_p_k_fld * this%omega0 * this%a0 * this % t_envelope( inj_time ) * &
        cos( this%omega0*inj_time + this%phase0)

        emf%e%f3(1, i1, i2, i3) = emf%e%f3(1, i1, i2, i3) + real( amp * this % pol_vector(1), p_k_fld )
        emf%e%f3(2, i1, i2, i3) = emf%e%f3(2, i1, i2, i3) + real( amp * this % pol_vector(2), p_k_fld )
        emf%e%f3(3, i1, i2, i3) = emf%e%f3(3, i1, i2, i3) + real( amp * this % pol_vector(3), p_k_fld )

      endif

    end select

  endif

end subroutine launch_zpulse_point
!-----------------------------------------------------------------------------------------

end module m_zpulse_point
