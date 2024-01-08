! OS-ZPULSE-SPECKLE.F03
! Change Log:
! Ver 1 (circa 2018):  This is documented in Han Wen et al, PPCF 
! Ver 1.1 : added transverse Gaussian envelope (per_center and per_w0)
!           and tilt (k_perp), October, 2019
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_zpulse_speckle

#include "memory/memory.h"

use m_zpulse_wall
use m_zpulse_std
use m_system
use m_parameters
use m_math

use m_fparser

use m_node_conf
use m_space
use m_grid_define
use m_vdf_define

use m_emf_define, only: t_emf_bound, t_emf
use m_emf_bound, only: is_open
use m_random_class
use m_random_hash
use m_random
use m_fft

implicit none

private

! speckle_types
integer, parameter :: p_rpp      = 1
integer, parameter :: p_cpp      = 2
integer, parameter :: p_isi      = 3
integer, parameter :: p_ssd      = 4
integer, parameter :: p_file     = 5
integer, parameter :: p_cpurp    = 6

! SSD specific
integer, parameter :: p_max_nfm  = 8     ! maximum number of sinusoidal phase modulation

! shape of the bandwidth
integer, parameter :: p_default_spectrum      = 0   ! different spectra for different speckle_types
integer, parameter :: p_gaussian_spectrum     = 1
integer, parameter :: p_lorentzian_spectrum   = 2
integer, parameter :: p_line_spectrum         = 3

integer, parameter :: per_ndim = max(1, p_x_dim-1)

type, extends(t_zpulse_wall) :: t_zpulse_speckle

  integer         :: bandwidth_type, speckle_type, rseed_phase, rseed_stud
  integer         :: nt
  logical         :: if_dump_speckle, stud_was_on
  real (p_double) :: chirp_period, pol_stepping, lon_focus
  real (p_double) :: dt_update_speckle, dt_update_pol
  real (p_double) :: stud_period, stud_jitter, duty_cycle
  real (p_double) :: stud_on_t, stud_rdc, stud_duration
  integer, dimension(per_ndim)         :: n_speckle, per_lnx, nfm, per_dir
  real (p_double), dimension(per_ndim) :: color_cycle, beamlet_delay, spatial_tilt, pdk
  real (p_double), dimension(per_ndim) :: ar1_rho, ar1_sigma, laser_bandwidth
  real (p_double), dimension(p_max_nfm, per_ndim ) :: fm, phasemod_amplitude, pm_bw

  real (p_double), dimension(:),     pointer     :: pp_2d, pmp_2d
  real (p_double), dimension(:,:),   pointer     :: buff_1d, pp_3d, per_wall_2d, pmp_3d
  real (p_double), dimension(:,:,:), pointer     :: buff_2d, per_wall_3d, per_mode_2d
  real (p_double), dimension(:,:,:,:), pointer   :: buff_3d, per_mode_3d
!  real (p_double) :: per_center(2)
!  real (p_double) :: per_w0(2)
  real (p_double) :: k_perp(2)
  !character(len = p_max_filename_len)      :: filename

  class (t_random), pointer :: rbuff => null()   ! rng for (complex) phase modulation
  class (t_random), pointer :: rstud => null()   ! rng for amplitude modulation
contains

  procedure :: read_input => read_input_zpulse_speckle
  procedure :: launch     => launch_zpulse_speckle
  procedure :: t_envelope => t_env_zpulse_speckle

  ! zpulse speckle has internal states to keep track of. need to override init and cleanup.
  procedure :: init       => init_zpulse_speckle
  procedure :: cleanup    => cleanup_zpulse_speckle
  procedure :: t_duration => t_duration_zpulse_speckle
end type t_zpulse_speckle

interface fastforward
  module procedure fastforward_buff1d
  module procedure fastforward_buff1d2
  module procedure fastforward_lorentzian_psd
end interface

interface ar1proc
  module procedure ar1proc_scaler
  module procedure ar1proc_1d
  module procedure ar1proc_2d
end interface

public :: t_zpulse_speckle

contains


!-----------------------------------------------------------------------------------------
! convert a string to lower case
!-----------------------------------------------------------------------------------------
function to_lower(str) result(res)
  implicit none
  character(*), intent(in) :: str

  character(len=len(str))  :: res
  integer :: i

  do i = 1, len(str)
    select case(str(i:i))
      case("A":"Z")
        res(i:i) = achar(iachar(str(i:i))+32)
      case default
        res(i:i) = str(i:i)
    end select
  end do
end function to_Lower

function ar1process(val, b, sigma, prng) result(res)
  implicit none
  real (p_double),               intent(inout) :: val
  real (p_double),               intent(in)    :: b, sigma
  class(t_random), pointer,      intent(in)    :: prng
  real (p_double) :: res

  res = b * val + sigma * prng%genrand_gaussian()

end function ar1process

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ar1proc_scaler( val, b, sigma, prng )
  implicit none
  real (p_double),               intent(inout) :: val
  real (p_double),               intent(in)    :: b, sigma
  class(t_random), pointer,      intent(in)    :: prng

  val = ar1process(val, b, sigma, prng)

end subroutine ar1proc_scaler
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ar1proc_1d( val, b, sigma, n_val, prng )
  implicit none
  real (p_double), dimension(1:), intent(inout) :: val
  real (p_double),                intent(in)    :: b, sigma
  integer,                        intent(in)    :: n_val
  class(t_random), pointer,       intent (in)   :: prng

  integer :: i1

  do i1 = 1, n_val
    val(i1) = ar1process(val(i1), b, sigma, prng)
  end do

end subroutine ar1proc_1d

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ar1proc_2d( val, b, sigma, n_val, prng )
  implicit none
  real (p_double), dimension(1:,1:), intent(inout) :: val
  real (p_double),                   intent(in)    :: b, sigma
  integer, dimension(1:2),           intent(in)    :: n_val
  class(t_random), pointer,          intent (in)   :: prng

  integer :: i1, i2
  do i1 = 1, n_val(1)
    do i2 = 1, n_val(2)
      val(i1, i2) = ar1process(val(i1, i2), b, sigma, prng)
    end do
  end do

end subroutine ar1proc_2d

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine fastforward_buff1d( buff_1d, n, tn, ts, b, sigma, prng )
  implicit none
  real(p_double), dimension(1:, 1:), intent(inout) :: buff_1d
  integer,                           intent(in)    :: n, tn, ts
  real(p_double),                    intent(in)    :: b, sigma
  class(t_random), pointer,          intent(in)    :: prng

  integer :: i, ind1, ind2
  do i = 0, tn-1
    ind1 = mod(ts + i, n) + 1
    ind2 = mod(ind1 + n - 1, n) + 1
    buff_1d(1, ind1) = ar1process( buff_1d(1, ind2), b, sigma, prng )
  enddo

end subroutine fastforward_buff1d
!-------------------------------------------------------------------------------

subroutine fastforward_buff1d2( buff_1d, n, tn, ts, b, sigma, prng )
  implicit none
  real(p_double), dimension(1:, 1:),     intent(inout) :: buff_1d
  integer,                               intent(in)    :: n, tn
  integer,        dimension(1:per_ndim), intent(in)    :: ts
  real(p_double), dimension(1:per_ndim), intent(in)    :: b, sigma
  class(t_random), pointer,              intent(in)    :: prng

  integer :: i, ind1, ind2
  do i = 0, tn-1
    ind1 = mod(ts(1) + i, n) + 1
    ind2 = mod(ind1 + n - 1, n) + 1
    buff_1d(1, ind1) = ar1process( buff_1d(1, ind2), b(1), sigma(1), prng )

    ind1 = mod(ts(2) + i, n) + 1
    ind2 = mod(ind1 + n - 1, n) + 1
    buff_1d(2, ind1) = ar1process( buff_1d(2, ind2), b(2), sigma(2), prng )
  enddo

end subroutine fastforward_buff1d2
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine fastforward_lorentzian_psd( this, tn )
  implicit none
  class( t_zpulse_speckle ), intent(inout) :: this
  integer,                   intent(in)    :: tn

  integer :: ti, i, j

  select case ( p_x_dim )

  case ( 1, 2 )
    if ( this%speckle_type == p_cpurp ) then

      do ti = 1, tn
        call ar1proc( this%pp_2d(:), this%ar1_rho(1), &
                        this%ar1_sigma(1), this%n_speckle(1), this%rbuff )
      enddo
      do i = 1, this%n_speckle(1)
        this%buff_2d(1, i, 1) = cos( this%phasemod_amplitude(1,1) * this%pp_2d(i) )
        this%buff_2d(2, i, 1) = sin( this%phasemod_amplitude(1,1) * this%pp_2d(i) )
      enddo

    else

      do ti = 1, tn
        do i = 1, 2
          call ar1proc( this%buff_2d(i, :, 1), this%ar1_rho(1), &
                          this%ar1_sigma(1), this%n_speckle(1), this%rbuff )
        enddo
      enddo

    endif

  case ( 3 )
    if ( this%speckle_type == p_cpurp ) then

      do ti = 1, tn
        call ar1proc( this%pp_3d, this%ar1_rho(1), &
                        this%ar1_sigma(1), this%n_speckle(1:2), this%rbuff )
      enddo
      do i = 1, this%n_speckle(1)
        do j = 1, this%n_speckle(2)
          this%buff_3d(1, i, j, 1) = cos( this%phasemod_amplitude(1,1) * this%pp_3d(i, j) )
          this%buff_3d(2, i, j, 1) = sin( this%phasemod_amplitude(1,1) * this%pp_3d(i, j) )
        enddo
      enddo

    else

      do ti = 1, tn
        do i = 1, 2
          call ar1proc( this%buff_3d(i, :, :, 1), this%ar1_rho(1), &
                          this%ar1_sigma(1), this%n_speckle(1:2), this%rbuff )
        enddo
      enddo

    endif
  end select
end subroutine fastforward_lorentzian_psd
!-------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine new_phase_plate( this )
  implicit none
  class( t_zpulse_speckle ), intent(inout) :: this

  ! local variables
  integer :: i, j

  select case ( p_x_dim )

  case ( 1, 2 )

    do i = 1, this%n_speckle(1)

      if (this%speckle_type == p_rpp) then
        this%pp_2d(i) =  (sign(1.0_p_double, this % rbuff % genrand_gaussian()) &
                            + 1.0_p_double) * pi_2
      else
        this%pp_2d(i) =  this % rbuff % genrand_gaussian() * PI
      endif

      this%buff_2d(1, i, 1) = cos(this%pp_2d(i))
      this%buff_2d(2, i, 1) = sin(this%pp_2d(i))

    enddo

  case ( 3 )

    do i = 1, this%n_speckle(1)
      do j = 1, this%n_speckle(2)

        if (this%speckle_type == p_rpp) then
          this%pp_3d(i,j) =  (sign(1.0_p_double, this % rbuff % genrand_gaussian()) &
                                 + 1.0_p_double) * pi_2
        else
          this%pp_3d(i,j) =  this % rbuff % genrand_gaussian() * PI
        endif

        this%buff_3d(1, i, j, 1) = cos(this%pp_3d(i, j))
        this%buff_3d(2, i, j, 1) = sin(this%pp_3d(i, j))

      enddo
    enddo

  end select

end subroutine new_phase_plate
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine copy_rpmssd_buffer_2d(buff_2d, buff_1d, nt, ns, am, ti, tiu)
  real(p_double), dimension(1:,1:,1:), intent(inout) :: buff_2d
  real(p_double), dimension(1:,1:),    intent(in)    :: buff_1d
  integer,                             intent(in)    :: nt, ti, tiu
  integer,                             intent(in)    :: ns
  real(p_double),                      intent(in)    :: am

  integer :: i, ind

  do i = 1, ns
    ind = mod(ti+nt-1+(i-1)*tiu, nt) + 1
    buff_2d(1, i, 1) = cos( am * buff_1d(1, ind) )
    buff_2d(2, i, 1) = sin( am * buff_1d(1, ind) )
  enddo

end subroutine copy_rpmssd_buffer_2d

!-----------------------------------------------------------------------------------------
subroutine copy_rpmssd_buffer_3d(buff_3d, buff_1d, nt, ns, am, ti, tiu)
  real(p_double), dimension(1:,1:,1:,1:), intent(inout) :: buff_3d
  real(p_double), dimension(1:,1:),       intent(in)    :: buff_1d
  integer,                                intent(in)    :: nt, ti
  integer,        dimension(1:),          intent(in)    :: tiu
  integer,        dimension(1:),          intent(in)    :: ns
  real(p_double), dimension(1:),          intent(in)    :: am

  integer :: i, j, tj, ind1, ind2
  real(p_double) :: ph

  do i = 1, ns(1)
    tj = ti + (i - 1) * tiu(1)
    ind1 = mod(tj + nt - 1, nt) + 1

    do j = 1, ns(2)
      ind2 = mod(ind1 + nt - 1 + (j - 1) * tiu(2), nt) + 1
      ph = am(1)*buff_1d(1, ind1) + am(2)*buff_1d(2, ind2)
      buff_3d(2, i, j, 1) = cos( ph )
      buff_3d(2, i, j, 1) = sin( ph )
    enddo

  enddo
end subroutine copy_rpmssd_buffer_3d

!-----------------------------------------------------------------------------------------
! generate a time series with gaussian power density spectrum
!-----------------------------------------------------------------------------------------
subroutine new_gaussian_psd_timeseries( this )
  implicit none
  class( t_zpulse_speckle ), intent(inout) :: this

  ! local variables
  integer :: i, j, k, nt
  real(p_double) :: df, total
  real(p_double), dimension(per_ndim)   :: c1, bw
  real(p_double), dimension(:), pointer :: profile

  if ( this%speckle_type == p_ssd .or. this%speckle_type == p_cpurp) then
    bw = this%pm_bw(1, 1:per_ndim)
  else
    bw = this%laser_bandwidth
  endif

  c1 = 0.5 * alog(2.0) / (bw/2.0)**2
  nt = this%nt
  df = 1.0 / (nt * this%dt_update_speckle)
  call fft_init( nt, .true. )
  call alloc( profile, (/ nt /) )

  do i = 1, nt
    profile(i) = exp( - c1(1) * ((i-(nt-1)*0.5)*df)**2 )
  enddo

  if ( this%speckle_type == p_ssd ) then

    do i = 1, nt
      this%buff_1d(1, i) =  this% rbuff % genrand_gaussian() * profile(i)
      this%buff_1d(2, i) =  this% rbuff % genrand_gaussian() * profile(i)
    enddo

    call ifft( this%buff_1d(1:2,:), nt )

    total = 0
    do i = 1, nt
      total = total + this%buff_1d(1, i)**2
    enddo
    this%buff_1d = this%buff_1d / sqrt(total / real(nt, p_double))

  endif


  select case ( p_x_dim )
  case ( 1, 2 )

    if ( this%speckle_type == p_ssd ) then
      call copy_rpmssd_buffer_2d(this%buff_2d, this%buff_1d, nt, this%n_speckle(1), &
                      this%phasemod_amplitude(1,1), 1, &
                      floor(this%beamlet_delay(1) / this%dt_update_speckle))
    else

      do k = 1, nt
        do i = 1, this%n_speckle(1)
          this%buff_2d(1, i, k) = this% rbuff % genrand_gaussian() * profile(k)
          this%buff_2d(2, i, k) = this% rbuff % genrand_gaussian() * profile(k)
        enddo
      enddo

      call ifft( this%buff_2d, (/ nt, this%n_speckle(1) /) )

      total = 0
      do k = 1, nt
        do i = 1, this%n_speckle(1)
          total = total + this%buff_2d(1, i, k)**2 + this%buff_2d(2, i, k)**2
        enddo
      enddo
      this%buff_2d = this%buff_2d / sqrt(total / real(nt*this%n_speckle(1), p_double))

      if ( this%speckle_type == p_cpurp ) then
        do k = 1, nt
          do i = 1, this%n_speckle(1)
            this%buff_2d(1, i, k) = cos( this%phasemod_amplitude(1,1) * this%buff_2d(1,i,k) )
            this%buff_2d(2, i, k) = sin( this%phasemod_amplitude(1,1) * this%buff_2d(2,i,k) )
          enddo
        enddo
      endif

    endif

  case ( 3 )

    if ( this%speckle_type == p_ssd ) then
      if ( this%laser_bandwidth(2) > 0 ) then

        do i = 1, nt
          profile(i) = exp( - c1(2) * ((i-(nt-1)*0.5)*df)**2 )
        enddo
        do i = 1, nt
          this%buff_1d(2, i) =  this% rbuff % genrand_gaussian() * profile(i)
          this%buff_1d(3, i) =  this% rbuff % genrand_gaussian() * profile(i)
        enddo

        call ifft( this%buff_1d(2:3,:), nt )

        total = 0
        do i = 1, nt
          total = total + this%buff_1d(2, i)**2
        enddo
        this%buff_1d(2,:) = this%buff_1d(2,:) / sqrt(total / real(nt, p_double))

      else

        this%buff_1d(2, :) = 0.0_p_double

      endif
      call copy_rpmssd_buffer_3d(this%buff_3d, this%buff_1d, nt, this%n_speckle(1:per_ndim), &
                      this%phasemod_amplitude(1,1:per_ndim), 1, &
                      floor(this%beamlet_delay(1:per_ndim) / this%dt_update_speckle))

    else

      do k = 1, nt
        do j = 1, this%n_speckle(2)
          do i = 1, this%n_speckle(1)
            this%buff_3d(1, i, j, k) =  this% rbuff % genrand_gaussian() * profile(k)
            this%buff_3d(2, i, j, k) =  this% rbuff % genrand_gaussian() * profile(k)
          enddo
        enddo
      enddo

      call ifft( this%buff_3d, (/ nt, this%n_speckle(1), this%n_speckle(2) /) )

      total = 0
      do k = 1, nt
        do j = 1, this%n_speckle(2)
          do i = 1, this%n_speckle(1)
            total = total + this%buff_3d(1, i, j, k)**2 + this%buff_3d(2, i, j, k)**2
          enddo
        enddo
      enddo
      this%buff_3d = this%buff_3d / sqrt(total / &
                      &  real(nt*this%n_speckle(1)*this%n_speckle(2), p_double))

      if ( this%speckle_type == p_cpurp ) then
        do k = 1, nt
          do j = 1, this%n_speckle(2)
            do i = 1, this%n_speckle(1)
              this%buff_3d(1, i, j, k) = cos( this%phasemod_amplitude(1,1) * &
                              &               this%buff_3d(1,i,j,k) )
              this%buff_3d(2, i, j, k) = sin( this%phasemod_amplitude(1,1) * &
                              &               this%buff_3d(2,i,j,k) )
            enddo
          enddo
        enddo
      endif

    endif
  end select

  call freemem( profile )
  call fft_cleanup()
end subroutine new_gaussian_psd_timeseries

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_speckle_pattern_2d( wall_2d, per_mode_2d, phase, ns, lnx )
  implicit none
  real(p_double), dimension(1:, 1:),     intent(inout) :: wall_2d
  real(p_double), dimension(1:, 1:, 1:), intent(in)    :: per_mode_2d
  real(p_double), dimension(1:, 1:),     intent(in)    :: phase
  integer,                               intent(in)    :: ns, lnx


  ! local variables
  integer ind1, hmode, i

  wall_2d(:,:) = 0.0_p_double
  hmode = int( ns / 2 )

  do ind1 = 1, lnx
    do i = 1, ns
      wall_2d(1,ind1) = wall_2d(1,ind1) + per_mode_2d(1,i,ind1) * phase(1, i) - &
                      &                   per_mode_2d(2,i,ind1) * phase(2, i)
      wall_2d(2,ind1) = wall_2d(2,ind1) + per_mode_2d(1,i,ind1) * phase(2, i) + &
                      &                   per_mode_2d(2,i,ind1) * phase(1, i)

      wall_2d(3,ind1) = wall_2d(3,ind1) + per_mode_2d(3,i,ind1) * phase(1, i) - &
                      &                   per_mode_2d(4,i,ind1) * phase(2, i)
      wall_2d(4,ind1) = wall_2d(4,ind1) + per_mode_2d(3,i,ind1) * phase(2, i) + &
                      &                   per_mode_2d(4,i,ind1) * phase(1, i)
    end do
  end do
end subroutine update_speckle_pattern_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_speckle_pattern_3d( wall_3d, per_mode_3d, phase, ns, lnx )
  implicit none
  real(p_double), dimension(1:, 1:, 1:),     intent(inout) :: wall_3d
  real(p_double), dimension(1:, 1:, 1:, 1:), intent(in)    :: per_mode_3d
  real(p_double), dimension(1:, 1:, 1:),     intent(in)    :: phase
  integer, dimension(1:per_ndim),            intent(in)    :: ns, lnx

  ! local variables
  integer :: l, m, ind1, ind2, i

  wall_3d(:,:,:) = 0.0_p_double

  do ind2 = 1, lnx(2)
  do ind1 = 1, lnx(1)
    i = 1

    do m = 1, ns(2)
    do l = 1, ns(1)
      wall_3d(1,ind1,ind2) = wall_3d(1,ind1,ind2) + per_mode_3d(1,i,ind1,ind2) * phase(1, l, m) - &
                      &      per_mode_3d(2,i,ind1,ind2) * phase(2, l, m)
      wall_3d(2,ind1,ind2) = wall_3d(2,ind1,ind2) + per_mode_3d(1,i,ind1,ind2) * phase(2, l, m) + &
                      &      per_mode_3d(2,i,ind1,ind2) * phase(1, l, m)

      wall_3d(3,ind1,ind2) = wall_3d(3,ind1,ind2) + per_mode_3d(3,i,ind1,ind2) * phase(1, l, m) - &
                      &      per_mode_3d(4,i,ind1,ind2) * phase(2, l, m)
      wall_3d(4,ind1,ind2) = wall_3d(4,ind1,ind2) + per_mode_3d(3,i,ind1,ind2) * phase(2, l, m) + &
                      &      per_mode_3d(4,i,ind1,ind2) * phase(1, l, m)

      i = i + 1
    end do
    end do

  end do
  end do
end subroutine update_speckle_pattern_3d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_ssd_buffer_2d( buff_2d, pp_2d, pmp_2d, phasemod_amplitude, pm_bw, nfm, ns, &
                        &        inj_time, beamlet_delay )
  implicit none
  real(p_double), dimension(1:,1:,1:), intent(inout) :: buff_2d
  real(p_double), dimension(1:),       intent(in)    :: pp_2d, pmp_2d, phasemod_amplitude, pm_bw
  integer,                             intent(in)    :: nfm, ns
  real(p_double),                      intent(in)    :: inj_time, beamlet_delay

  integer :: ind1, ind2
  real(p_double) :: ph

  do ind1 = 1, ns
    ph = pp_2d(ind1)
    do ind2 = 1, nfm
      ph = phasemod_amplitude(ind2) * sin(pm_bw(ind2) * (inj_time + &
                      & beamlet_delay * (ind1-1) ) + pmp_2d(ind2) ) + ph
    enddo
    buff_2d(1, ind1, 1) = cos(ph)
    buff_2d(2, ind1, 1) = sin(ph)
  end do
end subroutine update_ssd_buffer_2d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_ssd_buffer_3d( buff_3d, pp_3d, pmp_3d, phasemod_amplitude, pm_bw, nfm, ns, &
                        &        inj_time, beamlet_delay )
  implicit none
  real(p_double), dimension(1:,1:,1:,1:),  intent(inout) :: buff_3d
  real(p_double), dimension(1:,1:),        intent(in)    :: pmp_3d
  real(p_double), dimension(1:,1:),        intent(in)    :: pp_3d, phasemod_amplitude, pm_bw
  integer, dimension(1:per_ndim),          intent(in)    :: nfm, ns
  real(p_double),                          intent(in)    :: inj_time
  real(p_double), dimension(1:per_ndim),   intent(in)    :: beamlet_delay

  integer :: i, j, k, ind1
  real(p_double) :: ph

  do j = 1, ns(2)
    do i = 1, ns(1)
      ph = pp_3d(i,j)
      do k = 1, 2
        do ind1 = 1, nfm(k)
          ph = phasemod_amplitude(ind1,k) * sin(pm_bw(ind1,k) * (inj_time + beamlet_delay(k) * (i - 1) + pmp_3d(ind1, k))) + ph
        enddo
      enddo
      buff_3d(1, i, j, 1) = cos(ph)
      buff_3d(2, i, j, 1) = sin(ph)
    enddo
  end do
end subroutine update_ssd_buffer_3d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine new_modulation_phase_2d( pmp_2d, nfm, prng )
  implicit none
  real(p_double), dimension(1:), intent(inout) :: pmp_2d
  integer,                       intent(in)    :: nfm
  class(t_random), pointer,      intent (in)   :: prng

  integer :: j
  do j = 1, nfm
    pmp_2d(j) = prng%genrand_gaussian() * pi
  enddo
end subroutine new_modulation_phase_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine new_modulation_phase_3d( pmp_3d, nfm, prng )
  implicit none
  real(p_double), dimension(1:,1:), intent(inout) :: pmp_3d
  integer, dimension(1:2),          intent(in)    :: nfm
  class(t_random), pointer,         intent (in)   :: prng

  integer :: k, l
  do l = 1, 2
    do k = 1, nfm(l)
      pmp_3d(k, l) = prng%genrand_gaussian() * pi
    enddo
  enddo
end subroutine new_modulation_phase_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! **temporary**
!early termination if FFTW is not present and Gaussian bandwidth is chosen
subroutine check_for_fft_support()
!-----------------------------------------------------------------------------------------
#ifndef FFTW_ENABLED
  if ( mpi_node() == 0 ) then
    write(0,*)  ""
    write(0,*)  "   Error in zpulse_speckle parameters"
    write(0,*)  "   FFT library not enabled in compilation, "
    write(0,*)  "   'Gaussian' bandwidth_type not supported."
    write(0,*)  "   aborting..."
  endif
  stop
#endif
end subroutine check_for_fft_support
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine read_input_zpulse_speckle( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

  use m_input_file

  implicit none

  class( t_zpulse_speckle ), intent(inout) :: this
  class( t_input_file ),      intent(inout) :: input_file
  type( t_space ),           intent(in)    :: g_space
  class (t_emf_bound),       intent(in)    :: bnd_con
  logical, dimension(1:),    intent(in)    :: periodic
  class( t_grid ),           intent(in)    :: grid
  type( t_options ),         intent(in)    :: sim_options

  ! local variables
  real(p_double) :: a0
  real(p_double) :: omega0
  real(p_double) :: phase
  integer        :: pol_type
  real(p_double) :: pol, pol_stepping
  character(len=16) :: propagation
  integer         :: direction

  integer            :: chirp_order
  real(p_double), dimension(p_max_chirp_order) :: chirp_coefs
  real(p_double) :: chirp_period

  ! Temporal profile
  character(len=16) :: tenv_type

  real(p_double) :: tenv_fwhm

  real(p_double) :: tenv_rise, tenv_flat, tenv_fall
  real(p_double) :: tenv_duration, tenv_range
  character(len = p_max_expr_len) :: tenv_math_func

  real(p_double) :: duty_cycle, stud_jitter

  ! Beam profile
  !character(len=16) :: per_type
  real(p_double), dimension(2) :: spatial_tilt
  real(p_double), dimension(2) :: per_center

  real(p_double), dimension(2) :: per_w0, per_fwhm, per_focus
  real(p_double), dimension(2) :: theta

  integer, dimension(2) :: per_tem_mode

  ! real(p_double), dimension(2,2) :: per_w0_asym, per_fwhm_asym, per_focus_asym
  ! real(p_double), dimension(2) :: per_asym_trans

  integer         :: per_n, per_0clip
  real(p_double) :: per_kt, lon_focus

  real(p_double) :: launch_time
  logical        :: if_launch
  integer :: pos

  ! Speckle pattern
  logical               :: if_dump_speckle
  integer, dimension(per_ndim) :: n_speckle, nfm, color_cycle
  real(p_double), dimension(p_max_nfm, per_ndim) :: fm, phasemod_amplitude
  Character(len=16)  :: speckle_type, bandwidth_type
  real(p_double)     :: recurrence_time, stud_duration
  real(p_double)     :: dt_update_speckle, dt_update_pol
  integer            :: rseed_phase, rseed_stud
  !character(len = p_max_filename_len) :: filename
  real(p_double), dimension(per_ndim) :: laser_bandwidth, beamlet_delay

  namelist /nl_zpulse_speckle/ &
        a0, omega0,phase, pol_type, pol, propagation, direction, &
        speckle_type, laser_bandwidth, recurrence_time, &
        duty_cycle, stud_duration, stud_jitter, n_speckle, &
        bandwidth_type, color_cycle, nfm, fm, &
        phasemod_amplitude, dt_update_speckle, if_dump_speckle, &
        rseed_phase, rseed_stud, dt_update_pol, pol_stepping, &
        chirp_order, chirp_coefs, chirp_period, tenv_fwhm, &
        tenv_type, tenv_rise, tenv_flat, tenv_fall, lon_focus, &
        tenv_duration, tenv_range, tenv_math_func, spatial_tilt, &
        per_center, per_w0, theta, &
        !per_fwhm, &
        !per_w0_asym, per_fwhm_asym, per_focus_asym, per_asym_trans, &
        !per_n, per_0clip, per_kt, per_tem_mode, &
        launch_time, if_launch

  integer :: i, j, ierr, ndim
  integer, dimension(1)  :: maxfmloc
  real(p_double)     :: tot_t, lbw
  integer, parameter :: p_default_seed_phase   = 23 ! random seed
  integer, parameter :: p_default_seed_stud    = 71 ! random seed

  SCR_ROOT(" - Reading zpulse_speckle configuration...")

  ! set default values
  if_launch = .true.

  a0             = 1.0_p_double
  omega0         = 1.0_p_double
  phase          = 0.0_p_double
  pol_type       = 0
  pol            = 0.0_p_double
  pol_stepping   = 90.0_p_double
  dt_update_pol  = - huge(1.0_p_double)

  chirp_order    = 0
  chirp_coefs    = 0.0_p_double
  chirp_period   = huge(1.0_p_double)

  propagation    = "forward"
  direction      = 1

  tenv_type      = "polynomial"

  tenv_fwhm      = -1.0_p_double

  tenv_rise      = 0.0_p_double
  tenv_flat      = 0.0_p_double
  tenv_fall      = 0.0_p_double

  tenv_duration  = 0.0_p_double
  tenv_range     = - huge(1.0_p_double)
  tenv_math_func = "NO_FUNCTION_SUPPLIED!"

  duty_cycle     = 50.0_p_double
  stud_jitter    = 0.0_p_double
  stud_duration  = - huge(1.0_p_double)

  spatial_tilt   = 0.0_p_double

  lon_focus      = 0.0_p_double

  per_center     = - huge(1.0_p_double)

  per_w0         = 0.0_p_double
  per_fwhm       = 0.0_p_double
  per_focus      = - huge(1.0_p_double)
  per_tem_mode   = 0 ! tem mode, integer, default = 0,0

!  per_w0_asym    = 0.0_p_double
!  per_fwhm_asym  = 0.0_p_double
!  per_focus_asym = - huge(1.0_p_double)
!  per_asym_trans = 0.0_p_double

  per_n          = 0
  per_0clip      = 1
  per_kt         = 0.0_p_double

  theta          = 0.0_p_double


  launch_time           = 0.0_p_double

  n_speckle          = 1
  speckle_type       = "default"
  bandwidth_type     = 'default'
  laser_bandwidth    = 0.0_p_double
  phasemod_amplitude = 2 * PI  ! should be > pi
  dt_update_speckle  = - huge(1.0_p_double)
  rseed_phase        = p_default_seed_phase
  rseed_stud         = p_default_seed_stud
  color_cycle        = 1.0_p_double
  nfm                = 1
  fm                 = 1.0_p_double
  !filename     = "NO_FILE_SUPPLIED"
  recurrence_time    = - huge(1.0_p_double)
  if_dump_speckle    = .false.
  beamlet_delay      = - huge(1.0_p_double)

  ! read data from file
  call get_namelist( input_file, "nl_zpulse_speckle", ierr )

  if ( ierr /= 0 ) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        print *, "Error reading zpulse_speckle parameters"
      else
        print *, "Error: zpulse_speckle parameters missing"
      endif
      print *, "aborting..."
    endif
    stop
  endif

  read (input_file%nml_text, nml = nl_zpulse_speckle, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error reading zpulse_speckle parameters"
      write(0,*)  "   aborting..."
    endif
    stop
  endif

  this%a0         = a0
  this%omega0     = omega0
  this%phase0     = real( phase * pi_180, p_double )
  this%lon_focus  = lon_focus

  ! Polarization parameters
  if ((pol_type <-1) .or. (pol_type >1)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_speckle parameters"
      write(0,*)  "   pol must be in the range [-1,0,+1]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%pol_type  = pol_type
  this%pol = real( pol * pi_180, p_double )
  this%pol_stepping = pol_stepping

  ! Frequency chirp parameters
  if ((chirp_order < 0) .or. (chirp_order > p_max_chirp_order)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_speckle parameters"
      write(0,*)  "   chirp_order must be in the range [0,",p_max_chirp_order,"]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  if (chirp_period == 0) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_speckle parameters"
      write(0,*)  "   chirp_period cannot be zero"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%chirp_order     = chirp_order
  this%chirp_coefs     = chirp_coefs
  this%chirp_period    = chirp_period

  ! Propagation direction
  select case ( trim(propagation) )
  case ("forward")
    this%propagation = p_forward
  case ("backward")
    this%propagation = p_backward
  case default
    if ( mpi_node() == 0 ) then
      write(0,*)  ''
      write(0,*)  '   Error in zpulse_speckle parameters'
      write(0,*)  '   propagation must be either "forward" or "backward"'
      write(0,*)  '   aborting...'
    endif
    stop
  end select

  ! Launch direction
  if ((direction < 1) .or. (direction > p_x_dim)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_speckle parameters"
      write(0,*)  "   direction must be in the range [1 .. x_dim]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%direction = direction


  ! spatial tilt
  if ( p_x_dim > 1 ) then
    do i = 1, p_x_dim-1
      if ( spatial_tilt(i) /= 0.0_p_double ) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ""
          write(0,*)  "   Error in zpulse_speckle parameters"
          write(0,*)  "   spatial_tilt has not been implemented yet."
          write(0,*)  "   aborting..."
        endif
        stop
        !if (( spatial_tilt(i) < -45.0_p_double ) .or. (spatial_tilt(i) > 45.0_p_double )) then
        !  if ( mpi_node() == 0 ) then
        !    write(0,*)  ""
        !    write(0,*)  "   Error in zpulse_speckle parameters"
        !    write(0,*)  "   spatial_tilt must be in the range [ -45 , 45 ] deg"
        !    write(0,*)  "   aborting..."
        !  endif
        !  stop
        !endif
        !this%spatial_tilt(i) = tan( real( spatial_tilt(i) * pi_180, p_double ))
        !if ( this%propagation == p_forward ) this%spatial_tilt(i) = -this%spatial_tilt(i)
        !this%if_spatial_tilt = .true.
      else
        this%spatial_tilt(i) = 0.0_p_double
      endif
    enddo
  endif


  ! Temporal envelope parameters
  select case ( to_lower(trim(tenv_type)) )
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
        write(0,*)  '   Error in zpulse_speckle parameters'
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
        write(0,*)  '   Error in zpulse_speckle parameters'
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
        write(0,*)  '   Error in zpulse_speckle parameters'
        write(0,*)  '   Supplied tenv_math_func failed to compile'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

  case default
    if ( mpi_node() == 0 ) then
      write(0,*)  ''
      write(0,*)  '   Error in zpulse_speckle parameters'
      write(0,*)  '   lon_type must be either "polynomial", "gaussian",'
      write(0,*)  '   "sin2" or "math".'
      write(0,*)  '   aborting...'
    endif
    stop
  end select


  do i = 1, p_x_dim - 1
    if ( n_speckle(i) < 1 ) then
      if (mpi_node() == 0) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_speckle parameters'
        write(0,*)  '   n_speckle must be greater than zero'
        write(0,*)  '   aborting...'
      endif
        stop
    endif
  enddo
  this%n_speckle = n_speckle

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

  this%stud_duration = stud_duration
  if ( stud_duration > 0.0 ) then

    if ( duty_cycle < 0.0 ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_speckle parameters'
        write(0,*)  '   duty_cycle must not be smaller than zero'
        write(0,*)  '   aborting...'
      endif
      stop
    endif
    if ( (abs(stud_jitter) + duty_cycle) > 100. ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_speckle parameters'
        write(0,*)  '   (stud_jitter+duty_cycle) must be less than 100 (%)'
      endif
      stop
    endif
    if ( abs(stud_jitter) > duty_cycle ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_speckle parameters'
        write(0,*)  '   stud_jitter must be smaller than duty_cycle'
      endif
      stop
    endif
    this%duty_cycle = duty_cycle / 100.
    this%stud_jitter = abs(stud_jitter) / 100.

    select case ( this%tenv_type )
    case (p_polynomial, p_sin2)
      this%stud_period = (this%tenv_rise + this%tenv_flat + this%tenv_fall) / this%duty_cycle
    case (p_gaussian, p_func)
      this%stud_period = (this%tenv_range) / this%duty_cycle
    case default
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_speckle parameters'
        write(0,*)  '   when using STUD pulses,'
        write(0,*)  '   tenv_type must be one of the following'
        write(0,*)  '   polynomial, sin2, gaussian or math'
      endif
      stop
    end select

    ! change some default values for stud pulses
    if ( to_lower(trim(speckle_type)) == 'default' ) then
      speckle_type = 'cpurp'
    endif
    if ( to_lower(trim(bandwidth_type)) == 'default' ) then
      bandwidth_type = 'lorentzian'
    endif
    do i = 1, per_ndim
      if ( laser_bandwidth(i) < 0.0 ) then
        laser_bandwidth(i) = 1.0
      endif
    enddo
    if ( dt_update_speckle < 0 ) then
      dt_update_speckle = this%stud_period
    endif
  endif

  this%per_center(1:2)     = per_center(1:2)
  this%per_w0(1:2)         = per_w0(1:2)


  this%k_perp(1:2)         = sin(theta(1:2)/360.0*(2.0_p_double*pi))


  ! do not allow speckle zpulse without open boundaries
  if ( this%propagation == p_backward) then
    pos = p_upper
  else
    pos = p_lower
  endif

  if ( periodic(this%direction) .or. (.not. is_open(bnd_con, this%direction, pos)) ) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_speckle parameters"
      write(0,*)  "   Boundary condition must be open ('vpml' or 'lindman') in the speckle launching the pulse."
      write(0,*)  "   aborting..."
    endif
    stop
  endif

  this%launch_time    = launch_time
  this%if_launch = if_launch

  ! 1D is treated as a special 2D case
  if ( p_x_dim == 1 ) then
    n_speckle(:) = 1
  endif

  if ( product( n_speckle(1:per_ndim) ) == 1 ) then
    beamlet_delay = 0.0
  endif

  ! speckle parameters
  ndim = max(1, p_x_dim - 1)

  select case ( to_lower(trim(speckle_type)) )
  case ( 'rpp' )
    this%speckle_type = p_rpp
    this%nt    = 1
    !laser_bandwidth = 0
    this%dt_update_speckle = -1.0

  ! the random phases in both CPP and DPP can change contiguously
  case ( 'cpp', 'dpp', 'default' )
    this%speckle_type = p_cpp
    this%nt    = 1
    !laser_bandwidth = 0
    this%dt_update_speckle = -1.0
  end select

  lbw = max( laser_bandwidth(1), laser_bandwidth(per_ndim) )

  if ( lbw > 0 ) then

    if ( dt_update_speckle <= 0 ) then
      this%dt_update_speckle = 0.01 / lbw
    else
      this%dt_update_speckle = dt_update_speckle
    endif

    do i=1, per_ndim

      if ( nfm(i) < 1 .or. nfm(i) > p_max_nfm) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ""
          write(0,*)  "   Error in zpulse_speckle parameters"
          write(0,*)  "   The values of nfm must be integer between 1 and ", p_max_nfm
          write(0,*)  "   Please contact the development team if you need more frequency bands."
          write(0,*)  "   aborting..."
        endif
        stop
      endif
      this%nfm(i) = nfm(i)

      do j = 1, nfm(i)

        if ( phasemod_amplitude(j, i) > 0 ) then
          if ( mpi_node() == 0 .and. phasemod_amplitude(j, i) < pi ) then
            write(0,*)  ""
            write(0,*)  "   Warning: phasemod_amplitude(", j,',',i, ") smaller than PI"
            if ( fm(j, i) == 1.0 ) then
              write(0,*)  "   laser_bandwidth=", laser_bandwidth, " unachievable."
            endif
          endif
          this%phasemod_amplitude(j, i) = phasemod_amplitude(j, i)
          this%pm_bw(j, i) = 0.5 * laser_bandwidth(i) / phasemod_amplitude(j, i)
        else
          if ( mpi_node() == 0 ) then
            write(0,*)  ""
            write(0,*)  "   Error in zpulse_speckle parameters"
            write(0,*)  "   phasemod_amplitude must be greater than 0"
            write(0,*)  "   aborting..."
          endif
          stop
        endif

      enddo
    enddo
    this%phasemod_amplitude = phasemod_amplitude

    select case ( to_lower(trim(speckle_type)) )
    case ( 'isi', 'rprp' )
      this%speckle_type = p_isi
      ! /2pi for frequency to omega normalization
      laser_bandwidth = laser_bandwidth * 0.5 / pi / sqrt( real(per_ndim, p_k_fld) )

      select case ( to_lower(trim(bandwidth_type)) )
      ! We greatly appreciate the discussions with researchers from the Naval Research Laboratory (NRL)
      ! The Gaussian bandwidth model is closely related to the model for KrF laser at NRL
      case ( 'gaussian', 'default' )
        call check_for_fft_support()
        this%bandwidth_type = p_gaussian_spectrum
        if ( recurrence_time <= 0 ) then
          if ( mpi_node() == 0 ) then
            write(0,*)  ""
            write(0,*)  "   Error in zpulse_speckle parameters"
            write(0,*)  "   recurrence_time must be specified and"
            write(0,*)  "   must be greater than zero when using Gaussian (default) spectrum"
            write(0,*)  "   aborting..."
          endif
          stop
        endif
        ! 2 * how many times updating the speckle pattern, because the amplitude in ISI is a complex number
        this%nt = ceiling( recurrence_time / this%dt_update_speckle )
        this%dt_update_speckle = recurrence_time / this%nt
      case ( 'lorentzian' )
        this%bandwidth_type = p_lorentzian_spectrum
        laser_bandwidth = laser_bandwidth * sqrt(8.0)
        this%nt = 1  ! we use AR(1) process to generate random phases
      case default
        if ( mpi_node() == 0 ) then
          write(0,*)  ""
          write(0,*)  "   Error in zpulse_speckle parameters"
          write(0,*)  "   valid options for ISI speckle_type are:"
          write(0,*)  "   default, gaussian and lorentzian."
          write(0,*)  "   aborting..."
        endif
        stop
      end select

    case ( 'ssd', 'cprp' )
      this%speckle_type = p_ssd
      do i = 1, per_ndim
        if ( color_cycle(i) <= 0. ) then
          if ( mpi_node() == 0 ) then
            write(0,*)  ""
            write(0,*)  "   Error in zpulse_speckle parameters"
            write(0,*)  "   color_cycle must be greater than 0"
            write(0,*)  "   aborting..."
          endif
          stop
        endif
      enddo
      this%color_cycle = color_cycle

      select case ( to_lower(trim(bandwidth_type)) )
      case ( 'line', 'default' )  ! FM modulation
        this%bandwidth_type         = p_line_spectrum
        this%nt                     = 1

        do i=1, per_ndim
          this%pm_bw(1:nfm(i), i)   = this%pm_bw(1:nfm(i), i) * pi
          if ( beamlet_delay(i) >= 0 ) then
            this%beamlet_delay(i) = beamlet_delay(i)
          else
            maxfmloc =  maxloc(fm(1:nfm(i), i))
            if ( this%pm_bw(maxfmloc(1), i) > 0.0 ) then
              this%beamlet_delay(i) = color_cycle(i) * 2*pi / (this%n_speckle(i) * &
                              & this%pm_bw(maxfmloc(1), i))
            else
              this%beamlet_delay(i) = 0.0
            endif
          endif
        end do

      case ( 'gaussian' )
        call check_for_fft_support()
        this%bandwidth_type = p_gaussian_spectrum
        if ( recurrence_time <= 0 ) then
          if ( mpi_node() == 0 ) then
            write(0,*)  ""
            write(0,*)  "   Error in zpulse_speckle parameters"
            write(0,*)  "   recurrence_time must be specified and"
            write(0,*)  "   must be greater than zero when using Gaussian spectrum"
            write(0,*)  "   aborting..."
          endif
          stop
        endif
        tot_t = 0.0
        do i=1, per_ndim
          this%pm_bw(1, i)          = 0.5 * this%pm_bw(1, i)
          if ( beamlet_delay(i) > 0 ) then
            this%beamlet_delay(i)   = beamlet_delay(i)
          else
            if ( this%pm_bw(1, i) > 0.0 ) then
              this%beamlet_delay(i) = color_cycle(i) * 2*pi / (this%n_speckle(i) * this%pm_bw(1, i))
            else
              this%beamlet_delay(i) = 0.0_p_double
            endif
          endif
          tot_t = tot_t + this%beamlet_delay(i) * (this%n_speckle(i) - 1)
        enddo
        this%nt                     = ceiling( (tot_t + recurrence_time) / this%dt_update_speckle )
        this%dt_update_speckle      = (tot_t + recurrence_time) / this%nt

      case ( 'lorentzian' )
        this%bandwidth_type = p_lorentzian_spectrum
        tot_t = 0.0
        do i=1, per_ndim
          this%pm_bw(1, i)       = 0.5 * laser_bandwidth(i) / (phasemod_amplitude(1, i)**2)
          if ( beamlet_delay(i) > 0 ) then
            this%beamlet_delay(i) = beamlet_delay(i)
          else
            if ( this%pm_bw(1, i) > 0.0 ) then
              this%beamlet_delay(i) = color_cycle(i) * 2*pi / (this%n_speckle(i) * this%pm_bw(1, i))
            else
              this%beamlet_delay(i) = 0.0_p_double
            endif
          endif
          tot_t = tot_t + this%beamlet_delay(i) * this%n_speckle(i)
        enddo
        this%nt                  = ceiling( tot_t / this%dt_update_speckle )
        this%dt_update_speckle   = tot_t / this%nt

      case default
        if ( mpi_node() == 0 ) then
          write(0,*)  ""
          write(0,*)  "   Error in zpulse_speckle parameters"
          write(0,*)  "   Valid bandwidth_type in SSD are"
          write(0,*)  "   line, default, Gaussian or Lorentzian"
          write(0,*)  "   aborting..."
        endif
        stop
      end select

    case ( 'cpurp' )
      this%speckle_type = p_cpurp
      select case ( to_lower(trim(bandwidth_type)) )
      case ( 'gaussian', 'default' )
        call check_for_fft_support()
        this%bandwidth_type      = p_gaussian_spectrum
        this%pm_bw(1, 1) = laser_bandwidth(1) / (phasemod_amplitude(1, 1)) / sqrt(real(per_ndim, p_k_fld))
        if ( recurrence_time <= 0 ) then
          if ( mpi_node() == 0 ) then
            write(0,*)  ""
            write(0,*)  "   Error in zpulse_speckle parameters"
            write(0,*)  "   recurrence_time must be specified and"
            write(0,*)  "   must be greater than zero when using CPURP spectrum"
            write(0,*)  "   aborting..."
          endif
          stop
        endif

        ! 2 * how many times updating the speckle pattern, because the amplitude in ISI is a complex number
        this%nt                  = ceiling( recurrence_time / this%dt_update_speckle )
        this%dt_update_speckle   = recurrence_time / this%nt

      case ( 'lorentzian' )
        this%bandwidth_type      = p_lorentzian_spectrum
        this%pm_bw(1, :)         = laser_bandwidth / (phasemod_amplitude(1, 1)**2) / sqrt(real(per_ndim, p_k_fld))
        this%nt                  = 1  ! we use AR(1) process to generate random phases

      case default
        if ( mpi_node() == 0 ) then
          write(0,*)  ""
          write(0,*)  "   Error in zpulse_speckle parameters"
          write(0,*)  "   valid options for CPURP speckle_type are:"
          write(0,*)  "   default, gaussian and lorentzian."
          write(0,*)  "   aborting..."
        endif
        stop
      end select

    !case ( 'file' )
    !  this%speckle_type = p_file
    !  this%filename = filename
    case default
      if ( mpi_node() == 0 ) then
        write(0,*)  ""
        write(0,*)  "   Error in zpulse_speckle parameters"
        write(0,*)  "   valid options for speckle_type (bandwidth > 0) are:"
        write(0,*)  "   isi or rprp (random power random phase),"
        write(0,*)  "   ssd or cprp (constant power random phase), and "
        write(0,*)  "   cpurp (constant power uncorrelated random phase)."
        write(0,*)  "   valid options for speckle_type (bandwidth = 0) are:"
        write(0,*)  "   cpp, rpp or default (not set in inputdeck)"
        write(0,*)  "   aborting..."
      endif
      stop
    end select

  else

    if ( this%speckle_type /= p_rpp .and. this%speckle_type /= p_cpp ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ""
        write(0,*)  "   Error in zpulse_speckle parameters"
        write(0,*)  "   laser_bandwidth must be greater than 0 unless"
        write(0,*)  "   'rpp' or 'cpp' speckle_type is used"
        write(0,*)  "   aborting..."
      endif
      stop
    endif

  endif

  this%laser_bandwidth    = laser_bandwidth
  this%rseed_phase        = rseed_phase
  this%rseed_stud         = rseed_stud
  this%dt_update_pol      = dt_update_pol
  this%pol_stepping       = pol_stepping * real( pi_180, p_k_fld )
  j = 1
  do i = 1, p_x_dim
    if ( i /= this%direction ) then
      this%per_dir(j) = i
      j = j + 1
    endif
  enddo
  do i = 1, p_x_dim-1
    this%pdk(i) = 2*PI / xmax(g_space, this%per_dir(i))
  enddo

end subroutine read_input_zpulse_speckle

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine init_static_modes_12d( this, ti, b, g_x_range, nx_p_min, no_co )
  class( t_zpulse_speckle ),                intent(inout) :: this
  integer,                                  intent(in)    :: ti
  type( t_vdf ) ,                           intent(in)    :: b
  real(p_double), dimension(1:2,1:p_x_dim), intent(in)    :: g_x_range
  class( t_node_conf ),                     intent(in)    :: no_co
  integer,        dimension(1:),            intent(in)    :: nx_p_min

  integer          :: i, hmode, imode1, ind1
  real( p_double ) :: r, r2, ldxp, rmin, ph, fac, fnumber, curvature, curvature2
  real( p_double ) :: r_norm, r2_norm, env, env2

  if ( p_x_dim == 1 ) then
    ldxp = 0.0
    rmin = 0.0
  else
    ldxp = b%dx(this%per_dir(1))
    rmin =  g_x_range(p_lower,this%per_dir(1)) + (nx_p_min(this%per_dir(1))-1)*ldxp
    select case( this%interpolation )
    case( p_linear, p_cubic )
      ! Nothing to change here
    case( p_quadratic, p_quartic )
      rmin = rmin + ldxp*0.5_p_double
    end select
    fnumber = 0.5 * (g_x_range(p_upper,this%per_dir(1)) - g_x_range(p_lower,this%per_dir(1))) * &
             &    this%omega0 / pi / this%n_speckle(1)
  endif


  hmode = int(this%n_speckle(1) / 2)
  fac = 1.0_p_double / sqrt( real(this%n_speckle(1), p_double) )

  do ind1 = 1, this%per_lnx(1)
    i = 1
    r = rmin + (ind1-2) * ldxp
    r2 = r + ldxp/2.0_p_double

    if ( p_x_dim == 1 ) then
      curvature = 0.0
    else
      curvature = this%omega0 * (this%lon_focus + r * this%k_perp(1)) * 0.5 / (fnumber * fnumber)
      curvature2 = this%omega0 * (this%lon_focus + r2 * this%k_perp(1)) * 0.5 / (fnumber * fnumber)
    endif


    if(this%per_w0(1) .ne. 0) then
      r_norm = abs(r - this%per_center(1)) / this%per_w0(1)
      r2_norm = abs(r2 - this%per_center(1)) / this%per_w0(1)
      env = p_env(r_norm)
      env2 = p_env(r2_norm)
    else
      env = 1.
      env2 = 1.
    endif

    do imode1 = hmode - this%n_speckle(1) + 1, hmode
      ph                           = imode1 * ( this%pdk(1)*r + imode1 * curvature )
      this%per_mode_2d(1, i, ind1) = cos(ph) * fac * env
      this%per_mode_2d(2, i, ind1) = sin(ph) * fac * env

      ph                           = imode1 * ( this%pdk(1)*r2 + imode1 * curvature2 )
      this%per_mode_2d(3, i, ind1) = cos(ph) * fac * env2
      this%per_mode_2d(4, i, ind1) = sin(ph) * fac * env2
      i = i + 1
    enddo

  enddo
  call update_speckle_pattern_2d(this%per_wall_2d, this%per_mode_2d, this%buff_2d(:, :, ti), &
                  this%n_speckle(1), this%per_lnx(1))
end subroutine init_static_modes_12d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine init_static_modes_3d( this, ti, b, g_x_range, nx_p_min, no_co )
  class( t_zpulse_speckle ),                intent(inout) :: this
  integer,                                  intent(in)    :: ti
  type( t_vdf ),                            intent(in)    :: b
  real(p_double), dimension(1:2,1:p_x_dim), intent(in)    :: g_x_range
  class( t_node_conf ),                     intent(in)    :: no_co
  integer,        dimension(1:),            intent(in)    :: nx_p_min

  integer  :: i, hmode1, hmode2, imode1, imode2, ind1, ind2
  real( p_double ), dimension(1:2) :: r, r2, ldxp, rmin
  real( p_double ) :: ph, fac, fnumber, curvature, r_norm, r2_norm, k_cs, env, env2

  hmode1    = int(this%n_speckle(1) / 2)
  hmode2    = int(this%n_speckle(2) / 2)
  fnumber   = 0.5 * (g_x_range(p_upper,this%per_dir(1)) - g_x_range(p_lower,this%per_dir(1))) * &
              &    this%omega0 / pi / this%n_speckle(1)

  ldxp(1)   = b%dx(this%per_dir(1))
  ldxp(2)   = b%dx(this%per_dir(2))
  rmin(1)   = g_x_range(p_lower,this%per_dir(1)) + (nx_p_min(this%per_dir(1))-1)*ldxp(1)
  rmin(2)   = g_x_range(p_lower,this%per_dir(2)) + (nx_p_min(this%per_dir(2))-1)*ldxp(2)

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    rmin(1) = rmin(1) + ldxp(1)*0.5_p_double
    rmin(2) = rmin(2) + ldxp(2)*0.5_p_double
  end select

  fac       = 1.0_p_double / sqrt( real(product(this%n_speckle(1:2)), p_double) )

  do ind2 = 1, this%per_lnx(2)
    r(2)  = rmin(2) + (ind2-2) * ldxp(2)
    r2(2) = r(2) + ldxp(2)/2.0_p_double

    do ind1 = 1, this%per_lnx(1)
      r(1)  = rmin(1) + (ind1-2) * ldxp(1)
      r2(1) = r(1) + ldxp(1)/2.0_p_double

      i = 1

      if(this%per_w0(1) .ne. 0 .and. this%per_w0(2) .ne. 0) then
        r_norm = sqrt(((r(1) - this%per_center(1)) / this%per_w0(1))**2 + &
                    & ((r(2) - this%per_center(2)) / this%per_w0(2))**2)
        r2_norm = sqrt(((r2(1) - this%per_center(1)) / this%per_w0(1))**2 + &
                    & ((r2(2) - this%per_center(2)) / this%per_w0(2))**2)
        env = p_env(r_norm)
        env2 = p_env(r2_norm)
      else
        env = 1.
        env2 = 1.
      endif

      k_cs = sqrt(1. - this%k_perp(2)**2) * this%k_perp(1)
      curvature = this%omega0 * (this%lon_focus + k_cs * r(1) + &
                              &  this%k_perp(2) * r(2)) * 0.5 / (fnumber * fnumber)

      do imode2 = hmode2 - this%n_speckle(2) + 1, hmode2

        do imode1 = hmode1 - this%n_speckle(1) + 1, hmode1

          ph = imode1 * ( this%pdk(1)*r(1) + imode1 * curvature ) + &
            &   imode2 * ( this%pdk(2)*r2(2) + imode2 * curvature )
          this%per_mode_3d(1, i, ind1, ind2) = cos(ph) * fac * env
          this%per_mode_3d(2, i, ind1, ind2) = sin(ph) * fac * env

          ph = imode1 * ( this%pdk(1)*r2(1) + imode1 * curvature ) + &
            &   imode2 * ( this%pdk(2)*r(2) + imode2 * curvature )
          this%per_mode_3d(3, i, ind1, ind2) = cos(ph) * fac * env2
          this%per_mode_3d(4, i, ind1, ind2) = sin(ph) * fac * env2

          i = i + 1
        enddo

      enddo

  enddo
  enddo

  call update_speckle_pattern_3d(this%per_wall_3d, this%per_mode_3d, this%buff_3d(:, :, :, ti), &
                  this%n_speckle(1:2), this%per_lnx(1:2))
end subroutine init_static_modes_3d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine init_zpulse_speckle( this, restart, t, dt, b, g_space, nx_p_min, no_co, &
                                interpolation )

  implicit none

  class( t_zpulse_speckle ), intent(inout) :: this
  logical,                   intent(in)    :: restart
  real( p_double ),          intent(in)    :: t, dt
  type( t_vdf ) ,            intent(in)    :: b
  type( t_space ),           intent(in)    :: g_space
  class( t_node_conf ),      intent(in)    :: no_co
  integer, dimension(1:),    intent(in)    :: nx_p_min
  integer,                   intent(in)    :: interpolation

  integer                              :: active_wall, i, j, ind1, ind2, ti, tb
  real(p_double), dimension(2,p_x_dim) :: g_x_range
  real(p_double)                       :: inj_time
  real(p_double), dimension(per_ndim)  :: bw


  inj_time = t - this%launch_time

  if ( this%if_launch ) then
    if (this%propagation == p_forward) then
     active_wall = p_lower
    else
     active_wall = p_upper
    endif

    ti = 1
    tb = 1
    j = 1
    if ( p_x_dim == 1 ) then
      this%per_lnx(1) = 1
    else
      do i = 1, p_x_dim-1
        this%per_lnx(j) = b%nx_( this%per_dir(i) ) + 3
        j = j + 1
      enddo
    endif

    if ( no_co % on_edge(this%direction, active_wall) ) then

      ! not all rng's are necessary but we init them all just make things easier later
      allocate( t_random_hash :: this % rbuff )
      call this % rbuff % init( this % rseed_phase )
      allocate( t_random_hash :: this % rstud )
      call this % rstud % init( this % rseed_stud )

      call get_x_bnd( g_space, g_x_range )
      ! alloc buffer storing the time varying phases
      select case ( p_x_dim )

      case ( 1, 2 )
        if ( this%speckle_type == p_ssd ) then

          call alloc( this%buff_2d, (/ 2, this%n_speckle(1), 1 /) )

          select case( this%bandwidth_type )
          case ( p_gaussian_spectrum, p_lorentzian_spectrum )
            call alloc( this%buff_1d, (/ 2, this%nt /) )
          case ( p_line_spectrum )
            call alloc( this%pmp_2d, (/ this%nfm(1) /) )
          end select

        else

          call alloc( this%buff_2d, (/ 2, this%n_speckle(1), this%nt /) )

        endif

        call alloc( this%pp_2d, (/ this % n_speckle(1) /) )
        call alloc( this%per_wall_2d, (/ 4, this%per_lnx(1) /) )
        call alloc( this%per_mode_2d, (/ 4, this%n_speckle(1), this%per_lnx(1) /) )

      case ( 3 )
        if ( this%speckle_type == p_ssd ) then

          call alloc( this%buff_3d, (/ 2, this%n_speckle(1), this%n_speckle(2), 1 /) )

          select case( this%bandwidth_type )
          case ( p_gaussian_spectrum, p_lorentzian_spectrum )
            call alloc( this%buff_1d, (/ 3, this%nt /) )
          case ( p_line_spectrum )
            call alloc( this%pmp_3d, (/ max(this%nfm(1), this%nfm(2)), 2 /) )
          end select

        else

          call alloc( this%buff_3d, (/ 2, this%n_speckle(1), this%n_speckle(2), this%nt /) )

        endif
        call alloc( this%pp_3d, (/ this % n_speckle(1), this % n_speckle(2) /) )
        call alloc( this%per_wall_3d, (/ 4, this%per_lnx(1), this%per_lnx(2) /) )
        call alloc( this%per_mode_3d, (/ 4, this%n_speckle(1)*this%n_speckle(2), this%per_lnx(1), this%per_lnx(2) /) )
      end select

      select case ( this%speckle_type )

      case ( p_isi, p_cpurp )
        select case ( this%bandwidth_type )

        case ( p_gaussian_spectrum )
          call check_for_fft_support()
          call new_gaussian_psd_timeseries( this )
          if ( restart ) then
            ti = floor(inj_time/this%dt_update_speckle)
            ti = modulo(ti, this%nt) + 1
          endif

        case ( p_lorentzian_spectrum )
          if ( this%speckle_type == p_isi ) then
            bw = this%laser_bandwidth
          else
            bw = this%pm_bw(1, :)
          endif
          this%ar1_rho = exp(- NINT(this%dt_update_speckle / dt) * bw * dt)
          this%ar1_sigma = sqrt(0.5*(1-this%ar1_rho**2))
          call new_phase_plate( this )
          if ( restart .and. this%if_launch ) then
            call fastforward( this, floor(inj_time/this%dt_update_speckle) )
          endif
        end select

      case ( p_ssd )
        ! we would need phase plate for SSD anyway.
        call new_phase_plate( this )

        select case ( this%bandwidth_type )
        case ( p_line_spectrum )
          select case ( p_x_dim )

          case ( 1, 2 )
            call new_modulation_phase_2d( this%pmp_2d, this%nfm(1), this%rbuff )
            call update_ssd_buffer_2d( this%buff_2d, this%pp_2d, this%pmp_2d, this%phasemod_amplitude(:,1), this%pm_bw(:,1), &
                            &          this%nfm(1), this%n_speckle(1), inj_time, this%beamlet_delay(1) )
          case ( 3 )
            call new_modulation_phase_3d( this%pmp_3d, this%nfm(1:2), this%rbuff )
            call update_ssd_buffer_3d( this%buff_3d, this%pp_3d, this%pmp_3d, this%phasemod_amplitude(:,1:2), this%pm_bw(:,1:2), &
                            &          this%nfm(1:2), this%n_speckle(1:2), inj_time, this%beamlet_delay(1:2) )
          end select

        case ( p_gaussian_spectrum )
          call check_for_fft_support()
          call new_gaussian_psd_timeseries( this )
          if ( restart .and. this%if_launch ) then
            select case ( p_x_dim )

            case ( 1, 2 )
              call copy_rpmssd_buffer_2d(this%buff_2d, this%buff_1d, this%nt, this%n_speckle(1), &
                          this%phasemod_amplitude(1,1), floor(inj_time / this%dt_update_speckle), &
                          floor(this%beamlet_delay(1) / this%dt_update_speckle))
            case ( 3 )
              call copy_rpmssd_buffer_3d(this%buff_3d, this%buff_1d, this%nt, this%n_speckle(1:per_ndim), &
                          this%phasemod_amplitude(1,1:per_ndim), floor(inj_time / this%dt_update_speckle), &
                          floor(this%beamlet_delay(1:per_ndim) / this%dt_update_speckle))
            end select
          endif

        case ( p_lorentzian_spectrum )
          this%ar1_rho(1:per_ndim)   = exp(- NINT(this%dt_update_speckle / dt) * this%pm_bw(1, 1:per_ndim) * dt)
          this%ar1_sigma(1:per_ndim) = sqrt((1-this%ar1_rho(1:per_ndim)**2))
          this%buff_1d(1, 1)         = this%rbuff%genrand_gaussian()

          do i = 2, this%nt
            this%buff_1d(1, i) = ar1process( this%buff_1d(1, i-1), this%ar1_rho(1), &
                                              this%ar1_sigma(1), this%rbuff )
          enddo

          if ( p_x_dim == 3 ) then
            this%buff_1d(2, 1) = this%rbuff%genrand_gaussian()
            do i = 2, this%nt
              this%buff_1d(2, i) = ar1process( this%buff_1d(2, i-1), this%ar1_rho(2), &
                                                this%ar1_sigma(2), this%rbuff )
            enddo
          endif

          if ( restart .and. this%if_launch ) then
            tb    = floor(inj_time/this%dt_update_speckle)
            ind1  = floor(this%beamlet_delay(1) / this%dt_update_speckle)*this%n_speckle(1)
            if ( p_x_dim == 3 ) then
              ind2 = floor(this%beamlet_delay(2) / this%dt_update_speckle)*this%n_speckle(2)
              call fastforward( this%buff_1d, this%nt, tb, (/ ind1, ind2 /), &
                              this%ar1_rho(1:2), this%ar1_sigma(1:2), this%rbuff )
            else
              call fastforward( this%buff_1d, this%nt, tb, ind1, &
                              this%ar1_rho(1), this%ar1_sigma(1), this%rbuff )
            endif
          endif

          select case ( p_x_dim )
          case ( 1, 2 )
            call copy_rpmssd_buffer_2d(this%buff_2d, this%buff_1d, this%nt, this%n_speckle(1), &
                        this%phasemod_amplitude(1,1), tb, floor(this%beamlet_delay(1) / this%dt_update_speckle))
          case ( 3 )
            call copy_rpmssd_buffer_3d(this%buff_3d, this%buff_1d, this%nt, this%n_speckle(1:per_ndim), &
                        this%phasemod_amplitude(1,1:per_ndim), tb, floor(this%beamlet_delay(1:per_ndim) / this%dt_update_speckle))
          end select

        end select

      case ( p_cpp, p_rpp )
        call new_phase_plate( this )
      end select


!      select case ( this%per_type )
!      case (p_plane)
!        ! plane wave pulses do not require divergence correction
!        this%no_div_corr = .true.
!
!      case (p_gaussian)
!        ! this needs to be done locally because of chirp
!
!        ! get rayleigh range
!        !this%per_z0(1) = this%omega0 * this%per_w0(1)**2 /2
!
!      case (p_gaussian_asym)
!        ! no initialization required
!
!      end select

      select case ( p_x_dim )
      !case ( 1 )
      case ( 1, 2 )
        call init_static_modes_12d( this, ti, b, g_x_range, nx_p_min, no_co )
      case ( 3 )
        call init_static_modes_3d( this, ti, b, g_x_range, nx_p_min, no_co )
      end select
    endif

    this%stud_rdc = 1.0_p_double

    if ( this%stud_duration > 0. ) then
      ti = max(1, floor(inj_time / this%dt_update_speckle))
      do i = 1, ti
        call update_stud_params( this )
      enddo
    endif
  endif

  this % interpolation = interpolation

end subroutine init_zpulse_speckle
!-----------------------------------------------------------------------------------------

subroutine cleanup_zpulse_speckle( this )

  implicit none

  class( t_zpulse_speckle ), intent(inout) :: this

  ! buff_1d may be used in any p_x_dim when modeling gaussian rpm ssd
  if ( associated( this%buff_1d ) ) then
    call freemem( this%buff_1d )
  endif
  ! no cleanup required
  select case ( p_x_dim )
  !case ( 1 )
  case ( 1, 2 )
    call freemem( this%buff_2d     )
    call freemem( this%pp_2d       )
    call freemem( this%per_wall_2d )
    call freemem( this%per_mode_2d )
    if ( associated( this%pmp_2d ) ) then
      call freemem( this%pmp_2d )
    endif
  case ( 3 )
    call freemem( this%buff_3d     )
    call freemem( this%pp_3d       )
    call freemem( this%per_wall_3d )
    call freemem( this%per_mode_3d )
    if ( associated( this%pmp_3d ) ) then
      call freemem( this%pmp_3d )
    endif
  end select

end subroutine cleanup_zpulse_speckle

subroutine update_stud_params( this )
  implicit none
  class( t_zpulse_speckle ), intent(inout) :: this

  real(p_double) :: r

  r = (this%rstud%genrand_res53() - 0.5_p_double) * 2 * this%stud_jitter
  this%stud_on_t = this%stud_period * (this%duty_cycle + r)
  this%stud_rdc = 1.0_p_double / (this%duty_cycle + r)
end subroutine update_stud_params

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
function t_env_zpulse_speckle( this, t ) result(t_env)

  implicit none

  class( t_zpulse_speckle ), intent(inout) :: this ! intent must be inout because of
                                                ! the eval( f_parser ) function
  real(p_double), intent(in) :: t
  real(p_double) :: t_env

  ! local variables

  real(p_k_fparse), dimension(1) ::t_dbl
  real(p_double) :: te

  ! executable statements
  if ( this%stud_duration > 0.0 ) then
    te = mod( t, this%stud_period )

    if ( te < this%stud_on_t ) then
      if ( .not. this%stud_was_on ) then
        this%stud_was_on = .true.
      endif

    else

      if ( this%stud_was_on ) then
        call update_stud_params( this )
        this%stud_was_on = .false.
      endif
      t_env = 0
      return

    endif
  else
    te = t
  endif


  select case ( this%tenv_type )

  case (p_polynomial)
    t_env = fenv_poly( te, this%tenv_rise, this%tenv_flat, this%tenv_fall ) * this%stud_rdc

  case (p_sin2)
    t_env = fenv_sin2( te, this%tenv_rise, this%tenv_flat, this%tenv_fall ) * this%stud_rdc

  case (p_gaussian)
    t_env = exp( - 2*((te - 0.5*this % tenv_range)/this%tenv_duration)**2 ) * this%stud_rdc


  case (p_func)

    if ( te > this % tenv_range ) then
      t_env = 0
    else
      t_dbl(1) = te
      t_env = eval( this%lon_math_func, t_dbl ) * this%stud_rdc
    endif

  case default
    t_env = 0

  end select

end function t_env_zpulse_speckle
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Total duration of the laser pulse
!-----------------------------------------------------------------------------------------
function t_duration_zpulse_speckle( this ) result( t_duration )

  implicit none

  class( t_zpulse_speckle ), intent(in) :: this
  real(p_double) :: t_duration

  if ( this%stud_duration > 0.0 ) then

    t_duration = this%stud_duration

  else

    select case ( this%tenv_type )
    case (p_polynomial, p_sin2)
      t_duration = this%tenv_rise + this%tenv_flat + this%tenv_fall

    case (p_gaussian, p_func)
      t_duration = this%tenv_range

    case default
      t_duration = 0

    end select

  endif

end function t_duration_zpulse_speckle
!----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_speckle( this, emf, g_space, &
                               nx_p_min, g_nx, t, dt, no_co )

  implicit none

  class( t_zpulse_speckle ), intent(inout) :: this
  class( t_emf ) ,           intent(inout) :: emf

  type( t_space ),           intent(in) :: g_space
  integer, dimension(:),     intent(in) :: nx_p_min, g_nx
  real(p_double),            intent(in) :: t
  real(p_double),            intent(in) :: dt
  class( t_node_conf ),      intent(in) :: no_co

  real(p_double) :: phase0_temp, pol_temp

  if ( this%if_launch .and. ( t > this%launch_time + this % t_duration() ) ) then
    this%if_launch = .false.
  endif

  if ( this%if_launch .and. ( t >= this%launch_time ) ) then

    ! call the routine with the appropriate dimensions

    select case ( emf % b % x_dim_)
    case(1)
      call launch_em_speckle_1d( this, emf%b, t, dt, no_co )

    case(2)
      call launch_em_speckle_2d( this, emf%b, g_space, nx_p_min, t, dt, no_co )

    case(3)
      call launch_em_speckle_3d( this, emf%b, g_space, nx_p_min, t, dt, no_co )

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
        call launch_em_speckle_1d( this, emf%b, t, dt, no_co )

      case(2)
        call launch_em_speckle_2d( this, emf%b, g_space, nx_p_min, t, dt, no_co )

      case(3)
        call launch_em_speckle_3d( this, emf%b, g_space, nx_p_min, t, dt, no_co )

      end select

      ! Restore phase and polarization for first pulse
      this%phase0 = phase0_temp
      this%pol = pol_temp

    endif

  endif

end subroutine launch_zpulse_speckle
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       returns the value of the perpendicular
!       envelope at the requested position
!-----------------------------------------------------------------------------------------
function speckle_envelope( z_center, w0, phase0, wall_2d, x1, propagation, &
                                & chirp_order, chirp_coefs, chirp_period )

  implicit none

  integer, parameter :: rank = 2

  ! x1 is the longitudinal coordinate
  real(p_double),                intent(in) :: z_center, w0, phase0, x1, chirp_period
  real(p_double), dimension(1:), intent(in) :: wall_2d, chirp_coefs
  integer,                       intent(in) :: propagation, chirp_order

  real(p_double) :: speckle_envelope

  ! local variables
  real(p_double) :: k, ph, cosmode, sinmode

  integer :: i

  ! executable statements

  ! get wavenumber
  k = w0

  ! add chirp
  if ( propagation == p_forward ) then
    do i = 1, chirp_order
      k = k + chirp_coefs(i) * ((mod(x1 - z_center, chirp_period))**(i))
    enddo
  else
    do i = 1, chirp_order
      k = k + chirp_coefs(i) * ((mod(z_center - x1, chirp_period))**(i))
    enddo
  endif

  ! get envelope
  ph = k*(x1-z_center) + phase0
  cosmode = cos(ph)
  sinmode = sin(ph)

  speckle_envelope = cosmode * wall_2d(2) + sinmode * wall_2d(1)
end function speckle_envelope
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_and_get_buffer_index_12d( this, ti, dt )
  implicit none
  class( t_zpulse_speckle ), intent(inout) :: this
  integer,                   intent(inout) :: ti
  real( p_double ),          intent(in)    :: dt

  integer :: i, tiu, ind1
  ! update speckle pattern if needed
  if ( this%dt_update_speckle > 0 ) then
  if ( mod(this%inj_time + dt, this%dt_update_speckle) < dt ) then

    select case ( this%speckle_type )
    !case ( p_file )
    !  ti = floor(this%inj_time/this%dt_update_speckle) + 1


    case ( p_ssd )
      select case ( this%bandwidth_type )
      case ( p_line_spectrum )
        ti = 1
        call update_ssd_buffer_2d( this%buff_2d, this%pp_2d, this%pmp_2d, this%phasemod_amplitude(:,1), this%pm_bw(:,1), &
                        &          this%nfm(1), this%n_speckle(1), this%inj_time, this%beamlet_delay(1) )

      case ( p_gaussian_spectrum )
        ti = 1
        call copy_rpmssd_buffer_2d(this%buff_2d, this%buff_1d, this%nt, this%n_speckle(1), &
                      this%phasemod_amplitude(1,1), floor(this%inj_time/this%dt_update_speckle) + 1, &
                      floor(this%beamlet_delay(1) / this%dt_update_speckle))

      case ( p_lorentzian_spectrum )
        ti = floor(this%inj_time/this%dt_update_speckle) + 1
        tiu = floor(this%beamlet_delay(1) / this%dt_update_speckle)
        ind1 = ti + tiu*this%n_speckle(1)
        call copy_rpmssd_buffer_2d(this%buff_2d, this%buff_1d, this%nt, this%n_speckle(1), &
                      this%phasemod_amplitude(1,1), ti, tiu)
        call fastforward( this%buff_1d, this%nt, 1, ind1, this%ar1_rho(1), this%ar1_sigma(1), this%rbuff )
        ti = 1
      end select


    case ( p_isi )
      select case ( this%bandwidth_type )
      case ( p_gaussian_spectrum )
        ti = modulo(floor(this%inj_time/this%dt_update_speckle), this%nt) + 1
      case ( p_lorentzian_spectrum )
        ti = 1
        do i = 1, 2
          call ar1proc( this%buff_2d(i, :, 1), this%ar1_rho(1), &
                          this%ar1_sigma(1), this%n_speckle(1), this%rbuff )
        enddo
      end select


    case ( p_cpurp )
      select case ( this%bandwidth_type )
      case ( p_gaussian_spectrum )
        ti = modulo(floor(this%inj_time/this%dt_update_speckle), this%nt) + 1
      case ( p_lorentzian_spectrum )
        ti = 1
        call ar1proc( this%pp_2d(:), this%ar1_rho(1), &
                        this%ar1_sigma(1), this%n_speckle(1), this%rbuff )
        do i = 1, this%n_speckle(1)
          this%buff_2d(1, i, 1) = cos( this%phasemod_amplitude(1, 1) * this%pp_2d(i) )
          this%buff_2d(2, i, 1) = sin( this%phasemod_amplitude(1, 1) * this%pp_2d(i) )
        enddo
      end select

    end select
    call update_speckle_pattern_2d(this%per_wall_2d, this%per_mode_2d, this%buff_2d(:, :, ti), &
                        this%n_speckle(1), this%per_lnx(1))
  endif
  endif

end subroutine update_and_get_buffer_index_12d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_and_get_buffer_index_3d( this, ti, dt )
  implicit none
  class( t_zpulse_speckle ), intent(inout) :: this
  integer,                   intent(inout) :: ti
  real( p_double ),          intent(in)    :: dt

  integer, dimension(1:per_ndim) :: tiu, ind
  integer                      :: i, j

  ! update speckle pattern if needed
  if ( this%dt_update_speckle > 0 ) then
  if ( mod(this%inj_time + dt, this%dt_update_speckle) < dt ) then

    select case ( this%speckle_type )
    !case ( p_file )
    !  ti = floor(this%inj_time/this%dt_update_speckle) + 1


    case ( p_ssd )
      select case ( this%bandwidth_type )
      case ( p_line_spectrum )
        ti = 1
        call update_ssd_buffer_3d( this%buff_3d, this%pp_3d, this%pmp_3d, this%phasemod_amplitude(:,1:per_ndim), this%pm_bw(:,1:per_ndim), &
                        &          this%nfm(1:per_ndim), this%n_speckle(1:per_ndim), this%inj_time, this%beamlet_delay(1:per_ndim) )

      case ( p_gaussian_spectrum )
        ti = 1
        call copy_rpmssd_buffer_3d(this%buff_3d, this%buff_1d, this%nt, this%n_speckle(1:per_ndim), &
                      this%phasemod_amplitude(1,1:per_ndim), floor(this%inj_time/this%dt_update_speckle) + 1, &
                      floor(this%beamlet_delay(1:per_ndim) / this%dt_update_speckle))

      case ( p_lorentzian_spectrum )
        ti = floor(this%inj_time/this%dt_update_speckle) + 1
        tiu = floor(this%beamlet_delay / this%dt_update_speckle)
        ind = ti + tiu * this%n_speckle(1:per_ndim)
        call copy_rpmssd_buffer_3d(this%buff_3d, this%buff_1d, this%nt, this%n_speckle(1:per_ndim), &
                      this%phasemod_amplitude(1,1:per_ndim), ti, tiu)
        call fastforward( this%buff_1d, this%nt, 1, ind, this%ar1_rho(1:per_ndim), this%ar1_sigma(1:per_ndim), this%rbuff )
        ti = 1
      end select


    case ( p_isi )
      select case ( this%bandwidth_type )
      case ( p_gaussian_spectrum )
        ti = modulo(floor(this%inj_time/this%dt_update_speckle), this%nt) + 1
      case ( p_lorentzian_spectrum )
        ti = 1
        do i = 1, 2
          call ar1proc( this%buff_3d(i, :, :, 1), this%ar1_rho(1), &
                        this%ar1_sigma(1), this%n_speckle(1:2), this%rbuff )
        enddo
      end select


    case ( p_cpurp )
      select case ( this%bandwidth_type )
      case ( p_gaussian_spectrum )
        ti = modulo(floor(this%inj_time/this%dt_update_speckle), this%nt) + 1
      case ( p_lorentzian_spectrum )
        ti = 1
        call ar1proc( this%pp_3d, this%ar1_rho(1), &
                        this%ar1_sigma(1), this%n_speckle(1:2), this%rbuff )
        do i = 1, this%n_speckle(1)
          do j = 1, this%n_speckle(2)
            this%buff_3d(1, i, j, 1) = cos( this%phasemod_amplitude(1, 1) * this%pp_3d(i, j) )
            this%buff_3d(2, i, j, 1) = sin( this%phasemod_amplitude(1, 1) * this%pp_3d(i, j) )
          enddo
        enddo
      end select

    end select
    call update_speckle_pattern_3d(this%per_wall_3d, this%per_mode_3d, this%buff_3d(:, :, :, ti), &
                        this%n_speckle(1:per_ndim), this%per_lnx(1:per_ndim))
  endif
  endif

end subroutine update_and_get_buffer_index_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine launch_em_speckle_1d( this, b, t, dt, no_co )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_zpulse_speckle ), intent(inout) :: this
  type( t_vdf ) ,            intent(inout)  :: b           ! magnetic field

  real( p_double ),     intent(in) :: t           ! simulation time
  real( p_double ),     intent(in) :: dt         ! time step
  class( t_node_conf ), intent(in) :: no_co       ! node configuration

  ! local variables
  integer :: inpos, active_wall, ti
  real( p_double ) :: prop_sign
  real( p_double ) :: amp, z_center
  real(p_double)  :: ldx, rdtdx
  real(p_double) :: cos_pol, sin_pol

  ! executable statements

  if (this%propagation == p_forward) then
   active_wall = p_lower
   inpos = 1
   prop_sign = 1.0_p_double

  else
   active_wall = p_upper
   inpos = b %nx_( this%direction)
   prop_sign = -1.0_p_double

  endif
   this % inj_time = t - this%launch_time

  if ( no_co % on_edge(this%direction, active_wall) ) then

  ldx = b%dx(this%direction)
  rdtdx = real( dt/ldx, p_k_fld )

  ! update polarization if needed
  if ( this%dt_update_pol > 0 ) then
    if ( mod(this%inj_time + dt, this%dt_update_pol) < dt ) then
      this%pol = mod(this%pol+this%pol_stepping, 2*PI)
    endif
  endif

  call update_and_get_buffer_index_12d( this, ti, dt )

  z_center = this%lon_center()

  ldx = b%dx(this%direction)
  rdtdx = real( dt/ldx, p_k_fld )

  cos_pol = + cos( this%pol ) * prop_sign
  sin_pol = - sin( this%pol ) * prop_sign

  amp = 2.0_p_k_fld * this%omega0 * this%a0 * this % t_envelope( this % inj_time ) * &
             &  rdtdx * speckle_envelope( z_center, this%omega0, this%phase0, this%per_wall_2d(1:2, 1), &
                    & 0.0_p_double, this%propagation, this%chirp_order, this%chirp_coefs, this%chirp_period )


  b%f1(2, inpos) = b%f1(2, inpos) + real( amp * sin_pol, p_k_fld )
  b%f1(3, inpos) = b%f1(3, inpos) + real( amp * cos_pol, p_k_fld )

  endif

end subroutine launch_em_speckle_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! The trick to using the standard per_envelope routines is to use the
! lon_start /
!-----------------------------------------------------------------------------------------
subroutine launch_em_speckle_2d( this, b, g_space, nx_p_min, t, dt, no_co )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  class( t_zpulse_speckle ), intent(inout) :: this
  type( t_vdf ) ,            intent(inout) :: b           ! magnetic field

  type( t_space ),        intent(in) :: g_space     ! global space information
  integer, dimension(1:), intent(in) :: nx_p_min
  real( p_double ),       intent(in) :: t           ! simulation time
  real( p_double ),       intent(in) :: dt         ! time step
  class( t_node_conf ),   intent(in) :: no_co       ! node configuration

! local variables
  integer :: inpos, ip, iperp, ti
  real( p_double ) :: prop_sign, z_center
  real( p_double ) :: amp, lenv, z
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double)  :: ldx, rdtdx, b12, b3
  real(p_double) :: cos_pol, sin_pol
  integer :: active_wall

  real(p_double) :: r_local, r2_local
  real(p_double) :: rmin, k_paral


  call get_x_bnd( g_space, g_x_range )

  if (this%propagation == p_forward) then
   active_wall = p_lower
   inpos = 1
   prop_sign = 1.0_p_double

  else
   active_wall = p_upper
   inpos = b %nx_( this%direction)
   prop_sign = -1.0_p_double

  endif
   z =  g_x_range(active_wall,this%direction)
   this % inj_time = t - this%launch_time

  if ( no_co % on_edge(this%direction, active_wall) ) then

  ldx = b%dx(this%direction)
  rdtdx = real( dt/ldx, p_k_fld )


  ! perpendicular direction
  iperp = 3 - this%direction

  ! update polarization if needed
  if ( this%dt_update_pol > 0 ) then
    if ( mod(this%inj_time + dt, this%dt_update_pol) < dt ) then
      this%pol = mod(this%pol+this%pol_stepping, 2*PI)
    endif
  endif

  rmin = g_x_range(p_lower, iperp) + real(nx_p_min(iperp)-1, p_double) * b%dx(iperp)

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    z = z + ldx*0.5_p_double
    rmin = rmin + b%dx(iperp)*0.5_p_double
  end select

  call update_and_get_buffer_index_12d( this, ti, dt )

  cos_pol = + cos( this%pol ) * prop_sign
  sin_pol = - sin( this%pol ) * prop_sign
  amp = 2.0_p_k_fld * this%omega0 * this%a0

  ! inject column of cells
  lenv = amp * this % t_envelope( this % inj_time   )

  z_center = this%lon_center()
  k_paral = sqrt(1. - this%k_perp(1)**2)

  do ip = 0, this%per_lnx(1) - 1
    r_local = (rmin + real(ip-1, p_double) * real(b%dx(iperp), p_double)) * prop_sign
    r2_local = (r_local + real(b%dx(iperp), p_double) * 0.5) * prop_sign
    b12 = lenv   * rdtdx * speckle_envelope( z_center, this%omega0, this%phase0, this%per_wall_2d(1:2, ip+1), &
                    & z+r_local*this%k_perp(1), this%propagation, this%chirp_order, this%chirp_coefs, this%chirp_period ) * sin_pol
    b3  = lenv   * rdtdx * speckle_envelope( z_center, this%omega0, this%phase0, this%per_wall_2d(3:4, ip+1), &
                    & z+r2_local*this%k_perp(1), this%propagation, this%chirp_order, this%chirp_coefs, this%chirp_period ) * cos_pol

    if (this%direction == 1) then
      ! launching along x1
        b%f2(1, inpos, ip) = b%f2(1, inpos, ip) - real( b12, p_k_fld ) * this%k_perp(1)
        b%f2(2, inpos, ip) = b%f2(2, inpos, ip) + real( b12, p_k_fld ) * k_paral
        b%f2(3, inpos, ip) = b%f2(3, inpos, ip) + real( b3, p_k_fld )
    else
      ! launching along x2
        b%f2(1, ip, inpos) = b%f2(1, ip, inpos) + real( b12, p_k_fld ) * k_paral
        b%f2(2, inpos, ip) = b%f2(2, inpos, ip) + real( b12, p_k_fld ) * this%k_perp(1)
        b%f2(3, ip, inpos) = b%f2(3, ip, inpos) + real( b3, p_k_fld )
    endif

  enddo

  endif

end subroutine launch_em_speckle_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine launch_em_speckle_3d( this, b, g_space, nx_p_min, t, dt, no_co )

  implicit none

  integer, parameter :: rank = 3

  class( t_zpulse_speckle ), intent(inout) :: this
  type( t_vdf ) ,            intent(inout) :: b           ! magnetic field

  type( t_space ),        intent(in) :: g_space     ! global space information
  integer, dimension(1:), intent(in) :: nx_p_min
  real( p_double ),       intent(in) :: t           ! simulation time
  real( p_double ),       intent(in) :: dt         ! time step
  class( t_node_conf ),   intent(in) :: no_co       ! node configuration

  ! local variables
  integer :: ti, inj_node, inpos
  integer :: ip1, ip2
  real( p_double ) :: z, prop_sign
  real( p_double ) :: amp, lenv!, rmin1, rmin2, r1, r2, r1_2, r2_2, dr1_2, dr2_2
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(2) :: rmin, c, r_local
  real(p_double)  :: ldx, rdtdx, b0, b90, z_center
  real(p_double) :: cos_pol, sin_pol, cs

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

    z =  - g_x_range(p_upper,this%direction)
  endif

  if ( my_ngp( no_co, this%direction) == inj_node ) then

    ldx = b%dx(this%direction)
    rdtdx = real( dt/ldx, p_k_fld )

    ! update polarization if needed
    if ( this%dt_update_pol > 0 ) then
      if ( mod(this%inj_time + dt, this%dt_update_pol) < dt ) then
        this%pol = mod(this%pol+this%pol_stepping, 2*PI)
      endif
    endif

    call update_and_get_buffer_index_3d( this, ti, dt )

    cos_pol =  + cos( this%pol ) * prop_sign
    sin_pol =  - sin( this%pol ) * prop_sign
    amp = 2.0_p_k_fld * this%omega0 * this%a0

    ! inject column of cells
    lenv = amp * this % t_envelope( this % inj_time )
    z_center = this%lon_center()
    rmin(1) = g_x_range(p_lower, this%per_dir(1)) + real(nx_p_min(this%per_dir(1))-1, p_double) * b%dx(this%per_dir(1))
    rmin(2) = g_x_range(p_lower, this%per_dir(2)) + real(nx_p_min(this%per_dir(2))-1, p_double) * b%dx(this%per_dir(2))
    c(1:2) = sqrt(1. - this%k_perp(1:2)**2)
    cs = sqrt(1. - this%k_perp(2)**2) * this%k_perp(1)

    select case( this%interpolation )
    case( p_linear, p_cubic )
      ! Nothing to change here
    case( p_quadratic, p_quartic )
      z = z + ldx*0.5_p_double
      rmin(1) = rmin(1) + b%dx(this%per_dir(1))*0.5_p_double
      rmin(2) = rmin(2) + b%dx(this%per_dir(2))*0.5_p_double
    end select

    do ip2 = 0, this%per_lnx(2) - 1
      r_local(1) = (rmin(2) + real(ip2, p_double) * real(b%dx(this%per_dir(2)), p_double)) * this%k_perp(2)
      do ip1 = 0, this%per_lnx(1) - 1

        r_local(2) = z + r_local(1) + (rmin(1) + real(ip1, p_double) * real(b%dx(this%per_dir(1)), p_double)) * cs
        b0  = lenv   * rdtdx * speckle_envelope( z_center, this%omega0, this%phase0, this%per_wall_3d(1:2, ip1+1, ip2+1), &
                        & r_local(2), this%propagation, this%chirp_order, this%chirp_coefs, this%chirp_period ) * sin_pol
        b90 = lenv   * rdtdx * speckle_envelope( z_center, this%omega0, this%phase0, this%per_wall_3d(3:4, ip1+1, ip2+1), &
                        & r_local(2), this%propagation, this%chirp_order, this%chirp_coefs, this%chirp_period ) * cos_pol

        select case (this%direction)

        case (1)
          b%f3(1, inpos, ip1, ip2) = b%f3(1, inpos, ip1, ip2) - real( b0, p_k_fld ) * this%k_perp(1) - real( b90, p_k_fld ) * this%k_perp(2) * c(1)
          b%f3(2, inpos, ip1, ip2) = b%f3(2, inpos, ip1, ip2) + real( b0, p_k_fld ) * c(1) - real( b90, p_k_fld ) * this%k_perp(2) * this%k_perp(1)
          b%f3(3, inpos, ip1, ip2) = b%f3(3, inpos, ip1, ip2) + real( b90, p_k_fld ) * c(2)

        case (2)
          b%f3(1, ip1, inpos, ip2) = b%f3(1, ip1, inpos, ip2) + real( b0, p_k_fld ) * c(1) - real( b90, p_k_fld ) * this%k_perp(2) * this%k_perp(1)
          b%f3(2, ip1, inpos, ip2) = b%f3(2, ip1, inpos, ip2) + real( b0, p_k_fld ) * this%k_perp(1) - real( b90, p_k_fld ) * this%k_perp(2) * c(1)
          b%f3(3, ip1, inpos, ip2) = b%f3(3, ip1, inpos, ip2) + real( b90, p_k_fld ) * c(2)

        case (3)
          b%f3(1, ip1, ip2, inpos) = b%f3(1, ip1, ip2, inpos) + real( b0, p_k_fld ) * c(1) - real( b90, p_k_fld ) * this%k_perp(2) * this%k_perp(1)
          b%f3(2, ip1, ip2, inpos) = b%f3(2, ip1, ip2, inpos) + real( b90, p_k_fld ) * c(2)
          b%f3(3, ip1, ip2, inpos) = b%f3(3, ip1, ip2, inpos) + real( b0, p_k_fld ) * this%k_perp(1) - real( b90, p_k_fld ) * this%k_perp(2) * c(1)

        end select

      enddo
    enddo

  endif

end subroutine launch_em_speckle_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(p_double) function p_env(rad)
real (p_double) rad

! local variable
real (p_double) ab_rad

ab_rad = abs(rad)
if(ab_rad.lt. 1.0) then
    p_env = (1.0-2.5*ab_rad**2+2.5*ab_rad**4-ab_rad**5)
else
    p_env = 0.0
endif

end function p_env

end module m_zpulse_speckle

