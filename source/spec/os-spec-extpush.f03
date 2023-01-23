#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_extpush

#include "memory/memory.h"

use m_parameters

use m_species_define
use m_species_current

use m_vdf_define

use m_emf_define
use m_emf
use m_emf_interpolate

use m_logprof

implicit none

private

interface dudt_exact
  module procedure dudt_exact
end interface

interface dudt_exact_rr
  module procedure dudt_exact_rr
end interface

public :: dudt_exact, dudt_exact_rr

integer, parameter :: p_iter_max = 10

contains

!-----------------------------------------------------------------------------------------
subroutine dudt_exact( this, emf, dt, i0, i1, energy, time )
!-----------------------------------------------------------------------------------------
! Advance velocities and calculate time-centered total energy
! process particles in the range [i0,i1]
! This pusher follows the ref. https://arxiv.org/abs/2007.07556
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  real(p_double), intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_double), intent(inout) :: energy
  real(p_double), intent(in) :: time

  integer :: i, j, ptrcur, np, pp
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep
  real(p_k_part), dimension(4) :: u0, u, fu0, f2u0, f3u0, uk0, fuk0, uk, uo0, fuo0, uo
  real(p_k_part) :: t, cn_o, snc_o, cnc_o, ch_k, shc_k, chc_k, tmp, k2o2, c1, c2
  real(p_k_part) :: inv1, inv2, e2, b2, fac_rr, qm
  real(p_double) :: loc_ene
  real(p_k_part), dimension(p_cache_size) :: k, o, k2, o2
  integer, dimension(p_cache_size) :: ix1, ix2
  integer :: n1, n2

#ifdef __HAS_SPIN__
  real(p_k_part) :: k2w2, o2w2
  real(p_k_part), dimension(4) :: s, sk0, so0, fsk0, fso0, fu, fuk, fuo, fs0, f2s0, f3s0
  real(p_k_part) :: w2, w, hw, ho, hk, at, a, f0, ia_dfdt, fk0, fo0, inv3
  real(p_k_part) :: ch_ak, shc_ak, chc_ak, cn_ao, snc_ao, cnc_ao, cn_w, snc_w, cnc_w
  real(p_k_part) :: ch_1ak, cn_1ao, shc_1ak, snc_1ao
  real(p_k_part) :: c_k, c_o, d_k, d_o
#endif

  real(p_double), save :: rel_tol = 1.0d-12
  real(p_double), save :: weak_fld_tol = 1.0d-12

  ! executable statements
  qm = 1.0_p_k_part / this%rqm
  fac_rr = this%k_rr * this%q_real * dt * 0.5_p_k_part * qm
#ifdef __HAS_SPIN__
  a = this%anom_mag_moment
#endif

  !loop through all particles
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    call this % get_emf( emf, bp, ep, np, ptrcur, t=time )

    loc_ene = 0.0_p_double

    pp = ptrcur
    do i = 1, np
      ep(1,i) = ep(1,i) * qm
      ep(2,i) = ep(2,i) * qm
      ep(3,i) = ep(3,i) * qm
      bp(1,i) = bp(1,i) * qm
      bp(2,i) = bp(2,i) * qm
      bp(3,i) = bp(3,i) * qm
      pp = pp + 1
    enddo

    ! -------------------------------------------------------------------------
    ! half-dt push of radiation reaction
    ! -------------------------------------------------------------------------
    if ( this%rad_react ) then
      pp = ptrcur
      do i = 1, np

        u0(1:3) = this%p(1:3,pp)
        u0(4) = sqrt( 1.0_p_k_part + u0(1)*u0(1) + u0(2)*u0(2) + u0(3)*u0(3) )

        fu0  = get_fu( ep(:,i), bp(:,i), u0 )
        f2u0 = get_fu( ep(:,i), bp(:,i), fu0 )
        tmp  = dot4( u0, f2u0 )

        u0(4) = fac_rr / u0(4)
        this%p(1,pp) = u0(1) + u0(4) * ( f2u0(1) - tmp * u0(1) )
        this%p(2,pp) = u0(2) + u0(4) * ( f2u0(2) - tmp * u0(2) )
        this%p(3,pp) = u0(3) + u0(4) * ( f2u0(3) - tmp * u0(3) )

        pp = pp + 1
      enddo
    endif

    ! -------------------------------------------------------------------------
    ! categorize particles
    ! -------------------------------------------------------------------------
    n1 = 0; ix1 = 0 
    n2 = 0; ix2 = 0
    do i = 1, np

      e2 = ep(1,i) * ep(1,i) + ep(2,i) * ep(2,i) + ep(3,i) * ep(3,i)
      b2 = bp(1,i) * bp(1,i) + bp(2,i) * bp(2,i) + bp(3,i) * bp(3,i)      

      if ( e2 + b2 < weak_fld_tol ) then

        k2(i) = 0.0_p_k_part; k(i) = 0.0_p_k_part
        o2(i) = 0.0_p_k_part; o(i) = 0.0_p_k_part

      else

        ! Lorentz invariants
        inv1 = e2 - b2
        inv2 = ep(1,i) * bp(1,i) + ep(2,i) * bp(2,i) + ep(3,i) * bp(3,i)

        ! calculate eigenvalues
        tmp = sqrt( inv1 * inv1 + 4.0_p_k_part * inv2 * inv2 )
        k2(i) = 0.5_p_k_part * ( inv1 + tmp ); k(i) = sqrt(k2(i))
        o2(i) = 0.5_p_k_part * (-inv1 + tmp ); o(i) = sqrt(o2(i))

        if ( k2(i) < rel_tol .and. o2(i) < rel_tol ) then
          ! light-like
          n1 = n1 + 1; ix1(n1) = i
        else
          ! general case
          n2 = n2 + 1; ix2(n2) = i
        endif

      endif
    enddo

    ! -------------------------------------------------------------------------
    ! advance particles in light-like field (category 1)
    ! -------------------------------------------------------------------------
    do j = 1, n1

      i = ix1(j); pp = ptrcur + i - 1

      u0(1:3) = this%p(:,pp)
      u0(4) = sqrt( 1.0_p_k_part + u0(1)*u0(1) + u0(2)*u0(2) + u0(3)*u0(3) )

      fu0  = get_fu( ep(:,i), bp(:,i), u0 )
      f2u0 = get_fu( ep(:,i), bp(:,i), fu0 )
      f3u0 = get_fu( ep(:,i), bp(:,i), f2u0 )

      ! solve for the proper time step according to dt
      call get_t_light( dt, u0(4), fu0(4), f2u0(4), f3u0(4), c1, c2, t, this%iter_tol )
      u = u0 + t * fu0 + c1 * f2u0 + c2 * f3u0

      this%p(:,pp) = u(1:3)

#ifdef __HAS_SPIN__
      ! get spin in lab frame
      s(1:3) = this%s(:,pp)
      s(4) = ( s(1) * u0(1) + s(2) * u0(2) + s(3) * u0(3) ) / u0(4)

      fs0  = get_fu( ep(:,i), bp(:,i), s )
      f2s0 = get_fu( ep(:,i), bp(:,i), fs0 )
      f3s0 = get_fu( ep(:,i), bp(:,i), f2s0 )
      fu   = get_fu( ep(:,i), bp(:,i), u )

      f0 = dot4( u0, fs0 )
      ia_dfdt = dot4( u0, f2s0 )
      at = a * t
      w  = at + t ! (1+a)t
      w2 = w * w

      ! coefficient C and D
      c_k = -at * ( f0 + 0.5_p_double * ia_dfdt * at )
      d_k = -at * at * ( 0.5_p_double * f0 + 0.166666666666667_p_double * ia_dfdt * at )

      ! homogeneous solution
      this%s(:,pp) = s(1:3) + fs0(1:3) * w + 0.5_p_double * f2s0(1:3) * w2 + 0.166666666666667_p_double * f3s0(1:3) * w2 * w

      ! inhomogeneous solution
      this%s(:,pp) = this%s(:,pp) + c_k * u(1:3) + d_k * fu(1:3)
#endif

    enddo

    ! -------------------------------------------------------------------------
    ! advance particles in arbitrary field (category 2)
    ! -------------------------------------------------------------------------
    do j = 1, n2

      i = ix2(j); pp = ptrcur + i - 1

      u0(1:3) = this%p(:,pp)
      u0(4) = sqrt( 1.0_p_k_part + u0(1)*u0(1) + u0(2)*u0(2) + u0(3)*u0(3) )

      ! subspace decomposition
      k2o2 = k2(i) + o2(i)
      fu0  = get_fu( ep(:,i), bp(:,i), u0 )
      f2u0 = get_fu( ep(:,i), bp(:,i), fu0 )
      uk0  = ( f2u0 + o2(i) * u0 ) / k2o2
      uo0  = u0 - uk0

      fuk0 = get_fu( ep(:,i), bp(:,i), uk0 )
      fuo0 = fu0 - fuk0

      ! solve for the proper time step according to dt
      call get_t( dt, k(i), uk0(4), fuk0(4), o(i), uo0(4), fuo0(4), &
        ch_k, shc_k, chc_k, cn_o, snc_o, cnc_o, t, this%iter_tol )

      uk = uk0 * ch_k + fuk0 * shc_k * t
      uo = uo0 * cn_o + fuo0 * snc_o * t
      u = uk + uo

      this%p(1:3,pp) = u(1:3)

#ifdef __HAS_SPIN__

      ! get spin in lab frame
      s(1:3) = this%s(:,pp)
      s(4) = ( s(1) * u0(1) + s(2) * u0(2) + s(3) * u0(3) ) / u0(4)

      ! subspace decomposition
      fs0  = get_fu( ep(:,i), bp(:,i), s )
      f2s0 = get_fu( ep(:,i), bp(:,i), fs0 )
      sk0  = ( f2s0 + o2(i) * s ) / k2o2
      so0  = s - sk0

      fsk0 = get_fu( ep(:,i), bp(:,i), sk0 )
      fso0 = fs0 - fsk0 

      ! calculate F*u_kappa and F*u_omega
      fuk = get_fu( ep(:,i), bp(:,i), uk )
      fuo = get_fu( ep(:,i), bp(:,i), uo )

      w2 = o2(i) - k2(i) - dot4( fu0, fu0 )
      w  = sqrt( w2 )
      k2w2 = k2(i) + w2
      o2w2 = o2(i) - w2

      fk0 = dot4( uk0, fsk0 )
      fo0 = dot4( uo0, fso0 )
      inv3 = o2(i) * fk0 - k2(i) * fo0

      f0 = fk0 + fo0
      hw = inv3 - f0 * w2
      hk = inv3 + f0 * k2(i)
      ho = inv3 - f0 * o2(i)
      ia_dfdt = k2o2 * dot4( uk0, sk0 )
      at = a * t

      call hyper_func( k(i) * at, ch_ak, shc_ak, chc_ak )
      call trig_func(  o(i) * at, cn_ao, snc_ao, cnc_ao )
      call trig_func( w * at, cn_w, snc_w, cnc_w )

      ! calculate homogeneous solution
      ch_1ak = ch_k * ch_ak + shc_k * shc_ak * at * t * k2(i)
      cn_1ao = cn_o * cn_ao - snc_o * snc_ao * at * t * o2(i)
      shc_1ak = ( ch_k * shc_ak * a + shc_k * ch_ak ) / ( 1 + a )
      snc_1ao = ( cn_o * snc_ao * a + snc_o * cn_ao ) / ( 1 + a )

      this%s(:,pp) = sk0(1:3) * ch_1ak + fsk0(1:3) * shc_1ak * (t+at) + &
                     so0(1:3) * cn_1ao + fso0(1:3) * snc_1ao * (t+at)

      ! calculate the inhomogeneous solution
      c_k = ia_dfdt * ( cn_w - ch_ak )   + ( hw * snc_w - hk * shc_ak ) * at
      d_k = ia_dfdt * ( snc_w - shc_ak ) * at - inv3 * at**2 * ( cnc_w + chc_ak ) + f0 * ( cn_w - ch_ak )
      c_k = c_k / k2w2
      d_k = d_k / k2w2

      ! TODO: the singularity needs to be treated
      c_o = ia_dfdt * ( cn_ao - cn_w )   + ( ho * snc_ao - hw * snc_w ) * at
      d_o = ia_dfdt * ( snc_ao - snc_w ) * at - inv3 * at**2 * ( cnc_ao - cnc_w ) + f0 * ( cn_ao - cn_w )
      c_o = c_o / o2w2
      d_o = d_o / o2w2

      this%s(:,pp) = this%s(:,pp) + c_k * uk(1:3) + d_k * fuk(1:3) + &
                                    c_o * uo(1:3) + d_o * fuo(1:3)
#endif

    enddo

    ! -------------------------------------------------------------------------
    ! half-dt push of radiation reaction
    ! -------------------------------------------------------------------------
    if ( this%rad_react ) then
      pp = ptrcur
      do i = 1, np
        u0(1:3) = this%p(1:3,pp)
        u0(4) = sqrt( 1.0_p_k_part + u0(1) * u0(1) + u0(2) * u0(2) + u0(3) * u0(3) )

        fu0  = get_fu( ep(:,i), bp(:,i), u0 )
        f2u0 = get_fu( ep(:,i), bp(:,i), fu0 )
        tmp  = dot4( u0, f2u0 )

        u0(4) = fac_rr / u0(4)
        this%p(1,pp) = this%p(1,pp) + u0(4) * ( f2u0(1) - tmp * u0(1) )
        this%p(2,pp) = this%p(2,pp) + u0(4) * ( f2u0(2) - tmp * u0(2) )
        this%p(3,pp) = this%p(3,pp) + u0(4) * ( f2u0(3) - tmp * u0(3) )

        pp = pp + 1
      enddo
    endif

    pp = ptrcur
    do i = 1, np
      tmp = sqrt(1.0 + this%p(1,pp)*this%p(1,pp) + this%p(2,pp)*this%p(2,pp) + this%p(3,pp)*this%p(3,pp))
      loc_ene = loc_ene + this%q(pp) * (tmp - 1.0)
      pp = pp + 1
    enddo

    energy = energy + loc_ene

  enddo

end subroutine dudt_exact
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_t_light( dt, u0, fu0, f2u0, f3u0, c1, c2, t, tol )
!-----------------------------------------------------------------------------------------
! Using Newton-Raphson method to calculate proper time step for the light-like case
!-----------------------------------------------------------------------------------------

  implicit none

  real(p_double), intent(in) :: dt, tol
  real(p_k_part), intent(in) :: u0, fu0, f2u0, f3u0
  real(p_k_part), intent(out) :: t, c1, c2

  real(p_k_part) :: t1, t2, t3, fun, dfun
  integer :: iter

  t1 = dt / u0
  iter = 0

  do

    t2 = t1 * t1
    t3 = t2 * t1
    c1 = 0.5_p_double * t2
    c2 = 0.166666666666667_p_double * t3

    fun  = t1 * u0 + c1 * fu0 + c2 * f2u0 + 0.041666666666667_p_double * t3 * t1 * f3u0 - dt
    dfun = u0 + t1 * fu0 + c1 * f2u0 + c2 * f3u0

    t = t1 - fun / dfun
    iter = iter + 1

    if ( abs(t - t1) / t1 < tol .or. iter > p_iter_max ) then
      exit
    else
      t1 = t
    endif

  enddo

end subroutine get_t_light
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_t( dt, k, uk0, fuk0, o, uo0, fuo0, ch_k, shc_k, chc_k, cn_o, snc_o, cnc_o, t, tol )
!-----------------------------------------------------------------------------------------
! Using Newton-Raphson method to calculate proper time step for general case
!-----------------------------------------------------------------------------------------

  implicit none

  real(p_double), intent(in) :: dt, tol
  real(p_k_part), intent(in) :: k, uk0, fuk0
  real(p_k_part), intent(in) :: o, uo0, fuo0
  real(p_k_part), intent(out) :: ch_k, shc_k, chc_k, cn_o, snc_o, cnc_o, t

  real(p_k_part) :: fun, dfun, t0
  integer :: iter

  t0 = dt / ( uk0 + uo0 )
  iter = 0

  do

    call  trig_func( o * t0, cn_o, snc_o, cnc_o )
    call hyper_func( k * t0, ch_k, shc_k, chc_k )

    fun = ( uk0 * shc_k + fuk0 * chc_k * t0 ) * t0 + &
          ( uo0 * snc_o - fuo0 * cnc_o * t0 ) * t0 - dt

    dfun = uk0 * ch_k + fuk0 * shc_k * t0 + &
           uo0 * cn_o + fuo0 * snc_o * t0

    t = t0 - fun / dfun
    iter = iter + 1

    if ( abs(t - t0) / t0 < tol .or. iter > p_iter_max ) then
      exit
    else
      t0 = t
    endif

  enddo

end subroutine get_t
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dudt_exact_rr( this, emf, dt, i0, i1, energy, time )
!-----------------------------------------------------------------------------------------
! Advance velocities and calculate time-centered total energy
! process particles in the range [i0,i1]
! This pusher follows the ref. https://arxiv.org/abs/2007.07556
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  real(p_double), intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_double), intent(inout) :: energy
  real(p_double), intent(in) :: time

  integer :: i, j, ptrcur, np, pp
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep
  real(p_k_part), dimension(4) :: u0, u, fu0, f2u0, f3u0, uk0, fuk0, uk, uo0, fuo0, uo
  real(p_k_part) :: t, tmp, k2o2, amp_k, amp_o, amp
  real(p_k_part) :: ch_k, shc_k, chc_k, cn_o, snc_o, cnc_o
  real(p_k_part) :: c1, c2, c3, c4
  real(p_k_part) :: inv1, inv2, e2, b2, sigma0, alpha, qm, uk0_sq, uo0_sq, y0
  real(p_double) :: loc_ene
  real(p_k_part), dimension(p_cache_size) :: k, o, k2, o2
  integer, dimension(p_cache_size) :: ix1, ix2
  integer :: n1, n2

  real(p_double), save :: rel_tol = 1.0d-12
  real(p_double), save :: weak_fld_tol = 1.0d-12

  ! executable statements
  qm = 1.0_p_k_part / this%rqm
  sigma0 = this%k_rr * this%q_real * qm

  !loop through all particles
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    call this % get_emf( emf, bp, ep, np, ptrcur, t=time )

    loc_ene = 0.0_p_double

    pp = ptrcur
    do i = 1, np
      ep(1,i) = ep(1,i) * qm
      ep(2,i) = ep(2,i) * qm
      ep(3,i) = ep(3,i) * qm
      bp(1,i) = bp(1,i) * qm
      bp(2,i) = bp(2,i) * qm
      bp(3,i) = bp(3,i) * qm
      pp = pp + 1
    enddo

    ! -------------------------------------------------------------------------
    ! categorize particles
    ! -------------------------------------------------------------------------
    n1 = 0; ix1 = 0 
    n2 = 0; ix2 = 0
    do i = 1, np

      e2 = ep(1,i) * ep(1,i) + ep(2,i) * ep(2,i) + ep(3,i) * ep(3,i)
      b2 = bp(1,i) * bp(1,i) + bp(2,i) * bp(2,i) + bp(3,i) * bp(3,i)      

      if ( e2 + b2 < weak_fld_tol ) then

        k2(i) = 0.0_p_k_part; k(i) = 0.0_p_k_part
        o2(i) = 0.0_p_k_part; o(i) = 0.0_p_k_part

      else

        ! Lorentz invariants
        inv1 = e2 - b2
        inv2 = ep(1,i) * bp(1,i) + ep(2,i) * bp(2,i) + ep(3,i) * bp(3,i)

        ! calculate eigenvalues
        tmp = sqrt( inv1 * inv1 + 4.0_p_k_part * inv2 * inv2 )
        k2(i) = 0.5_p_k_part * ( inv1 + tmp ); k(i) = sqrt(k2(i))
        o2(i) = 0.5_p_k_part * (-inv1 + tmp ); o(i) = sqrt(o2(i))

        if ( k2(i) < rel_tol .and. o2(i) < rel_tol ) then
          ! light-like
          n1 = n1 + 1; ix1(n1) = i
        else
          ! general case
          n2 = n2 + 1; ix2(n2) = i
        endif

      endif
    enddo

    ! -------------------------------------------------------------------------
    ! advance particles in light-like field (category 1)
    ! -------------------------------------------------------------------------
    do j = 1, n1

      i = ix1(j); pp = ptrcur + i - 1

      u0(1:3) = this%p(:,pp)
      u0(4) = sqrt( 1.0_p_k_part + u0(1)*u0(1) + u0(2)*u0(2) + u0(3)*u0(3) )

      fu0  = get_fu( ep(:,i), bp(:,i), u0 )
      f2u0 = get_fu( ep(:,i), bp(:,i), fu0 )
      f3u0 = get_fu( ep(:,i), bp(:,i), f2u0 )
      y0 = dot4( u0, f2u0 )

      ! solve for the proper time step according to dt
      call get_t_rr_light( dt, k(i), o(i), u0(4), fu0(4), f2u0(4), f3u0(4), &
        y0, sigma0, t, amp, c1, c2, c3, c4, this%iter_tol )
      u = ( c1 * u0 + c2 * fu0 + c3 * f2u0 + c4 * f3u0 ) * amp

      this%p(:,pp) = u(1:3)

    enddo

    ! -------------------------------------------------------------------------
    ! advance particles in arbitrary field (category 2)
    ! -------------------------------------------------------------------------
    do j = 1, n2

      i = ix2(j); pp = ptrcur + i - 1

      u0(1:3) = this%p(:,pp)
      u0(4) = sqrt( 1.0_p_k_part + u0(1)*u0(1) + u0(2)*u0(2) + u0(3)*u0(3) )

      ! subspace decomposition
      k2o2 = k2(i) + o2(i)
      fu0  = get_fu( ep(:,i), bp(:,i), u0 )
      f2u0 = get_fu( ep(:,i), bp(:,i), fu0 )
      uk0  = ( f2u0 + o2(i) * u0 ) / k2o2
      uo0  = u0 - uk0

      uk0_sq = dot4( uk0, uk0 )
      uo0_sq = 1.0_p_k_part - uk0_sq

      fuk0 = get_fu( ep(:,i), bp(:,i), uk0 )
      fuo0 = fu0 - fuk0

      alpha = sigma0 * k2o2

      ! solve for the proper time step according to dt
      call get_t_rr( dt, k(i), uk0(4), fuk0(4), uk0_sq, o(i), uo0(4), fuo0(4), &
        uo0_sq, ch_k, shc_k, chc_k, cn_o, snc_o, cnc_o, t, amp_k, amp_o, alpha, this%iter_tol )

      uk = ( uk0 * ch_k + fuk0 * shc_k * t ) * amp_k
      uo = ( uo0 * cn_o + fuo0 * snc_o * t ) * amp_o
      u = uk + uo

      this%p(1:3,pp) = u(1:3)

    enddo

    pp = ptrcur
    do i = 1, np
      tmp = sqrt(1.0 + this%p(1,pp)*this%p(1,pp) + this%p(2,pp)*this%p(2,pp) + this%p(3,pp)*this%p(3,pp))
      loc_ene = loc_ene + this%q(pp) * (tmp - 1.0)
      pp = pp + 1
    enddo

    energy = energy + loc_ene

  enddo

end subroutine dudt_exact_rr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_t_rr_light( dt, k, o, u0, fu0, f2u0, f3u0, y0, sigma0, t, amp, c1, c2, c3, c4, tol )
!-----------------------------------------------------------------------------------------
! Using Newton-Raphson method to calculate proper time step for light-like case
!-----------------------------------------------------------------------------------------

  implicit none

  real(p_double), intent(in) :: dt, tol
  real(p_k_part), intent(in) :: k, o, u0, fu0, f2u0, f3u0, sigma0, y0
  real(p_k_part), intent(out) :: t, amp, c1, c2, c3, c4

  real(p_k_part) :: t1, t2, fun, dfun, y0s0t, s0t, k2, o2, o2k2_h, a1, a2, a3, a4
  integer :: iter

  k2 = k * k; o2 = o * o
  o2k2_h = 0.5_p_double * ( o2 - k2 )
  t1 = dt / u0
  iter = 0

  do

    t2    = t1 * t1
    s0t   = t1 * sigma0
    y0s0t = y0 * s0t

    a1 = ( 1.0_p_double + y0s0t ) * t1
    a2 = ( 0.5_p_double + 0.333333333333333_p_double * y0s0t ) * t2
    a3 = ( 0.5_p_double * s0t + 0.166666666666667_p_double * t2 + 0.125_p_double * y0s0t * t2 ) * t1
    a4 = ( s0t + 0.125_p_double * t2 + 0.1_p_double * y0s0t * t2 ) * 0.333333333333333_p_double * t2

    fun  = a1 * u0 + a2 * fu0 + a3 * f2u0 + a4 * f3u0 - dt

    amp = 1.0_p_double / sqrt( 1.0 + 2.0 * (o2k2_h - y0) * s0t )
    c1 = 1.0_p_double + o2k2_h * s0t
    c2 = c1 * t1
    c3 = s0t + 0.5_p_double * t2
    c4 = (s0t + 0.166666666666667_p_double * t2) * t1

    dfun = ( c1 * u0 + c2 * fu0 + c3 * f2u0 + c4 * f3u0 ) * amp

    t = t1 - fun / dfun
    iter = iter + 1

    if ( abs(t - t1) / t1 < tol .or. iter > p_iter_max ) then
      exit
    else
      t1 = t
    endif

  enddo

end subroutine get_t_rr_light
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_t_rr( dt, k, uk0, fuk0, uk0_sq, o, uo0, fuo0, uo0_sq, &
  ch_k, shc_k, chc_k, cn_o, snc_o, cnc_o, t, amp_k, amp_o, alpha, tol )
!-----------------------------------------------------------------------------------------
! Using Newton-Raphson method to calculate proper time step for general case
!-----------------------------------------------------------------------------------------

  implicit none

  real(p_double), intent(in) :: dt, tol
  real(p_k_part), intent(in) :: k, uk0, fuk0, o, uo0, fuo0, alpha, uk0_sq, uo0_sq
  real(p_k_part), intent(out) :: ch_k, shc_k, chc_k, cn_o, snc_o, cnc_o, t, amp_k, amp_o

  real(p_k_part) :: fun, dfun, t0, expn, tmp_k, tmp_o, shcc_k, sncc_o
  integer :: iter

  t0 = dt / ( uk0 + uo0 )
  iter = 0

  do

    call trig_func( o * t0, cn_o, snc_o, cnc_o, sncc_o )
    call hyper_func( k * t0, ch_k, shc_k, chc_k, shcc_k )

    tmp_k = ( uk0 * shc_k + fuk0 * chc_k * t0 ) * t0
    tmp_k = tmp_k + alpha * uo0_sq * ( tmp_k - ( uk0 * chc_k + fuk0 * shcc_k * t0 ) * t0 ) * t0

    tmp_o = ( uo0 * snc_o - fuo0 * cnc_o * t0 ) * t0
    tmp_o = tmp_o - alpha * uk0_sq * ( tmp_o + ( uo0 * cnc_o + fuo0 * sncc_o * t0 ) * t0 ) * t0

    fun = tmp_k + tmp_o - dt

    expn = exp( -2.0 * alpha * t0 )

    amp_k = 1.0_p_k_part / sqrt( uk0_sq + uo0_sq * expn )
    amp_o = 1.0_p_k_part / sqrt( uo0_sq + uk0_sq / expn )

    dfun = ( uk0 * ch_k + fuk0 * shc_k * t0 ) * amp_k + &
           ( uo0 * cn_o + fuo0 * snc_o * t0 ) * amp_o

    t = t0 - fun / dfun
    iter = iter + 1

    if ( abs(t - t0) / t0 < tol .or. iter > p_iter_max ) then
      exit
    else
      t0 = t
    endif

  enddo

end subroutine get_t_rr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine trig_func( x, cn, snc, cnc, sncc )
!-----------------------------------------------------------------------------------------
! Fast evaluation of the trigonometric functions cos(x), sinc(x), (cos(x)-1)/x^2 and 
! (sinc(x)-1)/x^2. When x is smaller than a threshold value, they are evaluated using
! Taylor expansion, otherwise using tan(x/2) for evaluation.
!-----------------------------------------------------------------------------------------
  implicit none
  real(p_k_part), intent(in) :: x
  real(p_k_part), intent(out) :: cn, snc, cnc
  real(p_k_part), intent(out), optional :: sncc

  real(p_double) :: xd, x2, x4, tmp1, tmp2, tmp3, cn_d, snc_d
  ! threshold making Taylor expansion accurate to machine precision
  real(p_double), save :: eps = 0.007367081148666_p_double

  xd = real(x,p_double)
  x2 = xd * xd

  if ( abs(xd) < eps ) then
    x4 = x2 * x2
    cn  = real( 1.0_p_double - 0.500000000000000_p_double * x2 + 0.041666666666667_p_double * x4, p_k_part )
    snc = real( 1.0_p_double - 0.166666666666667_p_double * x2 + 0.008333333333333_p_double * x4, p_k_part )
    cnc =real( -0.5_p_double + 0.041666666666667_p_double * x2 - 0.001388888888889_p_double * x4, p_k_part )
    if ( present(sncc) ) then
      sncc = real( -0.166666666666667_p_double + 0.008333333333333_p_double * x2 - 1.984126984126984d-4 * x4, p_k_part )
    endif
  else
    tmp1 = tan( 0.5_p_double * xd )
    tmp2 = tmp1 * tmp1
    tmp3 = 1.0_p_double / (1.0_p_double + tmp2)
    cn_d = ( 1.0_p_double - tmp2 ) * tmp3
    cn = real( cn_d, p_k_part )
    snc_d = 2.0_p_double * tmp1 * tmp3 / xd
    snc = real( snc_d, p_k_part )
    cnc = real( ( cn_d - 1.0_p_double ) / x2, p_k_part )
    if ( present(sncc) ) then
      sncc = real( ( snc_d - 1.0_p_double ) / x2, p_k_part )
    endif
  endif

end subroutine trig_func

!-----------------------------------------------------------------------------------------
subroutine hyper_func( x, ch, shc, chc, shcc )
!-----------------------------------------------------------------------------------------
! Fast evaluation of the hyperbolic functions cosh(x), sinhc(x), (cosh(x)-1)/x^2 and 
! (sinhc(x)-1)/x^2. When x is smaller than a threshold value, they are evaluated using
! Taylor expansion, otherwise using tanh(x/2) for evaluation.
!-----------------------------------------------------------------------------------------
  implicit none
  real(p_k_part), intent(in) :: x
  real(p_k_part), intent(out) :: ch, shc, chc
  real(p_k_part), intent(out), optional :: shcc

  real(p_double) :: xd, x2, x4, tmp1, tmp2, tmp3, ch_d, shc_d
  ! threshold making Taylor expansion accurate to machine precision
  real(p_double), save :: eps = 0.007367081148666_p_double

  xd = real(x,p_double)
  x2 = xd * xd

  if ( abs(xd) < eps ) then
    x4 = x2 * x2
    ch  = real( 1.0_p_double + 0.500000000000000_p_double * x2 + 0.041666666666667_p_double * x4, p_k_part )
    shc = real( 1.0_p_double + 0.166666666666667_p_double * x2 + 0.008333333333333_p_double * x4, p_k_part )
    chc = real( 0.5_p_double + 0.041666666666667_p_double * x2 + 0.001388888888889_p_double * x4, p_k_part )
    if ( present(shcc) ) then
      shcc = real( 0.166666666666667_p_double + 0.008333333333333_p_double * x2 + 1.984126984126984d-4 * x4, p_k_part )
    endif
  else
    tmp1 = tanh( 0.5_p_double * xd )
    tmp2 = tmp1 * tmp1
    tmp3 = 1.0_p_double / (1.0_p_double - tmp2)
    ch_d = ( 1.0_p_double + tmp2 ) * tmp3
    ch = real( ch_d, p_k_part )
    shc_d = 2.0_p_double * tmp1 * tmp3 / xd
    shc = real( shc_d, p_k_part )
    chc = real( ( ch_d - 1.0_p_double ) / x2, p_k_part )
    if ( present(shcc) ) then
      shcc = real( ( shc_d - 1.0_p_double ) / x2, p_k_part )
    endif
  endif

end subroutine hyper_func

!-----------------------------------------------------------------------------------------
function dot4( u, v )

  implicit none
  real(p_k_part), intent(in), dimension(4) :: u, v
  real(p_k_part) :: dot4

  dot4 = u(4) * v(4) - u(1) * v(1) - u(2) * v(2) - u(3) * v(3)

end function dot4

!-----------------------------------------------------------------------------------------
function get_fu( ep, bp, u )
!-----------------------------------------------------------------------------------------
! Calculate the product of field tensor F and a four-vector u
!-----------------------------------------------------------------------------------------
  implicit none
  real(p_k_part), intent(in), dimension(3) :: ep, bp
  real(p_k_part), intent(in), dimension(4) :: u
  real(p_k_part), dimension(4) :: get_fu

  get_fu(4) =                ep(1) * u(1) + ep(2) * u(2) + ep(3) * u(3)
  get_fu(1) = ep(1) * u(4) +                bp(3) * u(2) - bp(2) * u(3)
  get_fu(2) = ep(2) * u(4) - bp(3) * u(1) +                bp(1) * u(3)
  get_fu(3) = ep(3) * u(4) + bp(2) * u(1) - bp(1) * u(2)
end function get_fu

end module m_species_extpush
