!-----------------------------------------------------------------------------------------
! Particle pusher
!-----------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_push

#include "memory/memory.h"

use m_parameters

use m_species_define
use m_species_current

use m_vdf_define

use m_emf_define
use m_emf
use m_emf_interpolate

implicit none

private

interface dudt_boris
  module procedure dudt_boris
end interface

interface dudt_vay
  module procedure dudt_vay
end interface

interface dudt_fullrot
  module procedure dudt_fullrot
end interface

interface dudt_euler
  module procedure dudt_euler
end interface

interface dudt_cond_vay
  module procedure dudt_cond_vay
end interface

interface dudt_cary
  module procedure dudt_cary
end interface

#ifdef __HAS_SPIN__

interface dsdt_vieira
  module procedure dsdt_vieira
end interface

public :: dsdt_vieira
#endif

interface advance_deposit
  module procedure advance_deposit
end interface

interface ntrim
  module procedure ntrim
end interface

public :: advance_deposit

public :: dudt_boris, dudt_vay, dudt_fullrot, dudt_euler
public :: dudt_cond_vay, dudt_cary

public :: ntrim

contains

!-----------------------------------------------------------------------------------------
! Advance velocities and calculate time-centered total energy using the full calculation
! of the tangent of the velocity rotation angle.
!
! process particles in the range [i0,i1]
!-----------------------------------------------------------------------------------------
subroutine dudt_fullrot( this, emf, dt, i0, i1, energy, time )

  implicit none

  ! dummy variables
  class( t_species ),    intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  real(p_double),                 intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_double), intent(inout) :: energy
  real(p_double), intent(in) :: time

  ! local variables
  real(p_k_part) :: tem

  integer :: i, ptrcur, np, pp

  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp
  real(p_k_part), dimension(p_cache_size) ::  tem_gam
  real(p_k_part) :: gamma, u2, b2, t, r, s !, bnorm
  real(p_double) :: loc_ene

  ! executable statements
  !print *, '(*warn*) Using dudt_fullrot'

  ! get factor that includes timestep and charge-to-mass ratio
  ! note that the charge-to-mass ratio is used since the
  ! momentum p is not the total momentum of a particles but
  ! the momentum per unit restmass (which is the electron mass)

  tem = real( 0.5_p_double * (dt / this%rqm), p_k_part )

  !loop through all particles
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! modify bp & ep to include timestep and charge-to-mass ratio
    ! and perform half the electric field acceleration.
    ! Result is stored in UTEMP.

    call this % get_emf( emf, bp, ep, np, ptrcur, t=time )

    do i=1, np
      ep(1,i) = ep(1,i) * tem
      ep(2,i) = ep(2,i) * tem
      ep(3,i) = ep(3,i) * tem
    end do

    loc_ene = 0

    pp = ptrcur
    do i=1,np
      utemp(1,i) = this%p(1,pp) + ep(1,i)
      utemp(2,i) = this%p(2,pp) + ep(2,i)
      utemp(3,i) = this%p(3,pp) + ep(3,i)

      ! Get time centered gamma
      u2 =  (utemp(1,i)**2 + utemp(2,i)**2) + utemp(3,i)**2
      gamma = sqrt( u2 + 1 )

      ! accumulate time centered energy
      ! this is done in double precicion always
      loc_ene = loc_ene + this%q(pp) * u2 / (gamma + 1.0_p_double)

      tem_gam(i)= tem / gamma

      pp = pp + 1
    enddo

    ! accumulate global energy
    energy = energy + loc_ene

    do i=1,np
      bp(1,i) = bp(1,i)*tem_gam(i)
      bp(2,i) = bp(2,i)*tem_gam(i)
      bp(3,i) = bp(3,i)*tem_gam(i)
    end do

    do i=1,np
      b2 = (bp(1,i)**2 + bp(2,i)**2) + bp(3,i)**2

      ! Full tangent calculation
      !if ( b2 > 0 ) then
      ! bnorm = sqrt(b2)
      !  t = tan( bnorm ) / bnorm
      !else
      !  t = 1
      !endif

      ! rational approximation to tan(b)/b, approx. 20.3 digits accurate
      r = 2005892.32789228135999056d0 + b2*(-267559.850143894287032640d0 + b2*( 6870.7950060474599763728d0 - 35.8427314349709199363798d0*b2)) 
      s = 2005892.32789228136020125d0 + b2*(-936190.626107988089080768d0 + b2*(51482.0266564063006502261d0 + b2*(b2-625.608053837898687047447d0)))
      
      ! Avoid division by 0 if b = \pi/2
      if ( s /= 0 ) then
        t = r/s
      else
        ! Use first order taylor approximation
        t = 1
      endif

      bp(1,i) = bp(1,i)*t
      bp(2,i) = bp(2,i)*t
      bp(3,i) = bp(3,i)*t
    end do

    ! Get u' and store it in p
    pp = ptrcur
    do i=1,np
      this%p(1,pp) = utemp(1,i) + utemp(2,i) * bp(3,i) - utemp(3,i) * bp(2,i)
      this%p(2,pp) = utemp(2,i) + utemp(3,i) * bp(1,i) - utemp(1,i) * bp(3,i)
      this%p(3,pp) = utemp(3,i) + utemp(1,i) * bp(2,i) - utemp(2,i) * bp(1,i)
      pp = pp + 1
    end do

    do i=1,np
      s = 2.0_p_k_part / ( ((1.0_p_k_part + bp(1,i)**2) + bp(2,i)**2) + bp(3,i)**2)
      bp(1,i) = bp(1,i) * s
      bp(2,i) = bp(2,i) * s
      bp(3,i) = bp(3,i) * s
    end do

    ! Get u+ and store it in utemp
    pp = ptrcur
    do i=1,np
      utemp(1,i) = utemp(1,i) + this%p(2,pp) * bp(3,i) - this%p(3,pp) * bp(2,i)
      utemp(2,i) = utemp(2,i) + this%p(3,pp) * bp(1,i) - this%p(1,pp) * bp(3,i)
      utemp(3,i) = utemp(3,i) + this%p(1,pp) * bp(2,i) - this%p(2,pp) * bp(1,i)
      pp = pp + 1
    end do

    ! Perform second half of electric field acceleration.
    pp = ptrcur
    do i=1,np
      this%p(1,pp) = utemp(1,i) + ep(1,i)
      this%p(2,pp) = utemp(2,i) + ep(2,i)
      this%p(3,pp) = utemp(3,i) + ep(3,i)
      pp = pp + 1
    end do

  enddo

end subroutine dudt_fullrot
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dudt_euler( this, emf, dt, i0, i1, energy, t )

  implicit none

  class( t_species ),    intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  real(p_double),                 intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_double), intent(inout) :: energy
  real(p_double), intent(in) :: t

  real(p_k_part) :: tem
  integer :: i, ptrcur, np, pp
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp
  real(p_k_part), dimension(p_cache_size) :: gam_tem
  real(p_k_part) :: gamma, u2, b2, bnorm
  real(p_k_part) :: s, a, b, c, d
  real(p_k_part) :: r11, r21, r31, r12, r22, r32, r13, r23, r33
  
  real(p_double) :: loc_ene
  
  tem = real( 0.5_p_double * (dt / this%rqm), p_k_part )

  !loop through all particles
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! modify bp & ep to include timestep and charge-to-mass ratio
    ! and perform half the electric field acceleration.
    ! Result is stored in UTEMP.

    call this % get_emf( emf, bp, ep, np, ptrcur, t=t )

    do i=1, np
      ep(1,i) = ep(1,i) * tem
      ep(2,i) = ep(2,i) * tem
      ep(3,i) = ep(3,i) * tem
    end do

    loc_ene = 0

    pp = ptrcur
    do i=1,np
      utemp(1,i) = this%p(1,pp) + ep(1,i)
      utemp(2,i) = this%p(2,pp) + ep(2,i)
      utemp(3,i) = this%p(3,pp) + ep(3,i)

      ! Get time centered gamma
      u2 =  (utemp(1,i)**2 + utemp(2,i)**2) + utemp(3,i)**2
      gamma = sqrt( u2 + 1 )

      ! accumulate time centered energy
      ! this is done in double precicion always
      loc_ene = loc_ene + this%q(pp) * u2 / (gamma + 1.0_p_double)

      gam_tem(i)= dt / (this%rqm * gamma )

      pp = pp + 1
    enddo

    ! accumulate global energy
    energy = energy + loc_ene

    do i=1,np
      bp(1,i) = bp(1,i)*gam_tem(i)
      bp(2,i) = bp(2,i)*gam_tem(i)
      bp(3,i) = bp(3,i)*gam_tem(i)
    end do

    pp = ptrcur

    do i=1,np
      b2 = (bp(1,i)**2 + bp(2,i)**2) + bp(3,i)**2
      
      bnorm = sqrt(b2)

      s = sin( bnorm / 2 )
      a = cos( bnorm / 2 )

      if ( b2 > 0 ) then
        s = -s / bnorm  
      else
        s = -1
      endif

      b = bp(1,i) * s
      c = bp(2,i) * s
      d = bp(3,i) * s

      r11 = a*a+b*b-c*c-d*d;  r21=2*(b*c-a*d);      r31=2*(b*d+a*c)
      r12 = 2*(b*c+a*d);      r22=a*a+c*c-b*b-d*d;  r32=2*(c*d-a*b)
      r13 = 2*(b*d-a*c);      r23=2*(c*d+a*b);      r33=a*a+d*d-b*b-c*c
      
      ! get u+ and store in p
      this%p(1,pp) = r11 * utemp(1,i) + r21 * utemp(2,i) + r31 * utemp(3,i) 
      this%p(2,pp) = r12 * utemp(1,i) + r22 * utemp(2,i) + r32 * utemp(3,i) 
      this%p(3,pp) = r13 * utemp(1,i) + r23 * utemp(2,i) + r33 * utemp(3,i) 
      pp = pp + 1
    end do

    ! Perform second half of electric field acceleration.
    pp = ptrcur
    do i=1,np
      this%p(1,pp) = this%p(1,pp) + ep(1,i)
      this%p(2,pp) = this%p(2,pp) + ep(2,i)
      this%p(3,pp) = this%p(3,pp) + ep(3,i)
      pp = pp + 1
    end do

  enddo

end subroutine dudt_euler
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dudt_boris( this, emf, dt, i0, i1, energy, t )
!-----------------------------------------------------------------------------------------
! Advance velocities and calculate time-centered total energy
! process particles in the range [i0,i1]
! This routine is meant to be used for testing simd code and is written to use the same numerics
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species ),    intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  real(p_double),                 intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_double), intent(inout) :: energy
  real(p_double), intent(in) :: t

  real(p_k_part) :: tem
  integer :: i, ptrcur, np, pp
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp

#ifdef __HAS_SPIN__
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: p_old
  real(p_k_part), dimension(p_cache_size) :: gamma
#else
  real(p_k_part) :: gamma
#endif

  real(p_k_part), dimension(p_cache_size) ::  gam_tem, otsq
  real(p_k_part) :: u2
  real(p_double) :: loc_ene

  real(p_k_part), dimension(p_p_dim) :: fp, vtmp
  real(p_k_part) :: f2, gam, v_dot_e, fac, cor

  ! get factor that includes timestep and charge-to-mass ratio
  ! note that the charge-to-mass ratio is used since the
  ! momentum p is not the total momentum of a particles but
  ! the momentum per unit restmass (which is the electron mass)

  fac = this%k_rr * this%q_real * dt * 0.5_p_k_part / this%rqm
  tem = real( 0.5_p_double * (dt / this%rqm), p_k_part )

  !loop through all particles
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! modify bp & ep to include timestep and charge-to-mass ratio
    ! and perform half the electric field acceleration.
    ! Result is stored in UTEMP.

    call this % get_emf( emf, bp, ep, np, ptrcur, t=t )

    ! half push of radiation reaction
    if ( this%rad_react ) then
      pp = ptrcur
      do i = 1, np

        u2 = this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2
        gam = sqrt( 1.0_p_k_part + u2 )

        vtmp(1) = this%p(1,pp) / gam
        vtmp(2) = this%p(2,pp) / gam
        vtmp(3) = this%p(3,pp) / gam

        ! Lorentz force
        fp(1) = ep(1,i) + vtmp(2) * bp(3,i) - vtmp(3) * bp(2,i)
        fp(2) = ep(2,i) + vtmp(3) * bp(1,i) - vtmp(1) * bp(3,i)
        fp(3) = ep(3,i) + vtmp(1) * bp(2,i) - vtmp(2) * bp(1,i)

        f2 = fp(1)*fp(1) + fp(2)*fp(2) + fp(3)*fp(3)

        v_dot_e = vtmp(1)*ep(1,i) + vtmp(2)*ep(2,i) + vtmp(3)*ep(3,i)

        f2 = f2 - v_dot_e * v_dot_e
        utemp(1,i) = -gam * this%p(1,pp) * f2
        utemp(2,i) = -gam * this%p(2,pp) * f2
        utemp(3,i) = -gam * this%p(3,pp) * f2

        utemp(1,i) = utemp(1,i) + v_dot_e * ep(1,i)
        utemp(2,i) = utemp(2,i) + v_dot_e * ep(2,i)
        utemp(3,i) = utemp(3,i) + v_dot_e * ep(3,i)

        utemp(1,i) = utemp(1,i) + fp(2) * bp(3,i) - fp(3) * bp(2,i)
        utemp(2,i) = utemp(2,i) + fp(3) * bp(1,i) - fp(1) * bp(3,i)
        utemp(3,i) = utemp(3,i) + fp(1) * bp(2,i) - fp(2) * bp(1,i)

        this%p(1,pp) = this%p(1,pp) + fac * utemp(1,i)
        this%p(2,pp) = this%p(2,pp) + fac * utemp(2,i)
        this%p(3,pp) = this%p(3,pp) + fac * utemp(3,i)

        pp = pp + 1
      enddo
    endif

#ifdef __HAS_SPIN__
    ! store the old momentum for spin pusher
    pp = ptrcur
    do i=1, np
      p_old(1,i) = this%p(1,pp)
      p_old(2,i) = this%p(2,pp)
      p_old(3,i) = this%p(3,pp)
      pp = pp + 1
    end do
#endif

    do i=1, np
      ep(1,i) = ep(1,i) * tem
      ep(2,i) = ep(2,i) * tem
      ep(3,i) = ep(3,i) * tem
    end do

    ! Perform first half of electric field acceleration.
    ! and get time centered gamma

    loc_ene = 0

    pp = ptrcur
    do i=1,np
      utemp(1,i) = this%p(1,pp) + ep(1,i)
      utemp(2,i) = this%p(2,pp) + ep(2,i)
      utemp(3,i) = this%p(3,pp) + ep(3,i)

      ! Get time centered gamma
      u2 =  (utemp(1,i)**2 + utemp(2,i)**2) + utemp(3,i)**2

#ifdef __HAS_SPIN__

      gamma(i) = sqrt( u2 + 1 )
      ! accumulate time centered energy
      ! this is done in double precicion always
      loc_ene = loc_ene + this%q(pp) * u2 / (gamma(i) + 1.0_p_double)
      gam_tem(i)= tem / gamma(i)

#else
      gamma = sqrt( u2 + 1 )

      ! accumulate time centered energy
      ! this is done in double precicion always
      loc_ene = loc_ene + this%q(pp) * u2 / (gamma + 1.0_p_double)

      gam_tem(i)= tem / gamma

#endif

      pp = pp + 1
    enddo

    ! accumulate global energy
    energy = energy + loc_ene

    do i=1,np
      bp(1,i) = bp(1,i)*gam_tem(i)
      bp(2,i) = bp(2,i)*gam_tem(i)
      bp(3,i) = bp(3,i)*gam_tem(i)
    end do

    ! Get u'
    pp = ptrcur
    do i=1,np
      this%p(1,pp) = utemp(1,i) + utemp(2,i) * bp(3,i)
      this%p(2,pp) = utemp(2,i) + utemp(3,i) * bp(1,i)
      this%p(3,pp) = utemp(3,i) + utemp(1,i) * bp(2,i)
      pp = pp + 1
    end do

    pp = ptrcur
    do i=1,np
      this%p(1,pp) = this%p(1,pp) - utemp(3,i) * bp(2,i)
      this%p(2,pp) = this%p(2,pp) - utemp(1,i) * bp(3,i)
      this%p(3,pp) = this%p(3,pp) - utemp(2,i) * bp(1,i)
      pp = pp + 1
    end do

#ifdef __HAS_SPIN__
    ! push spin
    call this%dsdt( ep, bp, p_old, gamma, dt, ptrcur, np )
#endif

    do i=1,np
      otsq(i) = 2.0_p_k_part / ( ((1.0_p_k_part + bp(1,i)**2) + bp(2,i)**2) + bp(3,i)**2)
    end do

    do i=1,np
      bp(1,i) = bp(1,i) * otsq(i)
      bp(2,i) = bp(2,i) * otsq(i)
      bp(3,i) = bp(3,i) * otsq(i)
    end do

    ! Rotate from u- to u+ using u'
    pp = ptrcur
    do i=1,np
      utemp(1,i) = utemp(1,i) + this%p(2,pp) * bp(3,i)
      utemp(2,i) = utemp(2,i) + this%p(3,pp) * bp(1,i)
      utemp(3,i) = utemp(3,i) + this%p(1,pp) * bp(2,i)
      pp = pp + 1
    end do

    pp = ptrcur
    do i=1,np
      utemp(1,i) = utemp(1,i) - this%p(3,pp) * bp(2,i)
      utemp(2,i) = utemp(2,i) - this%p(1,pp) * bp(3,i)
      utemp(3,i) = utemp(3,i) - this%p(2,pp) * bp(1,i)
      pp = pp + 1
    end do

    ! Perform second half of electric field acceleration.
    pp = ptrcur
    do i=1,np
      this%p(1,pp) = utemp(1,i) + ep(1,i)
      this%p(2,pp) = utemp(2,i) + ep(2,i)
      this%p(3,pp) = utemp(3,i) + ep(3,i)
      pp = pp + 1
    end do


    ! second half push of radiation reaction
    if ( this%rad_react ) then
      pp = ptrcur
      do i = 1, np

        ! retrieve the fields
        cor = 1.0_p_k_part / tem
        ep(1,i) = ep(1,i) * cor
        ep(2,i) = ep(2,i) * cor
        ep(3,i) = ep(3,i) * cor

        cor = 1.0_p_k_part / ( gam_tem(i) * otsq(i) )
        bp(1,i) = bp(1,i) * cor
        bp(2,i) = bp(2,i) * cor
        bp(3,i) = bp(3,i) * cor

        u2 = this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2
        gam = sqrt( 1.0_p_k_part + u2 )

        vtmp(1) = this%p(1,pp) / gam
        vtmp(2) = this%p(2,pp) / gam
        vtmp(3) = this%p(3,pp) / gam

        ! Lorentz force
        fp(1) = ep(1,i) + vtmp(2) * bp(3,i) - vtmp(3) * bp(2,i)
        fp(2) = ep(2,i) + vtmp(3) * bp(1,i) - vtmp(1) * bp(3,i)
        fp(3) = ep(3,i) + vtmp(1) * bp(2,i) - vtmp(2) * bp(1,i)

        f2 = fp(1)*fp(1) + fp(2)*fp(2) + fp(3)*fp(3)

        v_dot_e = vtmp(1)*ep(1,i) + vtmp(2)*ep(2,i) + vtmp(3)*ep(3,i)

        f2 = f2 - v_dot_e * v_dot_e
        utemp(1,i) = -gam * this%p(1,pp) * f2
        utemp(2,i) = -gam * this%p(2,pp) * f2
        utemp(3,i) = -gam * this%p(3,pp) * f2

        utemp(1,i) = utemp(1,i) + v_dot_e * ep(1,i)
        utemp(2,i) = utemp(2,i) + v_dot_e * ep(2,i)
        utemp(3,i) = utemp(3,i) + v_dot_e * ep(3,i)

        utemp(1,i) = utemp(1,i) + fp(2) * bp(3,i) - fp(3) * bp(2,i)
        utemp(2,i) = utemp(2,i) + fp(3) * bp(1,i) - fp(1) * bp(3,i)
        utemp(3,i) = utemp(3,i) + fp(1) * bp(2,i) - fp(2) * bp(1,i)

        this%p(1,pp) = this%p(1,pp) + fac * utemp(1,i)
        this%p(2,pp) = this%p(2,pp) + fac * utemp(2,i)
        this%p(3,pp) = this%p(3,pp) + fac * utemp(3,i)

        pp = pp + 1
      enddo
    endif

  enddo

end subroutine dudt_boris
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dudt_vay( this, emf, dt, i0, i1, energy, t )
!-----------------------------------------------------------------------------------------
! Advance velocities and calculate time-centered total energy
! process particles in the range [i0,i1] using VAY pusher [Ref. xxx]
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species ),    intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  real(p_double),                 intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_double), intent(inout) :: energy
  real(p_double), intent(in) :: t

  real(p_k_part) :: tem, gamma_c
  real(p_k_part), dimension(p_p_dim) :: uc
  real(p_k_part) :: uc2
  integer :: i, ptrcur, np, pp
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp, temp_vec
  real(p_k_part), dimension(p_cache_size) :: rgam, ustar, sigma, spar, uptp, bpsq
  real( p_double ) :: loc_ene

  ! get factor that includes timestep and charge-to-mass ratio
  ! note that the charge-to-mass ratio is used since the
  ! momentum p is not the total momentum of a particles but
  ! the momentum per unit restmass (which is the electron mass)
  tem = real( 0.5_p_double * dt / this%rqm, p_k_part )

  !loop through all particles
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    call this % get_emf( emf, bp, ep, np, ptrcur, t=t )

    pp = ptrcur
    do i=1, np
      ep(1,i) = ep(1,i) * tem
      ep(2,i) = ep(2,i) * tem
      ep(3,i) = ep(3,i) * tem

      bp(1,i) = bp(1,i) * tem
      bp(2,i) = bp(2,i) * tem
      bp(3,i) = bp(3,i) * tem

      rgam(i) = 1.0_p_k_part / sqrt(1.0_p_k_part+this%p(1,pp)**2 &
                + this%p(2,pp)**2 + this%p(3,pp)**2)

      pp = pp + 1
    end do

    pp = ptrcur
    do i=1, np
      ! temp_vec = initial velocity
      temp_vec(1,i) = this%p(1,pp) * rgam(i)
      temp_vec(2,i) = this%p(2,pp) * rgam(i)
      temp_vec(3,i) = this%p(3,pp) * rgam(i)

      bpsq(i) = bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2

      pp = pp + 1
    end do

    ! perform first half-step
    pp = ptrcur
    do i=1, np
      utemp(1,i) = this%p(1,pp) + ep(1,i) + temp_vec(2,i)*bp(3,i)-temp_vec(3,i)*bp(2,i)
      utemp(2,i) = this%p(2,pp) + ep(2,i) + temp_vec(3,i)*bp(1,i)-temp_vec(1,i)*bp(3,i)
      utemp(3,i) = this%p(3,pp) + ep(3,i) + temp_vec(1,i)*bp(2,i)-temp_vec(2,i)*bp(1,i)

      pp = pp + 1
    end do


    do i=1, np
      ! temp_vec = u prime
      temp_vec(1,i) = utemp(1,i) + ep(1,i)
      temp_vec(2,i) = utemp(2,i) + ep(2,i)
      temp_vec(3,i) = utemp(3,i) + ep(3,i)
    end do

    do i=1, np
      ustar(i) = temp_vec(1,i)*bp(1,i) + temp_vec(2,i)*bp(2,i) + temp_vec(3,i)*bp(3,i)
      sigma(i) = 1.0_p_k_part + temp_vec(1,i)**2 + temp_vec(2,i)**2 + temp_vec(3,i)**2  - bpsq(i)
    end do

    do i=1, np
      rgam(i) = 1.0_p_k_part / sqrt( 0.5_p_k_part * ( sigma(i) + sqrt( sigma(i)**2 &
                + 4.0_p_k_part * ( bpsq(i) + ustar(i)**2 ) ) ) )
    end do

    do i=1, np
      ! bp => tau parameter
      bp(1,i) = bp(1,i) * rgam(i)
      bp(2,i) = bp(2,i) * rgam(i)
      bp(3,i) = bp(3,i) * rgam(i)
    end do

    do i=1, np
      spar(i) = 1.0_p_k_part / ( 1.0_p_k_part + bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2)
      uptp(i) = temp_vec(1,i)*bp(1,i) + temp_vec(2,i)*bp(2,i) + temp_vec(3,i)*bp(3,i)
    end do

    ! Store initial momentum for energy calculations
    ! to do: find a way to do this without storing the original momentum
    pp = ptrcur
    do i=1, np
      utemp(1,i) = this%p(1,pp)
      utemp(2,i) = this%p(2,pp)
      utemp(3,i) = this%p(3,pp)

      pp = pp + 1
    enddo

    ! perform second half-step
    pp = ptrcur
    do i=1, np
      this%p(1,pp) = spar(i) * ( temp_vec(1,i) + uptp(i)*bp(1,i) &
                      + temp_vec(2,i)*bp(3,i)-temp_vec(3,i)*bp(2,i) )

      this%p(2,pp) = spar(i) * ( temp_vec(2,i) + uptp(i)*bp(2,i) &
                      + temp_vec(3,i)*bp(1,i)-temp_vec(1,i)*bp(3,i) )

      this%p(3,pp) = spar(i) * ( temp_vec(3,i) + uptp(i)*bp(3,i) &
                      + temp_vec(1,i)*bp(2,i)-temp_vec(2,i)*bp(1,i) )

      pp = pp + 1
    end do


    ! do energy diagnostic
    loc_ene = 0
    pp = ptrcur
    do i=1,np

      ! get time centered velocities
      uc(1) = 0.5_p_k_part * ( utemp(1,i) + this%p(1,pp) )
      uc(2) = 0.5_p_k_part * ( utemp(2,i) + this%p(2,pp) )
      uc(3) = 0.5_p_k_part * ( utemp(3,i) + this%p(3,pp) )

      uc2 = uc(1)**2 + uc(2)**2 + uc(3)**2
      gamma_c = sqrt( uc2 + 1 )

      loc_ene = loc_ene + this%q(pp) * uc2 / ( gamma_c + 1 )

      pp = pp+1
    enddo

    ! accumulate global energy
    ! Only accumulating after every bunch has better roundoff properties
    energy = energy + loc_ene

  enddo

end subroutine dudt_vay
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dudt_cond_vay( this, emf, dt, i0, i1, energy, t )
!-----------------------------------------------------------------------------------------
! Advance velocities and calculate time-centered total energy
! process particles in the range [i0,i1]
! using VAY pusher [Ref. xxx] only for particles with gamma > 5
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species ),    intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  real(p_double),                 intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_double), intent(inout) :: energy
  real(p_double), intent(in) :: t

  real(p_k_part) :: tem, gamma_c
  real(p_k_part), dimension(p_p_dim) :: uc
  real(p_k_part) :: uc2
  integer :: i, j, ptrcur, np, pp, nsimd, nvay
  integer, dimension(2,p_cache_size) :: psimd, pvay ! particle indices for each pusher
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp, temp_vec
  real(p_k_part), dimension(p_cache_size) :: rgam, ustar, sigma, spar, uptp, bpsq
  real(p_k_part), dimension(p_cache_size) ::  gam_tem, otsq
  real(p_k_part) :: gamma, u2
  real( p_double ) :: loc_ene

  ! executable statements

  ! get factor that includes timestep and charge-to-mass ratio
  ! note that the charge-to-mass ratio is used since the
  ! momentum p is not the total momentum of a particles but
  ! the momentum per unit restmass (which is the electron mass)

  tem = real( 0.5_p_double * dt / this%rqm, p_k_part )

  !loop through all particles
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! do things common to both pushers
    call this % get_emf( emf, bp, ep, np, ptrcur, t=t )

    do i=1, np
      ep(1,i) = ep(1,i) * tem
      ep(2,i) = ep(2,i) * tem
      ep(3,i) = ep(3,i) * tem
    end do

    loc_ene = 0

    ! calculate gamma for each particle, populate index arrays for each pusher
    pp = ptrcur
    nsimd = 0; nvay = 0;
    do i=1, np
      gam_tem(i) = sqrt(1.0_p_k_part+this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2)

      if ( gam_tem(i) > 5.0_p_k_part ) then
        nvay = nvay + 1
        pvay(1,nvay) = i
        if (nvay==1) then
          pvay(2,nvay) = i - 1
        else
          pvay(2,nvay) = pvay(1,nvay) - pvay(1,nvay-1)
        endif

      else
        nsimd = nsimd + 1
        psimd(1,nsimd) = i
        if (nsimd==1) then
          psimd(2,nsimd) = i - 1
        else
          psimd(2,nsimd) = psimd(1,nsimd) - psimd(1,nsimd-1)
        endif

      endif

      pp = pp + 1
    enddo

    !-----------------------------------
    ! vay pusher
    pp = ptrcur
    do j=1, nvay
      i = pvay(1,j)
      pp = pp + pvay(2,j)

      bp(1,i) = bp(1,i) * tem
      bp(2,i) = bp(2,i) * tem
      bp(3,i) = bp(3,i) * tem

      rgam(i) = 1.0_p_k_part / gam_tem(i)
    end do

    pp = ptrcur
    do j=1, nvay
      i = pvay(1,j)
      pp = pp + pvay(2,j)

      ! temp_vec = initial velocity
      temp_vec(1,i) = this%p(1,pp) * rgam(i)
      temp_vec(2,i) = this%p(2,pp) * rgam(i)
      temp_vec(3,i) = this%p(3,pp) * rgam(i)

      bpsq(i) = bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2
    end do

    ! perform first half-step
    pp = ptrcur
    do j=1, nvay
      i = pvay(1,j)
      pp = pp + pvay(2,j)

      utemp(1,i) = this%p(1,pp) + ep(1,i) + temp_vec(2,i)*bp(3,i)-temp_vec(3,i)*bp(2,i)
      utemp(2,i) = this%p(2,pp) + ep(2,i) + temp_vec(3,i)*bp(1,i)-temp_vec(1,i)*bp(3,i)
      utemp(3,i) = this%p(3,pp) + ep(3,i) + temp_vec(1,i)*bp(2,i)-temp_vec(2,i)*bp(1,i)
    end do


    do j=1, nvay
      i = pvay(1,j)
      ! temp_vec = u prime
      temp_vec(1,i) = utemp(1,i) + ep(1,i)
      temp_vec(2,i) = utemp(2,i) + ep(2,i)
      temp_vec(3,i) = utemp(3,i) + ep(3,i)
    end do

    do j=1, nvay
      i = pvay(1,j)
      ustar(i) = temp_vec(1,i)*bp(1,i) + temp_vec(2,i)*bp(2,i) + temp_vec(3,i)*bp(3,i)
      sigma(i) = 1.0_p_k_part + temp_vec(1,i)**2 + temp_vec(2,i)**2 + temp_vec(3,i)**2  - bpsq(i)
    end do

    do j=1, nvay
      i = pvay(1,j)
      rgam(i) = 1.0_p_k_part / sqrt( 0.5_p_k_part * ( sigma(i) + sqrt( sigma(i)**2 &
                + 4.0_p_k_part * ( bpsq(i) + ustar(i)**2 ) ) ) )
    end do

    do j=1, nvay
      i = pvay(1,j)
      ! bp => tau parameter
      bp(1,i) = bp(1,i) * rgam(i)
      bp(2,i) = bp(2,i) * rgam(i)
      bp(3,i) = bp(3,i) * rgam(i)
    end do

    do j=1, nvay
      i = pvay(1,j)
      spar(i) = 1.0_p_k_part / ( 1.0_p_k_part + bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2)
      uptp(i) = temp_vec(1,i)*bp(1,i) + temp_vec(2,i)*bp(2,i) + temp_vec(3,i)*bp(3,i)
    end do

    ! Store initial momentum for energy calculations
    ! to do: find a way to do this without storing the original momentum
    pp = ptrcur
    do j=1, nvay
      i = pvay(1,j)
      pp = pp + pvay(2,j)

      utemp(1,i) = this%p(1,pp)
      utemp(2,i) = this%p(2,pp)
      utemp(3,i) = this%p(3,pp)
    enddo

    ! perform second half-step
    pp = ptrcur
    do j=1, nvay
      i = pvay(1,j)
      pp = pp + pvay(2,j)

      this%p(1,pp) = spar(i) * ( temp_vec(1,i) + uptp(i)*bp(1,i) &
                      + temp_vec(2,i)*bp(3,i)-temp_vec(3,i)*bp(2,i) )

      this%p(2,pp) = spar(i) * ( temp_vec(2,i) + uptp(i)*bp(2,i) &
                      + temp_vec(3,i)*bp(1,i)-temp_vec(1,i)*bp(3,i) )

      this%p(3,pp) = spar(i) * ( temp_vec(3,i) + uptp(i)*bp(3,i) &
                      + temp_vec(1,i)*bp(2,i)-temp_vec(2,i)*bp(1,i) )
    end do


    ! do energy diagnostic
    pp = ptrcur
    do j=1, nvay
      i = pvay(1,j)
      pp = pp + pvay(2,j)

      ! get time centered velocities
      uc(1) = 0.5_p_k_part * ( utemp(1,i) + this%p(1,pp) )
      uc(2) = 0.5_p_k_part * ( utemp(2,i) + this%p(2,pp) )
      uc(3) = 0.5_p_k_part * ( utemp(3,i) + this%p(3,pp) )

      uc2 = uc(1)**2 + uc(2)**2 + uc(3)**2
      gamma_c = sqrt( uc2 + 1 )

      loc_ene = loc_ene + this%q(pp) * uc2 / ( gamma_c + 1 )
    enddo



    !-----------------------------------
    ! regular simd pusher
    pp = ptrcur
    do j=1,nsimd
      i = psimd(1,j)
      pp = pp + psimd(2,j)

      utemp(1,i) = this%p(1,pp) + ep(1,i)
      utemp(2,i) = this%p(2,pp) + ep(2,i)
      utemp(3,i) = this%p(3,pp) + ep(3,i)

      ! Get time centered gamma
      u2 =  (utemp(1,i)**2 + utemp(2,i)**2) + utemp(3,i)**2
      gamma = sqrt( u2 + 1 )

      ! accumulate time centered energy
      ! this is done in double precicion always
      loc_ene = loc_ene + this%q(pp) * u2 / (gamma + 1.0_p_double)

      gam_tem(i)= tem / gamma
    enddo

    do j=1,nsimd
      i = psimd(1,j)
      bp(1,i) = bp(1,i)*gam_tem(i)
      bp(2,i) = bp(2,i)*gam_tem(i)
      bp(3,i) = bp(3,i)*gam_tem(i)
    end do

    pp = ptrcur
    do j=1,nsimd
      i = psimd(1,j)
      pp = pp + psimd(2,j)
      this%p(1,pp) = utemp(1,i) + utemp(2,i) * bp(3,i)
      this%p(2,pp) = utemp(2,i) + utemp(3,i) * bp(1,i)
      this%p(3,pp) = utemp(3,i) + utemp(1,i) * bp(2,i)
    end do

    pp = ptrcur
    do j=1,nsimd
      i = psimd(1,j)
      pp = pp + psimd(2,j)
      this%p(1,pp) = this%p(1,pp) - utemp(3,i) * bp(2,i)
      this%p(2,pp) = this%p(2,pp) - utemp(1,i) * bp(3,i)
      this%p(3,pp) = this%p(3,pp) - utemp(2,i) * bp(1,i)
    end do

    do j=1,nsimd
      i = psimd(1,j)
      otsq(i) = 2.0_p_k_part / ( ((1.0_p_k_part + bp(1,i)**2) + bp(2,i)**2) + bp(3,i)**2)
    end do

    do j=1,nsimd
      i = psimd(1,j)
      bp(1,i) = bp(1,i) * otsq(i)
      bp(2,i) = bp(2,i) * otsq(i)
      bp(3,i) = bp(3,i) * otsq(i)
    end do

    pp = ptrcur
    do j=1,nsimd
      i = psimd(1,j)
      pp = pp + psimd(2,j)
      utemp(1,i) = utemp(1,i) + this%p(2,pp) * bp(3,i)
      utemp(2,i) = utemp(2,i) + this%p(3,pp) * bp(1,i)
      utemp(3,i) = utemp(3,i) + this%p(1,pp) * bp(2,i)
    end do

    pp = ptrcur
    do j=1,nsimd
      i = psimd(1,j)
      pp = pp + psimd(2,j)
      utemp(1,i) = utemp(1,i) - this%p(3,pp) * bp(2,i)
      utemp(2,i) = utemp(2,i) - this%p(1,pp) * bp(3,i)
      utemp(3,i) = utemp(3,i) - this%p(2,pp) * bp(1,i)
    end do

    ! Perform second half of electric field acceleration.
    pp = ptrcur
    do j=1,nsimd
      i = psimd(1,j)
      pp = pp + psimd(2,j)
      this%p(1,pp) = utemp(1,i) + ep(1,i)
      this%p(2,pp) = utemp(2,i) + ep(2,i)
      this%p(3,pp) = utemp(3,i) + ep(3,i)
    end do

    ! accumulate global energy
    ! Only accumulating after every bunch has better roundoff properties
    energy = energy + loc_ene

  enddo

end subroutine dudt_cond_vay
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dudt_cary( this, emf, dt, i0, i1, energy, t )
!-----------------------------------------------------------------------------------------
! Advance velocities and calculate time-centered total energy
! process particles in the range [i0,i1]
! This routine is from A. V. Higuera and J. R. Cary, Phys. Plasmas 24, 052104 (2017)
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species ),    intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  real(p_double),                 intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_double), intent(inout) :: energy
  real(p_double), intent(in) :: t

  real(p_k_part) :: tem
  integer :: i, ptrcur, np, pp
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp
  real(p_k_part), dimension(p_cache_size) ::  gamma, otsq
  real(p_k_part) :: u2, gam_minus_sq, bpsq, bdotusq
  real(p_double) :: loc_ene

  ! executable statements
  ! print *, '(*warn*) Using dudt_cary'

  ! get factor that includes timestep and charge-to-mass ratio
  ! note that the charge-to-mass ratio is used since the
  ! momentum p is not the total momentum of a particles but
  ! the momentum per unit restmass (which is the electron mass)

  tem = real( 0.5_p_double * (dt / this%rqm), p_k_part )

  !loop through all particles
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! modify bp & ep to include timestep and charge-to-mass ratio
    ! and perform half the electric field acceleration.
    ! Result is stored in UTEMP.

    call this % get_emf( emf, bp, ep, np, ptrcur, t=t )

    do i=1, np
      ep(1,i) = ep(1,i) * tem
      ep(2,i) = ep(2,i) * tem
      ep(3,i) = ep(3,i) * tem

      bp(1,i) = bp(1,i) * tem
      bp(2,i) = bp(2,i) * tem
      bp(3,i) = bp(3,i) * tem
    end do

    loc_ene = 0

    pp = ptrcur
    do i=1,np
      utemp(1,i) = this%p(1,pp) + ep(1,i)
      utemp(2,i) = this%p(2,pp) + ep(2,i)
      utemp(3,i) = this%p(3,pp) + ep(3,i)

      ! Get time centered gamma
      u2 =  (utemp(1,i)**2 + utemp(2,i)**2) + utemp(3,i)**2
      gam_minus_sq = u2 + 1.0_p_double

      ! accumulate time centered energy
      ! this is done in double precicion always
      ! note, this is not done with unew and gammanew from Cary paper due to complexity
      loc_ene = loc_ene + this%q(pp) * u2 / (sqrt(gam_minus_sq) + 1.0_p_double)

      ! The Cary pusher only differs in how gamma is computed
      bpsq = bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2
      bdotusq = bp(1,i)*utemp(1,i) + bp(2,i)*utemp(2,i) + bp(3,i)*utemp(3,i)
      bdotusq = bdotusq**2

      gamma(i)= sqrt( 0.5_p_k_part * ( gam_minus_sq - bpsq + &
        sqrt( ( gam_minus_sq - bpsq )**2 + 4.0_p_k_part *( bpsq + bdotusq ) ) ) )

      pp = pp + 1
    enddo

    ! accumulate global energy
    energy = energy + loc_ene

    do i=1,np
      bp(1,i) = bp(1,i) / gamma(i)
      bp(2,i) = bp(2,i) / gamma(i)
      bp(3,i) = bp(3,i) / gamma(i)
    end do

    pp = ptrcur
    do i=1,np
      this%p(1,pp) = utemp(1,i) + utemp(2,i) * bp(3,i)
      this%p(2,pp) = utemp(2,i) + utemp(3,i) * bp(1,i)
      this%p(3,pp) = utemp(3,i) + utemp(1,i) * bp(2,i)
      pp = pp + 1
    end do

    pp = ptrcur
    do i=1,np
      this%p(1,pp) = this%p(1,pp) - utemp(3,i) * bp(2,i)
      this%p(2,pp) = this%p(2,pp) - utemp(1,i) * bp(3,i)
      this%p(3,pp) = this%p(3,pp) - utemp(2,i) * bp(1,i)
      pp = pp + 1
    end do

    do i=1,np
      otsq(i) = 2.0_p_k_part / ( ((1.0_p_k_part + bp(1,i)**2) + bp(2,i)**2) + bp(3,i)**2)
    end do

    do i=1,np
      bp(1,i) = bp(1,i) * otsq(i)
      bp(2,i) = bp(2,i) * otsq(i)
      bp(3,i) = bp(3,i) * otsq(i)
    end do

    pp = ptrcur
    do i=1,np
      utemp(1,i) = utemp(1,i) + this%p(2,pp) * bp(3,i)
      utemp(2,i) = utemp(2,i) + this%p(3,pp) * bp(1,i)
      utemp(3,i) = utemp(3,i) + this%p(1,pp) * bp(2,i)
      pp = pp + 1
    end do

    pp = ptrcur
    do i=1,np
      utemp(1,i) = utemp(1,i) - this%p(3,pp) * bp(2,i)
      utemp(2,i) = utemp(2,i) - this%p(1,pp) * bp(3,i)
      utemp(3,i) = utemp(3,i) - this%p(2,pp) * bp(1,i)
      pp = pp + 1
    end do

    ! Perform second half of electric field acceleration.
    pp = ptrcur
    do i=1,np
      this%p(1,pp) = utemp(1,i) + ep(1,i)
      this%p(2,pp) = utemp(2,i) + ep(2,i)
      this%p(3,pp) = utemp(3,i) + ep(3,i)
      pp = pp + 1
    end do

  enddo

end subroutine dudt_cary
!-----------------------------------------------------------------------------------------

#ifdef __HAS_SPIN__

!-----------------------------------------------------------------------------------------
subroutine dsdt_vieira( this, ep, bp, p_old, gamma, dt, ptrcur, chunk_size )
!-----------------------------------------------------------------------------------------
! Advance spin
! process particles in the range [i0,i1]
! This routine is from J. Vieira et al, PRAB 14, 071303 (2011)
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species ), intent(inout) :: this
  real(p_k_part), intent(in), dimension(:,:) :: bp, ep, p_old
  real(p_k_part), intent(in), dimension(:) :: gamma
  real(p_double), intent(in) :: dt
  integer, intent(in) :: ptrcur, chunk_size

  integer :: i, pp
  real(p_k_part), dimension(p_p_dim, chunk_size) :: omega, vtemp, stemp
  real(p_k_part) :: a, coef, vdotb, igam

  a = this%anom_mag_moment

  pp = ptrcur
  do i = 1, chunk_size

    igam = 1.0_p_k_part / gamma(i)

    ! calculate the time-centered velocity
    vtemp(1,i) = 0.5_p_k_part * ( p_old(1,i) + this%p(1,pp) ) * igam
    vtemp(2,i) = 0.5_p_k_part * ( p_old(2,i) + this%p(2,pp) ) * igam
    vtemp(3,i) = 0.5_p_k_part * ( p_old(3,i) + this%p(3,pp) ) * igam

    ! Now calculate the precession frequency Omega
    ! calculate contribution from B
    coef = a + 1.0_p_k_part * igam
    omega(1,i) = coef * bp(1,i) * gamma(i)
    omega(2,i) = coef * bp(2,i) * gamma(i)
    omega(3,i) = coef * bp(3,i) * gamma(i)

    ! calculate (v cross E) contribution
    coef = -1.0_p_k_part * ( a + 1.0_p_k_part / ( 1.0_p_k_part + gamma(i) ) )
    omega(1,i) = omega(1,i) + coef * ( vtemp(2,i) * ep(3,i) - vtemp(3,i) * ep(2,i) )
    omega(2,i) = omega(2,i) + coef * ( vtemp(3,i) * ep(1,i) - vtemp(1,i) * ep(3,i) )
    omega(3,i) = omega(3,i) + coef * ( vtemp(1,i) * ep(2,i) - vtemp(2,i) * ep(1,i) )

    ! calculate (v dot B) contribution
    ! note that omega should be multiplied by 2 because ep and bp are already
    ! divided by 2 in dudt. But rotation vector is obtained via d = omega * dt/2,
    ! therefore the factor 2 is canceled.
    vdotb = vtemp(1,i) * bp(1,i) + vtemp(2,i) * bp(2,i) + vtemp(3,i) * bp(3,i)
    coef = -1.0_p_k_part * ( a * gamma(i)**2 / ( 1.0_p_k_part + gamma(i) ) * vdotb )
    omega(1,i) = omega(1,i) + coef * vtemp(1,i)
    omega(2,i) = omega(2,i) + coef * vtemp(2,i)
    omega(3,i) = omega(3,i) + coef * vtemp(3,i)

    ! calculate s_prime
    stemp(1,i) = this%s(1,pp) + ( this%s(2,pp) * omega(3,i) - this%s(3,pp) * omega(2,i) )
    stemp(2,i) = this%s(2,pp) + ( this%s(3,pp) * omega(1,i) - this%s(1,pp) * omega(3,i) )
    stemp(3,i) = this%s(3,pp) + ( this%s(1,pp) * omega(2,i) - this%s(2,pp) * omega(1,i) )

    coef = 2.0_p_k_part / ( 1.0_p_k_part + omega(1,i)**2 + omega(2,i)**2 + omega(3,i)**2 )
    this%s(1,pp) = this%s(1,pp) + coef * ( stemp(2,i) * omega(3,i) - stemp(3,i) * omega(2,i) )
    this%s(2,pp) = this%s(2,pp) + coef * ( stemp(3,i) * omega(1,i) - stemp(1,i) * omega(3,i) )
    this%s(3,pp) = this%s(3,pp) + coef * ( stemp(1,i) * omega(2,i) - stemp(2,i) * omega(1,i) )

    pp = pp + 1
  enddo

end subroutine dsdt_vieira
!-----------------------------------------------------------------------------------------

#endif

!-----------------------------------------------------------------------------------------
function ntrim(x)
!-----------------------------------------------------------------------------------------
! Returns the integer shift (-1, 0 or +1) so that the coordinate remains in the [-0.5, 0.5[
! range. This is the fastest implementation (twice as fast as a sequence of ifs) because
! the two if structures compile as conditional moves and can be processed independently.
! This has no precision problem and is only 12% slower than the previous "int(x+1.5)-1"
! routine that would break for x = nearest( 0.5, -1.0 )
!-----------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine advance_deposit_1d( this, emf, jay, energy, gdt, i0, i1, t )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 1

  class( t_species ), intent(inout) :: this
  class( t_emf ), intent(in)  ::  emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: i0, i1
  real(p_double), intent(in) :: t

  integer :: np, ptrcur, pp, i
  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rgamma
  real(p_k_part) :: dt_dx1

  dt_dx1 = real( gdt/this%dx(1), p_k_part )

  ! update momenta
  if ( .not. this%free_stream ) call this%dudt( emf, gdt, i0, i1, energy, t )

  ! loop through all particles
  do ptrcur = i0, i1, p_cache_size

     ! check if last copy of table and set np
     if( ptrcur + p_cache_size > i1 ) then
         np = i1 - ptrcur + 1
     else
         np = p_cache_size
     endif

     pp = ptrcur
     do i=1,np

        rgamma(i) = 1.0_p_k_part / &
          sqrt( ( ( 1.0_p_k_part + this%p(1,pp)**2 ) + this%p(2,pp)**2 ) + this%p(3,pp)**2 )

        pp = pp + 1
     end do

     pp = ptrcur
     do i=1,np
       xbuf (1,i) = this%x (1,pp) + ( this%p(1,pp) * rgamma(i) ) * dt_dx1
       dxi(1,i) = ntrim( xbuf(1,i) )

       pp = pp + 1
     end do

     ! Deposit current
     call this % dep_current_1d( jay, dxi, xbuf, &
        this%ix(:,ptrcur:), this%x(:,ptrcur:), &
        this%q(ptrcur:), rgamma,     &
        this%p(:,ptrcur:),      &
        np, gdt )

     ! copy data from buffer to species data
     pp = ptrcur
     do i = 1, np
       this%x(1,pp)  = xbuf(1,i)     - dxi(1,i)
       this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)

       pp = pp + 1
     end do

  enddo

end subroutine advance_deposit_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine advance_deposit_2d( this, emf, jay, energy, gdt, i0, i1, t )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  class( t_species ),    intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  type( t_vdf ),        intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_double),   intent(in) :: gdt
  integer, intent(in) :: i0, i1
  real(p_double), intent(in) :: t

  integer :: i, pp, np, ptrcur
  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rgamma
  real(p_k_part) :: dt_dx1, dt_dx2
  
  dt_dx1 = real( gdt/this%dx(1), p_k_part )
  dt_dx2 = real( gdt/this%dx(2), p_k_part )

  ! update momenta
  if ( .not. this%free_stream ) call this%dudt( emf, gdt, i0, i1, energy, t )

  ! advance position of particles i0 to i1 in chunks of p_cache_size
  do ptrcur = i0, i1, p_cache_size

     ! check if last copy of table and set np
     if( ptrcur + p_cache_size > i1 ) then
         np = i1 - ptrcur + 1
     else
         np = p_cache_size
     endif

     ! this type of loop is actually faster than
     ! using a forall construct
     pp = ptrcur

     do i=1,np

        rgamma(i) = 1.0_p_k_part / &
          sqrt( ( ( 1.0_p_k_part + this%p(1,pp)**2 ) + this%p(2,pp)**2 ) + this%p(3,pp)**2 )

        pp = pp + 1
     end do

     pp = ptrcur
     do i=1,np
       xbuf(1,i) = this%x(1,pp) + ( this%p(1,pp) * rgamma(i) ) * dt_dx1
       xbuf(2,i) = this%x(2,pp) + ( this%p(2,pp) * rgamma(i) ) * dt_dx2

       dxi(1,i) = ntrim( xbuf(1,i) )
       dxi(2,i) = ntrim( xbuf(2,i) )

       pp = pp + 1
     end do

     ! Deposit current
     call this % dep_current_2d( jay, dxi, xbuf, &
                         this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                         this%q(ptrcur:), rgamma,     &
                         this%p(:,ptrcur:),      &
                         np, gdt )

    ! copy data from buffer to species data trimming positions
    pp = ptrcur
    do i = 1, np
      this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
      this%x(2,pp)  = xbuf(2,i) - dxi(2,i)
      this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
      this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

      pp = pp + 1
    end do

  enddo

end subroutine advance_deposit_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine advance_deposit_3d( this, emf, jay, energy, gdt, i0, i1, t )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 3

  class( t_species ),    intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  type( t_vdf), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_double),                     intent(in) :: gdt
  integer, intent(in) :: i0, i1
  real(p_double), intent(in) :: t

  real(p_k_part), dimension(p_max_dim) :: dt_dx
  integer :: np, ptrcur, pp, i
  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size) :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rgamma

  dt_dx(1) = real( gdt/this%dx(1), p_k_part )
  dt_dx(2) = real( gdt/this%dx(2), p_k_part )
  dt_dx(3) = real( gdt/this%dx(3), p_k_part )

  ! update momenta
  if ( .not. this%free_stream ) call this%dudt( emf, gdt, i0, i1, energy, t )

  ! loop through all the particles
  do ptrcur = i0, i1, p_cache_size

     ! check if last copy of table and set np
     if( ptrcur + p_cache_size > i1 ) then
         np = i1 - ptrcur + 1
     else
         np = p_cache_size
     endif

     ! calculate 1/gamma
     pp = ptrcur
     do i=1, np
        rgamma(i) = 1.0_p_k_part / &
          sqrt( ( ( 1.0_p_k_part + this%p(1,pp)**2 ) + this%p(2,pp)**2 ) + this%p(3,pp)**2 )

        pp = pp + 1
     end do

     ! advance particle position
     pp = ptrcur
     do i=1,np

       xbuf(1,i) = this%x(1,pp) + (this%p(1,pp) * rgamma(i)) * dt_dx(1)
       xbuf(2,i) = this%x(2,pp) + (this%p(2,pp) * rgamma(i)) * dt_dx(2)
       xbuf(3,i) = this%x(3,pp) + (this%p(3,pp) * rgamma(i)) * dt_dx(3)

       dxi(1,i) = ntrim( xbuf(1,i) )
       dxi(2,i) = ntrim( xbuf(2,i) )
       dxi(3,i) = ntrim( xbuf(3,i) )

       pp = pp + 1

     end do

     ! Deposit current
     call this % dep_current_3d( jay, dxi, xbuf, &
                       this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                       this%q(ptrcur:), np, gdt )

    ! copy data from buffer to species data trimming positions
     pp = ptrcur
     do i = 1, np
       this%x(1,pp)  = xbuf(1,i)  - dxi(1,i)
       this%x(2,pp)  = xbuf(2,i)  - dxi(2,i)
       this%x(3,pp)  = xbuf(3,i)  - dxi(3,i)
       this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
       this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)
       this%ix(3,pp) = this%ix(3,pp) + dxi(3,i)

       pp = pp + 1
     end do

  enddo

end subroutine advance_deposit_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine advance_deposit_cyl_2d( this, emf, jay, energy, gdt, i0, i1, t )

  implicit none

  integer, parameter :: rank = 2

  class( t_species ), intent(inout) :: this
  class( t_emf ), intent( in )  :: emf
  type( t_vdf ), intent(inout) :: jay
  real(p_double), intent(inout) :: energy
  real(p_double), intent(in) :: gdt
  integer, intent(in) :: i0, i1
  real(p_double), intent(in) :: t

  integer :: i, pp, np, ptrcur
  real(p_double), dimension(p_x_dim)   :: xmin_g
  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rgamma

  real(p_k_part) :: dt_dx1
  integer :: gix2, shift_ix2
  real(p_double) :: dr, rdr, tmp

  real(p_double) :: x2_new, x3_new, r_old, r_new
  real(p_double) :: r_shift

  ! executable statements

  xmin_g(1) = this%g_box( p_lower, 1 )
  xmin_g(2) = this%g_box( p_lower, 2 )

  shift_ix2 = this%my_nx_p(p_lower, 2) - 2
  dr  = this%dx(p_r_dim)
  rdr = 1.0_p_double/dr

  dt_dx1 = real( gdt / this%dx(1), p_k_part )

  ! advance momenta using EM fields
  ! the particle momenta will be further changed below
  if ( .not. this%free_stream ) call this%dudt( emf, gdt, i0, i1, energy, t )

  ! Since r = 0 is at the center of cell ix2 = 1 the radial position of particles will be
  ! - Positions defined with regard to the center of the cell (odd interpolation)
  !       r = ( (gix2-1) + x2 ) * dr
  ! - Positions defined with regard to the corner of the cell (even interpolation)
  !       r = ( (gix2-1) + x2 - 0.5 ) * dr

  if ( this%pos_type == p_cell_near ) then
    r_shift = 0.5_p_double
  else
    r_shift = 0.0_p_double
  endif

  ! advance position
  do ptrcur = i0, i1, p_cache_size

     ! check if last copy of table and set np
     if( ptrcur + p_cache_size > i1 ) then
         np = i1 - ptrcur + 1
     else
         np = p_cache_size
     endif

     ! this type of loop is actually faster than
     ! using a forall construct
     pp = ptrcur

     do i=1,np
        rgamma(i) = 1.0_p_k_part / &
          sqrt( ( ( 1.0_p_k_part + this%p(1,pp)**2 ) + this%p(2,pp)**2 ) + this%p(3,pp)**2 )

        pp = pp + 1
     end do

     ! advance particle position
     pp = ptrcur

     do i=1,np
       xbuf(1,i) = this%x(1,pp) + ( this%p(1,pp) * rgamma(i) ) * dt_dx1

       ! Convert radial "cell" position to "box" position in double precision
       gix2  = this%ix(2,pp)  + shift_ix2
       r_old = ( this%x(2,pp) + (gix2 - r_shift) ) * dr

       ! Particle is pushed in 3D like fashion
       x2_new    = r_old + ( this%p(2,pp) * rgamma(i) ) * gdt
       x3_new    =         ( this%p(3,pp) * rgamma(i) ) * gdt

       r_new     = sqrt( x2_new**2 + x3_new**2 )

       ! Convert new position to "cell" type position

       ! this is a protection against roundoff for cold plasmas
       if ( r_old == r_new ) then
         xbuf(2,i) = this%x(2,pp)
       else
         xbuf(2,i) = real( r_new * rdr - ( gix2 - r_shift ) , p_k_part)
       endif

       ! Correct p_r and p_\theta to conserve angular momentum
       ! there is a potential division by zero here
       tmp          = 1.0_p_double / r_new
       this%p(2,pp) = real( ( this%p(2,pp) * x2_new + this%p(3,pp)*x3_new ) * tmp, p_k_part )
       this%p(3,pp) = real( this%p(3,pp) * ( r_old * tmp ), p_k_part )

       dxi(1,i) = ntrim( xbuf(1,i) )
       dxi(2,i) = ntrim( xbuf(2,i) )

       pp = pp + 1
     end do

     ! Deposit current
     call this % dep_current_2d( jay, dxi, xbuf, &
                         this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                         this%q(ptrcur:), rgamma,     &
                         this%p(:,ptrcur:),      &
                         np, gdt )

    ! copy data from buffer to species data trimming positions
    pp = ptrcur
    do i = 1, np
      this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
      this%x(2,pp)  = xbuf(2,i) - dxi(2,i)
      this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
      this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

      pp = pp + 1
    end do

  enddo


end subroutine advance_deposit_cyl_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Push particles and deposit electric current
!-----------------------------------------------------------------------------------------
subroutine advance_deposit( this, emf, current, t, tstep, tid, n_threads )

  use m_time_step

  implicit none

  class( t_species ), intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  type( t_vdf ), intent(inout) :: current

  real(p_double), intent(in) :: t
  type( t_time_step ) :: tstep

  integer, intent(in) :: tid        ! local thread id
  integer, intent(in) :: n_threads  ! total number of threads

  ! local variables

  real(p_double) :: dtcycle, energy
  integer :: chunk, i0, i1

  ! executable statements

  ! call this%validate( "before advance deposit" )

  ! if before push start time return silently
  if ( t < this%push_start_time ) return

  ! initialize time centered energy diagnostic
  energy = 0.0_p_double

  dtcycle = real( dt(tstep), p_k_part )

  ! range of particles for each thread
  chunk = ( this%num_par + n_threads - 1 ) / n_threads
  i0    = tid * chunk + 1
  i1    = min( (tid+1) * chunk, this%num_par )

  ! Push particles. Boundary crossings will be checked at update_boundary
  select case ( this%coordinates )

    case default
      select case ( p_x_dim )
      case (1)
        call advance_deposit_1d( this, emf, current, energy, dtcycle, i0, i1, t )
      case (2)
        call advance_deposit_2d( this, emf, current, energy, dtcycle, i0, i1, t )
      case (3)
        call advance_deposit_3d( this, emf, current, energy, dtcycle, i0, i1, t )
      case default
        ERROR('Not implemented for x_dim = ',p_x_dim)
        call abort_program(p_err_invalid)
      end select
    case ( p_cylindrical_b )
      call advance_deposit_cyl_2d( this, emf, current, energy, dtcycle, i0, i1, t )
  end select

  this%energy(tid+1) = energy


  ! call this%validate( "after advance deposit", over=.true. )

end subroutine advance_deposit
!-----------------------------------------------------------------------------------------


end module m_species_push
