!-----------------------------------------------------------------------------------------
! Particle accelerator (for beam initialization)
!-----------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_accelerate

#include "memory/memory.h"

use m_parameters

use m_species_define
use m_species_current

use m_vdf_define

use m_emf_define
use m_emf
use m_emf_interpolate

use m_math

implicit none

private

interface accelerate_deposit
  module procedure accelerate_deposit
end interface

public :: accelerate_deposit

contains

!-----------------------------------------------------------------------------------------
! Accelerate particles using a fixed acceleration, advance positions only in the x1
! direction and deposit current
!-----------------------------------------------------------------------------------------
subroutine accelerate_deposit( this, jay, t, tstep, tid, n_threads )

  use m_time_step

  implicit none

  class( t_species ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: jay

  real(p_double), intent(in) :: t
  type( t_time_step ) :: tstep

  integer, intent(in) :: tid		! local thread id
  integer, intent(in) :: n_threads  ! total number of threads

  ! local variables

  real(p_k_part) :: dt_dx1, u_accel, u_frac, q_frac
  integer :: chunk, i0, i1


  ! Time centered energy diagnostic is not yet available, force recalculation later
  this%energy = p_ene_recalc

  dt_dx1 = real( dt(tstep)/this%dx(1), p_k_part )

  ! range of particles for each thread
  chunk = ( this%num_par + n_threads - 1 ) / n_threads
  i0    = tid * chunk + 1
  i1    = min( (tid+1) * chunk, this%num_par )

  ! ufl momentum acceleration in p1
  if ( this%udist%n_accelerate > 0 ) then
    ! calculate charge fraction range=(0., 1.]
    u_frac = 1.0
    select case ( this%udist%n_accelerate_type )
    case ( incr_linear )
      u_frac = real( grad_linear( n(tstep)+1, this%udist%n_accelerate ), p_k_part )

    case ( incr_regressive )
      u_frac = real( grad_regressive( n(tstep)+1, this%udist%n_accelerate ), p_k_part )

    case default
      ERROR('Not implemented')
      call abort_program( p_err_notimplemented )
    end select

  else
    u_frac = 1.0
  endif

  ! increase charge over n steps
  if ( this%udist%n_q_incr > 0 ) then
    ! calculate charge fraction range=(0., 1.]
    select case ( this%udist%n_q_incr_type )
    case ( incr_linear )
      q_frac = real( grad_linear( n(tstep)+1, this%udist%n_q_incr ), p_k_part )

    case ( incr_regressive )
      q_frac = real( grad_regressive( n(tstep)+1, this%udist%n_q_incr ), p_k_part )

    case default
      ERROR('Not implemented')
      call abort_program( p_err_notimplemented )
    end select
  else
    q_frac = 1.0
  endif

  ! Push particles. Boundary crossings will be checked at update_boundary
  select case ( p_x_dim )

  case (1)
    call accelerate_deposit_1d( this, jay, dt(tstep), i0, i1, u_frac, q_frac )

  case (2)
    call accelerate_deposit_2d( this, jay, dt(tstep), i0, i1, u_frac, q_frac )

  case (3)
    call accelerate_deposit_3d( this, jay, dt(tstep), i0, i1, u_frac, q_frac )

  end select


end subroutine accelerate_deposit
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine accelerate_deposit_1d( this, jay, dt, i0, i1, u_frac, q_frac )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 1

  class( t_species ), intent(inout) :: this

  type( t_vdf ),     intent(inout) :: jay
  real( p_double ),     intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_k_part), intent(in) :: u_frac, q_frac
  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rg
  real(p_k_part), dimension(p_cache_size)      :: q_reduced
  real(p_k_part) :: vdt_dx1, dt_dx1, u_accel
  integer :: i, ptrcur, np, pp

  dt_dx1 = real( dt/this%dx(1), p_k_part )
  u_accel = u_frac * this%udist%ufl(1)
  vdt_dx1 = u_accel * dt_dx1 / &
                    sqrt( 1.0_p_k_part + u_accel**2 )
  ! This is used to cancel currents in other directions
  do i = 1, p_cache_size
    rg(i) = 0.0
  end do

  ! Move particles only in the accel. direction
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    pp = ptrcur
    if(this%udist%use_particle_uacc) then
      do i=1,np

        u_accel = u_frac * this%p(1,pp)

        vdt_dx1 = u_accel * dt_dx1 / &
                          sqrt( 1.0_p_k_part + u_accel**2 )

        xbuf(1,i) = this%x(1,pp) + vdt_dx1

        dxi(1,i) = ntrim( xbuf(1,i) )

        pp = pp + 1
      end do
    else
      do i=1,np

        xbuf(1,i) = this%x(1,pp) + vdt_dx1

        dxi(1,i) = ntrim( xbuf(1,i) )

        pp = pp + 1
      end do
    endif

    q_reduced(1:np) = this%q(ptrcur:ptrcur+np-1) * q_frac

    ! Deposit current
    call this % dep_current_1d( jay, dxi, xbuf, &
                this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                q_reduced, rg,     &
                this%p(:,ptrcur:),      &
                np, dt )

    ! copy data from buffer to species data trimming positions
    pp = ptrcur
    do i = 1, np
      this%x(1,pp)  = xbuf(1,i)     - dxi(1,i)
      this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
      pp = pp + 1
    end do

  end do

end subroutine accelerate_deposit_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine accelerate_deposit_2d( this, jay, dt, i0, i1, u_frac, q_frac )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  class( t_species ), intent(inout) :: this

  type( t_vdf ),     intent(inout) :: jay
  real( p_double ),     intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_k_part), intent(in) :: u_frac, q_frac

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: rg
  real(p_k_part), dimension(p_cache_size)      :: q_reduced
  real(p_k_part) :: vdt_dx1, dt_dx1, u_accel
  integer :: i, ptrcur, np, pp

  dt_dx1 = real( dt/this%dx(1), p_k_part )
  u_accel = u_frac * this%udist%ufl(1)
  vdt_dx1 = u_accel * dt_dx1 / &
                    sqrt( 1.0_p_k_part + u_accel**2 )

  ! This is used to cancel currents in other directions
  do i = 1, p_cache_size
    rg(i) = 0.0
    dxi(2,i) = 0
  end do

  ! Move particles only in the accel. direction
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    pp = ptrcur
    if(this%udist%use_particle_uacc) then
      do i=1,np
        u_accel = u_frac * this%p(1,pp)
        vdt_dx1 = u_accel * dt_dx1 / &
                          sqrt( 1.0_p_k_part + u_accel**2 )
        xbuf(1,i) = this%x(1,pp) + vdt_dx1
        xbuf(2,i) = this%x(2,pp)

        dxi(1,i) = ntrim( xbuf(1,i) )

        pp = pp + 1
      end do
    else
      do i=1,np

        xbuf(1,i) = this%x(1,pp) + vdt_dx1
        xbuf(2,i) = this%x(2,pp)

        dxi(1,i) = ntrim( xbuf(1,i) )

        pp = pp + 1
      end do
    endif

    q_reduced(1:np) = this%q(ptrcur:ptrcur+np-1) * q_frac

    ! Deposit current
    call this % dep_current_2d( jay, dxi, xbuf, &
                this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                q_reduced, rg,     &
                this%p(:,ptrcur:),      &
                np, dt )

    ! copy data from buffer to species data trimming positions
    pp = ptrcur
    do i = 1, np
      this%x(1,pp)  = xbuf(1,i)     - dxi(1,i)
      this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
      pp = pp + 1
    end do

  end do

end subroutine accelerate_deposit_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine accelerate_deposit_3d( this, jay, dt, i0, i1, u_frac, q_frac )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 3

  class( t_species ), intent(inout) :: this

  type( t_vdf ),     intent(inout) :: jay
  real( p_double ),     intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_k_part), intent(in) :: u_frac, q_frac

  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size)        :: dxi
  real(p_k_part), dimension(p_cache_size)      :: q_reduced
  real(p_k_part) :: vdt_dx1, dt_dx1, u_accel
  integer :: i, ptrcur, np, pp


  dt_dx1 = real( dt/this%dx(1), p_k_part )
  u_accel = u_frac * this%udist%ufl(1)
  vdt_dx1 = u_accel * dt_dx1 / &
                    sqrt( 1.0_p_k_part + u_accel**2 )
  ! This is used to cancel currents in other directions
  do i = 1, p_cache_size
    dxi(2,i) = 0
    dxi(3,i) = 0
  enddo

  ! Move particles only in the accel. direction
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    pp = ptrcur
    if(this%udist%use_particle_uacc) then
      do i=1,np
        u_accel = u_frac * this%p(1,pp)

        vdt_dx1 = u_accel * dt_dx1 / &
                          sqrt( 1.0_p_k_part + u_accel**2 )

        xbuf(1,i) = this%x(1,pp) + vdt_dx1
        xbuf(2,i) = this%x(2,pp)
        xbuf(3,i) = this%x(3,pp)

        dxi(1,i) = ntrim( xbuf(1,i) )

        pp = pp + 1
      end do
    else
      do i=1,np

        xbuf(1,i) = this%x(1,pp) + vdt_dx1
        xbuf(2,i) = this%x(2,pp)
        xbuf(3,i) = this%x(3,pp)

        dxi(1,i) = ntrim( xbuf(1,i) )

        pp = pp + 1
      end do
    endif


    q_reduced(1:np) = this%q(ptrcur:ptrcur+np-1) * q_frac

    ! Deposit current
    call this % dep_current_3d( jay, dxi, xbuf, &
                this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                q_reduced, np, dt )

    ! copy data from buffer to species data trimming positions
    pp = ptrcur
    do i = 1, np
      this%x(1,pp)  = xbuf(1,i)     - dxi(1,i)
      this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
      pp = pp + 1
    end do

  end do

end subroutine accelerate_deposit_3d
!-----------------------------------------------------------------------------------------

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


end module m_species_accelerate
