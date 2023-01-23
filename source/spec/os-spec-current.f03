! m_species_current module
!
! Handles current deposition
!

! -------------------------------------------------------------------------------------------------
! The new current deposition routines expect the old cell index (ixold) and position (xold), the
! new position (xnew) still indexed to the old cell index, and the number of cells moved (dxi)
! which can be -1, 0, or 1.
!
! This has the advantage of removing roundoff problems in the current splitters where the total
! particle motion is required. Since there is no trimming of xnew before the current split the particle motion will never be
! 0 (which could happen before due to roundoff) when there is a cell edge crossing.
! -------------------------------------------------------------------------------------------------

! if __TEMPLATE__ is defined then read the template definition at the end of the file
#ifndef __TEMPLATE__


#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_current

use m_species_define

use m_system
use m_parameters

use m_vdf_define

use m_node_conf
use m_space

use m_utilities

use m_random

implicit none

! restrict access to things explicitly declared public
private

type :: t_vp1D
  real(p_k_part) :: x0
  real(p_k_part) :: x1
  real(p_k_part) :: q
  real(p_k_part) :: vy, vz
  integer        :: i
end type t_vp1D

type :: t_vp2D
!  sequence
  real(p_k_part) :: x0, y0
  real(p_k_part) :: x1, y1
  real(p_k_part) :: q
  real(p_k_part) :: vz
  integer        :: i, j
end type t_vp2D

type :: t_vp3D
  real(p_k_part) :: x0, y0, z0
  real(p_k_part) :: x1, y1, z1
  real(p_k_part) :: q
  integer        :: i, j, k
end type t_vp3D

public :: t_vp1D, t_vp2D, t_vp3D

! Turns on the SIMD like splitting routines
! (We should test if these versions outperform the normal ones.)
!
! Define the VEC_WIDTH macro to the appropriate value to enable these routines.
! These are the values that should be used for compatibility with the SIMD code.
!
!            single  double
! --------------------------
! SSE             4       2
! QPX             4       4
! AVX/AVX2        8       4
! AVX-512        16       8
!
! Please leave this commented out unless you know what you are doing
!#define VEC_WIDTH 4


! cell based positions

interface dep_current_1d_s1
   module procedure dep_current_1d_s1
end interface

interface dep_current_1d_s2
   module procedure dep_current_1d_s2
end interface

interface dep_current_1d_s3
   module procedure dep_current_1d_s3
end interface

interface dep_current_1d_s4
   module procedure dep_current_1d_s4
end interface

interface dep_current_2d_s1
   module procedure dep_current_2d_s1
end interface

interface dep_current_2d_s2
   module procedure dep_current_2d_s2
end interface

interface dep_current_2d_s3
   module procedure dep_current_2d_s3
end interface

interface getjr_2d_s4
   module procedure dep_current_2d_s4
end interface

interface dep_current_3d_s1
   module procedure dep_current_3d_s1
end interface

interface dep_current_3d_s2
   module procedure dep_current_3d_s2
end interface

interface dep_current_3d_s3
   module procedure dep_current_3d_s3
end interface

interface dep_current_3d_s4
   module procedure dep_current_3d_s4
end interface

interface deposit_current_1d
  module procedure deposit_current_1d_cell
end interface

interface deposit_current_2d
  module procedure deposit_current_2d_cell
end interface

interface deposit_current_3d
  module procedure deposit_current_3d_cell
end interface


! declare things that should be public

public :: deposit_current_1d, deposit_current_2d, deposit_current_3d

public :: dep_current_1d_s1, dep_current_1d_s2, dep_current_1d_s3, dep_current_1d_s4
public :: dep_current_2d_s1, dep_current_2d_s2, dep_current_2d_s3, dep_current_2d_s4
public :: dep_current_3d_s1, dep_current_3d_s2, dep_current_3d_s3, dep_current_3d_s4

! Auxiliary buffers for current deposition from boundary conditions / cathodes

real(p_k_part), pointer, dimension(:,:)   :: tmp_xold => null(), tmp_xnew => null()
real(p_k_part), pointer, dimension(:,:)   :: tmp_p => null()
real(p_k_part), pointer, dimension(:)     :: tmp_rg => null(), tmp_q => null()
integer, pointer, dimension(:,:)          :: tmp_dxi => null(), tmp_ix => null()

public :: tmp_xold, tmp_xnew, tmp_p, tmp_rg, tmp_q, tmp_dxi, tmp_ix

contains

!-------------------------------------------------------------------------------------------------
! Deposit current wrapper functions
!-------------------------------------------------------------------------------------------------


! -----------------------------------------------------------------------------
subroutine deposit_current_1d_cell( species, jay, dxi, xnew, ixold, xold, &
                                    q, rgamma, u, np, dt )
! -----------------------------------------------------------------------------

  implicit none

  ! dummy variables

  class( t_species ), intent(inout) :: species
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(inout) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(inout) :: xnew, xold
  real(p_k_part), dimension( : ), intent(inout) :: q, rgamma
  real(p_k_part), dimension(:,:), intent(inout) :: u
  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt

  ! local variables

  integer :: idx, lnp

  ! executable statements
  do idx = 1, np, p_cache_size

     if (idx+p_cache_size <= np) then
       lnp = p_cache_size
     else
       lnp = np - idx + 1
     endif

     select case ( species%interpolation )

        case (p_linear)

          call dep_current_1d_s1( jay, dxi(:,idx:), xnew(:,idx:), &
                                     ixold(:,idx:), xold(:,idx:), &
                                     q(idx:), rgamma(idx:), &
                                     u(:,idx:), lnp, dt )

        case (p_quadratic)

          call dep_current_1d_s2( jay, dxi(:,idx:), xnew(:,idx:), &
                                     ixold(:,idx:), xold(:,idx:), &
                                     q(idx:), rgamma(idx:), &
                                     u(:,idx:), lnp, dt )

        case (p_cubic)

          call dep_current_1d_s3( jay, dxi(:,idx:), xnew(:,idx:), &
                                     ixold(:,idx:), xold(:,idx:), &
                                     q(idx:), rgamma(idx:), &
                                     u(:,idx:), lnp, dt )

        case (p_quartic)

          call dep_current_1d_s4( jay, dxi(:,idx:), xnew(:,idx:), &
                                     ixold(:,idx:), xold(:,idx:), &
                                     q(idx:), rgamma(idx:), &
                                     u(:,idx:), lnp, dt )
     end select

  enddo

end subroutine deposit_current_1d_cell
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine deposit_current_2d_cell( species, jay, dxi, xnew, ixold, xold, &
                                    q, rgamma, u, np, dt )
! -----------------------------------------------------------------------------

  implicit none

  ! dummy variables

  class( t_species ), intent(inout) :: species
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(inout) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(inout) :: xnew, xold
  real(p_k_part), dimension( : ), intent(inout) :: q, rgamma
  real(p_k_part), dimension(:,:), intent(inout) :: u
  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt

  ! local variables

  integer :: idx, lnp

  ! executable statements
  do idx = 1, np, p_cache_size

     if (idx+p_cache_size <= np) then
       lnp = p_cache_size
     else
       lnp = np - idx + 1
     endif

     select case ( species%interpolation )

        case (p_linear)

           call dep_current_2d_s1( jay, dxi(:,idx:), xnew(:,idx:), &
                                      ixold(:,idx:), xold(:,idx:), &
                                      q(idx:), rgamma(idx:), &
                                      u(:,idx:), lnp, dt )

        case (p_quadratic)

           call dep_current_2d_s2( jay, dxi(:,idx:), xnew(:,idx:), &
                                      ixold(:,idx:), xold(:,idx:), &
                                      q(idx:), rgamma(idx:), &
                                      u(:,idx:), lnp, dt )

        case (p_cubic)

           call dep_current_2d_s3( jay, dxi(:,idx:), xnew(:,idx:), &
                                      ixold(:,idx:), xold(:,idx:), &
                                      q(idx:), rgamma(idx:), &
                                      u(:,idx:), lnp, dt )

        case (p_quartic)

           call dep_current_2d_s4( jay, dxi(:,idx:), xnew(:,idx:), &
                                      ixold(:,idx:), xold(:,idx:), &
                                      q(idx:), rgamma(idx:), &
                                      u(:,idx:), lnp, dt )

     end select
  enddo

end subroutine deposit_current_2d_cell
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine deposit_current_3d_cell( species, jay, dxi, xnew, ixold, xold, q, np, dt )
! -----------------------------------------------------------------------------

  implicit none

  ! dummy variables

  class( t_species ), intent(inout) :: species
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(inout) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(inout) :: xnew, xold
  real(p_k_part), dimension( : ), intent(inout) :: q

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt

  ! local variables

  integer :: idx, lnp

  ! executable statements


  ! particles must be processed in bunches of up to p_cache_size

  do idx = 1, np, p_cache_size

     if (idx+p_cache_size <= np) then
       lnp = p_cache_size
     else
       lnp = np - idx + 1
     endif

     select case ( species%interpolation )

        case (p_linear)

           call dep_current_3d_s1( jay, dxi(:,idx:), xnew(:,idx:), &
                                      ixold(:,idx:), xold(:,idx:), &
                                      q(idx:), lnp, dt )

        case (p_quadratic)

           call dep_current_3d_s2( jay, dxi(:,idx:), xnew(:,idx:), &
                                      ixold(:,idx:), xold(:,idx:), &
                                      q(idx:), lnp, dt )

        case (p_cubic)

           call dep_current_3d_s3( jay, dxi(:,idx:), xnew(:,idx:), &
                                      ixold(:,idx:), xold(:,idx:), &
                                      q(idx:), lnp, dt )

        case (p_quartic)

           call dep_current_3d_s4( jay, dxi(:,idx:), xnew(:,idx:), &
                                      ixold(:,idx:), xold(:,idx:), &
                                      q(idx:), lnp, dt )

     end select

  end do


end subroutine deposit_current_3d_cell
! -----------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! Longitudinal current weights
! ----------------------------------------------------------------------------------------

subroutine wl_s1( qnx, x0, x1, wl )

  implicit none

  real( p_k_fld ), intent(in) :: qnx, x0, x1
  real( p_k_fld ), intent(inout), dimension( 0:0 ) :: wl

  wl( 0) = qnx * ( x1 - x0 )

end subroutine wl_s1

subroutine wl_s2( qnx, x0, x1, wl )

  implicit none

  real( p_k_fld ), intent(in) :: qnx, x0, x1
  real( p_k_fld ), intent(inout), dimension( -1:0 ) :: wl

  real( p_k_fld ) :: d, s1_2, p1_2, n

  d     = x1 - x0
  s1_2  = 0.5_p_k_fld - x0
  p1_2  = 0.5_p_k_fld + x0

  wl(-1) =  s1_2 - 0.5_p_k_fld * d
  wl( 0) =  p1_2 + 0.5_p_k_fld * d

  n = qnx * d
  wl(-1) = n * wl(-1)
  wl( 0) = n * wl( 0)

end subroutine wl_s2

subroutine wl_s3( qnx, x0, x1, wl )

  implicit none

  real( p_k_fld ), intent(in) :: qnx, x0, x1
  real( p_k_fld ), intent(inout), dimension( -1:1 ) :: wl

  real( p_k_fld ) :: d, s1_2, p1_2, d_3, n

  real( p_k_fld ), parameter :: c1_3 = 1.0_p_k_fld / 3.0_p_k_fld

  d     =  x1 - x0
  s1_2  = 0.5_p_k_fld - x0
  p1_2  = 0.5_p_k_fld + x0
  d_3   = c1_3 * d

  wl(-1) =  0.5_p_k_fld * ( s1_2**2 - d * ( s1_2 - d_3 ) )
  wl( 0) =  ( 0.75_p_k_fld - x0**2 ) - d * ( x0 + d_3 )
  wl(+1) =  0.5_p_k_fld * ( p1_2**2 + d * ( p1_2 + d_3 ) )

  n = qnx * d
  wl(-1) = n * wl(-1)
  wl( 0) = n * wl( 0)
  wl(+1) = n * wl(+1)

end subroutine wl_s3

subroutine wl_s4( qnx, x0, x1, wl )

  implicit none

  real( p_k_fld ), intent(in) :: qnx, x0, x1
  real( p_k_fld ), intent(inout), dimension( -2:1 ) :: wl

  real( p_k_fld ) :: d, s, p, t, u, n
  real( p_k_fld ) :: s2, s3, p2, p3

  real( p_k_fld ), parameter :: c1_6 = 1.0_p_k_fld / 6.0_p_k_fld
  real( p_k_fld ), parameter :: c1_3 = 1.0_p_k_fld / 3.0_p_k_fld
  real( p_k_fld ), parameter :: c2_3 = 2.0_p_k_fld / 3.0_p_k_fld

  d = x1 - x0
  s  = 0.5_p_k_fld - x0
  p  = 0.5_p_k_fld + x0
  t  = x0 - c1_6
  u  = x0 + c1_6

  s2 = s * s;   s3 = s2 * s
  p2 = p * p;   p3 = p2 * p

  wl(-2) =  c1_6 * ( d * ( d * ( s - 0.25_p_k_fld*d ) - 1.5_p_k_fld * s2 ) + s3 )
  wl(-1) =    d * ( ( 0.5_p_k_fld * d ) * (t + 0.25_p_k_fld*d) + ( 0.75_p_k_fld*(t*t) - c1_3 ) ) + &
            ( 0.5_p_k_fld * p3 + ( c2_3 - p2 ) )
  wl( 0) =  - d * ( ( 0.5_p_k_fld * d ) * (u + 0.25_p_k_fld*d) + ( 0.75_p_k_fld*(u*u) - c1_3 ) ) + &
            ( 0.5_p_k_fld * s3 + ( c2_3 - s2 ) )
  wl( 1) =  c1_6 * ( d * ( d * ( p + 0.25_p_k_fld*d ) + 1.5_p_k_fld * p2 ) + p3 )

  n = qnx * d
  wl(-2) = n * wl(-2)
  wl(-1) = n * wl(-1)
  wl( 0) = n * wl( 0)
  wl( 1) = n * wl( 1)

end subroutine wl_s4

! ----------------------------------------------------------------------------------------
! Spiines
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! Linear
! ----------------------------------------------------------------------------------------
subroutine spline_s1( x, s )

  implicit none

  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(0:1), intent(inout) :: s

  s(0) = 0.5_p_k_fld - x
  s(1) = 0.5_p_k_fld + x

end subroutine spline_s1

! ----------------------------------------------------------------------------------------
! Quadratic
! ----------------------------------------------------------------------------------------
subroutine spline_s2( x, s )

  implicit none

  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(-1:1), intent(inout) :: s

  real( p_k_fld ) :: t0, t1

  t0 = 0.5_p_k_fld - x
  t1 = 0.5_p_k_fld + x

  s(-1) = 0.5_p_k_fld * t0**2
  s( 0) = 0.5_p_k_fld + t0*t1
  s( 1) = 0.5_p_k_fld * t1**2

end subroutine spline_s2

! ----------------------------------------------------------------------------------------
! Cubic
! ----------------------------------------------------------------------------------------

subroutine spline_s3( x, s )

  implicit none

  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(-1:2), intent(inout) :: s

  real( p_k_fld ) :: t0, t1, t2, t3
  real( p_k_fld ), parameter :: c1_6 = 1.0_p_k_fld / 6.0_p_k_fld
  real( p_k_fld ), parameter :: c2_3 = 2.0_p_k_fld / 3.0_p_k_fld

  t0 = 0.5_p_k_fld - x
  t1 = 0.5_p_k_fld + x

  t2 = t0 * t0
  t3 = t1 * t1

  t0 = t0 * t2
  t1 = t1 * t3

  s(-1) = c1_6 * t0
  s( 0) = ( c2_3 - t3 ) + 0.5_p_k_fld*t1
  s( 1) = ( c2_3 - t2 ) + 0.5_p_k_fld*t0
  s( 2) = c1_6 * t1

end subroutine spline_s3

! ----------------------------------------------------------------------------------------
! Quartic
! ----------------------------------------------------------------------------------------

subroutine spline_s4( x, s )

  implicit none

  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(-2:2), intent(inout) :: s

  real( p_k_fld ) :: t0, t1, t2

  real( p_k_fld ), parameter :: c1_6  = 1.0_p_k_fld / 6.0_p_k_fld
  real( p_k_fld ), parameter :: c1_24 = 1.0_p_k_fld / 24.0_p_k_fld
  real( p_k_fld ), parameter :: c11_24 = 11.0_p_k_fld / 24.0_p_k_fld

  t0 = 0.5_p_k_fld - x
  t1 = 0.5_p_k_fld + x
  t2 = t0 * t1

!  s(-2) = c1_24 * t0**4
  s(-2) = c1_24 * ((t0**2)*(t0**2))
  s(-1) = c1_6 * ( 0.25_p_k_fld + t0 * ( 1.0_p_k_fld + t0 * ( 1.5_p_k_fld + t2 ) ) )
  s( 0) = c11_24 + t2 * ( 0.5_p_k_fld + 0.25_p_k_fld * t2 )
  s( 1) = c1_6 * ( 0.25_p_k_fld + t1 * ( 1.0_p_k_fld + t1 * ( 1.5_p_k_fld + t2 ) ) )
!  s( 2) = c1_24 * t1**4
  s( 2) = c1_24 * ((t1**2)*(t1**2))

end subroutine spline_s4

!-------------------------------------------------------------------------------------------------
! Normal versions of the splitting routines
!-------------------------------------------------------------------------------------------------
#ifndef VEC_WIDTH

!-------------------------------------------------------------------------------------------------
! Splits particle trajectories so that all virtual particles have a motion starting and ending
! in the same cell. This routines also calculates vy, vz for each virtual particle.
!-------------------------------------------------------------------------------------------------
subroutine split_1d( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf1D, nsplit )

  implicit none

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer,                 intent(in) :: np       ! number of particles to process

  type(t_vp1D), dimension(:), intent(out) :: vpbuf1D
  integer,                 intent(out) :: nsplit  ! number of virtual particles created

  ! local variables
  real(p_k_part) :: xint, delta
  real(p_k_part) :: vz, vy
  integer        :: k,i

  ! create virtual particles
  k=0
  do i=1,np
     if (dxi(1,i) == 0) then
        k=k+1
        vpbuf1D(k)%x0 = xold(1,i)
        vpbuf1D(k)%x1 = xnew(1,i)
        vpbuf1D(k)%q = q(i)

        vpbuf1D(k)%vy = u(2,i)*rgamma(i)
        vpbuf1D(k)%vz = u(3,i)*rgamma(i)

        vpbuf1D(k)%i = ixold(1,i)
     else
        xint = 0.5_p_k_part * dxi(1,i)

        vy = u(2,i)*rgamma(i)
        vz = u(3,i)*rgamma(i)
        delta = (xint - xold(1,i)) / (xnew(1,i) - xold(1,i))

        k=k+1
        vpbuf1D(k)%x0 = xold(1,i)
        vpbuf1D(k)%x1 = xint
        vpbuf1D(k)%q = q(i)
        vpbuf1D(k)%vy = vy * delta
        vpbuf1D(k)%vz = vz * delta

        vpbuf1D(k)%i = ixold(1,i)

        k=k+1
        vpbuf1D(k)%x0 = - xint
        vpbuf1D(k)%x1 = xnew(1,i) - dxi(1,i)
        vpbuf1D(k)%q = q(i)
        vpbuf1D(k)%vy = vy * (1-delta)
        vpbuf1D(k)%vz = vz * (1-delta)

        vpbuf1D(k)%i = ixold(1,i) + dxi(1,i)
     endif
  enddo

  nsplit = k

end subroutine split_1d
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
! Splits particle trajectories so that all virtual particles have a motion starting and ending
! in the same cell. This routines also calculates vz for each virtual particle.
!-------------------------------------------------------------------------------------------------
subroutine split_2d( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit )

  implicit none

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer,                 intent(in) :: np       ! number of particles to process

  type(t_vp2D), dimension(:), intent(out) :: vpbuf2D
  integer,                 intent(out) :: nsplit  ! number of virtual particles created

  ! local variables
  real(p_k_part)                 :: xint,yint,xint2,yint2, delta
  real(p_k_part) :: vz, vzint
  integer                        :: k,i

  integer :: cross

  k=0

  ! create virtual particles that correspond to a motion that starts and ends
  ! in the same grid cell

  do i=1,np

    cross = abs(dxi(1,i)) + 2 * abs(dxi(2,i))

    select case( cross )
      case(0) ! no cross
         k=k+1
         vpbuf2D(k)%x0 = xold(1,i)
         vpbuf2D(k)%y0 = xold(2,i)
         vpbuf2D(k)%x1 = xnew(1,i)
         vpbuf2D(k)%y1 = xnew(2,i)
         vpbuf2D(k)%q = q(i)

         vpbuf2D(k)%vz = u(3,i)*rgamma(i)

         vpbuf2D(k)%i = ixold(1,i)
         vpbuf2D(k)%j = ixold(2,i)


      case(1) ! x cross only

         xint = 0.5_p_k_fld * dxi(1,i)
         delta = ( xint - xold(1,i) ) / ( xnew(1,i) - xold(1,i) )
         yint =  xold(2,i) + (xnew(2,i) - xold(2,i)) * delta

         vz = u(3,i)*rgamma(i)

         k=k+1

         vpbuf2D(k)%x0 = xold(1,i)
         vpbuf2D(k)%y0 = xold(2,i)

         vpbuf2D(k)%x1 = xint
         vpbuf2D(k)%y1 = yint

         vpbuf2D(k)%q = q(i)
         vpbuf2D(k)%vz = vz * delta

         vpbuf2D(k)%i = ixold(1,i)
         vpbuf2D(k)%j = ixold(2,i)

         k=k+1
         vpbuf2D(k)%x0 = -xint
         vpbuf2D(k)%y0 = yint

         vpbuf2D(k)%x1 = xnew(1,i) - dxi(1,i)
         vpbuf2D(k)%y1 = xnew(2,i)

         vpbuf2D(k)%q = q(i)
         vpbuf2D(k)%vz = vz * (1-delta)

         vpbuf2D(k)%i = ixold(1,i) + dxi(1,i)
         vpbuf2D(k)%j = ixold(2,i)

      case(2) ! y cross only

         yint = 0.5_p_k_fld * dxi(2,i)
         delta = ( yint - xold(2,i) ) / ( xnew(2,i) - xold(2,i))
         xint =  xold(1,i) + (xnew(1,i) - xold(1,i)) * delta

         vz = u(3,i)*rgamma(i)

         k=k+1
         vpbuf2D(k)%x0 = xold(1,i)
         vpbuf2D(k)%y0 = xold(2,i)

         vpbuf2D(k)%x1 = xint
         vpbuf2D(k)%y1 = yint
         vpbuf2D(k)%q = q(i)
         vpbuf2D(k)%vz = vz * delta

         vpbuf2D(k)%i = ixold(1,i)
         vpbuf2D(k)%j = ixold(2,i)

         k=k+1
         vpbuf2D(k)%x0 = xint
         vpbuf2D(k)%y0 = -yint

         vpbuf2D(k)%x1 = xnew(1,i)
         vpbuf2D(k)%y1 = xnew(2,i) - dxi(2,i)
         vpbuf2D(k)%q = q(i)
         vpbuf2D(k)%vz = vz * (1-delta)

         vpbuf2D(k)%i = ixold(1,i)
         vpbuf2D(k)%j = ixold(2,i) + dxi(2,i)

      case(3) ! x,y cross

         ! split in x direction first
         xint = 0.5_p_k_fld * dxi(1,i)
         delta = ( xint - xold(1,i) ) / ( xnew(1,i) - xold(1,i))
         yint =  xold(2,i) + ( xnew(2,i) - xold(2,i)) * delta


         vz = u(3,i)*rgamma(i)

         ! check if y intersection occured for 1st or 2nd split
         if ((yint >= -0.5_p_k_fld) .and. ( yint < 0.5_p_k_fld )) then

            ! no y cross on 1st vp
            k=k+1
            vpbuf2D(k)%x0 = xold(1,i)
            vpbuf2D(k)%y0 = xold(2,i)

            vpbuf2D(k)%x1 = xint
            vpbuf2D(k)%y1 = yint

            vpbuf2D(k)%q = q(i)
            vpbuf2D(k)%vz = vz * delta

            vzint = vz*(1-delta)

            vpbuf2D(k)%i = ixold(1,i)
            vpbuf2D(k)%j = ixold(2,i)

            ! y split 2nd vp
            k=k+1

            yint2 = 0.5_p_k_fld * dxi(2,i)
            delta = ( yint2 - yint ) / ( xnew(2,i) - yint )
            xint2 = -xint + ( xnew(1,i) - xint ) * delta


            vpbuf2D(k)%x0 = -xint
            vpbuf2D(k)%y0 =  yint

            vpbuf2D(k)%x1 = xint2
            vpbuf2D(k)%y1 = yint2

            vpbuf2D(k)%q = q(i)
            vpbuf2D(k)%vz = vzint * delta

            vpbuf2D(k)%i = ixold(1,i) + dxi(1,i)
            vpbuf2D(k)%j = ixold(2,i)

            k=k+1

            vpbuf2D(k)%x0 = xint2
            vpbuf2D(k)%y0 = -yint2

            vpbuf2D(k)%x1 = xnew(1,i) - dxi(1,i)
            vpbuf2D(k)%y1 = xnew(2,i) - dxi(2,i)
            vpbuf2D(k)%q = q(i)
            vpbuf2D(k)%vz = vzint * (1-delta)

            vpbuf2D(k)%i = ixold(1,i) + dxi(1,i)
            vpbuf2D(k)%j = ixold(2,i) + dxi(2,i)

         else

            vzint = vz * delta

            ! y split 1st vp
            yint2 = 0.5_p_k_fld * dxi(2,i)
            delta = ( yint2 - xold(2,i) ) / ( yint - xold(2,i))
            xint2 = xold(1,i) + (xint - xold(1,i)) * delta

            k=k+1
            vpbuf2D(k)%x0 = xold(1,i)
            vpbuf2D(k)%y0 = xold(2,i)

            vpbuf2D(k)%x1 = xint2
            vpbuf2D(k)%y1 = yint2

            vpbuf2D(k)%q = q(i)
            vpbuf2D(k)%vz = vzint*delta

            vpbuf2D(k)%i = ixold(1,i)
            vpbuf2D(k)%j = ixold(2,i)


            k=k+1
            vpbuf2D(k)%x0 = xint2
            vpbuf2D(k)%y0 = -yint2

            vpbuf2D(k)%x1 = xint
            vpbuf2D(k)%y1 = yint - dxi(2,i)

            vpbuf2D(k)%q = q(i)
            vpbuf2D(k)%vz = vzint*(1-delta)

            vpbuf2D(k)%i = ixold(1,i)
            vpbuf2D(k)%j = ixold(2,i) + dxi(2,i)

            ! no y cross on second vp
            k=k+1
            vpbuf2D(k)%x0 = - xint
            vpbuf2D(k)%y0 = yint - dxi(2,i)
            vpbuf2D(k)%x1 = xnew(1,i) - dxi(1,i)
            vpbuf2D(k)%y1 = xnew(2,i) - dxi(2,i)
            vpbuf2D(k)%q = q(i)
            vpbuf2D(k)%vz = vz - vzint

            vpbuf2D(k)%i = ixold(1,i) + dxi(1,i)
            vpbuf2D(k)%j = ixold(2,i) + dxi(2,i)

         endif

    end select

  enddo

  nsplit = k

end subroutine split_2d
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
! Splits particle trajectories so that all virtual particles have a motion starting and ending
! in the same cell.
!-------------------------------------------------------------------------------------------------
subroutine split_3d( dxi, xnew, ixold, xold, q, np, vpbuf3D, nsplit )

  implicit none

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold
  real(p_k_part), dimension( : ), intent(in) :: q

  integer,                 intent(in) :: np       ! number of particles to process

  type(t_vp3D), dimension(:), intent(out) :: vpbuf3D
  integer,                 intent(out) :: nsplit  ! number of virtual particles created

  ! local variables
  real(p_k_part) :: xint,yint,zint
  real(p_k_part) :: xint2,yint2,zint2
  real(p_k_part) :: delta
  integer        :: k,i, cross
  integer :: j, vpidx

  k=0

  ! create virtual particles

  do i=1,np

    ! zspi = 0
    ! yspi = 0

    cross    =  abs(dxi(1,i)) + 2 * abs(dxi(2,i)) + 4 * abs(dxi(3,i))

    select case( cross )
      case(0) ! no cross
         k=k+1
         vpbuf3D(k)%x0 = xold(1,i)
         vpbuf3D(k)%y0 = xold(2,i)
         vpbuf3D(k)%z0 = xold(3,i)

         vpbuf3D(k)%x1 = xnew(1,i)
         vpbuf3D(k)%y1 = xnew(2,i)
         vpbuf3D(k)%z1 = xnew(3,i)

         vpbuf3D(k)%q = q(i)

         vpbuf3D(k)%i = ixold(1,i)
         vpbuf3D(k)%j = ixold(2,i)
         vpbuf3D(k)%k = ixold(3,i)


      case(1) ! x cross only

         xint = 0.5_p_k_fld * dxi(1,i)
         delta = ( xint - xold(1,i) ) / (xnew(1,i) - xold(1,i))
         yint =  xold(2,i) + (xnew(2,i) - xold(2,i)) * delta
         zint =  xold(3,i) + (xnew(3,i) - xold(3,i)) * delta

         k=k+1

         vpbuf3D(k)%x0 = xold(1,i)
         vpbuf3D(k)%y0 = xold(2,i)
         vpbuf3D(k)%z0 = xold(3,i)

         vpbuf3D(k)%x1 = xint
         vpbuf3D(k)%y1 = yint
         vpbuf3D(k)%z1 = zint

         vpbuf3D(k)%q = q(i)

         vpbuf3D(k)%i = ixold(1,i)
         vpbuf3D(k)%j = ixold(2,i)
         vpbuf3D(k)%k = ixold(3,i)

         k=k+1
         vpbuf3D(k)%x0 = -xint
         vpbuf3D(k)%y0 = yint
         vpbuf3D(k)%z0 = zint
         vpbuf3D(k)%x1 = xnew(1,i) - dxi(1,i)
         vpbuf3D(k)%y1 = xnew(2,i)
         vpbuf3D(k)%z1 = xnew(3,i)

         vpbuf3D(k)%q = q(i)

         vpbuf3D(k)%i = ixold(1,i) + dxi(1,i)
         vpbuf3D(k)%j = ixold(2,i)
         vpbuf3D(k)%k = ixold(3,i)

      case(2) ! y cross only

         yint = 0.5_p_k_fld * dxi(2,i)
         delta = (yint - xold(2,i)) / (xnew(2,i) - xold(2,i))
         xint =  xold(1,i) + (xnew(1,i) - xold(1,i)) * delta
         zint =  xold(3,i) + (xnew(3,i) - xold(3,i)) * delta

         k=k+1

         vpbuf3D(k)%x0 = xold(1,i)
         vpbuf3D(k)%y0 = xold(2,i)
         vpbuf3D(k)%z0 = xold(3,i)

         vpbuf3D(k)%x1 = xint
         vpbuf3D(k)%y1 = yint
         vpbuf3D(k)%z1 = zint

         vpbuf3D(k)%q = q(i)

         vpbuf3D(k)%i = ixold(1,i)
         vpbuf3D(k)%j = ixold(2,i)
         vpbuf3D(k)%k = ixold(3,i)

         k=k+1
         vpbuf3D(k)%x0 = xint
         vpbuf3D(k)%y0 = -yint
         vpbuf3D(k)%z0 = zint
         vpbuf3D(k)%x1 = xnew(1,i)
         vpbuf3D(k)%y1 = xnew(2,i) - dxi(2,i)
         vpbuf3D(k)%z1 = xnew(3,i)

         vpbuf3D(k)%q = q(i)

         vpbuf3D(k)%i = ixold(1,i)
         vpbuf3D(k)%j = ixold(2,i) + dxi(2,i)
         vpbuf3D(k)%k = ixold(3,i)


      case(3) ! x,y cross

         xint = 0.5_p_k_fld * dxi(1,i)
         delta = ( xint - xold(1,i) ) / (xnew(1,i) - xold(1,i))
         yint = xold(2,i) + (xnew(2,i) - xold(2,i)) * delta
         zint = xold(3,i) + (xnew(3,i) - xold(3,i)) * delta ! no z cross

         if ((yint >= -0.5_p_k_fld) .and. ( yint < 0.5_p_k_fld )) then

            ! yspi = 1

            ! no y cross on 1st vp
            k=k+1
            vpbuf3D(k)%x0 = xold(1,i)
            vpbuf3D(k)%y0 = xold(2,i)
            vpbuf3D(k)%z0 = xold(3,i)

            vpbuf3D(k)%x1 = xint
            vpbuf3D(k)%y1 = yint
            vpbuf3D(k)%z1 = zint

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i)
            vpbuf3D(k)%j = ixold(2,i)
            vpbuf3D(k)%k = ixold(3,i)

               ! y split 2nd vp
               yint2 = 0.5_p_k_fld * dxi(2,i)

            delta = ( yint2 - yint ) / (xnew(2,i) - yint)
            xint2 =  -xint + (xnew(1,i) - xint) * delta
            zint2 =   zint + (xnew(3,i) - zint) * delta

            k=k+1

            vpbuf3D(k)%x0 = -xint
            vpbuf3D(k)%y0 = yint
            vpbuf3D(k)%z0 = zint

            vpbuf3D(k)%x1 = xint2
            vpbuf3D(k)%y1 = yint2
            vpbuf3D(k)%z1 = zint2

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i) + dxi(1,i)
            vpbuf3D(k)%j = ixold(2,i)
            vpbuf3D(k)%k = ixold(3,i)

            k=k+1
            vpbuf3D(k)%x0 = xint2
            vpbuf3D(k)%y0 = -yint2
            vpbuf3D(k)%z0 = zint2

            vpbuf3D(k)%x1 = xnew(1,i) - dxi(1,i)
            vpbuf3D(k)%y1 = xnew(2,i) - dxi(2,i)
            vpbuf3D(k)%z1 = xnew(3,i)

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i) + dxi(1,i)
            vpbuf3D(k)%j = ixold(2,i) + dxi(2,i)
            vpbuf3D(k)%k = ixold(3,i)

         else

            ! yspi = 2

            ! y split 1st vp
            yint2 = 0.5_p_k_fld * dxi(2,i)
            delta = (yint2 - xold(2,i)) / (yint - xold(2,i))
            xint2 =  xold(1,i) + ( xint - xold(1,i)) * delta
            zint2 =  xold(3,i) + ( zint - xold(3,i)) * delta

            k=k+1

            vpbuf3D(k)%x0 = xold(1,i)
            vpbuf3D(k)%y0 = xold(2,i)
            vpbuf3D(k)%z0 = xold(3,i)

            vpbuf3D(k)%x1 = xint2
            vpbuf3D(k)%y1 = yint2
            vpbuf3D(k)%z1 = zint2

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i)
            vpbuf3D(k)%j = ixold(2,i)
            vpbuf3D(k)%k = ixold(3,i)

            k=k+1

            vpbuf3D(k)%x0 = xint2
            vpbuf3D(k)%y0 = -yint2
            vpbuf3D(k)%z0 = zint2

            vpbuf3D(k)%x1 = xint
            vpbuf3D(k)%y1 = yint - dxi(2,i)
            vpbuf3D(k)%z1 = zint

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i)
            vpbuf3D(k)%j = ixold(2,i) + dxi(2,i)
            vpbuf3D(k)%k = ixold(3,i)

            k=k+1

            vpbuf3D(k)%x0 = -xint
            vpbuf3D(k)%y0 = yint - dxi(2,i)
            vpbuf3D(k)%z0 = zint

            vpbuf3D(k)%x1 = xnew(1,i)  - dxi(1,i)
            vpbuf3D(k)%y1 = xnew(2,i)  - dxi(2,i)
            vpbuf3D(k)%z1 = xnew(3,i)

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i) + dxi(1,i)
            vpbuf3D(k)%j = ixold(2,i) + dxi(2,i)
            vpbuf3D(k)%k = ixold(3,i)

         endif

      case(4) ! z cross only

         zint = 0.5_p_k_fld * dxi(3,i)
         delta = (zint - xold(3,i)) / (xnew(3,i) - xold(3,i))
         xint =  xold(1,i) + (xnew(1,i) - xold(1,i)) * delta
         yint =  xold(2,i) + (xnew(2,i) - xold(2,i)) * delta

         k=k+1

         vpbuf3D(k)%x0 = xold(1,i)
         vpbuf3D(k)%y0 = xold(2,i)
         vpbuf3D(k)%z0 = xold(3,i)

         vpbuf3D(k)%x1 = xint
         vpbuf3D(k)%y1 = yint
         vpbuf3D(k)%z1 = zint

         vpbuf3D(k)%q = q(i)

         vpbuf3D(k)%i = ixold(1,i)
         vpbuf3D(k)%j = ixold(2,i)
         vpbuf3D(k)%k = ixold(3,i)

         k=k+1
         vpbuf3D(k)%x0 = xint
         vpbuf3D(k)%y0 = yint
         vpbuf3D(k)%z0 = -zint

         vpbuf3D(k)%x1 = xnew(1,i)
         vpbuf3D(k)%y1 = xnew(2,i)
         vpbuf3D(k)%z1 = xnew(3,i) - dxi(3,i)

         vpbuf3D(k)%q = q(i)

         vpbuf3D(k)%i = ixold(1,i)
         vpbuf3D(k)%j = ixold(2,i)
         vpbuf3D(k)%k = ixold(3,i) + dxi(3,i)


      case (5) ! x, z split - 3 vp

         xint = 0.5_p_k_fld * dxi(1,i)
         delta = ( xint - xold(1,i) ) / (xnew(1,i) - xold(1,i))
         yint = xold(2,i) + (xnew(2,i) - xold(2,i)) * delta
         zint = xold(3,i) + (xnew(3,i) - xold(3,i)) * delta

         zint2 = 0.5_p_k_fld * dxi(3,i)

         if ((zint >= -0.5_p_k_fld) .and. ( zint < 0.5_p_k_fld )) then

            ! zspi = 1

            ! no z cross on 1st vp
            k=k+1
            vpbuf3D(k)%x0 = xold(1,i)
            vpbuf3D(k)%y0 = xold(2,i)
            vpbuf3D(k)%z0 = xold(3,i)

            vpbuf3D(k)%x1 = xint
            vpbuf3D(k)%y1 = yint
            vpbuf3D(k)%z1 = zint

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i)
            vpbuf3D(k)%j = ixold(2,i)
            vpbuf3D(k)%k = ixold(3,i)

               ! z split 2nd vp
            delta = (zint2 - zint) / (xnew(3,i) - zint)
            xint2 =  -xint + (xnew(1,i) - xint) * delta
            yint2 =   yint + (xnew(2,i) - yint) * delta

            k=k+1

            vpbuf3D(k)%x0 = -xint
            vpbuf3D(k)%y0 = yint
            vpbuf3D(k)%z0 = zint

            vpbuf3D(k)%x1 = xint2
            vpbuf3D(k)%y1 = yint2
            vpbuf3D(k)%z1 = zint2

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i) + dxi(1,i)
            vpbuf3D(k)%j = ixold(2,i)
            vpbuf3D(k)%k = ixold(3,i)

            k=k+1
            vpbuf3D(k)%x0 = xint2
            vpbuf3D(k)%y0 = yint2
            vpbuf3D(k)%z0 = -zint2

            vpbuf3D(k)%x1 = xnew(1,i) - dxi(1,i)
            vpbuf3D(k)%y1 = xnew(2,i)
            vpbuf3D(k)%z1 = xnew(3,i) - dxi(3,i)
            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i) + dxi(1,i)
            vpbuf3D(k)%j = ixold(2,i)
            vpbuf3D(k)%k = ixold(3,i) + dxi(3,i)

         else

            ! zspi = 2

            ! z split 1st vp
            delta = (zint2 - xold(3,i)) / (zint - xold(3,i))
            xint2 =  xold(1,i) + ( xint - xold(1,i)) * delta
            yint2 =  xold(2,i) + ( yint - xold(2,i)) * delta

            k=k+1

            vpbuf3D(k)%x0 = xold(1,i)
            vpbuf3D(k)%y0 = xold(2,i)
            vpbuf3D(k)%z0 = xold(3,i)

            vpbuf3D(k)%x1 = xint2
            vpbuf3D(k)%y1 = yint2
            vpbuf3D(k)%z1 = zint2

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i)
            vpbuf3D(k)%j = ixold(2,i)
            vpbuf3D(k)%k = ixold(3,i)

            k=k+1

            vpbuf3D(k)%x0 =  xint2
            vpbuf3D(k)%y0 =  yint2
            vpbuf3D(k)%z0 = -zint2

            vpbuf3D(k)%x1 = xint
            vpbuf3D(k)%y1 = yint
            vpbuf3D(k)%z1 = zint - dxi(3,i)

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i)
            vpbuf3D(k)%j = ixold(2,i)
            vpbuf3D(k)%k = ixold(3,i) + dxi(3,i)

            k=k+1

            vpbuf3D(k)%x0 = -xint
            vpbuf3D(k)%y0 = yint
            vpbuf3D(k)%z0 = zint - dxi(3,i)

            vpbuf3D(k)%x1 = xnew(1,i) - dxi(1,i)
            vpbuf3D(k)%y1 = xnew(2,i)
            vpbuf3D(k)%z1 = xnew(3,i) - dxi(3,i)

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i) + dxi(1,i)
            vpbuf3D(k)%j = ixold(2,i)
            vpbuf3D(k)%k = ixold(3,i) + dxi(3,i)

         endif

       case (6) ! y, z split - 3 vp

         yint  = 0.5_p_k_fld * dxi(2,i)
         delta = (yint - xold(2,i)) / (xnew(2,i) - xold(2,i))
         xint = xold(1,i) + (xnew(1,i) - xold(1,i)) * delta ! no x cross
         zint = xold(3,i) + (xnew(3,i) - xold(3,i)) * delta


         if ((zint >= -0.5_p_k_fld) .and. ( zint < 0.5_p_k_fld )) then

            ! zspi = 1

            ! no z cross on 1st vp
            k=k+1
            vpbuf3D(k)%x0 = xold(1,i)
            vpbuf3D(k)%y0 = xold(2,i)
            vpbuf3D(k)%z0 = xold(3,i)

            vpbuf3D(k)%x1 = xint
            vpbuf3D(k)%y1 = yint
            vpbuf3D(k)%z1 = zint

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i)
            vpbuf3D(k)%j = ixold(2,i)
            vpbuf3D(k)%k = ixold(3,i)

               ! z split 2nd vp
            zint2 = 0.5_p_k_fld * dxi(3,i)
            delta = (zint2 - zint) / (xnew(3,i) - zint)
            xint2 =  xint + (xnew(1,i) - xint) * delta
            yint2 = -yint + (xnew(2,i) - yint) * delta

            k=k+1

            vpbuf3D(k)%x0 = xint
            vpbuf3D(k)%y0 = -yint
            vpbuf3D(k)%z0 = zint

            vpbuf3D(k)%x1 = xint2
            vpbuf3D(k)%y1 = yint2
            vpbuf3D(k)%z1 = zint2

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i)
            vpbuf3D(k)%j = ixold(2,i) + dxi(2,i)
            vpbuf3D(k)%k = ixold(3,i)

            k=k+1
            vpbuf3D(k)%x0 = xint2
            vpbuf3D(k)%y0 = yint2
            vpbuf3D(k)%z0 = -zint2

            vpbuf3D(k)%x1 = xnew(1,i)
            vpbuf3D(k)%y1 = xnew(2,i) - dxi(2,i)
            vpbuf3D(k)%z1 = xnew(3,i) - dxi(3,i)

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i)
            vpbuf3D(k)%j = ixold(2,i) + dxi(2,i)
            vpbuf3D(k)%k = ixold(3,i) + dxi(3,i)

         else

            ! zspi = 2

            ! z split 1st vp
            zint2 = 0.5_p_k_fld * dxi(3,i)
            delta = (zint2 - xold(3,i)) / (zint - xold(3,i))

            xint2 =  xold(1,i) + ( xint - xold(1,i)) * delta
            yint2 =  xold(2,i) + ( yint - xold(2,i)) * delta

            k=k+1

            vpbuf3D(k)%x0 = xold(1,i)
            vpbuf3D(k)%y0 = xold(2,i)
            vpbuf3D(k)%z0 = xold(3,i)

            vpbuf3D(k)%x1 = xint2
            vpbuf3D(k)%y1 = yint2
            vpbuf3D(k)%z1 = zint2

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i)
            vpbuf3D(k)%j = ixold(2,i)
            vpbuf3D(k)%k = ixold(3,i)

            k=k+1

            vpbuf3D(k)%x0 = xint2
            vpbuf3D(k)%y0 = yint2
            vpbuf3D(k)%z0 = -zint2

            vpbuf3D(k)%x1 = xint
            vpbuf3D(k)%y1 = yint
            vpbuf3D(k)%z1 = zint - dxi(3,i)

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i)
            vpbuf3D(k)%j = ixold(2,i)
            vpbuf3D(k)%k = ixold(3,i) + dxi(3,i)

            k=k+1

            vpbuf3D(k)%x0 = xint
            vpbuf3D(k)%y0 = -yint
            vpbuf3D(k)%z0 = zint - dxi(3,i)

            vpbuf3D(k)%x1 = xnew(1,i)
            vpbuf3D(k)%y1 = xnew(2,i) - dxi(2,i)
            vpbuf3D(k)%z1 = xnew(3,i) - dxi(3,i)

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i)
            vpbuf3D(k)%j = ixold(2,i) + dxi(2,i)
            vpbuf3D(k)%k = ixold(3,i) + dxi(3,i)

         endif


      case(7) ! x,y,z cross

         ! do an x,y cross first
         xint = 0.5_p_k_fld * dxi(1,i)
         delta = ( xint - xold(1,i) ) / (xnew(1,i) - xold(1,i))
         yint = xold(2,i) + (xnew(2,i) - xold(2,i)) * delta
         zint = xold(3,i) + (xnew(3,i) - xold(3,i)) * delta ! no z cross

         if ((yint >= -0.5_p_k_fld) .and. ( yint < 0.5_p_k_fld )) then

            ! no y cross on 1st vp
            k=k+1
            vpbuf3D(k)%x0 = xold(1,i)
            vpbuf3D(k)%y0 = xold(2,i)
            vpbuf3D(k)%z0 = xold(3,i)

            vpbuf3D(k)%x1 = xint
            vpbuf3D(k)%y1 = yint
            vpbuf3D(k)%z1 = zint

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i)
            vpbuf3D(k)%j = ixold(2,i)
            vpbuf3D(k)%k = ixold(3,i)

               ! y split 2nd vp
               yint2 = 0.5_p_k_fld * dxi(2,i)

            delta = ( yint2 - yint ) / (xnew(2,i) - yint)
            xint2 =  -xint + (xnew(1,i) - xint) * delta
            zint2 =   zint + (xnew(3,i) - zint) * delta

            k=k+1

            vpbuf3D(k)%x0 = -xint
            vpbuf3D(k)%y0 = yint
            vpbuf3D(k)%z0 = zint

            vpbuf3D(k)%x1 = xint2
            vpbuf3D(k)%y1 = yint2
            vpbuf3D(k)%z1 = zint2

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i) + dxi(1,i)
            vpbuf3D(k)%j = ixold(2,i)
            vpbuf3D(k)%k = ixold(3,i)

            k=k+1
            vpbuf3D(k)%x0 = xint2
            vpbuf3D(k)%y0 = -yint2
            vpbuf3D(k)%z0 = zint2

            vpbuf3D(k)%x1 = xnew(1,i) - dxi(1,i)
            vpbuf3D(k)%y1 = xnew(2,i) - dxi(2,i)
            vpbuf3D(k)%z1 = xnew(3,i)

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i) + dxi(1,i)
            vpbuf3D(k)%j = ixold(2,i) + dxi(2,i)
            vpbuf3D(k)%k = ixold(3,i)

         else

            ! y split 1st vp
            yint2 = 0.5_p_k_fld * dxi(2,i)
            delta = (yint2 - xold(2,i)) / (yint - xold(2,i))
            xint2 =  xold(1,i) + ( xint - xold(1,i)) * delta
            zint2 =  xold(3,i) + ( zint - xold(3,i)) * delta

            k=k+1

            vpbuf3D(k)%x0 = xold(1,i)
            vpbuf3D(k)%y0 = xold(2,i)
            vpbuf3D(k)%z0 = xold(3,i)

            vpbuf3D(k)%x1 = xint2
            vpbuf3D(k)%y1 = yint2
            vpbuf3D(k)%z1 = zint2

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i)
            vpbuf3D(k)%j = ixold(2,i)
            vpbuf3D(k)%k = ixold(3,i)

            k=k+1

            vpbuf3D(k)%x0 = xint2
            vpbuf3D(k)%y0 = -yint2
            vpbuf3D(k)%z0 = zint2

            vpbuf3D(k)%x1 = xint
            vpbuf3D(k)%y1 = yint - dxi(2,i)
            vpbuf3D(k)%z1 = zint

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i)
            vpbuf3D(k)%j = ixold(2,i) + dxi(2,i)
            vpbuf3D(k)%k = ixold(3,i)

            k=k+1

            vpbuf3D(k)%x0 = -xint
            vpbuf3D(k)%y0 = yint - dxi(2,i)
            vpbuf3D(k)%z0 = zint

            vpbuf3D(k)%x1 = xnew(1,i)  - dxi(1,i)
            vpbuf3D(k)%y1 = xnew(2,i)  - dxi(2,i)
            vpbuf3D(k)%z1 = xnew(3,i)

            vpbuf3D(k)%q = q(i)

            vpbuf3D(k)%i = ixold(1,i) + dxi(1,i)
            vpbuf3D(k)%j = ixold(2,i) + dxi(2,i)
            vpbuf3D(k)%k = ixold(3,i)

         endif

         ! one of the 3 vp requires an additional z split
         zint = 0.5_p_k_fld * dxi(3,i)

         do vpidx = k-2, k
           if ( (vpbuf3D(vpidx)%z1 < -0.5_p_k_fld) .or. ( vpbuf3D(vpidx)%z1 >= 0.5_p_k_fld )) then

             delta = (zint - vpbuf3D(vpidx)%z0) / (vpbuf3D(vpidx)%z1 - vpbuf3D(vpidx)%z0)
             xint =  vpbuf3D(vpidx)%x0 + ( vpbuf3D(vpidx)%x1 - vpbuf3D(vpidx)%x0) * delta
             yint =  vpbuf3D(vpidx)%y0 + ( vpbuf3D(vpidx)%y1 - vpbuf3D(vpidx)%y0) * delta

             k = k+1

             ! store new vp
             vpbuf3D(k)%x0 =  xint;             vpbuf3D(k)%x1 = vpbuf3D(vpidx)%x1
             vpbuf3D(k)%y0 =  yint;             vpbuf3D(k)%y1 = vpbuf3D(vpidx)%y1
             vpbuf3D(k)%z0 = -zint;             vpbuf3D(k)%z1 = vpbuf3D(vpidx)%z1 - dxi(3,i)
             vpbuf3D(k)%q  = vpbuf3D(vpidx)%q

             vpbuf3D(k)%i = vpbuf3D(vpidx)%i
             vpbuf3D(k)%j = vpbuf3D(vpidx)%j
             vpbuf3D(k)%k = vpbuf3D(vpidx)%k + dxi(3,i)

             ! correct old vp
             vpbuf3D(vpidx)%x1 = xint
             vpbuf3D(vpidx)%y1 = yint
             vpbuf3D(vpidx)%z1 = zint

             ! correct remaining vp
             do j = vpidx+1, k-1
               vpbuf3D(j)%z0 = vpbuf3D(j)%z0 - dxi(3,i)
               vpbuf3D(j)%z1 = vpbuf3D(j)%z1 - dxi(3,i)
               vpbuf3D(j)%k  = vpbuf3D(j)%k  + dxi(3,i)
             enddo

             exit
           endif
         enddo


    end select

  enddo


  nsplit = k

end subroutine split_3d
!-------------------------------------------------------------------------------------------------

#else

!-------------------------------------------------------------------------------------------------
! SIMD compatibie versions of the splitting routines
!
! These routines will split the trajectories in the same sequence as the SIMD code splitter wouid
! for a vector width of p_vec_width (defined above)
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
! Splits particle trajectories so that all virtual particles have a motion starting and ending
! in the same cell. This routines also calculates vy, vz for each virtual particle.
!-------------------------------------------------------------------------------------------------
subroutine split_1d( dix, x1, ix, x0, q, rgamma, u, np, vpbuf1D, nsplit )

  implicit none

  integer, parameter :: p_vec_width = VEC_WIDTH

  integer, dimension(:,:), intent(inout) :: dix, ix
  real(p_k_part), dimension(:,:), intent(inout) :: x1, x0 ,u
  real(p_k_part), dimension( : ), intent(inout) :: q, rgamma

  integer,                 intent(in) :: np       ! number of particles to process

  type(t_vp1D), dimension(:), intent(inout) :: vpbuf1D
  integer,                 intent(inout) :: nsplit  ! number of virtual particles created

  ! local variables
  real(p_k_part) :: xint, delta
  real(p_k_part) :: vy, vz
  integer        :: i, ip, vp, s, k, v, nptotal

  ! sanity check
  ! If p_cache_size is not a muitipie of p_vec_width we will aiways be padding the particle
  ! buffer with dummy particles. Whiie this does not affect the calculations
  if ( mod( p_cache_size, p_vec_width ) /= 0 ) then
    write(0,*) 'p_cache_size is not a muitipie of p_vec_width(', p_vec_width,')'
  endif

  ! If number of particles to split is not a muitipie of p_vec_width fiii with dummy
  ! particles at the end
  if ( mod( np, p_vec_width ) /= 0 ) then
    nptotal = p_vec_width * ( np / p_vec_width + 1)
    do i = np + 1,  nptotal
      x1(1,i) = 0
      x0(1,i) = 0
      ix(1,i) = 1
      dix(1,i) = 0
      q(i) = 0; rgamma(i) = 1
      u(1,i) = 0; u(2,i) = 0; u(3,i) = 0
    enddo
  else
    nptotal = np
  endif


  ! create virtual particles
  vp = 0
  do ip = 1, nptotal, p_vec_width

    s = vp + p_vec_width

    ! split tracjectories
    do k = 1, p_vec_width

      ! particle index
      i = ip + k - 1

      ! virtual particle index
      v = vp + k

      vy = u(2,i)*rgamma(i)
      vz = u(3,i)*rgamma(i)

      if ( dix(1,i) == 0 ) then
        ! no split required, just copy particle
        vpbuf1D(v)%x0 = x0(1,i)
        vpbuf1D(v)%x1 = x1(1,i)
        vpbuf1D(v)%q  = q(i)

        vpbuf1D(v)%vy = vy
        vpbuf1D(v)%vz = vz

        vpbuf1D(v)%i  = ix(1,i)
      else
        ! crossing cell boundary
        xint = 0.5_p_k_fld * dix(1,i)
        delta = ( xint - x0(1,i) ) / ( x1(1,i) - x0(1,i) )

        ! 1st segment
        vpbuf1D(v)%x0 = x0(1,i); vpbuf1D(v)%x1 = xint
        vpbuf1D(v)%i = ix(1,i)

        vpbuf1D(v)%q = q(i)
        vpbuf1D(v)%vy = vy * delta
        vpbuf1D(v)%vz = vz * delta

        ! add extra particle (2nd segment)
        s = s + 1
        vpbuf1D(s)%x0 = -xint; vpbuf1D(s)%x1 = x1(1,i) - dix(1,i)

        vpbuf1D(s)%i  = ix(1,i) + dix(1,i)

        vpbuf1D(s)%q = q(i)
        vpbuf1D(s)%vy = vy * (1-delta)
        vpbuf1D(s)%vz = vz * (1-delta)
      endif

    enddo

    vp = s

  enddo

  nsplit = vp

end subroutine split_1d
!-------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Splits particle trajectories so that all virtual particles have a motion starting and
! ending in the same cell. This routine also calculates vz for each virtual particle.
!
! This version mimics the splitter used in SIMD code by processing the particles in
! bunches of p_vec_width particles. Setting p_vec_width to 1 will reproduce the originai
! (seriai) splitter.
!
! To avold a iot of dummy (0 charge, 0 motion) particles, p_cache_size should be a muitipie
! of p_vec_width
!-----------------------------------------------------------------------------------------


subroutine split_2d( dix, x1, ix, x0, q, rgamma, u, np, vpbuf2D, nsplit )

  implicit none

  integer, parameter :: p_vec_width = VEC_WIDTH

  integer, dimension(:,:), intent(inout) :: dix, ix
  real(p_k_part), dimension(:,:), intent(inout) :: x1, x0 ,u
  real(p_k_part), dimension( : ), intent(inout) :: q,rgamma

  integer,                 intent(in) :: np       ! number of particles to process

  type(t_vp2D), dimension(:), intent(inout) :: vpbuf2D
  integer,                 intent(inout) :: nsplit  ! number of virtual particles created

  ! local variables
  real(p_k_part) :: xint,yint,xint2,yint2, delta
  real(p_k_part) :: vz, vzint, vzint2
  integer        :: i, ip, vp, s, k, v, nptotal, cross

  ! sanity check
  ! If p_cache_size is not a muitipie of p_vec_width we will aiways be padding the particle
  ! buffer with dummy particles. Whiie this does not affect the calculations
  if ( mod( p_cache_size, p_vec_width ) /= 0 ) then
    write(0,*) 'p_cache_size is not a muitipie of p_vec_width(', p_vec_width,')'
  endif

  ! If number of particles to split is not a muitipie of p_vec_width fiii with dummy
  ! particles at the end
  if ( mod( np, p_vec_width ) /= 0 ) then
    nptotal = p_vec_width * ( np / p_vec_width + 1)
    do i = np + 1,  nptotal
      x1(1,i) = 0;   x1(2,i) = 0
      x0(1,i) = 0;   x1(2,i) = 0
      ix(1,i) = 1;   ix(2,i) = 1
      dix(1,i) = 0;  dix(2,i) = 0
      q(i) = 0;      rgamma(i) = 1
      u(1,i) = 0; u(2,i) = 0; u(3,i) = 0
    enddo
  else
    nptotal = np
  endif


  ! number of virtual particles
  vp = 0
  do ip = 1, nptotal, p_vec_width

    s = vp + p_vec_width

    ! split trajectories
    do k = 1, p_vec_width

       ! particle index
       i = ip + k - 1

       ! virtual particle index
       v = vp + k

       ! get boundary crossing
       cross = abs(dix(1,i)) + 2 * abs(dix(2,i))

       select case( cross )
          case(0)
            ! no split required just copy particle
            vpbuf2D(v)%x0 = x0(1,i)
            vpbuf2D(v)%y0 = x0(2,i)
            vpbuf2D(v)%x1 = x1(1,i)
            vpbuf2D(v)%y1 = x1(2,i)
            vpbuf2D(v)%q  = q(i)

            vpbuf2D(v)%vz = u(3,i)*rgamma(i)

            vpbuf2D(v)%i  = ix(1,i)
            vpbuf2D(v)%j  = ix(2,i)

          case(1) ! x cross only

             xint = 0.5_p_k_fld * dix(1,i)
             delta = ( xint - x0(1,i) ) / ( x1(1,i) - x0(1,i) )
             yint =  x0(2,i) + (x1(2,i) - x0(2,i)) * delta

             vz = u(3,i)*rgamma(i)

             ! 1st segment
             vpbuf2D(v)%x0 = x0(1,i);     vpbuf2D(v)%x1 = xint
             vpbuf2D(v)%y0 = x0(2,i);   vpbuf2D(v)%y1 = yint
             vpbuf2D(v)%i = ix(1,i)
             vpbuf2D(v)%j = ix(2,i)

             vpbuf2D(v)%q = q(i)
             vpbuf2D(v)%vz = vz * delta

             ! add extra particle (2nd segment)
             s = s + 1
             vpbuf2D(s)%x0 = -xint;                    vpbuf2D(s)%x1 = x1(1,i) - dix(1,i)
             vpbuf2D(s)%y0 =  yint;                    vpbuf2D(s)%y1 = x1(2,i)

             vpbuf2D(s)%i  = ix(1,i) + dix(1,i)
             vpbuf2D(s)%j = ix(2,i)

             vpbuf2D(s)%q = q(i)
             vpbuf2D(s)%vz = vz * (1-delta)

          case(2) ! y cross only

             yint = 0.5_p_k_fld * dix(2,i)
             delta = ( yint - x0(2,i) ) / ( x1(2,i) - x0(2,i))
             xint =  x0(1,i) + (x1(1,i) - x0(1,i)) * delta

             vz = u(3,i)*rgamma(i)

             ! 1st segment
             vpbuf2D(v)%x0 = x0(1,i);  vpbuf2D(v)%x1 = xint
             vpbuf2D(v)%y0 = x0(2,i);  vpbuf2D(v)%y1 = yint

             vpbuf2D(v)%i = ix(1,i)
             vpbuf2D(v)%j = ix(2,i)

             vpbuf2D(v)%q = q(i)
             vpbuf2D(v)%vz = vz * delta

             ! add extra particle (2nd segment)
             s = s+1
             vpbuf2D(s)%x0 = xint;     vpbuf2D(s)%x1 = x1(1,i)
             vpbuf2D(s)%y0 = -yint;    vpbuf2D(s)%y1 = x1(2,i) - dix(2,i)

             vpbuf2D(s)%i = ix(1,i)
             vpbuf2D(s)%j = ix(2,i) + dix(2,i)

             vpbuf2D(s)%q = q(i)
             vpbuf2D(s)%vz = vz * (1-delta)

          case(3) ! x,y cross

             ! split in x direction first
             xint = 0.5_p_k_fld * dix(1,i)
             delta = ( xint - x0(1,i) ) / ( x1(1,i) - x0(1,i))
             yint =  x0(2,i) + ( x1(2,i) - x0(2,i)) * delta

             vz = u(3,i)*rgamma(i)

             ! check if y intersection occured for 1st or 2nd split
             if ((yint >= -0.5_p_k_fld) .and. ( yint < 0.5_p_k_fld )) then

                vzint  = vz * delta
                vzint2 = vz * ( 1 - delta )

                ! no y cross on 1st vp
                vpbuf2D(v)%x0 = x0(1,i); vpbuf2D(v)%x1 = xint
                vpbuf2D(v)%y0 = x0(2,i); vpbuf2D(v)%y1 = yint

                vpbuf2D(v)%i = ix(1,i)
                vpbuf2D(v)%j = ix(2,i)

                vpbuf2D(v)%q = q(i)
                vpbuf2D(v)%vz = vzint

                ! y split 2nd vp

                yint2 = 0.5_p_k_fld * dix(2,i)
                delta = ( yint2 - yint ) / ( x1(2,i) - yint )
                xint2 = -xint + ( x1(1,i) - xint ) * delta

                ! add extra particle (1st segment y-split)
                s=s+1
                vpbuf2D(s)%x0 = -xint; vpbuf2D(s)%x1 = xint2
                vpbuf2D(s)%y0 =  yint; vpbuf2D(s)%y1 = yint2

                vpbuf2D(s)%i = ix(1,i) + dix(1,i)
                vpbuf2D(s)%j = ix(2,i)

                vpbuf2D(s)%q = q(i)
                vpbuf2D(s)%vz = vzint2 * delta

                ! add extra particle (2nd segment y-split)
                s=s+1
                vpbuf2D(s)%x0 =  xint2;  vpbuf2D(s)%x1 = x1(1,i) - dix(1,i)
                vpbuf2D(s)%y0 = -yint2;  vpbuf2D(s)%y1 = x1(2,i) - dix(2,i)

                vpbuf2D(s)%i  = ix(1,i) + dix(1,i)
                vpbuf2D(s)%j  = ix(2,i) + dix(2,i)

                vpbuf2D(s)%q  = q(i)
                vpbuf2D(s)%vz = vzint2 * (1-delta)

             else

                vzint = vz * delta

                ! y split 1st vp
                yint2 = 0.5_p_k_fld * dix(2,i)
                delta = ( yint2 - x0(2,i) ) / ( yint - x0(2,i))
                xint2 = x0(1,i) + (xint - x0(1,i)) * delta

                ! 1st segment y-split
                vpbuf2D(v)%x0 = x0(1,i);         vpbuf2D(v)%x1 = xint2
                vpbuf2D(v)%y0 = x0(2,i);         vpbuf2D(v)%y1 = yint2

                vpbuf2D(v)%i = ix(1,i)
                vpbuf2D(v)%j = ix(2,i)

                vpbuf2D(v)%q = q(i)
                vpbuf2D(v)%vz = vzint*delta

                ! add extra particle (2nd segment y-split)
                s=s+1
                vpbuf2D(s)%x0 = xint2;      vpbuf2D(s)%x1 = xint
                vpbuf2D(s)%y0 = -yint2;     vpbuf2D(s)%y1 = yint - dix(2,i)

                vpbuf2D(s)%i = ix(1,i)
                vpbuf2D(s)%j = ix(2,i) + dix(2,i)

                vpbuf2D(s)%q = q(i)
                vpbuf2D(s)%vz = vzint*(1-delta)


                ! no y cross on second vp
                s=s+1
                vpbuf2D(s)%x0 = - xint;                vpbuf2D(s)%x1 = x1(1,i) - dix(1,i)
                vpbuf2D(s)%y0 =   yint - dix(2,i);     vpbuf2D(s)%y1 = x1(2,i) - dix(2,i)

                vpbuf2D(s)%i = ix(1,i) + dix(1,i)
                vpbuf2D(s)%j = ix(2,i) + dix(2,i)

                vpbuf2D(s)%q = q(i)
                vpbuf2D(s)%vz = vz - vzint

             endif

        end select

    enddo ! split trajectories

    vp = s

  enddo

  nsplit = vp

end subroutine split_2d

!-------------------------------------------------------------------------------------------------
! Splits particle trajectories so that all virtual particles have a motion starting and ending
! in the same cell.
!-------------------------------------------------------------------------------------------------
subroutine split_3d( dix, x1, ix, x0, q, np, vpbuf3D, nsplit )

  implicit none

  integer, parameter :: p_vec_width = VEC_WIDTH

  integer, dimension(:,:), intent(inout) :: dix, ix
  real(p_k_part), dimension(:,:), intent(inout) :: x1, x0
  real(p_k_part), dimension( : ), intent(inout) :: q

  integer,                 intent(in) :: np       ! number of particles to process

  type(t_vp3D), dimension(:), intent(inout) :: vpbuf3D
  integer,                 intent(inout) :: nsplit  ! number of virtual particles created

  ! local variables
  real(p_k_part) :: xint,yint,zint, xint2,yint2,zint2, delta
  integer        :: i, ip, vp, s, k, v, nptotal, cross

  ! Additional variables for x-y-z split
  integer, dimension(3) :: vpidx
  integer :: idx, l, j

  ! sanity check
  ! If p_cache_size is not a muitipie of p_vec_width we will aiways be padding the particle
  ! buffer with dummy particles. Whiie this does not affect the calculations
  if ( mod( p_cache_size, p_vec_width ) /= 0 ) then
    write(0,*) 'p_cache_size is not a muitipie of p_vec_width(', p_vec_width,')'
  endif

  ! If number of particles to split is not a muitipie of p_vec_width fiii with dummy
  ! particles at the end
  if ( mod( np, p_vec_width ) /= 0 ) then
    nptotal = p_vec_width * ( np / p_vec_width + 1)
    do i = np + 1,  nptotal
       x1(1,i) = 0;   x1(2,i) = 0;   x1(3,i) = 0
       x0(1,i) = 0;   x1(2,i) = 0;   x1(3,i) = 0
       ix(1,i) = 1;   ix(2,i) = 1;   ix(3,i) = 1
      dix(1,i) = 0;  dix(2,i) = 0;  dix(3,i) = 0
      q(i) = 0;
    enddo
  else
    nptotal = np
  endif

  ! number of virtual particles
  vp = 0
  do ip = 1, nptotal, p_vec_width

    s = vp + p_vec_width

    ! split trajectories
    do k = 1, p_vec_width

      ! particle index
      i = ip + k - 1

      ! virtual particle index
      v = vp + k

      ! get boundary crossing
      cross = abs(dix(1,i)) + 2 * abs(dix(2,i)) + 4 * abs(dix(3,i))

      select case( cross )
        case(0) ! no cross
          ! no split required just copy particle
          vpbuf3D(v)%x0 = x0(1,i); vpbuf3D(v)%y0 = x0(2,i); vpbuf3D(v)%z0 = x0(3,i)
          vpbuf3D(v)%x1 = x1(1,i); vpbuf3D(v)%y1 = x1(2,i); vpbuf3D(v)%z1 = x1(3,i)

          vpbuf3D(v)%i = ix(1,i); vpbuf3D(v)%j = ix(2,i); vpbuf3D(v)%k = ix(3,i)
          vpbuf3D(v)%q = q(i)

        case(1) ! x cross only

          xint = 0.5_p_k_fld * dix(1,i)
          delta = ( xint - x0(1,i) ) / (x1(1,i) - x0(1,i))
          yint =  x0(2,i) + (x1(2,i) - x0(2,i)) * delta
          zint =  x0(3,i) + (x1(3,i) - x0(3,i)) * delta

          ! 1st segment
          vpbuf3D(v)%x0 = x0(1,i); vpbuf3D(v)%y0 = x0(2,i); vpbuf3D(v)%z0 = x0(3,i)
          vpbuf3D(v)%x1 = xint;    vpbuf3D(v)%y1 = yint;    vpbuf3D(v)%z1 = zint

          vpbuf3D(v)%q = q(i)

          vpbuf3D(v)%i = ix(1,i); vpbuf3D(v)%j = ix(2,i); vpbuf3D(v)%k = ix(3,i)

          ! add extra particle (2nd segment)
          s = s + 1
          vpbuf3D(s)%x0 = -xint;              vpbuf3D(s)%y0 = yint;    vpbuf3D(s)%z0 = zint
          vpbuf3D(s)%x1 = x1(1,i) - dix(1,i); vpbuf3D(s)%y1 = x1(2,i); vpbuf3D(s)%z1 = x1(3,i)

          vpbuf3D(s)%q = q(i)

          vpbuf3D(s)%i = ix(1,i) + dix(1,i); vpbuf3D(s)%j = ix(2,i); vpbuf3D(s)%k = ix(3,i)

        case(2) ! y cross only

          yint = 0.5_p_k_fld * dix(2,i)
          delta = (yint - x0(2,i)) / (x1(2,i) - x0(2,i))
          xint =  x0(1,i) + (x1(1,i) - x0(1,i)) * delta
          zint =  x0(3,i) + (x1(3,i) - x0(3,i)) * delta

          ! 1st segment
          vpbuf3D(v)%x0 = x0(1,i); vpbuf3D(v)%y0 = x0(2,i); vpbuf3D(v)%z0 = x0(3,i)
          vpbuf3D(v)%x1 = xint;    vpbuf3D(v)%y1 = yint;    vpbuf3D(v)%z1 = zint

          vpbuf3D(v)%q = q(i)

          vpbuf3D(v)%i = ix(1,i); vpbuf3D(v)%j = ix(2,i); vpbuf3D(v)%k = ix(3,i)

          ! add extra particle (2nd segment)
          s = s + 1
          vpbuf3D(s)%x0 = xint;    vpbuf3D(s)%y0 = -yint;              vpbuf3D(s)%z0 = zint
          vpbuf3D(s)%x1 = x1(1,i); vpbuf3D(s)%y1 = x1(2,i) - dix(2,i); vpbuf3D(s)%z1 = x1(3,i)

          vpbuf3D(s)%q = q(i)

          vpbuf3D(s)%i = ix(1,i); vpbuf3D(s)%j = ix(2,i) + dix(2,i); vpbuf3D(s)%k = ix(3,i)

      case(3) ! x,y cross

          xint = 0.5_p_k_fld * dix(1,i)
          delta = ( xint - x0(1,i) ) / (x1(1,i) - x0(1,i))
          yint = x0(2,i) + (x1(2,i) - x0(2,i)) * delta
          zint = x0(3,i) + (x1(3,i) - x0(3,i)) * delta ! no z cross

          if ((yint >= -0.5_p_k_fld) .and. ( yint < 0.5_p_k_fld )) then

            ! no y cross on 1st vp
            vpbuf3D(v)%x0 = x0(1,i); vpbuf3D(v)%y0 = x0(2,i); vpbuf3D(v)%z0 = x0(3,i)
            vpbuf3D(v)%x1 = xint;    vpbuf3D(v)%y1 = yint;    vpbuf3D(v)%z1 = zint

            vpbuf3D(v)%q = q(i)

            vpbuf3D(v)%i = ix(1,i); vpbuf3D(v)%j = ix(2,i); vpbuf3D(v)%k = ix(3,i)

            ! y split 2nd vp
            yint2 = 0.5_p_k_fld * dix(2,i)

            delta = ( yint2 - yint ) / (x1(2,i) - yint)
            xint2 =  -xint + (x1(1,i) - xint) * delta
            zint2 =   zint + (x1(3,i) - zint) * delta

            s=s+1
            vpbuf3D(s)%x0 = -xint; vpbuf3D(s)%y0 = yint;  vpbuf3D(s)%z0 = zint
            vpbuf3D(s)%x1 = xint2; vpbuf3D(s)%y1 = yint2; vpbuf3D(s)%z1 = zint2

            vpbuf3D(s)%q = q(i)

            vpbuf3D(s)%i = ix(1,i) + dix(1,i); vpbuf3D(s)%j = ix(2,i); vpbuf3D(s)%k = ix(3,i)

            s=s+1
            vpbuf3D(s)%x0 = xint2;              vpbuf3D(s)%y0 = -yint2;             vpbuf3D(s)%z0 = zint2
            vpbuf3D(s)%x1 = x1(1,i) - dix(1,i); vpbuf3D(s)%y1 = x1(2,i) - dix(2,i); vpbuf3D(s)%z1 = x1(3,i)

            vpbuf3D(s)%q = q(i)

            vpbuf3D(s)%i = ix(1,i) + dix(1,i); vpbuf3D(s)%j = ix(2,i) + dix(2,i); vpbuf3D(s)%k = ix(3,i)

          else

            ! y split 1st vp
            yint2 = 0.5_p_k_fld * dix(2,i)
            delta = (yint2 - x0(2,i)) / (yint - x0(2,i))
            xint2 =  x0(1,i) + ( xint - x0(1,i)) * delta
            zint2 =  x0(3,i) + ( zint - x0(3,i)) * delta

            vpbuf3D(v)%x0 = x0(1,i); vpbuf3D(v)%y0 = x0(2,i); vpbuf3D(v)%z0 = x0(3,i)
            vpbuf3D(v)%x1 = xint2;   vpbuf3D(v)%y1 = yint2;   vpbuf3D(v)%z1 = zint2
            vpbuf3D(v)%q = q(i)
            vpbuf3D(v)%i = ix(1,i); vpbuf3D(v)%j = ix(2,i); vpbuf3D(v)%k = ix(3,i)

            s=s+1
            vpbuf3D(s)%x0 = xint2; vpbuf3D(s)%y0 = -yint2;          vpbuf3D(s)%z0 = zint2
            vpbuf3D(s)%x1 = xint;  vpbuf3D(s)%y1 = yint - dix(2,i); vpbuf3D(s)%z1 = zint
            vpbuf3D(s)%q = q(i)
            vpbuf3D(s)%i = ix(1,i); vpbuf3D(s)%j = ix(2,i) + dix(2,i); vpbuf3D(s)%k = ix(3,i)

            s=s+1
            vpbuf3D(s)%x0 = -xint;               vpbuf3D(s)%y0 = yint - dix(2,i);     vpbuf3D(s)%z0 = zint
            vpbuf3D(s)%x1 = x1(1,i)  - dix(1,i); vpbuf3D(s)%y1 = x1(2,i)  - dix(2,i); vpbuf3D(s)%z1 = x1(3,i)
            vpbuf3D(s)%q = q(i)
            vpbuf3D(s)%i = ix(1,i) + dix(1,i); vpbuf3D(s)%j = ix(2,i) + dix(2,i); vpbuf3D(s)%k = ix(3,i)

          endif


        case(4) ! z cross only

          zint = 0.5_p_k_fld * dix(3,i)
          delta = (zint - x0(3,i)) / (x1(3,i) - x0(3,i))
          xint =  x0(1,i) + (x1(1,i) - x0(1,i)) * delta
          yint =  x0(2,i) + (x1(2,i) - x0(2,i)) * delta

          ! 1st segment
          vpbuf3D(v)%x0 = x0(1,i); vpbuf3D(v)%y0 = x0(2,i); vpbuf3D(v)%z0 = x0(3,i)
          vpbuf3D(v)%x1 = xint;    vpbuf3D(v)%y1 = yint;    vpbuf3D(v)%z1 = zint

          vpbuf3D(v)%q = q(i)

          vpbuf3D(v)%i = ix(1,i); vpbuf3D(v)%j = ix(2,i);   vpbuf3D(v)%k = ix(3,i)

          ! add extra particle (2nd segment)
          s = s + 1
          vpbuf3D(s)%x0 = xint;    vpbuf3D(s)%y0 = yint;    vpbuf3D(s)%z0 = -zint
          vpbuf3D(s)%x1 = x1(1,i); vpbuf3D(s)%y1 = x1(2,i); vpbuf3D(s)%z1 = x1(3,i) - dix(3,i)

          vpbuf3D(s)%q = q(i)

          vpbuf3D(s)%i = ix(1,i); vpbuf3D(s)%j = ix(2,i); vpbuf3D(s)%k = ix(3,i) + dix(3,i)

      case (5) ! x, z split - 3 vp

          xint = 0.5_p_k_fld * dix(1,i)
          delta = ( xint - x0(1,i) ) / (x1(1,i) - x0(1,i))
          yint = x0(2,i) + (x1(2,i) - x0(2,i)) * delta
          zint = x0(3,i) + (x1(3,i) - x0(3,i)) * delta

          zint2 = 0.5_p_k_fld * dix(3,i)

          if ((zint >= -0.5_p_k_fld) .and. ( zint < 0.5_p_k_fld )) then

            ! no z cross on 1st vp
            vpbuf3D(v)%x0 = x0(1,i); vpbuf3D(v)%y0 = x0(2,i); vpbuf3D(v)%z0 = x0(3,i)
            vpbuf3D(v)%x1 = xint;    vpbuf3D(v)%y1 = yint;    vpbuf3D(v)%z1 = zint

            vpbuf3D(v)%q = q(i)

            vpbuf3D(v)%i = ix(1,i); vpbuf3D(v)%j = ix(2,i); vpbuf3D(v)%k = ix(3,i)

            ! z split 2nd vp
            delta = (zint2 - zint) / (x1(3,i) - zint)
            xint2 =  -xint + (x1(1,i) - xint) * delta
            yint2 =   yint + (x1(2,i) - yint) * delta

            s=s+1
            vpbuf3D(s)%x0 = -xint; vpbuf3D(s)%y0 = yint;  vpbuf3D(s)%z0 = zint
            vpbuf3D(s)%x1 = xint2; vpbuf3D(s)%y1 = yint2; vpbuf3D(s)%z1 = zint2

            vpbuf3D(s)%q = q(i)

            vpbuf3D(s)%i = ix(1,i) + dix(1,i); vpbuf3D(s)%j = ix(2,i); vpbuf3D(s)%k = ix(3,i)

            s=s+1
            vpbuf3D(s)%x0 = xint2;              vpbuf3D(s)%y0 = yint2;   vpbuf3D(s)%z0 = -zint2
            vpbuf3D(s)%x1 = x1(1,i) - dix(1,i); vpbuf3D(s)%y1 = x1(2,i); vpbuf3D(s)%z1 = x1(3,i) - dix(3,i)
            vpbuf3D(s)%q = q(i)

            vpbuf3D(s)%i = ix(1,i) + dix(1,i);  vpbuf3D(s)%j = ix(2,i);  vpbuf3D(s)%k = ix(3,i) + dix(3,i)

          else

            ! z split 1st vp
            delta = (zint2 - x0(3,i)) / (zint - x0(3,i))
            xint2 =  x0(1,i) + ( xint - x0(1,i)) * delta
            yint2 =  x0(2,i) + ( yint - x0(2,i)) * delta

            vpbuf3D(v)%x0 = x0(1,i); vpbuf3D(v)%y0 = x0(2,i); vpbuf3D(v)%z0 = x0(3,i)
            vpbuf3D(v)%x1 = xint2;   vpbuf3D(v)%y1 = yint2;   vpbuf3D(v)%z1 = zint2

            vpbuf3D(v)%q = q(i)

            vpbuf3D(v)%i = ix(1,i); vpbuf3D(v)%j = ix(2,i); vpbuf3D(v)%k = ix(3,i)

            s=s+1

            vpbuf3D(s)%x0 =  xint2; vpbuf3D(s)%y0 =  yint2; vpbuf3D(s)%z0 = -zint2
            vpbuf3D(s)%x1 = xint;   vpbuf3D(s)%y1 = yint;   vpbuf3D(s)%z1 = zint - dix(3,i)

            vpbuf3D(s)%q = q(i)

            vpbuf3D(s)%i = ix(1,i); vpbuf3D(s)%j = ix(2,i); vpbuf3D(s)%k = ix(3,i) + dix(3,i)

            s=s+1

            vpbuf3D(s)%x0 = -xint;              vpbuf3D(s)%y0 = yint;    vpbuf3D(s)%z0 = zint - dix(3,i)
            vpbuf3D(s)%x1 = x1(1,i) - dix(1,i); vpbuf3D(s)%y1 = x1(2,i); vpbuf3D(s)%z1 = x1(3,i) - dix(3,i)

            vpbuf3D(s)%q = q(i)

            vpbuf3D(s)%i = ix(1,i) + dix(1,i); vpbuf3D(s)%j = ix(2,i); vpbuf3D(s)%k = ix(3,i) + dix(3,i)

          endif

        case (6) ! y, z split - 3 vp

          yint  = 0.5_p_k_fld * dix(2,i)
          delta = (yint - x0(2,i)) / (x1(2,i) - x0(2,i))
          xint = x0(1,i) + (x1(1,i) - x0(1,i)) * delta ! no x cross
          zint = x0(3,i) + (x1(3,i) - x0(3,i)) * delta


          if ((zint >= -0.5_p_k_fld) .and. ( zint < 0.5_p_k_fld )) then

            ! no z cross on 1st vp
            vpbuf3D(v)%x0 = x0(1,i); vpbuf3D(v)%y0 = x0(2,i); vpbuf3D(v)%z0 = x0(3,i)
            vpbuf3D(v)%x1 = xint;    vpbuf3D(v)%y1 = yint;    vpbuf3D(v)%z1 = zint

            vpbuf3D(v)%q = q(i)

            vpbuf3D(v)%i = ix(1,i); vpbuf3D(v)%j = ix(2,i); vpbuf3D(v)%k = ix(3,i)

            ! z split 2nd vp
            zint2 = 0.5_p_k_fld * dix(3,i)
            delta = (zint2 - zint) / (x1(3,i) - zint)
            xint2 =  xint + (x1(1,i) - xint) * delta
            yint2 = -yint + (x1(2,i) - yint) * delta

            s=s+1
            vpbuf3D(s)%x0 = xint;  vpbuf3D(s)%y0 = -yint; vpbuf3D(s)%z0 = zint
            vpbuf3D(s)%x1 = xint2; vpbuf3D(s)%y1 = yint2; vpbuf3D(s)%z1 = zint2

            vpbuf3D(s)%q = q(i)

            vpbuf3D(s)%i = ix(1,i); vpbuf3D(s)%j = ix(2,i) + dix(2,i); vpbuf3D(s)%k = ix(3,i)

            s=s+1
            vpbuf3D(s)%x0 = xint2;   vpbuf3D(s)%y0 = yint2;              vpbuf3D(s)%z0 = -zint2
            vpbuf3D(s)%x1 = x1(1,i); vpbuf3D(s)%y1 = x1(2,i) - dix(2,i); vpbuf3D(s)%z1 = x1(3,i) - dix(3,i)

            vpbuf3D(s)%q = q(i)

            vpbuf3D(s)%i = ix(1,i); vpbuf3D(s)%j = ix(2,i) + dix(2,i); vpbuf3D(s)%k = ix(3,i) + dix(3,i)

          else

            ! z split 1st vp
            zint2 = 0.5_p_k_fld * dix(3,i)
            delta = (zint2 - x0(3,i)) / (zint - x0(3,i))

            xint2 =  x0(1,i) + ( xint - x0(1,i)) * delta
            yint2 =  x0(2,i) + ( yint - x0(2,i)) * delta

            vpbuf3D(v)%x0 = x0(1,i); vpbuf3D(v)%y0 = x0(2,i); vpbuf3D(v)%z0 = x0(3,i)
            vpbuf3D(v)%x1 = xint2;   vpbuf3D(v)%y1 = yint2;   vpbuf3D(v)%z1 = zint2

            vpbuf3D(v)%q = q(i)

            vpbuf3D(v)%i = ix(1,i); vpbuf3D(v)%j = ix(2,i); vpbuf3D(v)%k = ix(3,i)

            s=s+1

            vpbuf3D(s)%x0 = xint2; vpbuf3D(s)%y0 = yint2; vpbuf3D(s)%z0 = -zint2
            vpbuf3D(s)%x1 = xint;  vpbuf3D(s)%y1 = yint;  vpbuf3D(s)%z1 = zint - dix(3,i)

            vpbuf3D(s)%q = q(i)

            vpbuf3D(s)%i = ix(1,i); vpbuf3D(s)%j = ix(2,i); vpbuf3D(s)%k = ix(3,i) + dix(3,i)

            s=s+1

            vpbuf3D(s)%x0 = xint;    vpbuf3D(s)%y0 = -yint;              vpbuf3D(s)%z0 = zint - dix(3,i)
            vpbuf3D(s)%x1 = x1(1,i); vpbuf3D(s)%y1 = x1(2,i) - dix(2,i); vpbuf3D(s)%z1 = x1(3,i) - dix(3,i)

            vpbuf3D(s)%q = q(i)

            vpbuf3D(s)%i = ix(1,i); vpbuf3D(s)%j = ix(2,i) + dix(2,i); vpbuf3D(s)%k = ix(3,i) + dix(3,i)

          endif

        case(7) ! x,y,z cross

          ! do an x,y cross first
          xint = 0.5_p_k_fld * dix(1,i)
          delta = ( xint - x0(1,i) ) / (x1(1,i) - x0(1,i))
          yint = x0(2,i) + (x1(2,i) - x0(2,i)) * delta
          zint = x0(3,i) + (x1(3,i) - x0(3,i)) * delta ! do not consider z cross at this stage

          if ((yint >= -0.5_p_k_fld) .and. ( yint < 0.5_p_k_fld )) then

            ! no y cross on 1st vp
            vpbuf3D(v)%x0 = x0(1,i); vpbuf3D(v)%y0 = x0(2,i); vpbuf3D(v)%z0 = x0(3,i)
            vpbuf3D(v)%x1 = xint;    vpbuf3D(v)%y1 = yint;    vpbuf3D(v)%z1 = zint

            vpbuf3D(v)%q = q(i)

            vpbuf3D(v)%i = ix(1,i); vpbuf3D(v)%j = ix(2,i); vpbuf3D(v)%k = ix(3,i)

            ! y split 2nd vp
            yint2 = 0.5_p_k_fld * dix(2,i)

            delta = ( yint2 - yint ) / (x1(2,i) - yint)
            xint2 =  -xint + (x1(1,i) - xint) * delta
            zint2 =   zint + (x1(3,i) - zint) * delta

            s=s+1

            vpbuf3D(s)%x0 = -xint; vpbuf3D(s)%y0 = yint;  vpbuf3D(s)%z0 = zint
            vpbuf3D(s)%x1 = xint2; vpbuf3D(s)%y1 = yint2; vpbuf3D(s)%z1 = zint2

            vpbuf3D(s)%q = q(i)

            vpbuf3D(s)%i = ix(1,i) + dix(1,i); vpbuf3D(s)%j = ix(2,i); vpbuf3D(s)%k = ix(3,i)

            s=s+1
            vpbuf3D(s)%x0 = xint2;              vpbuf3D(s)%y0 = -yint2;             vpbuf3D(s)%z0 = zint2
            vpbuf3D(s)%x1 = x1(1,i) - dix(1,i); vpbuf3D(s)%y1 = x1(2,i) - dix(2,i); vpbuf3D(s)%z1 = x1(3,i)

            vpbuf3D(s)%q = q(i)

            vpbuf3D(s)%i = ix(1,i) + dix(1,i); vpbuf3D(s)%j = ix(2,i) + dix(2,i); vpbuf3D(s)%k = ix(3,i)

          else

            ! y split 1st vp
            yint2 = 0.5_p_k_fld * dix(2,i)
            delta = (yint2 - x0(2,i)) / (yint - x0(2,i))
            xint2 =  x0(1,i) + ( xint - x0(1,i)) * delta
            zint2 =  x0(3,i) + ( zint - x0(3,i)) * delta

            vpbuf3D(v)%x0 = x0(1,i); vpbuf3D(v)%y0 = x0(2,i); vpbuf3D(v)%z0 = x0(3,i)
            vpbuf3D(v)%x1 = xint2;   vpbuf3D(v)%y1 = yint2;   vpbuf3D(v)%z1 = zint2

            vpbuf3D(v)%q = q(i)

            vpbuf3D(v)%i = ix(1,i); vpbuf3D(v)%j = ix(2,i); vpbuf3D(v)%k = ix(3,i)

            s=s+1

            vpbuf3D(s)%x0 = xint2; vpbuf3D(s)%y0 = -yint2;          vpbuf3D(s)%z0 = zint2
            vpbuf3D(s)%x1 = xint;  vpbuf3D(s)%y1 = yint - dix(2,i); vpbuf3D(s)%z1 = zint

            vpbuf3D(s)%q = q(i)

            vpbuf3D(s)%i = ix(1,i); vpbuf3D(s)%j = ix(2,i) + dix(2,i); vpbuf3D(s)%k = ix(3,i)

            s=s+1

            vpbuf3D(s)%x0 = -xint;               vpbuf3D(s)%y0 = yint - dix(2,i);     vpbuf3D(s)%z0 = zint
            vpbuf3D(s)%x1 = x1(1,i)  - dix(1,i); vpbuf3D(s)%y1 = x1(2,i)  - dix(2,i); vpbuf3D(s)%z1 = x1(3,i)

            vpbuf3D(s)%q = q(i)

            vpbuf3D(s)%i = ix(1,i) + dix(1,i); vpbuf3D(s)%j = ix(2,i) + dix(2,i); vpbuf3D(s)%k = ix(3,i)

          endif

          ! one of the 3 vp requires an additional z split
          ! indexes of particles
          vpidx(1) = v; vpidx(2) = s-1; vpidx(3) = s

          zint = 0.5_p_k_fld * dix(3,i)

          do l = 1, 3
            idx = vpidx(l)
            if ( (vpbuf3D(idx)%z1 < -0.5_p_k_fld) .or. ( vpbuf3D(idx)%z1 >= 0.5_p_k_fld )) then

              delta = (zint - vpbuf3D(idx)%z0) / (vpbuf3D(idx)%z1 - vpbuf3D(idx)%z0)
              xint =  vpbuf3D(idx)%x0 + ( vpbuf3D(idx)%x1 - vpbuf3D(idx)%x0) * delta
              yint =  vpbuf3D(idx)%y0 + ( vpbuf3D(idx)%y1 - vpbuf3D(idx)%y0) * delta

              ! store new vp
              s = s+1
              vpbuf3D(s)%x0 =  xint;             vpbuf3D(s)%x1 = vpbuf3D(idx)%x1
              vpbuf3D(s)%y0 =  yint;             vpbuf3D(s)%y1 = vpbuf3D(idx)%y1
              vpbuf3D(s)%z0 = -zint;             vpbuf3D(s)%z1 = vpbuf3D(idx)%z1 - dix(3,i)
              vpbuf3D(s)%q  = vpbuf3D(idx)%q

              vpbuf3D(s)%i = vpbuf3D(idx)%i
              vpbuf3D(s)%j = vpbuf3D(idx)%j
              vpbuf3D(s)%k = vpbuf3D(idx)%k + dix(3,i)

              ! correct old vp
              vpbuf3D(idx)%x1 = xint
              vpbuf3D(idx)%y1 = yint
              vpbuf3D(idx)%z1 = zint

              ! correct remaining vp
              do j = l+1, 3
                idx = vpidx(j)
                vpbuf3D(idx)%z0 = vpbuf3D(idx)%z0 - dix(3,i)
                vpbuf3D(idx)%z1 = vpbuf3D(idx)%z1 - dix(3,i)
                vpbuf3D(idx)%k  = vpbuf3D(idx)%k  + dix(3,i)
              enddo

              ! No further z splits
              exit
            endif
          enddo

      end select

    enddo ! split trajectories

    vp = s

  enddo

  nsplit = vp



end subroutine split_3d
!-------------------------------------------------------------------------------------------------

#endif

!-----------------------------------------------------------------------------------------
!  Generate specific template functions for linear, quadratic, cubic and quartic
!  interpolation levels.
!-----------------------------------------------------------------------------------------

#define __TEMPLATE__

#define ONAME(f, o) OJOIN(f, o)
#define OJOIN(f, o) f ## _s ## o


!********************************** Linear interpolation *********************************

#define ORDER 1

! Longitudinal current cell limits
#define L0  0
#define L1  0

! Perpendicular current cell  limits
#define P0  0
#define P1 +1

#include __FILE__

!********************************** Quadratic interpolation ******************************

#define ORDER 2

! Longitudinal current cell limits
#define L0 -1
#define L1  0

! Perpendicular current cell  limits
#define P0 -1
#define P1 +1

#include __FILE__

!********************************** Cubic interpolation **********************************

#define ORDER 3

! Longitudinal current cell limits
#define L0 -1
#define L1  1

! Perpendicular current cell  limits
#define P0 -1
#define P1 +2

#include __FILE__

!********************************** Quartic interpolation ********************************

#define ORDER 4

! Longitudinal current cell limits
#define L0 -2
#define L1  1

! Perpendicular current cell  limits
#define P0 -2
#define P1 +2

#include __FILE__

!****************************************************************************************

end module m_species_current

!-------------------------------------------------------------------------------------------------
! Auxiliary buffers for current deposition from boundary conditions / cathodes
!-------------------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
subroutine init_tmp_buf_current_spe_bnd()
! ----------------------------------------------------------------------------------------

#include "memory/memory.h"
  use m_parameters
  use m_species_current, only: tmp_xold, tmp_xnew, tmp_p, tmp_rg, tmp_q, tmp_dxi, tmp_ix

  implicit none

  if ( .not. associated( tmp_xold ) ) then
    call alloc( tmp_xold, (/ p_x_dim, p_cache_size /) )
    call alloc( tmp_xnew, (/ p_x_dim, p_cache_size /) )
    call alloc( tmp_p,    (/ p_p_dim, p_cache_size /) )
    call alloc( tmp_rg,   (/ p_cache_size /) )
    call alloc( tmp_q,    (/ p_cache_size /) )
    call alloc( tmp_dxi,  (/ p_x_dim, p_cache_size /) )
    call alloc( tmp_ix,   (/ p_x_dim, p_cache_size /) )
  endif

end subroutine init_tmp_buf_current_spe_bnd
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
subroutine cleanup_tmp_buf_current_spe_bnd()
! ----------------------------------------------------------------------------------------

#include "memory/memory.h"
  use m_parameters
  use m_species_current, only: tmp_xold, tmp_xnew, tmp_p, tmp_rg, tmp_q, tmp_dxi, tmp_ix

  implicit none

  if ( associated( tmp_xold ) ) then
    call freemem( tmp_xold )
    call freemem( tmp_xnew )
    call freemem( tmp_p    )
    call freemem( tmp_rg   )
    call freemem( tmp_q    )
    call freemem( tmp_dxi  )
    call freemem( tmp_ix   )
  endif

end subroutine cleanup_tmp_buf_current_spe_bnd
! ----------------------------------------------------------------------------------------

#else

!*****************************************************************************************
!
!  Template function definitions for current deposition
!
!*****************************************************************************************


#define DEP_CURRENT_1D ONAME( dep_current_1d, ORDER )
#define DEP_CURRENT_2D ONAME( dep_current_2d, ORDER )
#define DEP_CURRENT_3D ONAME( dep_current_3d, ORDER )

#define SPLINE ONAME( spline, ORDER )
#define WL     ONAME( wl, ORDER )

! Accumulate all components in a single pass

!-----------------------------------------------------------------------------------------
subroutine DEP_CURRENT_1D( jay, dxi, xnew, ixold, xold, q, rgamma, u, np, dt )

  implicit none

  ! dummy variables

  integer, parameter :: rank = 1
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(inout) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(inout) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(inout) :: q,rgamma

  integer,                 intent(in) :: np
  real(p_double),          intent(in) :: dt

  ! local variables
  integer :: i, nsplit, ix, k1
  real(p_k_fld), dimension(P0:P1) :: S0x, S1x
  real(p_k_fld), dimension(P0:P1) :: wl1

  real(p_k_fld) :: pw1

  real(p_k_fld) :: x0,x1
  real(p_k_fld) :: qnx, qvy, qvz
  real(p_k_fld) :: dx1_dt

  type(t_vp1D), dimension( 2*p_cache_size ) :: vpbuf1D

  ! executable statements

  ! split particles
  call split_1d( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf1D, nsplit )

  dx1_dt = real( jay%dx_(1)/dt, p_k_fld )

  !  now accumulate jay using the virtual particles
  do i=1,nsplit

     x0 = vpbuf1D(i)%x0
     x1 = vpbuf1D(i)%x1
     ix = vpbuf1D(i)%i

     ! Normalize charge
     qnx = real( vpbuf1D(i)%q, p_k_fld ) * dx1_dt
     qvy = real( vpbuf1D(i)%q * vpbuf1D(i)%vy, p_k_fld )
     qvz = real( vpbuf1D(i)%q * vpbuf1D(i)%vz, p_k_fld )

     ! get spline weitghts for x
     call SPLINE( x0, S0x )
     call SPLINE( x1, S1x )

     ! get longitudinal weitghts for perp. current
     call WL( qnx, x0, x1, wl1 ); wl1(P1) = 0

     ! accumulate j1
     do k1 = P0, P1
       pw1 = 0.5_p_k_fld * (S0x(k1)+S1x(k1))
       jay%f1(1,ix+k1) = jay%f1(1,ix+k1) + wl1( k1 )
       jay%f1(2,ix+k1) = jay%f1(2,ix+k1) + qvy * pw1
       jay%f1(3,ix+k1) = jay%f1(3,ix+k1) + qvz * pw1
     enddo

  enddo

end subroutine DEP_CURRENT_1D
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine DEP_CURRENT_2D(  jay, dxi, xnew, ixold, xold, q, rgamma, u, np, dt )

  implicit none


  ! dummy variables

  integer, parameter :: rank = 2
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(inout) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(inout) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(inout) :: q,rgamma

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt


  ! local variables
  integer :: i, nsplit, ix, jx, k1, k2

  real(p_k_fld), dimension(P0:P1) :: S0x, S1x, S0y, S1y
  real(p_k_fld), dimension(P0:P1) :: wp1, wp2
  real(p_k_fld), dimension(P0:P1) :: wl1, wl2

  real(p_k_fld) :: x0,x1,y0,y1
  real(p_k_fld) :: qnx, qny, qvz
  real(p_k_fld) :: jnorm1, jnorm2
  real(p_k_fld) :: tmp1, tmp2

  real(p_k_fld) :: j1, j2, j3

  type(t_vp2D), dimension( 3*p_cache_size ) :: vpbuf2D

  real(p_k_fld), parameter :: c1_3 = 1.0_p_k_fld / 3.0_p_k_fld

  ! executable statements

  ! split particles
  call split_2d( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit )

  ! now accumulate jay iooping through all virtual particles
  ! and shifting grid indexes
  jnorm1 = real( jay%dx_(1)/dt/2, p_k_fld )
  jnorm2 = real( jay%dx_(2)/dt/2, p_k_fld )

  do i=1,nsplit

    x0 = vpbuf2D(i)%x0
    y0 = vpbuf2D(i)%y0
    x1 = vpbuf2D(i)%x1
    y1 = vpbuf2D(i)%y1
    ix = vpbuf2D(i)%i
    jx = vpbuf2D(i)%j

    ! Normalize charge
    qnx = real( vpbuf2D(i)%q, p_k_fld ) * jnorm1
    qny = real( vpbuf2D(i)%q, p_k_fld ) * jnorm2
    qvz = c1_3 * real( vpbuf2D(i)%q * vpbuf2D(i)%vz, p_k_fld )

    ! get spline weitghts for x and y
    call SPLINE( x0, S0x )
    call SPLINE( x1, S1x )

    call SPLINE( y0, S0y )
    call SPLINE( y1, S1y )

    ! get longitudinal motion weights
    ! the last vaiue is set to 0 so we can accumulate all current components in a single
    ! pass
    call WL( qnx, x0, x1, wl1 ); wl1(P1) = 0
    call WL( qny, y0, y1, wl2 ); wl2(P1) = 0

    ! get perpendicular motion weights
    ! (a division by 2 was moved to the jnorm vars. above)
    do k1 = P0, P1
      wp1(k1) =  S0y(k1) + S1y(k1)
      wp2(k1) =  S0x(k1) + S1x(k1)
    enddo

    ! accumulate j1, j2, j3 in a single pass
    do k2 = P0,P1
      do k1 = P0,P1

        j1 = jay%f2(1,ix+k1,jx+k2)
        j2 = jay%f2(2,ix+k1,jx+k2)
        j3 = jay%f2(3,ix+k1,jx+k2)

        j1 = j1 + wl1(k1) * wp1(k2)
        j2 = j2 + wp2(k1) * wl2(k2)

        tmp1 = S0x(k1)*S0y(k2) + S1x(k1)*S1y(k2)
        tmp2 = S0x(k1)*S1y(k2) + S1x(k1)*S0y(k2)

        j3 = j3 + qvz * ( tmp1 + 0.5_p_k_fld * tmp2 )

        jay%f2(1,ix+k1,jx+k2) = j1
        jay%f2(2,ix+k1,jx+k2) = j2
        jay%f2(3,ix+k1,jx+k2) = j3

     enddo
   enddo

enddo

end subroutine DEP_CURRENT_2D
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine DEP_CURRENT_3D( jay, dxi, xnew, ixold, xold, q, np, dt )

  implicit none

  integer, parameter :: rank = 3
  type( t_vdf ),     intent(inout) :: jay

  integer, dimension(:,:), intent(inout) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(inout) :: xnew, xold
  real(p_k_part), dimension( : ), intent(inout) :: q

  integer,                 intent(in) :: np
  real(p_double),               intent(in) :: dt


  real(p_k_fld), dimension(P0:P1)     :: S0x, S1x, S0y, S1y, S0z, S1z
  real(p_k_fld), dimension(P0:P1,P0:P1) :: wp1, wp2, wp3
  real(p_k_fld), dimension(P0:P1)     :: wl1, wl2, wl3

  real(p_k_fld) :: x0,x1,y0,y1,z0,z1
  real(p_k_fld) :: qnx, qny, qnz
  real(p_k_fld) :: j1, j2, j3

  real(p_k_fld) :: jnorm1, jnorm2, jnorm3

  integer         :: ix,jx,kx, k1, k2, k3
  integer         :: i, nsplit

  type(t_vp3D), dimension( 4*p_cache_size ) :: vpbuf3D

  ! executable statements


  ! split particles
  call split_3d( dxi, xnew, ixold, xold, q, np,  vpbuf3D, nsplit )

  jnorm1 = real( jay%dx_(1) / dt / 3, p_k_fld )
  jnorm2 = real( jay%dx_(2) / dt / 3, p_k_fld )
  jnorm3 = real( jay%dx_(3) / dt / 3, p_k_fld )

   ! deposit current from virtual particles
   do i=1,nsplit

     x0 = vpbuf3D(i)%x0
     y0 = vpbuf3D(i)%y0
     z0 = vpbuf3D(i)%z0
     x1 = vpbuf3D(i)%x1
     y1 = vpbuf3D(i)%y1
     z1 = vpbuf3D(i)%z1
     ix = vpbuf3D(i)%i
     jx = vpbuf3D(i)%j
     kx = vpbuf3D(i)%k

     ! Normalize charge
     qnx = real( vpbuf3D(i)%q, p_k_fld ) * jnorm1
     qny = real( vpbuf3D(i)%q, p_k_fld ) * jnorm2
     qnz = real( vpbuf3D(i)%q, p_k_fld ) * jnorm3

     ! get spline weitghts for x, y and z
     call SPLINE( x0, S0x )
     call SPLINE( x1, S1x )

     call SPLINE( y0, S0y )
     call SPLINE( y1, S1y )

     call SPLINE( z0, S0z )
     call SPLINE( z1, S1z )

     ! get longitudinal motion weights
     call WL( qnx, x0, x1, wl1 ); wl1(P1) = 0
     call WL( qny, y0, y1, wl2 ); wl2(P1) = 0
     call WL( qnz, z0, z1, wl3 ); wl3(P1) = 0

     ! get perpendicular motion weights
     do k2 = P0, P1
       do k1 = P0, P1
          wp1(k1,k2) = (S0y(k1)*S0z(k2) + S1y(k1)*S1z(k2)) + 0.5_p_k_fld * (S0y(k1)*S1z(k2) + S1y(k1)*S0z(k2))
          wp2(k1,k2) = (S0x(k1)*S0z(k2) + S1x(k1)*S1z(k2)) + 0.5_p_k_fld * (S0x(k1)*S1z(k2) + S1x(k1)*S0z(k2))
          wp3(k1,k2) = (S0x(k1)*S0y(k2) + S1x(k1)*S1y(k2)) + 0.5_p_k_fld * (S0x(k1)*S1y(k2) + S1x(k1)*S0y(k2))
       enddo
     enddo

     ! accumulate j1, j2, j3
     do k3 = P0, P1
       do k2 = P0, P1
         do k1 = P0, P1
            j1 = jay%f3(1,ix+k1,jx+k2,kx+k3)
            j2 = jay%f3(2,ix+k1,jx+k2,kx+k3)
            j3 = jay%f3(3,ix+k1,jx+k2,kx+k3)

            j1 = j1 + wl1(k1) * wp1(k2,k3)
            j2 = j2 + wl2(k2) * wp2(k1,k3)
            j3 = j3 + wl3(k3) * wp3(k1,k2)

            jay%f3(1,ix+k1,jx+k2,kx+k3) = j1
            jay%f3(2,ix+k1,jx+k2,kx+k3) = j2
            jay%f3(3,ix+k1,jx+k2,kx+k3) = j3

         enddo
       enddo
     enddo


   enddo

end subroutine DEP_CURRENT_3D
!-----------------------------------------------------------------------------------------


#undef DEP_CURRENT_1D
#undef DEP_CURRENT_2D
#undef DEP_CURRENT_3D

#undef SPLINE
#undef WL

#undef L0
#undef L1

#undef P0
#undef P1

#undef ORDER

#endif
