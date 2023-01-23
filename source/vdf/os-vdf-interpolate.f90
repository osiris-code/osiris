!---------------------------------------------------------------------------------------------------
! Grid -> Particle interpolation
!
! The rountines in this module allow for the interpolation of a quantity stored in a vdf grid to be
! interpolated at given particle positions:
!  - Only one field component can be interpolated at a time.
!  - The routines also work for quantities not defined in the lower left corner of the cell (e.g.
!    EM field components)
!  - linear, quadratic, cubic and quartic interpolations supported
!  - The routines expect the positions to be defined as "cell" particle positions, namely cell
!    index and cell position, normalized to the cell size
!  - The routines expect the cell positions to be referred to the lower boundary of the cell for
!    odd valued interpolation levels, and nearest cell boundary for even valued interpolation
!    levels
!
! When interpolating multiple field components of the same vdf, if performance is critical then a
! specific routine should be written, similar to the get_emf* routines in os-emf-interpolate.f90
!---------------------------------------------------------------------------------------------------

#include "os-preprocess.fpp"
#include "os-config.h"


module m_vdf_interpolate

use m_parameters
use m_vdf_define

implicit none

private

integer, parameter :: p_buf_size = 1024

interface interpolate
  module procedure interpolate_vdf
end interface

public :: interpolate

contains

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine interpolate_vdf( vdf, fc, fpos, ix, x, np, int_type, fval )

  implicit none

  type( t_vdf ), intent(in) :: vdf
  integer, intent(in) :: fc
  integer, intent(in) :: fpos

  integer, dimension(:,:), intent(in) :: ix
  real( p_k_part ), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np
  integer, intent(in) :: int_type

  real( p_k_fld ), dimension(:), intent(out) :: fval

  select case ( vdf%x_dim_ ) 

    case (1)

      call interpolate_1d( vdf, fc, iand( fpos, 1 ), ix, x, np, int_type, fval )

    case (2)
      call interpolate_2d( vdf, fc, iand( fpos, 3 ), ix, x, np, int_type, fval )

    case (3)
      call interpolate_3d( vdf, fc, fpos, ix, x, np, int_type, fval )

  end select

end subroutine interpolate_vdf
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine interpolate_1D( vdf, fc, fpos, ix, x, np, int_type, fval )

  implicit none

  type( t_vdf ), intent(in) :: vdf
  integer, intent(in) :: fc
  integer, intent(in) :: fpos

  integer, dimension(:,:), intent(in) :: ix
  real( p_k_part ), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np
  integer, intent(in) :: int_type

  real( p_k_fld ), dimension(:), intent(out) :: fval

  real( p_k_fld ) :: hpth
  real( p_k_fld ), dimension( p_buf_size ) :: dx1
  integer, dimension( p_buf_size ) :: i1
  integer :: i, l, lnp, idx
  integer :: i1h, di1
  real( p_k_fld ), dimension( 5 ) :: w1

  if ( mod( int_type, 2 ) == 1 ) then
    hpth = 0.5
  else
    hpth = 0.0
  endif

  idx = 1
  do
    if ( idx > np ) exit

    if ( idx + p_buf_size - 1 > np ) then
      lnp = np - idx + 1
    else
      lnp = p_buf_size
    endif

    l = idx

    ! correct positions for field position
    if ( fpos == 0 ) then
      ! field position (0)
      do i = 1, lnp
        i1(i)  = ix(1,l)
        dx1(i) = x(1,l)
        l = l+1
      enddo
    else
      ! field position (0.5,0)
      do i = 1, lnp
        if ( x(1,l) < hpth ) then
          i1h = -1
        else
          i1h = 0
        endif

        i1(i)  = ix(1,l) + i1h
        dx1(i) = x(1,l) - 0.5 - i1h
        l = l+1
      enddo
    endif

    l = idx

    select case ( int_type )
      case( 1 )
        do i = 1, lnp
          call spline_n1( dx1(i), w1 )

          fval( l ) = 0.0
          do di1 = 0, 1
            fval( l ) = fval( l ) + vdf%f1( fc, i1(i) + di1 ) * w1( di1 + 1 )
          enddo

          l = l+1
        enddo

      case( 2 )
        do i = 1, lnp
          call spline_n2( dx1(i), w1 )

          fval( l ) = 0.0
          do di1 = -1, 1
            fval( l ) = fval( l ) + vdf%f1( fc, i1(i) + di1 ) * w1( di1 + 2 )
          enddo

          l = l+1
        enddo

      case( 3 )
        do i = 1, lnp
          call spline_n3( dx1(i), w1 )

          fval( l ) = 0.0
          do di1 = -1, 2
            fval( l ) = fval( l ) + vdf%f1( fc, i1(i) + di1 ) * w1( di1 + 2 )
          enddo

          l = l+1
        enddo

      case( 4 )
        do i = 1, lnp
          call spline_n4( dx1(i), w1 )

          fval( l ) = 0.0
          do di1 = -2, 2
            fval( l ) = fval( l ) + vdf%f1( fc, i1(i) + di1 ) * w1( di1 + 3 )
          enddo

          l = l+1
        enddo
    end select

    idx = idx + p_buf_size
  enddo

end subroutine interpolate_1D
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine interpolate_2D( vdf, fc, fpos, ix, x, np, int_type, fval )

  implicit none

  type( t_vdf ), intent(in) :: vdf
  integer, intent(in) :: fc
  integer, intent(in) :: fpos

  integer, dimension(:,:), intent(in) :: ix
  real( p_k_part ), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np
  integer, intent(in) :: int_type

  real( p_k_fld ), dimension(:), intent(out) :: fval

  real( p_k_fld ) :: hpth, line
  real( p_k_fld ), dimension( p_buf_size ) :: dx1, dx2
  integer, dimension( p_buf_size ) :: i1, i2
  integer :: i, l, lnp, idx
  integer :: i1h, i2h, di1, di2
  real( p_k_fld ), dimension( 5 ) :: w1, w2

  if ( mod( int_type, 2 ) == 1 ) then
    hpth = 0.5
  else
    hpth = 0.0
  endif

  idx = 1
  do
    if ( idx > np ) exit

    if ( idx + p_buf_size - 1 > np ) then
      lnp = np - idx + 1
    else
      lnp = p_buf_size
    endif

    l = idx

    ! correct positions for field position
    select case ( fpos )
      case(0) ! field position (0,0)
        do i = 1, lnp
          i1(i)  = ix(1,l)
          i2(i)  = ix(2,l)
          dx1(i) = x(1,l)
          dx2(i) = x(2,l)
          l = l+1
        enddo

      case(1) ! field position (0.5,0)
        do i = 1, lnp
          if ( x(1,l) < hpth ) then
            i1h = -1
          else
            i1h = 0
          endif

          i1(i)  = ix(1,l) + i1h
          i2(i)  = ix(2,l)

          dx1(i) = x(1,l) - 0.5 - i1h
          dx2(i) = x(2,l)
          l = l+1
        enddo

      case(2) ! field position (0.0,0.5)
        do i = 1, lnp
          if ( x(2,l) < hpth ) then
            i2h = -1
          else
            i2h = 0
          endif

          i1(i)  = ix(1,l)
          i2(i)  = ix(2,l) + i2h

          dx1(i) = x(1,l)
          dx2(i) = x(2,l) - 0.5 - i2h
          l = l+1
        enddo

      case(3)
        do i = 1, lnp
          if ( x(1,l) < hpth ) then
            i1h = -1
          else
            i1h = 0
          endif
          if ( x(2,l) < hpth ) then
            i2h = -1
          else
            i2h = 0
          endif

          i1(i)  = ix(1,l) + i1h
          i2(i)  = ix(2,l) + i2h

          dx1(i) = x(1,l) - 0.5 - i1h
          dx2(i) = x(2,l) - 0.5 - i2h
          l = l+1
        enddo

    end select

    l = idx

    select case ( int_type )
      case( 1 )
        do i = 1, lnp
          call spline_n1( dx1(i), w1 )
          call spline_n1( dx2(i), w2 )

          fval( l ) = 0.0
          do di2 = 0, 1
            line = 0.0
            do di1 = 0, 1
              line = line + vdf%f2( fc, i1(i) + di1, i2(i) + di2 ) * w1( di1 + 1 )
            enddo
            fval( l ) = fval( l ) + line * w2( di2 + 1 )
          enddo

          l = l+1
        enddo

      case( 2 )
        do i = 1, lnp
          call spline_n2( dx1(i), w1 )
          call spline_n2( dx2(i), w2 )

          fval( l ) = 0.0
          do di2 = -1, 1
            line = 0.0
            do di1 = -1, 1
              line = line + vdf%f2( fc, i1(i) + di1, i2(i) + di2 ) * w1( di1 + 2 )
            enddo
            fval( l ) = fval( l ) + line * w2( di2 + 2 )
          enddo

          l = l+1
        enddo

      case( 3 )
        do i = 1, lnp
          call spline_n3( dx1(i), w1 )
          call spline_n3( dx2(i), w2 )

          fval( l ) = 0.0
          do di2 = -1, 2
            line = 0.0
            do di1 = -1, 2
              line = line + vdf%f2( fc, i1(i) + di1, i2(i) + di2 ) * w1( di1 + 2 )
            enddo
            fval( l ) = fval( l ) + line * w2( di2 + 2 )
          enddo

          l = l+1
        enddo

      case( 4 )
        do i = 1, lnp
          call spline_n4( dx1(i), w1 )
          call spline_n4( dx2(i), w2 )

          fval( l ) = 0.0
          do di2 = -2, 2
            line = 0.0
            do di1 = -2, 2
              line = line + vdf%f2( fc, i1(i) + di1, i2(i) + di2 ) * w1( di1 + 3 )
            enddo
            fval( l ) = fval( l ) + line * w2( di2 + 3 )
          enddo

          l = l+1
        enddo
    end select

    idx = idx + p_buf_size
  enddo

end subroutine interpolate_2D
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine interpolate_3D( vdf, fc, fpos, ix, x, np, int_type, fval )

  implicit none

  type( t_vdf ), intent(in) :: vdf
  integer, intent(in) :: fc
  integer, intent(in) :: fpos

  integer, dimension(:,:), intent(in) :: ix
  real( p_k_part ), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np
  integer, intent(in) :: int_type

  real( p_k_fld ), dimension(:), intent(out) :: fval

  real( p_k_fld ) :: hpth, line, slice
  real( p_k_fld ), dimension( p_buf_size ) :: dx1, dx2, dx3
  integer, dimension( p_buf_size ) :: i1, i2, i3
  integer :: i, l, lnp, idx
  integer :: i1h, i2h, i3h, di1, di2, di3
  real( p_k_fld ), dimension( 5 ) :: w1, w2, w3

  if ( mod( int_type, 2 ) == 1 ) then
    hpth = 0.5
  else
    hpth = 0.0
  endif

  idx = 1
  do
    if ( idx > np ) exit

    if ( idx + p_buf_size - 1 > np ) then
      lnp = np - idx + 1
    else
      lnp = p_buf_size
    endif

    l = idx

    ! correct positions for field position
    select case ( fpos )
      case(0) ! field position (0,0,0)
        do i = 1, lnp
          i1(i)  = ix(1,l)
          i2(i)  = ix(2,l)
          i3(i)  = ix(3,l)
          dx1(i) = x(1,l)
          dx2(i) = x(2,l)
          dx3(i) = x(3,l)
          l = l+1
        enddo

      case(1) ! field position (0.5,0,0)
        do i = 1, lnp
          if ( x(1,l) < hpth ) then
            i1h = -1
          else
            i1h = 0
          endif

          i1(i)  = ix(1,l) + i1h
          i2(i)  = ix(2,l)
          i3(i)  = ix(3,l)

          dx1(i) = x(1,l) - 0.5 - i1h
          dx2(i) = x(2,l)
          dx3(i) = x(3,l)
          l = l+1
        enddo

      case(2) ! field position (0.0,0.5,0)
        do i = 1, lnp
          if ( x(2,l) < hpth ) then
            i2h = -1
          else
            i2h = 0
          endif

          i1(i)  = ix(1,l)
          i2(i)  = ix(2,l) + i2h
          i3(i)  = ix(3,l)

          dx1(i) = x(1,l)
          dx2(i) = x(2,l) - 0.5 - i2h
          dx3(i) = x(3,l)
          l = l+1
        enddo

      case(3) ! field position (0.5,0.5,0)
        do i = 1, lnp
          if ( x(1,l) < hpth ) then
            i1h = -1
          else
            i1h = 0
          endif
          if ( x(2,l) < hpth ) then
            i2h = -1
          else
            i2h = 0
          endif

          i1(i)  = ix(1,l) + i1h
          i2(i)  = ix(2,l) + i2h
          i3(i)  = ix(3,l)

          dx1(i) = x(1,l) - 0.5 - i1h
          dx2(i) = x(2,l) - 0.5 - i2h
          dx3(i) = x(3,l)
          l = l+1
        enddo

      case(4) ! field position (0,0,0.5)
        do i = 1, lnp
          if ( x(3,l) < hpth ) then
            i3h = -1
          else
            i3h = 0
          endif

          i1(i)  = ix(1,l)
          i2(i)  = ix(2,l)
          i3(i)  = ix(3,l) + i3h
          dx1(i) = x(1,l)
          dx2(i) = x(2,l)
          dx3(i) = x(3,l) - 0.5 - i3h
          l = l+1
        enddo

      case(5) ! field position (0.5,0,0.5)
        do i = 1, lnp
          if ( x(1,l) < hpth ) then
            i1h = -1
          else
            i1h = 0
          endif
          if ( x(3,l) < hpth ) then
            i3h = -1
          else
            i3h = 0
          endif

          i1(i)  = ix(1,l) + i1h
          i2(i)  = ix(2,l)
          i3(i)  = ix(3,l) + i3h

          dx1(i) = x(1,l) - 0.5 - i1h
          dx2(i) = x(2,l)
          dx3(i) = x(3,l) - 0.5 - i3h
          l = l+1
        enddo

      case(6) ! field position (0.0,0.5,0.5)
        do i = 1, lnp
          if ( x(2,l) < hpth ) then
            i2h = -1
          else
            i2h = 0
          endif
          if ( x(3,l) < hpth ) then
            i3h = -1
          else
            i3h = 0
          endif

          i1(i)  = ix(1,l)
          i2(i)  = ix(2,l) + i2h
          i3(i)  = ix(3,l) + i3h

          dx1(i) = x(1,l)
          dx2(i) = x(2,l) - 0.5 - i2h
          dx3(i) = x(3,l) - 0.5 - i3h
          l = l+1
        enddo

      case(7) ! field position (0.5,0.5,0.5)
        do i = 1, lnp
          if ( x(1,l) < hpth ) then
            i1h = -1
          else
            i1h = 0
          endif
          if ( x(2,l) < hpth ) then
            i2h = -1
          else
            i2h = 0
          endif
          if ( x(3,l) < hpth ) then
            i3h = -1
          else
            i3h = 0
          endif

          i1(i)  = ix(1,l) + i1h
          i2(i)  = ix(2,l) + i2h
          i3(i)  = ix(3,l) + i3h

          dx1(i) = x(1,l) - 0.5 - i1h
          dx2(i) = x(2,l) - 0.5 - i2h
          dx3(i) = x(3,l) - 0.5 - i3h
          l = l+1
        enddo
    end select

    l = idx

    select case ( int_type )
      case( 1 )
        do i = 1, lnp
          call spline_n1( dx1(i), w1 )
          call spline_n1( dx2(i), w2 )
          call spline_n1( dx3(i), w3 )

          fval( l ) = 0.0
          do di3 = 0, 1
            slice = 0.0
            do di2 = 0, 1
              line = 0.0
              do di1 = 0, 1
                line = line + vdf%f3(fc, i1(i)+di1, i2(i)+di2, i3(i)+di3) * w1( di1 + 1 )
              enddo
              slice = slice + line * w2( di2 + 1 )
            enddo
            fval( l ) = fval( l ) + slice * w3( di3 + 1 )
          enddo

          l = l+1
        enddo

      case( 2 )
        do i = 1, lnp
          call spline_n2( dx1(i), w1 )
          call spline_n2( dx2(i), w2 )
          call spline_n2( dx3(i), w3 )

          fval( l ) = 0.0
          do di3 = -1, 1
            slice = 0.0
            do di2 = -1, 1
              line = 0.0
              do di1 = -1, 1
                line = line + vdf%f3(fc, i1(i)+di1, i2(i)+di2, i3(i)+di3) * w1( di1 + 2 )
              enddo
              slice = slice + line * w2( di2 + 2 )
            enddo
            fval( l ) = fval( l ) + slice * w3( di3 + 2 )
          enddo

          l = l+1
        enddo

      case( 3 )
        do i = 1, lnp
          call spline_n3( dx1(i), w1 )
          call spline_n3( dx2(i), w2 )
          call spline_n3( dx3(i), w3 )

          fval( l ) = 0.0
          do di3 = -1, 2
            slice = 0.0
            do di2 = -1, 2
              line = 0.0
              do di1 = -1, 2
                line = line + vdf%f3(fc, i1(i)+di1, i2(i)+di2, i3(i)+di3) * w1( di1 + 2 )
              enddo
              slice = slice + line * w2( di2 + 2 )
            enddo
            fval( l ) = fval( l ) + slice * w3( di3 + 2 )
          enddo

          l = l+1
        enddo

      case( 4 )
        do i = 1, lnp
          call spline_n4( dx1(i), w1 )
          call spline_n4( dx2(i), w2 )
          call spline_n4( dx3(i), w3 )

          fval( l ) = 0.0
          do di3 = -2, 2
            slice = 0.0
            do di2 = -2, 2
              line = 0.0
              do di1 = -2, 2
                line = line + vdf%f3(fc, i1(i)+di1, i2(i)+di2, i3(i)+di3) * w1( di1 + 3 )
              enddo
              slice = slice + line * w2( di2 + 3 )
            enddo
            fval( l ) = fval( l ) + slice * w3( di3 + 3 )
          enddo

          l = l+1
        enddo
    end select

    idx = idx + p_buf_size
  enddo

end subroutine interpolate_3D
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! B-spline calculation routines
!   - These expect the position dx to be in the range [ -0.5, 0.5 [ for all interpolation levels
!---------------------------------------------------------------------------------------------------

subroutine spline_n1( dx, w )

  implicit none

  real( p_k_fld ), intent(in) :: dx
  real( p_k_fld ), dimension(:), intent(out) :: w

  w(1) = 0.5_p_k_fld - dx
  w(2) = 0.5_p_k_fld + dx

end subroutine spline_n1

subroutine spline_n2( dx, w )

  implicit none

  real( p_k_fld ), intent(in) :: dx
  real( p_k_fld ), dimension(:), intent(out) :: w

  w(1) = 0.5_p_k_fld*(0.5_p_k_fld - dx)**2
  w(2) = 0.75_p_k_fld - dx**2
  w(3) = 0.5_p_k_fld*(0.5_p_k_fld + dx)**2

end subroutine spline_n2

subroutine spline_n3( dx, w )

  implicit none

  real( p_k_fld ), intent(in) :: dx
  real( p_k_fld ), dimension(:), intent(out) :: w

  w(1) = -(-0.5 + dx)**3/6.
  w(2) = (4 - 6*(0.5 + dx)**2 + 3*(0.5 + dx)**3)/6.
  w(3) = (23 + 30*dx - 12*dx**2 - 24*dx**3)/48.
  w(4) = (0.5 + dx)**3/6.

end subroutine spline_n3

subroutine spline_n4( dx, w )

  implicit none

  real( p_k_fld ), intent(in) :: dx
  real( p_k_fld ), dimension(:), intent(out) :: w

  w(1) = (1 - 2*dx)**4/384.
  w(2) = (19 - 44*dx + 24*dx**2 + 16*dx**3 - 16*dx**4)/96.
  w(3) = 0.5989583333333334 - (5*dx**2)/8. + dx**4/4.
  w(4) = (19 + 44*dx + 24*dx**2 - 16*dx**3 - 16*dx**4)/96.
  w(5) = (1 + 2*dx)**4/384.

end subroutine spline_n4


end module m_vdf_interpolate
