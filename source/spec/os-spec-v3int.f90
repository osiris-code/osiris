!*****************************************************************************************
! m_species_v3dep module
! 
! Handles interpolation of a vector 3 quantity at the cell corner
!*****************************************************************************************

! if __TEMPLATE__ is defined then read the template definition at the end of the file
#ifndef __TEMPLATE__


#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_v3int

use m_species_define

use m_system
use m_parameters
use m_vdf_define

private

interface interpolate_v3
  module procedure interpolate_v3_spec
end interface

public :: interpolate_v3

! -----------------------------------------------------------------------------
contains

! ----------------------------------------------------------------------------------------
! Splines
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! Linear
! ----------------------------------------------------------------------------------------
subroutine spline_s1( x, s )
  
  implicit none
  
  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(0:1), intent(out) :: s
  
  s(0) = 0.5 - x
  s(1) = 0.5 + x

end subroutine spline_s1

! ----------------------------------------------------------------------------------------
! Quadratic
! ----------------------------------------------------------------------------------------
subroutine spline_s2( x, s )
  
  implicit none
  
  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(-1:1), intent(out) :: s
  
  real( p_k_fld ) :: t0, t1
  
  t0 = 0.5 - x
  t1 = 0.5 + x
  
  s(-1) = 0.5 * t0**2
  s( 0) = 0.5 + t0*t1
  s( 1) = 0.5 * t1**2

end subroutine spline_s2

! ----------------------------------------------------------------------------------------
! Cubic
! ----------------------------------------------------------------------------------------

subroutine spline_s3( x, s )
  
  implicit none
  
  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(-1:2), intent(out) :: s
  
  real( p_k_fld ) :: t0, t1, t2, t3
  
  t0 = 0.5 - x
  t1 = 0.5 + x
  
  t2 = t0 * t0
  t3 = t1 * t1
  
  t0 = t0 * t2
  t1 = t1 * t3 

  s(-1) = t0/6.
  s( 0) = (0.6666666666666666_p_k_fld - t3 ) + 0.5*t1
  s( 1) = (0.6666666666666666_p_k_fld - t2 ) + 0.5*t0
  s( 2) = t1/6.
  
end subroutine spline_s3

! ----------------------------------------------------------------------------------------
! Quartic
! ----------------------------------------------------------------------------------------

subroutine spline_s4( x, s )
  
  implicit none
  
  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(-2:2), intent(out) :: s

  real( p_k_fld ) :: t0, t1, t2
  
  t0 = 0.5 - x
  t1 = 0.5 + x
  t2 = t0 * t1
    
  s(-2) = t0**4/24.
  s(-1) = ( 0.25 + t0 * ( 1 + t0 * ( 1.5 + t2 ) ) ) / 6.
  s( 0) = 0.4583333333333333_p_k_fld + t2 * ( 0.5 + 0.25 * t2 )
  s( 1) = ( 0.25 + t1 * ( 1 + t1 * ( 1.5 + t2 ) ) ) / 6.
  s( 2) = t1**4/24.

end subroutine spline_s4

!*****************************************************************************************
!*                                       Interface                                       *
!*****************************************************************************************

!-----------------------------------------------------------------------------------------
subroutine interpolate_v3_spec( this, f, i1, i2, fp )

  implicit none

  ! dummy variables
  class( t_species ), intent(in) :: this
  type( t_vdf ),     intent(in) :: f
  integer, intent(in) :: i1, i2
  
  real( p_k_part ), dimension(:,:), intent(out) :: fp
  
  integer :: np
  
  ! number of particles to deposit
  np = i2 - i1 + 1
    
  ! Interpolate v3 grid quantity f at particle position using the appropriate
  ! interpolation level
  
  select case ( this%interpolation )
     
     case (p_linear)
        select case ( p_x_dim )
         case (1)
           call interpolate_v3_1d_s1( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (2)
           call interpolate_v3_2d_s1( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (3)
           call interpolate_v3_3d_s1( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
        end select

     case (p_quadratic)
        select case ( p_x_dim )
         case (1)
           call interpolate_v3_1d_s2( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (2)
           call interpolate_v3_2d_s2( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (3)
           call interpolate_v3_3d_s2( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
        end select

     case (p_cubic)
        select case ( p_x_dim )
         case (1)
           call interpolate_v3_1d_s3( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (2)
           call interpolate_v3_2d_s3( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (3)
           call interpolate_v3_3d_s3( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
        end select

     case (p_quartic)
        select case ( p_x_dim )
         case (1)
           call interpolate_v3_1d_s4( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (2)
           call interpolate_v3_2d_s4( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
         case (3)
           call interpolate_v3_3d_s4( f, this%ix(:,i1:), this%x(:,i1:), np, fp )
        end select
     
     case default
        ERROR('Not implemented')
        call abort_program( p_err_notimplemented )

  end select
    

end subroutine interpolate_v3_spec
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!  Generate specific template functions for linear, quadratic, cubic and quartic
!  interpolation levels.
!-----------------------------------------------------------------------------------------

#define __TEMPLATE__

!********************************** Linear interpolation ********************************

#define INTERPOLATE_V3_1D   interpolate_v3_1d_s1
#define INTERPOLATE_V3_2D   interpolate_v3_2d_s1
#define INTERPOLATE_V3_3D   interpolate_v3_3d_s1

#define SPLINE spline_s1

! Lower point
#define LP 0
! Upper point
#define UP 1

#include __FILE__

!********************************** Quadratic interpolation ********************************

#define INTERPOLATE_V3_1D   interpolate_v3_1d_s2
#define INTERPOLATE_V3_2D   interpolate_v3_2d_s2
#define INTERPOLATE_V3_3D   interpolate_v3_3d_s2

#define SPLINE spline_s2

! Lower point
#define LP -1
! Upper point
#define UP 1

#include __FILE__

!********************************** Cubic interpolation ********************************

#define INTERPOLATE_V3_1D   interpolate_v3_1d_s3
#define INTERPOLATE_V3_2D   interpolate_v3_2d_s3
#define INTERPOLATE_V3_3D   interpolate_v3_3d_s3

#define SPLINE spline_s3

! Lower point
#define LP -1
! Upper point
#define UP 2

#include __FILE__

!********************************** Quartic interpolation ********************************

#define INTERPOLATE_V3_1D   interpolate_v3_1d_s4
#define INTERPOLATE_V3_2D   interpolate_v3_2d_s4
#define INTERPOLATE_V3_3D   interpolate_v3_3d_s4

#define SPLINE spline_s4

! Lower point
#define LP -2
! Upper point
#define UP 2

#include __FILE__

!*****************************************************************************************


end module m_species_v3int

#else

!*****************************************************************************************
!
!  Template function definitions for v3 quantity interpolation
!
!*****************************************************************************************


!-----------------------------------------------------------------------------------------
subroutine INTERPOLATE_V3_1D( f, ix, x, np, fp )
  
  implicit none

  integer, parameter :: rank = 1
  
  type( t_vdf ), intent(in) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld), dimension(:,:), intent(out) :: fp
  
  ! local variables
  integer :: l, i1, k1
  real( p_k_fld ) :: dx1, f1, f2, f3
    
  real( p_k_fld ), dimension(LP:UP) :: w1

  ! quantity is defined in the lower corner of the cell

  do l = 1, np
    
    i1 = ix(1,l)
    dx1 = real( x(1,l), p_k_fld )
     
    ! get spline weights for x 
    call SPLINE( dx1, w1 )
     
    ! Interpolate vector3 quantity
    f1 = 0
    f2 = 0
    f3 = 0

    do k1 = LP, UP
      f1 = f1 + f%f1(1,i1+k1) * w1(k1)
      f2 = f2 + f%f1(2,i1+k1) * w1(k1)
      f3 = f3 + f%f1(3,i1+k1) * w1(k1)
    enddo
    
    fp( 1, l ) = f1
    fp( 2, l ) = f2
    fp( 3, l ) = f3

  enddo

end subroutine INTERPOLATE_V3_1D
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine INTERPOLATE_V3_2D( f, ix, x, np, fp )
  
  implicit none

  integer, parameter :: rank = 2
  
  type( t_vdf ), intent(in) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld), dimension(:,:), intent(out) :: fp
  
  ! local variables
  integer :: l, i1, i2, k1, k2
  real( p_k_fld ) :: dx1, dx2
  
  real(p_k_fld) :: f1, f2, f3, f1line, f2line, f3line
  
  real( p_k_fld ), dimension(LP:UP) :: w1, w2

  ! quantity is defined in the lower corner of the cell

  do l = 1, np
    
    ! order 1 vector3 interpolation
    ! generated automatically by z-2.1
    
    i1 = ix(1,l)
    i2 = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
     
    ! get spline weights for x and y
    call SPLINE( dx1, w1 )
    call SPLINE( dx2, w2 )
     
    ! Interpolate vector3 quantity
    f1 = 0
    f2 = 0
    f3 = 0

    do k2 = LP, UP
      f1line = 0
      f2line = 0
      f3line = 0

      do k1 = LP, UP
        f1line = f1line + f%f2(1,i1+k1,i2+k1) * w1(k1)
        f2line = f2line + f%f2(2,i1+k1,i2+k1) * w1(k1)
        f3line = f3line + f%f2(3,i1+k1,i2+k1) * w1(k1)
      enddo

      f1 = f1 + f1line * w2( k2 )
      f2 = f2 + f2line * w2( k2 )
      f3 = f3 + f3line * w2( k2 )
    enddo
     
    fp(1,l) = f1
    fp(2,l) = f2
    fp(3,l) = f3

  enddo

end subroutine INTERPOLATE_V3_2D
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine INTERPOLATE_V3_3D( f, ix, x, np, fp )
  
  implicit none

  integer, parameter :: rank = 3
  
  type( t_vdf ), intent(in) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  real(p_k_fld), dimension(:,:), intent(out) :: fp
  
  ! local variables
  integer :: l, i1, i2, i3, k1, k2, k3
  real( p_k_fld ) :: dx1, dx2, dx3
  
  real(p_k_fld) :: f1, f2, f3 
  real(p_k_fld) :: f1line, f2line, f3line
  real(p_k_fld) :: f1plane, f2plane, f3plane
  
  real( p_k_fld ), dimension(LP:UP) :: w1, w2, w3

  ! quantity is defined in the lower corner of the cell

  do l = 1, np
    
    i1 = ix(1,l)
    i2 = ix(2,l)
    i3 = ix(3,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    dx3 = real( x(3,l), p_k_fld )
     
    ! get spline weights for x,y and z
    call SPLINE( dx1, w1 )
    call SPLINE( dx2, w2 )
    call SPLINE( dx3, w3 )

    ! Interpolate vector3 quantity
    f1 = 0
    f2 = 0
    f3 = 0
    
    do k3 = LP, UP
      f1plane = 0
      f2plane = 0
      f3plane = 0

      do k2 = LP, UP
        f1line = 0
        f2line = 0
        f3line = 0
      
        do k1 = LP, UP
          f1line = f1line + f%f3(1,i1 + k1, i2 + k2, i3 + k3) * w1(k1)
          f2line = f2line + f%f3(2,i1 + k1, i2 + k2, i3 + k3) * w1(k1)
          f3line = f3line + f%f3(3,i1 + k1, i2 + k2, i3 + k3) * w1(k1)
        enddo
      
        f1plane = f1plane + f1line * w2(k2)
        f2plane = f2plane + f2line * w2(k2)
        f3plane = f3plane + f3line * w2(k2)
      enddo

      f1 = f1 + f1plane * w3(k3)
      f2 = f2 + f2plane * w3(k3)
      f3 = f3 + f3plane * w3(k3)
    enddo
     
     fp( 1, l ) = f1
     fp( 2, l ) = f2
     fp( 3, l ) = f3     

  enddo

end subroutine INTERPOLATE_V3_3D
!-----------------------------------------------------------------------------------------

! Clear template definitions

#undef INTERPOLATE_V3_1D
#undef INTERPOLATE_V3_2D
#undef INTERPOLATE_V3_3D

#undef SPLINE

#undef LP
#undef UP

#endif




