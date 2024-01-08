!*****************************************************************************************
! Field interpolation routines
!  
!*****************************************************************************************


! if __TEMPLATE__ is defined then read the template definition at the end of the file
#ifndef __TEMPLATE__


#include "os-preprocess.fpp"
#include "os-config.h"

module m_emf_interpolate

use m_fparser
use m_emf_define
use m_logprof
use m_vdf_define
use m_parameters

implicit none

private


interface get_emf
   module procedure get_emf_cell
end interface

public :: get_emf


contains

!-----------------------------------------------------------------------------------------
! Interpolate fields at particle positions for cell based positions. The fields pointed to by
! this%e_part and this%b_part already include smoothed and/or external fields
!-----------------------------------------------------------------------------------------
subroutine get_emf_cell( this, bp, ep, ix, x, np, interpolation )

  implicit none

  class( t_emf ), intent(in), target :: this
  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  
  integer, intent(in) :: np
  
  integer, intent(in) :: interpolation

  type( t_vdf ), pointer :: e, b
        
  e => this%e_part
  b => this%b_part

  if ( this%part_grid_center ) then
    
    ! Interpolate fields centered at the cell corner 
    
    select case (interpolation)
       case(p_linear)
          select case ( p_x_dim )
            case (1)
               call get_emf_1d_cs1( bp, ep, b, e, ix, x, np )
            case (2)
               call get_emf_2d_cs1( bp, ep, b, e, ix, x, np )
            case (3)
               call get_emf_3d_cs1( bp, ep, b, e, ix, x, np )
          end select

       case(p_quadratic)
          select case ( p_x_dim )
            case (1)
               call get_emf_1d_cs2( bp, ep, b, e, ix, x, np )
            case (2)
               call get_emf_2d_cs2( bp, ep, b, e, ix, x, np )
            case (3)
               call get_emf_3d_cs2( bp, ep, b, e, ix, x, np )
          end select

       case(p_cubic)
          select case ( p_x_dim )
            case (1)
               call get_emf_1d_cs3( bp, ep, b, e, ix, x, np )
            case (2)
               call get_emf_2d_cs3( bp, ep, b, e, ix, x, np )
            case (3)
               call get_emf_3d_cs3( bp, ep, b, e, ix, x, np )
          end select

       case(p_quartic)
          select case ( p_x_dim )
            case (1)
               call get_emf_1d_cs4( bp, ep, b, e, ix, x, np )
            case (2)
               call get_emf_2d_cs4( bp, ep, b, e, ix, x, np )
            case (3)
               call get_emf_3d_cs4( bp, ep, b, e, ix, x, np )
          end select

       case default
          ERROR('Not implemented yet')
          call abort_program( p_err_notimplemented )
    end select
  
  else

    ! Interpolate fields on a staggered Yee mesh 
  
    select case (interpolation)
       case(p_linear)
          select case ( p_x_dim )
            case (1)
               call get_emf_1d_s1( bp, ep, b, e, ix, x, np )
            case (2)
               call get_emf_2d_s1( bp, ep, b, e, ix, x, np )
            case (3)
               call get_emf_3d_s1( bp, ep, b, e, ix, x, np )
          end select

       case(p_quadratic)
          select case ( p_x_dim )
            case (1)
               call get_emf_1d_s2( bp, ep, b, e, ix, x, np )
            case (2)
               call get_emf_2d_s2( bp, ep, b, e, ix, x, np )
            case (3)
               call get_emf_3d_s2( bp, ep, b, e, ix, x, np )
          end select

       case(p_cubic)
          select case ( p_x_dim )
            case (1)
               call get_emf_1d_s3( bp, ep, b, e, ix, x, np )
            case (2)
               call get_emf_2d_s3( bp, ep, b, e, ix, x, np )
            case (3)
               call get_emf_3d_s3( bp, ep, b, e, ix, x, np )
          end select

       case(p_quartic)
          select case ( p_x_dim )
            case (1)
               call get_emf_1d_s4( bp, ep, b, e, ix, x, np )
            case (2)
               call get_emf_2d_s4( bp, ep, b, e, ix, x, np )
            case (3)
               call get_emf_3d_s4( bp, ep, b, e, ix, x, np )
          end select

       case default
          ERROR('Not implemented yet')
          call abort_program( p_err_notimplemented )
    end select
  
  endif
  
end subroutine get_emf_cell
!-----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! Splines
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

subroutine splineh_s1( x, h, s )
  
  implicit none
  
  real( p_k_fld ), intent(in) :: x
  integer, intent(in) :: h
  real( p_k_fld ), dimension(0:1), intent(inout) :: s
    
  s(0) = ( 1 - h ) - x
  s(1) = (     h ) + x

end subroutine splineh_s1

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

subroutine splineh_s2( x, h, s )
  
  implicit none
  
  real( p_k_fld ), intent(in) :: x
  integer, intent(in) :: h
  real( p_k_fld ), dimension(-1:1), intent(inout) :: s
  
  real( p_k_fld ) :: t0, t1
  
  t0 = ( 1 - h ) - x
  t1 = (     h ) + x
  
  s(-1) = 0.5_p_k_fld * t0**2
  s( 0) = 0.5_p_k_fld + t0*t1
  s( 1) = 0.5_p_k_fld * t1**2

end subroutine splineh_s2

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

subroutine splineh_s3( x, h, s )
  
  implicit none
  
  real( p_k_fld ), intent(in) :: x
  integer, intent(in) :: h
  real( p_k_fld ), dimension(-1:2), intent(inout) :: s
  
  real( p_k_fld ) :: t0, t1, t2, t3
  real( p_k_fld ), parameter :: c1_6 = 1.0_p_k_fld / 6.0_p_k_fld 
  real( p_k_fld ), parameter :: c2_3 = 2.0_p_k_fld / 3.0_p_k_fld 
  
  t0 = ( 1 - h ) - x
  t1 = (     h ) + x
  
  t2 = t0 * t0
  t3 = t1 * t1
  
  t0 = t0 * t2
  t1 = t1 * t3 

  s(-1) = c1_6 * t0
  s( 0) = ( c2_3 - t3 ) + 0.5_p_k_fld*t1
  s( 1) = ( c2_3 - t2 ) + 0.5_p_k_fld*t0
  s( 2) = c1_6 * t1
  
end subroutine splineh_s3

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

subroutine splineh_s4( x, h, s )
  
  implicit none
  
  real( p_k_fld ), intent(in) :: x
  integer, intent(in) :: h
  real( p_k_fld ), dimension(-2:2), intent(inout) :: s

  real( p_k_fld ) :: t0, t1, t2

  real( p_k_fld ), parameter :: c1_6  = 1.0_p_k_fld / 6.0_p_k_fld 
  real( p_k_fld ), parameter :: c1_24 = 1.0_p_k_fld / 24.0_p_k_fld 
  real( p_k_fld ), parameter :: c11_24 = 11.0_p_k_fld / 24.0_p_k_fld 
  
  t0 = ( 1 - h ) - x
  t1 = (     h ) + x
  t2 = t0 * t1
  
!  s(-2) = c1_24 * t0**4
  s(-2) = c1_24 * ((t0**2)*(t0**2))
  s(-1) = c1_6 * ( 0.25_p_k_fld + t0 * ( 1.0_p_k_fld + t0 * ( 1.5_p_k_fld + t2 ) ) )
  s( 0) = c11_24 + t2 * ( 0.5_p_k_fld + 0.25_p_k_fld * t2 )
  s( 1) = c1_6 * ( 0.25_p_k_fld + t1 * ( 1.0_p_k_fld + t1 * ( 1.5_p_k_fld + t2 ) ) )
!  s( 2) = c1_24 * t1**4
  s( 2) = c1_24 * ((t1**2)*(t1**2))

end subroutine splineh_s4


!-----------------------------------------------------------------------------------------
! Signbit function
!-----------------------------------------------------------------------------------------
function signbit(x)

  implicit none
  real(p_k_fld), intent(in) :: x
  
  integer :: signbit
  
  if ( x < 0 ) then
    signbit = 1
  else
    signbit = 0
  endif
  
end function signbit
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!  Generate specific template functions for linear, quadratic, cubic and quartic
!  interpolation levels.
!-----------------------------------------------------------------------------------------

#define __TEMPLATE__


!********************************** Linear interpolation ********************************

#define GET_EMF_1D   get_emf_1d_s1
#define GET_EMF_2D   get_emf_2d_s1
#define GET_EMF_3D   get_emf_3d_s1

#define GET_EMF_1D_CS   get_emf_1d_cs1
#define GET_EMF_2D_CS   get_emf_2d_cs1
#define GET_EMF_3D_CS   get_emf_3d_cs1

#define SPLINE spline_s1
#define SPLINEH splineh_s1

! Lower point
#define LP 0
! Upper point
#define UP 1

#include __FILE__

!********************************** Quadratic interpolation ********************************

#define GET_EMF_1D   get_emf_1d_s2
#define GET_EMF_2D   get_emf_2d_s2
#define GET_EMF_3D   get_emf_3d_s2

#define GET_EMF_1D_CS   get_emf_1d_cs2
#define GET_EMF_2D_CS   get_emf_2d_cs2
#define GET_EMF_3D_CS   get_emf_3d_cs2

#define SPLINE spline_s2
#define SPLINEH splineh_s2

! Lower point
#define LP -1
! Upper point
#define UP 1

#include __FILE__

!********************************** Cubic interpolation **********************************

#define GET_EMF_1D   get_emf_1d_s3
#define GET_EMF_2D   get_emf_2d_s3
#define GET_EMF_3D   get_emf_3d_s3

#define GET_EMF_1D_CS   get_emf_1d_cs3
#define GET_EMF_2D_CS   get_emf_2d_cs3
#define GET_EMF_3D_CS   get_emf_3d_cs3

#define SPLINE spline_s3
#define SPLINEH splineh_s3

! Lower point
#define LP -1
! Upper point
#define UP 2

#include __FILE__

!********************************** Quartic interpolation ********************************

#define GET_EMF_1D   get_emf_1d_s4
#define GET_EMF_2D   get_emf_2d_s4
#define GET_EMF_3D   get_emf_3d_s4

#define GET_EMF_1D_CS   get_emf_1d_cs4
#define GET_EMF_2D_CS   get_emf_2d_cs4
#define GET_EMF_3D_CS   get_emf_3d_cs4

#define SPLINE spline_s4
#define SPLINEH splineh_s4

! Lower point
#define LP -2
! Upper point
#define UP 2

#include __FILE__

!****************************************************************************************

end module m_emf_interpolate

#else

!*****************************************************************************************
!
!  Template function definitions for field interpolation
!
!*****************************************************************************************

!*****************************************************************************************
!* Interpolation routines for staggered Yee mesh
!*****************************************************************************************

!-----------------------------------------------------------------------------------------
subroutine GET_EMF_1D( bp, ep, b, e, ix, x, np )
!-----------------------------------------------------------------------------------------
!  calculates the values of the electric-magnetic field 
!  at the positions of the array x on a 1d grid, using ORDER spline 
!  interpolation
!---------------------------------------------------

  implicit none

  integer, parameter :: rank = 1

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  type( t_vdf ), intent(in) :: b, e

  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix

  integer, intent(in) :: np


  real(p_k_fld) :: dx1, f1, f2, f3
  real(p_k_fld), dimension(LP:UP) :: w1, w1h
  integer :: i1, i1h, l, k1, h


  do l = 1, np
     
     i1 = ix(1,l)
     dx1 = real( x(1,l), p_k_fld )
     
     h = signbit(dx1)
     i1h = i1 - h
     
     ! get spline weitghts for x
     call SPLINE( dx1, w1 )
     call SPLINEH( dx1, h, w1h )
     
     ! Interpolate Fields
     f1 = 0
     f2 = 0
     f3 = 0
     
     do k1 = LP, UP
       f1 = f1 + e%f1(1,i1h + k1) * w1h(k1)
       f2 = f2 + e%f1(2,i1  + k1) * w1(k1)
       f3 = f3 + e%f1(3,i1  + k1) * w1(k1)
     enddo
     
     ep( 1, l ) = f1
     ep( 2, l ) = f2
     ep( 3, l ) = f3

     f1 = 0
     f2 = 0
     f3 = 0
     
     do k1 = LP, UP
       f1 = f1 + b%f1(1,i1  + k1) * w1(k1)
       f2 = f2 + b%f1(2,i1h + k1) * w1h(k1)
       f3 = f3 + b%f1(3,i1h + k1) * w1h(k1)
     enddo
     
     bp( 1, l ) = f1
     bp( 2, l ) = f2
     bp( 3, l ) = f3
     
     ! end of automatic code
  enddo

end subroutine GET_EMF_1D
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine GET_EMF_2D( bp, ep, b, e, ix, x, np )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep
  type( t_vdf ), intent(in) :: b, e
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  integer :: i1, i2, i1h, i2h, l, k1, k2, h1, h2
  real(p_k_fld) :: dx1, dx2
  real(p_k_fld), dimension(LP:UP) :: w1, w1h, w2, w2h
  
  real(p_k_fld) :: f1, f2, f3, f1line, f2line, f3line
  
  
  do l = 1, np
          
     i1 = ix(1,l)
     i2 = ix(2,l)
     dx1 = real( x(1,l), p_k_fld )
     dx2 = real( x(2,l), p_k_fld )
 
     h1 = signbit(dx1)
     h2 = signbit(dx2)
     
     i1h = i1 - h1
     i2h = i2 - h2
          
     ! get spline weitghts for x and y
     call SPLINE( dx1, w1 )
     call SPLINEH( dx1, h1, w1h )

     call SPLINE( dx2, w2 )
     call SPLINEH( dx2, h2, w2h )
     
     ! Interpolate E - Field
     f1 = 0
     f2 = 0
     f3 = 0
     
     do k2 = LP, UP
       f1line = 0
       f2line = 0
       f3line = 0
       
       do k1 = LP, UP
         f1line = f1line + e%f2(1,i1h + k1, i2  + k2) * w1h(k1)
         f2line = f2line + e%f2(2,i1  + k1, i2h + k2) * w1(k1)
         f3line = f3line + e%f2(3,i1  + k1, i2  + k2) * w1(k1)
       enddo
       
       f1 = f1 + f1line * w2(k2)
       f2 = f2 + f2line * w2h(k2)
       f3 = f3 + f3line * w2(k2)
     enddo
     
     ep( 1, l ) = f1
     ep( 2, l ) = f2
     ep( 3, l ) = f3

     ! Interpolate B - Field

     f1 = 0
     f2 = 0
     f3 = 0
     
     do k2 = LP, UP
       f1line = 0
       f2line = 0
       f3line = 0
       
       do k1 = LP, UP
         f1line = f1line + b%f2(1,i1  + k1, i2h + k2) * w1(k1)
         f2line = f2line + b%f2(2,i1h + k1, i2  + k2) * w1h(k1)
         f3line = f3line + b%f2(3,i1h + k1, i2h + k2) * w1h(k1)
       enddo
       
       f1 = f1 + f1line * w2h(k2)
       f2 = f2 + f2line * w2(k2)
       f3 = f3 + f3line * w2h(k2)
     enddo
     
     bp( 1, l ) = f1
     bp( 2, l ) = f2
     bp( 3, l ) = f3

  enddo


end subroutine GET_EMF_2D
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine GET_EMF_3D( bp, ep, b, e, ix, x, np )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  integer, parameter :: rank = 3

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  type( t_vdf ), intent(in) :: b, e

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  ! local variables
  real(p_k_fld) :: dx1, dx2, dx3
  integer :: h1, h2, h3
  real(p_k_fld), dimension(LP:UP) :: w1, w1h, w2, w2h,  w3, w3h
  real(p_k_fld) :: f1, f2, f3 
  real(p_k_fld) :: f1line, f2line, f3line
  real(p_k_fld) :: f1plane, f2plane, f3plane

  integer :: i1, i2, i3, i1h, i2h, i3h, l
  integer :: k1, k2, k3

  do l = 1, np
     
     i1 = ix(1,l)
     i2 = ix(2,l)
     i3 = ix(3,l)
     dx1 = real( x(1,l), p_k_fld )
     dx2 = real( x(2,l), p_k_fld )
     dx3 = real( x(3,l), p_k_fld )

     h1 = signbit(dx1)
     h2 = signbit(dx2)
     h3 = signbit(dx3)

     i1h = i1 - h1
     i2h = i2 - h2
     i3h = i3 - h3
     
     ! get spline weights for x,y and z
     call SPLINE( dx1, w1 )
     call SPLINEH( dx1, h1, w1h )

     call SPLINE( dx2, w2 )
     call SPLINEH( dx2, h2, w2h )

     call SPLINE( dx3, w3 )
     call SPLINEH( dx3, h3, w3h )

     ! Interpolate E-Field
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
           f1line = f1line + e%f3(1,i1h + k1, i2  + k2, i3  + k3) * w1h(k1)
           f2line = f2line + e%f3(2,i1  + k1, i2h + k2, i3  + k3) * w1(k1)
           f3line = f3line + e%f3(3,i1  + k1, i2  + k2, i3h + k3) * w1(k1)
         enddo
      
         f1plane = f1plane + f1line * w2(k2)
         f2plane = f2plane + f2line * w2h(k2)
         f3plane = f3plane + f3line * w2(k2)
       enddo

       f1 = f1 + f1plane * w3(k3)
       f2 = f2 + f2plane * w3(k3)
       f3 = f3 + f3plane * w3h(k3)
     enddo
     
     ep( 1, l ) = f1
     ep( 2, l ) = f2
     ep( 3, l ) = f3

     ! Interpolate B-Field
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
           f1line = f1line + b%f3(1,i1  + k1, i2h + k2, i3h + k3) * w1(k1)
           f2line = f2line + b%f3(2,i1h + k1, i2  + k2, i3h + k3) * w1h(k1)
           f3line = f3line + b%f3(3,i1h + k1, i2h + k2, i3  + k3) * w1h(k1)
         enddo
      
         f1plane = f1plane + f1line * w2h(k2)
         f2plane = f2plane + f2line * w2(k2)
         f3plane = f3plane + f3line * w2h(k2)
       enddo

       f1 = f1 + f1plane * w3h(k3)
       f2 = f2 + f2plane * w3h(k3)
       f3 = f3 + f3plane * w3(k3)
     enddo
     
     bp( 1, l ) = f1
     bp( 2, l ) = f2
     bp( 3, l ) = f3
     
  enddo


end subroutine GET_EMF_3D
!-----------------------------------------------------------------------------------------

!*****************************************************************************************
!* Interpolation routines for field values defined on the cell corner
!*****************************************************************************************

!-----------------------------------------------------------------------------------------
subroutine GET_EMF_1D_CS( bp, ep, b, e, ix, x, np )
!-----------------------------------------------------------------------------------------
!  calculates the values of the electric-magnetic field 
!  at the positions of the array x on a 1d grid, using ORDER spline 
!  interpolation
!---------------------------------------------------

  implicit none

  integer, parameter :: rank = 1

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  type( t_vdf ), intent(in) :: b, e

  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, dimension(:,:), intent(in) :: ix

  integer, intent(in) :: np


  real(p_k_fld) :: dx1, f1, f2, f3
  real(p_k_fld), dimension(LP:UP) :: w1
  integer :: l, i1, k1


  do l = 1, np
     
     i1 = ix(1,l)
     dx1 = real( x(1,l), p_k_fld )
     
     ! get spline weitghts for x
     call SPLINE( dx1, w1 )
     
     ! Interpolate Fields
     f1 = 0
     f2 = 0
     f3 = 0
     
     do k1 = LP, UP
       f1 = f1 + e%f1(1,i1 + k1) * w1(k1)
       f2 = f2 + e%f1(2,i1 + k1) * w1(k1)
       f3 = f3 + e%f1(3,i1 + k1) * w1(k1)
     enddo
     
     ep( 1, l ) = f1
     ep( 2, l ) = f2
     ep( 3, l ) = f3

     f1 = 0
     f2 = 0
     f3 = 0
     
     do k1 = LP, UP
       f1 = f1 + b%f1(1,i1 + k1) * w1(k1)
       f2 = f2 + b%f1(2,i1 + k1) * w1(k1)
       f3 = f3 + b%f1(3,i1 + k1) * w1(k1)
     enddo
     
     bp( 1, l ) = f1
     bp( 2, l ) = f2
     bp( 3, l ) = f3
     
  enddo

end subroutine GET_EMF_1D_CS
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine GET_EMF_2D_CS( bp, ep, b, e, ix, x, np )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep
  type( t_vdf ), intent(in) :: b, e
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  integer :: i1, i2, l, k1, k2
  real(p_k_fld) :: dx1, dx2
  real(p_k_fld), dimension(LP:UP) :: w1, w2
  
  real(p_k_fld) :: f1, f2, f3, f1line, f2line, f3line
  
  do l = 1, np
          
     i1 = ix(1,l)
     i2 = ix(2,l)
     dx1 = real( x(1,l), p_k_fld )
     dx2 = real( x(2,l), p_k_fld )
     
     ! get spline weitghts for x and y
     call SPLINE( dx1, w1 )
     call SPLINE( dx2, w2 )

     
     ! Interpolate E - Field
     f1 = 0
     f2 = 0
     f3 = 0
     
     do k2 = LP, UP
       f1line = 0
       f2line = 0
       f3line = 0
       
       do k1 = LP, UP
         f1line = f1line + e%f2(1,i1 + k1, i2 + k2) * w1(k1)
         f2line = f2line + e%f2(2,i1 + k1, i2 + k2) * w1(k1)
         f3line = f3line + e%f2(3,i1 + k1, i2 + k2) * w1(k1)
       enddo
       
       f1 = f1 + f1line * w2(k2)
       f2 = f2 + f2line * w2(k2)
       f3 = f3 + f3line * w2(k2)
     enddo
     
     ep( 1, l ) = f1
     ep( 2, l ) = f2
     ep( 3, l ) = f3

     ! Interpolate B - Field

     f1 = 0
     f2 = 0
     f3 = 0
     
     do k2 = LP, UP
       f1line = 0
       f2line = 0
       f3line = 0
       
       do k1 = LP, UP
         f1line = f1line + b%f2(1,i1 + k1, i2 + k2) * w1(k1)
         f2line = f2line + b%f2(2,i1 + k1, i2 + k2) * w1(k1)
         f3line = f3line + b%f2(3,i1 + k1, i2 + k2) * w1(k1)
       enddo
       
       f1 = f1 + f1line * w2(k2)
       f2 = f2 + f2line * w2(k2)
       f3 = f3 + f3line * w2(k2)
     enddo
     
     bp( 1, l ) = f1
     bp( 2, l ) = f2
     bp( 3, l ) = f3

  enddo


end subroutine GET_EMF_2D_CS
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine GET_EMF_3D_CS( bp, ep, b, e, ix, x, np )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  integer, parameter :: rank = 3

  real(p_k_part), dimension(:,:), intent( out ) :: bp, ep

  type( t_vdf ), intent(in) :: b, e

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  ! local variables
  real(p_k_fld) :: dx1, dx2, dx3
  real(p_k_fld), dimension(LP:UP) :: w1, w2, w3
  real(p_k_fld) :: f1, f2, f3 
  real(p_k_fld) :: f1line, f2line, f3line
  real(p_k_fld) :: f1plane, f2plane, f3plane

  integer :: i1, i2, i3, l, k1, k2, k3

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

     ! Interpolate E-Field
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
           f1line = f1line + e%f3(1,i1 + k1, i2 + k2, i3 + k3) * w1(k1)
           f2line = f2line + e%f3(2,i1 + k1, i2 + k2, i3 + k3) * w1(k1)
           f3line = f3line + e%f3(3,i1 + k1, i2 + k2, i3 + k3) * w1(k1)
         enddo
      
         f1plane = f1plane + f1line * w2(k2)
         f2plane = f2plane + f2line * w2(k2)
         f3plane = f3plane + f3line * w2(k2)
       enddo

       f1 = f1 + f1plane * w3(k3)
       f2 = f2 + f2plane * w3(k3)
       f3 = f3 + f3plane * w3(k3)
     enddo
     
     ep( 1, l ) = f1
     ep( 2, l ) = f2
     ep( 3, l ) = f3

     ! Interpolate B-Field
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
           f1line = f1line + b%f3(1,i1 + k1, i2 + k2, i3 + k3) * w1(k1)
           f2line = f2line + b%f3(2,i1 + k1, i2 + k2, i3 + k3) * w1(k1)
           f3line = f3line + b%f3(3,i1 + k1, i2 + k2, i3 + k3) * w1(k1)
         enddo
      
         f1plane = f1plane + f1line * w2(k2)
         f2plane = f2plane + f2line * w2(k2)
         f3plane = f3plane + f3line * w2(k2)
       enddo

       f1 = f1 + f1plane * w3(k3)
       f2 = f2 + f2plane * w3(k3)
       f3 = f3 + f3plane * w3(k3)
     enddo
     
     bp( 1, l ) = f1
     bp( 2, l ) = f2
     bp( 3, l ) = f3
     
  enddo


end subroutine GET_EMF_3D_CS
!-----------------------------------------------------------------------------------------


! Clear template definitions

#undef GET_EMF_1D
#undef GET_EMF_2D
#undef GET_EMF_3D

#undef GET_EMF_1D_CS
#undef GET_EMF_2D_CS
#undef GET_EMF_3D_CS

#undef SPLINE
#undef SPLINEH

#undef LP
#undef UP

#endif

