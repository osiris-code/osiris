!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     module with general utility subroutines
!     that are utilized by other modules
!     each subroutine in this module should not depend on
!     any definitions external to that subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_utilities

use m_system

implicit none

private

interface swap_value
  module procedure swap_value_int_vector
  module procedure swap_value_single_vector
  module procedure swap_value_double_vector
  module procedure swap_value_logical_vector
  module procedure swap_value_int_array2d
  module procedure swap_value_single_array2d
  module procedure swap_value_double_array2d
end interface

interface heapsort
  module procedure heapsort_inplace
  module procedure heapsort_index
end interface

interface lorentz_transform
  module procedure lorentz_transform_single
  module procedure lorentz_transform_double
  module procedure lorentz_transform_0_single
  module procedure lorentz_transform_0_double
end interface

! declare things that should be public
public :: swap_value
public :: heapsort, lorentz_transform

contains

!---------------------------------------------------------------------------------------------------
! Swaps 2 integer values within a vector
!---------------------------------------------------------------------------------------------------
subroutine swap_value_int_vector( vector, i0, i1 )

  implicit none

  integer, dimension(:), intent(inout) :: vector
  integer, intent(in) :: i0, i1

  integer :: temp

  temp = vector(i0)
  vector(i0) = vector(i1)
  vector(i1) = temp

end subroutine swap_value_int_vector

!---------------------------------------------------------------------------------------------------
! Swaps 2 single values within a vector
!---------------------------------------------------------------------------------------------------
subroutine swap_value_single_vector( vector, i0, i1 )

  implicit none

  real(p_single), dimension(:), intent(inout) :: vector
  integer, intent(in) :: i0, i1

  real(p_single) :: temp

  temp = vector(i0)
  vector(i0) = vector(i1)
  vector(i1) = temp

end subroutine swap_value_single_vector

!---------------------------------------------------------------------------------------------------
! Swaps 2 double values within a vector
!---------------------------------------------------------------------------------------------------
subroutine swap_value_double_vector( vector, i0, i1 )

  implicit none

  real(p_double), dimension(:), intent(inout) :: vector
  integer, intent(in) :: i0, i1

  real(p_double) :: temp

  temp = vector(i0)
  vector(i0) = vector(i1)
  vector(i1) = temp

end subroutine swap_value_double_vector

!---------------------------------------------------------------------------------------------------
! Swaps 2 integer values within a vector
!---------------------------------------------------------------------------------------------------
subroutine swap_value_logical_vector( vector, i0, i1 )

  implicit none

  logical, dimension(:), intent(inout) :: vector
  integer, intent(in) :: i0, i1

  logical :: temp

  temp = vector(i0)
  vector(i0) = vector(i1)
  vector(i1) = temp

end subroutine swap_value_logical_vector


!---------------------------------------------------------------------------------------------------
! Swaps 2 integer rows/cols within a 2D array
!---------------------------------------------------------------------------------------------------
subroutine swap_value_int_array2d( array, diridx, i0, i1 )

  implicit none

  integer, dimension(:,:), intent(inout) :: array
  integer, intent(in) :: dirIdx, i0, i1

  integer :: i, temp

  if ( dirIdx == 2 ) then
    do i = lbound(array, 1), ubound(array, 1)
      temp = array(i, i0)
      array(i, i0 ) = array(i, i1)
      array(i, i1 ) = temp
    enddo
  else
    do i = lbound(array, 2), ubound(array, 2)
      temp = array(i0,i)
      array(i0,i) = array(i1,i)
      array(i1,i) = temp
    enddo
  endif

end subroutine swap_value_int_array2d

!---------------------------------------------------------------------------------------------------
! Swaps 2 single prec. rows/cols within a 2D array
!---------------------------------------------------------------------------------------------------
subroutine swap_value_single_array2d( array, dirIdx, i0, i1 )

  implicit none

  real(p_single), dimension(:,:), intent(inout) :: array
  integer, intent(in) :: dirIdx, i0, i1

  integer :: i
  real(p_single) :: temp

  if ( dirIdx == 2 ) then
    do i = lbound(array, 1), ubound(array, 1)
      temp = array(i, i0)
      array(i, i0 ) = array(i, i1)
      array(i, i1 ) = temp
    enddo
  else
    do i = lbound(array, 2), ubound(array, 2)
      temp = array(i0,i)
      array(i0,i) = array(i1,i)
      array(i1,i) = temp
    enddo
  endif

end subroutine swap_value_single_array2d


!---------------------------------------------------------------------------------------------------
! Swaps 2 double prec. rows/cols within a 2D array
!---------------------------------------------------------------------------------------------------
subroutine swap_value_double_array2d( array, dirIdx, i0, i1 )

  implicit none

  real(p_double), dimension(:,:), intent(inout) :: array
  integer, intent(in) :: dirIdx, i0, i1

  integer :: i
  real(p_double) :: temp

  if ( dirIdx == 2 ) then
    do i = lbound(array, 1), ubound(array, 1)
      temp = array(i, i0)
      array(i, i0 ) = array(i, i1)
      array(i, i1 ) = temp
    enddo
  else
    do i = lbound(array, 2), ubound(array, 2)
      temp = array(i0,i)
      array(i0,i) = array(i1,i)
      array(i1,i) = temp
    enddo
  endif

end subroutine swap_value_double_array2d



!---------------------------------------------------------------------------------------------------
! heapsort integer array in place
! Adapted from Numerical Recipes in C, 2nd edition
! this needs to be fixed to use 1: based arrays
!---------------------------------------------------------------------------------------------------
subroutine heapsort_inplace( list )

  implicit none

  integer, dimension(0:), intent(inout) :: list

  integer :: t
  integer :: n, i, parent, child

  n = size(list)

  if ( n<2 ) return

  i = n/2

  do
    if ( i > 0 ) then
      i = i-1
      t = list(i)
    else
      n = n-1
      if ( n == 0 ) return
      t = list(n)
      list(n) = list(0)
    endif

    parent = i
    child = i*2+1

    do
      if (child >= n) exit
      if (child +1 < n ) then
         if (list(child) < list(child+1)) child = child+1
      endif

      if ( t < list(child) ) then
        list(parent) = list(child)
        parent = child
        child = parent*2 + 1
      else
        exit
      endif
    enddo

    list(parent) = t
  enddo


end subroutine heapsort_inplace

!---------------------------------------------------------------------------------------------------
! index integer array using heapsort (store indexes of sorted array in and index array)
! this needs to be fixed to use 1: based arrays
!---------------------------------------------------------------------------------------------------
subroutine heapsort_index( list, idx )

  implicit none

  integer, dimension(0:), intent(in) :: list
  integer, dimension(0:), intent(inout) :: idx

  integer :: t, t_idx
  integer :: n, i, parent, child

  n = size(idx)

  if ( n<2 ) then
    idx(0) = 1
    return
  endif

  i = n/2

  do i = 0, n-1
    idx(i) = i
  enddo

  do
    if ( i > 0 ) then
      i = i-1
      t_idx = idx(i)
      t = list(idx(i))
    else
      n = n-1
      if ( n == 0 ) exit
      t_idx = idx(n)
      t = list(idx(n))
      idx(n) = idx(0)
    endif

    parent = i
    child = i*2+1

    do
      if (child >= n) exit
      if (child +1 < n ) then
         if (list(idx(child)) < list(idx(child+1))) child = child+1
      endif

      if ( t < list(idx(child)) ) then
        idx(parent) = idx(child)
        parent = child
        child = parent*2 + 1
      else
        exit
      endif
    enddo

    idx(parent) = t_idx
  enddo

  do i = 0, size(idx) -1
    idx(i) = idx(i) + 1
  enddo


end subroutine heapsort_index
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine lorentz_transform_single( u, g, a, a0, b )
!---------------------------------------------------
! Lorenz transforms the 4-vector {a0, a} into {b0, b}
! The origin of the b reference frame moves in the a
! reference frame with proper velocity u. g is the
! gamma factor corresponding to u (gamma * beta = u).
!---------------------------------------------------

  implicit none

! dummy variables

  real(p_single), dimension(:), intent(in)   :: u
  real(p_single), intent(in)                 :: g
  real(p_single), dimension(:), intent(in)   :: a
  real(p_single), intent(in)                 :: a0
  real(p_single), dimension(:), intent(out)  :: b

  real(p_single) :: a_u ! dot_product(a,u)
  real(p_single) :: eta ! dot_product(a,u) / (1 + g) - a0

  ! executable statements

  !  b = a + (1.0_p_single / (g + 1.0_p_single)) * (dot_product(a,u) * u) - u * a0
  a_u = a(1)*u(1) + a(2)*u(2) + a(3)*u(3)
  eta = a_u / ( g + 1.0_p_single ) - a0

  b(1) = a(1) + eta*u(1)
  b(2) = a(2) + eta*u(2)
  b(3) = a(3) + eta*u(3)

end subroutine lorentz_transform_single
!---------------------------------------------------

!---------------------------------------------------
subroutine lorentz_transform_double( u, g, a, a0, b )
!---------------------------------------------------
! Lorenz transforms the 4-vector {a0, a} into {b0, b}
! The origin of the b reference frame moves in the a
! reference frame with proper velocity u. g is the
! gamma factor corresponding to u (gamma * beta = u).
!---------------------------------------------------

  implicit none

! dummy variables

  real(p_double), dimension(:), intent(in)   :: u
  real(p_double), intent(in)                 :: g
  real(p_double), dimension(:), intent(in)   :: a
  real(p_double), intent(in)                 :: a0
  real(p_double), dimension(:), intent(out)  :: b

  real(p_double) :: a_u ! dot_product(a,u)
  real(p_double) :: eta ! dot_product(a,u) / (1 + g) - a0

  !  b = a + (1.0_p_double / (g + 1.0_p_double)) * (dot_product(a,u) * u) - u * a0
  a_u = a(1)*u(1) + a(2)*u(2) + a(3)*u(3)
  eta = a_u / ( g + 1.0_p_double ) - a0

  b(1) = a(1) + eta*u(1)
  b(2) = a(2) + eta*u(2)
  b(3) = a(3) + eta*u(3)

end subroutine lorentz_transform_double
!---------------------------------------------------

!---------------------------------------------------
subroutine lorentz_transform_0_single( u, g, a, a0, b, b0 )
!---------------------------------------------------
! Lorenz transforms the 4-vector {a0, a} into {b0, b}
! The origin of the b reference frame moves in the a
! reference frame with proper velocity u. g is the
! gamma factor corresponding to u (gamma * beta = u).
!---------------------------------------------------

  implicit none

! dummy variables

  real(p_single), dimension(:), intent(in)   :: u
  real(p_single), intent(in)                 :: g
  real(p_single), dimension(:), intent(in)   :: a
  real(p_single), intent(in)                 :: a0
  real(p_single), dimension(:), intent(out)  :: b
  real(p_single), intent(out)                :: b0

  real(p_single) :: a_u ! dot_product(a,u)
  real(p_single) :: eta ! dot_product(a,u) / (1 + g) - a0

  ! executable statements

  !  b = a + (1.0_p_single / (g + 1.0_p_single)) * (dot_product(a,u) * u) - u * a0
  a_u = a(1)*u(1) + a(2)*u(2) + a(3)*u(3)
  eta = a_u / ( g + 1.0_p_single ) - a0

  b(1) = a(1) + eta*u(1)
  b(2) = a(2) + eta*u(2)
  b(3) = a(3) + eta*u(3)

  ! b0 = g * a0 - dot_product(a,u)
  b0 = g * a0 - a_u

end subroutine lorentz_transform_0_single
!---------------------------------------------------


!---------------------------------------------------
subroutine lorentz_transform_0_double( u, g, a, a0, b, b0 )
!---------------------------------------------------
! Lorenz transforms the 4-vector {a0, a} into {b0, b}
! The origin of the b reference frame moves in the a
! reference frame with proper velocity u. g is the
! gamma factor corresponding to u (gamma * beta = u).
!---------------------------------------------------

  implicit none

! dummy variables

  real(p_double), dimension(:), intent(in)   :: u
  real(p_double), intent(in)                 :: g
  real(p_double), dimension(:), intent(in)   :: a
  real(p_double), intent(in)                 :: a0
  real(p_double), dimension(:), intent(out)  :: b
  real(p_double), intent(out)                :: b0

  real(p_double) :: a_u ! dot_product(a,u)
  real(p_double) :: eta ! dot_product(a,u) / (1 + g) - a0

  !  b = a + (1.0_p_double / (g + 1.0_p_double)) * (dot_product(a,u) * u) - u * a0
  a_u = a(1)*u(1) + a(2)*u(2) + a(3)*u(3)
  eta = a_u / ( g + 1.0_p_double ) - a0

  b(1) = a(1) + eta*u(1)
  b(2) = a(2) + eta*u(2)
  b(3) = a(3) + eta*u(3)

  ! b0 = g * a0 - dot_product(a,u)
  b0 = g * a0 - a_u

end subroutine lorentz_transform_0_double
!---------------------------------------------------


end module m_utilities


