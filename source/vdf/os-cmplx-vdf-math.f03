#include "os-preprocess.fpp"
#include "os-config.h"

module m_cmplx_vdf_math

#include "memory/memory.h"

use m_parameters

use m_cmplx_vdf
use m_vdf_define

! This is isn't necessary for compilation but it helps identify
! incompatibilities between the m_cmplx_vdf_math and the m_vdf_math
! modules
use m_vdf_math

private

interface add
  module procedure add_cmplx_cmplx_vdf
  module procedure add_real_cmplx_vdf
  module procedure add_mult_cmplx_vdf
  module procedure add_scalar_cmplx_vdf
  module procedure add_fc_cmplx_vdf
end interface

interface add_pow
  module procedure add_pow_int_cmplx_vdf
end interface

interface sub
  module procedure sub_cmplx_cmplx_vdf
  module procedure sub_real_cmplx_vdf
end interface

interface mult
  module procedure mult_cmplx_vdf
  module procedure mult_scalar_cmplx_vdf
end interface

interface sqrt_vdf
  module procedure sqrt_cmplx_vdf
end interface

public :: add, add_pow, sub, mult, sqrt_vdf

contains

!-------------------------------------------------------------------------------
! a = a + real b
!-------------------------------------------------------------------------------
subroutine add_real_cmplx_vdf( vdf_a, vdf_b )

  implicit none

  type(t_cmplx_vdf), intent(inout) :: vdf_a
  type(t_vdf), intent(in) :: vdf_b

  integer :: i

  !$omp parallel do
  do i = 1, size( vdf_a % zbuffer )
    vdf_a % zbuffer(i) = vdf_a % zbuffer(i) + vdf_b % buffer(i)
  enddo

end subroutine add_real_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! a = a + real b
!-------------------------------------------------------------------------------
subroutine add_cmplx_cmplx_vdf( vdf_a, vdf_b )

  implicit none

  type(t_cmplx_vdf), intent(inout) :: vdf_a
  type(t_cmplx_vdf), intent(in) :: vdf_b

  integer :: i

  !$omp parallel do
  do i = 1, size( vdf_a % zbuffer )
    vdf_a % zbuffer(i) = vdf_a % zbuffer(i) + vdf_b % zbuffer(i)
  enddo

end subroutine add_cmplx_cmplx_vdf
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! a = a - real b
!-------------------------------------------------------------------------------
subroutine sub_real_cmplx_vdf( vdf_a, vdf_b )

  implicit none

  type(t_cmplx_vdf), intent(inout) :: vdf_a
  type(t_vdf), intent(in) :: vdf_b

  integer :: i

  !$omp parallel do
  do i = 1, size( vdf_a % zbuffer )
    vdf_a % zbuffer(i) = vdf_a % zbuffer(i) - vdf_b % buffer(i)
  enddo

end subroutine sub_real_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! a = a - real b
!-------------------------------------------------------------------------------
subroutine sub_cmplx_cmplx_vdf( vdf_a, vdf_b )

  implicit none

  type(t_cmplx_vdf), intent(inout) :: vdf_a
  type(t_cmplx_vdf), intent(in) :: vdf_b

  integer :: i

  !$omp parallel do
  do i = 1, size( vdf_a % zbuffer )
    vdf_a % zbuffer(i) = vdf_a % zbuffer(i) - vdf_b % zbuffer(i)
  enddo

end subroutine sub_cmplx_cmplx_vdf
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!  vdf_a = a1*vdf_a + a2*vdf_b
!-------------------------------------------------------------------------------
subroutine add_mult_cmplx_vdf( vdf_a, vdf_b, a1, a2 )

  implicit none

  type(t_cmplx_vdf), intent(inout) :: vdf_a
  type(t_cmplx_vdf), intent(in) :: vdf_b

  complex(p_k_fld), intent(in) :: a1
  complex(p_k_fld), intent(in) :: a2

  integer :: i

  !$omp parallel do
  do i = 1, size( vdf_a % zbuffer )
    vdf_a % zbuffer(i) = a1 * vdf_a % zbuffer(i) + a2 * vdf_b % zbuffer(i)
  enddo

end subroutine add_mult_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! does a = a + val
!-------------------------------------------------------------------------------
subroutine add_scalar_cmplx_vdf( vdf_a, val )

  implicit none

  type(t_cmplx_vdf), intent(inout) :: vdf_a
  complex(p_k_fld), intent(in) :: val

  integer :: i

  !$omp parallel do
  do i = 1, size( vdf_a % buffer )
    vdf_a % zbuffer(i) = vdf_a % zbuffer(i) + val
  enddo

end subroutine add_scalar_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vdf_a(fc_a) = a1*vdf_a(fc_a) + a2*vdf_b(fc_b)
!-------------------------------------------------------------------------------
subroutine add_fc_cmplx_vdf( vdf_a, fc_a, vdf_b, fc_b, a1, a2 )

  implicit none

  type(t_cmplx_vdf), intent(inout) :: vdf_a
  integer,           intent(in)    :: fc_a
  type(t_cmplx_vdf), intent(in)    :: vdf_b
  integer,           intent(in)    :: fc_b

  real(p_k_fld), optional, intent(in) :: a1
  real(p_k_fld), optional, intent(in) :: a2

  real(p_k_fld) :: la1, la2

  if (present(a1)) then
    la1 = a1
  else
    la1 = 1.0_p_k_fld
  endif

  if (present(a2)) then
    la2 = a2
  else
    la2 = 1.0_p_k_fld
  endif

  select case (vdf_a%x_dim_)
    case(1)
      vdf_a%z1(fc_a,:) = la1*vdf_a%z1(fc_a,:) + la2*vdf_b%z1(fc_b,:)
    case(2)
      vdf_a%z2(fc_a,:,:) = la1*vdf_a%z2(fc_a,:,:) + la2*vdf_b%z2(fc_b,:,:)
    case(3)
      vdf_a%z3(fc_a,:,:,:) = la1*vdf_a%z3(fc_a,:,:,:) + la2*vdf_b%z3(fc_b,:,:,:)
  end select

end subroutine add_fc_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vdf_a(1) = vdf_a(1) + vdf_b(fc)^int_exp
!-------------------------------------------------------------------------------
subroutine add_pow_int_cmplx_vdf( vdf_a, vdf_b, int_exp, fc )

  implicit none

  type(t_cmplx_vdf), intent(inout) :: vdf_a
  type(t_cmplx_vdf),    intent(in) :: vdf_b
  integer,        intent(in) :: int_exp
  integer,        intent(in) :: fc

  integer :: j, k, l


  select case (vdf_a%x_dim_)
  case (1)
    !$omp parallel do
    do j = lbound( vdf_a%z1, 2 ), ubound( vdf_a%z1, 2 )
      vdf_a%z1( 1, j ) = vdf_a%z1( 1, j ) + &
                         vdf_b%z1( fc, j)**int_exp
    enddo

  case (2)
    !$omp parallel do private(j)
    do k = lbound( vdf_a%z2, 3 ), ubound( vdf_a%z2, 3 )
      do j = lbound( vdf_a%z2, 2 ), ubound( vdf_a%z2, 2 )
        vdf_a%z2( 1, j, k ) = vdf_a%z2( 1, j, k ) + &
                              vdf_b%z2( fc,j, k )**int_exp
      enddo
    enddo

  case (3)
    !$omp parallel do private(j,k)
    do l = lbound( vdf_a%z3, 4 ), ubound( vdf_a%z3, 4 )
      do k = lbound( vdf_a%z3, 3 ), ubound( vdf_a%z3, 3 )
        do j = lbound( vdf_a%z3, 2 ), ubound( vdf_a%z3, 2 )
          vdf_a%z3( 1, j, k, l ) = vdf_a%z3( 1, j, k, l ) + &
                                   vdf_b%z3( fc,j, k, l )**int_exp
        enddo
      enddo
    enddo

  end select

end subroutine add_pow_int_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! a = a * b
!-------------------------------------------------------------------------------
subroutine mult_cmplx_vdf( vdf_a, vdf_b )

  implicit none

  type(t_cmplx_vdf), intent(inout) :: vdf_a
  type(t_cmplx_vdf), intent(in) :: vdf_b

  integer :: i

  !$omp parallel do
  do i = 1, size( vdf_a % zbuffer )
    vdf_a % zbuffer(i) = vdf_a % zbuffer(i) * vdf_b % zbuffer(i)
  enddo

end subroutine mult_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! a = a * val
!-------------------------------------------------------------------------------
subroutine mult_scalar_cmplx_vdf( vdf_a, val )

  implicit none

  type(t_cmplx_vdf), intent(inout) :: vdf_a
  complex(p_k_fld), intent(in) :: val

  integer :: i

  !$omp parallel do
  do i = 1, size( vdf_a % zbuffer )
    vdf_a % zbuffer(i) = vdf_a % zbuffer(i) * val
  enddo

end subroutine mult_scalar_cmplx_vdf
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! a = sqrt(a)
!-------------------------------------------------------------------------------
subroutine sqrt_cmplx_vdf( vdf_a )

  implicit none

  type(t_cmplx_vdf), intent(inout) :: vdf_a

  integer :: i

  !$omp parallel do
  do i = 1, size( vdf_a % zbuffer )
    vdf_a % zbuffer(i) = sqrt( vdf_a % zbuffer(i) )
  enddo

end subroutine sqrt_cmplx_vdf
!-------------------------------------------------------------------------------


end module m_cmplx_vdf_math
