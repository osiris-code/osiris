!-----------------------------------------------------------------------------------------
! (P)erfect (E)lectric (C)onductor field boundaries
!-----------------------------------------------------------------------------------------
!
! These routines implement PEC field boundaries for the FDTD solver for 2 cases:
!   pec_bc_?d_?_odd  - Boundary at the cell edge (lower edge [1], right edge [nx])
!   pec_bc_?d_?_even - Boundary at the cell center (center [0], center [nx])
!-----------------------------------------------------------------------------------------

#include "os-preprocess.fpp"
#include "os-config.h"

module m_emf_pec

#include "memory/memory.h"

use m_parameters
use m_vdf_define

implicit none

private

interface pec_bc_e_1d_odd
  module procedure pec_bc_e_1d_odd
end interface

interface pec_bc_e_2d_odd
  module procedure pec_bc_e_2d_odd
end interface

interface pec_bc_e_3d_odd
  module procedure pec_bc_e_3d_odd
end interface

interface pec_bc_e_1d_even
  module procedure pec_bc_e_1d_even
end interface

interface pec_bc_e_2d_even
  module procedure pec_bc_e_2d_even
end interface

interface pec_bc_e_3d_even
  module procedure pec_bc_e_3d_even
end interface

interface pec_bc_b_1d_odd
  module procedure pec_bc_b_1d_odd
end interface

interface pec_bc_b_2d_odd
  module procedure pec_bc_b_2d_odd
end interface

interface pec_bc_b_3d_odd
  module procedure pec_bc_b_3d_odd
end interface

interface pec_bc_b_1d_even
  module procedure pec_bc_b_1d_even
end interface

interface pec_bc_b_2d_even
  module procedure pec_bc_b_2d_even
end interface

interface pec_bc_b_3d_even
  module procedure pec_bc_b_3d_even
end interface

public :: pec_bc_e_1d_odd, pec_bc_e_2d_odd, pec_bc_e_3d_odd
public :: pec_bc_b_1d_odd, pec_bc_b_2d_odd, pec_bc_b_3d_odd
public :: pec_bc_e_1d_even, pec_bc_e_2d_even, pec_bc_e_3d_even
public :: pec_bc_b_1d_even, pec_bc_b_2d_even, pec_bc_b_3d_even

contains

!-----------------------------------------------------------------------------------------
! pec_bc_?d_?_odd  - Boundary at the cell edge (lower edge [1], right edge [nx])
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pec_bc_e_1d_odd( e, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: e
  integer, intent(in) :: bnd_idx

  integer :: i1

  select case (bnd_idx)
     case ( p_lower )

       do i1 = lbound( e%f1, 2 ), 0
         e%f1( 1, i1 ) =  e%f1( 1, 1-i1 )

         e%f1( 2, i1 ) = -e%f1( 2, 2-i1 )

         e%f1( 3, i1 ) = -e%f1( 3, 2-i1 )

       enddo

       ! The tangential field components at the surface are zero
       e%f1(2,1) = 0
       e%f1(3,1) = 0

     case ( p_upper )

       do i1 = 2, ubound( e%f1, 2 ) - e%nx_(1)
         e%f1( 1, e%nx_(1) + i1) =  e%f1( 1, e%nx_(1) + 1 - i1 )

         e%f1( 2, e%nx_(1) + i1) = -e%f1( 2, e%nx_(1) + 2 - i1 )

         e%f1( 3, e%nx_(1) + i1) = -e%f1( 3, e%nx_(1) + 2 - i1 )

       enddo

       e%f1(1,e%nx_(1)+1) = e%f1(1,e%nx_(1))
       e%f1(2,e%nx_(1)+1) = 0
       e%f1(3,e%nx_(1)+1) = 0

  end select

end subroutine pec_bc_e_1d_odd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pec_bc_e_2d_odd( e, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: e
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary
             do i2 = lbound( e%f2, 3 ), ubound( e%f2, 3 )
               do i1 = lbound( e%f2, 2 ), 0
                  e%f2( 1, i1, i2 ) =  e%f2( 1, 1-i1, i2 )

                  e%f2( 2, i1, i2 ) = -e%f2( 2, 2-i1, i2 )

                  e%f2( 3, i1, i2 ) = -e%f2( 3, 2-i1, i2 )

               enddo
               e%f2(2,1,i2) = 0
               e%f2(3,1,i2) = 0
             enddo

          case(2) ! x2 boundary

             do i2 = lbound( e%f2, 3 ), 0
               do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
                  e%f2( 1, i1, i2 ) = -e%f2( 1, i1, 2-i2 )

                  e%f2( 2, i1, i2 ) =  e%f2( 2, i1, 1-i2 )

                  e%f2( 3, i1, i2 ) = -e%f2( 3, i1, 2-i2 )

               enddo
             enddo

             do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
               e%f2(1,i1,1) = 0
               e%f2(3,i1,1) = 0
             enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
             do i2 = lbound( e%f2, 3 ), ubound( e%f2, 3 )

               e%f2(1,e%nx_(1)+1, i2) = e%f2(1,e%nx_(1), i2)
               e%f2(2,e%nx_(1)+1, i2) = 0
               e%f2(3,e%nx_(1)+1, i2) = 0

               do i1 = 2, ubound( e%f2, 2 ) - e%nx_(1)
                 e%f2( 1, e%nx_(1) + i1, i2) =  e%f2( 1, e%nx_(1) + 1 - i1, i2)

                 e%f2( 2, e%nx_(1) + i1, i2) = -e%f2( 2, e%nx_(1) + 2 - i1, i2)

                 e%f2( 3, e%nx_(1) + i1, i2) = -e%f2( 3, e%nx_(1) + 2 - i1, i2)

               enddo
             enddo

          case(2) ! x2 boundary

             do i2 = 2, ubound( e%f2, 3 ) - e%nx_(2)
               do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
                 e%f2( 1, i1, e%nx_(2) + i2) = -e%f2( 1, i1, e%nx_(2) + 2 - i2)

                 e%f2( 2, i1, e%nx_(2) + i2) =  e%f2( 2, i1, e%nx_(2) + 1 - i2)

                 e%f2( 3, i1, e%nx_(2) + i2) = -e%f2( 3, i1, e%nx_(2) + 2 - i2)

               enddo
             enddo

             do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
               e%f2(1,i1,e%nx_(2)+1) = 0
               e%f2(2,i1,e%nx_(2)+1) = e%f2(2,i1,e%nx_(2))
               e%f2(3,i1,e%nx_(2)+1) = 0
             enddo

       end select
  end select

end subroutine pec_bc_e_2d_odd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pec_bc_e_3d_odd( e, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: e
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2, i3

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary
             do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
               do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
                 do i1 = lbound( e%f3, 2 ), 0
                    e%f3( 1, i1, i2, i3 ) =  e%f3( 1, 1-i1, i2, i3 )

                    e%f3( 2, i1, i2, i3 ) = -e%f3( 2, 2-i1, i2, i3 )

                    e%f3( 3, i1, i2, i3 ) = -e%f3( 3, 2-i1, i2, i3 )

                 enddo
                 e%f3(2, 1, i2, i3 ) = 0
                 e%f3(3, 1, i2, i3 ) = 0
               enddo
             enddo

          case(2) ! x2 boundary

             do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
               do i2 = lbound( e%f3, 3 ), 0
                 do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                    e%f3( 1, i1, i2, i3 ) = -e%f3( 1, i1, 2-i2, i3 )

                    e%f3( 2, i1, i2, i3 ) =  e%f3( 2, i1, 1-i2, i3 )

                    e%f3( 3, i1, i2, i3 ) = -e%f3( 3, i1, 2-i2, i3 )

                 enddo
               enddo

               do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                 e%f3(1, i1, 1, i3) = 0
                 e%f3(3, i1, 1, i3) = 0
               enddo
             enddo

          case(3) ! x3 boundary

             do i3 = lbound( e%f3, 4 ), 0

               do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
                 do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                    e%f3( 1, i1, i2, i3 ) = -e%f3( 1, i1, i2, 2-i3 )

                    e%f3( 2, i1, i2, i3 ) = -e%f3( 2, i1, i2, 2-i3 )

                    e%f3( 3, i1, i2, i3 ) =  e%f3( 3, i1, i2, 1-i3 )

                 enddo
               enddo
             enddo

             do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
               do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                 e%f3( 1, i1, i2, 1 ) = 0
                 e%f3( 2, i1, i2, 1 ) = 0
               enddo
             enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
             do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
               do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )

                 e%f3(1,e%nx_(1)+1, i2, i3 ) = e%f3(1,e%nx_(1), i2, i3 )
                 e%f3(2,e%nx_(1)+1, i2, i3 ) = 0
                 e%f3(3,e%nx_(1)+1, i2, i3 ) = 0

                 do i1 = 2, ubound( e%f3, 2 ) - e%nx_(1)
                   e%f3( 1, e%nx_(1) + i1, i2, i3 ) =  e%f3( 1, e%nx_(1) + 1 - i1, i2, i3 )

                   e%f3( 2, e%nx_(1) + i1, i2, i3 ) = -e%f3( 2, e%nx_(1) + 2 - i1, i2, i3 )

                   e%f3( 3, e%nx_(1) + i1, i2, i3 ) = -e%f3( 3, e%nx_(1) + 2 - i1, i2, i3 )

                 enddo
               enddo
             enddo

          case(2) ! x2 boundary

             do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
               do i2 = 2, ubound( e%f3, 3 ) - e%nx_(2)
                 do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                   e%f3( 1, i1, e%nx_(2) + i2, i3 ) = -e%f3( 1, i1, e%nx_(2) + 2 - i2, i3 )

                   e%f3( 2, i1, e%nx_(2) + i2, i3 ) =  e%f3( 2, i1, e%nx_(2) + 1 - i2, i3 )

                   e%f3( 3, i1, e%nx_(2) + i2, i3 ) = -e%f3( 3, i1, e%nx_(2) + 2 - i2, i3 )

                 enddo
               enddo

               do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                 e%f3(1,i1,e%nx_(2)+1,i3) = 0
                 e%f3(2,i1,e%nx_(2)+1,i3) = e%f3(2,i1,e%nx_(2),i3)
                 e%f3(3,i1,e%nx_(2)+1,i3) = 0
               enddo
             enddo

          case(3) ! x3 boundary

             do i3 = 2, ubound( e%f3, 4 ) - e%nx_(3)
               do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
                 do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                   e%f3( 1, i1, i2, e%nx_(3) + i3 ) = -e%f3( 1, i1, i2, e%nx_(3) + 2 - i3 )

                   e%f3( 2, i1, i2, e%nx_(3) + i3 ) = -e%f3( 2, i1, i2, e%nx_(3) + 2 - i3 )

                   e%f3( 3, i1, i2, e%nx_(3) + i3 ) =  e%f3( 3, i1, i2, e%nx_(3) + 1 - i3 )

                 enddo
               enddo
             enddo

             do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
               do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                 e%f3(1,i1,i2,e%nx_(3)+1) = 0
                 e%f3(2,i1,i2,e%nx_(3)+1) = 0
                 e%f3(3,i1,i2,e%nx_(3)+1) = e%f3(3,i1,i2,e%nx_(3))
               enddo
             enddo

       end select
  end select

end subroutine pec_bc_e_3d_odd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pec_bc_b_1d_odd( b, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: b
  integer, intent(in) :: bnd_idx

  integer :: i1

  select case (bnd_idx)
     case ( p_lower )

       do i1 = lbound( b%f1, 2 ), 0
         b%f1( 1, i1 ) = -b%f1( 1, 2-i1 )

         b%f1( 2, i1 ) =  b%f1( 2, 1-i1 )

         b%f1( 3, i1 ) =  b%f1( 3, 1-i1 )

       enddo

       ! The perp field components at the surface are zero
       b%f1(1,1) = 0

     case ( p_upper )

       do i1 = 2, ubound( b%f1, 2 ) - b%nx_(1)
         b%f1( 1, b%nx_(1) + i1) = -b%f1( 1, b%nx_(1) + 2 - i1 )

         b%f1( 2, b%nx_(1) + i1) =  b%f1( 2, b%nx_(1) + 1 - i1 )

         b%f1( 3, b%nx_(1) + i1) =  b%f1( 3, b%nx_(1) + 1 - i1 )

       enddo

       b%f1(1,b%nx_(1)+1) = 0

  end select

end subroutine pec_bc_b_1d_odd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pec_bc_b_2d_odd( b, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: b
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary

             do i2 = lbound( b%f2, 3 ), ubound( b%f2, 3 )
               do i1 = lbound( b%f2, 2 ), 0
                 b%f2( 1, i1, i2 ) = -b%f2( 1, 2-i1, i2 )

                 b%f2( 2, i1, i2 ) =  b%f2( 2, 1-i1, i2 )

                 b%f2( 3, i1, i2 ) =  b%f2( 3, 1-i1, i2 )

               enddo

               ! The perp field components at the surface are zero
               b%f2( 1, 1, i2) = 0
             enddo

          case(2) ! x2 boundary

             do i2 = lbound( b%f2, 3 ), 0
               do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )
                 b%f2( 1, i1, i2 ) =   b%f2( 1, i1, 1-i2 )

                 b%f2( 2, i1, i2 ) =  -b%f2( 2, i1, 2-i2 )

                 b%f2( 3, i1, i2 ) =   b%f2( 3, i1, 1-i2 )

               enddo
             enddo

             ! The perp field components at the surface are zero
             do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )
               b%f2( 2, i1, 1 ) = 0
             enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
             do i2 = lbound( b%f2, 3 ), ubound( b%f2, 3 )
               b%f2(2,b%nx_(1)+1, i2) = 0

               do i1 = 2, ubound( b%f2, 2 ) - b%nx_(1)
                 b%f2( 1, b%nx_(1) + i1, i2) = -b%f2( 1, b%nx_(1) + 2 - i1, i2)

                 b%f2( 2, b%nx_(1) + i1, i2) =  b%f2( 2, b%nx_(1) + 1 - i1, i2)

                 b%f2( 3, b%nx_(1) + i1, i2) =  b%f2( 3, b%nx_(1) + 1 - i1, i2)

               enddo
             enddo

          case(2) ! x2 boundary
             do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )
               b%f2( 2, i1, b%nx_(2) + 1) = 0
             enddo

             do i2 = 2, ubound( b%f2, 3 ) - b%nx_(2)
                 do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )
                 b%f2( 1, i1, b%nx_(2) + i2) =  b%f2( 1, i1, b%nx_(2) + 1 - i2)

                 b%f2( 2, i1, b%nx_(2) + i2) = -b%f2( 2, i1, b%nx_(2) + 2 - i2)

                 b%f2( 3, i1, b%nx_(2) + i2) =  b%f2( 3, i1, b%nx_(2) + 1 - i2)

               enddo
             enddo

       end select
  end select

end subroutine pec_bc_b_2d_odd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pec_bc_b_3d_odd( b, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: b
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2, i3

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary

             do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
               do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
                 do i1 = lbound( b%f3, 2 ), 0
                   b%f3( 1, i1, i2, i3 ) = -b%f3( 1, 2-i1, i2, i3 )

                   b%f3( 2, i1, i2, i3 ) =  b%f3( 2, 1-i1, i2, i3 )

                   b%f3( 3, i1, i2, i3 ) =  b%f3( 3, 1-i1, i2, i3 )

                 enddo

                 ! The perp field components at the surface are zero
                 b%f3( 1, 1, i2, i3 ) = 0
               enddo
             enddo

          case(2) ! x2 boundary

             do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
               do i2 = lbound( b%f3, 3 ), 0
                 do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
                   b%f3( 1, i1, i2, i3 ) =   b%f3( 1, i1, 1-i2, i3 )

                   b%f3( 2, i1, i2, i3 ) =  -b%f3( 2, i1, 2-i2, i3 )

                   b%f3( 3, i1, i2, i3 ) =   b%f3( 3, i1, 1-i2, i3 )

                 enddo
               enddo

               ! The perp field components at the surface are zero
               do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
                 b%f3( 2, i1, 1, i3 ) = 0
               enddo
             enddo

          case(3) ! x3 boundary

             do i3 = lbound( b%f3, 4 ), 0
               do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
                 do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
                   b%f3( 1, i1, i2, i3 ) =   b%f3( 1, i1, i2, 1-i3 )

                   b%f3( 2, i1, i2, i3 ) =   b%f3( 2, i1, i2, 1-i3 )

                   b%f3( 3, i1, i2, i3 ) =  -b%f3( 3, i1, i2, 2-i3 )

                 enddo
               enddo
             enddo

             ! The perp field components at the surface are zero
             do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
               do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
                 b%f3( 3, i1, i2, 1 ) = 0
               enddo
             enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
             do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
               do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
                 b%f3( 1, b%nx_(1)+1, i2, i3 ) = 0

                 do i1 = 2, ubound( b%f3, 2 ) - b%nx_(1)
                   b%f3( 1, b%nx_(1) + i1, i2, i3 ) = -b%f3( 1, b%nx_(1) + 2 - i1, i2, i3 )

                   b%f3( 2, b%nx_(1) + i1, i2, i3 ) =  b%f3( 2, b%nx_(1) + 1 - i1, i2, i3 )

                   b%f3( 3, b%nx_(1) + i1, i2, i3 ) =  b%f3( 3, b%nx_(1) + 1 - i1, i2, i3 )

                 enddo
               enddo
             enddo

          case(2) ! x2 boundary
             do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
               do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
                 b%f3( 2, i1, b%nx_(2) + 1, i3 ) = 0
               enddo

               do i2 = 2, ubound( b%f3, 3 ) - b%nx_(2)
                 do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
                   b%f3( 1, i1, b%nx_(2) + i2, i3 ) =  b%f3( 1, i1, b%nx_(2) + 1 - i2, i3 )

                   b%f3( 2, i1, b%nx_(2) + i2, i3 ) = -b%f3( 2, i1, b%nx_(2) + 2 - i2, i3 )

                   b%f3( 3, i1, b%nx_(2) + i2, i3 ) =  b%f3( 3, i1, b%nx_(2) + 1 - i2, i3 )

                 enddo
               enddo
             enddo

          case(3) ! x3 boundary

             do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
               do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
                 b%f3( 3, i1, i2, b%nx_(3) + 1 ) = 0
               enddo
             enddo

             do i3 = 2, ubound( b%f3, 4 ) - b%nx_(3)
               do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
                 do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
                   b%f3( 1, i1, i2, b%nx_(3) + i3 ) =  b%f3( 1, i1, i2, b%nx_(3) + 1 - i3 )

                   b%f3( 2, i1, i2, b%nx_(3) + i3 ) =  b%f3( 2, i1, i2, b%nx_(3) + 1 - i3 )

                   b%f3( 3, i1, i2, b%nx_(3) + i3 ) = -b%f3( 3, i1, i2, b%nx_(3) + 2 - i3 )

                 enddo
               enddo
             enddo

       end select
  end select

end subroutine pec_bc_b_3d_odd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! pec_bc_?d_?_even  - Boundary at the cell center (center [0], center [nx])
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pec_bc_e_1d_even( e, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: e
  integer, intent(in) :: bnd_idx

  integer :: i1

  select case (bnd_idx)
     case ( p_lower )

       do i1 = lbound( e%f1, 2 ), 0
         e%f1( 1, i1 ) =  e%f1( 1,  -i1 )

         e%f1( 2, i1 ) = -e%f1( 2, 1-i1 )

         e%f1( 3, i1 ) = -e%f1( 3, 1-i1 )

       enddo

     case ( p_upper )

       do i1 = 1, ubound( e%f1, 2 ) - e%nx_(1)
         e%f1( 1, e%nx_(1) + i1) =  e%f1( 1, e%nx_(1)     - i1 )

         e%f1( 2, e%nx_(1) + i1) = -e%f1( 2, e%nx_(1) + 1 - i1 )

         e%f1( 3, e%nx_(1) + i1) = -e%f1( 3, e%nx_(1) + 1 - i1 )

       enddo

       e%f1(2,e%nx_(1)+1) = 0
       e%f1(3,e%nx_(1)+1) = 0

  end select

end subroutine pec_bc_e_1d_even
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pec_bc_e_2d_even( e, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: e
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary
             do i2 = lbound( e%f2, 3 ), ubound( e%f2, 3 )
               do i1 = lbound( e%f2, 2 ), 0
                  e%f2( 1, i1, i2 ) =  e%f2( 1,  -i1, i2 )

                  e%f2( 2, i1, i2 ) = -e%f2( 2, 1-i1, i2 )

                  e%f2( 3, i1, i2 ) = -e%f2( 3, 1-i1, i2 )

               enddo
             enddo

          case(2) ! x2 boundary

             do i2 = lbound( e%f2, 3 ), 0
               do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
                  e%f2( 1, i1, i2 ) = -e%f2( 1, i1, 1-i2 )

                  e%f2( 2, i1, i2 ) =  e%f2( 2, i1,  -i2 )

                  e%f2( 3, i1, i2 ) = -e%f2( 3, i1, 1-i2 )

               enddo
             enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
             do i2 = lbound( e%f2, 3 ), ubound( e%f2, 3 )
               do i1 = 1, ubound( e%f2, 2 ) - e%nx_(1)
                 e%f2( 1, e%nx_(1) + i1, i2) =  e%f2( 1, e%nx_(1)     - i1, i2)

                 e%f2( 2, e%nx_(1) + i1, i2) = -e%f2( 2, e%nx_(1) + 1 - i1, i2)

                 e%f2( 3, e%nx_(1) + i1, i2) = -e%f2( 3, e%nx_(1) + 1 - i1, i2)

               enddo
             enddo

          case(2) ! x2 boundary

             do i2 = 1, ubound( e%f2, 3 ) - e%nx_(2)
               do i1 = lbound( e%f2, 2 ), ubound( e%f2, 2 )
                 e%f2( 1, i1, e%nx_(2) + i2) = -e%f2( 1, i1, e%nx_(2) + 1 - i2)

                 e%f2( 2, i1, e%nx_(2) + i2) =  e%f2( 2, i1, e%nx_(2)     - i2)

                 e%f2( 3, i1, e%nx_(2) + i2) = -e%f2( 3, i1, e%nx_(2) + 1 - i2)

               enddo
             enddo

       end select
  end select

end subroutine pec_bc_e_2d_even
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pec_bc_e_3d_even( e, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: e
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2, i3

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary
             do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
               do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
                 do i1 = lbound( e%f3, 2 ), 0
                    e%f3( 1, i1, i2, i3 ) =  e%f3( 1,  -i1, i2, i3 )

                    e%f3( 2, i1, i2, i3 ) = -e%f3( 2, 1-i1, i2, i3 )

                    e%f3( 3, i1, i2, i3 ) = -e%f3( 3, 1-i1, i2, i3 )

                 enddo
               enddo
             enddo

          case(2) ! x2 boundary

             do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
               do i2 = lbound( e%f3, 3 ), 0
                 do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                    e%f3( 1, i1, i2, i3 ) = -e%f3( 1, i1, 1-i2, i3 )

                    e%f3( 2, i1, i2, i3 ) =  e%f3( 2, i1,  -i2, i3 )

                    e%f3( 3, i1, i2, i3 ) = -e%f3( 3, i1, 1-i2, i3 )

                 enddo
               enddo
             enddo

          case(3) ! x3 boundary

             do i3 = lbound( e%f3, 4 ), 0

               do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
                 do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                    e%f3( 1, i1, i2, i3 ) = -e%f3( 1, i1, i2, 1-i3 )

                    e%f3( 2, i1, i2, i3 ) = -e%f3( 2, i1, i2, 1-i3 )

                    e%f3( 3, i1, i2, i3 ) =  e%f3( 3, i1, i2,  -i3 )

                 enddo
               enddo
             enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
             do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
               do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
                 do i1 = 1, ubound( e%f3, 2 ) - e%nx_(1)
                   e%f3( 1, e%nx_(1) + i1, i2, i3 ) =  e%f3( 1, e%nx_(1)     - i1, i2, i3 )

                   e%f3( 2, e%nx_(1) + i1, i2, i3 ) = -e%f3( 2, e%nx_(1) + 1 - i1, i2, i3 )

                   e%f3( 3, e%nx_(1) + i1, i2, i3 ) = -e%f3( 3, e%nx_(1) + 1 - i1, i2, i3 )

                 enddo
               enddo
             enddo

          case(2) ! x2 boundary

             do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
               do i2 = 1, ubound( e%f3, 3 ) - e%nx_(2)
                 do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                   e%f3( 1, i1, e%nx_(2) + i2, i3 ) = -e%f3( 1, i1, e%nx_(2) + 1 - i2, i3 )

                   e%f3( 2, i1, e%nx_(2) + i2, i3 ) =  e%f3( 2, i1, e%nx_(2)     - i2, i3 )

                   e%f3( 3, i1, e%nx_(2) + i2, i3 ) = -e%f3( 3, i1, e%nx_(2) + 1 - i2, i3 )

                 enddo
               enddo
             enddo

          case(3) ! x3 boundary

             do i3 = 1, ubound( e%f3, 4 ) - e%nx_(3)
               do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
                 do i1 = lbound( e%f3, 2 ), ubound( e%f3, 2 )
                   e%f3( 1, i1, i2, e%nx_(3) + i3 ) = -e%f3( 1, i1, i2, e%nx_(3) + 1 - i3 )

                   e%f3( 2, i1, i2, e%nx_(3) + i3 ) = -e%f3( 2, i1, i2, e%nx_(3) + 1 - i3 )

                   e%f3( 3, i1, i2, e%nx_(3) + i3 ) =  e%f3( 3, i1, i2, e%nx_(3)     - i3 )

                 enddo
               enddo
             enddo

       end select
  end select

end subroutine pec_bc_e_3d_even
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pec_bc_b_1d_even( b, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: b
  integer, intent(in) :: bnd_idx

  integer :: i1

  select case (bnd_idx)
     case ( p_lower )

       do i1 = lbound( b%f1, 2 ), 0
         b%f1( 1, i1 ) = -b%f1( 1, 1-i1 )

         b%f1( 2, i1 ) =  b%f1( 2,  -i1 )

         b%f1( 3, i1 ) =  b%f1( 3,  -i1 )

       enddo

     case ( p_upper )

       do i1 = 1, ubound( b%f1, 2 ) - b%nx_(1)
         b%f1( 1, b%nx_(1) + i1) = -b%f1( 1, b%nx_(1) + 1 - i1 )

         b%f1( 2, b%nx_(1) + i1) =  b%f1( 2, b%nx_(1)     - i1 )

         b%f1( 3, b%nx_(1) + i1) =  b%f1( 3, b%nx_(1)     - i1 )

       enddo

  end select

end subroutine pec_bc_b_1d_even
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pec_bc_b_2d_even( b, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: b
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary

             do i2 = lbound( b%f2, 3 ), ubound( b%f2, 3 )
               do i1 = lbound( b%f2, 2 ), 0
                 b%f2( 1, i1, i2 ) = -b%f2( 1, 1-i1, i2 )

                 b%f2( 2, i1, i2 ) =  b%f2( 2,  -i1, i2 )

                 b%f2( 3, i1, i2 ) =  b%f2( 3,  -i1, i2 )

               enddo
             enddo

          case(2) ! x2 boundary

             do i2 = lbound( b%f2, 3 ), 0
               do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )
                 b%f2( 1, i1, i2 ) =   b%f2( 1, i1,  -i2 )

                 b%f2( 2, i1, i2 ) =  -b%f2( 2, i1, 1-i2 )

                 b%f2( 3, i1, i2 ) =   b%f2( 3, i1,  -i2 )

               enddo
             enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
             do i2 = lbound( b%f2, 3 ), ubound( b%f2, 3 )
               do i1 = 1, ubound( b%f2, 2 ) - b%nx_(1)
                 b%f2( 1, b%nx_(1) + i1, i2) = -b%f2( 1, b%nx_(1) + 1 - i1, i2)

                 b%f2( 2, b%nx_(1) + i1, i2) =  b%f2( 2, b%nx_(1)     - i1, i2)

                 b%f2( 3, b%nx_(1) + i1, i2) =  b%f2( 3, b%nx_(1)     - i1, i2)

               enddo
             enddo

          case(2) ! x2 boundary

             do i2 = 1, ubound( b%f2, 3 ) - b%nx_(2)
                 do i1 = lbound( b%f2, 2 ), ubound( b%f2, 2 )
                 b%f2( 1, i1, b%nx_(2) + i2) =  b%f2( 1, i1, b%nx_(2)     - i2)

                 b%f2( 2, i1, b%nx_(2) + i2) = -b%f2( 2, i1, b%nx_(2) + 1 - i2)

                 b%f2( 3, i1, b%nx_(2) + i2) =  b%f2( 3, i1, b%nx_(2)     - i2)

               enddo
             enddo

       end select
  end select

end subroutine pec_bc_b_2d_even
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine pec_bc_b_3d_even( b, i_dim, bnd_idx )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_vdf ) , intent(inout) :: b
  integer, intent(in) :: i_dim, bnd_idx

  integer :: i1, i2, i3

  select case (bnd_idx)

     case ( p_lower )

       select case ( i_dim )
          case(1) ! x1 boundary

             do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
               do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
                 do i1 = lbound( b%f3, 2 ), 0
                   b%f3( 1, i1, i2, i3 ) = -b%f3( 1, 1-i1, i2, i3 )

                   b%f3( 2, i1, i2, i3 ) =  b%f3( 2,  -i1, i2, i3 )

                   b%f3( 3, i1, i2, i3 ) =  b%f3( 3,  -i1, i2, i3 )

                 enddo
               enddo
             enddo

          case(2) ! x2 boundary

             do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
               do i2 = lbound( b%f3, 3 ), 0
                 do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
                   b%f3( 1, i1, i2, i3 ) =   b%f3( 1, i1,  -i2, i3 )

                   b%f3( 2, i1, i2, i3 ) =  -b%f3( 2, i1, 1-i2, i3 )

                   b%f3( 3, i1, i2, i3 ) =   b%f3( 3, i1,  -i2, i3 )

                 enddo
               enddo
             enddo

          case(3) ! x3 boundary

             do i3 = lbound( b%f3, 4 ), 0
               do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
                 do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
                   b%f3( 1, i1, i2, i3 ) =   b%f3( 1, i1, i2,  -i3 )

                   b%f3( 2, i1, i2, i3 ) =   b%f3( 2, i1, i2,  -i3 )

                   b%f3( 3, i1, i2, i3 ) =  -b%f3( 3, i1, i2, 1-i3 )

                 enddo
               enddo
             enddo

       end select

     case ( p_upper )

       select case ( i_dim )
          case(1) ! x1 boundary
             do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
               do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
                 do i1 = 1, ubound( b%f3, 2 ) - b%nx_(1)
                   b%f3( 1, b%nx_(1) + i1, i2, i3 ) = -b%f3( 1, b%nx_(1) + 1 - i1, i2, i3 )

                   b%f3( 2, b%nx_(1) + i1, i2, i3 ) =  b%f3( 2, b%nx_(1)     - i1, i2, i3 )

                   b%f3( 3, b%nx_(1) + i1, i2, i3 ) =  b%f3( 3, b%nx_(1)     - i1, i2, i3 )

                 enddo
               enddo
             enddo

          case(2) ! x2 boundary

             do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
               do i2 = 1, ubound( b%f3, 3 ) - b%nx_(2)
                 do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
                   b%f3( 1, i1, b%nx_(2) + i2, i3 ) =  b%f3( 1, i1, b%nx_(2)     - i2, i3 )

                   b%f3( 2, i1, b%nx_(2) + i2, i3 ) = -b%f3( 2, i1, b%nx_(2) + 1 - i2, i3 )

                   b%f3( 3, i1, b%nx_(2) + i2, i3 ) =  b%f3( 3, i1, b%nx_(2)     - i2, i3 )

                 enddo
               enddo
             enddo

          case(3) ! x3 boundary

             do i3 = 1, ubound( b%f3, 4 ) - b%nx_(3)
               do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
                 do i1 = lbound( b%f3, 2 ), ubound( b%f3, 2 )
                   b%f3( 1, i1, i2, b%nx_(3) + i3 ) =  b%f3( 1, i1, i2, b%nx_(3)     - i3 )

                   b%f3( 2, i1, i2, b%nx_(3) + i3 ) =  b%f3( 2, i1, i2, b%nx_(3)     - i3 )

                   b%f3( 3, i1, i2, b%nx_(3) + i3 ) = -b%f3( 3, i1, i2, b%nx_(3) + 1 - i3 )

                 enddo
               enddo
             enddo

       end select
  end select

end subroutine pec_bc_b_3d_even
!-----------------------------------------------------------------------------------------

end module m_emf_pec
