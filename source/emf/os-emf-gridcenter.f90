!-----------------------------------------------------------------------------------------
! Centers field values on the corner of the grid
!
!-----------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_gridcenter

#include "memory/memory.h"

use m_parameters
use m_emf_define

implicit none

private

interface grid_center
  module procedure grid_center
end interface

public :: grid_center

contains

!-----------------------------------------------------------------------------------------
! Converts from a staggered Yee mesh to a grid with field values centered on the corner
! of the cell
!-----------------------------------------------------------------------------------------
subroutine grid_center( this )
!-----------------------------------------------------------------------------------------

  use m_vdf_define
  use m_vdf

  implicit none

  class( t_emf ), intent(inout) :: this
  integer :: i1, i2, i3

  ! temporary vdf buffer for centering calculations
  type( t_vdf ) :: buf

  call buf % new( this%e )

  select case ( p_x_dim )

    case (1)
      do i1 = lbound( buf%f1, 2 ) + 1, ubound( buf%f1, 2 )
        buf%f1( 1, i1 ) = 0.5 * ( this%e_part%f1( 1, i1  ) + &
                                  this%e_part%f1( 1, i1-1) )
        buf%f1( 2, i1 ) =         this%e_part%f1( 2, i1  )
        buf%f1( 3, i1 ) =         this%e_part%f1( 3, i1  )
      enddo

      this%e_part = buf

      do i1 = lbound( buf%f1, 2 ) + 1, ubound( buf%f1, 2 )
        buf%f1( 1, i1 ) =          this%b_part%f1( 1, i1  )
        buf%f1( 2, i1 ) =  0.5 * ( this%b_part%f1( 2, i1  ) + &
                                   this%b_part%f1( 2, i1-1) )
        buf%f1( 3, i1 ) =  0.5 * ( this%b_part%f1( 3, i1  ) + &
                                   this%b_part%f1( 3, i1-1) )
      enddo

      this%b_part = buf

    case (2)

      do i2 = lbound( buf%f2, 3 ) + 1, ubound( buf%f2, 3 )
        do i1 = lbound( buf%f2, 2 ) + 1, ubound( buf%f2, 2 )
          buf%f2( 1, i1, i2 ) = 0.5 * ( this%e_part%f2( 1, i1  , i2  ) + &
                                        this%e_part%f2( 1, i1-1, i2  ) )
          buf%f2( 2, i1, i2 ) = 0.5 * ( this%e_part%f2( 2, i1  , i2  ) + &
                                        this%e_part%f2( 2, i1  , i2-1) )
          buf%f2( 3, i1, i2 ) =         this%e_part%f2( 3, i1  , i2  )
        enddo
      enddo

      this%e_part = buf

      do i2 = lbound( buf%f2, 3 ) + 1, ubound( buf%f2, 3 )
        do i1 = lbound( buf%f2, 2 ) + 1, ubound( buf%f2, 2 )
          buf%f2( 1, i1, i2 ) =  0.5 * ( this%b_part%f2( 1, i1  , i2  ) + &
                                         this%b_part%f2( 1, i1  , i2-1) )

          buf%f2( 2, i1, i2) =  0.5 * ( this%b_part%f2( 2, i1  , i2  ) + &
                                         this%b_part%f2( 2, i1-1, i2  ) )

          buf%f2( 3, i1, i2 ) = 0.25 * ( this%b_part%f2( 3, i1  , i2  ) + &
                                         this%b_part%f2( 3, i1-1, i2  ) + &

                                         this%b_part%f2( 3, i1  , i2-1) + &
                                         this%b_part%f2( 3, i1-1, i2-1) )
        enddo
      enddo

      this%b_part = buf

    case (3)

      do i3 = lbound( buf%f3, 4 ) + 1, ubound( buf%f3, 4 )
        do i2 = lbound( buf%f3, 3 ) + 1, ubound( buf%f3, 3 )
          do i1 = lbound( buf%f3, 2 ) + 1, ubound( buf%f3, 2 )
            buf%f3( 1, i1, i2, i3 ) = 0.5 * ( this%e_part%f3( 1, i1  , i2  ,   i3 ) + &
                                              this%e_part%f3( 1, i1-1, i2  ,   i3 ) )
            buf%f3( 2, i1, i2, i3 ) = 0.5 * ( this%e_part%f3( 2, i1  , i2  ,   i3 ) + &
                                              this%e_part%f3( 2, i1  , i2-1,   i3 ) )
            buf%f3( 3, i1, i2, i3 ) = 0.5 * ( this%e_part%f3( 3, i1  , i2  ,   i3 ) + &
                                              this%e_part%f3( 3, i1  , i2  , i3-1 ) )
          enddo
        enddo
      enddo

      this%e_part = buf

      do i3 = lbound( buf%f3, 4 ) + 1, ubound( buf%f3, 4 )
        do i2 = lbound( buf%f3, 3 ) + 1, ubound( buf%f3, 3 )
          do i1 = lbound( buf%f3, 2 ) + 1, ubound( buf%f3, 2 )
            buf%f3( 1, i1, i2, i3 ) = 0.25 * ( this%b_part%f3( 1, i1  , i2  , i3   ) + &
                                               this%b_part%f3( 1, i1  , i2-1, i3   ) + &

                                               this%b_part%f3( 1, i1  , i2  , i3-1 ) + &
                                               this%b_part%f3( 1, i1  , i2-1, i3-1 ) )

            buf%f3( 2, i1, i2, i3 ) = 0.25 * ( this%b_part%f3( 2, i1  , i2  , i3   ) + &
                                               this%b_part%f3( 2, i1-1, i2  , i3   ) + &

                                               this%b_part%f3( 2, i1  , i2  , i3-1 ) + &
                                               this%b_part%f3( 2, i1-1, i2  , i3-1 ) )

            buf%f3( 3, i1, i2, i3 ) = 0.25 * ( this%b_part%f3( 3, i1  , i2  , i3   ) + &
                                               this%b_part%f3( 3, i1-1, i2  , i3   ) + &

                                               this%b_part%f3( 3, i1  , i2-1, i3   ) + &
                                               this%b_part%f3( 3, i1-1, i2-1, i3   ) )
          enddo
        enddo
      enddo

      this%b_part = buf

  end select

  call buf % cleanup()

end subroutine grid_center
!-----------------------------------------------------------------------------------------

end module m_emf_gridcenter
