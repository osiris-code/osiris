!-----------------------------------------------------------------------------------------
!  Poynting Flux (S) Calculations
!  S is defined as S = E x B and is calculated at the lower corner of the cell (in the
!  same position as the charge)
!-----------------------------------------------------------------------------------------


#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_poynting


  use m_emf_define
  use m_vdf_define

  use m_parameters

  implicit none

  private

  interface poynting
    module procedure poynting_emf
  end interface

  public :: poynting

contains


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine poynting_1d( emf, comp, S )

  implicit none

  class( t_emf ), intent(in) :: emf
  integer, intent(in) :: comp
  type( t_vdf ), intent(inout) :: S

  integer :: i1
  real( p_k_fld ) :: e1, e2, e3, b1, b2, b3

  select case ( comp )

    case (1)

       do i1 = 1, emf%e%nx_(1)

          !e1 = 0.5_p_k_fld * ( emf%e%f1(1, i1) + emf%e%f1(1, i1-1) )
          e2 = emf%e%f1(2, i1)
          e3 = emf%e%f1(3, i1)

          !b1 = emf%b%f1( 1, i1 )

          b2 = 0.5_p_k_fld * ( emf%b%f1( 2, i1 ) + emf%b%f1( 2, i1-1) )

          b3 = 0.5_p_k_fld * ( emf%b%f1( 3, i1 ) + emf%b%f1( 3, i1-1) )

          S%f1(1, i1) = e2*b3 - b2*e3   ! S1

          ! S%f1(1, i1) = e3*b1 - e1*b3   ! S2
          ! S%f1(1, i1) = e1*b2 - b1*e2   ! S3

       enddo

    case (2)

       do i1 = 1, emf%e%nx_(1)

          e1 = 0.5_p_k_fld * ( emf%e%f1(1, i1) + emf%e%f1(1, i1-1) )
          !e2 = emf%e%f1(2, i1)
          e3 = emf%e%f1(3, i1)

          b1 = emf%b%f1( 1, i1 )
          !b2 = 0.5_p_k_fld * ( emf%b%f1( 2, i1 ) + emf%b%f1( 2, i1-1) )
          b3 = 0.5_p_k_fld * ( emf%b%f1( 3, i1 ) + emf%b%f1( 3, i1-1) )

          ! S%f1(1, i1) = e2*b3 - b2*e3   ! S1
          S%f1(1, i1) = e3*b1 - e1*b3   ! S2
          ! S%f1(1, i1) = e1*b2 - b1*e2   ! S3

       enddo

    case (3)

       do i1 = 1, emf%e%nx_(1)

          e1 = 0.5_p_k_fld * ( emf%e%f1(1, i1) + emf%e%f1(1, i1-1) )
          e2 = emf%e%f1(2, i1)
          ! e3 = emf%e%f1(3, i1)

          b1 = emf%b%f1( 1, i1 )
          b2 = 0.5_p_k_fld * ( emf%b%f1( 2, i1 ) + emf%b%f1( 2, i1-1) )
          ! b3 = 0.5_p_k_fld * ( emf%b%f1( 3, i1 ) + emf%b%f1( 3, i1-1) )

          ! S%f1(1, i1) = e2*b3 - b2*e3   ! S1
          ! S%f1(1, i1) = e3*b1 - e1*b3   ! S2
          S%f1(1, i1) = e1*b2 - b1*e2   ! S3

       enddo

  end select


end subroutine poynting_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine poynting_2d( emf, comp, S )

  implicit none

  class( t_emf ), intent(in) :: emf
  integer, intent(in) :: comp
  type( t_vdf ), intent(inout) :: S

  integer :: i1, i2
  real( p_k_fld ) :: e1, e2, e3, b1, b2, b3

  select case ( comp )

    case (1)

       do i2 = 1, emf%e%nx_(2)
         do i1 = 1, emf%e%nx_(1)

            !e1 = 0.5_p_k_fld * ( emf%e%f2(1, i1, i2) + emf%e%f2(1, i1-1, i2) )
            e2 = 0.5_p_k_fld * ( emf%e%f2(2, i1, i2) + emf%e%f2(2, i1, i2-1) )
            e3 = emf%e%f2(3, i1, i2)

            !b1 = 0.5_p_k_fld * ( emf%b%f2( 1, i1, i2  ) + emf%b%f2( 1, i1, i2-1)  )

            b2 = 0.5_p_k_fld * ( emf%b%f2( 2, i1, i2  ) + emf%b%f2( 2, i1-1, i2) )

            b3 = 0.25_p_k_fld * ( emf%b%f2( 3, i1, i2  ) + emf%b%f2( 3, i1-1, i2) + &
                                  emf%b%f2( 3, i1, i2-1) + emf%b%f2( 3, i1-1, i2-1) )

            S%f2(1, i1, i2) = e2*b3 - b2*e3   ! S1

            ! S%f2(1, i1, i2) = e3*b1 - e1*b3   ! S2
            ! S%f2(1, i1, i2) = e1*b2 - b1*e2   ! S3

         enddo
       enddo

    case (2)

       do i2 = 1, emf%e%nx_(2)
         do i1 = 1, emf%e%nx_(1)

            e1 = 0.5_p_k_fld * ( emf%e%f2(1, i1, i2) + emf%e%f2(1, i1-1, i2) )
            !e2 = 0.5_p_k_fld * ( emf%e%f2(2, i1, i2) + emf%e%f2(2, i1, i2-1) )
            e3 = emf%e%f2(3, i1, i2)

            b1 = 0.5_p_k_fld * ( emf%b%f2( 1, i1, i2  ) + emf%b%f2( 1, i1, i2-1)  )

            !b2 = 0.5_p_k_fld * ( emf%b%f2( 2, i1, i2  ) + emf%b%f2( 2, i1-1, i2) )

            b3 = 0.25_p_k_fld * ( emf%b%f2( 3, i1, i2  ) + emf%b%f2( 3, i1-1, i2) + &
                                  emf%b%f2( 3, i1, i2-1) + emf%b%f2( 3, i1-1, i2-1) )

            ! S%f2(1, i1, i2) = e2*b3 - b2*e3   ! S1
            S%f2(1, i1, i2) = e3*b1 - e1*b3   ! S2
            ! S%f2(1, i1, i2) = e1*b2 - b1*e2   ! S3

         enddo
       enddo

    case (3)

       do i2 = 1, emf%e%nx_(2)
         do i1 = 1, emf%e%nx_(1)

            e1 = 0.5_p_k_fld * ( emf%e%f2(1, i1, i2) + emf%e%f2(1, i1-1, i2) )
            e2 = 0.5_p_k_fld * ( emf%e%f2(2, i1, i2) + emf%e%f2(2, i1, i2-1) )
            ! e3 = emf%e%f2(3, i1, i2)

            b1 = 0.5_p_k_fld * ( emf%b%f2( 1, i1, i2  ) + emf%b%f2( 1, i1, i2-1)  )

            b2 = 0.5_p_k_fld * ( emf%b%f2( 2, i1, i2  ) + emf%b%f2( 2, i1-1, i2) )

            ! b3 = 0.25_p_k_fld * ( emf%b%f2( 3, i1, i2  ) + emf%b%f2( 3, i1-1, i2) + &
            !                         emf%b%f2( 3, i1, i2-1) + emf%b%f2( 3, i1-1, i2-1) )

            ! S%f2(1, i1, i2) = e2*b3 - b2*e3   ! S1
            ! S%f2(1, i1, i2) = e3*b1 - e1*b3   ! S2
            S%f2(1, i1, i2) = e1*b2 - b1*e2   ! S3

         enddo
       enddo

  end select


end subroutine poynting_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine poynting_3d( emf, comp, S )

  implicit none

  class( t_emf ), intent(in) :: emf
  integer, intent(in) :: comp
  type( t_vdf ), intent(inout) :: S

  integer :: i1, i2, i3
  real( p_k_fld ) :: e1, e2, e3, b1, b2, b3

  select case ( comp )

    case (1)

       do i3 = 1, emf%e%nx_(3)
         do i2 = 1, emf%e%nx_(2)
           do i1 = 1, emf%e%nx_(1)

              !e1 = 0.5_p_k_fld * ( emf%e%f3(1, i1, i2, i3) + emf%e%f3(1, i1-1, i2, i3) )
              e2 = 0.5_p_k_fld * ( emf%e%f3(2, i1, i2, i3) + emf%e%f3(2, i1, i2-1, i3) )
              e3 = 0.5_p_k_fld * ( emf%e%f3(3, i1, i2, i3) + emf%e%f3(3, i1, i2, i3-1) )

              !b1 = 0.25_p_k_fld * ( emf%b%f3( 1, i1, i2, i3  ) + emf%b%f3( 1, i1, i2-1, i3) + &
              !                      emf%b%f3( 1, i1, i2, i3-1) + emf%b%f3( 1, i1, i2-1, i3-1) )

              b2 = 0.25_p_k_fld * ( emf%b%f3( 2, i1, i2, i3  ) + emf%b%f3( 2, i1-1, i2, i3) + &
                                    emf%b%f3( 2, i1, i2, i3-1) + emf%b%f3( 2, i1-1, i2, i3-1) )

              b3 = 0.25_p_k_fld * ( emf%b%f3( 3, i1, i2, i3  ) + emf%b%f3( 3, i1-1, i2, i3) + &
                                    emf%b%f3( 3, i1, i2-1, i3) + emf%b%f3( 3, i1-1, i2-1, i3) )

              S%f3(1, i1, i2, i3) = e2*b3 - b2*e3   ! S1
              ! S%f3(1, i1, i2, i3) = e3*b1 - e1*b3   ! S2
              ! S%f3(1, i1, i2, i3) = e1*b2 - b1*e2   ! S3

           enddo
         enddo
       enddo

    case (2)

       do i3 = 1, emf%e%nx_(3)
         do i2 = 1, emf%e%nx_(2)
           do i1 = 1, emf%e%nx_(1)

              e1 = 0.5_p_k_fld * ( emf%e%f3(1, i1, i2, i3) + emf%e%f3(1, i1-1, i2, i3) )
              ! e2 = 0.5_p_k_fld * ( emf%e%f3(2, i1, i2, i3) + emf%e%f3(2, i1, i2-1, i3) )
              e3 = 0.5_p_k_fld * ( emf%e%f3(3, i1, i2, i3) + emf%e%f3(3, i1, i2, i3-1) )

              b1 = 0.25_p_k_fld * ( emf%b%f3( 1, i1, i2, i3  ) + emf%b%f3( 1, i1, i2-1, i3) + &
                                    emf%b%f3( 1, i1, i2, i3-1) + emf%b%f3( 1, i1, i2-1, i3-1) )

              ! b2 = 0.25_p_k_fld * ( emf%b%f3( 2, i1, i2, i3  ) + emf%b%f3( 2, i1-1, i2, i3) + &
              !                       emf%b%f3( 2, i1, i2, i3-1) + emf%b%f3( 2, i1-1, i2, i3-1) )

              b3 = 0.25_p_k_fld * ( emf%b%f3( 3, i1, i2, i3  ) + emf%b%f3( 3, i1-1, i2, i3) + &
                                    emf%b%f3( 3, i1, i2-1, i3) + emf%b%f3( 3, i1-1, i2-1, i3) )

              ! S%f3(1, i1, i2, i3) = e2*b3 - b2*e3  ! S1
              S%f3(1, i1, i2, i3) = e3*b1 - e1*b3  ! S2
              ! S%f3(1, i1, i2, i3) = e1*b2 - b1*e2  ! S3

           enddo
         enddo
       enddo

    case (3)

       do i3 = 1, emf%e%nx_(3)
         do i2 = 1, emf%e%nx_(2)
           do i1 = 1, emf%e%nx_(1)

              e1 = 0.5_p_k_fld * ( emf%e%f3(1, i1, i2, i3) + emf%e%f3(1, i1-1, i2, i3) )
              e2 = 0.5_p_k_fld * ( emf%e%f3(2, i1, i2, i3) + emf%e%f3(2, i1, i2-1, i3) )
              ! e3 = 0.5_p_k_fld * ( emf%e%f3(3, i1, i2, i3) + emf%e%f3(3, i1, i2, i3-1) )

              b1 = 0.25_p_k_fld * ( emf%b%f3( 1, i1, i2, i3  ) + emf%b%f3( 1, i1, i2-1, i3) + &
                                    emf%b%f3( 1, i1, i2, i3-1) + emf%b%f3( 1, i1, i2-1, i3-1) )

              b2 = 0.25_p_k_fld * ( emf%b%f3( 2, i1, i2, i3  ) + emf%b%f3( 2, i1-1, i2, i3) + &
                                    emf%b%f3( 2, i1, i2, i3-1) + emf%b%f3( 2, i1-1, i2, i3-1) )

              ! b3 = 0.25_p_k_fld * ( emf%b%f3( 3, i1, i2, i3  ) + emf%b%f3( 3, i1-1, i2, i3) + &
              !                       emf%b%f3( 3, i1, i2-1, i3) + emf%b%f3( 3, i1-1, i2-1, i3) )

              !S%f3(1, i1, i2, i3) = e2*b3 - b2*e3  ! S1
              !S%f3(1, i1, i2, i3) = e3*b1 - e1*b3  ! S2
              S%f3(1, i1, i2, i3) = e1*b2 - b1*e2  ! S3

           enddo
         enddo
       enddo

  end select


end subroutine poynting_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine poynting_emf( emf, comp, S )

  implicit none

  class( t_emf ), intent(in) :: emf
  integer, intent(in) :: comp
  type( t_vdf ), intent(inout) :: S

  select case ( p_x_dim )
    case (1)
      call poynting_1d( emf, comp, S )
    case (2)
      call poynting_2d( emf, comp, S )
    case (3)
      call poynting_3d( emf, comp, S )
  end select

end subroutine poynting_emf
!-----------------------------------------------------------------------------------------


end module m_emf_poynting
