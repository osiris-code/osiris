!
! Numerical Cherenkonv Instability filter
!
! Based on B. B. Godfrey and J.-L. Vay, “Suppressing the numerical Cherenkov instability
! in FDTD PIC codes,” Journal of Computational Physics, vol. 267, pp. 1–6, Jun. 2014.
!
! This is a special filter that is applied in the x1 direction only, and only to fields
! seen by the particles

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_ncifilter

#include "memory/memory.h"

use m_system
use m_parameters
use m_vdf_define

implicit none

private

interface filter_nci
  module procedure filter_nci
end interface

public :: filter_nci

contains

!-----------------------------------------------------------------------------------------
subroutine filter_nci( e, b, dt, grid_center )

  implicit none

  type( t_vdf ), intent(inout) :: e, b
  real(p_double), intent(in) :: dt
  logical, intent(in) :: grid_center

  real(p_k_fld), dimension(4) :: s0, s1

  ! Get stencil values
  if ( grid_center ) then
    call get_stencil_gcenter( dt / e%dx_(1), s0, s1 )
  else
    call get_stencil_std( dt / e%dx_(1), s0, s1 )
  endif

  select case ( e%x_dim_ )
    case(1)
      call filter_nci_1d( e, b, s0, s1 )
    case(2)
      call filter_nci_2d( e, b, s0, s1 )
    case(3)
      call filter_nci_3d( e, b, s0, s1 )
  end select

end subroutine filter_nci

! ----------------------------------------------------------------------------------------
!
! Since calculating the coefficients for these filters is non-trivial, we approximate then
! using rational interpolation, based on precomputed values in the range dt_dz ]0,1[
!

! ----------------------------------------------------------------------------------------
! Get stencil values for standard (Yee mesh) interpolation
subroutine get_stencil_std( dt_dz, s0, s1 )

  implicit none

  real(p_double), intent(in) :: dt_dz
  real(p_k_fld), intent(inout), dimension(:) :: s0, s1

! Type: Standard



real(p_double), parameter, dimension(5,2,4) :: S0c = &
reshape( (/ 0.24608734804650872,2.9260848740735863,-5.06079044592636,0.8467390609562634,1.0431576814974235, &
1.0,11.890449463592006,-20.56752256495116,3.419799879497752,4.262407653667142, &
0.2201439658615669,0.323578408711018,-0.9298518396154316,0.4379701352570405,-0.004124703774763843, &
1.0,1.4698878776257118,-4.254394304931845,1.955532556644995,0.02993745437789668, &
0.14234586188779313,0.3530875573285242,-0.44487578670364747,0.23690453916100165,0.0005461439568704511, &
1.0,2.480830789310552,-3.3696501423806433,1.1586665617871381,0.2611678751122907, &
0.02732404938998069,-0.03028333261465197,0.08743928622469097,-0.07921251317134334,0.04057225234047862, &
1.0,-1.1082731082004649,0.7828517878776163,-0.21279057027998521,0.1187716702999842/), (/ 5, 2, 4 /) )

real(p_double), parameter, dimension(5,2,4) :: S1c = &
reshape( (/ 0.2463424226521974,0.08152548339290902,0.2660683750736488,0.035893801469176305,-0.01704278106197646, &
1.0,0.3309437187162323,1.081437617467492,0.14615101007188355,-0.06800856184645067, &
0.22319809633636536,0.062312486807837617,0.22435867089678882,0.024107788374300953,-0.02333347571350356, &
1.0,0.27917992363120187,1.0227154443405428,0.11282683668920913,-0.08983628263142367, &
0.15663481753118244,0.05399522779117138,0.15261027852236692,0.017482037790291914,-0.03315910178131886, &
1.0,0.3447175589519361,1.0919457441333944,0.15136713566551208,-0.10264996826382146, &
0.05139073363640456,-0.017896674471897255,0.0029349299604851005,0.0039622915500230425,-0.008821526982677687, &
1.0,-0.34824007651716693,0.7844178098197293,-0.17438686824115085,0.09392358532612638/), (/ 5, 2, 4 /) )

  real(p_double) :: x, n, d
  integer :: i, j

  x = dt_dz

  do j = 1, 4
    ! E2 filters
    n = S0c(4,1,j) + S0c(5,1,j) * x
    d = S0c(4,2,j) + S0c(5,2,j) * x

    do i = 3,1,-1
      n = S0c(i,1,j) + x*n
      d = S0c(i,2,j) + x*d
    enddo

    s0( j ) = n/d

    ! E1, B3 filters
      n = S1c(4,1,j) + S1c(5,1,j) * x
    d = S1c(4,2,j) + S1c(5,2,j) * x

    do i = 3,1,-1
      n = S1c(i,1,j) + x*n
      d = S1c(i,2,j) + x*d
    enddo

    s1( j ) = n/d
  enddo

end subroutine get_stencil_std

! ----------------------------------------------------------------------------------------
! Get stencil values for field values centered on the grid interpolation
subroutine get_stencil_gcenter( dt_dz, s0, s1 )

  implicit none

  real(p_double), intent(in) :: dt_dz
  real(p_k_fld), intent(inout), dimension(:) :: s0, s1

! Type: Grid Center

real(p_double), parameter, dimension(5,2,4) :: S0c = &
reshape( (/ 0.24929613438035986,-0.3520583738750654,0.43943904975326104,-0.6281382126400855,0.30018054608902855, &
            1.0,-1.4122120354564773,1.7629737524535756,-2.5204917519047894,1.204606206731683, &
            0.2372290592015112,-0.34309764076269955,0.4550841878390453,-0.5964892195620909,0.27451025140840796, &
            1.0,-1.4463579735495258,1.9273407143992771,-2.5441698970285347,1.1740911749442848, &
            0.18652734031519105,-0.26402043344082815,0.4741475298025892,-0.49692304315303754,0.21081964528179398, &
            1.0,-1.4161813371922989,2.6273840533922495,-2.9297466586447265,1.2358114427986138, &
            0.07387939622941167,0.1129789035519798,1.0708883312399595,-0.8734699673390198,2.0623555999839893, &
            1.0,1.5524067362944998,14.479184310198997,-5.909380097417019,10.146564352489142/), (/ 5, 2, 4 /) )

real(p_double), parameter, dimension(5,2,4) :: S1c = &
reshape( (/ 0.24626263786985514,-0.0012730665758197452,1.963114810196638,2.664899632601468,-0.7825556524396512, &
            1.0,-0.005179341229568304,7.976638542460345,10.818157513605255,-3.1332880242792425, &
            0.222259959685894,0.08743946193807703,1.6442937364397208,1.3011676005487618,-0.5232972337742652, &
            1.0,0.3933714979484943,7.461925704064187,5.864187141135722,-1.9371002887787667, &
            0.15233817877446523,24.836603045468713,31.770680985403462,-58.81484785453512,17.097835634752016, &
            1.0,163.31191607925186,205.53436032597693,-277.43441863231743,38.3545301442178, &
            0.042861119763267806,-0.009272215092752915,0.08314051009302832,-0.08732623070482323,-0.029403356281416364, &
            1.0,-0.2139832577549167,5.118083037615773,-2.1421112550846817,1.306484477199788/), (/ 5, 2, 4 /) )

  real(p_double) :: x, n, d
  integer :: i, j

  x = dt_dz

  do j = 1, 4
    ! E2 filters
    n = S0c(4,1,j) + S0c(5,1,j) * x
    d = S0c(4,2,j) + S0c(5,2,j) * x

    do i = 3,1,-1
      n = S0c(i,1,j) + x*n
      d = S0c(i,2,j) + x*d
    enddo

    s0( j ) = n/d

    ! E1, B3 filters
      n = S1c(4,1,j) + S1c(5,1,j) * x
    d = S1c(4,2,j) + S1c(5,2,j) * x

    do i = 3,1,-1
      n = S1c(i,1,j) + x*n
      d = S1c(i,2,j) + x*d
    enddo

    s1( j ) = n/d
  enddo

end subroutine get_stencil_gcenter

!-----------------------------------------------------------------------------------------
! Filtering is done using multiple passages of a 3 point stencil

!-----------------------------------------------------------------------------------------
! 1D filter  ! there is no nci in 1d
!-----------------------------------------------------------------------------------------
subroutine filter_nci_1d( e, b, s0, s1 )

   implicit none

   type(t_vdf), intent(inout) :: e, b
   real(p_k_fld), dimension(:), intent(in) :: s0, s1

   ! Point field values
   real(p_k_fld) :: af1, af2, af3
   real(p_k_fld) :: bf1, bf2, bf3

   real(p_k_fld) :: ac13, ac2
   real(p_k_fld) :: bc13, bc2

   integer :: i1, j, lb, ub

   ! limits for applying the filter
   lb = lbound( e%f1, 2 ) + 1
   ub = ubound( e%f1, 2 ) - 1

   ! Filter E field

   ! do 4 passes
   do j = 1, 4
     ! Get stencil coefficients
     ac13 = s1( j )
     ac2  = 1 - 2*s1( j )

     bc13 = s0( j )
     bc2  = 1 - 2*s0( j )

     ! read first 2 points
     af1 = e%f1( 1, lb-1 )
     bf1 = e%f1( 2, lb-1 )

     af2 = e%f1( 1, lb   )
     bf2 = e%f1( 2, lb   )

     do i1 = lb, ub
       ! read the next data point from memory
       af3 = e%f1( 1, i1 + 1 )
       bf3 = e%f1( 2, i1 + 1 )

       ! store filtered data
       e%f1( 1, i1 ) = ac13*af1 + ac2*af2 + ac13*af3
       e%f1( 2, i1 ) = bc13*bf1 + bc2*bf2 + bc13*bf3

       ! Shift left existing data
       af1 = af2; bf1 = bf2
       af2 = af3; bf2 = bf3
     enddo
   enddo

   ! Filter B field

   ! do 4 passes
   do j = 1, 4
     ! Get stencil coefficients
     ac13 = s1( j )
     ac2  = 1 - 2*s1( j )

     ! read first 2 points
     af1 = b%f1( 3, lb-1 )
     af2 = b%f1( 3, lb   )

     do i1 = lb, ub
       ! read the next data point from memory
       af3 = b%f1( 3, i1 + 1 )

       ! store filtered data
       b%f1( 3, i1 ) = ac13*af1 + ac2*af2 + ac13*af3

       ! Shift left existing data
       af1 = af2
       af2 = af3
     enddo
   enddo

end subroutine filter_nci_1d


!-----------------------------------------------------------------------------------------
! 2D filter
!-----------------------------------------------------------------------------------------
subroutine filter_nci_2d( e, b, s0, s1 )

   implicit none

   type(t_vdf), intent(inout) :: e, b
   real(p_k_fld), dimension(:), intent(in) :: s0, s1

   ! Point field values
   real(p_k_fld) :: af1, af2, af3
   real(p_k_fld) :: bf1, bf2, bf3
   real(p_k_fld) :: cf1, cf2, cf3

   real(p_k_fld) :: ac13, ac2
   real(p_k_fld) :: bc13, bc2

   integer :: i1, i2, j, lb, ub

   ! limits for applying the filter
   lb = lbound( e%f2, 2 ) + 1
   ub = ubound( e%f2, 2 ) - 1

   ! Filter E field
   do i2 = lbound( e%f2, 3 ), ubound( e%f2, 3 )
     ! do 4 passes
     do j = 1, 4
       ! Get stencil coefficients
       ac13 = s1( j )
       ac2  = 1 - 2*s1( j )

       bc13 = s0( j )
       bc2  = 1 - 2*s0( j )

       ! read first 2 points
       af1 = e%f2( 1, lb-1, i2 )
       bf1 = e%f2( 2, lb-1, i2 )
       cf1 = e%f2( 3, lb-1, i2 )

       af2 = e%f2( 1, lb  , i2 )
       bf2 = e%f2( 2, lb  , i2 )
       cf2 = e%f2( 3, lb  , i2 )

       do i1 = lb, ub
         ! read the next data point from memory
         af3 = e%f2( 1, i1 + 1, i2 )
         bf3 = e%f2( 2, i1 + 1, i2 )
         cf3 = e%f2( 3, i1 + 1, i2 )

         ! store filtered data
         e%f2( 1, i1, i2 ) = ac13*af1 + ac2*af2 + ac13*af3
         e%f2( 2, i1, i2 ) = bc13*bf1 + bc2*bf2 + bc13*bf3
             e%f2( 3, i1, i2 ) = bc13*cf1 + bc2*cf2 + bc13*cf3
         ! Shift left existing data
         af1 = af2; bf1 = bf2; cf1 = cf2;
         af2 = af3; bf2 = bf3; cf2 = cf3;
       enddo
     enddo
   enddo

   ! Filter B field
   do i2 = lbound( b%f2, 3 ), ubound( b%f2, 3 )
     ! do 4 passes
     do j = 1, 4
       ! Get stencil coefficients
       ac13 = s1( j )
       ac2  = 1 - 2*s1( j )

       ! read first 2 points
       af1 = b%f2( 1, lb-1, i2 )
       bf1 = b%f2( 2, lb-1, i2 )
       cf1 = b%f2( 3, lb-1, i2 )
       af2 = b%f2( 1, lb  , i2 )
       bf2 = b%f2( 2, lb  , i2 )
       cf2 = b%f2( 3, lb  , i2 )

       do i1 = lb, ub
         ! read the next data point from memory
         af3 = b%f2( 1, i1 + 1, i2 )
         bf3 = b%f2( 2, i1 + 1, i2 )
         cf3 = b%f2( 3, i1 + 1, i2 )

         ! store filtered data
         b%f2( 1, i1, i2 ) = ac13*af1 + ac2*af2 + ac13*af3
          b%f2( 2, i1, i2 ) = ac13*bf1 + ac2*bf2 + ac13*bf3
         b%f2( 3, i1, i2 ) = ac13*cf1 + ac2*cf2 + ac13*cf3
         ! Shift left existing data
         af1 = af2; bf1 = bf2; cf1 = cf2;
         af2 = af3; bf2 = bf3; cf2 = cf3;
       enddo
     enddo
   enddo

end subroutine filter_nci_2d

!-----------------------------------------------------------------------------------------
! 3D filter
!-----------------------------------------------------------------------------------------
subroutine filter_nci_3d( e, b, s0, s1 )

   implicit none

   type(t_vdf), intent(inout) :: e, b
   real(p_k_fld), dimension(:), intent(in) :: s0, s1

   ! Point field values
   real(p_k_fld) :: af1, af2, af3
   real(p_k_fld) :: bf1, bf2, bf3
   real(p_k_fld) :: cf1, cf2, cf3

   real(p_k_fld) :: ac13, ac2
   real(p_k_fld) :: bc13, bc2

   integer :: i1, i2, i3, j, lb, ub

   ! limits for applying the filter
   lb = lbound( e%f3, 2 ) + 1
   ub = ubound( e%f3, 2 ) - 1

   ! Filter E field
   do i3 = lbound( e%f3, 4 ), ubound( e%f3, 4 )
     do i2 = lbound( e%f3, 3 ), ubound( e%f3, 3 )
       ! do 4 passes
       do j = 1, 4
         ! Get stencil coefficients
         ac13 = s1( j )
         ac2  = 1 - 2*s1( j )

         bc13 = s0( j )
         bc2  = 1 - 2*s0( j )

         ! read first 2 points
         af1 = e%f3( 1, lb-1, i2, i3 )
         bf1 = e%f3( 2, lb-1, i2, i3 )
           cf1 = e%f3( 3, lb-1, i2, i3 )

         af2 = e%f3( 1, lb  , i2, i3 )
         bf2 = e%f3( 2, lb  , i2, i3 )
         cf2 = e%f3( 3, lb  , i2, i3 )

         do i1 = lb, ub
           ! read the next data point from memory
           af3 = e%f3( 1, i1 + 1, i2, i3 )
           bf3 = e%f3( 2, i1 + 1, i2, i3 )
           cf3 = e%f3( 3, i1 + 1, i2, i3 )

           ! store filtered data
           e%f3( 1, i1, i2, i3 ) = ac13*af1 + ac2*af2 + ac13*af3
           e%f3( 2, i1, i2, i3 ) = bc13*bf1 + bc2*bf2 + bc13*bf3
           e%f3( 3, i1, i2, i3 ) = bc13*cf1 + bc2*cf2 + bc13*cf3
           ! Shift left existing data
           af1 = af2; bf1 = bf2; cf1 = cf2;
           af2 = af3; bf2 = bf3; cf2 = cf3;
         enddo
       enddo
     enddo
   enddo

   ! Filter B field
   do i3 = lbound( b%f3, 4 ), ubound( b%f3, 4 )
     do i2 = lbound( b%f3, 3 ), ubound( b%f3, 3 )
       ! do 4 passes
       do j = 1, 4
         ! Get stencil coefficients
         ac13 = s1( j )
         ac2  = 1 - 2*s1( j )

         ! read first 2 points
         af1 = b%f3( 1, lb-1, i2, i3 )
         bf1 = b%f3( 2, lb-1, i2, i3 )
         cf1 = b%f3( 3, lb-1, i2, i3 )
         af2 = b%f3( 1, lb  , i2, i3 )
         bf2 = b%f3( 2, lb  , i2, i3 )
         cf2 = b%f3( 3, lb  , i2, i3 )

         do i1 = lb, ub
           ! read the next data point from memory
           af3 = b%f3( 1, i1 + 1, i2, i3 )
           bf3 = b%f3( 2, i1 + 1, i2, i3 )
           cf3 = b%f3( 3, i1 + 1, i2, i3 )

           ! store filtered data
           b%f3( 1, i1, i2, i3 ) = ac13*af1 + ac2*af2 + ac13*af3
            b%f3( 2, i1, i2, i3 ) = ac13*bf1 + ac2*bf2 + ac13*bf3
           b%f3( 3, i1, i2, i3 ) = ac13*cf1 + ac2*cf2 + ac13*cf3
           ! Shift left existing data
           af1 = af2; bf1 = bf2; cf1 = cf2;
           af2 = af3; bf2 = bf3; cf2 = cf3;
         enddo
       enddo
     enddo
   enddo

end subroutine filter_nci_3d



end module m_emf_ncifilter
