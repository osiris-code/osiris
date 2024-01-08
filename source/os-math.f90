! Mathematical functions module

module m_math

  use m_system

  implicit none
  
  private

  ! Mathematical constants
  
  real(p_double), parameter :: pi      = 3.14159265358979323846264_p_double    ! pi
  real(p_double), parameter :: pi_2    = 1.57079632679489661923132_p_double    ! pi/2
  real(p_double), parameter :: pi_180  = 0.0174532925199432957692369_p_double  ! pi/180
  real(p_double), parameter :: sqrt_2  = 1.41421356237309504880169_p_double    ! sqrt(2)

  interface fenv_poly
    module procedure fenv_poly
  end interface

  interface fenv_sin2
    module procedure fenv_sin2
  end interface

  interface besselj0
    module procedure besselj0_rsing
    !module procedure besselj0_rdoub    ! not implemented yet
  end interface

  interface besselj1
    module procedure besselj1_rsing
    !module procedure besselj1_rdoub    ! not implemented yet
  end interface

  interface besseljn
    module procedure besseljn_rsing
    !module procedure besseljn_rdoub    ! not implemented yet
  end interface

  interface besselJnZero
    module procedure besselJnZero  
  end interface

  interface besseli0
    !module procedure besseli0_rsing     ! not implemented yet
    module procedure besseli0_rdoub
  end interface

  interface hermite
    module procedure hermite_single  
    module procedure hermite_double
  end interface
  
  interface laguerre
    module procedure laguerre_single  
    module procedure laguerre_double
  end interface
  
  interface smooth_heaviside
    module procedure smooth_heaviside_single
    module procedure smooth_heaviside_double
  end interface
  
  interface grad_linear
    module procedure gradient_linear
  end interface

  interface grad_regressive
    module procedure gradient_regressive
  end interface

  public :: pi, pi_2, pi_180, sqrt_2
  public :: fenv_poly, fenv_sin2
  public :: besselj0, besselj1, besseljn, besselJnZero, besseli0
  public :: hermite, laguerre, smooth_heaviside, grad_linear, grad_regressive

!-----------------------------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function fenv_poly(x,rise,flat,fall)
!-----------------------------------------------------------------------------------------
!       polynomial envelope function
!-----------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  real(p_double) :: fenv_poly

  real(p_double), intent(in) :: x, rise, flat, fall

!       local variables

  real(p_double) :: length

!       executable statements

  length = rise + flat + fall

  if ( x < 0.0 ) then
   fenv_poly = 0.0
  else if ( x < rise ) then
   fenv_poly = fenv( x/rise )
  else if ( x < rise + flat ) then
   fenv_poly = 1.0
  else if ( x < length ) then
  fenv_poly = fenv( (length-x)/fall )
  else
  fenv_poly = 0.0
  endif

  contains

  ! function declaration of polynomial profile

  function fenv(tt)
  implicit none
  real(p_double) :: fenv
  real(p_double), intent(in) :: tt

  fenv = 10 * tt**3 - 15 * tt**4 + 6 * tt**5
  end function fenv

end function fenv_poly
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function fenv_sin2(x,rise,flat,fall)
!-----------------------------------------------------------------------------------------
!       sin^2 envelope function
!-----------------------------------------------------------------------------------------

  implicit none

!       dummy variables

  real(p_double) :: fenv_sin2

  real(p_double), intent(in) :: x, rise, flat, fall

!       local variables

  real(p_double) :: length

!       executable statements

  length = rise + flat + fall

  if ( x < 0.0 ) then
   fenv_sin2 = 0.0
  else if ( x < rise ) then
   fenv_sin2 = sin( real(pi_2,p_double) * x/rise )**2
  else if ( x < rise + flat ) then
   fenv_sin2 = 1.0
  else if ( x < length ) then
  fenv_sin2 = sin( real(pi_2,p_double)*(length-x)/fall ) ** 2
  else
  fenv_sin2 = 0.0
  endif

end function fenv_sin2
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
function besselJnZero(n,i)
!-----------------------------------------------------------------------------------------
! returns the ith zero for the Jn function
! Values were computed using Mathematica 4.2
! << NumericalMath`BesselZeros`
! BesselJZeros[5, 5, WorkingPrecision -> 16]
! if n or i fall outside the computed range the function
! returns 0.0
!-----------------------------------------------------------------------------------------
  implicit none
  
  ! dummy variables
  integer, intent(in) :: n, i
  real(p_single) :: besselJnZero
  
  ! local variables
  integer, parameter :: NMAX = 5, IMAX = 5
  
  
  real(p_single), dimension(0:NMAX, IMAX), parameter :: jzeros = reshape( &
       ! J0       J1       J2       J3       J4      J5
  
    (/ (/2.40483, 3.83171, 5.13562, 6.38016, 7.58834, 8.77148/), &        ! 1st zero
       (/5.52008, 7.01559, 8.41724, 9.76102, 11.0647, 12.3386/), &        ! 2nd zero
       (/8.65373, 10.1735, 11.6198, 13.0152, 14.3725, 15.7002/), &        ! 3rd zero
       (/11.7915, 13.3237, 16.2235, 16.2235, 17.6160, 18.9801/), &        ! 4th zero
       (/14.9309, 14.3725, 19.4094, 19.4094, 20.8269, 22.2178/) /), &     ! 5th zero
             (/NMAX+1, IMAX/) )
  
  if (n>NMAX .or. i>IMAX) then  
    besselJnZero = 0.0_p_single
  else
    besselJnZero =  jzeros(n,i)
  endif

end function besselJnZero
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function besselj0_rsing(x)
!-----------------------------------------------------------------------------------------
! Returns the Bessel function J0(x) for any real x.
! Adapted from 'Numerical Recipes in Fortran 90'
!-----------------------------------------------------------------------------------------

  implicit none
  
  ! dummy variables
  real(p_single), intent(in) :: x
  real(p_single) :: besselj0_rsing
  
  ! local variables
  real(p_single) :: ax,xx,z
  real(p_double) :: y  
  
  ! accumulate polynomials in double precision.
  real(p_double), dimension(5), parameter :: &
             p = (/1.0_p_double             ,-0.1098628627e-2_p_double,&
                   0.2734510407e-4_p_double ,-0.2073370639e-5_p_double,&
                   0.2093887211e-6_p_double/)
  real(p_double), dimension(5), parameter :: &
             q = (/-0.1562499995e-1_p_double, 0.1430488765e-3_p_double,&
                   -0.6911147651e-5_p_double, 0.7621095161e-6_p_double,&
                   -0.934945152e-7_p_double/)
  real(p_double), dimension(6), parameter :: &
             r = (/57568490574.0_p_double   ,-13362590354.0_p_double,&
                     651619640.7_p_double   ,-11214424.18_p_double,&
                     77392.33017_p_double   ,-184.9052456_p_double/)
  real(p_double), dimension(6), parameter :: &
             s = (/57568490411.0_p_double   , 1029532985.0_p_double,&
                     9494680.718_p_double   , 59272.64853_p_double,&
                     267.8532712_p_double   , 1.0_p_double/)
  
  real(p_double) :: poly1, poly2, y1
  integer :: i
  
  ! executable statements  
  
  y1=x**2
  if (y1 <= 8.0_p_single*tiny(x)) then !Underflow limit.
    besselj0_rsing = 1.0
  else if (abs(x) < 8.0) then !Direct rational function fit.

  
  y = y1

!  besselj0 = poly(y,r) / poly(y,s)

  poly1 = r(1)*y
  poly2 = s(1)*y
  do i = 2, 6  
      y = y*y1
      poly1 = poly1 + r(i)*y
    poly2 = poly2 + s(i)*y
  enddo
  besselj0_rsing = real(poly1/poly2,p_single)

  else !Fitting function (6.5.9).
  ax=abs(x)
  z=8.0_p_single/ax
  y1=z**2
  xx=ax-0.785398164_p_single

!  besselj0=sqrt(0.636619772_p_single/ax)*(cos(xx)*poly(y,p)-z*sin(xx)*poly(y,q))

  y = y1
  poly1 = p(1)*y
  poly2 = q(1)*y
  do i = 2, 5  
      y = y*y1
      poly1 = poly1 + p(i)*y
    poly2 = poly2 + q(i)*y
  enddo
  besselj0_rsing = sqrt(0.636619772_p_single/ax)*&
                   (cos(xx)*real(poly1,p_single)-z*sin(xx)*real(poly2,p_single))

  end if

end function besselj0_rsing
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
function besselj1_rsing(x)
!-----------------------------------------------------------------------------------------
! Returns the Bessel function J1(x) for any real x.
! Adapted from 'Numerical Recipes in Fortran 90'
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  real(p_single), intent(in) :: x
  real(p_single) :: besselj1_rsing

  ! local variables
  real(p_single) :: ax,xx,z
  real(p_double) :: y , y1
  
  ! accumulate polynomials in double precision.
  real(p_double), dimension(6), parameter :: &
              r = (/72362614232.0_p_double, -7895059235.0_p_double, 242396853.1_p_double,&
                     -2972611.439_p_double,   15704.48260_p_double,-30.16036606_p_double/)
  real(p_double), dimension(6), parameter :: &
              s = (/144725228442.0_p_double,2300535178.0_p_double,  18583304.74_p_double,&
                       99447.43394_p_double, 376.9991397_p_double,          1.0_p_double/)
  real(p_double), dimension(5), parameter :: &
              p = (/ 1.0_p_double, 0.183105e-2_p_double,-0.3516396496e-4_p_double,&
                     0.2457520174e-5_p_double,-0.240337019e-6_p_double/)
  real(p_double), dimension(5), parameter :: &
              q = (/0.04687499995_p_double, -0.2002690873e-3_p_double,0.8449199096e-5_p_double,&
                   -0.88228987e-6_p_double,  0.105787412e-6_p_double/)

  real(p_double) :: poly1, poly2
  integer :: i
  
  ! executable statements  

  y1=real(x**2,p_double)
  if (y1 <= 8.0_p_single*tiny(x)) then !Underflow limit.
    besselj1_rsing = 0.0_p_single
  else if (abs(x) < 8.0_p_single) then !Direct rational approximation.
  y1=x**2
  y = y1
  !besselj1_rsing=x*(poly(y,r)/poly(y,s))

  poly1 = r(1)*y
  poly2 = s(1)*y
  do i = 2, 6  
      y = y*y1
      poly1 = poly1 + r(i)*y
    poly2 = poly2 + s(i)*y
  enddo
  besselj1_rsing = x*real(poly1/poly2, p_single)

  else !Fitting function (6.5.9).
  ax=abs(x)
  z=8.0_p_single/ax
  y1=z**2
  xx=ax-2.356194491_p_single
  
  !besselj1=sqrt(0.636619772_p_single/ax)*(cos(xx)*poly(y,p)-z*sin(xx)*poly(y,q))*&
  ! sign(1.0_p_single,x)
    y = y1
  poly1 = p(1)*y
  poly2 = q(1)*y
  do i = 2, 5  
      y = y*y1
      poly1 = poly1 + p(i)*y
    poly2 = poly2 + q(i)*y
  enddo
    besselj1_rsing=sqrt(0.636619772_p_single/ax)* &
                   (cos(xx)*real(poly1, p_single)-z*sin(xx)*real(poly2,p_single))* &
                 sign(1.0_p_single,x)
  end if
  
end function besselj1_rsing
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function besseljn_rsing(n,x)
!-----------------------------------------------------------------------------------------
! Returns the Bessel function Jn(x) for any real x. Make the parameter IACC
! larger to increase accuracy.
! Adapted from 'Numerical Recipes in Fortran 90'
!-----------------------------------------------------------------------------------------

  implicit none
  
  ! dummy variables
  integer, intent(in) :: n
  real(p_single), intent(in) :: x
  real(p_single) :: besseljn_rsing
  
  ! local variables
  integer, parameter :: IACC=40,IEXP=maxexponent(x)/2
  integer :: j,jsum,m
  real(p_single) :: ax,bj,bjm,bjp,summ,tox
  
  ! executable statements
  if (n == 0) then
  besseljn_rsing = besselj0_rsing(x)
  return
  else if (n == 1) then
  besseljn_rsing = besselj1_rsing(x)
  return
  endif
  
  ax=abs(x)
  if (ax*ax <= 8.0_p_single*tiny(x)) then !Underflow limit.
    besseljn_rsing=0.0_p_single
  else if (ax > real(n,p_single)) then !Upwards recurrence from J0 and J1.
  tox=2.0_p_single/ax
  bjm=besselj0_rsing(ax)
  bj=besselj1_rsing(ax)
  do j=1,n-1
  bjp=j*tox*bj-bjm
  bjm=bj
  bj=bjp
  end do
  besseljn_rsing=bj
  else !Downwards recurrence from an even m
  tox=2.0_p_single/ax ! here computed.  
  m=2*((n+int(sqrt(real(IACC*n,p_single))))/2)
  besseljn_rsing=0.0_p_single
  jsum=0  
  !jsum will alternate between 0 and 1; when
  !it is 1, we accumulate in sum the
  !even terms in (5.5.16).
  summ=0.0
  bjp=0.0
  bj=1.0
  do j=m,1,-1 !The downward recurrence.
    bjm=j*tox*bj-bjp
    bjp=bj
    bj=bjm
    if (exponent(bj) > IEXP) then !Renormalize to prevent overflows.
    bj=scale(bj,-IEXP)
    bjp=scale(bjp,-IEXP)
    besseljn_rsing=scale(besseljn_rsing,-IEXP)
    summ=scale(summ,-IEXP)
    end if
    if (jsum /= 0) summ=summ+bj !Accumulate the sum.
    jsum=1-jsum !Change 0 to 1 or vice versa.
    if (j == n) besseljn_rsing=bjp !Save the unnormalized answer.
  end do
  summ=2.0_p_single*summ-bj !Compute (5.5.16)
  besseljn_rsing=besseljn_rsing/summ !and use it to normalize the answer.
  end if
  
  if (x < 0.0 .and. mod(n,2) == 1) besseljn_rsing=-besseljn_rsing
end function besseljn_rsing
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! Modified Bessel function of the first kind I0(x) for any real x
! Adapted from 'Numerical Recipes in Fortran 90'
!-----------------------------------------------------------------------------------------
function besseli0_rdoub( x )  

  implicit none
  
  real(p_double), intent(in) :: x
  real(p_double) :: besseli0_rdoub
  
  real(p_double) :: ax
  real(p_double) :: y ! accumulate polynomials in double precision

  ax = abs(x)
  if ( ax < 3.75 ) then
    y = (x / 3.75)**2
    besseli0_rdoub = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 + &
                     y*(0.2659732 + y*(0.360768d-1 + y*0.45813d-2)))))
  else
    y = 3.75 / ax
    besseli0_rdoub = (exp(ax)/sqrt(ax)) * (0.39894228 + y*(0.1328592d-1 + &
                     y*(0.225319d-2 + y*(-0.157565d-2 + y*(0.916281d-2 + &
                     y*(-0.2057706d-1 + y*(0.2635537d-1 + y*(-0.1647633d-1 + y*0.392377d-2))))))))
  endif

end function besseli0_rdoub
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Physicists' Hermite polynomials - double precision
!-----------------------------------------------------------------------------------------
function hermite_double( mode, u )

  implicit none

  real(p_double) :: hermite_double
  integer, intent(in) :: mode
  real(p_double), intent(in) :: u

  select case (mode)
  case(0)
    hermite_double = 1.0_p_double
  case(1)
    hermite_double = 2.0_p_double * u
  case(2)
    hermite_double = 4.0_p_double * u**2   - 2.0_p_double
  case(3)
    hermite_double = 8.0_p_double * u**3   - 12.0_p_double* u
  case(4)
    hermite_double = 16.0_p_double * u**4  - 48.0_p_double* u**2 + 12.0_p_double
  case(5)
    hermite_double = 32.0_p_double * u**5  - 160.0_p_double* u**3 + 120.0_p_double * u
  case(6)
    hermite_double = 64.0_p_double * u**6  - 480.0_p_double* u**4 + 720.0_p_double * u**2 - &
      120.0_p_double
  case(7)
    hermite_double = 128.0_p_double * u**7 - 1344.0_p_double* u**5 + 3360.0_p_double * u**3 - &
      1680.0_p_double * u
  case default
    hermite_double = 0
  end select

end function hermite_double
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Physicists' Hermite polynomials - single precision
!-----------------------------------------------------------------------------------------
function hermite_single( mode, u )

  implicit none

  real(p_single) :: hermite_single
  integer, intent(in) :: mode
  real(p_single), intent(in) :: u

  select case (mode)
  case(0)
    hermite_single = 1.0_p_single
  case(1)
    hermite_single = 2.0_p_single * u
  case(2)
    hermite_single = 4.0_p_single * u**2   - 2.0_p_single
  case(3)
    hermite_single = 8.0_p_single * u**3   - 12.0_p_single* u
  case(4)
    hermite_single = 16.0_p_single * u**4  - 48.0_p_single* u**2 + 12.0_p_single
  case(5)
    hermite_single = 32.0_p_single * u**5  - 160.0_p_single* u**3 + 120.0_p_single * u
  case(6)
    hermite_single = 64.0_p_single * u**6  - 480.0_p_single* u**4 + 720.0_p_single * u**2 - &
      120.0_p_single
  case(7)
    hermite_single = 128.0_p_single * u**7 - 1344.0_p_single* u**5 + 3360.0_p_single * u**3 - &
      1680.0_p_single * u
  case default
    hermite_single = 0
  end select

end function hermite_single
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Generalized Laguerre polynomials - double precision
!-----------------------------------------------------------------------------------------
function laguerre_double( p, l, u )

  implicit none

  real(p_double) :: laguerre_double
  integer, intent(in) :: p, l
  real(p_double), intent(in) :: u

  select case (p)
  case(0)
    laguerre_double = 1.0_p_double
  case(1)
    laguerre_double = 1.0_p_double + l - u
  case(2)
    laguerre_double = (2 + 3.0_p_double* l + l**2 - 4.0_p_double* u - &
        2.0_p_double* l*u + u**2)/2.0_p_double
  case(3)
    laguerre_double = (6.0_p_double + 11.0_p_double* l + 6.0_p_double* l**2 + l**3 - &
          18.0_p_double* u - 15.0_p_double* l*u - 3.0_p_double* l**2*u + &
          9.0_p_double* u**2 + 3.0_p_double* l*u**2 - u**3)/6.0_p_double
  case(4)
    laguerre_double =   (24.0_p_double + 50.0_p_double* l + 35.0_p_double* l**2 + &
          10.0_p_double* l**3 + l**4 - 96.0_p_double* u - 104.0_p_double* l*u - &
          36.0_p_double* l**2*u - 4.0_p_double* l**3*u + 72.0_p_double* u**2 + &
          42.0_p_double* l*u**2 + 6.0_p_double* l**2*u**2 - 16.0_p_double* u**3 - &
          4.0_p_double* l*u**3 + u**4)/24.0_p_double
  case(5)
    laguerre_double = (120.0_p_double + 274.0_p_double* l + 225.0_p_double* l**2 + &
          85.0_p_double* l**3 + 15.0_p_double* l**4 + l**5 - 600.0_p_double* u - &
          770.0_p_double* l*u - 355.0_p_double* l**2*u - 70.0_p_double* l**3*u - &
          5.0_p_double* l**4*u + 600.0_p_double* u**2 + 470.0_p_double* l*u**2 + &
          120.0_p_double* l**2*u**2 + 10.0_p_double* l**3*u**2 - &
          200.0_p_double* u**3 - 90.0_p_double* l*u**3 - 10.0_p_double* l**2*u**3 + &
          25.0_p_double* u**4 + 5.0_p_double* l*u**4 - u**5)/120.0_p_double
  case(6)
    laguerre_double = (720.0_p_double + 1764.0_p_double* l + 1624.0_p_double* l**2 + &
          735.0_p_double* l**3 + 175.0_p_double* l**4 + 21.0_p_double* l**5 + &
          l**6 - 4320.0_p_double* u - 6264.0_p_double* l*u - &
          3480.0_p_double* l**2*u - 930.0_p_double* l**3*u - &
          120.0_p_double* l**4*u - 6.0_p_double* l**5*u + 5400.0_p_double* u**2 + &
          5130.0_p_double* l*u**2 + 1785.0_p_double* l**2*u**2 + &
          270.0_p_double* l**3*u**2 + 15.0_p_double* l**4*u**2 - &
          2400.0_p_double* u**3 - 1480.0_p_double* l*u**3 - &
          300.0_p_double* l**2*u**3 - 20.0_p_double* l**3*u**3 + &
          450.0_p_double* u**4 + 165.0_p_double* l*u**4 + 15.0_p_double* l**2*u**4 - &
          36.0_p_double* u**5 - 6.0_p_double* l*u**5 + u**6)/720.0_p_double
  case(7)
    laguerre_double =  (5040.0_p_double + 13068.0_p_double* l + 13132.0_p_double* l**2 + &
          6769.0_p_double* l**3 + 1960.0_p_double* l**4 + 322.0_p_double* l**5 + &
          28.0_p_double* l**6 + l**7 - 35280.0_p_double* u - &
          56196.0_p_double* l*u - 35728.0_p_double* l**2*u - &
          11655.0_p_double* l**3*u - 2065.0_p_double* l**4*u - &
          189.0_p_double* l**5*u - 7.0_p_double* l**6*u + 52920.0_p_double* u**2 + &
          57834.0_p_double* l*u**2 + 24675.0_p_double* l**2*u**2 + &
          5145.0_p_double* l**3*u**2 + 525.0_p_double* l**4*u**2 + &
          21.0_p_double* l**5*u**2 - 29400.0_p_double* u**3 - &
          22330.0_p_double* l*u**3 - 6265.0_p_double* l**2*u**3 - &
          770.0_p_double* l**3*u**3 - 35.0_p_double* l**4*u**3 + &
          7350.0_p_double* u**4 + 3745.0_p_double* l*u**4 + &
          630.0_p_double* l**2*u**4 + 35.0_p_double* l**3*u**4 - &
          882.0_p_double* u**5 - 273.0_p_double* l*u**5 - 21.0_p_double* l**2*u**5 + &
          49.0_p_double* u**6 + 7.0_p_double* l*u**6 - u**7)/5040.0_p_double
  case default
    laguerre_double = 0
  end select

end function laguerre_double
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Generalized Laguerre polynomials - single precision
!-----------------------------------------------------------------------------------------
function laguerre_single( p, l, u )
  
  implicit none

  real(p_single) :: laguerre_single
  integer, intent(in) :: p, l
  real(p_single), intent(in) :: u
  
  select case (p)
  case(0)
    laguerre_single = 1.0_p_single
  case(1)
    laguerre_single = 1.0_p_single + l - u
  case(2)
    laguerre_single = (2 + 3.0_p_single* l + l**2 - 4.0_p_single* u - &
          2.0_p_single* l*u + u**2)/2.0_p_single
  case(3)
    laguerre_single = (6.0_p_single + 11.0_p_single* l + 6.0_p_single* l**2 + l**3 - &
          18.0_p_single* u - 15.0_p_single* l*u - 3.0_p_single* l**2*u + &
          9.0_p_single* u**2 + 3.0_p_single* l*u**2 - u**3)/6.0_p_single
  case(4)
    laguerre_single =   (24.0_p_single + 50.0_p_single* l + 35.0_p_single* l**2 + &
          10.0_p_single* l**3 + l**4 - 96.0_p_single* u - 104.0_p_single* l*u - &
          36.0_p_single* l**2*u - 4.0_p_single* l**3*u + 72.0_p_single* u**2 + &
          42.0_p_single* l*u**2 + 6.0_p_single* l**2*u**2 - 16.0_p_single* u**3 - &
          4.0_p_single* l*u**3 + u**4)/24.0_p_single
  case(5)
    laguerre_single = (120.0_p_single + 274.0_p_single* l + 225.0_p_single* l**2 + &
          85.0_p_single* l**3 + 15.0_p_single* l**4 + l**5 - 600.0_p_single* u - &
          770.0_p_single* l*u - 355.0_p_single* l**2*u - 70.0_p_single* l**3*u - &
          5.0_p_single* l**4*u + 600.0_p_single* u**2 + 470.0_p_single* l*u**2 + &
          120.0_p_single* l**2*u**2 + 10.0_p_single* l**3*u**2 - &
          200.0_p_single* u**3 - 90.0_p_single* l*u**3 - 10.0_p_single* l**2*u**3 + &
          25.0_p_single* u**4 + 5.0_p_single* l*u**4 - u**5)/120.0_p_single
  case(6)
    laguerre_single = (720.0_p_single + 1764.0_p_single* l + 1624.0_p_single* l**2 + &
          735.0_p_single* l**3 + 175.0_p_single* l**4 + 21.0_p_single* l**5 + &
          l**6 - 4320.0_p_single* u - 6264.0_p_single* l*u - &
          3480.0_p_single* l**2*u - 930.0_p_single* l**3*u - &
          120.0_p_single* l**4*u - 6.0_p_single* l**5*u + 5400.0_p_single* u**2 + &
          5130.0_p_single* l*u**2 + 1785.0_p_single* l**2*u**2 + &
          270.0_p_single* l**3*u**2 + 15.0_p_single* l**4*u**2 - &
          2400.0_p_single* u**3 - 1480.0_p_single* l*u**3 - &
          300.0_p_single* l**2*u**3 - 20.0_p_single* l**3*u**3 + &
          450.0_p_single* u**4 + 165.0_p_single* l*u**4 + 15.0_p_single* l**2*u**4 - &
          36.0_p_single* u**5 - 6.0_p_single* l*u**5 + u**6)/720.0_p_single
  case(7)
    laguerre_single =  (5040.0_p_single + 13068.0_p_single* l + 13132.0_p_single* l**2 + &
          6769.0_p_single* l**3 + 1960.0_p_single* l**4 + 322.0_p_single* l**5 + &
          28.0_p_single* l**6 + l**7 - 35280.0_p_single* u - &
          56196.0_p_single* l*u - 35728.0_p_single* l**2*u - &
          11655.0_p_single* l**3*u - 2065.0_p_single* l**4*u - &
          189.0_p_single* l**5*u - 7.0_p_single* l**6*u + 52920.0_p_single* u**2 + &
          57834.0_p_single* l*u**2 + 24675.0_p_single* l**2*u**2 + &
          5145.0_p_single* l**3*u**2 + 525.0_p_single* l**4*u**2 + &
          21.0_p_single* l**5*u**2 - 29400.0_p_single* u**3 - &
          22330.0_p_single* l*u**3 - 6265.0_p_single* l**2*u**3 - &
          770.0_p_single* l**3*u**3 - 35.0_p_single* l**4*u**3 + &
          7350.0_p_single* u**4 + 3745.0_p_single* l*u**4 + &
          630.0_p_single* l**2*u**4 + 35.0_p_single* l**3*u**4 - &
          882.0_p_single* u**5 - 273.0_p_single* l*u**5 - 21.0_p_single* l**2*u**5 + &
          49.0_p_single* u**6 + 7.0_p_single* l*u**6 - u**7)/5040.0_p_single
  case default
    laguerre_single = 0
  end select

end function laguerre_single
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! This polynomial goes from 0 to 1 in the interval [0,1] and has 0 derivative at the ends  
!-----------------------------------------------------------------------------------------
function smooth_heaviside_double(x)  
  implicit none
  
  real(p_double) :: smooth_heaviside_double
  real(p_double), intent(in) :: x
  
  real(p_double) :: lx
  
  if ( x < -1.0_p_double ) then  
    smooth_heaviside_double = 0.0_p_double
  else if ( x > 1.0_p_double ) then
    smooth_heaviside_double = 1.0_p_double
  else
    lx = (x+1.0_p_double)*0.5_p_double
    smooth_heaviside_double = 10.0_p_double * lx**3 - 15.0_p_double * lx**4 + 6.0_p_double * lx**5
  endif

end function smooth_heaviside_double

function smooth_heaviside_single(x)  
  implicit none
  
  real(p_single) :: smooth_heaviside_single
  real(p_single), intent(in) :: x
  
  real(p_single) :: lx
  
  if ( x < -1.0_p_single ) then  
    smooth_heaviside_single = 0.0_p_single
  else if ( x > 1.0_p_single ) then
    smooth_heaviside_single = 1.0_p_single
  else
    lx = (x+1.0_p_single)*0.5_p_single
    smooth_heaviside_single = 10.0_p_single * lx**3 - 15.0_p_single * lx**4 + 6.0_p_single * lx**5
  endif

end function smooth_heaviside_single

!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Linear gradient
! @Return Real value between 0. and 1.
!-----------------------------------------------------------------------------------------
function gradient_linear( n, nmax )

  implicit none

  real( p_double ) :: gradient_linear, x
  integer, intent(in) :: n, nmax

  x = real( n ) / real( nmax )
  gradient_linear = x

end function gradient_linear
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Regressive (first strong, then weaker) gradient
! @Return Real value between 0. and 1.
!-----------------------------------------------------------------------------------------
function gradient_regressive( n, nmax )

  implicit none

  real( p_double ) :: gradient_regressive, x
  integer, intent(in) :: n, nmax

  real( p_double ), parameter :: st = 3.0 ! steepness
  real( p_double ), parameter :: scaling = SQRT( st**2 + 1.0 ) / st ! scale to [0., 1.]

  x = real( n ) / real( nmax )
  gradient_regressive = st * x / SQRT( ( st * x )**2 + 1.0 ) * scaling
  !gradient_regressive = 2.0 * x / ( x**2 + 1.0 )

end function gradient_regressive
!-----------------------------------------------------------------------------------------


end module m_math

