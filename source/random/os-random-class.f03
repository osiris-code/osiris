!-------------------------------------------------------------------------------
! Random number generator module
!
!
! The folowing routines are available:
!
! init_genrand(seed) - Initialize the random number generator. seed is an 
!                      integer scalar
!
! genrand_int32()    - return a random integer in the range 
!                      [-2147483648, 2147483647] (all possible 32 bit integers)
!
! genrand_real1()    - return a random double precision real in the range
!                      [0,1] with 32 bit resolution
!
! genrand_real2()    - return a random double precision real in the range
!                      [0,1) with 32 bit resolution
!
! genrand_real3()    - return a random double precision real in the range
!                      (0,1) with 32 bit resolution
!
! genrand_res53()    - return a random double precision real in the range
!                      [0,1) with 53 bit resolution. (slower, requires two 
!                      random integers)
!
! genrand_gaussian() - return a random number with from a normal distribution
!                      with 0 mean and unit variance.
!
! genrand_half_max() - return a random number with from a half maxwellian distribution
!                      f(x) = x * exp( -x^2 / 2 ) 
!
!-------------------------------------------------------------------------------

! Enable this flag if, while testing, want to keep track of the number of random numbers
!  generated (when using the t_random_grid)
#define INC_COUNT call this%inc_count(ix) 


#include "os-config.h"
#include "os-preprocess.fpp"

module m_random_class

use m_parameters
use m_restart

use m_grid_define
use m_node_conf

private

integer, parameter :: p_maxlen_random_class_name = 64

type, abstract :: t_random

  ! Gaussian distribution variables
  real(p_double) :: gaussian_set = 0.0d0 
  logical        :: new_gaussian = .true.
  logical        :: enforce_random_constancy = .false.

contains

  ! The actual random number generator
  ! must be overriden by subclasses
  procedure(genrand_int32), deferred :: genrand_int32

  procedure(init_genrand_scalar), deferred :: init_genrand_scalar

  procedure :: init_genrand_scalar_grid

  ! The following is not used so we'll remove it
  !procedure :: init_genrand_rst
  
  generic   :: init => init_genrand_scalar,init_genrand_scalar_grid !, init_genrand_rst

  procedure(write_checkpoint), deferred :: write_checkpoint
  procedure(read_checkpoint), deferred :: read_checkpoint

  procedure(class_name), deferred :: class_name

  procedure :: genrand_real1
  procedure :: genrand_real2
  procedure :: genrand_real3
  procedure :: genrand_res53

  procedure :: harvest_real1_p4
  procedure :: harvest_real1_p8
  generic   :: harvest_real1 => harvest_real1_p4, harvest_real1_p8

  procedure :: harvest_real2_p4
  procedure :: harvest_real2_p8  
  generic   :: harvest_real2 => harvest_real2_p4, harvest_real2_p8

  procedure :: harvest_real3_p4
  procedure :: harvest_real3_p8 
  generic   :: harvest_real3 => harvest_real3_p4, harvest_real3_p8

  procedure :: genrand_gaussian
  procedure :: genrand_half_max
  procedure :: genrand_spherical

  procedure :: harvest_spherical_p4
  procedure :: harvest_spherical_p8
  generic   :: harvest_spherical => harvest_spherical_p4, harvest_spherical_p8

  procedure :: harvest_relmax_p4
  procedure :: harvest_relmax_p8
  generic   :: harvest_relmax => harvest_relmax_p4, harvest_relmax_p8

  procedure :: harvest_relmax_boost_p4
  procedure :: harvest_relmax_boost_p8
  generic   :: harvest_relmax_boost => harvest_relmax_boost_p4, &
                                       harvest_relmax_boost_p8
  
end type t_random

abstract interface
  function genrand_int32( this, ix)
    import :: t_random
    integer, dimension(:), optional ::  ix
    class(t_random), intent(inout) :: this
    integer :: genrand_int32
  end function genrand_int32
end interface

abstract interface
  subroutine init_genrand_scalar( this, s )
    import :: t_random
    class(t_random), intent(inout) :: this
    integer, intent(in)            :: s
  end subroutine init_genrand_scalar
end interface

abstract interface
  subroutine write_checkpoint( this, restart_handle )

    import :: t_random, t_restart_handle

    class( t_random ), intent(in) :: this
    class( t_restart_handle ), intent(inout) :: restart_handle

  end subroutine write_checkpoint
end interface

abstract interface
  subroutine read_checkpoint( this, restart_handle )

    import :: t_random, t_restart_handle

    class( t_random ), intent(inout) :: this
    class( t_restart_handle ), intent(in) :: restart_handle

  end subroutine read_checkpoint
end interface

abstract interface
  function class_name( this )
    import :: t_random, p_maxlen_random_class_name
    class( t_random ), intent(in) :: this
    character( len = p_maxlen_random_class_name ) :: class_name
  end function class_name
end interface


public :: t_random, p_maxlen_random_class_name

contains

!-------------------------------------------------------------------------------
! generates a random number on [0,1]-real-interval (32 bit resolution)
!-------------------------------------------------------------------------------
function genrand_real1( this, ix )

  implicit none
  
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional  ::  ix
  real(p_double) :: genrand_real1
  genrand_real1 = ( this%genrand_int32(ix) + 2147483648.0d0 ) * &
                               (1.0d0/4294967295.0d0) ! divided by 2**32 - 1  

end function genrand_real1

!-------------------------------------------------------------------------------
! generates a random number on [0,1]-real-interval (32 bit resolution)
!-------------------------------------------------------------------------------
subroutine harvest_real1_p8( this, rnd, ix )

  implicit none
  
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_double), intent(out) :: rnd
  rnd = ( this%genrand_int32(ix) + 2147483648.0d0 ) * &
                               (1.0d0/4294967295.0d0) ! divided by 2**32 - 1  

end subroutine harvest_real1_p8

!-------------------------------------------------------------------------------
! generates a random number on [0,1]-real-interval (24 bit resolution)
!-------------------------------------------------------------------------------
subroutine harvest_real1_p4( this, rnd, ix )

  implicit none
  
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_single), intent(out) :: rnd

  ! throw away excess precision
  rnd = real(ishft( this%genrand_int32(ix), -8), p_single) * &
         (1.0_p_single/16777215.0_p_single)

end subroutine harvest_real1_p4

!-------------------------------------------------------------------------------
! generates a random number on [0,1)-real-interval (32 bit resolution) 
!-------------------------------------------------------------------------------
function genrand_real2( this, ix )

  implicit none
  
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_double) :: genrand_real2
  
  genrand_real2 = ( this%genrand_int32(ix) + 2147483648.0d0 ) * &
                                (1.0d0/4294967296.0d0) ! divided by 2**32
  
end function genrand_real2
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! generates a random number on [0,1)-real-interval (32 bit resolution)
!-------------------------------------------------------------------------------
subroutine harvest_real2_p8( this, rnd, ix)

  implicit none
  
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_double), intent(out) :: rnd
  rnd = ( this%genrand_int32(ix) + 2147483648.0d0 ) * &
                                (1.0d0/4294967296.0d0) ! divided by 2**32
                              
end subroutine harvest_real2_p8

!-------------------------------------------------------------------------------
! generates a random number on [0,1)-real-interval (24 bit resolution)
!-------------------------------------------------------------------------------
subroutine harvest_real2_p4( this, rnd, ix )

  implicit none
  
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_single), intent(out) :: rnd

  rnd = real(ishft(this%genrand_int32(ix), -8), p_single) * &
         (1.0_p_single/16777216.0_p_single)

end subroutine harvest_real2_p4



!-------------------------------------------------------------------------------
! generates a random number on (0,1)-real-interval (32 bit resolution)
!-------------------------------------------------------------------------------
function genrand_real3( this, ix )

  implicit none
  
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_double) :: genrand_real3
  
  genrand_real3 = ( this%genrand_int32(ix) + 2147483649.0d0) * &
                                (1.0d0/4294967297.0d0 ) 
                                ! divided by 2**32 + 1 

end function genrand_real3

!-------------------------------------------------------------------------------
! generates a random number on (0,1)-real-interval (32 bit resolution)
!-------------------------------------------------------------------------------
subroutine harvest_real3_p8( this, rnd, ix )

  implicit none
  
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_double), intent(out) :: rnd
  rnd = ( this%genrand_int32(ix) + 2147483649.0d0) * &
                                (1.0d0/4294967297.0d0 ) 

end subroutine harvest_real3_p8

!-------------------------------------------------------------------------------
! generates a random number on (0,1)-real-interval (24 bit resolution)
!-------------------------------------------------------------------------------
subroutine harvest_real3_p4( this, rnd, ix )

  implicit none
  
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_single), intent(out) :: rnd
  
  rnd = ( real(ishft( this%genrand_int32(ix), -8), p_single) + 1.0) * &
         (1.0_p_single/16777218.0_p_single)

end subroutine harvest_real3_p4

!-------------------------------------------------------------------------------
! generates a random number on [0,1)-real-interval with 53 bit resolution
! (full double precision resolution)
! please note that this routine will give 0 when the two integer
! random numbers are  0 (0x00000000), and 1-eps when the two integer
! random numbers are -1 (0xffffffff), unlike the previous ones
!-------------------------------------------------------------------------------
function genrand_res53( this, ix ) 

  implicit none
  
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_double) :: genrand_res53
  integer :: a, b
  
  a = ishft(this%genrand_int32(ix), -5) ! convert a to 27 bit integer (discard LSBs)
  b = ishft(this%genrand_int32(ix), -6) ! convert a to 26 bit integer (discard LSBs)
  
  genrand_res53 = (a*67108864.0d0+b)*(1.0d0/9007199254740992.0d0)              
                 ! mult by 2**26            divided by 2**53

end function genrand_res53

!-------------------------------------------------------------------------------
! generates a random number with a normally distributed deviate with 0 mean
! and unit variance. This routine follows the method found in section 7.2
! (pg. 289) of Numerical Recipes in C, 2nd edition.
!
! This routine uses genrand_real3 that always returns numbers in the (0,1)
! interval.
!
! The module variables iset and gset need to saved at restart.
!-------------------------------------------------------------------------------
function genrand_gaussian( this, ix ) 

  implicit none
  
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix

  real(p_double) :: genrand_gaussian
  
  real(p_double) :: v1, v2, rsq, fac

  if ( this%new_gaussian ) then
    
    do
      ! pick 2 uniform numbers inside the square extending from -1 to 1
      ! in each direction
      v1 = ( this%genrand_int32(ix) + 0.5_p_double ) * &
             (1.0_p_double/2147483648.0_p_double ) 
      v2 = ( this%genrand_int32(ix) + 0.5_p_double ) * &
             (1.0_p_double/2147483648.0_p_double ) 
      
      ! check if they are inside the unit circle, and not (0,0) 
      ! otherwise try again
      rsq = v1**2 + v2**2
      if (( rsq < 1.0_p_double ) .and. ( rsq > 0.0d0 )) exit
    end do
    
    ! Use Box-Muller method to generate random deviates with 
    ! normal (gaussian) distribution
    fac = sqrt(-2.0*log(rsq)/rsq)
    
    ! store 1 value for future use
    this%gaussian_set = v1*fac 
    this%new_gaussian = .false.

    genrand_gaussian = v2*fac
  else

    ! Use previously generated value
    genrand_gaussian  = this%gaussian_set
    this%new_gaussian = .true.
  endif
  
  
end function genrand_gaussian
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Generates a random number obeying the distribution f(x) = x * exp( -x^2 / 2 )
!
! The result is in the range [2.157918643883382E-005, 6.66043688929654] ( which
! correspond to the max and min value of the integer random number generator )
!-------------------------------------------------------------------------------
function genrand_half_max( this, ix )

  implicit none

  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix

  real(p_double) :: genrand_half_max
  
  real(p_double) :: v1

  ! generate a random number between 2.328306435996595E-010 and 0.999999999767169
  v1 = ( this%genrand_int32(ix) + 2147483649.0d0) * (1.0d0/4294967297.0d0 )
  
  ! Use transformation methodo to get the required distribution
  ! D[ InverseFunc[sqrt( -2*log(x))], x] = x * exp( -x^2 / 2)
  genrand_half_max = sqrt( -2 * log( v1 ))

end function genrand_half_max

!-------------------------------------------------------------------------------
! random number theta with f[theta] = Sin[theta], theta element [0,Pi]
! used for spherical angle to create isotropic distribution
!-------------------------------------------------------------------------------
function genrand_spherical( this, ix )

  implicit none
  
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_double) :: genrand_spherical
  
  genrand_spherical = 2.0_p_double * asin( sqrt(this % genrand_real1(ix)) )

end function genrand_spherical

subroutine harvest_spherical_p8( this, rnd, ix )

  implicit none
  
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_double), intent(out) :: rnd

  call this % harvest_real1_p8(rnd,ix)
  rnd = 2.0_p_double * asin( sqrt(rnd) )

end subroutine harvest_spherical_p8

subroutine harvest_spherical_p4( this, rnd, ix )

  implicit none
  
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_single), intent(out) :: rnd

  call this % harvest_real1_p4(rnd,ix)
  rnd = 2.0_p_single * asin( sqrt(rnd) )

end subroutine harvest_spherical_p4

!-------------------------------------------------------------------------------
! Harvests n random numbers distributed according to the radial component of 
! the relativistic maxwellian (Maxwell-Juettner): 
!
!    f(u) ~ u^2 exp(-sqrt(1+u^2)/T)
!
! u_max is the cutoff for the distribution
!-------------------------------------------------------------------------------

subroutine harvest_relmax_p8( this, T, u_max, n, f, ix )
  
  implicit none
    
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_double), intent(in) :: T
  real(p_double), intent(in) :: u_max  
  integer, intent(in) :: n
  real(p_double), dimension(:), intent(inout) :: f

  real(p_double) :: s0, s1, s2, s3, s4, norm
  real(p_double) :: test_x, test_y, rnd, cap, juettner
  integer :: i

  ! Constant initializations
  s0 = 2*(T**2 + sqrt(T**2 + T**4))
  s1 = sqrt(s0)
  s2 = 1 + 5*T
  norm = exp( sqrt(1 + s0)/T ) / s0
  
  s3 = atan( 2 * s1/s2 ) 
  s4 = atan( 2 * (s1 - u_max) / s2 )
  
  do i = 1, n
     do
       call this%harvest_real1_p8(rnd, ix)
    
       test_x = s1 + 0.5 * s2 * tan( (rnd-1) * s3 - rnd * s4 )
    
       cap = s2**2 / ( s2**2 + 4*( test_x - s1)**2 )

       juettner = norm * test_x**2 * exp( -sqrt( 1 + test_x**2 ) / T )

       ! roundoff protection - in some rare combination of parameters the value 
       ! of juettner is slightly above 1 due to roundoff
    
       if ( juettner > 1 ) juettner = 1

       call this%harvest_real1_p8(test_y, ix)
       test_y = test_y * cap
    
       if( test_y < juettner) exit

     enddo

     f(i) = test_x 
  enddo 

end subroutine harvest_relmax_p8

subroutine harvest_relmax_p4( this, T, u_max, n, f, ix )
    
  implicit none
  
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_single), intent(in) :: T
  real(p_single), intent(in) :: u_max  
  integer, intent(in) :: n
  real(p_single), dimension(:), intent(inout) :: f

  real(p_single) :: s0, s1, s2, s3, s4, norm
  real(p_single) :: test_x, test_y, rnd, cap, juettner
  integer :: i

  ! Constant initializations
  s0 = 2*(T**2 + sqrt(T**2 + T**4))
  s1 = sqrt(s0)
  s2 = 1 + 5*T
  norm = exp( sqrt(1 + s0)/T ) / s0
  
  s3 = atan( 2 * s1/s2 ) 
  s4 = atan( 2 * (s1 - u_max) / s2 )
  
  do i = 1, n
     do
       call this%harvest_real1_p4(rnd, ix)
    
       test_x = s1 + 0.5 * s2 * tan( (rnd-1) * s3 - rnd * s4 )
    
       cap = s2**2 / ( s2**2 + 4*( test_x - s1)**2 )

       juettner = norm * test_x**2 * exp( -sqrt( 1 + test_x**2 ) / T )

       ! roundoff protection - in some rare combination of parameters the value 
       ! of juettner is slightly above 1 due to roundoff
    
       if ( juettner > 1 ) juettner = 1

       call this%harvest_real1_p4(test_y, ix)
       test_y = test_y * cap
    
       if( test_y < juettner) exit

     enddo

     f(i) = test_x 
  enddo 

end subroutine harvest_relmax_p4


!-----------------------------------------------------------------------------------------
! Harvests n random numbers distributed according to the boosted Maxwell Juettner
! distribution
!-----------------------------------------------------------------------------------------

subroutine harvest_relmax_boost_p8( this, T, u_boost, n, upar_out, uperp_out, ix )
  
  implicit none
    
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_double), intent(in) :: T
  real(p_double), intent(in) :: u_boost

  integer, intent(in) :: n
  real(p_double), dimension(:), intent(inout) :: upar_out, uperp_out
  
  real( p_double ) :: gamma_boost, beta_boost, gamma_boost_T, aux, tampa, uperpmax
  real( p_double ) :: rnd1, rnd2, rnd3
  real( p_double ) :: uparmin, deltaupar, gamma, boosted_juettner, upar, uperp
  
  integer :: i
  
  ! Constant initializations
  gamma_boost = sqrt( u_boost**2 + 1 )
  beta_boost = u_boost / gamma_boost

  gamma_boost_T = gamma_boost/T
  aux = (T**2+T*sqrt(4+T**2)) / 2
  
  ! this corresponds to the peak value of the distribution function for a given T
  ! it does not depend on beta_boost
  tampa = exp(( - sqrt(1+aux))/T)*sqrt(aux) 

  ! this is a way of helping the algorithm when dealing with sub-relativistic temperatures
  if (T >= 0.1) then
     uperpmax = 20 * T
  else
     uperpmax = sqrt( 32 * T)
  endif
  
  uparmin   = gamma_boost*( -uperpmax + beta_boost*sqrt( 1 + uperpmax**2) )
  deltaupar = 2 * gamma_boost * uperpmax
  
  ! loop over all values
  do i = 1, n
  
    ! Generate upar and uperp using the rejection method
    do
  
      call this % harvest_real1_p8( rnd1, ix )
      upar = uparmin + rnd1 * deltaupar
      
      call this % harvest_real1_p8( rnd2, ix )
      uperp = uperpmax*rnd2
      
      ! gamma factor of test values
      gamma = sqrt( 1 + upar**2 + uperp**2 )
      
      boosted_juettner = uperp * exp( - gamma_boost_T* ( gamma - beta_boost*upar ) )

      call this % harvest_real1_p8( rnd3, ix )
      if ( rnd3*tampa < boosted_juettner) exit

    enddo

    upar_out(i)  = upar
    uperp_out(i) = uperp

  enddo

end subroutine harvest_relmax_boost_p8


subroutine harvest_relmax_boost_p4( this, T, u_boost, n, upar_out, uperp_out, ix )
  
  implicit none
    
  class( t_random ), intent(inout) :: this
  integer, dimension(:), optional ::  ix
  real(p_single), intent(in) :: T
  real(p_single), intent(in) :: u_boost

  integer, intent(in) :: n
  real(p_single), dimension(:), intent(inout) :: upar_out, uperp_out
  
  ! Intermediate calculations are always done in double precision
  real( p_double ) :: gamma_boost, beta_boost, gamma_boost_T, aux, tampa, uperpmax
  real( p_double ) :: rnd1, rnd2, rnd3
  real( p_double ) :: uparmin, deltaupar, gamma, boosted_juettner, upar, uperp
  
  integer :: i
  
  ! Constant initializations
  gamma_boost = sqrt( u_boost**2 + 1.0_p_double )
  beta_boost = u_boost / gamma_boost

  gamma_boost_T = gamma_boost/T
  aux = (T**2+T*sqrt(4+T**2)) / 2
  
  ! this corresponds to the peak value of the distribution function for a given T
  ! it does not depend on beta_boost
  tampa = exp(- sqrt(1+aux)/T ) * sqrt(aux) 
 
  ! this is a way of helping the algorithm when dealing with sub-relativistic temperatures
  if (T >= 0.1) then
     uperpmax = 20 * T
  else
     uperpmax = sqrt( 32 * T)
  endif
  
  uparmin   = gamma_boost*( -uperpmax + beta_boost*sqrt( 1 + uperpmax**2) )
  deltaupar = 2 * gamma_boost * uperpmax
  
  ! loop over all values
  do i = 1, n
  
    ! Generate upar and uperp using the rejection method
    do
  
      call this % harvest_real1_p8( rnd1, ix )
      upar = uparmin + rnd1 * deltaupar
      
      call this % harvest_real1_p8( rnd2, ix )
      uperp = uperpmax*rnd2
      
      ! gamma factor of test values
      gamma = sqrt( 1 + upar**2 + uperp**2 )
      
      boosted_juettner = uperp * exp( - gamma_boost_T* ( gamma - beta_boost*upar ) )

      call this % harvest_real1_p8( rnd3, ix )
      if ( rnd3*tampa < boosted_juettner) exit

    enddo

    upar_out(i)  = real( upar, p_single )
    uperp_out(i) = real( uperp, p_single )

  enddo

end subroutine harvest_relmax_boost_p4


!-------------------------------------------------------------------------------
! Initializes random number generator seed with a scalar integer or from a
! checkpoint information
!
!    Commented out because not currently used
!-------------------------------------------------------------------------------
! subroutine init_genrand_rst( this, s, restart, restart_handle)
 
!   implicit none

!   class( t_random ), intent(inout) :: this
  
!   integer, intent(in) :: s
!   logical, intent(in) :: restart
!   type( t_restart_handle ), intent(in) :: restart_handle
  
!   if ( restart ) then
!     call this % read_checkpoint( restart_handle )
!   else
!     call this % init( s )
!   endif
  
! end subroutine init_genrand_rst

subroutine init_genrand_scalar_grid( this, s, grid, class_idx )
  use m_grid
  use m_restart

  implicit none

  class(t_random), intent(inout)  :: this
  integer, intent(in)             :: s
  class( t_grid ), intent(in)      ::  grid
  integer, intent(inout) :: class_idx

  call this%init_genrand_scalar( s)

end subroutine init_genrand_scalar_grid

end module m_random_class

