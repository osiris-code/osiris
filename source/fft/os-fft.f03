#include "os-config.h"
#include "os-preprocess.fpp"

! Only one-dimensional fft wrappers are currently implemented in this module
module m_fft

use, intrinsic :: iso_c_binding
use m_system

implicit none

! define data type of fft/ifft wrappers
#if defined(PRECISION_SINGLE)

integer, parameter :: type_real = C_FLOAT
integer, parameter :: type_complex = C_FLOAT_COMPLEX
#define _fftw_alloc_real(...) fftwf_alloc_real(__VA_ARGS__)
#define _fftw_alloc_complex(...) fftwf_alloc_complex(__VA_ARGS__)
#define _fftw_plan_dft_1d(...) fftwf_plan_dft_1d(__VA_ARGS__)
#define _fftw_plan_dft_r2c_1d(...) fftwf_plan_dft_r2c_1d(__VA_ARGS__)
#define _fftw_plan_dft_c2r_1d(...) fftwf_plan_dft_c2r_1d(__VA_ARGS__)
#define _fftw_execute_dft(...) fftwf_execute_dft(__VA_ARGS__)
#define _fftw_execute_dft_r2c(...) fftwf_execute_dft_r2c(__VA_ARGS__)
#define _fftw_execute_dft_c2r(...) fftwf_execute_dft_c2r(__VA_ARGS__)
#define _fftw_destroy_plan(...) fftwf_destroy_plan(__VA_ARGS__)
#define _fftw_free(...) fftwf_free(__VA_ARGS__)

#elif defined(PRECISION_DOUBLE)

integer, parameter :: type_real = C_DOUBLE
integer, parameter :: type_complex = C_DOUBLE_COMPLEX
#define _fftw_alloc_real(...) fftw_alloc_real(__VA_ARGS__)
#define _fftw_alloc_complex(...) fftw_alloc_complex(__VA_ARGS__)
#define _fftw_plan_dft_1d(...) fftw_plan_dft_1d(__VA_ARGS__)
#define _fftw_plan_dft_r2c_1d(...) fftw_plan_dft_r2c_1d(__VA_ARGS__)
#define _fftw_plan_dft_c2r_1d(...) fftw_plan_dft_c2r_1d(__VA_ARGS__)
#define _fftw_execute_dft(...) fftw_execute_dft(__VA_ARGS__)
#define _fftw_execute_dft_r2c(...) fftw_execute_dft_r2c(__VA_ARGS__)
#define _fftw_execute_dft_c2r(...) fftw_execute_dft_c2r(__VA_ARGS__)
#define _fftw_destroy_plan(...) fftw_destroy_plan(__VA_ARGS__)
#define _fftw_free(...) fftw_free(__VA_ARGS__)

#else

integer, parameter :: type_real = C_DOUBLE
integer, parameter :: type_complex = C_DOUBLE_COMPLEX
#define _fftw_alloc_real(...) fftw_alloc_real(__VA_ARGS__)
#define _fftw_alloc_complex(...) fftw_alloc_complex(__VA_ARGS__)
#define _fftw_plan_dft_1d(...) fftw_plan_dft_1d(__VA_ARGS__)
#define _fftw_plan_dft_r2c_1d(...) fftw_plan_dft_r2c_1d(__VA_ARGS__)
#define _fftw_plan_dft_c2r_1d(...) fftw_plan_dft_c2r_1d(__VA_ARGS__)
#define _fftw_execute_dft(...) fftw_execute_dft(__VA_ARGS__)
#define _fftw_execute_dft_r2c(...) fftw_execute_dft_r2c(__VA_ARGS__)
#define _fftw_execute_dft_c2r(...) fftw_execute_dft_c2r(__VA_ARGS__)
#define _fftw_destroy_plan(...) fftw_destroy_plan(__VA_ARGS__)
#define _fftw_free(...) fftw_free(__VA_ARGS__)

#endif


interface fft_init
  module procedure fftw_init_1d
  module procedure fftw_init2_1d
end interface

interface fft
  module procedure fftw_r2c_1d2
  module procedure fftw_r2c_1d3
end interface

interface ifft
  module procedure fftw_idft_1d0
  module procedure fftw_idft_1d1
  module procedure fftw_idft_1d2
  module procedure fftw_c2r_1d2
  module procedure fftw_c2r_1d3
end interface

interface fft_cleanup
  module procedure fftw_cleanup_1d
end interface


type(C_PTR) ::  plan_fft_x, plan_fft_y, plan_fft_z, &
                plan_ifft_x, plan_ifft_y, plan_ifft_z, &
                ptr_real_x, ptr_real_y, ptr_real_z, &
                ptr_complex_x, ptr_complex_y, ptr_complex_z

real(type_real), dimension(:), pointer :: fftw_real_x, fftw_real_y, fftw_real_z
complex(type_complex), dimension(:), pointer :: fftw_complex_x, fftw_complex_y, fftw_complex_z

private
public :: fft_init, fft, ifft, fft_cleanup, type_real, type_complex

save

contains


#ifdef FFTW_ENABLED

!-------------------------------------------------------------------------------
subroutine fftw_init_1d( n )
!-------------------------------------------------------------------------------
! Initialize FFT
!-------------------------------------------------------------------------------

  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  ! dummy variables
  integer, intent(in) :: n

  ! local variables
  integer :: fft_n

  fft_n = 2 * ( n/2 + 1 )

  ptr_real_x = _fftw_alloc_real( int( fft_n, C_SIZE_T ) )
  ptr_real_y = _fftw_alloc_real( int( fft_n, C_SIZE_T ) )
  ptr_real_z = _fftw_alloc_real( int( fft_n, C_SIZE_T ) )
  call c_f_pointer( ptr_real_x, fftw_real_x, [fft_n] )
  call c_f_pointer( ptr_real_y, fftw_real_y, [fft_n] )
  call c_f_pointer( ptr_real_z, fftw_real_z, [fft_n] )

  ptr_complex_x = _fftw_alloc_complex( int( fft_n, C_SIZE_T ) )
  ptr_complex_y = _fftw_alloc_complex( int( fft_n, C_SIZE_T ) )
  ptr_complex_z = _fftw_alloc_complex( int( fft_n, C_SIZE_T ) )
  call c_f_pointer( ptr_complex_x, fftw_complex_x, [fft_n] )
  call c_f_pointer( ptr_complex_y, fftw_complex_y, [fft_n] )
  call c_f_pointer( ptr_complex_z, fftw_complex_z, [fft_n] )


  plan_fft_x     = _fftw_plan_dft_r2c_1d( n, fftw_real_x, fftw_complex_x, FFTW_ESTIMATE )
  plan_fft_y     = _fftw_plan_dft_r2c_1d( n, fftw_real_y, fftw_complex_y, FFTW_ESTIMATE )
  plan_fft_z     = _fftw_plan_dft_r2c_1d( n, fftw_real_z, fftw_complex_z, FFTW_ESTIMATE )

  plan_ifft_x    = _fftw_plan_dft_c2r_1d( n, fftw_complex_x, fftw_real_x, FFTW_ESTIMATE )
  plan_ifft_y    = _fftw_plan_dft_c2r_1d( n, fftw_complex_y, fftw_real_y, FFTW_ESTIMATE )
  plan_ifft_z    = _fftw_plan_dft_c2r_1d( n, fftw_complex_z, fftw_real_z, FFTW_ESTIMATE )

end subroutine fftw_init_1d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_init2_1d( n, is_backward )
!-------------------------------------------------------------------------------
! Initialize FFT
!-------------------------------------------------------------------------------

  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  ! dummy variables
  integer, intent(in) :: n
  logical, intent(in) :: is_backward


  ptr_complex_x = _fftw_alloc_complex( int( n, C_SIZE_T ) )
  ptr_complex_y = _fftw_alloc_complex( int( n, C_SIZE_T ) )
  call c_f_pointer( ptr_complex_x, fftw_complex_x, [n] )
  call c_f_pointer( ptr_complex_y, fftw_complex_y, [n] )


  if ( is_backward .eqv. .true. ) then
    plan_ifft_x = _fftw_plan_dft_1d( n, fftw_complex_x, fftw_complex_y, FFTW_BACKWARD, FFTW_ESTIMATE )
  else
    plan_fft_x = _fftw_plan_dft_1d( n, fftw_complex_x, fftw_complex_y, FFTW_FORWARD, FFTW_ESTIMATE )
  endif

end subroutine fftw_init2_1d
!-------------------------------------------------------------------------------

subroutine fftw_idft_1d0( f, n )
!-------------------------------------------------------------------------------
! 1D IFFT (complex to complex) for 1D data
! this subroutine expect the frequency order to be [-n/2+1, ...0, ... n/2]
!-------------------------------------------------------------------------------


  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  real(p_double), intent(inout), dimension(:, :) :: f
  integer, intent(in) :: n

  ! local variables
  integer :: i
  !real(type_real) :: scale


  !scale = 1.0_type_real / real( n, type_real )


  fftw_complex_x = 0.0_type_complex

  ! duplicate to the buffer and do the fft shift
  do i = 1, n
    if ( i < (n+1) / 2 ) then
      fftw_complex_x(n/2+1+i) = cmplx(f(1, i), f(2, i))
    else
      fftw_complex_x(i-(n+1)/2+1) = cmplx(f(1, i), f(2, i))
    endif
  enddo

  call _fftw_execute_dft( plan_ifft_x, fftw_complex_x, fftw_complex_y )

  do i = 1, n
    f(1, i) = real(fftw_complex_y(i))!* scale
    f(2, i) = aimag(fftw_complex_y(i))! * scale
  enddo
    
end subroutine fftw_idft_1d0

subroutine fftw_idft_1d1( f, n )
!-------------------------------------------------------------------------------
! 2D IFFT (complex to complex) for 2D data on the last dimension
! this subroutine expect the frequency order to be [-nt/2+1, ...0, ... nt/2]
!-------------------------------------------------------------------------------


  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  real(p_double), intent(inout), dimension(:, :, :) :: f
  integer, dimension(2), intent(in) :: n

  ! local variables
  integer :: i, j, ind, nt, n1
  !real(type_real) :: scale


  !scale = 1.0_type_real / real( nt, type_real )

  nt = n(1)
  n1 = n(2)
  fftw_complex_x = 0.0_type_complex

  do j = 1, n1
    ! duplicate to the buffer and do the fft shift
    do i = 1, nt
      ind = mod( i+nt/2, nt ) + 1
      fftw_complex_x(ind) = cmplx(f(1, j, i), f(2, j, i))
    enddo

    call _fftw_execute_dft( plan_ifft_x, fftw_complex_x, fftw_complex_y )

    do i = 1, nt
      f(1, j, i) = real(fftw_complex_y(i)) !* scale
      f(2, j, i) = aimag(fftw_complex_y(i))  !* scale
    enddo
  enddo
   
end subroutine fftw_idft_1d1

subroutine fftw_idft_1d2( f, ns )
!-------------------------------------------------------------------------------
! 3D IFFT (complex to complex) for 3D data on the last dimension
! this subroutine expect the frequency order to be [-n/2+1, ...0, ... n/2]
!-------------------------------------------------------------------------------


  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  real(p_double), intent(inout), dimension(:, :, :, :) :: f
  integer, dimension(2), intent(in) :: ns

  ! local variables
  integer :: i, j, k, ind, n


  n = ns(1)


  fftw_complex_x = 0.0_type_complex

  do k = 1, ns(3)
    do j = 1, ns(2)

      ! duplicate to the buffer and do the fft shift
      do i = 1, n
        ind = mod( i+n/2, n ) + 1
        fftw_complex_x(ind) = cmplx(f(1, j, k, i), f(2, j, k, i))
      enddo

      call _fftw_execute_dft( plan_ifft_x, fftw_complex_x, fftw_complex_y )

      do i = 1, n
        f(1, j, k, i) = real(fftw_complex_y(i))
        f(2, j, k, i) = aimag(fftw_complex_y(i))
      enddo

    enddo
  enddo
    
end subroutine fftw_idft_1d2




!-------------------------------------------------------------------------------
subroutine fftw_r2c_1d2( f, n1, n2 )
!-------------------------------------------------------------------------------
! 1D FFT (real to complex) for 2D field
!-------------------------------------------------------------------------------


  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  real(type_real), intent(inout), dimension(:,:,:), pointer :: f
  integer, intent(in) :: n1, n2

  ! local variables
  integer :: i1, i2, half_n1
  real(type_real) :: scale

  half_n1 = n1/2 + 1

  scale = 1.0_type_real / real( n1, type_real )


  fftw_real_x = 0.0_type_real
  fftw_real_y = 0.0_type_real
  fftw_real_z = 0.0_type_real

  fftw_complex_x = 0.0_type_complex
  fftw_complex_y = 0.0_type_complex
  fftw_complex_z = 0.0_type_complex

  do i2 = 1, n2
    
    ! duplicate to the buffer
    fftw_real_x = f(1,:,i2)
    fftw_real_y = f(2,:,i2)
    fftw_real_z = f(3,:,i2)

    call _fftw_execute_dft_r2c( plan_fft_x, fftw_real_x, fftw_complex_x )
    call _fftw_execute_dft_r2c( plan_fft_y, fftw_real_y, fftw_complex_y )
    call _fftw_execute_dft_r2c( plan_fft_z, fftw_real_z, fftw_complex_z )

    do i1 = 1, half_n1

      f(1,2*i1-1,i2) = real( fftw_complex_x(i1) ) * scale
      f(2,2*i1-1,i2) = aimag( fftw_complex_x(i1) ) * scale
      f(3,2*i1-1,i2) = real( fftw_complex_y(i1) ) * scale
      f(1,2*i1  ,i2) = aimag( fftw_complex_y(i1) ) * scale
      f(2,2*i1  ,i2) = real( fftw_complex_z(i1) ) * scale
      f(3,2*i1  ,i2) = aimag( fftw_complex_z(i1) ) * scale

    enddo
  enddo
    
end subroutine fftw_r2c_1d2
!-------------------------------------------------------------------------------

subroutine fftw_c2r_1d2( f, n1, n2 )
!-------------------------------------------------------------------------------
! 1D inverse FFT (complex to real) for 2D field
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  real(type_real), intent(inout), dimension(:,:,:), pointer :: f
  integer, intent(in) :: n1, n2

  ! local variables
  integer :: i1, i2, half_n1
  real(type_real) :: scale

  scale = 1.0_type_real / real( n1, type_real )

  half_n1 = n1 / 2 + 1

  fftw_real_x = 0.0_type_real
  fftw_real_y = 0.0_type_real
  fftw_real_z = 0.0_type_real

  fftw_complex_x = 0.0_type_complex
  fftw_complex_y = 0.0_type_complex
  fftw_complex_z = 0.0_type_complex

  do i2 = 1, n2

    do i1 = 1, half_n1
      fftw_complex_x(i1) = cmplx( f(1,2*i1-1,i2), f(2,2*i1-1,i2) )
      fftw_complex_y(i1) = cmplx( f(3,2*i1-1,i2), f(1,2*i1  ,i2) )
      fftw_complex_z(i1) = cmplx( f(2,2*i1  ,i2), f(3,2*i1  ,i2) )
    enddo

    call _fftw_execute_dft_c2r( plan_ifft_x, fftw_complex_x, fftw_real_x )
    call _fftw_execute_dft_c2r( plan_ifft_y, fftw_complex_y, fftw_real_y )
    call _fftw_execute_dft_c2r( plan_ifft_z, fftw_complex_z, fftw_real_z )

    f(1,:,i2) = fftw_real_x
    f(2,:,i2) = fftw_real_y
    f(3,:,i2) = fftw_real_z

  enddo
    
end subroutine fftw_c2r_1d2
!-------------------------------------------------------------------------------

subroutine fftw_r2c_1d3( f, n1, n2, n3 )
!-------------------------------------------------------------------------------
! 1D FFT (real to complex) for 3D field
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  real(type_real), intent(inout), dimension(:,:,:,:), pointer :: f
  integer, intent(in) :: n1, n2, n3

  ! local variables
  integer :: i1, i2, i3, half_n1
  real(type_real) :: scale

  half_n1 = n1/2 + 1

  scale = 1.0_type_real / real( n1, type_real )


  fftw_real_x = 0.0_type_real
  fftw_real_y = 0.0_type_real
  fftw_real_z = 0.0_type_real

  fftw_complex_x = 0.0_type_complex
  fftw_complex_y = 0.0_type_complex
  fftw_complex_z = 0.0_type_complex

  do i3 = 1, n3
    do i2 = 1, n2
      
      ! duplicate to the buffer
      fftw_real_x = f(1,:,i2,i3)
      fftw_real_y = f(2,:,i2,i3)
      fftw_real_z = f(3,:,i2,i3)

      call _fftw_execute_dft_r2c( plan_fft_x, fftw_real_x, fftw_complex_x )
      call _fftw_execute_dft_r2c( plan_fft_y, fftw_real_y, fftw_complex_y )
      call _fftw_execute_dft_r2c( plan_fft_z, fftw_real_z, fftw_complex_z )

      do i1 = 1, half_n1

        f(1,2*i1-1,i2,i3) = real( fftw_complex_x(i1) ) * scale
        f(2,2*i1-1,i2,i3) = aimag( fftw_complex_x(i1) ) * scale
        f(3,2*i1-1,i2,i3) = real( fftw_complex_y(i1) ) * scale
        f(1,2*i1  ,i2,i3) = aimag( fftw_complex_y(i1) ) * scale
        f(2,2*i1  ,i2,i3) = real( fftw_complex_z(i1) ) * scale
        f(3,2*i1  ,i2,i3) = aimag( fftw_complex_z(i1) ) * scale

      enddo
    enddo
  enddo
    
end subroutine fftw_r2c_1d3
!-------------------------------------------------------------------------------

subroutine fftw_c2r_1d3( f, n1, n2, n3 )
!-------------------------------------------------------------------------------
! 1D inverse FFT (complex to real) for 3D field
!-------------------------------------------------------------------------------

  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  real(type_real), intent(inout), dimension(:,:,:,:), pointer :: f
  integer, intent(in) :: n1, n2, n3

  ! local variables
  integer :: i1, i2, i3, half_n1
  real(type_real) :: scale

  scale = 1.0_type_real / real( n1, type_real )

  half_n1 = n1 / 2 + 1

  fftw_real_x = 0.0_type_real
  fftw_real_y = 0.0_type_real
  fftw_real_z = 0.0_type_real

  fftw_complex_x = 0.0_type_complex
  fftw_complex_y = 0.0_type_complex
  fftw_complex_z = 0.0_type_complex

  do i3 = 1, n3
    do i2 = 1, n2

      do i1 = 1, half_n1
        fftw_complex_x(i1) = cmplx( f(1,2*i1-1,i2,i3), f(2,2*i1-1,i2,i3) )
        fftw_complex_y(i1) = cmplx( f(3,2*i1-1,i2,i3), f(1,2*i1  ,i2,i3) )
        fftw_complex_z(i1) = cmplx( f(2,2*i1  ,i2,i3), f(3,2*i1  ,i2,i3) )
      enddo

      call _fftw_execute_dft_c2r( plan_ifft_x, fftw_complex_x, fftw_real_x )
      call _fftw_execute_dft_c2r( plan_ifft_y, fftw_complex_y, fftw_real_y )
      call _fftw_execute_dft_c2r( plan_ifft_z, fftw_complex_z, fftw_real_z )

      f(1,:,i2,i3) = fftw_real_x
      f(2,:,i2,i3) = fftw_real_y
      f(3,:,i2,i3) = fftw_real_z

    enddo
  enddo
    
end subroutine fftw_c2r_1d3
!-------------------------------------------------------------------------------

subroutine fftw_cleanup_1d()

  use, intrinsic :: iso_c_binding

  implicit none

  include 'fftw3.f03'

  if (c_associated(plan_fft_x)) call _fftw_destroy_plan( plan_fft_x )
  if (c_associated(plan_fft_y)) call _fftw_destroy_plan( plan_fft_y )
  if (c_associated(plan_fft_z)) call _fftw_destroy_plan( plan_fft_z )
  if (c_associated(plan_ifft_x)) call _fftw_destroy_plan( plan_ifft_x )
  if (c_associated(plan_ifft_y)) call _fftw_destroy_plan( plan_ifft_y )
  if (c_associated(plan_ifft_z)) call _fftw_destroy_plan( plan_ifft_z )

  if (c_associated(ptr_real_x)) call _fftw_free( ptr_real_x )
  if (c_associated(ptr_real_y)) call _fftw_free( ptr_real_y )
  if (c_associated(ptr_real_z)) call _fftw_free( ptr_real_z )
  if (c_associated(ptr_complex_x)) call _fftw_free( ptr_complex_x )
  if (c_associated(ptr_complex_y)) call _fftw_free( ptr_complex_y )
  if (c_associated(ptr_complex_z)) call _fftw_free( ptr_complex_z )

end subroutine fftw_cleanup_1d


#else

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! NOTE: The following 6 routines are a stopgap fix to allow this code to compile 
! on systems without FFTW. It will be removed once a permanent FFT 
! library in introduced.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine err_out( name )
  implicit none
  character(*) :: name

  write(0,*) "An FFT routine was called but the Osiris was not compiled with FFT "
  write(0,*) "support via the FFTW library. " 
  write(0,*) "Please see the FFTW_DIR parameter in your system's configuration file. "
  write(0,*) "The routine called was: ", trim(name)
  call abort_program( ) 

end subroutine
!-------------------------------------------------------------------------------
subroutine fftw_init_1d( n )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none
  
  ! dummy variables
  integer, intent(in) :: n

  call err_out("fftw_init_1d")
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_init2_1d( n, is_backward )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  ! dummy variables
  integer, intent(in) :: n
  logical, intent(in) :: is_backward

  call err_out("fftw_init2_1d")
end subroutine fftw_init2_1d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_idft_1d0( f, n )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  real(p_double), intent(inout), dimension(:, :) :: f
  integer, intent(in) :: n
    
  call err_out("fftw_idft_1d0")
end subroutine fftw_idft_1d0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_idft_1d1( f, n )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  real(p_double), intent(inout), dimension(:, :, :) :: f
  integer, dimension(2), intent(in) :: n
    
  call err_out("fftw_idft_1d1")
end subroutine fftw_idft_1d1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_idft_1d2( f, n )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  real(p_double), intent(inout), dimension(:, :, :, :) :: f
  integer, dimension(2), intent(in) :: n
    
  call err_out("fftw_idft_1d2")
end subroutine fftw_idft_1d2
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_r2c_1d2( f, n1, n2 )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none
  
  ! dummy variables
  real(type_real), intent(inout), dimension(:,:,:), pointer :: f
  integer, intent(in) :: n1, n2

  call err_out("fftw_r2c_1d2")
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_c2r_1d2( f, n1, n2 )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none
  
  ! dummy variables
  real(type_real), intent(inout), dimension(:,:,:), pointer :: f
  integer, intent(in) :: n1, n2

  call err_out("fftw_c2r_1d2")
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_r2c_1d3( f, n1, n2, n3 )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none
  
  ! dummy variables
  real(type_real), intent(inout), dimension(:,:,:, :), pointer :: f
  integer, intent(in) :: n1, n2, n3  

  call err_out("fftw_r2c_1d3")
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_c2r_1d3( f, n1, n2, n3 )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none
  
  ! dummy variables
  real(type_real), intent(inout), dimension(:,:,:, :), pointer :: f 
  integer, intent(in) :: n1, n2, n3

  call err_out("fftw_c2r_1d3")
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_cleanup_1d( )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  call err_out("fftw_cleanup_1d")  
end subroutine
!-------------------------------------------------------------------------------

#endif

end module m_fft