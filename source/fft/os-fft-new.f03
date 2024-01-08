#include "os-config.h"
#include "os-preprocess.fpp"

! Only one-dimensional fft wrappers are currently implemented in this module
module m_fft_new

use, intrinsic :: iso_c_binding
use m_system

! use m_math, only : pi 

implicit none

real(p_double), parameter, private :: pi      = 3.14159265358979323846264_p_double    


! define data type of fft/ifft wrappers
#if defined(PRECISION_SINGLE)

integer, parameter :: type_real = C_FLOAT
integer, parameter :: type_complex = C_FLOAT_COMPLEX
#define _fftw_alloc_real(...) fftwf_alloc_real(__VA_ARGS__)
#define _fftw_alloc_complex(...) fftwf_alloc_complex(__VA_ARGS__)

#define _fftw_plan_dft_1d(...) fftwf_plan_dft_1d(__VA_ARGS__)
#define _fftw_plan_dft_r2c_1d(...) fftwf_plan_dft_r2c_1d(__VA_ARGS__)
#define _fftw_plan_dft_c2r_1d(...) fftwf_plan_dft_c2r_1d(__VA_ARGS__)

#define _fftw_plan_dft_2d(...) fftwf_plan_dft_2d(__VA_ARGS__)
#define _fftw_plan_dft_r2c_2d(...) fftwf_plan_dft_r2c_2d(__VA_ARGS__)
#define _fftw_plan_dft_c2r_2d(...) fftwf_plan_dft_c2r_2d(__VA_ARGS__)

#define _fftw_execute_dft(...) fftwf_execute_dft(__VA_ARGS__)
#define _fftw_execute_dft_r2c(...) fftwf_execute_dft_r2c(__VA_ARGS__)
#define _fftw_execute_dft_c2r(...) fftwf_execute_dft_c2r(__VA_ARGS__)
#define _fftw_destroy_plan(...) fftwf_destroy_plan(__VA_ARGS__)
#define _fftw_free(...) fftwf_free(__VA_ARGS__)

#else

integer, parameter :: type_real = C_DOUBLE
integer, parameter :: type_complex = C_DOUBLE_COMPLEX
#define _fftw_alloc_real(...) fftw_alloc_real(__VA_ARGS__)
#define _fftw_alloc_complex(...) fftw_alloc_complex(__VA_ARGS__)

#define _fftw_plan_dft_1d(...) fftw_plan_dft_1d(__VA_ARGS__)
#define _fftw_plan_dft_r2c_1d(...) fftw_plan_dft_r2c_1d(__VA_ARGS__)
#define _fftw_plan_dft_c2r_1d(...) fftw_plan_dft_c2r_1d(__VA_ARGS__)

#define _fftw_plan_dft_2d(...) fftw_plan_dft_2d(__VA_ARGS__)
#define _fftw_plan_dft_r2c_2d(...) fftw_plan_dft_r2c_2d(__VA_ARGS__)
#define _fftw_plan_dft_c2r_2d(...) fftw_plan_dft_c2r_2d(__VA_ARGS__)

#define _fftw_execute_dft(...) fftw_execute_dft(__VA_ARGS__)
#define _fftw_execute_dft_r2c(...) fftw_execute_dft_r2c(__VA_ARGS__)
#define _fftw_execute_dft_c2r(...) fftw_execute_dft_c2r(__VA_ARGS__)
#define _fftw_destroy_plan(...) fftw_destroy_plan(__VA_ARGS__)
#define _fftw_free(...) fftw_free(__VA_ARGS__)

#endif

type :: t_array_ptr 
  real(p_double), dimension(:), allocatable :: arr
end type 


!-------------------------------------------------------------------------------
type :: t_fft_manager 
  
  integer :: rank
  integer, dimension(:), allocatable :: shape
  integer :: size

  type(C_PTR) :: rbuf_cptr, cbuf1_cptr, cbuf2_cptr
  type(C_PTR) :: plan_fft_r2c, plan_ifft_c2r, plan_fft_c2c, plan_ifft_c2c

  real(type_real), dimension(:), pointer    :: rbuf
  complex(type_complex), dimension(:), pointer :: cbuf1
  complex(type_complex), dimension(:), pointer :: cbuf2

  real(type_real), dimension(:,:), pointer    :: rbuf_2d
  complex(type_complex), dimension(:,:), pointer :: cbuf1_2d
  complex(type_complex), dimension(:,:), pointer :: cbuf2_2d

  type(t_array_ptr), dimension(:), allocatable :: fft_freqs 

contains 

  procedure :: init => init_fft_manager 
  procedure :: cleanup => cleanup_fft_manager 
  procedure :: check_fftw_support => check_fftw_support_fft_manager

  procedure :: fft_c2c => fft_c2c_fft_manager
  procedure :: ifft_c2c => ifft_c2c_fft_manager

  procedure :: fft_2d_c2c => fft_2d_c2c_fft_manager
  procedure :: ifft_2d_c2c => ifft_2d_c2c_fft_manager

  ! to do 
  ! procedure :: fft_r2c => fft_r2c_fft_manager 
  ! procedure :: ifft_c2r => fft_c2r_fft_manager 

  procedure :: compute_fft_freqs => compute_fft_freqs_fft_manager

end type t_fft_manager
!-------------------------------------------------------------------------------


public :: t_fft_manager

contains 

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! only compile fftw functions if we linked fftw. otherwise 
! compile stubs so that the code still compiles (but will not run.)
#ifdef FFTW_ENABLED 
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
subroutine init_fft_manager( this, shape_in, rank )

  use, intrinsic :: iso_c_binding
  implicit none 

  include 'fftw3.f03'

  class( t_fft_manager ), intent(inout) :: this 
  integer, dimension(:), intent(in) :: shape_in
  integer, intent(in) :: rank 

  integer :: i, size
 
  this%rank = rank 

  allocate( this%shape( rank ) )

  do i = 1, rank 

    if( mod( shape_in(i), 2 ) .ne. 0 ) then 
      SCR_ROOT( "ERROR: attempted to initialize FFT buffers with odd n.")
      stop 
    endif 

    this%shape(i) = shape_in(i)

  enddo

  size = 1 
  do i = 1, rank   
    size = size * this%shape(i) 
  enddo
  this%size = size 

  this%rbuf_cptr = _fftw_alloc_real( int( size, C_SIZE_T ) )
  this%cbuf1_cptr = _fftw_alloc_complex( int( size, C_SIZE_T ) )
  this%cbuf2_cptr = _fftw_alloc_complex( int( size, C_SIZE_T ) )

  select case( rank ) 

  case( 1 ) 
  
    call c_f_pointer( this%rbuf_cptr, this%rbuf, this%shape )
    call c_f_pointer( this%cbuf1_cptr, this%cbuf1, this%shape )
    call c_f_pointer( this%cbuf2_cptr, this%cbuf2, this%shape )
   
    this%plan_fft_r2c     = _fftw_plan_dft_r2c_1d( this%shape(1), this%rbuf, this%cbuf1, FFTW_ESTIMATE )
    this%plan_ifft_c2r    = _fftw_plan_dft_c2r_1d( this%shape(1), this%cbuf1, this%rbuf, FFTW_ESTIMATE )
    this%plan_fft_c2c     = _fftw_plan_dft_1d( this%shape(1), this%cbuf1, this%cbuf2, FFTW_FORWARD, FFTW_ESTIMATE )
    this%plan_ifft_c2c    = _fftw_plan_dft_1d( this%shape(1), this%cbuf1, this%cbuf2, FFTW_BACKWARD, FFTW_ESTIMATE )

  case( 2 ) 
  
    call c_f_pointer( this%rbuf_cptr, this%rbuf_2d, this%shape )
    call c_f_pointer( this%cbuf1_cptr, this%cbuf1_2d, this%shape )
    call c_f_pointer( this%cbuf2_cptr, this%cbuf2_2d, this%shape )

    this%plan_fft_r2c     = _fftw_plan_dft_r2c_2d( this%shape(1), this%shape(2), this%rbuf_2d, this%cbuf1_2d, FFTW_ESTIMATE )
    this%plan_ifft_c2r    = _fftw_plan_dft_c2r_2d( this%shape(1), this%shape(2), this%cbuf1_2d, this%rbuf_2d, FFTW_ESTIMATE )
    this%plan_fft_c2c     = _fftw_plan_dft_2d( this%shape(1), this%shape(2), this%cbuf1_2d, this%cbuf2_2d, FFTW_FORWARD, FFTW_ESTIMATE )
    this%plan_ifft_c2c    = _fftw_plan_dft_2d( this%shape(1), this%shape(2), this%cbuf1_2d, this%cbuf2_2d, FFTW_BACKWARD, FFTW_ESTIMATE )

  case default 
    ! this%plan_fft_r2c     = _fftw_plan_dft_r2c( rank, shape, this%rbuf, this%cbuf1 FFTW_ESTIMATE )
    ! this%plan_ifft_c2r    = _fftw_plan_dft_c2r( rank, shape, this%cbuf1 this%rbuf, FFTW_ESTIMATE )
    ! this%plan_fft_c2c     = _fftw_plan_dft( rank, shape, this%cbuf1 this%cbuf2, FFTW_FORWARD, FFTW_ESTIMATE )
    ! this%plan_ifft_c2c    = _fftw_plan_dft( rank, shape, this%cbuf1 this%cbuf2, FFTW_BACKWARD, FFTW_ESTIMATE )

    SCR_ROOT( "ERROR: requsested FFT rank that has not been implemented.")
    SCR_ROOT( "FFT rank requested: ", rank )
    stop 

  end select 

end subroutine init_fft_manager
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
subroutine cleanup_fft_manager( this ) 
  use, intrinsic :: iso_c_binding

  implicit none

  include 'fftw3.f03'

  class( t_fft_manager ), intent(inout) :: this 

  integer :: i 

  if (c_associated(this%plan_fft_r2c)) call _fftw_destroy_plan( this%plan_fft_r2c )
  if (c_associated(this%plan_ifft_c2r)) call _fftw_destroy_plan( this%plan_ifft_c2r )
  if (c_associated(this%plan_fft_c2c)) call _fftw_destroy_plan( this%plan_fft_c2c )
  if (c_associated(this%plan_ifft_c2c)) call _fftw_destroy_plan( this%plan_ifft_c2c )

  call _fftw_free( this%rbuf_cptr ) 
  call _fftw_free( this%cbuf1_cptr )
  call _fftw_free( this%cbuf2_cptr )

  ! warning: if these are uncommented it shows that c_associated thinks they are
  ! still associated. this will not be a problem if the class is used properly 
  ! (i.e. not using it after calling cleanup), but could be a problem for 
  ! other parts of the code relying on c_associated to function properly. 
  ! print *, c_associated( this%rbuf_cptr )
  ! print *, c_associated( this%cbuf1cptr )
  ! print *, c_associated( this%cbuf2_cptr )

  if( allocated( this%fft_freqs )) then
    
    do i = 1, this%rank 
      deallocate( this%fft_freqs(i)%arr )
    enddo

    deallocate( this%fft_freqs )

  endif 

  deallocate( this%shape )

end subroutine cleanup_fft_manager
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
! do nothing since fftw is included. 
subroutine check_fftw_support_fft_manager( this ) 
  
  implicit none 
  class( t_fft_manager ), intent(in) :: this 

end subroutine check_fftw_support_fft_manager
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine compute_fft_freqs_fft_manager( this, spacing_arr ) 

  implicit none 
  class( t_fft_manager ), intent(inout) :: this 
  real(p_double), dimension(:), intent(in) :: spacing_arr

  real(p_double) :: dk
  integer :: i, dim, n

  if( .not. allocated( this%fft_freqs) ) then 
    allocate( this%fft_freqs( this%rank ) ) 
  endif 

  do dim = 1, this%rank 

    n = this%shape(dim)

    SCR_ROOT( dim, n ) 

    allocate( this%fft_freqs(dim)%arr(n) )
 
    dk = 2 * pi / (spacing_arr(dim) * n )

    do i = 1, n/2  
      this%fft_freqs(dim)%arr(i) = (i - 1) * dk
    enddo 

    do i = n/2 + 1, n 
      this%fft_freqs(dim)%arr(i) = (i - 1 - n) * dk 
    enddo

  enddo

end subroutine compute_fft_freqs_fft_manager
!-------------------------------------------------------------------------------



!--------------------- 1D Routines ---------------------------------------------

!-------------------------------------------------------------------------------
subroutine fft_c2c_fft_manager( this, x )

  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  class( t_fft_manager ), intent(inout) :: this 
  complex(p_double), intent(inout), dimension(:) :: x

  integer :: i

  ! duplicate to the buffer
  do i = 1, this%shape(1)
    this%cbuf1(i) = x(i)
  enddo

  call _fftw_execute_dft( this%plan_fft_c2c, this%cbuf1, this%cbuf2 )

  do i = 1, this%shape(1)
    x(i) = this%cbuf2(i) 
  enddo 

end subroutine fft_c2c_fft_manager
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
subroutine ifft_c2c_fft_manager( this, x )

  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  class( t_fft_manager ), intent(inout) :: this 
  complex(p_double), intent(inout), dimension(:) :: x

  real(p_double) :: scale
  integer :: i 

  scale = 1.0_p_double / this%size

  ! duplicate to the buffer 
  do i = 1, this%shape(1)
    this%cbuf1(i) = x(i)
  enddo

  call _fftw_execute_dft( this%plan_ifft_c2c, this%cbuf1, this%cbuf2 )

  do i = 1, this%shape(1) 
    x(i) = this%cbuf2(i) * scale 
  enddo

end subroutine ifft_c2c_fft_manager
!-------------------------------------------------------------------------------



!--------------------- 2D Routines ---------------------------------------------

!-------------------------------------------------------------------------------
subroutine fft_2d_c2c_fft_manager( this, x )

  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  class( t_fft_manager ), intent(inout) :: this 
  complex(p_double), intent(inout), dimension(:,:) :: x

  integer :: i,j 

  ! duplicate to the buffer

  do i = 1, this%shape(1)
    do j = 1, this%shape(2)
      this%cbuf1_2d(i,j) = x(i,j)
    enddo
  enddo

  call _fftw_execute_dft( this%plan_fft_c2c, this%cbuf1_2d, this%cbuf2_2d )

  do i = 1, this%shape(1)
    do j = 1, this%shape(2)
        x(i,j) = this%cbuf2_2d(i,j) 
    enddo
  enddo 

end subroutine fft_2d_c2c_fft_manager
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
subroutine ifft_2d_c2c_fft_manager( this, x )

  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  class( t_fft_manager ), intent(inout) :: this 
  complex(p_double), intent(inout), dimension(:,:) :: x

  real(p_double) :: scale
  integer :: i, j 

  scale = 1.0_p_double / this%size
  
  ! duplicate to the buffer 

  do i = 1, this%shape(1)
    do j = 1, this%shape(2)
      this%cbuf1_2d(i,j) = x(i,j)
    enddo
  enddo
  
  call _fftw_execute_dft( this%plan_ifft_c2c, this%cbuf1_2d, this%cbuf2_2d )

  do i = 1, this%shape(1)
    do j = 1, this%shape(2)
      x(i,j) = this%cbuf2_2d(i,j) * scale 
    enddo
  enddo

end subroutine ifft_2d_c2c_fft_manager
!-------------------------------------------------------------------------------







!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! stubs so that the code compiles without fft support 
#else 
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine init_fft_manager( this,shape_in, rank )

  implicit none 

  class( t_fft_manager ), intent(inout) :: this 
  integer, dimension(:), intent(in) :: shape_in
  integer, intent(in) :: rank 

  call this % check_fftw_support() 

end subroutine init_fft_manager
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine cleanup_fft_manager( this ) 
 
  implicit none
  class( t_fft_manager ), intent(inout) :: this 
  
  call this % check_fftw_support() 

end subroutine cleanup_fft_manager
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine compute_fft_freqs_fft_manager( this, spacing ) 

  implicit none 
  class( t_fft_manager ), intent(inout) :: this 
  real(p_double), dimension(:), intent(in) :: spacing 

  call this % check_fftw_support() 

end subroutine compute_fft_freqs_fft_manager
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine check_fftw_support_fft_manager( this ) 

  implicit none 
  class( t_fft_manager ), intent(in) :: this 

  SCR_ROOT( "ERROR: attempted to call t_fft_manager manager subroutine, ")
  SCR_ROOT( "but the code was not compiled with FFTW.")
  
  stop 

end subroutine check_fftw_support_fft_manager
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine fft_c2c_fft_manager( this, x )

  implicit none 
  class( t_fft_manager ), intent(inout) :: this 
  complex(p_double), intent(inout), dimension(:) :: x

  call this % check_fftw_support() 

end subroutine fft_c2c_fft_manager
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine ifft_c2c_fft_manager( this, x )

  implicit none
  class( t_fft_manager ), intent(inout) :: this 
  complex(p_double), intent(inout), dimension(:) :: x

  call this % check_fftw_support() 

end subroutine ifft_c2c_fft_manager
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine fft_2d_c2c_fft_manager( this, x )

  implicit none 
  class( t_fft_manager ), intent(inout) :: this 
  complex(p_double), intent(inout), dimension(:,:) :: x

  call this % check_fftw_support() 

end subroutine fft_2d_c2c_fft_manager
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine ifft_2d_c2c_fft_manager( this, x )

  implicit none
  class( t_fft_manager ), intent(inout) :: this 
  complex(p_double), intent(inout), dimension(:,:) :: x

  call this % check_fftw_support() 

end subroutine ifft_2d_c2c_fft_manager
!-------------------------------------------------------------------------------




!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! end check for FFTW_ENABLED 
#endif  
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


end module m_fft_new