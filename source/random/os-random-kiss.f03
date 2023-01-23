!-------------------------------------------------------------------------------
! KISS Pseudo Random Number Generator
! Marsaglia, George (2003) "Random Number Generators," Journal of Modern Applied
!  Statistical Methods: Vol. 2: Iss. 1, Article 2.
!
! Period > 2^124 ~ 10^37 
!
! Memory: 4 x 32 bit integers (16 bytes)
! Speed: ~ 520 M numbers/s
!-------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_random_kiss

use m_parameters
use m_restart
use m_random_class  

private  

! string to id restart data
character(len=*), parameter :: p_rand_rst_id = "random_kiss rst data - 0x0000"

type, extends( t_random ) :: t_random_kiss

  ! state data
  integer :: x = 123456789
  integer :: y = 362436000
  integer :: z = 521288629
  integer :: c = 7654321

  contains

  procedure :: genrand_int32         => kiss_int32
  procedure :: init_genrand_scalar   => init_kiss_scalar

  procedure :: write_checkpoint      => write_checkpoint_kiss
  procedure :: read_checkpoint       => read_checkpoint_kiss

  procedure :: class_name            => class_name_kiss

end type t_random_kiss

public :: t_random_kiss

contains

function class_name_kiss( this )
  
  implicit none

  class( t_random_kiss ), intent(in) :: this
  character( len = p_maxlen_random_class_name ) :: class_name_kiss

  class_name_kiss = "KISS"

end function class_name_kiss

!-------------------------------------------------------------------------------
! generates a integer random number on the [-2147483648, 2147483647] interval
! (all possible 32 bit integers)               -(2**31), (2**31)-1
!-------------------------------------------------------------------------------
function kiss_int32( this, ix )
  
  implicit none
  
  class( t_random_kiss ), intent(inout) :: this
  integer, dimension(:), optional  ::  ix
  integer :: kiss_int32

  integer, parameter :: p_int64 = selected_int_kind(10)

  integer(p_int64) :: t

  this%x = 69069 * this%x + 12345
  this%y = ieor( this%y, ishft(this%y, +13))
  this%y = ieor( this%y, ishft(this%y, -17))
  this%y = ieor( this%y, ishft(this%y, +5))

  t = 698769069_p_int64 * this%z + this%c
  this%z = t
  this%c = ishft( t, +32 )

  kiss_int32 = this%x + this%y + this%z

end function kiss_int32

!-------------------------------------------------------------------------------
! Initializes random number generator seed with a scalar integer 
!-------------------------------------------------------------------------------
subroutine init_kiss_scalar(this, s)
 
  implicit none
  
  class( t_random_kiss ), intent(inout) :: this
  integer, intent(in) :: s
  
  this%x = ieor( int(z'80000000'), s )
  this%y = 362436000
  this%z = 521288629
  this%c = 1812433253 * ieor( this%x, ishft(this%x, -30) )
    
end subroutine init_kiss_scalar


!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine write_checkpoint_kiss( this, restart_handle )

    implicit none

    class( t_random_kiss ), intent(in) :: this
    class( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for random number generator.'
    integer :: ierr

    ! Restart tag
    restart_io_wr( p_rand_rst_id, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    ! State data
    restart_io_wr( this%x, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%y, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%z, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%c, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    ! Write gaussian distribution data. This is done here for simplicity
    restart_io_wr( this%gaussian_set, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%new_gaussian, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

end subroutine write_checkpoint_kiss

!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine read_checkpoint_kiss( this, restart_handle )
    implicit none

    class( t_random_kiss ), intent(inout) :: this
    class( t_restart_handle ), intent(in) :: restart_handle

    integer :: ierr
    character(len=len(p_rand_rst_id)) :: rst_id
    character(len=*), parameter :: err_msg = 'error reading restart data for random number generator'

    restart_io_rd( rst_id, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    ! check if restart file is compatible
    if ( rst_id /= p_rand_rst_id) then
        ERROR('Corrupted restart file, or restart file ')
        ERROR('from incompatible binary (random number generator)')
        call abort_program(p_err_rstrd)
    endif

    ! State data
    restart_io_rd( this % x, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    restart_io_rd( this % y, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    restart_io_rd( this % z, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    restart_io_rd( this % c, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    ! Read gaussian distribution data. This is done here for simplicity
    restart_io_rd( this%gaussian_set, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    restart_io_rd( this%new_gaussian, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

end subroutine read_checkpoint_kiss


end module m_random_kiss