!-------------------------------------------------------------------------------
! Marsaglia MWC (multiply-with-carry) PRNG
!
! Period : ? 
!
! Memory: 2 x 32 bit integers (8 bytes)
! Speed: ~ 620 M numbers/s
!-------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_random_mwc

use m_parameters
use m_restart
use m_random_class  

private

! string to id restart data
character(len=*), parameter :: p_rand_rst_id = "random_mwc rst data - 0x0000"

type, extends( t_random ) :: t_random_mwc

    ! state data
    integer :: m_w = 12345
    integer :: m_z = 67890

contains

    procedure :: genrand_int32         => mwc_int32
    procedure :: init_genrand_scalar   => init_mwc_scalar
    
    procedure :: write_checkpoint      => write_checkpoint_mwc
    procedure :: read_checkpoint       => read_checkpoint_mwc

    procedure :: class_name            => class_name_mwc

end type t_random_mwc

public :: t_random_mwc

contains

function class_name_mwc( this )
  
  implicit none

  class( t_random_mwc ), intent(in) :: this
  character( len = p_maxlen_random_class_name ) :: class_name_mwc

  class_name_mwc = "Multiply With Carry Random Number Generator (MWC)"

end function class_name_mwc

!-------------------------------------------------------------------------------
! generates a integer random number on the [-2147483648, 2147483647] interval
! (all possible 32 bit integers)               -(2**31), (2**31)-1
!-------------------------------------------------------------------------------
function mwc_int32( this, ix )
  
  implicit none
  
  class( t_random_mwc ), intent(inout) :: this
  integer, dimension(:), optional  ::  ix
  integer :: mwc_int32

  integer, parameter :: MASK = int(z'ffff') ! least significant 16 bits 

  ! Reference C implementation
  !  m_z = 36969 * (m_z & 65535) + (m_z >> 16);
  !  m_w = 18000 * (m_w & 65535) + (m_w >> 16);
  !  return (m_z << 16) + m_w;  /* 32-bit result */
  
  this % m_z = 36969 * iand( this % m_z, MASK ) + ishft( this % m_z, -16 )
  this % m_w = 18000 * iand( this % m_w, MASK ) + ishft( this % m_w, -16 )
  
  mwc_int32 = ishft( this % m_z, 16 ) + this % m_w

end function mwc_int32

!-------------------------------------------------------------------------------
! Initializes random number generator seed with a scalar integer 
!-------------------------------------------------------------------------------
subroutine init_mwc_scalar(this, s)
 
  implicit none
  
  class( t_random_mwc ), intent(inout) :: this
  integer, intent(in) :: s
  
  this % m_z = ieor( int(z'80000000'), s )

  ! In the unlikely event that a user uses s = 0x80000000 set a default seed
  ! to prevent a bad sequence
  if( this % m_z == 0 ) this % m_z = 123456789

  this % m_w = 1812433253 * ieor( this % m_z, ishft(this % m_z, -30) )
    
end subroutine init_mwc_scalar


!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine write_checkpoint_mwc( this, restart_handle )

    implicit none

    class( t_random_mwc ), intent(in) :: this
    class( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for random number generator.'
    integer :: ierr

    restart_io_wr( p_rand_rst_id, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    ! State data
    restart_io_wr( this%m_w, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%m_z, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    ! Write gaussian distribution data. This is done here for simplicity
    restart_io_wr( this%gaussian_set, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%new_gaussian, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

end subroutine write_checkpoint_mwc

!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine read_checkpoint_mwc( this, restart_handle )
    implicit none

    class( t_random_mwc ), intent(inout) :: this
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
    restart_io_rd( this % m_w, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    restart_io_rd( this % m_z, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    ! Read gaussian distribution data. This is done here for simplicity
    restart_io_rd( this%gaussian_set, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    restart_io_rd( this%new_gaussian, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

end subroutine read_checkpoint_mwc


end module m_random_mwc