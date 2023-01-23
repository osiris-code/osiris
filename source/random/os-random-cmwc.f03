!-------------------------------------------------------------------------------
! CMCW Complimentary Multiply With Carry Random Number Generator
! Marsaglia, George (2003) "Random Number Generators," Journal of Modern Applied Statistical Methods: Vol. 2: Iss. 1, Article 2.
!
! Period > 2^131086 ~ 10^39461 
!
! Memory: 4096 x 32 bit integers + 1 x 64 bit integer ~ 16 kbytes
! Speed: ~ 700 M numbers/s
!-------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_random_cmwc

use m_parameters
use m_restart
use m_random_class  

private  

integer, parameter :: CMWC_CYCLE = 4096
integer, parameter :: CMWC_C_MAX = 809430660 

! string to id restart data
character(len=*), parameter :: p_rand_rst_id = "random_cmwc rst data - 0x0000"

type, extends( t_random ) :: t_random_cmwc

    ! state data
    integer, dimension(CMWC_CYCLE) :: Q
    integer :: i = CMWC_CYCLE

    integer(p_int64) :: c = 362436

contains

    procedure :: genrand_int32         => cmwc_int32
    procedure :: init_genrand_scalar   => init_cmwc_scalar
    
    procedure :: write_checkpoint      => write_checkpoint_cmwc
    procedure :: read_checkpoint       => read_checkpoint_cmwc

    procedure :: class_name            => class_name_cmwc

end type t_random_cmwc

public :: t_random_cmwc

contains

function class_name_cmwc( this )
  
  implicit none

  class( t_random_cmwc ), intent(in) :: this
  character( len = p_maxlen_random_class_name ) :: class_name_cmwc

  class_name_cmwc = "Complimentary Multiply With Carry Random Number Generator (CMWC)"

end function class_name_cmwc

!-------------------------------------------------------------------------------
! generates a integer random number on the [-2147483648, 2147483647] interval
! (all possible 32 bit integers)               -(2**31), (2**31)-1
!-------------------------------------------------------------------------------
function cmwc_int32( this, ix )
  
  implicit none
  
  class( t_random_cmwc ), intent(inout) :: this
  integer, dimension(:), optional  ::  ix
  integer :: cmwc_int32

    integer(p_int64) :: t, x

    integer(p_int64), parameter :: p_mask32 = int(z'FFFFFFFF',kind=p_int64)

 
    this%i = this%i + 1
    if ( this%i > CMWC_CYCLE ) this%i = 1

    ! The iand operation converts Q(i) to an unsigned int32 value
    t = 18782_p_int64 * iand( p_mask32, int( this%Q( this%i ), kind = p_int64 )) + this%c

    this%c = ishft( t, -32 )
    x = iand( p_mask32, t + this%c )

    if ( x < this%c ) then
        x = x + 1
        this%c = this%c + 1
    endif

    this%Q(this%i) = int( int(z'FFFFFFFE',kind=p_int64) - x )
    cmwc_int32 = this%Q(this%i)

end function cmwc_int32

!-------------------------------------------------------------------------------
! Initializes random number generator seed with a scalar integer 
!-------------------------------------------------------------------------------
subroutine init_cmwc_scalar(this, s)
 
  implicit none
  
  class( t_random_cmwc ), intent(inout) :: this
  integer, intent(in) :: s
  
  integer :: l, m_z, m_w

  ! set seed of simple MWC prng for initialization
  m_z = ieor( int(z'80000000'), s )
  m_w = 1812433253 * ieor( m_z, ishft( m_z, -30) )

  do l = 1, CMWC_CYCLE
      this % Q(l) = lrnd( m_z, m_w )
  enddo

  do
      this % c = lrnd( m_z, m_w )
      if ( this % c < CMWC_C_MAX ) exit
  enddo

  contains

  function lrnd( m_z, m_w )
      implicit none

      integer :: lrnd

      integer, parameter :: p_mask16 = int(z'ffff') ! least significant 16 bits 

      integer, intent(inout) :: m_z, m_w

      m_z = 36969 * iand( m_z, p_mask16 ) + ishft( m_z, -16 )
      m_w = 18000 * iand( m_w, p_mask16 ) + ishft( m_w, -16 )

      lrnd = ishft( m_z, 16 ) + m_w

  end function lrnd
    
end subroutine init_cmwc_scalar


!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine write_checkpoint_cmwc( this, restart_handle )

    implicit none

    class( t_random_cmwc ), intent(in) :: this
    class( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for random number generator.'
    integer :: ierr

    ! Restart tag
    restart_io_wr( p_rand_rst_id, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    ! State data
    restart_io_wr( this%Q, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%i, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%c, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    ! Write gaussian distribution data. This is done here for simplicity
    restart_io_wr( this%gaussian_set, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%new_gaussian, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

end subroutine write_checkpoint_cmwc

!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine read_checkpoint_cmwc( this, restart_handle )
    implicit none

    class( t_random_cmwc ), intent(inout) :: this
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
    restart_io_rd( this % Q, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    restart_io_rd( this % i, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    restart_io_rd( this % c, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    ! Read gaussian distribution data. This is done here for simplicity
    restart_io_rd( this%gaussian_set, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    restart_io_rd( this%new_gaussian, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

end subroutine read_checkpoint_cmwc


end module m_random_cmwc