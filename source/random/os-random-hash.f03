!-------------------------------------------------------------------------------
! Implements the Hashing RNG from Numerical Recipes Chapter 7.5
! Honestly Josh doesn't think it's very good, but it's standard and cheap
!
! Josh, 2016
! 
!-------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_random_hash

use m_parameters
use m_restart
use m_random_class  

private  

! string to id restart data
character(len=*), parameter :: p_rand_rst_id = "random_hash rst data - 0x0000"

type, extends( t_random ) :: t_random_hash

    ! state data
    integer :: seed
    integer :: index

contains

    procedure :: genrand_int32         => hash_int32
    procedure :: init_genrand_scalar   => init_hash_scalar
    
    procedure :: write_checkpoint      => write_checkpoint_hash
    procedure :: read_checkpoint       => read_checkpoint_hash

    procedure :: class_name            => class_name_hash
    
    procedure :: has_cycled
    
end type t_random_hash

public :: t_random_hash

contains

function has_cycled( this )
  
  implicit none
  
  class( t_random_hash ), intent(in) :: this
  logical :: has_cycled
  
  ! Technically this actually only means it has half-cycled, but that's hackish to begin with;
  ! and this way the test can be done only occasionally with a decent safety factor
  has_cycled = (this%index .lt. 0)
  
end function has_cycled


function class_name_hash( this )
  
  implicit none

  class( t_random_hash ), intent(in) :: this
  character( len = p_maxlen_random_class_name ) :: class_name_hash

  class_name_hash = "Hashing RNG"

end function class_name_hash


!-------------------------------------------------------------------------------
! generates a integer random number on the [-2147483648, 2147483647] interval
! (all possible 32 bit integers)               -(2**31), (2**31)-1
!-------------------------------------------------------------------------------
function hash_int32( this, ix )
  
  implicit none
  
  class( t_random_hash ), intent(inout) :: this
  integer, dimension(:), optional  ::  ix
  integer :: hash_int32

  integer, parameter :: NITER = 4
  !integer, dimension(NITER), parameter :: c1 = (/ -z'45569779', z'1e17d32c', z'03bcdc3c', z'0f33d1b2' /)
  !integer, dimension(NITER), parameter :: c2 = (/ z'4b0f3b58', -z'178b0f3d', z'6955c5a6', z'55a7ca46' /)

  integer, dimension(NITER), parameter :: c1 = (/ -1163302777,  504877868,   62708796,  255054258 /)
  integer, dimension(NITER), parameter :: c2 = (/  1259289432, -394989373, 1767228838, 1437059654 /)

  integer, parameter :: mask = int(z'ffff')


  integer :: i, ia, ib, iswap, itmph, itmpl, lword, irword
  
  lword = this%seed
  irword = this%index
  
  do i = 1, NITER
    iswap = irword
    ia = ieor(irword, c1(i))
    itmpl = iand(ia, mask)
    itmph = ishft(ia, -16)
    ib = itmpl*itmpl + not(itmph*itmph)
    ia = ishft(ib, -16)
    ib = ishft(ib, 16)
    irword = ieor(ieor((ia + ib), c2(i)) + itmpl*itmph, lword)
    lword = iswap
  enddo
  
  this%index = this%index + 1

  hash_int32 = irword


end function hash_int32

!-------------------------------------------------------------------------------
! Initializes random number generator seed with a scalar integer 
!-------------------------------------------------------------------------------
subroutine init_hash_scalar(this, s)
 
  implicit none
  
  class( t_random_hash ), intent(inout) :: this
  integer, intent(in) :: s
  
  this%seed = s
  this%index = 1
  
end subroutine init_hash_scalar


!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine write_checkpoint_hash( this, restart_handle )

    implicit none

    class( t_random_hash ), intent(in) :: this
    class( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for random number generator.'
    integer :: ierr

    restart_io_wr( p_rand_rst_id, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%seed, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%index, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    ! Write gaussian distribution data. This is done here for simplicity
    restart_io_wr( this%gaussian_set, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%new_gaussian, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

end subroutine write_checkpoint_hash

!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine read_checkpoint_hash( this, restart_handle )
    implicit none

    class( t_random_hash ), intent(inout) :: this
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

    restart_io_rd( this % seed, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    restart_io_rd( this % index, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    ! Read gaussian distribution data. This is done here for simplicity
    restart_io_rd( this%gaussian_set, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    restart_io_rd( this%new_gaussian, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

end subroutine read_checkpoint_hash


end module m_random_hash