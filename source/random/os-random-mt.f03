!-------------------------------------------------------------------------------
! Fortran implementation of the Mersenne Twister random number generator
! (period 2**19937-1). For details on the random number generator see:
!
! http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
!
! and references therein.
!
! Ported to fortran from the distribution mt19937ar.c
!
! Period = 2^19937-1  
!
! Memory: 625 x 32 bit integers (2.5 kb)
! Speed: ~ 1078 M numbers/s
!-------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_random_mt

use m_parameters
use m_restart
use m_random_class  

private

! period parameters
integer, parameter :: N = 624
integer, parameter :: M = 397

integer, parameter :: MATRIX_A   = int(z'9908b0df') ! constant vector a 
integer, parameter :: UPPER_MASK = int(z'80000000') ! most significant w-r bits 
integer, parameter :: LOWER_MASK = int(z'7fffffff') ! least significant r bits 


! default seed 
integer, parameter :: DEFAULT_SEED = 5489      

! mag01(x) = x * MATRIX_A  for x=0,1 
integer, parameter, dimension(0:1) :: mag01 = (/ 0, MATRIX_A /)
  

! string to id restart data
character(len=*), parameter :: p_rand_rst_id = "random_mt rst data - 0x0000"

type, extends( t_random ) :: t_random_mt

    ! state data
    integer, dimension(N) :: mt   ! the array for the state vector
    integer :: mti

contains

    procedure :: genrand_int32         => mt_int32
    procedure :: init_genrand_scalar   => init_mt_scalar
    
    procedure :: write_checkpoint      => write_checkpoint_mt
    procedure :: read_checkpoint       => read_checkpoint_mt

    procedure :: class_name            => class_name_mt

end type t_random_mt

public :: t_random_mt

contains

function class_name_mt( this )
  
  implicit none

  class( t_random_mt ), intent(in) :: this
  character( len = p_maxlen_random_class_name ) :: class_name_mt

  class_name_mt = "Mersenne Twister"

end function class_name_mt

!-------------------------------------------------------------------------------
! generates a integer random number on the [-2147483648, 2147483647] interval
! (all possible 32 bit integers)               -(2**31), (2**31)-1
!-------------------------------------------------------------------------------
function mt_int32( this, ix )
  
  implicit none
  
  class( t_random_mt ), intent(inout) :: this
  integer, dimension(:), optional  ::  ix
  integer :: mt_int32
  integer :: y
  
  integer :: kk
  
  if (this%mti > N) then ! generate N words at one time 
    
    do kk = 1, N-M
      y = ior(iand(this%mt(kk),UPPER_MASK),iand(this%mt(kk+1),LOWER_MASK))
      this%mt(kk) = ieor(ieor(this%mt(kk + M),ishft( y, -1)), mag01(iand(y,1)))
    enddo
 
    do kk = N-M+1, N-1
      y = ior(iand(this%mt(kk),UPPER_MASK),iand(this%mt(kk+1),LOWER_MASK))
      this%mt(kk) = ieor(ieor(this%mt(kk+(M-N)),ishft(y,-1)), mag01(iand(y,1)))
    enddo
  
    y = ior(iand(this%mt(N),UPPER_MASK),iand(this%mt(1),LOWER_MASK))
    this%mt(N) = ieor(ieor(this%mt(M),ishft(y,-1)), mag01(iand(y,1)))

    this % mti = 1
  endif
  
  y = this%mt( this%mti )
  this%mti = this%mti + 1

  ! Tempering 
  y = ieor(y, ishft( y, -11))
  y = ieor(y, iand( ishft( y, 7 ), int(z'9d2c5680')))
  y = ieor(y, iand( ishft( y, 15 ), int(z'efc60000')))
  y = ieor(y, ishft( y, -18))

  mt_int32 = y

end function mt_int32

!-------------------------------------------------------------------------------
! Initializes random number generator seed with a scalar integer 
!-------------------------------------------------------------------------------
subroutine init_mt_scalar(this, s)
 
  implicit none
  
  class( t_random_mt ), intent(inout) :: this
  integer, intent(in) :: s
  integer :: i
  
  this%mt(1) = ieor( s, int(z'80000000') )
  do i = 2, N 
    
    this%mt(i) = (1812433253 * ieor( this%mt(i-1), ishft(this%mt(i-1), -30)) + i-1)
    ! See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. 
    ! In the previous versions, MSBs of the seed affect   
    !only MSBs of the array mt[].                        

  enddo

  this%mti = N
    
end subroutine init_mt_scalar

!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine write_checkpoint_mt( this, restart_handle )

    implicit none

    class( t_random_mt ), intent(in) :: this
    class( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for random number generator.'
    integer :: ierr

    restart_io_wr( p_rand_rst_id, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    ! State data
    restart_io_wr( this%mti, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%mt, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    ! Write gaussian distribution data. This is done here for simplicity
    restart_io_wr( this%gaussian_set, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%new_gaussian, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

end subroutine write_checkpoint_mt

!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine read_checkpoint_mt( this, restart_handle )
    implicit none

    class( t_random_mt ), intent(inout) :: this
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
    restart_io_rd( this % mti, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    restart_io_rd( this % mt, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    ! Read gaussian distribution data. This is done here for simplicity
    restart_io_rd( this%gaussian_set, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

    restart_io_rd( this%new_gaussian, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

end subroutine read_checkpoint_mt


end module m_random_mt