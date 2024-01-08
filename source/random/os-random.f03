#include "os-config.h"
#include "os-preprocess.fpp"

module m_random

use m_grid_define

use m_random_class ! Virtual base class
use m_random_mt    ! Mersenne Twister
use m_random_r250  ! R250 PRNG
use m_random_cmwc  ! Marsaglia Cyclic Multiply with Carry
use m_random_mwc   ! Marsaglia Multiply with Carry
use m_random_kiss  ! Marsaglia KISS
use m_random_hash  ! Hash

use m_random_class_grid

private

! string to id restart data
character(len=*), parameter :: p_rand_rst_id = "random rst data - 0x0001"

integer, public, parameter :: p_random_mt      = 1
integer, public, parameter :: p_random_r250    = 2
integer, public, parameter :: p_random_cmwc    = 3
integer, public, parameter :: p_random_mwc     = 4
integer, public, parameter :: p_random_kiss    = 5
integer, public, parameter :: p_random_hash    = 6

! The default is currently set to the Mersenne Twister algorithm
integer, public, parameter :: p_random_default = p_random_mt

class(t_random), pointer, public :: rng => null()

interface new_random
    module procedure new_random
end interface new_random

interface init_genrand
    module procedure init_genrand
end interface

interface cleanup_genrand
    module procedure cleanup_genrand
end interface

interface random_type
    module procedure random_type_rnd
    module procedure random_type_str
end interface

interface enforce_random_constancy
    module procedure enforce_random_constancy
end interface

interface write_checkpoint_rnd
    module procedure write_checkpoint_rnd
end interface

interface list_algorithm_rnd
    module procedure list_algorithm_rnd
end interface

public :: t_random, new_random, random_type
public :: init_genrand, cleanup_genrand
public :: write_checkpoint_rnd, list_algorithm_rnd

contains

function random_type_rnd( this )
    implicit none

    integer :: random_type_rnd
    class( t_random ), intent(in),target :: this
    class( t_random ), pointer :: rnd_class

    rnd_class => this
    select type( this )
      class is ( t_random_grid ) 
        rnd_class => this%base_rng
    end select

    select type( rnd_class )
      
        class is ( t_random_mt )  
            random_type_rnd = p_random_mt

        class is ( t_random_r250 )  
            random_type_rnd = p_random_r250
               
        class is ( t_random_cmwc )
            random_type_rnd = p_random_cmwc

        class is ( t_random_mwc )
            random_type_rnd = p_random_mwc

        class is ( t_random_kiss )  
            random_type_rnd = p_random_kiss
     
        class is ( t_random_hash )  
            random_type_rnd = p_random_hash
     
        class default
            random_type_rnd = -1

    end select    

end function random_type_rnd

function random_type_str( str )
    implicit none

    integer :: random_type_str
    character(len=*), intent(in) :: str

    select case ( trim( str ) )
    case('mt')
        random_type_str = p_random_mt
    case('r250')
        random_type_str = p_random_r250
    case('cmwc')
        random_type_str = p_random_cmwc
    case('mwc')
        random_type_str = p_random_mwc
    case('kiss')
        random_type_str = p_random_kiss
    case('hash')
        random_type_str = p_random_hash
    case('default')
        random_type_str = p_random_default
    case default
        random_type_str = -1
    end select

end function random_type_str

function enforce_random_constancy(this)
    implicit none
    class( t_random ), intent(in) :: this
    logical :: enforce_random_constancy

    enforce_random_constancy = this%enforce_random_constancy

end function enforce_random_constancy

subroutine init_genrand( class_idx, enforce_rng_constancy, base_seed, seed, grid, restart, restart_handle )
    
    use m_restart
    use m_grid

    implicit none

    integer, intent(inout) :: class_idx
    integer, intent(in) :: base_seed, seed
    logical, intent(inout) :: enforce_rng_constancy
    class( t_grid ), intent(in)      ::  grid

    logical, intent(in) :: restart
    type( t_restart_handle ), intent(in) :: restart_handle   
    class(t_random_grid), pointer :: new_random_grid

    if ( enforce_rng_constancy ) then 
      ! create and initalize a grid (i.e. array) of rng objects.
      if ( restart ) then
        ! read the restart data which fills in the 'class_idx'
        !  of the base rng type. Then make one on these.
        call read_checkpoint_rnd( restart_handle, class_idx, enforce_rng_constancy)
        rng => new_random( class_idx )
        call rng % init_genrand_scalar( seed )
        
        ! now create the grid of rngs
        allocate( t_random_grid :: new_random_grid )
        new_random_grid % base_rng => rng
        rng => new_random_grid

        ! and read in the restart data for each of the rng in the grid.
        call rng % init_genrand_scalar_grid( base_seed, grid,  class_idx )
        call rng % read_checkpoint( restart_handle )
      else
        ! create an rng of the type requested by the user and save it as a template
        !   that will be duplicated for each spatial cell.
        rng => new_random( class_idx, enforce_rng_constancy )
        allocate( t_random_grid :: new_random_grid )
        new_random_grid % base_rng => rng
        rng => new_random_grid

        call rng % init_genrand_scalar_grid( base_seed, grid, class_idx )
      endif
    else
      if ( restart ) then
        
        call read_checkpoint_rnd( restart_handle, class_idx, enforce_rng_constancy)
        rng => new_random( class_idx )
        call rng % read_checkpoint( restart_handle )

      else
        ! Create random number object of the appropriate class
        rng => new_random( class_idx, enforce_rng_constancy )
        
        ! Initialize 
        call rng % init_genrand_scalar( seed )
      endif
    endif

end subroutine init_genrand

subroutine cleanup_genrand( )

    implicit none

    deallocate( rng )
    rng => null()

end subroutine cleanup_genrand

subroutine read_checkpoint_rnd( restart_handle, class_idx, enforce_rng_constancy )

    use m_restart
    use m_system
    use m_parameters, only : p_err_rstrd

  implicit none
  
  class( t_restart_handle ), intent(in) :: restart_handle
  integer, intent(inout) :: class_idx
  logical, intent(inout) :: enforce_rng_constancy

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
         
  restart_io_rd( class_idx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

  restart_io_rd( enforce_rng_constancy, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd ) 

end subroutine read_checkpoint_rnd


subroutine write_checkpoint_rnd( this, restart_handle )

    use m_restart
    use m_system
    use m_parameters, only : p_err_rstwrt

    implicit none

    class( t_random ), intent(in) :: this
    class( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for random number generator.'
    integer :: class_idx, ierr
    logical :: enforce_random_constancy

    ! 
    restart_io_wr( p_rand_rst_id, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    ! write integer identifying the prng type
    class_idx = random_type( this )    
    restart_io_wr( class_idx, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%enforce_random_constancy, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    ! Write object specific info
    call this % write_checkpoint( restart_handle )

end subroutine write_checkpoint_rnd

function new_random( class_idx_, enforce_rng_constancy_ )

    implicit none

    class(t_random), pointer :: new_random

    integer, intent(in), optional :: class_idx_
    logical, intent(in), optional :: enforce_rng_constancy_

    integer :: class_idx
    logical :: enforce_rng_constancy

    if ( present(class_idx_) ) then 
        class_idx = class_idx_
    else
        class_idx = p_random_default
    endif

    if ( present(enforce_rng_constancy_) ) then 
        enforce_rng_constancy = enforce_rng_constancy_
    else
        enforce_rng_constancy = .false.
    endif

    select case ( class_idx )
        case( p_random_mt )
            allocate( t_random_mt :: new_random ) 

        case( p_random_r250 )
            allocate( t_random_r250 :: new_random ) 

        case( p_random_cmwc )
            allocate( t_random_cmwc  :: new_random ) 

        case( p_random_mwc )
            allocate( t_random_mwc  :: new_random ) 

        case( p_random_kiss )
            allocate( t_random_kiss :: new_random ) 

        case( p_random_hash )
            allocate( t_random_hash :: new_random ) 

        case default
            write(0,*) '(*warning*) Invalid random number generator type,'&
                       &' using default Mersenne Twister.'
            allocate( t_random_mt :: new_random )
           
    end select

end function new_random

subroutine list_algorithm_rnd()
    implicit none

    print *, ' '
    print *, 'Global random number generator'
    print *, ' - Algorithm: ', trim( rng % class_name() )

end subroutine list_algorithm_rnd



end module m_random