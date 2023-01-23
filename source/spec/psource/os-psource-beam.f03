#include "os-config.h"
#include "os-preprocess.fpp"

module m_psource_beam

#include "memory/memory.h"

use m_species_define
use m_parameters
use m_fparser
use m_psource_std
use m_current_define
use m_node_conf

implicit none

private

interface transform_focal
  module procedure transform_focal_beam
end interface

integer, parameter :: p_gaussian  = 2  ! gaussian

type, extends(t_psource_std) :: t_psource_beam


! Focal plane distance after acceleration
real(p_k_part), dimension(p_p_dim)   :: focal_dist

! Params to rotate p2x2 and px3x3 phase space
real(p_k_part), dimension(p_p_dim)   :: alpha

! Beam uth parameter at Focal Plane
real(p_k_part), dimension(p_p_dim)   :: uth

real(p_k_part) :: gamma



contains

  ! Methods from abstract class
  procedure :: if_inject => if_inject_beam
  procedure :: read_input => read_input_beam
  procedure :: cleanup => cleanup_beam
  procedure :: inject => inject_beam

  ! Additional methods
  ! Allow access to density values from outside
  ! procedure :: get_den_value => get_den_value_beam
  ! procedure :: den_value => den_value_beam



end type t_psource_beam



public :: t_psource_beam



contains

!-------------------------------------------------------------------------------
! Read information from input file
!-------------------------------------------------------------------------------
subroutine read_input_beam( this, input_file, coordinates )
    use m_input_file
  
    implicit none
  
    class( t_psource_beam ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file
    integer, intent(in) :: coordinates
  
  
    ! minimum density to inject particles
    real(p_k_part)                      :: den_min
  
    ! global density (multiplies the selected density profile)
    real(p_k_part)               :: density
  
    ! variables needed for Gaussian profiles
    real(p_k_part), dimension(p_x_dim)   :: gauss_center
    real(p_k_part), dimension(p_x_dim)   :: gauss_sigma
    integer, dimension(  p_x_dim )       :: gauss_n_sigma
    real(p_k_part), dimension(2,p_x_dim) :: gauss_range
      
    ! This is only used for t_psource_constq but is included here for simplicity
    ! integer, dimension(p_max_dim) :: sample_rate

    real(p_double), dimension(p_p_dim) :: focal_dist, alpha
    real(p_double), dimension(p_p_dim) :: uth
    real(p_double) :: gamma

  
    ! namelists of input variables
    namelist /nl_profile/ den_min, density, gauss_n_sigma, gauss_range, &
                          gauss_center, gauss_sigma, focal_dist, alpha, &
                          uth, gamma
  
    integer  ::  j
    integer :: ierr
  
    den_min = 0
    density = 1.0
    gauss_center  = 0.0_p_k_part
    gauss_sigma   = huge( 1.0_p_k_part )
    gauss_n_sigma = 0
    gauss_range(p_lower,:) = -gauss_sigma
    gauss_range(p_upper,:) =  gauss_sigma

    focal_dist = 0.0
    alpha = 0.0

    uth = 0
    gamma = 1
  
    ! Get namelist text from input file
    call get_namelist( input_file, "nl_profile", ierr )
  
    if ( ierr == 0 ) then
      read (input_file%nml_text, nml = nl_profile, iostat = ierr)
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*) ""
          write(0,*) "   Error reading profile information"
          write(0,*) "   aborting..."
        endif
        stop
      endif
    else
       SCR_ROOT("   - profile parameters missing, using uniform density")
    endif
    
    ! minimum density for injection
    this%den_min      = den_min
    this%density = density
    this % if_mult = .true.
    do j = 1, p_x_dim
        this%type(j) = p_gaussian
  
        this%gauss_center(j) = gauss_center(j)
        this%gauss_sigma(j) = abs( gauss_sigma(j) )
        if ( gauss_n_sigma(j) > 0 ) then
            this%gauss_range(p_lower,j) = gauss_center(j)-abs(gauss_n_sigma(j)*gauss_sigma(j))
            this%gauss_range(p_upper,j) = gauss_center(j)+abs(gauss_n_sigma(j)*gauss_sigma(j))
        else
            this%gauss_range(:,j) = gauss_range(:,j)
        endif  
    enddo
  
    this%uth = uth
    this%alpha = alpha
    this%focal_dist = focal_dist
    this%gamma = gamma 
    call transform_focal(this)
end subroutine read_input_beam
!---------------------------------------------------------------------------------------------------




!---------------------------------------------------------------------------------------------------
! Returns .true. if this psource injects particles
!---------------------------------------------------------------------------------------------------
function if_inject_beam( this )

    implicit none
  
    logical :: if_inject_beam
    class( t_psource_beam ), intent( in ) :: this
  
    if_inject_beam = .true.
  
end function if_inject_beam
!---------------------------------------------------------------------------------------------------
  
!---------------------------------------------------------------------------------------------------
! Cleanup t_psource object
!---------------------------------------------------------------------------------------------------
subroutine cleanup_beam( this )

    implicit none
  
    class( t_psource_beam ), intent( inout ) :: this

    call this % t_psource_std % cleanup()
  
end subroutine cleanup_beam
!---------------------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! Transform the particle beam (3D/Quasi-3D/rz version)
!-----------------------------------------------------------------------------------------
subroutine transform_focal_beam( this )

    implicit none 

    class( t_psource_beam ), intent( inout ) :: this

    real(p_double)  :: sigma_mul, p1, dens_mul
    integer :: j
  
    p1 = sqrt( this % gamma ** 2 - 1 )
    dens_mul = 1.0_p_k_part
      
    do j = 2, p_x_dim
        ! focal plane distance defined after n_acceleration
        this%focal_dist(j) = this%focal_dist(j) 
        !transform params
        sigma_mul = sqrt(1 + this%focal_dist(j)**2 *this%uth(j)**2 /(p1**2 * this%gauss_sigma(j)**2 ))
        dens_mul  = dens_mul * sigma_mul
        this%gauss_sigma(j) = sigma_mul * this%gauss_sigma(j)
        this%gauss_range(:,j) = sigma_mul * this%gauss_range(:,j)
        this%uth(j) = this%uth(j)/sigma_mul
        if(this%focal_dist(j) /= 0.0) then
            this%alpha(j) = p1 / this%focal_dist(j) * (1.0_p_k_part- 1.0_p_k_part/sigma_mul**2)
        endif
    end do
    ! transform density
    if(p_x_dim <= 2) then
        !if using r-z; need to check for p_phi and scale it; angular action is conserved
        if(p_p_dim > p_x_dim) then
            this%uth(p_p_dim) =  this%uth(p_p_dim)/dens_mul 
        endif
        this%density = this%density / dens_mul**2
    else
        this%density = this%density / dens_mul
    endif
   
end subroutine transform_focal_beam
!-----------------------------------------------------------------------------------------
  


function inject_beam( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, &
                              send_msg, recv_msg ) result(num_inj)

    use m_species_udist
    use m_random
    implicit none

    class(t_psource_beam), intent(inout) :: this
    class(t_species), intent(inout), target :: species
    integer, dimension(:, :), intent(in) :: ig_xbnd_inj
    class( t_current ), intent(inout)   :: jay
    class( t_node_conf ), intent(in)     :: no_co
    type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
    type( t_part_idx ), intent(inout) :: node_cross
    type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

    integer :: num_inj, i, j
    real(p_double), dimension( p_x_dim ) :: x
    real(p_double) :: u1

    ! Call superclass to inject particles
    num_inj = this % t_psource_std % inject( species, ig_xbnd_inj, jay, no_co, &
                                             bnd_cross, node_cross, send_msg, recv_msg )
    
    do i = species%num_par+1, species%num_par + num_inj
        species % p(1,i) = this % uth(1) * rng % genrand_gaussian( )
        species % p(2,i) = this % uth(2) * rng % genrand_gaussian( )
        species % p(3,i) = this % uth(3) * rng % genrand_gaussian( )
    enddo    
    
    if ( any( this % focal_dist /= 0 ) ) then
        do i = species%num_par+1, species%num_par + num_inj
            ! 3D/rz geometry
            call species % get_position( i, x )
            do j = 1, p_x_dim
                species%p(j,i) = species%p(j,i) - x(j) * this%alpha(j)
            enddo
        enddo
    endif

    u1 = sqrt( this % gamma ** 2 - 1 )
    species%udist%ufl(1) = u1
    do i = species%num_par+1, species%num_par + num_inj
        species%p(1,i) = species%p(1,i) + u1
    enddo


end function inject_beam

end module m_psource_beam