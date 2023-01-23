! Initialize spin distribution

#include "os-config.h"
#include "os-preprocess.fpp"

#ifdef __HAS_SPIN__

module m_species_sdist

#include "memory/memory.h"

use m_system
use m_species_define
use m_node_conf
use m_parameters

implicit none

private

integer, parameter, private :: p_part_block = 65536
integer, parameter, private :: p_num_sdir_max = 32

interface set_spin
  module procedure set_spin
end interface

interface read_nml
  module procedure read_nml_spec_sdist
end interface

interface cleanup
  module procedure cleanup_sdist
end interface


public :: set_spin, read_nml, cleanup

contains

!-----------------------------------------------------------------------------------------
! Read spin distribution parameters
!-----------------------------------------------------------------------------------------
subroutine read_nml_spec_sdist( this, input_file )

  use m_input_file

  implicit none

  type( t_sdist ), intent(out) :: this
  class( t_input_file ), intent(inout) :: input_file

  character(len=20) :: sdist_type

  real(p_double), dimension( p_num_sdir_max, p_s_dim ) :: sdir
  real(p_double), dimension( p_num_sdir_max ) :: pdf
  real(p_double) :: sdir_norm, pdf_norm, cum_sum, gauss_theta
  real(p_double), dimension( p_s_dim ) :: gauss_mu

  integer :: n_dirs

  integer :: i, ierr

  ! namelists of input variables

  namelist /nl_sdist/ sdist_type, sdir, n_dirs, pdf, gauss_mu, gauss_theta

  ! default values

  sdist_type = "isotropic"

  sdir = 0.0_p_double; sdir(1,1) = 1.0_p_double
  pdf  = 0.0_p_double; pdf(1) = 1.0_p_double
  n_dirs = 1
  gauss_mu = (/-huge(1.0_p_double), 0.0_p_double, 0.0_p_double/)
  gauss_theta = 0.0_p_double

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_sdist", ierr )

  if ( ierr == 0 ) then
    read (input_file%nml_text, nml = nl_sdist, iostat = ierr)
    if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) ""
        write(0,*) "   Error reading spin distribution information"
        write(0,*) "   aborting..."
      endif
      stop
    endif
  else
     SCR_ROOT("   - spin distribution parameters missing, initializing particles isotropically")
  endif

  ! sdist_type
  select case ( trim( sdist_type ) )

    case ( "isotropic" )
      this%sdist_type = p_isotropic

    case ( "fixed dir" )
      this%sdist_type = p_fixed_dir

      if ( n_dirs > 0 ) then
        this%n_dirs = n_dirs
        call alloc( this%sdir, (/ n_dirs, p_s_dim /) )
        call alloc( this%pdf, (/ n_dirs /) )
      else
        ERROR( 'n_dirs must be positive, aborting...' )
        call abort_program( p_err_invalid )
      endif

      ! normalize spin
      do i = 1, this%n_dirs
        sdir_norm = sqrt( sdir(i,1)**2 + sdir(i,2)**2 + sdir(i,3)**2 )
        if ( sdir_norm == 0.0 ) then
          ERROR( 'Invalid spin direction, aborting...' )
          call abort_program( p_err_invalid )
        else
          this%sdir(i,1) = sdir(i,1) / sdir_norm
          this%sdir(i,2) = sdir(i,2) / sdir_norm
          this%sdir(i,3) = sdir(i,3) / sdir_norm
        endif
      enddo

      ! normalize probability
      if ( any( pdf(1:n_dirs) < 0.0_p_double ) ) then
        ERROR( 'Invalid probability of spin direction, aborting...' )
        call abort_program( p_err_invalid )
      endif

      pdf_norm = sum( pdf(1:n_dirs) )
      if ( pdf_norm == 0.0_p_double ) then
        ERROR( 'Invalid probability of spin direction, must be greater than 0, aborting...' )
        call abort_program( p_err_invalid )
      endif
      pdf(1:n_dirs) = pdf(1:n_dirs) / pdf_norm

      ! calculate integrated pdf
      cum_sum = 0.0
      do i = 1, n_dirs
        this%pdf(i) = pdf(i) + cum_sum
        cum_sum = cum_sum + pdf(i)
      enddo

    case ( "spherical-gaussian" )
      this%sdist_type = p_spherical_gauss

      if ( gauss_mu(1) == -huge(1.0_p_double) ) then
        ERROR( 'Invalid spin direction for spherical-gaussian distribution, aborting...' )
        call abort_program( p_err_invalid )
      endif

      if ( gauss_theta == 0.0_p_double ) then
        ERROR( 'Invalid gauss_theta for spherical-gaussian distribution, aborting...' )
        call abort_program( p_err_invalid )
      endif

      sdir_norm = sqrt( gauss_mu(1)**2 + gauss_mu(2)**2 + gauss_mu(3)**2 )

      this%gauss_mu = gauss_mu / sdir_norm
      this%gauss_theta = gauss_theta

    case default
      if ( mpi_node() == 0 ) then
        write(0,*) "   Error reading spin distribution parameters"
        write(0,*) "   invalid sdist_type, valid values are 'isotropic' and 'fixed dir'"
        write(0,*) "   aborting..."
      endif
      stop
  end select

end subroutine read_nml_spec_sdist
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Sets the spin of the specified particle range
!-----------------------------------------------------------------------------------------
subroutine set_spin( species, idx0, idx1 )

  use m_random

  implicit none

  class( t_species ), intent(inout) :: species
  integer, intent(in) :: idx0, idx1

  real(p_k_part) :: rnd, theta, phi, s1, s2, s3, s_theta, c_theta, s_phi, c_phi
  real(p_k_part) :: gam

  integer :: i, k, l, npart, n_dirs, idx_dir

  real(p_double), parameter :: pi      = 3.14159265358979323846264_p_double

  ! Process particles in chunks
  do k = idx0, idx1, p_cache_size

     if ( k + p_cache_size - 1 > idx1 ) then
       npart = idx1 - k + 1
     else
       npart = p_cache_size
     endif

     ! set spin
     select case ( species%sdist%sdist_type )

       case ( p_isotropic )

         do i = 1, npart
            l = k + i - 1

            ! get s1
            call rng % harvest_real3( rnd, species%ix(:,l) )
            ! theta is subject to 0.5*sin(theta) distribution on (0,pi)
            theta = acos( 1.0_p_double - 2.0_p_double * rnd )

            call rng % harvest_real3( rnd, species%ix(:,l) )
            phi = real( 2 * pi, p_k_part ) * rnd

            species%s(1,l) = cos(theta)
            species%s(2,l) = sin(theta) * cos(phi)
            species%s(3,l) = sin(theta) * sin(phi)
         enddo

       case ( p_fixed_dir )

          n_dirs = species%sdist%n_dirs

          do i = 1, npart
            l = k + i - 1

            call rng % harvest_real3( rnd, species%ix(:,l) )
            idx_dir = 1
            do while ( rnd > species%sdist%pdf(idx_dir) )
              idx_dir = idx_dir + 1
            enddo

            species%s(1,l) = species%sdist%sdir(idx_dir,1)
            species%s(2,l) = species%sdist%sdir(idx_dir,2)
            species%s(3,l) = species%sdist%sdir(idx_dir,3)
          enddo

       case ( p_spherical_gauss )

          c_theta = species%sdist%gauss_mu(1)
          s_theta = sqrt( 1.0_p_double - c_theta**2 )
          if ( s_theta /= 0.0_p_double ) then
            c_phi = species%sdist%gauss_mu(2) / s_theta
            s_phi = species%sdist%gauss_mu(3) / s_theta
          else
            c_phi = 1.0_p_double
            s_phi = 0.0_p_double
          endif

          do i = 1, npart
            l = k + i - 1

            theta = species%sdist%gauss_theta * rng%genrand_gaussian( species%ix(:,l) )
            call rng%harvest_real3( rnd, species%ix(:,l) )
            phi = real( 2 * pi, p_k_part ) * rnd

            s1 = cos(theta)
            s2 = sin(theta) * cos(phi)
            s3 = sin(theta) * sin(phi)

            species%s(1,l) = c_theta * s1 - s_theta * s2
            species%s(2,l) = s_theta * c_phi * s1 + c_theta * c_phi * s2 - s_phi * s3
            species%s(3,l) = s_theta * s_phi * s1 + c_theta * s_phi * s2 + c_phi * s3

          enddo

       case default
          ERROR('Not implemented')
          call abort_program( p_err_notimplemented )

     end select

     ! if using exact pusher, transform spin to lab frame
    if ( species%push_type == p_exact ) then
      do i = 1, npart
        l = k + i - 1
        gam = sqrt( 1.0_p_k_part + species%p(1,l)**2 + species%p(2,l)**2 + species%p(3,l)**2 )
        gam = 1.0_p_k_part / ( 1.0 + gam )
        gam = gam * ( species%s(1,l)*species%p(1,l) + species%s(2,l)*species%p(2,l) &
          + species%s(3,l)*species%p(3,l) )
        species%s(1,l) = species%s(1,l) + gam * species%p(1,l)
        species%s(2,l) = species%s(2,l) + gam * species%p(2,l)
        species%s(3,l) = species%s(3,l) + gam * species%p(3,l)
      enddo
    endif

  enddo

end subroutine set_spin
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleanup the sdist object
!-----------------------------------------------------------------------------------------
subroutine cleanup_sdist( this )

  implicit none

  type( t_sdist ), intent(inout) :: this

  ! do nothing, placeholder

end subroutine cleanup_sdist

end module m_species_sdist

#endif
