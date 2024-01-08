#include "os-config.h"
#include "os-preprocess.fpp"

module m_psource_std

#include "memory/memory.h"

use m_species_define
use m_parameters
use m_fparser
use m_current_define
use m_node_conf

implicit none

private

! density psource types

! Already defined in os-spec-define.f03
!integer, parameter :: p_none      = -1 ! no injection (always return 0)
!integer, parameter :: p_uniform   = 0  ! uniform
integer, parameter :: p_pw_linear = 1  ! piecewise-linear
integer, parameter :: p_gaussian  = 2  ! gaussian
integer, parameter :: p_channel   = 3  ! parabolic channel
integer, parameter :: p_sphere    = 4  ! sphere
integer, parameter :: p_func      = 5  ! math function

! used in cathode only
integer, parameter :: p_step     = 10
integer, parameter :: p_thruster = 11

! maximum number of points for the piecewise-linear profiles
integer, parameter :: p_num_x_max = 100

type, extends(t_psource) :: t_psource_std

  ! Minimum density to inject
  real(p_k_part) :: den_min

  ! type of psource ( defaults to none )
  integer, dimension(p_max_dim) :: type = p_none

  ! flag to define if the psource is the product of independent
  ! functions or not. determined automatically from the input deck.
  logical :: if_mult

  ! Variable for constant charge injection
  integer, dimension(p_max_dim) :: sample_rate = 32


  ! components needed for piecewise linear profiles

  ! number of points in psource ( for piecewise-linear profiles )
  integer :: num_x

  ! allocate later: x,fx( num_x, p_x_dim )
  real(p_k_part), dimension(:,:), pointer :: x  => NULL()
  real(p_k_part), dimension(:,:), pointer :: fx => NULL()

  ! components needed for Gaussian profiles

  ! gauss_center  : beam center in each direction direction
  ! gauss_sigma   : sigma of gaussian distribution
  ! gauss_range   : range over which the Gaussian psource is
  !                 considered to be significant and outside of
  !                 which the distribution is set to zero

  real(p_k_part), dimension(p_max_dim)   :: gauss_center
  real(p_k_part), dimension(p_max_dim)   :: gauss_sigma
  real(p_k_part), dimension(2,p_max_dim) :: gauss_range


  ! density multiplier
  real(p_k_part)               :: density


  ! Components needed parabolic channel profiles

  ! channel_dir       : direction for laser propagation
  ! channel_bottom    : density at the bottom of the channel
  ! channel_r0        : channel parabolic radius
  ! channel_depth     : channel depth
  ! channel_size      : dimension of the parabolic psource
  ! channel_center    : position of the center of the parabolic psource
  ! channel_wall      : channel wall size; if < 0 do a finite channel
  !                     else do a leaky channel with the given wall size
  ! channel_pos(2)    : position of channel entry and channel exit along
  !                     channel_dir

  integer                      :: channel_dir
  real(p_k_part)               :: channel_bottom
  real(p_k_part)               :: channel_r0
  real(p_k_part)               :: channel_depth
  real(p_k_part)               :: channel_size
  real(p_k_part),dimension(p_x_dim)   :: channel_center
  real(p_k_part)                      :: channel_wall
  real(p_k_part), dimension(2)        :: channel_pos

  ! Components needed for shpere profiles
  real(p_k_part), dimension(p_x_dim) :: sphere_center
  real(p_k_part)                     :: sphere_radius

  ! Components needed for math_func profiles
  character(len = p_max_expr_len) :: math_func_expr = " "
  type(t_fparser) :: math_func

contains

  ! Methods from abstract class
  procedure :: if_inject => if_inject_std
  procedure :: read_input => read_input_std
  procedure :: cleanup => cleanup_std
  procedure :: inject => inject_std

  ! Additional methods
  ! Allow access to density values from outside
  procedure :: get_den_value => get_den_value_std
  procedure :: den_value => den_value_std

end type t_psource_std

public :: t_psource_std

public :: p_channel, p_sphere

contains

!-------------------------------------------------------------------------------
! Read information from input file
!-------------------------------------------------------------------------------
subroutine read_input_std( this, input_file, coordinates )
  use m_input_file

  implicit none

  class( t_psource_std ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: coordinates

  ! number of points in piecewise-linear profiles
  integer  :: num_x
  ! components needed for piecewise linear profiles
  real(p_k_part), dimension(p_num_x_max,p_x_dim) :: x, fx

  ! variable used to choose different types of profiles
  character(20),   dimension(p_x_dim) :: profile_type

  ! minimum density to inject particles
  real(p_k_part)                      :: den_min

  ! global density (multiplies the selected density profile)
  real(p_k_part)               :: density

  ! variables needed for Gaussian profiles
  real(p_k_part), dimension(p_x_dim)   :: gauss_center
  real(p_k_part), dimension(p_x_dim)   :: gauss_sigma
  integer, dimension(  p_x_dim )        :: gauss_n_sigma
  real(p_k_part), dimension(2,p_x_dim) :: gauss_range

  ! variables needed for parabolic channel profiles
  integer                      :: channel_dir
  real(p_k_part)               :: channel_bottom
  real(p_k_part)               :: channel_r0
  real(p_k_part)               :: channel_depth
  real(p_k_part)               :: channel_size
  real(p_k_part), dimension(p_x_dim) :: channel_center
  real(p_k_part)               :: channel_wall
  real(p_k_part), dimension(2) :: channel_pos

  ! variables needed for sphere profiles
  real(p_k_part), dimension(p_x_dim) :: sphere_center
  real(p_k_part)                     :: sphere_radius

  ! variables needed for math func profiles
  character(len = p_max_expr_len) :: math_func_expr

  ! This is only used for t_psource_constq but is included here for simplicity
  integer, dimension(p_max_dim) :: sample_rate

  ! namelists of input variables
  namelist /nl_profile/ profile_type, den_min, &
                        density, num_x, x, fx, gauss_n_sigma, gauss_range, &
                        gauss_center, gauss_sigma, &
                        channel_dir, channel_r0, channel_depth, channel_size, &
                        channel_center, channel_wall, channel_pos, channel_bottom, &
                        sphere_center, sphere_radius, math_func_expr, &
                        sample_rate

  integer  ::  i, j
  integer :: ierr


  ! Default profile type (same as uniform)
  profile_type = "default"

  ! Minimum density for injection
  den_min = 0

  !default values for piecewise-linear profile
  num_x = -1
  density = 1.0

  x  = - huge( 1.0_p_k_part )
  fx = 0.0_p_k_part


  ! default values for gaussian psources
  gauss_center  = 0.0_p_k_part
  gauss_sigma   = huge( 1.0_p_k_part )
  gauss_n_sigma = 0
  gauss_range(p_lower,:) = -gauss_sigma
  gauss_range(p_upper,:) =  gauss_sigma

  ! default values for channel parameters
  channel_dir    = 1                          ! default along x1
  channel_r0     = -huge( 1.0_p_k_part )  ! no default r0
  channel_depth  = -huge( 1.0_p_k_part )  ! no default depth
  channel_size   = -huge( 1.0_p_k_part )  ! no default size
  channel_center = -huge( 1.0_p_k_part )  ! no default size
  channel_wall   = 0.0                     ! default wall size
  channel_pos    = -huge( 1.0_p_k_part )  ! no default pos
  channel_bottom = 1.0                     ! default bottom density

  ! default values for sphere parameters
  sphere_center = 0.0
  sphere_radius = 0.0

  ! default values for math function
  math_func_expr = "NO_FUNCTION_SUPPLIED!"

  ! Variables for fixed charge injection - not used in this psource
  sample_rate = 32

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

  ! set the default psource type
  if ( profile_type(1) == "default" ) then
    if (num_x > 0) then
       profile_type(1) = "piecewise-linear"
    else
      profile_type(1) = "uniform"
    endif
  endif

  ! minimum density for injection
  this%den_min          = den_min
  this%density = density

  ! This variable is only used by t_psource_constq but it was simpler to include it here
  ! When using constant charge injection this parameter will be checked in a subclass
  do i = 1, p_x_dim
    this%sample_rate(i) = sample_rate(i)
  enddo


  ! parse non mult options first

  this % if_mult = .false.

  select case ( trim( profile_type(1) ) )
      case ( "piecewise-linear" )
        this % if_mult = .true.
      case ( "gaussian" )
        this % if_mult = .true.
      case ( "channel" )
        this%type = p_channel
        if (p_x_dim == 1) then
           if ( mpi_node() == 0 ) then
              write(0,*) "(*error*) Channel psource is not available in 1D"
              write(0,*) "(*error*) aborting..."
          endif
          stop
        endif

        ! set channel variables
        this%channel_dir = channel_dir
        if ((this%channel_dir < 1) .or. (this%channel_dir > p_x_dim)) then
          if ( mpi_node() == 0 ) then
              write(0,*) "(*error*) Invalid channel direction"
              write(0,*) "(*error*) aborting..."
          endif
          stop
        endif

        this%channel_r0 = channel_r0
        if ( this%channel_r0 <= 0.0 ) then
          if ( mpi_node() == 0 ) then
             write(0,*) "(*error*) Invalid channel_r0 value, must be > 0.0"
             write(0,*) "(*error*) aborting..."
          endif
          stop
        endif

        this%channel_depth  = channel_depth
        if ( this%channel_depth < 0.0 ) then
          if ( mpi_node() == 0 ) then
              write(0,*) "(*error*) Invalid channel_depth value, must be >= 0.0"
              write(0,*) "(*error*) aborting..."
          endif
          stop
        endif

        this%channel_size   = channel_size
        if ( this%channel_size <= 0.0 ) then
          if ( mpi_node() == 0 ) then
             write(0,*) "(*error*) Invalid channel diameter (channel_size) value &
                         &must be > 0.0"
             write(0,*) "(*error*) aborting..."
          endif
          stop
        endif

        do j = 1, p_x_dim-1
          this%channel_center(j) = channel_center(j)
          if (this%channel_center(j) == -huge( 1.0_p_k_part )) then
            if ( mpi_node() == 0 ) then
              write(0,*) "(*error*) channel_center(",j,") was not specified"
              write(0,*) "(*error*) aborting..."
            endif
            stop
          endif
        enddo

        this%channel_wall   = channel_wall
        do j = 1, 2
          this%channel_pos(j) = channel_pos(j)
          if (this%channel_pos(j) == -huge( 1.0_p_k_part )) then
            if ( mpi_node() == 0 ) then
              write(0,*) "(*error*) channel_pos(",j,") was not specified"
              write(0,*) "(*error*) aborting..."
            endif
            stop
          endif
        enddo

        this%channel_bottom = channel_bottom

      case ( "sphere" )
        this%type = p_sphere

        ! set sphere variables
        this%sphere_center = sphere_center
        this%sphere_radius = sphere_radius

      case ( "math func" )
        this%type = p_func

        ! set math func variables
        this%math_func_expr = trim(math_func_expr)

      case ( "uniform" )

        this%type = p_uniform

        ! if not all directions set to uniform use a separable variable function
        do j = 2, p_x_dim
          if ( (trim( profile_type(j) ) /= "uniform") ) this%if_mult = .true.
        enddo

      case default
        if ( mpi_node() == 0 ) then
           write(0,*) '(*error*) Invalid profile_type(1), "', trim(profile_type(1)), '"'
           write(0,*) '(*error*) aborting...'
        endif
        stop
  end select

  ! allocate x, fx if using piecewise-linear
  this%num_x = num_x
  if ((this%if_mult) .and. ( num_x > 0 )) then
    call alloc( this%x, (/  num_x, p_x_dim /) )
    call alloc( this%fx, (/ num_x, p_x_dim /) )
  endif

  if ( this%if_mult ) then
    ! parse mult options
    do j=1, p_x_dim
      if ( profile_type(j) == "default" ) then
        if (num_x > 0) then
           profile_type(j) = "piecewise-linear"
        else
          profile_type(j) = "uniform"
        endif
      endif


      select case ( trim( profile_type(j) ) )
        case ( "uniform" )
          this%type(j) = p_uniform

        case ( "piecewise-linear" )
          if (this%num_x <= 0) then
            if ( mpi_node() == 0 ) then
               write(0,*) "(*error*) Invalid num_x. num_x must be > 0 when using &
                          &'piecewise-linear' density psources"
               write(0,*) "(*error*) aborting..."
            endif
            stop
          endif
          this%type(j) = p_pw_linear

          do i = 1, num_x-1
            if ( x(i,j) >= x(i+1,j) ) then
              if ( mpi_node() == 0 ) then
                 write(0,*) "(*error*) Invalid x value. All values must be larger than predecessor."
                 write(0,*) "(*error*) aborting..."
              endif
              stop
            endif
          enddo

          ! store the piecewise linear parameters for this direction
          do i=1, num_x
            this%x(i,j)  =  x(i,j)
            this%fx(i,j) = fx(i,j)
          enddo

        case ( "gaussian" )
          this%type(j) = p_gaussian

          this%gauss_center(j) = gauss_center(j)
          this%gauss_sigma(j) = abs( gauss_sigma(j) )
          if ( gauss_n_sigma(j) > 0 ) then
             this%gauss_range(p_lower,j) = gauss_center(j)-abs(gauss_n_sigma(j)*gauss_sigma(j))
             this%gauss_range(p_upper,j) = gauss_center(j)+abs(gauss_n_sigma(j)*gauss_sigma(j))
          else
             this%gauss_range(:,j) = gauss_range(:,j)
          endif

        case default
          if ( mpi_node() == 0 ) then
             write(0,*) "(*error*) Invalid profile_type(",j,") : '", trim(profile_type(j)), "'"
             write(0,*) "(*error*) You can only use 'piecewise-linear', 'gaussian', and 'uniform' &
                        &if you are using a separable function psource."
             write(0,*) "(*error*) aborting..."
          endif
          stop
      end select
    enddo

  else
    ! check that no other options were given
    do j = 2, p_x_dim
      if (profile_type(j) /= "default") then
          if ( mpi_node() == 0 ) then
             write(0,*) "(*error*) Invalid profile_type(",j,") : '", trim(profile_type(j)), "'"
             write(0,*) "(*error*) When using : '",trim(profile_type(1)),"' you cannot &
                         &specify any other options for the profile_type"
             write(0,*) "(*error*) aborting..."
          endif
          stop
      endif
    enddo

    ! Check math function
    if ( this%type(1) == p_func ) then
      select case (p_x_dim)
       case (1)
         call setup(this%math_func, trim(this%math_func_expr), &
                       (/'x1'/), ierr)
       case (2)
         call setup(this%math_func, trim(this%math_func_expr), &
                       (/'x1','x2'/), ierr)
       case (3)
         call setup(this%math_func, trim(this%math_func_expr), &
                       (/'x1','x2','x3'/), ierr)
      end select

      if (ierr /= 0) then
         if ( mpi_node() == 0 ) then
            write(0,*) "(*error*) Invalid function supplied : '", &
                        trim(this%math_func_expr), "'"
         endif
         stop
      endif

    endif
  endif


end subroutine read_input_std
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Cleanup t_psource object
!---------------------------------------------------------------------------------------------------
subroutine cleanup_std( this )

  implicit none

  class( t_psource_std ), intent( inout ) :: this

  call freemem( this%x )
  call freemem( this%fx )

  call cleanup( this%math_func )

  this%type = p_none

end subroutine cleanup_std
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Returns .true. if this psource injects particles
!---------------------------------------------------------------------------------------------------
function if_inject_std( this )

  implicit none

  logical :: if_inject_std
  class( t_psource_std ), intent( in ) :: this

  if_inject_std = (this%type(1) /= p_none)

end function if_inject_std
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Gets density values for a set of npxx points
!---------------------------------------------------------------------------------------------------
subroutine get_den_value_std( this, pxx, npxx, den_value )

   implicit none

   ! this needs to be inout because of the eval function
   class( t_psource_std ), intent(inout) :: this
   real(p_k_part), dimension(:,:), intent(in) :: pxx
   integer, intent(in) :: npxx
   real(p_k_part), dimension(:), intent(out) :: den_value

   ! local variables

   integer  ::  i, j, k
   real(p_k_part)  ::  r
   integer  ::  i1, i2
   integer :: x_dim

   ! executable statements
   if (this%type(1) == p_none) then
     den_value = 0.0_p_k_part
     return
   endif

   if (this%if_mult) then ! psource is a product of p_x_dim functions

     den_value = 1.0d0

     do i = 1, p_x_dim

       select case ( this%type(i) )

         case(p_pw_linear) ! piecewise linear -----------------------

            do k = 1, npxx

               if (pxx(i,k) <= this%x(1,i) ) then
                 den_value(k) = den_value(k) * this%fx(1,i)
               elseif (pxx(i,k) >= this%x(this%num_x,i) ) then
                 den_value(k) = den_value(k) * this%fx(this%num_x,i)
               else
                 j = 2
                 do while (pxx(i,k) > this%x(j,i))
                   j = j+1
                 enddo

                 den_value(k) = den_value(k) * ( this%fx(j-1,i) &
                     + ( this%fx(j,i) - this%fx(j-1,i) ) &
                     / ( this%x( j,i) - this%x( j-1,i) ) &
                     * ( pxx(i,k)     - this%x( j-1,i) ) )

               endif

            enddo

         case(p_gaussian) ! gaussian --------------------------------

           do k = 1, npxx

             if ((pxx(i,k) < this%gauss_range(p_lower,i)) .or.  &
                 (pxx(i,k) > this%gauss_range(p_upper,i)) ) then
                   den_value(k) = 0.0
             else
               den_value(k) = den_value(k) * &
                   exp(-(pxx(i,k)-this%gauss_center(i))**2/(2*this%gauss_sigma(i)**2 ))
             endif

           enddo

         case( p_uniform ) ! uniform -----

           ! nothing required

       end select !case ( this%type(i) )

     enddo !i = 1, p_x_dim

   else ! (this%if_mult)

     ! psource is an arbitrary function of pxx

     select case ( this%type(1) )
       case(p_uniform) ! uniform ----------------------------------

          den_value = 1.0d0

       case(p_func) ! math func -----------------------------------

          do k = 1, npxx
            den_value(k) = real( eval( this%math_func, real( pxx(:,k), p_k_fparse ) ), p_k_part )
          enddo

       case (p_sphere)  ! sphere

          do k = 1, npxx
             r = ( pxx(1,k) - this%sphere_center(1) )**2
             do i = 2, p_x_dim
               r = r + ( pxx(i,k) - this%sphere_center(i) )**2
             enddo
             r = sqrt(r)

             if (r <= this%sphere_radius) then
               den_value(k) = 1.0_p_k_part
             else
               den_value(k) = 0.0_p_k_part
             endif
          enddo

       case (p_channel) ! parabolic channel

         if ( p_x_dim == 3 ) then
            i1 = 3 - this%channel_dir
            i2 = 3
            if (this%channel_dir == 3) then
              i1 = 1
              i2 = 2
            endif
         endif

         do k = 1, npxx

           if ((pxx(this%channel_dir,k) >= this%channel_pos(1)) .and. &
               (pxx(this%channel_dir,k) <= this%channel_pos(2))) then

              ! get distance to the center of the channel

              select case (p_x_dim)
                case (1)
                  r = 0.0
                case (2)
                  r = abs( pxx(3 - this%channel_dir,k) - this%channel_center(1))
                case (3)
                  x_dim = p_x_dim-1
                  r = sqrt((pxx(i1,k) - this%channel_center(1))**2 + &
                           (pxx(i2,k) - this%channel_center(x_dim))**2)

              end select

              ! get parabolic psource

              if (r <= this%channel_size/2) then ! inside the channel
                den_value(k) = this%channel_bottom + &
                            this%channel_depth * &
                            (r/this%channel_r0)**2
              else                                 ! outside the channel

                 if (this%channel_wall <= 0.0) then ! finite channel
                   den_value(k) = this%channel_bottom + &
                            this%channel_depth * &
                            ((this%channel_size/2)/this%channel_r0)**2

                 else                              ! leaky channel

                   if (r < (this%channel_size/2 + this%channel_wall )) then
                      den_value(k) = (-r + this%channel_size/2 + this%channel_wall) / &
                                this%channel_wall * (this%channel_bottom + &
                                this%channel_depth * &
                              ((this%channel_size/2)/this%channel_r0)**2)
                   else
                     den_value(k) = 0.0d0
                   endif

                 endif

              endif

           else
             den_value(k) = 0.0d0
           endif

         enddo

       case default ! invalid psource or psource not specified
         den_value(1: npxx) = 1.0d0

     end select

   endif

   den_value(1: npxx) = den_value(1: npxx) * this % density

end subroutine get_den_value_std
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Return density value at specified position
!---------------------------------------------------------------------------------------------------
function den_value_std( this, pxx ) result(den_value)

  implicit none

  class( t_psource_std ), intent(inout) :: this
  real(p_k_part), dimension(:),   intent(in) :: pxx
  real(p_k_part) :: den_value

  real(p_k_part), dimension(p_x_dim,1) ::  lpxx
  real(p_k_part), dimension(1) :: lden


  lpxx(1:p_x_dim,1) = pxx(1:p_x_dim)

  call this%get_den_value( lpxx, 1, lden )
  den_value = lden(1)

end function den_value_std
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
! This routine injects particles into an area defined by grid cell indexes ig_xbnd_inj, using
! variable charge per particle.
!
! It allows for the loading of arbitrary density psources that do not require being a product of
! x_dim functions.
!
! Particles are injected in fixed initial positions inside the simulation cell. The required density
! psource is obtained by varying the individual particle charge.
!---------------------------------------------------------------------------------------------------

function inject_std( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, &
                             send_msg, recv_msg ) result(num_inj)

  use m_species_udist
#ifdef __HAS_SPIN__
  use m_species_sdist
#endif

  implicit none

  class(t_psource_std), intent(inout) :: this
  class(t_species), intent(inout), target :: species
  integer, dimension(:, :), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout)   :: jay
  class( t_node_conf ), intent(in)     :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: num_inj

  integer, parameter :: p_max_x_dim = 3

  ! cell size (p_x_dim)
  real(p_double), dimension(p_max_x_dim) :: dx
  real(p_double), dimension(p_max_x_dim) :: g_xmin

  ! half distance between particles
  real(p_k_part), dimension(p_max_x_dim) :: dxp_2

  integer :: i, i1, i2, i3, ipart

  real(p_k_part)  :: x1, x2

  ! volume that each particle occupies
  real(p_k_part) :: pvol

  ! number of particles per cell
  integer :: ppcell

  ! particle positions (global / inside cell )
  real(p_k_part), dimension(:,:), pointer :: ppos, ppos_cell

  ! particle charges
  real(p_k_part), dimension(:), pointer :: pcharge

  ! Set total number of injected particles to 0
  num_inj = 0

  ! check if num_par_x is > 0 for all directions
  ! if not return silently
  do i=1, p_x_dim
    if (species%num_par_x(i) <= 0) return
  enddo

  do i = 1, p_x_dim
    ! get cell size
    dx(i) = species%dx(i)

    ! get global mininum (this is shifted by +0.5 cells from global simulation values)
    g_xmin( i ) = species%g_box( p_lower , i )

    ! get half distance between particles
    dxp_2(i) = 0.5_p_k_part/species%num_par_x(i)
  enddo

  ! find total number particles per cell
  ppcell = species%num_par_x(1)
  do i = 2, p_x_dim
    ppcell = ppcell * species%num_par_x(i)
  enddo

  ! initialize temp buffers
  call alloc( ppos, (/p_x_dim, ppcell/) )
  call alloc( ppos_cell, (/p_x_dim, ppcell/) )
  call alloc( pcharge, (/ppcell/))

  ! find normalization factor
  pvol = sign( 1.0_p_k_part/ppcell, species%rqm )

  ! Get position of particles inside the cell
  ! The position will always be in the range [-0.5, +0.5 [ regardless of interpolation type

  i = 0
  select case ( p_x_dim )
    case (1)
       do i1 = 1, species%num_par_x(1)
          i = i + 1
          ppos_cell(1,i) = (2*i1 - 1 - species%num_par_x(1)) * dxp_2(1)
       enddo

    case (2)
       do i1 = 1, species%num_par_x(1)
          x1 = (2*i1 - 1 - species%num_par_x(1)) * dxp_2(1)
          do i2 = 1, species%num_par_x(2)
             i = i + 1
             ppos_cell(1,i) = x1
             ppos_cell(2,i) = (2*i2 - 1 - species%num_par_x(2)) * dxp_2(2)
          enddo
       enddo

    case (3)
       do i1 = 1, species%num_par_x(1)
          x1 = (2*i1 - 1 - species%num_par_x(1)) * dxp_2(1)
          do i2 = 1, species%num_par_x(2)
             x2 = (2*i2 - 1 - species%num_par_x(2)) * dxp_2(2)
             do i3 = 1, species%num_par_x(3)
                i = i + 1
                ppos_cell(1,i) = x1
                ppos_cell(2,i) = x2
                ppos_cell(3,i) = (2*i3 - 1 - species%num_par_x(3)) * dxp_2(3)
             enddo
          enddo
       enddo

  end select

  ! Inject particles
  ipart = species%num_par + 1

  select case ( p_x_dim )

     case (1) ! -----------------------------------------------------------

       ! loop through all cells and count number of particles to inject
       do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)
         ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )

         call this%get_den_value( ppos, ppcell, pcharge )

         do i = 1, ppcell
           if (pcharge(i) > this % den_min ) then
             num_inj = num_inj + 1
           endif
         enddo
       enddo

       if ( num_inj > 0 ) then

          ! check if the buffer size is sufficient and grow it if necessary
          if ( num_inj > species%num_par_max - species%num_par ) then
             call species%grow_buffer( species%num_par_max + num_inj + &
                                              p_spec_buf_block )
          endif

          ! loop through all the injection cells and
          ! inject particles, normalizing charge

          do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)
            ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )

            call this%get_den_value( ppos, ppcell, pcharge )

            do i=1, ppcell
              if (pcharge(i) > this % den_min) then
                  ! add particle
                  species%q(ipart)  = pcharge(i) * pvol
                  species%x(1, ipart) = ppos_cell(1,i)
                  species%ix(1, ipart) = i1 - species%my_nx_p( p_lower, 1 ) + 1
                  ipart = ipart + 1
              endif

            enddo
          enddo

       endif

     case (2) ! -----------------------------------------------------------


       ! loop through all cells and count number of particles to inject
       select case ( species%coordinates )

         case default  ! cartesian coordinates
            do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)
              ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )

              do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2)
                ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )

                call this%get_den_value( ppos, ppcell, pcharge )

                do i = 1, ppcell
                  if (pcharge(i) > this % den_min ) then
                    num_inj = num_inj + 1
                  endif
                enddo

              enddo
            enddo


         case ( p_cylindrical_b ) ! cyl. coord. -> B1 on axis
                                  ! don't inject on axis

            do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)
              ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )

              do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2)
                ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )

                call this%get_den_value( ppos, ppcell, pcharge )

                do i = 1, ppcell
                  if ((pcharge(i) > this % den_min) .and. &
                       (ppos(p_r_dim,i) > 0.0_p_k_part)) then

                    num_inj = num_inj + 1
                  endif
                enddo

              enddo
            enddo
       end select

       if ( num_inj > 0 ) then

          ! check if the buffer size is sufficient and grow it if necessary
          if ( num_inj > species%num_par_max - species%num_par ) then
             call species%grow_buffer( species%num_par_max + num_inj + p_spec_buf_block )
          endif


          ! loop through all the injection cells and
          ! inject particles, normalizing charge
          select case ( species%coordinates )
            case default  ! cartesian coordinates

              do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)
                ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )

                do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2)
                  ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )

                  call this%get_den_value( ppos, ppcell, pcharge )

                  do i=1, ppcell
                    if (pcharge(i) > this % den_min) then
                        ! add particle
                        species%x(1, ipart) = ppos_cell(1,i)
                        species%x(2, ipart) = ppos_cell(2,i)
                        species%ix(1, ipart) = i1 - species%my_nx_p( p_lower, 1 ) + 1
                        species%ix(2, ipart) = i2 - species%my_nx_p( p_lower, 2 ) + 1
                        species%q(ipart)  = pcharge(i) * pvol
                        ipart = ipart + 1
                     endif

                  enddo

                enddo
              enddo

            case ( p_cylindrical_b ) ! cyl. coord. -> B1 on axis
                                     ! don't inject on axis

              do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)
                ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )

                do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2)
                  ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )

                  call this%get_den_value( ppos, ppcell, pcharge )

                  do i=1, ppcell
                     if ((pcharge(i) > this % den_min) .and. &
                         (ppos(p_r_dim,i) > 0.0_p_k_part)) then
                        ! add particle
                        species%x(1, ipart) = ppos_cell(1,i)
                        species%x(2, ipart) = ppos_cell(2,i)
                        species%ix(1, ipart) = i1 - species%my_nx_p( p_lower, 1 ) + 1
                        species%ix(2, ipart) = i2 - species%my_nx_p( p_lower, 2 ) + 1
                        species%q(ipart)  = pcharge(i) * pvol * ppos( p_r_dim, i )
                        ipart = ipart + 1
                     endif
                  enddo
                enddo
              enddo

          end select

       endif


     case (3) ! -----------------------------------------------------------

       ! loop through all cells and count number of particles to inject
       do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)
         ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )

         do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2)
           ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )

           do i3 = ig_xbnd_inj(p_lower,3), ig_xbnd_inj(p_upper,3)
              ppos(3,:) = real( g_xmin(3) + (ppos_cell(3,:) + (i3-1))*dx(3), p_k_part )

              call this%get_den_value( ppos, ppcell, pcharge )

              do i = 1, ppcell
                if (pcharge(i) > this % den_min ) then
                  num_inj = num_inj + 1
                endif
              enddo
           enddo
         enddo
       enddo

       if ( num_inj > 0 ) then

          ! check if the buffer size is sufficient and grow it if necessary
          if ( num_inj > species%num_par_max - species%num_par ) then
             call species%grow_buffer( species%num_par_max + num_inj + p_spec_buf_block )
          endif

          ! loop through all the injection cells and
          ! inject particles

          do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)
            ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )

            do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2)
              ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )

              do i3 = ig_xbnd_inj(p_lower,3), ig_xbnd_inj(p_upper,3)
                ppos(3,:) = real( g_xmin(3) + (ppos_cell(3,:) + (i3-1))*dx(3), p_k_part )

                call this%get_den_value( ppos, ppcell, pcharge )

                do i=1, ppcell
                  if (pcharge(i) > this % den_min) then
                         ! add particle
                         species%q(ipart)  = pcharge(i) * pvol
                         species%x(1, ipart) = ppos_cell(1,i)
                         species%x(2, ipart) = ppos_cell(2,i)
                         species%x(3, ipart) = ppos_cell(3,i)
                         species%ix(1, ipart) = i1 - species%my_nx_p( p_lower, 1 ) + 1
                         species%ix(2, ipart) = i2 - species%my_nx_p( p_lower, 2 ) + 1
                         species%ix(3, ipart) = i3 - species%my_nx_p( p_lower, 3 ) + 1
                         ipart = ipart + 1
                   endif

                enddo
              enddo
            enddo
          enddo

       endif

  end select

  ! free temporary memory
  call freemem( ppos )
  call freemem( ppos_cell )
  call freemem( pcharge )

  ! Set momentum of injected particles
  if ( num_inj > 0 ) then
    ! Initialize particle momentum
    call set_momentum( species, species%num_par+1, species%num_par + num_inj )
#ifdef __HAS_SPIN__
    call set_spin( species, species%num_par+1, species%num_par + num_inj )
#endif
  endif

end function inject_std
!-------------------------------------------------------------------------------



end module m_psource_std
