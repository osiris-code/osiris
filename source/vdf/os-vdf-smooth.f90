!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     smoothing operation class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_vdf_smooth

#include "memory/memory.h"

  use m_parameters
  use m_utilities
  use m_vdf_define

  use m_math

  implicit none

  ! restrict access to things explicitly declared public
  private

  integer, parameter :: p_none        = 0
  integer, parameter :: p_custom      = 1
  integer, parameter :: p_binomial    = 2
  integer, parameter :: p_compensated = 3
  integer, parameter :: p_5pass       = 4
  integer, parameter :: p_tri6        = 5
  integer, parameter :: p_tri8        = 6

  integer, parameter :: p_laserkdx    = 10
  integer, parameter :: p_digital     = 11

  integer, parameter :: p_order_max = 15


  type :: t_smooth

    private

    ! Type of smoothing in each direction
    integer, dimension(p_x_dim)  ::  type

    ! Order of smoothing (number of neighboring cells included)
    integer, dimension(p_x_dim)  ::  order

    integer :: max_order

    ! Single pass filter kernel for each direction
    ! scoef( 2*maxval(order)+1, p_x_dim )
    real(p_k_fld), dimension(:,:), pointer :: scoef => NULL()

    ! Parameters for lasedkdx filter
    real(p_double) :: laserkdx_omega0       ! frequency for 0 atenuation
    integer        :: laserkdx_dir          ! laser propagation direction

    ! Parameters for digital filter
    real(p_double), dimension(p_x_dim) :: digital_A   ! size of Gibbs phenomenon wiggles in -db
    real(p_double), dimension(p_x_dim) :: digital_fc  ! cutoff frequency

  end type t_smooth


  interface cleanup
    module procedure cleanup_smooth
  end interface

  interface read_nml
    module procedure read_nml_smooth
  end interface

  interface setup
    module procedure setup_smooth
  end interface

  interface smooth
    module procedure smooth_vdf
   end interface

  interface if_smooth
    module procedure if_smooth
  end interface

  interface smooth_order
    module procedure smooth_order
  end interface

  ! declare things that should be public
  public :: t_smooth
  public :: cleanup, setup, read_nml
  public :: smooth, smooth_order, if_smooth

 contains


!-----------------------------------------------------------------------------------------
! Cleanup smooth object
!-----------------------------------------------------------------------------------------
subroutine cleanup_smooth( this )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_smooth ), intent(inout) :: this

  call freemem( this%scoef )

end subroutine cleanup_smooth
!-----------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Digital Filter
! Based on "Digital Filters", Robert Walraven, Proceedings of the Digital Equipment User's Society,
! Fall, 1984.
!---------------------------------------------------------------------------------------------------
function digital_filter( fc, A, n )

  implicit none

  real(p_double), intent(in) :: fc, A
  integer, intent(in) :: n

  real(p_double), dimension( 2*n+1 ) :: digital_filter

  real(p_double) :: alpha
  integer :: i

  if ( A <= 21. ) then
    alpha = 0.
  else if ( A >= 50. ) then
    alpha = 0.1102d0 * ( A - 8.7d0 )
  else
    alpha = 0.5842d0 * (( A - 21.d0) ** 0.4d0) + 0.07886d0*(A-21.d0)
  endif

  do i = 1, n
    digital_filter( n + i + 1 ) = besseli0( alpha*sqrt( 1.0d0 - (real(i)/n)**2) ) / &
                                  besseli0( alpha ) * &
                                  sin( i * pi * fc ) / (i * pi)
  enddo
  digital_filter( n + 1 ) = fc
  do i = 1, n
    digital_filter( n + 1 - i ) = digital_filter( n + 1 + i )
  enddo


end function digital_filter
!---------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Combines n 3 point kernels into a single 2n+1 point kernel
!-----------------------------------------------------------------------------------------
function combine_kernel3( kern3, n, nmax )

  implicit none

  real(p_double), dimension(:,:), intent(in) :: kern3
  integer, intent(in) :: n, nmax

  real(p_double), dimension( 2*nmax+1 ) :: combine_kernel3
  real(p_double) :: f1, f2
  integer :: k, i

  combine_kernel3 = 0.0
  combine_kernel3(nmax+1) = 1.0

  do k = 1, n
    f1 = 0.
    do i = (nmax+1) - (k-1), (nmax+1) + k
      combine_kernel3(i-1) = kern3(3,k) * combine_kernel3(i) + combine_kernel3(i-1)
      f2                   = kern3(1,k) * combine_kernel3(i)
      combine_kernel3(i)   = kern3(2,k) * combine_kernel3(i) + f1
      f1 = f2
    enddo
  enddo

end function combine_kernel3
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Read setup information from input file
!-----------------------------------------------------------------------------------------
subroutine read_nml_smooth( this, input_file, nml_name )
!-----------------------------------------------------------------------------------------

  use m_input_file

  implicit none

  type( t_smooth ), intent(out) :: this
  class( t_input_file ), intent(inout) :: input_file
  character(len = *), optional, intent(in) :: nml_name

  character(len = 20), dimension(p_x_dim)   :: type
  integer, dimension(p_x_dim)               :: order

  ! parameters for custom filter
  real(p_double), dimension(3,p_order_max,p_x_dim) ::  swfj

  ! parameters for laserkdx filter
  real(p_double) :: laserkdx_omega0
  integer        :: laserkdx_dir

  ! parameters for digital filter
  real(p_double), dimension(p_x_dim) :: digital_A
  real(p_double), dimension(p_x_dim) :: digital_fc


  character( len=20 ) :: name
  integer  :: i

  namelist /nl_smooth/ type, order, swfj, &
                       laserkdx_omega0, laserkdx_dir, &
                       digital_A, digital_fc

  integer :: ierr


  type = "none"
  order = 0
  swfj = 1

  laserkdx_omega0 = -1.0
  laserkdx_dir    = 1

  digital_A  = -huge(1.0d0)
  digital_fc = -huge(1.0d0)

  ! Get namelist text from input file
  if( present(nml_name) ) then
    name = trim(nml_name)
    call get_namelist( input_file, "nl_" // trim(name), ierr, "nl_smooth" )
  else
    call get_namelist( input_file, "nl_smooth", ierr )
  endif

  if (ierr == 0) then
    read (input_file%nml_text, nml = nl_smooth, iostat = ierr)
    if (ierr /= 0) then
      print *, "Error reading smooth parameters"
      print *, "aborting..."
      stop
    endif
  else
    if (disp_out(input_file)) then
      SCR_ROOT(" - no smoothing specified")
    endif
  endif

  ! Check for laserkdx filter first, since it overrides the settings in every direction
  if ( trim( type(1) ) == 'laserkdx' ) then
    this%type = p_laserkdx

    if ( (laserkdx_dir < 1) .or. (laserkdx_dir > p_x_dim) ) then
      print *, "Error reading smooth parameters."
      print *, "laserkdx_dir must be between 1 and ", p_x_dim
      print *, 'aborting...'
      stop
    endif
    this%laserkdx_dir = laserkdx_dir
    do i = 1, p_x_dim
      if ( i == laserkdx_dir ) then
        this%order(i) = 2
      else
        this%order(i) = 0
      endif
    enddo

    if (laserkdx_omega0 <= 0.0) then
      print *, "Error reading smooth parameters."
      print *, "Invalid laserkdx_omega0, or laserkdx_omega0 not specified"
      print *, "laserkdx_omega0 must be >= 0."
      print *, 'aborting...'
      stop
    endif
    this%laserkdx_omega0 = laserkdx_omega0
  endif

  ! first pass, check filtering type and get maximum filter order
  this%max_order = 0
  do i = 1, p_x_dim

    select case ( trim( type(i) ) )
       case( 'none' )
          this%order(i) = 0
          this%type(i) = p_none
       case( 'custom' )
          this%type(i) = p_custom
          if ( order(i) < 1 .or. order(i) > p_order_max ) then
             print *, 'Error reading smooth parameters.'
             print *, 'order must be >= 1 and <= p_order_max when using custom smoothing.'
             print *, 'aborting...'
             stop
          endif
          this%order(i) = order(i)

       case( 'binomial' )
          this%type(i) = p_binomial
          if ( order(i) < 1 .or. order(i) > p_order_max ) then
             print *, 'Error reading smooth parameters.'
             print *, 'order must be >= 1 and <= p_order_max when using custom smoothing.'
             print *, 'aborting...'
             stop
          endif
          this%order(i) = order(i)
       case( 'compensated' )
          this%type(i)  = p_compensated
          if ( order(i) < 2 ) then
            this%order(i) = 2
          else
            this%order(i) = order(i)
          endif

       case( '5pass' )
          this%type(i)  = p_5pass
          this%order(i) = 5
       case( 'tri6' )
          this%type(i)  = p_tri6
          this%order(i) = 3
       case( 'tri8' )
          this%type(i)  = p_tri8
          this%order(i) = 4

       case( 'laserkdx' )
          ! this was processed above
          continue

       case( 'digital' )
          this%type(i)  = p_digital
          if ( order(i) < 1 ) then
             print *, 'Error reading smooth parameters.'
             print *, 'order must be >= 1 when using digital smoothing.'
             print *, 'aborting...'
             stop
          endif
          this%order(i) = order(i)

          this%digital_A(i)  = digital_A(i)
          this%digital_fc(i) = digital_fc(i)

       case default
          print *, 'Error reading smooth parameters.'
          print *, 'Invalid smooth type : "', trim( trim( type(i) ) ), '"'
          print *, 'aborting...'
          stop
    end select

    if ( this%order(i) > this%max_order ) this%max_order = this%order(i)
  enddo

  ! allocate the filter coefficients array if needed
  if ( this%max_order > 0 ) then
    call alloc( this%scoef, (/ 2*this%max_order + 1, p_x_dim /) )
    this%scoef = 0.

    ! store custom values
    do i = 1, p_x_dim
      if ( this%type(i) == p_custom ) then
         this%scoef( :, i ) = combine_kernel3( swfj(:,:,i), this%order(i), this%max_order )
      endif
    enddo
  endif

end subroutine read_nml_smooth
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Setup the smooth object
!-----------------------------------------------------------------------------------------
subroutine setup_smooth( this, dx, gamma )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_smooth ), intent( inout )  ::  this
  real( p_double ), dimension(:), intent(in) :: dx
  real( p_double ), intent(in) ::  gamma

  real( p_double ), dimension( 3, p_order_max ) :: kern3
  real( p_double ) :: laserkdx, total, comp

  integer :: i, j

  ! laserkdx is handled separately because it sets filtering in all directions
  if ( this%type(1) == p_laserkdx ) then
    ! adjust coeficients for laserkdx type
    ! from (-1, 6, -1) to (-1, 1/|wcomp|, -1)
    ! wcomp = -W1/[ 1 + 2*W1*(1 + cos[laserkdx]) ]  <= S. F. Martins et al., CPC, 2010
    ! For the binomial smooth, W1 = 1/2 and wcomp = -1/( 4 + 2*cos[laserkdx] )
    ! laserkdx = k*dx/[g(1+beta)]

    laserkdx = this%laserkdx_omega0 * dx(this%laserkdx_dir) / (gamma + sqrt(gamma**2 - 1))

    kern3( :, 1 ) = (/ 1.0d0, 2.0d0, 1.0d0 /)
    kern3( :, 2 ) = (/ -1.0d0,  4 + 2 * cos( laserkdx ), -1.0d0 /)

    this%scoef( :, this%laserkdx_dir ) = combine_kernel3( kern3, 2, this%max_order )

  else

    do i = 1, p_x_dim

      select case ( this%type(i) )
        case ( p_binomial )
          do j = 1, this%order(i)
            kern3( :, j ) = (/ 1.0d0, 2.0d0, 1.0d0 /)
          enddo
          this%scoef( :, i ) = combine_kernel3( kern3, this%order(i), this%max_order )

        case ( p_compensated )
          do j = 1, this%order(i) - 1
            kern3( :, j ) = (/  1.0d0, 2.0d0, 1.0d0 /)
          enddo

          comp = ( 4.0 + 2.0 * ( this%order(i) - 1 ) ) / ( this%order(i) - 1 )

          kern3( :, this%order(i) ) = (/ -1.0d0, comp,-1.0d0 /)

          this%scoef( :, i ) = combine_kernel3( kern3, this%order(i), this%max_order )

        case ( p_5pass )
          kern3( :, 1 ) = (/  1.0d0, 2.0d0, 1.0d0 /)
          kern3( :, 2 ) = (/  1.0d0, 2.0d0, 1.0d0 /)
          kern3( :, 3 ) = (/  1.0d0, 2.0d0, 1.0d0 /)
          kern3( :, 4 ) = (/  1.0d0, 2.0d0, 1.0d0 /)
          kern3( :, 5 ) = (/ -5.0d0, 15.0d0,-5.0d0 /)
          this%scoef( :, i ) = combine_kernel3( kern3, 5, this%max_order )

        case ( p_tri6 )
          this%scoef( (this%max_order+1)-3 : (this%max_order+1)+3 , i ) = &
                               (/ 0.015625, -0.09375, 0.234375, 0.6875, &
                                  0.234375, -0.09375, 0.015625 /)

        case ( p_tri8 )
          this%scoef( (this%max_order+1)-4 : (this%max_order+1)+4 , i ) = &
                               (/ -0.00390625,   0.03125, -0.109375,     0.21875, 0.7265625, &
                                      0.21875, -0.109375,   0.03125, -0.00390625 /)

        case ( p_digital )
          this%scoef( (this%max_order+1)-this%order(i) : (this%max_order+1)+this%order(i) , i ) = &
                           digital_filter( this%digital_fc(i), this%digital_A(i), this%order(i))

      end select
    enddo

  endif

  ! normalize coefficients
  do i = 1, p_x_dim
    if ( this%order(i) > 0 ) then
       total = 0.
       do j = (this%max_order+1) - this%order(i), (this%max_order+1) + this%order(i)
         total = total + this%scoef( j, i )
       enddo

       do j = (this%max_order+1) - this%order(i), (this%max_order+1) + this%order(i)
         this%scoef( j, i ) = this%scoef( j, i )/total
       enddo
    endif
  enddo



end subroutine setup_smooth
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! gives the order of smoothing in each direction
!-----------------------------------------------------------------------------------------
function smooth_order( this )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, dimension(p_x_dim) :: smooth_order

  type( t_smooth ), intent( in )  ::  this

  smooth_order = this%order

end function smooth_order
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Returns true if smoothing in any direction
!-----------------------------------------------------------------------------------------
function if_smooth( this )
!-----------------------------------------------------------------------------------------

  implicit none

  logical :: if_smooth

  type( t_smooth ), intent( in )  ::  this

  integer :: i

  if_smooth = ( this%order(1) > 0 )
  do i=2, p_x_dim
    if_smooth = if_smooth .or. ( this%order(i) > 0 )
  enddo

end function if_smooth
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Filter vdf data
!-----------------------------------------------------------------------------------------
subroutine smooth_vdf( pf, this )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_vdf ),     intent(inout) :: pf
  type( t_smooth ), intent(in)    :: this

  if ( associated(this%scoef) ) then
    select case ( pf%x_dim_ )
     case (1)
       ! call smooth_1d( pf, this%order, this%scoef )
       call filter_1d( pf, this%order, this%max_order, this%scoef )
     case (2)
       ! call smooth_2d( pf, this%order, this%scoef )
       call filter_2d( pf, this%order, this%max_order, this%scoef )
     case (3)
       ! call smooth_3d( pf, this%order, this%scoef )
       call filter_3d( pf, this%order, this%max_order, this%scoef )

     case default
       ERROR("invalid dimensions on pf vdf, xdim = ", pf % x_dim_ )
       call abort_program(p_err_invalid)
    end select
  endif

end subroutine smooth_vdf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Filter 1D vdf
!-----------------------------------------------------------------------------------------
subroutine filter_1d( pf, order, max_order, scoef )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_vdf ), intent(inout)      :: pf
  integer, dimension(:), intent(in) :: order
  integer,  intent(in)              :: max_order
  real(p_k_fld), dimension(:,:), intent(in) :: scoef

  integer :: i1, ext_lb1, ext_ub1, int_lb1, int_ub1
  integer :: nsm1, k

  real(p_k_fld), dimension( pf%f_dim_, -max_order:max_order ) :: o
  real(p_k_fld), dimension( pf%f_dim_ ) :: f
  real(p_k_fld), dimension( -max_order:max_order ) :: s

  ! exterior boundaries of grid
  ext_lb1 = lbound(pf%f1,2)
  ext_ub1 = ubound(pf%f1,2)

  ! interior boundaries of grid
  int_lb1 = ext_lb1 + order(1) ; int_ub1 = ext_ub1 - order(1)

  nsm1 = order(1)

  if ( nsm1 > 0 ) then
    s = scoef(:,1)

    ! use a special version for f_dim = 3
    ! since it is the most used
    if ( pf%f_dim_ == 3 ) then
       do k = -nsm1, nsm1 - 1
         o( 1, k ) = pf%f1(1, int_lb1 + k )
         o( 2, k ) = pf%f1(2, int_lb1 + k )
         o( 3, k ) = pf%f1(3, int_lb1 + k )
       enddo

       do i1 = int_lb1, int_ub1
         o( 1, nsm1 ) = pf%f1(1, i1 + nsm1 )
         o( 2, nsm1 ) = pf%f1(2, i1 + nsm1 )
         o( 3, nsm1 ) = pf%f1(3, i1 + nsm1 )

         f(1) = o( 1, -nsm1 ) * s(-nsm1)
         f(2) = o( 2, -nsm1 ) * s(-nsm1)
         f(3) = o( 3, -nsm1 ) * s(-nsm1)
         do k = -nsm1 + 1, nsm1
           f(1) = f(1) + o(1,k) * s(k)
           f(2) = f(2) + o(2,k) * s(k)
           f(3) = f(3) + o(3,k) * s(k)
         enddo

         do k = -nsm1, nsm1-1
           o( 1, k ) = o(1,k+1)
           o( 2, k ) = o(2,k+1)
           o( 3, k ) = o(3,k+1)
         enddo

         pf%f1(1, i1 ) = f(1)
         pf%f1(2, i1 ) = f(2)
         pf%f1(3, i1 ) = f(3)

       enddo
    else
       ! any number of field components
       do k = -nsm1, nsm1 - 1
         o( :, k ) = pf%f1(:, int_lb1 + k )
       enddo

       do i1 = int_lb1, int_ub1
         o( :, nsm1 ) = pf%f1(:, i1 + nsm1 )

         f(:) = o( :, -nsm1 ) * s(-nsm1)
         do k = -nsm1 + 1, nsm1
           f(:) = f(:) + o(:,k) * s(k)
         enddo

         do k = -nsm1, nsm1-1
           o(:, k) = o(:,k+1)
         enddo

         pf%f1(:, i1 ) = f(:)
       enddo
    endif
  endif

end subroutine filter_1d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Filter 2D vdf
!-----------------------------------------------------------------------------------------
subroutine filter_2d( pf, order, max_order, scoef )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_vdf ), intent(inout)      :: pf
  integer, dimension(:), intent(in) :: order
  integer,  intent(in)              :: max_order
  real(p_k_fld), dimension(:,:), intent(in) :: scoef

  integer :: i1, ext_lb1, ext_ub1, int_lb1, int_ub1
  integer :: i2, ext_lb2, ext_ub2, int_lb2, int_ub2
  integer :: nsm1, nsm2, k

  real(p_k_fld), dimension( pf%f_dim_, -max_order:max_order ) :: o
  real(p_k_fld), dimension( pf%f_dim_ ) :: f
  real(p_k_fld), dimension( -max_order:max_order ) :: s

  ! exterior boundaries of grid
  ext_lb1 = lbound(pf%f2,2)
  ext_ub1 = ubound(pf%f2,2)
  ext_lb2 = lbound(pf%f2,3)
  ext_ub2 = ubound(pf%f2,3)

  ! interior boundaries of grid
  int_lb1 = ext_lb1 + order(1) ; int_ub1 = ext_ub1 - order(1)
  int_lb2 = ext_lb2 + order(2) ; int_ub2 = ext_ub2 - order(2)

  nsm1 = order(1)
  nsm2 = order(2)

  ! use a special version for f_dim = 3
  ! since it is the most used
  if ( pf%f_dim_ == 3 ) then
     if ( nsm1 > 0 ) then

       s = scoef(:,1)

       !$omp parallel do private(k, o, i1, f)
       do i2 = ext_lb2, ext_ub2

         do k = -nsm1, nsm1 - 1
           o( 1, k ) = pf%f2(1, int_lb1 + k, i2 )
           o( 2, k ) = pf%f2(2, int_lb1 + k, i2 )
           o( 3, k ) = pf%f2(3, int_lb1 + k, i2 )
         enddo

         do i1 = int_lb1, int_ub1
           o( 1, nsm1 ) = pf%f2(1, i1 + nsm1, i2 )
           o( 2, nsm1 ) = pf%f2(2, i1 + nsm1, i2 )
           o( 3, nsm1 ) = pf%f2(3, i1 + nsm1, i2 )

           f(1) = o( 1, -nsm1 ) * s(-nsm1)
           f(2) = o( 2, -nsm1 ) * s(-nsm1)
           f(3) = o( 3, -nsm1 ) * s(-nsm1)
           do k = -nsm1 + 1, nsm1
             f(1) = f(1) + o(1,k) * s(k)
             f(2) = f(2) + o(2,k) * s(k)
             f(3) = f(3) + o(3,k) * s(k)
           enddo

           do k = -nsm1, nsm1-1
             o( 1, k ) = o(1,k+1)
             o( 2, k ) = o(2,k+1)
             o( 3, k ) = o(3,k+1)
           enddo

           pf%f2(1, i1, i2 ) = f(1)
           pf%f2(2, i1, i2 ) = f(2)
           pf%f2(3, i1, i2 ) = f(3)

         enddo

       enddo
       !$omp end parallel do

     endif

     if ( nsm2 > 0 ) then

       s = scoef(:,2)

       !$omp parallel do private(k, o, i2, f)
       do i1 = ext_lb1, ext_ub1

         do k = -nsm2, nsm2 - 1
           o( 1, k ) = pf%f2(1, i1, int_lb2 + k )
           o( 2, k ) = pf%f2(2, i1, int_lb2 + k )
           o( 3, k ) = pf%f2(3, i1, int_lb2 + k )
         enddo

         do i2 = int_lb2, int_ub2
           o( 1, nsm2 ) = pf%f2(1, i1, i2 + nsm2 )
           o( 2, nsm2 ) = pf%f2(2, i1, i2 + nsm2 )
           o( 3, nsm2 ) = pf%f2(3, i1, i2 + nsm2 )

           f(1) = o( 1, -nsm2 ) * s(-nsm2)
           f(2) = o( 2, -nsm2 ) * s(-nsm2)
           f(3) = o( 3, -nsm2 ) * s(-nsm2)
           do k = -nsm2 + 1, nsm2
             f(1) = f(1) + o(1,k) * s(k)
             f(2) = f(2) + o(2,k) * s(k)
             f(3) = f(3) + o(3,k) * s(k)
           enddo

           do k = -nsm2, nsm2-1
             o( 1, k ) = o(1,k+1)
             o( 2, k ) = o(2,k+1)
             o( 3, k ) = o(3,k+1)
           enddo

           pf%f2(1, i1, i2 ) = f(1)
           pf%f2(2, i1, i2 ) = f(2)
           pf%f2(3, i1, i2 ) = f(3)

         enddo

       enddo
       !$omp end parallel do

     endif

  else ! any number of field components

     if ( nsm1 > 0 ) then

       s = scoef(:,1)

       !$omp parallel do private(k, o, i1, f)
       do i2 = ext_lb2, ext_ub2

         do k = -nsm1, nsm1 - 1
           o( :, k ) = pf%f2(:, int_lb1 + k, i2 )
         enddo

         do i1 = int_lb1, int_ub1
           o( :, nsm1 ) = pf%f2(:, i1 + nsm1, i2 )

           f(:) = o( :, -nsm1 ) * s(-nsm1)
           do k = -nsm1 + 1, nsm1
             f(:) = f(:) + o(:,k) * s(k)
           enddo

           do k = -nsm1, nsm1-1
             o( :, k ) = o(:,k+1)
           enddo

           pf%f2(:, i1, i2 ) = f(:)
         enddo

       enddo
       !$omp end parallel do

     endif

     if ( nsm2 > 0 ) then

       s = scoef(:,2)

       !$omp parallel do private(k, o, i2, f)
       do i1 = ext_lb1, ext_ub1

         do k = -nsm2, nsm2 - 1
           o( :, k ) = pf%f2(:, i1, int_lb2 + k )
         enddo

         do i2 = int_lb2, int_ub2
           o( :, nsm2 ) = pf%f2(:, i1, i2 + nsm2 )

           f(:) = o( :, -nsm2 ) * s( -nsm2 )
           do k = -nsm2 + 1, nsm2
             f(:) = f(:) + o(:,k) * s(k)
           enddo

           do k = -nsm2, nsm2-1
             o( :, k ) = o(:,k+1)
           enddo

           pf%f2(:, i1, i2 ) = f(:)

         enddo

       enddo
       !$omp end parallel do

     endif

  endif

end subroutine filter_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Filter 3D vdf
!-----------------------------------------------------------------------------------------
subroutine filter_3d( pf, order, max_order, scoef )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_vdf ), intent(inout)      :: pf
  integer, dimension(:), intent(in) :: order
  integer,  intent(in)              :: max_order
  real(p_k_fld), dimension(:,:), intent(in) :: scoef

  integer :: i1, ext_lb1, ext_ub1, int_lb1, int_ub1
  integer :: i2, ext_lb2, ext_ub2, int_lb2, int_ub2
  integer :: i3, ext_lb3, ext_ub3, int_lb3, int_ub3
  integer :: nsm1, nsm2, nsm3, k

  real(p_k_fld), dimension( pf%f_dim_, -max_order:max_order ) :: o
  real(p_k_fld), dimension( pf%f_dim_ ) :: f
  real(p_k_fld), dimension( -max_order:max_order ) :: s

  ! exterior boundaries of grid
  ext_lb1 = lbound(pf%f3,2);   ext_ub1 = ubound(pf%f3,2)
  ext_lb2 = lbound(pf%f3,3);   ext_ub2 = ubound(pf%f3,3)
  ext_lb3 = lbound(pf%f3,4);   ext_ub3 = ubound(pf%f3,4)

  ! interior boundaries of grid
  int_lb1 = ext_lb1 + order(1) ; int_ub1 = ext_ub1 - order(1)
  int_lb2 = ext_lb2 + order(2) ; int_ub2 = ext_ub2 - order(2)
  int_lb3 = ext_lb3 + order(3) ; int_ub3 = ext_ub3 - order(3)

  nsm1 = order(1)
  nsm2 = order(2)
  nsm3 = order(3)

  ! use a special version for f_dim = 3
  ! since it is the most used
  if ( pf%f_dim_ == 3 ) then
     if ( nsm1 > 0 ) then

       s = scoef(:,1)

       !$omp parallel do private(i2, k, o, i1, f)
       do i3 = ext_lb3, ext_ub3
         do i2 = ext_lb2, ext_ub2

           do k = -nsm1, nsm1 - 1
             o( 1, k ) = pf%f3(1, int_lb1 + k, i2, i3 )
             o( 2, k ) = pf%f3(2, int_lb1 + k, i2, i3 )
             o( 3, k ) = pf%f3(3, int_lb1 + k, i2, i3 )
           enddo

           do i1 = int_lb1, int_ub1
             o( 1, nsm1 ) = pf%f3(1, i1 + nsm1, i2, i3 )
             o( 2, nsm1 ) = pf%f3(2, i1 + nsm1, i2, i3 )
             o( 3, nsm1 ) = pf%f3(3, i1 + nsm1, i2, i3 )

             f(1) = o( 1, -nsm1 ) * s(-nsm1)
             f(2) = o( 2, -nsm1 ) * s(-nsm1)
             f(3) = o( 3, -nsm1 ) * s(-nsm1)
             do k = -nsm1 + 1, nsm1
               f(1) = f(1) + o(1,k) * s(k)
               f(2) = f(2) + o(2,k) * s(k)
               f(3) = f(3) + o(3,k) * s(k)
             enddo

             do k = -nsm1, nsm1-1
               o( 1, k ) = o(1,k+1)
               o( 2, k ) = o(2,k+1)
               o( 3, k ) = o(3,k+1)
             enddo

             pf%f3(1, i1, i2, i3 ) = f(1)
             pf%f3(2, i1, i2, i3 ) = f(2)
             pf%f3(3, i1, i2, i3 ) = f(3)

           enddo

         enddo

       enddo
       !$omp end parallel do

     endif

     if ( nsm2 > 0 ) then
       s = scoef(:,2)

       !$omp parallel do private(i1, k, o, i2, f)
       do i3 = ext_lb3, ext_ub3
         do i1 = ext_lb1, ext_ub1

           do k = -nsm2, nsm2 - 1
             o( 1, k ) = pf%f3(1, i1, int_lb2 + k, i3 )
             o( 2, k ) = pf%f3(2, i1, int_lb2 + k, i3 )
             o( 3, k ) = pf%f3(3, i1, int_lb2 + k, i3 )
           enddo

           do i2 = int_lb2, int_ub2
             o( 1, nsm2 ) = pf%f3(1, i1, i2 + nsm2, i3 )
             o( 2, nsm2 ) = pf%f3(2, i1, i2 + nsm2, i3 )
             o( 3, nsm2 ) = pf%f3(3, i1, i2 + nsm2, i3 )

             f(1) = o( 1, -nsm2 ) * s(-nsm2)
             f(2) = o( 2, -nsm2 ) * s(-nsm2)
             f(3) = o( 3, -nsm2 ) * s(-nsm2)
             do k = -nsm2 + 1, nsm2
               f(1) = f(1) + o(1,k) * s(k)
               f(2) = f(2) + o(2,k) * s(k)
               f(3) = f(3) + o(3,k) * s(k)
             enddo

             do k = -nsm2, nsm2-1
               o( 1, k ) = o(1,k+1)
               o( 2, k ) = o(2,k+1)
               o( 3, k ) = o(3,k+1)
             enddo

             pf%f3(1, i1, i2, i3 ) = f(1)
             pf%f3(2, i1, i2, i3 ) = f(2)
             pf%f3(3, i1, i2, i3 ) = f(3)

           enddo

         enddo
       enddo
       !$omp end parallel do

     endif

     if ( nsm3 > 0 ) then
       s = scoef(:,3)

       !$omp parallel do private(i1, k, o, f, i3)
       do i2 = ext_lb2, ext_ub2
         do i1 = ext_lb1, ext_ub1

           do k = -nsm3, nsm3 - 1
             o( 1, k ) = pf%f3(1, i1, i2, int_lb3 + k )
             o( 2, k ) = pf%f3(2, i1, i2, int_lb3 + k )
             o( 3, k ) = pf%f3(3, i1, i2, int_lb3 + k )
           enddo

           do i3 = int_lb3, int_ub3
             o( 1, nsm2 ) = pf%f3(1, i1, i2, i3 + nsm3 )
             o( 2, nsm2 ) = pf%f3(2, i1, i2, i3 + nsm3 )
             o( 3, nsm2 ) = pf%f3(3, i1, i2, i3 + nsm3 )

             f(1) = o( 1, -nsm3 ) * s(-nsm3)
             f(2) = o( 2, -nsm3 ) * s(-nsm3)
             f(3) = o( 3, -nsm3 ) * s(-nsm3)
             do k = -nsm3 + 1, nsm3
               f(1) = f(1) + o(1,k) * s(k)
               f(2) = f(2) + o(2,k) * s(k)
               f(3) = f(3) + o(3,k) * s(k)
             enddo

             do k = -nsm3, nsm3-1
               o( 1, k ) = o(1,k+1)
               o( 2, k ) = o(2,k+1)
               o( 3, k ) = o(3,k+1)
             enddo

             pf%f3(1, i1, i2, i3 ) = f(1)
             pf%f3(2, i1, i2, i3 ) = f(2)
             pf%f3(3, i1, i2, i3 ) = f(3)

           enddo

         enddo
       enddo
       !$omp end parallel do

     endif

  else ! any number of field components

     if ( nsm1 > 0 ) then
       s = scoef(:,1)

       !$omp parallel do private(i2, k, o, i1, f)
       do i3 = ext_lb3, ext_ub3
         do i2 = ext_lb2, ext_ub2

           do k = -nsm1, nsm1 - 1
             o( :, k ) = pf%f3(:, int_lb1 + k, i2, i3 )
           enddo

           do i1 = int_lb1, int_ub1
             o( :, nsm1 ) = pf%f3(:, i1 + nsm1, i2, i3 )

             f(:) = o( :, -nsm1 ) * s(-nsm1)
             do k = -nsm1 + 1, nsm1
               f(:) = f(:) + o(:,k) * s(k)
             enddo

             do k = -nsm1, nsm1-1
               o( :, k ) = o(:,k+1)
             enddo

             pf%f3(:, i1, i2, i3 ) = f(:)
           enddo
         enddo
       enddo
       !$omp end parallel do

     endif

     if ( nsm2 > 0 ) then
       s = scoef(:,2)

       !$omp parallel do private(i1, k, o, i2, f)
       do i3 = ext_lb3, ext_ub3
         do i1 = ext_lb1, ext_ub1

           do k = -nsm2, nsm2 - 1
             o( :, k ) = pf%f3(:, i1, int_lb2 + k, i3 )
           enddo

           do i2 = int_lb2, int_ub2
             o( :, nsm2 ) = pf%f3(:, i1, i2 + nsm2, i3 )

             f(:) = o( :, -nsm2 ) * s(-nsm2)
             do k = -nsm2 + 1, nsm2
               f(:) = f(:) + o(:,k) * s(k)
             enddo

             do k = -nsm2, nsm2-1
               o( :, k ) = o(:,k+1)
             enddo

             pf%f3(:, i1, i2, i3 ) = f(:)
           enddo

         enddo
       enddo
       !$omp end parallel do

     endif

     if ( nsm3 > 0 ) then
       s = scoef(:,3)

       !$omp parallel do private(i1, k, o, f, i3)
       do i2 = ext_lb2, ext_ub2
         do i1 = ext_lb1, ext_ub1

           do k = -nsm3, nsm3 - 1
             o( :, k ) = pf%f3(:, i1, i2, int_lb3 + k )
           enddo

           do i3 = int_lb3, int_ub3
             o( :, nsm2 ) = pf%f3(:, i1, i2, i3 + nsm3 )

             f(:) = o( :, -nsm3 ) * s(-nsm3)
             do k = -nsm3 + 1, nsm3
               f(:) = f(:) + o(:,k) * s(k)
             enddo
             do k = -nsm3, nsm3-1
               o( :, k - 1 ) = o(:,k)
             enddo

             pf%f3(:, i1, i2, i3 ) = f(:)
           enddo

         enddo
       enddo
       !$omp end parallel do

     endif

  endif

end subroutine filter_3d
!-----------------------------------------------------------------------------------------


end module m_vdf_smooth
