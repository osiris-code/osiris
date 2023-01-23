#include "os-config.h"
#include "os-preprocess.fpp"

module m_cmplx_vdf

#include "memory/memory.h"

use m_system
use m_parameters

use m_vdf_define

implicit none

private

type, extends( t_vdf ) :: t_cmplx_vdf

  ! pointers to 1D, 2D and 3D complex field data
  complex(p_k_fld), dimension(:,:), pointer :: z1 => null()
  complex(p_k_fld), dimension(:,:,:), pointer :: z2 => null()
  complex(p_k_fld), dimension(:,:,:,:), pointer :: z3 => null()

  ! pointer to the complex data buffer
  complex(p_k_fld), dimension(:), pointer :: zbuffer => null()

contains

  procedure :: f_dim => f_dim_cmplx_vdf
  procedure :: new_vdf => new_cmplx_vdf
  procedure :: cleanup => cleanup_cmplx_vdf
  procedure :: read_checkpoint => read_checkpoint_cmplx_vdf

  procedure :: copy_vdf    => copy_cmplx_vdf
  procedure :: copy_fc_vdf => copy_fc_cmplx_vdf

  procedure :: write       => write_real_cmplx_vdf
  procedure :: write_cmplx => write_cmplx_vdf

  procedure :: abs => abs_cmplx_vdf

  procedure :: init_cmplx_pointers

  ! overwrites routines which require to be aware of data layout changes
  procedure :: size_vdf     => size_cmplx_vdf
  procedure :: size_dim_vdf => size_dim_cmplx_vdf

  procedure :: lbound_vdf     => lbound_cmplx_vdf
  procedure :: lbound_dim_vdf => lbound_dim_cmplx_vdf

  procedure :: ubound_vdf     => ubound_cmplx_vdf
  procedure :: ubound_dim_vdf => ubound_dim_cmplx_vdf

  procedure :: check_nan => check_nan_cmplx_vdf

  procedure :: copy_scalar_r4_vdf => copy_r4_cmplx_vdf
  procedure :: copy_scalar_r8_vdf => copy_r8_cmplx_vdf
  procedure :: copy_z4_cmplx_vdf
  procedure :: copy_z8_cmplx_vdf
  ! extends assignment operator
  generic   :: assignment(=) => copy_z4_cmplx_vdf, copy_z8_cmplx_vdf

end type t_cmplx_vdf

integer, parameter :: p_cmplx_real = 0
integer, parameter :: p_cmplx_imag = 1
integer, parameter :: p_cmplx_abs  = 2

public :: t_cmplx_vdf, p_cmplx_real, p_cmplx_imag, p_cmplx_abs

contains

!-------------------------------------------------------------------------------
! Returns number of field dimension in object
!-------------------------------------------------------------------------------
function f_dim_cmplx_vdf( this )

  implicit none

  integer :: f_dim_cmplx_vdf

  class(t_cmplx_vdf), intent(in) :: this

  ! The superclass actually has twice the field dimensions to store
  ! the real and imaginary parts
  if ( this%f_dim_ > 0 ) then
    f_dim_cmplx_vdf = this%f_dim_ / 2
  else
    f_dim_cmplx_vdf = -1
  endif

end function f_dim_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Associates complex data pointers with data buffer
!-------------------------------------------------------------------------------
subroutine init_cmplx_pointers( this )

  use, intrinsic :: iso_c_binding

  implicit none

  class(t_cmplx_vdf), intent(inout) :: this

  type(c_ptr) :: tmp
  integer, dimension(1) :: shape
  integer :: f_dim, i
  integer, dimension( this%x_dim_ + 1) :: lb

  f_dim = this % f_dim_ / 2
  if ( 2 * f_dim /= this % f_dim_ ) then
    ERROR("(*error*) Superclass f_dim is not compatible with complex data")
    call abort_program( p_err_invalid )
  endif

  ! The size of the complex array will have only f_dim components
  shape(1) = f_dim
  do i = 1, this % x_dim_
    shape(1) = shape(1) * ( this % gc_num_(p_lower,i) + this % nx_(i) + &
                            this % gc_num_(p_upper,i) )
  enddo

  ! Associate complex pointers to superclass real pointers
  ! use lowest boundary of `f1, f2, f3` data arrays
  select case ( this % x_dim_ )
  case(1)
    lb = lbound( this % f1 )
    tmp = c_loc( this % f1( lb(1),  lb(2) ) )
    call c_f_pointer( tmp, this % zbuffer, shape )
    this % z1( 1 : f_dim, 1 - this%gc_num_(p_lower,1):this%nx_(1) + this%gc_num_(p_upper,1) ) => this % zbuffer

  case(2)
    lb = lbound( this % f2 )
    tmp = c_loc( this % f2( lb(1), lb(2), lb(3) ) )
    call c_f_pointer( tmp, this % zbuffer, shape )
    this % z2( 1 : f_dim, 1 - this%gc_num_(p_lower,1):this%nx_(1) + this%gc_num_(p_upper,1), &
                          1 - this%gc_num_(p_lower,2):this%nx_(2) + this%gc_num_(p_upper,2) ) => this % zbuffer

  case(3)
    lb = lbound( this % f3 )
    tmp = c_loc( this % f3( lb(1), lb(2), lb(3), lb(4) ) )
    call c_f_pointer( tmp, this % zbuffer, shape )
    this % z3( 1 : f_dim, 1 - this%gc_num_(p_lower,1):this%nx_(1) + this%gc_num_(p_upper,1), &
                          1 - this%gc_num_(p_lower,2):this%nx_(2) + this%gc_num_(p_upper,2), &
                          1 - this%gc_num_(p_lower,3):this%nx_(3) + this%gc_num_(p_upper,3) ) => this % zbuffer

  end select

end subroutine init_cmplx_pointers
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Complex VDF constructor
!-------------------------------------------------------------------------------
subroutine new_cmplx_vdf( this, x_dim, f_dim, nx, gc_num, dx, zero )

  implicit none

  class( t_cmplx_vdf ),   intent( inout ) ::  this
  integer, intent(in) :: x_dim, f_dim
  integer, dimension(:), intent(in) :: nx
  integer, dimension(:,:), intent(in) :: gc_num
  real(p_double), dimension(:), intent(in) :: dx
  logical, intent(in) :: zero

  ! call superclass constructor with twice the field dimensions to
  ! store the real and imaginary parts
  call this % t_vdf % new( x_dim, 2 * f_dim, nx, gc_num, dx, zero )

  ! Associate complex pointers to superclass real pointers
  call this % init_cmplx_pointers()

end subroutine new_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Read checkpoint information for complex vdf
!-------------------------------------------------------------------------------
subroutine read_checkpoint_cmplx_vdf( this, restart_handle )

  use m_restart

  implicit none

  class( t_cmplx_vdf ), intent(inout) ::  this
  type( t_restart_handle ), intent(in) :: restart_handle

  ! Call superclass method to read checkpoint data
  call this % t_vdf % read_checkpoint( restart_handle )

  ! Associate complex pointers to superclass real pointers
  call this % init_cmplx_pointers()

end subroutine read_checkpoint_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Cleanup Complex VDF data
!-------------------------------------------------------------------------------
subroutine cleanup_cmplx_vdf(this)

  implicit none

  class( t_cmplx_vdf ),   intent( inout ) ::  this

  ! Call superclass cleanup
  call this % t_vdf % cleanup()

  ! Nullify complex pointers
  this % z1 => null()
  this % z2 => null()
  this % z3 => null()

  this % zbuffer => null()

end subroutine cleanup_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Writes the real part of the selected field component to disk
! - This method overrides the superclass write method
!-------------------------------------------------------------------------------
subroutine write_real_cmplx_vdf( this, report, fc, g_space, grid, no_co  )

  use m_space
  use m_grid_define
  use m_node_conf
  use m_vdf_report

  implicit none

  class( t_cmplx_vdf ), intent(in) :: this     ! vdf object
  type( t_vdf_report ), intent(in) :: report
  integer,              intent(in) :: fc       ! field component to report
  type( t_space ),      intent(in) :: g_space    ! spatial information
  class( t_grid ),      intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co    ! node configuration of the object

  ! Write the real part of the selected field component
  call this % t_vdf % write( report, 2*(fc-1) + 1, g_space, grid, no_co  )

end subroutine write_real_cmplx_vdf
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Writes the selected field component to disk. The part parameter controls
! whether to write the real part, imaginary part or the magnitude.
!-------------------------------------------------------------------------------
subroutine write_cmplx_vdf( this, report, fc, part, g_space, grid, no_co  )

  use m_space
  use m_grid_define
  use m_node_conf
  use m_vdf_report

  implicit none

  class( t_cmplx_vdf ), intent(in) :: this     ! vdf object
  type( t_vdf_report ), intent(in) :: report
  integer,              intent(in) :: fc       ! field component to report
  integer,              intent(in) :: part     ! part of the complex data to save
  type( t_space ),      intent(in) :: g_space    ! spatial information
  class( t_grid ),       intent(in) :: grid
  class( t_node_conf ),  intent(in) :: no_co    ! node configuration of the object

  ! Temporary vdf to hold magnitude data
  type(t_vdf) :: mag

  select case (part)
  case(p_cmplx_real)
    ! Save the real part
    call this % t_vdf % write( report, 2*(fc-1) + 1, g_space, grid, no_co  )

  case(p_cmplx_imag)
    ! Save the imaginary part
    call this % t_vdf % write( report,  2*(fc-1) + 2, g_space, grid, no_co  )

  case(p_cmplx_abs)
    ! Save the complex magnitude
    call mag % new( this % t_vdf, f_dim = 1 )
    call this % abs( fc, mag )
    call mag % write( report, 1, g_space, grid, no_co  )
    call mag % cleanup()

  end select

end subroutine write_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine abs_cmplx_vdf( this, fc, mag )

  implicit none

  class( t_cmplx_vdf ), intent(in) :: this
  integer, intent(in) :: fc
  type( t_vdf ), intent(inout) :: mag

  integer :: i, j, k

  select case ( this % x_dim_ )
  case(1)
    do i = lbound( this % z1, 2), ubound( this % z1, 2)
      mag % f1(1,i) = abs( this % z1( fc, i ) )
    enddo
  case(2)
    do j = lbound( this % z2, 3), ubound( this % z2, 3)
      do i = lbound( this % z2, 2), ubound( this % z2, 2)
        mag % f2(1,i,j) = abs( this % z2( fc, i,j ) )
      enddo
    enddo
  case(3)
    do k = lbound( this % z3, 4), ubound( this % z3, 4)
      do j = lbound( this % z3, 3), ubound( this % z3, 3)
        do i = lbound( this % z3, 2), ubound( this % z3, 2)
          mag % f3(1,i,j,k) = abs( this % z3( fc, i,j,k) )
        enddo
      enddo
    enddo
  end select

end subroutine abs_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Copies the double precision real value to the vdf
!-------------------------------------------------------------------------------
subroutine copy_r8_cmplx_vdf( this, r8 )

  implicit none

  class(t_cmplx_vdf), intent(inout)  :: this
  real(p_double), intent(in) :: r8

  complex(p_k_fld) :: val

  val = cmplx( r8, kind = p_k_fld )
  this%zbuffer = val

end subroutine copy_r8_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Copies the single precision real value to the vdf
!-------------------------------------------------------------------------------
subroutine copy_r4_cmplx_vdf( this, r4 )

  implicit none

  class(t_cmplx_vdf), intent(inout)  :: this
  real(p_single), intent(in) :: r4

  complex(p_k_fld) :: val

  val = cmplx( r4, kind = p_k_fld )
  this%zbuffer = val

end subroutine copy_r4_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Copies the double precision complex value to the vdf
!-------------------------------------------------------------------------------
subroutine copy_z8_cmplx_vdf( this, z8 )

  implicit none

  class(t_cmplx_vdf), intent(inout)  :: this
  complex(p_double), intent(in) :: z8

  complex(p_k_fld) :: val

  val = cmplx( z8, kind = p_k_fld )
  this%zbuffer = val

end subroutine copy_z8_cmplx_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Copies the single precision complex value to the vdf
!-------------------------------------------------------------------------------
subroutine copy_z4_cmplx_vdf( this, z4 )

  implicit none

  class(t_cmplx_vdf), intent(inout)  :: this
  complex(p_single), intent(in) :: z4

  complex(p_k_fld) :: val

  val = cmplx( z4, kind = p_k_fld )
  this%zbuffer = val

end subroutine copy_z4_cmplx_vdf
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Copies the values of source to this
!-------------------------------------------------------------------------------
subroutine copy_cmplx_vdf( this, source )

  implicit none

  class(t_cmplx_vdf), intent(inout) :: this
  class(t_vdf), intent(in) :: source

  select type (source)
  class is (t_cmplx_vdf)
    ! Source will be a complex data vdf (or a subclass), we can use the
    ! superclass method to do a memcpy of the complete data buffer
    call this % t_vdf % copy( source )
  class is (t_vdf)
    ! If source is a real data vdf we need to convert the data
    this%zbuffer = this%buffer
  end select

end subroutine copy_cmplx_vdf
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Copy the selected field component from source into the 1st field component
! of this
!-------------------------------------------------------------------------------
subroutine copy_fc_cmplx_vdf( this, source, fc )

  implicit none

  class(t_cmplx_vdf), intent(inout) :: this
  class(t_vdf), intent(in) :: source
  integer, intent(in) :: fc

  select type (source)
  class is (t_vdf)
    ! If source is a real data vdf we need to convert the data
    select case (this%x_dim_)
    case(1)
      this%z1(1,:)     = source%f1(fc,:)
    case(2)
      this%z2(1,:,:)   = source%f2(fc,:,:)
    case(3)
      this%z3(1,:,:,:) = source%f3(fc,:,:,:)
    end select

  class is (t_cmplx_vdf)
    ! Source will be a complex data vdf (or a subclass), copy complex data
    select case (this%x_dim_)
    case(1)
      this%z1(1,:)     = source%z1(fc,:)
    case(2)
      this%z2(1,:,:)   = source%z2(fc,:,:)
    case(3)
      this%z3(1,:,:,:) = source%z3(fc,:,:,:)
    end select

  end select

end subroutine copy_fc_cmplx_vdf
!-------------------------------------------------------------------------------

!---------------------------------------------------
function size_cmplx_vdf( this )
!---------------------------------------------------
!---------------------------------------------------

  implicit none

  class( t_cmplx_vdf ), intent(in) :: this

  integer, dimension(this%x_dim_) :: size_cmplx_vdf

  select case (this%x_dim_)
    case(1)
      size_cmplx_vdf(1) = size( this%z1, 2 )
    case(2)
      size_cmplx_vdf(1) = size( this%z2, 2 )
      size_cmplx_vdf(2) = size( this%z2, 3 )
    case(3)
      size_cmplx_vdf(1) = size( this%z3, 2 )
      size_cmplx_vdf(2) = size( this%z3, 3 )
      size_cmplx_vdf(3) = size( this%z3, 4 )
  end select

end function size_cmplx_vdf
!---------------------------------------------------

!---------------------------------------------------
function size_dim_cmplx_vdf( this, dim )
!---------------------------------------------------
!---------------------------------------------------

  implicit none

  integer :: size_dim_cmplx_vdf
  class( t_cmplx_vdf ), intent(in) :: this
  integer, intent(in) :: dim

  select case (this%x_dim_)
    case(1)
      size_dim_cmplx_vdf = size( this%z1, dim )
    case(2)
      size_dim_cmplx_vdf = size( this%z2, dim )
    case(3)
      size_dim_cmplx_vdf = size( this%z3, dim )
    case default
      size_dim_cmplx_vdf = 0
  end select

end function size_dim_cmplx_vdf
!---------------------------------------------------

!---------------------------------------------------
function lbound_cmplx_vdf( this )
!---------------------------------------------------
!---------------------------------------------------

  implicit none

  class( t_cmplx_vdf ), intent(in) :: this
  integer, dimension( this%x_dim_ + 1) :: lbound_cmplx_vdf

  select case (this%x_dim_)
    case(1)
      lbound_cmplx_vdf = lbound( this%z1 )
    case(2)
      lbound_cmplx_vdf = lbound( this%z2 )
    case(3)
      lbound_cmplx_vdf = lbound( this%z3 )
  end select

end function lbound_cmplx_vdf
!---------------------------------------------------


!---------------------------------------------------
function lbound_dim_cmplx_vdf( this, dim )
!---------------------------------------------------
!---------------------------------------------------

  implicit none

  integer :: lbound_dim_cmplx_vdf
  class( t_cmplx_vdf ), intent(in) :: this
  integer, intent(in) :: dim

  select case (this%x_dim_)
    case(1)
      lbound_dim_cmplx_vdf = lbound( this%z1, dim )
    case(2)
      lbound_dim_cmplx_vdf = lbound( this%z2, dim )
    case(3)
      lbound_dim_cmplx_vdf = lbound( this%z3, dim )
    case default
      lbound_dim_cmplx_vdf = 0
  end select

end function lbound_dim_cmplx_vdf
!---------------------------------------------------

!---------------------------------------------------
function ubound_cmplx_vdf( this )
!---------------------------------------------------
!---------------------------------------------------

  implicit none

  class( t_cmplx_vdf ), intent(in) :: this
  integer, dimension( this%x_dim_ + 1 ) :: ubound_cmplx_vdf

  select case (this%x_dim_)
    case(1)
      ubound_cmplx_vdf = ubound( this%z1 )
    case(2)
      ubound_cmplx_vdf = ubound( this%z2 )
    case(3)
      ubound_cmplx_vdf = ubound( this%z3 )
  end select

end function ubound_cmplx_vdf
!---------------------------------------------------


!---------------------------------------------------
function ubound_dim_cmplx_vdf( this, dim )
!---------------------------------------------------
!---------------------------------------------------

  implicit none

  integer :: ubound_dim_cmplx_vdf
  class( t_cmplx_vdf ), intent(in) :: this
  integer, intent(in) :: dim

  select case (this%x_dim_)
    case(1)
      ubound_dim_cmplx_vdf = ubound( this%z1, dim )
    case(2)
      ubound_dim_cmplx_vdf = ubound( this%z2, dim )
    case(3)
      ubound_dim_cmplx_vdf = ubound( this%z3, dim )
    case default
      ubound_dim_cmplx_vdf = 0
  end select

end function ubound_dim_cmplx_vdf
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine check_nan_cmplx_vdf( this, msg )
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

  implicit none

  class( t_cmplx_vdf ), intent(in) :: this
  character( len = * ), intent(in), optional :: msg

  ERROR("(*error*) check_nan_vdf is not implemented for complex valued VDFs")
  call abort_program( p_err_invalid )

  ! TODO: add an implementation here

end subroutine check_nan_cmplx_vdf
!---------------------------------------------------------------------------------------------------

end module m_cmplx_vdf
