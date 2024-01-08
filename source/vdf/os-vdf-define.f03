!#define DEBUG_FILE 1
!-------------------------------------------------------------------------------
! VDF definition module
!
! This file contains the class definition for the following classes:
!
! t_vdf
! t_vdf_msg
!-------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_vdf_define

#include "memory/memory.h"

use m_system
use m_parameters

use m_diagnostic_utilities, only : p_diag_prec, t_diag_file
use m_diagfile, only : p_no_offset
use m_space, only : t_space
use m_grid_define, only : t_grid
use m_node_conf, only : t_node_conf

private

integer, parameter :: p_hc_x = 1
integer, parameter :: p_hc_y = 2
integer, parameter :: p_hc_z = 4

character(len=*), parameter, private :: p_vdf_rst_id = "vdf rst data - 0x0001"

type :: t_vdf

  ! pointers to 1D, 2D and 3D field data
  real(p_k_fld), dimension(:,:), pointer :: f1 => null()
  real(p_k_fld), dimension(:,:,:), pointer :: f2 => null()
  real(p_k_fld), dimension(:,:,:,:), pointer :: f3 => null()

  ! pointer to the data buffer
  real(p_k_fld), dimension(:), pointer :: buffer => null()

  ! spacial dimensions
  integer :: x_dim_ = -1

  ! field dimensions
  integer :: f_dim_ = -1

  ! grid size
  integer, dimension(p_max_dim) :: nx_ = 0

  ! guard cell information
  integer, dimension(2, p_max_dim) :: gc_num_ = 0

  ! cell size information
  ! this should not be needed, temporary fix
  real(p_double), dimension(p_max_dim) :: dx_ = 0.0_p_double

contains

  procedure :: f_dim => f_dim_vdf

  procedure :: new_vdf
  procedure :: new_copy_vdf
  generic   :: new => new_vdf, new_copy_vdf
  procedure :: cleanup => cleanup_vdf
  procedure :: zero => zero_vdf

  procedure :: write_checkpoint => write_checkpoint_vdf
  procedure :: read_checkpoint => read_checkpoint_vdf

  procedure :: copy_vdf
  procedure :: copy_fc_vdf
  procedure :: copy_component_vdf
  generic   :: copy => copy_vdf, copy_fc_vdf, copy_component_vdf

  ! this type bound procedure is in a separate file,
  ! only the interface is declared here
  procedure :: write => write_vdf

  ! type bound procedures querying VDF informations
  procedure :: dvol => dvol_vdf

  procedure :: size_vdf
  procedure :: size_dim_vdf
  generic   :: size => size_vdf, size_dim_vdf

  procedure :: offset_vdf
  procedure :: offset_dim_vdf
  generic   :: offset => offset_vdf, offset_dim_vdf

  procedure :: gc_num_vdf
  procedure :: gc_num_bnd_vdf
  generic   :: gc_num => gc_num_vdf, gc_num_bnd_vdf

  procedure :: dx_vdf
  procedure :: dx_dim_vdf
  generic   :: dx => dx_vdf, dx_dim_vdf

  procedure :: lbound_vdf
  procedure :: lbound_dim_vdf
  generic   :: lbound => lbound_vdf, lbound_dim_vdf

  procedure :: ubound_vdf
  procedure :: ubound_dim_vdf
  generic   :: ubound => ubound_vdf, ubound_dim_vdf

  procedure :: check_nan => check_nan_vdf

  procedure :: copy_scalar_r4_vdf
  procedure :: copy_scalar_r8_vdf
  generic   :: assignment(=) => copy_scalar_r4_vdf, copy_scalar_r8_vdf, copy_vdf

end type t_vdf

! this extra type necessary to make an array of pointers in fortran
type :: t_vdf_arr

  type( t_vdf ), pointer :: v => null()

end type t_vdf_arr

! ---------------------------- report ---------------------------------

integer, parameter :: p_n_report_type = 5, &
                      p_full  = 1, &
                      p_savg  = 2, &
                      p_senv  = 3, &
                      p_line  = 4, &
                      p_slice = 5


integer, parameter :: p_max_reports = 1024
integer, parameter :: p_max_reports_len = 32

type :: t_vdf_report_item

  ! type of report: full, spatial average, spatial envelope, line, slice
  integer :: type
  ! use time averaged data
  logical :: tavg
  character(len=64)  :: name

  ! line and slice parameters
  integer :: id
  integer :: direction
  integer, dimension(2) :: gipos

  type(t_vdf_report_item), pointer :: next => null()

end type

! structure describing vdf reports of a given quantity
! this structure does not change during the simulation
type :: t_vdf_report

   ! these are set just before saving the file
   character(len=1024) :: path = ''
   character(len=1024) :: filename = ''
   integer            :: n = 0
   real(p_double)     :: t = 0.0d0

   integer :: quant ! integer id of quantity to report

   integer, dimension(p_n_report_type) :: ndump

   character(len=256) :: fileLabel = ''
   character(len=1024) :: basePath  = ''

   character(len=64)  :: name       ! this should be a suitable variable name e.g. E1 (spaces ok)
   character(len=64)  :: label      ! this is a label to use in plots e.g. "E_1^2" (LaTeX ok)
   character(len=64)  :: units

   ! Setting default values for these causes an "internal compiler error" on PGI
   ! on os-spec-diagnostics
   character(len=64), dimension(3) :: xname
   character(len=64), dimension(3) :: xlabel
   character(len=64), dimension(3) :: xunits

   ! fraction of cells offset from lower-left corner
   real(p_double), dimension(3) :: offset_x = p_no_offset
   ! fraction of time step offset
   real(p_double) :: offset_t = p_no_offset

   character(len=64)  :: time_units = ''
   real(p_double)     :: dt = 0.0d0

   ! precision of diagnostics
   integer :: prec = p_diag_prec

   ! spatial average/envelope parameters
   integer, dimension(3) :: n_ave

   ! number of line/slices in each direction
   integer, dimension(3) :: line_idx = 0
   integer, dimension(3) :: slice_idx = 0

   ! time average data
   integer :: n_tavg = -1         ! number of timesteps to average before report
   integer :: ndump_tavg = -1     ! the minimum number of timesteps between time averaged
                                  ! data dumps
   type( t_vdf ) :: tavg_data     ! vdf holding time averaged data

   ! list of individual report items
   type(t_vdf_report_item), pointer :: list => null()

   ! next vdf report
   type(t_vdf_report), pointer :: next => null()

end type t_vdf_report

! ---------------------------------------------------------------------


interface
  subroutine write_vdf( this, report, fc, g_space, grid, no_co  )

  import t_vdf, t_vdf_report, t_space, t_grid, t_node_conf

  class( t_vdf ),       intent(in) :: this
  type( t_vdf_report ), intent(in) :: report
  integer,              intent(in) :: fc
  type( t_space ),      intent(in) :: g_space
  class( t_grid ),      intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co

  end subroutine
end interface

! define exports for use statements
public :: t_vdf, t_vdf_arr, t_vdf_report_item, t_vdf_report
public :: p_hc_x, p_hc_y, p_hc_z
public :: p_n_report_type, p_full, p_savg, p_senv, p_line, p_slice
public :: p_max_reports, p_max_reports_len

contains

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
function f_dim_vdf( this )

  implicit none

  integer :: f_dim_vdf

  class(t_vdf),intent(in) :: this

  f_dim_vdf = this%f_dim_

end function f_dim_vdf
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine cleanup_vdf(this)

  implicit none

  class( t_vdf ),   intent( inout ) ::  this

  ! calling freemem with unassociated pointers is ok
  select case (this%x_dim_)
    case(1)
      call freemem( this%f1 )
    case(2)
      call freemem( this%f2 )
    case(3)
      call freemem( this%f3 )
    case default
      continue
  end select

  this % f1 => null()
  this % f2 => null()
  this % f3 => null()

  this % buffer => null()

  this%x_dim_ = -1
  this%f_dim_ = -1

  this%nx_ = 0
  this%gc_num_ = 0
  this%dx_ = 0.0_p_double

end subroutine cleanup_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vdf constructor
!-------------------------------------------------------------------------------
subroutine new_vdf( this, x_dim, f_dim, nx, gc_num, dx, zero )

  use, intrinsic :: iso_c_binding

  implicit none

  class( t_vdf ),   intent( inout ) ::  this
  integer, intent(in) :: x_dim, f_dim
  integer, dimension(:), intent(in) :: nx
  integer, dimension(:,:), intent(in) :: gc_num
  real(p_double), dimension(:), intent(in) :: dx
  logical, intent(in) :: zero

  integer :: i
  integer, dimension(1+p_max_dim) :: lb, ub
  type(c_ptr) :: tmp
  integer, dimension(1) :: shape

  ! check if supplied parameters are ok
  if ( (x_dim < 1) .or. (x_dim > p_max_dim)) then
    ERROR("Invalid value for x_dim, ", x_dim, " must be between 1 and ", p_max_dim)
    call abort_program( p_err_invalid )
  endif

  if (f_dim < 1) then
    ERROR("Invalid value for f_dim, ", f_dim, " must be >= 1")
    call abort_program( p_err_invalid )
  endif

  if ( size(nx) < x_dim ) then
    ERROR("The size of nx, ",size(nx),"is smaller that x_dim, ",x_dim)
    call abort_program( p_err_invalid )
  endif

  if (( size(gc_num,1) /= 2 ) .or. ( size(gc_num,2) < x_dim )) then
    ERROR("The size of gc_num is invalid, must be gc_num(2, >= x_dim)")
    call abort_program( p_err_invalid )
  endif

  ! Everything ok, create new vdf object
  call this % cleanup()

  this%x_dim_ = x_dim
  this%f_dim_ = f_dim

  this%nx_(1:x_dim) = nx(1:x_dim)
  this%gc_num_(:, 1:x_dim) = gc_num(:, 1:x_dim)

  this%dx_(1:x_dim) = dx(1:x_dim)

  lb(1) = 1
  ub(1) = this%f_dim_

  do i = 1, x_dim
    lb(i+1) =  1-this%gc_num_(p_lower,i)
    ub(i+1) =  this%nx_(i)+this%gc_num_(p_upper,i)
  enddo

  ! c_loc needs to point to first value in array (IBM compiler)
  select case (x_dim)
  case(1)
    call alloc( this%f1, lb, ub )
    tmp = c_loc( this % f1( lb(1), lb(2) ) )

  case(2)
    call alloc( this%f2, lb, ub )
    tmp = c_loc( this % f2( lb(1), lb(2), lb(3) ) )

  case(3)
    call alloc( this%f3, lb, ub )
    tmp = c_loc( this % f3( lb(1), lb(2), lb(3), lb(4) ) )

  case default
    tmp = c_null_ptr

  end select

  ! Associate buffer pointer with allocated data
  shape(1) = f_dim
  do i = 1, x_dim
    shape(1) = shape(1) * (ub(i+1) - lb(i+1) + 1)
  enddo
  call c_f_pointer( tmp, this % buffer, shape )

  if (zero) then
    ! Set the grid values to 0
    call this%zero()
  endif

end subroutine new_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vdf copy constructor
! Creates a new vdf object with the same structure as the source.
! source    - Souce VDF to copy structure (and data) from
! f_dim     - Number of field components of the new vdf. Defaults to the same
!             number of components as the source
! zero      - Controls whether to zero the new vdf. Defaults to .false.
! copy      - Controls whether to copy the vdf data from source. Defaults to
!             .false.
! fc        - Controls which field component to copy from source. Defaults to
!             copying all components.
!-------------------------------------------------------------------------------
subroutine new_copy_vdf( this, source, f_dim, zero, copy, fc )

  implicit none

  class( t_vdf ),   intent( inout ) ::  this
  class( t_vdf ),   intent( in ) ::  source
  integer, intent(in), optional :: f_dim
  logical, intent(in), optional :: zero
  logical, intent(in), optional :: copy
  integer, intent(in), optional :: fc

  integer :: f_dim_
  logical :: zero_, copy_

  if (present(f_dim)) then
    if (f_dim < 1) then
      ERROR("Field components (f_dim) must be >= 1")
      call abort_program(p_err_invalid)
    endif
    f_dim_ = f_dim
  else
    f_dim_ = source%f_dim()
  endif

  if ( present(zero) ) then
    zero_ = zero
  else
    zero_ = .false.
  endif

  if ( present( copy ) ) then
    copy_ = copy
    if ( copy_ .and. (.not. present( fc )) ) then
      if ( source%f_dim() /= f_dim_) then
        write(0,*) "(*error*) vdf%new() was called with copy = .true. but target and source"
        write(0,*) "(*error*) field dimensions do not match (may need to set the fc parameter)."
        call abort_program(p_err_invalid)
      endif
    endif
  else
    copy_ = .false.
  endif

  ! create new vdf
  call this%new( source%x_dim_, f_dim_, source%nx_, source%gc_num_, source%dx_, zero_ )

  if ( copy_ ) then
    if ( zero_ ) &
      write(0,*) "(*warning*) vdf%new() called with both copy and zero set to true"

    ! copy values to new vdf
    if (.not. present(fc)) then
      ! Copy all data from source
      call this % copy( source )
    else
      ! If target vdf has more than 1 field component zero the array first
      if ( f_dim_ > 1 ) call this % zero()

      ! Copy only the selected field component into the 1st component
      call this % copy( source, fc )
    endif
  else
    ! otherwise zero the array
    call this % zero()
  endif

end subroutine new_copy_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Sets the value of the vdf to zero using memset
!-------------------------------------------------------------------------------
subroutine zero_vdf( this )

  use iso_c_binding
  use m_system

  implicit none

  class(t_vdf), intent(inout)  :: this

  real(p_k_fld) :: tmp
  type(c_ptr)   :: b
  integer(c_size_t) :: len

  b = c_loc( this%buffer(1) )
  len = c_sizeof( tmp ) * size( this%buffer )

  call zero( b, len )

end subroutine zero_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Write checkpoint information for vdf
!-------------------------------------------------------------------------------
subroutine write_checkpoint_vdf( this, restart_handle )

  use m_restart

  implicit none

  class( t_vdf ),   intent(in) ::  this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for vdf object.'
  integer :: ierr

  ! Save version id
  restart_io_wr( p_vdf_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! save number of dimensions and grid parameters
  restart_io_wr( this%x_dim_, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%f_dim_, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%nx_, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%gc_num_, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%dx_, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  if (this%x_dim_ > 0) then
    select case (this%x_dim_ )
      case(1)
      restart_io_wr( this%f1, restart_handle, ierr )
      case(2)
      restart_io_wr( this%f2, restart_handle, ierr )
      case(3)
      restart_io_wr( this%f3, restart_handle, ierr )
    end select
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
  endif

end subroutine write_checkpoint_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! read checkpoint information for vdf
!-------------------------------------------------------------------------------
subroutine read_checkpoint_vdf( this, restart_handle )

  use, intrinsic :: iso_c_binding

  use m_restart

  implicit none

  class( t_vdf ), intent(inout) ::  this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for vdf object.'
  character(len=len(p_vdf_rst_id)) :: rst_id
  integer, dimension(1+p_max_dim) :: lb, ub
  integer :: i, ierr

  type(c_ptr) :: tmp
  integer, dimension(1) :: shape

  ! clear vdf data and free memory if necessary
  DEBUG('before this%cleanup')
  call this%cleanup()

  ! check if restart file is compatible
  DEBUG('before reading restart')
  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  if ( rst_id /= p_vdf_rst_id) then
    ERROR('Corrupted restart file, or restart file ')
    ERROR('from incompatible binary (vdf)')
    call abort_program(p_err_rstrd)
  endif

  ! read number of dimensions and grid parameters
  DEBUG('before reading this%x_dim_')
  restart_io_rd( this%x_dim_, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  DEBUG('before reading this%f_dim_')
  restart_io_rd( this%f_dim_, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  DEBUG('before reading this%nx_')
  restart_io_rd( this%nx_, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  DEBUG('before reading this%gc_num_')
  restart_io_rd( this%gc_num_, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  DEBUG('before reading this%dx_')
  restart_io_rd( this%dx_, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  ! read vdf data if necessary
  if (this%x_dim_ > 0) then
    DEBUG('reading vdf data ...')

    lb(1) = 1
    ub(1) = this%f_dim_

    do i = 1, this%x_dim_
      lb(i+1) =  1-this%gc_num_(p_lower,i)
      ub(i+1) =  this%nx_(i)+this%gc_num_(p_upper,i)
    enddo

    ! c_loc needs to point to first value in array (IBM compiler)
    select case ( this%x_dim_ )
    case(1)
      call alloc( this%f1, lb, ub )
      tmp = c_loc( this % f1( lb(1), lb(2) ) )
      restart_io_rd( this%f1, restart_handle, ierr )

    case(2)
      call alloc( this%f2, lb, ub )
      tmp = c_loc( this % f2( lb(1), lb(2), lb(3) ) )
      restart_io_rd( this%f2, restart_handle, ierr )

    case(3)
      call alloc( this%f3, lb, ub )
      tmp = c_loc( this % f3( lb(1), lb(2), lb(3), lb(4) ) )
      restart_io_rd( this%f3, restart_handle, ierr )

    case default
      tmp = c_null_ptr

    end select

    CHECK_ERROR( ierr, err_msg, p_err_rstrd )

    ! Associate buffer pointer with allocated data
    shape(1) = this%f_dim_
    do i = 1, this%x_dim_
      shape(1) = shape(1) * (ub(i+1) - lb(i+1) + 1)
    enddo
    DEBUG('before c_f_pointer')
    call c_f_pointer( tmp, this % buffer, shape )

  endif

  DEBUG('end of reading vdf')

end subroutine read_checkpoint_vdf
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Copies the values of source to this
!-------------------------------------------------------------------------------
subroutine copy_vdf( this, source )

  implicit none

  class(t_vdf), intent(inout) :: this
  class(t_vdf), intent(in) :: source

  call memcpy( this%buffer, source%buffer, size(this%buffer))

end subroutine copy_vdf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Copy the selected field component from source into the 1st field component
! of this
!-------------------------------------------------------------------------------
subroutine copy_fc_vdf( this, source, fc )

  implicit none

  class(t_vdf), intent(inout) :: this
  class(t_vdf), intent(in) :: source
  integer, intent(in) :: fc

  select case (this%x_dim_)
  case(1)
    this%f1(1,:)     = source%f1(fc,:)
  case(2)
    this%f2(1,:,:)   = source%f2(fc,:,:)
  case(3)
    this%f3(1,:,:,:) = source%f3(fc,:,:,:)
  end select

end subroutine
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Copy source_component from source into dest_component of this 
!-------------------------------------------------------------------------------
subroutine copy_component_vdf( this, source, source_component, dest_component )

  implicit none

  class(t_vdf), intent(inout) :: this
  class(t_vdf), intent(in) :: source
  integer, intent(in) :: source_component, dest_component

  select case (this%x_dim_)
  case(1)
    this%f1(dest_component,:)     = source%f1(source_component,:)
  case(2)
    this%f2(dest_component,:,:)   = source%f2(source_component,:,:)
  case(3)
    this%f3(dest_component,:,:,:) = source%f3(source_component,:,:,:)
  end select

end subroutine
!-------------------------------------------------------------------------------


!---------------------------------------------------
 function dvol_vdf( this )
!---------------------------------------------------
!       gives the volume of the grid cells
!       (deprecated)
!---------------------------------------------------

   implicit none

!       dummy variables

   real(p_double) :: dvol_vdf

   class(t_vdf), intent(in) :: this

!       local variables

   integer :: i

!       executable statements

   if (this%x_dim_ > 0) then
     dvol_vdf = this%dx_(1)
     do i=2, this%x_dim_
       dvol_vdf = dvol_vdf * this%dx_(i)
     enddo
   else
     dvol_vdf = 0.0_p_double
   endif

 end function dvol_vdf
!---------------------------------------------------


!---------------------------------------------------
function size_vdf( this )
!---------------------------------------------------
!---------------------------------------------------

  implicit none

  class( t_vdf ), intent(in) :: this

  integer, dimension(this%x_dim_) :: size_vdf

  select case (this%x_dim_)
    case(1)
      size_vdf(1) = size( this%f1, 2 )
    case(2)
      size_vdf(1) = size( this%f2, 2 )
      size_vdf(2) = size( this%f2, 3 )
    case(3)
      size_vdf(1) = size( this%f3, 2 )
      size_vdf(2) = size( this%f3, 3 )
      size_vdf(3) = size( this%f3, 4 )
  end select

end function size_vdf
!---------------------------------------------------

!---------------------------------------------------
function offset_vdf( this )
!---------------------------------------------------
!---------------------------------------------------

  implicit none

  class( t_vdf ), intent(in) :: this
  integer, dimension(this%x_dim_) :: offset_vdf

  select case (this%x_dim_)
    case(1)
      offset_vdf(1) = this%gc_num_( p_lower, 1 )
    case(2)
      offset_vdf(1) = this%gc_num_( p_lower, 1 )
      offset_vdf(2) = this%gc_num_( p_lower, 2 )
    case(3)
      offset_vdf(1) = this%gc_num_( p_lower, 1 )
      offset_vdf(2) = this%gc_num_( p_lower, 2 )
      offset_vdf(3) = this%gc_num_( p_lower, 3 )
  end select

end function offset_vdf
!---------------------------------------------------


!---------------------------------------------------
function size_dim_vdf( this, dim )
!---------------------------------------------------
!---------------------------------------------------

  implicit none

  integer :: size_dim_vdf
  class( t_vdf ), intent(in) :: this
  integer, intent(in) :: dim

  select case (this%x_dim_)
    case(1)
      size_dim_vdf = size( this%f1, dim )
    case(2)
      size_dim_vdf = size( this%f2, dim )
    case(3)
      size_dim_vdf = size( this%f3, dim )
    case default
      size_dim_vdf = 0
  end select

end function size_dim_vdf
!---------------------------------------------------

!---------------------------------------------------
function offset_dim_vdf( this, dim ) result(offset)
    
    implicit none

    class( t_vdf ), intent(in) :: this
    integer, intent(in) :: dim
    integer :: offset

    offset = this%gc_num_( p_lower, dim )

end function offset_dim_vdf
!---------------------------------------------------
    
    


!---------------------------------------------------
function gc_num_vdf( this )
!---------------------------------------------------
! gives number of guard cells at the boundaries
!---------------------------------------------------

  implicit none

! dummy variables

  class(t_vdf), intent(in) :: this
  integer, dimension(2,this%x_dim_) :: gc_num_vdf

! local variables - none

! executable statements

  gc_num_vdf = this%gc_num_(:,1:this%x_dim_)

end function gc_num_vdf
!---------------------------------------------------

!---------------------------------------------------
function gc_num_bnd_vdf( this, bnd_idx, i_dim )
!---------------------------------------------------
! gives number of guard cells at the boundaries for
! a given boundary/direction
!---------------------------------------------------

  implicit none

! dummy variables

  class(t_vdf), intent(in) :: this
  integer, intent(in) :: bnd_idx, i_dim
  integer :: gc_num_bnd_vdf

! local variables - none

! executable statements

  gc_num_bnd_vdf = this%gc_num_(bnd_idx,i_dim)

end function gc_num_bnd_vdf
!---------------------------------------------------

!---------------------------------------------------
function dx_vdf( this )
!---------------------------------------------------
! gives the sizes of the grid cells
!---------------------------------------------------

  implicit none

  ! dummy variables

  class(t_vdf), intent(in) :: this
  real(p_double) , dimension(this%x_dim_)  :: dx_vdf

  dx_vdf = this%dx_(1:this%x_dim_)

end function dx_vdf
!---------------------------------------------------

!---------------------------------------------------
function dx_dim_vdf( this, dim )
!---------------------------------------------------
! gives the sizes of the grid cells
!---------------------------------------------------

  implicit none

  ! dummy variables

  class(t_vdf), intent(in) :: this
  integer, intent(in) :: dim
  real(p_double) :: dx_dim_vdf

  dx_dim_vdf = this%dx_(dim)

end function dx_dim_vdf
!---------------------------------------------------

!---------------------------------------------------
function lbound_vdf( this )
!---------------------------------------------------
!---------------------------------------------------

  implicit none

  class( t_vdf ), intent(in) :: this
  integer, dimension( this%x_dim_ + 1) :: lbound_vdf

  select case (this%x_dim_)
    case(1)
      lbound_vdf = lbound( this%f1 )
    case(2)
      lbound_vdf = lbound( this%f2 )
    case(3)
      lbound_vdf = lbound( this%f3 )
  end select

end function lbound_vdf
!---------------------------------------------------


!---------------------------------------------------
function lbound_dim_vdf( this, dim )
!---------------------------------------------------
!---------------------------------------------------

  implicit none

  integer :: lbound_dim_vdf
  class( t_vdf ), intent(in) :: this
  integer, intent(in) :: dim

  select case (this%x_dim_)
    case(1)
      lbound_dim_vdf = lbound( this%f1, dim )
    case(2)
      lbound_dim_vdf = lbound( this%f2, dim )
    case(3)
      lbound_dim_vdf = lbound( this%f3, dim )
    case default
      lbound_dim_vdf = 0
  end select

end function lbound_dim_vdf
!---------------------------------------------------

!---------------------------------------------------
function ubound_vdf( this )
!---------------------------------------------------
!---------------------------------------------------

  implicit none

  class( t_vdf ), intent(in) :: this
  integer, dimension( this%x_dim_ + 1 ) :: ubound_vdf

  select case (this%x_dim_)
    case(1)
      ubound_vdf = ubound( this%f1 )
    case(2)
      ubound_vdf = ubound( this%f2 )
    case(3)
      ubound_vdf = ubound( this%f3 )
  end select

end function ubound_vdf
!---------------------------------------------------


!---------------------------------------------------
function ubound_dim_vdf( this, dim )
!---------------------------------------------------
!---------------------------------------------------

  implicit none

  integer :: ubound_dim_vdf
  class( t_vdf ), intent(in) :: this
  integer, intent(in) :: dim

  select case (this%x_dim_)
    case(1)
      ubound_dim_vdf = ubound( this%f1, dim )
    case(2)
      ubound_dim_vdf = ubound( this%f2, dim )
    case(3)
      ubound_dim_vdf = ubound( this%f3, dim )
    case default
      ubound_dim_vdf = 0
  end select

end function ubound_dim_vdf
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine check_nan_vdf( this, msg )
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

  implicit none

  class( t_vdf ),        intent(in) :: this
  character( len = * ) , intent(in), optional :: msg

  integer :: i1, i2, i3, f

  select case ( this%x_dim_ )

    case(1)

      do i1 = lbound( this%f1, 2 ), ubound( this%f1, 2 )
        do f = lbound( this%f1, 1 ), ubound( this%f1, 1 )
           if ( isnan( this%f1(f,i1) ) .or. isinf( this%f1(f,i1) ) ) then

              if ( present(msg) ) write(0,'(A,I0,A,A)') '[', mpi_node(), '] ', trim(msg)

              write(0,'(A,I0,A)') '[', mpi_node(), '] (* error *) check_nan_vdf failed '
              write(0,*) "[", mpi_node(), "] Bad cell [", f, i1, "]"
              write(0,*) "[", mpi_node(), "] Value = ", this%f1(:,i1)
              write(0,*) "[", mpi_node(), "] lbound = ", lbound( this%f1 )
              write(0,*) "[", mpi_node(), "] ubound = ", ubound( this%f1 )
              call abort_program( p_err_invalid )
           endif
        enddo
      enddo

    case(2)

      do i2 = lbound( this%f2, 3 ), ubound( this%f2, 3 )
        do i1 = lbound( this%f2, 2 ), ubound( this%f2, 2 )
          do f = lbound( this%f2, 1 ), ubound( this%f2, 1 )
             if ( isnan( this%f2(f,i1,i2) ) .or. isinf( this%f2(f,i1,i2) ) ) then

                  if ( present(msg) ) write(0,'(A,I0,A,A)') '[', mpi_node(), '] ', trim(msg)

                write(0,'(A,I0,A)') '[', mpi_node(), '] (* error *) check_nan_vdf failed '
                write(0,*) "[", mpi_node(), "] Bad cell [", f, i1, i2, "]"
                write(0,*) "[", mpi_node(), "] Value = ", this%f2(:,i1,i2)
                write(0,*) "[", mpi_node(), "] lbound = ", lbound( this%f2 )
                write(0,*) "[", mpi_node(), "] ubound = ", ubound( this%f2 )
                call abort_program( p_err_invalid )
             endif
          enddo
        enddo
      enddo

    case(3)

      do i3 = lbound( this%f3, 4 ), ubound( this%f3, 4 )
        do i2 = lbound( this%f3, 3 ), ubound( this%f3, 3 )
          do i1 = lbound( this%f3, 2 ), ubound( this%f3, 2 )
            do f = lbound( this%f3, 1 ), ubound( this%f3, 1 )
               if ( isnan( this%f3(f,i1,i2,i3) ) .or. isinf( this%f3(f,i1,i2,i3) ) ) then

                  if ( present(msg) ) write(0,'(A,I0,A,A)') '[', mpi_node(), '] ', trim(msg)

                  write(0,'(A,I0,A,A)') '[', mpi_node(), '] (* error *) check_nan_vdf failed '
                  write(0,*) "[", mpi_node(), "] Bad cell [", f, i1, i2, i3, "]"
                  write(0,*) "[", mpi_node(), "] Value = ", this%f3(:,i1,i2,i3)
                  write(0,*) "[", mpi_node(), "] lbound = ", lbound( this%f3 )
                  write(0,*) "[", mpi_node(), "] ubound = ", ubound( this%f3 )
                  call abort_program( p_err_invalid )
               endif
            enddo
          enddo
        enddo
      enddo

  end select

end subroutine check_nan_vdf
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine copy_scalar_r8_vdf( this, r8 )
!---------------------------------------------------
!       copies the value of r8 to vdf_a
!---------------------------------------------------

  implicit none

  class(t_vdf), intent(inout)  :: this
  real(p_double), intent(in) :: r8

  real(p_k_fld) :: val
  integer :: i, bsize

  val = real(r8, p_k_fld)
  bsize = size( this % buffer )

  !DIR$ vector nontemporal
  !$omp parallel do
  do i = 1, bsize
    this % buffer(i) = val
  enddo

end subroutine copy_scalar_r8_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine copy_scalar_r4_vdf( this, r4 )
!---------------------------------------------------
!       copies the value of r4 to vdf_a
!---------------------------------------------------

  implicit none

  class(t_vdf), intent(inout)  :: this
  real(p_single), intent(in) :: r4

  real(p_k_fld) :: val
  integer :: i, bsize

  val = real(r4, p_k_fld)
  bsize = size( this % buffer )

  !DIR$ vector nontemporal
  !$omp parallel do
  do i = 1, bsize
    this % buffer(i) = val
  enddo

end subroutine copy_scalar_r4_vdf
!---------------------------------------------------

end module m_vdf_define
