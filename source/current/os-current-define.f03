!#define __DEBUG__ 1

module m_current_define

#include "os-config.h"
#include "os-preprocess.fpp"
#include "memory/memory.h"

use m_system
use m_parameters

use m_utilities
use m_logprof


use m_vdf_define, only : t_vdf, t_vdf_report
use m_vdf_smooth, only : t_smooth, setup, smooth_order, if_smooth
use m_vdf_comm
use m_vdf_memory
use m_vdf

use m_node_conf
use m_grid_define
use m_space
use m_restart
use m_time_step

use m_current_diag

implicit none

private

! string to id restart data
character(len=*), parameter :: p_current_rst_id = "electric current rst data - 0x0001"
integer, public :: smoothev = 0, jayboundev

!-----------------------------------------------------------------------------------------
! t_current class
!-----------------------------------------------------------------------------------------

type :: t_current
  
  ! current iteration that is currently held by this
  !   t_current object. Used with CUDA to keep alert
  !   if jay needs to be copied from the GPU to be
  !   output for diagnostics.
  integer :: n_current_jay = 0

  ! vdf holding current data
  type( t_vdf ), dimension(:), pointer :: pf => null()

  !  definition of smoothing
  type( t_smooth )  ::  curr_smooth

  ! diagnostic setup
  class( t_current_diag ), pointer :: diag => null()

  ! coordinates type
  integer :: coordinates

  ! Position on global grid
  integer, dimension(p_max_dim) :: gix_pos

  ! boundary conditions
  integer, dimension(2,p_max_dim) :: bc_type

  ! interpolation level (used in the boundary conditions)
  integer :: interpolation

  contains

    procedure :: allocate_objs   => allocate_objs_current
    procedure :: init            => init_current
    procedure :: cleanup         => cleanup_current
    procedure :: smooth          => smooth_current
    procedure :: read_input      => read_input_current
    procedure :: update_boundary => update_boundary_current
    procedure :: report          => report_current
    procedure :: normalize_cyl   => normalize_cyl_current
    procedure :: fill_data       => fill_data_current


    procedure :: list_algorithm => list_algorithm_current
    procedure :: get_diag_buffer_size => get_diag_buffer_size_curr

    procedure :: read_checkpoint => read_checkpoint_current
    procedure :: write_checkpoint => write_checkpoint_current

    procedure :: get_min_gc => get_min_gc_current
    procedure :: nx         => nx_current
    procedure :: dx         => dx_current
    procedure :: reshape    => reshape_current

end type t_current


interface
  subroutine report_current( this, space, grid, no_co, tstep, t  )

  import t_current, t_space, t_grid, t_node_conf, t_time_step, p_double

  class( t_current ),  intent(inout) :: this
  type( t_space ),       intent(in) :: space
  class( t_grid ),        intent(in) :: grid
  class( t_node_conf ),   intent(in) :: no_co
  type( t_time_step ),   intent(in) :: tstep
  real(p_double),        intent(in) :: t

  end subroutine
end interface

public :: t_current

contains

!-----------------------------------------------------------------------------------------
! Allocate sub-objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_current( this )

  implicit none

  class( t_current ), intent(inout) :: this

  ! Allocate default t_diag_emf object class
  if ( .not. associated( this%diag ) ) then
    allocate( this%diag )
  endif

end subroutine allocate_objs_current

!-----------------------------------------------------------------------------------------
! Initialize object
!-----------------------------------------------------------------------------------------
subroutine init_current( this, grid, gc_min, dx, no_co, if_move, &
                          bcemf_type, interpolation, restart, restart_handle, sim_options )

  use m_restart
  use m_grid_define, only: t_grid
  use m_node_conf, only: t_node_conf, n_threads
  use m_system
  use m_parameters

  use m_current_diag
  use m_logprof

  implicit none

  class( t_current ), intent( inout )  ::  this

  class( t_grid ), intent(in) :: grid
  integer, dimension(:,:), intent(in) :: gc_min
  real(p_double), dimension(:), intent(in) :: dx
  class( t_node_conf ), intent(in)  :: no_co
  logical, dimension(:), intent(in) :: if_move
  integer, dimension(:,:), intent(in) :: bcemf_type
  integer, intent(in) :: interpolation

  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

  type( t_options ), intent(in) :: sim_options


  ! local variables

  integer, dimension(2,p_x_dim) :: gc_num_std
  integer, dimension(2,p_x_dim) :: gc_num_smooth
  integer, dimension(2,p_x_dim) :: gc_num_new

  integer :: i, k, nt

  ! executable statements

  ! setup smoothing first so that it can change the number of ghost cells if neccessary
  ! note that the upper boundary needs to have at least one more guard cell than the
  ! smoothing level since the current in the first upper guard cell is required in the
  ! advance of the e-field (dedt in m_el_mag_fld)

  ! also currently smooth has no restart information so set it up before restart
  call setup( this%curr_smooth, dx, sim_options%gamma )

  if ( restart ) then

    call this % read_checkpoint( no_co, restart_handle )

  else

    this%n_current_jay = 0
    
    gc_num_std = gc_min

    ! add guard cells for moving window algorithm
    do i = 1, p_x_dim
      if ( if_move(i)) gc_num_std(p_lower,i) = gc_num_std(p_lower,i) + 1
    enddo

    gc_num_smooth(p_lower,:) = smooth_order(this%curr_smooth)
    gc_num_smooth(p_upper,:) = gc_num_smooth(p_lower,:)+1

    do i=1, p_x_dim
      gc_num_new(p_lower,i)=max(gc_num_std(p_lower,i),gc_num_smooth(p_lower,i))
      gc_num_new(p_upper,i)=max(gc_num_std(p_upper,i),gc_num_smooth(p_upper,i))
    enddo

    ! setup of of the vdf object for the field

    nt = n_threads( no_co )

    call alloc( this%pf, (/ nt  /) )

    ! V0 - Main process allocates all temporary current buffers
    ! do i = 1, nt
    !    call this%pf(i) % new( p_x_dim, p_f_dim, grid%my_nx(3,:), gc_num_new, dx, .false. )
    ! enddo

    ! V1 - Memory allocation is done by each OpenMP thread. This should cause the memory to
    !      be allocated from a heap local to each thread, which combined with cpu thread affinity,
    !      improves the performance on NUMA nodes (e.g. Cray XT5). We use a critical section to
    !      ensure that the actual memory allocation is not done in parallel (even if the system
    !      allows it, our internal memory accounting routines do not)
    !$omp parallel do
    do i = 1, nt
      !$omp critical
      call this%pf(i) % new( p_x_dim, p_f_dim, grid%my_nx(3,:), gc_num_new, dx, .false. )
      !$omp end critical
    enddo
    !$omp end parallel do


  endif

  ! setup data that is not on restart file
  this%coordinates = grid%coordinates
  do i = 1, p_x_dim
    this%gix_pos(i) = grid%my_nx( p_lower, i )
  enddo

  ! The boundary conditions for currents are set to match boundary conditions for EM
  ! fields.
  do i = 1, p_x_dim
    do k = 1, 2
      if ( ( bcemf_type( k, i ) == p_bc_pec ) .or. &
        ( bcemf_type( k, i ) == p_bc_pmc ) ) then
        this%bc_type( k, i ) = p_bc_specular
      else
        this%bc_type( k, i ) = p_bc_none
      endif
    enddo
  enddo

  ! Interpolation level
  this%interpolation = interpolation

  ! setup diagnostic data (no restart data)
  call this%diag%init( interpolation )

  ! Create timing events
  if(smoothev==0) then
    smoothev      = create_event('current smooth')
    jayboundev    = create_event('update current boundary')
  endif

end subroutine init_current
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleanup object
!-----------------------------------------------------------------------------------------
subroutine cleanup_current(this)

  use m_current_diag
  use m_vdf_smooth

  implicit none

  class( t_current ),   intent( inout ) ::  this

  integer :: i

  if (associated(this%pf)) then
    do i = 1, size( this%pf )
      call this % pf(i) % cleanup()
    enddo
    call freemem( this%pf )
  endif

  call cleanup( this%curr_smooth )

  call this % diag % cleanup( )

end subroutine cleanup_current

!-----------------------------------------------------------------------------------------
! Smooth (digital filter)
!-----------------------------------------------------------------------------------------
subroutine smooth_current( this )

  use m_logprof

  use m_vdf_define
  use m_vdf_smooth
  use m_vdf_comm
  use m_vdf_memory
  use m_vdf

  implicit none

  class( t_current ), intent(inout) :: this

  call begin_event(smoothev)
  call smooth( this%pf(1), this%curr_smooth )
  call end_event(smoothev)

end subroutine smooth_current

!-----------------------------------------------------------------------------------------
! Read input
!  - this is just a separator for the input file, no data is actually read
!-----------------------------------------------------------------------------------------
subroutine read_input_current( this, input_file )

  use m_input_file
  use m_vdf_smooth
  use m_current_diag

  implicit none

  class( t_current ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  integer :: ierr

  call get_namelist( input_file, "nl_current", ierr )

  call this%allocate_objs()

  call read_nml( this%curr_smooth, input_file )

  call this%diag%read_input( input_file )

end subroutine read_input_current

!-----------------------------------------------------------------------------------------
! Normalize deposited current for ring charges and take care of axial boundary.
! Note that current beyond axial boundary is reversed since r is considered to be < 0.
! This algorithm assumes that B1 is defined on the cylindrical axis.
!-----------------------------------------------------------------------------------------
subroutine normalize_cyl_current( this )

  use m_parameters
  use m_vdf_define, only: t_vdf

  implicit none

  class( t_current ), intent( inout ), target  ::  this

  real(p_k_fld) :: vol_fac_cor, vol_fac_mid, r_dr
  integer :: i1, i2, gshift_i2

  type( t_vdf ), pointer :: j

  !  ASSERT( this%coordinates == p_cylindrical_b )

  ! Normalize current for ring charges
  j => this%pf(1)

  r_dr      = real( 1.0_p_double / j%dx_( p_r_dim), p_k_fld )
  gshift_i2 = this%gix_pos(2) - 2

  do i2 = lbound(j%f2, 3), ubound(j%f2, 3)

    ! Get normalization factors 1/r at cell corner and middle
    vol_fac_cor = r_dr / (i2 + gshift_i2 - 0.5_p_k_fld)
    if ( i2 + gshift_i2 /= 0 ) then
      vol_fac_mid = r_dr / (i2 + gshift_i2)
    else
      vol_fac_mid = 0.0
    endif

    do i1 = lbound(j%f2, 2), ubound(j%f2, 2)
      j%f2(1,i1,i2) = j%f2(1,i1,i2) * vol_fac_cor
      j%f2(2,i1,i2) = j%f2(2,i1,i2) * vol_fac_mid
      j%f2(3,i1,i2) = j%f2(3,i1,i2) * vol_fac_cor
    enddo
  enddo

  ! process axial boundary
  ! This could be moved into update_boundary_current to overlap this calculation with
  ! communication, but this is cleaner and makes this routine interchangeable with
  ! (update_boundary and smooth).

  if ( gshift_i2 < 0 ) then

    ! note that pf%f2(1,:,1) and pf%f2(3,:,1) are off axis, but that pf%f2(2,:,1) is on axis
    ! first take care of the "physical" cell on axis - this is not an actual guard cell

    do i1 = lbound( j%f2, 2 ), ubound( j%f2, 2 )
      j%f2(1,i1,2) =  j%f2(1,i1,2) + j%f2(1,i1,1)
      j%f2(1,i1,1) =  j%f2(1,i1,2)

      j%f2(2,i1,1) =  0.0_p_k_fld

      j%f2(3,i1,2) =  j%f2(3,i1,2) + j%f2(3,i1,1)
      j%f2(3,i1,1) =  j%f2(3,i1,2)
    enddo

    do i2=1, j%gc_num_( p_lower, p_r_dim ) ! use bound. guards
      do i1 = lbound( j%f2, 2 ), ubound( j%f2, 2 )
        ! fold guard cells back into simul. space
        j%f2( 1, i1, i2+2 ) = j%f2( 1, i1, i2+2 ) + j%f2( 1, i1, 1-i2 )
        j%f2( 2, i1, i2+1 ) = j%f2( 2, i1, i2+1 ) - j%f2( 2, i1, 1-i2 )
        j%f2( 3, i1, i2+2 ) = j%f2( 3, i1, i2+2 ) - j%f2( 3, i1, 1-i2 )

        ! update values in guard cells
        j%f2( 1, i1, 1-i2 ) =   j%f2( 1, i1, i2+2 )
        j%f2( 2, i1, 1-i2 ) = - j%f2( 2, i1, i2+1 )
        j%f2( 3, i1, 1-i2 ) = - j%f2( 3, i1, i2+2 )
      enddo
    enddo

  endif

end subroutine normalize_cyl_current
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Printout the algorithm used by the electric current object
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_current( this )

  implicit none

  class( t_current ), intent(in) :: this

  write(*,*) ' '
  write(*,*) 'Electrical current :'

  ! smoothing
  if (if_smooth(this%curr_smooth)) then
    write(*,*) '- Using smoothing'
  else
    write(*,*) '- No smoothing done'
  endif


end subroutine list_algorithm_current

!-----------------------------------------------------------------------------------------
! Get required diagnostic buffer size
!-----------------------------------------------------------------------------------------
subroutine get_diag_buffer_size_curr( this, gnx, diag_buffer_size )

  implicit none

  class( t_current ), intent(in) :: this
  integer, dimension(:), intent(in) :: gnx
  integer, intent(inout) :: diag_buffer_size

  call this%diag % get_diag_buffer_size( gnx, diag_buffer_size )

end subroutine get_diag_buffer_size_curr
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Write checkpoint information
!-----------------------------------------------------------------------------------------
subroutine write_checkpoint_current( this, restart_handle )

  use m_restart

  implicit none


  class( t_current ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  integer :: ierr
  character(len=*), parameter :: err_msg = 'error writing restart data for electric current object.'

  restart_io_wr( p_current_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%n_current_jay, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  ! no need to write electric current data, only structure information
  restart_io_wr( this%pf(1)%x_dim_, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%pf(1)%f_dim_, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%pf(1)%nx_, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%pf(1)%gc_num_, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%pf(1)%dx_, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

end subroutine write_checkpoint_current

!-----------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------
subroutine read_checkpoint_current( this, no_co, restart_handle )

  use m_restart

  implicit none

  ! dummy variables
  class( t_current ), intent(inout) :: this
  class( t_node_conf ), intent(in) :: no_co
  type( t_restart_handle ), intent(in) :: restart_handle

  ! local variables
  character(len=len(p_current_rst_id)) :: rst_id
  integer :: lx_dim, f_dim, i, nt
  integer, dimension(p_max_dim) :: lnx
  integer, dimension(2, p_max_dim) :: lgc_num
  real(p_double), dimension(p_max_dim) :: ldx

  integer :: ierr
  character(len=*), parameter :: err_msg = 'error reading restart data for electric current object.'

  ! executable statements
  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  ! check if restart file is compatible
  if ( rst_id /= p_current_rst_id) then
    ERROR('Corrupted restart file, or restart file ')
    ERROR('from incompatible binary (electric current)')
    call abort_program(p_err_rstrd)
  endif

  restart_io_rd( this%n_current_jay, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  ! read electric current parameters and create electric current structure
  restart_io_rd( lx_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( f_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( lnx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( lgc_num, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( ldx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  ! setup of of the vdf object for the field
  nt = n_threads( no_co )
  call alloc( this%pf, (/ nt  /) )

  !$omp parallel do
  do i = 1, nt
    !$omp critical
    call this%pf(i) % new( lx_dim, f_dim, lnx, lgc_num, ldx, .false. )
    !$omp end critical
  enddo
  !$omp end parallel do

end subroutine read_checkpoint_current
!---------------------------------------------------

!---------------------------------------------------
function get_min_gc_current( this )
!---------------------------------------------------
! returns the number of guard cells required by this
! object. Does not require setup
!---------------------------------------------------

  implicit none

  class( t_current ), intent(in) :: this
  integer, dimension(2,p_x_dim) :: get_min_gc_current

  integer, dimension(2,p_x_dim) :: gc_num_temp
  integer :: i

  get_min_gc_current(1,:) = 1
  get_min_gc_current(2,:) = 2

  gc_num_temp(1,:) = smooth_order(this%curr_smooth)
  gc_num_temp(2,:) = gc_num_temp(1,:)+1

  do i = 1, p_x_dim
  get_min_gc_current(1,i)=max(get_min_gc_current(1,i),gc_num_temp(1,i))
  get_min_gc_current(2,i)=max(get_min_gc_current(2,i),gc_num_temp(2,i))
  enddo

end function get_min_gc_current
!---------------------------------------------------


!---------------------------------------------------
function nx_current( this )
!---------------------------------------------------
! returns the number of guard cells required by this
! object. Does not require setup
!---------------------------------------------------

  implicit none

  class( t_current ), intent(in) :: this
  integer, dimension(p_x_dim) :: nx_current

  nx_current = this%pf(1)%nx_(1:p_x_dim)

end function nx_current
!---------------------------------------------------

!---------------------------------------------------
function dx_current( this )
!---------------------------------------------------
! returns the number of guard cells required by this
! object. Does not require setup
!---------------------------------------------------

  implicit none

  class( t_current ), intent(in) :: this
  real( p_double ), dimension(p_x_dim) :: dx_current

  dx_current = this%pf(1)%dx_(1:p_x_dim)

end function dx_current
!---------------------------------------------------


!---------------------------------------------------
subroutine reshape_current( this, old_grid, new_grid, no_co, send_msg, recv_msg )
!---------------------------------------------------
! reshape phy_field object when node grids change
!---------------------------------------------------

    implicit none

    class( t_current ), intent(inout) :: this
    class( t_grid ), intent(in) :: old_grid, new_grid
    class( t_node_conf ), intent(in) :: no_co
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

    integer :: i
    type(t_vdf_report), pointer :: report

    ! reshape the vdf objects
    do i = 1, size( this%pf )
        call reshape_nocopy( this%pf(i), new_grid )
        call this % pf(i) % zero()
    enddo

    ! store new position on global grid
    do i = 1, this % pf(1) % x_dim_
        this%gix_pos(i) = new_grid%my_nx( 1, i )
    enddo

    ! Reshape tavg_data
    report => this%diag%reports
    do
        if ( .not. associated(report) ) exit
        if ( report%tavg_data%x_dim_ > 0 ) then
            call reshape_copy( report%tavg_data, old_grid, new_grid, no_co, send_msg, recv_msg )
        endif
        report => report%next
    enddo

end subroutine reshape_current
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Update boundary of electrical current
!---------------------------------------------------------------------------------------------------
subroutine update_boundary_current( this, no_co, send_msg, recv_msg )

  use m_vdf_comm
  use m_node_conf
  use m_logprof
  use m_parameters

  use m_current_boundary

  implicit none

  class( t_current ), intent(inout) :: this
  class( t_node_conf ),   intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i

  call begin_event( jayboundev )

  ! it is not necessary to use move_num(space) since the values for the current are recalculated
  ! at each time-step so no data needs to be shifted

  ! Take care of communication with other nodes and single-node periodic
  call update_boundary( this%pf(1), p_vdf_add, no_co, send_msg, recv_msg )

  ! Process specular boundaries
  select case ( p_x_dim )
    case (1)
      if ( this%bc_type( p_lower, 1 ) == p_bc_specular ) call spec_bc_1d( this%pf(1), p_lower, this%interpolation )
      if ( this%bc_type( p_upper, 1 ) == p_bc_specular ) call spec_bc_1d( this%pf(1), p_upper, this%interpolation )
    case (2)
      do i = 1, 2
        if ( this%bc_type( p_lower, i ) == p_bc_specular ) call spec_bc_2d( this%pf(1), i, p_lower, this%interpolation )
        if ( this%bc_type( p_upper, i ) == p_bc_specular ) call spec_bc_2d( this%pf(1), i, p_upper, this%interpolation )
      enddo
    case (3)
      do i = 1, 3
        if ( this%bc_type( p_lower, i ) == p_bc_specular ) call spec_bc_3d( this%pf(1), i, p_lower, this%interpolation )
        if ( this%bc_type( p_upper, i ) == p_bc_specular ) call spec_bc_3d( this%pf(1), i, p_upper, this%interpolation )
      enddo
  end select

  call end_event( jayboundev )

end subroutine update_boundary_current
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Used for CUDA (and similar) situations where data may have to be copied from another location.
!---------------------------------------------------------------------------------------------------
subroutine fill_data_current(this, tstep)
  implicit none
  class( t_current ), intent(inout) :: this
  type( t_time_step ),   intent(in) :: tstep
  
  ! If used, the general format for this function would be along the lines of:
  !
  ! if( n(tstep) /= this%n_current_jay) then
  !   ! get the data needed
  !   call ....
  !   ! update the status counter so that the data is only copied as needed
  !   this%n_current_jay = n(tstep)
  ! endif

  ! put a fatal error here in case any subclasses forgot to increment 'n_current_jay'
  !   this error is most likely from a subclass not having the line
  !      "this%n_current_jay=this%n_current_jay+1" 
  !   after advanced deposit. Ask Ricardo or Anton or Adam if needed.
  if( n(tstep) /= this%n_current_jay) then
    ERROR( "Error in fill_data_current.. see 'fill_data_current' is os-current-define.f03")
    call abort_program(p_err_rstrd)
  endif

end subroutine


end module m_current_define


!---------------------------------------------------
subroutine report_current( this, space, grid, no_co, tstep, t )

  use m_time_step
  use m_space
  use m_node_conf
  use m_current_diag
  use m_grid_define
  use m_parameters
  use m_vdf_define
  use m_vdf_report
  use m_logprof
  use m_current_define
  use m_vdf_math

  implicit none

  class( t_current ), intent(inout) :: this
  type( t_space ),       intent(in) :: space
  class( t_grid ),       intent(in) :: grid
  class( t_node_conf ),  intent(in) :: no_co
  type( t_time_step ),   intent(in) :: tstep
  real(p_double),        intent(in) :: t

  ! report attributes
  type( t_vdf_report ), pointer :: rep

  ! temporary vdf object for diagnostics
  type( t_vdf ),pointer :: pf
  type( t_vdf ) :: vdf_a

  ! executable statements

  pf => this%pf(1)
  call begin_event(diag_current_ev)

  rep => this%diag%reports
  do
    if ( .not. associated( rep ) ) exit

    if ( if_report( rep, tstep ) ) then

      call this%fill_data( tstep )
      
      select case ( rep%quant )
      case ( p_j1, p_j2, p_j3 )
        call report_vdf( rep, pf, rep%quant - p_j1 + 1, space, grid, no_co, tstep, t )

      case ( p_div_j )
        call vdf_a % new( pf, f_dim = 1 )
        call div( pf, vdf_a )
        call report_vdf( rep, vdf_a, 1, space, grid, no_co, tstep, t )
        call vdf_a % cleanup()

      case default
        ! must be a superclass diagnostic, ignore
        continue
      end select

    endif

    rep => rep%next
  enddo

  call end_event(diag_current_ev)


end subroutine report_current
!---------------------------------------------------

