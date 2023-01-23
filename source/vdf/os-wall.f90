!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     wall class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_wall

#include "memory/memory.h"

  use m_wall_define

  use m_vdf_define

  use m_parameters
  use m_utilities

  use m_vdf_report

  use m_grid_define
  use m_node_conf

  implicit none

!       restrict access to things explicitly declared public
  private

  ! string to id restart data
  character(len=*), parameter :: p_wall_rst_id = "wall rst data - 0x0003"

  interface cleanup
    module procedure cleanup_wall
  end interface

  interface restart_write
    module procedure restart_write_wall
  end interface

  interface new
    module procedure new_wall
  end interface

  interface assignment(=)
    module procedure copy_scalar_double_wall
    module procedure copy_scalar_single_wall
  end interface

  ! this is for debug purposes only
  interface report_wall
    module procedure report_diag_wall
  end interface

!       declare things that should be public
  public :: cleanup, new
  public :: restart_write
  public :: assignment(=)

  public :: report_wall

 contains

!---------------------------------------------------
subroutine cleanup_wall(this)
!---------------------------------------------------
! clear this wall
!---------------------------------------------------

  implicit none

  ! dummy variables

  type( t_wall ),   intent( inout ) ::  this

  ! executable statements

  this%idir  = -1
  this%ibnd  = -1
  this%range = 0
  if(associated(this%f1)) then
    call freemem( this%f1 )
  endif
  if(associated(this%f2)) then
    call freemem( this%f2 )
  endif
  if(associated(this%f3)) then
    call freemem( this%f3 )
  endif

  this%x_dim = -1
  this%f_dim = -1
  this%nx = 0
  this%gc_num = 0
  this%gc_num_comm = 0

end subroutine cleanup_wall
!---------------------------------------------------

!------------------------------------------------------------------------------
!> Creates a new wall object
!! @param   this            Wall object
!! @param   vdf_source      VDF parent object to get dimensions from
!! @param   idir            Wall normal direction
!! @param   ibnd            Wall position (p_lower, p_upper)
!! @param   range           Range of wall positions
!! @param   restart         Set flag to load object from checkpoint data
!! @param   restart_handle  Checkpoint file handle
!! @param   fdim_           (optional) Number of field components, defaults to
!!                          the same number of components as the VDF parent
!! @param   corner_         (optional) Number of corner cells, defaults to 0
subroutine new_wall( this, vdf_source, idir, ibnd, range, &
    restart, restart_handle , &
    f_dim, corner )

  use m_restart

  implicit none

  type( t_wall ),   intent( inout ) ::  this
  type( t_vdf ),    intent( in ) :: vdf_source
  integer,               intent(in) :: idir, ibnd
    integer, dimension(:), intent(in) :: range

  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

    integer, intent(in), optional :: f_dim
    integer, dimension(:,:), intent(in), optional :: corner

    integer :: i
    integer, dimension(2,p_max_dim) :: corner_
  integer, dimension(p_max_dim + 1) :: lb, ub

  call cleanup( this )

  if ( restart ) then
    call restart_read_wall( this, restart_handle )
  else
        ! Get spatial dimensions
  this%x_dim = vdf_source % x_dim_

        ! Get field dimensions
        if ( present(f_dim) ) then
            this % f_dim = f_dim
        else
            this%f_dim = vdf_source % f_dim()
        endif

        ! Get corner parameters
        if ( present(corner) ) then
            do i = 1, vdf_source % x_dim_
                corner_(p_lower,i) = corner(p_lower,i)
                corner_(p_upper,i) = corner(p_upper,i)
            enddo
        else
            corner_ = 0
        endif


  if (( idir < 1 ) .or. (idir > this%x_dim)) then
    ERROR("Invalid direction for wall, ", idir)
    call abort_program( p_err_invalid )
  endif

  if (( ibnd /= p_upper ) .and. ( ibnd /= p_lower )) then
    ERROR("Invalid position for wall, ibnd = ", ibnd)
    call abort_program( p_err_invalid )
  endif

  ! set wall specific data
  this%idir  = idir
  this%ibnd  = ibnd
  this%range = range

  lb(1) = 1
  ub(1) = this%f_dim

        do i = 1, this % x_dim
            if ( i == idir ) then
                this%gc_num_comm(p_lower,i) = 0
                this%gc_num_comm(p_upper,i) = 0
                this%gc_num(p_lower,i) = 0
                this%gc_num(p_upper,i) = 0
                this%nx(i) = range(p_upper) - range(p_lower) + 1
                lb(i+1) = range(1)
                ub(i+1) = range(2)
            else
                this%gc_num_comm(p_lower,i) = vdf_source%gc_num_(p_lower, i)
                this%gc_num_comm(p_upper,i) = vdf_source%gc_num_(p_upper, i)
                this%gc_num(p_lower,i) = max(vdf_source%gc_num_(p_lower, i), corner_(p_lower,i)-1)
                this%gc_num(p_upper,i) = max(vdf_source%gc_num_(p_upper, i), corner_(p_upper,i))
                this%nx(i) = vdf_source%nx_(i)
                lb(i+1) = 1 - this%gc_num(p_lower, i)
                ub(i+1) = this%nx(i) + this%gc_num(p_upper, i)
            endif
        enddo

  ! allocate data structure
  select case (this%x_dim)
    case(1)
        call alloc( this%f1, lb, ub )
            this%f1 = 0

    case(2)
        call alloc( this%f2, lb, ub )
            this%f2 = 0

    case(3)
    call alloc( this%f3, lb, ub )
            this%f3 = 0

    end select
  endif

end subroutine new_wall

!---------------------------------------------------
subroutine restart_write_wall( this, restart_handle )
!---------------------------------------------------
!       write object information into a restart file
!---------------------------------------------------

  use m_restart

  implicit none

  type( t_wall ),   intent(in) ::  this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for wall object'
  integer :: ierr

  ! Restart id tag
  restart_io_wr( p_wall_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! wall grid information
  restart_io_wr( this%idir, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%ibnd, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%range, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%x_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%f_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%gc_num, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%gc_num_comm, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  if (this%x_dim > 0) then

  select case (this%x_dim )

    case(1)
      restart_io_wr( lbound(this%f1), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      restart_io_wr( ubound(this%f1), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%f1, restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    case(2)
      restart_io_wr( lbound(this%f2), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      restart_io_wr( ubound(this%f2), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%f2, restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    case(3)
      restart_io_wr( lbound(this%f3), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

      restart_io_wr( ubound(this%f3), restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

    restart_io_wr( this%f3, restart_handle, ierr )
        CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  end select

  endif

end subroutine restart_write_wall
!---------------------------------------------------

!---------------------------------------------------
subroutine restart_read_wall( this, restart_handle )
!---------------------------------------------------
!       read object information from a restart file
!---------------------------------------------------

  use m_restart

  implicit none

  type( t_wall ), intent(inout) ::  this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for wall object'

  character(len=len(p_wall_rst_id)) :: rst_id
  integer, dimension(4) :: lb, ub
  integer :: ierr

  ! check if restart file is compatible
  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  if ( rst_id /= p_wall_rst_id) then
  ERROR('Corrupted restart file, or restart file ')
  ERROR('from incompatible binary (wall)')
  call abort_program(p_err_rstrd)
  endif

  ! read number of dimensions and grid parameters
  restart_io_rd( this%idir, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%ibnd, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%range, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%x_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%f_dim, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%nx, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%gc_num, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%gc_num_comm, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  ! read wall data if necessary
   if (this%x_dim > 0) then

     restart_io_rd( lb(1:this%x_dim+1), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstrd )

     restart_io_rd( ub(1:this%x_dim+1), restart_handle, ierr )
     CHECK_ERROR( ierr, err_msg, p_err_rstrd )

   select case ( this%x_dim )

     case(1)
     call alloc( this%f1, lb, ub )
     restart_io_rd( this%f1, restart_handle, ierr )

     case(2)
     call alloc( this%f2, lb, ub )
     restart_io_rd( this%f2, restart_handle, ierr )

     case(3)
     call alloc( this%f3, lb, ub )
     restart_io_rd( this%f3, restart_handle, ierr )

   end select
     CHECK_ERROR( ierr, err_msg, p_err_rstrd )

   endif

end subroutine restart_read_wall
!---------------------------------------------------

!---------------------------------------------------
subroutine copy_scalar_double_wall( wall_a, double_b )
!---------------------------------------------------
!       copies the value of double_b to wall_a
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_wall), intent(inout)  :: wall_a
  real(p_double), intent(in) :: double_b

!       local variables - none

!       executable statements

  select case (wall_a%x_dim)

  case(1)
    wall_a%f1 = real( double_b, p_k_fld )

  case(2)
    wall_a%f2 = real( double_b, p_k_fld )

  case(3)
    wall_a%f3 = real( double_b, p_k_fld )
  end select

end subroutine copy_scalar_double_wall
!---------------------------------------------------

!---------------------------------------------------
subroutine copy_scalar_single_wall( wall_a, single_b )
!---------------------------------------------------
!       copies the value of single_b to wall_a
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_wall), intent(inout)  :: wall_a
  real(p_single), intent(in) :: single_b

!       local variables - none

!       executable statements

  select case (wall_a%x_dim)

  case(1)
    wall_a%f1 = real( single_b, p_k_fld )

  case(2)
    wall_a%f2 = real( single_b, p_k_fld )

  case(3)
    wall_a%f3 = real( single_b, p_k_fld )
  end select

end subroutine copy_scalar_single_wall
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Write wall values to diagnostics file. Used for debug only
!---------------------------------------------------------------------------------------------------
subroutine report_diag_wall( wall_e, wall_b, i_wall, space, tstep, t, dx )

  use m_vdf
  use m_time_step
  use m_diagnostic_utilities
  use m_space, only : t_space, get_x_bnd, xmin, set_x_bnd

  implicit none

  ! dummy variables

  type( t_wall ),                 intent(in) :: wall_e, wall_b

  integer, intent(in) :: i_wall
  type( t_space ),               intent(in) :: space
  type( t_time_step ),          intent(in) :: tstep
  real(p_double),               intent(in) :: t
  real(p_double), dimension(:), intent(in) :: dx

!       local variables

  ! report attributes
  type( t_vdf_report ) :: wall_report
  type( t_vdf ) :: vdf_e, vdf_b

  integer :: i, wall_dir
  integer, parameter :: izero = ichar('0')

  character(*), dimension(2), parameter :: pos = (/'lower','upper'/)

  integer :: ndump_fac_all
  logical :: if_dump_b, if_dump_e

  real( p_double ), dimension( 2, p_x_dim ) :: x_bnd

  type( t_space ) :: l_space

  ! executable statements

  ! SFM: Change manually for diagnostics
  if_dump_b = .true.
  if_dump_e = .true.
  ndump_fac_all = 1

   ! self generated fields diagnostics
  if ( test_if_report( tstep, ndump_fac_all ) ) then

  ! create vdf from wall object

  call vdf_e % new( wall_e%x_dim, wall_e%f_dim, wall_e%nx, wall_e%gc_num, dx, .false.)
  call vdf_b % new( wall_b%x_dim, wall_b%f_dim, wall_b%nx, wall_b%gc_num, dx, .false.)

  l_space = space
  wall_dir = wall_e%idir

    ! Correct boundaries
  call get_x_bnd( l_space, x_bnd )

  select case (wall_e%x_dim)

    case(1)
    vdf_e%f1 = wall_e%f1
    vdf_b%f1 = wall_b%f1

    x_bnd( p_lower, 1 ) = xmin(l_space,1) + (lbound( wall_e%f1, 2 )-1) * dx(1)
    x_bnd( p_upper, 1 ) = xmin(l_space,1) + ubound( wall_e%f1, 2 ) * dx(1)

    case(2)
    vdf_e%f2 = wall_e%f2

    vdf_b%f2 = wall_b%f2

    do i = 1, 2
       x_bnd( p_lower, i ) = xmin(l_space,i) + (lbound( wall_e%f2, i+1 )-1) * dx(i)
       x_bnd( p_upper, i ) = xmin(l_space,i) + ubound( wall_e%f2, i+1 ) * dx(i)
        enddo

    case(3)
    vdf_e%f3 = wall_e%f3
    vdf_b%f3 = wall_b%f3

    do i = 1, 3
       x_bnd( p_lower, i ) = xmin(l_space,i) + (lbound( wall_e%f3, i+1 )-1) * dx(i)
       x_bnd( p_upper, i ) = xmin(l_space,i) + ubound( wall_e%f3, i+1 ) * dx(i)
        enddo

  end select

  call set_x_bnd( l_space, x_bnd )

  wall_report%xname  = (/'x1', 'x2', 'x3'/)

  wall_report%xlabel = (/'x_1', 'x_2', 'x_3'/)
  wall_report%xunits = (/'c / \omega_p', &
              'c / \omega_p', &
              'c / \omega_p'/)

  wall_report%prec = p_diag_prec

  wall_report%t = t
  wall_report%n = n( tstep )

  ! units are the same for e and b fields
  wall_report%units = 'm_e c \omega_p e^{-1}'

    ! all wall dumps

  do i=1, wall_b%f_dim

    if ( if_dump_b ) then

      wall_report%label = 'B_'//char(izero+i)//', x_'//char(izero+wall_dir)//' '//&
                          pos( wall_e%ibnd )//' wall ['//char(izero+i_wall)//']'
      wall_report%name  = 'B'//char(izero+i)

      wall_report%path = trim(path_mass) // 'FLD_WALL_x'//char(izero+wall_dir) // '_' // pos( wall_e%ibnd ) &
                // p_dir_sep // trim(wall_report%name)//p_dir_sep

      wall_report%filename = get_filename( n( tstep ) / ndump( tstep ), wall_report%name )

      call report_array( vdf_b, i, wall_report )

    endif

  enddo

  ! loop through all e-field components and report the required ones

  do i=1, wall_e%f_dim

    if ( if_dump_e ) then

      wall_report%label = 'E_'//char(izero+i)//', x_'//char(izero+wall_dir)//' '//&
                          pos( wall_e%ibnd )//' wall ['//char(izero+i_wall)//']'
      wall_report%name  = 'E'//char(izero+i)

  !    wall_report%path = trim(path_mass) // 'FLD_WALL' //char(izero+i_wall) &
  !               // p_dir_sep // trim(wall_report%name)//p_dir_sep

      wall_report%path = trim(path_mass) // 'FLD_WALL_x'//char(izero+wall_dir) // '_' // pos( wall_e%ibnd ) &
                // p_dir_sep // trim(wall_report%name)//p_dir_sep

      wall_report%filename = get_filename( n( tstep ) / ndump( tstep ), wall_report%name )

      call report_array( vdf_e, i, wall_report )

    endif

  enddo

    call vdf_e%cleanup()
    call vdf_b%cleanup()

  endif

contains

  subroutine report_array( this, fc, report )

  use m_vdf_reportfile

  implicit none

  type( t_vdf ),        intent(inout) :: this

  integer,              intent(in) :: fc

  type( t_vdf_report ), intent(in) :: report

  class(t_diag_file), allocatable :: diagFile

  call create_diag_file( diagFile )
  call diagFile % open( p_diag_create )

  diagFile % grid % name  = report%name
  diagFile % grid % label = report%label
  diagFile % grid % units = report%units

  select case( this % x_dim_ )
    case(1)
     call diagFile % add_dataset( report%name, this%f1(fc,:) )
     !call add_h5_dataset( diagFile%id, report%name, this%f1(fc,:), &
     !           units = report%units, long_name = report%label )
    case(2)
     call diagFile % add_dataset( report%name, this%f2(fc,:,:) )
     !call add_h5_dataset( diagFile%id, report%name, this%f2(fc,:,:), &
     !           units = report%units, long_name = report%label )
    case(3)

     if ( mod(fc,2) == 1 ) then
      this%f3(fc,:,:,:) = this%f3(fc,:,:,:) + this%f3(fc+1,:,:,:)
     endif

     call diagFile % add_dataset( report%name, this%f3(fc,:,:,:) )
     !call add_h5_dataset( diagFile%id, report%name, this%f3(fc,:,:,:), &
     !           units = report%units, long_name = report%label )
  end select

  ! close file
  call diagFile % close( )

  end subroutine report_array

end subroutine report_diag_wall
!---------------------------------------------------------------------------------------------------

end module m_wall
