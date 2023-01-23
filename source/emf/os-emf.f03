!#define DEBUG_FILE 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     electro_magnetic_fields class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For a description of the 'stencil' field solver see
! A.D. Greenwood et al. / Journal of Computational Physics 201 (2004) 665–684

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf

#include "memory/memory.h"

use m_system
use m_parameters

use m_emf_define
use m_emf_bound
use m_emf_diag
use m_emf_solver
use m_emf_interpolate
use m_emf_gridval
use m_emf_psi

use m_logprof
use m_space
use m_grid_define
use m_node_conf

use m_diagnostic_utilities

use m_fparser
use m_parameters

use m_vdf_define
use m_vdf_smooth
use m_vdf_comm
use m_vdf_math
use m_vdf_report
use m_vdf_memory
use m_vdf

use m_restart

implicit none

private

interface reshape_obj
  module procedure reshape_emf
end interface

!public :: t_emf, t_emf_bound
public :: is_open
public :: list_algorithm
public :: cleanup, reshape_obj


contains


!-----------------------------------------------------------------------------------------
subroutine reshape_emf( this, old_lb, new_lb, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

  implicit none

  class(t_emf), intent(inout), target :: this
  class( t_grid ), intent(in) :: old_lb, new_lb
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  type(t_vdf_report), pointer :: report

  ! reshape b and e fields
  call reshape_copy( this%b, old_lb, new_lb, no_co, send_msg, recv_msg )
  call reshape_copy( this%e, old_lb, new_lb, no_co, send_msg, recv_msg )

  ! if necessary reshape external fields
  if ( this%ext_fld /= p_extfld_none ) then

    ! Reshape the vdfs
    call reshape_nocopy( this%ext_b, new_lb )
    call reshape_nocopy( this%ext_e, new_lb )

    ! Recalculate the fields for static fields (dynamic fields will be recalculated at next
    ! iteration
    if ( this%ext_fld /= p_extfld_static ) then

      ! I need to pass the necessary parameters down here
      ERROR( 'Not implemented yet' )
      call abort_program( p_err_notimplemented )
    endif
  endif

  ! if necessary reshape particle interpolation fields
  ! these do not require copying since they are updated from e and b at each timestep
  if ( this%part_fld_alloc ) then
    call reshape_nocopy( this%b_part, new_lb )
    call reshape_nocopy( this%e_part, new_lb )
  endif

  ! reshape vdf objects if needed

  ! Reshape psi object
  if ( this%psi%data%x_dim_ > 0 ) call reshape_nocopy( this%psi%data, new_lb )

  ! Reshape tavg_data
  report => this%diag%reports
  do
    if ( .not. associated(report) ) exit
    if ( report%tavg_data%x_dim_ > 0 ) then
      call reshape_copy( report%tavg_data, old_lb, new_lb, no_co, send_msg, recv_msg )
    endif
    report => report%next
  enddo

  ! reshape boundary condition data
  call reshape_obj( this%bnd_con, old_lb, new_lb, no_co, send_msg, recv_msg )

  ! store new position on global grid
  this%gix_pos(1:p_x_dim) = new_lb%my_nx( 1, 1:p_x_dim )

end subroutine reshape_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine update_emf_int( this, emf_int )
!-----------------------------------------------------------------------------------------
! Update the spatially centered values
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_emf ), intent(in) :: this
  type( t_vdf ), intent(inout) :: emf_int

  integer :: i1, i2, i3

  ! Note that emf_int must be previously allocated

  select case ( p_x_dim )
    case(1)

      do i1 = lbound( this%e%f1, 2 ) + 1, ubound( this%e%f1, 2 )
         emf_int%f1( 1, i1 ) = 0.5*( this%e%f1( 1, i1 ) + this%e%f1( 1, i1-1 ) )
         emf_int%f1( 2, i1 ) = this%e%f1( 2, i1 )
         emf_int%f1( 3, i1 ) = this%e%f1( 3, i1 )
         emf_int%f1( 4, i1 ) = this%b%f1( 1, i1 )
         emf_int%f1( 5, i1 ) = 0.5*( this%b%f1( 2, i1 ) + this%b%f1( 2, i1-1 ) )
         emf_int%f1( 6, i1 ) = 0.5*( this%b%f1( 3, i1 ) + this%b%f1( 3, i1-1 ) )
      enddo

    case(2)

      do i2 = lbound( this%e%f2, 3 ) + 1, ubound( this%e%f2, 3 )
         do i1 = lbound( this%e%f2, 2 ) + 1, ubound( this%e%f2, 2 )
            emf_int%f2(1, i1, i2) = 0.5*( this%e%f2(1, i1,  i2) + &
                                               this%e%f2(1, i1-1,i2) )
            emf_int%f2(2, i1, i2) = 0.5*( this%e%f2(2, i1, i2  ) + &
                                               this%e%f2(2, i1, i2-1 ) )
            emf_int%f2(3, i1, i2) =       this%e%f2(3, i1, i2 )

            emf_int%f2(4, i1, i2) = 0.5*( this%b%f2(1, i1, i2 ) + &
                                               this%b%f2(1, i1, i2-1) )
            emf_int%f2(5, i1, i2) = 0.5*( this%b%f2(2, i1, i2 ) + &
                                               this%b%f2(2, i1-1, i2 ) )
            emf_int%f2(6, i1, i2) = 0.25*( this%b%f2(3, i1, i2 ) + &
                                                this%b%f2(3, i1, i2-1 ) + &
                                                this%b%f2(3, i1-1, i2 ) + &
                                                this%b%f2(3, i1-1, i2-1) )
         enddo
      enddo

    case(3)

      do i3 = lbound( this%e%f3, 4 ) + 1, ubound( this%e%f3, 4 )
         do i2 = lbound( this%e%f3, 3 ) + 1, ubound( this%e%f3, 3 )
            do i1 = lbound( this%e%f3, 2 ) + 1, ubound( this%e%f3, 2 )
               emf_int%f3(1, i1, i2, i3) = 0.5*( this%e%f3(1, i1, i2, i3) + &
                                                      this%e%f3(1, i1-1, i2, i3) )
               emf_int%f3(2, i1, i2, i3) = 0.5*( this%e%f3(2, i1, i2, i3) + &
                                                      this%e%f3(2, i1, i2-1, i3) )
               emf_int%f3(3, i1, i2, i3) = 0.5*( this%e%f3(3, i1, i2, i3) + &
                                                      this%e%f3(3, i1, i2, i3-1) )

               emf_int%f3(4, i1, i2, i3) = 0.25*( this%b%f3(1, i1, i2, i3) + &
                                                       this%b%f3(1, i1, i2, i3-1) + &
                                                       this%b%f3(1, i1, i2-1, i3) + &
                                                       this%b%f3(1, i1, i2-1, i3-1) )
               emf_int%f3(5, i1, i2, i3) = 0.25*( this%b%f3(2, i1, i2, i3) + &
                                                       this%b%f3(2, i1, i2, i3-1) + &
                                                       this%b%f3(2, i1-1, i2, i3) + &
                                                       this%b%f3(2, i1-1, i2, i3-1) )
               emf_int%f3(6, i1, i2, i3) = 0.25*( this%b%f3(3, i1, i2, i3) + &
                                                       this%b%f3(3, i1, i2-1, i3) + &
                                                       this%b%f3(3, i1-1, i2, i3) + &
                                                       this%b%f3(3, i1-1, i2-1, i3) )
            enddo
         enddo
      enddo

  end select

end subroutine update_emf_int
!-----------------------------------------------------------------------------------------






end module m_emf
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! These had to be moved outside of the module since they are external functions
! This is required to make them overridable type bound procedures while not having to
! define them in the same module as the class definition.
! their interface (as type bound procedures) is defined in emf-define
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine allocate_objs_emf( this )
!-----------------------------------------------------------------------------------------
!       Allocate any objects contained within the EMF (t_emf) object
!-----------------------------------------------------------------------------------------

    use m_emf_define
    use m_emf_solver, only : new_solver_emf

    implicit none

    class( t_emf ), intent(inout) :: this

    ! print *, '[emf] In allocate_objs_emf'

    if ( .not. associated( this%new_solver ) ) then
      this%new_solver => new_solver_emf
    endif

    ! Allocate default t_diag_emf object class
    if ( .not. associated( this%diag ) ) then
      allocate( this%diag )
    endif

    ! Allocate default t_emf_bound object class
    if ( .not. associated( this%bnd_con ) ) then
      allocate( this%bnd_con )
    endif

    if ( .not. associated( this%init_emf ) ) then
      allocate( this%init_emf )
    endif

    if ( .not. associated( this%ext_emf ) ) then
      allocate( this%ext_emf )
    endif

    ! If needed, 'allocate_objs' on any of classes we just made so that any subobjects are created
    call this%diag%allocate_objs()
    !call this%bnd_con%allocate_objs()

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_input_emf( this, input_file, periodic, if_move, grid, dx, dt, gamma )
!-----------------------------------------------------------------------------------------
!       read necessary information from inputdec
!-----------------------------------------------------------------------------------------

  use m_system
  use m_parameters

  use m_emf_define
  use m_emf_bound
  use m_input_file
  use m_fparser
  use m_vdf_smooth
  use m_emf_diag
  use m_emf_solver
  use m_grid_define

  implicit none

  class( t_emf ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ), intent(in)  :: grid
  real(p_double), dimension(:), intent(in) :: dx
  real(p_double), intent(in) :: dt
  real(p_double), intent(in) :: gamma 

  character(len=20) :: solver

  ! smooth
  character(len=20) :: smooth_type
  integer :: smooth_niter, smooth_nmax
 
  character(len=20) :: ext_fld

  ! Initial fields
  character(len=20), dimension(p_f_dim) :: type_init_b ! type of initial b-field
  character(len=20), dimension(p_f_dim) :: type_init_e ! type of initial e-field

  real(p_k_fld), dimension(p_f_dim)   :: init_b0       ! magnitude of initial b_field
  real(p_k_fld), dimension(p_f_dim)   :: init_e0       ! magnitude of initial e-field

  character(len = p_max_expr_len), dimension(p_f_dim) :: init_b_mfunc
  character(len = p_max_expr_len), dimension(p_f_dim) :: init_e_mfunc

  real(p_k_fld), dimension(p_f_dim) :: init_dipole_b_m
  real(p_k_fld), dimension(p_x_dim) :: init_dipole_b_x0
  real(p_k_fld)                     :: init_dipole_b_r0
  real(p_k_fld), dimension(p_f_dim) :: init_dipole_e_p
  real(p_k_fld), dimension(p_x_dim) :: init_dipole_e_x0
  real(p_k_fld)                     :: init_dipole_e_r0

  ! external fields
  character(len=20), dimension(p_f_dim) :: type_ext_b ! type of external b-field
  character(len=20), dimension(p_f_dim) :: type_ext_e ! type of external e-field

  real(p_k_fld), dimension(p_f_dim)   :: ext_b0       ! magnitude of external b_field
  real(p_k_fld), dimension(p_f_dim)   :: ext_e0       ! magnitude of external e-field

  character(len = p_max_expr_len), dimension(p_f_dim) :: ext_b_mfunc
  character(len = p_max_expr_len), dimension(p_f_dim) :: ext_e_mfunc

  real(p_k_fld), dimension(p_f_dim) :: ext_dipole_b_m
  real(p_k_fld), dimension(p_x_dim) :: ext_dipole_b_x0
  real(p_k_fld)                     :: ext_dipole_b_r0
  real(p_k_fld), dimension(p_f_dim) :: ext_dipole_e_p
  real(p_k_fld), dimension(p_x_dim) :: ext_dipole_e_x0
  real(p_k_fld)                     :: ext_dipole_e_r0

  namelist /nl_el_mag_fld/ solver, smooth_type, smooth_niter, smooth_nmax, &
                           type_init_b, type_init_e, init_b0, init_e0, &
                           init_b_mfunc, init_e_mfunc, &
                           init_dipole_b_m, init_dipole_b_x0, init_dipole_b_r0, &
                           init_dipole_e_p, init_dipole_e_x0, init_dipole_e_r0, &
                           ext_fld, type_ext_b, type_ext_e, ext_b0, ext_e0, &
                           ext_b_mfunc, ext_e_mfunc, &
                           ext_dipole_b_m, ext_dipole_b_x0, ext_dipole_b_r0, &
                           ext_dipole_e_p, ext_dipole_e_x0, ext_dipole_e_r0

  integer :: ierr, i

  solver = "yee"

  smooth_type = "none"
  smooth_niter = 1
  smooth_nmax = -1

  ! Initial fields
  type_init_b = "uniform"
  type_init_e = "uniform"
  init_b0     = 0.0_p_k_fld
  init_e0     = 0.0_p_k_fld

  init_b_mfunc = "NO_FUNCTION_SUPPLIED!"
  init_e_mfunc = "NO_FUNCTION_SUPPLIED!"

  init_dipole_b_m  = 0.0_p_k_fld
  init_dipole_b_x0 = 0.0_p_k_fld
  init_dipole_b_r0 = 0.001_p_k_fld

  init_dipole_e_p  = 0.0_p_k_fld
  init_dipole_e_x0 = 0.0_p_k_fld
  init_dipole_e_r0 = 0.001_p_k_fld

  ! external fields
  ext_fld    = "none"
  type_ext_b = "uniform"
  type_ext_e = "uniform"
  ext_b0     = 0.0_p_k_fld
  ext_e0     = 0.0_p_k_fld
  ext_b_mfunc = "NO_FUNCTION_SUPPLIED!"
  ext_e_mfunc = "NO_FUNCTION_SUPPLIED!"

  ext_dipole_b_m  = 0.0_p_k_fld
  ext_dipole_b_x0 = 0.0_p_k_fld
  ext_dipole_b_r0 = 0.001_p_k_fld

  ext_dipole_e_p  = 0.0_p_k_fld
  ext_dipole_e_x0 = 0.0_p_k_fld
  ext_dipole_e_r0 = 0.001_p_k_fld

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_el_mag_fld", ierr )

  if (ierr == 0) then
    read (input_file%nml_text, nml = nl_el_mag_fld, iostat = ierr)
    if (ierr /= 0) then
      print *, "Error reading emf parameters"
      print *, "aborting..."
      stop
    endif
  else
    if (ierr < 0) then
      print *, "Error reading emf parameters"
      print *, "aborting..."
      stop
    else
      if (disp_out(input_file)) then
        SCR_ROOT(" - emf parameters missing, using default")
      endif
    endif
  endif

  call this%allocate_objs()

  ! Create new solver
  call this % new_solver( solver )
  if (.not. associated(this%solver)) then
         print *, "Error reading emf parameters"
         print *, "aborting..."
         stop
      endif

  ! smooth type
  select case(trim(smooth_type))
    case ( "none" )
      this%smooth_type = p_emfsmooth_none
    case ( "stand" )
      this%smooth_type = p_emfsmooth_stand
    case ( "nci" )
      this%smooth_type = p_emfsmooth_nci
    case ( "local" )

      ! Check zero niter
      if (smooth_niter == 0) then
        print *, "Error reading emf parameters"
        print *, "'smooth_niter' must be defined and > 0 for 'local' smooth"
        print *, "aborting..."
      stop

      else

        this%smooth_type = p_emfsmooth_local
        this%smooth_niter = smooth_niter
        this%smooth_nmax = smooth_nmax

      endif

    case default
      print *, "Error reading emf parameters"
      print *, "Invalid value for the smooth parameter"
      print *, "Available smooth types are 'none', 'stand', 'nci' and 'local'"
      print *, "aborting..."
      stop
  end select

  ! Initial fields
   do i=1, p_f_dim
    select case(trim(type_init_b(i)))
      case( "uniform" )
        this%init_emf%type_b(i) = p_emf_uniform
      case( "math func" )
        this%init_emf%type_b(i) = p_emf_math
        this%init_emf%mfunc_expr_b(i) = trim(init_b_mfunc(i))
      case( "dipole" )
        this%init_emf%type_b(i) = p_emf_dipole
        this%init_emf%dipole_b_m  = init_dipole_b_m
        this%init_emf%dipole_b_x0 = init_dipole_b_x0
        this%init_emf%dipole_b_r0 = init_dipole_b_r0
      case default
        print *, "Error reading emf parameters"
        print *, "Invalid value for the type_init_b parameter"
        print *, "Available initial B-field types are 'uniform', 'math func' and 'dipole'"
        print *, "aborting..."
        stop
    end select

    select case(trim(type_init_e(i)))
      case( "uniform" )
        this%init_emf%type_e(i) = p_emf_uniform
      case( "math func" )
        this%init_emf%type_e(i) = p_emf_math
        this%init_emf%mfunc_expr_e(i) = trim(init_e_mfunc(i))
      case( "dipole" )
        this%init_emf%type_e(i)   = p_emf_dipole
        this%init_emf%dipole_e_p  = init_dipole_e_p
        this%init_emf%dipole_e_x0 = init_dipole_e_x0
        this%init_emf%dipole_e_r0 = init_dipole_e_r0
      case default
        print *, "Error reading emf parameters"
        print *, "Invalid value for the type_init_e parameter"
        print *, "Available initial E-field types are 'uniform', 'math func' and 'dipole'"
        print *, "aborting..."
        stop
    end select
  end do

  this%init_emf%uniform_b0        = init_b0
  this%init_emf%uniform_e0        = init_e0


  ! external fields
  select case(trim(ext_fld))
    case ( "none" )
      this%ext_fld = p_extfld_none
    case ( "static" )
      this%ext_fld = p_extfld_static
    case ( "dynamic" )
      this%ext_fld = p_extfld_dynamic
    case default
      print *, "Error reading emf parameters"
      print *, "Invalid value for the ext_fld parameter"
      print *, "Available external field types are 'none', 'static' and 'dynamic'"
      print *, "aborting..."
      stop
  end select

  do i=1, p_f_dim
    select case(trim(type_ext_b(i)))
      case( "none" )
        this%ext_emf%type_b(i) = p_emf_none
      case( "uniform" )
        this%ext_emf%type_b(i) = p_emf_uniform
      case( "math func" )
        this%ext_emf%type_b(i) = p_emf_math
        this%ext_emf%mfunc_expr_b(i) = trim(ext_b_mfunc(i))
      case( "dipole" )
        this%ext_emf%type_b(i) = p_emf_dipole
        this%ext_emf%dipole_b_m  = ext_dipole_b_m
        this%ext_emf%dipole_b_x0 = ext_dipole_b_x0
        this%ext_emf%dipole_b_r0 = ext_dipole_b_r0
      case default
        print *, "Error reading emf parameters"
        print *, "Invalid value for the type_ext_b parameter"
        print *, "Available external B-field types are 'none', 'uniform', 'math func' and 'dipole'"
        print *, "aborting..."
        stop
    end select

    select case(trim(type_ext_e(i)))
      case( "none" )
        this%ext_emf%type_e(i) = p_emf_none
      case( "uniform" )
        this%ext_emf%type_e(i) = p_emf_uniform
      case( "math func" )
        this%ext_emf%type_e(i) = p_emf_math
        this%ext_emf%mfunc_expr_e(i) = trim(ext_e_mfunc(i))
      case( "dipole" )
        this%ext_emf%type_e(i)   = p_emf_dipole
        this%ext_emf%dipole_e_p  = ext_dipole_e_p
        this%ext_emf%dipole_e_x0 = ext_dipole_e_x0
        this%ext_emf%dipole_e_r0 = ext_dipole_e_r0
      case default
        print *, "Error reading emf parameters"
        print *, "Invalid value for the type_ext_e parameter"
        print *, "Available external E-field types are 'none', 'uniform', 'math func' and 'dipole'"
        print *, "aborting..."
        stop
    end select
  end do

  this%ext_emf%uniform_b0        = ext_b0
  this%ext_emf%uniform_e0        = ext_e0

  ! read boundary condtion information
  call this % bnd_con % read_input(input_file, periodic, if_move)

  ! Do not allow PML with local smoothing, for now...
  if (this%smooth_type == p_emfsmooth_local) then
    do i=1, p_x_dim
      if (.not. periodic(i)) then
        if ( (type(this%bnd_con,p_lower,i) == p_bc_vpml) &
              .or. (type(this%bnd_con,p_upper,i) == p_bc_vpml) ) then
          print *, "PML boundary conditions cannot be used with 'local' EMF smoothing."
          print *, "Aborting..."
          stop
        endif
      endif
    end do
  endif

  ! Read field solver parameters, if necessary
  call this % solver % read_input( input_file, dx, dt )

  ! read smoothing information
  call read_nml( this%smooth, input_file )

  ! read diagnostics information
  call this % diag % read_input( input_file, gamma )

end subroutine read_input_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       sets up this data structure from the given information
!-----------------------------------------------------------------------------------------
subroutine init_emf( this, part_grid_center, part_interpolation, &
                        g_space, grid, gc_min, dx, tstep, tmax, &
                        no_co, send_msg, recv_msg, restart, &
                        restart_handle, sim_options )
!-----------------------------------------------------------------------------------------

  use m_system
  use m_parameters
  
  use m_emf_define
  use m_space,     only : t_space, if_move
  use m_node_conf, only : t_node_conf
  use m_restart,   only : t_restart_handle
  use m_grid_define, only : t_grid

  use m_vdf_smooth
  use m_vdf_memory
  use m_vdf_define, only : t_vdf
  use m_vdf_comm, only : t_vdf_msg

  use m_emf_gridval
  use m_emf_diag

  use m_emf_bound
  use m_logprof
  use m_time_step, only : t_time_step, dt, ndump

  implicit none

  class( t_emf ), intent( inout ), target  ::  this

  logical, intent(in) :: part_grid_center
  integer, intent(in) :: part_interpolation
  type( t_space ),     intent(in) :: g_space
  class( t_grid ), intent(in) :: grid
  integer, dimension(:,:), intent(in) :: gc_min
  real(p_double), dimension(:), intent(in) :: dx
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: tmax 
  class( t_node_conf ), intent(in), target :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  type( t_options ), intent(in) :: sim_options

  integer, dimension(2,p_x_dim) :: gc_num_std, gc_num_solver
  integer :: i

  ! Currently smooth information is not set in restart file so set it up first
  call setup( this%smooth, dx, sim_options%gamma )

  ! Initial field values must be setup before the rest of the object
  call this%init_emf%setup( .false., part_interpolation, restart, restart_handle )

  if ( restart ) then

    call this % restart_read( restart_handle )

  else

    this%n_current_emf = 0
    
    ! Store if field values for particle interpolation should be centered on
    ! the cell corner
    this%part_grid_center = part_grid_center

    ! Minimum number of guard cells for particle interpolation
    gc_num_std = gc_min

    ! When using spatially centered values for field interpolation we need an extra
    ! guard cell at the bottom
    if ( this%part_grid_center ) then
      gc_num_std( p_lower, : ) = gc_num_std( p_lower, : ) + 1
    endif

    ! check guard cells are enough for field solver
    call this % solver % min_gc( gc_num_solver )
        do i = 1, p_x_dim
      if ( gc_num_std( p_lower, i ) < gc_num_solver( p_lower, i ) ) then
        gc_num_std( p_lower, i ) = gc_num_solver( p_lower, i )
      endif
      if ( gc_num_std( p_upper, i ) < gc_num_solver( p_upper, i ) ) then
        gc_num_std( p_upper, i ) = gc_num_solver( p_upper, i )
      endif
        enddo

    !  modify guard cell numbers if space moves
    do i = 1, p_x_dim
       if ( if_move(g_space, i) ) then
          gc_num_std(p_lower,i) = gc_num_std(p_lower,i) + 1
       endif
    enddo

    ! Number of guard cells for spatial field smoothing
    if (if_smooth( this%smooth )) then

      ! if no smooth type defined abort
      if (this%smooth_type == p_emfsmooth_none .or. &
          this%smooth_type == p_emfsmooth_nci) then
        print *, "Error setting up emf"
        print *, "smooth_type = 'none' or 'nci', but a smooth section was found"
        print *, "aborting..."
        stop
      endif

      ! add extra guard cells for smoothing
      ! In some situations it would be better to keep the minimum of guard cells
      ! and do an extra communication after the smooth
      gc_num_std(p_lower,:) = gc_num_std(p_lower,:) + smooth_order(this%smooth)
      gc_num_std(p_upper,:) = gc_num_std(p_upper,:) + smooth_order(this%smooth)

    else

      ! if no smooth is specified break if smooth_type /= none
      if (this%smooth_type /= p_emfsmooth_none .and. &
          this%smooth_type /= p_emfsmooth_nci) then
        print *, "Error setting up emf"
        print *, "smooth_type /= none, but no smooth specified"
        print *, "aborting..."
        stop
      endif

      if ( this%smooth_type == p_emfsmooth_nci) then
        ! add extra guard cells for nci smoothing
        gc_num_std(p_lower,:) = gc_num_std(p_lower,:) + 4
        gc_num_std(p_upper,:) = gc_num_std(p_upper,:) + 4
      endif
    endif

    ! create b and e objects
    call this%b%new( p_x_dim, p_f_dim, grid%my_nx(3,:), gc_num_std, dx, .false. )
    call this%e%new( p_x_dim, p_f_dim, grid%my_nx(3,:), gc_num_std, dx, .false. )

    ! set initial values for e and b
    call this%init_emf%set_fld_values( this%e, this%b, g_space, grid%my_nx( p_lower, : ), 0.0d0 )

  endif

  ! setup arrays for external fields if necessary
  if ( this%ext_fld /= p_extfld_none ) then
    ! Only print the first time through this function if using tiles
    if (fsolverev==0) then
      SCR_ROOT(' - setting up external fields...')
    endif

    ! allocate the arrays for the external fields
    call this%ext_b%new(this%b)
    call this%ext_e%new(this%e)

    ! setup the external field object
    call this%ext_emf%setup( (this%ext_fld==p_extfld_dynamic), part_interpolation, restart, restart_handle )

    ! setup the values of static external fields
    ! values of dynamic external fields will be updated at every time step in update_particle_fld
    if ( this%ext_fld == p_extfld_static ) then
         call this%ext_emf%set_fld_values( this%ext_e, this%ext_b, g_space, grid%my_nx( p_lower, : ), 0.0d0 )
    endif

  endif

  ! setup boundary conditions
  call this%bnd_con%init( this%b, this%e, dt(tstep), part_interpolation, &
              g_space, no_co, grid, restart, restart_handle )

  ! setup data that is not on restart file

  ! setup cyl. coordinate data
  this%coordinates = grid%coordinates

  do i = 1, p_x_dim
    this%gix_pos(i) = grid%my_nx( p_lower, i )
  enddo

  ! Perform additional field solver setup
  call this % solver % init( this % gix_pos, this % coordinates )

  ! setup arrays for spatial smoothed / external fields if necessary
  if ( ( this%smooth_type == p_emfsmooth_stand ) .or. &
       ( this%smooth_type == p_emfsmooth_nci ) .or. &
       ( this%ext_fld /= p_extfld_none ) .or. &
       ( this%part_grid_center ) ) then

    this%part_fld_alloc = .true.

    ! allocate memory for particle field interpolation arrays
    call alloc( this%b_part )
    call alloc( this%e_part )

    call this%b_part%new(this%b)
    call this%e_part%new(this%e)

  else

    ! Just point to the main emf data
    this%b_part => this%b
    this%e_part => this%e

    this%part_fld_alloc = .false.
  endif

  if (fsolverev==0) then
    fsolverev   = create_event('field solver')
    fsmoothev   = create_event('field smooth')
    getpsi_ev   = create_event('psi calculation')
  endif

  ! Setup diagnostics
  call this % diag % init( this%ext_fld == p_extfld_none, this%part_fld_alloc, part_interpolation )

end subroutine init_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine cleanup_emf( this )
!-----------------------------------------------------------------------------------------

  use m_emf_define
  use m_vdf
  use m_vdf_memory
  use m_vdf_smooth
  use m_emf_bound
  use m_emf_psi
  use m_emf_diag
  use m_emf_gridval

  implicit none

  class( t_emf ), intent( inout )  ::  this

  ! If e_part and b_part are using memory (and not just pointing to e and b) clean
  ! them up
  if ( this%part_fld_alloc ) then
    call this%b_part%cleanup()
    call this%e_part%cleanup()

    call freemem( this%b_part )
    call freemem( this%e_part )
  endif
  nullify( this%b_part, this%e_part )

  call this%b%cleanup()
  call this%e%cleanup()

  call cleanup( this%smooth )

  call this%init_emf%cleanup()
  deallocate( this%init_emf )

  call this%ext_emf%cleanup()
  deallocate( this%ext_emf )

  if ( this%ext_fld /= p_extfld_none ) then
    call this%ext_b%cleanup()
    call this%ext_e%cleanup()
  endif

  call this % diag % cleanup()
  deallocate( this%diag )

  call cleanup( this%psi )

  ! cleanup data from boundary conditions
  call this%bnd_con%cleanup()
  deallocate( this%bnd_con )

  ! cleanup solver data if necessary
  call this % solver % cleanup()
  deallocate( this%solver )

end subroutine cleanup_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_checkpoint_emf( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       write object information into a restart file
!-----------------------------------------------------------------------------------------

  use m_parameters
  use m_emf_define
  use m_restart

  use m_emf_gridval
  use m_vdf
  use m_emf_bound

  implicit none

  class( t_emf ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  ! local variables
  character(len=*), parameter :: err_msg = 'error writing restart data for emf object.'
  integer :: ierr

  ! Write restart data for initial fields (this must happen before the rest of the object)
  call this%init_emf%write_checkpoint( restart_handle )

  restart_io_wr( p_emf_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%n_current_emf, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%part_grid_center, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  call this%b%write_checkpoint( restart_handle )
  call this%e%write_checkpoint( restart_handle )

  restart_io_wr( this%ext_fld, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! write restart data for external fields
  if ( this%ext_fld /= p_extfld_none ) then
    call this%ext_emf%write_checkpoint( restart_handle )
  endif

  ! write restart data for boundary conditions
  call this%bnd_con%write_checkpoint( restart_handle )

end subroutine write_checkpoint_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_read_emf( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       read object information from a restart file
!-----------------------------------------------------------------------------------------
  use m_parameters
  use m_emf_define
  use m_restart

  use m_emf_gridval
  use m_vdf
  use m_emf_bound

  implicit none

  class( t_emf ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for emf object.'
  character(len=len(p_emf_rst_id)) :: rst_id
  integer :: ierr

  ! init_emf checkpoint is called in this % init, even though it's called as normally
  ! in write_checkpoint_emf

  DEBUG('restart read emf')
  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  ! check if restart file is compatible
  if ( rst_id /= p_emf_rst_id) then
    ERROR('Corrupted restart file, or restart file from incompatible binary (emf)')
    call abort_program(p_err_rstrd)
  endif

  restart_io_rd( this%n_current_emf, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  DEBUG('before reading this%part_grid_center')
  restart_io_rd( this%part_grid_center, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  DEBUG('before reading this%b')
  call this%b%read_checkpoint( restart_handle )
  DEBUG('before reading this%e')
  call this%e%read_checkpoint( restart_handle )

  DEBUG('before reading this%ext_fld')
  restart_io_rd( this%ext_fld, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

end subroutine restart_read_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine move_window_emf( this, g_space, nx_p_min, need_fld_val )
!-----------------------------------------------------------------------------------------
!       move boundaries of the electro-magnetic field
!-----------------------------------------------------------------------------------------
  use m_parameters
  use m_vdf
  use m_vdf_define
  use m_emf_define
  use m_space

  use m_emf_gridval

  implicit none

  class( t_emf ), intent( inout )  ::  this
  type( t_space ),     intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  logical, intent(in) :: need_fld_val

  integer, dimension( p_max_dim ) :: lb

  if (nx_move( g_space, 1 ) > 0) then

    ! move window for e and b field
    call move_window( this%b, g_space )
    call move_window( this%e, g_space )


    if ( need_fld_val ) then
      lb(1) = this%b%nx_(1) - nx_move( g_space, 1 )
      lb(2:p_x_dim ) = 1 - this%b%gc_num_( p_lower, 2:p_x_dim )

      call this%init_emf%set_fld_values( this%e, this%b, g_space, nx_p_min, 0.0d0, &
                                lbin = lb )
    endif

    ! move window for external fields: this is only required for static external fields,
    ! dynamic fields are recalculated at every iteration
    if ( this%ext_fld == p_extfld_static  ) then
       ! shift local data
       call move_window( this%ext_b, g_space )
       call move_window( this%ext_e, g_space )

       lb(1) = this%b%nx_(1) + this%b%gc_num_( p_upper, 1 ) - nx_move( g_space, 1 )
       lb(2:p_x_dim ) = 1 - this%b%gc_num_( p_lower, 2:p_x_dim )

       call this%ext_emf%set_fld_values( this%ext_e, this%ext_b, g_space, nx_p_min, 0.0d0, &
                               lbin = lb )
    endif

    ! move window for boundary condition data
    ! (namely lindman boundaries)
    call this%bnd_con%move_window( g_space )

  endif

end subroutine move_window_emf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!       update boundaries of the electro-magnetic field
!-----------------------------------------------------------------------------------------
subroutine update_boundary_emf( this, g_space, no_co, send_msg, recv_msg )

  use m_parameters
  use m_emf_define
  use m_space
  use m_node_conf
  use m_vdf_comm, only : t_vdf_msg

  implicit none

  class( t_emf ), intent( inout )  ::  this

  type( t_space ),     intent(in) :: g_space
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  ! update boundaries on emf fields
  call this%bnd_con%update_boundary( this%b, this%e,  g_space, no_co, send_msg, recv_msg )

  ! no update is required for external fields since these are constant or calculated at every time
  ! step

end subroutine update_boundary_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine update_particle_fld_emf( this, n, t, dt, space, nx_p_min )
!-----------------------------------------------------------------------------------------
! Update fields to be used by particles
!-----------------------------------------------------------------------------------------
  use m_parameters
  use m_logprof
  use m_emf_define
  use m_space
  use m_emf_gridcenter
  use m_emf_ncifilter
  use m_emf_gridval
  use m_vdf
  use m_vdf_define
  use m_vdf_smooth
  use m_vdf_math
  use m_vdf

  implicit none

  class( t_emf ), intent( inout )  ::  this
  integer, intent(in) :: n
  real( p_double ), intent(in) :: t, dt
  type( t_space ), intent(in) :: space
  integer, dimension(:), intent(in) :: nx_p_min

  ! Smoothing and external fields

  if ( this%ext_fld == p_extfld_dynamic ) then
    call this%ext_emf%set_fld_values( this%ext_e, this%ext_b, space, nx_p_min, t )
  endif

  ! do local smooth of fields if necessary
  if ( this%smooth_type == p_emfsmooth_local ) then
    ! Stop smoothing after given iteration
    if (this%smooth_nmax > 0 .and. n >= this%smooth_nmax) then
      this%smooth_type = p_emfsmooth_none
    else if ( mod( n, this%smooth_niter ) == 0) then

      ! smooth e and b fields
      call begin_event( fsmoothev )

      call smooth( this%e, this%smooth )
      call smooth( this%b, this%smooth )

      call end_event( fsmoothev )
    endif
  endif

  ! Update fields for particle interpolation
  if ( this%part_fld_alloc ) then

    ! Copy current emf field values to e_part and b_part
    this%e_part = this%e
    this%b_part = this%b

    ! smooth e and b fields using standard smoothing
    select case ( this%smooth_type )
    case( p_emfsmooth_stand )
      call begin_event( fsmoothev )

      call smooth( this%e_part, this%smooth )
      call smooth( this%b_part, this%smooth )

      call end_event( fsmoothev )

    case( p_emfsmooth_nci )
      call begin_event( fsmoothev )
      call filter_nci( this%e_part, this%b_part, dt, this%part_grid_center )
      call end_event( fsmoothev )

    end select

    ! Add external fields (this could be optimized in the situations where the external fields
    ! are constant - avoid reading ext_e0, ext_b0)
    if ( this%ext_fld /= p_extfld_none ) then
      call add( this%e_part, this%ext_e )
      call add( this%b_part, this%ext_b )
    endif

    ! Grid center fields if necessary
    if ( this%part_grid_center ) then
      call grid_center( this )
    endif

  endif

end subroutine update_particle_fld_emf
!-----------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Used for CUDA (and similar) situations where data may have to be copied from another location.
!---------------------------------------------------------------------------------------------------
subroutine fill_data_emf(this, tstep)
  use m_time_step
  use m_emf_define
  use m_parameters
  implicit none
  class( t_emf ), intent(inout)     :: this
  type( t_time_step ),   intent(in) :: tstep

  ! If used, the general format for this function would be along the lines of:
  !
  ! if( n(tstep) /= this%n_current_emf) then
  !   ! get the data needed
  !   call ....
  !   ! update the status counter so that the data is only copied as needed
  !   this%n_current_emf = n(tstep)
  ! endif

  ! put a fatal error here in case any subclasses forgot to increment 'n_current_emf'
  !   this error is most likely from a subclass not having the line
  !      "this%n_current_emf=this%n_current_emf+1" 
  !   after the field solver (e.g. 'advance_emf'). Ask Ricardo or Anton or Adam if needed.
  if( n(tstep) /= this%n_current_emf) then
    ERROR( "Error in fill_data_emf.. see 'fill_data_emf' is os-emf.f03")
    call abort_program(p_err_rstrd)
  endif

end subroutine
!-----------------------------------------------------------------------------------------
