!  update_pgc_part - updates grad |a|^2 (n), grad |a|^2 (n+1/2), |a|^2 (n-1/2)
!  Must be called after advance a, before calling this routine


! m_species_pgc module
!
! Handles particle advance / deposit for the ponderomotive guiding center algorithm
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_meta

#include "memory/memory.h"

use m_emf_define
use m_vdf
use m_vdf_memory

use m_node_conf

use m_space
use m_define_current
use m_define_grid

use m_system
use m_parameters

use m_file_system

use m_vdf_comm

use m_restart

implicit none

private

! string to id restart data
character(len=*), parameter :: p_emf_meta_rst_id = "emf_meta rst data - 0x0001"


!-------------------------------------------------------------------------------
! t_meta_gridval class definition
!-------------------------------------------------------------------------------
! Parameters for defining metamaterial on the grid
!-------------------------------------------------------------------------------
type :: t_meta_gridval

  ! Type of grid values to set (math, uniform, ...)
  integer :: type_meta_params

  ! variables for uniform fields
  real(p_k_fld), dimension(p_f_dim) :: uniform_wpe = 0.0_p_k_fld
  real(p_k_fld), dimension(p_f_dim) :: uniform_w0e = 0.0_p_k_fld
  real(p_k_fld), dimension(p_f_dim) :: uniform_Ge = 0.0_p_k_fld
  real(p_k_fld), dimension(p_f_dim) :: uniform_wpm = 0.0_p_k_fld
  real(p_k_fld), dimension(p_f_dim) :: uniform_w0m = 0.0_p_k_fld
  real(p_k_fld), dimension(p_f_dim) :: uniform_Gm = 0.0_p_k_fld

  ! variables for math function fields
  character(len = p_max_expr_len),    dimension(p_f_dim) :: mfunc_expr_wpe = " "
  type(t_fparser),                    dimension(p_f_dim) :: mfunc_wpe
  character(len = p_max_expr_len),    dimension(p_f_dim) :: mfunc_expr_w0e = " "
  type(t_fparser),                    dimension(p_f_dim) :: mfunc_w0e
  character(len = p_max_expr_len),    dimension(p_f_dim) :: mfunc_expr_Ge = " "
  type(t_fparser),                    dimension(p_f_dim) :: mfunc_Ge
  character(len = p_max_expr_len),    dimension(p_f_dim) :: mfunc_expr_wpm = " "
  type(t_fparser),                    dimension(p_f_dim) :: mfunc_wpm
  character(len = p_max_expr_len),    dimension(p_f_dim) :: mfunc_expr_w0m = " "
  type(t_fparser),                    dimension(p_f_dim) :: mfunc_w0m
  character(len = p_max_expr_len),    dimension(p_f_dim) :: mfunc_expr_Gm = " "
  type(t_fparser),                    dimension(p_f_dim) :: mfunc_Gm

  ! math function defining metamaterial domain
  character(len = p_max_expr_len) :: mfunc_expr_meta_domain = " "
  type(t_fparser)                 :: mfunc_meta_domain

  ! metamaterial grid domain
  type( t_vdf ) :: meta_grid_domain ! if meta_grid_domain(i,j,k)==1, then we are inside the metamaterial
								    ! if meta_grid_domain(i,j,k)==0, then we are outside the metamaterial


end type t_meta_gridval

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! t_emf_hard_point_source class definition
!-------------------------------------------------------------------------------
! Parameters for defining hard point source
!-------------------------------------------------------------------------------
type :: t_emf_hard_point_source

  integer, dimension(p_x_dim) :: inj_pos_i  ! hard_point_source injection index
  integer                     :: inj_fc     ! field component to inject hard point source

  ! variables for math function fields
  character(len = p_max_expr_len) :: mfunc_expr_hard_point_source = " "
  type(t_fparser)                 :: mfunc_hard_point_source

end type t_emf_hard_point_source
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! t_emf_meta class definition
!-------------------------------------------------------------------------------

type, extends( t_emf ) :: t_emf_meta

  ! magnetic/electric polarization fields and currents
  type( t_vdf ) :: Pm, Jpm, Pe, Jpe

  ! incident b and e fields for TFSF simulations
  type( t_vdf ) :: b_incident, e_incident

  ! metamaterial parameter grid values
  type( t_vdf ) :: wpesq, w0esq, Ge, wpmsq, w0msq, Gm

  ! metamaterial parameters
  type( t_meta_gridval ) :: metamaterial

  ! cell indices of lower left and upper right corners of TFSF box
  logical :: use_tfsf = .false.
  integer, dimension(2, p_x_dim) :: tfsf_box_i

  ! Hard Point Source
  logical :: use_hard_point_source = .false.
  type( t_emf_hard_point_source )  ::  hard_point_source

contains

    procedure :: read_input => read_input_pgc
    procedure :: init       => init_pgc
    procedure :: advance    => advance_pgc
    procedure :: report     => report_diag_pgc
    procedure :: cleanup => cleanup_pgc

    procedure :: write_checkpoint => write_checkpoint_pgc

	! Diagnostics
	procedure :: avail_report_quants => avail_report_quants_emf_pgc
	procedure :: init_report_quants => init_report_quants_emf_pgc

    ! local methods
    procedure :: reset_chi => reset_chi_pgc
    procedure, nopass :: get_emf_pgc_2D

end type t_emf_meta

character(len=31), dimension(5), parameter :: p_report_quants_meta = &
   'e1_inc    ', 'e2_inc    ', 'e3_inc    ', 'b1_inc    ', 'b2_inc    ', 'b3_inc    ', &
   'vbound    ', &
   'wpesq1    ', 'wpesq2    ', 'wpesq3    ', 'wpmsq1    ', 'wpmsq2    ', 'wpmsq3    ', &
   'w0esq1    ', 'w0esq2    ', 'w0esq3    ', 'w0msq1    ', 'w0msq2    ', 'w0msq3    ', &
   'Ge1       ', 'Ge2       ', 'Ge3       ', 'Gm1       ', 'Gm2       ', 'Gm3       ', &
   'Pe1       ', 'Pe2       ', 'Pe3       ', 'Pm1       ', 'Pm2       ', 'Pm3       '/)

! Meta specific diagnostics

integer, parameter :: p_e1_inc = 1 , p_e2_inc = 2, p_e3_inc = 3
integer, parameter :: p_b1_inc = 4 , p_b2_inc = 5, p_b3_inc = 6
integer, parameter :: p_virtualbound = 7
integer, parameter :: p_wpesq1 = 8, p_wpesq2 = 9, p_wpesq3 = 10
integer, parameter :: p_wpmsq1 = 11, p_wpmsq2 = 12, p_wpmsq3 = 13
integer, parameter :: p_w0esq1 = 14, p_w0esq2 = 15, p_w0esq3 = 16
integer, parameter :: p_w0msq1 = 17, p_w0msq2 = 18, p_w0msq3 = 19
integer, parameter :: p_Ge1 = 20, p_Ge2 = 21, p_Ge3 = 22
integer, parameter :: p_Gm1 = 23, p_Gm2 = 24, p_Gm3 = 25
integer, parameter :: p_Pe1 = 26, p_Pe2 = 27, p_Pe3 = 28
integer, parameter :: p_Pm1 = 29, p_Pm2 = 30, p_Pm3 = 31


!-----------------------------------------------------------------------------------------
contains

!-----------------------------------------------------------------------------------------
! Read input file data for Metamaterial algorithm
!   - this routine also calls the superclass (t_emf) read_nml routine to read that section
!     first
!-----------------------------------------------------------------------------------------
subroutine read_input_meta( this, input_file, periodic, if_move, grid )

   implicit none

   class( t_emf_pgc ), intent( inout )  ::  this
   class( t_input_file ), intent(inout) :: input_file
   logical, dimension(:), intent(in) :: periodic, if_move
   class( t_grid ), intent(in) :: grid

   real(p_k_fld) :: w0 , tau , omega , lon_center , per_center , a0 , per_focus
   logical :: free_stream

   integer :: ierr

   namelist /nl_pgc/ w0 , tau , omega , lon_center , per_focus , &
                     per_center , a0 , free_stream

   print *, '[pgc] in read_input_pgc'

   ! Add PGC specific report quantities
   call this % init_report_quants( )

   ! first read the superclass (emf) section
   call this % t_emf % read_input( input_file, periodic, if_move )

   ! Now read this class input section

   ! default values
   w0          = 1.0
   tau         = 1.0
   omega       = 1.0
   lon_center  = 0.0
   per_focus   = 0.0
   per_center  = 0.0
   a0          = 0.0
   free_stream = .false.

   ! read values from file
   call get_namelist( input_file, "nl_pgc", ierr )

   if (ierr == 0) then
	 read (input_file%nml_text, nml = nl_pgc, iostat = ierr)
	 if (ierr /= 0) then
	   print *, "Error reading pgc parameters"
	   print *, "aborting..."
	   stop
	 endif
   else
	 SCR_ROOT(" - no pgc parameters specified")
     stop
   endif

   this%w0          = w0
   this%tau         = tau
   this%omega       = omega
   this%lon_center  = lon_center
   this%per_focus   = per_focus
   this%per_center  = per_center
   this%a0          = a0
   this%free_stream = free_stream

end subroutine read_input_meta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize emf_meta object
!   - this routine also calls the superclass (t_emf) setup routine first
!-----------------------------------------------------------------------------------------
subroutine init_meta( this, part_grid_center, part_interpolation, &
                      g_space, grid, gc_min, dx, dt, &
					  no_co, restart, restart_handle, sim_options )

  use m_emf

  implicit none

  class( t_emf_pgc ), intent( inout ), target  ::  this

  logical, intent(in) :: part_grid_center
  integer, intent(in) :: part_interpolation
  type( t_space ),     intent(in) :: g_space
  class( t_grid ), intent(in) :: grid
  integer, dimension(:,:), intent(in) :: gc_min
  real( p_double ),   intent(in) :: dt
  real( p_double ), dimension(:), intent(in) :: dx
  class( t_node_conf ), intent(in) :: no_co
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  type( t_options ), intent(in) :: sim_options


  real(p_k_fld), dimension(p_x_dim) ::  g_xmin

  print *, '[meta] in init_meta'

  ! setup superclass data
  call this % t_emf % init( part_grid_center, part_interpolation, &
              g_space, grid, gc_min, dx, dt, &
			  no_co, restart, restart_handle, sim_options )

  ! diagnostics have an additional setup to handle Meta specific reports
  call setup_diag_meta( this )

  !executable statements

  ! Set the minimum number of guard cells for the Meta field solver
  do i = 1, p_x_dim
	if ( gc_num_std( p_lower, i ) < 1 ) gc_num_std( p_lower, i ) = 1
	if ( gc_num_std( p_upper, i ) < 2 ) gc_num_std( p_upper, i ) = 2
  enddo


  if ( restart ) then

	call restart_read_meta( this, restart_handle )

  else

	! create 1D incident b and e fields if using tfsf solver
	if ( this%use_tfsf ) then
	   call new(this%b_incident, 1, p_f_dim, grid%my_nx(3,:), gc_num_std, dx )
	   call new(this%e_incident, 1, p_f_dim, grid%my_nx(3,:), gc_num_std, dx )

	   this%b_incident = 0.0_p_k_fld
	   this%e_incident = 0.0_p_k_fld
	endif

	! create electric/magnetic polarization field and current grids if using meta solver
	call new(this%Jpm, p_x_dim, p_f_dim, grid%my_nx(3,:), gc_num_std, dx)
	call new(this%Jpe, p_x_dim, p_f_dim, grid%my_nx(3,:), gc_num_std, dx)
	call new(this%Pm, p_x_dim, p_f_dim, grid%my_nx(3,:), gc_num_std, dx)
	call new(this%Pe, p_x_dim, p_f_dim, grid%my_nx(3,:), gc_num_std, dx)

	! set initial values
	this%Jpm = 0.0_p_k_fld
	this%Jpe = 0.0_p_k_fld
	this%Pm = 0.0_p_k_fld
	this%Pe = 0.0_p_k_fld

  endif

  ! setup hard point source
  if ( this%use_hard_point_source ) then
	SCR_ROOT('setting up hard point source...')

	! setup the hard point source object
	call setup( this%hard_point_source )

	SCR_ROOT('hard point source ok')
  endif

  ! setup arrays for metamaterial parameters if using meta solver
  SCR_ROOT('setting up metamaterial parameters...')

  ! allocate the arrays for metamaterial parameters
  call new(this%wpesq, this%e)
  call new(this%w0esq, this%e)
  call new(this%Ge, this%e)
  call new(this%wpmsq, this%b)
  call new(this%w0msq, this%b)
  call new(this%Gm, this%b)

  call new(this%metamaterial%meta_grid_domain, this%e)

  ! setup the metamaterial object
  call setup( this%metamaterial )

  ! setup metamaterial parameters on the grid
  SCR_ROOT('Setting metamaterial values')
  call set_meta_values( this%metamaterial, this%wpesq, this%w0esq, this%Ge, this%wpmsq, this%w0msq, this%Gm, g_space, nx_p_min )

  ! correct metamaterial boundaries
  SCR_ROOT('Correcting metamaterial boundaries')
  call correct_meta_bounds( this%metamaterial, this%wpesq, this%w0esq, this%Ge, this%wpmsq, this%w0msq, this%Gm)


  SCR_ROOT('metamaterial ok')

  ! store information in bnd_con data structure that TF/SF method is being used
  if ( this%use_tfsf ) then
     this%bnd_con%if_tfsf = .true.
     call setup_incident_bound( this%bnd_con, this%b_incident, this%e_incident, dt, &
			  g_space, no_co, grid%coordinates, restart, restart_handle )
  endif

end subroutine init_meta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleanup emf_meta object
!-----------------------------------------------------------------------------------------
subroutine cleanup_meta( this )
!-----------------------------------------------------------------------------------------

  use m_emf, only : cleanup

  implicit none

  class(t_emf_meta) , intent(inout) :: this

  ! cleanup superclass
  call this % t_emf % cleanup()

  ! cleanup local data structures
  call cleanup( this%Jpm )
  call cleanup( this%Jpe )
  call cleanup( this%Pm )
  call cleanup( this%Pe )

  call cleanup( this%metamaterial )
  call cleanup( this%wpesq )
  call cleanup( this%w0esq )
  call cleanup( this%Ge )
  call cleanup( this%wpmsq )
  call cleanup( this%w0msq )
  call cleanup( this%Gm )

  ! cleanup hard point source data
  call cleanup(this%hard_point_source)

end subroutine cleanup_meta
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine write_checkpoint_meta( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       write object information into a restart file
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_emf_meta ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for emf_meta object.'
  integer :: ierr

  print *, '[meta] in write_checkpoint_meta'

  ! write superclass checkpoint data first
  call this % t_emf % write_checkpoint( restart_handle )

  ! write checkpoint id
  restart_io_wr( p_emf_pgc_meta_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! Write restart data for ...
 ! call restart_write( this%a_nm, restart_handle )
 ! call restart_write( this%a_n, restart_handle )
 ! call restart_write( this%a_np, restart_handle )

end subroutine write_checkpoint_pgc
!-----------------------------------------------------------------------------------------


subroutine restart_read_meta( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       read object information from a restart file
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_emf_pgc ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=len(p_emf_meta_rst_id)) :: rst_id
  integer :: ierr

  restart_io_rd( rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then
	ERROR('error reading restart data for emf object.')
	call abort_program(p_err_rstrd)
  endif

  ! check if restart file is compatible
  if ( rst_id /= p_emf_meta_rst_id) then
 ERROR('Corrupted restart file, or restart file')
 ERROR('from incompatible binary (emf_meta)')
	call abort_program(p_err_rstrd)
  endif

  ! read restart data
  ! call restart_read( this%a_nm, restart_handle )
  ! call restart_read( this%a_n,  restart_handle )
  ! call restart_read( this%a_np, restart_handle )


end subroutine restart_read_pgc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! These allow the list of report quantities to be expanded / replaced by subclasses
!-----------------------------------------------------------------------------------------
function avail_report_quants_emf_meta( this )

  class( t_emf_meta ), intent(in) :: this
  integer :: avail_report_quants_emf_meta

  integer :: n

  ! get the reports available on the superclass
  n = this % t_emf % avail_report_quants( )

  ! add the additional local reports
  n = n + size( p_report_quants_meta )

  avail_report_quants_emf_meta = n

end function avail_report_quants_emf_meta


subroutine init_report_quants_emf_meta( this )

  class( t_emf_meta ), intent(inout) :: this
  integer :: n

  if ( .not. associated( this % diag % report_quants ) ) then
    allocate( this % diag % report_quants( this % avail_report_quants( ) ) )
  endif

  n = this % t_emf % avail_report_quants(  )

  ! add the additional reports
  this % diag % report_quants( n+1 : n+size(p_report_quants_meta) ) = p_report_quants_meta

end subroutine init_report_quants_emf_meta

!-----------------------------------------------------------------------------------------
subroutine setup_diag_meta( this )

  use m_logprof

  implicit none

  class( t_emf_meta ), intent(inout) :: this
  type(t_vdf_report), pointer :: report
  integer, parameter :: izero = ichar('0')

  integer :: quant

  ! Normal reports
  report => this%diag%reports
  do
    if ( .not. associated( report ) ) exit

    quant = report%quant - this % t_emf % avail_report_quants( )

	select case ( quant )

      case ( p_Pe1, p_Pe2, p_Pe3 )
        report%label = 'Pe_'//char(izero + 1 + report%quant - p_Pe1)
        report%units = 'm_e c \omega_p e^{-1}'

      case ( p_Pm1, p_Pm2, p_Pm3 )
        report%label = 'Pm_'//char(izero + 1 + report%quant - p_Pm1)
        report%units = 'm_e c \omega_p e^{-1}'

      case ( p_e1_inc, p_e2_inc, p_e3_inc )
        report%label = 'E_'//char(izero + 1 + report%quant - p_e1_inc)//'^{inc}'
        report%units = 'm_e c \omega_p e^{-1}'

      case ( p_b1_inc, p_b2_inc, p_b3_inc )
        report%label = 'B_'//char(izero + 1 + report%quant - p_b1_inc)//'^{inc}'
        report%units = 'm_e c \omega_p e^{-1}'

      case ( p_virtualbound )
        report%label    = 'Fields@virtualbound'
        report%units    = 'm_e c \omega_p e^{-1}'

      case ( p_wpesq1, p_wpesq2, p_wpesq3 )
        report%label = '\omega_{pe}_{'//char(izero + 1 + report%quant - p_wpesq1) //'}^2'
        report%units = '\omega_p^2'

      case ( p_wpmsq1, p_wpmsq2, p_wpmsq3 )
        report%label = '\omega_{pm}_{'//char(izero + 1 + report%quant - p_wpmsq1) //'}^2'
        report%units = '\omega_p^2'

      case ( p_w0esq1, p_w0esq2, p_w0esq3 )
        report%label = '\omega_{0e}_{'//char(izero + 1 + report%quant - p_w0esq1) //'}^2'
        report%units = '\omega_p^2'

      case ( p_w0msq1, p_w0msq2, p_w0msq3 )
        report%label = '\omega_{0m}_{'//char(izero + 1 + report%quant - p_w0msq1) //'}^2'
        report%units = '\omega_p^2'

      case ( p_Ge1, p_Ge2, p_Ge3 )
        report%label = '\Gamma_{e}_'//char(izero + 1 + report%quant - p_Ge1)
        report%units = '\omega_p'

      case ( p_Gm1, p_Gm2, p_Gm3 )
        report%label = '\Gamma_{m}_'//char(izero + 1 + report%quant - p_Gm1)
        report%units = '\omega_p'

	  case default
		! must be a superclass diagnostic, ignore
		! print *, 'In setup_diag_meta, unknown report = ', report % name
		continue

	end select

    ! process next report
    report => report % next
  enddo

end subroutine setup_diag_meta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       report on electro-magnetic field - diagnostic
!-----------------------------------------------------------------------------------------
subroutine report_diag_meta( this, g_space, grid, no_co, tstep, t )
!-----------------------------------------------------------------------------------------

  use m_time_step

  ! (*debug*)
  ! use m_vpml

  use m_vdf_report

  implicit none

  class( t_emf_pgc ),                 intent(inout) :: this

  type( t_space ),     intent(in) :: g_space
  class( t_grid ),      intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double),      intent(in) :: t

  type( t_vdf_report ), pointer :: rep

  integer :: quant

  print *, '[pgc] in report_diag_pgc'

  ! run superclass diagnostics
  call this % t_emf % report( g_space, grid, no_co, tstep, t )

  print *, '[pgc] after calling superclass report'

  ! run PGC specific diagnostics
  rep => this%diag%reports
  do
    if ( .not. associated( rep ) ) exit

    if ( if_report( rep, tstep ) ) then

	   quant = rep%quant - this % t_emf % avail_report_quants( )

	   select case ( quant )

         case ( p_Pe1, p_Pe2, p_Pe3 )
           call report_vdf( rep, emf%Pe, rep%quant - p_Pe1 + 1, g_space, grid, no_co, tstep, t )

         case ( p_Pm1, p_Pm2, p_Pm3 )
           call report_vdf( rep, emf%Pm, rep%quant - p_Pm1 + 1, g_space, grid, no_co, tstep, t )

         case ( p_e1_inc, p_e2_inc, p_e3_inc )
           call report_vdf( rep, emf%e_incident, rep%quant - p_e1_inc + 1, g_space, grid, no_co, tstep, t )

         case ( p_b1_inc, p_b2_inc, p_b3_inc )
           call report_vdf( rep, emf%b_incident, rep%quant - p_b1_inc + 1, g_space, grid, no_co, tstep, t )

         case ( p_virtualbound )
           ! virtual boundary around tfsf box
           call report_4vdfs( rep, emf%e, emf%b, emf%e_incident, emf%b_incident, emf%tfsf_box_i, g_space, grid, no_co, tstep, t )

         case ( p_wpesq1, p_wpesq2, p_wpesq3 )
           call report_vdf( rep, emf%wpesq, rep%quant - p_wpesq1 + 1, g_space, grid, no_co, tstep, t  )

         case ( p_wpmsq1, p_wpmsq2, p_wpmsq3 )
           call report_vdf( rep, emf%wpmsq, rep%quant - p_wpmsq1 + 1, g_space, grid, no_co, tstep, t  )

         case ( p_w0esq1, p_w0esq2, p_w0esq3 )
           call report_vdf( rep, emf%w0esq, rep%quant - p_w0esq1 + 1, g_space, grid, no_co, tstep, t  )

         case ( p_w0msq1, p_w0msq2, p_w0msq3 )
           call report_vdf( rep, emf%w0msq, rep%quant - p_w0msq1 + 1, g_space, grid, no_co, tstep, t  )

         case ( p_Ge1, p_Ge2, p_Ge3 )
           call report_vdf( rep, emf%Ge, rep%quant - p_Ge1 + 1, g_space, grid, no_co, tstep, t  )

         case ( p_Gm1, p_Gm2, p_Gm3 )
           call report_vdf( rep, emf%Gm, rep%quant - p_Gm1 + 1, g_space, grid, no_co, tstep, t  )

		 case default
		   ! must be a superclass diagnostic, ignore
		   continue

		end select

    endif

    rep => rep%next
  enddo

  ! (*debug*)
  ! call report( emf%bnd_con%vpml_all, g_space, grid, no_co, tstep, t, emf%e%dx )

end subroutine report_diag_meta
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine move_window_meta( this, g_space, nx_p_min, need_fld_val )
!-----------------------------------------------------------------------------------------
!       move boundaries of the electro-magnetic field
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_emf_meta ), intent( inout )  ::  this
  type( t_space ),     intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  logical, intent(in) :: need_fld_val

  integer, dimension( p_max_dim ) :: lb

  ! call superclass method
  call this % t_emf % move_window( this, g_space, nx_p_min, need_fld_val )

  ! move window for metamaterial grids
  if (nx_move( g_space, 1 ) > 0) then

	 call move_window( this%Jpm, g_space )
	 call move_window( this%Jpe, g_space )
	 call move_window( this%Pm, g_space )
	 call move_window( this%Pe, g_space )

	 call move_window( this%wpesq, g_space )
	 call move_window( this%w0esq, g_space )
	 call move_window( this%Ge, g_space )
	 call move_window( this%wpmsq, g_space )
	 call move_window( this%w0msq, g_space )
	 call move_window( this%Gm, g_space )

	 ! initialize the field values where required
	 ! note that we do not use values from another node, we simply recalculate them (which
	 ! ends up being faster that having the extra communication)

	 lb(1) = this%b%nx(1) + this%b%gc_num( p_upper, 1 ) - nx_move( g_space, 1 )
	 lb(2:p_x_dim ) = 1 - this%b%gc_num( p_lower, 2:p_x_dim )

	 call set_meta_values( this%metamaterial, this%wpesq, this%w0esq, this%Ge, &
											  this%wpmsq, this%w0msq, this%Gm, g_space, nx_p_min, lbin = lb )

	 call correct_meta_bounds( this%metamaterial, this%wpesq, this%w0esq, this%Ge, &
												  this%wpmsq, this%w0msq, this%Gm, lbin = lb)

  endif

end subroutine move_window_meta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       update boundaries of the electro-magnetic field
!-----------------------------------------------------------------------------------------
subroutine update_boundary_meta( this, g_space, no_co )

  implicit none

  class( t_emf_meta ), intent( inout )  ::  this

  type( t_space ),     intent(in) :: g_space
  class( t_node_conf ), intent(in) :: no_co

  ! call superclass method
  call this % t_emf % update_boundary( g_space, no_co )

  if ( this%use_tfsf ) then
     call update_boundary_incident( this%bnd_con, this%b_incident, this%e_incident,  g_space, no_co )
  endif

end subroutine update_boundary_meta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function get_min_gc_meta( this )
!-----------------------------------------------------------------------------------------
! returns the number of guard cells required by this
! object. Does not require setup
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_emf_meta ), intent(in) :: this
  integer, dimension(2,p_x_dim) :: get_min_gc_meta

  integer, dimension(2,p_x_dim) :: gc_num_temp
  integer :: i

  ! The "meta" solver uses 1|2 guard cells
  get_min_gc_emf(p_lower,:) = 1
  get_min_gc_emf(p_upper,:) = 2

  ! account for smoothing
  gc_num_temp(p_lower,:) = smooth_order(this%smooth)
  gc_num_temp(p_upper,:) = gc_num_temp(p_lower,:)+1

  do i = 1, p_x_dim
	get_min_gc_meta(p_lower,i)=max(get_min_gc_meta(p_lower,i),gc_num_temp(p_lower,i))
	get_min_gc_meta(p_upper,i)=max(get_min_gc_meta(p_upper,i),gc_num_temp(p_upper,i))
  enddo

end function get_min_gc_meta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine update_hard_point_source( this, grid, no_co, t )
!-----------------------------------------------------------------------------------------
! Update fields to be used by particles
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_emf_meta ), intent( inout )  ::  this
  class( t_grid ), intent(in)      :: grid
  class( t_node_conf ), intent(in) :: no_co
  real( p_double ), intent(in) :: t

  ! local variables
  integer, dimension(3, grid%x_dim) :: nodeIdx
  integer, dimension(p_x_dim)       :: inj_pos_i

  if ( this%use_hard_point_source) then

     nodeIdx = nx_p( grid, no_co, no_co%my_aid() )
     inj_pos_i(1:p_x_dim) = this%hard_point_source%inj_pos_i(1:p_x_dim)

     if ( HasInjectionPoint( inj_pos_i, nodeIdx ) ) then
        call inject_hard_point_source( this%hard_point_source, this%e, this%b, nodeIdx, t )
     endif

  endif

end subroutine update_hard_point_source
!-----------------------------------------------------------------------------------------


end module m_emf_meta
