!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Electric current diagnostics class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module m_current_diag

#include "os-config.h"
#include "os-preprocess.fpp"

  use m_system
  use m_parameters

  use m_space
  use m_grid_define
  use m_node_conf

  use m_logprof
  use m_utilities
  use m_diagnostic_utilities
  use stringutil

  use m_vdf_define
  use m_vdf_report
  use m_vdf_average
  use m_vdf_math
  use m_vdf

  implicit none

  private

  integer, public :: diag_current_ev = 0

  character(len=10), dimension(4), parameter :: p_report_quants = &
   (/ 'j1   ','j2   ','j3   ','div_j'/)

  integer, parameter, public :: p_j1 = 1, p_j2 = 2, p_j3 = 3, p_div_j = 4

!-----------------------------------------------------------------------------------------
! t_current_diag class
!   Diagnostics for t_current objects
!-----------------------------------------------------------------------------------------

type :: t_current_diag

  ! Normal diagnostics
  character( len = p_max_diag_len ), dimension(:), pointer :: report_quants => null()

  type(t_vdf_report), pointer :: reports => null()

  contains

  procedure :: avail_report_quants  => avail_report_quants_current
  procedure :: init_report_quants   => init_report_quants_current
  procedure :: init                 => init_diag_current
  procedure :: read_input           => read_input_diag_current
  procedure :: get_diag_buffer_size => get_diag_buffer_size
  procedure :: cleanup              => cleanup_diag_current

end type t_current_diag

public :: t_current_diag

contains

!-----------------------------------------------------------------------------------------
! Return available quantities for diagnostics
!-----------------------------------------------------------------------------------------
function avail_report_quants_current( this )

  implicit none

  class( t_current_diag ), intent(in) :: this

  integer :: avail_report_quants_current

  avail_report_quants_current = size( p_report_quants )

end function avail_report_quants_current

!-----------------------------------------------------------------------------------------
! Initialize available report quantities
!-----------------------------------------------------------------------------------------
subroutine init_report_quants_current( this )

  class( t_current_diag ), intent(inout) :: this

  ! print *, '[diag-emf] In init_report_quants_current'

  if ( .not. associated( this % report_quants ) ) then
    allocate( this % report_quants( this % avail_report_quants( ) ) )
  endif

  ! A subclass may have a larger set of diagnostics
  this % report_quants( 1 : size( p_report_quants ) ) = p_report_quants

end subroutine init_report_quants_current

!-----------------------------------------------------------------------------------------
! Initialize object
!-----------------------------------------------------------------------------------------
subroutine init_diag_current( this, interpolation )

  use m_system
  use m_parameters

  use m_space
  use m_grid_define
  use m_node_conf

  use m_logprof
  use m_utilities
  use m_diagnostic_utilities
  use stringutil

  use m_vdf_define
  use m_vdf_report
  use m_vdf_average
  use m_vdf_math
  use m_vdf

  implicit none

  class( t_current_diag ), intent(inout) :: this
  integer, intent(in) :: interpolation

  type( t_vdf_report ), pointer :: report
  integer, parameter :: izero = iachar('0')
  real(p_double) :: interp_offset

  select case( interpolation )
  case( p_linear, p_cubic )
    interp_offset = 0.0_p_double
  case( p_quadratic, p_quartic )
    interp_offset = 0.5_p_double
  case default
    interp_offset = 0.0_p_double
    ERROR('Interpolation value not supported')
    call abort_program(p_err_invalid)
  end select

  ! Normal reports
  report => this%reports
  do
    if ( .not. associated( report ) ) exit

    report%xname  = (/'x1', 'x2', 'x3'/)
    report%xlabel = (/'x_1', 'x_2', 'x_3'/)
    report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

    ! these are just dummy values for now
    report%time_units = '1 / \omega_p'
    report%dt         = 1.0

    ! current is 1/2 time step behind
    report%offset_t = -0.5_p_double

    report%fileLabel = ''
    report%basePath  = trim(path_mass) // 'FLD' // p_dir_sep

    select case ( report%quant )
    case ( p_j1, p_j2, p_j3 )
      report%label = 'j_'//char(izero + 1 + report%quant - p_j1)
      select case ( p_x_dim )
      case (1)
        report%units = 'e \omega_p'
      case (2)
        report%units = 'e \omega_p^2 / c'
      case (3)
        report%units = 'e \omega_p^3 / c^2'
      end select
      report%offset_x = (/ 0.0_p_double, 0.0_p_double, 0.0_p_double /) + interp_offset
      report%offset_x(report%quant-p_j1+1) = report%offset_x(report%quant-p_j1+1) + 0.5_p_double

    case ( p_div_j )
      report%label = '\bf{\nabla}\cdot\bf{j}'
      if ( p_x_dim == 1 ) then
        report%units  = 'e \omega_p / c'
      else
        report%units  = 'e \omega_p^'//char(izero+p_x_dim)// &
                        '/ c^'//char(izero+p_x_dim)
      endif
      report%offset_x = (/ 0.0_p_double, 0.0_p_double, 0.0_p_double /) + interp_offset

    end select

    report => report%next
  enddo

  if (diag_current_ev==0) diag_current_ev = create_event('electric current diagnostics')

end subroutine init_diag_current
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Read input
!-----------------------------------------------------------------------------------------
subroutine read_input_diag_current( this, input_file )

  use m_system
  use m_parameters
  use m_input_file

  use m_space
  use m_grid_define
  use m_node_conf

  use m_logprof
  use m_utilities
  use m_diagnostic_utilities
  use stringutil

  use m_vdf_define
  use m_vdf_report
  use m_vdf_average
  use m_vdf_math
  use m_vdf

  implicit none

  class( t_current_diag ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  integer                     :: ndump_fac
  integer                     :: ndump_fac_ave
  integer                     :: ndump_fac_lineout

  integer                     :: prec
  integer, dimension(p_x_dim) :: n_ave
  integer                     :: n_tavg

  character( len = p_max_reports_len ), dimension( p_max_reports ) :: reports

  namelist /nl_diag_current/ ndump_fac, ndump_fac_ave, ndump_fac_lineout, &
                             prec, n_ave, n_tavg, reports

  integer, dimension( p_n_report_type ) :: ndump_fac_all
  integer :: ierr

  ! executable statements

  ndump_fac   =  0
  ndump_fac_ave   =  0
  ndump_fac_lineout = 0

  prec           =  p_single
  n_ave          = -1
  n_tavg         = -1

  reports      =  "-"

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_diag_current", ierr )

  if (ierr == 0) then
  read (input_file%nml_text, nml = nl_diag_current, iostat = ierr)
  if (ierr /= 0) then
    print *, "Error reading diag_current parameters"
    print *, "aborting..."
    stop
  endif
  else
    if (disp_out(input_file)) then
      SCR_ROOT(" - no diagnostics specified")
    endif
  endif

  ndump_fac_all(p_full)  = ndump_fac
  ndump_fac_all(p_savg)  = ndump_fac_ave
  ndump_fac_all(p_senv)  = ndump_fac_ave
  ndump_fac_all(p_line)  = ndump_fac_lineout
  ndump_fac_all(p_slice) = ndump_fac_lineout

  call this % init_report_quants()

  ! process normal reports
  call new( this%reports, reports, this % report_quants, &
            ndump_fac_all, n_ave, n_tavg, prec, &
            p_x_dim, ierr )

  if ( ierr /= 0 ) then
     print *, "(*error*) Invalid report"
   print *, "(*error*) aborting..."
   stop
  endif


end subroutine read_input_diag_current
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Returns the required diagnostic buffer size
!-----------------------------------------------------------------------------------------
subroutine get_diag_buffer_size( this, gnx, diag_buffer_size )

  implicit none

  class( t_current_diag ), intent(in) :: this
  integer, dimension(:), intent(in) :: gnx
  integer, intent(inout) :: diag_buffer_size

  integer, dimension(3) :: ln_avg

  integer :: i, bsize

  ln_avg = n_avg( this%reports, p_x_dim )

  if (ln_avg(1) > 0) then
   bsize = gnx(1) / ln_avg(1)
   do i = 2, p_x_dim
     bsize = bsize * ( gnx(i) / ln_avg(i) )
   enddo

   ! we need 2 buffers (why?)
   bsize = 2*bsize
   if ( bsize > diag_buffer_size ) diag_buffer_size = bsize
  endif

end subroutine get_diag_buffer_size


!-----------------------------------------------------------------------------------------
! Cleanup object
!-----------------------------------------------------------------------------------------
subroutine cleanup_diag_current( this )

  implicit none

  class( t_current_diag ),   intent(inout) ::  this

  call cleanup( this%reports )

end subroutine cleanup_diag_current

end module m_current_diag
