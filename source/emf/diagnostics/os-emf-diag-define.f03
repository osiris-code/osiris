#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_diag

  use stringutil
  use m_diagnostic_utilities
  use m_vdf_define
  use m_parameters

  implicit none

  private

  integer, public :: diag_emf_ev = 0

!-------------------------------------------------------------------------------
! Diagnostics class definition
!-------------------------------------------------------------------------------

character(len=10), dimension(34), parameter :: p_report_quants = &
   (/ 'e1        ', 'e2        ', 'e3        ', 'b1        ', 'b2        ', 'b3        ', &
      'ext_e1    ', 'ext_e2    ', 'ext_e3    ', 'ext_b1    ', 'ext_b2    ', 'ext_b3    ', &
      'part_e1   ', 'part_e2   ', 'part_e3   ', 'part_b1   ', 'part_b2   ', 'part_b3   ', &
      'ene_e1    ', 'ene_e2    ', 'ene_e3    ', 'ene_b1    ', 'ene_b2    ', 'ene_b3    ', &
      'ene_e     ', 'ene_b     ', 'ene_emf   ', &
      'div_e     ', 'div_b     ', 'psi       ', 'chargecons', &
      's1        ', 's2        ', 's3        '/)

integer, parameter, public :: p_e1 = 1, p_e2 = 2, p_e3 = 3
integer, parameter, public :: p_b1 = 4, p_b2 = 5, p_b3 = 6
integer, parameter, public :: p_ext_e1 = 7, p_ext_e2 = 8, p_ext_e3 = 9
integer, parameter, public :: p_ext_b1 = 10, p_ext_b2 = 11, p_ext_b3 = 12
integer, parameter, public :: p_part_e1 = 13, p_part_e2 = 14, p_part_e3 = 15
integer, parameter, public :: p_part_b1 = 16, p_part_b2 = 17, p_part_b3 = 18
integer, parameter, public :: p_ene_e1 = 19, p_ene_e2 = 20, p_ene_e3 = 21
integer, parameter, public :: p_ene_b1 = 22, p_ene_b2 = 23, p_ene_b3 = 24
integer, parameter, public :: p_ene_e = 25, p_ene_b = 26, p_ene_emf = 27
integer, parameter, public :: p_div_e = 28, p_div_b = 29
integer, parameter, public :: p_psi = 30, p_charge_cons = 31, p_s1 = 32, p_s2 = 33, p_s3 = 34


integer, parameter :: p_max_diag_len = 14


!-----------------------------------------------------------------------------------------
type :: t_diag_emf

  ! Normal diagnostics
  character( len = p_max_diag_len ), dimension(:), pointer :: report_quants => null()

  type(t_vdf_report), pointer :: reports

  ! Total integrated field energy
  integer :: ndump_fac_ene_int

  ! Charge conservation test
  integer :: ndump_fac_charge_cons

  ! precision of diagnostics
  integer :: prec = p_diag_prec

contains

  procedure :: restart_write       => restart_write_diag_emf
  procedure :: restart_read        => restart_read_diag_emf  
  procedure :: allocate_objs       => allocate_objs_diag_emf
  procedure :: avail_report_quants => avail_report_quants_emf
  procedure :: init_report_quants  => init_report_quants_emf
  procedure :: init                => init_diag_emf
  procedure :: read_input          => read_input_diag_emf
  procedure :: cleanup             => cleanup_diag_emf
  procedure :: quant_offset        => quant_offset_diag_emf

end type t_diag_emf
!-----------------------------------------------------------------------------------------

  public :: t_diag_emf

 contains

!-----------------------------------------------------------------------------------------
function quant_offset_diag_emf( this )

  class( t_diag_emf ), intent(in) :: this
  integer :: quant_offset_diag_emf

  ! This is only relevant for subclasses extending the number of available
  ! diagnostic quantities
  quant_offset_diag_emf = 0

end function quant_offset_diag_emf
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine restart_write_diag_emf( this, restart_handle )

  use m_restart

  implicit none

  class ( t_diag_emf ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  ! Do nothing for now (boosted diag to go here)

end subroutine restart_write_diag_emf
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine restart_read_diag_emf( this, restart_handle )

  use m_restart

  implicit none

  class ( t_diag_emf ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  ! Do nothing for now (boosted diag to go here)

end subroutine restart_read_diag_emf
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! Initialize individual diagnostics
!-----------------------------------------------------------------------------------------
subroutine init_diag_emf( this, no_ext_fld, part_fld_alloc, interpolation )

  use m_logprof

  implicit none

  class( t_diag_emf ), intent(inout) :: this
  logical, intent(in) :: no_ext_fld
  logical, intent(in) :: part_fld_alloc
  integer, intent(in) :: interpolation

  type(t_vdf_report), pointer :: report
  integer, parameter :: izero = ichar('0')
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

    ! Set default units for all diagnostics
    report%xname  = (/'x1', 'x2', 'x3'/)
    report%xlabel = (/'x_1', 'x_2', 'x_3'/)
    report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

    ! these are just dummy values for now
    report%time_units = '1 / \omega_p'
    report%dt         = 1.0

    ! all field quantities are aligned in time
    report%offset_t = 0.0_p_double

    report%fileLabel = ''
    report%basePath  = trim(path_mass) // 'FLD' // p_dir_sep

    ! Set specific diagnostic parameters
    select case ( report%quant )
      case ( p_e1, p_e2, p_e3 )
        report%label = 'E_'//char(izero + 1 + report%quant - p_e1)
        report%units = 'm_e c \omega_p e^{-1}'
        report%offset_x = (/ 0.0_p_double, 0.0_p_double, 0.0_p_double /) + interp_offset
        report%offset_x(report%quant-p_e1+1) = report%offset_x(report%quant-p_e1+1) + 0.5_p_double
      case ( p_b1, p_b2, p_b3 )
        report%label = 'B_'//char(izero + 1 + report%quant - p_b1)
        report%units = 'm_e c \omega_p e^{-1}'
        report%offset_x = (/ 0.5_p_double, 0.5_p_double, 0.5_p_double /) + interp_offset
        report%offset_x(report%quant-p_b1+1) = report%offset_x(report%quant-p_b1+1) - 0.5_p_double
      case ( p_ext_e1, p_ext_e2, p_ext_e3 )
        if ( no_ext_fld ) then
SCR_ROOT('(*warning*) External field diagnostic requested, but external fields not in use')
          report%ndump = 0
        else
          report%label = 'E_'//char(izero + 1 + report%quant - p_ext_e1)//'^{ext}'
          report%units = 'm_e c \omega_p e^{-1}'
          report%offset_x = (/ 0.0_p_double, 0.0_p_double, 0.0_p_double /) + interp_offset
          report%offset_x(report%quant-p_ext_e1+1) = report%offset_x(report%quant-p_ext_e1+1) + 0.5_p_double
        endif

      case ( p_ext_b1, p_ext_b2, p_ext_b3 )
        if ( no_ext_fld ) then
          SCR_ROOT('(*warning*) External field diagnostic requested, but external fields')
          SCR_ROOT('            not in use.')
          report%ndump = 0
        else
          report%label = 'B_'//char(izero + 1 + report%quant - p_ext_b1)//'^{ext}'
          report%units = 'm_e c \omega_p e^{-1}'
          report%offset_x = (/ 0.5_p_double, 0.5_p_double, 0.5_p_double /) + interp_offset
          report%offset_x(report%quant-p_ext_b1+1) = report%offset_x(report%quant-p_ext_b1+1) - 0.5_p_double
        endif

      case ( p_part_e1, p_part_e2, p_part_e3 )
        if ( .not. part_fld_alloc ) then
          SCR_ROOT('(*warning*) Particle fields diagnostic requested, but external/smoothed')
          SCR_ROOT('            fields are not in use. Use main field diagnostics instead.')
          report%ndump = 0
        else
          report%label = 'E_'//char(izero + 1 + report%quant - p_part_e1)//'^{part}'
          report%units = 'm_e c \omega_p e^{-1}'
          report%offset_x = (/ 0.0_p_double, 0.0_p_double, 0.0_p_double /) + interp_offset
          report%offset_x(report%quant-p_part_e1+1) = report%offset_x(report%quant-p_part_e1+1) + 0.5_p_double
        endif

      case ( p_part_b1, p_part_b2, p_part_b3 )
        if ( .not. part_fld_alloc ) then
          SCR_ROOT('(*warning*) Particle fields diagnostic requested, but external/smoothed')
          SCR_ROOT('            fields are not in use.')
          report%ndump = 0
        else
          report%label = 'B_'//char(izero + 1 + report%quant - p_part_b1)//'^{part}'
          report%units = 'm_e c \omega_p e^{-1}'
          report%offset_x = (/ 0.5_p_double, 0.5_p_double, 0.5_p_double /) + interp_offset
          report%offset_x(report%quant-p_part_b1+1) = report%offset_x(report%quant-p_part_b1+1) - 0.5_p_double
        endif

      case ( p_ene_e1, p_ene_e2, p_ene_e3 )
        report%label = 'E_'//char(izero + 1 + report%quant - p_ene_e1)//'^2'
        report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}'
        report%offset_x = (/ 0.0_p_double, 0.0_p_double, 0.0_p_double /) + interp_offset
        report%offset_x(report%quant-p_ene_e1+1) = report%offset_x(report%quant-p_ene_e1+1) + 0.5_p_double

      case ( p_ene_b1, p_ene_b2, p_ene_b3 )
        report%label = 'B_'//char(izero + 1 + report%quant - p_ene_b1)//'^2'
        report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}'
        report%offset_x = (/ 0.5_p_double, 0.5_p_double, 0.5_p_double /) + interp_offset
        report%offset_x(report%quant-p_ene_b1+1) = report%offset_x(report%quant-p_ene_b1+1) - 0.5_p_double

      case ( p_ene_e )
        report%label = 'E^2'
        report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}'

      case ( p_ene_b )
        report%label = 'B^2'
        report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}'

      case ( p_ene_emf )
        report%label = 'E^2 + B^2'
        report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}'

      case ( p_div_e )
        report%label = '\bf{\nabla}\cdot\bf{E}'
        report%units = 'm_e c \omega_p e^{-1}'
        report%offset_x = (/ 0.0_p_double, 0.0_p_double, 0.0_p_double /) + interp_offset

      case ( p_div_b )
        report%label = '\bf{\nabla}\cdot\bf{B}'
        report%units = 'm_e c \omega_p e^{-1}'
        report%offset_x = (/ 0.5_p_double, 0.5_p_double, 0.5_p_double /) + interp_offset

      case ( p_charge_cons )

        report%label     = '\bf{\nabla}\cdot\bf{E} - \rho'
        report%units     = 'm_e c \omega_p e^{-1}'
        report%offset_x = (/ 0.0_p_double, 0.0_p_double, 0.0_p_double /) + interp_offset

        ! this diagnostic overrides some default parameters
        report%basePath  = trim(path_mass) // 'CHARGECONS' // p_dir_sep
        report%prec      = 8
        report%ndump     = this%ndump_fac_charge_cons

      case ( p_psi )
        report%label     = '\Psi_x'
        report%units     = 'a.u.'
        ! aligned with e1
        report%offset_x = (/ 0.5_p_double, 0.0_p_double, 0.0_p_double /) + interp_offset

      case ( p_s1, p_s2, p_s3 )
        report%label = 'S_'//char(izero + 1 + report%quant - p_s1)
        report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}'
        report%offset_x = (/ 0.0_p_double, 0.0_p_double, 0.0_p_double /) + interp_offset

      case default
        ! unknown diagnostic belonging to subclass, will be handled later
        ! print *, 'In setup_diag_emf, unknown report = ', report % name
        continue

    end select

    ! process next report
    report => report % next
  enddo

  if (diag_emf_ev==0) then
    diag_emf_ev = create_event('EMF diagnostics')
  endif

end subroutine init_diag_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_input_diag_emf( this, input_file, gamma )
!-----------------------------------------------------------------------------------------
!       read necessary information from inputdec
!-----------------------------------------------------------------------------------------

  use m_input_file
  use m_vdf_report

  implicit none

  class( t_diag_emf ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  real(p_double), intent(in) :: gamma

  integer                     :: ndump_fac
  integer                     :: ndump_fac_ave
  integer                     :: ndump_fac_lineout

  integer                     :: ndump_fac_ene_int
  integer                     :: ndump_fac_charge_cons

  integer                     :: prec

  integer                     :: n_tavg
  integer, dimension(p_x_dim) :: n_ave


  character( len = p_max_reports_len ), dimension( p_max_reports ) :: reports

  namelist /nl_diag_emf/ ndump_fac, ndump_fac_ave, ndump_fac_lineout, &
                         ndump_fac_ene_int, ndump_fac_charge_cons, &
                         prec, n_tavg, n_ave, reports

  integer, dimension( p_n_report_type ) :: ndump_fac_all
  integer :: ierr

  ndump_fac         = 0
  ndump_fac_ave     = 0
  ndump_fac_lineout = 0

  ndump_fac_ene_int =  0
  ndump_fac_charge_cons = 0

  prec    =  p_single
  n_ave         = -1
  n_tavg        = -1

  reports =  "-"

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_diag_emf", ierr )

  if (ierr == 0) then
    read (input_file%nml_text, nml = nl_diag_emf, iostat = ierr)
    if (ierr /= 0) then
      print *, "Error reading diag_emf parameters"
      print *, "aborting..."
      stop
    endif
  else
    if (ierr < 0) then
      print *, "Error reading diag_emf parameters"
      print *, "aborting..."
      stop
    else
      if (disp_out(input_file)) then
        SCR_ROOT(" - no diagnostics specified")
      endif
    endif
  endif

  ndump_fac_all(p_full)  = ndump_fac
  ndump_fac_all(p_savg)  = ndump_fac_ave
  ndump_fac_all(p_senv)  = ndump_fac_ave
  ndump_fac_all(p_line)  = ndump_fac_lineout
  ndump_fac_all(p_slice) = ndump_fac_lineout

  this%ndump_fac_ene_int     =  ndump_fac_ene_int
  this%ndump_fac_charge_cons = ndump_fac_charge_cons

  ! get list of available diagnostics
  call this % init_report_quants( )

  ! Initialize reports
  call new( this % reports, reports, this  % report_quants, &
            ndump_fac_all, n_ave, n_tavg, prec, &
            p_x_dim, ierr )
  if ( ierr /= 0 ) then
     print *, "(*error*) Invalid report"
     print *, "(*error*) aborting..."
     stop
  endif

end subroutine read_input_diag_emf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine allocate_objs_diag_emf( this ) 

  implicit none

  class( t_diag_emf ), intent( inout )  ::  this

  ! Do nothing for now (boosted diag to go here)

end subroutine allocate_objs_diag_emf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine cleanup_diag_emf(this)

  use m_vdf_report

  implicit none

  class( t_diag_emf ), intent( inout )  ::  this

  ! cleanup reports
  call cleanup( this%reports )

  deallocate( this % report_quants )

end subroutine cleanup_diag_emf
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
function avail_report_quants_emf( this )

  class( t_diag_emf ), intent(in) :: this

  integer :: avail_report_quants_emf

  avail_report_quants_emf = size( p_report_quants )

end function avail_report_quants_emf


subroutine init_report_quants_emf( this )

  class( t_diag_emf ), intent(inout) :: this

  ! print *, '[diag-emf] In init_report_quants_emf'

  if ( .not. associated( this % report_quants ) ) then
    allocate( this % report_quants( this % avail_report_quants( ) ) )
  endif

  ! A subclass may have a larger set of diagnostics
  this % report_quants( 1 : size( p_report_quants ) ) = p_report_quants

end subroutine init_report_quants_emf


end module m_emf_diag
