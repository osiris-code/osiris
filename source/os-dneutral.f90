!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     neutral diagnostics class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

#ifdef __HAS_IONIZATION__

module m_diag_neutral

  use m_system
  use m_parameters

  use m_vdf_define
  use m_vdf_report
  use m_vdf_math
  use m_vdf

  use m_node_conf

  use m_grid_define

  use m_space

  use m_diagnostic_utilities
  use stringutil

  implicit none

  private

  character(len=*), dimension(2), parameter :: p_report_quants = &
     (/ 'ion_charge', 'neut_den  ' /)

  integer, parameter :: p_ion_charge = 1, &
                        p_neut_den = 2, &
                        p_ion_rate = 3

  type :: t_diag_neutral

    !  diagnostics
    type(t_vdf_report), pointer :: reports

  end type t_diag_neutral

  interface read_nml
    module procedure read_nml_diag_neutral
  end interface

  interface report
    module procedure report_diag_neutral
  end interface

  interface cleanup
    module procedure cleanup_diag_neutral
  end interface

  interface setup
    module procedure setup_diag_neutral
  end interface

!       declare things that should be public
  public :: t_diag_neutral, setup, cleanup
  public :: report, read_nml
  public :: p_ion_charge, p_neut_den

 contains

!---------------------------------------------------
subroutine read_nml_diag_neutral( this, input_file)

  use m_input_file

  implicit none

  type( t_diag_neutral ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  integer                     :: ndump_fac
  integer                     :: ndump_fac_ave
  integer                     :: ndump_fac_lineout
  integer                   :: prec
  integer, dimension(p_x_dim) :: n_ave
  integer :: n_tavg

  character( len = p_max_reports_len ), dimension( p_max_reports ) :: reports

  namelist /nl_diag_neutral/ ndump_fac, ndump_fac_ave, ndump_fac_lineout, &
                             prec, n_ave, n_tavg, reports

  integer, dimension( p_n_report_type ) :: ndump_fac_all
  integer :: ierr

  ! executable statements

  ndump_fac         = 0
  ndump_fac_ave     = 0
  ndump_fac_lineout = 0

  n_ave         =  1
  n_tavg = -1
  prec = p_single
  reports =  "-"

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_diag_neutral", ierr )

  if (ierr == 0) then
    read (input_file%nml_text, nml = nl_diag_neutral, iostat = ierr)
    if (ierr /= 0) then
      print *, "Error reading diag_neutral parameters"
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

  ! process reports
  call new( this%reports, reports, p_report_quants, &
            ndump_fac_all, n_ave, n_tavg, prec, &
            p_x_dim, ierr )
  if ( ierr /= 0 ) then
     print *, "(*error*) Invalid report"
	 print *, "(*error*) aborting..."
	 stop
  endif

end subroutine read_nml_diag_neutral
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine setup_diag_neutral( this, neut_name, interpolation )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  type( t_diag_neutral ), intent(inout) :: this
  character(len=*), intent(in) :: neut_name
  integer, intent(in) :: interpolation

  ! local variables
  type(t_vdf_report), pointer :: report
  real(p_double) :: interp_offset

  ! Process reports
  report => this%reports

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

  do
    if ( .not. associated( report ) ) exit

    report%xname  = (/'x1', 'x2', 'x3'/)

    report%xlabel = (/'x_1', 'x_2', 'x_3'/)
    report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

    ! these are just dummy values for now
    report%time_units = '1 / \omega_p'
    report%dt         = 1.0

    report%fileLabel = replace_blanks(trim(neut_name))
    report%basePath  = trim(path_mass) // 'ION' // p_dir_sep // &
                       trim(report%fileLabel) // p_dir_sep

    select case ( report%quant )
      case ( p_ion_charge )
        report%label = 'Ion Charge'
        report%units = 'n_0'
        report%offset_t = 0.0_p_double

      case ( p_neut_den )
        report%label = 'Neutral Density'
        report%units = 'n_0'
        report%offset_t = 0.0_p_double

      case default
        print *, '(*error*) Invalid quantity for diagnostic, ', trim(report%name)
        call abort_program(-1)

    end select

    report%offset_x = (/ 0.0_p_double, 0.0_p_double, 0.0_p_double /) + interp_offset

    ! process next report
    report => report % next
  enddo

end subroutine setup_diag_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine report_diag_neutral( this, multi_ion, ion_idx, neut_idx, g_space, &
                                grid, no_co, tstep, t )
!-----------------------------------------------------------------------------------------

  use m_time_step

  implicit none

  ! dummy variables

  type( t_diag_neutral ), intent(inout) :: this

  type(t_vdf),            intent(in) :: multi_ion

  integer,                intent(in) :: ion_idx, neut_idx
  type( t_space ),        intent(in) :: g_space
  class( t_grid ),         intent(in) :: grid
  class( t_node_conf ),    intent(in) :: no_co
  type( t_time_step ),    intent(in) :: tstep
  real(p_double),         intent(in) :: t

  ! local variables
  type( t_vdf_report ), pointer :: rep
  type( t_vdf ) :: neut_den, neut_default

  rep => this%reports
  do
    if ( .not. associated( rep ) ) exit

    if ( if_report( rep, tstep ) ) then
       select case ( rep%quant )
         case ( p_ion_charge )
           call report_vdf( rep, multi_ion, ion_idx, g_space, grid, no_co, tstep, t )

         case ( p_neut_den )

           ! get real neutral density
           call neut_den % new( multi_ion, copy = .true., f_dim = 1, fc = 1 )
           call neut_default % new( multi_ion, copy = .true., f_dim = 1, fc = neut_idx )

           call mult(neut_den, neut_default)

           ! save result
           call report_vdf( rep, neut_den, 1, g_space, grid, no_co, tstep, t )

           ! free used memory
           call neut_default % cleanup()
           call neut_den % cleanup( )

       end select
    endif

    rep => rep%next
  enddo

end subroutine report_diag_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine cleanup_diag_neutral( this )
!-----------------------------------------------------------------------------------------

  implicit none

  type( t_diag_neutral ),       intent(inout) :: this

  call cleanup( this%reports )

end subroutine cleanup_diag_neutral
!-----------------------------------------------------------------------------------------

end module m_diag_neutral

#endif
