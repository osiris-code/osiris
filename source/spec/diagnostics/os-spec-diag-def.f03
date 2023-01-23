
! This file is meant to be included through an #include statement

!-----------------------------------------------------------------------------------------
! Species diagnostic classes definitions
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Particle tracks diagnostics and particle track class
!-----------------------------------------------------------------------------------------

type t_track

  integer :: npoints = 0                ! number of points in memory
  integer :: savedpoints = 0            ! number of points saved to disk
  integer, dimension(2) :: tag = 0      ! tag of particle being tracked
  integer :: part_idx = -1              ! index of particle being tracked (-1) means
                                        ! that particle is not on local node

  real(p_k_part), dimension(:,:), pointer :: data => null()  ! track points data
  integer, dimension(:), pointer :: n => null()               ! track points iterations

end type t_track


type t_track_set

  ! maximum number of points to hold
  integer :: maxpoints = -1

  ! number of iterations between track points
  integer :: niter = 1

  ! file holding tags of particles to follow
  character(len = p_max_filename_len) :: file_tags = ''

  ! filename to write tracks to
  character(len = p_max_filename_len) :: file_path = ''
  character(len = p_max_filename_len) :: file_name = ''

  ! total number of tracks
  integer :: ntracks = 0
  ! actual tracks
  type( t_track ), dimension(:), pointer :: tracks  => NULL()

  ! Additional field data
  integer :: nfields = 0

  ! switches to decide which field components to write in tracks
  logical, dimension(p_f_dim)  ::  ifdmp_tracks_efl
  logical, dimension(p_f_dim)  ::  ifdmp_tracks_bfl

  ! Include write psi diagnostic at particle positions
  logical                      ::  ifdmp_tracks_psi

  ! buffers used when storing extra field data in tracks
  real(p_k_part), dimension(:,:), pointer :: field_buffer => null()
  real(p_k_part), dimension(:,:), pointer :: pos_buffer => null()
  integer,        dimension(:,:), pointer :: ipos_buffer => null()

  ! set the initial buffer size to 0, this will force the buffers to be allocated
  ! as soon as the number of particles present on this node is > 0
  integer :: size_buffer = 0

  ! list of tags present / missing
  ! each item is the index of the track on the tracks array
  integer, dimension(:), pointer :: present => null()
  integer, dimension(:), pointer :: missing => null()

  ! list sizes (for simplicity)
  integer :: npresent = 0, nmissing = 0

end type t_track_set


!-----------------------------------------------------------------------------------------
! Phasespace diagnostics and phasespace list Class
!-----------------------------------------------------------------------------------------

integer, parameter, public :: p_max_phasespace_dims = 3
integer, parameter, public :: p_max_n_ene_bins = 32

! Maximum lenght of phasespace name
integer, parameter, public :: p_max_phasespace_namelen = 32


type t_phasespace
  integer           :: ndims = -1
  character(len=p_max_phasespace_namelen) :: name = "-"

  ! quantity to use as the axis for the phasespace
  ! 1 -> position, 2 -> momenta, 3 -> gamma  4 -> log10(gamma)
  integer, dimension( p_max_phasespace_dims ) :: x_or_p = 0

  ! coordinate to use as the axis for the phasespace
  ! has no meaning if x_or_p = 3, 4
  integer, dimension( p_max_phasespace_dims ) :: xp_dim = 0
  type( t_time_avg ) :: tavg
  integer :: ps_type = -1

  type( t_phasespace ), pointer :: next => null()
end type t_phasespace

type t_phasespace_list
  type( t_phasespace ), pointer :: head => null()
  type( t_phasespace ), pointer :: tail => null()

  contains

  procedure :: add_to_list => add_phasespace_to_list
end type

type t_phasespace_diagnostics

  integer :: ndump_fac

  integer :: ndump_fac_tavg

  ! number of time steps to average over for time averaged
  integer :: n_tavg

  ! physical range for phasespace data dumps
  real(p_diag_prec), dimension(p_x_dim) :: xmin, xmax
  real(p_diag_prec), dimension(p_p_dim) :: pmin, pmax
  real(p_diag_prec), dimension(p_p_dim) :: lmin, lmax

  ! switch for autorange of p for phasespaces
  logical, dimension(p_p_dim) :: if_p_auto

  ! switch for autorange of l for phasespaces
  logical, dimension(p_p_dim) :: if_l_auto

  ! physical range for 1D phasespace data dumps
  real(p_diag_prec) :: gammamin, gammamax
  real(p_diag_prec) :: kemin, kemax

  ! switch for autorange of 1D phasespace data dumps
  logical :: if_gamma_auto, if_ke_auto

  ! resolutions for phasespace data dumps
  integer, dimension(p_x_dim) :: nx, nx_3D
  integer, dimension(p_p_dim) :: np, np_3D
  integer, dimension(p_p_dim) :: nl, nl_3D

  ! resolution for 1D phasespace data dumps
  integer :: ngamma, nke

  ! energy binned diagnostics parameters
  integer :: n_ene_bins
  real(p_diag_prec), dimension(p_max_n_ene_bins) :: ene_bins

  ! normal phasespaces
  class( t_phasespace_list ), pointer :: phasespace_list => null()

  ! energy binned phasespaces
  class( t_phasespace_list ), pointer :: pha_ene_bin_list => null()

  ! cell average phasespaces
  class( t_phasespace_list ), pointer :: pha_cell_avg_list => null()

  ! time averaged phasespaces
  class( t_phasespace_list ), pointer :: pha_time_avg_list => null()

  contains

  procedure :: allocate_objs => allocate_objs_phasespace_diag
  procedure :: init_list => init_phasespace_list
  procedure :: get_params => get_phasespace_parameters
  procedure :: size => size_phasespace

end type

type t_phasespace_params
  integer :: ndims

  integer, dimension(3) :: n
  real(p_k_part), dimension(3) :: min, max
  logical, dimension(3) :: move_c

  character( len = 80 ) :: name, long_name, unit
  character( len = 80 ), dimension(3) :: xname, xlabel, xunits
  real(p_double), dimension(3) :: offset_x = p_no_offset
  real(p_double) :: offset_t = p_no_offset
  ! This offset is for phasespace axes
  real(p_double), dimension(3) :: offset_t_ax = p_no_offset
end type t_phasespace_params

!-----------------------------------------------------------------------------------------
! Species Diagnostics Class
!-----------------------------------------------------------------------------------------
integer, parameter :: p_max_diag_len = 18

type :: t_diag_species

  
  ! Normal diagnostics
  character( len = p_max_diag_len ), dimension(:), pointer :: report_quants => null()

  ! density reports
  type(t_vdf_report), pointer :: reports => null()

  ! cell average reports
  type(t_vdf_report), pointer :: rep_cell_avg => null()

  ! udist reports
  type(t_vdf_report), pointer :: rep_udist => null()

  ! frequency of species dianostic data dumps
  integer :: ndump_fac_ene         ! energy
  integer :: ndump_fac_heatflux    ! Heat flux
  integer :: ndump_fac_temp        ! Temperature
  integer :: ndump_fac_raw         ! Particle data dumps

  ! parameters for raw data dump
  real(p_diag_prec) :: raw_gamma_limit
  real(p_diag_prec) :: raw_fraction
  ! if raw_if_pos_ref_box is true in one direction, the raw diagnostic position of this
  ! direction will refer the the simulation box edge instead of global coordinates
  logical, dimension(p_x_dim) :: raw_if_pos_ref_box
  ! parameters for function parser
  character(len = p_max_expr_len) :: raw_math_expr = " "
  type(t_fparser) :: raw_func

  ! particle tracking data
  integer :: ndump_fac_tracks      ! frequency of track diagnostic writes
  integer :: n_start_tracks = -1   ! iteration to start writing tracks
  type( t_track_set ) :: tracks

  ! phasespaces
  class( t_phasespace_diagnostics ), pointer :: phasespaces => null()

  contains

  procedure :: allocate_objs       => allocate_objs_diag_species
  procedure :: avail_report_quants => avail_report_quants_species
  procedure :: init_report_quants  => init_report_quants_species
  procedure :: init                => init_diag_species
  procedure :: read_input          => read_input_diag_species
  procedure :: raw_ref_box         => raw_ref_box_species

end type t_diag_species

interface
subroutine allocate_objs_diag_species( this )
  import t_diag_species
  class( t_diag_species ), intent( inout ) :: this
end subroutine
end interface

interface
subroutine init_diag_species( this, spec_name, n_x_dim, ndump_fac, interpolation, &
                              restart, restart_handle )
  import t_diag_species, t_restart_handle
  class ( t_diag_species ), intent( inout )  ::  this
  character( len=* ), intent(in) :: spec_name
  integer, intent(in) :: n_x_dim
  integer, intent(in) :: ndump_fac
  integer, intent(in) :: interpolation
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
end subroutine
end interface

interface
subroutine read_input_diag_species( this, input_file, gamma )
  import t_diag_species, t_input_file, p_double
  class ( t_diag_species ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  real(p_double), intent(in) :: gamma
end subroutine
end interface

interface
function avail_report_quants_species( this )
  import t_diag_species
  class( t_diag_species ), intent(in) :: this
  integer :: avail_report_quants_species
end function
end interface

interface
subroutine init_report_quants_species( this )
  import t_diag_species
  class( t_diag_species ), intent(inout) :: this
end subroutine
end interface

interface
pure function raw_ref_box_species( this, dim )
  import t_diag_species
  class( t_diag_species ), intent(in) :: this
  integer, intent(in) :: dim
  logical :: raw_ref_box_species
end function
end interface

!------- phasespace routines --------
interface
  subroutine add_phasespace_to_list( list, ndims, x_or_p, xp_dim, ps_type )
    import t_phasespace_list
    class( t_phasespace_list ), intent( inout ) :: list
    integer, intent(in) :: ndims
    integer, dimension(:), intent(in) :: x_or_p, xp_dim
    integer, intent(in) :: ps_type
  end subroutine
end interface

interface
  subroutine allocate_objs_phasespace_diag( this )
    import t_phasespace_diagnostics
    class( t_phasespace_diagnostics ), intent(inout) :: this
  end subroutine
end interface

interface
  subroutine init_phasespace_list( this, list, phasespaces, msg, time_average )
    import t_phasespace_diagnostics, t_phasespace_list
    class( t_phasespace_diagnostics ), intent(in) :: this
    class( t_phasespace_list ), intent(inout) :: list
    character(len=*), dimension(:), intent(in) :: phasespaces
    character(len=*), intent(in) :: msg
    logical, intent(in), optional :: time_average
  end subroutine
end interface

interface
  subroutine get_phasespace_parameters( this, phasespace, g_space, params )
    import t_phasespace_diagnostics, t_phasespace, t_space, t_phasespace_params
    class( t_phasespace_diagnostics ), intent(in) :: this
    type( t_phasespace ), intent(in) :: phasespace
    type( t_space ), intent(in) :: g_space
    type (t_phasespace_params), intent(out) :: params
  end subroutine
end interface

interface
  function size_phasespace( this, phasespace )
    import t_phasespace_diagnostics, t_phasespace
    class( t_phasespace_diagnostics ),  intent(in) :: this
    type( t_phasespace ),    intent(in) :: phasespace
    integer :: size_phasespace
  end function
end interface
