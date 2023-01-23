!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     cathode class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module m_cathode

#include "os-config.h"
#include "os-preprocess.fpp"
#include "memory/memory.h"

  use m_parameters
  use m_species_define, only : t_species
  use m_fparser, only : t_fparser

  implicit none

  private

  integer, parameter :: p_gaussian  = 2  ! gaussian
  integer, parameter :: p_channel   = 3  ! parabolic channel

  type :: t_cathode

    ! allow access to type components only to module procedures
    !private

    ! cathode number - internal id
    integer :: cathode_id

    ! direction to inject particles
    integer :: dir

    ! wall to inject particles from
    integer :: wall

    ! use moving injector plane
    logical :: mov_inj

    ! position of the next particle to be injected
    real(p_k_part) :: nppos
    integer :: inppos ! this is the global cell index of the particle

    ! time profile of the injected electrons
    real(p_k_part) :: t_start, t_rise, t_fall, t_flat ! input parameters
    real(p_k_part) :: t_total                ! total length of the bunch

    ! deposit current of injected particles
    logical :: deposit_current

    ! minimum density to inject particles
    real(p_double) :: den_min

    ! transverse density profile of the injected beam
    integer :: prof_type                  ! density profile type
    real(p_double) :: density            ! total(peak) density of the beam

    real(p_double), dimension(2) :: center

    real(p_double) :: gauss_width              ! total gauss_width to be considered
    real(p_double) :: gauss_w0                 ! gaussian waist to be considered

    ! channel parameters
    real(p_double)               :: channel_bottom
    real(p_double)               :: channel_r0
    real(p_double)               :: channel_depth
    real(p_double)               :: channel_size
    real(p_double)               :: channel_wall

    ! species to receive produced particles
    class( t_species ), pointer :: species => null()
    real(p_double) :: rgamma
    real(p_double), dimension(p_p_dim)  :: vel   ! fluid velocity (not momenta) of the injected
                                                 ! particles normalized to the cell size
    integer            :: sp_id

    ! positions of particles inside the cell
    ! (this avoids calculating/allocating it at every time)
    real(p_k_part), dimension(:,:), pointer :: ppos_cell => null()

  contains

    procedure :: read_input          => read_input_cathode
    procedure :: inject              => inject_cathode
    procedure :: inject_lowroundoff  => inject_cathode_lowroundoff
    procedure :: init                => init_cathode
    procedure :: cleanup             => cleanup_cathode
    procedure :: move_window         => move_window_cathode
    procedure :: write_checkpoint    => write_checkpoint_cathode

  end type t_cathode

  ! string to id restart data
  character(len=*), parameter :: p_cath_rst_id = "cathode rst data - 0x0003"

  integer, parameter, dimension(3) :: perp1_idx_3d = (/2,1,1/)
  integer, parameter, dimension(3) :: perp2_idx_3d = (/3,3,2/)


  ! Types of density profiles

  ! These are already defined in m_species_define
  !integer, parameter :: p_uniform   = 0  ! uniform
  !integer, parameter :: p_gaussian  = 2  ! gaussian
  !integer, parameter :: p_channel   = 3  ! parabolic channel

  integer, parameter :: p_step     = 10
  integer, parameter :: p_thruster = 11


  ! declare things that should be public

  public :: t_cathode

  interface alloc
    module procedure alloc_cathode
    module procedure alloc_1d_cathode
    module procedure alloc_bound_1d_cathode
  end interface

  interface freemem
    module procedure freemem_cathode
    module procedure freemem_1d_cathode
  end interface

  public :: alloc, freemem

 contains

!---------------------------------------------------

!-----------------------------------------------------------------------------------------
! Read information from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_cathode( this, input_file, pspecies, spname, periodic, if_move, &
                             grid, dt, sim_options )

  use m_species_define
  use m_input_file
  use m_fparser
  use m_grid_define

  implicit none

  class( t_cathode ), intent(out)         :: this
  class( t_input_file ), intent(inout)    :: input_file

  class( t_species ), intent(inout), target :: pspecies

  character(len=*), intent(in)           :: spname
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ),               intent(in) :: grid
  real( p_double ), intent(in) :: dt

  type( t_options ), intent(in) :: sim_options

  integer :: dir   ! direction to inject particles 1, ..., p_x_dim
  integer :: wall  ! wall from which to inject particles

  logical :: deposit_current ! deposit current of injected particles

  ! time profile information
  real(p_double) :: t_start, t_rise, t_fall, t_flat

  ! moving injector plane
  logical :: mov_inj

  ! transverse density profile of the injected beam
  character(20)   :: prof_type              ! density profile type
  real(p_double) :: density                ! total(peak) density of the beam
  real(p_double), dimension(2) :: center       ! central position of the beam
  real(p_double) :: gauss_width                  ! total width to be considered
  real(p_double) :: gauss_w0                     ! gaussian waist to be considered

  real(p_double)               :: channel_bottom
  real(p_double)               :: channel_r0
  real(p_double)               :: channel_depth
  real(p_double)               :: channel_size
  real(p_double)               :: channel_wall

  real(p_double)               :: den_min

  namelist /nl_cathode/ dir, wall, t_start, t_rise, t_fall, t_flat, mov_inj, &
            prof_type, density, center, gauss_width, gauss_w0, &
            channel_bottom, channel_r0, channel_depth, &
            channel_size, channel_wall, deposit_current, den_min

  integer :: ierr, i

!       executable statements

  dir = 1  ! default inject particles along the x direction
  wall = 1 ! default inject particles from lower (left) wall

  deposit_current = .true. ! default deposit current of injected particles

  t_start = 0.0 ! by default, no delay for particle injection
  t_rise = -1.0 ! no default for rise time
  t_fall = -1.0 ! no default for fall time
  t_flat = -1.0 ! no default for flat time

  mov_inj = .false. ! standard wall injection

  prof_type = "uniform" ! defaults to uniform transverse density
  density = 1.0_p_double         ! defaults to 1.0 density
  center = -huge(1.0_p_double)       ! no default center
  gauss_width = -1.0_p_double          ! no default width
  gauss_w0 = -1.0_p_double             ! no default gauss_w0

  channel_r0     = -huge( 1.0_p_double )  ! no default r0
  channel_depth  = -huge( 1.0_p_double )  ! no default depth
  channel_size   = -huge( 1.0_p_double )  ! no default size
  channel_wall   = 0.0                     ! default wall size
  channel_bottom = 1.0

  den_min = 0.0

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_cathode", ierr )

  if (ierr /= 0) then
  if ( ierr < 0 ) then
    print *, "Error reading cathode parameters"
  else
    print *, "Error: cathode parameters missing"
  endif
    print *, "aborting..."
  stop
  endif

  read (input_file%nml_text, nml = nl_cathode, iostat = ierr)
  if (ierr /= 0) then
  print *, "Error reading cathode parameters"
  print *, "aborting..."
  stop
  endif

  ! check the direction and the wall parameters

  this%dir  = dir
  this%wall = wall

  if (this%dir<1 .or. this%dir>p_x_dim) then
  write(0,*) "(*error*) invalid direction for cathode. dir must be in the range [1..ND]"
  write(0,*) "(*error*) unable to setup cathode"
  write(0,*) "(*error*) bailing out..."
  call abort_program()
  endif

  if (this%wall<1 .or. this%wall>2) then
  write(0,*) "(*error*) invalid wall for cathode. wall must be 1 or 2"
  write(0,*) "(*error*) unable to setup cathode"
  write(0,*) "(*error*) bailing out..."
  call abort_program()
  endif

  this%den_min = den_min

  this%deposit_current = deposit_current

  ! check time profile parameters

  this%t_start  = t_start
  this%t_rise  = t_rise
  this%t_fall  = t_fall
  this%t_flat  = t_flat

  if (this%t_rise<0.0 .or. this%t_fall<0.0 .or. this%t_flat<0.0) then
  write(0,*) "(*error*) invalid time profile for cathode."
  write(0,*) "(*error*) t_fall, t_flat and t_rise must be specified and >= 0.0"
  write(0,*) "(*error*) bailing out..."
  call abort_program()
  endif

  this%mov_inj = mov_inj

  ! clear transverse profile for 1D runs
  if ( p_x_dim == 1 ) then
  prof_type = "uniform"
  endif

  this%density = density
  this%gauss_width = gauss_width
  this%gauss_w0 = gauss_w0
  this%center = center

  ! check transverse profile parameters
  select case ( trim( prof_type ) )
    case ( "uniform" )
    this%prof_type = p_uniform

    ! When using a uniform profile clear the center variable to avoid floating point
    ! exceptions if the user did not set this variable
        this%center = 0.0

    case ( "gaussian" )
    this%prof_type = p_gaussian
    this%gauss_width = gauss_width
    this%gauss_w0 = gauss_w0

    if (this%gauss_width <= 0.0) then
      write(0,*) "(*error*) invalid transverse profile for cathode."
      write(0,*) "(*error*) 'gauss_width' must be specified and > 0.0 when using ", &
           "the gaussian profile."
      write(0,*) "(*error*) bailing out..."
      stop
    endif

    if (this%gauss_w0 <= 0.0) then
      write(0,*) "(*error*) invalid transverse profile for cathode."
      write(0,*) "(*error*) 'gauss_w0' must be specified and > 0.0 when using the gaussian profile"
      write(0,*) "(*error*) bailing out..."
      stop
    endif

    case ( "step" )
    this%prof_type = p_step
        this%gauss_width = gauss_width

    if (this%gauss_width <= 0.0) then
      write(0,*) "(*error*) invalid transverse profile for cathode."
      write(0,*) "(*error*) 'gauss_width' must be specified and > 0.0 when using ", &
           "the step profile."
      write(0,*) "(*error*) bailing out..."
      stop
    endif

    case ( "thruster" )
    this%prof_type = p_thruster
    this%gauss_width = gauss_width
    this%gauss_w0 = gauss_w0

    if (this%gauss_width <= 0.0) then
      write(0,*) "(*error*) invalid transverse profile for cathode."
      write(0,*) "(*error*) 'gauss_width' must be specified and > 0.0 when using ", &
           "the thruster profile."
      write(0,*) "(*error*) bailing out..."
      stop
    endif

    if (this%gauss_w0 <= 0.0) then
      write(0,*) "(*error*) invalid transverse profile for cathode."
      write(0,*) "(*error*) 'gauss_w0' must be specified and > 0.0 when using the thruster profile"
      write(0,*) "(*error*) bailing out..."
      stop
    endif
    case ( "channel" )
    this%prof_type = p_channel

    ! set channel variables
    this%channel_r0 = channel_r0
    if ( this%channel_r0 <= 0.0 ) then
      write(0,*) "(*error*) Invalid channel_r0 value, must be > 0.0"
      write(0,*) "(*error*) aborting..."
      stop
    endif

    this%channel_depth  = channel_depth
    if ( this%channel_depth < 0.0 ) then
      write(0,*) "(*error*) Invalid channel_depth value, must be >= 0.0"
      write(0,*) "(*error*) aborting..."
      stop
    endif

    this%channel_size   = channel_size
    if ( this%channel_size <= 0.0 ) then
      write(0,*) "(*error*) Invalid channel diameter (channel_size) value"
      write(0,*) "(*error*) must be > 0.0"
      write(0,*) "(*error*) aborting..."
      stop
    endif

    this%channel_wall   = channel_wall
    this%channel_bottom = channel_bottom

    case default
    write(0,*) "(*error*) invalid transverse profile for cathode."
    write(0,*) "(*error*) bailing out..."
    call abort_program()
  end select


  if (this%prof_type /= p_uniform ) then
  if ( this%center(1) <= -huge(1.0_p_double) ) then
    write(0,*) "(*error*) invalid transverse profile for cathode."
    write(0,*) "(*error*) 'center' must be specified when using any profile other ", &
                          "than uniform."
    write(0,*) "(*error*) bailing out..."
    stop
  endif

  if (( this%center(2) <= -huge(1.0_p_double)) .and. ( p_x_dim > 2 )) then
    write(0,*) "(*error*) invalid transverse profile for cathode."
    write(0,*) "(*error*) 'center' must be specified for the two transverse directions ", &
               "when using any profile other than uniform."
    write(0,*) "(*error*) bailing out..."
    stop
  endif
  endif

  ! Read species
  if (disp_out(input_file)) then
    SCR_ROOT("    - reading associated species configuration...")
  endif
  call pspecies%read_input( input_file, spname, periodic, if_move, grid, dt, &
                  .false., sim_options )

  ! Validate velocity distribution
  if ( pspecies%udist%ufl_type /= p_uniform ) then
    write(0,*) "(*error*) Invalid fluid momentum for cathode species."
    write(0,*) "(*error*) Only uniform fluid momenta are supported."
    write(0,*) "(*error*) bailing out..."
  stop
  endif

  if ( (wall == p_lower .and. pspecies%udist%ufl(dir) <= 0 ) .or. &
       (wall == p_upper .and. pspecies%udist%ufl(dir) >= 0 ) ) then
   write(0,*) "(*error*) Invalid fluid momentum for cathode species."
   write(0,*) "(*error*) Injection fluid momentum must point into the box."
   write(0,*) "(*error*) bailing out..."
   stop
  endif

  ! Associate species with cathode
  this%species => pspecies

end subroutine read_input_cathode
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize object
!-----------------------------------------------------------------------------------------
subroutine init_cathode( this, cathode_id, restart, restart_handle, t, dt, coordinates )
!---------------------------------------------------
!       sets up this data structure from the given information
!---------------------------------------------------

   use m_restart
   use m_species_current

   implicit none

   ! dummy variables

   class( t_cathode ), intent(inout) :: this
   integer,     intent(in) :: cathode_id
   logical, intent(in) :: restart
   type( t_restart_handle ), intent(in) :: restart_handle

   real(p_double) :: t, dt
   integer :: coordinates

   ! local variables
   integer, dimension(2)            :: iperp

   integer, dimension(p_x_dim) :: lnum_par_x ! number of particles per cell

   ! Variables required for particle initialization inside simulation box
   real(p_double)  :: wall_disp
   integer         :: iwall_disp
   real(p_double), dimension(p_p_dim) :: u          ! fluid generalized momenta

   integer :: i, j, k
   logical :: inject_node

   ! executable statements

   ! get particle spacing
   lnum_par_x(1:p_x_dim) = this%species%num_par_x(1:p_x_dim)

   ! setup particle positions inside the cell
   select case (p_x_dim)
   case (1)
     continue
   case (2)
     iperp(1) = 3-this%dir

     call alloc( this%ppos_cell, (/ 1, lnum_par_x(iperp(1)) /) )

     do i = 1, lnum_par_x(iperp(1))
     this%ppos_cell(1,i) = (2*i - 1 - lnum_par_x(iperp(1)) ) * 0.5_p_k_part/lnum_par_x(iperp(1))
     enddo

   case (3)
     iperp(1) = perp1_idx_3d(this%dir)
     iperp(2) = perp2_idx_3d(this%dir)
     call alloc( this%ppos_cell, (/ 2, lnum_par_x(iperp(1))*lnum_par_x(iperp(2)) /))

     k = 1
     do i = 1, lnum_par_x(iperp(1))
     do j = 1, lnum_par_x(iperp(2))
       this%ppos_cell(1,k) = (2*i - 1 - lnum_par_x(iperp(1)) ) * 0.5_p_k_part/lnum_par_x(iperp(1))
       this%ppos_cell(2,k) = (2*j - 1 - lnum_par_x(iperp(2)) ) * 0.5_p_k_part/lnum_par_x(iperp(2))
       k = k + 1
     enddo
     enddo
   end select


   if ( restart ) then

      call restart_read_cathode( this, restart_handle )

   else

    ! set the internal cathode id
    this%cathode_id = cathode_id

    ! check if associated species is ok
    if (.not. associated(this%species)) then
        ERROR('species is not associated in setup_cathode')
    call abort_program(p_err_invalid)
    endif

    ! get the species id for restart information
    this%sp_id = this%species%sp_id

    ! get particle spacing in the injection direction
    lnum_par_x(1:p_x_dim) = this%species%num_par_x(1:p_x_dim)

    select case(this%wall)
    case (p_lower)
      this%inppos = 0
      this%nppos  = 0.5_p_double - 0.5_p_double /lnum_par_x(this%dir)
    case (p_upper)
      this%inppos = this%species%g_nx( this%dir ) + 1
      this%nppos  = 0.5_p_double /lnum_par_x(this%dir) - 0.5_p_double
    end select

    ! get fluid velocity from particle momenta for initial particle injection
    u = this%species%udist%ufl
    this%rgamma = 1.0_p_double/sqrt(1.0_p_double + u(1)**2 + u(2)**2 + u(3)**2)
    this%vel = u*this%rgamma / real( this%species%dx(this%dir), p_double )

    ! correct particle injection position with t_start (that can be positive or negative)
    if (this%t_start /= 0.0_p_double) then

      wall_disp   = this%vel(this%dir) * this%t_start
      iwall_disp  = floor( wall_disp )
      wall_disp   = wall_disp - iwall_disp ! wall_disp will be in ]-1,1[

      this%inppos = this%inppos - iwall_disp
      this%nppos  = this%nppos  - wall_disp

        call ntrim_pos( this%nppos, this%inppos )

    endif

    ! check if initial particle injection
    if (this%t_start < 0.0_p_double) then

        ! adjust multi-node injection
    select case(this%wall)
      case (p_lower)
      inject_node = ( this%inppos >= this%species%my_nx_p( p_lower, this%dir ) )

          case (p_upper)
      inject_node = ( this%inppos <= this%species%my_nx_p( p_upper, this%dir ) )

          case default
            inject_node = .false.
    end select

        ! inject only if necessary
        if ( inject_node ) then
      ! correct for zero time advance and inject particles
      this%nppos  = this%nppos  - this%vel(this%dir) * real( dt, p_double )

      call ntrim_pos( this%nppos, this%inppos )

      ! no need to deposit current
      call create_particles_cathode( this, t, dt, coordinates )
    endif

      endif

      this%t_total = this%t_start + this%t_rise + this%t_flat + this%t_fall

   endif

   if ( this%deposit_current ) call this % species % bnd_con % init_tmp_buf_current()

end subroutine init_cathode
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleanup dynamically allocated memory
!-----------------------------------------------------------------------------------------
subroutine cleanup_cathode( this )

  implicit none

  class( t_cathode ), intent(inout) :: this

  call freemem( this%ppos_cell )

end subroutine cleanup_cathode
!-----------------------------------------------------------------------------------------

!---------------------------------------------------
subroutine write_checkpoint_cathode( this, restart_handle )
!---------------------------------------------------
!       write object information into a restart file
!---------------------------------------------------

  use m_restart

  implicit none

!       dummy variables

  class( t_cathode ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

!       local variables
  character(len=*), parameter :: err_msg = 'error writing restart data for cathode object.'
  integer :: ierr

!       executable statements



  restart_io_wr( p_cath_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%cathode_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%dir, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%wall, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%deposit_current, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%inppos, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%nppos, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%vel, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%rgamma, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%t_start, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%t_rise, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%t_fall, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%t_flat, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%t_total, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%prof_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%density, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%center, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%gauss_width, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%gauss_w0, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%sp_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%channel_bottom, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%channel_r0, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%channel_depth, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%channel_size, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_wr( this%channel_wall, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )



end subroutine write_checkpoint_cathode
!---------------------------------------------------

!---------------------------------------------------
subroutine restart_read_cathode( this, restart_handle )
!---------------------------------------------------
!       read object information from a restart file
!---------------------------------------------------

  use m_restart

  implicit none

!       dummy variables

  type( t_cathode ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

!       local variables

  character(len=len(p_cath_rst_id)) :: rst_id
  character(len=*), parameter :: err_msg = 'error reading restart data for cathode object.'
  integer :: ierr

!       executable statements



  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  ! check if restart file is compatible
  if ( rst_id /= p_cath_rst_id) then
  ERROR('Corrupted restart file, or restart file ')
  ERROR('from incompatible binary (cathode)')
  call abort_program(p_err_rstrd)
  endif

  restart_io_rd( this%cathode_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%dir, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%wall, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%deposit_current, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%inppos, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%nppos, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%vel, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%rgamma, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%t_start, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%t_rise, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%t_fall, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%t_flat, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%t_total, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%prof_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%density, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%center, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%gauss_width, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%gauss_w0, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%sp_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%channel_bottom, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%channel_r0, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%channel_depth, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%channel_size, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%channel_wall, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )



end subroutine restart_read_cathode
!---------------------------------------------------



!---------------------------------------------------
function sp_id_cathode( this )
!---------------------------------------------------
!       returns the species id of the species
!       associated with this cathode
!---------------------------------------------------

  implicit none

!       dummy variables

  integer  :: sp_id_cathode

  type(t_cathode), intent(in) :: this

!       local variables

!       executable statements

  sp_id_cathode = this%sp_id

end function sp_id_cathode
!---------------------------------------------------


!-----------------------------------------------------------------------------------------
! Deposit current from injected cathode particles
!-----------------------------------------------------------------------------------------
subroutine deposit_current_cathode( this, ip0, ip1, jay, gdt )

   use m_vdf_define
   use m_species_current

   implicit none

   class(t_cathode), intent(inout)    :: this        ! cathode object
   integer, intent(in)               :: ip0, ip1    ! first and last particle injected
   type( t_vdf ),      intent(inout) :: jay         ! electrical current object
   real( p_double ),   intent(in)    :: gdt         ! time step

   ! local variables
   integer :: i, k, ip, npar
   real(p_k_part), dimension( p_x_dim ) :: dt_dx
   real(p_k_part), dimension( p_p_dim ) :: c

   ! Deposit current for cell type positions
   ! Due to the details of the current deposition routines we need to specify the original
   ! cell and position plus the new position of every motion. To avoid extra calculations
   ! we use the current particle position as the original data, and the original particle
   ! position as the new position. We then deposit the current using the opposite charge value
   ! to get the correct current. The sign of out of the plane components of the momentum has
   ! to be reversed to correct these current components (1D and 2D only).

   do i = 1, p_x_dim
   dt_dx(i) = gdt / this%species%dx(i)
   enddo

   ! Correct out of the plane components
   do i = 1, p_x_dim
     c(i) = 1.0_p_k_part
   enddo

   do i = p_x_dim+1, p_p_dim
     c(i) = -1.0_p_k_part
   enddo

   do ip = ip0, ip1, p_cache_size

    ! check if last copy of table and set np
    if( ip + p_cache_size > ip1 ) then
      npar = ip1 - ip + 1
    else
      npar = p_cache_size
    endif

    ! Get corrected particle momenta, corrected charge, and 1 / gamma
    do k = 1, npar
    tmp_p( 1, k ) = this%species%p(1, ip + k - 1 )
    tmp_p( 2, k ) = c(2) * this%species%p(2, ip + k - 1 )
    tmp_p( 3, k ) = c(3) * this%species%p(3, ip + k - 1 )

    tmp_q( k ) = -this%species%q( ip + k - 1 )

    tmp_rg( k ) = 1.0_p_k_part/sqrt(1.0_p_k_part + tmp_p(1,k)**2 + tmp_p(2,k)**2 + tmp_p(3,k)**2)
    enddo

    ! Get previous position and deposit current
    select case ( p_x_dim )
    case(1)

       do k = 1, npar
       tmp_xold(1,k)  = this%species%x(1,ip + k - 1) - tmp_p(1,k)*tmp_rg(k)*dt_dx(1)
       tmp_dxi(1,k) = ntrim( tmp_xold(1,k) )
       enddo

       call deposit_current_1d( this%species, jay, &
            tmp_dxi, tmp_xold, &
            this%species%ix(:, ip:), this%species%x(:, ip:), &
            tmp_q, tmp_rg, tmp_p, npar,  gdt)


    case(2)
       do k = 1, npar
       tmp_xold(1,k)  = this%species%x(1,ip + k - 1) - tmp_p(1,k)*tmp_rg(k)*dt_dx(1)
       tmp_xold(2,k)  = this%species%x(2,ip + k - 1) - tmp_p(2,k)*tmp_rg(k)*dt_dx(2)

       tmp_dxi(1,k) = ntrim( tmp_xold(1,k) )
       tmp_dxi(2,k) = ntrim( tmp_xold(2,k) )
       enddo

       call deposit_current_2d( this%species, jay, &
            tmp_dxi, tmp_xold, &
            this%species%ix(:, ip:), this%species%x(:, ip:), &
            tmp_q, tmp_rg, tmp_p, npar, gdt)

    case(3)
       do k = 1, npar
       tmp_xold(1,k)  = this%species%x(1,ip + k - 1) - tmp_p(1,k)*tmp_rg(k)*dt_dx(1)
       tmp_xold(2,k)  = this%species%x(2,ip + k - 1) - tmp_p(2,k)*tmp_rg(k)*dt_dx(2)
       tmp_xold(3,k)  = this%species%x(3,ip + k - 1) - tmp_p(3,k)*tmp_rg(k)*dt_dx(3)

       tmp_dxi(1,k) = ntrim( tmp_xold(1,k) )
       tmp_dxi(2,k) = ntrim( tmp_xold(2,k) )
       tmp_dxi(3,k) = ntrim( tmp_xold(3,k) )
       enddo

       call deposit_current_3d( this%species, jay, &
            tmp_dxi, tmp_xold, &
            this%species%ix(:, ip:), this%species%x(:, ip:), tmp_q, &
            npar, gdt)
    end select

   enddo


end subroutine deposit_current_cathode
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Inject particles from cathode
!-----------------------------------------------------------------------------------------
subroutine inject_cathode( this, jay, t, gdt, no_co, coords)

  use m_vdf_define
  use m_node_conf

  implicit none

  class(t_cathode), intent(inout)     :: this        ! cathode object
  type( t_vdf ),       intent(inout) :: jay         ! electrical current object
  real( p_double ),   intent(in)    :: t            ! simulation time
  real( p_double ),   intent(in)    :: gdt          ! time step
  class( t_node_conf ), intent(in)    :: no_co       ! node configuration
  integer :: coords ! coordinate system being used

  ! local variables
  integer                             :: np         ! number of particles before injection

  ! parallel run variables
  integer :: inj_node

  ! check if correct node for launching particles
  select case (this%wall)
  case (1)
    inj_node = 1
  case (2)
    inj_node = nx( no_co, this%dir )
  case default
    inj_node = -1
  end select

  if ( my_ngp( no_co, this%dir) == inj_node ) then

  ! check if injection is finished
  if (t >= this%t_start .and. t <= this%t_total) then

      ! get number of particles currently in the buffer
      np = this%species%num_par

    call create_particles_cathode( this, t, gdt, coords )

    ! deposit current corresponding to the injection

    if (this%species%num_par > np .and. this%deposit_current) then

    call deposit_current_cathode( this, np + 1, this%species%num_par, jay, gdt )

    endif

  endif ! (t <= this%t_total)

  endif ! injecting node


  ! call validate( this%species, "after cathode injection" )

end subroutine inject_cathode
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Inject particles from cathode, using low roundoff current deposition
!-----------------------------------------------------------------------------------------
subroutine inject_cathode_lowroundoff( this, jay, jay_tmp, t, gdt, no_co, coords)

  use m_vdf_define
  use m_vdf_math
  use m_node_conf

  implicit none

  class(t_cathode), intent(inout)     :: this        ! cathode object
  type( t_vdf ),       intent(inout) :: jay, jay_tmp         ! electrical current object
  real( p_double ),   intent(in)    :: t            ! simulation time
  real( p_double ),   intent(in)    :: gdt          ! time step
  class( t_node_conf ), intent(in)    :: no_co       ! node configuration
  integer :: coords ! coordinate system being used


  ! local variables
  integer                             :: np         ! number of particles before injection

  ! parallel run variables
  integer :: inj_node

  integer, dimension( 2, 3 ) :: range
  integer :: i_dim, i1, i2, i3



  ! check if correct node for launching particles

  select case (this%wall)
  case (1)
    inj_node = 1
  case (2)
    inj_node = nx( no_co, this%dir )
  case default
    inj_node = -1
  end select

  if ( my_ngp( no_co, this%dir) == inj_node ) then

  ! check if injection is finished
  if (t >= this%t_start .and. t <= this%t_total) then

      ! get number of particles currently in the buffer
      np = this%species%num_par

    call create_particles_cathode( this, t, gdt, coords )

    ! deposit current corresponding to the injection

    if ((this%species%num_par > np) .and. (this%deposit_current)) then

        ! reference implementation
        ! jay_tmp = 0.0_p_k_fld
        ! call deposit_current_cathode( this, np + 1, this%species%num_par, jay_tmp, gdt )
        ! call add( jay, jay_tmp )
        ! ~ reference implementation

    ! optimized implementation - only zero / add relevant cells
    do i_dim = 1, p_x_dim
      if ( i_dim == this%dir ) then
        select case ( this%wall )
              case(1)
                range(1,i_dim) = 1 - this%species%interpolation
                range(2,i_dim) = 1 + this%species%interpolation
              case(2)
                range(1,i_dim) = jay%nx_(i_dim) + 1 - this%species%interpolation
                range(2,i_dim) = jay%nx_(i_dim) + 1 + this%species%interpolation
        end select
      else
        range(1,i_dim) = jay%lbound( i_dim + 1 )
        range(2,i_dim) = jay%ubound( i_dim + 1 )
      endif
    enddo

    i_dim = this%dir

    ! zero the temporary grid
    select case( p_x_dim )
      case(1)
        do i1 = range(1,1),range(2,1)
          jay_tmp%f1(1:3,i1) = 0
        enddo
      case(2)
        do i2 = range(1,2),range(2,2)
        do i1 = range(1,1),range(2,1)
        jay_tmp%f2(1:3,i1, i2) = 0
        enddo
        enddo
      case(3)
        do i3 = range(1,3),range(2,3)
        do i2 = range(1,2),range(2,2)
        do i1 = range(1,1),range(2,1)
          jay_tmp%f3(1:3,i1, i2, i3) = 0
        enddo
        enddo
        enddo
    end select

    ! Deposit current on temporary grid
    call deposit_current_cathode( this, np + 1, this%species%num_par, jay_tmp, gdt )

    ! add the values to electric current grid in the relevant range
    call add( jay, jay_tmp, range )

    endif

  endif ! (t <= this%t_total)

  endif ! injecting node


 ! call validate( this%species, "after cathode injection (low roundoff)" )


end subroutine inject_cathode_lowroundoff
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Simulation window move
!-----------------------------------------------------------------------------------------
subroutine move_window_cathode( this, g_space )

  use m_space

  implicit none

  class(t_cathode), intent(inout)     :: this
  type( t_space ), intent(in) :: g_space

  this%inppos = this%inppos - nx_move( g_space, this%dir )

end subroutine move_window_cathode
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Advance cathode particle position, put cathode particles that are inside simulation space
! in the species buffer and set their momenta
!-----------------------------------------------------------------------------------------
subroutine create_particles_cathode( this, t, dt, coords )

  use m_species_define
  use m_species_udist
  use m_fparser

  implicit none

  ! dummy variables
  class(t_cathode), intent(inout)     :: this        ! cathode object
  real( p_double ),   intent(in)    :: t           ! simulation time
  real( p_double ),   intent(in)    :: dt          ! time step
  integer :: coords                                 ! coordinate system being used


  ! local variables
  real(p_k_part)                      :: dpx        ! particle spacing along the injection direction
  real(p_k_part)                      :: lden_min   ! den_min for the species to inject
  integer                             :: np         ! number of particles to inject
                          ! along the injection direction

  integer :: i, npar
  real(p_k_part)                       :: qnorm     ! normalization factor for injected charge
  real(p_k_part), dimension(p_max_dim) :: xnewpart  ! position,
  integer, dimension(p_max_dim)        :: ixnewpart ! momenta and
  real(p_k_part)                       :: qnewpart  ! charge of the new particle

  real(p_k_part)                         :: den
  real(p_k_fparse), dimension(p_max_dim) :: eval_var

  integer :: ptrcur

  integer, dimension(2)               :: iperp
  real(p_k_part), dimension(p_x_dim)  :: gpos_min, dx, dxp
  real(p_k_part), dimension(2)        :: x_cell_perp
  real(p_k_part) :: inj_time, inj_pos

  integer :: j, j1, j2, k

  integer :: npi                        ! total number of particles to inject
  integer, dimension( 2, 3 ) :: range   ! grid range in which to inject

  lden_min = this%den_min
  dx(1:p_x_dim) = real( this%species%dx(1:p_x_dim), p_k_part )

  ! advance the position of the next particle to inject
  this%nppos  = this%nppos + this%vel(this%dir)*real( dt, p_k_part )
  call ntrim_pos( this%nppos, this%inppos )

  ! get number of particles in the buffer
  np = this%species%num_par

  ! get normalization factor for injected charge
  qnorm = sign( 1.0_p_k_part, this%species%rqm )
  do i=1, p_x_dim
  qnorm = qnorm/this%species%num_par_x(i)
  enddo

  ! Get perpendicular directions information
  select case ( p_x_dim )
    case(1)
      continue
    case(2)
      iperp(1) = 3 - this%dir
    case(3)
      iperp(1) = perp1_idx_3d(this%dir)
      iperp(2) = perp2_idx_3d(this%dir)
    case default
      iperp = 0
  end select

  npi = 1
  do i = 1, p_x_dim - 1
    range( p_lower, i ) = this%species%my_nx_p(p_lower,iperp(i))
    range( p_upper, i ) = this%species%my_nx_p(p_upper,iperp(i))
    gpos_min( i ) = this%species%g_box( p_lower, iperp(i) )
    dxp( i ) = dx(iperp(i))

    npi = npi * this%species%num_par_x(iperp(i))
  enddo

  ! Get particle spacing along injection direction
  ! This is now done in cell normalized units
  dpx = 1.0 / this%species%num_par_x(this%dir)

  do

  ! Get longitudinal beam density
  select case (this%wall)
    case (p_lower)
    if (this%inppos < this%species%my_nx_p( p_lower, this%dir) ) exit

    inj_time = real(t, p_k_part) - this%t_start - &
               (this%inppos - this%nppos - 1)/this%vel(this%dir)

    case (p_upper)
    if (this%inppos > this%species%my_nx_p( p_upper, this%dir)) exit

    inj_time = real(t, p_k_part) - this%t_start - &
               ( this%inppos - this%species%my_nx_p( p_upper, this%dir) + &
                                    this%nppos) / this%vel(this%dir)
  end select

  if ( this%inppos >= this%species%my_nx_p( p_lower, this%dir) .and. &
       this%inppos <= this%species%my_nx_p( p_upper, this%dir) ) then

     ! the total_xmoved term corrects for motion of injecting wall (moving window)
     ! after injection started
     inj_time = inj_time + real( this%species%total_xmoved( this%dir ), p_k_part ) - &
               this%t_start

     ! multiply by global density
     den = this%density * t_env( inj_time, this%t_rise, this%t_flat, this%t_fall )

     ! inject particles
     xnewpart (this%dir ) = this%nppos
     ixnewpart(this%dir)  = this%inppos - this%species%my_nx_p( p_lower, this%dir ) + 1

     select case (p_x_dim)

     case (1) ! Inject in 1D

      qnewpart = den

      if (abs(qnewpart) > lden_min) then

        ! add particle to the buffer
        call this%species%create_particle( ixnewpart, xnewpart, qnorm*qnewpart )

      endif

     case (2) ! Inject in 2D

      do j = range(p_lower, 1), range( p_upper, 1)
         x_cell_perp(1) = gpos_min(1) + (j-1)*dxp(1)

         do i = 1, npi
         xnewpart(iperp(1)) = x_cell_perp(1) + this%ppos_cell(1,i)*dxp(1)
         ! set the charge of the new particle
         qnewpart = den * tr_profile( this, real( xnewpart(iperp(1)) - this%center(1), p_k_part ) )

         if (abs(qnewpart) > lden_min) then
           ! correct for cylindrical coords.
           if (coords == p_cylindrical_b) then
          qnewpart = qnewpart*abs(xnewpart(p_r_dim))
           endif

           ixnewpart(iperp(1)) = j - this%species%my_nx_p( 1, iperp(1) ) + 1
           xnewpart(iperp(1))  = this%ppos_cell(1,i)

           ! add particle to the buffer
           call this%species%create_particle( ixnewpart, xnewpart, qnorm*qnewpart )

         endif
         enddo
      enddo

     case(3)

      do j1 = range(p_lower, 1), range( p_upper, 1)
        x_cell_perp(1) = gpos_min(1) + (j1-1)*dxp(1)

        do j2 = range(p_lower, 2), range( p_upper, 2)
         x_cell_perp(2) = gpos_min(2) + (j2-1)*dxp(2)

         do i = 1, npi
           xnewpart(iperp(1)) = x_cell_perp(1) + this%ppos_cell(1,i)*dxp(1)
           xnewpart(iperp(2)) = x_cell_perp(2) + this%ppos_cell(2,i)*dxp(2)

           ! set the charge of the new particle
           qnewpart = den* tr_profile( this,  real( &
                sqrt((xnewpart(iperp(1))-this%center(1))**2 + &
                   (xnewpart(iperp(2))-this%center(2))**2), p_k_part ))

           if (abs(qnewpart) > lden_min) then

           ixnewpart(iperp(1)) = j1 - this%species%my_nx_p(p_lower,iperp(1)) + 1
           xnewpart (iperp(1)) = this%ppos_cell(1,i)
           ixnewpart(iperp(2)) = j2 - this%species%my_nx_p(p_lower,iperp(2)) + 1
           xnewpart (iperp(2)) = this%ppos_cell(2,i)

           ! add particle to the buffer
           call this%species%create_particle( ixnewpart, xnewpart, qnorm*qnewpart )
           endif
         enddo

        enddo
      enddo

     end select

    endif

  ! process next set of particles
  select case (this%wall)
    case (p_lower)
    this%nppos  = this%nppos - dpx

    ! Check injection halt for when using the moving injector plane
    if( this%mov_inj) then

      inj_pos = this%species%g_box( p_lower, this%dir ) + (this%inppos-1)*dx(this%dir)
      if (inj_pos < int(-this%t_start - t) ) exit

    endif

    case (p_upper)
    this%nppos = this%nppos + dpx

    ! Check injection halt for when using the moving injector plane
    if( this%mov_inj ) then

      inj_pos = this%species%g_box( p_lower, this%dir ) + (this%inppos-1)*dx(this%dir)
          if( inj_pos > int(this%species%g_box(p_upper, this%dir) + this%t_start + t) ) exit

    endif

  end select

  call ntrim_pos( this%nppos, this%inppos )

  enddo

  ! Set momentum of injected particles
   if (this%species%num_par > np) then

   ptrcur = np+1
   npar = this%species%num_par-np

   ! set momentum of injected particles
   call set_momentum( this%species, ptrcur, ptrcur + npar - 1 )

   ! process transversal momentum parser
   eval_var(1) = real( t, p_k_fparse )

  endif

  ! ------------------------------------------

  contains

  ! function declaration of time envelope

  function t_env( t, rise, flat, fall )
  implicit none
  real(p_k_part) :: t_env
  real(p_k_part), intent(in) :: t, rise, flat, fall

  if (t < 0.0d0 .or. t > rise+flat+fall) then
    t_env = 0.0
  else if (t < rise) then
    t_env = poly_env(t/rise)
  else if ( t < rise+flat ) then
    t_env = 1.0
  else
    t_env = poly_env((rise+flat+fall-t)/fall)
  endif

  end function t_env


  ! gaussian like polynomial for time envelope

  function poly_env(x)
  implicit none
  real(p_k_part) :: poly_env
  real(p_k_part), intent(in) :: x

  poly_env = 10 * x**3 - 15 * x**4 + 6 * x**5
  end function poly_env


  ! transverse profile functions

  function tr_profile( cathode, r )
  implicit none
  real(p_k_part) :: tr_profile
  type( t_cathode ), intent(in) :: cathode
  real(p_k_part)  , intent(in) :: r
!            real(p_k_part) :: x

  select case (cathode%prof_type)
    case(p_uniform)
    tr_profile = 1.0_p_k_part
    case(p_step)
    if ( abs(r) <= cathode%gauss_width/2.0 ) then
       tr_profile = 1.0_p_k_part
    else
       tr_profile = 0.0_p_k_part
    endif
    case(p_gaussian)
    if ( abs(r) <= cathode%gauss_width/2.0) then
!                  x = r
      tr_profile = exp(-2.0*(abs(r)/cathode%gauss_w0)**2)
    else
      tr_profile = 0.0_p_k_part
    endif
   case(p_thruster)
    if (( abs(r) <= cathode%gauss_width/2.0_p_k_part ) .and. &
      ( abs(r) > cathode%gauss_w0/2.0_p_k_part )) then
      tr_profile = 1.0_p_k_part
    else
      tr_profile = 0.0_p_k_part
    endif
    case (p_channel)

       if (abs(r) <= this%channel_size/2) then ! inside the channel
      tr_profile = this%channel_bottom + &
            this%channel_depth * &
            (abs(r)/this%channel_r0)**2
    else                                 ! outside the channel

       if (this%channel_wall <= 0.0) then ! finite channel
       tr_profile = this%channel_bottom + &
            this%channel_depth * &
            ((this%channel_size/2)/this%channel_r0)**2

       else                              ! leaky channel

       if (abs(r) < (this%channel_size/2 + this%channel_wall )) then
        tr_profile = (-abs(r) + this%channel_size/2 + this%channel_wall) / &
              this%channel_wall * (this%channel_bottom + &
              this%channel_depth * &
            ((this%channel_size/2)/this%channel_r0)**2)
       else
         tr_profile = 0.0d0
       endif

       endif

    endif

    case default

        tr_profile = 1.0_p_k_part

  end select

  end function tr_profile

end subroutine create_particles_cathode
!-----------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Returns the integer shift (-1, 0 or +1) so that the coordinate remains in the [-0.5, 0.5[
! range. This is the fastest implementation (twice as fast as a sequence of ifs) because
! the two if structures compile as conditional moves and can be processed independently.
!---------------------------------------------------------------------------------------------------
function ntrim(x)
!---------------------------------------------------------------------------------------------------
  implicit none

  real(p_k_part), intent(in) :: x
  integer :: ntrim, a, b

  if ( x < -.5 ) then
  a = -1
  else
    a = 0
  endif

  if ( x >= .5 ) then
  b = +1
  else
    b = 0
  endif

  ntrim = a+b

end function ntrim
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine ntrim_pos( x, ix )

  implicit none

  real(p_k_part), intent(inout) :: x
  integer, intent(inout) :: ix

  integer :: dx

  dx = ntrim(x)
  x  = x - dx
  ix = ix + dx

end subroutine ntrim_pos
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines

#define __TYPE__ type( t_cathode )
#define __TYPE_STR__ "t_cathode"
#define FNAME( a )  a ## _cathode
#define __MAX_DIM__ 1
#include "memory/mem-template.h"

!---------------------------------------------------------------------------------------------------

!---------------------------------------------------
      end module m_cathode
