! #define DEBUG_GLOBAL
! #define DEBUG_FILE 1

module m_simulation

#include "os-config.h"
#include "os-preprocess.fpp"

use m_system
use m_parameters

use m_node_conf
use m_grid_define
use m_grid
use m_grid_parallel

use m_time_step
use m_restart
use m_system

use m_space
use m_time
use m_emf_define
use m_particles_define
use m_bnd

use m_zpulse
use m_current_define
use m_antenna_array

use m_input_file

implicit none

private

type :: t_simulation

   type( t_time_step )       ::  tstep
   type( t_restart   )       ::  restart
   type( t_options   )       ::  options

   type( t_space        ), pointer  ::  g_space     => null()
   type( t_time         )           ::  time

   class( t_input_file  ), pointer  ::  input_file  => null()
   class( t_node_conf   ), pointer  ::  no_co       => null()
   class( t_grid        ), pointer  ::  grid        => null()
   class( t_zpulse_list ), pointer  ::  zpulse_list => null()
   class( t_current     ), pointer  ::  jay         => null()
   class( t_emf         ), pointer  ::  emf         => null()
   class( t_particles   ), pointer  ::  part        => null()
   class(t_antenna_array), pointer  ::  antenna_array => null()

   ! Boundary information
   class( t_bnd ), pointer :: bnd => null()

contains

   procedure :: iter               => iter_sim
   procedure :: iter_finished      => iter_finished_sim

   procedure :: allocate_objs      => allocate_objs_sim
   procedure :: init               => init_sim
   procedure :: cleanup            => cleanup_sim
   procedure :: read_input         => read_input_sim

   procedure :: write_checkpoint   => write_checkpoint_sim
   procedure :: list_algorithm     => list_algorithm_sim
   procedure :: test_input         => test_input_sim

   procedure :: set_int_load       => set_int_load_sim
   procedure :: dlb                => dlb_sim
   procedure :: report_load        => report_load_sim
   procedure :: min_gc_sim         => min_gc_sim
end type t_simulation

interface
  subroutine dlb_sim( this )
    import t_simulation
    class(t_simulation), intent(inout) :: this
  end subroutine
end interface

interface
  subroutine report_load_sim( this )
    import t_simulation
    class(t_simulation), intent(inout) :: this
  end subroutine
end interface

! variables to time program execution
integer :: loopev

integer :: dynlbev               ! load-balance (total)
integer :: dynlb_set_int_load_ev ! load-balance, get current int load
integer :: dynlb_new_lb_ev       ! load-balance, new lb configuration
integer :: dynlb_reshape_ev      ! load-balance, reshape simulation
integer :: dynlb_reshape_part_ev ! load-balance, reshape simulation (particles)
integer :: dynlb_reshape_grid_ev ! load-balance, reshape simulation (grids)

integer :: restart_write_ev      ! writing restart information
integer :: report_load_ev        ! report load diagnostic

! exported symbols
public :: t_simulation
public :: loopev, dynlbev, dynlb_set_int_load_ev, dynlb_new_lb_ev, dynlb_reshape_ev
public :: dynlb_reshape_part_ev, dynlb_reshape_grid_ev, restart_write_ev, report_load_ev

contains

!-----------------------------------------------------------------------------------------
! Read simulation parameters from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_sim( sim, input_file )

  implicit none

  class ( t_simulation ), intent(inout) :: sim
  class( t_input_file ), intent(inout) :: input_file

  real( p_double ), dimension( p_x_dim ) :: ldx

  SCR_ROOT('Reading parallel node configuration... ')
  call sim % no_co % read_input( input_file, p_x_dim  )

  SCR_ROOT('Reading grid configuration...')
  call sim % grid % read_input( input_file, sim % no_co % nx  )

  SCR_ROOT('Reading tstep configuration... ')
  call read_nml(  sim%tstep, input_file    )

  SCR_ROOT('Reading restart configuration... ')
  call read_nml(  sim%restart, input_file, sim%options%restart  )

  SCR_ROOT('Reading g_space configuration... ')
  call read_nml(  sim%g_space, input_file, periodic(sim%no_co), sim%grid%coordinates  )

  SCR_ROOT('Reading time configuration... ')
  call read_nml(  sim%time, input_file     )

  SCR_ROOT('Reading emf configuration... ')
  ! determine cell size
  ldx = (xmax( sim%g_space ) - xmin( sim%g_space ))/ sim%grid%g_nx( 1:p_x_dim )
  call sim % emf % read_input( input_file, periodic(sim%no_co), if_move(sim%g_space), &
                                sim%grid, ldx, dt(sim%tstep), sim%options%gamma )

  SCR_ROOT('Reading part configuration... ')
  call sim%part%read_input( input_file,  periodic(sim%no_co), &
                  if_move(sim%g_space), sim%grid, dt(sim%tstep), &
                  sim%options )

  SCR_ROOT('Reading zpulses ... ')
  call sim % zpulse_list % read_input( input_file,  sim%g_space, sim%emf%bnd_con, &
                 periodic(sim%no_co), sim%grid, sim%options )

  SCR_ROOT('Reading current configuration...')
  call sim % jay % read_input( input_file )

  SCR_ROOT('Reading antenna configuration...')
  call sim % antenna_array % read_input( input_file )

end subroutine read_input_sim


!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_sim( sim )

  implicit none

  class( t_simulation ), intent(inout) :: sim


  ! Allocate default space object class
  if ( .not. associated( sim%g_space ) ) then
    allocate( sim%g_space )
  endif

  ! Allocate default grid object class
  if ( .not. associated( sim%grid ) ) then
    allocate( sim%grid )
  endif

  ! Allocate default emf object class
  if ( .not. associated( sim%emf ) ) then
    allocate( sim%emf )
  endif

  ! Allocate default part object class
  if ( .not. associated( sim%part ) ) then
    allocate( sim%part )
  endif

  ! Allocate default current object class
  if ( .not. associated( sim%jay ) ) then
    allocate( sim%jay )
  endif

  ! Allocate default antenna array object class
  if ( .not. associated( sim%antenna_array) ) then
    allocate( sim%antenna_array )
  endif

  ! Allocate default zpulse object class
  if ( .not. associated( sim%zpulse_list) ) then
    allocate( sim%zpulse_list )
  endif

  ! Allocate default no_co object class
  if ( .not. associated( sim%no_co) ) then
    allocate( sim%no_co )
  endif

  ! Allocate boundary class
  if ( .not. associated( sim%bnd) ) then
    allocate( sim%bnd )
  endif

  ! If needed, 'allocate_objs' on the classes we just made so that any subobjects are created
  !   sim%emf, sim%part, sim%jay do object allocation in their 'read_input' proceedures
  !call sim%emf%allocate_objs()
  !call sim%part%allocate_objs()
  !call sim%jay%allocate_objs()

end subroutine allocate_objs_sim

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine init_sim( sim, grid_gc_min )

  use m_grid_parallel
  use m_random
  use m_logprof
  use m_particles
  use m_emf
  use m_diagnostic_utilities, only : init_dutil

  implicit none

  ! dummy variables

  class ( t_simulation ), intent(inout) ::  sim



  real( p_double ), dimension( p_x_dim ) :: ldx

  integer :: i_dim
  integer, dimension( 2, p_x_dim ) :: min_gc, min_gc_tmp
  integer, dimension( 2, p_x_dim ), intent(in) :: grid_gc_min


  logical :: restart
  type( t_restart_handle   )    ::  restart_handle
  integer :: diag_buffer_size

  ! initialize timing structure
  loopev                = create_event('loop')
  dynlbev               = create_event('dynamic load balance (total)')
  dynlb_set_int_load_ev = create_event('dynlb get int load')
  dynlb_new_lb_ev       = create_event('dynlb get new lb')
  dynlb_reshape_ev      = create_event('dynlb reshape simulation')
  dynlb_reshape_part_ev = create_event('dynlb reshape simulation (particles)')
  dynlb_reshape_grid_ev = create_event('dynlb reshape simulation (grids)')
  restart_write_ev      = create_event('writing restart information')
  report_load_ev        = create_event('load diagnostics')

  ! use restart files if restart is set in the input file
  restart = if_restart_read(sim%restart)

  SCR_ROOT('')
  if ( restart ) then
    SCR_ROOT('**********************************************')
    SCR_ROOT('* Restarting simulation from checkpoint data *')
    SCR_ROOT('**********************************************')
  else
    SCR_ROOT('Initializing simulation:')
  endif

  SCR_ROOT(" Setting up no_co")
  call sim % no_co % init()

  if ( restart ) then
    ! open restart file (this needs to happen after setup no_co)
    call restart_read_open( sim%restart, comm(sim%no_co), sim%no_co%ngp_id(), file_id_rst, restart_handle )

    ! read restart information for no_co
    ! This only checks for inconsistencies in the input file vs. restart files, all setup
    ! is completed above
    call sim % no_co % restart_read( restart_handle )
  endif
  
  ! setup the simulation grid
  SCR_ROOT(" Setting up grid")
  
  call sim%grid%init( sim%no_co, grid_gc_min, x_bnd( sim%g_space ),  &
              restart, restart_handle )
  
  ! if not restarting generate initial load balance partition
  if ( .not. restart ) then

     if ( sim % grid % needs_int_load( 0 ) ) then
       SCR_ROOT(' - generating initial integral load...')
       ! create even partition for calculating particle (and grid) loads
       call sim % grid % parallel_partition( sim%no_co, -2 )
       call sim % set_int_load( sim%grid, 0 )
     endif

     SCR_ROOT(' - setting parallel grid boundaries...')
     call sim % grid % parallel_partition( sim%no_co, 0 )
     call clear_int_load( sim%grid )

  endif
  
  ! Initialize I/O object
  ! Hi Ricardo Look here! (#199)
  ! (The following line was added)
  call sim % grid % io % init( sim%grid, sim%no_co )

  ! initialization of the random number seed
  SCR_ROOT(" Setting up random number seed")
  call init_genrand( sim%options%random_class, sim%options%enforce_rng_constancy, &
                     sim%options%random_seed, sim%no_co%ngp_id() + sim%options%random_seed, &
                     sim%grid, restart, restart_handle  )

  ! determine cell size
  ldx = (xmax( sim%g_space ) - xmin( sim%g_space ))/ sim%grid%g_nx( 1:p_x_dim )

  SCR_ROOT(" Setting up tstep")
  call setup( sim%tstep, -1 , restart, restart_handle )

  SCR_ROOT(" Setting up restart")
  call setup( sim%restart, restart, restart_handle )

  SCR_ROOT(" Setting up g_space")
  call setup( sim%g_space, ldx, sim%grid%coordinates , restart, restart_handle )

  SCR_ROOT(" Setting up time")
  call setup( sim%time, tmin(sim%time) - dt(sim%tstep), &
              dt(sim%tstep) , restart, restart_handle )

  SCR_ROOT(" Setting up boundary buffers")
  call sim % bnd % init()

  SCR_ROOT(" Setting up emf")

  min_gc = get_min_gc( sim%part )

  call sim%emf%init( get_grid_center( sim%part ), interpolation( sim%part ), &
                sim%g_space, sim%grid, min_gc, ldx, sim%tstep, tmax(sim%time), sim%no_co, &
                sim%bnd%send_vdf, sim%bnd%recv_vdf, restart, restart_handle, sim%options )

  SCR_ROOT(" Setting up current")

  ! The field solver may require than particles ( but the cells for emf smoothing don't
  ! need to be taken into account )
  min_gc_tmp = sim%emf%min_gc()

  do i_dim = 1, p_x_dim
    if ( min_gc_tmp( 1, i_dim ) > min_gc( 1, i_dim ) ) min_gc( 1, i_dim ) = min_gc_tmp( 1, i_dim )
    if ( min_gc_tmp( 2, i_dim ) > min_gc( 2, i_dim ) ) min_gc( 2, i_dim ) = min_gc_tmp( 2, i_dim )
  enddo

  call sim%jay%init( sim%grid, min_gc, ldx, &
              sim%no_co, if_move(sim%g_space), &
              sim % emf % bc_type(), interpolation( sim%part ), &
              restart, restart_handle, sim%options )

  SCR_ROOT(" Setting up particles")
  call sim%part%init( sim%g_space, sim%jay, sim%emf, &
              sim%grid, sim%no_co, sim%bnd, ndump(sim%tstep) , &
              restart, restart_handle, t(sim%time), sim%tstep, tmax(sim%time), sim%options  )

  SCR_ROOT(" Setting up zpulses")
  call sim % zpulse_list % init( restart, t(sim%time), dt(sim%tstep), sim%emf, sim%g_space, &
                                 sim%grid%my_nx(p_lower, :), sim%no_co, &
                                 interpolation( sim%part ) )

  SCR_ROOT(" Setting up antenna array")
  call sim % antenna_array % init( restart, restart_handle )


  SCR_ROOT(" Setting up diagnostics utilities")
  diag_buffer_size = 0
  call sim % part % get_diag_buffer_size( sim%grid%g_nx, diag_buffer_size )
  call sim % emf  % get_diag_buffer_size( sim%grid%g_nx, diag_buffer_size )
  call sim % jay  % get_diag_buffer_size( sim%grid%g_nx, diag_buffer_size )
  call init_dutil( sim%options, diag_buffer_size )

  if ( restart ) then
    call restart_read_close( sim%restart, restart_handle )
  endif

  SCR_ROOT('')
  SCR_ROOT('Initialization complete!')

end subroutine init_sim

!-----------------------------------------------------------------------------------------
! Cleanup simulation objects
!-----------------------------------------------------------------------------------------
subroutine cleanup_sim( sim )

  use m_vdf_comm
  use m_wall_comm
  use m_particles
  use m_diagnostic_utilities

  implicit none

  class( t_simulation ), intent(inout) :: sim

  call sim % grid % cleanup()

  ! cleanup particles
  call sim % part % cleanup( )

  ! cleanup current object
  call sim%jay%cleanup()

  ! cleanup emf object
  call sim%emf%cleanup()

  ! cleanup laser pulses
  call sim % zpulse_list % cleanup( )

  ! cleanup antennas
  call sim% antenna_array % cleanup( )

  ! cleanup vdf and spec comm buffers
  call sim % bnd % cleanup()

  ! cleanup diagnostic utilities
  call cleanup_dutil()

  ! Cleanup MPI system
  call sim % no_co % cleanup()

  ! deallocate polymorphic objects
  deallocate( sim%input_file )
  deallocate( sim%no_co )
  deallocate( sim%grid )
  deallocate( sim%emf )
  deallocate( sim%part )
  deallocate( sim%jay )
  deallocate( sim%zpulse_list )
  deallocate( sim%bnd )

end subroutine cleanup_sim


!-----------------------------------------------------------------------------------------
! Called at the end of an simulation iteration (i.e. after rach invokation of 'iter_sim')
!-----------------------------------------------------------------------------------------
subroutine iter_finished_sim( sim )

  use m_emf_define
  use m_current_define
  use m_particles_define
  use m_species_define

  implicit none
  
  class ( t_simulation ), intent(inout) :: sim

  class(t_species), pointer :: species

  ! 
  !   Increment these values in one central palce versus expllictly adding
  !         n_current_xxx = n_current_xxx + 1
  !    to all the relevent areas in all the various subclasses... we may
  !    need to go down the that route, however, if we need more fine grained 
  !    inter-iteration control of these flags.
  !

  ! signal that jay is up-to-date and pret-a-porter
  sim%jay%n_current_jay = sim%jay%n_current_jay + 1

  ! signal that the fields are also good to go.
  sim%emf%n_current_emf = sim%emf%n_current_emf + 1

  ! traverse all species object and also signal that their particle data is up-tp-date.
  species => sim%part%species
  do
    if (.not. associated(species)) exit
    species % n_current_spec = species % n_current_spec + 1
    species => species % next
  enddo

end subroutine

!-----------------------------------------------------------------------------------------
! Simulation iteration
!-----------------------------------------------------------------------------------------
subroutine iter_sim( sim )

  use m_emf
  use m_particles
  use m_logprof
  use m_particles_charge
  use m_space
  use m_vdf_define

  implicit none

  class ( t_simulation ), intent(inout) :: sim

  ! set the initial timer for the loop
  call begin_event(loopev)

  ! advance time step counter and time
  DEBUG('Before advance time')
  call advance( sim%tstep, 1 )
  call update( sim%time, dt(sim%tstep), n(sim%tstep) )

  if ( root(sim%no_co) ) then
     print '(A,I10,A,F0.3)', 'n = ', n(sim%tstep), ', t = ', t(sim%time)
  endif

  ! Electric current diagnostics  must occur before dynamic load
  ! balancing, since it destroys the current values
  DEBUG('Before report current')
  call sim%jay%report(sim%g_space, sim%grid, sim%no_co, sim%tstep, t(sim%time) )

  DEBUG('Before zpulse launch')
  call sim%zpulse_list%launch( sim%emf, sim%g_space, &
               sim%grid%my_nx(p_lower, :), sim%grid%g_nx, &
               t(sim%time), dt(sim%tstep), sim%no_co )

  ! ********************* move window ******************************

  ! test if simulation window needs to be moved - use
  ! global space to avoid numerical inconsistencies
  ! between move on different nodes
  DEBUG('Before move window sim%g_space')
  call move_boundary( sim%g_space, sim % emf % dx(), dt(sim%tstep) )

  ! move boundaries of other objects as required
  DEBUG('Before move window sim%emf')
  call sim%emf % move_window( sim%g_space, sim%grid%my_nx( p_lower, : ), &
                    sim%no_co%on_edge( 1, p_upper ) )

  DEBUG('Before move window sim%part')
  call sim%part % move_window( sim%g_space, sim%grid, &
                     sim%jay, sim%no_co, sim%bnd )

  ! jay does not need to be moved since it is recalculated
  ! from scratch during each interation (in advance_deposit)
  ! use antenna if necessary
  DEBUG('Before antenna')
  if (antenna_on(sim%antenna_array)) then
     call sim % antenna_array % launch( sim%emf, dt(sim%tstep), t(sim%time), &
                                        sim%g_space, sim%grid%my_nx( p_lower, :) )
  end if
  ! update boundaries of em-fields and calculate
  ! time averaged fields if necessary
  DEBUG('Before update boundary sim%emf')
  call sim%emf%update_boundary( sim%g_space, sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )

  ! Update the fields used for particle interpolation. This includes spatial smoothing,
  ! external fields and subcycling
  DEBUG('Before update_particle_fld sim%emf')
  call sim%emf%update_particle_fld( n(sim%tstep), t(sim%time), dt(sim%tstep), sim%g_space, &
                            sim%grid%my_nx( p_lower, : ))

  ! diagnostic of em-field and particles
  ! currents can only be diagnosed after advance_deposit and are taken care of at the
  ! beggining of the loop
  DEBUG('Before report emf')
  if ( sim % emf % if_charge_cons( sim%tstep ) ) then
    call update_charge( sim%part, sim%g_space, sim%grid, sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )
    ERROR( "Charge conservation diagnostic is temporarily disabled" )
    call abort_program()
    ! call update_charge_cons( sim%emf, sim%part%charge%current, sim%bnd%send_vdf, sim%bnd%recv_vdf )
  endif

  call sim % emf % report( sim%g_space, sim%grid,  sim%no_co, sim%tstep, t(sim%time), sim%bnd%send_vdf, sim%bnd%recv_vdf )

  DEBUG('Before report part')
  call sim % part % report( sim%emf, sim%g_space, sim%grid, sim%no_co, &
                    sim%tstep, t(sim%time), sim%bnd%send_vdf, sim%bnd%recv_vdf )

  ! move particles a timestep further and get new current
  DEBUG('Before advance_deposit')
  call sim % part % advance_deposit( sim%emf, sim%jay, sim%tstep, t(sim%time), &
                                     sim%no_co, sim%options )

  ! take care of boundaries for particles
  DEBUG('Before update_boundary sim%part')
  call sim % part % update_boundary( sim%jay, sim%g_space, sim%no_co, dt(sim%tstep), sim%bnd )

  ! report time centered energy
  call sim % emf % report_energy( sim%no_co, sim%tstep, t(sim%time))
  call sim % part % report_energy( sim%no_co, sim%tstep, t(sim%time))

  ! --------------------------------------------------
  ! begin injection of particles (cathode, ionization)
  ! --------------------------------------------------

  ! injection of particles must be done outside the
  ! move_boundary(part) <...> update_boundary(part) section
  ! of the code or else some paralell partitions will fail on
  ! moving window runs

  ! note that particles injected in this section will only
  ! be pushed by the EMF in the next iteration

  ! also note that particle injection here should contribute
  ! (if necessary) to the current (e.g. cathode)

  DEBUG('Before cathode')
  call cathode(sim%part, sim%jay%pf(1), t(sim%time), dt(sim%tstep), &
               sim%no_co, sim%grid%coordinates)

  DEBUG('Before ionization')
  ! --------------------------------------------------
  ! end injection of particles
  call ionize_neutral(sim%part, sim%emf, dt(sim%tstep), sim%grid%coordinates)

  ! --------------------------------------------------

  ! --------------------------------------------------
  ! begin particle sorting / collisions
  ! --------------------------------------------------

  ! sorting and collisions are done here, outside the
  ! move_boundary() <...> update_boundary() section
  ! of the loop to avoid some issues related to the moving
  ! window algorithm.
  ! Outside this section all the data (grids/particles)
  ! correctly placed in the local space
  ! sorting and collisions are done together since the
  ! collision algorithm requires sorting the particles

  call sort_collide( sim%part, n(sim%tstep), t(sim%time), n_threads(sim%no_co) )
  ! --------------------------------------------------
  ! end particle sorting / collisions
  ! --------------------------------------------------

  ! normalize the current using the volume elements for cylindrical coordinates
  ! this also takes care of axial boundary
  if ( sim%grid%coordinates == p_cylindrical_b ) then
    DEBUG('Before normalize sim%jay')
    call sim%jay%normalize_cyl( )
  endif

  ! take care of boundaries for currents
  DEBUG('Before update_boundary sim%jay')
  call  sim%jay%update_boundary( sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )

  ! smooth currents
  ! - smoothing has to be after update boundary
  DEBUG('Before smooth sim%jay')
  call sim % jay % smooth()

  ! calculate new electro-magnetic fields
  DEBUG('Before advance sim%emf')
  call sim % emf % advance( sim%jay, dt(sim%tstep), sim%bnd%send_vdf, sim%bnd%recv_vdf )

  ! set the final time for the loop
  call end_event(loopev)

end subroutine iter_sim


!-----------------------------------------------------------------------------------------
! Get the minimum number of guard cells in all directions / boundaries
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_sim( sim )

  use m_particles
  use m_emf
  use m_random

  implicit none

  class( t_simulation ), intent(in) :: sim

  call list_algorithm( sim%options )
  call list_algorithm_rnd()
  call sim % emf % list_algorithm( )
  call sim % part % list_algorithm( )
  call sim%jay%list_algorithm( )

end subroutine list_algorithm_sim


!---------------------------------------------------------------------------------------
! Write checkpoint data
!---------------------------------------------------------------------------------------
subroutine write_checkpoint_sim( sim, restart_handle )

  use m_random
  use m_particles

  implicit none

  class ( t_simulation ), intent(in) ::  sim
  type ( t_restart_handle ), intent(inout)    ::  restart_handle

  call sim % no_co % write_checkpoint( restart_handle )
  call sim % grid % write_checkpoint( restart_handle )
  call write_checkpoint_rnd( rng, restart_handle )
  call restart_write( sim%tstep   , restart_handle )
  call restart_write( sim%restart , restart_handle )
  call restart_write( sim%g_space , restart_handle )
  call restart_write( sim%time    , restart_handle )

  call sim % emf % write_checkpoint( restart_handle )
  call sim % jay % write_checkpoint( restart_handle )
  call sim % part % write_checkpoint( restart_handle )

  call sim % antenna_array % write_checkpoint( restart_handle )

end subroutine write_checkpoint_sim

!-----------------------------------------------------------------------------------------
! Check the input information for errors
!-----------------------------------------------------------------------------------------
subroutine test_input_sim( sim )

  use m_emf
  use m_grid

  implicit none

  class( t_simulation ), intent(inout) :: sim

  ! test the courant condition for EM field solver
  ! (this is more stringent than the stability limit for particles so we only need to test this
  call sim % emf % test_stability( dt(sim%tstep), (xmax(sim%g_space) - xmin(sim%g_space))/sim%grid%g_nx(1:p_x_dim))

  ! test the partition
  call sim%grid%test_partition( sim%no_co, sim% min_gc_sim() )
end subroutine test_input_sim

!-----------------------------------------------------------------------------------------
! Sets the integral of load on all directions for node partition / dynamic load balance
! - The result is saved on a grid object, which may not be the main grid object
!-----------------------------------------------------------------------------------------
subroutine set_int_load_sim( sim, grid, n_ )

  use m_grid_parallel
  use m_particles
  use m_logprof

  implicit none

  class( t_simulation ), intent(inout) :: sim
  class( t_grid ), intent(inout) :: grid
  integer, intent(in) :: n_

  call begin_event(dynlb_set_int_load_ev)

  ! initialize int_load array
  call init_int_load( grid )

  ! add particle load if required
  ! density data will be used for iteration 0
  if ( use_part_load( grid ) ) then
    call add_load( sim % part, grid, n_ )
  endif

  ! add grid / expression load if required
  call add_load( grid, n_ )

  ! gather results from all nodes
  call grid % gather_load( sim % no_co )

  ! convert to integral
  call conv_load_2_int_load( grid )

  call end_event(dynlb_set_int_load_ev)

end subroutine set_int_load_sim
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! Get the minimum number of guard cells in all directions / boundaries
!-----------------------------------------------------------------------------------------
function min_gc_sim( sim )

  use m_particles
  use m_emf

  implicit none

  class ( t_simulation ), intent(in) ::  sim

  integer, dimension(2,p_x_dim) :: min_gc_sim

  integer :: i
  integer, dimension(2, p_x_dim) :: gc

  ! get the required number of guard cells required by the
  ! deposition algorithm chosen
  min_gc_sim = get_min_gc( sim%part )

  ! get the required number of guard cells required by the
  ! field solver / smoothing algorithm chosen
  ! gc = get_min_gc( sim%emf )
  gc = sim%emf%min_gc()
  do i = 1, p_x_dim
    if ( gc( p_lower, i ) > min_gc_sim( p_lower, i ) ) min_gc_sim( p_lower, i ) = gc( p_lower, i )
    if ( gc( p_upper, i ) > min_gc_sim( p_upper, i ) ) min_gc_sim( p_upper, i ) = gc( p_upper, i )
  enddo

  ! get the required number of guard cells required by the
  ! current smoothing algorithm chosen
  gc = sim%jay%get_min_gc()
  do i = 1, p_x_dim
    if ( gc( p_lower, i ) > min_gc_sim( p_lower, i ) ) min_gc_sim( p_lower, i ) = gc( p_lower, i )
    if ( gc( p_upper, i ) > min_gc_sim( p_upper, i ) ) min_gc_sim( p_upper, i ) = gc( p_upper, i )
  enddo

  ! account for moving window
  do i = 1, p_x_dim
    if ( if_move( sim%g_space, i) ) then
       min_gc_sim(p_lower, i) = min_gc_sim(p_lower, i) + 1
    endif
  enddo

end function min_gc_sim


end module m_simulation
