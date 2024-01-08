!-----------------------------------------------------------------------------------------
! Particles definition module
!
! This file contains the class definition for the following classes:
!
!  t_particles
!  t_particles_charge
!-----------------------------------------------------------------------------------------

#include "os-preprocess.fpp"
#include "os-config.h"

module m_particles_define

#include "memory/memory.h"

use m_system
use m_parameters

use m_species_define
use m_cathode

use m_diagnostic_utilities, only: p_diag_prec
use m_bnd

#ifdef __HAS_COLLISIONS__
  use m_species_collisions
#endif

#ifdef __HAS_IONIZATION__
  use m_neutral, only : t_neutral
#endif

use m_emf_define,  only : t_emf
use m_vdf_define,  only : t_vdf, t_vdf_report
use m_time_step,   only : t_time_step
use m_input_file, only : t_input_file
use m_node_conf,   only : t_node_conf
use m_grid_define, only : t_grid
use m_restart,     only : t_restart_handle
use m_space,       only : t_space
use m_current_define, only : t_current
use m_vdf_comm,    only : t_vdf_msg

implicit none

private

! string to id restart data
character(len=*), parameter :: p_part_rst_id = "part rst data - 0x0005"

integer, parameter :: p_charge = 1, &
                      p_charge_htc = 2, &
                      p_dcharge_dt = 3

character(len=*), dimension(3), parameter :: &
   p_report_quants = (/ 'charge    ', &
                        'charge_htc', &
                        'dcharge_dt' /)

type :: t_particles_charge

  ! global charge vdfs
  type( t_vdf ) :: charge0, charge1

  type( t_vdf ), pointer :: current => null(), previous => null()

  integer :: n_last_update = - 1 ! Iteration of the last deposit
  integer :: n_prev_update = - 1 ! Iteration of the previous deposit

  logical :: low_roundoff  = .false.
  logical :: keep_previous = .false.

end type t_particles_charge

!-----------------------------------------------------------------------------------------
! t_particles
!  - Particle class definition - holds all simulation particles
!-----------------------------------------------------------------------------------------

type :: t_particles

  ! the following refer to simple charged particles
  ! that are initialized from density profiles at
  ! setup time

  ! number of species
  integer :: num_species = 0

  ! current iteration - this is updated after advance deposit
  integer :: n_current = 0

  ! interpolation type for all species
  integer :: interpolation

  integer :: pos_type

  ! Field values are centered at the corner of the cell (full momentum conserving)
  logical :: grid_center = .false.

  ! linked list of species
  class( t_species ), pointer :: species  => NULL()

  ! number of cathode objects
  integer :: num_cathode = 0

  ! pointer to array of cathode objects
  type( t_cathode ), dimension(:), pointer :: cathode  => NULL()

#ifdef __HAS_IONIZATION__
  ! number of neutral objects
  integer  :: num_neutral = 0

  ! pointer to array of neutral objects
  class( t_neutral ), dimension(:), pointer :: neutral  => NULL()
#endif

  ! low roundoff multi-species current deposition
  logical :: low_jay_roundoff = .false.
  type( t_vdf ) :: jay_tmp

  ! global charge
  type( t_particles_charge ) :: charge

  type(t_vdf_report), pointer :: reports

  ! precision of diagnostic
  integer :: prec = p_diag_prec

  ! Total kinetic energy diagnostic
  integer :: ndump_fac_ene = 0

#ifdef __HAS_COLLISIONS__
  ! collisions
  type( t_collisions ) :: coll
#endif

contains

  procedure :: allocate_objs   => allocate_objs_particles

  ! these type bound procedures are in separate files,
  ! only the interfaces are declared here
  procedure :: read_input      => read_input_particles
  procedure :: init            => init_particles
  procedure :: update_boundary => update_boundary_particles
  procedure :: report          => report_particles
  procedure :: report_energy   => report_energy_particles
  procedure :: get_diag_buffer_size => get_diag_buffer_size_part
  procedure :: advance_deposit => advance_deposit_particles
  procedure :: cleanup         => cleanup_particles
  procedure :: write_checkpoint => write_checkpoint_particles
  procedure :: restart_read     => restart_read_particles

  procedure :: validate_names => validate_names_part
  procedure :: list_algorithm => list_algorithm_part
  procedure :: move_window => move_window_part
  procedure :: get_spec_by_idx => get_spec_by_idx

end type t_particles

! ----------------------------------------------  type bound procedure (method) interfaces

interface
  subroutine read_input_particles( this, input_file, periodic, if_move, grid, dt, &
                               sim_options )

  import t_particles, t_input_file, p_double, t_options, t_grid

  class( t_particles ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ),               intent(in) :: grid
  real(p_double), intent(in) :: dt
  type( t_options ), intent(in) :: sim_options

  end subroutine
end interface

interface
  subroutine init_particles( this, g_space, jay, emf, &
                            grid, no_co, bnd, ndump_fac, &
                            restart, restart_handle, t, tstep, tmax, sim_options )

  import t_particles, t_space, t_emf, t_vdf, t_grid, t_node_conf, t_restart_handle, &
         p_double, t_options, t_current, t_bnd, t_time_step

  class( t_particles ), intent(inout) :: this
  type( t_space ),     intent(in) :: g_space
  class( t_emf) ,intent(inout) :: emf
  class( t_current ) ,intent(inout) :: jay
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  class(t_bnd), intent(inout) :: bnd
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  real(p_double), intent(in) :: t
  type(t_time_step), intent(in) :: tstep
  real(p_double), intent(in) :: tmax
  type( t_options ), intent(in) :: sim_options

  end subroutine
end interface

interface
  subroutine write_checkpoint_particles( this, restart_handle )

  import t_particles, t_restart_handle

  class( t_particles ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  end subroutine
end interface

interface
  subroutine restart_read_particles( this, restart_handle )

  import t_particles, t_restart_handle

  class( t_particles ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  end subroutine
end interface

interface
  subroutine cleanup_particles( this )

  import t_particles

  class( t_particles ), intent(inout) :: this

  end subroutine
end interface

interface
  subroutine update_boundary_particles( this, jay, g_space, no_co, dt, bnd )

  import t_particles, t_vdf, t_space, t_node_conf, p_double, t_current, t_bnd

  class( t_particles ), intent(inout) :: this
  class( t_current ),       intent(inout) :: jay

  type( t_space ),     intent(in) :: g_space
  class( t_node_conf ), intent(in) :: no_co
  real(p_double),     intent(in) :: dt
  class( t_bnd ), intent(inout) :: bnd

  end subroutine

end interface

interface
  subroutine advance_deposit_particles( this, emf, jay, tstep, t, no_co, options )

  import t_particles, t_emf, t_current, p_double, t_time_step, t_node_conf, t_options

  class( t_particles ), intent(inout) :: this
  class( t_emf ), intent( inout )  ::  emf
  class( t_current ), intent(inout) :: jay

  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t
  class( t_node_conf ), intent(in) :: no_co
  type( t_options ), intent(in) :: options

  end subroutine

end interface

interface
  subroutine report_particles( this, emf, g_space, grid, no_co, tstep, t, send_msg, recv_msg )

    import t_particles, t_emf, t_space, t_grid, t_node_conf, t_time_step, p_double, &
           t_vdf_msg

    class( t_particles ), intent(inout) :: this
    class( t_emf ),       intent(inout) :: emf

    type( t_space ),        intent(in) :: g_space
    class( t_grid ),         intent(in) :: grid
    class( t_node_conf ),    intent(in) :: no_co
    type( t_time_step ),    intent(in) :: tstep
    real(p_double),         intent(in) :: t
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  end subroutine
end interface

interface
  subroutine report_energy_particles( this, no_co, tstep, t )

    import t_particles, t_node_conf, t_time_step, p_double

    class( t_particles ), intent(in) :: this
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t

  end subroutine
end interface

interface
  subroutine get_diag_buffer_size_part( this, gnx, diag_buffer_size )

  import t_particles

  class( t_particles ), intent(in) :: this
  integer, dimension(:), intent(in) :: gnx
  integer, intent(inout) :: diag_buffer_size

  end subroutine
end interface

!-----------------------------------------------------------------------------------------

public :: t_particles_charge, t_particles

public :: p_part_rst_id, p_charge, p_charge_htc, p_dcharge_dt, p_report_quants

integer, public :: pushev = 0, partboundev, cathodeev, neutralev
integer, public :: diag_part_ev, reduce_current_ev

contains

!-----------------------------------------------------------------------------------------
! Allocate particle species, cathodes, etc.
!  - This is now done outside of read_input to allow overriding by subclasses
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_particles( this )

  use m_species_memory
  use m_cathode
  use m_neutral

  implicit none

  class( t_particles ), intent(inout) :: this

  integer :: i
  class(t_species), pointer :: spec, tail


  if ( this%num_species > 0 ) then

    allocate(spec)
    this % species => spec
    do i = 2, this%num_species
      tail => spec
      allocate(spec)
      tail % next => spec
    enddo

    if ( this%num_cathode > 0 ) then
      allocate( this%cathode( this%num_cathode ) )
    endif

#ifdef __HAS_IONIZATION__
    if ( this%num_neutral > 0 ) then
      allocate( this%neutral( this%num_neutral ) )
    endif
#endif

  endif

end subroutine allocate_objs_particles
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Printout the algorithm used by each species pusher
!-------------------------------------------------------------------------------
subroutine list_algorithm_part( this )

#ifdef __HAS_IONIZATION__
  use m_neutral, only : list_algorithm
#endif

  implicit none
  class( t_particles ), intent(in) :: this
  class( t_species ), pointer :: species
  integer :: i

  print *, ' '
  print '(A)', 'Particles:'

  print '(A,I0)', '- Number of species: ', this%num_species

  if ( this%num_species > 0 ) then

  ! particle shape
  select case ( this%interpolation )
    case (p_linear)
    print '(A)', '- Linear (1st order) interpolation'
    case (p_quadratic)
    print '(A)', '- Quadratic (2nd order) interpolation'
    case (p_cubic)
    print '(A)', '- Cubic (3rd order) interpolation'
    case (p_quartic)
    print '(A)', '- Quartic (4th order) interpolation'
  end select

  endif

  if ( this%grid_center ) then
    print '(A)', '- Using spatially centered grid values for interpolation'
  else
    print '(A)', '- Using staggered Yee grid values for interpolation'
  endif

  ! Current deposition
  if ( this%low_jay_roundoff ) then
    print *, '- Depositing current using low roundoff algorithm'
  endif

  species => this % species
  do
    if (.not. associated(species)) exit
    call species % list_algorithm()
    species => species % next
  enddo

#ifdef __HAS_COLLISIONS__

  ! report on the collision algorithm if required
  if ( this%coll%n_collide > 0 ) then
     call list_algorithm_collisions( this%coll )
  endif

#endif

#ifdef __HAS_IONIZATION__
  do i = 1, this%num_neutral
    call list_algorithm(this%neutral(i))
  enddo
#endif

end subroutine list_algorithm_part
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! moves the boundaries for particle data
!-----------------------------------------------------------------------------------------
subroutine move_window_part( this, g_space, grid, jay, no_co, bnd )

  use m_species
  use m_neutral
  use m_node_conf

  implicit none

  class( t_particles ), intent(inout) :: this

  type( t_space ),     intent(in) :: g_space
  class( t_grid ),  intent(in) :: grid
  class(t_current), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  class(t_bnd), intent(inout) :: bnd

  class( t_species ), pointer :: species
  integer :: i

  species => this % species
  do
    if (.not. associated(species)) exit
    call  move_window( species, g_space, grid, jay, no_co, bnd%bnd_cross, bnd%node_cross, bnd%send_spec, bnd%recv_spec )
    species => species % next
  enddo

  do i = 1, this%num_cathode
   call  this % cathode( i ) % move_window( g_space )
  enddo

#ifdef __HAS_IONIZATION__

  do i = 1, this%num_neutral
   call this%neutral(i)%move_window(grid%my_nx(p_lower, : ), g_space , &
                                    no_co%on_edge( 1, p_upper ))
  enddo

#endif

end subroutine move_window_part
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Check that all species have unique names
!-----------------------------------------------------------------------------------------
subroutine validate_names_part( this )

  use m_species_define
  use stringutil

  implicit none

  class( t_particles ), intent(inout) :: this

  class( t_species ), pointer :: sp1, sp2
  integer :: i1, i2

  sp1 => this % species; i1 = 1
  do
    if (.not. associated(sp1)) exit
    sp2 => sp1 % next; i2 = i1 + 1
    do
      if (.not. associated(sp2)) exit
      if (lowercase(replace_blanks(sp1 % name)) == &
          lowercase(replace_blanks(sp2 % name))) then
        if (mpi_node() == 0) then
          print '(A)', 'Error reading species parameters, invalid species name.'
          print '(A,I0,A,A,A)', 'Name for species (',i1,') - ', sp1 %name, '"'
          print '(A,I0,A,A,A)', 'Conflicts with name for species (',i2,') - ', sp2 %name, '"'
          print '(A)', 'aborting...'
        endif
      endif
      sp2 => sp2 % next; i2 = i2 + 1
    enddo
    sp1 => sp1 % next; i1 = i1 + 1
  enddo


end subroutine validate_names_part
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Utitlity function to return the species object at idx in the species list.
!   TODO: decide if there should be error checking a a thrown runtime error if the idx is invalid
!-----------------------------------------------------------------------------------------
function get_spec_by_idx( this, idx )
  implicit none
  class( t_species ), pointer      :: get_spec_by_idx
  class( t_particles ), intent(in) :: this
  integer                          :: idx

  integer                          :: temp
  class( t_species ),pointer       :: species

  temp = idx
  species => this % species
  do
    if (.not. associated(species)) exit
    temp = temp - 1
    if (temp .le. 0) exit
    species => species % next
  enddo

  if (temp .gt. 0) then
    get_spec_by_idx => null()
  else
    get_spec_by_idx => species
  endif

end function get_spec_by_idx

end module m_particles_define
