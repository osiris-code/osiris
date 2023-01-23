!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     particle class 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#define __DEBUG__ 1

#include "os-config.h"
#include "os-preprocess.fpp"

module m_particles

#include "memory/memory.h"

use m_system
use m_parameters

use m_particles_define
use m_particles_charge

use m_emf_define

use m_vdf_define
use m_vdf_average
use m_vdf_report
use m_vdf_math
use m_vdf

use m_species_define
use m_species_memory
use m_species_diagnostics
use m_species_loadbalance
use m_species_charge
use m_species

#ifdef __HAS_COLLISIONS__
use m_species_collisions
#endif

#ifdef __HAS_IONIZATION__
use m_neutral
#endif

use m_space
use m_grid_define
use m_node_conf
use m_restart

use stringutil

use m_logprof
use m_diagnostic_utilities
use m_diagfile, only: p_time_length

use m_math

implicit none

private

!  interface restart_read
!  module procedure restart_read_particles
!  end interface

!  interface num_species
!  module procedure num_species_particles
!  end interface

interface cathode
module procedure inject_cathode_part
end interface

interface ionize_neutral
module procedure ionize_neutral_part
end interface

interface validate
module procedure validate_part
end interface

interface sort_collide
module procedure sort_collide_part
end interface

interface reshape_obj
module procedure reshape_part
end interface

interface get_min_gc
module procedure get_min_gc_part
end interface

interface get_grid_center
module procedure get_grid_center_part
end interface

interface report_global_load
module procedure report_global_load_part
end interface

interface report_node_load
module procedure report_node_load_part
end interface

interface report_grid_load
module procedure report_grid_load_part
end interface

interface add_load
module procedure add_load_particles
end interface

interface num_par
module procedure num_par
end interface

interface interpolation
module procedure interpolation_par
end interface

! declare things that should be public
public :: setup, cleanup
public :: restart_read, restart_write
public :: get_min_gc
public :: validate, sort_collide

public :: report_global_load, report_node_load, report_grid_load

public :: cathode, get_grid_center
public :: ionize_neutral

public :: add_load, num_par, reshape_obj

public ::  interpolation

contains

!-----------------------------------------------------------------------------------------
! result is the number of species
!-----------------------------------------------------------------------------------------
function num_species_particles(this)

    implicit none

    integer :: num_species_particles

    class( t_particles ), intent( in )  ::  this

    num_species_particles = this%num_species

end function num_species_particles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Inject particles from cathodes
!-----------------------------------------------------------------------------------------
subroutine inject_cathode_part( this, jay, t, dt, no_co, coordinates )

    implicit none

    class( t_particles ), intent(inout) :: this
    type( t_vdf ),       intent(inout) :: jay
    real( p_double ),   intent(in)    :: t
    real( p_double ),   intent(in)    :: dt
    integer, intent(in) :: coordinates
    class( t_node_conf ), intent(in)  ::  no_co

    ! local variables

    integer :: i

    ! executable statements

    call begin_event(cathodeev)

    if ( this%low_jay_roundoff ) then

    ! this version has a lower roundoff error because the current from each cathode
    ! is deposited onto a clean grid, and the grids are then added
    do i=1, this%num_cathode
        call this % cathode(i) % inject_lowroundoff( jay, this%jay_tmp, &
        t, dt, no_co, coordinates)
    enddo

    else

    do i=1, this%num_cathode
        call this % cathode(i) % inject( jay, t, dt, no_co, coordinates)
    enddo

    endif

    call end_event(cathodeev)

end subroutine inject_cathode_part
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine ionize_neutral_part(this, emf, dt, coordinates)
!-----------------------------------------------------------------------------------------
! take care of neutrals (ionization)
!-----------------------------------------------------------------------------------------

    implicit none

    ! dummy variables

    class( t_particles ), intent(inout) :: this
    class( t_emf ), intent(in) :: emf
    real( p_double ), intent(in) :: dt
    integer, intent(in) :: coordinates

    ! local variables
    integer :: n
    ! executable statements

    call begin_event(neutralev)

#ifdef __HAS_IONIZATION__

    if (this%num_neutral > 0) then
    do n = 1, this%num_neutral
        call this% neutral(n) % ionize(this%species, emf, dt, coordinates)
    end do
    endif

#endif

    call end_event(neutralev)

end subroutine ionize_neutral_part
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! sorts/collides all particles
!-----------------------------------------------------------------------------------------
subroutine sort_collide_part( this, n, t, n_threads )
!-----------------------------------------------------------------------------------------

    implicit none

    class( t_particles ), intent(inout) :: this
    integer, intent(in) :: n
    real(p_double), intent(in) :: t
    integer, intent(in) :: n_threads

    class( t_species ), pointer :: species

    ! executable statements

#ifdef __HAS_COLLISIONS__

    if ( if_collide( this%coll, n ) ) then
    ! collide cycle: collide all species - will be sorted there
    call collide_particles( this%coll, this%species, this%num_species, t, &
    this%species%dt, n, n_threads)
    else
    ! we are not in a collide cycle: just sort all particles
    species => this % species
    do
        if (.not. associated(species)) exit
        call sort( species, n, t )
        species => species % next
    enddo
    endif

#else

    species => this % species
    do
        if (.not. associated(species)) exit
        call sort( species, n, t )
        species => species % next
    enddo

#endif

end subroutine sort_collide_part
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine validate_part( this, msg, over  )
!-----------------------------------------------------------------------------------------
! checks if all the particles have valid values:
!  - Check for NaNs and Infs on positions and momenta
!  - Check particle positions against node boundaries
!  - Check for 0, NaNs and Infs on charge
! Use for debugging purposes only.
!-----------------------------------------------------------------------------------------

    implicit none

    class( t_particles ), intent(in) :: this
    character( len = * ), intent(in) :: msg
    logical, intent(in), optional :: over

    class( t_species ), pointer :: species

    if ( present(over) ) then
        species => this % species
        do
            if (.not. associated(species)) exit
            call species%validate( msg, over )
            species => species % next
        enddo
    else
        species => this % species
        do
            if (.not. associated(species)) exit
            call species%validate( msg )
            species => species % next
        enddo
    endif

end subroutine validate_part
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Return the minimum number of guard cells required by the particle deposition
! scheme
!-------------------------------------------------------------------------------
function get_min_gc_part( this )

    implicit none

    class( t_particles ), intent(in) :: this
    integer, dimension(2, p_x_dim) :: get_min_gc_part
    class( t_species ), pointer :: species

    select case ( this%interpolation )
    case( p_linear )
        get_min_gc_part(p_lower, :) = 1
        get_min_gc_part(p_upper, :) = 2
        
    case( p_quadratic )
        get_min_gc_part(p_lower, :) = 2
        get_min_gc_part(p_upper, :) = 3
        
    case( p_cubic )
        get_min_gc_part(p_lower, :) = 3
        get_min_gc_part(p_upper, :) = 4
        
    case( p_quartic )
        get_min_gc_part(p_lower, :) = 4
        get_min_gc_part(p_upper, :) = 5
        
    end select

    ! If any of the species is using the radiation cooling pusher 1 extra gc is required
    species => this % species
    do
        if (.not. associated(species)) exit
        if ( species%push_type == p_radcool ) then
            get_min_gc_part(p_lower, :) = get_min_gc_part(p_lower, :) + 1
            get_min_gc_part(p_upper, :) = get_min_gc_part(p_upper, :) + 1
            exit
        endif
        species => species % next
    enddo

end function get_min_gc_part
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Return true if the field values are expected to be centered on the cell corner
!-----------------------------------------------------------------------------------------
function get_grid_center_part( this )

    implicit none

    class( t_particles ), intent(in) :: this
    logical :: get_grid_center_part

    get_grid_center_part = this%grid_center

end function get_grid_center_part
!-----------------------------------------------------------------------------------------

!*********************** Load Balance

!---------------------------------------------------------------------------------------------------
! Reports global particle load information (total, min, max and avg. number of particles
! per node)
!---------------------------------------------------------------------------------------------------
subroutine report_global_load_part( this, n, no_co )

    !use mpi

    implicit none

    class( t_particles ), intent(in) :: this
    integer, intent(in)             :: n
    class( t_node_conf ), intent(in) :: no_co

    integer( p_int64 ), dimension(2) :: npart, tmp
    integer( p_int64 ), dimension(1) :: parts

    character(80)   :: path, full_name
    class( t_species ), pointer :: species
    integer :: ierr

    parts(1) = 0

    ! get total particles on local node
    species => this % species
    do
        if (.not. associated(species)) exit
        parts(1) = parts(1) + species%num_par
        species => species % next
    enddo

    if ( no_num(no_co) > 1 ) then
        
        ! get maximum and minimum number of particles on all nodes
        npart(1) =  parts(1)
        npart(2) = -parts(1)
        
        call MPI_REDUCE( npart, tmp, 2, MPI_INTEGER8, MPI_MAX, 0, comm(no_co), ierr )
        
        npart(1) =  tmp(1)
        npart(2) = -tmp(2)
        
        ! get total number of particles (this needs to be done in 64 bits because the total number
        ! of particles can be larger than the 32 bit signed limit, 2^32-1 = 2.147e9)
        call MPI_REDUCE( parts, tmp, 1, MPI_INTEGER8, MPI_SUM, 0, comm(no_co), ierr )
        parts(1) = tmp(1)
        
    else
        
        npart(1) = parts(1)
        npart(2) = parts(1)
        
    endif

    if ( root( no_co ) ) then
        
        path = trim(path_hist) // 'LOAD' // p_dir_sep // 'GLOBAL'
        full_name = trim(path) // p_dir_sep // 'global_particle_load'
        
        if ( n == 0 ) then
            ! create directory
            call mkdir( path, ierr )
            ! create file
            open (unit=file_id_part_load, file=full_name, status = 'replace' , &
            form='formatted')
            ! print header
            write( file_id_part_load, '(A6, 4(1X,A18))' ) &
            'Iter', 'Total Parts.', 'Min', 'Max', 'Avg'
            write( file_id_part_load, '(A)' ) &
            '----------------------------------------------------------------------------------'
        else
            open (unit=file_id_part_load, file=full_name, position = 'append', &
            form='formatted')
        endif
        
        ! write particle load information to disk
        write( file_id_part_load, '(I6, 4(1X,I18))' ) &
        n, parts(1), npart(2), npart(1), parts(1) / no_num( no_co )
        
        close(file_id_part_load)
    endif

end subroutine report_global_load_part
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Reports number of particles per node for all nodes (small grid file)
!---------------------------------------------------------------------------------------------------
subroutine report_node_load_part( this, n, ndump, t, no_co )

    implicit none

    class( t_particles ), intent(in) :: this
    integer, intent(in)             :: n, ndump
    real( p_double ), intent(in)    :: t
    class( t_node_conf ), intent(in) :: no_co

    real( p_double ), dimension(:), pointer     :: node_load
    real( p_double ), dimension(:), pointer     :: f1
    real( p_double ), dimension(:,:), pointer   :: f2
    real( p_double ), dimension(:,:,:), pointer :: f3

    class( t_diag_file ), allocatable :: diagFile

    integer, dimension(p_max_dim) :: lnx

    class( t_species ), pointer :: species
    integer :: i, npart, nnodes, ierr
    real( p_double ) :: npart_dbl

    node_load => null()
    f1 => null()
    f2 => null()
    f3 => null()

    ! get total number of particles on node
    npart = 0
    species => this % species
    do
        if (.not. associated(species)) exit
        npart = npart + species%num_par
        species => species % next
    enddo
    npart_dbl = npart

    ! gather data
    nnodes =  no_num( no_co )
    call alloc( node_load, (/nnodes/) )

    if ( nnodes > 1 ) then
        call mpi_gather( npart_dbl, 1, MPI_DOUBLE_PRECISION, &
        node_load, 1, MPI_DOUBLE_PRECISION, &
        0, comm(no_co), ierr )
        if ( ierr /= MPI_SUCCESS ) then
            ERROR('MPI Error')
            call abort_program( p_err_mpi )
        endif
    else
        node_load = npart_dbl
    endif

    ! save data to disk
    if ( root( no_co ) ) then
        
        lnx(1:p_x_dim) = nx( no_co )
        
        ! Open the output file
        call create_diag_file( diagFile )
        diagFile%ftype   = p_diag_grid
        
        diagFile%filepath  = trim(path_mass) // 'LOAD' // p_dir_sep // 'NODE' // p_dir_sep
        diagFile%filename  = trim(get_filename(n/ndump, 'node_load'))
        diagFile%name      = 'Particles per node'
        diagFile%iter%n         = n
        diagFile%iter%t         = t
        diagFile%iter%time_units = '1 / \omega_p'
        
        
        diagFile%grid % ndims = p_x_dim
        diagFile%grid % name  = "load"
        diagFile%grid % label = "particles per node"
        diagFile%grid % units = "particles"
        
        do i = 1, p_x_dim
            diagFile % grid % axis(i) % min = 0
            diagFile % grid% axis(i) % max = lnx(i)
            write( diagFile % grid% axis(i) % name, '(A,I0)' ) 'x', i
            write( diagFile % grid% axis(i) % label, '(A,I0)' ) 'x_', i
            diagFile % grid% axis(i) % units = 'cell'
        enddo
        
        call diagFile % open( p_diag_create )
        
        ! Move node load values to the proper position and write node load values
        select case (p_x_dim)
        case(1)
            call alloc( f1, lnx )
            f1 = -1.0
            do i = 1, nnodes
                f1( no_co%ngp( i, 1 ) ) = node_load(i)
            enddo
            call diagFile % add_dataset( "load", f1 )
            call freemem( f1 )
            
        case(2)
            call alloc( f2, lnx )
            f2 = -1.0
            do i = 1, nnodes
                f2( no_co%ngp( i, 1 ), no_co%ngp( i, 2 ) ) = node_load(i)
            enddo
            call diagFile % add_dataset( "load", f2 )
            call freemem( f2 )
            
        case(3)
            call alloc( f3, lnx )
            f3 = -1.0
            do i = 1, nnodes
                f3( no_co%ngp( i, 1 ), no_co%ngp( i, 2 ), no_co%ngp( i, 3 ) ) = node_load(i)
            enddo
            call diagFile % add_dataset( "load", f3 )
            call freemem( f3 )
        end select
        
        ! close output file
        call diagFile % close( )
        deallocate( diagFile )
    endif

    ! free load array
    call freemem( node_load )

end subroutine report_node_load_part
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Reports number of particles per node for each grid cell
!---------------------------------------------------------------------------------------------------
subroutine report_grid_load_part( this, n, ndump, t, g_space, grid, no_co )

    implicit none

    class( t_particles ), intent(in) :: this
    integer, intent(in)             :: n, ndump
    real( p_double ), intent(in)    :: t
    type( t_space ), intent(in)     :: g_space
    class( t_grid ), intent(in)      :: grid
    class( t_node_conf ), intent(in) :: no_co

    ! local variables
    type( t_vdf ) :: load
    type( t_vdf_report ) :: load_report
    class( t_species ), pointer :: species

    integer, dimension(2,3) :: gc_num
    real( p_double ), dimension(3) :: dx

    ! for the old version of p_lower_cell positions we may have particles at nx+1 if a physical
    ! boundary exists at that edge. This is just a simple hack to aacount for that situation
    gc_num( p_lower, : ) = 0
    gc_num( p_upper, : ) = 1

    ! create vdf to hold the load
    dx = 1.0
    call load % new( p_x_dim, 1,  grid%my_nx( 3, : ), gc_num, dx, .true. )

    ! get the load on each cell
    species => this % species
    do
        if (.not. associated(species)) exit
        call deposit_cell_load( species, load )
        species => species % next
    enddo

    ! write the data
    load_report%name = 'cell_load'
    load_report%ndump = ndump

    load_report%xname  = (/'x1', 'x2', 'x3'/)
    load_report%xlabel = (/'x_1', 'x_2', 'x_3'/)
    load_report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

    load_report%time_units = '1 / \omega_p'
    load_report%dt         = 1.0

    load_report%fileLabel = ''
    load_report%path  = trim(path_mass) // 'LOAD' // p_dir_sep // 'CELL'

    load_report%label = 'Particles per cell'
    load_report%units = 'particles'

    load_report%n = n
    load_report%t = t

    load_report%path = trim(path_mass) // 'LOAD' // p_dir_sep // 'GRID'  // p_dir_sep
    load_report%filename = 'cell_load-'// idx_string( n/ndump, p_time_length )

    call load % write( load_report, 1, g_space, grid, no_co )

    ! cleanup the vdf
    call load % cleanup()

end subroutine report_grid_load_part
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Sets the particle load.
! For iteration n = 0 the load is set for the whole simulation volume (since this is called
! when no node partitions exist yet), for n > 0 the load is set for the
! local volume/particles only
!---------------------------------------------------------------------------------------------------
subroutine add_load_particles( this, grid, n )

    use m_grid_parallel

    implicit none

    class( t_particles ), intent(inout) :: this
    class( t_grid ), intent(inout) :: grid
    integer, intent(in) :: n


    class( t_species ), pointer :: species

    ! loop through all species and add to int_load array
    species => this % species
    if ( n > 0 ) then
        do
            if (.not. associated(species)) exit
            call add_particle_load( species, grid  )
            species => species % next
        enddo
    else
        do
            if (.not. associated(species)) exit
            call add_density_load( species, grid  )
            species => species % next
        enddo
    endif

end subroutine add_load_particles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reshape_part( this, old_lb, new_lb, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------
! redistribute particles to new simulation partition
!-----------------------------------------------------------------------------------------

    use m_grid_parallel
    use m_vdf_comm

    implicit none

    class( t_particles ), intent(inout) :: this
    class( t_grid ), intent(in) :: new_lb, old_lb
    class( t_node_conf ), intent(in) :: no_co
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

    type( t_msg_patt ) :: msg_patt
    class( t_species ), pointer :: species
    integer :: n
    integer, dimension(2,p_max_dim) :: gc_num
    type(t_vdf_report), pointer :: report

    ! reshape species
    if ( this%num_species > 0 ) then
        ! get message pattern (guard cells are not required since no particles
        ! should be in guard cells). This is the same for all species
        ! so we only need to do it once.
        gc_num = 0
        call new_msg_patt( msg_patt, old_lb, new_lb, no_co, gc_num )
        
        ! loop through all species
        species => this % species
        do
            if (.not. associated(species)) exit
            call species % reshape( old_lb, new_lb, msg_patt, no_co, send_msg, recv_msg )
            species => species % next
        enddo
        
        ! clear message pattern data
        call msg_patt % cleanup()
    endif

    #ifdef __HAS_IONIZATION__

    ! reshape neutral objects ( which are actually grids )
    if ( this%num_neutral > 0 ) then
        do n =1, this%num_neutral
            call this%neutral(n) % reshape_obj(old_lb, new_lb, no_co, send_msg, recv_msg )
        enddo
    endif

    #endif

    ! reshape local vdf objects if needed
    if ( this%jay_tmp%x_dim_ > 0 ) call reshape_nocopy( this%jay_tmp, new_lb )

    call reshape( this%charge, old_lb, new_lb, no_co, send_msg, recv_msg )

    ! Reshape tavg_data
    report => this%reports
    do
        if ( .not. associated(report) ) exit
        if ( report%tavg_data%x_dim_ > 0 ) then
            call reshape_copy( report%tavg_data, old_lb, new_lb, no_co, send_msg, recv_msg )
        endif
        report => report%next
    enddo

end subroutine reshape_part
!-----------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
! Return total number of particles on all species
!--------------------------------------------------------------------------------------------------
function num_par( this )

    implicit none

    integer :: num_par
    class( t_particles ), intent(in) :: this

    class( t_species ), pointer :: species

    num_par = 0
    species => this % species
    do
        if (.not. associated(species)) exit
        num_par = num_par + species%num_par
        species => species % next
    enddo

end function num_par
!--------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Return interpolation level used by the particles
!-----------------------------------------------------------------------------------------
function interpolation_par( this )
!-----------------------------------------------------------------------------------------

    implicit none

    integer :: interpolation_par
    class( t_particles ), intent(in) :: this

    interpolation_par = this%interpolation

end function interpolation_par
!-----------------------------------------------------------------------------------------

end module m_particles

!-----------------------------------------------------------------------------------------
! read necessary information from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_particles( this, input_file, periodic, if_move, grid, dt, &
    sim_options )
!-----------------------------------------------------------------------------------------

    use m_system
    use m_parameters

    use m_particles_define, only : t_particles, p_report_quants
    use m_input_file
    use m_vdf_define,       only : p_max_reports_len, p_max_reports, p_n_report_type, &
        p_full, p_savg, p_senv, p_line, p_slice

    use m_species_define
    use m_species

    use m_diagnostic_utilities,  only : p_diag_prec
    use m_neutral
    use m_vdf_report, only : new

    use m_species_collisions, only : read_nml
    use m_grid_define, only : t_grid

    use stringutil

    implicit none

    ! dummy variables

    class( t_particles ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file
    logical, dimension(:), intent(in) :: periodic, if_move
    class( t_grid ),             intent(in) :: grid
    real(p_double), intent(in) :: dt
    type( t_options ), intent(in) :: sim_options

    integer  :: num_species, i
    class( t_species ), pointer :: species
    integer  :: num_cathode
    integer  :: num_neutral, num_neutral_mov_ions
    logical  :: low_jay_roundoff

    character( len = p_max_reports_len ), dimension( p_max_reports ) :: reports
    integer  :: ndump_fac, ndump_fac_ave, ndump_fac_lineout, ndump_fac_ene
    integer, dimension(p_x_dim) :: n_ave
    integer                     :: prec, n_tavg

    character(len=20) :: interpolation
    logical :: grid_center

    namelist /nl_particles/ num_species, num_cathode, num_neutral, &
    num_neutral_mov_ions, low_jay_roundoff, &
    ndump_fac, ndump_fac_ave, ndump_fac_lineout, n_ave, prec, &
    reports, n_tavg, interpolation, grid_center, ndump_fac_ene

    integer, dimension( p_n_report_type ) :: ndump_fac_all
    integer :: ierr

    character(len = p_max_spname_len) :: spname

    ! executable statements

    num_species = 0
    num_cathode = 0
    num_neutral = 0
    num_neutral_mov_ions = 0

    low_jay_roundoff = .false.

    ! Diagnostics
    reports = "-"
    ndump_fac = 0
    ndump_fac_ave = 0
    ndump_fac_lineout = 0
    ndump_fac_ene = 0

    n_ave = -1
    n_tavg = -1
    prec = p_diag_prec

    interpolation = "quadratic"
    grid_center = .false.

    ! Get namelist text from input file
    call get_namelist( input_file, "nl_particles", ierr )

    if ( ierr /= 0 ) then
        if (ierr < 0) then
            print *, "Error reading particles parameters"
        else
            print *, "Error: particles parameters missing"
        endif
        print *, "aborting..."
        stop
    endif

    read (input_file%nml_text, nml = nl_particles, iostat = ierr)
    if (ierr /= 0) then
        print *, "Error reading particles parameters"
        print *, "aborting..."
        stop
    endif

    ! process interpolation scheme
    select case ( trim( interpolation ))
    ! case ( "ngp" ) ! not implemented
    !  this%interpolation = p_ngp
    case ( "linear" )
        this%interpolation = p_linear
    case ( "quadratic" )
        this%interpolation = p_quadratic
    case ( "cubic" )
        this%interpolation = p_cubic
    case ( "quartic" )
        this%interpolation = p_quartic

    case default
        print *, '   Error reading species parameters'
        print *, '   invalid interpolation: "', trim( interpolation ), '"'
        print *, '   Valid values are "linear", "quadratic", "cubic" and "quartic".'
        print *, '   aborting...'
        stop
    end select

    this%grid_center = grid_center

    ! process position type
    select case ( this%interpolation )
    case ( p_linear, p_cubic )
        this%pos_type = p_cell_low
    case ( p_quadratic, p_quartic )
        this%pos_type = p_cell_near
    end select

    this%low_jay_roundoff = low_jay_roundoff


    ! Particle kinetic energy diagnostics
    this%ndump_fac_ene = ndump_fac_ene

    ! global charge diagnostics
    ndump_fac_all(p_full)  = ndump_fac
    ndump_fac_all(p_savg)  = ndump_fac_ave
    ndump_fac_all(p_senv)  = ndump_fac_ave
    ndump_fac_all(p_line)  = ndump_fac_lineout
    ndump_fac_all(p_slice) = ndump_fac_lineout

    ! process normal reports
    call new( this%reports, reports, p_report_quants, &
    ndump_fac_all, n_ave, n_tavg, prec, &
    p_x_dim, ierr )
    if ( ierr /= 0 ) then
        print *, "(*error*) Invalid report"
        print *, "(*error*) aborting..."
        stop
    endif

    ! Allocate sub-objects
    this%num_species = num_species + num_cathode + num_neutral + 2*num_neutral_mov_ions
    this%num_cathode = num_cathode
    this%num_neutral = num_neutral + num_neutral_mov_ions
    call this % allocate_objs()

    species => this % species

    ! Read the input files for standard species
    if ( num_species>0 ) then
        do i=1, num_species
            
            write(spname, '(A,I0)') 'species ',i
            if ( mpi_node() == 0 .and. disp_out(input_file) ) then
                print '(A,I0,A)', " - species (",i,") configuration..."
            endif
            call species % read_input( input_file, spname, periodic, if_move, grid, &
                dt, .true., sim_options )
            
            species => species % next
        enddo
    endif

    ! Take care of cathodes
    if ( num_cathode > 0 ) then

        do i=1, num_cathode
            write(spname, '(A,I0)') 'cathode ',i
            if ( mpi_node() == 0 .and. disp_out(input_file) ) then
                print '(A,I0,A)', " - cathode (",i,") configuration..."
            endif
            
            ! read the cathode configuration, and the associated species
            ! configuration
            call this % cathode(i) % read_input( input_file, species, spname,  periodic, if_move, &
            grid, dt, sim_options )
            
            species => species % next
            
        end do

    endif

#ifdef __HAS_IONIZATION__

    ! Take care of neutrals

    if ( num_neutral > 0 ) then

        if (sim_options%omega_p0 <= 0.0) then
            print *, "(*error*) When using ionization (num_neutral > 0) omega_p0 or n0 must be"
            print *, "(*error*) set in the simulation section at the beggining of the input file"
            print *, "(*error*) bailing out..."
            call abort_program()
        endif

        do i=1, num_neutral
            write(spname, '(A,I0)') 'neutral ',i
            
            if ( mpi_node() == 0 .and. disp_out(input_file) ) then
                print '(A,I0,A)', " - neutral (",i,") configuration..."
            endif
            
            call this%neutral( i ) % read_input( input_file, species, spname,  .false., &
            periodic, if_move, grid, dt, &
            sim_options )
            
            species => species % next
        end do

    endif

    ! Take care of neutrals with moving ions

    if ( num_neutral_mov_ions > 0 ) then

        if (sim_options%omega_p0 <= 0.0) then
            print *, "(*error*) When using ionization (num_neutral_mov_ions > 0) omega_p0 or n0 must be"
            print *, "(*error*) set in the simulation section at the beggining of the input file"
            print *, "(*error*) bailing out..."
            call abort_program()
        endif

        do i=1, num_neutral_mov_ions
            
            write(spname, '(A,I0)') "neutral_mov_ions ", i
            
            if (disp_out(input_file)) then
                SCR_ROOT(" - neutral with moving ions (",i,") configuration...")
            endif
            
            call this%neutral( num_neutral + i )%read_input( input_file, species, spname,  .true., &
            periodic, if_move, grid, dt, &
            sim_options )
            
            species => species % next
            
        end do

    endif

#else

    if (( num_neutral > 0 ) .or. ( num_neutral_mov_ions > 0 )) then
        print *, "(*error*) Ionization is not supported in this version"
        stop
    endif

#endif

    ! validate species names
    call this % validate_names()

#ifdef __HAS_COLLISIONS__
    ! Read collision data
    call read_nml( this%coll, input_file, this%species, this%num_species, sim_options )
#endif

end subroutine read_input_particles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_particles( this, g_space, jay, emf, &
    grid, no_co, bnd, ndump_fac, &
    restart, restart_handle, t, tstep, tmax, sim_options )
!-----------------------------------------------------------------------------------------
! sets up this data structure from the given information
!-----------------------------------------------------------------------------------------

    use m_system
    use m_parameters

    use m_particles_define
    use m_particles_charge
    use m_species_collisions
    use m_species_define
    use m_species, only : sortev, sort_genidx_ev, sort_rearrange_ev, init_buffers_spec

    use m_vdf_define, only : t_vdf, t_vdf_report
    use m_grid_define, only : t_grid
    use m_node_conf, only : t_node_conf
    use m_restart, only : t_restart_handle
    use m_particles, only : restart_read

    use m_space,       only : t_space
    use m_emf_define,  only : t_emf
    use m_current_define,  only : t_current
    use m_bnd

    use m_logprof
    use m_time_step

    implicit none

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

    type( t_vdf_report ), pointer :: report
    integer, parameter :: izero = iachar('0')
    logical :: keep_previous_charge

    class(t_species), pointer :: species
    integer :: i, cathode_id, neutral_id
    real(p_double) :: interp_offset

    ! check that restart values match input deck
    if ( restart ) then
        call this % restart_read( restart_handle )
    else
        this % n_current = 0
    endif

    select case( this%interpolation )
    case( p_linear, p_cubic )
        interp_offset = 0.0_p_double
    case( p_quadratic, p_quartic )
        interp_offset = 0.5_p_double
    case default
        interp_offset = 0.0_p_double
        ERROR('Interpolation value not supported')
        call abort_program(p_err_invalid)
    end select

    ! setup diagnostics
    keep_previous_charge = .false.

    report => this%reports
    do
        if ( .not. associated( report ) ) exit
        
        report%xname  = (/'x1', 'x2', 'x3'/)
        report%xlabel = (/'x_1', 'x_2', 'x_3'/)
        report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)
        
        ! these are just dummy values for now
        report%time_units = '1 / \omega_p'
        report%dt         = 1.0
        
        report%fileLabel = ''
        report%basePath  = trim(path_mass) // 'FLD' // p_dir_sep
        
        select case ( report%quant )
        case ( p_charge )
            report%label = '\rho'
            if ( p_x_dim == 1 ) then
                report%units  = 'e \omega_p / c'
            else
                report%units  = 'e \omega_p^'//char(izero+p_x_dim)// &
                '/ c^'//char(izero+p_x_dim)
            endif
            report%offset_t = 0.0_p_double
            
        case ( p_charge_htc )
            report%label = '\rho_{n-1/2}'
            if ( p_x_dim == 1 ) then
                report%units  = 'e \omega_p / c'
            else
                report%units  = 'e \omega_p^'//char(izero+p_x_dim)// &
                '/ c^'//char(izero+p_x_dim)
            endif
            report%offset_t = -0.5_p_double
            keep_previous_charge = .true.
            
        case ( p_dcharge_dt )
            report%label = 'd\rho/dt'
            if ( p_x_dim == 1 ) then
                report%units  = 'e \omega_p / c'
            else
                report%units  = 'e \omega_p^'//char(izero+p_x_dim)// &
                '/ c^'//char(izero+p_x_dim)
            endif
            report%offset_t = -0.5_p_double
            keep_previous_charge = .true.
            
        end select

        report%offset_x = (/ 0.0_p_double, 0.0_p_double, 0.0_p_double /) + interp_offset
        
        report => report%next
    enddo

    ! Setup charge
    call setup( this%charge, keep_previous_charge, restart, restart_handle )
    ! setup species
    if ( this%num_species > 0 ) then
        
        ! Initialize buffers for communications and current deposition
        call init_buffers_spec( )
        
        species => this % species; i = 1
        do
            if (.not. associated(species)) exit
            call species%init( i, this%interpolation, this%grid_center, grid, &
            g_space, emf, jay, no_co, bnd%send_vdf, bnd%recv_vdf, bnd%bnd_cross, &
            bnd%node_cross, bnd%send_spec, bnd%recv_spec, ndump_fac, &
            restart, restart_handle, sim_options, tstep, tmax )
            species => species % next; i = i+1
        enddo
        
        ! setup additional buffers for low roundoff current deposition
        if (( this%num_species > 1 ) .and. (this%low_jay_roundoff)) then
            call this % jay_tmp % new( jay % pf(1) )
        else
            this%low_jay_roundoff = .false.
        endif
    endif

    ! setup cathodes
    if ( this%num_cathode > 0 ) then
        do cathode_id=1, this%num_cathode
            call this%cathode(cathode_id) % init( cathode_id, restart, restart_handle, t, &
            dt(tstep), grid%coordinates )
        enddo
    endif

    #ifdef __HAS_IONIZATION__

    ! setup neutrals
    if ( this%num_neutral > 0 ) then
        do neutral_id=1, this%num_neutral
            call this%neutral(neutral_id)%init( neutral_id, emf, grid%my_nx( p_lower, :), g_space, &
            restart, restart_handle, sim_options )
        enddo
    endif

    #endif

    #ifdef __HAS_COLLISIONS__

    ! setup collisions
    call setup( this%coll, grid%g_nx, grid%my_nx( 3, : ), restart, restart_handle )

    #endif

    ! setup events in this file
    if (pushev==0) then
        pushev            = create_event('advance deposit')
        reduce_current_ev = create_event('reduce current')
        
        partboundev       = create_event('update particle boundary')
        cathodeev         = create_event('cathode injection')
        neutralev         = create_event('neutral injection')
        diag_part_ev      = create_event('particle diagnostics')

        species => this % species
        
        ! setup events in the species class
        sortev            = create_event('particle sort (total)')
        sort_genidx_ev    = create_event('particle sort, gen. idx')
        sort_rearrange_ev = create_event('particle sort, rearrange particles')
    endif

end subroutine init_particles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! updates the particle data at the boundaries
!-----------------------------------------------------------------------------------------
subroutine update_boundary_particles( this, jay, g_space, no_co, dt, bnd )

    use m_system
    use m_particles_define
    use m_species
    use m_neutral
    use m_vdf_define
    use m_space
    use m_node_conf
    use m_logprof
    use m_current_define
    use m_bnd

    use m_species_define

    implicit none

    class( t_particles ), intent(inout) :: this
    class( t_current ), intent(inout) :: jay

    type( t_space ),     intent(in) :: g_space
    class( t_node_conf ), intent(in) :: no_co
    real(p_double),     intent(in) :: dt
    class(t_bnd),     intent(inout) :: bnd

    class(t_species), pointer :: species
    integer :: i

    call begin_event(partboundev)

    species => this % species
    do
        if (.not. associated(species)) exit
        call species % update_boundary( jay, no_co, dt, bnd%bnd_cross, bnd%node_cross, &
        bnd%send_spec, bnd%recv_spec )
        species => species % next
    enddo

    #ifdef __HAS_IONIZATION__

    do i=1, this%num_neutral
        call this%neutral(i) % update_boundary( nx_move(g_space), no_co, bnd%send_vdf, &
        bnd%recv_vdf )
    enddo

    #endif

    call end_event(partboundev)

end subroutine update_boundary_particles
!-----------------------------------------------------------------------------------------

#if 1

!-----------------------------------------------------------------------------------------
! Advance/deposit particles using multiple tasks per node (OpenMP)
!  - Low jay roundoff is currently not supported
!-----------------------------------------------------------------------------------------
subroutine advance_deposit_particles( this, emf, jay, tstep, t, no_co, options )

    use m_system
    use m_parameters

    use m_species_define
    use m_particles_define
    use m_particles
    use m_emf_define
    use m_time_step

    use m_vdf_define
    use m_vdf_math
    use m_vdf

    use m_current_define

    use m_logprof

    use m_node_conf

    implicit none

    class( t_particles ), intent(inout) :: this
    class( t_emf ), intent( inout )  ::  emf
    class( t_current ), intent(inout) :: jay

    real(p_double), intent(in) :: t
    type( t_time_step ), intent(in) :: tstep
    class( t_node_conf ), intent(in) :: no_co
    type( t_options ), intent(in) :: options

    class(t_species), pointer :: species

    integer :: nt, tid

    ! executable statements
#ifdef __DEBUG__
    call validate( this, 'before advance_deposit_particles' )

    if ( .not. associated(jay%pf) ) then
        ERROR('jay % pf is not associated')
        call abort_program()
    endif

    if ( n_threads( no_co ) /= size(jay%pf)) then
        ERROR('jay%pf size is invalid : ', size(jay%pf) )
        call abort_program()
    endif

    ! Sync nodes
    call no_co % barrier()
#endif

    nt = n_threads( no_co )
    if ( nt > 1 ) then
        
        call begin_event( pushev )

        ! Multi-threaded pusher       
        
        !$omp parallel do private(species)
        do tid = 1, nt
            ! Zero the current grid for this thread
            call jay % pf( tid ) % zero()
            
            species => this % species
            do
                if (.not. associated(species)) exit
                call species % push( emf, jay % pf(tid), t, tstep, tid-1, nt, options )
                species => species % next
            enddo
        enddo
        !$omp end parallel do
        
        call end_event( pushev )
        
        ! Add current from all threads - this is also done with OpenMP parallelism
        call begin_event( reduce_current_ev )
        call reduce( jay % pf )

        call end_event( reduce_current_ev )
        
    else
        
        call begin_event( pushev )

        ! Single-threaded pusher
                
        species => this % species
        
        call jay % pf(1) % zero()
        
        if ( this%low_jay_roundoff ) then
            do
                if (.not. associated(species)) exit
                call this % jay_tmp % zero()
                call species % push( emf, this%jay_tmp, t, tstep, 0, 1, options )
                call add( jay % pf(1), this%jay_tmp )
                species => species % next
            enddo
        else
            do
                if (.not. associated(species)) exit
                call species % push( emf, jay % pf(1), t, tstep, 0, 1, options )
                species => species % next
            enddo
        endif
        
        call end_event( pushev )
        
    endif

    ! Advance iteration information
    this%n_current = this%n_current + 1

#ifdef __DEBUG__
    call no_co % barrier()
    call validate( this, 'after advance_deposit_particles', over = .true. )
#endif

end subroutine advance_deposit_particles
!---------------------------------------------------------------------------------------------------

#else

#warning advance_deposit_particles() using non OpenMP version

!-----------------------------------------------------------------------------------------
! Advance/deposit particles using a single thread per node (no OpenMP support)
!-----------------------------------------------------------------------------------------
subroutine advance_deposit_particles( this, emf, jay, tstep, t, no_co, options )
    
    use m_system
    
    use m_species_define
    use m_particles_define
    use m_emf_define
    use m_time_step
    
    use m_vdf_define
    use m_vdf_math
    use m_vdf
    
    use m_current_define
    
    use m_logprof
    
    use m_node_conf
    
    implicit none
    
    ! dummy variables
    
    class( t_particles ), intent(inout) :: this
    class( t_emf ), intent( inout )  ::  emf
    class( t_current ), intent(inout) :: jay
    
    real(p_double), intent(in) :: t
    type( t_time_step ), intent(in) :: tstep
    class( t_node_conf ), intent(in) :: no_co
    type( t_options ), intent(in) :: options
    
    class(t_species), pointer :: species
    
    call begin_event( pushev )
    
    call jay % pf(1) % zero()
    
    species => this % species
    
    if ( this%low_jay_roundoff ) then
        do
            if (.not. associated(species)) exit
            call this % jay_tmp % zero()
            call species % push( emf, this%jay_tmp, t, tstep, 0, 1, options )
            call add( jay % pf(1), this%jay_tmp )
            species => species % next
        enddo
    else
        do
            if (.not. associated(species)) exit
            call species % push( emf, jay % pf(1), t, tstep, 0, 1, options )
            species => species % next
        enddo
    endif
    
    call end_event( pushev )
    
    ! Advance iteration information
    this%n_current = this%n_current + 1
    
end subroutine advance_deposit_particles
!---------------------------------------------------------------------------------------------------


#endif

!-----------------------------------------------------------------------------------------
! report particle data
!-----------------------------------------------------------------------------------------
subroutine report_particles( this, emf, g_space, grid, no_co, tstep, t, send_msg, recv_msg )
    
    use m_system
    use m_species_define
    use m_particles_define
    use m_particles_charge
    use m_particles
    use m_emf_define
    use m_space
    use m_grid_define
    use m_node_conf
    use m_time_step
    use m_vdf_define
    use m_vdf_report
    use m_vdf_comm, only : t_vdf_msg
    use m_logprof
    use m_neutral
    
    implicit none
    
    class( t_particles ), intent(inout) :: this
    class( t_emf ),       intent(inout) :: emf
    
    type( t_space ),      intent(in) :: g_space
    class( t_grid ),      intent(in) :: grid
    class( t_node_conf ),  intent(in) :: no_co
    type( t_time_step ),  intent(in) :: tstep
    real(p_double),       intent(in) :: t
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
    
    ! local variables
    type( t_vdf ) :: tmp
    integer, parameter :: izero = ichar('0')
    integer :: i
    class(t_species), pointer :: species
    logical :: needs_update_charge
    type( t_vdf_report ), pointer :: rep
    ! executable statements
    call begin_event(diag_part_ev)
    ! Species diagnostics (phasespaces, energy, temperature, heat flux, RAW, tracks, etc.)
    species => this % species
    do
        
        if (.not. associated(species)) exit
        call  species % report( emf, g_space, grid, no_co, tstep, t, send_msg, recv_msg )
        species => species % next
    enddo
    
    #ifdef __HAS_IONIZATION__
    do i=1, this%num_neutral
        call this%neutral(i) % report(g_space, grid, no_co, tstep, t )
    enddo
    #endif
    
    ! Check if charge deposit is required
    needs_update_charge = .false.
    rep => this%reports
    do
        if ( .not. associated( rep ) .or. needs_update_charge ) exit
        
        needs_update_charge = if_report( rep, tstep ) .or. &
        ( ( rep%quant == p_charge_htc .or. rep%quant == p_dcharge_dt ) .and. &
        if_report( rep, tstep, iter = +1 ) )
        
        rep => rep%next
    enddo
    
    ! Deposit charge if required
    if ( needs_update_charge ) then
        call update_charge( this, g_space, grid, no_co, send_msg, recv_msg )
    endif
    
    ! Do selected reports
    rep => this%reports
    do
        if ( .not. associated( rep ) ) exit
        
        if ( if_report( rep, tstep ) ) then
            
            ! deposit present charge
            
            select case ( rep%quant )
            case ( p_charge )
                call report_vdf( rep, this%charge%current, 1, g_space, grid, no_co, tstep, t )
                
            case ( p_charge_htc )
                call get_charge_htc( this, tmp )
                call report_vdf( rep, tmp, 1, g_space, grid, no_co, tstep, t )
                
            case ( p_dcharge_dt )
                call get_dcharge_dt( this, dt( tstep ), tmp )
                call report_vdf( rep, tmp, 1, g_space, grid, no_co, tstep, t )
                
            end select
        endif
        
        rep => rep%next
    enddo
    
    call tmp % cleanup( )
    
    call end_event(diag_part_ev)
    
end subroutine report_particles
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine report_energy_particles( this, no_co, tstep, t )
    
    use m_system
    use m_parameters, only: file_id_parene, p_err_diagfile, path_hist
    use m_species_define
    use m_particles_define
    use m_particles_charge
    use m_particles
    use m_emf_define
    use m_space
    use m_grid_define
    use m_node_conf
    use m_time_step
    use m_vdf_define
    use m_vdf_report
    use m_vdf_comm, only : t_vdf_msg
    use m_logprof
    use m_neutral
    
    implicit none
    
    class( t_particles ), intent(in) :: this
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t
    
    integer :: i, ierr
    integer, parameter :: p_header_length = 23
    real(p_double), dimension(:), pointer :: ene_all
    character(len = p_header_length), dimension(:), allocatable :: header
    character(len = 64) :: fmt
    character(len = 80) :: path, full_name
    
    character(len=*), parameter :: p_err_msg = "Unable to create directory for particle energy diagnostic"
    
    class( t_species ), pointer :: species
    
    call begin_event( diag_part_ev )
    
    ! Report global particle energy
    if ( test_if_report( tstep, this%ndump_fac_ene ) ) then
        
        call alloc( ene_all, (/ this%num_species + 1 /) )
        
        ! Get normalized kinetic energy for all species and calculate total part. energy
        ene_all(1) = 0
        species => this % species
        i = 1
        do
            if (.not. associated(species)) exit
            ene_all(i+1) = species % get_energy()
            ene_all(1) = ene_all(1) + ene_all(i+1)
            species => species % next
            i = i + 1
        enddo
        
        
        ! Gather results on root node
        call reduce_array( no_co, ene_all )
        
        ! Write output to file
        if ( root(no_co) ) then
            
            ! prepare path and file names
            
            path  = trim(path_hist)
            
            full_name = trim(path) // 'par_ene'
            
            ! Open file and position at the last record
            if ( t == 0.0_p_double ) then
                
                call mkdir( path, ierr )
                CHECK_ERROR( ierr, p_err_msg, p_err_diagfile )
                
                open (unit=file_id_parene, file=full_name, status = 'REPLACE' , form='formatted')
                
                ! Write Header
                
                write(file_id_parene, '(A)' ) '! Total and per species particle kinetic energy'
                
                allocate( header( this%num_species ) )
                
                write( fmt, '(A,I0,A,I0,A)' ) '( A6, (1X,A15), ', this%num_species + 1, &
                                                '(1X, A', p_header_length, ') )'
                
                i = 1
                species => this % species
                do
                    if (.not. associated(species)) exit
                    header(i) = trim(species%name)
                    header(i) = adjustr(header(i))
                    species => species % next
                    i = i+1
                enddo
                
                write(file_id_parene, fmt ) 'Iter','Time    ', 'Total', header
                deallocate( header )
            else
                
                open (unit=file_id_parene, file=full_name, position = 'append', &
                form='formatted')
            endif
            
            !  Writes the particle energies to file
            write( fmt, '(A,I0,A)' ) '( I6, (1X,g15.8), ', this%num_species + 1, '(1X, es23.16) )'
            write(file_id_parene, fmt) n(tstep), t, ene_all
            
            ! Close file
            close(file_id_parene)
            
        endif
        
        call freemem( ene_all )
        
    endif
    
    ! Report individual species energy
    species => this % species
    do
        if (.not. associated(species)) exit
        call species % report_energy( no_co, tstep, t )
        species => species % next
    enddo
    
    call end_event( diag_part_ev )
    
end subroutine report_energy_particles
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get diagnostics buffer size requirements
!-----------------------------------------------------------------------------------------
subroutine get_diag_buffer_size_part( this, gnx, diag_buffer_size )
    
    use m_particles_define
    use m_species_define
    
    implicit none
    
    class( t_particles ), intent(in) :: this
    integer, dimension(:), intent(in) :: gnx
    integer, intent(inout) :: diag_buffer_size
    
    class(t_species), pointer :: species
    
    species => this % species
    do
        if (.not. associated(species)) exit
        call species % get_diag_buffer_size( gnx, diag_buffer_size )
        species => species % next
    enddo
    
end subroutine get_diag_buffer_size_part
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleanup particle species, cathodes, etc.
!-----------------------------------------------------------------------------------------
subroutine cleanup_particles( this )
    
    use m_particles_define
    use m_particles_charge
    use m_species
    use m_species_define
    use m_species_collisions
    use m_species_current
    use m_cathode
    use m_vdf
    use m_vdf_report, only : cleanup
    
    implicit none
    
    class( t_particles ), intent(inout) :: this
    
    integer :: i
    class( t_species ), pointer :: species, next

    ! cleanup temp current buffers
    call this%jay_tmp%cleanup()
    if (associated(this%species)) call this % species % bnd_con % cleanup_tmp_buf_current()
    
    ! free species data
    if ( this%num_species > 0 ) then
        
        ! cleanup buffers for advance_deposit and get_jr
        call cleanup_buffers_spec()
        
        ! cleanup species
        species => this % species
        do
            if (.not. associated(species)) exit
            next => species % next
            call species % cleanup()
            deallocate(species)
            species => next
        enddo
        
        this % species => null()
    endif
    
    ! free cathode data
    if ( this%num_cathode > 0 ) then
        
        ! cleanup cathode
        do i=1, this%num_cathode
            call  this % cathode(i) % cleanup()
        enddo
        
        ! free cathode buffer
        ! call freemem( this%cathode )
        deallocate( this%cathode )
    endif
    
#ifdef __HAS_IONIZATION__
    ! free neutral data
    if ( this%num_neutral > 0 ) then
        
        ! cleanup neutral
        do i=1, this%num_neutral
            call  this%neutral(i)%cleanup()
        enddo
        
        ! free species buffer
        ! call freemem( this%neutral )
        deallocate( this%neutral )
    endif
#endif
    
    ! cleanup charge buffers
    call cleanup( this%charge )
    
    call cleanup( this%reports )
    
#ifdef __HAS_COLLISIONS__
    
    ! cleanup collision data
    call cleanup( this%coll )
    
#endif
    
end subroutine cleanup_particles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! write checkpoint information
!-----------------------------------------------------------------------------------------
subroutine write_checkpoint_particles( this, restart_handle )
    
    use m_species_define
    use m_particles_define
    use m_particles_charge
    use m_restart
    use m_parameters
    
    #ifdef __HAS_IONIZATION__
    use m_neutral
    #endif
    
    #ifdef __HAS_COLLISIONS__
    use m_species_collisions
    #endif
    
    implicit none
    
    class( t_particles ), intent(in) :: this
    type( t_restart_handle ), intent(inout) :: restart_handle
    
    character(len=*), parameter :: err_msg = 'error writing checkpoint data for particles object.'
    class( t_species ), pointer :: species
    integer :: i, ierr
    
    restart_io_wr( p_part_rst_id, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
    
    ! write local data to file
    restart_io_wr( this%interpolation, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
    
    restart_io_wr( this%num_species, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
    
    restart_io_wr( this%num_cathode, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
    
    restart_io_wr( this%n_current, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
    
    #ifdef __HAS_IONIZATION__
    restart_io_wr( this%num_neutral, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstwrt )
    #endif
    
    call restart_write( this%charge, restart_handle )
    
    ! write contained objects data
    species => this % species
    do
        if (.not. associated(species)) exit
        call species % restart_write( restart_handle )
        species => species % next
    enddo
    
    do i=1, this%num_cathode
        call this % cathode(i) % write_checkpoint( restart_handle )
    enddo
    
#ifdef __HAS_IONIZATION__
    do i=1, this%num_neutral
        call this%neutral(i)%restart_write( restart_handle )
    enddo
#endif
    
#ifdef __HAS_COLLISIONS__
    call restart_write( this%coll, restart_handle )
#endif
    
end subroutine write_checkpoint_particles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! read object information from a restart file
!-----------------------------------------------------------------------------------------
subroutine restart_read_particles( this, restart_handle )

    use m_species_define
    use m_particles_define
    use m_particles_charge
    use m_restart
    use m_parameters

    #ifdef __HAS_IONIZATION__
    use m_neutral
    #endif
    
    #ifdef __HAS_COLLISIONS__
    use m_species_collisions
    #endif

    implicit none

    class( t_particles ), intent(inout) :: this
    type( t_restart_handle ), intent(in) :: restart_handle

    ! local variables
    character(len=*), parameter :: err_msg = 'error reading restart data for particles object.'
    integer :: interpolation, num_species, num_cathode, num_neutral
    integer :: ierr

    character(len=len(p_part_rst_id)) :: rst_id

    restart_io_rd( rst_id, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd )

    ! check if restart file is compatible
    if ( rst_id /= p_part_rst_id) then
        ERROR('Corrupted restart file, or restart file ')
        ERROR('from incompatible binary (part)')
        ERROR('rst_id = ', rst_id)
        call abort_program(p_err_rstrd)
    endif

    restart_io_rd( interpolation, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd )

    restart_io_rd( num_species, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd )

    restart_io_rd( num_cathode, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd )

    restart_io_rd( this%n_current, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd )

#ifdef __HAS_IONIZATION__

    restart_io_rd( num_neutral, restart_handle, ierr )
    CHECK_ERROR( ierr, err_msg, p_err_rstrd )

    if (this%num_neutral /= num_neutral) then
    ERROR('The number of neutrals specified in the input deck is different')
    ERROR('input : ', this%num_neutral, ' rst: ', num_neutral )
    ERROR('from the number of neutrals in the restart file')
    call abort_program(p_err_invalid)
    endif

#endif

    if ( this%interpolation /= interpolation ) then
        ERROR('The interpolation level in the input deck is different')
        ERROR('from the interpolation level in the restart file')
        call abort_program(p_err_invalid)
    endif

    if (this%num_species /= num_species) then
        ERROR('The number of species specified in the input deck is different')
        ERROR('from the number of species in the restart file')
        call abort_program(p_err_invalid)
    endif

    if (this%num_cathode /= num_cathode) then
        ERROR('The number of cathodes specified in the input deck is different')
        ERROR('from the number of cathodes in the restart file')
        call abort_program(p_err_invalid)
    endif

end subroutine restart_read_particles
!-----------------------------------------------------------------------------------------