!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     electro-magnetic field boundary condition class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_bound

  use m_emf_define

  use m_logprof
  use m_lindman

  use m_space
  use m_node_conf
  use m_grid_define

  use m_parameters

  use m_vdf_define
  use m_vdf_comm

  use m_vpml

  implicit none

  private

  integer, parameter :: p_vpml_size_default = 8

  interface report_bcemf
    module procedure report_bcemf
  end interface

  interface type
    module procedure type_emf_bound
    module procedure type_emf_bound_bnd
  end interface

  interface reshape_obj
    module procedure reshape_bcemf
  end interface

  ! used by wall antennae
  interface is_open
    module procedure is_open_bcemf
  end interface

  interface bnd_type_name
    module procedure bnd_type_name
  end interface

  integer :: emfboundev = 0

  public :: t_emf_bound, emfboundev, p_vpml_size_default

  public :: read_nml
  public :: report_bcemf, bnd_type_name
  public :: type, is_open, reshape_obj

contains

!-----------------------------------------------------------------------------------------
! Get the name of a boundary type
!-----------------------------------------------------------------------------------------
function bnd_type_name( type )

  implicit none

  integer, intent(in) :: type
  character(len = 39) :: bnd_type_name

  select case ( type )
    case( p_bc_other_node )
      bnd_type_name = "Other node"
    case( p_bc_periodic )
      bnd_type_name = "Periodic"

    case( p_bc_pec )
      bnd_type_name = "Conducting (perfect electric conductor)"
    case( p_bc_pmc )
      bnd_type_name = "Reflecting (perfect magnetic conductor)"

    case( p_bc_cyl_axis )
      bnd_type_name = "Cyl. axis"
    case( p_bc_lindman )
      bnd_type_name = "Lindman"
    case( p_bc_vpml )
      bnd_type_name = "Perfectly matched layer"
    case( p_bc_move_c )
      bnd_type_name = "Moving window"
    case default
      bnd_type_name = "unknown"
  end select

end function bnd_type_name
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       report on electro-magnetic field - diagnostic
!-----------------------------------------------------------------------------------------
subroutine report_bcemf( this, g_space, grid, no_co, tstep, t, dx )

  use m_time_step

  implicit none

  class (t_emf_bound), intent(inout) :: this
  type( t_space ),      intent(in) :: g_space
  class( t_grid ),        intent(in) :: grid
  class( t_node_conf ),  intent(in) :: no_co
  type( t_time_step ),  intent(in) :: tstep
  real(p_double),      intent(in) :: t
  real(p_double), dimension(:), intent(in) :: dx


  ! report VPML
    call report( this%vpml_all, g_space, no_co, tstep, t, dx )

end subroutine report_bcemf

!-----------------------------------------------------------------------------------------
!       gives the boundary type variable for this boundary condition
!-----------------------------------------------------------------------------------------
function type_emf_bound( this )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, dimension(2,p_x_dim)  ::  type_emf_bound

  class (t_emf_bound),intent(in) :: this

  type_emf_bound(:,1:p_x_dim) = this%type(:,1:p_x_dim)

end function type_emf_bound
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       gives the boundary type variable for

!       the specified boundary

!-----------------------------------------------------------------------------------------
function type_emf_bound_bnd( this , i_bnd, i_dir )
!-----------------------------------------------------------------------------------------

  implicit none

  integer ::  type_emf_bound_bnd

  class (t_emf_bound),intent(in) :: this
  integer , intent(in) :: i_bnd,i_dir

  type_emf_bound_bnd = this%type(i_bnd,i_dir)

end function type_emf_bound_bnd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
logical function is_open_bcemf(this,idir,iside)
!-----------------------------------------------------------------------------------------
implicit none

class (t_emf_bound), intent(in) :: this
integer, intent(in) :: idir,iside

  is_open_bcemf= ((this%type(iside,idir) == p_bc_lindman) .or. &
                  (this%type(iside,idir) == p_bc_vpml))

end function is_open_bcemf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reshape_bcemf( this, old_lb, new_lb, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

   implicit none

   class (t_emf_bound), intent(inout) :: this
   class(t_grid), intent(in) :: old_lb, new_lb
   class( t_node_conf ), intent(in) :: no_co
   type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

   integer :: i_bnd

   do i_bnd = p_lower, p_upper
     if ( is_active( this%lindman(i_bnd) ) ) then
       call reshape_obj( this%lindman(i_bnd), old_lb, new_lb, no_co, send_msg, recv_msg )
     endif
   enddo

   if ( is_active( this%vpml_all ) ) then
     call reshape_obj( this%vpml_all, old_lb, new_lb, no_co, send_msg, recv_msg )

   endif

end subroutine reshape_bcemf
!-----------------------------------------------------------------------------------------
end module m_emf_bound

!---------------------------------------------------
subroutine read_input_emf_bound( this, input_file, periodic, if_move )

  use m_emf_bound, only: t_emf_bound, p_vpml_size_default
  use m_input_file

  use m_parameters

  implicit none

  ! dummy variables
  class (t_emf_bound), intent( inout ) :: this
  class( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic, if_move

  ! local variables
  character(20), dimension(2, p_x_dim) :: type

  integer :: vpml_bnd_size
  logical, dimension(2,p_x_dim) :: vpml_diffuse

  namelist /nl_emf_bound/ type, vpml_bnd_size, vpml_diffuse

  integer :: ierr, i, j, have_lind, have_cond

  ! executable statements
  type = "---"
  vpml_bnd_size = p_vpml_size_default
  vpml_diffuse = .false.

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_emf_bound", ierr )

  if ( ierr /= 0 ) then
    if (ierr < 0) then
      print *, "Error reading emf_bound parameters"
    else

      print *, "Error: emf_bound parameters missing"
    endif
    print *, "aborting..."
    stop
  endif

  read (input_file%nml_text, nml = nl_emf_bound, iostat = ierr)
  if (ierr /= 0) then
    print *, "Error reading emf_bound parameters"
    print *, "aborting..."
    stop
  endif

  ! Read boundary condition types
  do i = 1, p_x_dim
    do j = 1, 2
      select case ( trim( type(j,i) ) )
          case ( "axial" )
            this%type(j,i) = p_bc_cyl_axis
          case ( "lindman", "lindmann" )
            this%type(j,i) = p_bc_lindman
          case ( "vpml", "open" )
            this%type(j,i) = p_bc_vpml

          case ( "pec", "conducting" )

            this%type(j,i) = p_bc_pec
          case ( "pmc", "reflecting" )
            this%type(j,i) = p_bc_pmc

          case ( "---" )
            this%type(j,i) = p_bc_invalid
          case default
            if ( mpi_node() == 0 ) then
               print "(A,I0,A,I0)", "    Error reading emf_bound parameters, dir = ", i, ", wall = ", j
               print *, '   Unknown boundary condition type, "', trim(type(i,j)), &
                        '", valid values are "conducting", "axial", "lindman", "vpml", "open", "pec", "pmc" and "reflecting"'
               print *, '   aborting...'
            endif
            stop
      end select
    enddo
  enddo

  ! Verify boundary conditions
  do i = 1, p_x_dim
    ! periodic boundaries and moving window will override local settings
    if ((.not. periodic(i)) .and. (.not. if_move(i))) then

       do j = 1, 2
         if (this%type(j,i) == p_bc_invalid) then
           print *, ""
           print *, "   Error reading emf_bound parameters, dir = ", i, ", wall = ", j
           print *, "   Invalid or no boundary conditions specified."
           print *, "   aborting..."
           stop

         elseif (this%type(j,i) == p_bc_periodic) then
           print *, ""
           print *, "   Error reading emf_bound parameters, dir = ", i
           print *, "   Periodic boundaries cannot be specified unless global periodic"
           print *, "   boundaries were set in the node_conf section"
           print *, "   aborting..."
           stop

         endif
       enddo
    endif
  enddo

  ! check on Lindman boundaries
  have_lind = 0
  do i = 1, p_x_dim
     ! periodic boundaries and moving window override local settings
     if ((.not. periodic(i)) .and. (.not. if_move(i)) .and. &
         ((this%type(p_lower,i) == p_bc_lindman) .or. &
          (this%type(p_upper,i) == p_bc_lindman))) then

        if (have_lind > 0) then

           print *, ""
           print *, "   Error reading emf_bound parameters, dir = ", i
           print *, "   Lindman boundaries cannot be specified for more than"
           print *, "   one direction (Lindman bound. cond. already specified"
           print *, "   for dir =",have_lind, ")."
           print *, "   aborting..."
           stop

        else

           have_lind = i
        endif

     endif
  enddo

  ! check on cylindrical boundaries
  do i = 1, p_x_dim
     do j = 1, 2
       if ((this%type(j,i) == p_bc_cyl_axis) .and. &
           ((i /= p_r_dim) .or. ( j /= p_lower ))) then
           print *,        ""
           print "(A,I0)", "   Error reading emf_bound parameters, dir = ", i
           print "(A)",    "   Cylindrical axial boundaries can only be specified for"
           print "(A,I0)", "   the lower boundary of dimension ",p_r_dim
           print "(A)",    "   aborting..."
           stop

       endif
     enddo
  enddo

  ! check on VPML boundaries
  this%vpml_bnd_size = vpml_bnd_size
  this%vpml_diffuse(:,1:p_x_dim)   =  vpml_diffuse(:,1:p_x_dim)

  have_cond = -1
  do i=1, p_x_dim

     ! check if VPML is contacting with conducting / reflecting
     if ( (this%type( p_lower, i ) == p_bc_pec) .or. &
          (this%type( p_upper, i ) == p_bc_pec) .or. &
          (this%type( p_lower, i ) == p_bc_pmc) .or. &
          (this%type( p_upper, i ) == p_bc_pmc) ) then

        have_cond = i

     endif

     if ( (this%type( p_lower, i ) == p_bc_vpml) &
         .or. (this%type( p_upper, i ) == p_bc_vpml) ) then

       if (( have_cond > 0 ) .and. ( have_cond /= i )) then
           print *, ""
           print *, "   Error reading emf_bound parameters, dir = ", i
           print *, "   Conducting or reflecting boundaries contacting with VPML are not implemented yet."
           print *, "   aborting..."
           stop

       endif

     endif
  enddo
end subroutine read_input_emf_bound
!---------------------------------------------------

!---------------------------------------------------
subroutine init_emf_bound( this, b, e, dt, part_interpolation, &
                            space, no_co, grid, restart, restart_handle )

   use m_emf_bound

   use m_restart
   use m_emf_define

   use m_logprof
   use m_lindman

   use m_space
   use m_node_conf
   use m_grid_define

   use m_parameters

   use m_vdf_define
   use m_vdf_comm

   use m_vpml

  implicit none

  ! dummy variables

  class (t_emf_bound), intent( inout )  ::  this

  type( t_vdf ),       intent(in) :: b, e
  real( p_double ), intent(in) :: dt
  integer, intent(in) :: part_interpolation
  type( t_space ),     intent(in) :: space
  class( t_node_conf ), intent(in) :: no_co
  class( t_grid ), intent(in) :: grid

  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

  ! local variables

  logical, dimension(p_x_dim) :: ifpr_l, if_move_l
  logical, dimension(2,p_x_dim) :: if_vpml
  integer :: i_bnd, i_x_dim, previous, n_vpml

  ! executable statements

  if ( restart ) then

     call this%restart_read( restart_handle )

  else

     ifpr_l = no_co%ifpr()
     if_move_l = if_move(space)

     ! overwrite b.c. flags if motion of simulation space
     do i_x_dim = 1, p_x_dim
        if ( if_move_l(i_x_dim) ) then
           this%type( p_lower , i_x_dim ) = p_bc_move_c
           this%type( p_upper , i_x_dim ) = p_bc_move_c
        endif
     enddo

     ! Store the global type of boundary condition used
     this%global_type = this%type
     do i_x_dim = 1, p_x_dim
        if ( periodic( no_co, i_x_dim) ) this%global_type( 1:2 , i_x_dim ) = p_bc_periodic
     enddo

     ! overwrite b.c.-flags if single node periodicity
     do i_x_dim = 1, p_x_dim
        if ( ifpr_l(i_x_dim) ) then
           this%type( : , i_x_dim ) = p_bc_periodic
        endif
     enddo

     ! overwrite b.c. flags for boundaries to other nodes
     do i_x_dim = 1, p_x_dim
        do i_bnd=1, 2
           if ( neighbor( no_co, i_bnd, i_x_dim ) > 0 ) then
              this%type( i_bnd, i_x_dim ) = p_bc_other_node
           endif
        enddo
     enddo

     ! do additional cylindrical coordinates setup

     if ( grid % coordinates == p_cylindrical_b ) then

        ! If node has the axial boundary set the bc accordingly
        if ( no_co%on_edge( p_r_dim, p_lower )) then
           this%type(p_lower, p_r_dim) = p_bc_cyl_axis
        endif

        ! check if we have periodic boundaries along radial direction
        if ((this%type(p_lower, p_r_dim) == p_bc_periodic) .or. &
            (this%type(p_upper, p_r_dim) == p_bc_periodic)) then
           ERROR('no periodic b.c. for radial direction')
           call abort_program( p_err_invalid )
        endif

     endif

     ! When using PEC or PMC boundaries adjust to the correct one according to the

     ! particle interpolation
     do i_x_dim = 1, p_x_dim
       do i_bnd=1, 2
          select case ( this%type(i_bnd,i_x_dim) )
             case ( p_bc_pec )
               if ( mod( part_interpolation, 2 ) == 1 ) then
                 this%type(i_bnd,i_x_dim) = p_bc_pec_odd
               else
                 this%type(i_bnd,i_x_dim) = p_bc_pec_even
               endif
             case ( p_bc_pmc )
               if ( mod( part_interpolation, 2 ) == 1 ) then
                 this%type(i_bnd,i_x_dim) = p_bc_pmc_odd
               else
                 this%type(i_bnd,i_x_dim) = p_bc_pmc_even
               endif
             case default
               continue
          end select
       enddo
     enddo

  endif

  ! setup lindman boundaries
  previous = -1

  do i_x_dim=1, p_x_dim
    do i_bnd = p_lower, p_upper
      if ( this%type( i_bnd, i_x_dim ) == p_bc_lindman ) then
        ! check if lindman boundaries were already defined in another direction
        ! (this is redundant because it was already checked when reading the input file)
        if ( previous > 0 .and. previous /= i_x_dim ) then
           ERROR("Lindman boundaries are only allowed in 1 direction")
           call abort_program( p_err_invalid )
        endif

        ! setup lindman boundary object
        call setup( this%lindman(i_bnd), i_bnd, i_x_dim, b, e, dt, restart, restart_handle )
        previous = i_x_dim
      endif
    enddo

  enddo

  ! setup VPML boundaries
  n_vpml = 0
  if_vpml = .false.
  do i_x_dim=1, p_x_dim
     if ( this%type( p_lower, i_x_dim ) == p_bc_vpml ) then
       if_vpml(p_lower,i_x_dim) = .true.
       n_vpml = n_vpml + 1
     endif

     if ( this%type( p_upper, i_x_dim ) == p_bc_vpml ) then
       if_vpml(p_upper,i_x_dim) = .true.
       n_vpml = n_vpml + 1
     endif
  enddo

  ! setup if at least one VPML
  if (n_vpml > 0 ) then
    call setup( this%vpml_all, e, b, grid, n_vpml, this%vpml_bnd_size, this%vpml_diffuse, &
                if_vpml, if_move_l, dt, restart, restart_handle )
  endif

  ! setup events
  if (emfboundev==0) emfboundev = create_event('update emf boundary')

end subroutine init_emf_bound
!---------------------------------------------------

!---------------------------------------------------
subroutine cleanup_emf_bound( this )
!---------------------------------------------------
!       release memory used by ptr-components of this variable
!---------------------------------------------------
  use m_emf_bound
  use m_emf_define
  use m_lindman
  use m_parameters
  use m_vpml

  implicit none

  class (t_emf_bound), intent( inout ) :: this

  ! cleanup lindman data
  call cleanup( this%lindman(p_lower) )
  call cleanup( this%lindman(p_upper) )

  ! cleanup VPML
  call cleanup( this%vpml_all )

end subroutine cleanup_emf_bound
!---------------------------------------------------

!---------------------------------------------------
subroutine write_checkpoint_emf_bound( this, restart_handle )
!---------------------------------------------------
!       write object information into a restart file
!---------------------------------------------------

  use m_parameters
  use m_emf_define, only : t_emf_bound, p_bcemf_rst_id
  use m_restart
  use m_lindman, only : restart_write
  use m_vpml, only : restart_write, is_active

  implicit none

  class (t_emf_bound), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for emf_bound object.'
  integer :: i, bnd, ierr

  restart_io_wr( p_bcemf_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%global_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! write lindman boundary data if needed
  do i = 1, p_x_dim
    do bnd = p_lower, p_upper
      if ( this%type( bnd, i ) == p_bc_lindman ) then
        call restart_write(this%lindman(bnd), restart_handle )
      endif
    enddo
  enddo

#ifdef __HAS_PML__

  ! write vpml boundary data if needed
  if(is_active(this%vpml_all)) then
    call restart_write(this%vpml_all, restart_handle )
  endif

#endif

end subroutine write_checkpoint_emf_bound
!---------------------------------------------------

!---------------------------------------------------
subroutine restart_read_emf_bound( this, restart_handle )
!---------------------------------------------------
!       read object information from a restart file
!---------------------------------------------------

  use m_parameters
  use m_emf_define, only : t_emf_bound, p_bcemf_rst_id
  use m_restart

  implicit none

  class (t_emf_bound), intent(inout)  ::  this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for emf_bound object.'
  character(len=len(p_bcemf_rst_id)) :: rst_id
  integer :: ierr

  ! check if restart file is compatible
  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  if ( rst_id /= p_bcemf_rst_id) then
    ERROR('Corrupted restart file, or restart file ')
    ERROR('from incompatible binary (emf_bound)')
    ERROR('rst_id = "', rst_id, '"' )
    call abort_program(p_err_rstrd)
  endif

  restart_io_rd( this%type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%global_type, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

end subroutine restart_read_emf_bound
!---------------------------------------------------

!---------------------------------------------------
!       move boundary of electromagnetic field
!---------------------------------------------------
subroutine move_window_emf_bound( this, space )

 use m_restart
 use m_emf_define
 use m_logprof
 use m_lindman

 use m_space
 use m_node_conf
 use m_grid_define

 use m_parameters

 use m_vdf_define
 use m_vdf_comm

 use m_vpml

  implicit none

  class (t_emf_bound), intent(inout) :: this
  type( t_space ),  intent(in) :: space

  integer :: bnd


  ! move lindman boundary data
  do bnd = p_lower, p_upper
    if( is_active(this%lindman(bnd)) ) call move_window(this%lindman(bnd), space )
  enddo


  ! move vpml boundary data
  if ( is_active(this%vpml_all) ) then
     call move_window( this%vpml_all, space )
  endif

end subroutine move_window_emf_bound
!---------------------------------------------------

!---------------------------------------------------
subroutine update_boundary_emf_bound( this, b, e, g_space, no_co, send_msg, recv_msg )
!---------------------------------------------------
!       update boundary of electromagnetic field
!---------------------------------------------------
   use m_parameters
   use m_restart
   use m_emf_define

   use m_logprof
   use m_lindman

   use m_space
   use m_node_conf
   use m_grid_define

   use m_parameters

   use m_vdf_define
   use m_vdf_comm
   use m_emf_bound

   use m_vpml

  implicit none

!       dummy variables and functions

  class (t_emf_bound), intent(inout) :: this
  type( t_vdf ),       intent(inout) :: b, e

  type( t_space ),     intent(in) :: g_space
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i, my_node
  integer, dimension(p_max_dim) :: lnx_move
  integer :: lneighbor, uneighbor

  call begin_event(emfboundev)

  lnx_move( 1 : x_dim( g_space ) ) = nx_move( g_space )

  ! update to the lindman boundaries must occur before updating em fields to
  ! get the corner value right

  call update_boundary(this%lindman(p_lower), no_co, lnx_move, send_msg, recv_msg )
  call update_boundary(this%lindman(p_upper), no_co, lnx_move, send_msg, recv_msg )

  ! Update VPML walls if necessary
  if (is_active(this%vpml_all)) then
    call update_boundary(this%vpml_all, no_co, lnx_move, send_msg, recv_msg )
  endif

  ! Update other boundary types

  my_node = no_co % my_aid()

  do i = 1, e%x_dim_

    if ( my_node == neighbor(no_co,p_lower,i) ) then

       !single node periodic
       call update_periodic_vdf( e, i, p_vdf_replace )
       call update_periodic_vdf( b, i, p_vdf_replace )

    else

       ! communication with other nodes
       lneighbor = neighbor( no_co, p_lower, i )
       uneighbor = neighbor( no_co, p_upper, i )

       ! post receives
       if ( lneighbor > 0 ) call irecv( e, p_vdf_replace, i, p_lower, no_co, lnx_move(i), recv_msg, b )
       if ( uneighbor > 0 ) call irecv( e, p_vdf_replace, i, p_upper, no_co, lnx_move(i), recv_msg, b )

       ! post sends
       ! Note: sends MUST be posted in the opposite order of receives to account for a 2 node

       ! periodic partition, where 1 node sends 2 messages to the same node
       if ( uneighbor > 0 ) call isend( e, p_vdf_replace, i, p_upper, no_co, lnx_move(i), send_msg, b )
       if ( lneighbor > 0 ) call isend( e, p_vdf_replace, i, p_lower, no_co, lnx_move(i), send_msg, b )

       ! Process lindman boundaries if needed
       if (this%type( p_lower, i ) == p_bc_lindman) then
          call update_emfbound( this%lindman( p_lower ), b, e, i)
       endif
       if (this%type( p_upper, i ) == p_bc_lindman) then
          call update_emfbound( this%lindman( p_upper ), b, e, i)
       endif

       ! wait for receives and unpack data
       if ( lneighbor > 0 ) call irecv_wait( e, p_lower, recv_msg, b )
       if ( uneighbor > 0 ) call irecv_wait( e, p_upper, recv_msg, b )

       ! wait for sends to complete
       if ( uneighbor > 0 ) call isend_wait( e, p_upper, send_msg )
       if ( lneighbor > 0 ) call isend_wait( e, p_lower, send_msg )

    endif

  enddo

  call end_event(emfboundev)

  ! call e%check_nan(" Checking E field after update boundary")
  ! call b%check_nan(" Checking B field after update boundary")

end subroutine update_boundary_emf_bound
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update physical emf boundaries after each B field advance. Used for:
!  - Open boundaries (Lindman or PML)
!  - Conducting (to be removed)
!  - Perfect electric conductor (PEC)

!-----------------------------------------------------------------------------------------
subroutine update_boundary_bfld_emf_bound( this, e, b, step )

  use m_emf_pec
  use m_emf_pmc
   use m_restart
   use m_emf_define

   use m_logprof
   use m_lindman

   use m_space
   use m_node_conf
   use m_grid_define

   use m_parameters

   use m_vdf_define
   use m_vdf_comm

   use m_vpml

  implicit none

  class( t_emf_bound ), intent(inout) :: this
  class( t_vdf ), intent(inout) :: e, b
  integer,             intent(in   ) :: step

  ! local variables
  integer :: bnd, i_dim

  ! executable statements
  select case ( b%x_dim_ )

   case (1)

     i_dim = 1

     do bnd = p_lower, p_upper
       select case (this%type( bnd, i_dim ))
         case (p_bc_lindman)
           call update_b_1d(this%lindman(bnd), e)

         case ( p_bc_pec_odd )
           call pec_bc_b_1d_odd( b, bnd )

         case ( p_bc_pec_even )
           call pec_bc_b_1d_even( b, bnd )

         case ( p_bc_pmc_odd )
           call pmc_bc_b_1d_odd( b, bnd )

         case ( p_bc_pmc_even )
           call pmc_bc_b_1d_even( b, bnd )

         case default

       end select
     enddo

     ! advance vpml E if necessary

     if ( (step == 1) .and. is_active(this%vpml_all) ) then
         call update_e_1d( this%vpml_all, b )

     endif

   case (2)

     do i_dim = 1, 2
       do bnd = p_lower, p_upper
         select case (this%type( bnd, i_dim ))
           case (p_bc_lindman)
             call update_b_2d(this%lindman(bnd), e)

           case ( p_bc_pec_odd )
             call pec_bc_b_2d_odd( b, i_dim, bnd )

           case ( p_bc_pmc_odd )
             call pmc_bc_b_2d_odd( b, i_dim, bnd )

           case ( p_bc_pec_even )
             call pec_bc_b_2d_even( b, i_dim, bnd )

           case ( p_bc_pmc_even )
             call pmc_bc_b_2d_even( b, i_dim, bnd )

           case default
         end select
       enddo
     enddo


     ! advance vpml E if necessary (cylindrical require all emf for radii correction)
     if ( (step == 1) .and. is_active(this%vpml_all) ) then
         call update_e_2d( this%vpml_all, b )
     endif

   case (3)
     do i_dim = 1, 3
       do bnd = p_lower, p_upper
         select case (this%type( bnd, i_dim ))
           case (p_bc_lindman)
             call update_b_3d(this%lindman(bnd), e)

           case ( p_bc_pec_odd )
             call pec_bc_b_3d_odd( b, i_dim, bnd )

           case ( p_bc_pmc_odd )
             call pmc_bc_b_3d_odd( b, i_dim, bnd )

           case ( p_bc_pec_even )
             call pec_bc_b_3d_even( b, i_dim, bnd )

           case ( p_bc_pmc_even )
             call pmc_bc_b_3d_even( b, i_dim, bnd )

           case default
         end select
       enddo
     enddo

     ! advance vpml E if necessary

     if ( (step == 1) .and. is_active(this%vpml_all) ) then
         call update_e_3d( this%vpml_all, b )

     endif

  end select

end subroutine update_boundary_bfld_emf_bound
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update physical emf boundaries after each E field advance. Used for:
!  - Open boundaries (Lindman or PML)
!  - Conducting (to be removed)
!  - Perfect electric conductor (PEC)

!-----------------------------------------------------------------------------------------
subroutine update_boundary_efld_emf_bound( this, e, b )

  use m_emf_pec
  use m_emf_pmc
     use m_restart
   use m_emf_define

   use m_logprof
   use m_lindman

   use m_space
   use m_node_conf
   use m_grid_define

   use m_parameters

   use m_vdf_define
   use m_vdf_comm
   use m_vpml

  implicit none

  class( t_emf_bound ), intent(inout) :: this
  class( t_vdf ), intent(inout) :: e, b

  ! local variables
  integer :: bnd, i_dim

  ! executable statements

  select case ( e%x_dim_ )

   case (1)

     i_dim = 1

     ! process lower boundary
     do bnd = p_lower, p_upper
       select case (this%type( bnd, i_dim ))
         case (p_bc_lindman)
           call update_e_1d(this%lindman(bnd), e)

         case ( p_bc_pec_odd )
           call pec_bc_e_1d_odd( e, bnd )

         case ( p_bc_pmc_odd )
           call pmc_bc_e_1d_odd( e, bnd )

         case ( p_bc_pec_even )
           call pec_bc_e_1d_even( e, bnd )

         case ( p_bc_pmc_even )
           call pmc_bc_e_1d_even( e, bnd )

         case default
       end select
     enddo

     ! advance vpml B if necessary

     if ( is_active(this%vpml_all) ) then
         call update_b_1d( this%vpml_all, e )

     endif

   case (2)

     do i_dim = 1, 2
       do bnd = p_lower, p_upper
         select case (this%type( bnd, i_dim ))
           case (p_bc_lindman)
             call update_e_2d(this%lindman(bnd), e)

           case ( p_bc_pec_odd )
             call pec_bc_e_2d_odd( e, i_dim, bnd )

           case ( p_bc_pmc_odd )
             call pmc_bc_e_2d_odd( e, i_dim, bnd )

           case ( p_bc_pec_even )
             call pec_bc_e_2d_even( e, i_dim, bnd )

           case ( p_bc_pmc_even)
             call pmc_bc_e_2d_even( e, i_dim, bnd )

           case default
         end select
       enddo

     enddo

     ! advance vpml B if necessary (cylindrical require all emf for radii correction)
     if ( is_active(this%vpml_all) ) then
         call update_b_2d( this%vpml_all, e )

     endif

   case (3)

     do i_dim = 1, 3
       do bnd = p_lower, p_upper
         select case (this%type( bnd, i_dim ))
           case (p_bc_lindman)
             call update_e_3d(this%lindman(bnd), e)

           case ( p_bc_pec_odd )
             call pec_bc_e_3d_odd( e, i_dim, bnd )

           case ( p_bc_pmc_odd )
             call pmc_bc_e_3d_odd( e, i_dim, bnd )

           case ( p_bc_pec_even )
             call pec_bc_e_3d_even( e, i_dim, bnd )

           case ( p_bc_pmc_even )
             call pmc_bc_e_3d_even( e, i_dim, bnd )

           case default
         end select
       enddo

     enddo

     ! advance vpml B if necessary

     if ( is_active(this%vpml_all) ) then
         call update_b_3d( this%vpml_all, e )

     endif

  end select

end subroutine update_boundary_efld_emf_bound
!-----------------------------------------------------------------------------------------

subroutine list_algorithm_emf_bound( this )

  use m_emf_define
  use m_parameters
  use m_emf_bound

  implicit none

  class( t_emf_bound ), intent(in) :: this

  integer :: i

  print *, ' '
  print *, '- EM Field Boundary condtions:'
  do i = 1, p_x_dim
    print '(A,I0,A,A,A,A,A)', '     x', i, ' : [', &
                trim(bnd_type_name( this%global_type( p_lower, i ))), ', ', &
                trim(bnd_type_name( this%global_type( p_upper, i ))), ']'
  enddo

end subroutine list_algorithm_emf_bound
