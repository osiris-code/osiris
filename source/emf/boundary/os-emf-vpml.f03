!
! Warning: if floating point exception are enabled then, under some compilers,

!          especially in when running in single precision, it is possible that this code
!          will raise an 'underflow' exception. This is to be expected as waves are being
!          damped in the PML region and reach very low amplitudes. When floating point
!          exceptions are not enabled, these values will be correctly rounded down to 0
!

#include "os-config.h"
#include "os-preprocess.fpp"

#ifdef __HAS_PML__

module m_vpml

#include "memory/memory.h"

  use m_restart
  use m_emf_define

  use m_space
  use m_node_conf

  use m_grid_define
  use m_grid

  use m_parameters

  use m_vdf_define
  use m_vdf_comm, only : t_vdf_msg

  use m_wall_define
  use m_wall_comm
  use m_wall

  implicit none

  private

  ! string to id restart data
  character(len=*), parameter :: p_vpml_rst_id = "vpml rst data - 0x0001"

  ! Attenuation exponent n
  ! Eq (38) of Vay's paper: JCP 165, 511 (2000)

  ! This value seems to be the most efficient
  real (p_k_fld), parameter :: p_att_nexp = 4.0_p_k_fld

  ! parameter for correction to avoid noise when particles reach boundary
  real (p_k_fld), parameter :: p_difrate = 0.1_p_k_fld

  interface setup
  module procedure setup_vpml
  end interface

  interface cleanup
    module procedure cleanup_vpml
  end interface

  interface move_window
    module procedure move_window_vpml
  end interface

  interface update_boundary
    module procedure update_boundary_vpml
  end interface

  interface restart_read
    module procedure restart_read_vpml
  end interface

  interface restart_write
    module procedure restart_write_vpml
  end interface

  interface is_active
    module procedure is_active_vpml
  end interface

  interface update_e_1d
    module procedure update_e_vpml_1d
  end interface

  interface update_e_2d
    module procedure update_e_vpml_2d
  end interface

  interface update_e_3d
    module procedure update_e_vpml_3d
  end interface

  interface update_b_1d
    module procedure update_b_vpml_1d
  end interface

  interface update_b_2d
    module procedure update_b_vpml_2d
  end interface

  interface update_b_3d
    module procedure update_b_vpml_3d
  end interface

  interface report
    module procedure report_vpml
  end interface

  interface reshape_obj
    module procedure reshape_vpml
  end interface

  interface update_beam
    module procedure update_beam_vpml
  end interface

  interface update_interface_x1_2d
    module procedure update_interface_x1_2d
  end interface

  interface update_interface_x2_2d
    module procedure update_interface_x2_2d
  end interface

  interface update_e_x1_2d
    module procedure update_e_x1_2d
  end interface

  interface update_cyl_e_x1_2d
    module procedure update_cyl_e_x1_2d
  end interface

  interface update_e_x2_2d
    module procedure update_e_x2_2d
  end interface

  interface update_b_x1_2d
    module procedure update_b_x1_2d
  end interface

  interface update_cyl_b_x1_2d
    module procedure update_cyl_b_x1_2d
  end interface

  interface update_b_x2_2d
    module procedure update_b_x2_2d
  end interface

  public :: t_vpml, setup, cleanup, is_active
  public :: update_e_1d, update_b_1d
  public :: update_e_2d, update_b_2d
  public :: update_e_3d, update_b_3d
  public :: report
  public :: restart_write, restart_read

  !      ! load balancing and new moving window routines
  public :: reshape_obj
  public :: move_window, update_boundary

  ! Update initial values from beam fields
  public :: update_beam

  public :: update_interface_x1_2d
  public :: update_interface_x2_2d
  public :: update_e_x1_2d, update_cyl_e_x1_2d
  public :: update_e_x2_2d
  public :: update_b_x1_2d, update_cyl_b_x1_2d
  public :: update_b_x2_2d

contains

!---------------------------------------------------
subroutine cleanup_vpml(this)
!---------------------------------------------------
   implicit none

   type (t_vpml), intent(inout) :: this

   integer :: i_vpml

   call freemem( this%dir_vpml )
   call freemem( this%loc_vpml )

   do i_vpml = 1, this%n_vpml
     call cleanup( this%wall_array_e(i_vpml) )
     call cleanup( this%wall_array_b(i_vpml) )
   enddo

  call freemem( this%wall_array_e )
  call freemem( this%wall_array_b )

  call freemem( this%pos_corner )

  call freemem( this%coef_e_low )
  call freemem( this%coef_b_low )
  call freemem( this%coef_e_up )
  call freemem( this%coef_b_up )

  this%n_vpml = 0

end subroutine cleanup_vpml
!---------------------------------------------------

!---------------------------------------------------
subroutine restart_write_vpml(this, restart_handle)
!---------------------------------------------------
  implicit none

  type (t_vpml), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  integer :: ierr, i

  restart_io_wr( p_vpml_rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then

    ERROR('error writing restart data for vpml object.')
    call abort_program(p_err_rstwrt)
  endif

  restart_io_wr( this%n_vpml, restart_handle, ierr )

  do i=1, this%n_vpml

    call restart_write(this%wall_array_e(i), restart_handle)
    call restart_write(this%wall_array_b(i), restart_handle)
  enddo

end subroutine restart_write_vpml
!---------------------------------------------------

!---------------------------------------------------
subroutine restart_read_vpml(this, restart_handle)
!---------------------------------------------------

   implicit none
   type (t_vpml), intent(inout) :: this
   type( t_restart_handle ), intent(in) :: restart_handle

   character(len=len(p_vpml_rst_id)) :: rst_id
   integer :: ierr

   ! check if restart file is compatible
   restart_io_rd( rst_id, restart_handle, ierr )
   if ( ierr/=0 ) then

     ERROR('error reading restart data for vpml object.')
     call abort_program(p_err_rstrd)
   endif

   if ( rst_id /= p_vpml_rst_id) then
     ERROR('Corrupted restart file ')
     ERROR('from incompatible binary (vpml)')
     call abort_program(p_err_rstrd)
   endif

   restart_io_rd( this%n_vpml, restart_handle, ierr )
   if ( ierr/=0 ) then

     ERROR('error reading restart data for vpml object.')
     call abort_program(p_err_rstrd)
   endif

end subroutine restart_read_vpml
!---------------------------------------------------

!---------------------------------------------------
!    setup vpml boundaries
!---------------------------------------------------
subroutine setup_vpml(this, e, b, grid, n_vpml, vpml_bnd_size, if_diffuse, if_vpml, if_move, dt, &
                      restart, restart_handle )

implicit none

! dummy variables

type (t_vpml), intent(inout) :: this
type ( t_vdf ), intent(in) :: e, b
class(t_grid), intent(in) :: grid
integer, intent(in) :: n_vpml
logical, dimension(:,:), intent(in) :: if_vpml
integer, intent(in) :: vpml_bnd_size

logical, dimension(:,:), intent(in) :: if_diffuse
logical, dimension(:), intent(in) :: if_move
real (p_double), intent(in) :: dt
logical, intent(in) :: restart
type( t_restart_handle ), intent(in) :: restart_handle

! local variables
integer :: j, i_vpml, i_dir, vpml_dims, adj_vpml
integer, dimension(2) :: range
integer, dimension(2, p_x_dim) :: nx_corner

  ! Read restart info
  ! - this currently only has the n_vpml value
  if ( restart ) then
    call restart_read( this, restart_handle )

    if ( this%n_vpml /= n_vpml ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "The number of vpml boundaries requested in the input deck ", &
                    " does not match the restart file."
        call abort_program(p_err_invalid)
      endif
    endif
  else
    this%n_vpml = n_vpml
  endif

  ! Store simulation coordinates
  this % coordinates = grid % coordinates
  if ( this % coordinates == p_cylindrical_b ) then
    ! Position of lower axial boundary
    this % gir_pos = grid%my_nx( p_lower, 2 )
  endif


  ! allocate vpml parameter arrays
  call alloc(this%dir_vpml, (/n_vpml/) )
  call alloc(this%loc_vpml, (/n_vpml/) )

  ! fill vpml parameter arrays
  this%dir_vpml = 0
  this%loc_vpml = 0
  this%bnd_size = 0

  i_vpml = 0
  do j=1, p_x_dim

    if ( if_vpml(p_lower,j) ) then
      i_vpml = i_vpml+1

      this%dir_vpml(i_vpml) = j
      this%loc_vpml(i_vpml) = p_lower
    endif

    if ( if_vpml(p_upper,j) ) then
      i_vpml = i_vpml+1

      this%dir_vpml(i_vpml) = j
      this%loc_vpml(i_vpml) = p_upper
    endif

  enddo

  ! initialize walls for each vpml boundary
  call alloc( this%wall_array_e, (/ n_vpml /) )
  call alloc( this%wall_array_b, (/ n_vpml /) )

  ! allocate corner position
  call alloc(this%pos_corner, (/ n_vpml, p_x_dim, 2 /) )

  ! PML wall thickness (same for all PML walls)
  this%bnd_size = vpml_bnd_size
  this%if_diffuse = if_diffuse

  this%pos_corner = -1
  do i_vpml = 1, n_vpml

    if ( this%loc_vpml(i_vpml) == p_lower ) then
      range(1) = -this%bnd_size + 2
    else
      range(1) = e%nx_(this%dir_vpml(i_vpml))+1
    endif

    range(2) = range(1) + this%bnd_size-1

    ! Setup corners
    ! since VPML has always 2 cells of overlap with the grid,
    ! NO corner => nx_corner = 2
    ! nx_corner is general for all directions (new_wall will take care of the rest)

    ! Corner priority: X1 > X2 > X3 (no corners in X3 walls)

    nx_corner = 2

    !  modify guard cell numbers if space moves
    do j = 1, p_x_dim
      if ( if_move(j) ) then
          nx_corner(p_lower,j) = nx_corner(p_lower,j) + 1
      endif
    enddo

    if (n_vpml > 1 .and. p_x_dim > 1) then

      ! X2 edges
      if ( this%dir_vpml(i_vpml) == 1 ) then

        ! bottom corner
        adj_vpml = vpml_exists(this, 2, p_lower)
        if ( adj_vpml > 0 ) then
          nx_corner(p_lower, 2) = this%bnd_size
          this%pos_corner( i_vpml, 2, p_lower) = adj_vpml
          this%pos_corner( adj_vpml, 1, this%loc_vpml(i_vpml) ) = i_vpml
        endif

        ! top corner
        adj_vpml = vpml_exists(this, 2, p_upper)
        if ( adj_vpml > 0 ) then

          nx_corner(p_upper, 2) = this%bnd_size
          this%pos_corner(i_vpml, 2, p_upper) = adj_vpml
          this%pos_corner( adj_vpml, 1, this%loc_vpml(i_vpml) ) = i_vpml

        endif

      endif

      ! Take care of X3 edges
      if (p_x_dim == 3) then

        i_dir = this%dir_vpml(i_vpml)

        if ( (i_dir == 1) .or. (i_dir == 2) ) then

          ! bottom

          adj_vpml = vpml_exists(this, p_x_dim, p_lower)
          if ( adj_vpml > 0 ) then
            nx_corner(p_lower, p_x_dim) = this%bnd_size
            this%pos_corner(i_vpml, p_x_dim, p_lower) = adj_vpml
            this%pos_corner( adj_vpml, i_dir, this%loc_vpml(i_vpml) ) = i_vpml
          endif

          ! top
          adj_vpml = vpml_exists(this, p_x_dim, p_upper)
          if ( adj_vpml > 0 ) then

            nx_corner(p_upper, p_x_dim) = this%bnd_size
            this%pos_corner(i_vpml, p_x_dim, p_upper) = adj_vpml
            this%pos_corner( adj_vpml, i_dir, this%loc_vpml(i_vpml) ) = i_vpml

          endif

        endif

      endif

    endif

    ! Dimensions of VPML fields
    ! Case could be avoided with vpml_dims = 3 + p_x_dim*(p_x_dim-1)/2
    select case ( p_x_dim )
      case (1)
        vpml_dims = 3  ! X, Y, Z

      case (2)
        vpml_dims = 4  ! X, Y, ZX, ZY

      case (3)
        vpml_dims = 6  ! XY, XZ, YX, YZ, ZX, ZY
    end select

    ! Create walls
    call new(this%wall_array_e(i_vpml), e, this%dir_vpml(i_vpml), &
            this%loc_vpml(i_vpml), range, &
            restart, restart_handle, &
            f_dim = vpml_dims, corner = nx_corner )
    call new(this%wall_array_b(i_vpml), b, this%dir_vpml(i_vpml), &
            this%loc_vpml(i_vpml), range, &
            restart, restart_handle, &
            f_dim = vpml_dims, corner = nx_corner )

  end do

  do j=1, p_x_dim
    this%dtdx(j) = real( dt/e%dx_(j), p_k_fld )
  enddo

  ! setup attenuation coefficients
  call setup_att_coeff(this)

end subroutine setup_vpml
!---------------------------------------------------

!---------------------------------------------------
subroutine move_window_vpml( this, space )
!---------------------------------------------------
  implicit none

  type (t_vpml), intent(inout) :: this
  type (t_space), intent(in) :: space ! local or global, only nx_move is needed

  ! local variables
  integer :: i_vpml

  ! only update boundaries if boundary is active
  do i_vpml=1, this%n_vpml
    call move_window(this%wall_array_e(i_vpml), space )
    call move_window(this%wall_array_b(i_vpml), space )
  enddo

end subroutine move_window_vpml
!---------------------------------------------------


!------------------------------------------------------------------------------
!> Updates the VPML wall objects on communication boundaries
!! @param   this            VPML object
!! @param   no_co           Parallel universe configuration
!! @param   nx_move         Number of cells to move in this iteration (used
!!                          when moving window is active)
!! @param   send_msg        Send message buffers to use
!! @param   recv_msg        Receive message buffers to use
subroutine update_boundary_vpml(this, no_co, nx_move, send_msg, recv_msg)

  implicit none

  type (t_vpml), intent(inout) :: this
  class( t_node_conf ),intent(in) :: no_co
  integer, dimension(:), intent(in):: nx_move
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

    integer :: i

    do i=1, this%n_vpml

        call update_boundary(this%wall_array_e(i), p_vdf_replace, no_co, nx_move, send_msg, recv_msg )
        call update_boundary(this%wall_array_b(i), p_vdf_replace, no_co, nx_move, send_msg, recv_msg )

    enddo

end subroutine update_boundary_vpml


!---------------------------------------------------
  function is_active_vpml(this)
!---------------------------------------------------
    implicit none
    logical :: is_active_vpml
    type (t_vpml), intent(in) :: this

    is_active_vpml = ( this%n_vpml > 0 )

    end function is_active_vpml
!---------------------------------------------------

!---------------------------------------------------
  function vpml_exists(this, dir, loc)
! checks if a neighbor vpml wall exists in the direction
! and location specified. Returns the vpml index.
!---------------------------------------------------
    implicit none
    integer :: vpml_exists
    type (t_vpml), intent(in) :: this
    integer, intent(in) :: dir
    integer, intent(in) :: loc

    integer :: i

    vpml_exists = -1
    do i=1, this%n_vpml

      if ( (this%dir_vpml(i) == dir) .and. (this%loc_vpml(i) == loc) ) then
        vpml_exists = i
      endif

    enddo

    end function vpml_exists
!---------------------------------------------------

!---------------------------------------------------
subroutine update_e_vpml_1d(this, b)
!---------------------------------------------------
  implicit none

  type (t_vpml), intent(inout) :: this
  type(t_vdf),intent(inout) :: b

  ! local variables
  integer :: i_vpml

  ! loop through vpml boundaries and call corresponding push
  do i_vpml=1, this%n_vpml

    call update_interface_x1_1d(this%wall_array_b(i_vpml), b, this%loc_vpml(i_vpml))

    ! Left wall
    if (this%loc_vpml(i_vpml) == p_lower) then

      call update_e_x1_1d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                          this%coef_e_low(1,:,:), p_lower )

    ! Right wall
    else

      call update_e_x1_1d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                             this%coef_e_up(1,:,:), p_upper )
    endif
  enddo

end subroutine update_e_vpml_1d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_e_vpml_2d(this, b)
!---------------------------------------------------
  implicit none

  type ( t_vpml ), intent(inout) :: this
  class( t_vdf ), intent(inout) :: b

  ! local variables
  integer, dimension(2) :: vpml_range, update_range
  integer :: i_vpml

  ! loop through vpml boundaries to update interfaces
   do i_vpml=1, this%n_vpml

    select case (this%dir_vpml(i_vpml))

      ! x1 direction
      case (1)

        call update_interface_x1_2d(this%wall_array_b, i_vpml, this%pos_corner(i_vpml,2,:), &
                                    b, this%loc_vpml(i_vpml) )

      ! x2 direction
      case (2)

        call update_interface_x2_2d(this%wall_array_b, i_vpml, &
                                    b, this%loc_vpml(i_vpml) )

    end select

  enddo

  ! loop through vpml boundaries to update fields
  do i_vpml=1, this%n_vpml

    ! E and B must have the same range
    vpml_range = range(this%wall_array_e(i_vpml))

    select case (this%dir_vpml(i_vpml))

      ! x1 direction
      case (1)

        ! Left wall
        if (this%loc_vpml(i_vpml) == p_lower) then
          update_range(1) = vpml_range(1)+1
          update_range(2) = vpml_range(2)-1

         select case ( this % coordinates )

           case default

            call update_e_x1_2d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                                this%dtdx(2), update_range, vpml_range(1), &
                                this%coef_e_low, this%coef_e_up, &
                                this%pos_corner(i_vpml,2,:), p_lower, &
                                this%if_diffuse(p_lower, 1))

           case ( p_cylindrical_b )

            call update_cyl_e_x1_2d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                                this%dtdx(2), update_range, vpml_range(1), &
                                this%coef_e_low, this%coef_e_up, &
                                this%pos_corner(i_vpml,2,:), p_lower, this % gir_pos)

         end select

        ! Right wall
        else
          update_range(1) = vpml_range(1)+1
          update_range(2) = vpml_range(2)

          select case ( this % coordinates )

            case default

              call update_e_x1_2d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                                  this%dtdx(2), update_range, vpml_range(1), &
                                  this%coef_e_low, this%coef_e_up, &
                                  this%pos_corner(i_vpml,2,:), p_upper, &
                                  this%if_diffuse(p_upper, 1) )

            case ( p_cylindrical_b )

              call update_cyl_e_x1_2d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                                  this%dtdx(2), update_range, vpml_range(1), &
                                  this%coef_e_low, this%coef_e_up, &
                                  this%pos_corner(i_vpml,2,:), p_upper, this % gir_pos )

          end select

        endif

      ! x2 direction

      ! for now, cylindrical use the standard solver: solver in cyl. coord differs by dr/2r factor
      ! which is always small, since only upper boundary pml is possible
      case (2)

        ! Bottom wall
        if (this%loc_vpml(i_vpml) == p_lower) then
          update_range(1) = vpml_range(1)+1
          update_range(2) = vpml_range(2)-1

          call update_e_x2_2d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                              this%coef_e_low(2,:,:), this%dtdx(1), &
                              update_range, vpml_range(1)-1, &
                              this%if_diffuse(p_lower, 2) )

        ! Top wall
        else
          update_range(1) = vpml_range(1)+1
          update_range(2) = vpml_range(2)

          call update_e_x2_2d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                              this%coef_e_up(2,:,:), this%dtdx(1), &
                              update_range, vpml_range(1), &
                              this%if_diffuse(p_upper, 2) )
        endif

    end select

  enddo

end subroutine update_e_vpml_2d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_e_vpml_3d(this, b)
!---------------------------------------------------
  implicit none

  type (t_vpml), intent(inout) :: this
  type(t_vdf),intent(inout) :: b

  ! local variables
  integer, dimension(2) :: vpml_range, update_range
  integer :: i_vpml

  ! loop through vpml boundaries to update interfaces
  do i_vpml=1, this%n_vpml

    ! E and B must have the same range
    vpml_range = range(this%wall_array_e(i_vpml))

    select case (this%dir_vpml(i_vpml))

      ! x1 direction
      case (1)

        call update_interface_x1_3d(this%wall_array_b, i_vpml, this%pos_corner, &
                                    b, this%loc_vpml(i_vpml) )

      ! x2 direction
      case (2)

        call update_interface_x2_3d(this%wall_array_b, i_vpml, this%pos_corner, &
                                    b, this%loc_vpml(i_vpml) )

      ! x3 direction
      case (3)

         call update_interface_x3_3d(this%wall_array_b, i_vpml, &
                                      b, this%loc_vpml(i_vpml) )

    end select

  enddo

  ! loop through vpml boundaries to update fields
  do i_vpml=1, this%n_vpml

    ! E and B must have the same range
    vpml_range = range(this%wall_array_e(i_vpml))

    select case (this%dir_vpml(i_vpml))

      ! x1 direction
      case (1)

        ! Left wall
        if (this%loc_vpml(i_vpml) == p_lower) then

          update_range(1) = vpml_range(1)+1
          update_range(2) = vpml_range(2)-1

          call update_e_x1_3d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                              this%dtdx, update_range, vpml_range(1), &
                              this%coef_e_low, this%coef_e_up, &
                              this%pos_corner(i_vpml,:,:), p_lower)

        ! Right wall
        else

          update_range(1) = vpml_range(1)+1
          update_range(2) = vpml_range(2)

          call update_e_x1_3d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                              this%dtdx, update_range, vpml_range(1), &
                              this%coef_e_low, this%coef_e_up, &
                              this%pos_corner(i_vpml,:,:), p_upper)
        endif

      ! x2 direction
      case (2)

        ! Front wall
        if (this%loc_vpml(i_vpml) == p_lower) then

          update_range(1) = vpml_range(1)+1
          update_range(2) = vpml_range(2)-1

          call update_e_x2_3d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                              this%dtdx, update_range, vpml_range(1), &
                              this%coef_e_low, this%coef_e_up, &
                              this%pos_corner(i_vpml,:,:), p_lower)

        ! Back wall
        else

          update_range(1) = vpml_range(1)+1
          update_range(2) = vpml_range(2)

          call update_e_x2_3d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                              this%dtdx, update_range, vpml_range(1), &
                              this%coef_e_low, this%coef_e_up, &
                              this%pos_corner(i_vpml,:,:), p_upper)

        endif

      ! x3 direction
      case (3)

        ! Bottom wall
        if (this%loc_vpml(i_vpml) == p_lower) then

          update_range(1) = vpml_range(1)+1
          update_range(2) = vpml_range(2)-1

          call update_e_x3_3d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                              this%coef_e_low(3,:,:), this%dtdx, &
                              update_range, vpml_range(1)-1 )

        ! Top wall
        else

          update_range(1) = vpml_range(1)+1
          update_range(2) = vpml_range(2)

          call update_e_x3_3d(this%wall_array_e(i_vpml), this%wall_array_b(i_vpml), &
                              this%coef_e_up(3,:,:), this%dtdx, &
                              update_range, vpml_range(1) )
        endif

    end select

  enddo

end subroutine update_e_vpml_3d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_e_x1_1d(vpml_e, vpml_b, coef_e, loc_vpml)
!
! Updates electric field in the left vay boundary
! Field components in 1D are
! 1 => X

! 2 => Y

! 3 => Z
!---------------------------------------------------
  implicit none

  type (t_wall), intent(inout) :: vpml_e
  type(t_wall), intent(in) :: vpml_b
  real(p_k_fld), dimension(:,:), intent(in) :: coef_e

  integer, intent(in)  :: loc_vpml

  ! local variables
  integer, dimension(2) :: upd_range, vpml_range
  integer :: i1, coef_i, coef_idx

  vpml_range = range(vpml_e)

  ! Left wall
  if (loc_vpml == p_lower) then

    upd_range(1) = vpml_range(1)+1
    upd_range(2) = vpml_range(2)-1

    coef_idx = vpml_range(1)-1

  ! Right wall
  else

    upd_range(1) = vpml_range(1)+1
    upd_range(2) = vpml_range(2)

    coef_idx = vpml_range(1)
  endif

  do i1=upd_range(1), upd_range(2)

    ! coefficients are always defined from 1:max_bnd_size
    coef_i = i1-coef_idx

    ! X component
    ! vpml_e%f1(1,i1) = vpml_e%f1(1,i1)

    ! Y component
    vpml_e%f1(2,i1) = coef_e(1,coef_i) *   vpml_e%f1(2,i1) &
                         - coef_e(2,coef_i) * vpml_b%f1(3,i1  ) &
                         + coef_e(3,coef_i) * vpml_b%f1(3,i1-1)

    ! Z component
    vpml_e%f1(3,i1) = coef_e(1,coef_i) *   vpml_e%f1(3,i1) &
                         + coef_e(2,coef_i) * vpml_b%f1(2,i1  ) &
                         - coef_e(3,coef_i) * vpml_b%f1(2,i1-1)

  enddo

end subroutine update_e_x1_1d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_e_x1_2d(vpml_e, vpml_b, dtdx, range, coef_idx, &
                          coef_e_low, coef_e_up, pos_corner, pos_vpml, if_diffuse )
!
! Updates electric field in the left vay boundary
! Field components in 2D are
! 1 => X (XY = XZ)
! 2 => Y (YX = YZ)
! 3 => ZX
! 4 => ZY
!---------------------------------------------------
  implicit none

  type (t_wall), intent(inout) :: vpml_e
  type(t_wall), intent(inout) :: vpml_b
  real(p_k_fld), intent(in) :: dtdx

  integer, dimension(:), intent(in) :: range
  integer, intent(in)  :: coef_idx
  integer, dimension(:), intent(in) :: pos_corner

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_e_low

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_e_up

  integer, intent(in) :: pos_vpml
  logical, intent(in) :: if_diffuse

  ! local variables
  integer :: gc1, gc2, i1, i2, coef_i, coef_std
  real(p_k_fld), dimension(:,:), pointer :: coef_e

  real(p_k_fld), dimension(3) :: coef_e_cnr

  logical :: if_lower_cnr, if_upper_cnr, do_corner
  real(p_k_fld) :: difrate_m, Bzx_temp

  gc1 = vpml_e%gc_num(p_lower,2)
  gc2 = vpml_e%gc_num(p_upper,2)

  ! select coefficients for invariant Field components
  if (pos_vpml == p_lower) then
    coef_e => coef_e_low(1,:,:)
    coef_std = coef_idx - 1  ! shift coefficients for lower boundary
  else
    coef_e => coef_e_up(1,:,:)
    coef_std = coef_idx
  endif

  ! check corner positions
  if_lower_cnr = (pos_corner(p_lower) > 0)
  if_upper_cnr = (pos_corner(p_upper) > 0)

  difrate_m = 1.0_p_k_fld - p_difrate

  ! loop vertically through wall
  do i2=2-gc1, vpml_e%nx(2)+gc2

    ! check if in a corner
    do_corner = .false.

    ! original
    if ( if_lower_cnr .and. (i2 <= 0) ) then

      coef_i = i2+gc1
      coef_e_cnr = coef_e_low(2,:,coef_i)
      do_corner = .true.
    else if ( if_upper_cnr .and. (i2 > vpml_e%nx(2)+1 ) ) then

      coef_i = i2 - vpml_e%nx(2)-1
      coef_e_cnr = coef_e_up(2,:,coef_i)
      do_corner = .true.
    endif

    ! loop horizontally through wall
    do i1=range(1), range(2)

      ! correction to avoid noise when particles reach boundary
      if (if_diffuse) then
        Bzx_temp = vpml_b%f2(3,i1,i2)
        vpml_b%f2(3,i1,i2) = difrate_m * vpml_b%f2(3,i1,i2) + p_difrate*vpml_b%f2(4,i1,i2)
        vpml_b%f2(4,i1,i2) = difrate_m * vpml_b%f2(4,i1,i2) + p_difrate*Bzx_temp
      endif

      ! Use equations for corner/edge
      if (do_corner) then
        ! X component
        vpml_e%f2(1,i1,i2) = coef_e_cnr(1) *   vpml_e%f2(1,i1,i2) + &
                             coef_e_cnr(2) * ( vpml_b%f2(3,i1,i2  ) + vpml_b%f2(4,i1,i2  ) ) - &
                             coef_e_cnr(3) * ( vpml_b%f2(3,i1,i2-1) + vpml_b%f2(4,i1,i2-1) )

        ! ZY component
        vpml_e%f2(4,i1,i2) = coef_e_cnr(1) * vpml_e%f2(4,i1,i2  ) - &
                             coef_e_cnr(2) * vpml_b%f2(1,i1,i2  ) + &
                             coef_e_cnr(3) * vpml_b%f2(1,i1,i2-1)

      ! standard solver for wall
      else
        ! X component
        vpml_e%f2(1,i1,i2) = vpml_e%f2(1,i1,i2) &
                             + dtdx*( vpml_b%f2(3,i1,i2)   + vpml_b%f2(4,i1,i2) &
                                    - vpml_b%f2(3,i1,i2-1) - vpml_b%f2(4,i1,i2-1) )

        ! ZY component
        vpml_e%f2(4,i1,i2) = vpml_e%f2(4,i1,i2) - dtdx*( vpml_b%f2(1,i1,i2) - vpml_b%f2(1,i1,i2-1) )

      endif

      ! Update invariant field components
      coef_i = i1 - coef_std

      ! Y component
      vpml_e%f2(2,i1,i2) = coef_e(1,coef_i) *   vpml_e%f2(2,i1,i2) &
                           - coef_e(2,coef_i) * ( vpml_b%f2(3,i1,i2)   + vpml_b%f2(4,i1,i2)   ) &
                           + coef_e(3,coef_i) * ( vpml_b%f2(3,i1-1,i2) + vpml_b%f2(4,i1-1,i2) )

      ! ZX component
      vpml_e%f2(3,i1,i2) = coef_e(1,coef_i) * vpml_e%f2(3,i1,i2) &
                           + coef_e(2,coef_i) * vpml_b%f2(2,i1,i2) &
                           - coef_e(3,coef_i) * vpml_b%f2(2,i1-1,i2)

    enddo
  enddo

end subroutine update_e_x1_2d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_cyl_e_x1_2d(vpml_e, vpml_b, dtdr, range, coef_idx, &
                              coef_e_low, coef_e_up, pos_corner, pos_vpml, gir_pos )
!
! Updates electric field in the left vay boundary
! Field components in 2D are
! 1 => X (XY = XZ) <= Only this component is different in Cyl (and not for the corner)
! 2 => Y (YX = YZ)
! 3 => ZX
! 4 => ZY
!---------------------------------------------------
  implicit none

  type (t_wall), intent(inout) :: vpml_e
  type(t_wall), intent(in) :: vpml_b
  real(p_k_fld), intent(in) :: dtdr

  integer, dimension(:), intent(in) :: range
  integer, intent(in)  :: coef_idx
  integer, dimension(:), intent(in) :: pos_corner

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_e_low

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_e_up

  integer, intent(in) :: pos_vpml
  integer, intent(in) :: gir_pos

  ! local variables
  integer :: gc1, gc2, i1, i2, coef_i, coef_std
  real(p_k_fld), dimension(:,:), pointer :: coef_e

  real(p_k_fld), dimension(3) :: coef_e_cnr

  logical :: if_upper_cnr, do_corner

  real(p_k_fld) :: tmp, rcp, rcm
  integer :: gshift_i2

  gc1 = vpml_e%gc_num(p_lower,2)
  gc2 = vpml_e%gc_num(p_upper,2)

  ! select coefficients for invariant Field components
  if (pos_vpml == p_lower) then
    coef_e => coef_e_low(1,:,:)
    coef_std = coef_idx - 1  ! shift coefficients for lower boundary
  else
    coef_e => coef_e_up(1,:,:)
    coef_std = coef_idx
  endif

  ! check corner positions
  if_upper_cnr = (pos_corner(p_upper) > 0)

  ! Get global cell shift
  gshift_i2 = gir_pos - 2

  ! loop vertically through wall
  do i2=2-gc1, vpml_e%nx(2)+gc2

    ! check if in a corner
    do_corner = .false.
    if ( if_upper_cnr .and. (i2 > vpml_e%nx(2)+1 ) ) then

      coef_i = i2 - vpml_e%nx(2)-1
      coef_e_cnr = coef_e_up(2,:,coef_i)
      do_corner = .true.
    endif

    tmp = dtdr / (( i2 + gshift_i2 ) - 0.5)  ! (1/r)*(dt/dr) at the corner of cell i2
    rcp = ( i2 + gshift_i2 )                 ! r at the center of cell i2
    rcm = ( i2 + gshift_i2 - 1 )             ! r at the center of cell i2-1

    ! loop horizontally through wall
    do i1=range(1), range(2)

      ! Use equations for corner/edge
      if (do_corner) then

        ! X component
        vpml_e%f2(1,i1,i2) = coef_e_cnr(1) *   vpml_e%f2(1,i1,i2) &
                             + coef_e_cnr(2) * ( vpml_b%f2(3,i1,i2  ) + vpml_b%f2(4,i1,i2  ) ) &
                             - coef_e_cnr(3) * ( vpml_b%f2(3,i1,i2-1) + vpml_b%f2(4,i1,i2-1) )

        ! ZY component
        vpml_e%f2(4,i1,i2) = coef_e_cnr(1) *   vpml_e%f2(4,i1,i2  ) &
                             - coef_e_cnr(2) * vpml_b%f2(1,i1,i2  ) &
                             + coef_e_cnr(3) * vpml_b%f2(1,i1,i2-1)

      ! standard solver for wall
      else
        ! X component
        vpml_e%f2(1,i1,i2) = vpml_e%f2(1,i1,i2) + tmp * &
                              ( rcp * (vpml_b%f2(3,i1,i2  ) + vpml_b%f2(4,i1,i2  )) - &
                                rcm * (vpml_b%f2(3,i1,i2-1) + vpml_b%f2(4,i1,i2-1)))

        ! ZY component
        vpml_e%f2(4,i1,i2) = vpml_e%f2(4,i1,i2) - dtdr*( vpml_b%f2(1,i1,i2) - vpml_b%f2(1,i1,i2-1) )

      endif

      ! Update invariant field components
      coef_i = i1 - coef_std

      ! Y component
      vpml_e%f2(2,i1,i2) = coef_e(1,coef_i) *   vpml_e%f2(2,i1,i2) &
                           - coef_e(2,coef_i) * ( vpml_b%f2(3,i1,i2)   + vpml_b%f2(4,i1,i2)   ) &
                           + coef_e(3,coef_i) * ( vpml_b%f2(3,i1-1,i2) + vpml_b%f2(4,i1-1,i2) )

      ! ZX component
      vpml_e%f2(3,i1,i2) = coef_e(1,coef_i) * vpml_e%f2(3,i1,i2) &
                           + coef_e(2,coef_i) * vpml_b%f2(2,i1,i2) &
                           - coef_e(3,coef_i) * vpml_b%f2(2,i1-1,i2)

    enddo
  enddo

end subroutine update_cyl_e_x1_2d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_e_x1_3d(vpml_e, vpml_b, dtdx, range, coef_idx, &
                          coef_e_low, coef_e_up, pos_corner, pos_vpml )
!
! Updates electric field in the left vay boundary
! Field components in 2D are
! 1 => XY
! 2 => XZ
!
! 3 => YX
! 4 => YZ
!
! 5 => ZX
! 6 => ZY
!---------------------------------------------------
  implicit none

  type (t_wall), intent(inout) :: vpml_e
  type(t_wall), intent(in) :: vpml_b
  real(p_k_fld), dimension(:), intent(in) :: dtdx
  integer, dimension(:), intent(in) :: range
  integer, intent(in)  :: coef_idx
  integer, dimension(:,:), intent(in) :: pos_corner

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_e_low

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_e_up

  integer, intent(in) :: pos_vpml

  ! local variables
  integer :: gc2a, gc2b, gc3a, gc3b, i1, i2, i3, coef_i, coef_std
  real(p_k_fld), dimension(:,:), pointer :: coef_e

  real(p_k_fld), dimension(3) :: coef_e_cnr_x2, coef_e_cnr_x3

  logical, dimension(2) :: if_cnr_x2, if_cnr_x3
  logical :: do_corner_x2, do_corner_x3

  integer :: rank = 3

  gc2a = vpml_e%gc_num(p_lower,2)
  gc2b = vpml_e%gc_num(p_upper,2)

  gc3a = vpml_e%gc_num(p_lower,rank)
  gc3b = vpml_e%gc_num(p_upper,rank)

  ! select coefficients for invariant Field components
  if (pos_vpml == p_lower) then
    coef_e => coef_e_low(1,:,:)
    coef_std = coef_idx - 1  ! shift coefficients for lower boundary
  else
    coef_e => coef_e_up(1,:,:)
    coef_std = coef_idx
  endif

  ! check corner positions
  ! SFM: should these be global variables of the VPML?
  if_cnr_x2(p_lower) = (pos_corner(2,p_lower) > 0)
  if_cnr_x2(p_upper) = (pos_corner(2,p_upper) > 0)

  if_cnr_x3(p_lower) = (pos_corner(rank,p_lower) > 0)
  if_cnr_x3(p_upper) = (pos_corner(rank,p_upper) > 0)

  do i3=2-gc3a, vpml_e%nx(rank)+gc3b

    ! check if in a X3 edge
    do_corner_x3 = .false.
    if ( if_cnr_x3(p_lower) .and. (i3 <= 0) ) then

      coef_i = i3+gc3a
      coef_e_cnr_x3 = coef_e_low(rank,:,coef_i)
      do_corner_x3 = .true.
    else if ( if_cnr_x3(p_upper) .and. (i3 > vpml_e%nx(rank)+1 ) ) then

      coef_i = i3 - vpml_e%nx(rank)-1
      coef_e_cnr_x3 = coef_e_up(rank,:,coef_i)
      do_corner_x3 = .true.
    endif

    do i2=2-gc2a, vpml_e%nx(2)+gc2b

      ! check if in a X2 edge
      do_corner_x2 = .false.
      if ( if_cnr_x2(p_lower) .and. (i2 <= 0) ) then

        coef_i = i2+gc2a
        coef_e_cnr_x2 = coef_e_low(2,:,coef_i)
        do_corner_x2 = .true.
      else if ( if_cnr_x2(p_upper) .and. (i2 > vpml_e%nx(2)+1 ) ) then

        coef_i = i2 - vpml_e%nx(2) -1
        coef_e_cnr_x2 = coef_e_up(2,:,coef_i)
        do_corner_x2 = .true.
      endif

      do i1=range(1), range(2)

        ! Use equations for corner X2
        if (do_corner_x2) then

         ! XY component
         vpml_e%f3(1,i1,i2,i3) = coef_e_cnr_x2(1) * vpml_e%f3(1, i1, i2, i3) &
                                 + coef_e_cnr_x2(2) * ( vpml_b%f3(5, i1, i2  , i3) &
                                                      + vpml_b%f3(6, i1, i2  , i3) ) &
                                 - coef_e_cnr_x2(3) * ( vpml_b%f3(5, i1, i2-1, i3) &
                                                      + vpml_b%f3(6, i1, i2-1, i3) )

         ! ZY component
         vpml_e%f3(6,i1,i2,i3) = coef_e_cnr_x2(1) * vpml_e%f3(6, i1, i2, i3) &
                                 - coef_e_cnr_x2(2) * ( vpml_b%f3(1, i1, i2  , i3) &
                                                      + vpml_b%f3(2, i1, i2  , i3) ) &
                                 + coef_e_cnr_x2(3) * ( vpml_b%f3(1, i1, i2-1, i3) &
                                                      + vpml_b%f3(2, i1, i2-1, i3) )

        ! standard solver for wall

        else

         ! XY component
         vpml_e%f3(1,i1,i2,i3) = vpml_e%f3(1, i1, i2, i3) &
                                 + dtdx(2)*( vpml_b%f3(5, i1, i2  , i3) &
                                           + vpml_b%f3(6, i1, i2  , i3) &
                                           - vpml_b%f3(5, i1, i2-1, i3) &
                                           - vpml_b%f3(6, i1, i2-1, i3) )

         ! ZY component => standard solver
         vpml_e%f3(6,i1,i2,i3) = vpml_e%f3(6, i1, i2, i3) &
                                 - dtdx(2)*( vpml_b%f3(1, i1, i2  , i3) &
                                           + vpml_b%f3(2, i1, i2  , i3) &
                                           - vpml_b%f3(1, i1, i2-1, i3) &
                                           - vpml_b%f3(2, i1, i2-1, i3) )

        endif

        ! Use equations for corner X3
        if (do_corner_x3) then

          ! XZ component
         vpml_e%f3(2,i1,i2,i3) = coef_e_cnr_x3(1) * vpml_e%f3(2, i1, i2, i3) &
                                 - coef_e_cnr_x3(2) * ( vpml_b%f3(3, i1, i2, i3  ) &
                                                      + vpml_b%f3(4, i1, i2, i3  ) ) &
                                 + coef_e_cnr_x3(3) * ( vpml_b%f3(3, i1, i2, i3-1) &
                                                      + vpml_b%f3(4, i1, i2, i3-1) )

         ! YZ component
         vpml_e%f3(4,i1,i2,i3) = coef_e_cnr_x3(1) * vpml_e%f3(4, i1, i2, i3) &
                                 + coef_e_cnr_x3(2) * ( vpml_b%f3(1, i1, i2, i3  ) &
                                                      + vpml_b%f3(2, i1, i2, i3  ) ) &
                                 - coef_e_cnr_x3(3) * ( vpml_b%f3(1, i1, i2, i3-1) &
                                                      + vpml_b%f3(2, i1, i2, i3-1) )

        ! standard solver for wall
        else

         ! XZ component
         vpml_e%f3(2,i1,i2,i3) = vpml_e%f3(2, i1, i2, i3) &
                                 - dtdx(rank)*( vpml_b%f3(3, i1, i2, i3  ) &
                                           + vpml_b%f3(4, i1, i2, i3  ) &
                                           - vpml_b%f3(3, i1, i2, i3-1) &
                                           - vpml_b%f3(4, i1, i2, i3-1) )

         ! YZ component
         vpml_e%f3(4,i1,i2,i3) = vpml_e%f3(4, i1, i2, i3) &
                                 + dtdx(rank)*( vpml_b%f3(1, i1, i2, i3  ) &
                                           + vpml_b%f3(2, i1, i2, i3  ) &
                                           - vpml_b%f3(1, i1, i2, i3-1) &
                                           - vpml_b%f3(2, i1, i2, i3-1) )

        endif

        ! Update invariant field components
        coef_i = i1-coef_std

        ! YX component

        vpml_e%f3(3,i1,i2,i3) = coef_e(1,coef_i) * vpml_e%f3(3, i1, i2, i3) &
                                - coef_e(2,coef_i) * ( vpml_b%f3(5, i1  , i2, i3 ) &
                                                     + vpml_b%f3(6, i1  , i2, i3 ) ) &
                                + coef_e(3,coef_i) * ( vpml_b%f3(5, i1-1, i2, i3 ) &
                                                     + vpml_b%f3(6, i1-1, i2, i3 ) )

        ! ZX component

        vpml_e%f3(5,i1,i2,i3) = coef_e(1,coef_i) * vpml_e%f3(5, i1, i2, i3) &
                                + coef_e(2,coef_i) * ( vpml_b%f3(3, i1  , i2, i3 ) &
                                                     + vpml_b%f3(4, i1  , i2, i3 ) ) &
                                - coef_e(3,coef_i) * ( vpml_b%f3(3, i1-1, i2, i3 ) &
                                                     + vpml_b%f3(4, i1-1, i2, i3 ) )

      enddo
    enddo
  enddo

end subroutine update_e_x1_3d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_e_x2_2d(vpml_e, vpml_b, coef_e, dtdx, range, coef_idx, if_diffuse)
!
! Updates electric field in the left vay boundary
! Field components in 2D are
! 1 => X (XY = XZ)
! 2 => Y (YX = YZ)
! 3 => ZX
! 4 => ZY
!---------------------------------------------------
  implicit none

  type (t_wall), intent(inout) :: vpml_e
  type(t_wall), intent(inout) :: vpml_b
  real(p_k_fld), dimension(:,:), intent(in) :: coef_e

  real(p_k_fld), intent(in) :: dtdx

  integer, dimension(:), intent(in) :: range
  integer, intent(in)  :: coef_idx
  logical, intent(in) :: if_diffuse

  ! local variables
  integer :: gc1, gc2, i1, i2, coef_i
  real(p_k_fld) :: difrate_m, Bzx_temp

 !  logical :: diffuse

  gc1 = vpml_e%gc_num(p_lower,1)
  gc2 = vpml_e%gc_num(p_upper,1)

  if (if_diffuse) then
    do i2= range(1), range(2)
      do i1=2-gc1, vpml_e%nx(1)+gc2
        difrate_m = 1.0_p_k_fld - p_difrate
        ! correction to avoid noise when particles reach boundary
        Bzx_temp = vpml_b%f2(3,i1,i2)
        vpml_b%f2(3,i1,i2) = difrate_m * vpml_b%f2(3,i1,i2) + p_difrate*vpml_b%f2(4,i1,i2)
        vpml_b%f2(4,i1,i2) = difrate_m * vpml_b%f2(4,i1,i2) + p_difrate*Bzx_temp
      enddo
    enddo
  endif

  ! On the lower boundary also do the lowest row for components Y and ZX
  if ( range(1) > lbound( vpml_b%f2, 3 ) ) then
    i2 = range(1) - 1

    do i1=2-gc1, vpml_e%nx(1)+gc2

      ! Y component => standard solver
      vpml_e%f2(2,i1,i2) = vpml_e%f2(2,i1,i2) &
                           - dtdx*( vpml_b%f2(3,i1  ,i2) + vpml_b%f2(4,i1  ,i2) &
                                  - vpml_b%f2(3,i1-1,i2) - vpml_b%f2(4,i1-1,i2) )

      ! ZX component => standard solver
      vpml_e%f2(3,i1,i2) = vpml_e%f2(3,i1,i2) + dtdx*( vpml_b%f2(2,i1,i2) - vpml_b%f2(2,i1-1,i2) )

    enddo

  endif

  do i2= range(1), range(2)
    ! coefficients are always defined from 1:max_bnd_size
    coef_i = i2-coef_idx

    do i1=2-gc1, vpml_e%nx(1)+gc2

      ! X component
      vpml_e%f2(1,i1,i2) = coef_e(1,coef_i) *   vpml_e%f2(1,i1,i2) + &
                           coef_e(2,coef_i) * ( vpml_b%f2(3,i1,i2  ) + vpml_b%f2(4,i1,i2  ) ) - &
                           coef_e(3,coef_i) * ( vpml_b%f2(3,i1,i2-1) + vpml_b%f2(4,i1,i2-1) )

      ! Y component => standard solver
      vpml_e%f2(2,i1,i2) = vpml_e%f2(2,i1,i2) &
                           - dtdx*( vpml_b%f2(3,i1  ,i2) + vpml_b%f2(4,i1  ,i2) &
                                  - vpml_b%f2(3,i1-1,i2) - vpml_b%f2(4,i1-1,i2) )

      ! ZX component => standard solver
      vpml_e%f2(3,i1,i2) = vpml_e%f2(3,i1,i2) + dtdx*( vpml_b%f2(2,i1,i2) - vpml_b%f2(2,i1-1,i2) )

      ! ZY component
      vpml_e%f2(4,i1,i2) = coef_e(1,coef_i) * vpml_e%f2(4,i1,i2  ) - &
                           coef_e(2,coef_i) * vpml_b%f2(1,i1,i2  ) + &
                           coef_e(3,coef_i) * vpml_b%f2(1,i1,i2-1)

    enddo
  enddo

end subroutine update_e_x2_2d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_e_x2_3d(vpml_e, vpml_b, dtdx, range, coef_idx, &
                          coef_e_low, coef_e_up, pos_corner, pos_vpml )
!
! Updates electric field in the left vay boundary
! Field components in 2D are
! 1 => XY
! 2 => XZ
!
! 3 => YX
! 4 => YZ
!
! 5 => ZX
! 6 => ZY
!---------------------------------------------------
  implicit none

  type (t_wall), intent(inout) :: vpml_e
  type(t_wall), intent(in) :: vpml_b
  real(p_k_fld), dimension(:), intent(in) :: dtdx
  integer, dimension(:), intent(in) :: range
  integer, intent(in)  :: coef_idx
  integer, dimension(:,:), intent(in) :: pos_corner

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_e_low

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_e_up

  integer, intent(in) :: pos_vpml

  ! local variables
  integer :: gc1a, gc1b, gc3a, gc3b, i1, i2, i3, coef_i, coef_std
  real(p_k_fld), dimension(:,:), pointer :: coef_e

  real(p_k_fld), dimension(3) :: coef_e_cnr_x3

  logical, dimension(2) :: if_cnr_x3
  logical :: do_corner_x3

  integer :: pi3

  pi3 = 3

  gc1a = vpml_e%gc_num(p_lower,1)
  gc1b = vpml_e%gc_num(p_upper,1)

  gc3a = vpml_e%gc_num(p_lower,3)
  gc3b = vpml_e%gc_num(p_upper,3)

  ! select coefficients for invariant Field components
  if (pos_vpml == p_lower) then
    coef_e => coef_e_low(2,:,:)
    coef_std = coef_idx - 1  ! shift coefficients for lower boundary
  else
    coef_e => coef_e_up(2,:,:)
    coef_std = coef_idx
  endif

  ! check corner positions
  ! SFM: should these be global variables of the VPML?
  if_cnr_x3(p_lower) = (pos_corner(pi3,p_lower) > 0)
  if_cnr_x3(p_upper) = (pos_corner(pi3,p_upper) > 0)

  do i3=2-gc3a, vpml_e%nx(3)+gc3b

    ! check if in a X3 edge
    do_corner_x3 = .false.
    if ( if_cnr_x3(p_lower) .and. (i3 <= 0) ) then

      coef_i = i3+gc3a
      coef_e_cnr_x3 = coef_e_low(3,:,coef_i)
      do_corner_x3 = .true.
    else if ( if_cnr_x3(p_upper) .and. (i3 > vpml_e%nx(3)+1 ) ) then

      coef_i = i3 - vpml_e%nx(3)-1
      coef_e_cnr_x3 = coef_e_up(3,:,coef_i)
      do_corner_x3 = .true.
    endif

    do i2=range(1), range(2)

      ! coefficient index for invariant field components
      coef_i = i2-coef_std

      do i1=2-gc1a, vpml_e%nx(1)+gc1b

        ! Use equations for corner X3
        if (do_corner_x3) then

          ! XZ component
          vpml_e%f3(2,i1,i2,i3) = coef_e_cnr_x3(1) * vpml_e%f3(2, i1, i2, i3) &
                                  - coef_e_cnr_x3(2) * ( vpml_b%f3(3, i1, i2, i3  ) &
                                                       + vpml_b%f3(4, i1, i2, i3  ) ) &
                                  + coef_e_cnr_x3(3) * ( vpml_b%f3(3, i1, i2, i3-1) &
                                                       + vpml_b%f3(4, i1, i2, i3-1) )

          ! YZ component
          vpml_e%f3(4,i1,i2,i3) = coef_e_cnr_x3(1) * vpml_e%f3(4, i1, i2, i3) &
                                  + coef_e_cnr_x3(2) * ( vpml_b%f3(1, i1, i2, i3  ) &
                                                       + vpml_b%f3(2, i1, i2, i3  ) ) &
                                  - coef_e_cnr_x3(3) * ( vpml_b%f3(1, i1, i2, i3-1) &
                                                       + vpml_b%f3(2, i1, i2, i3-1) )

        ! standard solver for wall
        else

          ! XZ component
          vpml_e%f3(2,i1,i2,i3) = vpml_e%f3(2, i1, i2, i3) &
                                  - dtdx(pi3)*( vpml_b%f3(3, i1, i2, i3  ) &
                                            + vpml_b%f3(4, i1, i2, i3  ) &
                                            - vpml_b%f3(3, i1, i2, i3-1) &
                                            - vpml_b%f3(4, i1, i2, i3-1) )

          ! YZ component
          vpml_e%f3(4,i1,i2,i3) = vpml_e%f3(4, i1, i2, i3) &
                                  + dtdx(pi3)*( vpml_b%f3(1, i1, i2, i3  ) &
                                            + vpml_b%f3(2, i1, i2, i3  ) &
                                            - vpml_b%f3(1, i1, i2, i3-1) &
                                            - vpml_b%f3(2, i1, i2, i3-1) )

        endif

        ! Update invariant field components

        ! XY component
        vpml_e%f3(1,i1,i2,i3) = coef_e(1,coef_i) * vpml_e%f3(1, i1, i2, i3) &
                                + coef_e(2,coef_i) * ( vpml_b%f3(5, i1, i2  , i3 ) &
                                                     + vpml_b%f3(6, i1, i2  , i3 ) ) &
                                - coef_e(3,coef_i) * ( vpml_b%f3(5, i1, i2-1, i3 ) &
                                                     + vpml_b%f3(6, i1, i2-1, i3 ) )

        ! YX component => standard solver
        vpml_e%f3(3,i1,i2,i3) = vpml_e%f3(3, i1, i2, i3) &
                                - dtdx(1)*( vpml_b%f3(5, i1  , i2, i3) &
                                          + vpml_b%f3(6, i1  , i2, i3) &
                                          - vpml_b%f3(5, i1-1, i2, i3) &
                                          - vpml_b%f3(6, i1-1, i2, i3) )

        ! ZX component => standard solver
        vpml_e%f3(5,i1,i2,i3) = vpml_e%f3(5, i1, i2, i3) &
                                + dtdx(1)*( vpml_b%f3(3, i1  , i2, i3) &
                                          + vpml_b%f3(4, i1  , i2, i3) &
                                          - vpml_b%f3(3, i1-1, i2, i3) &
                                          - vpml_b%f3(4, i1-1, i2, i3) )

        ! ZY component
        vpml_e%f3(6,i1,i2,i3) =   coef_e(1,coef_i) * vpml_e%f3(6, i1, i2, i3) &
                                - coef_e(2,coef_i) * ( vpml_b%f3(1, i1, i2  , i3 ) &
                                                     + vpml_b%f3(2, i1, i2  , i3 ) ) &
                                + coef_e(3,coef_i) * ( vpml_b%f3(1, i1, i2-1, i3 ) &
                                                     + vpml_b%f3(2, i1, i2-1, i3 ) )

      enddo
    enddo
  enddo

end subroutine update_e_x2_3d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_e_x3_3d(vpml_e, vpml_b, coef_e, dtdx, range, coef_idx)
!
! Updates electric field in the left vay boundary
! Field components in 3D are
! 1 => XY
! 2 => XZ
!
! 3 => YX
! 4 => YZ
!
! 5 => ZX
! 6 => ZY
!---------------------------------------------------
  implicit none

  type (t_wall), intent(inout) :: vpml_e
  type(t_wall), intent(in) :: vpml_b
  real(p_k_fld), dimension(:,:), intent(in) :: coef_e

  real(p_k_fld), dimension(:), intent(in) :: dtdx
  integer, dimension(:), intent(in) :: range
  integer, intent(in)  :: coef_idx

  ! local variables
  integer :: gc1a, gc1b, gc2a, gc2b, i1, i2, i3, coef_i

  gc1a = vpml_e%gc_num(p_lower,1)
  gc1b = vpml_e%gc_num(p_upper,1)

  gc2a = vpml_e%gc_num(p_lower,2)
  gc2b = vpml_e%gc_num(p_upper,2)

  do i3=range(1), range(2)

    coef_i = i3-coef_idx

    do i2=2-gc2a, vpml_e%nx(2)+gc2b
      do i1=2-gc1a, vpml_e%nx(1)+gc1b

        ! XY component => standard solver
        vpml_e%f3(1,i1,i2,i3) = vpml_e%f3(1, i1, i2, i3) &
                                + dtdx(2)*( vpml_b%f3(5, i1, i2  , i3  ) &
                                          + vpml_b%f3(6, i1, i2  , i3  ) &
                                          - vpml_b%f3(5, i1, i2-1, i3) &
                                          - vpml_b%f3(6, i1, i2-1, i3) )

        ! XZ component
        vpml_e%f3(2,i1,i2,i3) = coef_e(1,coef_i) * vpml_e%f3(2, i1, i2, i3) &
                                - coef_e(2,coef_i) * ( vpml_b%f3(3, i1, i2, i3   ) &
                                                     + vpml_b%f3(4, i1, i2, i3   ) ) &
                                + coef_e(3,coef_i) * ( vpml_b%f3(3, i1, i2, i3-1 ) &
                                                     + vpml_b%f3(4, i1, i2, i3-1 ) )

        ! YX component => standard solver
        vpml_e%f3(3,i1,i2,i3) = vpml_e%f3(3, i1, i2, i3) &
                                - dtdx(1)*( vpml_b%f3(5, i1  , i2, i3) &
                                          + vpml_b%f3(6, i1  , i2, i3) &
                                          - vpml_b%f3(5, i1-1, i2, i3) &
                                          - vpml_b%f3(6, i1-1, i2, i3) )

        ! YZ component
        vpml_e%f3(4,i1,i2,i3) = coef_e(1,coef_i) * vpml_e%f3(4, i1, i2, i3) &
                                + coef_e(2,coef_i) * ( vpml_b%f3(1, i1, i2, i3   ) &
                                                     + vpml_b%f3(2, i1, i2, i3   ) ) &
                                - coef_e(3,coef_i) * ( vpml_b%f3(1, i1, i2, i3-1 ) &
                                                     + vpml_b%f3(2, i1, i2, i3-1 ) )

        ! ZX component => standard solver
        vpml_e%f3(5,i1,i2,i3) = vpml_e%f3(5, i1, i2, i3) &
                                + dtdx(1)*( vpml_b%f3(3, i1  , i2, i3) &
                                          + vpml_b%f3(4, i1  , i2, i3) &
                                          - vpml_b%f3(3, i1-1, i2, i3) &
                                          - vpml_b%f3(4, i1-1, i2, i3) )

        ! ZY component => standard solver
        vpml_e%f3(6,i1,i2,i3) = vpml_e%f3(6, i1, i2, i3) &
                                - dtdx(2)*( vpml_b%f3(1, i1, i2  , i3  ) &
                                          + vpml_b%f3(2, i1, i2  , i3  ) &
                                          - vpml_b%f3(1, i1, i2-1, i3) &
                                          - vpml_b%f3(2, i1, i2-1, i3) )

      enddo
    enddo
  enddo

end subroutine update_e_x3_3d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_b_vpml_1d(this, e)
!---------------------------------------------------
  implicit none

  type (t_vpml), intent(inout) :: this
  type(t_vdf),intent(inout) :: e

  ! local variables
  integer :: i_vpml

  ! loop through vpml boundaries and call corresponding push
  do i_vpml=1, this%n_vpml

    call update_interface_x1_1d(this%wall_array_e(i_vpml), e, this%loc_vpml(i_vpml))

    ! Left wall
    if (this%loc_vpml(i_vpml) == p_lower) then

      call update_b_x1_1d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                             this%coef_b_low(1,:,:), p_lower )

    ! Right wall
    else

      call update_b_x1_1d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                             this%coef_b_up(1,:,:), p_upper )

    endif
  enddo

end subroutine update_b_vpml_1d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_b_vpml_2d(this, e)
!---------------------------------------------------
  implicit none

  type (t_vpml), intent(inout) :: this
  class( t_vdf ), intent(inout) :: e

  ! local variables
  integer, dimension(2) :: vpml_range, update_range
  integer :: i_vpml

  ! loop through vpml boundaries and call corresponding push
  do i_vpml=1, this%n_vpml

    select case (this%dir_vpml(i_vpml))

      ! x1 direction
      case (1)

        call update_interface_x1_2d(this%wall_array_e, i_vpml, this%pos_corner(i_vpml,2,:), &
                                    e, this%loc_vpml(i_vpml) )

      ! x2 direction
      case (2)

        call update_interface_x2_2d(this%wall_array_e, i_vpml, &
                                    e, this%loc_vpml(i_vpml) )

    end select

  enddo

  ! loop through vpml boundaries and call corresponding push
  do i_vpml=1, this%n_vpml

    ! E and B must have the same range
    vpml_range = range(this%wall_array_e(i_vpml))

    select case (this%dir_vpml(i_vpml))

      ! x1 direction
      case (1)

        ! Left wall
        if (this%loc_vpml(i_vpml) == p_lower) then
          update_range(1) = vpml_range(1)
          update_range(2) = vpml_range(2)-1

         select case ( this%coordinates )

           case default

             call update_b_x1_2d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                                 this%dtdx(2), update_range, vpml_range(1), &

                                 this%coef_b_low, this%coef_b_up, &
                                 this%pos_corner(i_vpml,2,:), p_lower )

           case ( p_cylindrical_b )

             call update_cyl_b_x1_2d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                                 this%dtdx(2), update_range, vpml_range(1), &
                                 this%coef_b_low, this%coef_b_up, &
                                 this%pos_corner(i_vpml,2,:), p_lower, this % gir_pos )

           end select

        ! Right wall
        else
          update_range(1) = vpml_range(1)+1
          update_range(2) = vpml_range(2)-1

         select case ( this%coordinates )

           case default

             call update_b_x1_2d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                                    this%dtdx(2), update_range, vpml_range(1), &
                                    this%coef_b_low, this%coef_b_up, &
                                    this%pos_corner(i_vpml,2,:), p_upper )

           case ( p_cylindrical_b )

             call update_cyl_b_x1_2d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                                    this%dtdx(2), update_range, vpml_range(1), &
                                    this%coef_b_low, this%coef_b_up, &
                                    this%pos_corner(i_vpml,2,:), p_upper, this % gir_pos )

            end select

        endif

      ! x2 direction
      case (2)

        ! Bottom wall
        if (this%loc_vpml(i_vpml) == p_lower) then
          update_range(1) = vpml_range(1)
          update_range(2) = vpml_range(2)-1

          call update_b_x2_2d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                                 this%coef_b_low(2,:,:), this%dtdx(1), &
                                 update_range, vpml_range(1)-2 )

        ! Top wall
        else
          update_range(1) = vpml_range(1)+1
          update_range(2) = vpml_range(2)-1

          call update_b_x2_2d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                                 this%coef_b_up(2,:,:), this%dtdx(1), &
                                 update_range, vpml_range(1) )

        endif

    end select

  enddo

end subroutine update_b_vpml_2d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_b_vpml_3d(this, e)
!---------------------------------------------------
  implicit none

  type (t_vpml), intent(inout) :: this
  type(t_vdf),intent(inout) :: e

  ! local variables
  integer, dimension(2) :: vpml_range, update_range
  integer :: i_vpml

  ! loop through vpml boundaries to update interfaces
  do i_vpml=1, this%n_vpml

    ! E and B must have the same range
    vpml_range = range(this%wall_array_e(i_vpml))

    select case (this%dir_vpml(i_vpml))

      ! x1 direction
      case (1)

        call update_interface_x1_3d(this%wall_array_e, i_vpml, this%pos_corner, &
                                    e, this%loc_vpml(i_vpml) )

      ! x2 direction
      case (2)

        call update_interface_x2_3d(this%wall_array_e, i_vpml, this%pos_corner, &
                                    e, this%loc_vpml(i_vpml) )

      ! x3 direction
      case (3)

        call update_interface_x3_3d(this%wall_array_e, i_vpml, &
                                      e, this%loc_vpml(i_vpml))

    end select

  enddo

  ! loop through vpml boundaries to update fields
  do i_vpml=1, this%n_vpml

    ! E and B must have the same range
    vpml_range = range(this%wall_array_e(i_vpml))

    select case (this%dir_vpml(i_vpml))

      ! x1 direction
      case (1)

        ! Left wall
        if (this%loc_vpml(i_vpml) == p_lower) then

          update_range(1) = vpml_range(1)
          update_range(2) = vpml_range(2)-1

          call update_b_x1_3d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                              this%dtdx, update_range, vpml_range(1), &
                              this%coef_b_low, this%coef_b_up, &
                              this%pos_corner(i_vpml,:,:), p_lower)

        ! Right wall
        else

          update_range(1) = vpml_range(1)+1
          update_range(2) = vpml_range(2)-1

          call update_b_x1_3d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                              this%dtdx, update_range, vpml_range(1), &
                              this%coef_b_low, this%coef_b_up, &
                              this%pos_corner(i_vpml,:,:), p_upper)

        endif

      ! x2 direction
      case (2)

        ! Front wall
        if (this%loc_vpml(i_vpml) == p_lower) then

          update_range(1) = vpml_range(1)
          update_range(2) = vpml_range(2)-1

          call update_b_x2_3d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                              this%dtdx, update_range, vpml_range(1), &
                              this%coef_b_low, this%coef_b_up, &
                              this%pos_corner(i_vpml,:,:), p_lower)

        ! Back wall
        else

          update_range(1) = vpml_range(1)+1
          update_range(2) = vpml_range(2)-1

          call update_b_x2_3d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                              this%dtdx, update_range, vpml_range(1), &
                              this%coef_b_low, this%coef_b_up, &
                              this%pos_corner(i_vpml,:,:), p_upper)

        endif

      ! x3 direction
      case (3)

        ! Bottom wall
        if (this%loc_vpml(i_vpml) == p_lower) then

          update_range(1) = vpml_range(1)
          update_range(2) = vpml_range(2)-1

          call update_b_x3_3d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                                 this%coef_b_low(3,:,:), this%dtdx, &
                                 update_range, vpml_range(1)-2 )

        ! Top wall
        else

          update_range(1) = vpml_range(1)+1
          update_range(2) = vpml_range(2)-1

          call update_b_x3_3d(this%wall_array_b(i_vpml), this%wall_array_e(i_vpml), &
                                 this%coef_b_up(3,:,:), this%dtdx, &
                                 update_range, vpml_range(1) )

        endif

    end select

  enddo

end subroutine update_b_vpml_3d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_b_x1_1d(vpml_b, vpml_e, coef_b, loc_vpml)
!
! Updates magnetic field in the vay boundary
! Field components in 2D are
! 1 => X (XY = XZ)
! 2 => Y (YX = YZ)
! 3 => ZX
! 4 => ZY
!---------------------------------------------------
  implicit none

  type (t_wall), intent(inout) :: vpml_b
  type(t_wall), intent(in) :: vpml_e
  real(p_k_fld), dimension(:,:), intent(in) :: coef_b
  integer, intent(in)  :: loc_vpml

  ! local variables
   integer, dimension(2) :: vpml_range, upd_range
  integer :: i1, coef_i, coef_idx

  vpml_range = range(vpml_b)

  ! Left wall
  if (loc_vpml == p_lower) then

    upd_range(1) = vpml_range(1)
    upd_range(2) = vpml_range(2)-1

    coef_idx = vpml_range(1)-2

  ! Right wall
  else

    upd_range(1) = vpml_range(1)+1
    upd_range(2) = vpml_range(2)-1

    coef_idx = vpml_range(1)

  endif

  ! Update vpml
  do i1=upd_range(1), upd_range(2)

    ! coefficients are always defined from 1:max_bnd_size
    coef_i = i1-coef_idx

    ! X component
    ! vpml_b%f1(1,i1) = vpml_b%f1(1,i1)

    ! Y component
    vpml_b%f1(2,i1) = coef_b(1,coef_i) *   vpml_b%f1(2,i1) &
                         + coef_b(2,coef_i) * vpml_e%f1(3,i1+1) &
                         - coef_b(3,coef_i) * vpml_e%f1(3,i1  )

    ! Z component
    vpml_b%f1(3,i1) = coef_b(1,coef_i) *   vpml_b%f1(3,i1) &
                         - coef_b(2,coef_i) * vpml_e%f1(2,i1+1) &
                         + coef_b(3,coef_i) * vpml_e%f1(2,i1  )

  enddo

end subroutine update_b_x1_1d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_b_x1_2d(vpml_b, vpml_e, dtdx, range, coef_idx, &
                          coef_b_low, coef_b_up, pos_corner, pos_vpml )
!
! Updates magnetic field in the vay boundary
! Field components in 2D are
! 1 => X (XY = XZ)
! 2 => Y (YX = YZ)
! 3 => ZX
! 4 => ZY
!---------------------------------------------------
  implicit none

  type (t_wall), intent(inout) :: vpml_b
  type(t_wall), intent(in) :: vpml_e
  real(p_k_fld), intent(in) :: dtdx

  integer, dimension(:), intent(in) :: range
  integer, intent(in)  :: coef_idx

  integer, dimension(:), intent(in) :: pos_corner

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_b_low

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_b_up

  integer, intent(in) :: pos_vpml

  ! local variables
  integer :: gc1, gc2, i1, i2, coef_i, coef_std
  real(p_k_fld), dimension(:,:), pointer :: coef_b

  real(p_k_fld), dimension(3) :: coef_b_cnr

  logical :: if_lower_cnr, if_upper_cnr, do_corner

  gc1 = vpml_b%gc_num(p_lower,2)
  gc2 = vpml_b%gc_num(p_upper,2)

  ! update invariant components
  if (pos_vpml == p_lower) then
    coef_b => coef_b_low(1,:,:)
    coef_std = coef_idx - 2
  else
    coef_b => coef_b_up(1,:,:)
    coef_std = coef_idx
  endif

  ! check corner positions
  if_lower_cnr = (pos_corner(p_lower) > 0)
  if_upper_cnr = (pos_corner(p_upper) > 0)

  ! loop vertically through wall
  do i2=1-gc1, vpml_b%nx(2)+gc2-1

    do_corner = .false.

    if ( if_lower_cnr .and. (i2 <= 0) ) then
      coef_i = i2+gc1+1
      coef_b_cnr = coef_b_low(2,:,coef_i)
      do_corner = .true.
    else if ( if_upper_cnr .and. (i2 > vpml_b%nx(2)+1 ) ) then

      coef_i = i2 - vpml_b%nx(2)-1
      coef_b_cnr = coef_b_up(2,:,coef_i)
      do_corner = .true.
    endif

    ! loop horizontally through wall
    do i1=range(1), range(2)

      ! Use perpendicular equations for corner
      if (do_corner) then

        ! X component => perp solver
        vpml_b%f2(1,i1,i2) = coef_b_cnr(1) *   vpml_b%f2(1,i1,i2) - &
                             coef_b_cnr(2) * ( vpml_e%f2(3,i1,i2+1) + vpml_e%f2(4,i1,i2+1) ) + &
                             coef_b_cnr(3) * ( vpml_e%f2(3,i1,i2  ) + vpml_e%f2(4,i1,i2  ) )

        ! ZY component => perp solver
        vpml_b%f2(4,i1,i2) = coef_b_cnr(1) * vpml_b%f2(4,i1,i2) + &
                             coef_b_cnr(2) * vpml_e%f2(1,i1,i2+1) - &
                             coef_b_cnr(3) * vpml_e%f2(1,i1,i2  )

      ! standard solver for wall
      else

        ! X component => standard solver
        vpml_b%f2(1,i1,i2) = vpml_b%f2(1,i1,i2) &
                             - dtdx*( vpml_e%f2(3,i1,i2+1) + vpml_e%f2(4,i1,i2+1) &
                                    - vpml_e%f2(3,i1,i2  ) - vpml_e%f2(4,i1,i2  ) )

        ! ZY component => standard solver
        vpml_b%f2(4,i1,i2) = vpml_b%f2(4,i1,i2) + dtdx*( vpml_e%f2(1,i1,i2+1) - vpml_e%f2(1,i1,i2) )

      endif

      coef_i = i1 - coef_std

      ! Y component
      vpml_b%f2(2,i1,i2) = coef_b(1,coef_i) *   vpml_b%f2(2,i1,i2) + &
                           coef_b(2,coef_i) * ( vpml_e%f2(3,i1+1,i2) + vpml_e%f2(4,i1+1,i2) ) - &
                           coef_b(3,coef_i) * ( vpml_e%f2(3,i1  ,i2) + vpml_e%f2(4,i1  ,i2) )

      ! ZX component
      vpml_b%f2(3,i1,i2) = coef_b(1,coef_i) * vpml_b%f2(3,i1  ,i2) - &
                           coef_b(2,coef_i) * vpml_e%f2(2,i1+1,i2) + &
                           coef_b(3,coef_i) * vpml_e%f2(2,i1  ,i2)

    enddo
  enddo

end subroutine update_b_x1_2d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_cyl_b_x1_2d(vpml_b, vpml_e, dtdr, range, coef_idx, &
                          coef_b_low, coef_b_up, pos_corner, pos_vpml, gir_pos )
!
! Updates magnetic field in the vay boundary
! Field components in 2D are
! 1 => X (XY = XZ) <= Only this component is different in Cyl (and not for the corner)
! 2 => Y (YX = YZ)
! 3 => ZX
! 4 => ZY
!
! Only component 1 is different from the cartesian case
!---------------------------------------------------
  implicit none

  type (t_wall), intent(inout) :: vpml_b
  type(t_wall), intent(in) :: vpml_e
  real(p_k_fld), intent(in) :: dtdr

  integer, dimension(:), intent(in) :: range
  integer, intent(in)  :: coef_idx

  integer, dimension(:), intent(in) :: pos_corner

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_b_low

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_b_up

  integer, intent(in) :: pos_vpml
  integer, intent(in) :: gir_pos

  ! local variables
  integer :: gc1, gc2, i1, i2, coef_i, coef_std
  real(p_k_fld), dimension(:,:), pointer :: coef_b

  real(p_k_fld), dimension(3) :: coef_b_cnr

  logical :: if_upper_cnr, do_corner

  real(p_k_fld) :: tmp, rm, rp
  integer :: gshift_i2, i2_0

  gc1 = vpml_b%gc_num(p_lower,2)
  gc2 = vpml_b%gc_num(p_upper,2)

  ! update invariant components
  if (pos_vpml == p_lower) then
    coef_b => coef_b_low(1,:,:)
    coef_std = coef_idx - 2
  else
    coef_b => coef_b_up(1,:,:)
    coef_std = coef_idx
  endif

  ! check corner positions
  if_upper_cnr = (pos_corner(p_upper) > 0)

  ! Get global cell shift
  gshift_i2 = gir_pos - 2

  if ( gshift_i2 < 0 ) then
    ! This node contains the cylindrical axis so start the solver in cell 2
    i2_0 = 2
  else
    i2_0 = 0
  endif

  ! loop vertically through wall
  do i2 = i2_0, vpml_b%nx(2)+gc2-1

    do_corner = .false.
    if ( if_upper_cnr .and. (i2 > vpml_b%nx(2)+1 ) ) then

      coef_i = i2 - vpml_b%nx(2)-1
      coef_b_cnr = coef_b_up(2,:,coef_i)
      do_corner = .true.
    endif

    rp   = i2 + gshift_i2 + 0.5     ! r at the corner of cell i2 + 1
    rm   = i2 + gshift_i2 - 0.5     ! r at the corner of cell i2
    tmp  = dtdr/(i2 + gshift_i2)    ! (1/r)*(dt/dr) at the center of cell i2

    ! loop horizontally through wall
    do i1=range(1), range(2)

      ! Use perpendicular equations for corner
      if (do_corner) then

        ! X component => perp solver
        vpml_b%f2(1,i1,i2) = coef_b_cnr(1) *   vpml_b%f2(1,i1,i2) &
                             - coef_b_cnr(2) * ( vpml_e%f2(3,i1,i2+1) + vpml_e%f2(4,i1,i2+1) ) &
                             + coef_b_cnr(3) * ( vpml_e%f2(3,i1,i2  ) + vpml_e%f2(4,i1,i2  ) )

        ! ZY component => perp solver
        vpml_b%f2(4,i1,i2) = coef_b_cnr(1) *   vpml_b%f2(4,i1,i2) &
                             + coef_b_cnr(2) * vpml_e%f2(1,i1,i2+1) &
                             - coef_b_cnr(3) * vpml_e%f2(1,i1,i2  )

      ! standard solver for wall
      else

        ! X component => standard solver
        vpml_b%f2(1,i1,i2) = vpml_b%f2(1,i1,i2) - tmp * &

                              ( rp * (vpml_e%f2(3,i1,i2+1) + vpml_e%f2(4,i1,i2+1)) - &
                                rm * (vpml_e%f2(3,i1,i2  ) + vpml_e%f2(4,i1,i2  )))

        ! ZY component => standard solver
        vpml_b%f2(4,i1,i2) = vpml_b%f2(4,i1,i2) + dtdr*( vpml_e%f2(1,i1,i2+1) - vpml_e%f2(1,i1,i2) )

      endif

      coef_i = i1 - coef_std

      ! Y component
      vpml_b%f2(2,i1,i2) = coef_b(1,coef_i) *   vpml_b%f2(2,i1,i2) &
                           + coef_b(2,coef_i) * ( vpml_e%f2(3,i1+1,i2) + vpml_e%f2(4,i1+1,i2) ) &
                           - coef_b(3,coef_i) * ( vpml_e%f2(3,i1  ,i2) + vpml_e%f2(4,i1  ,i2) )

      ! ZX component
      vpml_b%f2(3,i1,i2) = coef_b(1,coef_i) * vpml_b%f2(3,i1  ,i2) &
                           - coef_b(2,coef_i) * vpml_e%f2(2,i1+1,i2) &
                           + coef_b(3,coef_i) * vpml_e%f2(2,i1  ,i2)

    enddo
  enddo

end subroutine update_cyl_b_x1_2d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_b_x1_3d(vpml_b, vpml_e, dtdx, range, coef_idx, &
                          coef_b_low, coef_b_up, pos_corner, pos_vpml )
!
! Updates magnetic field in the left vay boundary
! Field components in 2D are
! 1 => XY
! 2 => XZ
!
! 3 => YX
! 4 => YZ
!
! 5 => ZX
! 6 => ZY
!---------------------------------------------------
  implicit none

  type (t_wall), intent(inout) :: vpml_b
  type(t_wall), intent(in) :: vpml_e
  real(p_k_fld), dimension(:), intent(in) :: dtdx

  integer, dimension(:), intent(in) :: range
  integer, intent(in)  :: coef_idx

  integer, dimension(:,:), intent(in) :: pos_corner

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_b_low

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_b_up

  integer, intent(in) :: pos_vpml

  ! local variables
  integer :: gc2a, gc2b, gc3a, gc3b, i1, i2, i3, coef_i, coef_std
  real(p_k_fld), dimension(:,:), pointer :: coef_b

  real(p_k_fld), dimension(3) :: coef_b_cnr_x2, coef_b_cnr_x3

  logical, dimension(2) :: if_cnr_x2, if_cnr_x3
  logical :: do_corner_x2, do_corner_x3

  integer :: pi3

  pi3 = 3

  gc2a = vpml_b%gc_num(p_lower,2)
  gc2b = vpml_b%gc_num(p_upper,2)

  gc3a = vpml_b%gc_num(p_lower,3)
  gc3b = vpml_b%gc_num(p_upper,3)

  ! update invariant components
  if (pos_vpml == p_lower) then
    coef_b => coef_b_low(1,:,:)
    coef_std = coef_idx - 2
  else
    coef_b => coef_b_up(1,:,:)
    coef_std = coef_idx
  endif

  ! check corner positions
  ! SFM: should these be global variables of the VPML?
  if_cnr_x2(p_lower) = (pos_corner(2,p_lower) > 0)
  if_cnr_x2(p_upper) = (pos_corner(2,p_upper) > 0)

  if_cnr_x3(p_lower) = (pos_corner(pi3,p_lower) > 0)
  if_cnr_x3(p_upper) = (pos_corner(pi3,p_upper) > 0)

  do i3=1-gc3a, vpml_b%nx(3)+gc3b-1

    ! check if in a X3 edge
    do_corner_x3 = .false.
    if ( if_cnr_x3(p_lower) .and. (i3 <= 0) ) then

      coef_i = i3+gc3a+1
      coef_b_cnr_x3 = coef_b_low(3,:,coef_i)
      do_corner_x3 = .true.
    else if ( if_cnr_x3(p_upper) .and. (i3 > vpml_b%nx(3)+1 ) ) then

      coef_i = i3 - vpml_b%nx(3)-1
      coef_b_cnr_x3 = coef_b_up(3,:,coef_i)
      do_corner_x3 = .true.
    endif

    do i2=1-gc2a, vpml_b%nx(2)+gc2b-1

      ! check if in a X2 edge
      do_corner_x2 = .false.
      if ( if_cnr_x2(p_lower) .and. (i2 <= 0) ) then

        coef_i = i2+gc2a+1
        coef_b_cnr_x2 = coef_b_low(2,:,coef_i)
        do_corner_x2 = .true.
      else if ( if_cnr_x2(p_upper) .and. (i2 > vpml_b%nx(2)+1 ) ) then

        coef_i = i2 - vpml_b%nx(2)-1
        coef_b_cnr_x2 = coef_b_up(2,:,coef_i)
        do_corner_x2 = .true.
      endif

      do i1=range(1), range(2)

        ! Use equations for corner X2
        if (do_corner_x2) then

          ! XY component
          vpml_b%f3(1,i1,i2,i3) = coef_b_cnr_x2(1) * vpml_b%f3(1, i1, i2, i3) &
                                  - coef_b_cnr_x2(2) * ( vpml_e%f3(5, i1, i2+1, i3) &
                                                       + vpml_e%f3(6, i1, i2+1, i3) ) &
                                  + coef_b_cnr_x2(3) * ( vpml_e%f3(5, i1, i2  , i3) &
                                                       + vpml_e%f3(6, i1, i2  , i3) )

          ! ZY component => standard solver
          vpml_b%f3(6,i1,i2,i3) = coef_b_cnr_x2(1) * vpml_b%f3(6, i1, i2, i3) &
                                  + coef_b_cnr_x2(2) * ( vpml_e%f3(1, i1, i2+1, i3) &
                                                       + vpml_e%f3(2, i1, i2+1, i3) ) &
                                  - coef_b_cnr_x2(3) * ( vpml_e%f3(1, i1, i2  , i3) &
                                                       + vpml_e%f3(2, i1, i2  , i3) )

        ! standard solver for wall

        else

          ! XY component
          vpml_b%f3(1,i1,i2,i3) = vpml_b%f3(1, i1, i2, i3) &
                                  - dtdx(2)*( vpml_e%f3(5, i1, i2+1, i3) &
                                            + vpml_e%f3(6, i1, i2+1, i3) &
                                            - vpml_e%f3(5, i1, i2  , i3) &
                                            - vpml_e%f3(6, i1, i2  , i3) )

          ! ZY component => standard solver
          vpml_b%f3(6,i1,i2,i3) = vpml_b%f3(6, i1, i2, i3) &
                                  + dtdx(2)*( vpml_e%f3(1, i1, i2+1, i3) &
                                            + vpml_e%f3(2, i1, i2+1, i3) &
                                            - vpml_e%f3(1, i1, i2  , i3) &
                                            - vpml_e%f3(2, i1, i2  , i3) )

        endif

        ! Use equations for corner X3
        if (do_corner_x3) then

          ! XZ component
          vpml_b%f3(2,i1,i2,i3) = coef_b_cnr_x3(1) * vpml_b%f3(2, i1, i2, i3) &
                                  + coef_b_cnr_x3(2) * ( vpml_e%f3(3, i1, i2, i3+1) &
                                                       + vpml_e%f3(4, i1, i2, i3+1) ) &
                                  - coef_b_cnr_x3(3) * ( vpml_e%f3(3, i1, i2, i3  ) &
                                                       + vpml_e%f3(4, i1, i2, i3  ) )

          ! YZ component
          vpml_b%f3(4,i1,i2,i3) = coef_b_cnr_x3(1) * vpml_b%f3(4, i1, i2, i3) &
                                  - coef_b_cnr_x3(2) * ( vpml_e%f3(1, i1, i2, i3+1) &
                                                       + vpml_e%f3(2, i1, i2, i3+1) ) &
                                  + coef_b_cnr_x3(3) * ( vpml_e%f3(1, i1, i2, i3  ) &
                                                       + vpml_e%f3(2, i1, i2, i3  ) )

        ! standard solver for wall

        else

          ! XZ component
          vpml_b%f3(2,i1,i2,i3) = vpml_b%f3(2, i1, i2, i3) &
                                  + dtdx(pi3)*( vpml_e%f3(3, i1, i2, i3+1) &
                                            + vpml_e%f3(4, i1, i2, i3+1) &
                                            - vpml_e%f3(3, i1, i2, i3  ) &
                                            - vpml_e%f3(4, i1, i2, i3  ) )

          ! YZ component
          vpml_b%f3(4,i1,i2,i3) = vpml_b%f3(4, i1, i2, i3) &
                                  - dtdx(pi3)*( vpml_e%f3(1, i1, i2, i3+1) &
                                            + vpml_e%f3(2, i1, i2, i3+1) &
                                            - vpml_e%f3(1, i1, i2, i3  ) &
                                            - vpml_e%f3(2, i1, i2, i3  ) )

        endif

        ! Update invariant field components
        coef_i = i1-coef_std

        ! YX component

        vpml_b%f3(3,i1,i2,i3) = coef_b(1,coef_i) * vpml_b%f3(3, i1, i2, i3) &
                                + coef_b(2,coef_i) * ( vpml_e%f3(5, i1+1, i2, i3 ) &
                                                     + vpml_e%f3(6, i1+1, i2, i3 ) ) &
                                - coef_b(3,coef_i) * ( vpml_e%f3(5, i1  , i2, i3 ) &
                                                     + vpml_e%f3(6, i1  , i2, i3 ) )

        ! ZX component

        vpml_b%f3(5,i1,i2,i3) = coef_b(1,coef_i) * vpml_b%f3(5, i1, i2, i3) &
                                - coef_b(2,coef_i) * ( vpml_e%f3(3, i1+1, i2, i3 ) &
                                                     + vpml_e%f3(4, i1+1, i2, i3 ) ) &
                                + coef_b(3,coef_i) * ( vpml_e%f3(3, i1  , i2, i3 ) &
                                                     + vpml_e%f3(4, i1  , i2, i3 ) )
      enddo
    enddo
  enddo

end subroutine update_b_x1_3d
!---------------------------------------------------

!---------------------------------------------------
subroutine update_b_x2_2d(vpml_b, vpml_e, coef_b, dtdx, range, coef_idx)
!
! Updates magnetic field in the vay boundary
! Field components in 2D are
! 1 => X (XY = XZ)
! 2 => Y (YX = YZ)
! 3 => ZX
! 4 => ZY
!---------------------------------------------------
  implicit none

  type (t_wall), intent(inout) :: vpml_b
  type(t_wall), intent(in) :: vpml_e
  real(p_k_fld), dimension(:,:), intent(in) :: coef_b
  real(p_k_fld), intent(in) :: dtdx

  integer, dimension(:), intent(in) :: range
  integer, intent(in)  :: coef_idx

  ! local variables
  integer :: i1, i2, coef_i

  do i2=range(1), range(2)

    ! coefficients are always defined from 1:max_bnd_size
    coef_i = i2-coef_idx

    do i1 = lbound( vpml_b%f2, 2 ), ubound( vpml_b%f2, 2 ) - 1

      ! X component
      vpml_b%f2(1,i1,i2) = coef_b(1,coef_i) *   vpml_b%f2(1,i1,i2) - &
                           coef_b(2,coef_i) * ( vpml_e%f2(3,i1,i2+1) + vpml_e%f2(4,i1,i2+1) ) + &
                           coef_b(3,coef_i) * ( vpml_e%f2(3,i1,i2  ) + vpml_e%f2(4,i1,i2  ) )

      ! Y component => standard solver
      vpml_b%f2(2,i1,i2) = vpml_b%f2(2,i1,i2) &
                           + dtdx*( vpml_e%f2(3,i1+1,i2) + vpml_e%f2(4,i1+1,i2) &
                                  - vpml_e%f2(3,i1  ,i2) - vpml_e%f2(4,i1  ,i2) )

      ! ZX component => standard solver
      vpml_b%f2(3,i1,i2) = vpml_b%f2(3,i1,i2) - dtdx*( vpml_e%f2(2,i1+1,i2) - vpml_e%f2(2,i1,i2) )

      ! ZY component
      vpml_b%f2(4,i1,i2) = coef_b(1,coef_i) * vpml_b%f2(4,i1,i2  ) + &
                           coef_b(2,coef_i) * vpml_e%f2(1,i1,i2+1) - &
                           coef_b(3,coef_i) * vpml_e%f2(1,i1,i2  )

    enddo
  enddo

  ! On the upper boundary also do the highest row for components Y and ZX
  if ( range(2) < ubound( vpml_b%f2, 3 ) ) then
    i2 = range(2) + 1
    do i1 = lbound( vpml_b%f2, 2 ), ubound( vpml_b%f2, 2 ) - 1

      ! Y component => standard solver
      vpml_b%f2(2,i1,i2) = vpml_b%f2(2,i1,i2) &
                           + dtdx*( vpml_e%f2(3,i1+1,i2) + vpml_e%f2(4,i1+1,i2) &
                                  - vpml_e%f2(3,i1  ,i2) - vpml_e%f2(4,i1  ,i2) )

      ! ZX component => standard solver
      vpml_b%f2(3,i1,i2) = vpml_b%f2(3,i1,i2) - dtdx*( vpml_e%f2(2,i1+1,i2) - vpml_e%f2(2,i1,i2) )

    enddo

  endif

end subroutine update_b_x2_2d
!---------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine update_b_x2_3d(vpml_b, vpml_e, dtdx, range, coef_idx, &
                          coef_b_low, coef_b_up, pos_corner, pos_vpml )
!
! Updates magnetic field in the left vay boundary
! Field components in 2D are
! 1 => XY
! 2 => XZ
!
! 3 => YX
! 4 => YZ
!
! 5 => ZX
! 6 => ZY
!-----------------------------------------------------------------------------------------
  implicit none

  type (t_wall), intent(inout) :: vpml_b
  type(t_wall), intent(in) :: vpml_e
  real(p_k_fld), dimension(:), intent(in) :: dtdx

  integer, dimension(:), intent(in) :: range
  integer, intent(in)  :: coef_idx

  integer, dimension(:,:), intent(in) :: pos_corner

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_b_low

  real(p_k_fld), dimension(:,:,:), intent(in), target :: coef_b_up

  integer, intent(in) :: pos_vpml

  ! local variables
  integer :: gc1a, gc1b, gc3a, gc3b, i1, i2, i3, coef_i, coef_std
  real(p_k_fld), dimension(:,:), pointer :: coef_b

  real(p_k_fld), dimension(3) :: coef_b_cnr_x3

  logical, dimension(2) :: if_cnr_x3
  logical :: do_corner_x3

  integer :: pi3

  pi3 = 3

  gc1a = vpml_e%gc_num(p_lower,1)
  gc1b = vpml_e%gc_num(p_upper,1)

  gc3a = vpml_e%gc_num(p_lower,3)
  gc3b = vpml_e%gc_num(p_upper,3)

  ! update invariant components
  if (pos_vpml == p_lower) then
    coef_b => coef_b_low(2,:,:)
    coef_std = coef_idx - 2
  else
    coef_b => coef_b_up(2,:,:)
    coef_std = coef_idx
  endif

  ! check corner positions
  ! SFM: should these be global variables of the VPML?
  if_cnr_x3(p_lower) = (pos_corner(pi3,p_lower) > 0)
  if_cnr_x3(p_upper) = (pos_corner(pi3,p_upper) > 0)

  do i3=1-gc3a, vpml_e%nx(3)+gc3b-1

    ! check if in a X3 edge
    do_corner_x3 = .false.
    if ( if_cnr_x3(p_lower) .and. (i3 <= 0) ) then

      coef_i = i3+gc3a+1
      coef_b_cnr_x3 = coef_b_low(3,:,coef_i)
      do_corner_x3 = .true.
    else if ( if_cnr_x3(p_upper) .and. (i3 > vpml_b%nx(3)+1 ) ) then

      coef_i = i3 - vpml_b%nx(3)-1
      coef_b_cnr_x3 = coef_b_up(3,:,coef_i)
      do_corner_x3 = .true.
    endif

    do i2=range(1), range(2)

      ! coefficient index for invariant field components
      coef_i = i2-coef_std

      do i1=1-gc1a, vpml_e%nx(1)+gc1b-1

        ! Use equations for corner X3
        if (do_corner_x3) then

          ! XZ component => standard solver
          vpml_b%f3(2,i1,i2,i3) = coef_b_cnr_x3(1) * vpml_b%f3(2, i1, i2, i3) &
                                  + coef_b_cnr_x3(2) * ( vpml_e%f3(3, i1, i2, i3+1) &
                                                       + vpml_e%f3(4, i1, i2, i3+1) ) &
                                  - coef_b_cnr_x3(3) * ( vpml_e%f3(3, i1, i2, i3  ) &
                                                       + vpml_e%f3(4, i1, i2, i3  ) )

          ! YZ component => standard solver
          vpml_b%f3(4,i1,i2,i3) = coef_b_cnr_x3(1) * vpml_b%f3(4, i1, i2, i3) &
                                  - coef_b_cnr_x3(2) * ( vpml_e%f3(1, i1, i2, i3+1) &
                                                       + vpml_e%f3(2, i1, i2, i3+1) ) &
                                  + coef_b_cnr_x3(3) * ( vpml_e%f3(1, i1, i2, i3  ) &
                                                       + vpml_e%f3(2, i1, i2, i3  ) )

        ! standard solver for wall

        else

          ! XZ component
          vpml_b%f3(2,i1,i2,i3) = vpml_b%f3(2, i1, i2, i3) &
                                  + dtdx(pi3)*( vpml_e%f3(3, i1, i2, i3+1) &
                                            + vpml_e%f3(4, i1, i2, i3+1) &
                                            - vpml_e%f3(3, i1, i2, i3  ) &
                                            - vpml_e%f3(4, i1, i2, i3  ) )

          ! YZ component
          vpml_b%f3(4,i1,i2,i3) = vpml_b%f3(4, i1, i2, i3) &
                                  - dtdx(pi3)*( vpml_e%f3(1, i1, i2, i3+1) &
                                            + vpml_e%f3(2, i1, i2, i3+1) &
                                            - vpml_e%f3(1, i1, i2, i3  ) &
                                            - vpml_e%f3(2, i1, i2, i3  ) )

        endif

        ! Update invariant field components

        ! XY component
        vpml_b%f3(1,i1,i2,i3) = coef_b(1,coef_i) * vpml_b%f3(1, i1, i2, i3) &
                                - coef_b(2,coef_i) * ( vpml_e%f3(5, i1, i2+1, i3 ) &
                                                     + vpml_e%f3(6, i1, i2+1, i3 ) ) &
                                + coef_b(3,coef_i) * ( vpml_e%f3(5, i1, i2  , i3 ) &
                                                     + vpml_e%f3(6, i1, i2  , i3 ) )

        ! YX component => standard solver
        vpml_b%f3(3,i1,i2,i3) = vpml_b%f3(3, i1, i2, i3) &
                                + dtdx(1)*( vpml_e%f3(5, i1+1, i2, i3) &
                                          + vpml_e%f3(6, i1+1, i2, i3) &
                                          - vpml_e%f3(5, i1  , i2, i3) &
                                          - vpml_e%f3(6, i1  , i2, i3) )

        ! ZX component => standard solver
        vpml_b%f3(5,i1,i2,i3) = vpml_b%f3(5, i1, i2, i3) &
                                - dtdx(1)*( vpml_e%f3(3, i1+1, i2, i3) &
                                          + vpml_e%f3(4, i1+1, i2, i3) &
                                          - vpml_e%f3(3, i1  , i2, i3) &
                                          - vpml_e%f3(4, i1  , i2, i3) )

        ! ZY component
        vpml_b%f3(6,i1,i2,i3) = coef_b(1,coef_i) * vpml_b%f3(6, i1, i2, i3) &
                                + coef_b(2,coef_i) * ( vpml_e%f3(1, i1, i2+1, i3 ) &
                                                     + vpml_e%f3(2, i1, i2+1, i3 ) ) &
                                - coef_b(3,coef_i) * ( vpml_e%f3(1, i1, i2  , i3 ) &
                                                     + vpml_e%f3(2, i1, i2  , i3 ) )

      enddo
    enddo
  enddo

end subroutine update_b_x2_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine update_b_x3_3d(vpml_b, vpml_e, coef_b, dtdx, range, coef_idx)
! Updates magnetic field in the left vay boundary
! Field components in 3D are
! 1 => XY
! 2 => XZ
!
! 3 => YX
! 4 => YZ
!
! 5 => ZX
! 6 => ZY
!-----------------------------------------------------------------------------------------
  implicit none

  type (t_wall), intent(inout) :: vpml_b
  type(t_wall), intent(in) :: vpml_e
  real(p_k_fld), dimension(:,:), intent(in) :: coef_b
  real(p_k_fld), dimension(:), intent(in) :: dtdx

  integer, dimension(:), intent(in) :: range
  integer, intent(in)  :: coef_idx

  ! local variables
  integer :: gc1a, gc1b, gc2a, gc2b, i1, i2, i3, coef_i

  gc1a = vpml_e%gc_num(p_lower,1)
  gc1b = vpml_e%gc_num(p_upper,1)

  gc2a = vpml_e%gc_num(p_lower,2)
  gc2b = vpml_e%gc_num(p_upper,2)

  do i3=range(1), range(2)
    coef_i = i3-coef_idx

    do i2=1-gc2a, vpml_e%nx(2)+gc2b-1
      do i1=1-gc1a, vpml_e%nx(1)+gc1b-1

        ! XY component => standard solver
        vpml_b%f3(1,i1,i2,i3) = vpml_b%f3(1, i1, i2, i3) &
                                - dtdx(2)*( vpml_e%f3(5, i1, i2+1, i3  ) &
                                          + vpml_e%f3(6, i1, i2+1, i3  ) &
                                          - vpml_e%f3(5, i1, i2  , i3) &
                                          - vpml_e%f3(6, i1, i2  , i3) )

        ! XZ component
        vpml_b%f3(2,i1,i2,i3) = coef_b(1,coef_i) * vpml_b%f3(2, i1, i2, i3) &
                                + coef_b(2,coef_i) * ( vpml_e%f3(3, i1, i2, i3+1 ) &
                                                     + vpml_e%f3(4, i1, i2, i3+1 ) ) &
                                - coef_b(3,coef_i) * ( vpml_e%f3(3, i1, i2, i3   ) &
                                                     + vpml_e%f3(4, i1, i2, i3   ) )

        ! YX component => standard solver
        vpml_b%f3(3,i1,i2,i3) = vpml_b%f3(3, i1, i2, i3) &
                                + dtdx(1)*( vpml_e%f3(5, i1+1, i2, i3) &
                                          + vpml_e%f3(6, i1+1, i2, i3) &
                                          - vpml_e%f3(5, i1  , i2, i3) &
                                          - vpml_e%f3(6, i1  , i2, i3) )

        ! YZ component
        vpml_b%f3(4,i1,i2,i3) = coef_b(1,coef_i) * vpml_b%f3(4, i1, i2, i3) &
                                - coef_b(2,coef_i) * ( vpml_e%f3(1, i1, i2, i3+1 ) &
                                                     + vpml_e%f3(2, i1, i2, i3+1 ) ) &
                                + coef_b(3,coef_i) * ( vpml_e%f3(1, i1, i2, i3   ) &
                                                     + vpml_e%f3(2, i1, i2, i3   ) )

        ! ZX component => standard solver
        vpml_b%f3(5,i1,i2,i3) = vpml_b%f3(5, i1, i2, i3) &
                                - dtdx(1)*( vpml_e%f3(3, i1+1, i2, i3) &
                                          + vpml_e%f3(4, i1+1, i2, i3) &
                                          - vpml_e%f3(3, i1  , i2, i3) &
                                          - vpml_e%f3(4, i1  , i2, i3) )

        ! ZY component => standard solver
        vpml_b%f3(6,i1,i2,i3) = vpml_b%f3(6, i1, i2, i3) &
                                + dtdx(2)*( vpml_e%f3(1, i1, i2+1, i3  ) &
                                          + vpml_e%f3(2, i1, i2+1, i3  ) &
                                          - vpml_e%f3(1, i1, i2  , i3) &
                                          - vpml_e%f3(2, i1, i2  , i3) )

      enddo
    enddo
  enddo

end subroutine update_b_x3_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine update_interface_x1_1d(vpml_wall, field, loc_vpml)
!
! Updates field at the interface
! Field components in 1D are
! 1 => X
! 2 => Y
! 3 => Z
!-----------------------------------------------------------------------------------------
  implicit none

  type (t_wall), intent(inout) :: vpml_wall
  type(t_vdf),intent(inout) :: field
  integer, intent(in) :: loc_vpml

  ! local variables
  integer, dimension(2) :: vpml_range

  integer :: box_to_vpml, vpml_to_box

  ! Get vpml update ranges
  vpml_range = range(vpml_wall)

  ! Left wall
  if (loc_vpml == p_lower) then
    box_to_vpml = vpml_range(p_upper)
    vpml_to_box = vpml_range(p_upper)-1

  ! Right wall
  else
    box_to_vpml = vpml_range(p_lower)
    vpml_to_box = vpml_range(p_lower)+1
  endif

  ! Copy box values to VPML boundary from simulation box
  vpml_wall%f1(1, box_to_vpml) = field%f1(1, box_to_vpml)      ! X component

  vpml_wall%f1(2, box_to_vpml) = field%f1(2, box_to_vpml)      ! Y component
  vpml_wall%f1(3 ,box_to_vpml) = field%f1(3, box_to_vpml)      ! Z component

  ! Copy VPML boundary values to box
  field%f1(1, vpml_to_box) = vpml_wall%f1(1, vpml_to_box)
  field%f1(2, vpml_to_box) = vpml_wall%f1(2, vpml_to_box)
  field%f1(3, vpml_to_box) = vpml_wall%f1(3, vpml_to_box)

end subroutine update_interface_x1_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine update_interface_x1_2d(vpml_wall, i_vpml, pos_corner, field, loc_vpml )
!
! Updates field at the interface
! Field components in 2D are
! 1 => X (XY = XZ)
! 2 => Y (YX = YZ)
! 3 => ZX
! 4 => ZY
!-----------------------------------------------------------------------------------------
  implicit none

  type (t_wall), dimension(:), intent(inout) :: vpml_wall
  integer, intent(in) :: i_vpml                               ! vpml to update
  integer, dimension(:), intent(in) :: pos_corner    ! contact vpmls
  type(t_vdf),intent(inout) :: field
  integer, intent(in) :: loc_vpml

  ! local variables
  integer :: box_to_vpml, vpml_to_box ! cells to update
  integer :: i1, i2, adj_vpml, bnd_idx
  integer, dimension(2) :: vpml_range, bnd_vpml, i2_range

  ! Get vpml update ranges
  vpml_range = range(vpml_wall(i_vpml))

  ! Left wall
  if ( loc_vpml == p_lower ) then
      box_to_vpml = vpml_range(p_upper)
      vpml_to_box = vpml_range(p_upper)-1
  else
      box_to_vpml = vpml_range(p_lower)
      vpml_to_box = vpml_range(p_lower)+1
  endif

  ! Copy box values to VPML boundary
  do i2=0, field%nx_(2)+2
    vpml_wall(i_vpml)%f2(1, box_to_vpml, i2) = field%f2(1, box_to_vpml, i2)      ! X component

    vpml_wall(i_vpml)%f2(2, box_to_vpml, i2) = field%f2(2, box_to_vpml, i2)      ! Y component
    vpml_wall(i_vpml)%f2(3 ,box_to_vpml, i2) = field%f2(3, box_to_vpml, i2)      ! ZX component
    vpml_wall(i_vpml)%f2(4, box_to_vpml, i2) = 0.0_p_k_fld                       ! ZY component (damped)
  enddo

  ! Copy VPML boundary values to box
  do i2=0, field%nx_(2)+2
    field%f2(1, vpml_to_box, i2) = vpml_wall(i_vpml)%f2(1, vpml_to_box, i2)
    field%f2(2, vpml_to_box, i2) = vpml_wall(i_vpml)%f2(2, vpml_to_box, i2)
    field%f2(3, vpml_to_box, i2) = vpml_wall(i_vpml)%f2(3, vpml_to_box, i2) + &
                                   vpml_wall(i_vpml)%f2(4, vpml_to_box, i2)
  enddo

  ! update contact interfaces if necessary
  do bnd_idx=1,2

    ! adjacent vpml
    adj_vpml = pos_corner(bnd_idx)

    if (adj_vpml > 0) then

      bnd_vpml = range(vpml_wall(adj_vpml))

      ! bottom corner
      if (bnd_vpml(p_lower) < 0) then
         i2_range(p_lower) = bnd_vpml(p_lower)
         i2_range(p_upper) = bnd_vpml(p_upper)-1
      ! top corner
      else
         i2_range(p_lower) = bnd_vpml(p_lower)+1
         i2_range(p_upper) = bnd_vpml(p_upper)
      endif

      select case ( loc_vpml )
        case ( p_lower )
          do i2=i2_range(p_lower), i2_range(p_upper)
            vpml_wall(adj_vpml)%f2(1, 0, i2) = vpml_wall(i_vpml)%f2(1, 0, i2)

            vpml_wall(adj_vpml)%f2(2, 0, i2) = vpml_wall(i_vpml)%f2(2, 0, i2)

            vpml_wall(adj_vpml)%f2(3, 0, i2) = vpml_wall(i_vpml)%f2(3, 0, i2)
            vpml_wall(adj_vpml)%f2(4, 0, i2) = vpml_wall(i_vpml)%f2(4, 0, i2)

            vpml_wall(i_vpml)%f2(1, 1, i2) = vpml_wall(adj_vpml)%f2(1, 1, i2)

            vpml_wall(i_vpml)%f2(2, 1, i2) = vpml_wall(adj_vpml)%f2(2, 1, i2)

            vpml_wall(i_vpml)%f2(3, 1, i2) = vpml_wall(adj_vpml)%f2(3, 1, i2)
            vpml_wall(i_vpml)%f2(4, 1, i2) = vpml_wall(adj_vpml)%f2(4, 1, i2)

          enddo

        case ( p_upper )
          i1 = field%nx_(1)+1
          do i2=i2_range(p_lower), i2_range(p_upper)
            vpml_wall(i_vpml)%f2(1, i1, i2) = vpml_wall(adj_vpml)%f2(1, i1, i2)

            vpml_wall(i_vpml)%f2(2, i1, i2) = vpml_wall(adj_vpml)%f2(2, i1, i2)

            vpml_wall(i_vpml)%f2(3, i1, i2) = vpml_wall(adj_vpml)%f2(3, i1, i2)
            vpml_wall(i_vpml)%f2(4, i1, i2) = vpml_wall(adj_vpml)%f2(4, i1, i2)

            vpml_wall(adj_vpml)%f2(1, i1+1, i2) = vpml_wall(i_vpml)%f2(1, i1+1, i2)

            vpml_wall(adj_vpml)%f2(2, i1+1, i2) = vpml_wall(i_vpml)%f2(2, i1+1, i2)

            vpml_wall(adj_vpml)%f2(3, i1+1, i2) = vpml_wall(i_vpml)%f2(3, i1+1, i2)
            vpml_wall(adj_vpml)%f2(4, i1+1, i2) = vpml_wall(i_vpml)%f2(4, i1+1, i2)

          enddo
      end select

    endif
  enddo

end subroutine update_interface_x1_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine update_interface_x1_3d(vpml_wall, i_vpml, pos_corner, field, loc_vpml)
!
! Updates field at the interface
! Field components in 2D are
! 1 => XY
! 2 => XZ
!
! 3 => YX
! 4 => YZ
!
! 5 => ZX
! 6 => ZY
!-----------------------------------------------------------------------------------------
  implicit none

  type (t_wall), dimension(:), intent(inout) :: vpml_wall
  integer, intent(in) :: i_vpml                               ! vpml to update
  integer, dimension(:, :, :), intent(in) :: pos_corner    ! contact vpmls
  type(t_vdf),intent(inout) :: field
  integer, intent(in) :: loc_vpml

  ! local variables
  integer :: box_to_vpml, vpml_to_box
  integer :: i2_min, i2_max, i3_min, i3_max
  integer :: i1, i2, i3, adj_vpml, bnd_idx, dir_idx
  integer, dimension(2) :: vpml_range, bnd_vpml
  integer, dimension(3,2) :: i_range

  integer :: pi3
  pi3 = 3

  ! Get vpml update ranges
  call get_range(vpml_wall(i_vpml), vpml_range)

  ! Left wall
  if (loc_vpml == p_lower) then
    box_to_vpml = vpml_range(p_upper)
    vpml_to_box = vpml_range(p_upper)-1

  ! Right wall
  else
   box_to_vpml = vpml_range(p_lower)
   vpml_to_box = vpml_range(p_lower)+1

  endif

  ! Copy box values to VPML boundary
  i2_min = 1
  i2_max = field%nx_(2)+1
  i3_min = 1
  i3_max = field%nx_(3)+1

  do i3= i3_min, i3_max
    do i2 = i2_min, i2_max

      vpml_wall(i_vpml)%f3(1, box_to_vpml, i2, i3) = field%f3(1, box_to_vpml, i2, i3)    ! XY component

      vpml_wall(i_vpml)%f3(2, box_to_vpml, i2, i3) = 0.0_p_k_fld                         ! XZ component

      vpml_wall(i_vpml)%f3(3, box_to_vpml, i2, i3) = field%f3(2, box_to_vpml, i2, i3)    ! YX component
      vpml_wall(i_vpml)%f3(4, box_to_vpml, i2, i3) = 0.0_p_k_fld                         ! YZ component

      vpml_wall(i_vpml)%f3(5, box_to_vpml, i2, i3) = field%f3(3, box_to_vpml, i2, i3)    ! ZX component
      vpml_wall(i_vpml)%f3(6, box_to_vpml, i2, i3) = 0.0_p_k_fld                         ! ZY component

    enddo
  enddo

  ! Copy VPML boundary values to box
  ! (the boundaries for the loop need to be set as variables to go around a compiler bug in the
  ! intel ifort 11.1 series)
  i2_min = 0
  i2_max = field%nx_(2)+2
  i3_min = 0
  i3_max = field%nx_(3)+2

  ! Copy VPML boundary values to box
  do i3= i3_min, i3_max
    do i2 = i2_min, i2_max

      field%f3(1, vpml_to_box, i2, i3) = vpml_wall(i_vpml)%f3(1, vpml_to_box, i2, i3) + &
                                         vpml_wall(i_vpml)%f3(2, vpml_to_box, i2, i3)

      field%f3(2, vpml_to_box, i2, i3) = vpml_wall(i_vpml)%f3(3, vpml_to_box, i2, i3) + &
                                         vpml_wall(i_vpml)%f3(4, vpml_to_box, i2, i3)

      field%f3(3, vpml_to_box, i2, i3) = vpml_wall(i_vpml)%f3(5, vpml_to_box, i2, i3) + &
                                         vpml_wall(i_vpml)%f3(6, vpml_to_box, i2, i3)

    enddo
  enddo

  i2_min = 1
  i2_max = vpml_wall(i_vpml)%nx(2) + 1

  i3_min = 1 - vpml_wall(i_vpml)%gc_num(p_lower,3)
  i3_max = vpml_wall(i_vpml)%nx(3) + vpml_wall(i_vpml)%gc_num(p_upper,3)

  ! update contact interfaces if necessary
  do dir_idx=2,3

    ! set default ranges
    ! one direction will be changed in the next loop
    i_range(2,p_lower) = i2_min
    i_range(2,p_upper) = i2_max
    i_range(3,p_lower) = i3_min
    i_range(3,p_upper) = i3_max

    do bnd_idx=1,2

      ! adjacent vpml
      adj_vpml = pos_corner(i_vpml, dir_idx, bnd_idx)

      if (adj_vpml > 0) then

        call get_range(vpml_wall(adj_vpml), bnd_vpml)

        ! bottom corner
        if (bnd_vpml(p_lower) < 0) then
           i_range(dir_idx,p_lower) = bnd_vpml(p_lower)
           i_range(dir_idx,p_upper) = bnd_vpml(p_upper)-1

        ! top corner
        else
           i_range(dir_idx,p_lower) = bnd_vpml(p_lower)+1
           i_range(dir_idx,p_upper) = bnd_vpml(p_upper)
        endif

        select case ( loc_vpml )
          case ( p_lower )
            do i3=i_range(3,p_lower), i_range(3,p_upper)
              do i2=i_range(2, p_lower), i_range(2, p_upper)
                vpml_wall(adj_vpml)%f3(1, 0, i2, i3 ) = vpml_wall(i_vpml)%f3(1, 0, i2, i3 )

                vpml_wall(adj_vpml)%f3(2, 0, i2, i3 ) = vpml_wall(i_vpml)%f3(2, 0, i2, i3 )

                vpml_wall(adj_vpml)%f3(3, 0, i2, i3 ) = vpml_wall(i_vpml)%f3(3, 0, i2, i3 )
                vpml_wall(adj_vpml)%f3(4, 0, i2, i3 ) = vpml_wall(i_vpml)%f3(4, 0, i2, i3 )

                vpml_wall(adj_vpml)%f3(5, 0, i2, i3 ) = vpml_wall(i_vpml)%f3(5, 0, i2, i3 )
                vpml_wall(adj_vpml)%f3(6, 0, i2, i3 ) = vpml_wall(i_vpml)%f3(6, 0, i2, i3 )

                vpml_wall(i_vpml)%f3(1, 1, i2, i3 ) = vpml_wall(adj_vpml)%f3(1, 1, i2, i3 )

                vpml_wall(i_vpml)%f3(2, 1, i2, i3 ) = vpml_wall(adj_vpml)%f3(2, 1, i2, i3 )

                vpml_wall(i_vpml)%f3(3, 1, i2, i3 ) = vpml_wall(adj_vpml)%f3(3, 1, i2, i3 )
                vpml_wall(i_vpml)%f3(4, 1, i2, i3 ) = vpml_wall(adj_vpml)%f3(4, 1, i2, i3 )

                vpml_wall(i_vpml)%f3(5, 1, i2, i3 ) = vpml_wall(adj_vpml)%f3(5, 1, i2, i3 )
                vpml_wall(i_vpml)%f3(6, 1, i2, i3 ) = vpml_wall(adj_vpml)%f3(6, 1, i2, i3 )

              enddo
            enddo

          case ( p_upper )
            i1 = field%nx_(1)+1
            do i3=i_range(3,p_lower), i_range(3,p_upper)
              do i2=i_range(2, p_lower), i_range(2, p_upper)
                vpml_wall(i_vpml)%f3(1, i1, i2, i3 ) = vpml_wall(adj_vpml)%f3(1, i1, i2, i3 )

                vpml_wall(i_vpml)%f3(2, i1, i2, i3 ) = vpml_wall(adj_vpml)%f3(2, i1, i2, i3 )

                vpml_wall(i_vpml)%f3(3, i1, i2, i3 ) = vpml_wall(adj_vpml)%f3(3, i1, i2, i3 )
                vpml_wall(i_vpml)%f3(4, i1, i2, i3 ) = vpml_wall(adj_vpml)%f3(4, i1, i2, i3 )

                vpml_wall(i_vpml)%f3(5, i1, i2, i3 ) = vpml_wall(adj_vpml)%f3(5, i1, i2, i3 )
                vpml_wall(i_vpml)%f3(6, i1, i2, i3 ) = vpml_wall(adj_vpml)%f3(6, i1, i2, i3 )

                vpml_wall(adj_vpml)%f3(1, i1+1, i2, i3 ) = vpml_wall(i_vpml)%f3(1, i1+1, i2, i3 )

                vpml_wall(adj_vpml)%f3(2, i1+1, i2, i3 ) = vpml_wall(i_vpml)%f3(2, i1+1, i2, i3 )

                vpml_wall(adj_vpml)%f3(3, i1+1, i2, i3 ) = vpml_wall(i_vpml)%f3(3, i1+1, i2, i3 )
                vpml_wall(adj_vpml)%f3(4, i1+1, i2, i3 ) = vpml_wall(i_vpml)%f3(4, i1+1, i2, i3 )

                vpml_wall(adj_vpml)%f3(5, i1+1, i2, i3 ) = vpml_wall(i_vpml)%f3(5, i1+1, i2, i3 )
                vpml_wall(adj_vpml)%f3(6, i1+1, i2, i3 ) = vpml_wall(i_vpml)%f3(6, i1+1, i2, i3 )

              enddo
            enddo
        end select

      endif

    enddo

  enddo

end subroutine update_interface_x1_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine update_interface_x2_2d(vpml_wall, i_vpml, field, loc_vpml)
!
! Updates field at the interface
! Field components in 2D are
! 1 => X (XY = XZ)
! 2 => Y (YX = YZ)
! 3 => ZX
! 4 => ZY
!-----------------------------------------------------------------------------------------
  implicit none

  type (t_wall), dimension(:), intent(inout) :: vpml_wall
  integer, intent(in) :: i_vpml                               ! vpml to update
  type(t_vdf),intent(inout) :: field
  integer, intent(in) :: loc_vpml

  ! local variables
  integer :: box_to_vpml, vpml_to_box ! cells to update

  integer :: i1
  integer, dimension(2) :: vpml_range, bnd_field

  ! Get vpml update ranges
  vpml_range = range(vpml_wall(i_vpml))

  ! Left wall
  if (loc_vpml== p_lower) then
    box_to_vpml = vpml_range(2)
    vpml_to_box = vpml_range(2)-1

  ! Right wall
  else
    box_to_vpml = vpml_range(1)
    vpml_to_box = vpml_range(1)+1
  endif

  bnd_field(p_lower) = 0
  bnd_field(p_upper) = field%nx_(1) + 2

  ! Copy box values to VPML boundary from simulation box
  do i1=0, field%nx_(1)+2
    vpml_wall(i_vpml)%f2(1, i1, box_to_vpml) = field%f2(1, i1, box_to_vpml)    ! X component

    vpml_wall(i_vpml)%f2(2, i1, box_to_vpml) = field%f2(2, i1, box_to_vpml)    ! Y component
    vpml_wall(i_vpml)%f2(3, i1, box_to_vpml) = field%f2(3, i1, box_to_vpml)    ! ZX component
    vpml_wall(i_vpml)%f2(4, i1, box_to_vpml) = 0.0_p_k_fld                     ! ZY component (damped)
  enddo

  ! Copy VPML boundary values to box
  do i1=0, field%nx_(1)+2
    field%f2(1, i1, vpml_to_box) = vpml_wall(i_vpml)%f2(1, i1, vpml_to_box)
    field%f2(2, i1, vpml_to_box) = vpml_wall(i_vpml)%f2(2, i1, vpml_to_box)
    field%f2(3, i1, vpml_to_box) = vpml_wall(i_vpml)%f2(3, i1, vpml_to_box) + &
                                   vpml_wall(i_vpml)%f2(4, i1, vpml_to_box)
  enddo

  ! Corners are all handled in update_interface_x1_2d

end subroutine update_interface_x2_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine update_interface_x2_3d(vpml_wall, i_vpml, pos_corner, field, loc_vpml)
!
! Updates field at the interface
! Field components in 2D are
! 1 => XY
! 2 => XZ
!
! 3 => YX
! 4 => YZ
!
! 5 => ZX
! 6 => ZY
!-----------------------------------------------------------------------------------------
  implicit none

  type (t_wall), dimension(:), intent(inout) :: vpml_wall
  integer, intent(in) :: i_vpml                               ! vpml to update
  integer, dimension(:,:,:), intent(in) :: pos_corner    ! contact vpmls
  type(t_vdf),intent(inout) :: field
  integer, intent(in) :: loc_vpml

  ! local variables
  integer :: box_to_vpml, vpml_to_box
  integer :: i1_min, i1_max, i3_min, i3_max, i1, i2, i3
  integer :: adj_vpml, bnd_idx, dir_idx
  integer, dimension(2) :: vpml_range, bnd_vpml, i_range
  integer, dimension(p_x_dim,2) :: bnd_field

  ! Get vpml update ranges
  vpml_range = range(vpml_wall(i_vpml))

  ! Front wall
  if (loc_vpml == p_lower) then
    box_to_vpml = vpml_range(2)
    vpml_to_box = vpml_range(2)-1

  ! Back wall
  else
    box_to_vpml = vpml_range(1)
    vpml_to_box = vpml_range(1)+1

  endif

  bnd_field(1:p_x_dim,p_lower) = 0
  bnd_field(1:p_x_dim,p_upper) = field%nx_(1:p_x_dim) + 2

  ! Copy box values to VPML boundary
  i1_min = 1
  i1_max = field%nx_(1)+1

  i3_min = 1
  i3_max = field%nx_(3)+1

  do i3=i3_min, i3_max

    do i1=i1_min, i1_max

      vpml_wall(i_vpml)%f3(1, i1, box_to_vpml, i3) = field%f3(1, i1, box_to_vpml, i3)    ! XY component

      vpml_wall(i_vpml)%f3(2, i1, box_to_vpml, i3) = 0.0_p_k_fld                         ! XZ component

      vpml_wall(i_vpml)%f3(3, i1, box_to_vpml, i3) = field%f3(2, i1, box_to_vpml, i3)    ! YX component
      vpml_wall(i_vpml)%f3(4, i1, box_to_vpml, i3) = 0.0_p_k_fld                         ! YZ component

      vpml_wall(i_vpml)%f3(5, i1, box_to_vpml, i3) = field%f3(3, i1, box_to_vpml, i3)    ! ZX component
      vpml_wall(i_vpml)%f3(6, i1, box_to_vpml, i3) = 0.0_p_k_fld                         ! ZY component

    enddo
  enddo

  ! Copy VPML boundary values to box

  do i3=i3_min, i3_max

    do i1=i1_min, i1_max

      field%f3(1, i1, vpml_to_box, i3) =   vpml_wall(i_vpml)%f3(1, i1, vpml_to_box, i3) &
                                         + vpml_wall(i_vpml)%f3(2, i1, vpml_to_box, i3)

      field%f3(2, i1, vpml_to_box, i3) =   vpml_wall(i_vpml)%f3(3, i1, vpml_to_box, i3) &
                                         + vpml_wall(i_vpml)%f3(4, i1, vpml_to_box, i3)

      field%f3(3, i1, vpml_to_box, i3) =   vpml_wall(i_vpml)%f3(5, i1, vpml_to_box, i3) &
                                         + vpml_wall(i_vpml)%f3(6, i1, vpml_to_box, i3)

    enddo
  enddo

  ! Take care of x2 / x3 corners
  ! x1 / x2 corners were handled in update_interface_x1_3d

  i1_min = 1
  i1_max = vpml_wall(i_vpml)%nx(1)+1

  dir_idx = 3
  do bnd_idx=1,2

    ! adjacent vpml
    adj_vpml = pos_corner(i_vpml,dir_idx, bnd_idx)

    if (adj_vpml > 0) then

      bnd_vpml = range(vpml_wall(adj_vpml))

      ! bottom corner
      if (bnd_vpml(p_lower) < 0) then
         i_range(p_lower) = bnd_vpml(p_lower)
         i_range(p_upper) = bnd_vpml(p_upper)-1

      ! top corner
      else
         i_range(p_lower) = bnd_vpml(p_lower)+1
         i_range(p_upper) = bnd_vpml(p_upper)
      endif

      select case ( loc_vpml )
        case ( p_lower )
          do i3=i_range(p_lower), i_range(p_upper)
            do i1=i1_min, i1_max
              vpml_wall(adj_vpml)%f3(1, i1, 0, i3 ) = vpml_wall(i_vpml)%f3(1, i1, 0, i3 )

              vpml_wall(adj_vpml)%f3(2, i1, 0, i3 ) = vpml_wall(i_vpml)%f3(2, i1, 0, i3 )

              vpml_wall(adj_vpml)%f3(3, i1, 0, i3 ) = vpml_wall(i_vpml)%f3(3, i1, 0, i3 )
              vpml_wall(adj_vpml)%f3(4, i1, 0, i3 ) = vpml_wall(i_vpml)%f3(4, i1, 0, i3 )

              vpml_wall(adj_vpml)%f3(5, i1, 0, i3 ) = vpml_wall(i_vpml)%f3(5, i1, 0, i3 )
              vpml_wall(adj_vpml)%f3(6, i1, 0, i3 ) = vpml_wall(i_vpml)%f3(6, i1, 0, i3 )

            enddo

            do i1=i1_min, i1_max
              vpml_wall(i_vpml)%f3(1, i1, 1, i3 ) = vpml_wall(adj_vpml)%f3(1, i1, 1, i3 )

              vpml_wall(i_vpml)%f3(2, i1, 1, i3 ) = vpml_wall(adj_vpml)%f3(2, i1, 1, i3 )

              vpml_wall(i_vpml)%f3(3, i1, 1, i3 ) = vpml_wall(adj_vpml)%f3(3, i1, 1, i3 )
              vpml_wall(i_vpml)%f3(4, i1, 1, i3 ) = vpml_wall(adj_vpml)%f3(4, i1, 1, i3 )

              vpml_wall(i_vpml)%f3(5, i1, 1, i3 ) = vpml_wall(adj_vpml)%f3(5, i1, 1, i3 )
              vpml_wall(i_vpml)%f3(6, i1, 1, i3 ) = vpml_wall(adj_vpml)%f3(6, i1, 1, i3 )

            enddo
          enddo

        case ( p_upper )
          i2 = field%nx_(2)+1
          do i3=i_range(p_lower), i_range(p_upper)
            do i1=i1_min, i1_max
              vpml_wall(i_vpml)%f3(1, i1, i2, i3 ) = vpml_wall(adj_vpml)%f3(1, i1, i2, i3 )

              vpml_wall(i_vpml)%f3(2, i1, i2, i3 ) = vpml_wall(adj_vpml)%f3(2, i1, i2, i3 )

              vpml_wall(i_vpml)%f3(3, i1, i2, i3 ) = vpml_wall(adj_vpml)%f3(3, i1, i2, i3 )
              vpml_wall(i_vpml)%f3(4, i1, i2, i3 ) = vpml_wall(adj_vpml)%f3(4, i1, i2, i3 )

              vpml_wall(i_vpml)%f3(5, i1, i2, i3 ) = vpml_wall(adj_vpml)%f3(5, i1, i2, i3 )
              vpml_wall(i_vpml)%f3(6, i1, i2, i3 ) = vpml_wall(adj_vpml)%f3(6, i1, i2, i3 )

            enddo

            do i1=i1_min, i1_max
              vpml_wall(adj_vpml)%f3(1, i1, i2+1, i3 ) = vpml_wall(i_vpml)%f3(1, i1, i2+1, i3 )

              vpml_wall(adj_vpml)%f3(2, i1, i2+1, i3 ) = vpml_wall(i_vpml)%f3(2, i1, i2+1, i3 )

              vpml_wall(adj_vpml)%f3(3, i1, i2+1, i3 ) = vpml_wall(i_vpml)%f3(3, i1, i2+1, i3 )
              vpml_wall(adj_vpml)%f3(4, i1, i2+1, i3 ) = vpml_wall(i_vpml)%f3(4, i1, i2+1, i3 )

              vpml_wall(adj_vpml)%f3(5, i1, i2+1, i3 ) = vpml_wall(i_vpml)%f3(5, i1, i2+1, i3 )
              vpml_wall(adj_vpml)%f3(6, i1, i2+1, i3 ) = vpml_wall(i_vpml)%f3(6, i1, i2+1, i3 )

            enddo
          enddo
      end select

    endif

  enddo

end subroutine update_interface_x2_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine update_interface_x3_3d(vpml_wall, i_vpml, field, loc_vpml)
!
! Updates field at the interface
! Field components in 2D are
! 1 => XY
! 2 => XZ
!
! 3 => YX
! 4 => YZ
!
! 5 => ZX
! 6 => ZY
!-----------------------------------------------------------------------------------------
  implicit none

  type (t_wall), dimension(:), intent(inout) :: vpml_wall
  integer, intent(in) :: i_vpml                               ! vpml to update
  type(t_vdf),intent(inout) :: field
  integer, intent(in) :: loc_vpml

  ! local variables
  integer :: box_to_vpml, vpml_to_box

  integer :: i1_min, i1_max, i2_min, i2_max, i1, i2

  integer, dimension(2) :: vpml_range
  integer, dimension(p_x_dim,2) :: bnd_field

  ! Get vpml update ranges
  vpml_range = range(vpml_wall(i_vpml))

  ! Bottom wall
  if (loc_vpml == p_lower) then
    box_to_vpml = vpml_range(2)
    vpml_to_box = vpml_range(2)-1

  ! Top wall
  else
    box_to_vpml = vpml_range(1)
    vpml_to_box = vpml_range(1)+1

  endif

  bnd_field(1:p_x_dim,p_lower) = 0
  bnd_field(1:p_x_dim,p_upper) = field%nx_(1:p_x_dim) + 2

  i1_min = 1
  i1_max = field%nx_(1)+1

  i2_min = 1
  i2_max = field%nx_(2)+1

  ! Copy box values to VPML boundary

  do i2=i2_min, i2_max

    do i1=i1_min, i1_max

      vpml_wall(i_vpml)%f3(1, i1, i2, box_to_vpml) = field%f3(1, i1, i2, box_to_vpml)    ! XY component

      vpml_wall(i_vpml)%f3(2, i1, i2, box_to_vpml) = 0.0_p_k_fld                         ! XZ component

      vpml_wall(i_vpml)%f3(3, i1, i2, box_to_vpml) = field%f3(2, i1, i2, box_to_vpml)    ! YX component
      vpml_wall(i_vpml)%f3(4, i1, i2, box_to_vpml) = 0.0_p_k_fld                         ! YZ component

      vpml_wall(i_vpml)%f3(5, i1, i2, box_to_vpml) = field%f3(3, i1, i2, box_to_vpml)    ! ZX component
      vpml_wall(i_vpml)%f3(6, i1, i2, box_to_vpml) = 0.0_p_k_fld                         ! ZY component

    enddo
  enddo

  ! Copy VPML boundary values to box

  do i2=i2_min, i2_max

    do i1=i1_min, i1_max

      field%f3(1, i1, i2, vpml_to_box) = vpml_wall(i_vpml)%f3(1, i1, i2, vpml_to_box) &
                                         + vpml_wall(i_vpml)%f3(2, i1, i2, vpml_to_box)

      field%f3(2, i1, i2, vpml_to_box) = vpml_wall(i_vpml)%f3(3, i1, i2, vpml_to_box) &
                                         + vpml_wall(i_vpml)%f3(4, i1, i2, vpml_to_box)

      field%f3(3, i1, i2, vpml_to_box) = vpml_wall(i_vpml)%f3(5, i1, i2, vpml_to_box) &
                                         + vpml_wall(i_vpml)%f3(6, i1, i2, vpml_to_box)

    enddo
  enddo

  ! all corner values are handled either in

  ! update_interface_x1_3d or update_interface_x2_3d

end subroutine update_interface_x3_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reshape_vpml( this, old_lb, new_lb, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------
   implicit none

   type(t_vpml), intent(inout) :: this
   class( t_grid ), intent(in) :: old_lb, new_lb
   class( t_node_conf ), intent(in) :: no_co
   type(t_vdf_msg), dimension(:), intent(inout) :: send_msg, recv_msg

   integer :: i

   ! reshape wall objects
   do i=1, this%n_vpml
     call reshape_copy( this%wall_array_e(i), old_lb, new_lb, no_co, send_msg, recv_msg )
     call reshape_copy( this%wall_array_b(i), old_lb, new_lb, no_co, send_msg, recv_msg )
   enddo

end subroutine reshape_vpml
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Report E and B in VPML walls (debug only)
!-----------------------------------------------------------------------------------------
subroutine report_vpml( this, g_space, no_co, tstep, t, dx )
!-----------------------------------------------------------------------------------------

  use m_time_step

  implicit none

  type( t_vpml ),   intent(inout) :: this
  type( t_space ),     intent(in) :: g_space
  class( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double),      intent(in) :: t
  real(p_double), dimension(:), intent(in) :: dx

  integer :: i

  do i=1, this%n_vpml
    call report_wall( this%wall_array_e(i), this%wall_array_b(i), i, g_space, &
               tstep, t, dx )
  enddo

end subroutine report_vpml
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Compute coefficients for field push
! Coeficients of E and B in eq (27) of Vay's paper: JCP 165, 511 (2000)
! (p_x_dim, coefficient idx [1-3], vector idx)

!-----------------------------------------------------------------------------------------
subroutine setup_att_coeff(this)
!-----------------------------------------------------------------------------------------
  implicit none

  type (t_vpml), intent(inout) :: this

  ! local variables
  integer :: thickness, i, j
  real (p_k_fld) :: tj, tj12, tj1, beta1e, beta1b, dtdx, dxdtdif

  real (p_k_fld) :: beta2e, beta2b

  real (p_k_fld) :: s1, s2

  ! Allocate coefficient arrays
  ! see object definition for component details
  thickness = this%bnd_size
  call alloc(this%coef_e_low, (/ p_x_dim,3,thickness /))
  call alloc(this%coef_b_low, (/ p_x_dim,3,thickness /))
  call alloc(this%coef_e_up,  (/ p_x_dim,3,thickness /))
  call alloc(this%coef_b_up,  (/ p_x_dim,3,thickness /))

  ! Initialize
  this%coef_e_low = 0
  this%coef_b_low = 0

  this%coef_e_up  = 0
  this%coef_b_up  = 0

  ! This gives optimal absorption according to eq (7.61) in Taflove
  s1 = - (0.8_p_k_fld * ( p_att_nexp + 1 )) / 2.0_p_k_fld
  s2 = 1.0_p_k_fld / (thickness-3)

  ! This was the previous behavior (which is what you get from the above formulas
  ! for thickness = 8 and p_att_nexp = 4 )
  ! s1 = -2.0
  ! s2 = 0.2

  do i=1, p_x_dim

    ! use dx in corresponding direction
    dtdx = this%dtdx( i )

    ! (dx-dt)/(dx+dt) = (1-dt/dx)/(1+dt/dx)
    dxdtdif = (1.0_p_k_fld-dtdx)/(1.0_p_k_fld+dtdx)

    ! j=0 at the interface [after Eq. (38)]
    do j=0, thickness-1

      ! Attenuation coefficients for LOWER boundaries [Eqs (38-40)]:
      ! tj = exp[-2 * (abs(j)/5)^n]

      ! Check field positionings to obtain this (B is the outside point):
      ! j =  -0.5 => interface
      ! j =  -1.5 => one cell inside wall
      tj   = exp( s1 * ( s2 * ( j - 0.5_p_k_fld ) )**p_att_nexp )
      tj12 = exp( s1 * ( s2 * ( j ) )**p_att_nexp )
      tj1  = exp( s1 * ( s2 * ( j + 0.5_p_k_fld ) )**p_att_nexp )

      ! beta coefficients
      beta1e = dtdx*( 1.0_p_k_fld + dxdtdif * (1.0_p_k_fld - tj12) )
      beta2e = dtdx

      beta1b = dtdx*( 1.0_p_k_fld + dxdtdif * (1.0_p_k_fld - tj1 ) )
      beta2b = dtdx

      this%coef_e_low(i,1,thickness-j) = 1.0_p_k_fld - beta1e + tj12 * beta2e
      this%coef_e_low(i,2,thickness-j) = tj * beta1e
      this%coef_e_low(i,3,thickness-j) = beta2e

      this%coef_b_low(i,1,thickness-j) = 1.0_p_k_fld - beta1b + tj1 * beta2b
      this%coef_b_low(i,2,thickness-j) = tj12 * beta1b
      this%coef_b_low(i,3,thickness-j) = beta2b

      !Attenuation coefficients for UPPER boundaries
      ! E is the outside point
      ! j =  0 => interface
      ! j =  1 => one cell inside wall
      tj   = exp( s1 * ( s2 *  j )**p_att_nexp )
      tj12 = exp( s1 * ( s2 * (j + 0.5_p_k_fld) )**p_att_nexp )
      tj1  = exp( s1 * ( s2 * (j + 1.0_p_k_fld) )**p_att_nexp )

      ! beta coefficients
      beta1e = dtdx*( 1.0_p_k_fld + dxdtdif * (1.0_p_k_fld - tj12) )
      beta1b = dtdx*( 1.0_p_k_fld + dxdtdif * (1.0_p_k_fld - tj1 ) )

      this%coef_e_up(i,1,j+1) = 1.0_p_k_fld - beta1e + tj12 * dtdx
      this%coef_e_up(i,2,j+1) = dtdx
      this%coef_e_up(i,3,j+1) = tj * beta1e

      this%coef_b_up(i,1,j+1) = 1.0_p_k_fld - beta1b + tj1 * dtdx
      this%coef_b_up(i,2,j+1) = dtdx
      this%coef_b_up(i,3,j+1) = tj12 * beta1b

    enddo

  enddo

end subroutine setup_att_coeff
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_beam_vpml( this, ef, bf, no_co, send_msg, recv_msg )

  implicit none

  type (t_vpml), intent(inout) :: this
  type (t_vdf), intent(inout) :: ef, bf
  class( t_node_conf ),intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i
  integer, dimension( p_max_dim ) :: nx_move

  select case (p_x_dim)
    case(1)
      ! Nothing required in 1D
      continue

    case(2)
      ! Only x2 updates are required
      do i = 1, this%n_vpml
        if ( this%dir_vpml(i) == 2 ) then
           if (this%loc_vpml(i) == p_lower) then
              call update_beam_vpml_2d_lower( this%wall_array_b(i), this%wall_array_e(i), &
                                              this%coef_e_low(2,:,:), this%dtdx(1), &
                                              ef, bf )
           else
              call update_beam_vpml_2d_upper( this%wall_array_b(i), this%wall_array_e(i), &

                                              this%coef_e_up(2,:,:), this%dtdx(1), &
                                              ef, bf )
           endif
        endif
      enddo

    case(3)

      do i = 1, this%n_vpml
        select case (this%dir_vpml(i))

          case (1)
          ! we don't initialize the value in boundaries in x1 direction
          continue

          case (2)
            if (this%loc_vpml(i) == p_lower) then
              call update_beam_vpml_3d_x2_lower( this%wall_array_b(i), this%wall_array_e(i), &
                                                 this%coef_e_low(:,:,:), this%coef_e_up(:,:,:), &
                                                 this%coef_b_low(:,:,:), this%coef_b_up(:,:,:), &
                                                 this%dtdx(1), this%pos_corner(2,:,:), ef, bf)
            else
              call update_beam_vpml_3d_x2_upper( this%wall_array_b(i), this%wall_array_e(i), &
                                                 this%coef_e_low(:,:,:), this%coef_e_up(:,:,:), &
                                                 this%coef_b_low(:,:,:), this%coef_b_up(:,:,:), &
                                                 this%dtdx(1), this%pos_corner(2,:,:), ef, bf)
            endif

          case (3)
            if (this%loc_vpml(i) == p_lower) then
              call update_beam_vpml_3d_x3_lower( this%wall_array_b(i), this%wall_array_e(i), &
                                                 this%coef_e_low(:,:,:), this%coef_e_up(:,:,:), &
                                                 this%coef_b_low(:,:,:), this%coef_b_up(:,:,:), &
                                                 this%dtdx(1), ef, bf)
            else
              call update_beam_vpml_3d_x3_upper( this%wall_array_b(i), this%wall_array_e(i), &
                                                 this%coef_e_low(:,:,:), this%coef_e_up(:,:,:), &
                                                 this%coef_b_low(:,:,:), this%coef_b_up(:,:,:), &
                                                 this%dtdx(1), ef, bf)
            endif

        end select

      enddo

  end select

  ! Fix values in boundaries with other nodes
  nx_move = 0
  call update_boundary_vpml( this, no_co, nx_move, send_msg, recv_msg )

end subroutine update_beam_vpml

!-----------------------------------------------------------------------------------------
! Update beam fields, 2D, lower wall
!-----------------------------------------------------------------------------------------
subroutine update_beam_vpml_2d_lower( vpml_b, vpml_e, coef, dtdx, ef, bf )

  implicit none

  type (t_wall), intent(inout) :: vpml_b, vpml_e
  real(p_k_fld), dimension(:,:), intent(in) :: coef
  real(p_k_fld), intent(in) :: dtdx

  type( t_vdf ), intent(in) :: ef, bf

  integer, dimension(2)  :: range
  real( p_k_fld ) :: e0, b0
  integer :: i1, i2, coef_i

  call get_range( vpml_b, range )

  ! Update B field
  ! The magnetic field needs to be advanced half time step
  do i1 = 1, bf%nx_(1)
    b0 = bf%f2(3, i1, 1 ) - 0.5 * dtdx * ( ef%f2( 2, i1+1, 1 ) - ef%f2( 2, i1, 1 ))

    vpml_b%f2(3,i1,range(2)  ) = vpml_b%f2(3,i1,range(2)) + b0
    vpml_b%f2(3,i1,range(2)-1) = vpml_b%f2(3,i1,range(2)-1) + b0

    do i2 = range(2)-2, range(1), -1
       coef_i = i2 - range(1) + 2
       b0 = coef( 2,coef_i) / coef( 3, coef_i) * b0
       vpml_b%f2(3,i1,i2) = vpml_b%f2(3,i1,i2) + b0

    enddo

  enddo

  ! Update E field
  do i1 = 1, ef%nx_(1)
    e0 = ef%f2(2, i1, 1 )

    vpml_e%f2(2,i1,range(2))   = vpml_e%f2(2,i1,range(2))   + e0
    vpml_e%f2(2,i1,range(2)-1) = vpml_e%f2(2,i1,range(2)-1) + e0

    do i2 = range(2)-2, range(1), -1
       coef_i = i2 - range(1) + 2
       e0 = coef( 2,coef_i) / coef( 3, coef_i) * e0
       vpml_e%f2(2,i1,i2) = vpml_e%f2(2,i1,i2) + e0

    enddo

  enddo

end subroutine update_beam_vpml_2d_lower
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update beam fields, 2D, upper wall
! - This works for both cartesian and cylindrical geometries
! - Note that we are taking the field values on the guard cell; for cartesian geometry this
!   is the same value as the last cell inside the box, but for cylindrical geometry this
!   will be slightly lower 1 / r -> 1 / ( r + dr )

!-----------------------------------------------------------------------------------------
subroutine update_beam_vpml_2d_upper( vpml_b, vpml_e, coef, dtdx, ef, bf )

  implicit none

  type (t_wall), intent(inout) :: vpml_b, vpml_e
  real(p_k_fld), dimension(:,:), intent(in) :: coef
  real(p_k_fld), intent(in) :: dtdx

  type( t_vdf ), intent(in) :: ef, bf

  integer, dimension(2)  :: range
  real( p_k_fld ) :: e0, b0
  integer :: i1, i2, coef_i

  ! Update B field
  call get_range( vpml_b, range )

  ! The magnetic field needs to be advanced half time step
  do i1 = 1, bf%nx_(1)
    b0 = bf%f2(3, i1, bf%nx_(2)+1 ) - 0.5 * dtdx * ( ef%f2( 2, i1+1, bf%nx_(2)+1 ) - ef%f2( 2, i1, bf%nx_(2)+1))
    vpml_b%f2(3,i1,range(1)) = vpml_b%f2(3,i1,range(1)) + b0
    do i2 = range(1)+1, range(2)
       coef_i = i2 - range(1)

       b0 = b0 * coef( 3,coef_i) / coef( 2, coef_i)
       vpml_b%f2(3,i1,i2) =  vpml_b%f2(3,i1,i2) + b0
    enddo
  enddo

  ! Update E field
  do i1 = 1, ef%nx_(1)
    e0 = ef%f2(2, i1, ef%nx_(2)+1 )
    vpml_e%f2(2,i1,range(1)) = vpml_e%f2(2,i1,range(1)) + e0
    do i2 = range(1)+1, range(2)
       coef_i = i2 - range(1)

       e0 = e0 * coef( 3,coef_i) / coef( 2, coef_i)
       vpml_e%f2(2,i1,i2) = vpml_e%f2(2,i1,i2) + e0
    enddo
  enddo

end subroutine update_beam_vpml_2d_upper
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Module Overview: 

! The modules update_beam_vpml_3d (x2/x3, lower/uppper) aim to update (initialize)
! the 3D fields setup inside the VPML walls, after the beam fields are initialized in the box.
! It shows that the current implementation works well for long-beam simulations with 
! disruption parameter D < 10

! For more information about these modules, please contact wenlong.zhang@tecnico.ulisboa.pt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update beam fields, 3D, lower wall, x2 direction (x2-VPML, lower wall)
!-----------------------------------------------------------------------------------------
subroutine update_beam_vpml_3d_x2_lower( vpml_b, vpml_e, coef_e_low, coef_e_up, coef_b_low, coef_b_up, dtdx, pos_corner, ef, bf )

  implicit none

  type (t_wall), intent(inout) :: vpml_b, vpml_e
  real(p_k_fld), dimension(:,:,:), intent(in) ::  coef_e_low, coef_e_up, coef_b_low, coef_b_up
  real(p_k_fld), intent(in) :: dtdx
  type(t_vdf), intent(in) :: ef, bf

  integer, dimension(2) :: range
  logical, dimension(2) :: if_cnr_x3
  integer, dimension(:,:), intent(in) :: pos_corner
  real(p_k_fld), dimension(3) :: coef_e_cnr_x3, coef_b_cnr_x3
  real(p_k_fld) :: e0_3, e0_5, b0_3, b0_5

  integer :: i1, i2, i3, coef_i
  integer :: gc1a, gc1b, gc2a, gc2b, gc3a, gc3b, i2_connecting, i3_lowercorner, i3_uppercorner, coef_i_cnr

  ! Update B fields
  call get_range( vpml_b, range )

  gc1a = vpml_e%gc_num(p_lower, 1)
  gc1b = vpml_e%gc_num(p_upper, 1)

  gc2a = vpml_e%gc_num(p_lower, 2)
  gc2b = vpml_e%gc_num(p_upper, 2)

  gc3a = vpml_e%gc_num(p_lower, 3)
  gc3b = vpml_e%gc_num(p_upper, 3)

  if_cnr_x3(p_lower) = (pos_corner(3, p_lower) > 0)
  if_cnr_x3(p_upper) = (pos_corner(3, p_upper) > 0)

  ! Boundary indices of the lower and upper corners
  i3_lowercorner = 0
  i3_uppercorner = vpml_b%nx(3) + 2


  ! x2-lowerVPML: field setup in the conventional VPML grids (not in the corner regions)

  ! Bzx
  do i1 = 1, ef%nx_(1)
    do i3 = 1, ef%nx_(3) + 1

       ! Index connecting the simulation box and the VPML wall region
       i2_connecting = 0

       ! Advance the B-field with one half time step
       b0_5 = bf%f3(3, i1, i2_connecting, i3 ) - 0.5 * dtdx * ( ef%f3( 2, i1+1, i2_connecting, i3 ) - ef%f3( 2, i1, i2_connecting, i3 ))

       ! Assign the outermost VPML grids
       vpml_b%f3(5, i1, range(2), i3  ) = vpml_b%f3(5, i1, range(2), i3  ) + b0_5
       vpml_b%f3(5, i1, range(2)-1, i3) = vpml_b%f3(5, i1, range(2)-1, i3) + b0_5

       ! Set the damped profile
       do i2 = range(2) - 2, range(1), -1
          coef_i = i2 - range(1) + 2

          b0_5 = b0_5 * coef_e_low(2, 2, coef_i) / coef_e_low(2, 3, coef_i)
          vpml_b%f3(5, i1, i2, i3) = vpml_b%f3(5, i1, i2, i3) + b0_5

       enddo

    enddo
  enddo    


  ! Ezx
  do i1 = 1, ef%nx_(1)
    do i3 = lbound( ef%f3 , 4 ), ubound( ef%f3, 4 )

       i2_connecting = 0

       e0_5 = ef%f3(3, i1, i2_connecting, i3 )
       vpml_e%f3(5, i1, range(2), i3  ) = vpml_e%f3(5, i1, range(2), i3  ) + e0_5
       vpml_e%f3(5, i1, range(2)-1, i3) = vpml_e%f3(5, i1, range(2)-1, i3) + e0_5

       do i2 = range(2) - 2, range(1), -1
          coef_i = i2 - range(1) + 2

          e0_5 = e0_5 * coef_b_low(2, 2, coef_i) / coef_b_low(2, 3, coef_i)
          vpml_e%f3(5, i1, i2, i3) = vpml_e%f3(5, i1, i2, i3) + e0_5

       enddo

    enddo
  enddo


  ! Eyx
  do i1 = 1, ef%nx_(1)
    do i3 = lbound( ef%f3 , 4 ), ubound( ef%f3, 4 )

       i2_connecting = 0

       e0_3 = ef%f3(2, i1, i2_connecting, i3 )
       vpml_e%f3(3, i1, range(2), i3  ) = vpml_e%f3(3, i1, range(2), i3  ) + e0_3
       vpml_e%f3(3, i1, range(2)-1, i3) = vpml_e%f3(3, i1, range(2)-1, i3) + e0_3

       do i2 = range(2) - 2, range(1), -1
          coef_i = i2 - range(1) + 2

          e0_3 = e0_3 * coef_b_low(2, 2, coef_i) / coef_b_low(2, 3, coef_i)
          vpml_e%f3(3, i1, i2, i3) = vpml_e%f3(3, i1, i2, i3) + e0_3

       enddo

    enddo
  enddo


  ! x2-lowerVPML-Lower-Corner

  if (i3_lowercorner < lbound( ef%f3 , 4 )) then
      print *, 'Warning: There is a gap of index between the x2-lowerVPML-Lower-Corner and the lower boundary of simulation fields in the box'
      ERROR('Error: x2-lowerVPML-Lower-Corner boundary index < lower boundary of simulation fields in the box')
      ! If there is such an index gap, the initialization of the fields (at index of i3_lowercorner) is needed
  endif


  ! Set Eyx in the corner

  ! Condition to check if we are in the lower corner
  if ( (1 - gc3a) < lbound( ef%f3 , 4 ) ) then

     do i1 = 1, ef%nx_(1)
       do i2 = range(2), range(1), -1
         do i3 = i3_lowercorner - 1, 1 - gc3a, -1

            coef_i_cnr = i3 + gc3a + 1
            coef_e_cnr_x3 = coef_e_low(3, :, coef_i_cnr)
            coef_b_cnr_x3 = coef_b_low(3, :, coef_i_cnr)

            if ( (coef_e_cnr_x3(2) == 0.0) .or. (coef_e_cnr_x3(3) == 0.0) .or. (coef_b_cnr_x3(2) == 0.0) .or. (coef_b_cnr_x3(3) == 0.0)) then
               ERROR('Error: Coefficients for x2-lowerVPML-Lower-Corner are 0, please check the index')
            endif

            vpml_e%f3(3, i1, i2, i3) = coef_b_cnr_x3(2) / coef_b_cnr_x3(3) * vpml_e%f3(3, i1, i2, i3+1)

         enddo
       enddo
     enddo

  endif


  ! Set Ezx in the corner

  if ( (1 - gc3a) < lbound( ef%f3 , 4 ) ) then

     ! x3-damp
     do i1 = 1, ef%nx_(1)
       do i3 = i3_lowercorner - 1, 1 - gc3a, -1
          coef_i = i3 + gc3a + 1

          vpml_e%f3(5, i1, range(2), i3  ) = coef_b_low(3, 2, coef_i) / coef_b_low(3, 3, coef_i) * vpml_e%f3(5, i1, range(2), i3+1  )
          vpml_e%f3(5, i1, range(2)-1, i3) = coef_b_low(3, 2, coef_i) / coef_b_low(3, 3, coef_i) * vpml_e%f3(5, i1, range(2)-1, i3+1)

       enddo
     enddo

     ! x2-damp
     do i1 = 1, ef%nx_(1)
       do i3 = i3_lowercorner, 1 - gc3a, -1

          do i2 = range(2) - 2, range(1), -1
             coef_i = i2 - range(1) + 2

             vpml_e%f3(5, i1, i2, i3) = coef_b_low(2, 2, coef_i) / coef_b_low(2, 3, coef_i) * vpml_e%f3(5, i1, i2+1, i3)

          enddo

       enddo
     enddo

  endif


  ! x2-lowerVPML-Upper-Corner

  if (i3_uppercorner > ubound( ef%f3 , 4 )) then
     print *, 'Warning: There is an index gap between the x2-lowerVPML-Upper-Corner and the upper boundary of simulation fields in the box'
     ERROR('Error: x2-lowerVPML-Upper-Corner boundary index > upper boundary of simulation fields in the box')
     !  If there is such an index gap, the initialization of the fields (at index of i3_uppercorner) is needed
  endif

  ! Set Eyx in the corner

  if ( (vpml_e%nx(3) + gc3b) > ubound( ef%f3 , 4 ) ) then

     do i1 = 1, ef%nx_(1)
       do i2 = range(2), range(1), -1
         do i3 = i3_uppercorner + 1, vpml_e%nx(3) + gc3b

            coef_i_cnr = i3 - vpml_e%nx(3) - 2
            coef_e_cnr_x3 = coef_e_up(3, :, coef_i_cnr)
            coef_b_cnr_x3 = coef_b_up(3, :, coef_i_cnr)

            if ( (coef_e_cnr_x3(2) == 0.0) .or. (coef_e_cnr_x3(3) == 0.0) .or. (coef_b_cnr_x3(2) == 0.0) .or. (coef_b_cnr_x3(3) == 0.0)) then
               ERROR('Error: Coefficients for x2-lowerVPML-Upper-Corner are 0, please check the index')
            endif

            vpml_e%f3(3, i1, i2, i3) = coef_b_cnr_x3(3) / coef_b_cnr_x3(2) * vpml_e%f3(3, i1, i2, i3-1)

         enddo
       enddo
     enddo

  endif


  ! Set Ezx in the corner

  if ( (vpml_e%nx(3) + gc3b) > ubound( ef%f3 , 4 ) ) then

     ! x3-damp
     do i1 = 1, ef%nx_(1)
       do i3 = i3_uppercorner + 1, vpml_e%nx(3) + gc3b
          coef_i = i3 - vpml_e%nx(3) - 1

          vpml_e%f3(5, i1, range(2), i3  ) = coef_b_up(3, 3, coef_i) / coef_b_up(3, 2, coef_i) * vpml_e%f3(5, i1, range(2), i3-1  )
          vpml_e%f3(5, i1, range(2)-1, i3) = coef_b_up(3, 3, coef_i) / coef_b_up(3, 2, coef_i) * vpml_e%f3(5, i1, range(2)-1, i3-1)

       enddo
     enddo

     ! x2-damp
     do i1 = 1, ef%nx_(1)
       do i3 = i3_uppercorner, vpml_e%nx(3) + gc3b

          do i2 = range(2) - 2, range(1), -1
             coef_i = i2 - range(1) + 2

             vpml_e%f3(5, i1, i2, i3) = coef_b_low(2, 2, coef_i) / coef_b_low(2, 3, coef_i) * vpml_e%f3(5, i1, i2+1, i3)

          enddo

       enddo
     enddo

  endif

end subroutine update_beam_vpml_3d_x2_lower
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update beam fields, 3D, upper wall, x2 direction (x2-VPML, upper wall)
!-----------------------------------------------------------------------------------------
subroutine update_beam_vpml_3d_x2_upper( vpml_b, vpml_e, coef_e_low, coef_e_up, coef_b_low, coef_b_up, dtdx, pos_corner, ef, bf )

  implicit none

  type (t_wall), intent(inout) :: vpml_b, vpml_e
  real(p_k_fld), dimension(:,:,:), intent(in) :: coef_e_low, coef_e_up, coef_b_low, coef_b_up
  real(p_k_fld), intent(in) :: dtdx

  type( t_vdf ), intent(in) :: ef, bf

  integer, dimension(2)  :: range
  logical, dimension(2) :: if_cnr_x3
  integer, dimension(:,:), intent(in) :: pos_corner
  real(p_k_fld), dimension(3) :: coef_e_cnr_x3, coef_b_cnr_x3
  real( p_k_fld ) :: e0_3, e0_5, b0_3, b0_5

  integer :: i1, i2, i3, coef_i
  integer :: gc1a, gc1b, gc2a, gc2b, gc3a, gc3b, i2_connecting, i3_lowercorner, i3_uppercorner, coef_i_cnr

  ! Update B field
  call get_range( vpml_b, range )

  gc1a = vpml_e%gc_num(p_lower, 1)
  gc1b = vpml_e%gc_num(p_upper, 1)

  gc2a = vpml_e%gc_num(p_lower, 2)
  gc2b = vpml_e%gc_num(p_upper, 2)

  gc3a = vpml_e%gc_num(p_lower, 3)
  gc3b = vpml_e%gc_num(p_upper, 3)

  if_cnr_x3(p_lower) = (pos_corner(3, p_lower) > 0)
  if_cnr_x3(p_upper) = (pos_corner(3, p_upper) > 0)

  ! Boundary indices of the lower and upper corners
  i3_lowercorner = 0
  i3_uppercorner = vpml_b%nx(3) + 2

  ! x2-upperVPML: field setup in the conventional VPML grids (not in the corner regions)

  ! Bzx
  do i1 = 1, ef%nx_(1)
    do i3 = 1, ef%nx_(3) + 1

       i2_connecting = bf%nx_(2) + 1

       b0_5 = bf%f3(3, i1, i2_connecting, i3) - 0.5 * dtdx * ( ef%f3( 2, i1+1, i2_connecting, i3 ) - ef%f3( 2, i1, i2_connecting, i3))
       vpml_b%f3(5, i1, range(1), i3  ) = vpml_b%f3(5, i1, range(1), i3  ) + b0_5
       vpml_b%f3(5, i1, range(1)+1, i3) = vpml_b%f3(5, i1, range(1)+1, i3) + b0_5

       do i2 = range(1) + 2, range(2)
          coef_i = i2 - range(1)

          b0_5 = b0_5 * coef_e_up(2, 3, coef_i) / coef_e_up(2, 2, coef_i)
          vpml_b%f3(5, i1, i2, i3) = vpml_b%f3(5, i1, i2, i3) + b0_5

       enddo

    enddo
  enddo

  ! Ezx
  do i1 = 1, ef%nx_(1)
    do i3 = lbound( ef%f3 , 4 ), ubound( ef%f3, 4 )

       i2_connecting = bf%nx_(2) + 1
       e0_5 = ef%f3(3, i1, i2_connecting ,i3)
       vpml_e%f3(5, i1, range(1), i3  ) = vpml_e%f3(5, i1, range(1), i3  ) + e0_5
       vpml_e%f3(5, i1, range(1)+1, i3) = vpml_e%f3(5, i1, range(1)+1, i3) + e0_5

       do i2 = range(1) + 2, range(2)
          coef_i = i2 - range(1) - 1

          e0_5 = e0_5 * coef_b_up(2, 3, coef_i) / coef_b_up(2, 2, coef_i)
          vpml_e%f3(5, i1, i2, i3) = vpml_e%f3(5, i1, i2, i3) + e0_5

       enddo

    enddo
  enddo

  ! Eyx
  do i1 = 1, ef%nx_(1)
    do i3 = lbound( ef%f3 , 4 ), ubound( ef%f3, 4 )

       i2_connecting = bf%nx_(2) + 1
       e0_3 = ef%f3(2, i1, i2_connecting, i3)
       vpml_e%f3(3, i1, range(1), i3  ) = vpml_e%f3(3, i1, range(1), i3  ) + e0_3
       vpml_e%f3(3, i1, range(1)+1, i3) = vpml_e%f3(3, i1, range(1)+1, i3) + e0_3

       do i2 = range(1) + 2, range(2)
          coef_i = i2 - range(1)

          e0_3 = e0_3 * coef_b_up(2, 3, coef_i) / coef_b_up(2, 2, coef_i)
          vpml_e%f3(3, i1, i2, i3) = vpml_e%f3(3, i1, i2, i3) + e0_3

       enddo

    enddo
  enddo


  ! x2-upperVPML-Lower-Corner

  if (i3_lowercorner < lbound( ef%f3 , 4 )) then
      print *, 'Warning: There is a gap of index between the x2-upperVPML-Lower-Corner and the lower boundary of simulation fields in the box'
      ERROR('Error: x2-upperVPML-Lower-Corner boundary index < lower boundary of simulation fields in the box')
      !  If there is such an index gap, the initialization of the fields (at index of i3_lowercorner) is needed
  endif

  ! Set Eyx in the corner
  if ( (1 - gc3a) < lbound( ef%f3 , 4 ) ) then

     do i1 = 1, ef%nx_(1)
       do i2 = range(1), range(2)
         do i3 = i3_lowercorner - 1, 1 - gc3a, -1

            coef_i_cnr = i3 + gc3a + 1
            coef_e_cnr_x3 = coef_e_low(3, :, coef_i_cnr)
            coef_b_cnr_x3 = coef_b_low(3, :, coef_i_cnr)

            if ( (coef_e_cnr_x3(2) == 0.0) .or. (coef_e_cnr_x3(3) == 0.0) .or. (coef_b_cnr_x3(2) == 0.0) .or. (coef_b_cnr_x3(3) == 0.0)) then
               ERROR('Error: Coefficients for x2-upperVPML-Lower-Corner are 0, please check the index')
            endif

            vpml_e%f3(3, i1, i2, i3) = coef_b_cnr_x3(2) / coef_b_cnr_x3(3) * vpml_e%f3(3, i1, i2, i3+1)

         enddo
       enddo
     enddo

  endif

  ! Set Ezx in the corner
  if ( (1 - gc3a) < lbound( ef%f3 , 4 ) ) then

     ! x3-damp
     do i1 = 1, ef%nx_(1)
       do i3 = i3_lowercorner - 1, 1 - gc3a, -1
          coef_i = i3 + gc3a + 1

          vpml_e%f3(5, i1, range(1), i3  ) = coef_b_low(3, 2, coef_i) / coef_b_low(3, 3, coef_i) * vpml_e%f3(5, i1, range(1), i3+1  )
          vpml_e%f3(5, i1, range(1)+1, i3) = coef_b_low(3, 2, coef_i) / coef_b_low(3, 3, coef_i) * vpml_e%f3(5, i1, range(1)+1, i3+1)

       enddo
     enddo

     ! x2-damp
     do i1 = 1, ef%nx_(1)
       do i3 = i3_lowercorner, 1 - gc3a, -1

          do i2 = range(1) + 2, range(2)
             coef_i = i2 - range(1) - 1

             vpml_e%f3(5, i1, i2, i3) = coef_b_up(2, 3, coef_i) / coef_b_up(2, 2, coef_i) * vpml_e%f3(5, i1, i2-1, i3)

          enddo

       enddo
     enddo

  endif
  

  ! x2-upperVPML-Upper-Corner

  if (i3_uppercorner > ubound( ef%f3 , 4 )) then
     print *, 'Warning: There is an index gap between the x2-upperVPML-Upper-Corner and the upper boundary of simulation fields in the box'
     ERROR('Error: x2-upperVPML-Upper-Corner boundary index > upper boundary of simulation fields in the box')
  endif

  ! Set Eyx in the corner
  if ( (vpml_e%nx(3) + gc3b) > ubound( ef%f3 , 4 ) ) then

     do i1 = 1, ef%nx_(1)
       do i2 = range(1), range(2)
         do i3 = i3_uppercorner + 1, vpml_e%nx(3) + gc3b

            coef_i_cnr = i3 - vpml_e%nx(3) - 2
            coef_e_cnr_x3 = coef_e_up(3, :, coef_i_cnr)
            coef_b_cnr_x3 = coef_b_up(3, :, coef_i_cnr)

            if ( (coef_e_cnr_x3(2) == 0.0) .or. (coef_e_cnr_x3(3) == 0.0) .or. (coef_b_cnr_x3(2) == 0.0) .or. (coef_b_cnr_x3(3) == 0.0)) then
               ERROR('Error: Coefficients for x2-upperVPML-Upper-Corner are 0, please check the index')
            endif

            vpml_e%f3(3, i1, i2, i3) = coef_b_cnr_x3(3) / coef_b_cnr_x3(2) * vpml_e%f3(3, i1, i2, i3-1)

         enddo
       enddo
     enddo

  endif
  

  ! Set Ezx in the corner
  if ( (vpml_e%nx(3) + gc3b) > ubound( ef%f3 , 4 ) ) then

     ! x3-damp
     do i1 = 1, ef%nx_(1)
       do i3 = i3_uppercorner + 1, vpml_e%nx(3) + gc3b
          coef_i = i3 - vpml_e%nx(3) - 1

          vpml_e%f3(5, i1, range(1), i3  ) = coef_b_up(3, 3, coef_i) / coef_b_up(3, 2, coef_i) * vpml_e%f3(5, i1, range(1), i3-1  )
          vpml_e%f3(5, i1, range(1)+1, i3) = coef_b_up(3, 3, coef_i) / coef_b_up(3, 2, coef_i) * vpml_e%f3(5, i1, range(1)+1, i3-1)

       enddo
     enddo

     ! x2-damp
     do i1 = 1, ef%nx_(1)
       do i3 = i3_uppercorner, vpml_e%nx(3) + gc3b

          do i2 = range(1) + 2, range(2)
             coef_i = i2 - range(1) - 1

             vpml_e%f3(5, i1, i2, i3) = coef_b_up(2, 3, coef_i) / coef_b_up(2, 2, coef_i) * vpml_e%f3(5, i1, i2-1, i3)

          enddo

       enddo
     enddo

  endif


end subroutine update_beam_vpml_3d_x2_upper  
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update beam fields, 3D, lower wall, x3 direction (x3-VPML, lower wall)
!-----------------------------------------------------------------------------------------
subroutine update_beam_vpml_3d_x3_lower( vpml_b, vpml_e, coef_e_low, coef_e_up, coef_b_low, coef_b_up, dtdx, ef, bf )

  implicit none

  type (t_wall), intent(inout) :: vpml_b, vpml_e
  real(p_k_fld), dimension(:,:,:), intent(in) :: coef_e_low, coef_e_up, coef_b_low, coef_b_up
  real(p_k_fld), intent(in) :: dtdx
  type( t_vdf ), intent(in) :: ef, bf

  integer, dimension(2)  :: range
  real( p_k_fld ) :: e0_3, e0_5, b0_3, b0_5
  integer :: i1, i2, i3, coef_i

  integer :: i3_connecting

  call get_range( vpml_b, range )

  ! x3-lowerVPML: field setup in the VPML grids. There are no corners in x3-VPML

  ! Byx
  do i1 = 1, bf%nx_(1)
    do i2 = lbound( ef%f3 , 3 ), ubound( ef%f3, 3 )

       i3_connecting = 0

       b0_3 = bf%f3(2, i1, i2, i3_connecting ) + 0.5 * dtdx * ( ef%f3( 3, i1+1, i2, i3_connecting ) - ef%f3( 3, i1, i2, i3_connecting ))
       vpml_b%f3(3, i1, i2, range(2)  ) = vpml_b%f3(3, i1, i2, range(2)  ) + b0_3
       vpml_b%f3(3, i1, i2, range(2)-1) = vpml_b%f3(3, i1, i2, range(2)-1) + b0_3

       do i3 = range(2) - 2, range(1), -1
          coef_i = i3 - range(1) + 2

          b0_3 = b0_3 * coef_e_low(3, 2, coef_i) / coef_e_low(3, 3, coef_i)
          vpml_b%f3(3, i1, i2, i3) = vpml_b%f3(3, i1, i2, i3) + b0_3

       enddo

    enddo
  enddo  

  ! Eyx
  do i1 = 1, ef%nx_(1)
    do i2 = lbound( ef%f3 , 3 ), ubound( ef%f3, 3 )

       i3_connecting = 0
       e0_3 = ef%f3(2, i1, i2, i3_connecting )
       vpml_e%f3(3, i1, i2, range(2)  ) = vpml_e%f3(3, i1, i2, range(2)  ) + e0_3
       vpml_e%f3(3, i1, i2, range(2)-1) = vpml_e%f3(3, i1, i2, range(2)-1) + e0_3

       do i3 = range(2) - 2, range(1), -1
          coef_i = i3 - range(1) + 2

          e0_3 = e0_3 * coef_b_low(3, 2, coef_i) / coef_b_low(3, 3, coef_i)
          vpml_e%f3(3, i1, i2, i3) = vpml_e%f3(3, i1, i2, i3) + e0_3

       enddo

    enddo
  enddo


  ! Ezx
  do i1 = 1, bf%nx_(1)
    do i2 = lbound( ef%f3 , 3 ), ubound( ef%f3, 3 )

       i3_connecting = 0

       e0_5 = ef%f3(3, i1, i2, i3_connecting )
       vpml_e%f3(5, i1, i2, range(2)  ) = vpml_e%f3(5, i1, i2, range(2)  ) + e0_5
       vpml_e%f3(5, i1, i2, range(2)-1) = vpml_e%f3(5, i1, i2, range(2)-1) + e0_5

       do i3 = range(2) - 2, range(1), -1
          coef_i = i3 - range(1) + 2

          e0_5 = e0_5 * coef_b_low(3, 2, coef_i) / coef_b_low(3, 3, coef_i)
          vpml_e%f3(5, i1, i2, i3) = vpml_e%f3(5, i1, i2, i3) + e0_5

       enddo

    enddo
  enddo

end subroutine update_beam_vpml_3d_x3_lower
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update beam fields, 3D, upper wall, x3 direction (x3-VPML, upper wall)
!-----------------------------------------------------------------------------------------
subroutine update_beam_vpml_3d_x3_upper( vpml_b, vpml_e, coef_e_low, coef_e_up, coef_b_low, coef_b_up, dtdx, ef, bf )

  implicit none

  type (t_wall), intent(inout) :: vpml_b, vpml_e
  real(p_k_fld), dimension(:,:,:), intent(in) :: coef_e_low, coef_e_up, coef_b_low, coef_b_up
  real(p_k_fld), intent(in) :: dtdx

  type( t_vdf ), intent(in) :: ef, bf

  integer, dimension(2)  :: range
  real( p_k_fld ) :: e0_3, e0_5, b0_3, b0_5
  integer :: i1, i2, i3, coef_i

  integer :: i3_connecting

  ! Update B field
  call get_range( vpml_b, range )

  ! x3-upperVPML: field setup in the VPML grids. There are no corners in x3-VPML

  ! Byx
  do i1 = 1, bf%nx_(1)
    do i2 =  lbound( ef%f3 , 3 ), ubound( ef%f3, 3 )

       i3_connecting = bf%nx_(3)+1

       b0_3 = bf%f3(2, i1, i2, i3_connecting) + 0.5 * dtdx * ( ef%f3( 3, i1+1, i2, i3_connecting ) - ef%f3( 3, i1, i2, i3_connecting))
       vpml_b%f3(3, i1, i2, range(1)  ) = vpml_b%f3(3, i1, i2, range(1)  ) + b0_3
       vpml_b%f3(3, i1, i2, range(1)+1) = vpml_b%f3(3, i1, i2, range(1)+1) + b0_3

       do i3 = range(1) + 2, range(2)
          coef_i = i3 - range(1)

          b0_3 = b0_3 * coef_e_up(3, 3, coef_i) / coef_e_up(3, 2, coef_i)
          vpml_b%f3(3, i1, i2, i3) =  vpml_b%f3(3, i1, i2, i3) + b0_3

       enddo

    enddo
  enddo

  ! Eyx
  do i1 = 1, ef%nx_(1)
    do i2 = lbound( ef%f3 , 3 ), ubound( ef%f3, 3 )

       i3_connecting = bf%nx_(3)+1
       e0_3 = ef%f3(2, i1, i2, i3_connecting)
       vpml_e%f3(3, i1, i2, range(1)  ) = vpml_e%f3(3, i1, i2, range(1)  ) + e0_3
       vpml_e%f3(3, i1, i2, range(1)+1) = vpml_e%f3(3, i1, i2, range(1)+1) + e0_3

       do i3 = range(1) + 2, range(2)
          coef_i = i3 - range(1) - 1

          e0_3 = e0_3 * coef_b_up(3, 3, coef_i) / coef_b_up(3, 2, coef_i)
          vpml_e%f3(3, i1, i2, i3) = vpml_e%f3(3, i1, i2, i3) + e0_3

       enddo

    enddo
  enddo

  ! Ezx
  do i1 = 1, bf%nx_(1)
    do i2 =  lbound( ef%f3 , 3 ), ubound( ef%f3, 3 )

       i3_connecting = bf%nx_(3)+1

       e0_5 = ef%f3(3, i1, i2, i3_connecting)
       vpml_e%f3(5, i1, i2, range(1)  ) = vpml_e%f3(5, i1, i2, range(1)  ) + e0_5
       vpml_e%f3(5, i1, i2, range(1)+1) = vpml_e%f3(5, i1, i2, range(1)+1) + e0_5

       do i3 = range(1) + 2, range(2)
          coef_i = i3 - range(1)

          e0_5 = e0_5 * coef_b_up(3, 3, coef_i) / coef_b_up(3, 2, coef_i)
          vpml_e%f3(5, i1, i2, i3) = vpml_e%f3(5, i1, i2, i3) + e0_5

       enddo

    enddo
  enddo
  
end subroutine update_beam_vpml_3d_x3_upper
!-----------------------------------------------------------------------------------------

end module m_vpml

#endif
