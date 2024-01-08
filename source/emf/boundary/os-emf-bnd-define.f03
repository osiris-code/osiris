
! string to id restart data
character(len=*), parameter :: p_bcemf_rst_id = "bcemf rst data - 0x0002"

!-------------------------------------------------------------------------------
! t_lindman class definition
!-------------------------------------------------------------------------------

type:: t_lindman

    integer :: idir = -1                     ! direction of lindman bc (1->p_x_dim)
    integer :: orientation                   ! upper or lower wall
    integer, dimension(p_f_dim-1)  :: perp   ! perpendicular field directions
    integer, dimension(2) :: range                   ! range of buffers

    type (t_wall) :: e_buffer,b_buffer               ! e and b buffers for storing the
                                                    ! fields at previous timesteps

    real (p_k_fld), dimension(p_max_dim) :: dtdxi  ! (dt/dx/2)
    real (p_k_fld) :: disp_corr                      ! dispersion correction
                                                    ! (1-dx/(dt*c))/(1+dx/(dt*c))

    integer , dimension(p_max_dim) :: imin, imax

    integer , dimension(2) :: bound

end type t_lindman

!-------------------------------------------------------------------------------
! t_vpml class definition
!-------------------------------------------------------------------------------

type:: t_vpml

    ! Simulation coordinates (cartesian/r-z cylindrical)    
    integer :: coordinates

    ! Position of lower axial boundary for r-z cylindrical simulations
    integer :: gir_pos

    ! Boundaries are sorted as follows
    ! x1_lower, x1_upper, x2_lower, x2_upper, x3_lower, x3_upper

    ! Number of vpml boundaries
    integer :: n_vpml = 0

    ! Current vpml when looping over update boundary with tiles
    integer :: c

    ! Perp. direction of the boundary
    integer, dimension(:), pointer  :: dir_vpml => null()

    ! Boundary location (p_lower, p_upper)
    integer, dimension(:), pointer  :: loc_vpml => null()

    ! Number of cells in vpml boundaries
    integer :: bnd_size

    logical, dimension(2, p_x_dim) :: if_diffuse

    ! dt/dx
    real(p_k_fld), dimension(p_x_dim) :: dtdx

    ! Contact walls for each corner
    integer, dimension(:,:,:), pointer :: pos_corner => null()

    ! Coeficients of E and B in eq (27) of Vay's paper: JCP 165, 511 (2000)
    ! [
    !  direction  => wall direction (p_x_dim elements)
    !  coeff. idx => index of coefficients in eq. (27) of Vay's paper (3 elements)
    !  array pos. => position in the vpml boundary (bnd_size elements)
    ! ]

    real(p_k_fld), dimension(:,:,:), pointer  :: coef_e_low  => null()

    real(p_k_fld), dimension(:,:,:), pointer  :: coef_b_low  => null()

    real(p_k_fld), dimension(:,:,:), pointer  :: coef_e_up   => null()

    real(p_k_fld), dimension(:,:,:), pointer  :: coef_b_up   => null()

    ! 2 walls per vpml boundary
    type (t_wall), pointer, dimension(:) :: wall_array_e => null()

    type (t_wall), pointer, dimension(:) :: wall_array_b => null()

end type t_vpml

!-------------------------------------------------------------------------------
! t_emf_bound class definition
!-------------------------------------------------------------------------------

type :: t_emf_bound

    ! type of boundary conditions for each boundary

    !indexes: (lower-upper bound, direction)
    integer, dimension(2,p_max_dim) :: type

    integer, dimension(2,p_max_dim) :: global_type

    ! lindman boundary data
    type (t_lindman), dimension(2):: lindman

    ! pml boundary data (all boundaries)
    type (t_vpml) :: vpml_all
    integer :: vpml_bnd_size
    logical, dimension(2,p_x_dim) :: vpml_diffuse

contains

    procedure :: read_input         => read_input_emf_bound
    procedure :: init               => init_emf_bound 
    procedure :: cleanup            => cleanup_emf_bound
    procedure :: write_checkpoint   => write_checkpoint_emf_bound
    procedure :: restart_read       => restart_read_emf_bound
    procedure :: update_boundary    => update_boundary_emf_bound
    procedure :: update_boundary_b  => update_boundary_bfld_emf_bound
    procedure :: update_boundary_e  => update_boundary_efld_emf_bound
    procedure :: move_window        => move_window_emf_bound
    procedure :: list_algorithm     => list_algorithm_emf_bound

end type t_emf_bound

interface
    subroutine read_input_emf_bound( this, input_file, periodic, if_move )
        import t_emf_bound, t_input_file
        class (t_emf_bound), intent( inout ) :: this
        class( t_input_file ), intent(inout) :: input_file
        logical, dimension(:), intent(in) :: periodic, if_move
    end subroutine
end interface

interface
    subroutine init_emf_bound( this, b, e, dt, part_interpolation, &
                               space, no_co, grid, restart, restart_handle )
        import t_emf_bound, t_vdf, t_space, t_node_conf, t_grid, t_restart_handle, p_double
        class (t_emf_bound), intent( inout )  ::  this
        type( t_vdf ),       intent(in) :: b, e
        real( p_double ), intent(in) :: dt
        integer, intent(in) :: part_interpolation
        type( t_space ),     intent(in) :: space
        class( t_node_conf ), intent(in) :: no_co
        class( t_grid ), intent(in) :: grid
        logical, intent(in) :: restart
        type( t_restart_handle ), intent(in) :: restart_handle
    end subroutine
end interface

interface
    subroutine cleanup_emf_bound( this )
        import t_emf_bound
        class (t_emf_bound), intent( inout ) :: this
    end subroutine
end interface

interface
    subroutine write_checkpoint_emf_bound( this, restart_handle )
        import t_emf_bound, t_restart_handle
        class (t_emf_bound), intent(in) :: this
        type( t_restart_handle ), intent(inout) :: restart_handle
    end subroutine
end interface

interface
    subroutine restart_read_emf_bound( this, restart_handle )
        import t_emf_bound, t_restart_handle
        class (t_emf_bound), intent(inout)  ::  this
        type( t_restart_handle ), intent(in) :: restart_handle
    end subroutine
end interface

interface
    subroutine update_boundary_emf_bound( this, b, e, g_space, no_co, send_msg, recv_msg )
        import t_emf_bound, t_vdf, t_space, t_node_conf, t_vdf_msg
        class( t_emf_bound ), intent(inout) :: this
        type( t_vdf ),       intent(inout) :: b, e
        type( t_space ),     intent(in) :: g_space
        class( t_node_conf ), intent(in) :: no_co
        type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
    end subroutine
end interface

interface
    subroutine update_boundary_bfld_emf_bound( this, e, b, step )
        import t_emf_bound, t_vdf
        class( t_emf_bound ), intent(inout) :: this
        class( t_vdf ), intent(inout) :: e, b
        integer, intent(in) :: step
    end subroutine
end interface

interface
    subroutine update_boundary_efld_emf_bound( this, e, b )
        import t_emf_bound, t_vdf
        class( t_emf_bound ), intent(inout) :: this
        class( t_vdf ), intent(inout) :: e, b
    end subroutine
end interface

interface
    subroutine move_window_emf_bound( this, space )
        import t_emf_bound, t_space
        class( t_emf_bound ), intent(inout) :: this
        type( t_space ),  intent(in) :: space
    end subroutine
end interface

interface
    subroutine list_algorithm_emf_bound( this )
        import t_emf_bound
        class( t_emf_bound ), intent(in) :: this
    end subroutine list_algorithm_emf_bound
end interface

public :: p_bcemf_rst_id, t_lindman, t_vpml, t_emf_bound
