!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Grid class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!#define DEBUG_FILE 1
#include "os-preprocess.fpp"

module m_grid

  #include "memory/memory.h"

  use m_system
  use m_parameters

  use m_grid_define

  use m_node_conf

  use m_restart
  use m_math

  use m_utilities
  use m_logprof

  use m_fparser
  implicit none

  private

  interface my_nx_p
    module procedure my_nx_p_grid
  end interface

  interface g_nx
    module procedure g_nx_grid
    module procedure g_nx_grid_dir
  end interface

  interface my_nx_p_min
    module procedure my_nx_p_min_grid
    module procedure my_nx_p_min_grid_dir
  end interface

  interface my_nx
    module procedure my_nx_grid
  end interface

  interface nx_p
    module procedure nx_p_grid
  end interface

  interface equal_my_nx
    module procedure equal_my_nx_grid
  end interface

  interface operator(==)
    module procedure equal_grid
  end interface

  interface operator(/=)
    module procedure different_grid
  end interface

  interface coordinates
    module procedure coordinates_grid
  end interface

  interface local_vol
    module procedure local_grid_vol
  end interface

! declare things that should be public
  public :: coordinates, nx_p, g_nx
  public :: restart_write
  public :: my_nx, my_nx_p, my_nx_p_min
  public :: local_vol
  public :: equal_my_nx

  public :: p_density, p_particles, p_global, p_node, p_grid

  public :: operator(==), operator(/=)


contains



!-----------------------------------------------------------------------------------------
!  Gets the node limits for the supplied node with respect to global simulation.
!  The function returns an integer array in the form (bound, x_dim) where bound is
!  1:3 (lower/upper/#cells) and x_dim the  direction we are interested in
!-----------------------------------------------------------------------------------------
function nx_p_grid( this, no_co, node )

  implicit none


  class( t_grid ), intent(in)     :: this
  class( t_node_conf ), intent(in)        :: no_co
  integer, intent(in) :: node

  integer, dimension(3,this%x_dim) :: nx_p_grid

  integer :: i, idx

  idx = 0
  do i = 1, this%x_dim
    nx_p_grid( :, i ) = this%nx_node( :, idx + no_co%ngp( node, i))
    idx = idx + this%nnodes(i)
  enddo

end function nx_p_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!  gets the node limits for the local node with respect to
!  global simulation The function returns an integer array in the form
!  (bound, x_dim) where bound is 1:3 (lower/upper/#cells) and x_dim the
!  direction we are interested in
!-----------------------------------------------------------------------------------------
function my_nx_p_grid( this )

  implicit none


  class( t_grid ), intent(in)     :: this

  integer, dimension(3,this%x_dim) :: my_nx_p_grid

  my_nx_p_grid = this%my_nx(:,1:this%x_dim)

end function my_nx_p_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!  gets the global number of cells for a grid
!-----------------------------------------------------------------------------------------
function g_nx_grid( this )

  implicit none

  class( t_grid ), intent(in)      :: this

  integer, dimension(1:this%x_dim) :: g_nx_grid

  g_nx_grid = this%g_nx(1:this%x_dim)

end function g_nx_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!  gets the global number of cells for a grid for a given direction
!-----------------------------------------------------------------------------------------
function g_nx_grid_dir( this, dir )

  implicit none

  class( t_grid ), intent(in) :: this
  integer,         intent(in) :: dir

  integer :: g_nx_grid_dir

  g_nx_grid_dir = this%g_nx(dir)

end function g_nx_grid_dir

!-----------------------------------------------------------------------------------------
!  gets the lower node limits for the local node with respect to
!  global simulation The function returns an integer array in the form
!  (x_dim) where x_dim is the direction we are interested in
!-----------------------------------------------------------------------------------------
function my_nx_p_min_grid( this )

  implicit none

  class( t_grid ), intent(in)     :: this

  integer, dimension(this%x_dim)    :: my_nx_p_min_grid

  my_nx_p_min_grid = this%my_nx(1,1:this%x_dim)

end function my_nx_p_min_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!  gets the lower node limits for the local node with respect to
!  global simulation for the requested direction
!-----------------------------------------------------------------------------------------
function my_nx_p_min_grid_dir( this, dir )

  implicit none

  class( t_grid ), intent(in)     :: this
  integer, intent(in) :: dir

  integer    :: my_nx_p_min_grid_dir


  my_nx_p_min_grid_dir = this%my_nx(1,dir)


end function my_nx_p_min_grid_dir
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!  gets the node size for the local node.
!-----------------------------------------------------------------------------------------
function my_nx_grid( this )

  implicit none

  class( t_grid ), intent(in)     :: this

  integer, dimension(this%x_dim) :: my_nx_grid

  my_nx_grid = this%my_nx(3,1:this%x_dim)

end function my_nx_grid
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!  Gets the maximum number of cells along specified direction for every partition
!-----------------------------------------------------------------------------------------
function get_max_nx_grid( this, dir )

  implicit none

  class( t_grid ), intent(in)     :: this
  integer, intent(in) :: dir
  integer :: get_max_nx_grid

  get_max_nx_grid = this%max_nx( dir )

end function get_max_nx_grid
!-----------------------------------------------------------------------------------------





!-----------------------------------------------------------------------------------------
! check if grids are the same for local node in the two grid objects
!-----------------------------------------------------------------------------------------
function equal_my_nx_grid( lb_a, lb_b )

  implicit none

  class( t_grid ), intent(in)     :: lb_a, lb_b

  logical :: equal_my_nx_grid

  integer :: i

  equal_my_nx_grid = .true.

  do i = 1, lb_a%x_dim

    if (( lb_a%my_nx( 1, i ) /= lb_b%my_nx( 1, i ) ) .or. &
        ( lb_a%my_nx( 2, i ) /= lb_b%my_nx( 2, i ) )) then
      equal_my_nx_grid = .false.
      exit
    endif

  enddo

end function equal_my_nx_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! check if the 2 grid objects are equal
!-----------------------------------------------------------------------------------------
function equal_grid(lb_a, lb_b)

  implicit none

  logical :: equal_grid
  class( t_grid ), intent(in) :: lb_a, lb_b

  integer :: i, tnodes

  equal_grid = .true.

  ! check if dimensions match
  if ( lb_a%x_dim /= lb_b%x_dim ) then
  equal_grid = .false.
  return
  endif

  ! check if parallel partition matches
  tnodes = 1
  do i = 1, lb_a%x_dim
    if ( lb_a%nnodes(i) /= lb_b%nnodes(i) ) then
       equal_grid = .false.
       return
    endif
    tnodes = tnodes * lb_a%nnodes(i)
  enddo

  do i = 1, tnodes
    if (( lb_a%nx_node(1,i) /= lb_b%nx_node(1,i)) .or. &
        ( lb_a%nx_node(2,i) /= lb_b%nx_node(2,i)) .or. &
        ( lb_a%nx_node(3,i) /= lb_b%nx_node(3,i))) then
       equal_grid = .false.
       exit
    endif
  enddo

end function equal_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! check if the 2 grid objects are different
!-----------------------------------------------------------------------------------------
function different_grid(lb_a, lb_b)

  implicit none

  logical :: different_grid
  class( t_grid ), intent(in) :: lb_a, lb_b


  different_grid = .not. (lb_a == lb_b)

end function different_grid
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Returns coordinate system used
!-----------------------------------------------------------------------------------------
function coordinates_grid( this )

  implicit none

  integer :: coordinates_grid

  class( t_grid ), intent( in )  ::  this

  coordinates_grid = this%coordinates

end function coordinates_grid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Return local grid volume
!-----------------------------------------------------------------------------------------
function local_grid_vol( this )

  implicit none

  integer :: local_grid_vol
  class( t_grid ), intent( in )  ::  this

  integer :: i

  local_grid_vol = 1
  do i = 1, p_x_dim
    local_grid_vol = local_grid_vol * this%my_nx(3,i)
  enddo

end function local_grid_vol
!-----------------------------------------------------------------------------------------

end module m_grid
