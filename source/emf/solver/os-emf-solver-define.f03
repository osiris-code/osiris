! This file is meant to be included by os-emf-define.f03

integer, parameter, public :: p_maxlen_emf_solver_name = 32

type, abstract, public :: t_emf_solver

contains
  procedure(init_emf_solver), deferred :: init
  procedure(advance_emf_solver), deferred  :: advance
  procedure(min_gc_emf_solver), deferred  :: min_gc
  procedure(name_emf_solver), deferred  :: name
  procedure(read_input_emf_solver), deferred  :: read_input
  procedure(test_stability_emf_solver), deferred :: test_stability
  procedure :: cleanup => cleanup_emf_solver
end type

abstract interface
  subroutine init_emf_solver( this, gix_pos, coords )
    import t_emf_solver
    class( t_emf_solver ), intent(inout) :: this
    integer, intent(in), dimension(:) :: gix_pos
    integer, intent(in) :: coords
  end subroutine
end interface

abstract interface
  subroutine advance_emf_solver( this, e, b, jay, dt, bnd_con )
    import t_emf_solver, t_vdf, p_double, t_emf_bound
    class( t_emf_solver ), intent(inout) :: this
    class( t_vdf ), intent(inout) :: e, b, jay
    real(p_double), intent(in) :: dt
    class( t_emf_bound ), intent(inout) :: bnd_con
  end subroutine
end interface

abstract interface
  subroutine min_gc_emf_solver( this, min_gc )
    import t_emf_solver, p_x_dim
    class( t_emf_solver ), intent(in) :: this
    integer, dimension(:,:), intent(inout) :: min_gc
  end subroutine
end interface

abstract interface
  function name_emf_solver( this )
    import t_emf_solver, p_maxlen_emf_solver_name
    class( t_emf_solver ), intent(in) :: this
    character( len = p_maxlen_emf_solver_name ) :: name_emf_solver
  end function name_emf_solver
end interface

abstract interface
  subroutine read_input_emf_solver( this, input_file, dx, dt )
    import t_emf_solver, t_input_file, p_double
    class( t_emf_solver ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file
    real(p_double), dimension(:), intent(in) :: dx
    real(p_double), intent(in) :: dt
  end subroutine read_input_emf_solver
end interface

abstract interface
  subroutine test_stability_emf_solver( this, dt, dx )
    import t_emf_solver, p_double
    class( t_emf_solver ), intent(in) :: this
    real(p_double), intent(in)  ::  dt
    real(p_double), dimension(:), intent(in)  ::  dx
  end subroutine test_stability_emf_solver
end interface  

! Most solvers don't use the cleanup feature, so a simple subroutine is defined
! in os-emf-solver that doesn't do anything
interface
  subroutine cleanup_emf_solver( this )
    import t_emf_solver
    class( t_emf_solver ), intent(inout) :: this
  end subroutine cleanup_emf_solver
end interface
