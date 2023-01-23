module m_emf_solver

private

integer, parameter, public :: p_emf_yee = 0
! integer, parameter, public :: p_emf_4order = 1
! integer, parameter, public :: p_emf_ck = 2
integer, parameter, public :: p_emf_lehe = 3
integer, parameter, public :: p_emf_ndfx = 4
integer, parameter, public :: p_emf_stencil = 5
integer, parameter, public :: p_emf_fei = 6

interface new_solver_emf
    module procedure new_solver_emf
end interface

public :: new_solver_emf

contains

subroutine new_solver_emf( this, name )

    use m_emf_define, only : t_emf
    use m_emf_solver_yee
    use m_emf_solver_lehe
    use m_emf_solver_ndfx
    use m_emf_solver_stencil
    use m_emf_solver_fei
    use m_system

    implicit none

    class( t_emf ), intent(inout) :: this
    character( len = * ), intent(in) :: name

    select case( trim(name) )
    case( "yee" )
        allocate( t_emf_solver_yee :: this%solver )
    case( "lehe" )
        allocate( t_emf_solver_lehe :: this%solver )
    case( "ndfx" )
        allocate( t_emf_solver_ndfx :: this%solver )
    case( "stencil" )
        allocate( t_emf_solver_stencil :: this%solver )
    case( "custom", "fei" )
        allocate( t_emf_solver_fei :: this%solver )
    case default
        if ( mpi_node() == 0 ) then
            write(0,*) "Invalid EMF solver type: '", trim(name), "'"
        endif
        this%solver => null()
    end select

end subroutine new_solver_emf


end module m_emf_solver


function solver_type_class_emf( this )

    use m_emf_solver
    use m_emf_define
    use m_emf_solver_yee
    use m_emf_solver_lehe
    use m_emf_solver_ndfx
    use m_emf_solver_stencil
    use m_emf_solver_fei

    implicit none

    integer :: solver_type_class_emf
    class( t_emf ), intent(in) :: this

    select type( solver => this%solver )   
    class is ( t_emf_solver_yee )  
        solver_type_class_emf = p_emf_yee
    class is ( t_emf_solver_lehe )
        solver_type_class_emf = p_emf_lehe
    class is ( t_emf_solver_ndfx )  
        solver_type_class_emf = p_emf_ndfx
    class is ( t_emf_solver_stencil )  
        solver_type_class_emf = p_emf_stencil
    class is ( t_emf_solver_fei )  
        solver_type_class_emf = p_emf_fei
    class default
        ! Unknown class
        solver_type_class_emf = -1
    end select 

end function solver_type_class_emf

function solver_type_name_emf( this, name )

    use m_emf_solver
    use m_emf_define

    implicit none

    integer :: solver_type_name_emf
    class( t_emf ), intent(in) :: this
    character( len = * ), intent(in) :: name

    select case( trim(name) )
    case( "yee" )
        solver_type_name_emf = p_emf_yee
    case( "lehe" )
        solver_type_name_emf = p_emf_lehe
    case( "ndfx" )
        solver_type_name_emf = p_emf_ndfx
    case( "stencil" )
        solver_type_name_emf = p_emf_stencil
    case( "custom", "fei" )
        solver_type_name_emf = p_emf_fei
    case default    
        ! Unknown class
        solver_type_name_emf = -1
    end select 

end function solver_type_name_emf

! Most solvers don't have anything to cleanup
subroutine cleanup_emf_solver( this )
    use m_emf_define
    implicit none
    class( t_emf_solver ), intent(inout) :: this
end subroutine cleanup_emf_solver
