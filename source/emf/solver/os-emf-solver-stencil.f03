!-----------------------------------------------------------------------------------------
! Stencil based field solver
!   A. Greenwood, et. al. "On the elimination of numerical Cerenkov radiation in PIC
!   simulations", Journal of Computational Physics, vol. 201, no. 2, pp. 665-684,
!   Dec. 2004
!-----------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_solver_stencil

use m_parameters
use m_restart
use m_emf_define
use m_vdf_define
use m_input_file
    
private  

type, extends( t_emf_solver ), public :: t_emf_solver_stencil

    integer :: gir_pos
    procedure(dedt_stencil), pointer :: dedt => null()
    procedure(dbdt_stencil), pointer :: dbdt => null()

    ! Stencil parameters
    real(p_double) :: k1, k2

contains

    procedure :: init => init_stencil
    procedure :: advance => advance_stencil
    procedure :: read_input => read_input_stencil
    procedure :: min_gc => min_gc_stencil
    procedure :: name => name_stencil
    procedure :: test_stability => test_stability_stencil

end type t_emf_solver_stencil

interface
    subroutine dedt_stencil( this, e, b, jay, dt )
        import t_emf_solver_stencil, t_vdf, p_double
        class( t_emf_solver_stencil ), intent(in) :: this
        type( t_vdf ), intent(inout) :: e
        type( t_vdf ), intent(in) :: b, jay
        real(p_double), intent(in) :: dt
    end subroutine
end interface

interface
    subroutine dbdt_stencil( this, b, e, dt )
        import t_emf_solver_stencil, t_vdf, p_double
        class( t_emf_solver_stencil ), intent(in) :: this
        type( t_vdf ), intent(inout) :: b
        type( t_vdf ), intent(in) :: e
        real(p_double), intent(in) :: dt
    end subroutine
end interface

contains

function name_stencil( this )
  
    implicit none
  
    class( t_emf_solver_stencil ), intent(in) :: this
    character( len = p_maxlen_emf_solver_name ) :: name_stencil
  
    name_stencil = "Stencil"
  
end function name_stencil

subroutine dbdt_2d_stencil( this, b, e, dt )

    implicit none
    
    class( t_emf_solver_stencil ), intent(in) :: this
    type( t_vdf ), intent( inout ) :: b
    type( t_vdf ),    intent(in) :: e
    real(p_double), intent(in) :: dt

    real(p_k_fld) :: dtdx1, dtdx2
    real(p_k_fld) :: A0, A1, A2
    integer :: i1, i2

    dtdx1 = real(dt/ b%dx_(1), p_k_fld)
    dtdx2 = real(dt/ b%dx_(2), p_k_fld)

    A0 = real( 1.0_p_double - this%k1 - this%k2, p_k_fld )
    A1 = real( this%k1 / 3.0_p_double, p_k_fld )
    A2 = real( this%k2 / 6.0_p_double, p_k_fld )

    !$omp parallel do private(i1)
    do i2 = -2, b%nx_(2)+3
        do i1 = -2, b%nx_(1)+3

            !B1
            b%f2( 1, i1, i2 ) = b%f2( 1, i1, i2 ) &
                                - dtdx2 * ( A0 * (e%f2( 3, i1,   i2+1 ) - e%f2( 3, i1,   i2   )) + &
                                            A1 * (e%f2( 3, i1,   i2+2 ) - e%f2( 3, i1,   i2-1 )) + &
                                            A2 * (e%f2( 3, i1+1, i2+2 ) - e%f2( 3, i1+1, i2-1 ) + &
                                                e%f2( 3, i1-1, i2+2 ) - e%f2( 3, i1-1, i2-1 )))

            !B2
            b%f2( 2, i1, i2 ) = b%f2( 2, i1, i2 ) &
                                + dtdx1 * ( A0 * ( e%f2( 3, i1+1, i2 )   - e%f2( 3, i1,   i2 )) + &
                                            A1 * ( e%f2( 3, i1+2, i2 )   - e%f2( 3, i1-1, i2 )) + &
                                            A2 * ( e%f2( 3, i1+2, i2+1 ) - e%f2( 3, i1-1, i2+1 ) + &
                                                    e%f2( 3, i1+2, i2-1 ) - e%f2( 3, i1-1, i2-1 )))

            !B3
            b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) &
                                - dtdx1 * ( A0 * (e%f2( 2, i1+1, i2 )   - e%f2( 2, i1,   i2 )) + &
                                            A1 * (e%f2( 2, i1+2, i2 )   - e%f2( 2, i1-1, i2 )) + &
                                            A2 * (e%f2( 2, i1+2, i2+1 ) - e%f2( 2, i1-1, i2+1 ) + &
                                                e%f2( 2, i1+2, i2-1 ) - e%f2( 2, i1-1, i2-1 ))) &
                                + dtdx2 * ( A0 * (e%f2( 1, i1,   i2+1 ) - e%f2( 1, i1,   i2 )) + &
                                            A1 * (e%f2( 1, i1,   i2+2 ) - e%f2( 1, i1,   i2-1 )) + &
                                            A2 * (e%f2( 1, i1+1, i2+2 ) - e%f2( 1, i1+1, i2-1 ) + &
                                                e%f2( 1, i1-1, i2+2 ) - e%f2( 1, i1-1, i2-1 ) ))

        enddo
    enddo
    !$omp end parallel do

end subroutine dbdt_2d_stencil
    
subroutine dbdt_3d_stencil( this, b, e, dt )

    implicit none

    class( t_emf_solver_stencil ), intent(in) :: this
    type( t_vdf ),   intent( inout ) :: b
    type( t_vdf ),   intent(in) :: e
    real(p_double), intent(in) :: dt

    real(p_k_fld) :: dtdx1, dtdx2, dtdx3
    real(p_k_fld) :: A0, A1, A2
    integer :: i1, i2, i3

    dtdx1 = real( dt/ b%dx_(1), p_k_fld )
    dtdx2 = real( dt/ b%dx_(2), p_k_fld )
    dtdx3 = real( dt/ b%dx_(3), p_k_fld )

    A0 = real( 1.0_p_double - this%k1 - this%k2, p_k_fld )
    A1 = real( this%k1 / 3.0_p_double, p_k_fld )
    A2 = real( this%k2 / 6.0_p_double, p_k_fld )


    !$omp parallel do private(i2,i1)
    do i3 = -2, b%nx_(3)+3
        do i2 = -2, b%nx_(2)+3
            do i1 = -2, b%nx_(1)+3

                b%f3( 1, i1, i2, i3 ) =  b%f3( 1, i1, i2, i3 ) &
                - dtdx2 * ( A0 * ( e%f3( 3, i1  , i2+1, i3   ) - e%f3( 3, i1  , i2  , i3   )) + &
                            A1 * ( e%f3( 3, i1  , i2+2, i3   ) - e%f3( 3, i1  , i2-1, i3   )) + &
                            A2 * ( e%f3( 3, i1  , i2+2, i3+1 ) - e%f3( 3, i1  , i2-1, i3+1 ) + &
                                    e%f3( 3, i1  , i2+2, i3-1 ) - e%f3( 3, i1  , i2-1, i3-1 )))  &
                + dtdx3 * ( A0 * ( e%f3( 2, i1  , i2  , i3+1 ) - e%f3( 2, i1  , i2  , i3   )) + &
                            A1 * ( e%f3( 2, i1  , i2  , i3+2 ) - e%f3( 2, i1  , i2  , i3-1 )) + &
                            A2 * ( e%f3( 2, i1  , i2+1, i3+2 ) - e%f3( 2, i1  , i2+1, i3-1 ) + &
                                    e%f3( 2, i1  , i2-1, i3+2 ) - e%f3( 2, i1  , i2-1, i3-1 )))

                b%f3( 2, i1, i2, i3 ) =  b%f3( 2, i1, i2, i3 ) &
                - dtdx3 * ( A0 * ( e%f3( 1, i1  , i2  , i3+1 ) - e%f3( 1, i1  , i2  , i3   )) + &
                            A1 * ( e%f3( 1, i1  , i2  , i3+2 ) - e%f3( 1, i1  , i2  , i3-1 )) + &
                            A2 * ( e%f3( 1, i1+1, i2  , i3+2 ) - e%f3( 1, i1+1, i2  , i3-1 ) + &
                                    e%f3( 1, i1-1, i2  , i3+2 ) - e%f3( 1, i1-1, i2  , i3-1 ))) &
                + dtdx1 * ( A0 * ( e%f3( 3, i1+1, i2  , i3   ) - e%f3( 3, i1  , i2  , i3   )) + &
                            A1 * ( e%f3( 3, i1+2, i2  , i3   ) - e%f3( 3, i1-1, i2  , i3   )) + &
                            A2 * ( e%f3( 3, i1+2, i2  , i3+1 ) - e%f3( 3, i1-1, i2  , i3+1 ) + &
                                    e%f3( 3, i1+2, i2  , i3-1 ) - e%f3( 3, i1-1, i2  , i3-1 )))

            b%f3( 3, i1, i2, i3 ) =  b%f3( 3, i1, i2, i3 ) &
                - dtdx1 * ( A0 * ( e%f3( 2, i1+1, i2  , i3   ) - e%f3( 2, i1  , i2  , i3   )) + &
                            A1 * ( e%f3( 2, i1+2, i2  , i3   ) - e%f3( 2, i1-1, i2  , i3   )) + &
                            A2 * ( e%f3( 2, i1+2, i2+1, i3   ) - e%f3( 2, i1-1, i2+1, i3   ) + &
                                    e%f3( 2, i1+2, i2-1, i3   ) - e%f3( 2, i1-1, i2-1, i3   ))) &
                + dtdx2 * ( A0 * ( e%f3( 1, i1  , i2+1, i3   ) - e%f3( 1, i1  , i2  , i3   )) + &
                            A1 * ( e%f3( 1, i1  , i2+2, i3   ) - e%f3( 1, i1  , i2-1, i3   )) + &
                            A2 * ( e%f3( 1, i1+1, i2+2, i3   ) - e%f3( 1, i1+1, i2-1, i3   ) + &
                                    e%f3( 1, i1-1, i2+2, i3   ) - e%f3( 1, i1-1, i2-1, i3   )))

            enddo
        enddo
    enddo
    !$omp end parallel do

end subroutine dbdt_3d_stencil
    
subroutine dedt_2d_stencil( this, e, b, jay, dt )

    implicit none
    class( t_emf_solver_stencil ), intent(in) :: this
    type( t_vdf ), intent(inout) :: e
    type( t_vdf ), intent(in) :: b, jay
    real(p_double), intent(in) :: dt

!       local variables
    integer :: i1, i2
    real(p_k_fld) :: dtdx1, dtdx2, dtif
    real(p_k_fld) :: A0, A1, A2

!       executable statements

    dtdx1 = real( dt/e%dx_(1), p_k_fld )
    dtdx2 = real( dt/e%dx_(2), p_k_fld )
    dtif = real( dt, p_k_fld )

    A0 = real( 1.0_p_double - this%k1 - this%k2, p_k_fld )
    A1 = real( this%k1 / 3.0_p_double, p_k_fld )
    A2 = real( this%k2 / 6.0_p_double, p_k_fld )

    !$omp parallel do private(i1)
    do i2 = 0, e%nx_(2)+2
    do i1 = 0, e%nx_(1)+2
        ! E1
        e%f2(1, i1, i2) = e%f2(1, i1, i2) &
                        - dtif * jay%f2(1, i1, i2) &
                        + dtdx2 * ( A0 * (b%f2(3, i1,   i2)   - b%f2(3, i1,   i2-1)) + &
                                    A1 * (b%f2(3, i1,   i2+1) - b%f2(3, i1,   i2-2)) + &
                                    A2 * (b%f2(3, i1+1, i2+1) - b%f2(3, i1+1, i2-2) + &
                                            b%f2(3, i1-1, i2+1) - b%f2(3, i1-1, i2-2)))

        ! E2
        e%f2(2, i1, i2) = e%f2(2, i1, i2) &
                        - dtif * jay%f2(2, i1, i2) &
                        - dtdx1 * ( A0 * (b%f2(3, i1,   i2)   - b%f2(3, i1-1, i2)) + &
                                    A1 * (b%f2(3, i1+1, i2)   - b%f2(3, i1-2, i2)) + &
                                    A2 * (b%f2(3, i1+1, i2+1) - b%f2(3, i1-2, i2+1) + &
                                            b%f2(3, i1+1, i2-1) - b%f2(3, i1-2, i2-1)) )

        ! E3
        e%f2(3, i1, i2) = e%f2(3, i1, i2) &
                        - dtif * jay%f2(3, i1, i2) &
                        + dtdx1 * ( A0 * (b%f2(2, i1,   i2)   - b%f2(2, i1-1, i2)) + &
                                    A1 * (b%f2(2, i1+1, i2)   - b%f2(2, i1-2, i2)) + &
                                    A2 * (b%f2(2, i1+1, i2+1) - b%f2(2, i1-2, i2+1) + &
                                            b%f2(2, i1+1, i2-1) - b%f2(2, i1-2, i2-1))) &
                        - dtdx2 * ( A0 * (b%f2(1, i1,   i2)   - b%f2(1, i1,   i2-1)) + &
                                    A1 * (b%f2(1, i1,   i2+1) - b%f2(1, i1,   i2-2)) + &
                                    A2 * (b%f2(1, i1+1, i2+1) - b%f2(1, i1+1, i2-2) + &
                                            b%f2(1, i1-1, i2+1) - b%f2(1, i1-1, i2-2)) )
    enddo
    enddo
    !$omp end parallel do

end subroutine dedt_2d_stencil
    
subroutine dedt_3d_stencil( this, e, b, jay, dt )

    implicit none

    class( t_emf_solver_stencil ), intent(in) :: this
    type( t_vdf ), intent(inout) :: e
    type( t_vdf ), intent(in) :: b, jay
    real(p_double), intent(in) :: dt

    integer :: i1, i2, i3
    real(p_k_fld) :: dtdx1, dtdx2, dtdx3, dtif
    real(p_k_fld) :: A0, A1, A2

    dtdx1 = real( dt/e%dx_(1), p_k_fld )
    dtdx2 = real( dt/e%dx_(2), p_k_fld )
    dtdx3 = real( dt/e%dx_(3), p_k_fld )
    dtif = real( dt, p_k_fld )

    A0 = real( 1.0_p_double - this%k1 - this%k2 , p_k_fld )
    A1 = real( this%k1 / 3.0_p_double, p_k_fld )
    A2 = real( this%k2 / 6.0_p_double, p_k_fld )

    !$omp parallel do private(i2,i1)
    do i3 = 0, e%nx_(3)+2
        do i2 = 0, e%nx_(2)+2
            do i1 = 0, e%nx_(1)+2

            e%f3( 1, i1, i2, i3 ) = e%f3( 1, i1, i2, i3 ) &
            - dtif  * jay%f3( 1, i1, i2, i3 )   &
            + dtdx2 * ( A0 * ( b%f3( 3, i1,   i2,   i3   ) - b%f3( 3, i1,   i2-1, i3   )) + &
                        A1 * ( b%f3( 3, i1,   i2+1, i3   ) - b%f3( 3, i1,   i2-2, i3   )) + &
                        A2 * ( b%f3( 3, i1,   i2+1, i3+1 ) - b%f3( 3, i1,   i2-2, i3+1 ) + &
                                b%f3( 3, i1,   i2+1, i3-1 ) - b%f3( 3, i1,   i2-2, i3-1 ))) &
            - dtdx3 * ( A0 * ( b%f3( 2, i1,   i2,   i3   ) - b%f3( 2, i1,   i2,   i3-1 )) + &
                        A1 * ( b%f3( 2, i1,   i2,   i3+1 ) - b%f3( 2, i1,   i2,   i3-2 )) + &
                        A2 * ( b%f3( 2, i1,   i2+1, i3+1 ) - b%f3( 2, i1,   i2+1, i3-2 ) + &
                                b%f3( 2, i1,   i2-1, i3+1 ) - b%f3( 2, i1,   i2-1, i3-2 )))

            e%f3( 2, i1, i2, i3 ) = e%f3( 2, i1, i2, i3 ) &
            - dtif  * jay%f3( 2, i1, i2, i3 ) &
            + dtdx3 * ( A0 * ( b%f3( 1, i1,   i2,   i3   ) - b%f3( 1, i1,   i2,   i3-1 )) + &
                        A1 * ( b%f3( 1, i1,   i2,   i3+1 ) - b%f3( 1, i1,   i2,   i3-2 )) + &
                        A2 * ( b%f3( 1, i1+1, i2,   i3+1 ) - b%f3( 1, i1+1, i2,   i3-2 ) + &
                                b%f3( 1, i1-1, i2,   i3+1 ) - b%f3( 1, i1-1, i2,   i3-2 ))) &
            - dtdx1 * ( A0 * ( b%f3( 3, i1,   i2,   i3   ) - b%f3( 3, i1-1, i2,   i3   )) + &
                        A1 * ( b%f3( 3, i1+1, i2,   i3   ) - b%f3( 3, i1-2, i2,   i3   )) + &
                        A2 * ( b%f3( 3, i1+1, i2,   i3+1 ) - b%f3( 3, i1-2, i2,   i3+1 ) + &
                                b%f3( 3, i1+1, i2,   i3-1 ) - b%f3( 3, i1-2, i2,   i3-1 )))

            e%f3( 3, i1, i2, i3 ) = e%f3( 3, i1, i2, i3 ) &
            - dtif  * jay%f3( 3, i1, i2, i3 )   &
            + dtdx1 * ( A0 * ( b%f3( 2, i1,   i2,   i3   ) - b%f3( 2, i1-1, i2,   i3   )) + &
                        A1 * ( b%f3( 2, i1+1, i2,   i3   ) - b%f3( 2, i1-2, i2,   i3   )) + &
                        A2 * ( b%f3( 2, i1+1, i2+1, i3   ) - b%f3( 2, i1-2, i2+1, i3   ) + &
                                b%f3( 2, i1+1, i2-1, i3   ) - b%f3( 2, i1-2, i2-1, i3   ))) &
            - dtdx2 * ( A0 * ( b%f3( 1, i1,   i2,   i3   ) - b%f3( 1, i1,   i2-1, i3   )) + &
                        A1 * ( b%f3( 1, i1,   i2+1, i3   ) - b%f3( 1, i1,   i2-2, i3   )) + &
                        A2 * ( b%f3( 1, i1+1, i2+1, i3   ) - b%f3( 1, i1+1, i2-2, i3   ) + &
                                b%f3( 1, i1-1, i2+1, i3   ) - b%f3( 1, i1-1, i2-2, i3   )))

            enddo
        enddo
    enddo
    !$omp end parallel do

end subroutine dedt_3d_stencil

subroutine test_stability_stencil( this, dt, dx )
    
    implicit none

    class( t_emf_solver_stencil ), intent(in) :: this
    real(p_double), intent(in)  ::  dt
    real(p_double), dimension(:), intent(in)  ::  dx

    real(p_double) :: cour
    integer :: i

    cour = 0.0
    do i = 1, p_x_dim
        cour = cour + 1.0/(dx(i))**2
    enddo
    cour = sqrt(1.0/cour)

    if ( this%k2 == 0 ) then
        cour = cour * 1.0_p_double/(1.0_p_double - 4.0_p_double * this%k1 / 3.0_p_double)
    else
        cour = cour * 1.0_p_double/(1.0_p_double - 8.0_p_double * this%k1 / 3.0_p_double)
    endif


    if (dt > cour) then
        if ( mpi_node() == 0 ) then
           print *, '(*error*) Stencil EMF solver stability condition violated, aborting'
           print *, 'dx = ', dx(1:p_x_dim)
           print *, 'dt = ', dt, ' max(dt) = ', cour
           print *, 'k1 = ', this%k1, ' k2 = ', this%k2 
        endif
        call abort_program(p_err_invalid)
    endif

end subroutine test_stability_stencil

subroutine min_gc_stencil( this, min_gc )

    implicit none

    class( t_emf_solver_stencil ), intent(in) :: this
    integer, dimension(:,:), intent(inout) :: min_gc

    min_gc(p_lower,:) = 4
    min_gc(p_upper,:) = 5

end subroutine

subroutine init_stencil( this, gix_pos, coords )

    implicit none

    class( t_emf_solver_stencil ), intent(inout) :: this
    integer, intent(in), dimension(:) :: gix_pos
    integer, intent(in) :: coords

    select case (p_x_dim)
    case(1)
        if ( mpi_node() == 0 ) then
            write(0,*) "(*error*) Stencil EMF solver not implemented in 1D"
            write(0,*) "(*error*) aborting..."
        endif
        call abort_program( p_err_notimplemented )
    case(2)
        if ( coords == p_cylindrical_b ) then
            if ( mpi_node() == 0 ) then
                write(0,*) "(*error*) Stencil EMF solver not implemented in 2D cylindrical geometry"
                write(0,*) "(*error*) aborting..."
            endif
            call abort_program( p_err_notimplemented )
        else
            this % dedt => dedt_2d_stencil
            this % dbdt => dbdt_2d_stencil
        endif
    case(3)
        this % dedt => dedt_3d_stencil
        this % dbdt => dbdt_3d_stencil
    end select

end subroutine init_stencil

subroutine read_input_stencil( this, input_file, dx, dt )

    implicit none

    class( t_emf_solver_stencil ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file
    real(p_double), dimension(:), intent(in) :: dx
    real(p_double), intent(in) :: dt

    real(p_double) :: k1, k2
    integer :: ierr

    namelist /nl_emf_solver/ k1, k2

    k1 = 0
    k2 = 0

    call get_namelist( input_file, "nl_emf_solver", ierr )

    if ( ierr /= 0 ) then
        if ( mpi_node() == 0 ) then
            if (ierr < 0) then
                write(0,*) "Error reading EMF Stencil solver parameters"
            else
                write(0,*) "Error: EMF Stencil solver parameters missing"
            endif
            write(0,*) "aborting..."
        endif
        stop
    endif

    read (input_file%nml_text, nml = nl_emf_solver, iostat = ierr)
    if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
            write(0,*) "Error reading EMF Stencil solver parameters"
            write(0,*) "aborting..."
        endif
        stop
    endif

    if ( k1 > 0.0_p_k_fld ) then
        if ( mpi_node() == 0 ) then
            write(0,*) "Error reading EMF Stencil solver parameters"
            write(0,*) "k1 must be <= 0.0"
            write(0,*) "aborting..."
        endif
        stop
    endif

    if (( k2 /= 0.0_p_k_fld ) .and. (( k2 /= 2 * k1 ))) then
        if ( mpi_node() == 0 ) then
            write(0,*) "Error reading EMF Stencil solver parameters"
            write(0,*) "k2 must be 0.0 or 2 * k1"
            write(0,*) "aborting..."
        endif
        stop
    endif

    this%k1 = k1
    this%k2 = k2

end subroutine read_input_stencil

subroutine advance_stencil( this, e, b, jay, dt, bnd_con )
  
    implicit none
  
    class( t_emf_solver_stencil ), intent(inout) :: this
    class( t_vdf ), intent( inout )  ::  e, b, jay
    class( t_emf_bound ), intent(inout) :: bnd_con
    real(p_double), intent(in) :: dt
  
    real(p_double) :: dt_b, dt_e
  
    dt_b = dt / 2.0_p_double
    dt_e = dt
  
    ! Advance B half time step
    call this % dbdt( b, e, dt_b )
    call bnd_con%update_boundary_b( e, b, step = 1 )
  
    ! Advance E one full time step
    call this % dedt( e, b, jay, dt_e )
    call bnd_con%update_boundary_e( e, b )
  
    ! Advance B another half time step
    call this % dbdt( b, e, dt_b )
    call bnd_con%update_boundary_b( e, b, step = 2 )
    
end subroutine advance_stencil

end module