!-----------------------------------------------------------------------------------------
! NDFX Solver
! A. Pukhov, "Three-dimensional electromagnetic relativistic particle-in-cell code VLPL
!    (Virtual Laser Plasma Lab)", J Plasma Phys, vol. 61, pp. 425-433, 1999.
!-----------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_solver_ndfx

use m_parameters
use m_emf_define
use m_vdf_define
use m_input_file

private  

type, extends( t_emf_solver ), public :: t_emf_solver_ndfx

    integer :: gir_pos
    procedure(dedt_ndfx), pointer :: dedt => null()
    procedure(dbdt_ndfx), pointer :: dbdt => null()

contains

    procedure :: init => init_ndfx
    procedure :: advance => advance_ndfx
    procedure :: min_gc => min_gc_ndfx
    procedure :: name => name_ndfx
    procedure :: read_input => read_input_ndfx
    procedure :: test_stability => test_stability_ndfx

end type t_emf_solver_ndfx

interface
    subroutine dedt_ndfx( this, e, b, jay, dt )
        import t_emf_solver_ndfx, t_vdf, p_double
        class( t_emf_solver_ndfx ), intent(in) :: this
        type( t_vdf ), intent(inout) :: e
        type( t_vdf ), intent(in) :: b, jay
        real(p_double), intent(in) :: dt
    end subroutine
end interface

interface
    subroutine dbdt_ndfx( this, b, e, dt )
        import t_emf_solver_ndfx, t_vdf, p_double
        class( t_emf_solver_ndfx ), intent(in) :: this
        type( t_vdf ), intent(inout) :: b
        type( t_vdf ), intent(in) :: e
        real(p_double), intent(in) :: dt
    end subroutine
end interface


contains

function name_ndfx( this )
  
    implicit none
  
    class( t_emf_solver_ndfx ), intent(in) :: this
    character( len = p_maxlen_emf_solver_name ) :: name_ndfx
  
    name_ndfx = "NDFX (Pukhov)"
  
end function name_ndfx

subroutine dbdt_2d_ndfx( this, b, e, dt )
    
    implicit none

    class( t_emf_solver_ndfx ), intent(in) :: this
    type( t_vdf ), intent( inout ) :: b
    type( t_vdf ),    intent(in) :: e
    real(p_double),               intent(in) :: dt

    real(p_k_fld) :: dtdx1, dtdx2, Ax, Ay, Bx, By
    integer :: i1, i2

    dtdx1 = real(dt/b%dx_(1), p_k_fld)
    dtdx2 = real(dt/b%dx_(2), p_k_fld)

    Bx = 0.75_p_k_fld
    By = real( 1.0_p_k_fld - 0.25_p_k_fld * ( e%dx_(1)/e%dx_(2) )**2, p_k_fld )
    Ax = 0.5_p_k_fld * ( 1.0_p_k_fld - Bx )
    Ay = 0.5_p_k_fld * ( 1.0_p_k_fld - By )

    !$omp parallel do private(i1)
    do i2 = -1, b%nx_(2)+2
        do i1 = -1, b%nx_(1)+2

            !B1
            b%f2( 1, i1, i2 ) = b%f2( 1, i1, i2 ) &
                                - dtdx2 * ( e%f2( 3, i1, i2+1 ) - e%f2( 3, i1, i2 ))

            !B2
            b%f2( 2, i1, i2 ) = b%f2( 2, i1, i2 ) &
                                + dtdx1 * ( e%f2( 3, i1+1, i2 ) - e%f2( 3, i1, i2 ))

            !B3
            b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) &
                                - dtdx1 * ( By * ( e%f2( 2, i1+1, i2   ) - e%f2( 2, i1, i2   ))  &
                                          + Ay * ( e%f2( 2, i1+1, i2+1 ) - e%f2( 2, i1, i2+1 )   &
                                                 + e%f2( 2, i1+1, i2-1 ) - e%f2( 2, i1, i2-1 ))) &
                                + dtdx2 * ( Bx * ( e%f2( 1, i1  , i2+1 ) - e%f2( 1, i1  , i2 ))  &
                                          + Ax * ( e%f2( 1, i1+1, i2+1 ) - e%f2( 1, i1+1, i2 )   &
                                                 + e%f2( 1, i1-1, i2+1 ) - e%f2( 1, i1-1, i2 )))
        enddo
    enddo
    !$omp end parallel do

end subroutine dbdt_2d_ndfx
    
subroutine dbdt_3d_ndfx( this, b, e, dt )

    implicit none

    class( t_emf_solver_ndfx ), intent(in) :: this
    type( t_vdf ),   intent( inout ) :: b
    type( t_vdf ),   intent(in) :: e
    real(p_double), intent(in) :: dt

    real(p_k_fld) :: dtdx1, dtdx2, dtdx3
    real(p_k_fld) :: k, Ax, Ay, Az, Bx, By, Bz
    integer :: i1, i2, i3

    dtdx1 = real( dt/ b%dx_(1), p_k_fld )
    dtdx2 = real( dt/ b%dx_(2), p_k_fld )
    dtdx3 = real( dt/ b%dx_(3), p_k_fld )

    k = 0.264156081 ! 1.25 * [sqrt(3) - 1] / [2*sqrt(3)]
    Bx = 0.735843918 ! 1 - k
    By = real( 1.0_p_k_fld - k * ( e%dx_(1)/e%dx_(2) )**2, p_k_fld )
    Bz = real( 1.0_p_k_fld - k * ( e%dx_(1)/e%dx_(3) )**2, p_k_fld )
    Ax = 0.5_p_k_fld * (1.0_p_k_fld - Bx)
    Ay = 0.5_p_k_fld * (1.0_p_k_fld - By)
    Az = 0.5_p_k_fld * (1.0_p_k_fld - Bz)

    !$omp parallel do private(i2,i1)
    do i3 = -1, b%nx_(3)+2
        do i2 = -1, b%nx_(2)+2
            do i1 = -1, b%nx_(1)+2

                b%f3( 1, i1, i2, i3 ) =  b%f3( 1, i1, i2, i3 ) &
                - dtdx2 * ( Bz * ( e%f3( 3, i1  , i2+1, i3   ) - e%f3( 3, i1  , i2  , i3   ))  &
                          + Az * ( e%f3( 3, i1  , i2+1, i3+1 ) - e%f3( 3, i1  , i2  , i3+1 )   &
                                + e%f3( 3, i1  , i2+1, i3-1 ) - e%f3( 3, i1  , i2  , i3-1 ))) &
                + dtdx3 * ( By * ( e%f3( 2, i1  , i2  , i3+1 ) - e%f3( 2, i1  , i2  , i3   ))  &
                          + Ay * ( e%f3( 2, i1  , i2+1, i3+1 ) - e%f3( 2, i1  , i2+1, i3   )   &
                                + e%f3( 2, i1  , i2-1, i3+1 ) - e%f3( 2, i1  , i2-1, i3   )))

            b%f3( 2, i1, i2, i3 ) =  b%f3( 2, i1, i2, i3 ) &
                - dtdx3 * ( Bx * ( e%f3( 1, i1  , i2  , i3+1 ) - e%f3( 1, i1  , i2  , i3   ))  &
                          + Ax * ( e%f3( 1, i1+1, i2  , i3+1 ) - e%f3( 1, i1+1, i2  , i3   )   &
                                + e%f3( 1, i1-1, i2  , i3+1 ) - e%f3( 1, i1-1, i2  , i3   ))) &
                + dtdx1 * ( Bz * ( e%f3( 3, i1+1, i2  , i3   ) - e%f3( 3, i1  , i2  , i3   ))  &
                          + Az * ( e%f3( 3, i1+1, i2  , i3+1 ) - e%f3( 3, i1  , i2  , i3+1 )   &
                                + e%f3( 3, i1+1, i2  , i3-1 ) - e%f3( 3, i1  , i2  , i3-1 )))

            b%f3( 3, i1, i2, i3 ) =  b%f3( 3, i1, i2, i3 ) &
                - dtdx1 * ( By * ( e%f3( 2, i1+1, i2  , i3   ) - e%f3( 2, i1  , i2  , i3   ))  &
                          + Ay * ( e%f3( 2, i1+1, i2+1, i3   ) - e%f3( 2, i1  , i2+1, i3   )   &
                                + e%f3( 2, i1+1, i2-1, i3   ) - e%f3( 2, i1  , i2-1, i3   ))) &
                + dtdx2 * ( Bx * ( e%f3( 1, i1  , i2+1, i3   ) - e%f3( 1, i1  , i2  , i3   ))  &
                          + Ax * ( e%f3( 1, i1+1, i2+1, i3   ) - e%f3( 1, i1+1, i2  , i3   )   &
                                + e%f3( 1, i1-1, i2+1, i3   ) - e%f3( 1, i1-1, i2  , i3   )))

            enddo
        enddo
    enddo
    !$omp end parallel do

end subroutine dbdt_3d_ndfx
    
subroutine dedt_2d_ndfx( this, e, b, jay, dt )

    implicit none

    class( t_emf_solver_ndfx ), intent(in) :: this
    type( t_vdf ), intent(inout) :: e
    type( t_vdf ), intent(in) :: b, jay
    real(p_double), intent(in) :: dt

    real(p_k_fld) :: dtdx1, dtdx2, dtif, Ax, Ay, Bx, By
    integer :: i1, i2

    dtdx1 = real( dt/e%dx_(1), p_k_fld )
    dtdx2 = real( dt/e%dx_(2), p_k_fld )
    dtif = real( dt, p_k_fld )

    Bx = 0.75_p_k_fld
    By = real( 1.0_p_k_fld - 0.25_p_k_fld * ( e%dx_(1)/e%dx_(2) )**2, p_k_fld)
    Ax = 0.5_p_k_fld * ( 1.0_p_k_fld - Bx )
    Ay = 0.5_p_k_fld * ( 1.0_p_k_fld - By )

    !$omp parallel do private(i1)
    do i2 = 0, e%nx_(2)+1
        do i1 = 0, e%nx_(1)+1

            ! E1
            e%f2(1, i1, i2) = e%f2(1, i1, i2) &
                            - dtif * jay%f2(1, i1, i2) &
                            + dtdx2 * ( b%f2(3, i1, i2) - b%f2(3, i1, i2-1) )

            ! E2
            e%f2(2, i1, i2) = e%f2(2, i1, i2) &
                            - dtif * jay%f2(2, i1, i2) &
                            - dtdx1 * ( b%f2(3, i1, i2) - b%f2(3, i1-1, i2) )

            ! E3
            e%f2(3, i1, i2) = e%f2(3, i1, i2) &
                            - dtif * jay%f2(3, i1, i2) &
                            + dtdx1 * ( Bx * ( b%f2(2, i1, i2  ) - b%f2(2, i1-1, i2  ))  &
                                        + Ax * ( b%f2(2, i1, i2+1) - b%f2(2, i1-1, i2+1)   &
                                                + b%f2(2, i1, i2-1) - b%f2(2, i1-1, i2-1))) &
                            - dtdx2 * ( By * ( b%f2(1, i1  , i2) - b%f2(1, i1  , i2-1))  &
                                        + Ay * ( b%f2(2, i1+1, i2) - b%f2(2, i1+1, i2-1)   &
                                                + b%f2(2, i1-1, i2) - b%f2(2, i1-1, i2-1)))

        enddo
    enddo
    !$omp end parallel do

end subroutine dedt_2d_ndfx
    
subroutine dedt_3d_ndfx( this, e, b, jay, dt )

    implicit none

    class( t_emf_solver_ndfx ), intent(in) :: this
    type( t_vdf ), intent(inout) :: e
    type( t_vdf ), intent(in) :: b, jay
    real(p_double), intent(in) :: dt

    integer :: i1, i2, i3
    real(p_k_fld) :: dtdx1, dtdx2, dtdx3, dtif
    real(p_k_fld) :: k, Ax, Ay, Az, Bx, By, Bz

    dtdx1 = real( dt/e%dx_(1), p_k_fld )
    dtdx2 = real( dt/e%dx_(2), p_k_fld )
    dtdx3 = real( dt/e%dx_(3), p_k_fld )
    dtif = real( dt, p_k_fld )

    k = 0.264156081 ! 1.25 * [sqrt(3) - 1] / [2*sqrt(3)]
    Bx = 0.735843918 ! 1 - k
    By = real( 1.0_p_k_fld - k * ( e%dx_(1)/e%dx_(2) )**2, p_k_fld )
    Bz = real( 1.0_p_k_fld - k * ( e%dx_(1)/e%dx_(3) )**2, p_k_fld )
    Ax = 0.5_p_k_fld * (1.0_p_k_fld - Bx)
    Ay = 0.5_p_k_fld * (1.0_p_k_fld - By)
    Az = 0.5_p_k_fld * (1.0_p_k_fld - Bz)

    !$omp parallel do private(i2,i1)
    do i3 = 0, e%nx_(3)+1
        do i2 = 0, e%nx_(2)+1
            do i1 = 0, e%nx_(1)+1

            e%f3( 1, i1, i2, i3 ) = e%f3( 1, i1, i2, i3 ) &
                - dtif  * jay%f3( 1, i1, i2, i3 )   &
                + dtdx2 * ( Bz * ( b%f3( 3, i1  , i2  , i3   ) - b%f3( 3, i1  , i2-1, i3   ))  &
                        + Az * ( b%f3( 3, i1  , i2  , i3+1 ) - b%f3( 3, i1  , i2-1, i3+1 )   &
                                + b%f3( 3, i1  , i2  , i3-1 ) - b%f3( 3, i1  , i2-1, i3-1 ))) &
                - dtdx3 * ( By * ( b%f3( 2, i1  , i2  , i3   ) - b%f3( 2, i1  , i2  , i3-1 ))  &
                        + Ay * ( b%f3( 2, i1  , i2+1, i3   ) - b%f3( 2, i1  , i2+1, i3-1 )   &
                                + b%f3( 2, i1  , i2-1, i3   ) - b%f3( 2, i1  , i2-1, i3-1 )))

            e%f3( 2, i1, i2, i3 ) = e%f3( 2, i1, i2, i3 ) &
                - dtif  * jay%f3( 2, i1, i2, i3 ) &
                + dtdx3 * ( Bx * ( b%f3( 1, i1  , i2  , i3   ) - b%f3( 1, i1  , i2  , i3-1 ))  &
                        + Ax * ( b%f3( 1, i1+1, i2  , i3   ) - b%f3( 1, i1+1, i2  , i3-1 )   &
                                + b%f3( 1, i1-1, i2  , i3   ) - b%f3( 1, i1-1, i2  , i3-1 ))) &
                - dtdx1 * ( Bz * ( b%f3( 3, i1,   i2  , i3   ) - b%f3( 3, i1-1, i2  , i3   ))  &
                        + Az * ( b%f3( 3, i1,   i2  , i3+1 ) - b%f3( 3, i1-1, i2  , i3+1 )   &
                                + b%f3( 3, i1,   i2  , i3-1 ) - b%f3( 3, i1-1, i2  , i3-1 )))

            e%f3( 3, i1, i2, i3 ) = e%f3( 3, i1, i2, i3 ) &
                - dtif  * jay%f3( 3, i1, i2, i3 )   &
                + dtdx1 * ( By * ( b%f3( 2, i1  , i2  , i3   ) - b%f3( 2, i1-1, i2  , i3   ))  &
                        + Ay * ( b%f3( 2, i1  , i2+1, i3   ) - b%f3( 2, i1-1, i2+1, i3   )   &
                                + b%f3( 2, i1  , i2-1, i3   ) - b%f3( 2, i1-1, i2-1, i3   ))) &
                - dtdx2 * ( Bx * ( b%f3( 1, i1  , i2  , i3   ) - b%f3( 1, i1  , i2-1, i3   ))  &
                        + Ax * ( b%f3( 1, i1+1, i2  , i3   ) - b%f3( 1, i1+1, i2-1, i3   )   &
                                + b%f3( 1, i1-1, i2  , i3   ) - b%f3( 1, i1-1, i2-1, i3   )))
            enddo
        enddo
    enddo
    !$omp end parallel do

end subroutine dedt_3d_ndfx

subroutine test_stability_ndfx( this, dt, dx )
    
    implicit none

    class( t_emf_solver_ndfx ), intent(in) :: this
    real(p_double), intent(in)  ::  dt
    real(p_double), dimension(:), intent(in)  ::  dx

    real(p_double) :: cour
    integer :: i

    cour = dx(1)
    do i = 2, p_x_dim
        if (dx(i) < dx(1)) then
            if ( mpi_node() == 0 ) then
                print *, '(*error*) NDFX solver stability condition violated, aborting'
                print *, 'dx   = ', dx(1:p_x_dim)
                print *, 'Please ensure that dx1 < dx2 (and dx1 < dx3).'
            endif
            call abort_program(p_err_invalid)
        endif
    enddo

    if (dt > cour) then
        if ( mpi_node() == 0 ) then
           print *, '(*error*) NDFX solver stability condition violated, aborting'
           print *, 'dx   = ', dx(1:p_x_dim)
           print *, 'dt   = ', dt, ' max(dt) = ', cour
        endif
        call abort_program(p_err_invalid)
    endif

end subroutine test_stability_ndfx

subroutine min_gc_ndfx( this, min_gc )

    implicit none

    class( t_emf_solver_ndfx ), intent(in) :: this
    integer, dimension(:,:), intent(inout) :: min_gc

    min_gc(p_lower,:) = 3
    min_gc(p_upper,:) = 3

end subroutine


subroutine init_ndfx( this, gix_pos, coords )

    implicit none

    class( t_emf_solver_ndfx ), intent(inout) :: this
    integer, intent(in), dimension(:) :: gix_pos
    integer, intent(in) :: coords

    select case (p_x_dim)
    case(1)
        if ( mpi_node() == 0 ) then
            write(0,*) "(*error*) NDFX EMF solver not implemented in 1D"
            write(0,*) "(*error*) aborting..."
        endif
    case(2)
        if ( coords == p_cylindrical_b ) then
            if ( mpi_node() == 0 ) then
                write(0,*) "(*error*) NDFX EMF solver not implemented in 2D cylindrical geometry"
                write(0,*) "(*error*) aborting..."
            endif
            call abort_program( p_err_notimplemented )
        else
            this % dedt => dedt_2d_ndfx
            this % dbdt => dbdt_2d_ndfx
        endif
    case(3)
        this % dedt => dedt_3d_ndfx
        this % dbdt => dbdt_3d_ndfx
    end select

end subroutine init_ndfx

subroutine read_input_ndfx( this, input_file, dx, dt )

    implicit none

    class( t_emf_solver_ndfx ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file
    real(p_double), dimension(:), intent(in) :: dx
    real(p_double), intent(in) :: dt

    real :: dummy_
    integer :: ierr

    namelist /nl_emf_solver/ dummy_

    dummy_ = huge(1.0)

    call get_namelist( input_file, "nl_emf_solver", ierr )

    if ( ierr == 0 ) then
        read (input_file%nml_text, nml = nl_emf_solver, iostat = ierr)
        if (ierr /= 0) then
            if ( mpi_node() == 0 ) then
                write(0,*) ""
                write(0,*) "   Error reading NDFX (Pukhov) EMF solver parameters"
                write(0,*) "   aborting..."
            endif
            stop
        endif
    endif

    if ( dummy_ /= huge(1.0) ) then
        if ( mpi_node() == 0 ) then
            write(0,*) ""
            write(0,*) "   Error reading NDFX (Pukhov) EMF solver parameters"
            write(0,*) "   This solver does not accept any input parameters."
            write(0,*) "   aborting..."
        endif
        stop
    endif

end subroutine read_input_ndfx

subroutine advance_ndfx( this, e, b, jay, dt, bnd_con )
  
    implicit none
  
    class( t_emf_solver_ndfx ), intent(inout) :: this
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
    
end subroutine advance_ndfx

end module