!-----------------------------------------------------------------------------------------
! Lehe Solver
! R. Lehe, et. al., "Numerical growth of emittance in simulations of laser-wakefield
!    acceleration", Physical Review Special Topics-Accelerators and Beams, vol. 16, no. 2,
!    p. 021301, Feb. 2013.
!-----------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_solver_lehe

use m_parameters
use m_emf_define
use m_vdf_define
use m_input_file

use m_math, only: pi

private  

type, extends( t_emf_solver ), public :: t_emf_solver_lehe

    integer :: gir_pos
    procedure(dedt_lehe), pointer :: dedt => null()
    procedure(dbdt_lehe), pointer :: dbdt => null()

contains

    procedure :: init => init_lehe
    procedure :: advance => advance_lehe
    procedure :: min_gc => min_gc_lehe
    procedure :: name => name_lehe
    procedure :: read_input => read_input_lehe
    procedure :: test_stability => test_stability_lehe

end type t_emf_solver_lehe

interface
    subroutine dedt_lehe( this, e, b, jay, dt )
        import t_emf_solver_lehe, t_vdf, p_double
        class( t_emf_solver_lehe ), intent(in) :: this
        type( t_vdf ), intent(inout) :: e
        type( t_vdf ), intent(in) :: b, jay
        real(p_double), intent(in) :: dt
    end subroutine
end interface

interface
    subroutine dbdt_lehe( this, b, e, dt )
        import t_emf_solver_lehe, t_vdf, p_double
        class( t_emf_solver_lehe ), intent(in) :: this
        type( t_vdf ), intent(inout) :: b
        type( t_vdf ), intent(in) :: e
        real(p_double), intent(in) :: dt
    end subroutine
end interface

contains

function name_lehe( this )
  
    implicit none
  
    class( t_emf_solver_lehe ), intent(in) :: this
    character( len = p_maxlen_emf_solver_name ) :: name_lehe
  
    name_lehe = "Lehe"
  
end function name_lehe

subroutine dbdt_2d_lehe( this, b, e, dt  )

    implicit none

    class( t_emf_solver_lehe ), intent(in) :: this
    type( t_vdf ), intent( inout ) :: b
    type( t_vdf ), intent(in) :: e
    real(p_double), intent(in) :: dt

    real(p_k_fld) :: Bxy, Byx, Dx, Dy, Ax, Ay
    real(p_k_fld) :: dtdx1, dtdx2
    integer :: i1, i2

    dtdx1 = real(dt/b%dx_(1), p_k_fld) !*c*c !c=1
    dtdx2 = real(dt/b%dx_(2), p_k_fld) !*c*c !c=1

    Bxy = real( 0.125_p_double * ( b%dx_(1) / b%dx_(2) )**2, p_k_fld )
    Byx = 0.125_p_k_fld
    Dx = real( 0.25_p_double * ( 1.0_p_double - &
        ( sin(pi*dtdx1) / (2.0*dtdx1) )**2 ) , p_k_fld )
    Dy = 0.0_p_k_fld
    Ax = real(1.0_p_double - 2.0_p_double*Bxy - &
        3.0_p_double*Dx, p_k_fld )
    Ay = real(1.0_p_double - 2.0_p_double*Byx - &
        3.0_p_double*Dy, p_k_fld )

    !$omp parallel do private(i1)
    do i2 = -1, b%nx_(2)+2
        do i1 = -1, b%nx_(1)+2

            ! B1
            b%f2( 1, i1, i2 ) = b%f2( 1, i1, i2 ) &
                            - dtdx2 * ( Ay * ( e%f2( 3, i1  , i2+1 ) - e%f2( 3, i1  , i2 )) &
                                    + Byx * ( e%f2( 3, i1+1, i2+1 ) - e%f2( 3, i1+1, i2 ) &
                                            + e%f2( 3, i1-1, i2+1 ) - e%f2( 3, i1-1, i2 )) &
                                    + Dy * ( e%f2( 3, i1  , i2+2 ) - e%f2( 3, i1  , i2-1 )) )

            !B2
            b%f2( 2, i1, i2 ) = b%f2( 2, i1, i2 ) &
                            + dtdx1 * ( Ax * ( e%f2( 3, i1+1, i2   ) - e%f2( 3, i1, i2   )) &
                                    + Bxy * ( e%f2( 3, i1+1, i2+1 ) - e%f2( 3, i1, i2+1 ) &
                                            + e%f2( 3, i1+1, i2-1 ) - e%f2( 3, i1, i2-1 )) &
                                    + Dx * ( e%f2( 3, i1+2, i2   ) - e%f2( 3, i1-1, i2   )) )

            !B3
            b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) &
                        - dtdx1 * ( Ax * ( e%f2( 2, i1+1, i2   ) - e%f2( 2, i1, i2   ))  &
                                + Bxy * ( e%f2( 2, i1+1, i2+1 ) - e%f2( 2, i1, i2+1 )   &
                                        + e%f2( 2, i1+1, i2-1 ) - e%f2( 2, i1, i2-1 ))  &
                                + Dx * ( e%f2( 2, i1+2, i2   ) - e%f2( 2, i1-1, i2 ))) &
                        + dtdx2 * ( Ay * ( e%f2( 1, i1  , i2+1 ) - e%f2( 1, i1  , i2 ))  &
                                + Byx * ( e%f2( 1,i1+1, i2+1 ) - e%f2( 1, i1+1, i2 )   &
                                        + e%f2( 1, i1-1, i2+1 ) - e%f2( 1, i1-1, i2 ))  &
                                + Dy * ( e%f2( 1, i1  , i2+2 ) - e%f2( 1, i1  , i2-1 )) )

        enddo
    enddo
    !$omp end parallel do

end subroutine dbdt_2d_lehe

subroutine dbdt_3d_lehe( this, b, e, dt )

    implicit none

    class( t_emf_solver_lehe ), intent(in) :: this
    type( t_vdf ), intent( inout ) :: b
    type( t_vdf ), intent(in) :: e
    real(p_double), intent(in)    :: dt

    real(p_k_fld) :: Bxy, Byx, Bxz, Bzx, Byz, Bzy, Dx, Dy, Dz, Ax, Ay, Az

    real(p_k_fld) ::  dtdx1,  dtdx2,  dtdx3
    integer :: i1, i2, i3


    dtdx1 = real( dt/b%dx_(1), p_k_fld )
    dtdx2 = real( dt/b%dx_(2), p_k_fld )
    dtdx3 = real( dt/b%dx_(3), p_k_fld )


    Bxy = real( 0.125_p_double * ( b%dx_(1) / b%dx_(2) )**2, p_k_fld )
    Byx = 0.125_p_k_fld
    Bxz = real( 0.125_p_double * ( b%dx_(1) / b%dx_(3) )**2, p_k_fld )
    Bzx = 0.125_p_k_fld
    Bzy = 0.0_p_k_fld
    Byz = 0.0_p_k_fld
    Dx = real( 0.25_p_double * ( 1.0_p_double - &
        ( sin(pi*dtdx1) / (2.0*dtdx1) )**2 ) , p_k_fld )
    Dy = 0.0_p_k_fld
    Dz = 0.0_p_k_fld
    Ax = real(1.0_p_double - 2.0_p_double*Bxy - 2.0_p_double*Bxz &
        - 3.0_p_double*Dx, p_k_fld )
    Ay = real(1.0_p_double - 2.0_p_double*Byx - 2.0_p_double*Byz &
        - 3.0_p_double*Dy, p_k_fld )
    Az = real(1.0_p_double - 2.0_p_double*Bzx - 2.0_p_double*Bzy &
        - 3.0_p_double*Dz, p_k_fld )

    !$omp parallel do private(i2,i1)
    do i3 = -1, b%nx_(3)+2
        do i2 = -1, b%nx_(2)+2
            do i1 = -1, b%nx_(1)+2

            !B1
            b%f3( 1, i1, i2, i3 ) = b%f3( 1, i1, i2, i3 ) &
                    - dtdx2 * ( Ay * ( e%f3( 3, i1  , i2+1, i3   ) - e%f3( 3, i1  , i2  , i3   )) &
                            + Byx * ( e%f3( 3, i1+1, i2+1, i3   ) - e%f3( 3, i1+1, i2  , i3   )  &
                                    + e%f3( 3, i1-1, i2+1, i3   ) - e%f3( 3, i1-1, i2  , i3   ))  &
                            + Byz * ( e%f3( 3, i1  , i2+1, i3+1 ) - e%f3( 3, i1  , i2  , i3+1 )  &
                                    + e%f3( 3, i1  , i2+1, i3-1 ) - e%f3( 3, i1  , i2  , i3-1 )) &
                            + Dy * ( e%f3( 3, i1  , i2+2, i3   ) - e%f3( 3, i1  , i2-1  , i3   )) )    &
                    + dtdx3 * ( Az * ( e%f3( 2, i1  , i2  , i3+1 ) - e%f3( 2, i1  , i2  , i3   )) &
                            + Bzy * ( e%f3( 2, i1  , i2+1, i3+1 ) - e%f3( 2, i1  , i2+1, i3   )  &
                                    + e%f3( 2, i1  , i2-1, i3+1 ) - e%f3( 2, i1  , i2-1, i3   ))  &
                            + Bzx * ( e%f3( 2, i1+1, i2  , i3+1 ) - e%f3( 2, i1+1, i2  , i3   )  &
                                    + e%f3( 2, i1-1, i2  , i3+1 ) - e%f3( 2, i1-1, i2  , i3   )) &
                            + Dz * ( e%f3( 2, i1  , i2  , i3+2 ) - e%f3( 2, i1  , i2  , i3-1 )) )

            !B2
            b%f3( 2, i1, i2, i3 ) = b%f3( 2, i1, i2, i3 ) &
                    + dtdx1 * ( Ax * ( e%f3( 3, i1+1, i2  , i3   ) - e%f3( 3, i1  , i2  , i3   )) &
                            + Bxy * ( e%f3( 3, i1+1, i2+1, i3   ) - e%f3( 3, i1  , i2+1, i3   )  &
                                    + e%f3( 3, i1+1, i2-1, i3   ) - e%f3( 3, i1  , i2-1, i3   ))  &
                            + Bxz * ( e%f3( 3, i1+1, i2  , i3+1 ) - e%f3( 3, i1  , i2  , i3+1 )  &
                                    + e%f3( 3, i1+1, i2  , i3-1 ) - e%f3( 3, i1  , i2  , i3-1 )) &
                            + Dx * ( e%f3( 3, i1+2, i2  , i3   ) - e%f3( 3, i1-1, i2  , i3   )) ) &
                    - dtdx3 * ( Az * ( e%f3( 1, i1  , i2  , i3+1 ) - e%f3( 1, i1  , i2  , i3   )) &
                            + Bzx * ( e%f3( 1, i1+1, i2  , i3+1 ) - e%f3( 1, i1+1, i2  , i3   )  &
                                    + e%f3( 1, i1-1, i2  , i3+1 ) - e%f3( 1, i1-1, i2  , i3   ))  &
                            + Bzy * ( e%f3( 1, i1  , i2+1, i3+1 ) - e%f3( 1, i1  , i2+1, i3   )  &
                                    + e%f3( 1, i1  , i2-1, i3+1 ) - e%f3( 1, i1  , i2-1, i3   )) &
                            + Dz * ( e%f3( 1, i1  , i2  , i3+2 ) - e%f3( 1, i1  , i2  , i3-1  )) )

            !B3
            b%f3( 3, i1, i2, i3 ) = b%f3( 3, i1, i2, i3 ) &
                    - dtdx1 * ( Ax * ( e%f3( 2, i1+1, i2  , i3   ) - e%f3( 2, i1  , i2  , i3   )) &
                            + Bxy * ( e%f3( 2, i1+1, i2+1, i3   ) - e%f3( 2, i1  , i2+1, i3   )  &
                                    + e%f3( 2, i1+1, i2-1, i3   ) - e%f3( 2, i1  , i2-1, i3   ))  &
                            + Bxz * ( e%f3( 2, i1+1, i2  , i3+1 ) - e%f3( 2, i1  , i2  , i3+1 )  &
                                    + e%f3( 2, i1+1, i2  , i3-1 ) - e%f3( 2, i1  , i2  , i3-1 )) &
                            + Dx * ( e%f3( 2, i1+2, i2  , i3   ) - e%f3( 2, i1-1, i2  , i3   )) ) &
                    + dtdx2 * ( Ay * ( e%f3( 1, i1  , i2+1, i3   ) - e%f3( 1, i1  , i2  , i3   )) &
                            + Byx * ( e%f3( 1, i1+1, i2+1, i3   ) - e%f3( 1, i1+1, i2  , i3   )  &
                                    + e%f3( 1, i1-1, i2+1, i3   ) - e%f3( 1, i1-1, i2  , i3   ))  &
                            + Byz * ( e%f3( 1, i1  , i2+1, i3+1 ) - e%f3( 1, i1  , i2  , i3+1 )  &
                                    + e%f3( 1, i1  , i2+1, i3-1 ) - e%f3( 1, i1  , i2  , i3-1 )) &
                            + Dy * ( e%f3( 1, i1  , i2+2, i3   ) - e%f3( 1, i1  , i2-1, i3   )) )

            enddo
        enddo
    enddo
    !$omp end parallel do

end subroutine dbdt_3d_lehe

subroutine dedt_2d_lehe( this, e, b, jay, dt )

    implicit none

    class( t_emf_solver_lehe ), intent(in) :: this
    type( t_vdf ), intent(inout) :: e
    type( t_vdf ),    intent(in) :: b, jay
    real(p_double),               intent(in) :: dt

    integer :: i1, i2
    real(p_k_fld) :: dtdx1, dtdx2, dtif

    dtdx1 = real( dt/e%dx_(1), p_k_fld )
    dtdx2 = real( dt/e%dx_(2), p_k_fld )
    dtif = real( dt, p_k_fld )

    !$omp parallel do private(i1)
    do i2 = 0, e%nx_(2)+2
        do i1 = 0, e%nx_(1)+2
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
                            + dtdx1 * ( b%f2(2, i1, i2) - b%f2(2, i1-1, i2)) &
                            - dtdx2 * ( b%f2(1, i1, i2) - b%f2(1, i1, i2-1))
        enddo
    enddo
    !$omp end parallel do

end subroutine dedt_2d_lehe

subroutine dedt_3d_lehe( this, e, b, jay, dt )

    implicit none

    class( t_emf_solver_lehe ), intent(in) :: this
    type( t_vdf ), intent(inout) :: e
    type( t_vdf ), intent(in) :: b, jay

    real(p_double), intent(in) :: dt

    real(p_k_fld) :: dtdx1, dtdx2, dtdx3, dtif
    integer :: i1, i2, i3

    dtdx1 = real( dt/e%dx_(1), p_k_fld )
    dtdx2 = real( dt/e%dx_(2), p_k_fld )
    dtdx3 = real( dt/e%dx_(3), p_k_fld )
    dtif = real( dt, p_k_fld )

    !$omp parallel do private(i2,i1)
    do i3 = 0, e%nx_(3)+2
        do i2 = 0, e%nx_(2)+2
            do i1 = 0, e%nx_(1)+2

            e%f3( 1, i1, i2, i3 ) = e%f3( 1, i1, i2, i3 ) &
            - dtif  * jay%f3( 1, i1, i2, i3 )   &
            + dtdx2 * ( b%f3( 3, i1, i2, i3 ) - b%f3( 3, i1, i2-1, i3 ) ) &
            - dtdx3 * ( b%f3( 2, i1, i2, i3 ) - b%f3( 2, i1, i2, i3-1 ) )

            e%f3( 2, i1, i2, i3 ) = e%f3( 2, i1, i2, i3 ) &
            - dtif  * jay%f3( 2, i1, i2, i3 ) &
            + dtdx3 * ( b%f3( 1, i1, i2, i3 ) - b%f3( 1, i1, i2, i3-1 ) ) &
            - dtdx1 * ( b%f3( 3, i1, i2, i3 ) - b%f3( 3, i1-1, i2, i3 ) )

            e%f3( 3, i1, i2, i3 ) = e%f3( 3, i1, i2, i3 ) &
            - dtif  * jay%f3( 3, i1, i2, i3 )   &
            + dtdx1 * ( b%f3( 2, i1, i2, i3 ) - b%f3( 2, i1-1, i2, i3 ) ) &
            - dtdx2 * ( b%f3( 1, i1, i2, i3 ) - b%f3( 1, i1, i2-1, i3 ) )

            enddo
        enddo
    enddo
    !$omp end parallel do

end subroutine dedt_3d_lehe

subroutine test_stability_lehe( this, dt, dx )
    
    implicit none

    class( t_emf_solver_lehe ), intent(in) :: this
    real(p_double), intent(in)  ::  dt
    real(p_double), dimension(:), intent(in)  ::  dx

    real(p_double) :: cour
    
    cour = dx(1)

    if (dt > cour) then
        if ( mpi_node() == 0 ) then
           print *, '(*error*) Lehe solver stability condition violated, aborting'
           print *, 'dx   = ', dx(1:p_x_dim)
           print *, 'dt   = ', dt, ' max(dt) = ', cour
        endif
        call abort_program(p_err_invalid)
    endif

end subroutine test_stability_lehe

subroutine min_gc_lehe( this, min_gc )

    implicit none

    class( t_emf_solver_lehe ), intent(in) :: this
    integer, dimension(:,:), intent(inout) :: min_gc

    min_gc(p_lower,:) = 3
    min_gc(p_upper,:) = 4

end subroutine


subroutine init_lehe( this, gix_pos, coords )

    implicit none

    class( t_emf_solver_lehe ), intent(inout) :: this
    integer, intent(in), dimension(:) :: gix_pos
    integer, intent(in) :: coords

    select case (p_x_dim)
    case(1)
        if ( mpi_node() == 0 ) then
            write(0,*) "(*error*) Lehe EMF solver not implemented in 1D"
            write(0,*) "(*error*) aborting..."
        endif
    case(2)
        if ( coords == p_cylindrical_b ) then
            if ( mpi_node() == 0 ) then
                write(0,*) "(*error*) Lehe EMF solver not implemented in 2D cylindrical geometry"
                write(0,*) "(*error*) aborting..."
            endif
            call abort_program( p_err_notimplemented )
        else
            this % dedt => dedt_2d_lehe
            this % dbdt => dbdt_2d_lehe
        endif
    case(3)
        this % dedt => dedt_3d_lehe
        this % dbdt => dbdt_3d_lehe
    end select

end subroutine init_lehe

subroutine read_input_lehe( this, input_file, dx, dt )

    implicit none

    class( t_emf_solver_lehe ), intent(inout) :: this
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
                write(0,*) "   Error reading Lehe EMF solver parameters"
                write(0,*) "   aborting..."
            endif
            stop
        endif
    endif

    if ( dummy_ /= huge(1.0) ) then
        if ( mpi_node() == 0 ) then
            write(0,*) ""
            write(0,*) "   Error reading Lehe EMF solver parameters"
            write(0,*) "   This solver does not accept any input parameters."
            write(0,*) "   aborting..."
        endif
        stop
    endif

end subroutine read_input_lehe

subroutine advance_lehe( this, e, b, jay, dt, bnd_con )
  
    implicit none
  
    class( t_emf_solver_lehe ), intent(inout) :: this
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
    
end subroutine advance_lehe

end module 