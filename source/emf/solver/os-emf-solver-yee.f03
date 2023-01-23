!-----------------------------------------------------------------------------------------
! Yee solver
! K. YEE, NUMERICAL SOLUTION OF INITIAL BOUNDARY VALUE PROBLEMS INVOLVING MAXWELLS
!   EQUATIONS IN ISOTROPIC MEDIA, IEEE Transactions on Antenna Propagation, vol. 14,
!   no. 3, pp. 302-307, 1966.
!-----------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_solver_yee

use m_parameters
use m_emf_define
use m_vdf_define
use m_input_file

private  

type, extends( t_emf_solver ), public :: t_emf_solver_yee

    integer :: gir_pos
    procedure(dedt_yee), pointer :: dedt => null()
    procedure(dbdt_yee), pointer :: dbdt => null()

contains

    procedure :: init => init_yee
    procedure :: advance => advance_yee
    procedure :: min_gc => min_gc_yee
    procedure :: name => name_yee
    procedure :: read_input => read_input_yee
    procedure :: test_stability => test_stability_yee

end type t_emf_solver_yee

interface
    subroutine dedt_yee( this, e, b, jay, dt )
        import t_emf_solver_yee, t_vdf, p_double
        class( t_emf_solver_yee ), intent(in) :: this
        type( t_vdf ), intent(inout) :: e
        type( t_vdf ), intent(in) :: b, jay
        real(p_double), intent(in) :: dt
    end subroutine
end interface

interface
    subroutine dbdt_yee( this, b, e, dt )
        import t_emf_solver_yee, t_vdf, p_double
        class( t_emf_solver_yee ), intent(in) :: this
        type( t_vdf ), intent(inout) :: b
        type( t_vdf ), intent(in) :: e
        real(p_double), intent(in) :: dt
    end subroutine
end interface


contains

function name_yee( this )
  
    implicit none
  
    class( t_emf_solver_yee ), intent(in) :: this
    character( len = p_maxlen_emf_solver_name ) :: name_yee
  
    name_yee = "Yee"
  
end function name_yee

subroutine dbdt_1d_yee( this, b, e, dt )

    implicit none

    class( t_emf_solver_yee ), intent(in) :: this
    type( t_vdf ), intent( inout ) :: b
    type( t_vdf ), intent(in) :: e
    real(p_double), intent(in) :: dt

    real(p_k_fld) :: dtdx1
    integer :: i1

    dtdx1 = real( dt/b%dx_(1), p_k_fld ) !*c*c !c=1

    !$omp parallel do
    do i1 = 0, b%nx_(1)+1
        !B1
        ! b%f1( 1, i1 ) = b%f1( 1, i1 )

        !B2
        b%f1( 2, i1 ) = b%f1( 2, i1 ) &
                            + dtdx1 * ( e%f1( 3, i1+1 ) - e%f1( 3, i1 ))

        !B3
        b%f1( 3, i1 ) = b%f1( 3, i1 ) &
                            - dtdx1 * ( e%f1( 2, i1+1 ) - e%f1( 2, i1 ))
    enddo
    !$omp end parallel do

end subroutine dbdt_1d_yee
  
subroutine dbdt_2d_yee( this, b, e, dt )

    implicit none

    class( t_emf_solver_yee ), intent(in) :: this
    type( t_vdf ), intent( inout ) :: b
    type( t_vdf ),    intent(in) :: e
    real(p_double),               intent(in) :: dt

    real(p_k_fld) :: dtdx1, dtdx2
    integer :: i1, i2

    dtdx1 = real(dt/b%dx_(1), p_k_fld) !*c*c !c=1
    dtdx2 = real(dt/b%dx_(2), p_k_fld) !*c*c !c=1

    !$omp parallel do private(i1)
    do i2 = 0, b%nx_(2)+1
        do i1 = 0, b%nx_(1)+1

        !B1
        b%f2( 1, i1, i2 ) = b%f2( 1, i1, i2 ) &
                            - dtdx2 * ( e%f2( 3, i1, i2+1 ) - e%f2( 3, i1, i2 ))

        !B2
        b%f2( 2, i1, i2 ) = b%f2( 2, i1, i2 ) &
                            + dtdx1 * ( e%f2( 3, i1+1, i2 ) - e%f2( 3, i1, i2 ))

        !B3
        b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) &
                            - dtdx1 * ( e%f2( 2, i1+1, i2 ) - e%f2( 2, i1, i2 )) &
                            + dtdx2 * ( e%f2( 1, i1, i2+1 ) - e%f2( 1, i1, i2 ))

        enddo
    enddo
    !$omp end parallel do

end subroutine dbdt_2d_yee

subroutine dbdt_cyl_2d_yee( this, b, e, dt )

    implicit none

    class( t_emf_solver_yee ), intent(in) :: this
    type( t_vdf ), intent( inout ) :: b
    type( t_vdf ), intent( in ) :: e
    real(p_double),  intent(in) :: dt

    real(p_k_fld) :: dtdz, dtdr
    real(p_k_fld) :: tmp, rm, rp
    integer :: i1, i2, i2_0, gshift_i2

    dtdz = real(dt/b%dx_(1), p_k_fld)
    dtdr = real(dt/b%dx_(2), p_k_fld)

    gshift_i2 = this % gir_pos - 2

    if ( gshift_i2 < 0 ) then
        ! This node contains the cylindrical axis so start the solver in cell 2
        i2_0 = 2
    else
        i2_0 = 0
    endif

    do i2 = i2_0, b%nx_(2)+1

        rp   = i2 + gshift_i2 + 0.5     ! position of the upper edge of the cell (normalized to dr)
        rm   = i2 + gshift_i2 - 0.5     ! position of the lower edge of the cell (normalized to dr)
        tmp  = dtdr/(i2 + gshift_i2)    ! (dt/dr) / position of the middle of the cell (normalized to dr)

        do i1 = 0, b%nx_(1)+1
        !B1
        b%f2( 1, i1, i2 ) = b%f2(1, i1, i2) - tmp * ( rp*e%f2(3, i1, i2+1) - rm*e%f2(3, i1, i2))

        !B2
        b%f2( 2, i1, i2 ) = b%f2(2, i1, i2) + dtdz * ( e%f2( 3, i1+1, i2 ) - e%f2( 3, i1, i2 ))

        !B3
        b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) - &
                                    dtdz * ( e%f2( 2, i1+1, i2 ) - e%f2( 2, i1, i2 )) + &
                                    dtdr * ( e%f2( 1, i1, i2+1 ) - e%f2( 1, i1, i2 ))
        enddo
    enddo

    if ( gshift_i2 < 0 ) then

        ! guard cell 1 (i2 = 0)
        do i1 = 0, b%nx_(1)+1
        b%f2( 1, i1, 0 ) =   b%f2( 1, i1, 2 )
        b%f2( 2, i1, 0 ) = - b%f2( 2, i1, 3 )
        b%f2( 3, i1, 0 ) = - b%f2( 3, i1, 2 )
        enddo

        ! axial cell (i2 = 1)
        do i1 = 0, b%nx_(1)+1
        ! rmid = 0 in this cell so the rot E component is calculated differently
        b%f2( 1, i1, 1 ) =   b%f2(1, i1, 1) - 4 * dtdr * e%f2(3, i1, 2)

        ! B2 is not defined on axis so just reflect the corresponding value inside the box
        b%f2( 2, i1, 1 ) = - b%f2( 2, i1, 2 )

        ! B3 is zero on axis
        b%f2( 3, i1, 1 ) = 0.0_p_k_fld
        enddo

    endif

end subroutine dbdt_cyl_2d_yee

  
subroutine dbdt_3d_yee( this, b, e, dt )

    implicit none

    class( t_emf_solver_yee ), intent(in) :: this
    type( t_vdf ), intent( inout ) :: b
    type( t_vdf ), intent(in) :: e
    real(p_double), intent(in) :: dt

    real(p_k_fld) ::  dtdx1,  dtdx2,  dtdx3
    integer :: i1, i2, i3

    dtdx1 = real( dt/b%dx_(1), p_k_fld )
    dtdx2 = real( dt/b%dx_(2), p_k_fld )
    dtdx3 = real( dt/b%dx_(3), p_k_fld )

    !$omp parallel do private(i2,i1)
    do i3 = 0, b%nx_(3)+1
        do i2 = 0, b%nx_(2)+1
            do i1 = 0, b%nx_(1)+1

                !B1
                b%f3( 1, i1, i2, i3 ) = b%f3( 1, i1, i2, i3 ) &
                                    - dtdx2 * ( e%f3( 3, i1, i2+1, i3 ) - e%f3( 3, i1, i2, i3 )) &
                                    + dtdx3 * ( e%f3( 2, i1, i2, i3+1 ) - e%f3( 2, i1, i2, i3 ))

                !B2
                b%f3( 2, i1, i2, i3 ) = b%f3( 2, i1, i2, i3 ) &
                                    + dtdx1 * ( e%f3( 3, i1+1, i2, i3 ) - e%f3( 3, i1, i2, i3 )) &
                                    - dtdx3 * ( e%f3( 1, i1, i2, i3+1 ) - e%f3( 1, i1, i2, i3 ))
                !B3
                b%f3( 3, i1, i2, i3 ) = b%f3( 3, i1, i2, i3 ) &
                                    - dtdx1 * ( e%f3( 2, i1+1, i2, i3 ) - e%f3( 2, i1, i2, i3 )) &
                                    + dtdx2 * ( e%f3( 1, i1, i2+1, i3 ) - e%f3( 1, i1, i2, i3 ))

            enddo
        enddo
    enddo
    !$omp end parallel do
  
end subroutine dbdt_3d_yee

subroutine dedt_1d_yee( this, e, b, jay, dt )

    implicit none
    
    class( t_emf_solver_yee ), intent(in) :: this
    type( t_vdf ), intent(inout) :: e
    type( t_vdf ),    intent(in) :: b, jay

    real(p_double),               intent(in) :: dt
    real(p_k_fld) :: dtdx1, dtif
    integer :: i1

    dtdx1 = real( dt/e%dx_(1), p_k_fld )
    dtif = real( dt, p_k_fld )

    !$omp parallel do
    do i1 = 1, e%nx_(1)+1
        ! E1
        e%f1(1, i1) = e%f1(1, i1) - dtif * jay%f1(1, i1)

        ! E2
        e%f1(2, i1) = e%f1(2, i1) &
                        - dtif * jay%f1(2, i1) &
                        - dtdx1 * ( b%f1(3, i1) - b%f1(3, i1-1) )

        ! E3
        e%f1(3, i1) = e%f1(3, i1) &
                        - dtif * jay%f1(3, i1) &
                        + dtdx1 * ( b%f1(2, i1) - b%f1(2, i1-1))
    enddo
    !$omp end parallel do

end subroutine dedt_1d_yee

subroutine dedt_2d_yee( this, e, b, jay, dt )

    implicit none

    class( t_emf_solver_yee ), intent(in) :: this
    type( t_vdf ), intent(inout) :: e
    type( t_vdf ),    intent(in) :: b, jay

    real(p_double),               intent(in) :: dt

    integer :: i1, i2
    real(p_k_fld) :: dtdx1, dtdx2, dtif

    dtdx1 = real( dt/e%dx_(1), p_k_fld )
    dtdx2 = real( dt/e%dx_(2), p_k_fld )
    dtif = real( dt, p_k_fld )

    !$omp parallel do private(i1)
    do i2 = 1, e%nx_(2)+1
        do i1 = 1, e%nx_(1)+1
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

end subroutine dedt_2d_yee

subroutine dedt_cyl_2d_yee( this, e, b, jay, dt )

    implicit none

    class( t_emf_solver_yee ), intent(in) :: this
    type( t_vdf ), intent(inout) :: e
    type( t_vdf ), intent(in) :: b, jay
    real(p_double), intent(in) :: dt

    real(p_k_fld) :: dtdz, dtdr, dtif
    real(p_k_fld) :: tmp, rcp, rcm
    integer :: i1, i2, gshift_i2

    dtdz = real(dt/b%dx_(1), p_k_fld)
    dtdr = real(dt/b%dx_(2), p_k_fld)
    dtif = real( dt, p_k_fld )

    gshift_i2 = this % gir_pos - 2

    do i2 = 1, e%nx_(2)+1

        tmp = dtdr / (( i2 + gshift_i2 ) - 0.5)
        rcp = ( i2 + gshift_i2 )
        rcm = ( i2 + gshift_i2 - 1 )

        do i1 = 1, e%nx_(1)+2
            ! E1
            e%f2(1, i1, i2) = e%f2(1, i1, i2) &
                            - dtif * jay%f2(1, i1, i2) &
                            + tmp * ( rcp * b%f2(3, i1, i2) - rcm * b%f2(3, i1, i2-1) )

            ! E2
            e%f2(2, i1, i2) = e%f2(2, i1, i2) &
                            - dtif * jay%f2(2, i1, i2) &
                            - dtdz * ( b%f2(3, i1, i2) - b%f2(3, i1-1, i2) )

            ! E3
            e%f2(3, i1, i2) = e%f2(3, i1, i2) &
                            - dtif * jay%f2(3, i1, i2) &
                            + dtdz * ( b%f2(2, i1, i2) - b%f2(2, i1-1, i2)) &
                            - dtdr * ( b%f2(1, i1, i2) - b%f2(1, i1, i2-1))
        enddo
    enddo

    ! Correct values on axial boundary. (there is a little redundancy between this and update_boundary
    ! it needs to be checked)
    if ( gshift_i2 < 0 ) then
        ! axial cell (i2 = 1)
        do i1 = 1, e%nx_(1)+2
            e%f2(1, i1, 1) =   e%f2(1, i1, 2)
            e%f2(2, i1, 1) =   0
            e%f2(3, i1, 1) = - e%f2(3, i1, 2)
        enddo
    endif

end subroutine dedt_cyl_2d_yee

subroutine dedt_3d_yee( this, e, b, jay, dt )

    implicit none

    class( t_emf_solver_yee ), intent(in) :: this
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
    do i3 = 1, e%nx_(3)+1
        do i2 = 1, e%nx_(2)+1
            do i1 = 1, e%nx_(1)+1

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

end subroutine dedt_3d_yee

subroutine test_stability_yee( this, dt, dx )
    
    implicit none

    class( t_emf_solver_yee ), intent(in) :: this
    real(p_double), intent(in)  ::  dt
    real(p_double), dimension(:), intent(in)  ::  dx

    real(p_double) :: cour
    integer :: i

    cour = 0.0
    do i = 1, p_x_dim
        cour = cour + 1.0/(dx(i))**2
    enddo
    cour = sqrt(1.0/cour)

    if (dt > cour) then
        if ( mpi_node() == 0 ) then
           print *, '(*error*) Yee EMF solver stability condition violated, aborting'
           print *, 'dx   = ', dx(1:p_x_dim)
           print *, 'dt   = ', dt, ' max(dt) = ', cour
        endif
        call abort_program(p_err_invalid)
    endif

end subroutine test_stability_yee

subroutine min_gc_yee( this, min_gc )

    implicit none

    class( t_emf_solver_yee ), intent(in) :: this
    integer, dimension(:,:), intent(inout) :: min_gc

    min_gc(p_lower,:) = 1
    min_gc(p_upper,:) = 2

end subroutine

subroutine init_yee( this, gix_pos, coords )

    implicit none

    class( t_emf_solver_yee ), intent(inout) :: this
    integer, intent(in), dimension(:) :: gix_pos
    integer, intent(in) :: coords

    select case (p_x_dim)
    case(1)
        this % dedt => dedt_1d_yee
        this % dbdt => dbdt_1d_yee
    case(2)
        if ( coords == p_cylindrical_b ) then
            ! Cylindrical (r-z) coordinates solvers
            this % dedt => dedt_cyl_2d_yee
            this % dbdt => dbdt_cyl_2d_yee
            ! Radial position of lower cell
            this % gir_pos = gix_pos(2)
        else
            this % dedt => dedt_2d_yee
            this % dbdt => dbdt_2d_yee
        endif
    case(3)
        this % dedt => dedt_3d_yee
        this % dbdt => dbdt_3d_yee
    end select

end subroutine init_yee

subroutine read_input_yee( this, input_file, dx, dt )

    implicit none

    class( t_emf_solver_yee ), intent(inout) :: this
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
                write(0,*) "   Error reading Yee EMF solver parameters"
                write(0,*) "   aborting..."
            endif
            stop
        endif
    endif

    if ( dummy_ /= huge(1.0) ) then
        if ( mpi_node() == 0 ) then
            write(0,*) ""
            write(0,*) "   Error reading Yee EMF solver parameters"
            write(0,*) "   This solver does not accept any input parameters."
            write(0,*) "   aborting..."
        endif
        stop
    endif

end subroutine read_input_yee

subroutine advance_yee( this, e, b, jay, dt, bnd_con )
  
    implicit none
  
    class( t_emf_solver_yee ), intent(inout) :: this
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
    
end subroutine advance_yee

end module m_emf_solver_yee