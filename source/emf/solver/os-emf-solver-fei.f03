!-----------------------------------------------------------------------------------------
! Customized (Fei) Solver
!-----------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_solver_fei

#include "memory/memory.h"

use m_parameters
use m_emf_define
use m_vdf_define
use m_input_file

use m_fft
use m_math, only : pi

private  

integer, parameter, public :: p_max_n_coef = 128
integer, parameter :: p_max_lst_len = 32
integer, parameter :: p_err_integral = -1

integer, parameter, public :: p_emf_fei_std = 0
integer, parameter, public :: p_emf_fei_bump = 1
integer, parameter, public :: p_emf_fei_xu = 2
integer, parameter, public :: p_emf_fei_dual = 3
integer, parameter, public :: p_emf_fei_coef = 99

type, extends( t_emf_solver ), public :: t_emf_solver_fei

    integer :: gir_pos
    procedure(dedt_fei), pointer :: dedt => null()
    procedure(dbdt_fei), pointer :: dbdt => null()

    ! Solver parameters
    integer :: type
    real(p_k_fld), allocatable, dimension(:) :: coef_e, coef_b
    real(p_k_fld), allocatable, dimension(:,:) :: taper_coef_e, taper_coef_b
    integer :: n_coef, solver_ord, weight_n
    real(p_double) :: ku, kl, dk, weight_w

    ! Current filtering
    integer :: n_damp_cell, taper_order
    real(p_k_fld) :: filter_limit, filter_width
    logical :: correct_current, filter_current, taper_bnd
    procedure(filter_corr_current_fei), pointer :: filter_corr_current => null()
    real(p_k_fld), dimension(:,:,:), pointer :: f_buffer2 => null()
    real(p_k_fld), dimension(:,:,:,:), pointer :: f_buffer3 => null()

contains

    procedure :: init => init_fei
    procedure :: advance => advance_fei
    procedure :: min_gc => min_gc_fei
    procedure :: name => name_fei
    procedure :: read_input => read_input_fei
    procedure :: test_stability => test_stability_fei
    procedure :: cleanup => cleanup_fei

end type t_emf_solver_fei

interface
    subroutine dedt_fei( this, e, b, jay, dt, bnd_type )
        import t_emf_solver_fei, t_vdf, p_double, p_max_dim
        class( t_emf_solver_fei ), intent(in) :: this
        type( t_vdf ), intent(inout) :: e
        type( t_vdf ), intent(in) :: b, jay
        real(p_double), intent(in) :: dt
        integer, intent(in), dimension(2,p_max_dim) :: bnd_type
    end subroutine
end interface

interface
    subroutine dbdt_fei( this, b, e, dt, bnd_type )
        import t_emf_solver_fei, t_vdf, p_double, p_max_dim
        class( t_emf_solver_fei ), intent(in) :: this
        type( t_vdf ), intent(inout) :: b
        type( t_vdf ), intent(in) :: e
        real(p_double), intent(in) :: dt
        integer, intent(in), dimension(2,p_max_dim) :: bnd_type
    end subroutine
end interface

interface
    subroutine filter_corr_current_fei( this, jay )
        import t_emf_solver_fei, t_vdf
        class( t_emf_solver_fei ), intent(inout) :: this
        type( t_vdf ), intent(inout) :: jay
    end subroutine
end interface

interface
    function kernel( x, arglst )
        import p_double
        real(p_double), intent(in) :: x
        real(p_double), intent(in), dimension(:) :: arglst
        real(p_double) :: kernel
    end function
end interface

public :: std_coef, gen_solver_coef, kernel_xu, kernel_dual_e, kernel_dual_b, factorial

contains

subroutine dbdt_2d_fei( this, b, e, dt, bnd_type )
    
    implicit none

    class( t_emf_solver_fei ), intent(in) :: this
    type(t_vdf), intent(inout)                :: b
    type(t_vdf), intent(in)                   :: e
    real(p_double), intent(in)                :: dt
    integer, intent(in), dimension(2,p_max_dim) :: bnd_type

    real(p_k_fld) :: dtdx1, dtdx2, tmp
    integer :: i1, i2, j, i1_lb, i1_ub, offset, ord_hf

    dtdx1 = real( dt/b%dx(1), p_k_fld )
    dtdx2 = real( dt/b%dx(2), p_k_fld )

    i1_lb = - this % n_coef + 1
    i1_ub = b%nx_(1) + this % n_coef

    if ( this % taper_bnd ) then

        offset = this%taper_order / 2

        if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
             bnd_type(p_lower, 1) /= p_bc_other_node ) then
            i1_lb = offset
        endif

        if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
            bnd_type(p_upper, 1) /= p_bc_other_node ) then
            i1_ub = b%nx_(1) - offset + 1
        endif

    endif

    !$omp parallel do private(i1,tmp,j)
    do i2 = 0, b%nx_(2)+1
        do i1 = i1_lb, i1_ub

            !B1
            b%f2(1, i1, i2) = b%f2(1, i1, i2) &
                                - dtdx2 * ( e%f2(3, i1, i2+1) - e%f2(3, i1, i2))

            !B2
            tmp = 0.0_p_k_fld
            do j = 1, this % n_coef
                tmp = tmp + this % coef_e(j) * ( e%f2(3, i1+j, i2) - e%f2(3, i1-j+1, i2) )
            enddo
            b%f2(2, i1, i2) = b%f2(2, i1, i2) + dtdx1 * tmp

            !B3
            tmp = 0.0_p_k_fld
            do j = 1, this % n_coef
                tmp = tmp + this % coef_e(j) * ( e%f2(2, i1+j, i2) - e%f2(2, i1-j+1, i2) )
            enddo
            b%f2(3, i1, i2) = b%f2(3, i1, i2) - dtdx1 * tmp &
                            + dtdx2 * ( e%f2(1, i1, i2+1) - e%f2(1, i1, i2) )


        enddo
    enddo
    !$omp end parallel do

    if ( this % taper_bnd ) then

        ! advance lower boundary
        if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
             bnd_type(p_lower, 1) /= p_bc_other_node ) then

            !$omp parallel do private(i1,ord_hf,tmp,j)
            do i2 = 0, b%nx_(2)+1

                ord_hf = 0

                do i1 = 0, i1_lb-1

                    ord_hf = ord_hf + 1

                    !B1
                    b%f2(1, i1, i2) = b%f2(1, i1, i2) &
                                        - dtdx2 * ( e%f2(3, i1, i2+1) - e%f2(3, i1, i2))

                    !B2
                    tmp = 0.0_p_k_fld
                    do j = 1, ord_hf
                        tmp = tmp + this % taper_coef_e(j,ord_hf) * ( e%f2(3, i1+j, i2) - e%f2(3, i1-j+1, i2) )
                    enddo
                    b%f2(2, i1, i2) = b%f2(2, i1, i2) + dtdx1 * tmp

                    !B3
                    tmp = 0.0_p_k_fld
                    do j = 1, ord_hf
                        tmp = tmp + this % taper_coef_e(j,ord_hf) * ( e%f2(2, i1+j, i2) - e%f2(2, i1-j+1, i2) )
                    enddo
                    b%f2(3, i1, i2) = b%f2(3, i1, i2) - dtdx1 * tmp &
                                    + dtdx2 * ( e%f2(1, i1, i2+1) - e%f2(1, i1, i2) )

                enddo
            enddo
            !$omp end parallel do

        endif

        ! advance upper boundary
        if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
             bnd_type(p_upper, 1) /= p_bc_other_node ) then

            !$omp parallel do private(i1,ord_hf,tmp,j)
            do i2 = 0, b%nx_(2)+1

                ord_hf = offset + 1

                do i1 = i1_ub+1, b%nx_(1)+1

                    ord_hf = ord_hf - 1

                    !B1
                    b%f2(1, i1, i2) = b%f2(1, i1, i2) &
                                        - dtdx2 * ( e%f2(3, i1, i2+1) - e%f2(3, i1, i2))

                    !B2
                    tmp = 0.0_p_k_fld
                    do j = 1, ord_hf
                        tmp = tmp + this % taper_coef_e(j,ord_hf) * ( e%f2(3, i1+j, i2) - e%f2(3, i1-j+1, i2) )
                    enddo
                    b%f2(2, i1, i2) = b%f2(2, i1, i2) + dtdx1 * tmp

                    !B3
                    tmp = 0.0_p_k_fld
                    do j = 1, ord_hf
                        tmp = tmp + this % taper_coef_e(j,ord_hf) * ( e%f2(2, i1+j, i2) - e%f2(2, i1-j+1, i2) )
                    enddo
                    b%f2(3, i1, i2) = b%f2(3, i1, i2) - dtdx1 * tmp &
                                    + dtdx2 * ( e%f2(1, i1, i2+1) - e%f2(1, i1, i2) )

                enddo
            enddo
            !$omp end parallel do

        endif

    endif

end subroutine dbdt_2d_fei

subroutine dedt_2d_fei( this, e, b, jay, dt, bnd_type )
    
    implicit none

    class( t_emf_solver_fei ), intent(in) :: this
    type(t_vdf), intent( in )                 :: b, jay
    type(t_vdf), intent(inout)                :: e
    real(p_double), intent(in)                :: dt
    integer, intent(in), dimension(2,p_max_dim) :: bnd_type

    real(p_k_fld) :: dtdx1, dtdx2, dtif, tmp
    integer :: i1, i2, j, i1_lb, i1_ub, offset, ord_hf

    dtdx1 = real( dt/e%dx(1), p_k_fld )
    dtdx2 = real( dt/e%dx(2), p_k_fld )
    dtif = real( dt, p_k_fld )

    i1_lb = 2 - this % n_coef
    i1_ub = e%nx_(1) + this % n_coef

    if ( this % taper_bnd ) then

        offset = this%taper_order / 2

        if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
             bnd_type(p_lower, 1) /= p_bc_other_node ) then
            i1_lb = offset + 1
        endif

        if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
            bnd_type(p_upper, 1) /= p_bc_other_node ) then
            i1_ub = e%nx_(1) - offset + 1
        endif

    endif

    !$omp parallel do private(i1,tmp,j)
    do i2 = 1, e%nx_(2)+1
        do i1 = i1_lb, i1_ub
            ! E1
            e%f2(1, i1, i2) = e%f2(1, i1, i2) &
                            - dtif * jay%f2(1, i1, i2) &
                            + dtdx2 * ( b%f2(3, i1, i2) - b%f2(3, i1, i2-1) )

            ! E2
            tmp = 0.0_p_k_fld
            do j = 1, this % n_coef
                tmp = tmp + this % coef_b(j) * ( b%f2(3, i1+j-1, i2) - b%f2(3, i1-j, i2) )
            enddo
            e%f2(2, i1, i2) = e%f2(2, i1, i2) - dtif * jay%f2(2, i1, i2) &
                            - dtdx1 * tmp

            ! E3
            tmp = 0.0_p_k_fld
            do j = 1, this % n_coef
                tmp = tmp + this % coef_b(j) * ( b%f2(2, i1+j-1, i2) - b%f2(2, i1-j, i2) )
            enddo
            e%f2(3, i1, i2) = e%f2(3, i1, i2) - dtif * jay%f2(3, i1, i2) &
                            - dtdx2 * ( b%f2(1, i1, i2) - b%f2(1, i1, i2-1) ) &
                            + dtdx1 * tmp

        enddo
    enddo
    !$omp end parallel do

    if ( this % taper_bnd ) then

        ! advance lower boundary
        if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
             bnd_type(p_lower, 1) /= p_bc_other_node ) then

            !$omp parallel do private(i1,ord_hf,tmp,j)
            do i2 = 1, e%nx_(2)+1

                ord_hf = 0

                do i1 = 1, i1_lb-1

                    ord_hf = ord_hf + 1

                    ! E1
                    e%f2(1, i1, i2) = e%f2(1, i1, i2) &
                                    - dtif * jay%f2(1, i1, i2) &
                                    + dtdx2 * ( b%f2(3, i1, i2) - b%f2(3, i1, i2-1) )

                    ! E2
                    tmp = 0.0_p_k_fld
                    do j = 1, ord_hf
                        tmp = tmp + this % taper_coef_b(j,ord_hf) * ( b%f2(3, i1+j-1, i2) - b%f2(3, i1-j, i2) )
                    enddo
                    e%f2(2, i1, i2) = e%f2(2, i1, i2) - dtif * jay%f2(2, i1, i2) &
                                    - dtdx1 * tmp

                    ! E3
                    tmp = 0.0_p_k_fld
                    do j = 1, ord_hf
                        tmp = tmp + this % taper_coef_b(j,ord_hf) * ( b%f2(2, i1+j-1, i2) - b%f2(2, i1-j, i2) )
                    enddo
                    e%f2(3, i1, i2) = e%f2(3, i1, i2) - dtif * jay%f2(3, i1, i2) &
                                    - dtdx2 * ( b%f2(1, i1, i2) - b%f2(1, i1, i2-1) ) &
                                    + dtdx1 * tmp

                enddo
            enddo
            !$omp end parallel do

        endif

        ! advance upper boundary
        if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
             bnd_type(p_upper, 1) /= p_bc_other_node ) then

            !$omp parallel do private(i1,ord_hf,tmp,j)
            do i2 = 1, e%nx_(2)+1

                ord_hf = offset + 1

                do i1 = i1_ub+1, e%nx_(1)+1

                    ord_hf = ord_hf - 1

                    ! E1
                    e%f2(1, i1, i2) = e%f2(1, i1, i2) &
                                    - dtif * jay%f2(1, i1, i2) &
                                    + dtdx2 * ( b%f2(3, i1, i2) - b%f2(3, i1, i2-1) )

                    ! E2
                    tmp = 0.0_p_k_fld
                    do j = 1, ord_hf
                        tmp = tmp + this % taper_coef_b(j,ord_hf) * ( b%f2(3, i1+j-1, i2) - b%f2(3, i1-j, i2) )
                    enddo
                    e%f2(2, i1, i2) = e%f2(2, i1, i2) - dtif * jay%f2(2, i1, i2) &
                                    - dtdx1 * tmp

                    ! E3
                    tmp = 0.0_p_k_fld
                    do j = 1, ord_hf
                        tmp = tmp + this % taper_coef_b(j,ord_hf) * ( b%f2(2, i1+j-1, i2) - b%f2(2, i1-j, i2) )
                    enddo
                    e%f2(3, i1, i2) = e%f2(3, i1, i2) - dtif * jay%f2(3, i1, i2) &
                                    - dtdx2 * ( b%f2(1, i1, i2) - b%f2(1, i1, i2-1) ) &
                                    + dtdx1 * tmp

                enddo
            enddo
            !$omp end parallel do

        endif

    endif

end subroutine dedt_2d_fei

subroutine dbdt_3d_fei( this, b, e, dt, bnd_type )

    implicit none

    class( t_emf_solver_fei ), intent(in) :: this
    type( t_vdf ), intent( inout )            :: b
    type( t_vdf ), intent(in)                 :: e
    real(p_double), intent(in)                :: dt
    integer, intent(in), dimension(2,p_max_dim) :: bnd_type
    
    real(p_k_fld) :: dtdx1, dtdx2, dtdx3
    integer :: i1, i2, i3, j, i1_ub, i1_lb, offset, ord_hf
    real(p_k_fld) :: tmp

    dtdx1 = real(dt/b%dx(1), p_k_fld)
    dtdx2 = real(dt/b%dx(2), p_k_fld)
    dtdx3 = real(dt/b%dx(3), p_k_fld)

    i1_lb = - this % n_coef + 1
    i1_ub = b%nx_(1) + this % n_coef

    if ( this%taper_bnd ) then

        offset = this%taper_order / 2

        if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
             bnd_type(p_lower, 1) /= p_bc_other_node ) then
            i1_lb = offset
        endif

        if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
             bnd_type(p_upper, 1) /= p_bc_other_node ) then
            i1_ub = b%nx_(1) - offset + 1
        endif

    endif

    !$omp parallel do private(i2,i1,tmp,j)
    do i3 = 0, b%nx_(3)+1
        do i2 = 0, b%nx_(2)+1
            do i1 = i1_lb, i1_ub

                !B1
                b%f3(1, i1, i2, i3) = b%f3(1, i1, i2, i3) &
                                    - dtdx2 * ( e%f3(3, i1, i2+1, i3) - e%f3(3, i1, i2, i3)) &
                                    + dtdx3 * ( e%f3(2, i1, i2, i3+1) - e%f3(2, i1, i2, i3))

                !B2
                tmp = 0.0_p_k_fld
                do j = 1, this % n_coef
                    tmp = tmp + this % coef_e(j) * ( e%f3(3, i1+j, i2, i3) - e%f3(3, i1-j+1, i2, i3) )
                enddo
                b%f3(2, i1, i2, i3) = b%f3( 2, i1, i2, i3 ) + dtdx1 * tmp &
                                    - dtdx3 * ( e%f3(1, i1, i2, i3+1) - e%f3(1, i1, i2, i3))

                !B3
                tmp = 0.0_p_k_fld
                do j = 1, this % n_coef
                    tmp = tmp + this % coef_e(j) * ( e%f3(2, i1+j, i2, i3) - e%f3(2, i1-j+1, i2, i3) )
                enddo
                b%f3(3, i1, i2, i3) = b%f3(3, i1, i2, i3) - dtdx1 * tmp &
                                    + dtdx2 * ( e%f3(1, i1, i2+1, i3) - e%f3(1, i1, i2, i3) )


            enddo
        enddo
    enddo
    !$omp end parallel do

    if ( this%taper_bnd ) then

        ! advance lower boundary
        if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
             bnd_type(p_lower, 1) /= p_bc_other_node ) then

            !$omp parallel do private(i2,i1,ord_hf,tmp,j)
            do i3 = 0, b%nx_(3)+1
                do i2 = 0, b%nx_(2)+1

                    ord_hf = 0

                    do i1 = 0, i1_lb-1

                        ord_hf = ord_hf + 1

                        !B1
                        b%f3(1, i1, i2, i3) = b%f3(1, i1, i2, i3) &
                                            - dtdx2 * ( e%f3(3, i1, i2+1, i3) - e%f3(3, i1, i2, i3)) &
                                            + dtdx3 * ( e%f3(2, i1, i2, i3+1) - e%f3(2, i1, i2, i3))

                        !B2
                        tmp = 0.0_p_k_fld
                        do j = 1, ord_hf
                            tmp = tmp + this % taper_coef_e(j,ord_hf) * ( e%f3(3, i1+j, i2, i3) - e%f3(3, i1-j+1, i2, i3) )
                        enddo
                        b%f3(2, i1, i2, i3) = b%f3( 2, i1, i2, i3 ) + dtdx1 * tmp &
                                            - dtdx3 * ( e%f3(1, i1, i2, i3+1) - e%f3(1, i1, i2, i3))

                        !B3
                        tmp = 0.0_p_k_fld
                        do j = 1, ord_hf
                            tmp = tmp + this % taper_coef_e(j,ord_hf) * ( e%f3(2, i1+j, i2, i3) - e%f3(2, i1-j+1, i2, i3) )
                        enddo
                        b%f3(3, i1, i2, i3) = b%f3(3, i1, i2, i3) - dtdx1 * tmp &
                                            + dtdx2 * ( e%f3(1, i1, i2+1, i3) - e%f3(1, i1, i2, i3) )


                    enddo
                enddo
            enddo
            !$omp end parallel do

        endif

        ! advance upper boundary
        if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
             bnd_type(p_upper, 1) /= p_bc_other_node ) then

            !$omp parallel do private(i2,i1,ord_hf,tmp,j)
            do i3 = 0, b%nx_(3)+1
                do i2 = 0, b%nx_(2)+1

                    ord_hf = offset + 1

                    do i1 = i1_ub+1, b%nx_(1)+1

                        ord_hf = ord_hf - 1

                        !B1
                        b%f3(1, i1, i2, i3) = b%f3(1, i1, i2, i3) &
                                            - dtdx2 * ( e%f3(3, i1, i2+1, i3) - e%f3(3, i1, i2, i3)) &
                                            + dtdx3 * ( e%f3(2, i1, i2, i3+1) - e%f3(2, i1, i2, i3))

                        !B2
                        tmp = 0.0_p_k_fld
                        do j = 1, ord_hf
                            tmp = tmp + this % taper_coef_e(j,ord_hf) * ( e%f3(3, i1+j, i2, i3) - e%f3(3, i1-j+1, i2, i3) )
                        enddo
                        b%f3(2, i1, i2, i3) = b%f3( 2, i1, i2, i3 ) + dtdx1 * tmp &
                                            - dtdx3 * ( e%f3(1, i1, i2, i3+1) - e%f3(1, i1, i2, i3))

                        !B3
                        tmp = 0.0_p_k_fld
                        do j = 1, ord_hf
                            tmp = tmp + this % taper_coef_e(j,ord_hf) * ( e%f3(2, i1+j, i2, i3) - e%f3(2, i1-j+1, i2, i3) )
                        enddo
                        b%f3(3, i1, i2, i3) = b%f3(3, i1, i2, i3) - dtdx1 * tmp &
                                            + dtdx2 * ( e%f3(1, i1, i2+1, i3) - e%f3(1, i1, i2, i3) )


                    enddo
                enddo
            enddo
            !$omp end parallel do

        endif

    endif

end subroutine dbdt_3d_fei

subroutine dedt_3d_fei( this, e, b, jay, dt, bnd_type )

    implicit none
 
    class( t_emf_solver_fei ), intent(in) :: this
    type(t_vdf), intent(inout)                :: e
    type(t_vdf), intent(in)                   :: b, jay
    real(p_double), intent(in)                :: dt
    integer, intent(in), dimension(2,p_max_dim) :: bnd_type

    integer :: i1, i2, i3, j, i1_lb, i1_ub, offset, ord_hf
    real(p_k_fld) :: dtdx1, dtdx2, dtdx3, dtif
    real(p_k_fld) :: tmp

    dtdx1 = real( dt/e%dx(1), p_k_fld )
    dtdx2 = real( dt/e%dx(2), p_k_fld )
    dtdx3 = real( dt/e%dx(3), p_k_fld )
    dtif = real( dt, p_k_fld )

    i1_lb = 2 - this % n_coef
    i1_ub = e%nx_(1) + this % n_coef

    if ( this % taper_bnd ) then

        offset = this%taper_order / 2

        if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
             bnd_type(p_lower, 1) /= p_bc_other_node ) then
            i1_lb = offset + 1
        endif

        if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
            bnd_type(p_upper, 1) /= p_bc_other_node ) then
            i1_ub = e%nx_(1) - offset + 1
        endif

    endif

    !$omp parallel do private(i2,i1,tmp,j)
    do i3 = 1, e%nx_(3)+1
        do i2 = 1, e%nx_(2)+1
            do i1 = i1_lb, i1_ub
                ! E1
                e%f3(1, i1, i2, i3) = e%f3(1, i1, i2, i3) &
                                    - dtif * jay%f3(1, i1, i2, i3) &
                                    + dtdx2 * ( b%f3(3, i1, i2, i3) - b%f3(3, i1, i2-1, i3) ) &
                                    - dtdx3 * ( b%f3(2, i1, i2, i3) - b%f3(2, i1, i2, i3-1) )

                ! E2
                tmp = 0.0_p_k_fld
                do j = 1, this %n_coef
                    tmp = tmp + this % coef_b(j) * ( b%f3(3, i1+j-1, i2, i3) - b%f3(3, i1-j, i2, i3) )
                enddo
                e%f3(2, i1, i2, i3) = e%f3(2, i1, i2, i3) - dtif * jay%f3(2, i1, i2, i3) &
                                + dtdx3 * ( b%f3( 1, i1, i2, i3 ) - b%f3( 1, i1, i2, i3-1 ) ) &
                                - dtdx1 * tmp

                ! E3
                tmp = 0.0_p_k_fld
                do j = 1, this % n_coef
                    tmp = tmp + this % coef_b(j) * ( b%f3(2, i1+j-1, i2, i3) - b%f3(2, i1-j, i2, i3) )
                enddo
                e%f3(3, i1, i2, i3) = e%f3(3, i1, i2, i3) - dtif * jay%f3(3, i1, i2, i3) &
                                - dtdx2 * ( b%f3(1, i1, i2, i3) - b%f3(1, i1, i2-1, i3) ) &
                                + dtdx1 * tmp

            enddo
        enddo
    enddo
    !$omp end parallel do

    if ( this%taper_bnd ) then

        ! advance lower boundary
        if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
             bnd_type(p_lower, 1) /= p_bc_other_node ) then

            !$omp parallel do private(i2,i1,ord_hf,tmp,j)
            do i3 = 1, e%nx_(3)+1
                do i2 = 1, e%nx_(2)+1

                    ord_hf = 0

                    do i1 = 1, i1_ub-1

                        ord_hf = ord_hf + 1
                        ! E1
                        e%f3(1, i1, i2, i3) = e%f3(1, i1, i2, i3) &
                                            - dtif * jay%f3(1, i1, i2, i3) &
                                            + dtdx2 * ( b%f3(3, i1, i2, i3) - b%f3(3, i1, i2-1, i3) ) &
                                            - dtdx3 * ( b%f3(2, i1, i2, i3) - b%f3(2, i1, i2, i3-1) )

                        ! E2
                        tmp = 0.0_p_k_fld
                        do j = 1, ord_hf
                            tmp = tmp + this % taper_coef_b(j,ord_hf) * ( b%f3(3, i1+j-1, i2, i3) - b%f3(3, i1-j, i2, i3) )
                        enddo
                        e%f3(2, i1, i2, i3) = e%f3(2, i1, i2, i3) - dtif * jay%f3(2, i1, i2, i3) &
                                        + dtdx3 * ( b%f3( 1, i1, i2, i3 ) - b%f3( 1, i1, i2, i3-1 ) ) &
                                        - dtdx1 * tmp

                        ! E3
                        tmp = 0.0_p_k_fld
                        do j = 1, ord_hf
                            tmp = tmp + this % taper_coef_b(j,ord_hf) * ( b%f3(2, i1+j-1, i2, i3) - b%f3(2, i1-j, i2, i3) )
                        enddo
                        e%f3(3, i1, i2, i3) = e%f3(3, i1, i2, i3) - dtif * jay%f3(3, i1, i2, i3) &
                                        - dtdx2 * ( b%f3(1, i1, i2, i3) - b%f3(1, i1, i2-1, i3) ) &
                                        + dtdx1 * tmp

                    enddo
                enddo
            enddo
            !$omp end parallel do

        endif

        ! advance upper boundary
        if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
             bnd_type(p_upper, 1) /= p_bc_other_node ) then

            !$omp parallel do private(i2,i1,ord_hf,tmp,j)
            do i3 = 1, e%nx_(3)+1
                do i2 = 1, e%nx_(2)+1

                    ord_hf = offset + 1

                    do i1 = i1_ub+1, e%nx_(1)+1

                        ord_hf = ord_hf - 1
                        ! E1
                        e%f3(1, i1, i2, i3) = e%f3(1, i1, i2, i3) &
                                            - dtif * jay%f3(1, i1, i2, i3) &
                                            + dtdx2 * ( b%f3(3, i1, i2, i3) - b%f3(3, i1, i2-1, i3) ) &
                                            - dtdx3 * ( b%f3(2, i1, i2, i3) - b%f3(2, i1, i2, i3-1) )

                        ! E2
                        tmp = 0.0_p_k_fld
                        do j = 1, ord_hf
                            tmp = tmp + this % taper_coef_b(j,ord_hf) * ( b%f3(3, i1+j-1, i2, i3) - b%f3(3, i1-j, i2, i3) )
                        enddo
                        e%f3(2, i1, i2, i3) = e%f3(2, i1, i2, i3) - dtif * jay%f3(2, i1, i2, i3) &
                                        + dtdx3 * ( b%f3( 1, i1, i2, i3 ) - b%f3( 1, i1, i2, i3-1 ) ) &
                                        - dtdx1 * tmp

                        ! E3
                        tmp = 0.0_p_k_fld
                        do j = 1, ord_hf
                            tmp = tmp + this % taper_coef_b(j,ord_hf) * ( b%f3(2, i1+j-1, i2, i3) - b%f3(2, i1-j, i2, i3) )
                        enddo
                        e%f3(3, i1, i2, i3) = e%f3(3, i1, i2, i3) - dtif * jay%f3(3, i1, i2, i3) &
                                        - dtdx2 * ( b%f3(1, i1, i2, i3) - b%f3(1, i1, i2-1, i3) ) &
                                        + dtdx1 * tmp

                    enddo
                enddo
            enddo
            !$omp end parallel do

        endif

    endif

end subroutine dedt_3d_fei


subroutine dbdt_cyl_2d_fei( this, b, e, dt, bnd_type )

    implicit none

    class( t_emf_solver_fei ), intent(in) :: this
    type( t_vdf ), intent( inout ) :: b
    type( t_vdf ), intent( in )    :: e
    real(p_double),  intent(in)    :: dt
    integer, intent(in), dimension(2,p_max_dim) :: bnd_type

    real(p_k_fld) :: dtdz, dtdr
    real(p_k_fld) :: tmp, rm, rp, tmp_
    integer :: i1, i2, j, i2_0, gshift_i2, i1_lb, i1_ub, offset, ord_hf

    dtdz = real(dt/b%dx_(1), p_k_fld)
    dtdr = real(dt/b%dx_(2), p_k_fld)

    gshift_i2 = this % gir_pos - 2

    i1_lb = - this % n_coef + 1
    i1_ub = b%nx_(1) + this % n_coef

    if ( this%taper_bnd ) then

        offset = this%taper_order / 2

        if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
             bnd_type(p_lower, 1) /= p_bc_other_node ) then
            i1_lb = offset
        endif

        if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
             bnd_type(p_upper, 1) /= p_bc_other_node ) then
            i1_ub = b%nx_(1) - offset + 1
        endif

    endif

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

        do i1 = i1_lb, i1_ub
            !B1
            b%f2( 1, i1, i2 ) = b%f2(1, i1, i2) - tmp * ( rp*e%f2(3, i1, i2+1) - rm*e%f2(3, i1, i2))

            !B2
            tmp_ = 0.0_p_k_fld
            do j = 1, this % n_coef
                tmp_ = tmp_ + this % coef_e(j) * ( e%f2( 3, i1+j, i2 ) - e%f2( 3, i1-j+1, i2 ) )
            enddo
            b%f2( 2, i1, i2 ) = b%f2(2, i1, i2) + dtdz * tmp_

            !B3
            tmp_ = 0.0_p_k_fld
            do j = 1, this % n_coef
                tmp_ = tmp_ + this % coef_e(j) * ( e%f2( 2, i1+j, i2 ) - e%f2( 2, i1-j+1, i2 ))
            enddo
            b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) - dtdz * tmp_ + &
                                        dtdr * ( e%f2( 1, i1, i2+1 ) - e%f2( 1, i1, i2 ))
        enddo
    enddo

    if ( this%taper_bnd ) then

        ! advance lower boundary
        if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
             bnd_type(p_lower, 1) /= p_bc_other_node ) then

            do i2 = i2_0, b%nx_(2)+1

                rp   = i2 + gshift_i2 + 0.5     ! position of the upper edge of the cell (normalized to dr)
                rm   = i2 + gshift_i2 - 0.5     ! position of the lower edge of the cell (normalized to dr)
                tmp  = dtdr/(i2 + gshift_i2)    ! (dt/dr) / position of the middle of the cell (normalized to dr)

                ord_hf = 0

                do i1 = 0, i1_lb-1

                    ord_hf = ord_hf + 1
                    !B1
                    b%f2( 1, i1, i2 ) = b%f2(1, i1, i2) - tmp * ( rp*e%f2(3, i1, i2+1) - rm*e%f2(3, i1, i2))

                    !B2
                    tmp_ = 0.0_p_k_fld
                    do j = 1, ord_hf
                        tmp_ = tmp_ + this % taper_coef_e(j,ord_hf) * ( e%f2( 3, i1+j, i2 ) - e%f2( 3, i1-j+1, i2 ) )
                    enddo
                    b%f2( 2, i1, i2 ) = b%f2(2, i1, i2) + dtdz * tmp_

                    !B3
                    tmp_ = 0.0_p_k_fld
                    do j = 1, ord_hf
                        tmp_ = tmp_ + this % taper_coef_e(j,ord_hf) * ( e%f2( 2, i1+j, i2 ) - e%f2( 2, i1-j+1, i2 ))
                    enddo
                    b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) - dtdz * tmp_ + &
                                                dtdr * ( e%f2( 1, i1, i2+1 ) - e%f2( 1, i1, i2 ))
                enddo
            enddo

        endif

        ! advance upper boundary
        if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
             bnd_type(p_upper, 1) /= p_bc_other_node ) then

            do i2 = i2_0, b%nx_(2)+1

                rp   = i2 + gshift_i2 + 0.5     ! position of the upper edge of the cell (normalized to dr)
                rm   = i2 + gshift_i2 - 0.5     ! position of the lower edge of the cell (normalized to dr)
                tmp  = dtdr/(i2 + gshift_i2)    ! (dt/dr) / position of the middle of the cell (normalized to dr)

                ord_hf = offset + 1

                do i1 = i1_ub+1, b%nx_(1)+1

                    ord_hf = ord_hf - 1
                    !B1
                    b%f2( 1, i1, i2 ) = b%f2(1, i1, i2) - tmp * ( rp*e%f2(3, i1, i2+1) - rm*e%f2(3, i1, i2))

                    !B2
                    tmp_ = 0.0_p_k_fld
                    do j = 1, ord_hf
                        tmp_ = tmp_ + this % taper_coef_e(j,ord_hf) * ( e%f2( 3, i1+j, i2 ) - e%f2( 3, i1-j+1, i2 ) )
                    enddo
                    b%f2( 2, i1, i2 ) = b%f2(2, i1, i2) + dtdz * tmp_

                    !B3
                    tmp_ = 0.0_p_k_fld
                    do j = 1, ord_hf
                        tmp_ = tmp_ + this % taper_coef_e(j,ord_hf) * ( e%f2( 2, i1+j, i2 ) - e%f2( 2, i1-j+1, i2 ))
                    enddo
                    b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) - dtdz * tmp_ + &
                                                dtdr * ( e%f2( 1, i1, i2+1 ) - e%f2( 1, i1, i2 ))
                enddo
            enddo

        endif

    endif

    if ( gshift_i2 < 0 ) then

    ! guard cell 1 (i2 = 0)
    do i1 = i1_lb, i1_ub
        b%f2( 1, i1, 0 ) =   b%f2( 1, i1, 2 )
        b%f2( 2, i1, 0 ) = - b%f2( 2, i1, 3 )
        b%f2( 3, i1, 0 ) = - b%f2( 3, i1, 2 )
    enddo

    ! axial cell (i2 = 1)
    do i1 = i1_lb, i1_ub
        ! rmid = 0 in this cell so the rot E component is calculated differently
        b%f2( 1, i1, 1 ) =   b%f2(1, i1, 1) - 4 * dtdr * e%f2(3, i1, 2)

        ! B2 is not defined on axis so just reflect the corresponding value inside the box
        b%f2( 2, i1, 1 ) = - b%f2( 2, i1, 2 )

        ! B3 is zero on axis
        b%f2( 3, i1, 1 ) = 0.0_p_k_fld
    enddo

    endif

end subroutine dbdt_cyl_2d_fei
    
subroutine dedt_cyl_2d_fei( this, e, b, jay, dt, bnd_type )

    implicit none

    class( t_emf_solver_fei ), intent(in) :: this
    type( t_vdf ), intent(inout) :: e
    type( t_vdf ),    intent(in) :: b, jay
    real(p_double),   intent(in) :: dt
    integer, intent(in), dimension(2,p_max_dim) :: bnd_type

    real(p_k_fld) :: dtdz, dtdr, dtif
    real(p_k_fld) :: tmp, rcp, rcm, tmp_
    integer :: i1, i2, j, i1_lb, i1_ub, gshift_i2, offset, ord_hf

    dtdz = real(dt/b%dx_(1), p_k_fld)
    dtdr = real(dt/b%dx_(2), p_k_fld)
    dtif = real( dt, p_k_fld )

    gshift_i2 = this % gir_pos - 2

    i1_lb = 2 - this % n_coef
    i1_ub = e%nx_(1) + this % n_coef

    if ( this%taper_bnd ) then

        offset = this%taper_order / 2

        if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
             bnd_type(p_lower, 1) /= p_bc_other_node ) then
            i1_lb = offset + 1
        endif

        if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
             bnd_type(p_upper, 1) /= p_bc_other_node ) then
            i1_ub = e%nx_(1) - offset + 1
        endif

    endif

    do i2 = 1, e%nx_(2)+1

        tmp = dtdr / (( i2 + gshift_i2 ) - 0.5)
        rcp = ( i2 + gshift_i2 )
        rcm = ( i2 + gshift_i2 - 1 )

        do i1 = i1_lb, i1_ub
            ! E1
            e%f2(1, i1, i2) = e%f2(1, i1, i2) &
                            - dtif * jay%f2(1, i1, i2) &
                            + tmp * ( rcp * b%f2(3, i1, i2) - rcm * b%f2(3, i1, i2-1) )

            ! E2
            tmp_ = 0.0_p_k_fld
            do j = 1, this % n_coef
                tmp_ = tmp_ + this % coef_b(j) * ( b%f2(3, i1+j-1, i2) - b%f2(3, i1-j, i2) )
            enddo
            e%f2(2, i1, i2) = e%f2(2, i1, i2) &
                            - dtif * jay%f2(2, i1, i2) - dtdz * tmp_

            ! E3
            tmp_ = 0.0_p_k_fld
            do j = 1, this % n_coef
                tmp_ = tmp_ + this % coef_b(j) * ( b%f2(2, i1+j-1, i2) - b%f2(2, i1-j, i2))
            enddo
            e%f2(3, i1, i2) = e%f2(3, i1, i2) &
                            - dtif * jay%f2(3, i1, i2) + dtdz * tmp_ &
                            - dtdr * ( b%f2(1, i1, i2) - b%f2(1, i1, i2-1))
        enddo
    enddo

    if ( this%taper_bnd ) then

        ! advance lower boundary
        if ( bnd_type(p_lower, 1) /= p_bc_periodic .and. &
             bnd_type(p_lower, 1) /= p_bc_other_node ) then

        do i2 = 1, e%nx_(2)+1

            tmp = dtdr / (( i2 + gshift_i2 ) - 0.5)
            rcp = ( i2 + gshift_i2 )
            rcm = ( i2 + gshift_i2 - 1 )

            ord_hf = 0

            do i1 = 1, i1_lb-1

                ord_hf = ord_hf + 1
                ! E1
                e%f2(1, i1, i2) = e%f2(1, i1, i2) &
                                - dtif * jay%f2(1, i1, i2) &
                                + tmp * ( rcp * b%f2(3, i1, i2) - rcm * b%f2(3, i1, i2-1) )

                ! E2
                tmp_ = 0.0_p_k_fld
                do j = 1, ord_hf
                    tmp_ = tmp_ + this % taper_coef_b(j,ord_hf) * ( b%f2(3, i1+j-1, i2) - b%f2(3, i1-j, i2) )
                enddo
                e%f2(2, i1, i2) = e%f2(2, i1, i2) &
                                - dtif * jay%f2(2, i1, i2) - dtdz * tmp_

                ! E3
                tmp_ = 0.0_p_k_fld
                do j = 1, ord_hf
                    tmp_ = tmp_ + this % taper_coef_b(j,ord_hf) * ( b%f2(2, i1+j-1, i2) - b%f2(2, i1-j, i2))
                enddo
                e%f2(3, i1, i2) = e%f2(3, i1, i2) &
                                - dtif * jay%f2(3, i1, i2) + dtdz * tmp_ &
                                - dtdr * ( b%f2(1, i1, i2) - b%f2(1, i1, i2-1))
            enddo
        enddo

        endif

        ! advance upper boundary
        if ( bnd_type(p_upper, 1) /= p_bc_periodic .and. &
             bnd_type(p_upper, 1) /= p_bc_other_node ) then

        do i2 = 1, e%nx_(2)+1

            tmp = dtdr / (( i2 + gshift_i2 ) - 0.5)
            rcp = ( i2 + gshift_i2 )
            rcm = ( i2 + gshift_i2 - 1 )

            ord_hf = offset + 1

            do i1 = i1_ub+1, e%nx_(1)+1

                ord_hf = ord_hf - 1
                ! E1
                e%f2(1, i1, i2) = e%f2(1, i1, i2) &
                                - dtif * jay%f2(1, i1, i2) &
                                + tmp * ( rcp * b%f2(3, i1, i2) - rcm * b%f2(3, i1, i2-1) )

                ! E2
                tmp_ = 0.0_p_k_fld
                do j = 1, ord_hf
                    tmp_ = tmp_ + this % taper_coef_b(j,ord_hf) * ( b%f2(3, i1+j-1, i2) - b%f2(3, i1-j, i2) )
                enddo
                e%f2(2, i1, i2) = e%f2(2, i1, i2) &
                                - dtif * jay%f2(2, i1, i2) - dtdz * tmp_

                ! E3
                tmp_ = 0.0_p_k_fld
                do j = 1, ord_hf
                    tmp_ = tmp_ + this % taper_coef_b(j,ord_hf) * ( b%f2(2, i1+j-1, i2) - b%f2(2, i1-j, i2))
                enddo
                e%f2(3, i1, i2) = e%f2(3, i1, i2) &
                                - dtif * jay%f2(3, i1, i2) + dtdz * tmp_ &
                                - dtdr * ( b%f2(1, i1, i2) - b%f2(1, i1, i2-1))
            enddo
        enddo

        endif

    endif

    ! Correct values on axial boundary.
    if ( gshift_i2 < 0 ) then
        ! axial cell (i2 = 1)
        do i1 = i1_lb, i1_ub
            e%f2(1, i1, 1) =   e%f2(1, i1, 2)
            e%f2(2, i1, 1) =   0
            e%f2(3, i1, 1) = - e%f2(3, i1, 2)
        enddo
    endif

end subroutine dedt_cyl_2d_fei

!-------------------------------------------------------------------------------
!       generate solver coefficients for standard FDTD algorithm
!-------------------------------------------------------------------------------
subroutine std_coef( ord, coef )

implicit none

    integer, intent(in) :: ord
    real(p_k_fld), dimension(:), intent(out) :: coef

    integer :: n
    integer :: k

    n = ord / 2

    do k = 1, n
    coef(k) = real( &
        (-1.0_p_double)**(k+1) * 16.0_p_double**(1-n) * factorial(2*n-1)**2 &
        / ( (2*k-1)**2 * factorial(n+k-1) * factorial(n-k) * factorial(n-1)**2 ), &
        p_k_fld )
    enddo

end subroutine std_coef
    
!-------------------------------------------------------------------------------
!       generate solver coefficients for bumped FDTD algorithm
!-------------------------------------------------------------------------------
subroutine bump_coef( ord, n_coef, kl, ku, dk, coef )

    implicit none

    integer, intent(in) :: ord, n_coef
    real(p_double), intent(in) :: kl, ku, dk
    real(p_k_fld), dimension(:), intent(out) :: coef

    integer :: i, j, k
    real(p_double), dimension(n_coef) :: a, c_p
    real(p_k_fld),  dimension(n_coef) :: c_p_
    real(p_double), dimension(ord/2, n_coef) :: m
    real(p_double), dimension(n_coef+ord/2, n_coef+ord/2) :: m_left, m_left_inv
    real(p_double), dimension(n_coef+ord/2) :: m_right, sol


    do i = 1, n_coef
        k = 2*i - 1
        a(i) = 8.0_p_double * dk * ( cos(k*pi*ku) - cos(k*pi*kl) ) &
                / ( k * (k**2 * (ku-kl)**2 - 4.0_p_double) )
    enddo

    do j = 1, n_coef
        do i = 1, ord/2
            m(i,j) = real(2*j-1, p_double)**(2*i-1) / factorial(2*i-1)
        enddo
    enddo

    m_left = 0.0_p_double
    m_right = 0.0_p_double; m_right(n_coef+1) = 1.0_p_double
    c_p = 0.0_p_double; c_p_ = 0.0_p_k_fld
    call std_coef( ord, c_p_ )
    c_p = real( c_p_, p_double )

    do i = 1, n_coef
        m_left(i,i) = 0.5_p_double * pi**(-2)
        m_right(i) = 0.5_p_double * pi**(-2) * ( a(i) + c_p(i) )
    enddo
    m_left(1:n_coef, n_coef+1:n_coef+ord/2) = transpose(m)
    m_left(n_coef+1:n_coef+ord/2, 1:n_coef) = m

    call precond( m_left, m_right, n_coef + ord/2 )

    call inverse( m_left, m_left_inv, n_coef + ord/2 )

    sol = 0.0_p_double
    do i = 1, n_coef + ord/2
        do j = 1, n_coef + ord/2
            sol(i) = sol(i) + m_left_inv(i,j) * m_right(j)
        enddo
    enddo

    coef(1:n_coef) = real( sol(1:n_coef), p_k_fld )

end subroutine bump_coef

!-------------------------------------------------------------------------------
!       generate solver coefficients for arbitrary [k] operator
!-------------------------------------------------------------------------------
subroutine gen_solver_coef( fun, arglst, ord, n_coef, coef, w, n )

    implicit none

    procedure(kernel), intent(in), pointer :: fun
    real(p_double), dimension(:), intent(in) :: arglst
    integer, intent(in) :: ord, n_coef, n
    real(p_double), intent(in) :: w
    real(p_k_fld), dimension(:), intent(out) :: coef

    integer :: i, j, iter
    real(p_double), dimension(n_coef+ord/2, n_coef+ord/2) :: m_left, m_left_inv
    real(p_double), dimension(n_coef+ord/2) :: m_right, sol

    real(p_double), dimension(p_max_lst_len) :: lst
    real(p_double), save :: lower = 0.0_p_double, upper = 0.5_p_double, tol = 1.0d-12
    procedure(kernel), pointer :: fun_ptr
    logical :: print_err_integral

    m_left = 0.0_p_double
    lst    = 0.0_p_double
    print_err_integral = .true.
    fun_ptr => kernel_sine
    do j = 1, n_coef
        do i = 1, n_coef
            lst(1:2) = (/ real(i, p_double), real(j, p_double) /)
            iter = 1
            m_left(i,j) = (2.0_p_double/pi**2) * integral( fun_ptr, lst, upper, lower, tol, w, n, iter )
            if ( iter == p_err_integral .and. print_err_integral .and. mpi_node() == 0 ) then
                write(0,*) "   (*error*) in gen_solver_coef for Customized (Fei) EMF solver"
                write(0,*) "   Integral did not properly converge, possible Courant limit violation"
                print_err_integral = .false.
            endif
        enddo
    enddo

    do j = 1, n_coef
        do i = 1, ord/2
            m_left(i+n_coef,j) = real(2*j-1, p_double)**(2*i-1) / factorial(2*i-1)
            m_left(j,i+n_coef) = m_left(i+n_coef,j)
        enddo
    enddo

    m_right           = 0.0_p_double
    m_right(n_coef+1) = 1.0_p_double
    lst               = 0.0_p_double
    do i = 1, n_coef
        lst(1) = real(i, p_double)
        lst( 2:size(arglst)+1 ) = arglst
        iter = 1
        m_right(i) = (2.0_p_double/pi) * integral( fun, lst, upper, lower, tol, w, n, iter )
        if ( iter == p_err_integral .and. print_err_integral .and. mpi_node() == 0 ) then
            write(0,*) "   (*error*) in gen_solver_coef for Customized (Fei) EMF solver"
            write(0,*) "   Integral did not properly converge, possible Courant limit violation"
            print_err_integral = .false.
        endif
    enddo

    call precond( m_left, m_right, n_coef + ord/2 )

    call inverse( m_left, m_left_inv, n_coef + ord/2 )

    sol = 0.0_p_double
    do i = 1, n_coef + ord/2
        do j = 1, n_coef + ord/2
            sol(i) = sol(i) + m_left_inv(i,j) * m_right(j)
        enddo
    enddo

    coef(1:n_coef) = real( sol(1:n_coef), p_k_fld )

end subroutine gen_solver_coef

function kernel_sine( x, arglst )

    implicit none

    real(p_double), intent(in) :: x
    real(p_double), intent(in), dimension(:) :: arglst
    real(p_double) :: kernel_sine

    real(p_double) :: u, v

    u = ( 2.0_p_double * arglst(1) - 1.0_p_double ) * pi
    v = ( 2.0_p_double * arglst(2) - 1.0_p_double ) * pi

    kernel_sine = sin(u*x) * sin(v*x)

end function kernel_sine

!-------------------------------------------------------------------------------
! Integral kernel for bump-type solver
! Not in use currently
!-------------------------------------------------------------------------------
! function kernel_bump( x, arglst )

!     implicit none

!     real(p_double), intent(in) :: x
!     real(p_double), intent(in), dimension(:) :: arglst
!     real(p_double) :: kernel_bump

!     real(p_double) :: u, kl, ku, dk, mid
!     real(p_double), save :: log2 = 0.693147180559945_p_double

!     u = ( 2.0_p_double * arglst(1) - 1.0_p_double ) * pi
!     kl = arglst(2)
!     ku = arglst(3)
!     dk = arglst(4)
!     mid = 0.5_p_double * ( kl + ku )

!     kernel_bump = x + ( dk * sin( pi * (x-kl) / (ku-kl) )**2 ) * &
!         exp( -log2 * (2.0_p_double * (x-mid) / (ku-kl))**20 )
!     kernel_bump = kernel_bump * sin(u*x)

! end function kernel_bump

!-------------------------------------------------------------------------------
! Integral kernel for Xu-type solver
!-------------------------------------------------------------------------------
function kernel_xu( x, arglst )

    implicit none

    real(p_double), intent(in) :: x
    real(p_double), intent(in), dimension(:) :: arglst
    real(p_double) :: kernel_xu

    real(p_double) :: u, v

    u = ( 2.0_p_double * arglst(1) - 1.0_p_double ) * pi
    v = arglst(2) * pi

    kernel_xu = sin(u*x) * sin(v*x) / v

end function kernel_xu

!-------------------------------------------------------------------------------
! Integral kernels for dual type solver
!-------------------------------------------------------------------------------
function kernel_dual_e( x, arglst )

    implicit none

    real(p_double), intent(in) :: x
    real(p_double), intent(in), dimension(:) :: arglst
    real(p_double) :: kernel_dual_e

    real(p_double) :: u, v, omega, k_cons

    u = ( 2.0_p_double * arglst(1) - 1.0_p_double ) * pi
    v = arglst(2) * pi

    k_cons = sin(v*x) / v
    omega = asin(v*k_cons) / v
    kernel_dual_e = sin(u*x) * k_cons / cos(v*omega)

end function kernel_dual_e

function kernel_dual_b( x, arglst )

    implicit none

    real(p_double), intent(in) :: x
    real(p_double), intent(in), dimension(:) :: arglst
    real(p_double) :: kernel_dual_b

    real(p_double) :: u, v, omega, k_cons

    u = ( 2.0_p_double * arglst(1) - 1.0_p_double ) * pi
    v = arglst(2) * pi

    k_cons = sin(v*x) / v
    omega = asin(v*k_cons) / v
    kernel_dual_b = sin(u*x) * k_cons * cos(v*omega)

end function kernel_dual_b

!-------------------------------------------------------------------------------
! Test algorithm stability for resolution and time-step
!-------------------------------------------------------------------------------
subroutine test_stability_fei( this, dt, dx )

    implicit none

    class( t_emf_solver_fei ), intent(in) :: this
    real(p_double), intent(in)  ::  dt
    real(p_double), dimension(:), intent(in)  ::  dx

    real(p_double) :: cour_new, cour_old, tmp
    real(p_k_fld), dimension(p_max_n_coef) :: coef_e, coef_b
    procedure(kernel), pointer :: kernel_ptr
    integer :: i, j, cnt
    integer, parameter :: cnt_max = 100

    SCR_ROOT( 'Calculating Courant stability condition...' )

    cnt = 0
    cour_new = 0.0_p_double
    cour_old = 0.0_p_double
    tmp = 0.0_p_double

    coef_e(1:this%n_coef) = this%coef_e(1:this%n_coef)
    coef_b(1:this%n_coef) = this%coef_b(1:this%n_coef)

    select case ( this%type )

    case ( p_emf_fei_std, p_emf_fei_bump, p_emf_fei_coef )

        do j = 1, this%n_coef
            cour_new = cour_new + abs( coef_e(j) )
            tmp  = tmp  + abs( coef_b(j) )
        enddo

        cour_new = cour_new * tmp / (dx(1))**2

        do i = 2, p_x_dim
            cour_new = cour_new + 1.0_p_double / (dx(i))**2
        enddo

        cour_new = sqrt( 1.0_p_double / cour_new )

    case ( p_emf_fei_xu )

        ! set the CFL of Yee as the initial guess
        do i = 1, p_x_dim
            cour_old = cour_old + 1.0_p_double / (dx(i))**2
        enddo
        cour_old = sqrt( 1.0_p_double / cour_old )

        ! iterate to get CFL
        kernel_ptr => kernel_xu
        cnt = 0
        do

            ! calculate the solver stencil coefficients
            call gen_solver_coef( kernel_ptr, (/ cour_old/dx(1) /), this%solver_ord, &
                                    this%n_coef, coef_e, this%weight_w, this%weight_n )
            coef_b(1:this%n_coef) = coef_e(1:this%n_coef)

            ! calculate the updated CFL
            tmp = 0.0_p_double
            cour_new = 0.0_p_double
            do j = 1, this%n_coef
                cour_new = cour_new + abs( coef_e(j) )
                tmp  = tmp  + abs( coef_b(j) )
            enddo

            cour_new = cour_new * tmp / (dx(1))**2

            do i = 2, p_x_dim
                cour_new = cour_new + 1.0_p_double / (dx(i))**2
            enddo

            cour_new = sqrt( 1.0_p_double / cour_new )
            cnt = cnt + 1

            if ( abs( cour_new - cour_old ) / cour_old < 1.0d-6 ) then
                exit
            elseif ( cnt > cnt_max ) then
                SCR_ROOT( '(*error*) Calculation of (Fei) EMF solver stability condition may not converge, aborting' )
                call abort_program(p_err_invalid)
            else
                cour_old = cour_new
                cour_new = 0.0_p_double
            endif

        enddo

    case ( p_emf_fei_dual )

        ! set the CFL of Yee as the initial guess
        do i = 1, p_x_dim
            cour_old = cour_old + 1.0_p_double / (dx(i))**2
        enddo
        cour_old = sqrt( 1.0_p_double / cour_old )

        ! iterate to get CFL
        cnt = 0
        do

            ! calculate the solver stencil coefficients
            kernel_ptr => kernel_dual_e
            call gen_solver_coef( kernel_ptr, (/ cour_old/dx(1) /), this%solver_ord, &
                                    this%n_coef, coef_e, this%weight_w, this%weight_n )
            kernel_ptr => kernel_dual_b
            call gen_solver_coef( kernel_ptr, (/ cour_old/dx(1) /), this%solver_ord, &
                                    this%n_coef, coef_b, this%weight_w, this%weight_n )

            ! calculate the updated CFL
            tmp = 0.0_p_double
            cour_new = 0.0_p_double
            do j = 1, this%n_coef
                cour_new = cour_new + abs( coef_e(j) )
                tmp = tmp  + abs( coef_b(j) )
            enddo

            cour_new = cour_new * tmp / (dx(1))**2

            do i = 2, p_x_dim
                cour_new = cour_new + 1.0_p_double / (dx(i))**2
            enddo

            cour_new = sqrt( 1.0_p_double / cour_new )
            cnt = cnt + 1

            if ( abs( cour_new - cour_old ) / cour_old < 1.0d-6 ) then
                exit
            elseif ( cnt > cnt_max ) then
                SCR_ROOT( '(*error*) Calculation of (Fei) EMF solver stability condition may not converge, aborting' )
                call abort_program(p_err_invalid)
            else
                cour_old = cour_new
                cour_new = 0.0_p_double
            endif

        enddo

    end select

    if (dt > cour_new) then
        if ( mpi_node() == 0 ) then
            print *, '(*error*) Customized (Fei) EMF solver stability condition violated, aborting'
            print *, 'Please reduce the time step and re-try.'
            print *, 'Note that CFL will change as the time step changes.'
            print *, 'dx   = ', dx(1:p_x_dim)
            print *, 'dt   = ', dt, ' dt(CFL) = ', cour_new
        endif
        call abort_program(p_err_invalid)
    endif

end subroutine test_stability_fei

!-------------------------------------------------------------------------------
! Cleanup solver buffers if necessary
!-------------------------------------------------------------------------------
subroutine cleanup_fei( this )

    implicit none

    class( t_emf_solver_fei ), intent(inout) :: this

    if ( associated(this%f_buffer2) ) call freemem( this%f_buffer2 )
    if ( associated(this%f_buffer3) ) call freemem( this%f_buffer3 )

    if ( allocated(this%coef_e) ) then
        deallocate( this%coef_e, this%coef_b )
    endif
    if ( allocated(this%taper_coef_e) ) then
        deallocate( this%taper_coef_e, this%taper_coef_b )
    endif

    call fft_cleanup()

end subroutine cleanup_fei
    
    
function factorial(n)

    implicit none

    integer, intent(in) :: n
    real(p_double) :: factorial

    integer :: j

    factorial = 1.0_p_double
    do j = 2,n
    factorial = factorial * real(j, p_double)
    enddo

end function factorial
    
    
!-------------------------------------------------------------------------------
! Inverse matrix based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-------------------------------------------------------------------------------
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! c(n,n) - inverse matrix of A
!-------------------------------------------------------------------------------
subroutine inverse( a, c, n )

    implicit none

    integer, intent(in) :: n
    real(p_double), intent(in), dimension(n,n)  :: a
    real(p_double), intent(out), dimension(n,n) :: c

    real(p_double), dimension(n,n)              :: L, U, a_tmp
    real(p_double), dimension(n)                :: b, d, x
    real(p_double)                              :: coeff
    integer                                     :: i, j, k

    ! step 0: initialization for matrices L and U and b
    L = 0.0_p_double
    U = 0.0_p_double
    b = 0.0_p_double

    ! duplicate matrix a due to the calculation will destroy the input matrix
    a_tmp = a

    ! step 1: forward elimination
    do k = 1, n-1
        do i = k+1, n
            coeff = a_tmp(i,k)/a_tmp(k,k)
            L(i,k) = coeff
            do j = k+1, n
                a_tmp(i,j) = a_tmp(i,j) - coeff*a_tmp(k,j)
            enddo
        enddo
    enddo

    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i = 1, n
        L(i,i) = 1.0_p_double
    enddo
    ! U matrix is the upper triangular part of A
    do j = 1, n
        do i = 1, j
            U(i,j) = a_tmp(i,j)
        enddo
    enddo

    ! Step 3: compute columns of the inverse matrix C
    do k = 1, n
        b(k) = 1.0_p_double
        d(1) = b(1)
        ! Step 3a: Solve Ld=b using the forward substitution
        do i = 2,n
            d(i) = b(i)
            do j = 1, i-1
            d(i) = d(i) - L(i,j)*d(j)
            enddo
        enddo
        ! Step 3b: Solve Ux=d using the back substitution
        x(n) = d(n)/U(n,n)
        do i = n-1, 1, -1
            x(i) = d(i)
            do j = n, i+1, -1
            x(i) = x(i) - U(i,j)*x(j)
            enddo
            x(i) = x(i)/u(i,i)
        enddo
        ! Step 3c: fill the solutions x(n) into column k of C
        do i = 1, n
            c(i,k) = x(i)
        enddo
        b(k) = 0.0_p_double
    enddo

end subroutine inverse
    
!-------------------------------------------------------------------------------
! precondition the linear equations Ax=b of dimension n to improve
! the condition number
!-------------------------------------------------------------------------------
subroutine precond( A, b, n )

    implicit none

    real(p_double), intent(inout), dimension(:,:) :: A
    real(p_double), intent(inout), dimension(:) :: b
    integer, intent(in) :: n

    real(p_double) :: d = 0.0_p_double
    integer :: i

    do i = 1, n
    d = maxval( abs( A(i,:) ) )**(-1)
    b(i) = b(i) * d
    A(i,:) = A(i,:) * d
    enddo

end subroutine precond

!-------------------------------------------------------------------------------
! numerical integral using adaptive quadrature and trapezoidal rule
!-------------------------------------------------------------------------------
recursive function integral( fun, arglst, lower, upper, abs_tol, w, n, iter ) result(res)

    implicit none

    procedure(kernel), intent(in), pointer :: fun
    real(p_double), intent(in), dimension(:) :: arglst
    real(p_double), intent(in) :: lower, upper, abs_tol, w
    integer, intent(in) :: n
    integer, intent(inout) :: iter
    real(p_double) :: res

    integer, parameter :: num_pts = 128, iter_max = 10
    integer :: i, iter1, iter2
    real(p_double) :: xi, dx_fine, dx_crs, int_fine, int_crs

    if ( iter > iter_max ) then
        res = 0.0_p_double
        iter = p_err_integral
    else

        dx_crs   = ( upper - lower ) / num_pts
        dx_fine  = 0.5_p_double * dx_crs
        int_fine = 0.0_p_double
        int_crs  = 0.0_p_double

        int_crs  = 0.5 * dx_crs * ( &
            fun( lower, arglst ) * weight( lower, w, n ) + &
            fun( upper, arglst ) * weight( upper, w, n ) )
        int_fine = 0.5 * int_crs
        do i = 1, num_pts-1
            xi       = lower + i * dx_crs
            int_crs  = int_crs + dx_crs * fun( xi, arglst ) * weight( xi, w, n )
            int_fine = int_fine + dx_fine * ( &
                fun( xi, arglst ) * weight( xi, w, n ) + &
                fun( xi-dx_fine, arglst ) * weight( xi-dx_fine, w, n ) )
        enddo
        int_fine = int_fine + dx_fine * fun( xi+dx_fine, arglst ) * weight( xi+dx_fine, w, n )

        if ( abs( int_crs - int_fine ) > abs_tol ) then
            iter1 = iter + 1; iter2 = iter + 1
            res = integral( fun, arglst, lower, 0.5_p_double*(lower+upper), 0.5_p_double*abs_tol, w, n, iter1 ) + &
                  integral( fun, arglst, 0.5_p_double*(lower+upper), upper, 0.5_p_double*abs_tol, w, n, iter2 )
            if ( iter1 == p_err_integral .or. iter2 == p_err_integral ) iter = p_err_integral
        else
            res = int_fine
        endif

    endif

end function integral

function weight( x, w, n )

    implicit none
    real(p_double), intent(in) :: x, w
    integer, intent(in) :: n
    real(p_double) :: weight

    real(p_double), save :: log2 = 0.693147180559945_p_double

    weight = exp( -log2*(x/w)**n )

end function weight

!-------------------------------------------------------------------------------
! Filter / correct current 
!-------------------------------------------------------------------------------
subroutine filter_correct_current_2d( this, jay )
  
    implicit none
  
    class( t_emf_solver_fei ), intent(inout) :: this
    type( t_vdf ), intent(inout) :: jay
  
    real(p_k_fld) :: filter_limit, filter_width
    real(p_k_fld) :: corr1, corr2
    real(p_k_fld) :: dm1, k1, k1_half
    integer :: gc_l1, gc_u1, gc_l2, gc_u2, n1, n2, m1, m2, dc, fft_m1
    integer, dimension(2,2) :: gc_num
    integer, dimension(3) :: bshape
    integer :: i, i1, i2, j
  
    filter_width = this%filter_width
    filter_limit = this%filter_limit
    gc_num = jay%gc_num()
    gc_l1 = gc_num(p_lower,1)
    gc_u1 = gc_num(p_upper,1)
    gc_l2 = gc_num(p_lower,2)
    gc_u2 = gc_num(p_upper,2)
    n1 = jay%nx_(1)
    n2 = jay%nx_(2)
    m1 = n1 + gc_l1 + gc_u1
    m2 = n2 + gc_l2 + gc_u2
    dm1 = 2.0_p_k_fld * pi / real(m1, p_k_fld)
    dc = this%n_damp_cell
  
    fft_m1 = 2 * ( m1/2 + 1 )
  
    if ( .not. associated( this%f_buffer2 ) ) then
      SCR_ROOT( ' - using FFTW subroutine' )
      call alloc( this%f_buffer2, (/ 3, fft_m1, m2 /) )
      call fft_init( m1 )
    else
      ! check if buffer should be resized (e.g., from dlb)
      bshape = shape(this%f_buffer2)
      if ( bshape(2)/=fft_m1 .or. bshape(3)/=m2 ) then
        call freemem( this%f_buffer2 )
        call alloc( this%f_buffer2, (/ 3, fft_m1, m2 /) )
        ! re-init fft if necessary
        if ( bshape(2)/=fft_m1 ) then
          call fft_cleanup()
          call fft_init( m1 )
        endif
      endif
    endif
  
    this%f_buffer2 = 0.0_p_k_fld
  
    ! damp the current
    do i2 = 1, m2
      do i1 = 1, m1
        do i = 1, 3
  
            if ( i1 <= dc ) then
              corr1 = cos( 0.5_p_k_fld * pi * (i1-dc-1) / dc )**2
            elseif ( i1 > m1-dc ) then
              corr1 = cos( 0.5_p_k_fld * pi * (i1-m1+dc) / dc )**2
            else
              corr1 = 1.0_p_k_fld
            endif
            this%f_buffer2(i,i1,i2) = jay%f2( i, i1-gc_l1, i2-gc_l2 ) * corr1
  
        enddo
      enddo
    enddo
  
    call fft( this%f_buffer2, m1, m2 )
  
    ! filter and corret current
    do i2 = 1, m2
      do i1 = 2, fft_m1/2
        k1 = dm1 * real(i1-1, p_k_fld)
        k1_half = k1 * 0.5_p_k_fld
  
        ! current filter factor
        corr1 = 1.0_p_k_fld
        if ( ( k1 > filter_limit*pi ) .and. &
              ( this%filter_current ) ) then
  
          if ( k1 < ( filter_limit+filter_width) * pi ) then
            corr1 = cos( 0.5_p_k_fld * pi * (k1 - filter_limit*pi) / (filter_width*pi) )**2
          else
            corr1 = 0.0_p_k_fld
          endif
  
        endif
  
        ! current corretion factor
        if ( this%correct_current ) then
  
          corr2 = 0.0_p_k_fld
          do j = 1, this%n_coef
            corr2 = corr2 + this%coef_b(j) * sin( real(2*j-1, p_k_fld) * k1_half )
          enddo
          corr2 = sin(k1_half) / corr2
  
        else
  
          corr2 = 1.0_p_k_fld
  
        endif
  
        this%f_buffer2(1,2*i1-1,i2) = this%f_buffer2(1,2*i1-1,i2) * corr1 * corr2
        this%f_buffer2(2,2*i1-1,i2) = this%f_buffer2(2,2*i1-1,i2) * corr1 * corr2
        this%f_buffer2(3,2*i1-1,i2) = this%f_buffer2(3,2*i1-1,i2) * corr1
        this%f_buffer2(1,2*i1  ,i2) = this%f_buffer2(1,2*i1  ,i2) * corr1
        this%f_buffer2(2,2*i1  ,i2) = this%f_buffer2(2,2*i1  ,i2) * corr1
        this%f_buffer2(3,2*i1  ,i2) = this%f_buffer2(3,2*i1  ,i2) * corr1
  
      enddo
    enddo
  
  
    call ifft( this%f_buffer2, m1, m2 )
  
    ! direct copy
    do i2 = 1, m2
      do i1 = 1,m1
        do i = 1, 3
          jay%f2( i, i1-gc_l1, i2-gc_l2 ) = this%f_buffer2(i,i1,i2)
        enddo
      enddo
    enddo
  
end subroutine filter_correct_current_2d
  
subroutine filter_correct_current_3d( this, jay )
  
    implicit none
  
    class( t_emf_solver_fei ), intent(inout) :: this
    type( t_vdf ), intent(inout) :: jay
  
    real(p_k_fld) :: filter_limit, filter_width
    real(p_k_fld) :: corr1, corr2
    real(p_k_fld) :: dm1, k1, k1_half
    integer :: gc_l1, gc_u1, gc_l2, gc_u2, gc_l3, gc_u3, n1, n2, n3, m1, m2, m3, dc, fft_m1
    integer, dimension(2,3) :: gc_num
    integer, dimension(4) :: bshape
    integer :: i, i1, i2, i3, j
  
    filter_width = this%filter_width
    filter_limit = this%filter_limit
    gc_num = jay%gc_num()
    gc_l1 = gc_num(p_lower,1)
    gc_u1 = gc_num(p_upper,1)
    gc_l2 = gc_num(p_lower,2)
    gc_u2 = gc_num(p_upper,2)
    gc_l3 = gc_num(p_lower,3)
    gc_u3 = gc_num(p_upper,3)
    n1 = jay%nx_(1)
    n2 = jay%nx_(2)
    n3 = jay%nx_(3)
    m1 = n1 + gc_l1 + gc_u1
    m2 = n2 + gc_l2 + gc_u2
    m3 = n3 + gc_l3 + gc_u3
    dm1 = 2.0_p_k_fld * pi / real(m1, p_k_fld)
    dc = this%n_damp_cell

    fft_m1 = 2 * ( m1/2 + 1 )
  
    if ( .not. associated( this%f_buffer3 ) ) then
      SCR_ROOT( ' - using FFTW subroutine' )
      call alloc( this%f_buffer3, (/ 3, fft_m1, m2, m3 /) )
      call fft_init( m1 )
    else
      ! check if buffer should be resized (e.g., from dlb)
      bshape = shape(this%f_buffer3)
      if ( bshape(2)/=fft_m1 .or. bshape(3)/=m2 .or. bshape(4)/=m3 ) then
        call freemem( this%f_buffer3 )
        call alloc( this%f_buffer3, (/ 3, fft_m1, m2, m3 /) )
        ! re-init fft if necessary
        if ( bshape(2)/=fft_m1 ) then
          call fft_cleanup()
          call fft_init( m1 )
        endif
      endif
    endif
  
    this%f_buffer3 = 0.0_p_k_fld
  
    ! damp the current
    do i3 = 1, m3
      do i2 = 1, m2
        do i1 = 1, m1
          do i = 1, 3
  
              if ( i1 <= dc ) then
                corr1 = cos( 0.5_p_k_fld * pi * (i1-dc-1) / dc )**2
              elseif ( i1 > m1-dc ) then
                corr1 = cos( 0.5_p_k_fld * pi * (i1-m1+dc) / dc )**2
              else
                corr1 = 1.0_p_k_fld
              endif
              this%f_buffer3(i,i1,i2,i3) = jay%f3( i, i1-gc_l1, i2-gc_l2, i3-gc_l3 ) * corr1
  
          enddo
        enddo
      enddo
    enddo
  
    call fft( this%f_buffer3, m1, m2, m3 )
  
    ! filter and corret current
    do i3 = 1, m3
      do i2 = 1, m2
        do i1 = 2, fft_m1/2
          k1 = dm1 * real(i1-1, p_k_fld)
          k1_half = k1 * 0.5_p_k_fld
  
          ! current filter factor
          corr1 = 1.0_p_k_fld
          if ( ( k1 > filter_limit*pi ) .and. ( this%filter_current ) ) then
  
            if ( k1 < ( filter_limit+filter_width) * pi ) then
              corr1 = cos( 0.5_p_k_fld * pi * (k1 - filter_limit*pi) / (filter_width*pi) )**2
            else
              corr1 = 0.0_p_k_fld
            endif
  
          endif
  
          ! current corretion factor
          if ( this%correct_current ) then
  
            corr2 = 0.0_p_k_fld
            do j = 1, this%n_coef
              corr2 = corr2 + this%coef_b(j) * sin( real(2*j-1, p_k_fld) * k1_half )
            enddo
            corr2 = sin(k1_half) / corr2
  
          else
  
            corr2 = 1.0_p_k_fld
  
          endif
  
          this%f_buffer3(1,2*i1-1,i2,i3) = this%f_buffer3(1,2*i1-1,i2,i3) * corr1 * corr2
          this%f_buffer3(2,2*i1-1,i2,i3) = this%f_buffer3(2,2*i1-1,i2,i3) * corr1 * corr2
          this%f_buffer3(3,2*i1-1,i2,i3) = this%f_buffer3(3,2*i1-1,i2,i3) * corr1
          this%f_buffer3(1,2*i1  ,i2,i3) = this%f_buffer3(1,2*i1  ,i2,i3) * corr1
          this%f_buffer3(2,2*i1  ,i2,i3) = this%f_buffer3(2,2*i1  ,i2,i3) * corr1
          this%f_buffer3(3,2*i1  ,i2,i3) = this%f_buffer3(3,2*i1  ,i2,i3) * corr1
  
        enddo
      enddo
    enddo
  
  
    call ifft( this%f_buffer3, m1, m2, m3 )
  
    ! direct copy
    do i3 = 1, m3
      do i2 = 1, m2
        do i1 = 1,m1
          do i = 1, 3
            jay%f3( i, i1-gc_l1, i2-gc_l2, i3-gc_l3 ) = this%f_buffer3(i,i1,i2,i3)
          enddo
        enddo
      enddo
    enddo
  
end subroutine filter_correct_current_3d

!-------------------------------------------------------------------------------
! Returns solver name
!-------------------------------------------------------------------------------
function name_fei( this )
  
    implicit none
  
    class( t_emf_solver_fei ), intent(in) :: this
    character( len = p_maxlen_emf_solver_name ) :: name_fei
  
    name_fei = "Customized (Fei)"
  
end function name_fei

!-------------------------------------------------------------------------------
! Returns minimum number of guard cells for solver
!-------------------------------------------------------------------------------
subroutine min_gc_fei( this, min_gc )

    implicit none

    class( t_emf_solver_fei ), intent(in) :: this
    integer, dimension(:,:), intent(inout) :: min_gc

    integer :: i

    ! Old method, relies on only one process in x1
    ! if ( this%taper_bnd ) then
    !     offset = this%taper_order / 2
    !     min_gc(p_lower,1) = max( 1, 2*this%n_coef - 1 - offset )
    !     min_gc(p_upper,1) = max( 2, 2*this%n_coef - offset )
    !     min_gc(p_lower,1) = max( min_gc(p_lower, 1), this%n_damp_cell )
    !     min_gc(p_upper,1) = max( min_gc(p_upper, 1), this%n_damp_cell )
    ! else
    !     min_gc(p_lower,1) = max( 2*this%n_coef - 1, this%n_damp_cell )
    !     min_gc(p_upper,1) = max( 2*this%n_coef, this%n_damp_cell )
    ! endif

    ! Safe assumption that we have more than one process in x1
    min_gc(p_lower,1) = max( 2*this%n_coef - 1, this%n_damp_cell )
    min_gc(p_upper,1) = max( 2*this%n_coef, this%n_damp_cell )


    do i = 2, p_x_dim
        min_gc(p_lower,i) = 1
        min_gc(p_upper,i) = 2
    enddo
  
end subroutine min_gc_fei

!-------------------------------------------------------------------------------
! Initializes the solver
!-------------------------------------------------------------------------------
subroutine init_fei( this, gix_pos, coords )

    implicit none

    class( t_emf_solver_fei ), intent(inout) :: this
    integer, intent(in), dimension(:) :: gix_pos
    integer, intent(in) :: coords

    integer :: i

    select case (p_x_dim)
    case(1)
        if ( mpi_node() == 0 ) then
            write(0,*) "(*error*) Customized (Fei) EMF solver not implemented in 1D"
            write(0,*) "(*error*) aborting..."
        endif
        call abort_program( p_err_notimplemented )
    case(2)
        if ( coords == p_cylindrical_b ) then
            ! Cylindrical (r-z) coordinates solvers
            this % dedt => dedt_cyl_2d_fei
            this % dbdt => dbdt_cyl_2d_fei
            ! Radial position of lower cell
            this % gir_pos = gix_pos(2)
        else
            this % dedt => dedt_2d_fei
            this % dbdt => dbdt_2d_fei
        endif
        this % filter_corr_current => filter_correct_current_2d
    case(3)
        this % dedt => dedt_3d_fei
        this % dbdt => dbdt_3d_fei
        this % filter_corr_current => filter_correct_current_3d
    end select

    ! generate tapered solver coefficients
    if ( this%taper_bnd ) then
        allocate( this%taper_coef_b( this%taper_order/2, this%taper_order/2 ) )
        allocate( this%taper_coef_e( this%taper_order/2, this%taper_order/2 ) )

        ! currently only standard solver is available
        do i = 1, this%taper_order/2
            call std_coef( 2*i, this%taper_coef_b(:,i) )
            call std_coef( 2*i, this%taper_coef_e(:,i) )
        enddo
    endif

end subroutine init_fei

subroutine read_input_fei( this, input_file, dx, dt )

    implicit none

    class( t_emf_solver_fei ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file
    real(p_double), dimension(:), intent(in) :: dx
    real(p_double), intent(in) :: dt

    character(len=20)                      :: type
    real(p_k_fld), dimension(p_max_n_coef) :: coef, coef_e, coef_b
    integer                                :: n_coef, solver_ord, weight_n
    real(p_double)                         :: kl, ku, dk, weight_w, dtdx1

    integer :: n_damp_cell, taper_order
    real(p_k_fld) :: filter_width, filter_limit
    logical :: correct_current, filter_current, taper_bnd

    procedure(kernel), pointer :: kernel_ptr
  
    integer :: i, ierr
  
    namelist /nl_emf_solver/ type, n_coef, coef_e, coef_b, solver_ord, kl, ku, dk, &
    filter_limit, filter_width, correct_current, filter_current, n_damp_cell, &
    weight_w, weight_n, taper_bnd, taper_order

    type = "standard"
    solver_ord = 2
    n_coef = 1
    coef_e = 0.0_p_k_fld; coef_e(1) = 1.0_p_k_fld
    coef_b = 0.0_p_k_fld; coef_b(1) = 1.0_p_k_fld
    kl = 0.0_p_double
    ku = 0.0_p_double
    dk = 0.0_p_double
    weight_n = 10
    weight_w = huge(1.0_p_double)

    filter_limit = 1.0_p_k_fld
    filter_width = 0.0_p_k_fld
    n_damp_cell  = 0
    correct_current = .false.
    filter_current = .false.  

    dtdx1 = dt / dx(1)

    taper_bnd = .false.
    taper_order = 2

    call get_namelist( input_file, "nl_emf_solver", ierr )

    if ( ierr /= 0 ) then
        if ( mpi_node() == 0 ) then
            if (ierr < 0) then
                write(0,*) "Error reading Customized (Fei) EMF solver parameters"
            else
                write(0,*) "Error: Customized (Fei) EMF solver parameters missing"
            endif
            write(0,*) "aborting..."
        endif
        stop
    endif

    read (input_file%nml_text, nml = nl_emf_solver, iostat = ierr)
    if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
            write(0,*) "Error reading Customized (Fei) EMF solver parameters"
            write(0,*) "aborting..."
        endif
        stop
    endif

    ! solver type
    if ( mod( solver_ord, 2 ) /= 0 .or. solver_ord <= 0 ) then
        if ( mpi_node() == 0 ) then
            write(0,*) "Error reading Customized (Fei) EMF solver parameters"
            print *, "Solver order must be a positive even number"
            print *, "aborting..."
        endif
        stop
    endif

    select case( trim(type) )
    case ( "standard" )

        if (disp_out(input_file)) then
            SCR_ROOT( "- using Fei solver: standard type" )
        endif
        n_coef = solver_ord/2
        call std_coef( solver_ord, coef )
        this%type = p_emf_fei_std

    case ( "bump" )

        if ( n_coef*2 <= solver_ord ) then
            if ( mpi_node() == 0 ) then
                write(0,*) "Error reading Customized (Fei) EMF solver parameters"
                print *, "Number of solver coefficients must be larger than half of solver order"
                print *, "aborting..."
            endif
            stop
        endif
            
        if (disp_out(input_file)) then
            SCR_ROOT( "- using Fei solver: bump type" )
        endif
        call bump_coef( solver_ord, n_coef, kl, ku, dk, coef )
        this%type = p_emf_fei_bump
        ! old method, not currently in use
        ! kernel_ptr => kernel_bump
        ! call gen_solver_coef( kernel_ptr, (/kl, ku, dk/), solver_ord, n_coef, coef, weight_w, weight_n )

    case ( "xu" )
        
        if ( n_coef*2 <= solver_ord ) then
            if ( mpi_node() == 0 ) then
                write(0,*) "Error reading Customized (Fei) EMF solver parameters"
                print *, "Number of solver coefficients must be larger than half of solver order"
                print *, "aborting..."
            endif
            stop
        endif

        if (disp_out(input_file)) then
            SCR_ROOT( "- using Fei solver: Xu type" )
        endif
        kernel_ptr => kernel_xu
        this%type = p_emf_fei_xu
        call gen_solver_coef( kernel_ptr, (/ dtdx1 /), solver_ord, n_coef, coef, weight_w, weight_n )

    case ( "dual" )
        
        if ( n_coef*2 <= solver_ord ) then
            if ( mpi_node() == 0 ) then
                write(0,*) "Error reading Customized (Fei) EMF solver parameters"
                print *, "Number of solver coefficients must be larger than half of solver order"
                print *, "aborting..."
            endif
            stop
        endif

        if (disp_out(input_file)) then
            SCR_ROOT( "- using Fei solver: dual type" )
        endif
        this%type = p_emf_fei_dual
        kernel_ptr => kernel_dual_e
        call gen_solver_coef( kernel_ptr, (/ dtdx1 /), solver_ord, n_coef, coef_e, weight_w, weight_n )
        kernel_ptr => kernel_dual_b
        call gen_solver_coef( kernel_ptr, (/ dtdx1 /), solver_ord, n_coef, coef_b, weight_w, weight_n )
    
    case ( "customized-coef" )

        if ( n_coef <= 0 ) then
            if ( mpi_node() == 0 ) then
                write(0,*) "Error reading Customized (Fei) EMF solver parameters"
                write(0,*) "Number of solver coefficients must be greater than zero"
                write(0,*) "aborting..."
            endif
            stop
        endif

        if (disp_out(input_file)) then
            SCR_ROOT( "- using Fei solver: customized-coefficient type" )
        endif
        this%type = p_emf_fei_coef

    case default

        if ( mpi_node() == 0 ) then
            write(0,*) "Error reading Customized (Fei) EMF solver parameters"
            write(0,*) "Invalid value for the solver type"
            write(0,*) "Available solver types are 'standard', 'bump', 'xu', 'dual' and 'customized-coef'"
            write(0,*) "aborting..."
        endif
        stop

    end select

    if ( mpi_node() == 0 .and. disp_out(input_file) ) then
        print *, "  Solver stencil coefficients are:"
        select case ( this%type )
        case ( p_emf_fei_std, p_emf_fei_bump, p_emf_fei_xu )
            do i = 1, n_coef
                print *, "  C(", i, ") = ", coef(i)
            enddo
        case ( p_emf_fei_dual, p_emf_fei_coef )
            do i = 1, n_coef
                print *, "  C_e(", i, ") = ", coef_e(i)
            enddo
            do i = 1, n_coef
                print *, "  C_b(", i, ") = ", coef_b(i)
            enddo
        end select
    endif

    ! store namelist parameters
    allocate( this%coef_e(n_coef), this%coef_b(n_coef) )
    this%n_coef = n_coef
    select case ( this%type )
    case ( p_emf_fei_std, p_emf_fei_bump, p_emf_fei_xu )
        this%coef_e = coef(1:n_coef)
        this%coef_b = coef(1:n_coef)
    case ( p_emf_fei_dual, p_emf_fei_coef )
        this%coef_e = coef_e(1:n_coef)
        this%coef_b = coef_b(1:n_coef)
    end select
    this%solver_ord = solver_ord
    this%ku         = ku
    this%kl         = kl
    this%dk         = dk
    this%weight_w   = weight_w
    this%weight_n   = weight_n

    ! Current filtering
    this%n_damp_cell = n_damp_cell
    this%correct_current = correct_current
    this%filter_current = filter_current

    ! check validity of spectral filter parameters
    if ( filter_current ) then
        if ( filter_width < 0.0 .or. filter_width > 1.0 .or. &
             filter_limit < 0.0 .or. filter_limit > 1.0 ) then

            if ( mpi_node() == 0 ) then
                write(0,*) "Error reading current spectral filter parameters"
                write(0,*) "filter_limit and filter_width must be set between 0 and 1"
                write(0,*) "aborting..."
            endif
            stop

        else

            this%filter_width = filter_width
            this%filter_limit = filter_limit

        endif
    endif

    ! tapered solver coef. at boundaries
    this%taper_bnd = taper_bnd
    if ( this%taper_bnd ) then
        this%taper_order = taper_order
        if ( taper_order < 2 .or. mod(taper_order,2) /= 0 ) then
            if ( mpi_node() == 0 ) then
                write(0,*) "Invalid parameter taper_order."
                write(0,*) "aborting..."
            endif
        endif
    endif

end subroutine read_input_fei

subroutine advance_fei( this, e, b, jay, dt, bnd_con )
  
    implicit none
  
    class( t_emf_solver_fei ), intent(inout) :: this
    class( t_vdf ), intent( inout )  ::  e, b, jay
    class( t_emf_bound ), intent(inout) :: bnd_con
    real(p_double), intent(in) :: dt
  
    real(p_double) :: dt_b, dt_e

    ! Filter current if necessary
    if ( this%filter_current .or. this%correct_current ) then
        call this % filter_corr_current( jay )
    endif
    
    ! Advance fields
    dt_b = dt / 2.0_p_double
    dt_e = dt
  
    ! Advance B half time step
    call this % dbdt( b, e, dt_b, bnd_con%type )
    call bnd_con%update_boundary_b( e, b, step = 1 )
  
    ! Advance E one full time step
    call this % dedt( e, b, jay, dt_e, bnd_con%type )
    call bnd_con%update_boundary_e( e, b )
  
    ! Advance B another half time step
    call this % dbdt( b, e, dt_b, bnd_con%type )
    call bnd_con%update_boundary_b( e, b, step = 2 )
    
end subroutine advance_fei

end module m_emf_solver_fei