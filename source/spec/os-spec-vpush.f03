!-----------------------------------------------------------------------------------------
! Particle pusher using hardware (vector) acceleration
!-----------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_vpush
    
#include "memory/memory.h"

use m_parameters
use m_species_define
use m_emf_define
use m_vdf_define

implicit none

private

interface vadvance_deposit
    module procedure vadvance_deposit
end interface

public :: vadvance_deposit

contains

!---------------------------------------------------------------------------------------------------
! Interfaces to hardware optimized pushers
!---------------------------------------------------------------------------------------------------

#ifdef SIMD

!---------------------------------------------------------------------------------------------------
!> Call 1D vector advance/deposit routine based on interpolation level
!!
!! @param   this    species object
!! @param   emf     EM fields
!! @param   jay     Electric current
!! @param   energy  Total species energy
!! @param   gdt     Time step
!! @param   i0      Index of 1st particle to push
!! @param   i1      Index of last particle to push
subroutine vadvance_deposit_1d( this, emf, jay, energy, gdt, i0, i1 )
    
    implicit none
    
    integer, parameter :: rank = 1
   
    class( t_species ),    intent(inout) :: this
    class( t_emf ), intent( in ), target  ::  emf
    type( t_vdf ),        intent(inout) :: jay
    real(p_double), intent(inout) :: energy
    real(p_double),   intent(in) :: gdt
    integer, intent(in) :: i0, i1
    
    integer, dimension(rank) :: size_b, offset_b, size_jay, offset_jay
    
    size_b(1)     = size( emf%b_part%f1, 2)
    offset_b(1)   = emf%b_part%gc_num_(p_lower,1)
    size_jay(1)   = size( jay%f1, 2)
    offset_jay(1) = jay%gc_num_(p_lower,1)
    
    ! executable statements
    select case (this%interpolation)
        
    case( p_linear )
        
        call vadvance_deposit_1d_s1( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
            emf%e_part%f1, emf%b_part%f1, size_b, offset_b, &
            jay%f1, size_jay, offset_jay, &
            this%dx, gdt, energy )
        
    case( p_quadratic )
        
        call vadvance_deposit_1d_s2( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
            emf%e_part%f1, emf%b_part%f1, size_b, offset_b, &
            jay%f1, size_jay, offset_jay, &
            this%dx, gdt, energy )
        
        
    case( p_cubic )
        
        call vadvance_deposit_1d_s3( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
            emf%e_part%f1, emf%b_part%f1, size_b, offset_b, &
            jay%f1, size_jay, offset_jay, &
            this%dx, gdt, energy )
        
    case( p_quartic )
        
        call vadvance_deposit_1d_s4( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
            emf%e_part%f1, emf%b_part%f1, size_b, offset_b, &
            jay%f1, size_jay, offset_jay, &
            this%dx, gdt, energy )
        
    case default
        
        ERROR('Not implemented yet')
        call abort_program( p_err_notimplemented )
        
    end select
    
    
end subroutine vadvance_deposit_1d

!---------------------------------------------------------------------------------------------------
!> Call 2D vector advance/deposit routine based on interpolation level
!!
!! @param   this    species object
!! @param   emf     EM fields
!! @param   jay     Electric current
!! @param   energy  Total species energy
!! @param   gdt     Time step
!! @param   i0      Index of 1st particle to push
!! @param   i1      Index of last particle to push
subroutine vadvance_deposit_2d( this, emf, jay, energy, gdt, i0, i1 )
    
    implicit none
    
    integer, parameter :: rank = 2
    
    class( t_species ),    intent(inout) :: this
    class( t_emf ), intent( in ), target  ::  emf
    type( t_vdf ),        intent(inout) :: jay
    real(p_double), intent(inout) :: energy
    real(p_double),   intent(in) :: gdt
    integer, intent(in) :: i0, i1
    
    integer, dimension(rank) :: size_b, offset_b, size_jay, offset_jay
    integer :: i
    
    do i = 1, rank
        size_b(i)     = size( emf%b_part%f2, i + 1)
        offset_b(i)   = emf%b_part%gc_num_(p_lower,i)
        size_jay(i)   = size( jay%f2, i + 1)
        offset_jay(i) = jay%gc_num_(p_lower,i)
    enddo

    select case (this%interpolation)
        
    case( p_linear )
        
        call vadvance_deposit_2d_s1( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
            emf%e_part%f2, emf%b_part%f2, size_b, offset_b, &
            jay%f2, size_jay, offset_jay, &
            this%dx, gdt, energy )
        
    case( p_quadratic )
        
        call vadvance_deposit_2d_s2( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
            emf%e_part%f2, emf%b_part%f2, size_b, offset_b, &
            jay%f2, size_jay, offset_jay, &
            this%dx, gdt, energy )
        
        
    case( p_cubic )
        
        call vadvance_deposit_2d_s3( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
            emf%e_part%f2, emf%b_part%f2, size_b, offset_b, &
            jay%f2, size_jay, offset_jay, &
            this%dx, gdt, energy )
        
    case( p_quartic )
        
        call vadvance_deposit_2d_s4( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
            emf%e_part%f2, emf%b_part%f2, size_b, offset_b, &
            jay%f2, size_jay, offset_jay, &
            this%dx, gdt, energy )
        
    case default
        
        ERROR('Not implemented yet')
        call abort_program( p_err_notimplemented )
        
    end select
    
end subroutine vadvance_deposit_2d

!---------------------------------------------------------------------------------------------------
!> Call 2D cylindrical vector advance/deposit routine based on interpolation level
!!
!! @param   this    species object
!! @param   emf     EM fields
!! @param   jay     Electric current
!! @param   energy  Total species energy
!! @param   gdt     Time step
!! @param   i0      Index of 1st particle to push
!! @param   i1      Index of last particle to push
subroutine vadvance_deposit_2d_cyl( this, emf, jay, energy, gdt, i0, i1 )
    
    implicit none
    
    integer, parameter :: rank = 2
    
    class( t_species ),    intent(inout) :: this
    class( t_emf ), intent( in ), target  ::  emf
    type( t_vdf ),        intent(inout) :: jay
    real(p_double), intent(inout) :: energy
    real(p_double),   intent(in) :: gdt
    integer, intent(in) :: i0, i1

    integer, dimension(rank) :: size_b, offset_b, size_jay, offset_jay
    integer :: i
    
    do i = 1, rank
        size_b(i)     = size( emf%b_part%f2, i + 1)
        offset_b(i)   = emf%b_part%gc_num_(p_lower,i)
        size_jay(i)   = size( jay%f2, i + 1)
        offset_jay(i) = jay%gc_num_(p_lower,i)
    enddo
    
    select case (this%interpolation)
       
    case( p_linear )
        
        call vadvance_deposit_2d_cyl_s1( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
            this%my_nx_p(p_lower, 2), &
            emf%e_part%f2, emf%b_part%f2, size_b, offset_b, &
            jay%f2, size_jay, offset_jay, &
            this%dx, gdt, energy )
        
    case( p_quadratic )
        
        call vadvance_deposit_2d_cyl_s2( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
            this%my_nx_p(p_lower, 2), &
            emf%e_part%f2, emf%b_part%f2, size_b, offset_b, &
            jay%f2, size_jay, offset_jay, &
            this%dx, gdt, energy )
        
        
    case( p_cubic )
        
        call vadvance_deposit_2d_cyl_s3( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
            this%my_nx_p(p_lower, 2), &
            emf%e_part%f2, emf%b_part%f2, size_b, offset_b, &
            jay%f2, size_jay, offset_jay, &
            this%dx, gdt, energy )
        
    case( p_quartic )
        
        call vadvance_deposit_2d_cyl_s4( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
            this%my_nx_p(p_lower, 2), &
            emf%e_part%f2, emf%b_part%f2, size_b, offset_b, &
            jay%f2, size_jay, offset_jay, &
            this%dx, gdt, energy )
        
    case default
        
        ERROR('Not implemented yet')
        call abort_program( p_err_notimplemented )
        
    end select
    
end subroutine vadvance_deposit_2d_cyl

!---------------------------------------------------------------------------------------------------
!> Call 3D vector advance/deposit routine based on interpolation level
!!
!! @param   this    species object
!! @param   emf     EM fields
!! @param   jay     Electric current
!! @param   energy  Total species energy
!! @param   gdt     Time step
!! @param   i0      Index of 1st particle to push
!! @param   i1      Index of last particle to push
subroutine vadvance_deposit_3d( this, emf, jay, energy, gdt, i0, i1 )
    
    implicit none
    
    integer, parameter :: rank = 3
    
    class( t_species ),    intent(inout) :: this
    class( t_emf ), intent( in )  ::  emf
    type( t_vdf ),        intent(inout) :: jay
    real(p_double), intent(inout) :: energy
    real(p_double),   intent(in) :: gdt
    integer, intent(in) :: i0, i1
    
    integer, dimension(rank) :: size_b, offset_b, size_jay, offset_jay
    integer :: i
    
    do i = 1, rank
        size_b(i)     = size( emf%b_part%f3, i + 1)
        offset_b(i)   = emf%b_part%gc_num_(p_lower,i)
        size_jay(i)   = size( jay%f3, i + 1)
        offset_jay(i) = jay%gc_num_(p_lower,i)
    enddo
    
    select case (this%interpolation)

    case( p_linear )
        
        call vadvance_deposit_3d_s1( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
        emf%e_part%f3, emf%b_part%f3, size_b, offset_b, &
        jay%f3, size_jay, offset_jay, &
        this%dx, gdt, energy )
        
    case( p_quadratic )
        
        call vadvance_deposit_3d_s2( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
        emf%e_part%f3, emf%b_part%f3, size_b, offset_b, &
        jay%f3, size_jay, offset_jay, &
        this%dx, gdt, energy )
        
        
    case( p_cubic )
        
        call vadvance_deposit_3d_s3( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
        emf%e_part%f3, emf%b_part%f3, size_b, offset_b, &
        jay%f3, size_jay, offset_jay, &
        this%dx, gdt, energy )
        
    case( p_quartic )
        
        call vadvance_deposit_3d_s4( this%ix, this%x, this%p, this%q, i0, i1, this%rqm, &
        emf%e_part%f3, emf%b_part%f3, size_b, offset_b, &
        jay%f3, size_jay, offset_jay, &
        this%dx, gdt, energy )
        
    case default
        ERROR('Not implemented yet')
        call abort_program( p_err_notimplemented )
        
    end select
    
end subroutine vadvance_deposit_3d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!> Push particles and deposit current using vector (SIMD) code
!!
!! @param   this        species object
!! @param   emf         EM fields
!! @param   current     Electric current
!! @param   t           Simulation time
!! @param   tstep       Time step object
!! @param   tid         Thread id of the current thread 
!! @param   n_threads   Total number of threads in the group
subroutine vadvance_deposit( this, emf, current, t, tstep, tid, n_threads )
    
    use m_time_step
    use m_emf
    
    implicit none
    
    class( t_species ), intent(inout) :: this
    class( t_emf ), intent( in )  ::  emf
    type( t_vdf ), intent(inout) :: current
    
    real(p_double), intent(in) :: t
    type( t_time_step ) :: tstep
    
    integer, intent(in) :: tid        ! local thread id
    integer, intent(in) :: n_threads  ! total number of threads
    
    
    real(p_double) :: dtcycle, energy
    integer :: chunk, i0, i1
        
    ! sanity checks - these should be done at the setup stage...
    if ( this%num_pistons > 0 ) then
        write(0,*) '(*error*) Pistons are not supported by the vector code, please use the standard fortran pusher'
        call abort_program( p_err_invalid )
    endif
    
    ! Sanity check
    ! Buffer size should be a multiple of vector size.
    if ( mod( this%num_par_max, p_vecwidth ) /= 0 ) then
        write(0,*) '(*error*) The particle buffer size is not a multiple of vector width = ', p_vecwidth
        call abort_program()
    endif
    
    dtcycle = real( dt(tstep), p_k_part )
    
    ! Range of particles for each thread
    ! The chunk size must be a multiple of p_vecwidth so that all threads but the last
    ! fill the vectors exactly and only the last one will do any padding. The code below
    ! also handles the situation where num_par < n_threads * p_vecwidth correctly.
    if ( n_threads > 1 ) then
        
        chunk = ( this%num_par + n_threads - 1 ) / n_threads
        ! round up to the nearest multiple of vector width
        chunk = ( ( chunk + p_vecwidth - 1 ) / p_vecwidth ) * p_vecwidth
        i0    = tid * chunk + 1
        i1    = min( (tid+1) * chunk, this%num_par )
        
        ! Sanity check
        if ( mod( chunk, p_vecwidth ) /= 0 ) then
            write(0,'(A,I0)') '(*error*) The chunk size is not a multiple of vector width = ', p_vecwidth
            call abort_program()
        endif
        
    else
        i0 = 1
        i1 = this%num_par
    endif
    
    ! initialize time centered energy diagnostic
    energy = 0.0_p_double

    if (i1 >= i0) then
        
        select case ( this%coordinates )
            
        case default
            select case ( p_x_dim )
                
            case (1)
                call vadvance_deposit_1d( this, emf, current, energy, dtcycle, i0, i1 )
                
            case (2)
                call vadvance_deposit_2d( this, emf, current, energy, dtcycle, i0, i1 )
                
            case (3)
                call vadvance_deposit_3d( this, emf, current, energy, dtcycle, i0, i1 )
                
            case default
                ERROR('space_dim(jay) has the value:',p_x_dim)
                call abort_program(p_err_invalid)
                
            end select
            
        case ( p_cylindrical_b )
            call vadvance_deposit_2d_cyl( this, emf, current, energy, dtcycle, i0, i1 )
            
        end select
        
    endif
    
    this%energy(tid+1) = energy
    
end subroutine vadvance_deposit
!---------------------------------------------------------------------------------------------------

#else

!---------------------------------------------------------------------------------------------------
! Push particles and deposit electric current using vector code
!---------------------------------------------------------------------------------------------------
subroutine vadvance_deposit( this, emf, jay, t, tstep, tid, n_threads )
    
    use m_time_step
    
    implicit none
    
    class( t_species ), intent(inout) :: this
    class( t_emf ), intent( in )  ::  emf
    type( t_vdf ), dimension(:), intent(inout) :: jay
    
    real(p_double), intent(in) :: t
    type( t_time_step ) :: tstep
    
    integer, intent(in) :: tid        ! local thread id
    integer, intent(in) :: n_threads  ! total number of threads
    
    WARNING( 'Vector code not available' )
    
end subroutine vadvance_deposit

#endif
    
end module m_species_vpush
        