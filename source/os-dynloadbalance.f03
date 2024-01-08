!#define __DEBUG__ 1

module m_dynamic_loadbalance
    
#include "os-config.h"
#include "os-preprocess.fpp"

use m_parameters
use m_grid_define
use m_grid
use m_node_conf

implicit none

public :: test_reshape_vdf, sim_imbalance

contains

!-----------------------------------------------------------------------------------------
!> Test the loadbalance of vdfs. The vdf_ref is used as a template vdf; it must use the
!! old_lb grid (i.e. it must not have been reshaped yet)
!-----------------------------------------------------------------------------------------
subroutine test_reshape_vdf( vdf_ref, old_lb, new_lb, no_co )
    
    use m_vdf_comm
    use m_vdf
    use m_vdf_define
    
    implicit none
   
    type (t_vdf), intent(in) :: vdf_ref
    class (t_grid ), intent(in) :: old_lb, new_lb
    class (t_node_conf), intent(in) :: no_co
    type (t_vdf_msg), dimension(2) :: send_msg, recv_msg
    
    ! local variables
    type (t_vdf ) :: vdf
    integer :: i1, i2, i3, k
    integer :: gi1, gi2, gi3
    real( p_k_fld ) :: v1
    
    call vdf % new( vdf_ref )
    
    if ( vdf%x_dim_ /= 3 ) then
        ERROR( 'test_reshape_vdf is only implemented for 3D' )
        call abort_program()
    endif
    
    vdf = huge(1.0_p_k_fld)
    
    ! Set reference values
    do i3 = lbound( vdf%f3, 4 ), ubound( vdf%f3, 4 )
        do i2 = lbound( vdf%f3, 3 ), ubound( vdf%f3, 3 )
            do i1 = lbound( vdf%f3, 2 ), ubound( vdf%f3, 2 )
                
                gi1 = i1 + old_lb%my_nx(p_lower,1) - 1
                gi2 = i2 + old_lb%my_nx(p_lower,2) - 1
                gi3 = i3 + old_lb%my_nx(p_lower,3) - 1
                
                v1 = gi1 + &
                gi2 * old_lb%g_nx(1) + &
                gi3 * old_lb%g_nx(1) * old_lb%g_nx(2)
                
                vdf%f3( :, i1, i2, i3 ) = v1
                
            enddo
        enddo
    enddo
    
    ! reshape the vdf
    call reshape_copy( vdf, old_lb, new_lb, no_co, send_msg, recv_msg )
    
    ! Test the operation
    do i3 = lbound( vdf%f3, 4 ), ubound( vdf%f3, 4 )
        do i2 = lbound( vdf%f3, 3 ), ubound( vdf%f3, 3 )
            do i1 = lbound( vdf%f3, 2 ), ubound( vdf%f3, 2 )
                
                gi1 = i1 + new_lb%my_nx(p_lower,1) - 1
                gi2 = i2 + new_lb%my_nx(p_lower,2) - 1
                gi3 = i3 + new_lb%my_nx(p_lower,3) - 1
                
                ! if doing a periodic run, the values in the guard cells may correspond to values in the
                ! opposite side of the box
                if ( periodic( no_co, 1 ) ) then
                    if ( gi1 < 1 ) gi1 = gi1 + new_lb%g_nx(1)
                    if ( gi1 > new_lb%g_nx(1) ) gi1 = gi1 - new_lb%g_nx(1)
                endif
                
                if ( periodic( no_co, 2 ) ) then
                    if ( gi2 < 1 ) gi2 = gi2 + new_lb%g_nx(2)
                    if ( gi2 > new_lb%g_nx(2) ) gi2 = gi2 - new_lb%g_nx(2)
                endif
                
                if ( periodic( no_co, 3 ) ) then
                    if ( gi3 < 1 ) gi3 = gi3 + new_lb%g_nx(3)
                    if ( gi3 > new_lb%g_nx(3) ) gi3 = gi3 - new_lb%g_nx(3)
                endif
                
                v1 = gi1 + &
                gi2 * old_lb%g_nx(1) + &
                gi3 * old_lb%g_nx(1) * old_lb%g_nx(2)
                
                do k = 1, vdf%f_dim()
                    if ( vdf%f3( k, i1, i2, i3 ) /= v1 ) then
                        SCR_MPINODE('reshape_vdf_copy failed, bad value at position ', i1, ', ', i2, ', ', i3, ' fc: ', k)
                        SCR_MPINODE('expected : ', v1, ' got : ', vdf%f3( k, i1, i2, i3 ) )
                        print *, '[', mpi_node(), '] Global cell : ', i1 + new_lb%my_nx(p_lower,1) - 1, ' , ', &
                        i2 + new_lb%my_nx(p_lower,2) - 1, ' , ', &
                        i3 + new_lb%my_nx(p_lower,3) - 1
                        call abort_program()
                    endif
                enddo
            enddo
        enddo
    enddo
    
    ! cleanup the temporary vdf
    call vdf % cleanup()

    ! cleanup message buffers
    call cleanup( send_msg )
    call cleanup( recv_msg )
    
end subroutine test_reshape_vdf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Gets the current load imbalance for the simulation
!-----------------------------------------------------------------------------------------
function sim_imbalance( grid, particles, no_co )
    
    use m_particles
    use m_particles_define
    
    implicit none
    
    real( p_single ) :: sim_imbalance
    
    class( t_grid ), intent(in) :: grid
    class( t_particles ), intent(in) :: particles
    class( t_node_conf ), intent(in) :: no_co
    
    real( p_single ) :: local_load, max_load, avg_load
    
    ! Get local load
    local_load = real( num_par( particles ), p_single ) + grid%cell_weight * local_vol( grid )
    
    ! Get maximum load
    max_load = local_load
    call reduce( no_co, max_load, p_max, all = .true. )
    
    ! Get average load
    avg_load = local_load
    call reduce( no_co, avg_load, p_sum, all = .true. )
    avg_load = avg_load / no_num( no_co )
    
    ! Calculate simulation imbalance
    sim_imbalance = max_load / avg_load
    
end function sim_imbalance
!-----------------------------------------------------------------------------------------
    
end module m_dynamic_loadbalance

!-----------------------------------------------------------------------------------------
!  dynamically load balance the simulation
!-----------------------------------------------------------------------------------------
subroutine dlb_sim( sim )
    
    use m_particles
    use m_emf
    use m_logprof
    use m_dynamic_loadbalance
    use m_simulation
    use m_grid_define
    use m_node_conf
    use m_time_step
    use m_grid_parallel

    use m_vdf_comm

    implicit none
    
    class ( t_simulation ), intent(inout) ::  sim
    
    type ( t_grid ) :: new_lb
    logical :: do_dyn_lb
        
    if (if_dynamic_lb(sim%grid, n(sim%tstep), sim%no_co)) then
        
        call begin_event(dynlbev)
        
        ! Check if the simulation imbalance is over the set threshold. (The imbalance is defined as
        ! the ratio between the maximum load and the average load)
        if ( sim%grid%max_imbalance >= 1 ) then
            do_dyn_lb = sim_imbalance( sim%grid, sim%part, sim%no_co ) > sim%grid%max_imbalance
        else
            do_dyn_lb = .true.
        endif
        
        if ( do_dyn_lb ) then
            
            if ( root( sim%no_co ) ) then
                print *, ''
                print *, ' Dynamically load balancing the simulation.'
                print *, ''
            endif
            
            ! copy data to new_lb
            call new_lb % copy( sim%grid )
            
            ! determine the current integral particle load
            call sim % set_int_load( new_lb, n(sim%tstep) )
            
            ! determine optimum distribution
            call begin_event(dynlb_new_lb_ev)
            call new_lb % parallel_partition( sim%no_co, n(sim%tstep) )
            call end_event(dynlb_new_lb_ev)

            ! clear the memory used by the integral particle load
            call clear_int_load( new_lb )

            ! Check if partitions are different
            if ( diff_partition(sim%grid, new_lb ) ) then
                ! We should also check if there is a significant global
                ! load balance advantage and only change it in that case               
                
                call begin_event(dynlb_reshape_ev)
                
                call begin_event(dynlb_reshape_part_ev)
                
#ifdef __DEBUG__
                ! (* debug *)
                call sim % no_co % barrier(  )
                SCR_ROOT(" (*debug*)")
                
                ! Check if reshaping a vdf with the same shape as the e field works ok
                !SCR_ROOT(" Verifying test vdf reshape...")
                !call test_reshape_vdf( sim%emf%e, sim%grid, new_lb, sim%no_co )
                
                ! check particles are ok for debug purposes
                SCR_ROOT(" Verifying particle data before reshape...")
                call validate( sim%part, "Veryfying dynamic load balance (A)" )
#endif
                
                ! move particles to proper process
                call reshape_obj( sim%part, sim%grid, new_lb, sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )
                
                call end_event(dynlb_reshape_part_ev)
                
                call begin_event(dynlb_reshape_grid_ev)
                
                ! reshape jay vdf (communication may be necessary)
                call sim%jay%reshape( sim%grid, new_lb, sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )
                
                ! reshape E and B and copy values to proper process
                call reshape_obj( sim%emf, sim%grid, new_lb, sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )
                call end_event(dynlb_reshape_grid_ev)
                
                call end_event(dynlb_reshape_ev)
                
                ! store the new load balance
                ! note that this must happen even if the local nx_p remains constant
                call sim%grid % copy ( new_lb )

                ! Recreate grid io object
                call sim % grid % io % init( sim % grid, sim % no_co )
                                
                ! synchronize all nodes
                call sim % no_co  % barrier( )
                
#ifdef __DEBUG__
            
                ! check particles are ok for debug purposes
                SCR_ROOT(" Verifying particle data after reshape...")
                call validate( sim%part, "Veryfying dynamic load balance (B)" )
                
                ! check all nodes have the same grid object
                SCR_ROOT(" Verifying consistency of grid object...")
                call check_consistency( sim%grid, sim%no_co )
                
                ! synchronize all nodes for debug purposes
                call  sim%no_co  % barrier( )
                SCR_ROOT(" No errors found!")
                SCR_ROOT("")
#endif           
            
            else
                if ( root( sim%no_co ) ) then
                    print *, ' (no partition change)'
                    print *, ''
                endif
            endif

            ! Cleanup temp. objects
            call new_lb % cleanup()

        endif
        
        call end_event(dynlbev) 
    endif
    
end subroutine dlb_sim
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Report on the load of the simulation. Currently this accounts for particles only
!-----------------------------------------------------------------------------------------
subroutine report_load_sim( sim )
    
    use m_grid_define
    use m_particles
    use m_dynamic_loadbalance
    use m_simulation
    use m_time_step
    use m_time
    use m_logprof
    
    implicit none
    
    class ( t_simulation ), intent(inout) ::  sim
    
    call begin_event(report_load_ev)
    
    ! global load : total, max, min, avg parts per node
    if ( sim%grid%if_report( n(sim%tstep), p_global ) ) then
        call report_global_load( sim%part, n(sim%tstep), sim%no_co )
    endif
    
    ! node load : number of particles per node for all nodes
    if ( sim%grid%if_report( n(sim%tstep), p_node ) ) then
        call report_node_load( sim%part, n(sim%tstep), sim%grid%ndump( p_node ), t(sim%time), &
        sim%no_co )
    endif
    
    ! grid load : number of particles per cell for all grid points
    if ( sim%grid%if_report( n(sim%tstep), p_grid ) ) then
        call report_grid_load( sim%part, n(sim%tstep), sim%grid%ndump( p_grid ), t(sim%time), &
        sim%g_space, sim%grid, sim%no_co )
    endif
    
    call end_event(report_load_ev)
    
end subroutine report_load_sim
!-----------------------------------------------------------------------------------------
