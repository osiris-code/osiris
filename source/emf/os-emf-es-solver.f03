!-----------------------------------------------------------------------------------------
! Open boundary electrostatic field solver
!
! This solver calculates the field in each grid point as the total field contributions
! of all other grid points.
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! The E field VDF is expected to be a staggered Yee grid
!-----------------------------------------------------------------------------------------

!
! Electrostatic field solver test
!
!        if ( test_if_report( n, ndump, this%ndump_fac_charge ) ) then
!           print *, '--------------ES Solver test--------------'
!           charge_report%label = 'E_1'
!           charge_report%name  = 'E1'
!           charge_report%units  = 'm_e c \omega_p e^{-1}'
!
!           call es_solver( this%charge0, ef, grid, no_co )
!
!           charge_report%path = trim(path_mass) // 'ESF' // p_dir_sep // 'E1-static' // &
!                                  p_dir_sep
!
!           charge_report%filename = get_filename( n / ndump, 'e1' )
!
!           call report( ef, 1, g_space, grid, no_co, charge_report )
!
!           charge_report%label = 'E_2'
!           charge_report%name  = 'E2'
!
!           charge_report%path = trim(path_mass) // 'ESF' // p_dir_sep // 'E2-static' // &
!                                  p_dir_sep
!
!           charge_report%filename = get_filename( n / ndump, 'e2' )
!
!           call report( ef, 2, g_space, grid, no_co, charge_report )
!
!           call cleanup( ef )
!        endif

#include "os-preprocess.fpp"
#include "os-config.h"

module m_emf_es_solver

#include "memory/memory.h"

use m_vdf_define
use m_vdf

use m_utilities
use m_node_conf

use m_grid_define
use m_grid

use m_parameters
use m_system

implicit none

private

interface es_solver
  module procedure es_solver
end interface

interface es_solver_beam
  module procedure es_solver_beam
end interface

public :: es_solver, es_solver_beam

contains

!-----------------------------------------------------------------------------------------
subroutine es_solver( rho, ef, grid, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

  use m_vdf_comm

  implicit none

  type(t_vdf), intent(in) :: rho
  type(t_vdf), intent(inout) :: ef
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  ! size of the grid on each node
  integer, dimension(:,:,:), pointer :: lnx_p

  ! get size of the grid on each node
  call alloc(lnx_p, (/ 2, rho%x_dim_, no_num(no_co) /))

  call grid % get_node_limits( no_co, lnx_p )

  ! Allocate the ef grid
  call ef % new( rho, f_dim = 3, copy = .false. )

  select case (rho%x_dim_)
    case(1)
       ERROR('1D ES field solver not implemented yet')
       call abort_program( p_err_notimplemented )
    case(2)
       call es_solver_2d( rho, ef, lnx_p, no_co )

    case(3)
       ERROR('3D ES field solver not implemented yet')
       call abort_program( p_err_notimplemented )

  end select

  call freemem( lnx_p )

  ! update fields in guard cells
  call update_boundary( ef, p_vdf_replace, no_co, send_msg, recv_msg )

end subroutine es_solver
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine es_solver_2d( rho, ef, nx_p, no_co )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  type(t_vdf), intent(in) :: rho
  type(t_vdf), intent(inout) :: ef
  integer, dimension(:,:,:), intent(in) :: nx_p
  class( t_node_conf ), intent(in) :: no_co

  integer :: ri1, ri2, ei1, ei2

  real( p_double ) :: r_2, dx1_e1, dx1_e2, dx2_e1, dx2_e2, cell_vol
  real( p_double ), pointer, dimension(:,:) :: charge
  integer, dimension(rank) :: lnx0 ! position of local data on global grid

  integer, dimension(rank) :: lb, ub

  integer :: i, ierr, count

  integer( p_int64 ) :: t0, t1

  t0 = timer_ticks()

  ! find position of local data on global grid
  do i = 1, rank
    lnx0( i ) = nx_p( p_lower, i, no_co%my_aid() ) - 1
  enddo

  ! set electric field to 0
  call ef % zero()

  ! loop over every nodes
  do i = 1, no_num(no_co)

    lb(1) = nx_p( p_lower, 1, i )
    lb(2) = nx_p( p_lower, 2, i )

    ub(1) = nx_p( p_upper, 1, i )
    ub(2) = nx_p( p_upper, 2, i )

    count =  (ub(1)-lb(1)+1) * (ub(2)-lb(2)+1)

    call alloc( charge, lb, ub )

    ! copy local charge of node being processed to all nodes
    if ( no_co%my_aid() == i ) then
      do ri2 = 1, rho%nx_(2)
        do ri1 = 1, rho%nx_(1)
          charge( lnx0(1) + ri1, lnx0(2) + ri2 ) = rho%f2(1,ri1,ri2)
        enddo
      enddo
    endif

    if ( no_num(no_co) > 1 ) then
      call MPI_BCAST( charge, count, MPI_DOUBLE_PRECISION, i-1, comm(no_co), ierr)
    endif

    ! calculate field generated by current data
    do ri2 = nx_p( p_lower, 2, i ), nx_p( p_upper, 2, i )
      do ri1 = nx_p( p_lower, 1, i ), nx_p( p_upper, 1, i )
         if ( charge(ri1,ri2) /= 0 ) then
            do ei2 = 1, ef%nx_(2)
              dx2_e1 = real(lnx0(2)+ei2-ri2, p_double )
              dx2_e2 = dx2_e1 + 0.5_p_double

              dx2_e1 = dx2_e1 * ef%dx_(2)
              dx2_e2 = dx2_e2 * ef%dx_(2)

              do ei1 = 1, ef%nx_(1)
                 dx1_e2 = real(lnx0(1)+ei1-ri1, p_double )
                 dx1_e1 = dx1_e2 + 0.5_p_double

                 dx1_e1 = dx1_e1 * ef%dx_(1)
                 dx1_e2 = dx1_e2 * ef%dx_(1)

                 r_2 = dx1_e1**2 + dx2_e1**2
                 ef%f2(1, ei1, ei2) = ef%f2(1, ei1, ei2) + charge(ri1,ri2)*dx1_e1/r_2

                 r_2 = dx1_e2**2 + dx2_e2**2
                 ef%f2(2, ei1, ei2) = ef%f2(2, ei1, ei2) + charge(ri1,ri2)*dx2_e2/r_2
              enddo
            enddo
         endif
      enddo
    enddo

    ! free temp memory
    call freemem( charge )

  enddo

  ! Normalize results to cell volume (charge is charge density) and clear E3

  cell_vol = ef%dx_(1) * ef%dx_(2)
  do ei2 = 1, ef%nx_(2)
    do ei1 = 1, ef%nx_(1)
       ef%f2(1, ei1, ei2) = ef%f2(1, ei1, ei2) * cell_vol
       ef%f2(2, ei1, ei2) = ef%f2(2, ei1, ei2) * cell_vol
       ef%f2(3, ei1, ei2) = 0
    enddo
  enddo

  t1 = timer_ticks()

  if ( root(no_co) ) then
    print *, 'Time for 2D ES field solve : ', timer_interval_seconds( t0, t1 )
  endif

end subroutine es_solver_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine es_solver_beam( rho, ef, bf, u1, grid, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

  use m_vdf_comm

  implicit none

  type(t_vdf), intent(in) :: rho
  type(t_vdf), intent(inout) :: ef, bf
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  real(p_k_part), intent(in) :: u1
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  ! size of the grid on each node
  integer, dimension(:,:,:), pointer :: lnx_p

  ! get size of the grid on each node
  call alloc(lnx_p, (/ 2, rho%x_dim_, no_num(no_co) /))

  call grid % get_node_limits( no_co, lnx_p )

  ! Allocate the ef grid
  call ef % new( rho, f_dim = 3, copy = .false. )

  select case (rho%x_dim_)
    case(1)
       ERROR('1D ES field solver (beam) not implemented yet')
       call abort_program( p_err_notimplemented )
    case(2)

       if ( coordinates(grid) == p_cylindrical_b ) then
         call es_solver_2dcyl_beam( rho, ef, bf, u1, lnx_p, no_co, send_msg, recv_msg )
       else
         call es_solver_2d_beam( rho, ef, bf, u1, lnx_p, no_co, send_msg, recv_msg )
       endif

    case(3)
       call es_solver_3d_beam( rho, ef, bf, u1, lnx_p, no_co, send_msg, recv_msg )

  end select

  call freemem( lnx_p )

  ! update fields in guard cells
  call update_boundary( ef, p_vdf_replace, no_co, send_msg, recv_msg )

end subroutine es_solver_beam
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! This version assumes the charge is a line charge along x1 with dimension equal to
! cell size.
!-----------------------------------------------------------------------------------------
subroutine es_solver_2d_beam_full( rho, ef, bf, u1, nx_p, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

  use m_vdf_comm

  implicit none

  integer, parameter :: rank = 2

  type(t_vdf), intent(in) :: rho
  type(t_vdf), intent(inout) :: ef, bf
  real( p_k_part ), intent(in) :: u1
  integer, dimension(:,:,:), intent(in) :: nx_p
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: ri1, ri2, ei1, ei2

  real( p_double ) ::  gamma, x1, x2, cell_vol, dx1, dx2, e1, e2, a
  real( p_double ), pointer, dimension(:,:) :: charge
  integer, dimension(rank) :: lnx0 ! position of local data on global grid

  integer, dimension(rank) :: lb, ub

  integer :: i, ierr, count

  integer( p_int64 ) :: t0, t1

  t0 = timer_ticks()

  ! find position of local data on global grid
  do i = 1, rank
    lnx0( i ) = nx_p( p_lower, i, no_co%my_aid() ) - 1
  enddo

  ! set electric field to 0
  call ef % zero()

  gamma = sqrt( u1**2 + 1.0_p_double )

  dx1 = ef%dx_(1) * gamma
  dx2 = ef%dx_(2)
  a = dx1/2

  ! loop over every nodes
  do i = 1, no_num(no_co)

    lb(1) = nx_p( p_lower, 1, i )
    lb(2) = nx_p( p_lower, 2, i )

    ub(1) = nx_p( p_upper, 1, i )
    ub(2) = nx_p( p_upper, 2, i )

    count =  (ub(1)-lb(1)+1) * (ub(2)-lb(2)+1)

    call alloc( charge, lb, ub )

    ! copy local charge of node being processed to all nodes
    if ( no_co%my_aid() == i ) then
      do ri2 = 1, rho%nx_(2)
        do ri1 = 1, rho%nx_(1)
          charge( lnx0(1) + ri1, lnx0(2) + ri2 ) = rho%f2(1,ri1,ri2)
        enddo
      enddo
    endif

    if ( no_num(no_co) > 1 ) then
      call MPI_BCAST( charge, count, MPI_DOUBLE_PRECISION, i-1, comm(no_co), ierr)
    endif

    ! calculate field generated by current data
    do ri2 = nx_p( p_lower, 2, i ), nx_p( p_upper, 2, i )
      do ri1 = nx_p( p_lower, 1, i ), nx_p( p_upper, 1, i )
         if ( charge(ri1,ri2) /= 0 ) then
            do ei2 = 1, ef%nx_(2)
              x2 = real(lnx0(2)+ei2-ri2, p_double ) * dx2

              do ei1 = 1, ef%nx_(1)
                 x1 = real(lnx0(1)+ei1-ri1, p_double ) * dx1

                 if ( x1 == 0.0 ) then
                   e1 = 0
                 else
                   e1 = 0.5_p_double * log( ( (a + x1)**2 + x2**2 ) / ((a - x1)**2 + x2**2 ) )
                 endif

                 if ( x2 == 0.0 ) then
                   e2 = 0
                 else
                   if ( x1 == 0.0 ) then
                     e2 = 2 * atan( a / x2 )
                   else
                     e2 = atan( ( 2* a * x2 ) / ( x1**2 + x2**2 - a**2 ) )
                   endif
                 endif

                 ef%f2(1, ei1, ei2) = ef%f2(1, ei1, ei2) + e1 * charge(ri1,ri2)
                 ef%f2(2, ei1, ei2) = ef%f2(2, ei1, ei2) + e2 * charge(ri1,ri2)
              enddo
            enddo
         endif
      enddo
    enddo

    ! free temp memory
    call freemem( charge )

  enddo

  ! Normalize results to cell volume (charge is charge density) and boost values to lab frame

  ! This is the same as cell_vol = ef%dx(2)
  cell_vol = ef%dx_(1) * ef%dx_(2) / (2 * a) / (2 * 3.141592653589793d0)
  do ei2 = 1, ef%nx_(2)
    do ei1 = 1, ef%nx_(1)

       ! E2 in the beam ref. frame
       e2 = ef%f2(2, ei1, ei2) * cell_vol

       ef%f2(1, ei1, ei2) = ef%f2(1, ei1, ei2) * cell_vol
       ef%f2(2, ei1, ei2) = e2 * gamma
       ef%f2(3, ei1, ei2) = 0

       bf%f2(1, ei1, ei2) = 0
       bf%f2(2, ei1, ei2) = 0

       ! u1 = \gamma v1
       bf%f2(3, ei1, ei2) = u1 * e2
    enddo
  enddo

  ! Fix values in guard cells
  call update_boundary( ef, p_vdf_replace, no_co, send_msg, recv_msg )
  call update_boundary( bf, p_vdf_replace, no_co, send_msg, recv_msg )

  ! Interpolate values to proper grid positions
  do ei2 = 1, ef%nx_(2)
    do ei1 = 1, ef%nx_(1)
       ef%f2(1, ei1, ei2) = 0.5_p_double * ( ef%f2(1, ei1, ei2) + ef%f2(1, ei1+1, ei2) )
       ef%f2(2, ei1, ei2) = 0.5_p_double * ( ef%f2(2, ei1, ei2) + ef%f2(2, ei1, ei2+1) )

       bf%f2(3, ei1, ei2) = 0.25_p_double * ( bf%f2(3, ei1, ei2  ) + bf%f2(3, ei1+1, ei2  ) + &
                                              bf%f2(3, ei1, ei2+1) + bf%f2(3, ei1+1, ei2+1) )
    enddo
  enddo

  ! Fix values in guard cells
  call update_boundary( ef, p_vdf_replace, no_co, send_msg, recv_msg )
  call update_boundary( bf, p_vdf_replace, no_co, send_msg, recv_msg )

  t1 = timer_ticks()

  if ( root(no_co) ) then
    print *, '(*info*) Time for es_solver_2d_beam : ', timer_interval_seconds( t0, t1 )
  endif

end subroutine es_solver_2d_beam_full
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! This solver assumes that the beam has a sufficiently high gamma such that each cell can
! be considered an infinite plane and that only transverse fields are relevant.
!-----------------------------------------------------------------------------------------
subroutine es_solver_2d_beam( rho, ef, bf, u1, nx_p, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

  use m_vdf_comm

  implicit none

  integer, parameter :: rank = 2

  type(t_vdf), intent(in) :: rho
  type(t_vdf), intent(inout) :: ef, bf
  real( p_k_part ), intent(in) :: u1
  integer, dimension(:,:,:), intent(in) :: nx_p
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: ri1, ri2, ei1, ei2

  real( p_double ) ::  gamma, dx1, dx2, charge_norm
  real( p_double ), pointer, dimension(:,:) :: charge

  integer, dimension(rank) :: lnx0 ! position of local data on global grid
  integer, dimension(rank) :: lb, ub

  integer :: i, ierr, count

  integer( p_int64 ) :: t0, t1

  integer :: color, key, myrank, newcomm, newgroup, newsize, oldgroup, k

  ! This is needed for compatibility with MPI F90+ interface
  integer, dimension(1) :: ranks1, ranks2

  t0 = timer_ticks()

  call ef % zero()

  gamma = sqrt( u1**2 + 1.0_p_double )
  dx1 = ef%dx_(1) * gamma
  dx2 = ef%dx_(2)

  ! field for an infinite plane is E = \sigma / 2
  charge_norm = dx2 / gamma / 2

  ! find position of local data on global grid
  do i = 1, rank
    lnx0( i ) = nx_p( p_lower, i, no_co%my_aid() ) - 1
  enddo

  if ( no_num(no_co) == 1 ) then

    ! serial code
    do ri2 = 1, rho%nx_(2)
      do ri1 = 1, rho%nx_(1)
         if ( rho%f2(1,ri1,ri2) /= 0 ) then
           ei1 = ri1
           do ei2 = lbound(ef%f2,3), ubound(ef%f2,3)
             ef%f2(2,ei1,ei2) = ef%f2(2,ei1,ei2) + rho%f2(1,ri1,ri2) * &
                                        sign( charge_norm, 0.5_p_double + (ei2 - ri2) )
           enddo
         endif
      enddo
    enddo

  else

     ! Create new communicators will all nodes that lie in the same x1
     color = my_ngp( no_co, 1 )
     key = 0
     call MPI_COMM_SPLIT( comm( no_co ), color, key, newcomm, ierr )
     call MPI_COMM_RANK( newcomm, myrank, ierr )
     call MPI_COMM_SIZE( newcomm, newsize, ierr )

     ! Get groups associated with each communicator to allow rank translation
     call MPI_COMM_GROUP( comm(no_co), oldgroup, ierr )
     call MPI_COMM_GROUP( newcomm, newgroup, ierr )

     ! loop over all nodes in new partition
     do i = 1, newsize

       ! Find rank of process in old (global) group
       ranks1(1) = i-1
       call MPI_GROUP_TRANSLATE_RANKS( newgroup, 1, ranks1, oldgroup, ranks2, ierr )
       k = ranks2(1)+1

       ! find position of charge grid on global grid
       lb(1) = nx_p( p_lower, 1, k )
       lb(2) = nx_p( p_lower, 2, k )

       ub(1) = nx_p( p_upper, 1, k )
       ub(2) = nx_p( p_upper, 2, k )

       ! copy local charge of node being processed to all nodes
       call alloc( charge, lb, ub )
       if ( myrank + 1 == i ) then
         do ri2 = 1, rho%nx_(2)
           do ri1 = 1, rho%nx_(1)
             charge( lnx0(1) + ri1, lnx0(2) + ri2 ) = rho%f2(1,ri1,ri2)
           enddo
         enddo
       endif

       count =  (ub(1)-lb(1)+1) * (ub(2)-lb(2)+1)
       call MPI_BCAST( charge, count, MPI_DOUBLE_PRECISION, i-1, newcomm, ierr)

       ! calculate field generated by current data
       do ri2 = lb(2), ub(2)
         do ri1 = lb(1), ub(1)
            if ( charge(ri1,ri2) /= 0 ) then
              ei1 = ri1 - lnx0(1)

              do ei2 = lbound(ef%f2,3), ubound(ef%f2,3)
                ef%f2(2,ei1,ei2) = ef%f2(2,ei1,ei2) + charge( ri1, ri2 ) * &
                                sign( charge_norm, 0.5_p_double + (lnx0(2) + ei2 - ri2) )
              enddo
            endif
         enddo
       enddo

       ! free temp memory
       call freemem( charge )

     enddo

     ! Free new communicators
     call MPI_COMM_FREE( newcomm, ierr )

  endif

  call update_boundary( ef, p_vdf_replace, no_co, send_msg, recv_msg )

  ! boost values and grid center B3
  do ei2 = lbound(ef%f2,3), ubound(ef%f2,3)
    do ei1 = 1, ef%nx_(1)

       bf%f2(1, ei1, ei2) = 0
       bf%f2(2, ei1, ei2) = 0
       ! u1 = \gamma v1
       bf%f2(3, ei1, ei2) = 0.5_p_double * u1 * ( ef%f2(2, ei1, ei2) + ef%f2(2, ei1+1, ei2) )

    enddo
  enddo

  call update_boundary( bf, p_vdf_replace, no_co, send_msg, recv_msg )

  do ei2 =     lbound(ef%f2,3), ubound(ef%f2,3)
    do ei1 = lbound(ef%f2,2), ubound(ef%f2,2)
       ef%f2(1, ei1, ei2) = 0
       ef%f2(2, ei1, ei2) = ef%f2(2, ei1, ei2) * gamma
       ef%f2(3, ei1, ei2) = 0
    enddo
  enddo

  t1 = timer_ticks()

  if ( root(no_co) ) then
    print *, '(*info*) Time for es_solver_2d_beam : ', timer_interval_seconds( t0, t1 )
  endif

end subroutine es_solver_2d_beam
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cylindrical geometry
!-----------------------------------------------------------------------------------------
subroutine es_solver_2dcyl_beam( rho, ef, bf, u1, nx_p, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

  use m_vdf_comm

  implicit none

  integer, parameter :: rank = 2

  type(t_vdf), intent(in) :: rho
  type(t_vdf), intent(inout) :: ef, bf
  real( p_k_part ), intent(in) :: u1
  integer, dimension(:,:,:), intent(in) :: nx_p
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: ri1, ri2, ei1, ei2

  real( p_double ) ::  gamma, dx1, dx2, charge_norm
  real( p_double ), pointer, dimension(:,:) :: charge
  integer, dimension(rank) :: lnx0 ! position of local data on global grid

  integer, dimension(rank) :: lb, ub

  integer :: i, ierr, count

  integer( p_int64 ) :: t0, t1
  integer :: color, key, myrank, newcomm, newgroup, newsize, oldgroup, k
  integer :: ei2_start

  integer, dimension(1) :: ranks1, ranks2

  t0 = timer_ticks()

  call ef%zero()

  gamma = sqrt( u1**2 + 1.0_p_double )
  dx1 = ef%dx_(1) * gamma
  dx2 = ef%dx_(2)

  ! field for an infinite cylinder is Er = \lambda / r
  charge_norm = rho%dx_(2) / gamma

  ! find position of local data on global grid
  do i = 1, rank
    lnx0( i ) = nx_p( p_lower, i, no_co%my_aid() ) - 1
  enddo

  if ( no_num(no_co) == 1 ) then

    ! In cylindrical geometry start in cell i2 = 2 ( cell i2 = 1 is below the axis )
    do ri2 = 2, rho%nx_(2)
      do ri1 = 1, rho%nx_(1)
         if ( rho%f2(1,ri1,ri2) /= 0 ) then
           ei1 = ri1
           ! do ei2 = ri2, ef%nx(2)
           do ei2 = ri2, ubound(ef%f2,3)
             ef%f2(2,ei1,ei2) = ef%f2(2,ei1,ei2) + rho%f2(1,ri1,ri2) * charge_norm  / ( (ei2-1) * dx2 )
           enddo
         endif
      enddo
    enddo

  else

    ! Create new communicators will all nodes that lie in the same x1
     color = my_ngp( no_co, 1 )
     key = 0
     call MPI_COMM_SPLIT( comm( no_co ), color, key, newcomm, ierr )
     call MPI_COMM_RANK( newcomm, myrank, ierr )
     call MPI_COMM_SIZE( newcomm, newsize, ierr )

     ! Get groups associated with each communicator to allow rank translation
     call MPI_COMM_GROUP( comm(no_co), oldgroup, ierr )
     call MPI_COMM_GROUP( newcomm, newgroup, ierr )

     ! loop over all nodes in new partition
     do i = 1, newsize

       ! Find rank of process in old (global) group
       ranks1(1) = i-1
       call MPI_GROUP_TRANSLATE_RANKS( newgroup, 1, ranks1, oldgroup, ranks2, ierr )
       k = ranks2(1)+1

       ! find position of charge grid on global grid
       lb(1) = nx_p( p_lower, 1, k )
       lb(2) = nx_p( p_lower, 2, k )

       ub(1) = nx_p( p_upper, 1, k )
       ub(2) = nx_p( p_upper, 2, k )

       ! copy local charge of node being processed to all nodes
       call alloc( charge, lb, ub )
       if ( myrank + 1 == i ) then
         do ri2 = 1, rho%nx_(2)
           do ri1 = 1, rho%nx_(1)
             charge( lnx0(1) + ri1, lnx0(2) + ri2 ) = rho%f2(1,ri1,ri2)
           enddo
         enddo
       endif

       count =  (ub(1)-lb(1)+1) * (ub(2)-lb(2)+1)
       call MPI_BCAST( charge, count, MPI_DOUBLE_PRECISION, i-1, newcomm, ierr)

       ! don't process charge below the axis
       if ( lb(2) == 1 ) lb(2) = 2

       ! calculate field generated by current data
       do ri2 = lb(2), ub(2)
         do ri1 = lb(1), ub(1)
            if ( charge(ri1,ri2) /= 0 ) then
              ei1 = ri1 - lnx0(1)
              ei2_start = ri2 - lnx0(2)
              if ( ei2_start < 1 ) ei2_start = 1
              do ei2 = ei2_start, ubound(ef%f2,3)
                 ef%f2(2,ei1,ei2) = ef%f2(2,ei1,ei2) + charge(ri1,ri2) * charge_norm  / ( (lnx0(2) + ei2-1) * dx2 )
              enddo
            endif
         enddo
       enddo

       ! free temp memory
       call freemem( charge )

     enddo

     ! Free new communicators
     call MPI_COMM_FREE( newcomm, ierr )

  endif

  call update_boundary( ef, p_vdf_replace, no_co, send_msg, recv_msg )

  ! boost values
!  do ei2 = 1, ef%nx(2)
  do ei2 = lbound(ef%f2,3), ubound(ef%f2,3)
    do ei1 = 1, ef%nx_(1)

       bf%f2(1, ei1, ei2) = 0
       bf%f2(2, ei1, ei2) = 0
       ! u1 = \gamma v1
       bf%f2(3, ei1, ei2) = 0.5_p_double * u1 * ( ef%f2(2, ei1, ei2) + ef%f2(2, ei1+1, ei2) )

       ef%f2(1, ei1, ei2) = 0
       ef%f2(2, ei1, ei2) = ef%f2(2, ei1, ei2) * gamma
       ef%f2(3, ei1, ei2) = 0
    enddo
  enddo

  ! correct values on axial boundary
  if ( lnx0(2) == 0 ) then
    do ei1 = 1, ef%nx_(1)
      ef%f2(2, ei1, 1) = 0
      ef%f2(2, ei1, 0) = -ef%f2(2, ei1, 2)
    enddo
  endif

  ! Fix values in guard cells
  call update_boundary( ef, p_vdf_replace, no_co, send_msg, recv_msg )
  call update_boundary( bf, p_vdf_replace, no_co, send_msg, recv_msg )

  t1 = timer_ticks()

  if ( root(no_co) ) then
    print *, '(*info*) Time for es_solver_2dcyl_beam : ', timer_interval_seconds( t0, t1 )
  endif

end subroutine es_solver_2dcyl_beam
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! 3D cartesian
!-----------------------------------------------------------------------------------------
subroutine es_solver_3d_beam( rho, ef, bf, u1, nx_p, no_co, send_msg, recv_msg )

  use m_vdf_comm

  implicit none

  integer, parameter :: rank = 3

  type(t_vdf), intent(in) :: rho
  type(t_vdf), intent(inout) :: ef, bf
  real( p_k_part ), intent(in) :: u1
  integer, dimension(:,:,:), intent(in) :: nx_p
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: ri1, ri2, ri3, ei1, ei2, ei3

  real( p_double ) :: dx2c, dx2s, dx2s_h, dx3c, dx3s, dx3s_h, r_2
  real( p_double ) ::  gamma, dx1, dx2, dx3, charge_norm

  real( p_double ), pointer, dimension(:,:,:) :: charge
  integer, dimension(rank) :: lnx0 ! position of local data on global grid
  integer, dimension(rank) :: lb, ub

  integer :: color, key, newcomm, myrank, newsize, k, oldgroup, newgroup

  integer :: i, ierr, count
  integer, dimension(1) :: ranks1, ranks2

  real( p_double ), parameter :: M_PI    = 3.14159265358979323846264338327950288d0

  integer( p_int64 ) :: t0, t1

  t0 = timer_ticks()

  call ef%zero()

  ! evaluate gamma for this species
  gamma=sqrt( 1.0d0 + u1**2 )
  dx1 = ef%dx_(1) * gamma
  dx2 = ef%dx_(2)
  dx3 = ef%dx_(3)

  ! Field for an infinit wire is \lambda / (2 \pi r)
  charge_norm = dx2 * dx3 / gamma / 2 / M_PI

  if ( no_num(no_co) == 1 ) then

    ! calculate E1 and E2 field generated by current data
      do ri3 = 1, rho%nx_(3)
          do ri2 = 1, rho%nx_(2)
          do ri1 = 1, rho%nx_(1)
            if ( rho%f3(1,ri1,ri2,ri3) /= 0 ) then

              ei1 = ri1

                do ei3 = lbound( ef%f3 , 4 ), ubound( ef%f3, 4 )

                  dx3c = real( ei3-ri3, p_double )
                  dx3s = dx3c * dx3
                  dx3s_h = (dx3c + 0.5_p_double) * dx3

                   do ei2 = lbound( ef%f3 , 3 ), ubound( ef%f3 , 3 )

                        dx2c = real(ei2-ri2, p_double )
                        dx2s = dx2c * dx2
                          dx2s_h = (dx2c + 0.5_p_double) * dx2

                          r_2 = dx2s_h**2 + dx3s**2
                          ef%f3(2,ei1,ei2,ei3) = ef%f3(2,ei1,ei2,ei3) + rho%f3(1,ri1,ri2,ri3) * charge_norm * dx2s_h / r_2

                          r_2 = dx2s**2 + dx3s_h**2
                        ef%f3(3,ei1,ei2,ei3) = ef%f3(3,ei1,ei2,ei3) + rho%f3(1,ri1,ri2,ri3) * charge_norm * dx3s_h / r_2
                    enddo
                   enddo
            endif
          enddo
        enddo
      enddo

  else

    ! find position of local data on global grid
    lnx0( 1 ) = nx_p( p_lower, 1, no_co%my_aid() ) - 1
    lnx0( 2 ) = nx_p( p_lower, 2, no_co%my_aid() ) - 1
    lnx0( 3 ) = nx_p( p_lower, 3, no_co%my_aid() ) - 1

    ! Create new communicators will all nodes that lie in the same x1
    color = my_ngp( no_co, 1 )
    key = 0
    call MPI_COMM_SPLIT( comm( no_co ), color, key, newcomm, ierr )
    call MPI_COMM_RANK( newcomm, myrank, ierr )
    call MPI_COMM_SIZE( newcomm, newsize, ierr )

    ! Get groups associated with each communicator to allow rank translation
    call MPI_COMM_GROUP( comm(no_co), oldgroup, ierr )
    call MPI_COMM_GROUP( newcomm, newgroup, ierr )

    ! loop over every nodes
    do i = 1, newsize

      ranks1(1) = i-1
      call MPI_GROUP_TRANSLATE_RANKS( newgroup, 1, ranks1, oldgroup, ranks2, ierr )
      k = ranks2(1)+1

      lb(1) = nx_p( p_lower, 1, k )
      lb(2) = nx_p( p_lower, 2, k )
      lb(3) = nx_p( p_lower, 3, k )

      ub(1) = nx_p( p_upper, 1, k )
      ub(2) = nx_p( p_upper, 2, k )
      ub(3) = nx_p( p_upper, 3, k )

      count =  (ub(1)-lb(1)+1) * (ub(2)-lb(2)+1) * (ub(3)-lb(3)+1)

      call alloc( charge, lb, ub )

      ! copy local charge of node being processed to all nodes, CHARGE is the charge density
      if ( no_co%my_aid() == k ) then
          do ri3 = 1, rho%nx_(3)
            do ri2 = 1, rho%nx_(2)
            do ri1 = 1, rho%nx_(1)
                      charge( lnx0(1) + ri1, lnx0(2) + ri2 , lnx0(3) + ri3 ) = rho%f3(1,ri1,ri2,ri3)
            enddo
            enddo
          enddo
      endif

      call MPI_BCAST( charge, count, MPI_DOUBLE_PRECISION, i-1, newcomm, ierr)

      ! calculate E1 and E2 field generated by current data
      do ri3 = nx_p( p_lower, 3, k), nx_p( p_upper, 3, k )
          do ri2 = nx_p( p_lower, 2, k ), nx_p( p_upper, 2, k )
          do ri1 = nx_p( p_lower, 1, k ), nx_p( p_upper, 1, k )
            if ( charge(ri1,ri2,ri3) /= 0 ) then

              ei1=ri1-lnx0(1)

                do ei3 = lbound( ef%f3 , 4 ), ubound( ef%f3, 4 )

                  dx3c = real(lnx0(3)+ei3-ri3, p_double )
                  dx3s = dx3c * dx3
                  dx3s_h = (dx3c + 0.5_p_double) * dx3

                   do ei2 = lbound( ef%f3 , 3 ), ubound( ef%f3 , 3 )

                        dx2c = real(lnx0(2)+ei2-ri2, p_double )
                        dx2s = dx2c * dx2
                          dx2s_h = (dx2c + 0.5_p_double) * dx2

                          r_2 = dx2s_h**2 + dx3s**2
                          ef%f3(2,ei1,ei2,ei3) = ef%f3(2,ei1,ei2,ei3) + charge(ri1,ri2,ri3) * charge_norm * dx2s_h / r_2

                          r_2 = dx2s**2 + dx3s_h**2
                        ef%f3(3,ei1,ei2,ei3) = ef%f3(3,ei1,ei2,ei3) + charge(ri1,ri2,ri3) * charge_norm * dx3s_h / r_2

                    enddo
                   enddo
            endif
          enddo
        enddo
      enddo

      ! free temp memory
      call freemem( charge )
    enddo

     ! Free new communicators
    call MPI_COMM_FREE( newcomm, ierr )

  endif

  ! update E fields in guard cells
  call update_boundary( ef, p_vdf_replace, no_co, send_msg, recv_msg )

  ! boost and grid center B2, B3
  do ei3 = lbound(ef%f3,4), ubound(ef%f3,4)
    do ei2 = lbound(ef%f3,3), ubound(ef%f3,3)
      do ei1 = 1, ef%nx_(1)

         bf%f3(1, ei1, ei2, ei3) = 0
         ! u1 = \gamma v1
         bf%f3(2, ei1, ei2, ei3) = - 0.5_p_double * u1 * (ef%f3(3,ei1,ei2,ei3) + ef%f3(3,ei1+1,ei2,ei3))
         bf%f3(3, ei1, ei2, ei3) =   0.5_p_double * u1 * (ef%f3(2,ei1,ei2,ei3) + ef%f3(2,ei1+1,ei2,ei3))
      enddo
    enddo
  enddo

  ! update B fields in guard cells
  call update_boundary( bf, p_vdf_replace, no_co, send_msg, recv_msg )

  ! Boost E2, E3
  do ei3 = lbound(ef%f3,4), ubound(ef%f3,4)
    do ei2 = lbound(ef%f3,3), ubound(ef%f3,3)
      do ei1 = lbound(ef%f3,2), ubound(ef%f3,2)
         ef%f3(1, ei1, ei2, ei3) = 0
         ef%f3(2, ei1, ei2, ei3) = ef%f3(2, ei1, ei2, ei3) * gamma
         ef%f3(3, ei1, ei2, ei3) = ef%f3(3, ei1, ei2, ei3) * gamma
      enddo
    enddo
  enddo

  t1 = timer_ticks()

  if ( root(no_co) ) then
    print *, '(*info*) Time for 3D beam field calculations was ', timer_interval_seconds( t0, t1 ), ' s'
  endif

end subroutine es_solver_3d_beam
!-----------------------------------------------------------------------------------------

end module m_emf_es_solver
