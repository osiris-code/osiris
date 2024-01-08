#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_initfields

#include "memory/memory.h"

use m_system
use m_parameters

use m_species_define, only: t_species
use m_emf_define,     only: t_emf
use m_grid_define,    only: t_grid
use m_node_conf,      only: t_node_conf
use m_vdf_define,     only: t_vdf
use m_vdf_comm,       only: update_boundary, t_vdf_msg
use m_vdf_math,       only: add
use m_species_charge, only: deposit_rho
use m_emf_es_solver,  only: es_solver, es_solver_beam

use m_vpml,           only: update_beam

implicit none

private

interface init_fields
  module procedure init_fields_species
end interface

public :: init_fields

contains

!-----------------------------------------------------------------------------------------
! Initialize fields from initial charge / momentum distribution
!-----------------------------------------------------------------------------------------
subroutine init_fields_species( this, emf, no_co, grid, send_msg, recv_msg )
  
  implicit none

  class( t_species ), intent(in) :: this
  class( t_emf ), intent(inout) :: emf
  class( t_node_conf ), intent(in) :: no_co
  class( t_grid ), intent(in) :: grid
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  type( t_vdf ) :: rho, ef, bf

  SCR_ROOT( ' - Initializing EM fields for species "', trim(this%name),'"' )

  ! Get charge density from species

  ! Create new charge grid with same dimensions as E field grid and 1 field component
  call rho%new( emf%e, f_dim = 1 )

  ! zero the array
  call rho%zero( )

  ! Deposit the charge density

  call deposit_rho( this, rho )

  ! Get contributions from neighboring parallel nodes
  call update_boundary( rho, p_vdf_add, no_co, send_msg, recv_msg )

  ! Find E field from initial charge distribution

  ! Create temporary E field grid
  call ef%new( emf%e )

  if ( this%udist%ufl(1) == 0.0 ) then

     ! Call ES solver to find E field from charge density
     call es_solver( rho, ef, grid, no_co, send_msg, recv_msg )

     ! Add ES field to global simulation field
     call add( emf%e, ef )

  else

     call bf % new( emf%b )

     ! Call ES solver to find E field from beam charge density
     call es_solver_beam( rho, ef, bf, this%udist%ufl(1), grid, no_co, send_msg, recv_msg )

     ! Add static fields to global simulation field
     call add( emf%e, ef )
     call add( emf%b, bf )

     ! Update PML boundaries, if any
     call update_beam( emf%bnd_con%vpml_all, ef, bf, no_co, send_msg, recv_msg )

     ! cleanup
     call bf % cleanup()

  endif

  ! cleanup
  call ef % cleanup()
  call rho % cleanup()

end subroutine

end module m_species_initfields
