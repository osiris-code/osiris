!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     boundary class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_bnd

#include "memory/memory.h"

use m_system
use m_parameters
use m_species_define, only : t_part_idx, t_spec_msg
use m_vdf_comm
use m_species_comm

implicit none

private

! In each t_sim_tiles
type :: t_bnd

  ! indexes of particles crossing lower/upper boundary
  type( t_part_idx ), dimension(2) :: bnd_cross
  ! indexes of particles leaving the node
  type( t_part_idx ) :: node_cross
  ! Replaces module variables for vdf comm
  type( t_vdf_msg ), dimension(2) :: send_vdf, recv_vdf
  ! Replaces module variables for part comm
  type( t_spec_msg ), dimension(2) :: send_spec, recv_spec

  contains

  procedure :: init => init_bnd
  procedure :: cleanup => cleanup_bnd

end type t_bnd

public :: t_bnd

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine init_bnd( this, buffer_size )

  implicit none

  class( t_bnd ), intent(inout) :: this
  integer, intent(in), optional :: buffer_size

  if (present(buffer_size)) then
    call init_spec_comm( this%send_spec, this%recv_spec, buffer_size )
  else
    call init_spec_comm( this%send_spec, this%recv_spec )
  endif

end subroutine init_bnd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine cleanup_bnd( this )

  implicit none

  class( t_bnd ), intent(inout) :: this

  integer :: i

  call freemem( this%node_cross % idx )
  do i = p_lower, p_upper
    call freemem( this%bnd_cross( i ) % idx )
  enddo

  call cleanup( this%send_spec )
  call cleanup( this%recv_spec )
  call cleanup( this%send_vdf )
  call cleanup( this%recv_vdf )


end subroutine cleanup_bnd
!-----------------------------------------------------------------------------------------

end module m_bnd
