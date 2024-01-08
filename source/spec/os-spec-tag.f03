#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_tag

#include "memory/memory.h"

use m_species_define

#ifdef __HAS_TRACKS__
use m_species_tracks
#endif

implicit none

private

interface set_tags
  module procedure set_tags_range
  module procedure set_tags_single
end interface

public :: set_tags

contains

!-------------------------------------------------------------------------------
subroutine set_tags_range( this, idx0, idx1 )
!-------------------------------------------------------------------------------
! Sets tags of particle range
!-------------------------------------------------------------------------------

  implicit none

  class( t_species ), intent(inout) :: this
  integer, intent(in) :: idx0, idx1

  integer :: i,tag

  ! version 1 of tags

  tag = this%num_created

  do i = idx0, idx1
    tag = tag + 1
    this%tag(1,i) = this%ngp_id
    this%tag(2,i) = tag
  enddo

#ifdef __HAS_TRACKS__

  ! if the tracking diagnostic is on and tags are in the missing list
  ! move them to the present list and mark the particle as being tracked
  if ( this%diag%ndump_fac_tracks > 0)  then
     call new_particles( this%diag%tracks, this, idx0, idx1 )
  endif

#endif

end subroutine set_tags_range
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine set_tags_single( this, idx )
!-------------------------------------------------------------------------------
! Sets tag of single particle
!-------------------------------------------------------------------------------

  implicit none

  class( t_species ), intent(inout) :: this
  integer, intent(in) :: idx

  integer :: tag

  ! version 1 of tags

  tag = this%num_created + 1
  this%tag(1,idx) = this%ngp_id
  this%tag(2,idx) = tag

#ifdef __HAS_TRACKS__

  ! if the tracking diagnostic is on and tags are in the missing list
  ! move them to the present list and mark the particle as being tracked
  if ( this%diag%ndump_fac_tracks > 0)  then
     call new_particles( this%diag%tracks, this%tag(:,idx), idx )
  endif

#endif

end subroutine set_tags_single
!-------------------------------------------------------------------------------

end module m_species_tag