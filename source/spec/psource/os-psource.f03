#include "os-config.h"
#include "os-preprocess.fpp"

module m_profile

#include "memory/memory.h"

use m_psource_std
use m_psource_constq
use m_psource_beam
use m_psource_file

use m_species_define

implicit none

private

interface new_source
    module procedure new_source
end interface new_source

public :: new_source

contains

function new_source( type, input_file, coordinates )

  use m_input_file

  implicit none

  class( t_psource ), pointer :: new_source
  character(len=*), intent(in) :: type
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: coordinates
  ! create t_profile object of the selected kind
  select case ( trim(type) )
  case('standard','profile')
    allocate( t_psource_std :: new_source )
  case('constq')
    allocate( t_psource_constq :: new_source )
  case('beamfocus')
    allocate( t_psource_beam :: new_source )
  case('file')
    allocate( t_psource_file :: new_source )
  case default
    new_source => null()
  end select

  ! If successfull read input parameters and finish initalization
  if ( associated(new_source) ) then
    call new_source % read_input( input_file, coordinates )
  endif

end function new_source

end module m_profile