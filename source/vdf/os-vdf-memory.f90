module m_vdf_memory

use m_vdf_define, only : t_vdf

#include "memory/memory.h"

integer, private, parameter :: p_stderr = 0
integer, private, parameter :: p_int64  = selected_int_kind(10)

private

interface alloc
  module procedure alloc_vdf
  module procedure alloc_1d_vdf
  module procedure alloc_bound_1d_vdf
end interface

interface freemem
  module procedure freemem_vdf
  module procedure freemem_1d_vdf
end interface

public :: alloc, freemem

contains

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines

#define __TYPE__ type( t_vdf )
#define __TYPE_STR__ "t_vdf"
#define FNAME( a )  a ## _vdf
#define __MAX_DIM__ 1
#include "memory/mem-template.h"

!---------------------------------------------------------------------------------------------------

end module m_vdf_memory
