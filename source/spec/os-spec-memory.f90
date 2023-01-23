module m_species_memory

use m_species_define, only : t_species, t_part_idx

#include "memory/memory.h"

private

integer, private, parameter :: p_stderr = 0
integer, private, parameter :: p_int64  = selected_int_kind(10)

interface alloc
  module procedure alloc_species
  module procedure alloc_1d_species
  module procedure alloc_bound_1d_species
  module procedure alloc_part_idx
  module procedure alloc_1d_part_idx
  module procedure alloc_bound_1d_part_idx
end interface

interface freemem
  module procedure freemem_species
  module procedure freemem_1d_species
  module procedure freemem_part_idx
  module procedure freemem_1d_part_idx
end interface

public :: alloc, freemem

contains

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines for type t_species

#define __TYPE__ class( t_species )
#define __TYPE_STR__ "t_species"
#define FNAME( a )  a ## _species
#define __MAX_DIM__ 1
#include "memory/mem-template.h"

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines for type t_part_idx

#define __TYPE__ type( t_part_idx )
#define __TYPE_STR__ "t_part_idx"
#define FNAME( a )  a ## _part_idx
#define __MAX_DIM__ 1
#include "memory/mem-template.h"

!---------------------------------------------------------------------------------------------------

end module m_species_memory

