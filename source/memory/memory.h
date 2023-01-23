! Include file for the memory module
! Files must include this file using #include "memory/memory.h" rather than just using the module
! the #include must be placed where the module use statement would usually be i.e.
!
! module module2
!
!  use module1
!  #include "memory/memory.h"
!

#ifndef _MEMORY_H
#define _MEMORY_H

#define alloc(...) alloc(__VA_ARGS__,__FILE__,__LINE__)
#define freemem(...)  freemem(__VA_ARGS__,__FILE__,__LINE__)

#endif

use memory
