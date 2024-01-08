! Generated automatically by memory.py on  2019-03-06 18:30:41
! Dimensions:  [1, 2, 3, 4]
! Types:  real(p_single), real(p_double), complex(p_single), complex(p_double), integer, integer(p_int64), integer(p_byte)


! Use posix_memalign for allocation instead of the default Fortran allocate routines
!   unless we are compling in Windows.
#define USE_POSIX_MEMALIGN

#if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64) || defined(__KNL__)
  #undef USE_POSIX_MEMALIGN
#endif

! Memory alignment to use with posix_memalign, if not set defaults to 64. Note that this
! has no effect on the default Fortran allocate routines
#ifndef _ALIGNMENT
#define _ALIGNMENT 64
#endif


! Issue a debug message for every allocation / deallocation
!#define DEBUG_MEM

! Maintain a list of all allocated memory blocks, including size and allocation file/line
!#define LOG_MEM

! if the sizeof intrinsic is not supported the routine will log all memory allocations as
! blocks of 1 byte. The total memory will be wrong, but memory leaks will still be properly
! identified

#ifdef __NO_SIZEOF__
#define sizeof( a ) 1
#endif

module memory

implicit none

! kind parameters
integer, private, parameter :: p_single = kind(1.0e0)
integer, private, parameter :: p_double = kind(1.0d0)
integer, private, parameter :: p_byte   = selected_int_kind(2)
integer, private, parameter :: p_int64  = selected_int_kind(10)

! stderr unit
integer, private, parameter :: p_stderr = 0

interface status_mem
  module procedure status_mem_local
  module procedure status_mem_parallel
end interface

! initialize memory system
interface init_mem
  module procedure init_mem
end interface init_mem

! shut down memory system
interface finalize_mem
  module procedure finalize_mem
end interface

interface log_allocation
  module procedure log_allocation
end interface

interface log_deallocation
  module procedure log_deallocation
end interface

#ifdef USE_POSIX_MEMALIGN

interface
  integer(c_int) function posix_memalign_f( memptr, alignment, size ) bind(c, name="posix_memalign")
    use iso_c_binding
    implicit none
    type(c_ptr), intent(out) :: memptr
    integer(c_size_t), intent(in), value :: alignment
    integer(c_size_t), intent(in), value :: size
  end function
end interface

interface
  subroutine free_f( ptr ) bind(c, name="free")
    use iso_c_binding
    implicit none
    type(c_ptr), intent(in), value :: ptr
  end subroutine
end interface

#endif

! Memory accounting
integer(p_int64), private :: total_mem
integer, private :: n_alloc


#ifdef LOG_MEM

type, private :: t_mem_block

  type( t_mem_block ), pointer :: next => null()

  integer( p_int64 ) :: addr, bsize

  character( len = 128 ) :: msg
  character( len = 128 ) :: fname
  integer :: fline

end type

type( t_mem_block ), private,  pointer :: mem_block_list => null(), &
                                mem_block_tail => null()

#endif


! Allocate memory
interface alloc
    module procedure alloc_1d_r4
    module procedure alloc_2d_r4
    module procedure alloc_3d_r4
    module procedure alloc_4d_r4
    module procedure alloc_5d_r4
    module procedure alloc_1d_r8
    module procedure alloc_2d_r8
    module procedure alloc_3d_r8
    module procedure alloc_4d_r8
    module procedure alloc_5d_r8
    module procedure alloc_1d_c4
    module procedure alloc_2d_c4
    module procedure alloc_3d_c4
    module procedure alloc_4d_c4
    module procedure alloc_1d_c8
    module procedure alloc_2d_c8
    module procedure alloc_3d_c8
    module procedure alloc_4d_c8
    module procedure alloc_1d_int4
    module procedure alloc_2d_int4
    module procedure alloc_3d_int4
    module procedure alloc_4d_int4
    module procedure alloc_5d_int4
    module procedure alloc_1d_int8
    module procedure alloc_2d_int8
    module procedure alloc_3d_int8
    module procedure alloc_4d_int8
    module procedure alloc_1d_byte
    module procedure alloc_2d_byte
    module procedure alloc_3d_byte
    module procedure alloc_4d_byte
    module procedure alloc_bound_1d_r4
    module procedure alloc_bound_2d_r4
    module procedure alloc_bound_3d_r4
    module procedure alloc_bound_4d_r4
    module procedure alloc_bound_5d_r4
    module procedure alloc_bound_1d_r8
    module procedure alloc_bound_2d_r8
    module procedure alloc_bound_3d_r8
    module procedure alloc_bound_4d_r8
    module procedure alloc_bound_5d_r8
    module procedure alloc_bound_1d_c4
    module procedure alloc_bound_2d_c4
    module procedure alloc_bound_3d_c4
    module procedure alloc_bound_4d_c4
    module procedure alloc_bound_1d_c8
    module procedure alloc_bound_2d_c8
    module procedure alloc_bound_3d_c8
    module procedure alloc_bound_4d_c8
    module procedure alloc_bound_1d_int4
    module procedure alloc_bound_2d_int4
    module procedure alloc_bound_3d_int4
    module procedure alloc_bound_4d_int4
    module procedure alloc_bound_5d_int4
    module procedure alloc_bound_1d_int8
    module procedure alloc_bound_2d_int8
    module procedure alloc_bound_3d_int8
    module procedure alloc_bound_4d_int8
    module procedure alloc_bound_1d_byte
    module procedure alloc_bound_2d_byte
    module procedure alloc_bound_3d_byte
    module procedure alloc_bound_4d_byte
    module procedure alloc_bound_5d_byte
end interface alloc

! Free memory
interface freemem
    module procedure freemem_1d_r4
    module procedure freemem_2d_r4
    module procedure freemem_3d_r4
    module procedure freemem_4d_r4
    module procedure freemem_5d_r4
    module procedure freemem_1d_r8
    module procedure freemem_2d_r8
    module procedure freemem_3d_r8
    module procedure freemem_4d_r8
    module procedure freemem_5d_r8
    module procedure freemem_1d_c4
    module procedure freemem_2d_c4
    module procedure freemem_3d_c4
    module procedure freemem_4d_c4
    module procedure freemem_5d_c4
    module procedure freemem_1d_c8
    module procedure freemem_2d_c8
    module procedure freemem_3d_c8
    module procedure freemem_4d_c8
    module procedure freemem_5d_c8
    module procedure freemem_1d_int4
    module procedure freemem_2d_int4
    module procedure freemem_3d_int4
    module procedure freemem_4d_int4
    module procedure freemem_5d_int4
    module procedure freemem_1d_int8
    module procedure freemem_2d_int8
    module procedure freemem_3d_int8
    module procedure freemem_4d_int8
    module procedure freemem_1d_byte
    module procedure freemem_2d_byte
    module procedure freemem_3d_byte
    module procedure freemem_4d_byte
    module procedure freemem_5d_byte
end interface freemem


!---------------------------------------------------------------------------------------------------

contains

function strmem( mem_size )

  implicit none

  integer(p_int64), intent(in) :: mem_size
  character(len = 64) :: strmem

  if ( mem_size > 1073741824 ) then
    write( strmem, "(F0.1,A3)" ) real(mem_size)/1073741824., "GB"
  else if ( mem_size > 1048576 ) then
    write( strmem, "(F0.1,A3)" ) real(mem_size)/1048576., "MB"
  else if ( mem_size > 1024 ) then
    write( strmem, "(F0.1,A3)" ) real(mem_size)/1024.,"kB"
  else
    if ( mem_size == 1 ) then
      write( strmem, "(I0,A5)" ) mem_size,"byte"
    else
      write( strmem, "(I0,A6)" ) mem_size,"bytes"
    endif
  endif

  strmem = adjustl( strmem )

end function strmem


function strfileline( file, line )

  implicit none

  character(len = *), intent(in) :: file
  integer, intent(in) :: line
  character(len = len_trim(file) + 15 ) :: strfileline

  write( strfileline, '(A,A,I0)' ) trim(file), ':', line

end function strfileline

function strrange( dim, a, b )

  implicit none
  integer, intent(in) :: dim
  integer, intent(in), dimension(:) :: a
  integer, intent(in), dimension(:), optional :: b

  character( len = 128 ) :: strrange

  integer :: i

  if ( present(b) ) then
    write( strrange, '(I0,A,I0)' ) a(1), ':', b(1)
    do i = 2, dim
    write( strrange, '(A,A,I0,A,I0)' ) trim(strrange), ',', a(i), ':', b(i)
    enddo
  else
    write( strrange, '(I0)' ) a(1)
    do i = 2, dim
    write( strrange, '(A,A,I0)' ) trim(strrange), ',', a(i)
    enddo
  endif

end function strrange

function strblocks( n )

  implicit none
  integer, intent(in) :: n
  character(len=32) :: strblocks

  if ( n == 1 ) then
    strblocks = '1 block'
  else
    write( strblocks, '(I0,A)' ) n, ' blocks'
  endif

end function strblocks

subroutine log_allocation( addr, msg, s, file, line )

  implicit none

  integer( p_int64 ), intent(in) :: addr, s
  character(len=*), intent(in) :: msg

  character(len=*), intent(in) :: file
  integer, intent(in) :: line

#ifdef LOG_MEM
  type( t_mem_block ), pointer :: mem_block
#endif

  !$omp atomic
  total_mem = total_mem + s
  !$omp atomic
  n_alloc = n_alloc + 1

#ifdef LOG_MEM
  allocate( mem_block )
  mem_block%next    => null()
  mem_block%addr    = addr
  mem_block%bsize   = s
  mem_block%msg     = msg
  mem_block%fname   = file
  mem_block%fline   = line

  ! use critical for tiles module, which may allocate in parallel
  !$omp critical (log_mem)
  if ( .not. associated( mem_block_list ) ) then
    mem_block_list => mem_block
  else
    mem_block_tail%next => mem_block
  endif
  mem_block_tail => mem_block
  !$omp end critical (log_mem)
#endif

#ifdef DEBUG_MEM
  write( p_stderr, '(A,Z0)' ) '(*debug*) Allocated ' // trim(msg) // ', ' // &
           trim(strmem(s)) // ' at ' // &
           trim(strfileline(file,line)) // ' addr : 0x', addr
  write( p_stderr, '(A)'  ) '(*debug*) Total allocated memory: '// &
           trim(strblocks(n_alloc)) // ', ' // &
           trim(strmem( total_mem ))

#endif

end subroutine log_allocation

subroutine log_deallocation( addr, msg, s, file, line )

  implicit none

  integer( p_int64 ), intent(in) :: addr, s
  character(len=*), intent(in) :: msg

  character(len=*), intent(in) :: file
  integer, intent(in) :: line

#ifdef LOG_MEM
  type( t_mem_block ), pointer :: mem_block, prev
  integer :: ierr
#endif

  !$omp atomic
  total_mem = total_mem - s
  !$omp atomic
  n_alloc = n_alloc - 1

#ifdef LOG_MEM
  ! use critical for tiles module, which may deallocate in parallel
  !$omp critical (log_mem)
  ! Search for previous allocation
  prev => null()
  mem_block => mem_block_list
  ierr = -1
  do
    if ( .not. associated( mem_block ) ) exit

    if ( mem_block%addr == addr ) then
      ierr = 0
      exit
    endif
    prev => mem_block
    mem_block => mem_block%next
  enddo

  ! If found delete the entry, otherwise abort
  if ( ierr == 0 ) then
    if ( associated( prev ) ) then
      prev % next => mem_block%next
    else
      mem_block_list => mem_block%next
    endif

    if ( associated( mem_block_tail, mem_block ) ) mem_block_tail => prev

    deallocate( mem_block )
  else
    write( p_stderr, '(A,Z0,A)' )'(*error*) The pointer at 0x', addr, &
              ' was not allocated by the memory module.'
    write( p_stderr, '(A)' ) '(*error*) Deallocation failed ' // &
           trim(msg) // ', ' //&
           trim(strmem(s)) // ' at ' // &
           trim(strfileline(file,line))
    call exit(ierr)
  endif
  !$omp end critical (log_mem)
#endif

#ifdef DEBUG_MEM
  write( p_stderr, '(A,Z0)' ) '(*debug*) Dellocated ' // trim(msg) // ', ' //&
           trim(strmem(s)) // ' at ' // &
           trim(strfileline(file,line))// ' addr : 0x', addr
  write( p_stderr, '(A)'  ) '(*debug*) Remaining allocated memory: '// &
           trim(strblocks(n_alloc)) // ', ' // &
           trim(strmem( total_mem ))
#endif


end subroutine log_deallocation

#ifdef LOG_MEM

subroutine list_blocks_mem()

  implicit none

  character(len=128) :: tmp1, tmp2
  type( t_mem_block ), pointer :: mem_block
  integer :: i

  if ( associated(mem_block_list)) print *, 'Allocated memory blocks:'

  mem_block => mem_block_list
  i = 0
  do
    if ( .not. associated( mem_block ) ) exit
    i = i+1
    tmp1 = strmem(mem_block%bsize)
    tmp2 = strfileline(mem_block%fname,mem_block%fline)
    print '(I4,A,Z16.16,A)', i, ' - 0x',mem_block%addr, &
                ', ' // trim(mem_block%msg) // ', ' //&
                trim(tmp1) // ' at ' // &
                trim(tmp2)

    mem_block => mem_block % next
  enddo

end subroutine list_blocks_mem

subroutine delete_mem_block_list()

  implicit none

  type( t_mem_block ), pointer :: mem_block

  mem_block => mem_block_list
  do
    if ( .not. associated( mem_block_list ) ) exit

    mem_block => mem_block_list%next
    deallocate( mem_block_list )
    mem_block_list => mem_block
  enddo

end subroutine delete_mem_block_list

#endif

subroutine status_mem_local( )

  implicit none

  character(len=64) :: tmp


  print *, 'Dynamic memory status '
  print *, 'Number of allocated memory blocks: ', n_alloc
  tmp = strmem(total_mem)
  print *, 'Total Allocated Memory: ', trim(tmp)

#ifdef LOG_MEM
  ! list allocated memory blocks
  print *, ''
  call list_blocks_mem()
#endif

end subroutine status_mem_local

subroutine status_mem_parallel( comm, vmrss, vmsize )

  use mpi

  implicit none

  integer, intent(in) :: comm
  logical, intent(in), optional :: vmrss
  logical, intent(in), optional :: vmsize

  character(len=64), dimension(2) :: tmp
  integer(p_int64) :: global_total_mem, global_max_mem
  integer :: global_n_alloc, global_max_alloc

  integer :: id, ierr

  if ( comm /= MPI_COMM_NULL ) then
    ! get total / max  allocated blocks
    call MPI_REDUCE( n_alloc, global_n_alloc, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
    call MPI_REDUCE( n_alloc, global_max_alloc, 1, MPI_INTEGER, MPI_MAX, 0, comm, ierr )

    ! get total / max  allocated ram
    call MPI_REDUCE( total_mem, global_total_mem, 1, MPI_INTEGER8, MPI_SUM, 0, comm, ierr)
    call MPI_REDUCE( total_mem, global_max_mem, 1, MPI_INTEGER8, MPI_MAX, 0, comm, ierr )

    ! Report results on rank 0
    call MPI_COMM_RANK( comm, id, ierr )
    if ( id == 0 ) then
      print *, ''
      print "(A)", 'Global dynamic memory status:'
      print "(A,I0,A,I0)", '- Number of allocated memory blocks (total/node max): ', global_n_alloc, " / ", global_max_alloc
      tmp(1) = strmem( global_total_mem )
      tmp(2) = strmem( global_max_mem )
      print "(A,A,A,A)", '- Total Allocated Memory (total/node max): ', trim(tmp(1)), " / ", trim(tmp(2))
    endif

  else
    ! Just report local node info
    print "(A)", 'Dynamic memory status:'
    print "(A,I0)", '- Number of allocated memory blocks: ', n_alloc
    tmp(1) = strmem(total_mem)
    print "(A,A)", 'Total Allocated Memory: ', trim(tmp(1))

  endif

  if ( present(vmrss) .and. vmrss ) then
    call status_vmrss( comm )
  endif

  if ( present(vmsize) .and. vmsize ) then
    call status_vmsize( comm )
  endif

end subroutine status_mem_parallel

subroutine status_vmrss( comm )

  use mpi
  use m_system

  implicit none

  integer, intent(in) :: comm

  integer( p_int64 ) :: vmrss, vmrss_max, vmrss_sum
  character(len=64), dimension(2) :: tmp
  integer :: id, ierr


  call get_vmrss( vmrss )

  if ( comm /= MPI_COMM_NULL ) then

    call MPI_REDUCE( vmrss, vmrss_sum, 1, MPI_INTEGER8, MPI_SUM, 0, comm, ierr)
    call MPI_REDUCE( vmrss, vmrss_max, 1, MPI_INTEGER8, MPI_MAX, 0, comm, ierr )

    tmp(1) = strmem( vmrss_sum )
    tmp(2) = strmem( vmrss_max )

    ! Report results on rank 0
    call MPI_COMM_RANK( comm, id, ierr )
    if ( id == 0 ) then
      print "(A,A,A,A)", '- Total Virtual Resident Memory Set (total/node max): ', trim(tmp(1)), " / ", trim(tmp(2))
    endif

  else

    tmp(1) = strmem(vmrss)
    print "(A,A)", 'Total Virtual Resident Memory Set: ', trim(tmp(1))

  endif


end subroutine status_vmrss

subroutine status_vmsize( comm )

  use mpi
  use m_system

  implicit none

  integer, intent(in) :: comm

  integer( p_int64 ) :: vmsize, vmsize_max, vmsize_sum
  character(len=64), dimension(2) :: tmp
  integer :: id, ierr


  call get_vmsize( vmsize )

  if ( comm /= MPI_COMM_NULL ) then

    call MPI_REDUCE( vmsize, vmsize_sum, 1, MPI_INTEGER8, MPI_SUM, 0, comm, ierr)
    call MPI_REDUCE( vmsize, vmsize_max, 1, MPI_INTEGER8, MPI_MAX, 0, comm, ierr )

    tmp(1) = strmem( vmsize_sum )
    tmp(2) = strmem( vmsize_max )

    ! Report results on rank 0
    call MPI_COMM_RANK( comm, id, ierr )
    if ( id == 0 ) then
      print "(A,A,A,A)", '- Total Virtual Memory Size (total/node max): ', trim(tmp(1)), " / ", trim(tmp(2))
    endif

  else

    tmp(1) = strmem(vmsize)
    print "(A,A)", 'Total Virtual Memory Size: ', trim(tmp(1))

  endif

end subroutine status_vmsize


subroutine init_mem( )

  implicit none

  total_mem = 0
  n_alloc   = 0

end subroutine init_mem

subroutine finalize_mem( )

  implicit none

  character(len=64) :: tmp

  if ( n_alloc /= 0 .or. total_mem /= 0 ) then

    write( p_stderr, * ) 'Allocated dynamic memory remaining:'
    write( p_stderr, * ) 'Number of allocated memory blocks: ', n_alloc
    tmp = strmem(total_mem)
    write( p_stderr, * ) 'Total Allocated Memory: ', trim(tmp)

  endif

#ifdef LOG_MEM
  if ( associated( mem_block_list ) ) then
     call list_blocks_mem()
     call delete_mem_block_list()
  endif
#endif

end subroutine finalize_mem

!---------------------------------------------------------------------------------------------------

subroutine alloc_1d_r4( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_single), dimension(:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 1
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'single' // '(' // &
             trim(strrange(1,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 1
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'single' // '(' // &
         trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single' // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_r4


subroutine alloc_bound_1d_r4( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_single), dimension(:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_single), dimension(:), pointer  :: ptmp
  real(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 1
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'single' // '(' // &
             trim(strrange(1,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 1
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single' // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single' // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_r4


subroutine freemem_1d_r4( p, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_single), dimension(:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'single' // '(:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'single' // '(' // &
        trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_1d_r4

!---------------------------------------------------------------------------------------------------

subroutine alloc_2d_r4( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_single), dimension(:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(2) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 2
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'single' // '(' // &
             trim(strrange(2,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 2
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'single' // '(' // &
         trim(strrange(2,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single' // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_2d_r4


subroutine alloc_bound_2d_r4( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_single), dimension(:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_single), dimension(:), pointer  :: ptmp
  real(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 2
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'single' // '(' // &
             trim(strrange(2,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 2
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single' // '(' // &
           trim(strrange(2,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single' // '('// trim(strrange(2,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_2d_r4


subroutine freemem_2d_r4( p, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_single), dimension(:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'single' // '(:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'single' // '(' // &
        trim(strrange(2,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_2d_r4

!---------------------------------------------------------------------------------------------------

subroutine alloc_3d_r4( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_single), dimension(:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(3) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 3
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'single' // '(' // &
             trim(strrange(3,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 3
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'single' // '(' // &
         trim(strrange(3,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single' // '('// trim(strrange(3,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_3d_r4


subroutine alloc_bound_3d_r4( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_single), dimension(:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_single), dimension(:), pointer  :: ptmp
  real(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 3
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'single' // '(' // &
             trim(strrange(3,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 3
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single' // '(' // &
           trim(strrange(3,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single' // '('// trim(strrange(3,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_3d_r4


subroutine freemem_3d_r4( p, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_single), dimension(:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'single' // '(:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'single' // '(' // &
        trim(strrange(3,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_3d_r4

!---------------------------------------------------------------------------------------------------

subroutine alloc_4d_r4( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_single), dimension(:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(4) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 4
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'single' // '(' // &
             trim(strrange(4,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 4
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3),n(4) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'single' // '(' // &
         trim(strrange(4,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single' // '('// trim(strrange(4,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_4d_r4


subroutine alloc_bound_4d_r4( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_single), dimension(:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_single), dimension(:), pointer  :: ptmp
  real(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 4
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'single' // '(' // &
             trim(strrange(4,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 4
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single' // '(' // &
           trim(strrange(4,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single' // '('// trim(strrange(4,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_4d_r4


subroutine freemem_4d_r4( p, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_single), dimension(:,:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'single' // '(:,:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3),lbound(p,4))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'single' // '(' // &
        trim(strrange(4,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_4d_r4

!---------------------------------------------------------------------------------------------------

subroutine alloc_5d_r4( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_single), dimension(:,:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(5) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 5
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'single' // '(' // &
             trim(strrange(5,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 5
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3),n(4),n(5) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'single' // '(' // &
         trim(strrange(5,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single' // '('// trim(strrange(5,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_5d_r4


subroutine alloc_bound_5d_r4( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_single), dimension(:,:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_single), dimension(:), pointer  :: ptmp
  real(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 5
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'single' // '(' // &
             trim(strrange(5,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 5
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single' // '(' // &
           trim(strrange(5,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single' // '('// trim(strrange(5,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_5d_r4


subroutine freemem_5d_r4( p, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_single), dimension(:,:,:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'single' // '(:,:,:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3),lbound(p,4),lbound(p,5))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'single' // '(' // &
        trim(strrange(5,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_5d_r4

!---------------------------------------------------------------------------------------------------

subroutine alloc_1d_r8( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_double), dimension(:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 1
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'double' // '(' // &
             trim(strrange(1,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 1
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'double' // '(' // &
         trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double' // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_r8


subroutine alloc_bound_1d_r8( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_double), dimension(:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_double), dimension(:), pointer  :: ptmp
  real(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 1
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'double' // '(' // &
             trim(strrange(1,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 1
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double' // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double' // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_r8


subroutine freemem_1d_r8( p, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_double), dimension(:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'double' // '(:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'double' // '(' // &
        trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_1d_r8

!---------------------------------------------------------------------------------------------------

subroutine alloc_2d_r8( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_double), dimension(:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(2) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 2
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'double' // '(' // &
             trim(strrange(2,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 2
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'double' // '(' // &
         trim(strrange(2,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double' // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_2d_r8


subroutine alloc_bound_2d_r8( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_double), dimension(:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_double), dimension(:), pointer  :: ptmp
  real(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 2
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'double' // '(' // &
             trim(strrange(2,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 2
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double' // '(' // &
           trim(strrange(2,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double' // '('// trim(strrange(2,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_2d_r8


subroutine freemem_2d_r8( p, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_double), dimension(:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'double' // '(:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'double' // '(' // &
        trim(strrange(2,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_2d_r8

!---------------------------------------------------------------------------------------------------

subroutine alloc_3d_r8( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_double), dimension(:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(3) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 3
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'double' // '(' // &
             trim(strrange(3,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 3
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'double' // '(' // &
         trim(strrange(3,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double' // '('// trim(strrange(3,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_3d_r8


subroutine alloc_bound_3d_r8( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_double), dimension(:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_double), dimension(:), pointer  :: ptmp
  real(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 3
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'double' // '(' // &
             trim(strrange(3,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 3
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double' // '(' // &
           trim(strrange(3,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double' // '('// trim(strrange(3,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_3d_r8


subroutine freemem_3d_r8( p, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_double), dimension(:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'double' // '(:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'double' // '(' // &
        trim(strrange(3,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_3d_r8

!---------------------------------------------------------------------------------------------------

subroutine alloc_4d_r8( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_double), dimension(:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(4) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 4
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'double' // '(' // &
             trim(strrange(4,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 4
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3),n(4) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'double' // '(' // &
         trim(strrange(4,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double' // '('// trim(strrange(4,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_4d_r8


subroutine alloc_bound_4d_r8( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_double), dimension(:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_double), dimension(:), pointer  :: ptmp
  real(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 4
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'double' // '(' // &
             trim(strrange(4,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 4
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double' // '(' // &
           trim(strrange(4,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double' // '('// trim(strrange(4,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_4d_r8


subroutine freemem_4d_r8( p, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_double), dimension(:,:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'double' // '(:,:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3),lbound(p,4))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'double' // '(' // &
        trim(strrange(4,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_4d_r8

!---------------------------------------------------------------------------------------------------

subroutine alloc_5d_r8( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_double), dimension(:,:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(5) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 5
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'double' // '(' // &
             trim(strrange(5,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 5
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3),n(4),n(5) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'double' // '(' // &
         trim(strrange(5,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double' // '('// trim(strrange(5,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_5d_r8


subroutine alloc_bound_5d_r8( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_double), dimension(:,:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  real(p_double), dimension(:), pointer  :: ptmp
  real(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 5
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'double' // '(' // &
             trim(strrange(5,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 5
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double' // '(' // &
           trim(strrange(5,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double' // '('// trim(strrange(5,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_5d_r8


subroutine freemem_5d_r8( p, file, line, stat )

  use iso_c_binding

  implicit none

  real(p_double), dimension(:,:,:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'double' // '(:,:,:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3),lbound(p,4),lbound(p,5))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'double' // '(' // &
        trim(strrange(5,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_5d_r8

!---------------------------------------------------------------------------------------------------

subroutine alloc_1d_c4( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_single), dimension(:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 1
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'single complex' // '(' // &
             trim(strrange(1,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 1
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'single complex' // '(' // &
         trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_c4


subroutine alloc_bound_1d_c4( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_single), dimension(:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_single), dimension(:), pointer  :: ptmp
  complex(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 1
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'single complex' // '(' // &
             trim(strrange(1,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 1
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single complex' // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_c4


subroutine freemem_1d_c4( p, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_single), dimension(:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'single complex' // '(:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'single complex' // '(' // &
        trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_1d_c4

!---------------------------------------------------------------------------------------------------

subroutine alloc_2d_c4( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_single), dimension(:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(2) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 2
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'single complex' // '(' // &
             trim(strrange(2,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 2
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'single complex' // '(' // &
         trim(strrange(2,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_2d_c4


subroutine alloc_bound_2d_c4( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_single), dimension(:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_single), dimension(:), pointer  :: ptmp
  complex(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 2
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'single complex' // '(' // &
             trim(strrange(2,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 2
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single complex' // '(' // &
           trim(strrange(2,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(2,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_2d_c4


subroutine freemem_2d_c4( p, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_single), dimension(:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'single complex' // '(:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'single complex' // '(' // &
        trim(strrange(2,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_2d_c4

!---------------------------------------------------------------------------------------------------

subroutine alloc_3d_c4( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_single), dimension(:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(3) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 3
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'single complex' // '(' // &
             trim(strrange(3,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 3
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'single complex' // '(' // &
         trim(strrange(3,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(3,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_3d_c4


subroutine alloc_bound_3d_c4( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_single), dimension(:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_single), dimension(:), pointer  :: ptmp
  complex(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 3
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'single complex' // '(' // &
             trim(strrange(3,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 3
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single complex' // '(' // &
           trim(strrange(3,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(3,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_3d_c4


subroutine freemem_3d_c4( p, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_single), dimension(:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'single complex' // '(:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'single complex' // '(' // &
        trim(strrange(3,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_3d_c4

!---------------------------------------------------------------------------------------------------

subroutine alloc_4d_c4( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_single), dimension(:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(4) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 4
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'single complex' // '(' // &
             trim(strrange(4,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 4
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3),n(4) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'single complex' // '(' // &
         trim(strrange(4,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(4,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_4d_c4


subroutine alloc_bound_4d_c4( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_single), dimension(:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_single), dimension(:), pointer  :: ptmp
  complex(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 4
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'single complex' // '(' // &
             trim(strrange(4,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 4
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single complex' // '(' // &
           trim(strrange(4,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(4,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_4d_c4


subroutine freemem_4d_c4( p, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_single), dimension(:,:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'single complex' // '(:,:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3),lbound(p,4))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'single complex' // '(' // &
        trim(strrange(4,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_4d_c4

!---------------------------------------------------------------------------------------------------

subroutine alloc_5d_c4( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_single), dimension(:,:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(5) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 5
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'single complex' // '(' // &
             trim(strrange(5,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 5
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3),n(4),n(5) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'single complex' // '(' // &
         trim(strrange(5,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(5,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_5d_c4


subroutine alloc_bound_5d_c4( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_single), dimension(:,:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_single), dimension(:), pointer  :: ptmp
  complex(p_single) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 5
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'single complex' // '(' // &
             trim(strrange(5,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 5
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'single complex' // '(' // &
           trim(strrange(5,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'single complex' // '('// trim(strrange(5,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_5d_c4


subroutine freemem_5d_c4( p, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_single), dimension(:,:,:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'single complex' // '(:,:,:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3),lbound(p,4),lbound(p,5))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'single complex' // '(' // &
        trim(strrange(5,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_5d_c4

!---------------------------------------------------------------------------------------------------

subroutine alloc_1d_c8( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_double), dimension(:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 1
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'double complex' // '(' // &
             trim(strrange(1,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 1
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'double complex' // '(' // &
         trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_c8


subroutine alloc_bound_1d_c8( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_double), dimension(:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_double), dimension(:), pointer  :: ptmp
  complex(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 1
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'double complex' // '(' // &
             trim(strrange(1,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 1
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double complex' // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_c8


subroutine freemem_1d_c8( p, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_double), dimension(:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'double complex' // '(:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'double complex' // '(' // &
        trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_1d_c8

!---------------------------------------------------------------------------------------------------

subroutine alloc_2d_c8( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_double), dimension(:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(2) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 2
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'double complex' // '(' // &
             trim(strrange(2,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 2
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'double complex' // '(' // &
         trim(strrange(2,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_2d_c8


subroutine alloc_bound_2d_c8( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_double), dimension(:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_double), dimension(:), pointer  :: ptmp
  complex(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 2
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'double complex' // '(' // &
             trim(strrange(2,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 2
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double complex' // '(' // &
           trim(strrange(2,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(2,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_2d_c8


subroutine freemem_2d_c8( p, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_double), dimension(:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'double complex' // '(:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'double complex' // '(' // &
        trim(strrange(2,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_2d_c8

!---------------------------------------------------------------------------------------------------

subroutine alloc_3d_c8( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_double), dimension(:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(3) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 3
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'double complex' // '(' // &
             trim(strrange(3,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 3
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'double complex' // '(' // &
         trim(strrange(3,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(3,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_3d_c8


subroutine alloc_bound_3d_c8( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_double), dimension(:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_double), dimension(:), pointer  :: ptmp
  complex(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 3
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'double complex' // '(' // &
             trim(strrange(3,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 3
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double complex' // '(' // &
           trim(strrange(3,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(3,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_3d_c8


subroutine freemem_3d_c8( p, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_double), dimension(:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'double complex' // '(:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'double complex' // '(' // &
        trim(strrange(3,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_3d_c8

!---------------------------------------------------------------------------------------------------

subroutine alloc_4d_c8( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_double), dimension(:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(4) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 4
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'double complex' // '(' // &
             trim(strrange(4,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 4
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3),n(4) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'double complex' // '(' // &
         trim(strrange(4,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(4,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_4d_c8


subroutine alloc_bound_4d_c8( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_double), dimension(:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_double), dimension(:), pointer  :: ptmp
  complex(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 4
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'double complex' // '(' // &
             trim(strrange(4,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 4
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double complex' // '(' // &
           trim(strrange(4,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(4,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_4d_c8


subroutine freemem_4d_c8( p, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_double), dimension(:,:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'double complex' // '(:,:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3),lbound(p,4))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'double complex' // '(' // &
        trim(strrange(4,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_4d_c8

!---------------------------------------------------------------------------------------------------

subroutine alloc_5d_c8( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_double), dimension(:,:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(5) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 5
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'double complex' // '(' // &
             trim(strrange(5,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 5
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3),n(4),n(5) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'double complex' // '(' // &
         trim(strrange(5,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(5,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_5d_c8


subroutine alloc_bound_5d_c8( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_double), dimension(:,:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  complex(p_double), dimension(:), pointer  :: ptmp
  complex(p_double) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 5
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'double complex' // '(' // &
             trim(strrange(5,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 5
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'double complex' // '(' // &
           trim(strrange(5,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'double complex' // '('// trim(strrange(5,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_5d_c8


subroutine freemem_5d_c8( p, file, line, stat )

  use iso_c_binding

  implicit none

  complex(p_double), dimension(:,:,:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'double complex' // '(:,:,:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3),lbound(p,4),lbound(p,5))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'double complex' // '(' // &
        trim(strrange(5,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_5d_c8

!---------------------------------------------------------------------------------------------------

subroutine alloc_1d_int4( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  integer, dimension(:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 1
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'integer' // '(' // &
             trim(strrange(1,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 1
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'integer' // '(' // &
         trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'integer' // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_int4


subroutine alloc_bound_1d_int4( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  integer, dimension(:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer, dimension(:), pointer  :: ptmp
  integer :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 1
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'integer' // '(' // &
             trim(strrange(1,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 1
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'integer' // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'integer' // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_int4


subroutine freemem_1d_int4( p, file, line, stat )

  use iso_c_binding

  implicit none

  integer, dimension(:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'integer' // '(:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'integer' // '(' // &
        trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_1d_int4

!---------------------------------------------------------------------------------------------------

subroutine alloc_2d_int4( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  integer, dimension(:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(2) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 2
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'integer' // '(' // &
             trim(strrange(2,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 2
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'integer' // '(' // &
         trim(strrange(2,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'integer' // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_2d_int4


subroutine alloc_bound_2d_int4( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  integer, dimension(:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer, dimension(:), pointer  :: ptmp
  integer :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 2
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'integer' // '(' // &
             trim(strrange(2,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 2
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'integer' // '(' // &
           trim(strrange(2,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'integer' // '('// trim(strrange(2,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_2d_int4


subroutine freemem_2d_int4( p, file, line, stat )

  use iso_c_binding

  implicit none

  integer, dimension(:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'integer' // '(:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'integer' // '(' // &
        trim(strrange(2,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_2d_int4

!---------------------------------------------------------------------------------------------------

subroutine alloc_3d_int4( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  integer, dimension(:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(3) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 3
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'integer' // '(' // &
             trim(strrange(3,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 3
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'integer' // '(' // &
         trim(strrange(3,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'integer' // '('// trim(strrange(3,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_3d_int4


subroutine alloc_bound_3d_int4( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  integer, dimension(:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer, dimension(:), pointer  :: ptmp
  integer :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 3
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'integer' // '(' // &
             trim(strrange(3,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 3
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'integer' // '(' // &
           trim(strrange(3,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'integer' // '('// trim(strrange(3,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_3d_int4


subroutine freemem_3d_int4( p, file, line, stat )

  use iso_c_binding

  implicit none

  integer, dimension(:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'integer' // '(:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'integer' // '(' // &
        trim(strrange(3,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_3d_int4

!---------------------------------------------------------------------------------------------------

subroutine alloc_4d_int4( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  integer, dimension(:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(4) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 4
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'integer' // '(' // &
             trim(strrange(4,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 4
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3),n(4) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'integer' // '(' // &
         trim(strrange(4,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'integer' // '('// trim(strrange(4,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_4d_int4


subroutine alloc_bound_4d_int4( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  integer, dimension(:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer, dimension(:), pointer  :: ptmp
  integer :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 4
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'integer' // '(' // &
             trim(strrange(4,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 4
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'integer' // '(' // &
           trim(strrange(4,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'integer' // '('// trim(strrange(4,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_4d_int4


subroutine freemem_4d_int4( p, file, line, stat )

  use iso_c_binding

  implicit none

  integer, dimension(:,:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'integer' // '(:,:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3),lbound(p,4))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'integer' // '(' // &
        trim(strrange(4,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_4d_int4

!---------------------------------------------------------------------------------------------------

subroutine alloc_5d_int4( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  integer, dimension(:,:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(5) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 5
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'integer' // '(' // &
             trim(strrange(5,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 5
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3),n(4),n(5) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'int' // '(' // &
         trim(strrange(5,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'integer' // '('// trim(strrange(5,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_5d_int4


subroutine alloc_bound_5d_int4( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  integer, dimension(:,:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer, dimension(:), pointer  :: ptmp
  integer :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 5
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'integer' // '(' // &
             trim(strrange(5,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 5
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'int' // '(' // &
           trim(strrange(5,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'integer' // '('// trim(strrange(5,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_5d_int4


subroutine freemem_5d_int4( p, file, line, stat )

  use iso_c_binding

  implicit none

  integer, dimension(:,:,:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'integer' // '(:,:,:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3),lbound(p,4),lbound(p,5))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'integer' // '(' // &
        trim(strrange(5,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_5d_int4

!---------------------------------------------------------------------------------------------------

subroutine alloc_1d_int8( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_int64), dimension(:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_int64) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 1
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             '64 bit integer' // '(' // &
             trim(strrange(1,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 1
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         '64 bit integer' // '(' // &
         trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, '64 bit integer' // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_int8


subroutine alloc_bound_1d_int8( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_int64), dimension(:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_int64), dimension(:), pointer  :: ptmp
  integer(p_int64) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 1
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             '64 bit integer' // '(' // &
             trim(strrange(1,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 1
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           '64 bit integer' // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, '64 bit integer' // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_int8


subroutine freemem_1d_int8( p, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_int64), dimension(:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, '64 bit integer' // '(:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            '64 bit integer' // '(' // &
        trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_1d_int8

!---------------------------------------------------------------------------------------------------

subroutine alloc_2d_int8( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_int64), dimension(:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(2) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_int64) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 2
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             '64 bit integer' // '(' // &
             trim(strrange(2,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 2
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         '64 bit integer' // '(' // &
         trim(strrange(2,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, '64 bit integer' // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_2d_int8


subroutine alloc_bound_2d_int8( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_int64), dimension(:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_int64), dimension(:), pointer  :: ptmp
  integer(p_int64) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 2
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             '64 bit integer' // '(' // &
             trim(strrange(2,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 2
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           '64 bit integer' // '(' // &
           trim(strrange(2,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, '64 bit integer' // '('// trim(strrange(2,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_2d_int8


subroutine freemem_2d_int8( p, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_int64), dimension(:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, '64 bit integer' // '(:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            '64 bit integer' // '(' // &
        trim(strrange(2,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_2d_int8

!---------------------------------------------------------------------------------------------------

subroutine alloc_3d_int8( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_int64), dimension(:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(3) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_int64) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 3
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             '64 bit integer' // '(' // &
             trim(strrange(3,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 3
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         '64 bit integer' // '(' // &
         trim(strrange(3,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, '64 bit integer' // '('// trim(strrange(3,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_3d_int8


subroutine alloc_bound_3d_int8( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_int64), dimension(:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_int64), dimension(:), pointer  :: ptmp
  integer(p_int64) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 3
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             '64 bit integer' // '(' // &
             trim(strrange(3,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 3
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           '64 bit integer' // '(' // &
           trim(strrange(3,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, '64 bit integer' // '('// trim(strrange(3,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_3d_int8


subroutine freemem_3d_int8( p, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_int64), dimension(:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, '64 bit integer' // '(:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            '64 bit integer' // '(' // &
        trim(strrange(3,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_3d_int8

!---------------------------------------------------------------------------------------------------

subroutine alloc_4d_int8( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_int64), dimension(:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(4) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_int64) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 4
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             '64 bit integer' // '(' // &
             trim(strrange(4,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 4
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3),n(4) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         '64 bit integer' // '(' // &
         trim(strrange(4,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, '64 bit integer' // '('// trim(strrange(4,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_4d_int8


subroutine alloc_bound_4d_int8( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_int64), dimension(:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_int64), dimension(:), pointer  :: ptmp
  integer(p_int64) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 4
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             '64 bit integer' // '(' // &
             trim(strrange(4,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 4
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           '64 bit integer' // '(' // &
           trim(strrange(4,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, '64 bit integer' // '('// trim(strrange(4,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_4d_int8


subroutine freemem_4d_int8( p, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_int64), dimension(:,:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, '64 bit integer' // '(:,:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3),lbound(p,4))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            '64 bit integer' // '(' // &
        trim(strrange(4,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_4d_int8

!---------------------------------------------------------------------------------------------------

subroutine alloc_1d_byte( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_byte), dimension(:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_byte) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 1
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'byte' // '(' // &
             trim(strrange(1,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 1
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'byte' // '(' // &
         trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'byte' // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_byte


subroutine alloc_bound_1d_byte( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_byte), dimension(:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_byte), dimension(:), pointer  :: ptmp
  integer(p_byte) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 1
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'byte' // '(' // &
             trim(strrange(1,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 1
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'byte' // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'byte' // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_byte


subroutine freemem_1d_byte( p, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_byte), dimension(:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'byte' // '(:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'byte' // '(' // &
        trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_1d_byte

!---------------------------------------------------------------------------------------------------

subroutine alloc_2d_byte( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_byte), dimension(:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(2) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_byte) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 2
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'byte' // '(' // &
             trim(strrange(2,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 2
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'byte' // '(' // &
         trim(strrange(2,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'byte' // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_2d_byte


subroutine alloc_bound_2d_byte( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_byte), dimension(:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_byte), dimension(:), pointer  :: ptmp
  integer(p_byte) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 2
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'byte' // '(' // &
             trim(strrange(2,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 2
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'byte' // '(' // &
           trim(strrange(2,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'byte' // '('// trim(strrange(2,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_2d_byte


subroutine freemem_2d_byte( p, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_byte), dimension(:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'byte' // '(:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'byte' // '(' // &
        trim(strrange(2,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_2d_byte

!---------------------------------------------------------------------------------------------------

subroutine alloc_3d_byte( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_byte), dimension(:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(3) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_byte) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 3
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'byte' // '(' // &
             trim(strrange(3,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 3
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'byte' // '(' // &
         trim(strrange(3,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'byte' // '('// trim(strrange(3,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_3d_byte


subroutine alloc_bound_3d_byte( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_byte), dimension(:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_byte), dimension(:), pointer  :: ptmp
  integer(p_byte) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 3
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'byte' // '(' // &
             trim(strrange(3,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 3
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'byte' // '(' // &
           trim(strrange(3,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'byte' // '('// trim(strrange(3,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_3d_byte


subroutine freemem_3d_byte( p, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_byte), dimension(:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'byte' // '(:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'byte' // '(' // &
        trim(strrange(3,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_3d_byte

!---------------------------------------------------------------------------------------------------

subroutine alloc_4d_byte( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_byte), dimension(:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(4) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_byte) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 4
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'byte' // '(' // &
             trim(strrange(4,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 4
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3),n(4) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'byte' // '(' // &
         trim(strrange(4,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'byte' // '('// trim(strrange(4,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_4d_byte


subroutine alloc_bound_4d_byte( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_byte), dimension(:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_byte), dimension(:), pointer  :: ptmp
  integer(p_byte) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 4
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'byte' // '(' // &
             trim(strrange(4,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 4
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'byte' // '(' // &
           trim(strrange(4,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'byte' // '('// trim(strrange(4,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_4d_byte


subroutine freemem_4d_byte( p, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_byte), dimension(:,:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'byte' // '(:,:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3),lbound(p,4))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'byte' // '(' // &
        trim(strrange(4,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_4d_byte

!---------------------------------------------------------------------------------------------------

subroutine alloc_5d_byte( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_byte), dimension(:,:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: n

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims   = -1001

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(5) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_byte) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 5
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             'byte' // '(' // &
             trim(strrange(5,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory

#ifdef USE_POSIX_MEMALIGN
  cbufsize = c_sizeof( tmp )
  do i = 1, 5
    cbufsize = cbufsize * n(i)
    shape(i) = n(i)
  enddo
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( n(1),n(2),n(3),n(4),n(5) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         'byte' // '(' // &
         trim(strrange(5,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, p, shape )
#endif

  ! log allocation
  addr  = loc(p)
  bsize = sizeof( p )
  call log_allocation(  addr, 'byte' // '('// trim(strrange(5,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_5d_byte


subroutine alloc_bound_5d_byte( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_byte), dimension(:,:,:,:,:), pointer    :: p
  integer, dimension(:), intent(in)     :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002

#ifdef USE_POSIX_MEMALIGN
  integer, dimension(1) :: shape
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  integer(c_size_t) :: cbufsize
  type(c_ptr) :: cbuf
  integer(p_byte), dimension(:), pointer  :: ptmp
  integer(p_byte) :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 5
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             'byte' // '(' // &
             trim(strrange(5,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  shape(1) = 1
  do i = 1, 5
    shape(1) = shape(1) * ( ub(i) - lb(i) + 1 )
  enddo
  cbufsize = c_sizeof( tmp ) * shape(1)
  ierr = posix_memalign_f( cbuf, alignment, cbufsize )
#else
  allocate( p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5) ), stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           'byte' // '(' // &
           trim(strrange(5,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif

#ifdef USE_POSIX_MEMALIGN
  ! Associate with fortran pointer
  call c_f_pointer( cbuf, ptmp, shape )
  p( lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5) ) => ptmp
#endif

  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation(  addr, 'byte' // '('// trim(strrange(5,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_5d_byte


subroutine freemem_5d_byte( p, file, line, stat )

  use iso_c_binding

  implicit none

  integer(p_byte), dimension(:,:,:,:,:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, 'byte' // '(:,:,:,:,:)', &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f(c_loc(p(lbound(p,1),lbound(p,2),lbound(p,3),lbound(p,4),lbound(p,5))))
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            'byte' // '(' // &
        trim(strrange(5,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_5d_byte


end module memory

