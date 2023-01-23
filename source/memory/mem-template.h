! Generated automatically by memory.py on  2019-03-06 18:36:58

#ifndef _MEMORY_H
#error memory.h file must be included in this source file (#include "memory.h")
#endif

#ifndef __TYPE__
#error The macro __TYPE__ must be defined before including this file.
#endif

#ifndef __TYPE_STR__
#error The macro __TYPE_STR__ must be defined before including this file.
#endif

#ifndef FNAME
#error The macro FNAME must be defined before including this file.
#endif

#ifndef __MAX_DIM__
#define __MAX_DIM__ 2
#endif


subroutine FNAME(alloc)( p, file, line, stat )

  use iso_c_binding

  implicit none

  __TYPE__, pointer    :: p

  character(len=*), intent(in) :: file
  integer, intent(in)       :: line
  integer, intent(out), optional     :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

#ifdef USE_POSIX_MEMALIGN
  integer(c_size_t), parameter :: alignment = _ALIGNMENT
  type(c_ptr) :: cbuf
  __TYPE__ :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Allocate memory
#ifdef USE_POSIX_MEMALIGN
  ierr = posix_memalign_f( cbuf, alignment, c_sizeof(tmp) )
#else
  allocate( p, stat = ierr )
#endif

  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         __TYPE_STR__ // ', ' // &
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
  call c_f_pointer( cbuf, p )
#endif

  ! log allocation
  addr  = loc( p )
  bsize = sizeof( p )
  call log_allocation(  addr, __TYPE_STR__ , &
                        bsize, &
                        file, line )

end subroutine FNAME(alloc)

#if __MAX_DIM__ > 0

subroutine FNAME(alloc_1d)( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  __TYPE__, dimension(:), pointer    :: p
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
  __TYPE__ :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 1
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             __TYPE_STR__ // '(' // &
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
         __TYPE_STR__ // '(' // &
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
  call log_allocation(  addr, __TYPE_STR__ // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine FNAME(alloc_1d)

#endif
#if __MAX_DIM__ > 1

subroutine FNAME(alloc_2d)( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  __TYPE__, dimension(:,:), pointer    :: p
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
  __TYPE__ :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 2
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             __TYPE_STR__ // '(' // &
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
         __TYPE_STR__ // '(' // &
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
  call log_allocation(  addr, __TYPE_STR__ // '('// trim(strrange(2,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine FNAME(alloc_2d)

#endif
#if __MAX_DIM__ > 0

subroutine FNAME(alloc_bound_1d)( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  __TYPE__, dimension(:), pointer    :: p
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
  __TYPE__, dimension(:), pointer  :: ptmp
  __TYPE__ :: tmp
#endif

  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 1
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             __TYPE_STR__ // '(' // &
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
           __TYPE_STR__ // '(' // &
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
  call log_allocation(  addr, __TYPE_STR__ // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine FNAME(alloc_bound_1d)

#endif

subroutine FNAME(freemem)( p, file, line, stat )

  use iso_c_binding

  implicit none

  __TYPE__, pointer :: p

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
   call log_deallocation( addr, __TYPE_STR__, &
                        bsize, file, line )

   ! deallocate memory
#ifdef USE_POSIX_MEMALIGN
   call free_f( c_loc(p) )
   ! free function does not return a value
   ierr = 0
#else
   deallocate( p, stat = ierr )
#endif

   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            __TYPE_STR__ // ', ' // &
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

end subroutine FNAME(freemem)

#if __MAX_DIM__ > 0

subroutine FNAME(freemem_1d)( p, file, line, stat )

  use iso_c_binding

  implicit none

  __TYPE__, dimension(:), pointer :: p

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
     call log_deallocation( addr, __TYPE_STR__ // '(:)', &
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
            __TYPE_STR__ // '(' // &
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

end subroutine FNAME(freemem_1d)

#endif
#if __MAX_DIM__ > 1

subroutine FNAME(freemem_2d)( p, file, line, stat )

  use iso_c_binding

  implicit none

  __TYPE__, dimension(:,:), pointer :: p

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
     call log_deallocation( addr, __TYPE_STR__ // '(:,:)', &
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
            __TYPE_STR__ // '(' // &
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

end subroutine FNAME(freemem_2d)

#endif


#undef __MAX_DIM__
#undef __TYPE_STR__
#undef __TYPE__
#undef FNAME

