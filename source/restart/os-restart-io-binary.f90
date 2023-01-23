!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     restart I/O fuctions for fortran binary output backend
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! $URL: $
! $Id: $

#ifndef __TEMPLATE__

#define __TEMPLATE__

! string
#define __TYPE__ character(len = *)
#undef __ARRAY__
#define __STRING__
#define FNAME(a) a ## _string
#include __FILE__

! string 1d
#define __TYPE__ character(len = *), dimension(:)
#define __ARRAY__
#define __STRING__
#define FNAME(a) a ## _string1d
#include __FILE__

! string 2d
#define __TYPE__ character(len = *), dimension(:,:)
#define __ARRAY__
#define __STRING__
#define FNAME(a) a ## _string2d
#include __FILE__

! logical
#define __TYPE__ logical
#undef __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _logical
#include __FILE__

! logical 1d
#define __TYPE__ logical, dimension(:)
#define __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _logical1d
#include __FILE__

! integer
#define __TYPE__ integer
#undef __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _integer
#include __FILE__

! int64
#define __TYPE__ integer(p_int64)
#undef __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _int64
#include __FILE__

! int64 1d
#define __TYPE__ integer(p_int64), dimension(:)
#define __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _int641d
#include __FILE__

! integer 1d
#define __TYPE__ integer, dimension(:)
#define __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _integer1d
#include __FILE__

! integer 2d
#define __TYPE__ integer, dimension(:,:)
#define __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _integer2d
#include __FILE__

! integer 3d
#define __TYPE__ integer, dimension(:,:,:)
#define __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _integer3d
#include __FILE__

! integer 4d
#define __TYPE__ integer, dimension(:,:,:,:)
#define __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _integer4d
#include __FILE__

! single
#define __TYPE__ real(p_single)
#undef __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _single
#include __FILE__

! single 1d
#define __TYPE__ real(p_single), dimension(:)
#define __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _single1d
#include __FILE__

! single 2d
#define __TYPE__ real(p_single), dimension(:,:)
#define __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _single2d
#include __FILE__

! single 3d
#define __TYPE__ real(p_single), dimension(:,:,:)
#define __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _single3d
#include __FILE__

! single 4d
#define __TYPE__ real(p_single), dimension(:,:,:,:)
#define __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _single4d
#include __FILE__

! double
#define __TYPE__ real(p_double)
#undef __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _double
#include __FILE__

! double 1d
#define __TYPE__ real(p_double), dimension(:)
#define __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _double1d
#include __FILE__

! double 2d
#define __TYPE__ real(p_double), dimension(:,:)
#define __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _double2d
#include __FILE__

! double 3d
#define __TYPE__ real(p_double), dimension(:,:,:)
#define __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _double3d
#include __FILE__

! double 4d
#define __TYPE__ real(p_double), dimension(:,:,:,:)
#define __ARRAY__
#undef __STRING__
#define FNAME(a) a ## _double4d
#include __FILE__


#else 


!---------------------------------------------------------------------------------------------------
!  Template Function definitions
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine FNAME(restart_io_write)( varname, var, restart_handle, ierr )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  character(len=*), intent(in) :: varname
  __TYPE__, intent(in) :: var
  type(t_restart_handle), intent(inout) :: restart_handle
  integer, intent(out) :: ierr

!       local variables
  integer(8)         :: datasize

  if ( restart_handle%if_acc_size ) then
    #ifdef __ARRAY__
     #ifdef __STRING__
       datasize=len(var)*size(var)
     #else
       datasize=kind(var)*size(var)
     #endif
    #else
     #ifdef __STRING__
       datasize=len(var)
     #else
       datasize=kind(var)
     #endif
    #endif
    restart_handle%data_size = restart_handle%data_size + datasize + 8
    ierr=0  
  else
    write(unit=restart_handle%file_id, iostat=ierr) var
  endif

end subroutine FNAME(restart_io_write)
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine FNAME(restart_io_read)( var, restart_handle, ierr )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  __TYPE__, intent(out) :: var
  type(t_restart_handle), intent(in) :: restart_handle
  integer, intent(out) :: ierr

    read(unit=restart_handle%file_id, iostat=ierr) var

end subroutine FNAME(restart_io_read)
!---------------------------------------------------------------------------------------------------

#undef FNAME
#undef __TYPE__

#endif
