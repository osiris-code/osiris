!-------------------------------------------------------------------------------------------------
! m_logprof module
!
! This module handles logging and profiling of events. It can optionally be used with the MPE
! library and/or with the PAPI library, just set the compiler macros __USE_MPE__ and or
! __USE_PAPI__ and link with the appropriate libs.
!
!-------------------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_logprof

use m_parameters
use m_system
use m_diagfile, only : p_time_length

implicit none

! restrict access to things explicitly declared public
private

integer, parameter :: p_max_events = 64
integer, parameter :: p_max_event_name = 128

!---------------------------------------------------------------------------------------------------
! PAPI event data
!---------------------------------------------------------------------------------------------------

#ifdef __USE_PAPI__

#warning Using PAPI library for profiling

! Number of hardware counters being used
integer, parameter :: p_nhwc = 3

type t_papi_event

  ! Timing information
  integer*8 :: rcy		= 0	! real cycles
  integer*8 :: rus		= 0	! real microseconds
  integer*8 :: ucy		= 0	! user cycles
  integer*8 :: uus		= 0	! user microseconds

  ! Hardware counters
  integer*8, dimension( p_nhwc ) :: hwc

end type t_papi_event

#endif
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! MPE event data
!---------------------------------------------------------------------------------------------------
#ifdef __USE_MPE__
integer, parameter :: num_colors = 16

character(len=*), dimension(num_colors), parameter :: colornames = &
   (/ "white      ", "black      ", "red        ", "yellow     ", "green      ", &
	  "cyan       ", "blue       ", "magenta    ", "aquamarine ", "forestgreen", &
	  "orange     ", "maroon     ", "brown      ", "pink       ", "coral      ", &
	  "gray       " /)

integer, external :: MPE_Log_get_event_number, MPE_Describe_state
integer, external :: MPE_Log_event
#endif
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! t_event class definition
!---------------------------------------------------------------------------------------------------
type :: t_event

  character(len = p_max_event_name) :: name

#ifdef __USE_MPE__
  ! Extra data for MPE
  integer :: idbegin, idend
#endif

#ifdef __USE_PAPI__
  type( t_papi_event ) :: evdelta, evtotal
#else
  integer( p_int64 ) :: timer_begin
  real( p_double )   :: total_time
#endif

end type t_event
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Module variables
!---------------------------------------------------------------------------------------------------

! Total number of events created
integer, save :: num_events = 0

! Events array
type( t_event ), dimension(p_max_events), save :: events

!---------------------------------------------------------------------------------------------------

interface create_event
  module procedure create_event
end interface

interface begin_event
  module procedure begin_event
end interface

interface end_event
  module procedure end_event
end interface

interface total_event_time
  module procedure total_event_time
end interface

interface list_total_event_times
  module procedure list_total_event_times
end interface



public :: create_event, begin_event, end_event
public :: total_event_time, list_total_event_times

#ifdef __USE_PAPI__

interface init_papi
  module procedure init_papi
end interface

public :: init_papi

#endif


 contains

!---------------------------------------------------------------------------------------------------
! Create an event for profiling
!---------------------------------------------------------------------------------------------------
function create_event( name )

  implicit none

  integer :: create_event
  character(len=*), intent(in)  ::  name

#ifdef __USE_MPE__
  integer :: ierr
#endif

  num_events = num_events+1

  if ( num_events > p_max_events ) then
    if (mpi_node() == 0 ) then
      write(0, "('(*error*) Unable to create event ', A,', too many events.')" ) trim(name)
      write(0, "('(*error*) The maximum number of events is currently set at ',I0)" ) p_max_events
    endif
    call abort_program(p_err_invalid)
  endif

  create_event = num_events

  events(num_events)%name = name

#ifdef __USE_MPE__
  events(num_events)%idbegin =  MPE_Log_get_event_number()
  events(num_events)%idend   =  MPE_Log_get_event_number()
  ierr = MPE_Describe_state( events(num_events)%idbegin, &
							 events(num_events)%idend,   &
							 events(num_events)%name,    &
							 trim( colornames(mod(num_events+1,16) ) ) )
#endif

#ifdef __USE_PAPI__
  call zero_papi_event( events(num_events)%evtotal )
#else
  events(num_events)%total_time = 0.0D0
#endif

end function create_event
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Start an event
!---------------------------------------------------------------------------------------------------
subroutine begin_event( event_num )

  implicit none

  integer, intent(in) :: event_num

#ifdef __USE_MPE__
  integer :: ierr

  if( if_profiling() ) then
    ierr = MPE_Log_Event(events(event_num)%idbegin, 0, "start")
#else
  if( if_profiling() ) then
#endif

#ifdef __USE_PAPI__
   ! Start the PAPI event
   call start_papi_event( events(event_num)%evdelta )
#else
   ! Read the current timer ticks
   events(event_num)%timer_begin = timer_ticks()
#endif
  endif

end subroutine begin_event
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Finish an event
!---------------------------------------------------------------------------------------------------
subroutine end_event( event_num )

  implicit none

  integer, intent(in) :: event_num

#ifdef __USE_MPE__
  integer :: ierr
#endif

  integer( p_int64 ) :: timer_end

  if( if_profiling() ) then
! Terminate the event
#ifdef __USE_PAPI__
  call finish_papi_event( events(event_num)%evdelta )
#else
  timer_end = timer_ticks()
#endif

#ifdef __USE_MPE__
  ierr = MPE_Log_Event(events(event_num)%idend, 0, "end")
#endif

! Accumulate time
#ifdef __USE_PAPI__
  call accum_papi_event( events(event_num)%evtotal, events(event_num)%evdelta )
#else
  events(event_num)%total_time = events(event_num)%total_time + &
                                 timer_interval_seconds(events(event_num)%timer_begin, timer_end)
#endif
  endif

end subroutine end_event
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Return the total event time for a given event
!---------------------------------------------------------------------------------------------------
function total_event_time( event_num )

  implicit none

  real( p_double ) :: total_event_time

  integer, intent(in)  ::  event_num

#ifdef __USE_PAPI__
  total_event_time = 1.0d-6 * events(event_num)%evtotal%uus
#else
  total_event_time = events(event_num)%total_time
#endif

end function total_event_time
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! List the avg, min and max times for all events
!---------------------------------------------------------------------------------------------------
subroutine list_total_event_times( n, comm, label )

  !use mpi

  implicit none

  ! dummy variables
  integer, intent(in) :: n, comm
  character(len=*), optional, intent(in) :: label

  ! local variables
  integer :: i, rank, size, ierr
  real( p_double ) , dimension( num_events ) :: time, time_min, time_max, time_avg
  character(len=128) :: fname
  character(len=2) :: tim_len

#ifdef __USE_PAPI__
  type( t_papi_event ), dimension( num_events ) :: total

  ! stop papi counters temporarily
  call papi_stop_counters()
#endif


  fname = trim(path_time)//p_dir_sep//'timings-'
  if ( present( label ) ) then
    fname = trim(fname) // trim(label)
  else
    write( tim_len, '(i2)' ) p_time_length
    write( fname, '(A,i'//trim(tim_len)//'.'//trim(tim_len)//')' ) trim(fname), n
  endif

  if ( comm == MPI_COMM_NULL ) then


	rank = 0

	! Create output directory if necessary
	call mkdir( trim(path_time), ierr )
	if ( ierr /= 0 ) then
	  write( 0, * ) 'Error creating TIMINGS ', strerror( ierr )
	  return
	endif

    open ( unit = file_id_prof , file = trim(fname), form = 'formatted' )

	write( file_id_prof, '(A39,1X,A19)' ) &
			     ' Event', 'Total [s]'
    ! Serial run
	write( file_id_prof, '(A)' ) &
				 '-----------------------------------------------------------'

#ifdef __USE_PAPI__
	do i = 1, num_events
	  write( file_id_prof, '(A39,1X,E19.3)' ) trim( events(i)%name ), 1.0d-6*events(i)%evtotal%rus
	enddo
#else
	do i = 1, num_events
	  write( file_id_prof, '(A39,1X,E19.3)' ) trim( events(i)%name ), events(i)%total_time
	enddo
#endif

    close(file_id_prof)

  else

    call mpi_comm_rank( comm, rank, ierr )
    call mpi_comm_size( comm, size, ierr )

    ! Parallel run
#ifdef __USE_PAPI__
    do i = 1, num_events
      time(i) = 1.0d-6 * events(i)%evtotal%rus
    enddo
#else
    do i = 1, num_events
      time(i) = events(i)%total_time
    enddo
#endif

    ! Get min, max and avg times
    call mpi_reduce( time, time_min, num_events, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr )
    call mpi_reduce( time, time_max, num_events, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr )

    call mpi_reduce( time, time_avg, num_events, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr )

    if ( rank == 0 ) then

       do i = 1, num_events
         time_avg(i) = time_avg(i) / size
       enddo

	   ! Create output directory if necessary
	   call mkdir( trim(path_time), ierr )
	   if ( ierr /= 0 ) then
		 write( 0, * ) 'Error creating TIMINGS ', strerror( ierr )
		 return
	   endif

	   open ( unit = file_id_prof , file = trim(fname), form = 'formatted' )

       write( file_id_prof, * ) ' Iterations = ', n
       write( file_id_prof, * ) ' '
       write( file_id_prof, '(A39,3(1X,A19))' ) &
                    'Event', 'Avg [s]', 'Min [s]', 'Max [s]'
       write( file_id_prof, '(A)' ) &
                    '------------------------------------------------------------------------------'

       do i = 1, num_events
         write( file_id_prof, '(A39,3(1X,E19.6))' ) trim( events(i)%name ), &
             time_avg(i), time_min(i), time_max(i)
       enddo

       close(file_id_prof)

    endif

  endif

#ifdef __USE_PAPI__

  ! Gather papi data on root node
  call gather_papi_event( events(1:num_events), total, comm )

  ! Open file for papi report
  if ( rank == 0 ) then

	 open ( unit = file_id_prof , file = trim(fname)//'-papi', form = 'formatted' )

	 write( file_id_prof, * ) ' Iterations = ', n
	 write( file_id_prof, * ) ' '
	 write( file_id_prof, '(A39,7(1X,A19))' ) &
				  'Event', 'real cycles', 'real usec', 'user cycles', 'user usec', &
				  'TOT_INS', 'TOT_FP', 'L2_DCM'
	 write( file_id_prof, '(A)' ) &
		'----------------------------------------------------------------------------------' // &
		'----------------------------------------------------------------------------------'

     ! write event data and close file

	 do i = 1, num_events
	   write( file_id_prof, '(A39,7(1X,I19))' ) trim( events(i)%name ), &
		   total(i)%rcy, total(i)%rus, total(i)%ucy, total(i)%uus, &
		   total(i)%hwc(1), total(i)%hwc(2), total(i)%hwc(3)
	 enddo

	 close(file_id_prof)

  endif

  ! restart papi counters
  call papi_start_counters()

#endif

end subroutine list_total_event_times
!---------------------------------------------------------------------------------------------------

#ifdef __USE_PAPI__

! PAPI specific code

!-----------------------------------------------------------------------------------------
! Initialize the PAPI hardware counters that will be used by PAPI events
!-----------------------------------------------------------------------------------------
subroutine init_papi( )

  implicit none

  integer :: ierr

  call init_papi_counters( ierr )

  if ( ierr /= 0 ) then
    write(0,*) 'Unable to initialize the required PAPI hardware counters'
    call abort_program()
  endif

end subroutine init_papi
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleanup event data
!-----------------------------------------------------------------------------------------
subroutine zero_papi_event( event )

  implicit none

  type( t_papi_event ), intent(inout) :: event

  event%rcy = 0
  event%ucy = 0
  event%rus = 0
  event%uus = 0

  event%hwc(1) = 0
  event%hwc(2) = 0
  event%hwc(3) = 0

end subroutine zero_papi_event
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Start the PAPI event
!-----------------------------------------------------------------------------------------
subroutine start_papi_event( event )

  implicit none

  type( t_papi_event ), intent(inout) :: event

  ! this is done with PAPI read counters because of nested events
  call get_papi_info( event%rcy, event%ucy, event%rus, event%uus, event%hwc )

end subroutine start_papi_event
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Finish the PAPI event
!-----------------------------------------------------------------------------------------
subroutine finish_papi_event( event )

  implicit none

  type( t_papi_event ), intent(inout) :: event

  type( t_papi_event ) :: finish

  call get_papi_info( finish%rcy, finish%ucy, finish%rus, finish%uus, finish%hwc )

  event%rcy = finish%rcy - event%rcy
  event%ucy = finish%ucy - event%ucy
  event%rus = finish%rus - event%rus
  event%uus = finish%uus - event%uus

  event%hwc(1) = finish%hwc(1) - event%hwc(1)
  event%hwc(2) = finish%hwc(2) - event%hwc(2)
  event%hwc(3) = finish%hwc(3) - event%hwc(3)

end subroutine finish_papi_event
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Accumulate PAPI events
!-----------------------------------------------------------------------------------------
subroutine accum_papi_event( total, delta )

  implicit none

  type( t_papi_event ), intent(inout) :: total
  type( t_papi_event ), intent(in)    :: delta

  total%rcy = total%rcy + delta%rcy
  total%ucy = total%ucy + delta%ucy
  total%rus = total%rus + delta%rus
  total%uus = total%uus + delta%uus

  total%hwc(1) = total%hwc(1) + delta%hwc(1)
  total%hwc(2) = total%hwc(2) + delta%hwc(2)
  total%hwc(3) = total%hwc(3) + delta%hwc(3)

end subroutine accum_papi_event
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Global papi measurement
!-----------------------------------------------------------------------------------------
subroutine gather_papi_event( events, total, comm )

  !use mpi

  implicit none

  type( t_event ), dimension(:), intent(in)  :: events
  type( t_papi_event ), dimension(:), intent(out) :: total
  integer, intent(in) :: comm

  integer*8, dimension(:), allocatable :: buffer_in, buffer_out

  integer :: i, k, ierr, rank, nevents

  nevents = size(events)

  if ( comm == MPI_COMM_NULL ) then
    ! single node run, just copy results to output
    do i = 1, nevents
      total(i)%rcy = events(i)%evtotal%rcy
      total(i)%ucy = events(i)%evtotal%ucy
      total(i)%rus = events(i)%evtotal%rus
      total(i)%uus = events(i)%evtotal%uus

      total(i)%hwc(1) = events(i)%evtotal%hwc(1)
      total(i)%hwc(2) = events(i)%evtotal%hwc(2)
      total(i)%hwc(3) = events(i)%evtotal%hwc(3)
    enddo

  else

	 ! Allocate communication buffers
	 allocate( buffer_in( 4*nevents ), buffer_out( 4*nevents ), stat = ierr )
	 if ( ierr /= 0 ) then
	   ERROR('Allocation failed')
	   call abort_program()
	 endif

	 ! Get local rank
	 call MPI_COMM_RANK( comm, rank, ierr )
	 if ( ierr /= 0 ) then
	   ERROR('MPI Error')
	   call abort_program()
	 endif

	 ! Gather event data in communications buffer
	 k = 1
	 do i = 1, nevents
	   buffer_in(k)   = events(i)%evtotal%rcy
	   buffer_in(k+1) = events(i)%evtotal%ucy
	   buffer_in(k+2) = events(i)%evtotal%rus
	   buffer_in(k+3) = events(i)%evtotal%uus
	   k = k + 4
	 enddo

	 ! For times we want the maximum values because this dominates for the work group
	 ! Since events are timed with the local clock at the beggining and end of each event, only the
	 ! elapsed time is recorded, and we shouldn't need to synchronize the clocks.
	 call MPI_REDUCE( buffer_in, buffer_out, 4 * nevents, MPI_INTEGER8, MPI_MAX, 0, comm, ierr )
	 if ( ierr /= 0 ) then
	   ERROR('MPI Error')
	   call abort_program()
	 endif

	 ! copy the values to the output variable on root node
	 if ( rank == 0 ) then
	   k = 1
	   do i = 1, nevents
		 total(i)%rcy = buffer_out(k)
		 total(i)%ucy = buffer_out(k+1)
		 total(i)%rus = buffer_out(k+2)
		 total(i)%uus = buffer_out(k+3)
		 k = k + 4
	   enddo
	 endif

	 ! Gather hardware counters event data in communications buffer
	 k = 1
	 do i = 1, nevents
	   buffer_in(k)   = events(i)%evtotal%hwc(1)
	   buffer_in(k+1) = events(i)%evtotal%hwc(2)
	   buffer_in(k+2) = events(i)%evtotal%hwc(3)
	   k = k + 3
	 enddo

	 ! For hardware counters we want the total values
	 call MPI_REDUCE( buffer_in, buffer_out, 3 * nevents, MPI_INTEGER8, MPI_SUM, 0, comm, ierr )
	 if ( ierr /= 0 ) then
	   ERROR('MPI Error')
	   call abort_program()
	 endif

	 ! copy the values to the output variable on root node
	 if ( rank == 0 ) then
	   k = 1
	   do i = 1, nevents
		 total(i)%hwc(1) = buffer_out(k)
		 total(i)%hwc(2) = buffer_out(k+1)
		 total(i)%hwc(3) = buffer_out(k+2)
		 k = k + 3
	   enddo
	 endif

	 ! Free communication buffers
	 deallocate( buffer_in, buffer_out , stat = ierr )
	 if ( ierr /= 0 ) then
	   ERROR('Deallocation failed')
	   call abort_program()
	 endif

  endif

end subroutine gather_papi_event
!-----------------------------------------------------------------------------------------


#endif


end module m_logprof
