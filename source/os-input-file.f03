
#include "os-config.h"
#include "os-preprocess.fpp"

!-----------------------------------------------------------------------------------------
!> @brief Input file module.
!! @details Contains classes and functions for reading in the input file.
!<----------------------------------------------------------------------------------------
module m_input_file

use m_parameters
use m_system

implicit none

private


!> @param Maximum size for input file after parsing
integer, parameter :: p_max_input_file_size = 262144    ! 256 k

!> @param Maximum size for a namelist section after parsing
integer, parameter :: p_max_nml_size = 65536            ! 64 k

!> @param Maximum size for an input file line before parsing
integer, parameter :: p_max_nml_line = 1024          ! 1k

!> @param Maximum size for the name of the file sections
integer, parameter, public :: p_max_nml_section_name = 32

! 
!-----------------------------------------------------------------------------------------
!> @class t_input_file
!! @brief Input file buffer
!! @details Stores the input file and current position.
!<----------------------------------------------------------------------------------------
type t_input_file

  character, dimension(p_max_input_file_size) :: cbuf !< C buffer
  integer :: data_size = 0 !< Data size
  integer :: pos = -1 !< Current position in buffer

  character(len=p_max_nml_size) :: nml_text !< Fortran buffer

  logical :: display_out = .true. !< Whether or not to display output when reading/setting up objects

contains
  procedure :: get_section_name

end type

! Read input file
interface open_input
  module procedure open_input_file
  module procedure directly_set_input
end interface

! get namelist text
interface get_namelist
  module procedure get_namelist
end interface

! Functions for displaying output during read
interface disp_out
  module procedure disp_out
end interface

interface no_disp_out
  module procedure no_disp_out
end interface

public :: path_mass, path_rest, path_hist, path_time
public :: t_input_file

public :: open_input, get_namelist
public :: disp_out, no_disp_out, p_max_input_file_size

contains

!-----------------------------------------------------------------------------------------
!> @brief Directly impose an input deck.
!! @details Initialize the input file synthetically by directly imposing an input deck
!! (rather than reading from a file).
!! @param[in,out] input_file        Input file type
!! @param[in,out] synthetic_input   Synthetic input to be stored in input_file
!! @param[in]     comm              Integer communicator
!! @param[in]     involve_mpi       Whether or not to use MPI
!<----------------------------------------------------------------------------------------
subroutine directly_set_input( input_file, synthetic_input, comm, involve_mpi)
  use mpi

  use stringutil

  implicit none

  class( t_input_file ), intent(inout) :: input_file
  character(len=*), target, intent(inout) :: synthetic_input
  integer, intent(in) :: comm
  logical, intent(in) :: involve_mpi

  integer :: mpi_rank, mpi_size, ierr, i

  ! check if it fits in the buffer
  if ( len(synthetic_input) > p_max_input_file_size ) then
    ! we could grow the buffer instead, but this will never happen
    write(0,*) "Input file buffer too small, aborting."
    stop
  endif

  if( involve_mpi ) then
    if ( comm == MPI_COMM_NULL ) then
      mpi_rank = 0
      mpi_size = 1
    else
      call mpi_comm_rank( comm, mpi_rank, ierr )
      if ( ierr /= 0 ) then
        write(0,*) 'MPI Error'
        stop
      endif
      call mpi_comm_size( comm, mpi_size, ierr )
      if ( ierr /= 0 ) then
        write(0,*) 'MPI Error'
        stop
      endif
    endif
  endif

    
  ! Just Copy over the synthetic input:
  !      To all nodes       if involve_mpi == False
  !      Just to the root   if involve_mpi == True and mpi_rank == 0
  if( .not. involve_mpi .or. mpi_rank == 0) then
    do i = 1, len(synthetic_input)
      input_file%cbuf(i) = synthetic_input(i:i)
    enddo
  endif

  if( involve_mpi ) then
    ! When running in parallel send character buffer to all nodes
    if ( mpi_size > 1 ) then
    ! broadcast buffer size
    call mpi_bcast( input_file%data_size, 1, MPI_INTEGER, 0, comm, ierr )
      if ( ierr /= 0 ) then
        write(0,*) '(*error*) MPI Error'
        stop
      endif

    ! broadcast data
    call mpi_bcast( input_file%cbuf, input_file%data_size, MPI_CHARACTER, 0, &
              comm, ierr )
      if ( ierr /= 0 ) then
        write(0,*) '(*error*) MPI Error'
        stop
      endif
    endif
  endif

  ! reset buffer pointer
  input_file%pos = 0
  input_file%data_size = len(synthetic_input)
end subroutine

!-----------------------------------------------------------------------------------------
!> @brief Open the input file
!! @details Read in the entire input file, then broadcast to other MPI ranks
!! @param[in,out] input_file  Input file object
!! @param[in]     fname       File name
!! @param[in]     comm        MPI communicator
!<----------------------------------------------------------------------------------------
subroutine open_input_file( input_file, fname, comm )

  use mpi

  use stringutil

  implicit none

  class( t_input_file ), intent(inout) :: input_file
  character(len=*), intent(in) :: fname
  integer, intent(in) :: comm

  integer :: mpi_rank, mpi_size, ierr
  integer :: pos, line, nml_start
  logical :: in_namelist
  character(len = p_max_nml_line) :: buffer, linebuffer
  character(len=6) :: nml_size
  character :: quote_char

  ! get parallel universe information
  if ( comm == MPI_COMM_NULL ) then
    mpi_rank = 0
    mpi_size = 1
  else
    call mpi_comm_rank( comm, mpi_rank, ierr )
    if ( ierr /= 0 ) then
      write(0,*) 'MPI Error'
      stop
    endif
    call mpi_comm_size( comm, mpi_size, ierr )
    if ( ierr /= 0 ) then
      write(0,*) 'MPI Error'
      stop
    endif
  endif

  ! only root node reads input file
  if ( mpi_rank == 0 ) then

  ! Open input file
  open (unit = file_id_tem, file = fname, &
      access = 'sequential', form = 'formatted', &
      status = 'old', iostat = ierr, action = "read")

  if ( ierr /= 0 ) then
    write(0,*) '(*error*) Unable to read input file : file not found or unreadable'
      stop
  endif

  line = 0
  quote_char = ''
  input_file%pos = 0
  in_namelist = .false.

  do
    line = line + 1

    read(unit = file_id_tem, iostat = ierr, fmt = '(A)') linebuffer

    if (ierr < 0) exit  ! eof

      call check_io_err( ierr, line )

    buffer = trim_input(linebuffer, quote_char)

    do
    if (len_trim(buffer) < 1) exit

    if (.not. in_namelist) then ! outside namelist
      nml_start = input_file%pos
      call store( input_file, "      ,")

      ! get section name
      pos = scan(buffer, "{", .false.)

      if (pos == 0) then ! { is in the next line

      call store( input_file, "&nl_"//trim(buffer)//" " )

      line = line + 1
      read(unit = file_id_tem, iostat = ierr, fmt = '(A)') linebuffer
            call check_io_err( ierr, line )

      buffer = trim_input(linebuffer, quote_char)
      if (buffer(1:1) /= '{') then
        write(0,*) "error processing input file, '{' must follow section name"
        write(0,*) "text : ", trim(linebuffer)
        write(0,*) "line : ", line
        stop
      endif
      pos = 1
      else ! { is in the same line
      call store( input_file, "&nl_"//trim(buffer(1:pos-1))//" " )
      endif

      in_namelist = .true.
      buffer = trim(buffer(pos+1:))

    else ! inside namelist

      pos = find_char_unquote(buffer, "}")

      if (pos == 0) then
      call store( input_file, trim(buffer) )
      exit
      else
      in_namelist = .false.
      if (pos > 1) call store( input_file, buffer(1:pos-1) )
      call store( input_file, "/ " )
      buffer = trim(buffer(pos+1:))

          ! Check if nml_size is over the limit of p_max_nml_size because
          ! in the read_nml sections the temp string will be created with this size
      if ( input_file%pos - nml_start > p_max_nml_size ) then
        write(0,*) "(*error*) error processing input file, the section is too long"
        write(0,*) "text : ", trim(linebuffer)
        write(0,*) "line : ", line
        stop
      else
        write(nml_size,'(i6)') input_file%pos - nml_start
      endif

      ! write the namelist section size
      call store_pos( input_file, nml_size, nml_start )
      endif

    endif
    enddo

  enddo

  ! check if all namelists were properly closed
  if ( in_namelist ) then
     write(0,*) "(*error*) error processing input file, end of file reached with unterminated input section."
     stop
  endif

    ! get the total size of the processed data
    input_file%data_size = input_file%pos

    ! close the input file
    close( file_id_tem )

  endif

  ! When running in parallel send character buffer to all nodes
  if ( mpi_size > 1 ) then
  ! broadcast buffer size
  call mpi_bcast( input_file%data_size, 1, MPI_INTEGER, 0, comm, ierr )
    if ( ierr /= 0 ) then
      write(0,*) '(*error*) MPI Error'
      stop
    endif

  ! broadcast data
  call mpi_bcast( input_file%cbuf, input_file%data_size, MPI_CHARACTER, 0, &
            comm, ierr )
    if ( ierr /= 0 ) then
      write(0,*) '(*error*) MPI Error'
      stop
    endif
  endif

  ! reset buffer pointer
  input_file%pos = 0

contains

!-----------------------------------------------------------------------------------------
!> @brief Check for I/O error
!! @param[in]     ierr  Integer error value
!! @param[in]     line  Line
!<----------------------------------------------------------------------------------------
subroutine check_io_err( ierr, line )

   implicit none

   integer, intent(in) :: ierr, line

   if ( ierr /= 0 ) then
      write(0,*) "(*error*) Error reading input file, line ", line
      if ( ierr < 0 ) then
         write(0,*) "(*error*) Unexpected end-of-file or end-of-record"
      else
         write(0,*) "(*error*) I/O error #", ierr
      endif

      stop
   endif

end subroutine check_io_err

!-----------------------------------------------------------------------------------------
!> @brief Stores a string at the end of the character buffer
!! @param[in,out] input_file  Input file object
!! @param[in]     str         String to store
!<----------------------------------------------------------------------------------------
subroutine store( input_file, str )

  implicit none

  class( t_input_file ), intent(inout) :: input_file
  character(len=*), intent(in) :: str

  integer :: i, l

  l = len(str)

  ! check if it fits in the buffer
  if ( input_file%pos + l > p_max_input_file_size ) then
  ! we could grow the buffer instead, but this will never happen
  write(0,*) "Input file buffer too small, aborting."
  stop
  endif

  do i = 1, len(str)
    input_file%pos = input_file%pos + 1
    input_file%cbuf(input_file%pos) = str(i:i)
  enddo

end subroutine store

!-----------------------------------------------------------------------------------------
!> @brief Stores a string at a specific position of the character buffer
!! @param[in,out] input_file  Input file object
!! @param[in]     str         String to store
!! @param[in]     pos         Position to store
!<----------------------------------------------------------------------------------------
subroutine store_pos( input_file, str, pos )

  implicit none

  class( t_input_file ), intent(inout) :: input_file
  character(len=*), intent(in) :: str
  integer, intent(in) :: pos

  integer :: i, l

  l = len(str)

  ! check if it fits in the buffer
  if ( pos + l > p_max_input_file_size ) then
  ! we could grow the buffer instead, but this will never happen
  write(0,*) "Input file buffer too small, aborting."
  stop
  endif

  do i = 1, len(str)
    input_file%cbuf(pos+i) = str(i:i)
  enddo

end subroutine store_pos

!-----------------------------------------------------------------------------------------
!> @brief Remove blanks outside of quotes and strip comments from input string
!! @param[in]     in_str      Input string
!! @param[in,out] quote_char  Character for quotes
!! @return        trim_input  Trimmed string
!<----------------------------------------------------------------------------------------
function trim_input( in_str, quote_char )

  use stringutil

  implicit none

  character ( len = * ), intent(in) :: in_str
  character, intent(inout) :: quote_char
  character (len = len(in_str)) :: trim_input

!> @param Comment character
  character, parameter :: p_comment_char = '!'

  integer :: i, j

  trim_input = ""

  j = 1
  do i = 1, len(in_str)
  if (quote_char /= '') then
    trim_input(j:j) = in_str(i:i)
    j = j +1
    if (in_str(i:i) == quote_char) quote_char = ''
  elseif (in_str(i:i) == p_comment_char) then
    exit
  elseif (in_str(i:i) == '"' .or. in_str(i:i) == "'") then
    quote_char = in_str(i:i)
    trim_input(j:j) = in_str(i:i)
    j = j +1
  elseif ((in_str(i:i) /= p_space) .and. &
      (in_str(i:i) /= p_tab))  then
    trim_input(j:j) = in_str(i:i)
    j = j +1
  endif
  enddo

end function trim_input

end subroutine open_input_file
! ---------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!> @brief Gets a string from a character buffer
!! @param[in]     cbuf      Character buffer
!! @param[in]     start     Start index
!! @param[in]     strlen    Length of string
!! @return        get_str   Desired string
!<----------------------------------------------------------------------------------------
function get_str( cbuf, start, strlen )

  implicit none

  character, dimension(:), intent(in) :: cbuf
  integer, intent(in) :: start, strlen

  character(len=strlen) :: get_str

  integer :: i

  do i = 1, strlen
    get_str(i:i) = cbuf(start + i - 1)
  enddo

end function get_str
! ---------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> @brief Get a section name
!! @param[in]     this            Input file object
!! @return        nml_name_file   Namelist section name
!<----------------------------------------------------------------------------------------
function get_section_name( this ) result( nml_name_file )

  implicit none

  class(t_input_file), intent(in) :: this

  character(len=p_max_nml_section_name) :: nml_name_file
  integer :: i,p

  ! get the namelist name in the buffer
  i = this%pos + 9
  p = 1
  nml_name_file = ""
  do
    if ( i > this%data_size ) exit ! protection against EOF
    if ( this%cbuf(i) == " " ) exit
    nml_name_file(p:p) = this%cbuf(i)
    p = p + 1
    i = i + 1
  enddo

end function get_section_name
! ---------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> @brief Returns the text of the requested namelist if present
!! @param[in,out] this        Input file object
!! @param[in]     nml_name    Namelist name
!! @param[out]    stat        Custom-made error status
!! @param[in]     nml_output  Optional, used to replace nml_name if present
!<----------------------------------------------------------------------------------------
subroutine get_namelist( this, nml_name, stat, nml_output )

  implicit none

  class(t_input_file), intent(inout) :: this
  character(len=*), intent(in) :: nml_name
  integer, intent(out) :: stat
  character(len=*), intent(in), optional :: nml_output

  character(len=6) :: nml_size
  character(len = p_max_nml_section_name ) :: nml_name_file
  integer :: nml_text_len
  integer :: len_name, len_output, subst_pos

  if ( this%pos + 6 > this%data_size ) then
    ! eof reached, data is not here
    stat = 1
    return
  endif

  nml_name_file = this % get_section_name()

  ! get the namelist size
  nml_size = get_str( this%cbuf, this%pos + 1, 6 )
  read( nml_size, '(i6)' ) nml_text_len

  ! debug
  !print *, '>> Found namelist "',trim( nml_name_file ),'"'
  !print *, '>> matches: "',trim( nml_name ) == trim( nml_name_file ),'"'

  ! check the namelist name
  if ( trim( nml_name ) /= trim( nml_name_file ) ) then
    ! Another section is here
    stat = 1
    return
  endif

  ! We have the proper section lets get the data
  stat = 0
  this%nml_text = get_str( this%cbuf, this%pos + 8, nml_text_len - 7 )

  ! debug
  ! print *, 'Namelist text >'//trim( nml_text )//'<'

  ! replaces found namelist by desired namelist if 'nml_output' is present
  if( present( nml_output ) ) then
    len_output = len_trim( nml_output )
    len_name = len_trim( nml_name )

    subst_pos = index( this%nml_text, nml_name(:len_name) )

    this%nml_text = this%nml_text(:subst_pos-1) // &
                          nml_output(:len_output) // &
                          this%nml_text(subst_pos+len_name:)
  endif

  ! advance buffer pointer
  this%pos = this%pos + nml_text_len


end subroutine get_namelist
! ---------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> @brief Whether or not to display output while reading input
!! @details When initializing the simulation with tiles, for example, a global read of
!! the input file is performed first, followed by a read for each tile.  Since some output
!! is generated upon reading, this flag is set to false after the global read so that the
!! output is not displayed for each tile.
!! @param[in]     this      Input file object
!! @return        disp_out  Whether or not to display output
!<----------------------------------------------------------------------------------------
pure function disp_out( this )

  implicit none

  class( t_input_file ), intent(in) :: this

  logical :: disp_out

  disp_out = this%display_out

end function disp_out
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> @brief Cancel all subsequent output display
!! @details As referred to in the description of disp_out(), output may not be desired
!! when the input is read in some cases.  This function turns off subsequent output
!! display.
!! @param[in,out] this  Input file object
!<----------------------------------------------------------------------------------------
subroutine no_disp_out( this )

  implicit none

  class( t_input_file ), intent(inout) :: this

  this%display_out = .false.

end subroutine no_disp_out
!-------------------------------------------------------------------------------

end module m_input_file
