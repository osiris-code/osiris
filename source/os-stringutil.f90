module stringutil

implicit none

character, parameter :: p_null = achar(0)
character, parameter :: p_tab = achar(9)
character, parameter :: p_space = ' '

interface lowercase
  module procedure lowercase
end interface

interface is_blank
  module procedure is_blank
end interface

interface is_digit
  module procedure is_digit
end interface

interface strtoint
  module procedure strtoint
end interface

interface strtodouble
  module procedure strtodouble
end interface

interface scanall
  module procedure scanall
end interface

interface parsenumber
  module procedure parsenumber
end interface

interface tostring
  module procedure tostring_int
  module procedure tostring_r4
  module procedure tostring_r8
end interface

interface idx_string
  module procedure idx_string
end interface

interface replace_blanks
  module procedure replace_blanks
end interface

interface time2seconds
  module procedure time2seconds
end interface

contains

!-------------------------------------------------------------------------------
! Converts the time string to seconds. timestr must be in the form  
! 's','m:s' or 'h:m:s'. If conversion fails the routine returns -1.0
!-------------------------------------------------------------------------------
function time2seconds( timestr )

  implicit none
  
  character(len=*), intent(in) :: timestr
  integer :: time2seconds
  
  integer, dimension(3) :: t
  
  integer :: i, l, p, ierr
  
  t(1) = 0
  t(2) = 0
  t(3) = 0
  
  p = 1
  l = len(timestr)
  
  do i = 1, 3
    p = scan( timestr(p:l), ':', back = .true. )
    read( timestr(p+1:l), '(i2)', iostat = ierr ) t(i)
    if ( ierr /= 0 ) then
      time2seconds = -1
      return
    endif
    if ( p == 0 ) then
      exit
    else
      l = p - 1
      p = 1
    endif
  enddo
 
  time2seconds = t(1) + 60*t(2) + 3600*t(3)

end function time2seconds
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
function lowercase(str)
!-------------------------------------------------------------------------------
! return the lowercase character corresponding to c
!-------------------------------------------------------------------------------
  implicit none
  
  character( len = * ), intent(in) :: str
  character( len = len(str) )  :: lowercase
  integer, parameter :: ichar_lcshift = iachar('a') - iachar('A')
  
  integer :: i
  character :: c
  
  do i = 1, len(str)
    c = str(i:i) 
    if (((c >= 'A') .and. (c <= 'Z'))) then
       lowercase(i:i) = achar( ichar(c) + ichar_lcshift )
    else
       lowercase(i:i) = c 
    endif
  enddo

end function lowercase
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
function is_blank(c)
!-------------------------------------------------------------------------------
! test if character is blank (white space) 
!-------------------------------------------------------------------------------
  implicit none
  
  character, intent(in) :: c
  logical :: is_blank
  
  if ((c == p_space) .or. (c == p_tab)) then
     is_blank = .true. 
  else
     is_blank = .false. 
  endif
  
end function is_blank
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function is_alpha(c)
!-------------------------------------------------------------------------------
! test if character is a letter 
!-------------------------------------------------------------------------------
  implicit none
  
  character, intent(in) :: c
  logical :: is_alpha

!  if (((c >= 'a') .and. (c <= 'z')) .or. &
!     ((c >= 'A') .and. (c <= 'Z'))) then
!    is_alpha = .true. 
!  else
!    is_alpha = .false. 
!  endif

  is_alpha = ((c >= 'a') .and. (c <= 'z')) .or. ((c >= 'A') .and. (c <= 'Z'))
  
end function is_alpha
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function is_digit(c)
!-------------------------------------------------------------------------------
! test if character is a digit 
!-------------------------------------------------------------------------------
  implicit none
  
  character, intent(in) :: c
  logical :: is_digit
  
!  if ((c >= '0') .and. (c <= '9')) then
!    is_digit = .true. 
!  else
!    is_digit = .false. 
!  endif
  
  is_digit = (c >= '0') .and. (c <= '9')

end function is_digit
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function strtoint(str, ierr)
!-------------------------------------------------------------------------------
! convert string to integer
!-------------------------------------------------------------------------------
  implicit none
  
  character(len=*), intent(in) :: str
  integer, intent(out) :: ierr
  integer :: strtoint
  
  read( str, *, iostat = ierr ) strtoint

end function strtoint
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function strtodouble(str, ierr)
!-------------------------------------------------------------------------------
! convert string to integer
!-------------------------------------------------------------------------------
  implicit none
  
  character(len=*), intent(in) :: str
  integer, intent(out) :: ierr
  real( kind(1.0d0) ) :: strtodouble
  
  read( str, *, iostat = ierr ) strtodouble
  
end function strtodouble
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine scanall( str, c, pos, nfound )
!-------------------------------------------------------------------------------
! scans for all occurences of character c in string str
!-------------------------------------------------------------------------------

  implicit none
  
  character(len=*), intent(in) :: str
  character, intent(in) :: c
  
  integer, dimension(:), intent(out) :: pos
  integer, intent(out) :: nfound
  
  integer :: idx, npos
  
  nfound = 0
  pos = 0
  
  idx = 1
  
  do
    if ( idx > len( str ) ) exit
    
    npos = scan( str(idx:), c) 
    if ( npos == 0 ) exit
    
    nfound = nfound+1
    if ( nfound <= size(pos) ) pos( nfound ) = idx + npos - 1

    idx = idx + npos 
  enddo


end subroutine scanall
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine parsenumber( str, numberType, numberLen )
!-------------------------------------------------------------------------------
! parses the supplied string looking for a valid numeric value. If found
! returns the length of the string holding the number and the number type (1 - integer,
! 2 - float), otherwise returns length = 0, and type = -1
!-------------------------------------------------------------------------------

  implicit none
  
  character(len=*), intent(in) :: str
  integer, intent(out) :: numberLen, numberType

  character(len=len_trim(str)) :: value_text
  character :: c, c2
  integer :: pos
  logical :: has_decpoint, has_value, has_exponent, finished

  pos = 0
  has_decpoint = .false.
  has_value    = .false.
  has_exponent = .false.
  finished     = .false. 
  value_text = ''
  
  ! check for signed number
  c = str(1:1)
  if ( c == '+' .or. c== '-' ) pos = pos+1
  
  do
    pos = pos + 1
    if ( pos > len( str ) ) then
      finished = .true.
      exit
    endif
    c = str(pos:pos)
    
    ! check for decimal point
    if ( c == '.' ) then
      if ( .not. has_decpoint ) then
         has_decpoint = .true.    
      else
         ! error, two decimal points found
         finished = .true. 
      endif
    
    else if ( is_digit(c) ) then
      has_value = .true.
    
    else if ( c=='e' .or. c=='E' .or. c=='d' .or. c=='D' ) then 
      if ( (.not. has_value) .or. has_exponent ) then
         ! of exponent found before mantissa, or exponent symbol found inside exponent
         finished = .true.    
      endif      
      has_exponent = .true.
      has_decpoint = .true. ! decimal points are not allowed in the exponent  
    
      ! check if next character is + or -
      pos = pos+1
      if ( pos > len( str ) ) then
        finished = .true.
        exit
      endif
      c2 = str( pos:pos )
      if ((c2 == '+') .or. (c2 == '-')) then
        value_text = trim(value_text)//c//c2
        cycle
      else
        pos = pos - 1
      endif
      
      has_value = .false. ! an integer must follow 
    else 
      ! invalid character found
      ! finish processing
      
      finished = .true.
    endif
    
    
    if ( finished ) then
      pos = pos-1 
      exit
    else
      value_text = trim(value_text)//c
    endif

  enddo

  if ( .not. has_value ) then
    numberLen = 0
    numberType = -1
    return
  endif

  numberLen = len_trim(value_text)
  if ( has_exponent .or. has_decpoint ) then
    numberType = 2
  else
    numberType = 1
  endif
  
  
!  print *, ">", trim(value_text), "< len = ", numberLen, ' type = ', numberType

end subroutine parseNumber
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
pure function tostring_int( i )
!-------------------------------------------------------------------------------
! Converts the given integer to a left justified string
!-------------------------------------------------------------------------------
  
  implicit none
  
  integer, intent(in) :: i

  character(len=15) :: tostring_int

  write( tostring_int, * ) i  
  
  ! discard leading blanks
  tostring_int = adjustl(tostring_int)
    
end function tostring_int
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
pure function tostring_r4( r4 )
!-------------------------------------------------------------------------------
! Converts the given real to a left justified string
!-------------------------------------------------------------------------------
  
  implicit none
  
  real( kind(1.0e0) ), intent(in) :: r4

  character(len=24) :: tostring_r4

  write( tostring_r4, * ) r4  
  
  ! discard leading blanks
  tostring_r4 = adjustl(tostring_r4)
    
end function tostring_r4
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
pure function tostring_r8( r8 )
!-------------------------------------------------------------------------------
! Converts the given real to a left justified string
!-------------------------------------------------------------------------------
  
  implicit none
  
  real( kind(1.0d0) ), intent(in) :: r8

  character(len=24) :: tostring_r8

  write( tostring_r8, * ) r8  
  
  ! discard leading blanks
  tostring_r8 = adjustl(tostring_r8)
    
end function tostring_r8
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function idx_string( i, digits )
!-------------------------------------------------------------------------------
! this function generates a string representing the integer i using the number
! of digits specified (i.e. adding leading zeroes) 
!-------------------------------------------------------------------------------
  implicit none
  
  ! arguments and return value
  
  integer, intent(in) :: i, digits  
  character(len = digits) :: idx_string 
  
  ! local variables
  character(len = 18) :: d
  
  d = trim(tostring(digits))
  write( idx_string, '(i'//d//'.'//d//')' ) i
                       
end function idx_string
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function find_char_unquote( in_str, in_char )
!-------------------------------------------------------------------------------
! find the position of the first occurence of the character in_char in the 
! string in_str that is not inside quotes (single or double)
!-------------------------------------------------------------------------------
    
  implicit none
  
  character ( len = * ), intent(in) :: in_str
  character, intent(in) :: in_char
  integer :: find_char_unquote 

  integer :: i
  logical :: in_quotes
  character :: quote_char
  
  in_quotes = .false.
  find_char_unquote = 0
   do i = 1, len(in_str)
    if (in_quotes) then
      if (in_str(i:i) .eq. quote_char) in_quotes = .false.
    elseif (in_str(i:i) == '"' .or. in_str(i:i) == "'") then
      in_quotes = .true.
      quote_char= in_str(i:i)
    else if (in_str(i:i) .eq. in_char) then
      find_char_unquote = i
      return
    endif
    
  enddo

end function find_char_unquote
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function replace_blanks( str, subst_char )
!-------------------------------------------------------------------------------
!  returns a version of string str where blanks (spaces and tabs) are replaced
!  by a given character (defaults to underscore)
!-------------------------------------------------------------------------------

  implicit none

  !       dummy variables
  
  character(len=*), intent(in) :: str
  character, intent(in), optional :: subst_char
  character(len = len(str)) :: replace_blanks
  
  !       local variables
  character :: schar
  integer :: i
                
  !       executable statements
  
  if (present(subst_char)) then
    schar = subst_char
  else
    schar = '_'
  endif
  
  do i=1, len(str)
    if( is_blank(str(i:i)) ) then
       replace_blanks(i:i) = schar
    else
       replace_blanks(i:i) = str(i:i)
    endif
  enddo

end function replace_blanks
!-------------------------------------------------------------------------------


end module stringutil

