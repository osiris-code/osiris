!
! Function parsing class
! 

! There is a known issue with variable names: if a variable name ends in digit+'e' (or 'd') it 
! i.e. 'a2e' may be interpreted as a numeric constant and evaluating something like var +/- number
! i.e. 'a2e - 1' may give the wrong result.

module m_fparser

  use stringutil
  use m_system
  
  implicit none

! restrict access to things explicitly declared public
  private

  ! precision to use for the calculations
  !integer, parameter, public :: p_k_fparse = p_single
  integer, parameter, public :: p_k_fparse = p_double
  
  
  ! predefined sizes
  integer, parameter :: p_max_expr_len    = 1024 ! maximum length an expression can have
  integer, parameter :: MAX_FUNC_LEN    = 16   ! maximum length a single function name can have
  integer, parameter :: MAX_VAR_LEN     = 128  ! maximum length a variable name can have
  integer, parameter :: MAX_VARS        = 16   ! maximum number of variables 
  integer, parameter :: MAX_NESTED_FUNC = 16   ! maximum number of nested functions allowed
  integer, parameter :: MAX_OP_LEN      = 2    ! maximum length an operator can have

  ! predefined bytecode and stack sizes 
  integer, parameter :: MAX_CODE_SIZE = 256    ! code will be limited to 256 expressions
  integer, parameter :: MAX_STACK_SIZE = 256   ! maximum number of values in the calculation stack 
  integer, parameter :: MAX_DATA_SIZE = 64     ! maximum number of values in the data stack

  ! error codes
  integer, parameter :: err_eof     = -1
  integer, parameter :: err_mult_op = -2
  integer, parameter :: err_inv_num = -3
  integer, parameter :: err_parent  = -4
  integer, parameter :: err_syntax  = -5
  integer, parameter :: err_fun_par = -6
  integer, parameter :: err_param   = -7
  
  
  ! List of valid operators:
  ! - operators are specified in groups of operations with the same priority, separated by
  !   blank ('  '), elements.
  ! - operator groups are ordered with increasing priority.
  ! - Note that if two operators within the same group begin with the same character and
  !   have different lengths (e.g. '<=' and '<'), the longest one MUST be specified first
  !   because of the implementation of the scan_ops function  
  !
  ! NOTE: Operators must not have alphanumeric characters or dots
    
  character(len = *), dimension(19), parameter :: oper = &
     (/ '||', '  ', &
        '&&', '  ', &
        '==', '!=', '  ', &
        '<=', '>=', '< ', '> ', '  ', &
        '+ ', '- ', '  ', &
        '* ', '/ ', '  ', &
        '^ ' /)
  
  ! NOTE2: operator ^ is used for a ^ int(b) i.e. b is converted to an integer before
  !   evaluating. If you want b to not be an integer use the function pow(a,b) instead but
  !   keep in mind that if a < 0 pow(a,b) returns NaN
  
  ! Operator op codes
  ! These must be in sequence with the same order as the operator list array
  
  integer, parameter :: op_or   =  1, &
                        op_and  =  2, &
                        op_eq   =  3, op_ne   =  4, &
                        op_le   =  5, op_ge   =  6, op_lt   =  7, op_gt   =  8, &
                        op_add  =  9, op_sub  = 10, &
                        op_mult = 11, op_div  = 12, &
                        op_pow  = 13
  
  ! Total number of operators
  integer, parameter :: n_oper = 13

  integer, parameter :: bcode_data = 0

  ! valid functions
  integer, parameter :: n_functions = 30
  
  type :: t_math_func
    character (len = MAX_FUNC_LEN) :: name
    integer                        :: n_params
  end type t_math_func

  ! note that if two functions begin with the same string
  ! the longest one MUST be specified first because of the
  ! implementation of the function if_function, e.g. log10 and log, 
  ! or atan2 and atan
  
  type (t_math_func), dimension(n_functions) , parameter :: functions = &
                                (/ t_math_func('abs             ', 1) , &   ! 1
                                   t_math_func('sinh            ', 1) , &   ! 2
                                   t_math_func('sin             ', 1) , &   ! 3
                                   t_math_func('cosh            ', 1) , &   ! 4
                                   t_math_func('cos             ', 1) , &   ! 5
                                   t_math_func('tanh            ', 1) , &   ! 6
                                   t_math_func('tan             ', 1) , &   ! 7
                                   t_math_func('exp             ', 1) , &   ! 8
                                   t_math_func('log10           ', 1) , &   ! 9
                                   t_math_func('log             ', 1) , &   ! 10
                                   t_math_func('asin            ', 1) , &   ! 11
                                   t_math_func('acos            ', 1) , &   ! 12
                                   t_math_func('atan2           ', 2) , &   ! 13
                                   t_math_func('atan            ', 1) , &   ! 14
                                   t_math_func('sqrt            ', 1) , &   ! 15
                                   t_math_func('not             ', 1) , &   ! 16
                                   t_math_func('neg             ', 1) , &   ! 17
                                   t_math_func('if              ', 3) , &   ! 18
                                   t_math_func('pow             ', 2) , &   ! 19
                                   t_math_func('int             ', 1) , &   ! 20
                                   t_math_func('nint            ', 1) , &   ! 21
                                   t_math_func('ceiling         ', 1) , &   ! 22
                                   t_math_func('floor           ', 1) , &   ! 23
                                   t_math_func('modulo          ', 2) , &   ! 24
                                   t_math_func('rect            ', 1) , &   ! 25
                                   t_math_func('step            ', 1) , &   ! 26
                                   t_math_func('min3            ', 3) , &   ! 27
                                   t_math_func('min             ', 2) , &   ! 28
                                   t_math_func('max3            ', 3) , &   ! 29
                                   t_math_func('max             ', 2)  /)   ! 30

  integer, parameter :: func_abs     = 1, &
                        func_sinh    = 2, &
                        func_sin     = 3, &
                        func_cosh    = 4, &
                        func_cos     = 5, &
                        func_tanh    = 6, &
                        func_tan     = 7, &
                        func_exp     = 8, &
                        func_log10   = 9, &
                        func_log     = 10, &
                        func_asin    = 11, &
                        func_acos    = 12, &
                        func_atan2   = 13, &
                        func_atan    = 14, &
                        func_sqrt    = 15, &
                        func_not     = 16, &
                        func_neg     = 17, &
                        func_if      = 18, &
                        func_pow     = 19, &
                        func_int     = 20, &
                        func_nint    = 21, &
                        func_ceiling = 22, &
                        func_floor   = 23, &
                        func_modulo  = 24, &
                        func_rect    = 25, &
                        func_step    = 26, &
                        func_min3    = 27, &
                        func_min     = 28, &
                        func_max3    = 29, &
                        func_max     = 30

     
  ! beginning of variable space code
  integer, parameter :: bcode_var = n_oper + n_functions
  
  ! class definition
  type :: t_fparser
  
     integer :: code_size        ! actual size of compiled code
     integer :: stack_size       ! actual size of calculation stack
     integer :: data_size        ! actual size of data stack
     
     integer :: stack_ptr        ! pointer to the current value on the stack
                                 ! used at compile time to determine required stack size
                                 ! since we are not allocating the stack dynamically it 
                                 ! serves no purpose
     
     ! compiled code
     integer, dimension(MAX_CODE_SIZE) :: bytecode  
     
     ! stack for temporary values during calculation
     real(p_k_fparse), dimension(MAX_STACK_SIZE):: stack     
     
     ! internal stack for compiled number constants
     real(p_k_fparse), dimension(MAX_DATA_SIZE):: data_stack 
     
     ! Expression text
     character (len = p_max_expr_len) :: func
     
     ! Number of variables
     integer :: nvars
     
     ! Variable names
     character (len = MAX_VAR_LEN), dimension(MAX_VARS) :: var_list
  
  end type t_fparser
  
  interface setup
    module procedure setup_fparser
  end interface
  
  interface cleanup
    module procedure cleanup_fparser
  end interface

  interface eval
    module procedure eval_fparser
  end interface
  
  interface functext
    module procedure functext_parser
  end interface
  
! declare things that should be public 
  
  public :: t_fparser, setup, cleanup, eval, p_max_expr_len, functext

  contains
  
  function functext_parser( this )
    implicit none
     type( t_fparser ),                   intent(in) :: this
    character( len = len_trim( this%func ) ) :: functext_parser
    
    functext_parser = trim( this%func )
  end function functext_parser
  
!-----------------------------------------------------------------------------------------
!   class constructor
!-----------------------------------------------------------------------------------------
  subroutine setup_fparser( this, str_expr, var_list, ierr )
!-----------------------------------------------------------------------------------------
  
    implicit none
    
    ! dummy variables
    
    type( t_fparser ),                   intent(inout) :: this
    character ( len = * ),               intent (in)   :: str_expr
    character ( len = * ), dimension(:), intent (in)   :: var_list
    integer,                             intent (out)  :: ierr 
        
    ! executable statements
    
    
    ierr = 0
    
    ! initialize local data
    
    this%func = remove_spaces(str_expr) 
    this%nvars = size(var_list) 
    this%var_list(1:this%nvars) = var_list
    
    ! not necessary
    if (this%nvars < MAX_VARS) this%var_list(this%nvars+1:) = '~'
    
    ! check expression syntax
    
    call check_syntax_fparser( this , ierr )
    if (ierr < 0) return   
    
    ! compile the function

    call compile_fparser( this )

  end subroutine setup_fparser
!-----------------------------------------------------------------------------------------

subroutine cleanup_fparser( this )

  implicit none
  
  type( t_fparser ),                   intent(inout) :: this
  
  this%func = ''
  this%nvars = 0
  this%var_list = ''
  
  this%code_size = 0
  
end subroutine cleanup_fparser

!-----------------------------------------------------------------------------------------
!   process bytecode and evaluate function
!-----------------------------------------------------------------------------------------
  function eval_fparser( this, var_values )
!-----------------------------------------------------------------------------------------
  
    implicit none
    
    ! dummy variables
    
    real( p_k_fparse )                                  :: eval_fparser
    type( t_fparser ),                   intent(inout) :: this
    real( p_k_fparse ), dimension(:),     intent(in)    :: var_values
    
    ! local variables
    
    integer :: ip, & ! instruction pointer
               sp, & ! stack pointer
               dp    ! data pointer
    
    integer :: exponent   ! workarround for bug in a**b
    
    ! executable statements

    sp = 0; dp = 1

    do ip = 1, this%code_size
      select case (this%bytecode(ip))
        
        case (bcode_data)             ! compiled value, grab it from the data stack
          sp = sp + 1; this%stack(sp) = this%data_stack(dp); dp = dp + 1
          
        case ( op_eq )
          if ( this%stack(sp-1) == this%stack(sp) ) then
            this%stack(sp-1) = 1.0
          else 
            this%stack(sp-1) = 0.0
          endif 
          sp = sp-1
        
        case ( op_ne )
          if ( this%stack(sp-1) /= this%stack(sp) ) then
            this%stack(sp-1) = 1.0
          else 
            this%stack(sp-1) = 0.0
          endif 
          sp = sp-1
        
        case ( op_ge )
          if ( this%stack(sp-1) >= this%stack(sp) ) then
            this%stack(sp-1) = 1.0
          else 
            this%stack(sp-1) = 0.0
          endif 
          sp = sp-1
        
        case ( op_le )
          if ( this%stack(sp-1) <= this%stack(sp) ) then
            this%stack(sp-1) = 1.0
          else 
            this%stack(sp-1) = 0.0
          endif 
          sp = sp-1
        
        case ( op_gt )
          if ( this%stack(sp-1) > this%stack(sp) ) then
            this%stack(sp-1) = 1.0
          else 
            this%stack(sp-1) = 0.0
          endif 
          sp = sp-1
          
        case ( op_lt )
          if ( this%stack(sp-1) < this%stack(sp) ) then
            this%stack(sp-1) = 1.0
          else 
            this%stack(sp-1) = 0.0
          endif 
          sp = sp-1
          
        case ( op_and )
          if ( this%stack(sp-1) /= 0.0 .and. this%stack(sp) /= 0.0 ) then
            this%stack(sp-1) = 1.0
          else 
            this%stack(sp-1) = 0.0
          endif 
          sp = sp-1
        
        case ( op_or )
          if ( this%stack(sp-1) /= 0.0 .or. this%stack(sp) /= 0.0 ) then
            this%stack(sp-1) = 1.0
          else 
            this%stack(sp-1) = 0.0
          endif 
          sp = sp-1
          
        case ( op_add )
           this%stack(sp-1) = this%stack(sp-1) + this%stack(sp); sp = sp -1
           
        case ( op_sub )
           this%stack(sp-1) = this%stack(sp-1) - this%stack(sp); sp = sp -1
           
        case ( op_mult )
            this%stack(sp-1) = this%stack(sp-1) * this%stack(sp); sp = sp -1
           
        case ( op_div )
            this%stack(sp-1) = this%stack(sp-1) / this%stack(sp); sp = sp -1
            
        case ( op_pow )
            exponent = int(this%stack(sp))
            sp = sp -1
            this%stack(sp) = this%stack(sp) ** exponent
            
        case ( n_oper + func_abs )
            this%stack(sp) = abs( this%stack(sp) ) 
        
        case ( n_oper + func_sin )
            this%stack(sp) = sin( this%stack(sp) )

        case ( n_oper + func_sinh )
            this%stack(sp) = sinh( this%stack(sp) )
        
        case ( n_oper + func_cos )
            this%stack(sp) = cos( this%stack(sp) )

        case ( n_oper + func_cosh )
            this%stack(sp) = cosh( this%stack(sp) )
        
        case ( n_oper + func_tan )
            this%stack(sp) = tan( this%stack(sp) )

        case ( n_oper + func_tanh )
            this%stack(sp) = tanh( this%stack(sp) )
            
        case ( n_oper + func_exp )
            this%stack(sp) = exp( this%stack(sp) )

        case ( n_oper + func_log10 )
            this%stack(sp) = log10( this%stack(sp) )
            
        case ( n_oper + func_log )
            this%stack(sp) = log( this%stack(sp) )
                        
        case ( n_oper + func_asin )
            this%stack(sp) = asin( this%stack(sp) )
             
        case ( n_oper + func_acos  )
            this%stack(sp) = acos( this%stack(sp) )
            
        case ( n_oper + func_atan2 )
            this%stack(sp-1) = atan2(this%stack(sp-1),this%stack(sp)); sp = sp -1
            
        case ( n_oper + func_atan )
             this%stack(sp) = atan( this%stack(sp) )
             
        case ( n_oper + func_sqrt )
             this%stack(sp) = sqrt( this%stack(sp) )
             
        case ( n_oper + func_not )
             if (this%stack(sp) == 0.0) then
               this%stack(sp) = 1.0
             else
               this%stack(sp) = 0.0
             endif
        case ( n_oper + func_neg )
             this%stack(sp) = - this%stack(sp) 
        
        case ( n_oper + func_if )
             if (this%stack(sp-2) /= 0.0) then
                this%stack(sp-2) = this%stack(sp-1)
             else
                this%stack(sp-2) = this%stack(sp)
             endif
             sp = sp - 2 
             
        case ( n_oper + func_pow )
            this%stack(sp-1) = this%stack(sp-1) ** this%stack(sp); sp = sp -1
 
        case ( n_oper + func_int )
            this%stack(sp) = int( this%stack(sp) )

        case ( n_oper + func_nint )
            this%stack(sp) = nint( this%stack(sp) )

        case ( n_oper + func_ceiling )
            this%stack(sp) = ceiling( this%stack(sp) )

        case ( n_oper + func_floor )
            this%stack(sp) = floor( this%stack(sp) )

        case ( n_oper + func_modulo )
            this%stack(sp-1) = modulo(this%stack(sp-1),this%stack(sp)); sp = sp -1

        case ( n_oper + func_rect )
            if ( abs(this%stack(sp)) <= 0.5 ) then
              this%stack(sp) = 1.0
            else
              this%stack(sp) = 0.0
            endif

        case ( n_oper + func_step )
            if ( this%stack(sp) >= 0.0 ) then
              this%stack(sp) = 1.0
            else
              this%stack(sp) = 0.0
            endif

        case ( n_oper + func_min3 )
            this%stack(sp-2) = &
              min(this%stack(sp-2),this%stack(sp-1),this%stack(sp))
            sp = sp -2

        case ( n_oper + func_min )
            this%stack(sp-1) = &
              min(this%stack(sp-1),this%stack(sp)); sp = sp -1

        case ( n_oper + func_max3 )
            this%stack(sp-2) = &
              max(this%stack(sp-2),this%stack(sp-1),this%stack(sp))
            sp = sp -2

        case ( n_oper + func_max )
            this%stack(sp-1) = &
              max(this%stack(sp-1),this%stack(sp)); sp = sp -1

        case default ! this is only reached for variables
           sp = sp + 1
           this%stack(sp) = var_values(this%bytecode(ip) - bcode_var )
      end select
    enddo
    
    
    eval_fparser = this%stack(1)

  end function eval_fparser
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!   compile the supplied function
!-----------------------------------------------------------------------------------------
  subroutine compile_fparser( this )
!-----------------------------------------------------------------------------------------

    implicit none
    
    ! dummy variables
    
    type( t_fparser ),                   intent(inout) :: this
        
    ! local variables
                           
    this%code_size = 0
    this%stack_size = 0
    this%data_size = 0
    
    this%stack_ptr = 0
    
    call compile_substr_fparser( this, trim(this%func) )
    
    ! (*debug*)
    !print *, "compiled code size = ", this%code_size
    !print *, "required stack space = ", this%stack_size
    !print *, "required data size = ", this%data_size
    !print *, "code = ", this%bytecode(1:this%code_size)
        
  end subroutine compile_fparser
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! recursively compile substring
!-----------------------------------------------------------------------------------------
recursive subroutine compile_substr_fparser( this, sub_str )
!-----------------------------------------------------------------------------------------    
  implicit none
  
  ! dummy variables
  
  type( t_fparser ),                   intent(inout) :: this
  character ( len = * ),                 intent (in) :: sub_str 
      
  ! local variables
  
  integer :: sub_str_len            ! length of the sub_str begin processed
  integer :: fidx, f_len, n_params  ! function index, function name length, number of params
  integer :: varidx, var_len
  
  integer :: i, j, & ! string position index
             p, &    ! parameter count
             k       ! parenthesis count
  
  integer :: op_code, &      ! offset for calculating op code
             n_op, &         ! number of operators in group
             op_len, &       ! length of operator found
             list_idx, &     ! position in the operator list
             op_idx          ! position of the operator in the operator group
  
  integer :: ierr
  
  ! print *, '>', trim(sub_str), '<'
  
  sub_str_len = len(sub_str)
        
  ! check for special cases of substring
    
  if (sub_str(1:1) == '+') then                      ! case 1: str = '+...'
    call compile_substr_fparser( this, sub_str(2:) ) ! no compilation required
    return
   
  elseif (if_enclosed_parent(sub_str)) then          ! case 2: str = '(...)'
    call compile_substr_fparser( this, &             ! no compilation required
               sub_str(2:sub_str_len-1) )
    return
  
  elseif (is_alpha(sub_str(1:1))) then
    fidx = is_function( sub_str , f_len, n_params)
    if (fidx > 0) then                               ! case 3: str = 'fcn(...)'
      if (if_enclosed_parent(sub_str(f_len+1:))) then
        i = f_len + 2                                ! compile the parameters first
        
        do p = 2, functions(fidx)%n_params           ! parse multiple parameters
          j = comma_pos(sub_str(i:))
          
          call compile_substr_fparser( this, &         
               sub_str(i:i+j-2) )  
          i = i + j 
        enddo
        
        call compile_substr_fparser( this, &         
               sub_str(i:sub_str_len-1) )
        
        ! add the compiled function to the code
        call add_function( this, fidx )

        ! decrease stack pointer if n_params > 1
        this%stack_ptr = this%stack_ptr - functions(fidx)%n_params + 1
        
        return
      endif
    endif
    
  elseif (sub_str(1:1) == '-') then
    if (if_enclosed_parent(sub_str(2:))) then        ! case 4: str = '-(...)'
                                                     ! the same as neg(...)
        call compile_substr_fparser( this, &         ! compile the parameters first
              sub_str(3:sub_str_len-1) )

        ! add the compiled function to the code
        call add_function( this, func_neg )

        return
      
    elseif (is_alpha(sub_str(2:2))) then
        
      fidx = is_function( sub_str(2:) , f_len, n_params)
    
      if (fidx > 0) then                              ! case 5: str = '-fcn(...)'
        if (if_enclosed_parent(sub_str(2+f_len:))) then
          i = f_len + 3                               ! compile the parameters first

          do p = 2, functions(fidx)%n_params          ! parse multiple parameters
             j = index(sub_str(i:),',')
             call compile_substr_fparser( this, &         
                sub_str(i:i+j-2) )  
             i = i + j 
          enddo
          
          call compile_substr_fparser( this, &         
                 sub_str(i:sub_str_len-1) )
        
          ! add the compiled function to the code followed by neg()
          call add_function( this, fidx )
          call add_function( this, func_neg )
          
          return
        endif
      endif          
        
    endif
  endif
  
  ! parse operators
  ! only base level (k == 0) will be processed
  ! the remaining levels are processed recursively
  ! Operators are evaluated from left to right i.e. a - b - c ==> (a - b) - c
  
  ! Loop over operator groups
  op_code = 0
  list_idx = 1
  do while ( list_idx <= size( oper ) )
    
    n_op = same_priority_ops( list_idx )
    
    k = 0
    
    ! Process string from end to beginning
    do i=sub_str_len,2,-1  
      
      ! skip sections inside parenthesis
      select case (sub_str(i:i))
        case (')')
          k = k+1
        case ('(')
          k = k-1
        case default
          continue
      end select
      
      if (k == 0) then 
         
         op_idx = scan_ops( sub_str, i, list_idx, n_op, op_len )  
         
         if ( op_idx > 0 ) then
           
		   if ( num_signed_exp( i, sub_str ) ) then  
			   ! number with signed exponent e.g. 1.0e-5
			   ! this is a more complicated test because a variable name can end in
			   ! e or d
			 
			   ! ignore
			   continue
		   elseif ( (sub_str(1:1) == '-') .and. ( op_code >= op_div ) ) then
			  ! This is the only monadic operator we have
			  ! case 6: str = - ... op ... with op with higher priority than * or / (i.e. ^)
			  ! compile as neg(... op ...)
			
			  ! compile the left operand
			  call compile_substr_fparser( this, sub_str(2:i-1) )
			
			  ! compile the right operand
			  call compile_substr_fparser( this, sub_str(i+op_len:) )
			
			  ! add the compiled operator to the code
			  call add_operator( this, op_code + op_idx )
			
			  ! decrease stack pointer
			  ! because the two values are converted into a single one
			  this%stack_ptr = this%stack_ptr - 1
			
			  call add_function( this, func_neg )
			
			  return
			
		   else                                           ! case 7: str = ... op ...
			
			  ! compile the left operand
			  call compile_substr_fparser( this, sub_str(1:i-1) )
			
			  ! compile the right operand
			  call compile_substr_fparser( this, sub_str(i+op_len:) )
			
			  ! add the compiled operator to the code
			  call add_operator( this, op_code + op_idx )
			
			  ! decrease stack pointer
			  ! because the two values are converted into a single one
			  this%stack_ptr = this%stack_ptr - 1
			
			  return
		   endif 
         
         endif
      endif
      
    enddo

	! Process next group of operators
	list_idx = list_idx + n_op + 1
	op_code = op_code + n_op
  enddo
  
  ! parse numbers and variables
  
  i = 1
  
  ! check for var and -var
  if (sub_str(1:1) == '-') then  
    varidx = is_variable_fparser( this, sub_str(2:), var_len)
    if (varidx > 0) i = i + 1 
  else   
    varidx = is_variable_fparser( this, sub_str(1:), var_len) 
  endif

  ! add the compiled value
  if (varidx > 0) then ! variable
    call add_variable( this, varidx )
  else                                              ! number
    call add_constant( this, strtodouble(sub_str, ierr) )
  endif
  
  ! increase stack pointer because one value will be added to the stack
  this%stack_ptr = this%stack_ptr + 1

  ! find maximum stack size
  if (this%stack_ptr > this%stack_size) this%stack_size = this%stack_ptr
  
  ! if -var found add a neg() function
  if (i>1) call add_function( this, func_neg )

  contains 
  
  !-------------------------------------------------
  !  number with signed exponent e.g. 1.0e-5
  ! this is a more complicated test because a variable name can end in
  !  e or d
  !-------------------------------------------------
  function num_signed_exp( i, sub_str )
    
    implicit none
    
    integer, intent(in) :: i
    character(len=*), intent(in) :: sub_str
    logical :: num_signed_exp
    
    num_signed_exp = .false.
    if ( i > 2 ) then
       num_signed_exp = ((sub_str(i:i)=='-') .or. (sub_str(i:i)=='+')) .and. &
                        (scan( sub_str(i-1:i-1), 'edED' ) > 0) .and. &
                        (scan( sub_str(i-2:i-2), '0123456789.') > 0)
    endif
  end function num_signed_exp
  !-------------------------------------------------

  !-------------------------------------------------
  ! Check if str begins with any of the strings in oplist, and returns which one
  !-------------------------------------------------
  function scan_ops( in_str, i, list_idx, n_op, op_len )
     
    implicit none
    
    character(len = *), intent(in) :: in_str
    integer, intent(in) :: i, list_idx, n_op
    integer, intent(out) :: op_len 
    
    integer :: scan_ops, j, l
    
    do scan_ops = 1, n_op
      j = list_idx + scan_ops - 1
      op_len = len_trim( oper( j ) )
      
      l = min( op_len, len(in_str) - i + 1)
	  if (in_str(i:i+l-1) == oper(j)) return
    enddo

    op_len = 0
    scan_ops = -1
    
  end function scan_ops
  !-------------------------------------------------
  
  !-------------------------------------------------
  ! count operators with same priority
  !-------------------------------------------------
  function same_priority_ops( i )
    
    implicit none
    
    integer, intent(in) :: i
    integer :: same_priority_ops
    
    same_priority_ops = 0
    do while ( i + same_priority_ops <= size( oper ) )
      if ( oper(i + same_priority_ops) == '  ' ) exit
      same_priority_ops = same_priority_ops + 1  
    enddo
  
  end function same_priority_ops
  !-------------------------------------------------

  
end subroutine compile_substr_fparser
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! add bytecode for a variable to the code stack
!-----------------------------------------------------------------------------------------
  subroutine add_variable( this, varidx ) 
!-----------------------------------------------------------------------------------------
    implicit none

    ! dummy variables
    
    type( t_fparser ), intent(inout) :: this
    integer,           intent(in) :: varidx
    
    ! local variables
    
    ! print *, "add variable = ", varidx
    
    this%code_size = this%code_size + 1
    this%bytecode(this%code_size) = bcode_var + varidx

  end subroutine add_variable
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! add bytecode for a function to the code stack
!-----------------------------------------------------------------------------------------
  subroutine add_function( this, fidx ) 
!-----------------------------------------------------------------------------------------
    implicit none

    ! dummy variables
    
    type( t_fparser ), intent(inout) :: this
    integer,           intent(in)    :: fidx
    
    ! local variables
    
    ! print *, "add func = '", trim(functions(fidx)%name),"'"
    
    this%code_size = this%code_size + 1
    this%bytecode(this%code_size)  = n_oper + fidx

  end subroutine add_function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! add bytecode for an operator to the code stack
!-----------------------------------------------------------------------------------------
  subroutine add_operator( this, opidx ) 
!-----------------------------------------------------------------------------------------
    implicit none

    ! dummy variables
    
    type( t_fparser ), intent(inout) :: this
    integer,           intent(in)    :: opidx
    
    ! local variables
    
    ! print *, "add op = ",opidx
    
    this%code_size = this%code_size + 1
    this%bytecode(this%code_size) = opidx

  end subroutine add_operator
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! add bytecode for a number to the code stack
! and stores number in the intern_float stack
!-----------------------------------------------------------------------------------------
  subroutine add_constant( this, float ) 
!-----------------------------------------------------------------------------------------
    implicit none

    ! dummy variables
    
    type( t_fparser ), intent(inout) :: this
    real(p_double),      intent(in) :: float
    
    ! local variables
    
    !print *, 'add const = ', float
    
    this%code_size = this%code_size + 1
    this%data_size = this%data_size + 1
    this%data_stack(this%data_size) = float
    this%bytecode(this%code_size) = bcode_data

  end subroutine add_constant
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! test if current string is completely enclosed in parenthesis
!-----------------------------------------------------------------------------------------
  function if_enclosed_parent(in_str)
!-----------------------------------------------------------------------------------------
    implicit none
    
    ! dummy variables
    
    character (len = *),           intent(in) :: in_str
    logical :: if_enclosed_parent
        
    ! local variables     
    integer :: str_len 
    integer :: j,k
    
    ! executable statements
    
    str_len = len(in_str)
    if_enclosed_parent = .false.
    
    if ((in_str(1:1) == '(') .and. (in_str(str_len:str_len) == ')')) then
      k = 0
      do j = 2, str_len-1
        select case (in_str(j:j))
          case ('(') 
            k = k+1
          case (')')
            k = k-1
          case default
            continue
        end select
        if (k < 0) return
      enddo
      if (k == 0) if_enclosed_parent = .true.  
    endif
    
  end function if_enclosed_parent
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! returns the position of the first ',' at the first parenthesis level
!-----------------------------------------------------------------------------------------
  function comma_pos(in_str)
!-----------------------------------------------------------------------------------------
    implicit none
    
    ! dummy variables
    
    character (len = *),           intent(in) :: in_str
    integer :: comma_pos
        
    ! local variables     
    integer :: str_len, k
    
    ! executable statements
    
    str_len = len(in_str)

    k = 0
    do comma_pos = 1, str_len
        select case (in_str(comma_pos:comma_pos))
          case ('(') 
            k = k+1
          case (')')
            k = k-1
          case (',')
            if (k == 0) return
          case default
            continue
        end select
        if (k < 0) return
    enddo

    comma_pos = -1
    
  end function comma_pos
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!   check syntax on the supplied function
!-----------------------------------------------------------------------------------------
  subroutine check_syntax_fparser( this, ierr )
!-----------------------------------------------------------------------------------------

    implicit none
    
    ! dummy variables
    
    type( t_fparser ),                   intent(inout) :: this
    integer,                             intent (out)  :: ierr 
        
    !local variables    
    
    integer, dimension(MAX_NESTED_FUNC)  :: parent_count  ! parethesis count at the current level
    integer, dimension(MAX_NESTED_FUNC)  :: param_count   ! parameters count at the current level
    integer, dimension(MAX_NESTED_FUNC)  :: n_params      ! number of parameters for the function
                                                          ! at the current level
    
    integer :: level             ! function level being processed
    
    integer :: i, &              ! current position in the string being parsed
               fun_len, &        ! length of the string being parsed
               r_len, &          ! length of the number processed
               op_len, &         ! length of the operator processed
               f_len, &          ! lengh of the function processed
               var_len, &        ! length of the variable processed
               tmp_nparam      
               
    character(len=1) :: c        ! character being processed
    
    character(len=len_trim(this%func)) :: func ! local copy of the string to process 
    
    integer :: numberType        !  type of number being processed

    ! executable statements
    
    ! using this local variables ensures the length is correct
    func = trim(this%func)
    fun_len = len(func)
    
    ! iniatilize 
    parent_count = 0
    i = 1 
    ierr = 0
    
    level = 1
        
parse:  do
      if (i > fun_len) then 
        ierr = err_eof; call syntax_error(func,i,ierr); return
      endif
      
      c = func(i:i)
      ! check for leading + or -     
      if (c == '-' .or. c == '+') then
        i = i+1
        if (i > fun_len) then
          ierr = err_eof; return
        endif
        ! an operator must not follow
        if (is_operator(func(i:)) > 0) then
          ierr = err_mult_op; call syntax_error(func,i,ierr); return
        endif
        c = func(i:i)
      endif
      
      ! check for functions
      if (is_function(func(i:), f_len , tmp_nparam ) > 0) then
        i = i + f_len
        if (i > fun_len) then
          ierr = err_eof; return
        endif
        c = func(i:i)
        if ( c /= '(' ) then
          ierr = err_fun_par; call syntax_error(func,i,ierr); return
        endif        
        
        i = i + 1
        level = level + 1
        parent_count(level) = 0
        param_count(level) = 0
        n_params(level) = tmp_nparam
        cycle parse
      endif
      
      ! check for opening parenthesis
      if ( c == '(' ) then
         ! if opening ( increase parent_count and go to the beginning of the loop
         parent_count(level) = parent_count(level) + 1
         i = i+1
         cycle parse
      endif

      
      ! check for numbers and variables
      if (scan(c, '0123456789.') > 0) then 
        ! constant
        
        call parsenumber( func(i:), numberType, r_len)
        if (numberType < 0) then
          ierr = err_inv_num; call syntax_error(func,i,ierr); return 
        endif
        i = i + r_len
        if (i > fun_len) exit
        c = this%func(i:i)
      elseif (is_variable_fparser(this, func(i:), var_len) > 0) then ! variable
        ! variable
        
        i = i + var_len
        if (i > fun_len) exit
        c = this%func(i:i)
        
      endif
      
      ! check for closing parenthesis
      do while (c == ')')
        if ( i <= 1 ) then
          ierr = err_parent; return
        endif

        ! decrement parethesis count
        parent_count(level) = parent_count(level) - 1

        
        ! check if closing too many )
        if (parent_count(level) < 0) then
          ! check if closing function parenthesis
           if (level > 1) then
             !check if number of parameters is correct
             if (param_count(level) /= n_params(level)-1) then
               ierr = err_param; call syntax_error(func,i,ierr); return
             endif
             
             level = level - 1
             i = i+ 1
             if (i > fun_len) exit
             cycle parse
           endif
          ierr = err_parent;call syntax_error(func,i,ierr); return
        endif
        
        if (this%func(i-1:i-1) == '(') then
          ierr = err_parent; call syntax_error(func,i,ierr);return
        endif
        
        i = i+1
        if (i > fun_len) exit
        c = this%func(i:i)
      enddo
      
      ! check for new parameter to function
      if (c == ',') then
        if ( i <= 1 ) then
          ierr = err_syntax; call syntax_error(func,i,ierr); return
        endif
        
        if ( (func(i-1:i-1) == '(') .or. (func(i-1:i-1) == ',') ) then
          ierr = err_param; call syntax_error(func,i,ierr); return
        endif
        
        if (parent_count(level) /= 0) then
          ierr = err_parent; call syntax_error(func,i,ierr); return
        endif
        
        param_count(level) = param_count(level) + 1
        if (param_count(level) >= n_params(level)) then
          ierr = err_param; call syntax_error(func,i,ierr); return
        endif  
        i = i + 1 
        cycle parse
      endif
      
      ! if we reach this point a valid operator must follow
      if (i > fun_len) exit
      if (is_operator(func(i:), op_len) > 0) then
        i = i + op_len
        if (i > fun_len) then
          ierr = err_eof; call syntax_error(func,i,ierr); return
        endif
        ! an operator must not follow
        if ( is_operator(func(i:)) > 0 ) then
          ierr = err_mult_op; call syntax_error(func,i,ierr); return
        endif
      else
        ierr = err_syntax; call syntax_error(func,i,ierr); return
      endif
      
    enddo parse
     
    if ((parent_count(1) > 0) .or. (level > 1)) then
      ierr = err_parent
      call syntax_error(func,i,ierr);
    endif
    
  end subroutine check_syntax_fparser
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!   issue the appropriate error message
!-----------------------------------------------------------------------------------------
  subroutine syntax_error(in_str, pos, ierr, nparam )
!-----------------------------------------------------------------------------------------

     character (len=*) , intent (in) :: in_str
     integer , intent(in) :: pos, ierr
     integer , intent(in), optional :: nparam
     
!     character (len=pos-1) :: err_point 
     character (len = 64)  :: err_msg
     integer :: local_nparam
     
     if (present(nparam)) then
       local_nparam = nparam
     else
       local_nparam = 1
     endif
     
!     err_point = " "
     
     write(0,'(A)') "Error compiling function:"
     write(0,'(A)') in_str
     write(0,'(A,A)') repeat('-',pos-1), '^'
     
     select case (ierr)
       case ( err_eof )
         err_msg =  "End of line reached"
       case ( err_mult_op )
         err_msg =  "Consecutive operators" 
       case ( err_inv_num )
         err_msg =  "Invalid number format"
       case ( err_parent )
         err_msg =  "Parenthesis mismatch"
       case ( err_syntax )
         err_msg =  "Syntax error"
       case ( err_fun_par )
         err_msg = "An opening parenthesis '(' must follow a function"
       case ( err_param )
         write(err_msg,'(A,I0,A)') "Invalid number of parameters, ",local_nparam," expected"
       case default
         err_msg =  "Syntax error"
     end select
     
     write(0,'(A,A,I0)') trim(err_msg), ", col = ", pos


  end subroutine syntax_error
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! checks if in_str begins with a variable
!-----------------------------------------------------------------------------------------
  function is_variable_fparser( this, in_str , var_len)
!-----------------------------------------------------------------------------------------

    implicit none
    
    ! dummy variables

    type( t_fparser ),                   intent(in) :: this
    character ( len = * ),               intent(in) :: in_str
    integer                                         :: is_variable_fparser 
    integer, optional,                   intent(out):: var_len
    
    ! local variables
    integer :: i, l
        
    ! executable statements
    
    i = 1
            
    do while (i <=  this%nvars)
      l = min(len_trim(this%var_list(i)), len(in_str))
      if (in_str(1:l) == trim(this%var_list(i))) then
        is_variable_fparser = i
        if (present(var_len)) var_len = len_trim(this%var_list(i))
        return
      endif
      i = i+1
    enddo
    
    is_variable_fparser = -1

  end function is_variable_fparser
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! checks if in_str begins with an operator and return the opcode
!-----------------------------------------------------------------------------------------
  function is_operator( in_str , op_len)
!-----------------------------------------------------------------------------------------

    implicit none
    
    ! dummy variables
    
    character ( len = * ),               intent(in) :: in_str
    integer                                         :: is_operator 
    integer, optional,                   intent(out):: op_len
    
    ! local variables
    integer :: i, l
        
    ! executable statements
    
    i = 1
    is_operator = 0
    
    do while (i <=  size( oper ))
      if ( oper(i) /= '  ' ) then
		is_operator = is_operator + 1
		l = min(len_trim(oper(i)), len(in_str))
		if (in_str(1:l) == oper(i)) then
		  if (present(op_len)) op_len = len_trim(oper(i))
		  return
		endif
      endif
	  i = i+1
    enddo
    
    is_operator = -1

  end function is_operator
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! checks if in_str begins with a function
!-----------------------------------------------------------------------------------------
  function is_function( in_str , f_len, n_params)
!-----------------------------------------------------------------------------------------

    implicit none
    
    ! dummy variables
    
    character(len=*), intent(in) :: in_str
    integer :: is_function
    integer, optional,intent(out):: f_len
    integer, optional,intent(out):: n_params
    
    ! local variables
    integer :: i, l
        
    ! executable statements
    
    i = 1
    
    do while (i <=  n_functions)
      l = min(len_trim(functions(i)%name), len(in_str))
      if (in_str(1:l) == functions(i)%name) then
        is_function = i
        if (present(f_len)) f_len = len_trim(functions(i)%name)
        if (present(n_params)) n_params = functions(i)%n_params
        return
      endif
      i = i+1
    enddo
    
    is_function = -1

  end function is_function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! return in_str with white spaces (spaces and tabs) removed
!-----------------------------------------------------------------------------------------
  function remove_spaces( in_str )
!-----------------------------------------------------------------------------------------
      
    character ( len = * ), intent(in) :: in_str
    character (len = len_trim(in_str))     :: remove_spaces 

    integer :: i, j
    remove_spaces = ""
    j = 1
    do i = 1, len(in_str)
      ! if (in_str(i:i) .ne. " ") then 
      if (.not. is_blank(in_str(i:i))) then 
        remove_spaces(j:j) = in_str(i:i)
        j = j +1
      endif
    enddo

  end function remove_spaces
!-----------------------------------------------------------------------------------------


end module m_fparser
