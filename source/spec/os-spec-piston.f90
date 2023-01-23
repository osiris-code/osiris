!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     piston class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_piston

#include "memory/memory.h"

  use m_system

  use m_parameters

  use m_space
  use m_vdf
  
  use m_utilities
  use m_random

  use m_species_define
  use m_species_current 
  
  use m_fparser

  

  implicit none

!       restrict access to things explicitly declared public
  private

  ! piston parameters
  character(1), dimension(2), parameter :: p_sp_piston_transverse_symbols=(/"A","B"/)
  integer, parameter     :: p_sp_piston_uniform = 1
  integer, parameter     :: p_sp_piston_fparser = 2
  character(20), parameter           :: p_sp_piston_uniform_str = 'uniform'
  character(20), parameter           :: p_sp_piston_fparser_str = 'fparser'
  integer, parameter     :: p_sp_piston_up = 1
  integer, parameter     :: p_sp_piston_down = 2
        
  ! buffers for particle functions
  
  !buffers for piston
  ! indexes of prticles that crossed
  integer, pointer, dimension(:)   :: piston_sp_idx => null()
  ! number of particles in piston_sp_idx
  integer                              :: num_idx
  !for current deposit of piston
  real(p_k_part), pointer, dimension(:,:) :: pst_curr_xold => null(), &
                                            pst_curr_xnew => null()
  real(p_k_part), pointer, dimension(:,:) :: pst_curr_udeposit => null()
  
  interface read_nml_piston
    module procedure read_nml_pst
  end interface

  interface check_piston
    module procedure check_pst
  end interface

  interface update_piston
    module procedure update_pst
  end interface

  interface cross_particles_piston
    module procedure cross_particles_pst
  end interface

  interface init_buffers_piston
    module procedure init_buffers_pst
  end interface

  interface cleanup_buffers_piston
    module procedure cleanup_buffers_pst
  end interface

!       declare things that should be public
  public :: t_piston
  public :: read_nml_piston
  public :: check_piston
  public :: cross_particles_piston
  public :: update_piston
  public :: init_buffers_piston, cleanup_buffers_piston

interface alloc
  module procedure alloc_piston
  module procedure alloc_1d_piston
  module procedure alloc_bound_1d_piston
end interface

interface freemem
  module procedure freemem_piston
  module procedure freemem_1d_piston
end interface

  public :: alloc, freemem
 
 contains 


!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine read_nml_pst( this, input_file )

  use m_input_file

  implicit none

  type( t_piston ),  intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  integer         :: dim
  integer         :: updown
  real(p_k_part)         :: v
  real(p_k_part)         :: start_pos
  real(p_k_part)         :: start_time
  real(p_k_part)         :: stop_time
  real(p_k_part)         :: opacity_factor
  character(20)           :: profile_type
  character(p_max_expr_len) :: profile_expression


  namelist /nl_piston/ dim, updown, v, start_pos, start_time, stop_time, &
    opacity_factor, profile_type, profile_expression

  integer :: ierr

!       executable statements

!         set defaults for variables that can be omitted in
!         the input file
  dim = 1
  updown = p_sp_piston_up
  v = 0.0_p_k_part
  start_pos = 0.0_p_k_part
  start_time = 0.0_p_k_part
  stop_time = 0.0_p_k_part
  opacity_factor = 1.0_p_k_part
  profile_type = "uniform"
  profile_expression = ""

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_piston", ierr )
  if ( ierr /= 0 ) then
	print *, "   No piston parameters specified."
	print *, "   aborting..."
	stop
  endif

  read (input_file%nml_text, nml = nl_piston, iostat = ierr)
  if (ierr /= 0) then
	print *, "   Error reading piston parameters."
	print *, "   aborting..."
	stop
  endif
  
  ! setup of piston values

  if ( v <= 0.0_p_k_part) then
	print *, ""
	print *, "   Error reading piston parameters - v."
	print *, "   v should be greater than zero."
	print *, "   aborting..."
	stop          
  endif

  if ( updown == p_sp_piston_down ) v = -v

  this%dim = dim
  this%updown = updown
  this%u = v
  this%start_pos = start_pos
  this%start_time = start_time
  this%stop_time = stop_time
  this%opacity_factor = opacity_factor
  
  
  if ( start_time >= stop_time) then
	print *, ""
	print *, "   Error reading piston parameters - start_/stop_ time."
	print *, "   Stop_time must be greater than start_time. "
	print *, "   aborting..."
	stop          
  endif

  if ( dim > p_x_dim) then
	print *, ""
	print *, "   Error reading piston parameters - dim."
	print *, "   dim > p_x_dim. "
	print *, "   aborting..."
	stop          
  endif

  if ( updown /= p_sp_piston_up .AND. updown /= p_sp_piston_down) then
	print *, ""
	print *, "   Error reading piston parameters - updown."
	print *, "   updown must be either of 1, 2. "
	print *, "   aborting..."
	stop          
  endif
  select case ( profile_type )
	case ( p_sp_piston_uniform_str )
	  this%profile_type = p_sp_piston_uniform
	case ( p_sp_piston_fparser_str )
	  if( p_x_dim == 1) then
		print *, ""
		print *, "   Error reading piston parameters - profile_type"
		print *, "   It does not make sense to use f_parser for a 1D run."
		print *, "   aborting..."
		stop          
	  endif
	  this%profile_type = p_sp_piston_fparser

	  ! set up profile
	  call setup( this%profile_parser, profile_expression, &
				   p_sp_piston_transverse_symbols(1:p_x_dim-1), ierr )
	  if (ierr /= 0) then
		ERROR('Error seting up opacity profile for piston.')
		call abort_program(p_err_invalid)
	  endif
	case default
	  ERROR('Invalid profile_type_piston specified')
	  stop
  end select
  
  ! calculate constannt values

  this%squ = this%u**2 !u^2
  this%sqg = 1.0_p_k_part + this%u**2 !gamma^2
  this%g = sqrt( this%sqg ) ! gamma
  this%v = this%u / this%g  ! v

  !this%sqv = this%v**2

end subroutine read_nml_pst
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine init_buffers_pst( )
!---------------------------------------------------------------------------------------------------
!       initializes cache buffers for piston pusher
!       routines
!---------------------------------------------------------------------------------------------------
  implicit none
  
  ! setup buffers for advance_deposit routines
  
  ! setup buffers for piston
  call alloc( piston_sp_idx, (/ p_cache_size /))

  call alloc( pst_curr_xold, (/ p_x_dim, p_cache_size /))

  call alloc( pst_curr_xnew, (/ p_x_dim, p_cache_size /))

  call alloc( pst_curr_udeposit, (/ p_p_dim, p_cache_size /))

end subroutine init_buffers_pst
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine cleanup_buffers_pst( )
!---------------------------------------------------------------------------------------------------

  implicit none

  ! setup buffers for advance_deposit routines
  
  ! setup buffers for piston
  call freemem( piston_sp_idx)

  call freemem( pst_curr_xold)

  call freemem( pst_curr_xnew)

  call freemem( pst_curr_udeposit)

end subroutine cleanup_buffers_pst
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
        subroutine cross_particles_pst( piston_all, np, &
                       &  x, xbuf, ptrcur )
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
 
          implicit none

!       dummy variables

!          class(t_species),              intent(inout) :: this
          type(t_piston),dimension(:), pointer  :: piston_all
          real(p_k_part), dimension(:,:), intent(in)       :: x, xbuf
          integer, intent(in)                       :: ptrcur
          integer, intent(in)                       :: np

!       local variables

          integer                            :: i, pst_id, num_pistons, dim
          logical                                    :: piston_crossed


!       executable statements

          num_pistons=size(piston_all)
          
          ! generate list of particles to be processed by piston
          num_idx = 0
          do i = 1, np ! all particles
            do pst_id=1, num_pistons !all piston
              dim = piston_all(pst_id)%dim
              if ( piston_all(pst_id)%inbox ) then
                if( piston_all(pst_id)%updown == p_sp_piston_up ) then
                  piston_crossed = x(dim,ptrcur+i) > piston_all(pst_id)%pos_before .AND.   &
      &                            xbuf(dim,i)     <= piston_all(pst_id)%pos_after
                else
                  piston_crossed = x(dim,ptrcur+i) < piston_all(pst_id)%pos_before .AND.   &
      &                            xbuf(dim,i)     >= piston_all(pst_id)%pos_after
                endif !updown
                if( piston_crossed ) then
                  ! add particle to list of particles to be processed
                  num_idx = num_idx + 1
                  piston_sp_idx( num_idx ) = i
                  exit !of loop all piston
                endif !crossed
              endif !inbox
            enddo !all piston
          enddo ! all particles

        end subroutine cross_particles_pst
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
        subroutine update_pst( this, t, dt, l_space )
!---------------------------------------------------------------------------------------------------
!       updates pos_before, pos_after, inbox for piston
!---------------------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

          implicit none

!       dummy variables

          type(t_piston),                      intent(inout) :: this
          real(p_double),                     intent(in) ::  t
          real(p_k_part), intent(in) :: dt
          type( t_space ),                     intent(in) :: l_space

!       local variables

          real(p_k_part), dimension(p_x_dim)   :: l_xmin, l_xmax
          real(p_k_part)                       :: t_before, t_after

!       executable statements

          l_xmin   = real( xmin(l_space), p_k_part )
          l_xmax   = real( xmax(l_space), p_k_part )
          
          t_before = real( t - dt, p_k_part )
          t_after  = real( t, p_k_part )

          if ( this%start_time <= t_after .AND. this%stop_time >= t_after ) then ! piston active during this timestep ?

            ! update pos_before and pos_after
            this%pos_before = this%start_pos + this%v * ( t_before - this%start_time )
            this%pos_after  = this%start_pos + this%v * ( t_after  - this%start_time )

            ! update inbox, adjust boundaries by 3 dt so we include guard cells 
            this%inbox = this%pos_after >= l_xmin(this%dim) - 3 * dt .AND. &
     &                   this%pos_after <  l_xmax(this%dim) + 3 * dt

          else
            this%inbox=.false.
          endif ! piston active in this dt

       end subroutine update_pst
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
        subroutine check_pst( this, dt )
!---------------------------------------------------------------------------------------------------
!       updates pos_before, pos_after, inbox for piston
!---------------------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

          implicit none

!       dummy variables

          type(t_piston), dimension(:), intent(in) :: this
          real(p_double),              intent(in) :: dt

!       local variables

          integer                       :: pst_id_a, pst_id_b, dima, dimb, updowna, updownb
          real(p_k_part)                       :: edge_a, edge_b, va, vb, tsa, tsb, intervall, t_cross

!       executable statements

          do pst_id_a = 1, size(this)

            dima = this(pst_id_a)%dim
            va = this(pst_id_a)%v
            tsa = this(pst_id_a)%start_time
            updowna = this(pst_id_a)%updown
            edge_a = this(pst_id_a)%start_pos

            do pst_id_b = pst_id_a + 1, size(this)

              dimb = this(pst_id_b)%dim
              vb = this(pst_id_b)%v
              tsb = this(pst_id_b)%start_time
              updownb = this(pst_id_b)%updown
              edge_b = this(pst_id_b)%start_pos
              
              if( dima == dimb .AND. updowna /= updownb) then
                t_cross = (edge_b - edge_a + va*tsa - vb*tsb) / (va-vb)
                intervall = abs(real(2 * dt, p_k_part)/(va -vb))
                if( overlap( t_cross - intervall, t_cross + intervall, this(pst_id_a)%start_time, this(pst_id_a)%stop_time ) .AND. &
              &     overlap( t_cross - intervall, t_cross + intervall, this(pst_id_b)%start_time, this(pst_id_b)%stop_time ) ) then
                  print *, "Two or more piston of the same species and oposite direction"
                  print *, "get too close to each other (<c * dt)"
                  print *, "Abort..."
                  stop
                endif
              endif !oposit rirections

            enddo ! remaining piston
          
          enddo ! all piston
       end subroutine check_pst
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
        function overlap( al, ah, bl, bh )
!---------------------------------------------------------------------------------------------------
!       returns true if [al,ah] and [bl,bh] intersect
!---------------------------------------------------------------------------------------------------

!         <use-statements for this subroutine>

          implicit none

!       dummy variables

          real(p_k_part), intent(in) :: al,ah,bl,bh
          logical        :: overlap

          overlap = ( bl <= ah ) .AND. ( al <= bh )

       end function overlap
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines

#define __TYPE__ type( t_piston )
#define __TYPE_STR__ "t_piston"
#define FNAME( a )  a ## _piston
#define __MAX_DIM__ 1
#include "memory/mem-template.h"

!---------------------------------------------------------------------------------------------------

      end module m_piston
