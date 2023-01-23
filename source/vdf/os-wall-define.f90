#include "os-config.h"
#include "os-preprocess.fpp"

module m_wall_define

#include "memory/memory.h"

  use m_parameters

  private

  type :: t_wall
        
  ! perp. direction of the wall surface
  integer :: idir

  ! boundary of the wall (lower=1,upper=2)
  integer :: ibnd

  ! range of wall at the boundary i.e. bounds perpendicular
  ! to the wall
  integer, dimension(2) :: range

  ! pointers to 1D, 2D and 3D field data 
  real(p_k_fld), dimension(:,:), pointer :: f1 => null()
  real(p_k_fld), dimension(:,:,:), pointer :: f2 => null()
  real(p_k_fld), dimension(:,:,:,:), pointer :: f3 => null()
  
    ! spatial dimensions
  integer :: x_dim = -1

  ! field dimensions
  integer :: f_dim = -1
  
  ! grid size
  integer, dimension(p_max_dim) :: nx = 0

    ! guard cell information - full
  integer, dimension(2, p_max_dim) :: gc_num = 0
  
        ! guard cell information - comms only, may be smaller than gc_num
    integer, dimension(2, p_max_dim) :: gc_num_comm = 0

  end type t_wall

  interface idir
  module procedure idir_wall
  end interface

  interface ibnd
  module procedure ibnd_wall 
  end interface

  interface range
  module procedure  range_wall
  end interface

  interface get_range
  module procedure  get_range_wall
  end interface

  interface nx
  module procedure  nx_wall
  end interface

  interface field_comp
  module procedure field_comp_wall
  end interface

  interface space_dim
  module procedure space_dim_wall
  end interface

  interface alloc
    module procedure alloc_wall
    module procedure alloc_1d_wall
    module procedure alloc_bound_1d_wall
  end interface

  interface freemem
    module procedure freemem_wall
    module procedure freemem_1d_wall
  end interface

interface lbound
  module procedure lbound_wall
  module procedure lbound_dim_wall
end interface

interface ubound
  module procedure ubound_wall
  module procedure ubound_dim_wall
end interface
 
interface check_nan
    module procedure check_nan_wall
end interface
 
 public :: t_wall, idir, ibnd, range, get_range, nx, field_comp, space_dim
 public :: alloc, freemem, lbound, ubound, check_nan

contains

! -----------------------------------------------------------------------------
!> Checks if any of the values are NaN or Inf
!! @param   this        The wall object to check
!! @param   msg         Optional additional message to display if any invalid 
!!                      values are found
subroutine check_nan_wall( this, msg )

    implicit none

    type(t_wall), intent(in) :: this
    character(len=*), intent(in), optional :: msg

    integer :: i1, i2, i3, f

    select case ( this%x_dim )
    case(1)
        do i1 = lbound( this%f1, 2 ), ubound( this%f1, 2 )
            do f = lbound( this%f1, 1 ), ubound( this%f1, 1 )
                if ( isnan( this%f1(f,i1) ) .or. isinf( this%f1(f,i1) ) ) then

                    if ( present(msg) ) write(0,'(A,I0,A,A)') '[', mpi_node(), '] ', trim(msg)

                    write(0,'(A,I0,A)') '[', mpi_node(), '] (* error *) check_nan_wall failed '
                    write(0,*) "[", mpi_node(), "] Bad cell [", f, i1, "]"
                    write(0,*) "[", mpi_node(), "] Value = ", this%f1(:,i1)
                    write(0,*) "[", mpi_node(), "] lbound = ", lbound( this%f1 )
                    write(0,*) "[", mpi_node(), "] ubound = ", ubound( this%f1 )
                    call abort_program( p_err_invalid )
                endif
            enddo
        enddo

    case(2)
        do i2 = lbound( this%f2, 3 ), ubound( this%f2, 3 )
            do i1 = lbound( this%f2, 2 ), ubound( this%f2, 2 )
                do f = lbound( this%f2, 1 ), ubound( this%f2, 1 )
                    if ( isnan( this%f2(f,i1,i2) ) .or. isinf( this%f2(f,i1,i2) ) ) then

                        if ( present(msg) ) write(0,'(A,I0,A,A)') '[', mpi_node(), '] ', trim(msg)

                        write(0,'(A,I0,A)') '[', mpi_node(), '] (* error *) check_nan_wall failed '
                        write(0,*) "[", mpi_node(), "] Bad cell [", f, i1, i2, "]"
                        write(0,*) "[", mpi_node(), "] Value = ", this%f2(:,i1,i2)
                        write(0,*) "[", mpi_node(), "] lbound = ", lbound( this%f2 )
                        write(0,*) "[", mpi_node(), "] ubound = ", ubound( this%f2 )
                        call abort_program( p_err_invalid )
                    endif
                enddo
            enddo
        enddo

    case(3)
        do i3 = lbound( this%f3, 4 ), ubound( this%f3, 4 )
            do i2 = lbound( this%f3, 3 ), ubound( this%f3, 3 )
                do i1 = lbound( this%f3, 2 ), ubound( this%f3, 2 )
                    do f = lbound( this%f3, 1 ), ubound( this%f3, 1 )
                        if ( isnan( this%f3(f,i1,i2,i3) ) .or. isinf( this%f3(f,i1,i2,i3) ) ) then

                            if ( present(msg) ) write(0,'(A,I0,A,A)') '[', mpi_node(), '] ', trim(msg)

                            write(0,'(A,I0,A,A)') '[', mpi_node(), '] (* error *) check_nan_wall failed '
                            write(0,*) "[", mpi_node(), "] Bad cell [", f, i1, i2, i3, "]"
                            write(0,*) "[", mpi_node(), "] Value = ", this%f3(:,i1,i2,i3)
                            write(0,*) "[", mpi_node(), "] lbound = ", lbound( this%f3 )
                            write(0,*) "[", mpi_node(), "] ubound = ", ubound( this%f3 )
                            call abort_program( p_err_invalid )
                        endif
                    enddo
                enddo
            enddo
        enddo

    end select
  
end subroutine

!---------------------------------------------------
pure function idir_wall( this )
!---------------------------------------------------
!       gives a the direction of the wall
!---------------------------------------------------

  implicit none

!       dummy variables

  integer :: idir_wall

  type(t_wall), intent(in) :: this

!       local variables 

!       executable statements

  idir_wall = this%idir

end function idir_wall
!---------------------------------------------------


!---------------------------------------------------
pure function ibnd_wall( this )
!---------------------------------------------------
!       gives a the boundary of the wall (upper or lower)
!---------------------------------------------------

  implicit none

!       dummy variables

  integer :: ibnd_wall

  type(t_wall), intent(in) :: this

!       local variables 

!       executable statements

  ibnd_wall = this%ibnd

end function ibnd_wall
!---------------------------------------------------


!---------------------------------------------------
function range_wall( this )
!---------------------------------------------------
!       gives a the range of the wall
!---------------------------------------------------

  implicit none

!       dummy variables

  integer, dimension(2) :: range_wall

  type(t_wall), intent(in) :: this

!       local variables 

!       executable statements

  range_wall = this%range

end function range_wall
!---------------------------------------------------

!---------------------------------------------------
subroutine get_range_wall( this, range )
!---------------------------------------------------
!       gives a the range of the wall
!---------------------------------------------------

  implicit none

  type(t_wall), intent(in) :: this
  integer, dimension(:), intent(out) :: range

  range(1) = this%range(1)
  range(2) = this%range(2)

end subroutine get_range_wall
!---------------------------------------------------


!---------------------------------------------------
function nx_wall( this )
!---------------------------------------------------
!       gives number of guard cells at the boundaries
!---------------------------------------------------

  implicit none

!       dummy variables

  type(t_wall), intent(in) :: this
 integer, dimension(this%x_dim) :: nx_wall


!       local variables 

!       executable statements

  nx_wall = this%nx(1:this%x_dim)

end function nx_wall
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
!       gives the number of components the field has
!---------------------------------------------------------------------------------------------------
function field_comp_wall( this )

  implicit none

  integer :: field_comp_wall

  type(t_wall), intent(in) :: this
  
  field_comp_wall = this%f_dim

end function field_comp_wall
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! gives the dimensionality of the space the wall is defined in
!---------------------------------------------------------------------------------------------------
function space_dim_wall( this )

  implicit none

  integer :: space_dim_wall

  type(t_wall),intent(in) :: this

  space_dim_wall = this%x_dim

end function space_dim_wall
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------
pure function lbound_wall( this )
!---------------------------------------------------
!---------------------------------------------------
  
  implicit none
  
  type( t_wall ), intent(in) :: this
  integer, dimension( this%x_dim + 1) :: lbound_wall
  
  select case (this%x_dim)
    case(1)
      lbound_wall = lbound( this%f1 )
    case(2)
      lbound_wall = lbound( this%f2 )
    case(3)
      lbound_wall = lbound( this%f3 )
    case default
      lbound_wall = 0
  end select
  
end function lbound_wall
!---------------------------------------------------


!---------------------------------------------------
pure function lbound_dim_wall( this, dim )
!---------------------------------------------------
!---------------------------------------------------
  
  implicit none
  
  integer :: lbound_dim_wall
  type( t_wall ), intent(in) :: this
  integer, intent(in) :: dim
  
  select case (this%x_dim)
    case(1)
      lbound_dim_wall = lbound( this%f1, dim )
    case(2)
      lbound_dim_wall = lbound( this%f2, dim )
    case(3)
      lbound_dim_wall = lbound( this%f3, dim )
    case default
      lbound_dim_wall = 0
  end select
  
end function lbound_dim_wall
!---------------------------------------------------


!---------------------------------------------------
pure function ubound_wall( this )
!---------------------------------------------------
!---------------------------------------------------
  
  implicit none
  
  type( t_wall ), intent(in) :: this
  integer, dimension( this%x_dim + 1) :: ubound_wall
  
  select case (this%x_dim)
    case(1)
      ubound_wall = ubound( this%f1 )
    case(2)
      ubound_wall = ubound( this%f2 )
    case(3)
      ubound_wall = ubound( this%f3 )
    case default
      ubound_wall = 0
  end select
  
end function ubound_wall
!---------------------------------------------------


!---------------------------------------------------
pure function ubound_dim_wall( this, dim )
!---------------------------------------------------
!---------------------------------------------------
  
  implicit none
  
  integer :: ubound_dim_wall
  type( t_wall ), intent(in) :: this
  integer, intent(in) :: dim
  
  select case (this%x_dim)
    case(1)
      ubound_dim_wall = ubound( this%f1, dim )
    case(2)
      ubound_dim_wall = ubound( this%f2, dim )
    case(3)
      ubound_dim_wall = ubound( this%f3, dim )
    case default
      ubound_dim_wall = 0
  end select
  
end function ubound_dim_wall
!---------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines

#define __TYPE__ type( t_wall )
#define __TYPE_STR__ "t_wall"
#define FNAME( a )  a ## _wall
#define __MAX_DIM__ 1
#include "memory/mem-template.h"

!---------------------------------------------------------------------------------------------------



end module m_wall_define

