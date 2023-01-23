#include "os-config.h"
#include "os-preprocess.fpp"

module m_psource_constq

#include "memory/memory.h"

use m_species_define
use m_parameters
use m_psource_std
use m_current_define
use m_node_conf

implicit none

private

type, extends( t_psource_std ) :: t_psource_constq

contains

  procedure :: inject => inject_constq
  procedure :: read_input => read_input_constq

end type t_psource_constq

public :: t_psource_constq

contains

!---------------------------------------------------------------------------------------------------
! This routine injects particles into an area defined by grid cell indexes ig_xbnd_inj, using
! fixed charge per particle.
!
! It allows for the loading of density psources that can be represented as separable functions, i.e.
! n(x,y) = nx(x) * ny(y) (2D) or n(x,y,z) = nx(x) * ny(y) * nz(z) (3D).
!
! The required density psource is obtained by placing more particles in regions of higher density.
!---------------------------------------------------------------------------------------------------
function inject_constq( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, send_msg, &
                        recv_msg ) result(num_inj)

  use m_species_udist

  implicit none

  class(t_psource_constq), intent(inout) :: this
  class(t_species), intent(inout), target :: species
  integer, dimension(:, :), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout)   :: jay
  class( t_node_conf ), intent(in)     :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: num_inj

  integer, parameter :: p_max_x_dim = 3

  ! cell size (p_x_dim)
  real(p_double), dimension(p_max_x_dim) :: dx
  real(p_double), dimension(p_max_x_dim) :: g_xmin

  ! half distance between particles
  real(p_k_part), dimension(p_max_x_dim) :: dxp_2

  integer :: i, i1, i2, i3, ipart

  ! volume that each particle occupies
  real(p_k_part) :: pvol

  ! particle positions (global / inside cell )
  real(p_k_part), dimension(:,:), pointer :: ppos, ppos_cell

  ! particle charges
  real(p_k_part), dimension(:), pointer :: pcharge

  ! for constant q
  integer, dimension(p_x_dim) :: gppe, pnst
  integer :: gnp, noff
  integer, dimension(p_x_dim) :: sample_rate
  integer :: srn
  real(p_k_part), dimension(p_x_dim) :: x, qlow, qx
  real(p_k_part) :: const_q
  real(p_k_part), dimension(1) :: norm
  real(p_double) :: den_ltot
  integer, dimension(p_max_dim) :: part_nx
  real(p_k_part), parameter :: epslon = 10_p_k_part/huge(1.0_p_k_part)

  ! Set total number of injected particles to 0
  num_inj = 0

  ! This is already checked in the input file
  if ( species%coordinates == p_cylindrical_b ) then
    ERROR('Cylindrical coordinates are not yet supported with fixed charge injection, aborting.')
    call abort_program(p_err_notimplemented)
  endif

  do i = 1, p_x_dim
    ! get cell size
    dx(i) = species%dx(i)

    ! get global mininum (this is shifted by +0.5 cells from global simulation values)
    g_xmin( i ) = species%g_box( p_lower , i )

    ! get half distance between sampling grid points
    sample_rate(i) = this % sample_rate(i)

    dxp_2(i) = 0.5_p_k_part/sample_rate(i)
  enddo

  ! find total number sampling grid points per cell
  gppe(1) = species%g_nx(1) * sample_rate(1) + 1
  gnp = gppe(1)
  ! a number of huge(p_int64) particles can not be crammed into one node
  if(huge(1) .gt. species%tot_par_x(1)) then
    srn = species%tot_par_x(1)
  else
    srn = huge(1)
  endif
  do i = 2, p_x_dim
    gppe(i) = species%g_nx(i) * sample_rate(i) + 1
    gnp = gnp + gppe(i)
    ! decide the buffer size,
    if (srn .lt. species%tot_par_x(i)) then
      if(huge(1) .gt. species%tot_par_x(i)) then
        srn = species%tot_par_x(i)
      else
        srn = huge(1)
      endif
    endif
  enddo

  ! initialize temp buffers
  call alloc( ppos, (/p_x_dim, gnp/) )
  call alloc( ppos_cell, (/p_x_dim, srn /) )
  call alloc( pcharge, (/ gnp /) )

  den_ltot = 0.0
  num_inj = 0
  ! Get position of sampling points for the full box
  ! Inject particles
  noff = 0
  pcharge(1) = 0
  i = 0
  ! find non-zero normalization factor
  if (p_x_dim .ne. 1) then
    do
      do i1 = 1, p_x_dim
        x(i1) = real( g_xmin(i1) + dx(i1) * dxp_2(i1) * &
          &  ( 2 * mod(gppe(i1)/2+i, gppe(i1)) + 1))
        ppos(i1,1) = x(i1)
      enddo
      call this % get_den_value( ppos(:,1:1), 1, norm )
      if ( (norm(1) > 0) .or. (i>product(gppe(1:p_x_dim))) ) exit
      i = i + 1
    enddo
  else
    norm(1) = 1.0
  endif

  ! flatten the axes data
  noff = 0
  i = 0
  do i1 = 1, p_x_dim
    do i = 1, gppe(i1)
      ppos(:,noff+i) = x(:)
      ! overwrite in one direction
      ppos(i1,noff+i) = real( g_xmin(i1) + (2 * i - 2) * dxp_2(i1) * dx(i1), p_k_part)
    enddo
    noff = noff + gppe(i1)
  enddo
  call this % get_den_value( ppos, gnp, pcharge )
  do i = 1, gnp
    if (pcharge(i) <= this%den_min) then
      pcharge(i) = 0
    endif
  enddo
  ! numerically integrate the density function
  i = 1
  noff = 0
  do i1 = 1, p_x_dim
    qlow(i1) = 0
    i2 = sample_rate(i1) * (species%my_nx_p(p_lower, i1) - 1)
    do while ( i <= noff + i2 )
      qlow(i1) = qlow(i1) + 0.5 * (pcharge(i)+pcharge(i+1))
      i = i + 1
    enddo
    ! accumulated density at lower boundaries of this node
    pnst(i1) = i
    den_ltot = qlow(i1)
    do while ( i < noff + gppe(i1) )
      den_ltot = den_ltot + 0.5 * (pcharge(i)+pcharge(i+1))
      i = i + 1
    enddo
    i = i + 1    ! skip the last points (bonudaries of axes)
    qx(i1) = den_ltot / species%tot_par_x(i1)
    noff = noff + gppe(i1)
    qlow(i1) = qlow(i1) - (floor(qlow(i1) / qx(i1) + 0.5) - 0.5) * qx(i1)
  enddo
  const_q = sign( product(qx(1:p_x_dim))/ &
    & ((norm(1)**(p_x_dim-1))*product(sample_rate(1:p_x_dim))), species%rqm )
  ! find particle positions in this node
  noff = 0
  do i1 = 1, p_x_dim
    part_nx(i1) = 0
    noff = pnst(i1) - 1
    i = 1
    i2 = (ig_xbnd_inj(p_upper,i1)-ig_xbnd_inj(p_lower,i1) + 1)*sample_rate(i1)
    do while (.true.)
      if (qlow(i1) < qx(i1)) then
        if ( i <= i2 ) then
          qlow(i1) = qlow(i1) + 0.5 * (pcharge(i+noff) + pcharge(i+noff+1))
        else
          exit
        endif
        i = i + 1
      else
        part_nx(i1) = part_nx(i1) + 1
        qlow(i1) = qlow(i1) - qx(i1)
        ! here pvol is used as tmp storing the density gradient.
        pvol = pcharge(i+noff)-pcharge(i+noff-1)
        if (pvol <= epslon) then  ! piece-wise constant
          ppos_cell(i1,part_nx(i1)) = dxp_2(i1) * 2_p_k_part * &
          & (i - qlow(i1)/pcharge(i+noff) - 1)
        else  ! assuming linear gradient within this piece
          ppos_cell(i1,part_nx(i1)) = dxp_2(i1) * 2_p_k_part * &
          & (i - (sqrt(pcharge(i+noff-1)**2 + &
            & 2*qlow(i1)*pvol) - pcharge(i+noff-1))/pvol - 1)
        endif
      endif
    enddo
  enddo

  num_inj = product(part_nx(1:p_x_dim))
  ipart = species%num_par + 1

  if ( num_inj > species%num_par_max - ipart ) then
    call species % grow_buffer( species%num_par_max + num_inj + p_spec_buf_block )
  endif


  if ( num_inj > 0 ) then

    ! loop through all the injection cells and
    ! inject particles, normalizing charge
    select case (p_x_dim)
    case(1)
      do i1 = 1, part_nx(1)
        species%x(1, ipart ) = ppos_cell(1,i1) - FLOOR(ppos_cell(1,i1)) - 0.5_p_k_part
        species%q(ipart)     = const_q
        species%ix(1, ipart) = FLOOR(ppos_cell(1,i1)) + 1
        ipart = ipart + 1
      enddo
    case(2)
      do i1 = 1, part_nx(1)
        do i2 = 1, part_nx(2)
          species%x(1, ipart ) = ppos_cell(1,i1) - FLOOR(ppos_cell(1,i1)) - 0.5_p_k_part
          species%x(2, ipart)  = ppos_cell(2,i2) - FLOOR(ppos_cell(2,i2)) - 0.5_p_k_part
          species%q(ipart)     = const_q
          species%ix(1, ipart) = FLOOR(ppos_cell(1,i1)) + 1
          species%ix(2, ipart) = FLOOR(ppos_cell(2,i2)) + 1
          ipart = ipart + 1
        enddo
      enddo
    case(3)
      do i1 = 1, part_nx(1)
        do i2 = 1, part_nx(2)
          do i3 = 1, part_nx(3)
            species%x(1, ipart ) = ppos_cell(1,i1) - FLOOR(ppos_cell(1,i1)) - 0.5_p_k_part
            species%x(2, ipart)  = ppos_cell(2,i2) - FLOOR(ppos_cell(2,i2)) - 0.5_p_k_part
            species%x(3, ipart)  = ppos_cell(3,i3) - FLOOR(ppos_cell(3,i3)) - 0.5_p_k_part
            species%q(ipart)     = const_q
            species%ix(1, ipart) = FLOOR(ppos_cell(1,i1)) + 1
            species%ix(2, ipart) = FLOOR(ppos_cell(2,i2)) + 1
            species%ix(3, ipart) = FLOOR(ppos_cell(3,i3)) + 1
            ipart = ipart + 1
          enddo
        enddo
      enddo
    end select
  endif

  ! free temporary memory
  call freemem( ppos )
  call freemem( ppos_cell )
  call freemem( pcharge )

  ! If any particles were injected set tags and add momentum
  if ( num_inj > 0 ) then

    ! Initialize particle momentum
    call set_momentum( species, species%num_par+1, species%num_par + num_inj )

  endif


end function inject_constq
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Read information from input file
!---------------------------------------------------------------------------------------------------
subroutine read_input_constq( this, input_file, coordinates )

  use m_input_file

  implicit none

  class( t_psource_constq ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: coordinates

  integer :: i

  if ( coordinates == p_cylindrical_b ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "   Error reading spec parameters"
      write(0,*) "   Using constant charge particles is not (yet) supported for cylindrical geometry"
      write(0,*) "   aborting..."
    endif
    stop
  endif

  ! No changes to the input file from standard injection
  call this % t_psource_std % read_input( input_file, coordinates )

  ! Additional parameter validation
  do i = 1, p_x_dim
    if ( this%sample_rate(i) < 1 ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "   Error reading spec parameters"
        write(0,*) "   Sampling rate for density function must be >= 1 in all directions:"
        write(0,"(A,I0,A,I0)") "    -> sample_rate(",i,") = ", this%sample_rate(i)
        write(0,*) "   aborting..."
      endif
      stop
    endif
  enddo

  ! Constant charge injection works only for separable density functions
  if ( .not. this % if_mult ) then

    ! Test for cases known not to work
    select case ( this%type(1) )
    case (p_channel)
      if ( p_x_dim > 2 ) then
        if ( mpi_node() == 0 ) then
          write(0,*) "   Error reading spec parameters"
          write(0,*) "   'channel' density psource cannot be used with constant charge injection in 3D"
          write(0,*) "   aborting..."
        endif
        stop
      endif

    case (p_sphere)
      if ( p_x_dim > 1 ) then
        if ( mpi_node() == 0 ) then
          write(0,*) "   Error reading spec parameters"
          write(0,*) "   'sphere' density psource cannot be used with constant charge injection in 2D and 3D"
          write(0,*) "   aborting..."
        endif
        stop
      endif

    case default
      if ( p_x_dim > 1 ) then
        if ( mpi_node() == 0 ) then
          write(0,*) "(*warning*) User chose constant charge injection but density psource may not"
          write(0,*) "(*warning*) be a separable function, please check injected density."
        endif
      endif

    end select

  endif

end subroutine read_input_constq
!---------------------------------------------------------------------------------------------------

end module m_psource_constq
