#include "os-config.h"
#include "os-preprocess.fpp"

module m_vdf_average

use m_system
use m_parameters

use m_grid_define
use m_grid

use m_vdf_define
use m_space
use m_node_conf
use m_diagnostic_utilities

implicit none

private

interface report_ave
  module procedure report_ave_vdf
end interface

public :: report_ave


contains

!--------------------------------------------------------------------------------------------------
! This routine determines the range of cells in the averaged grid that use information from
! the global grid data that belongs to the specified range. It will also calculate the size of the
! required buffer. Returns -1 in local_ave_range, lsize if no data is required
!--------------------------------------------------------------------------------------------------
subroutine aveDataRange( global_nx, local_nx, rank, ave_nx, ave_nx_range, &
                         local_ave_nx, lsize )
!--------------------------------------------------------------------------------------------------

  implicit none

  integer, dimension(:), intent(in)   :: global_nx    ! global grid size
  integer, dimension(:,:), intent(in) :: local_nx     ! position of local data on global grid
  integer, intent(in) :: rank                         ! rank of the grids

  integer, dimension(:), intent(in) :: ave_nx         ! averaged grid size
  integer, dimension(:,:), intent(in) :: ave_nx_range ! range of averaged grid being processed

  integer, dimension(:,:), intent(out) :: local_ave_nx
  integer, intent(out) :: lsize


  real(p_double) :: dx_ave, rdx_ave
  integer :: i, lbound, ubound

  lsize = 1

  ! find normalized cell sizes
  do i = 1, rank
    dx_ave = real( global_nx(i), p_double ) / ave_nx(i)
    rdx_ave = 1 / dx_ave

    lbound = int( (local_nx( p_lower, i )-1) * rdx_ave ) + 1

    ! proceed diferently if the division is even or not
    if ( modulo( global_nx(i), ave_nx(i) ) == 0 ) then

      ubound = int( (local_nx( p_upper, i )-1) * rdx_ave ) + 1

    else
      ! the last cell needs special treatment
      if ( local_nx( p_upper, i ) < global_nx(i) ) then
        ubound = int( local_nx( p_upper, i ) * rdx_ave ) + 1
      else
        ubound = ave_nx(i)
      endif
    endif

    ! check if inside range
    if (((lbound<=ave_nx_range(p_upper,i)) .and. (ubound>=ave_nx_range(p_lower,i)))) then
       local_ave_nx( p_lower,i) = max( lbound, ave_nx_range(p_lower,i) )
       local_ave_nx( p_upper,i) = min( ubound, ave_nx_range(p_upper,i) )
       lsize = lsize * ( local_ave_nx( p_upper,i ) - local_ave_nx( p_lower,i ) + 1 )
    else
       local_ave_nx = -1
       lsize = -1
       return
    endif

  enddo

end subroutine aveDataRange
!--------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------
subroutine average1D( fld, fcomp, buffer, global_nx, nx_p, ave_nx, local_ave_nx, envelope )
!--------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 1

  real( p_k_fld ), dimension(:,:), intent(in) :: fld        ! field data
  integer, intent(in) :: fcomp                                ! field component to average
  real( p_diag_prec ), dimension(:), intent(inout) :: buffer    ! buffer to save data to
  integer, dimension(:), intent(in) :: global_nx            ! global field size
  integer, dimension(:,:), intent(in) :: nx_p                ! local field data position on global grid
  integer, dimension(:), intent(in) :: ave_nx                ! averaged grid size
  integer, dimension(:,:), intent(in) :: local_ave_nx        ! local averaged data position on averaged grid
  logical, intent(in) :: envelope

  integer :: i, lbound, ubound
  real(p_double), dimension(rank) :: dx_ave, rdx_ave
  real(p_double) :: delta
  integer, dimension(rank) :: mod_ave
  integer, dimension(2, rank) :: grange
  integer :: gi1, ai1, ai1b, idx
  real( p_diag_prec ) :: tmp


  ! find local grid limits for depositions and cell normalizations
  do i = 1, rank
    dx_ave(i) = real( global_nx(i), p_double ) / ave_nx(i)
    rdx_ave(i) = 1 / dx_ave(i)

    mod_ave(i) = modulo( global_nx(i), ave_nx(i) )

    if ( mod_ave(i)  == 0 ) then
      lbound = int( (local_ave_nx( p_lower, i ) - 1)*dx_ave(i) ) + 1
      ubound = int( local_ave_nx( p_upper, i ) *dx_ave(i) )
    else
      lbound = int( (local_ave_nx( p_lower, i ) - 1)*dx_ave(i) ) + 1
      if ( local_ave_nx( p_upper,i ) < ave_nx(i) ) then
        ubound = int( local_ave_nx( p_upper, i )*dx_ave(i) ) + 1
      else
        ubound = global_nx(i)
      endif
    endif
    grange( p_lower, i ) = max( lbound, nx_p( p_lower, i ) )
    grange( p_upper, i ) = min( ubound, nx_p( p_upper, i ) )
  enddo

  if ( envelope ) then

    ! envelope data

    if (  mod_ave(1)  == 0 ) then
       ! each column only deposits to 1 row
       do gi1 = grange( p_lower, 1 ), grange( p_upper, 1 )
          ai1 = int( (gi1-1) * rdx_ave(1) ) + 1
          idx = (ai1 - local_ave_nx(p_lower,1)) + 1

          tmp = abs( real ( fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1), p_diag_prec ) )
          if ( tmp > buffer(idx) ) then
            buffer(idx) = tmp
          endif
       enddo
    else
       do gi1 = grange( p_lower, 1 ), grange( p_upper, 1 )
          ai1 = int( (gi1-1) * rdx_ave(1) ) + 1

          if ( gi1 < global_nx(1) ) then
            ai1b = int( gi1 * rdx_ave(1) ) + 1
          else
            ai1b = ai1
          endif

          idx = (ai1 - local_ave_nx(p_lower,1)) + 1
          tmp = abs( real ( fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1), p_diag_prec ) )

          if ( ai1b == ai1 ) then
             if ( tmp > buffer(idx) ) then
               buffer(idx) = tmp
             endif
          else
             ! check if they belong to the current chunk
             ! (this is only needed for the last dimension)
             if ( ai1 >= local_ave_nx(p_lower,1) ) then
               if ( tmp > buffer(idx) ) then
                 buffer(idx) = tmp
               endif
             endif

             if ( ai1b <= local_ave_nx(p_upper,1) ) then
               if ( tmp > buffer(idx+1) ) then
                 buffer(idx+1) = tmp
               endif
             endif
          endif
       enddo
    endif

  else

    ! average data

    if (  mod_ave(1)  == 0 ) then
       ! each column only deposits to 1 row
       do gi1 = grange( p_lower, 1 ), grange( p_upper, 1 )
          ai1 = int( (gi1-1) * rdx_ave(1) ) + 1
          idx = (ai1 - local_ave_nx(p_lower,1)) + 1
          buffer(idx) = buffer(idx) + real ( rdx_ave(1) * &
                                        fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1), p_diag_prec )
       enddo
    else
       do gi1 = grange( p_lower, 1 ), grange( p_upper, 1 )
          ai1 = int( (gi1-1) * rdx_ave(1) ) + 1

          if ( gi1 < global_nx(1) ) then
            ai1b = int( gi1 * rdx_ave(1) ) + 1
          else
            ai1b = ai1
          endif

          idx = (ai1 - local_ave_nx(p_lower,1)) + 1
          if ( ai1b == ai1 ) then
             buffer(idx) = buffer(idx) + real ( rdx_ave(1) * &
                                           fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1), p_diag_prec )
          else
             ! the local global cell contributes to 2 averaged cells

             ! get fraction of global grid cell in averaged grid cell
             delta = (ai1b - 1) * dx_ave(1) - gi1 + 1

             ! check if they belong to the current chunk
             ! (this is only needed for the last dimension)
             if ( ai1 >= local_ave_nx(p_lower,1) ) then
               buffer(idx) = buffer(idx) + real ( rdx_ave(1) * delta * &
                                             fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1), p_diag_prec )
             endif

             if ( ai1b <= local_ave_nx(p_upper,1) ) then
               buffer(idx+1) = buffer(idx+1) + real ( rdx_ave(1) * (1-delta) * &
                                             fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1), p_diag_prec )
             endif
          endif


       enddo
    endif
  endif



end subroutine average1D
!--------------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------------
! This routine deposits the averaged grid values on a buffer for communication / adding to
! main data
!--------------------------------------------------------------------------------------------------
subroutine average2D( fld, fcomp, buffer, global_nx, nx_p, ave_nx, local_ave_nx, envelope )
!--------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  real( p_k_fld ), dimension(:,:,:), intent(in) :: fld        ! field data
  integer, intent(in) :: fcomp                                ! field component to average
  real( p_diag_prec ), dimension(:), intent(inout) :: buffer    ! buffer to save data to
  integer, dimension(:), intent(in) :: global_nx            ! global field size
  integer, dimension(:,:), intent(in) :: nx_p                ! local field data position on global grid
  integer, dimension(:), intent(in) :: ave_nx                ! averaged grid size
  integer, dimension(:,:), intent(in) :: local_ave_nx        ! local averaged data position on averaged grid
  logical, intent(in) :: envelope

  integer :: i, lbound, ubound
  real(p_double), dimension(rank) :: dx_ave, rdx_ave
  real(p_double) :: delta
  integer, dimension(rank) :: mod_ave
  integer, dimension(2, rank) :: grange
  integer :: Sy, gi2, ai2, ai2b, idx


  ! find local grid limits for depositions and cell normalizations
  do i = 1, rank
    dx_ave(i) = real( global_nx(i), p_double ) / ave_nx(i)
    rdx_ave(i) = 1 / dx_ave(i)

    mod_ave(i) = modulo( global_nx(i), ave_nx(i) )

    if ( mod_ave(i)  == 0 ) then
      lbound = int( (local_ave_nx( p_lower, i ) - 1)*dx_ave(i) ) + 1
      ubound = int( local_ave_nx( p_upper, i ) *dx_ave(i) )
    else
      lbound = int( (local_ave_nx( p_lower, i ) - 1)*dx_ave(i) ) + 1
      if ( local_ave_nx( p_upper,i ) < ave_nx(i) ) then
        ubound = int( local_ave_nx( p_upper, i )*dx_ave(i) ) + 1
      else
        ubound = global_nx(i)
      endif
    endif
    grange( p_lower, i ) = max( lbound, nx_p( p_lower, i ) )
    grange( p_upper, i ) = min( ubound, nx_p( p_upper, i ) )
  enddo

  Sy = local_ave_nx(p_upper,1) - local_ave_nx(p_lower,1) + 1

  if (  mod_ave(2)  == 0 ) then
     ! each column only deposits to 1 row
     do gi2 = grange( p_lower, 2 ), grange( p_upper, 2 )
        ai2 = int( (gi2-1) * rdx_ave(2) ) + 1
        idx = (ai2 - local_ave_nx(p_lower,2)) * Sy + 1
        call depline2D( buffer, idx, local_ave_nx(p_lower,1), fld, nx_p, grange(:,1), gi2, &
                                mod_ave, rdx_ave, rdx_ave(2), envelope )
     enddo
  else
     do gi2 = grange( p_lower, 2 ), grange( p_upper, 2 )
        ai2 = int( (gi2-1) * rdx_ave(2) ) + 1

        if ( gi2 < global_nx(2) ) then
          ai2b = int( gi2 * rdx_ave(2) ) + 1
        else
          ai2b = ai2
        endif

        idx = (ai2 - local_ave_nx(p_lower,2)) * Sy + 1
        if ( ai2b == ai2 ) then
           call depline2D( buffer, idx, local_ave_nx(p_lower,1), fld, nx_p, grange(:,1), gi2, &
                                   mod_ave, rdx_ave, rdx_ave(2), envelope )
        else
           ! the local global cell contributes to 2 averaged cells

           ! get fraction of global grid cell in averaged grid cell
           delta = (ai2b - 1) * dx_ave(2) - gi2 + 1

           ! check if they belong to the current chunk
           ! (this is only needed for the last dimension)
           if ( ai2 >= local_ave_nx(p_lower,2) ) then
             call depline2D( buffer, idx, local_ave_nx(p_lower,1), fld, nx_p, grange(:,1), gi2, &
                                     mod_ave, rdx_ave, rdx_ave(2) * delta, envelope )
           endif

           if ( ai2b <= local_ave_nx(p_upper,2) ) then
             idx = (ai2b - local_ave_nx(p_lower,2)) * Sy + 1
             call depline2D( buffer, idx, local_ave_nx(p_lower,1), fld, nx_p, grange(:,1), gi2, &
                                     mod_ave, rdx_ave, rdx_ave(2) * (1-delta), envelope )
           endif
        endif


     enddo
  endif


! --------

contains

  ! -------------------------------------------------
  subroutine depline2D( buffer, idx, ai1_min, fld, nx_p, gi1_range, gi2, &
                                mod_ave, rdx_ave, norm2, envelope )

     implicit none

     real( p_diag_prec ), dimension(:), intent(inout) :: buffer
     integer, intent(in) :: idx
     integer, intent(in) :: ai1_min

     real( p_k_fld ), dimension(:,:,:), intent(in) :: fld
     integer, dimension(:,:), intent(in) :: nx_p
     integer, dimension(:), intent(in) :: gi1_range
     integer, intent(in) :: gi2

     integer, dimension(:), intent(in) :: mod_ave
     real( p_double ), dimension(:), intent(in) :: rdx_ave
     real( p_double ), intent(in) :: norm2
     logical, intent(in) :: envelope


     real( p_double ) :: lnorm, delta
     integer :: gi1, ai1, ai1b, lidx
     real( p_diag_prec ) :: tmp


     if ( envelope ) then

       ! envelope data
       if ( mod_ave(1) == 0 ) then

         do gi1 = gi1_range(p_lower), gi1_range(p_upper)
            ai1 = int( (gi1-1) * rdx_ave(1) ) + 1
            lidx = idx + ai1 - ai1_min
            tmp = abs( real ( fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1, &
                                          gi2 - nx_p( p_lower, 2 ) + 1 ), p_diag_prec ) )
            if ( tmp > buffer(lidx) ) then
              buffer(lidx) = tmp
            endif
         enddo

       else

         do gi1 = gi1_range(p_lower), gi1_range(p_upper)
            ai1 = int( (gi1-1) * rdx_ave(1) ) + 1
            tmp = abs( real ( fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1, &
                                          gi2 - nx_p( p_lower, 2 ) + 1 ), p_diag_prec ) )

            if ( gi1 < global_nx(1) ) then
              ai1b = int( gi1 * rdx_ave(1) ) + 1
            else
              ai1b = ai1
            endif

            lidx = idx + ai1 - ai1_min
            if ( ai1b == ai1 ) then
              if ( tmp > buffer(lidx) ) then
                buffer(lidx) = tmp
              endif
            else
              if ( tmp > buffer(lidx) ) then
                buffer(lidx) = tmp
              endif

              if ( tmp > buffer(lidx+1) ) then
                buffer(lidx+1) = tmp
              endif
            endif
         enddo

       endif
     else

       ! average data
       lnorm = norm2 * rdx_ave(1)
       if ( mod_ave(1) == 0 ) then

         do gi1 = gi1_range(p_lower), gi1_range(p_upper)
            ai1 = int( (gi1-1) * rdx_ave(1) ) + 1
            lidx = idx + ai1 - ai1_min

            buffer(lidx) = buffer(lidx) + real ( lnorm * &
                                          fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1, &
                                                      gi2 - nx_p( p_lower, 2 ) + 1 ), p_diag_prec )
         enddo

       else

         do gi1 = gi1_range(p_lower), gi1_range(p_upper)
            ai1 = int( (gi1-1) * rdx_ave(1) ) + 1

            if ( gi1 < global_nx(1) ) then
              ai1b = int( gi1 * rdx_ave(1) ) + 1
            else
              ai1b = ai1
            endif

            lidx = idx + ai1 - ai1_min
            if ( ai1b == ai1 ) then
              buffer(lidx) = buffer(lidx) + real ( lnorm * &
                                            fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1, &
                                                        gi2 - nx_p( p_lower, 2 ) + 1 ), p_diag_prec )
            else
              delta = (ai1b - 1) * dx_ave(1) - gi1 + 1
              buffer(lidx) = buffer(lidx) + real ( lnorm * delta * &
                                            fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1, &
                                                        gi2 - nx_p( p_lower, 2 ) + 1 ), p_diag_prec )

              buffer(lidx+1) = buffer(lidx+1) + real ( lnorm * (1-delta) * &
                                            fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1, &
                                                        gi2 - nx_p( p_lower, 2 ) + 1 ), p_diag_prec )

            endif
         enddo

       endif
     endif

  end subroutine depline2D
  ! -------------------------------------------------

end subroutine average2D
!--------------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------------
! This routine deposits the averaged grid values on a buffer for communication / adding to
! main data
!--------------------------------------------------------------------------------------------------
subroutine average3D( fld, fcomp, buffer, global_nx, nx_p, ave_nx, local_ave_nx, envelope )
!--------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 3

  real( p_k_fld ), dimension(:,:,:,:), intent(in) :: fld        ! field data
  integer, intent(in) :: fcomp                                ! field component to average
  real( p_diag_prec ), dimension(:), intent(inout) :: buffer    ! buffer to save data to
  integer, dimension(:), intent(in) :: global_nx            ! global field size
  integer, dimension(:,:), intent(in) :: nx_p                ! local field data position on global grid
  integer, dimension(:), intent(in) :: ave_nx                ! averaged grid size
  integer, dimension(:,:), intent(in) :: local_ave_nx        ! local averaged data position on averaged grid
  logical, intent(in) :: envelope

  integer :: i, lbound, ubound
  real(p_double), dimension(rank) :: dx_ave, rdx_ave
  real(p_double) :: delta
  integer, dimension(rank) :: mod_ave
  integer, dimension(2, rank) :: grange
  integer :: Sy, Sz, gi3, ai3, ai3b, idx


  ! find local grid limits for depositions and cell normalizations
  do i = 1, rank
    dx_ave(i) = real( global_nx(i), p_double ) / ave_nx(i)
    rdx_ave(i) = 1 / dx_ave(i)

    mod_ave(i) = modulo( global_nx(i), ave_nx(i) )

    if ( mod_ave(i)  == 0 ) then
      lbound = int( (local_ave_nx( p_lower, i ) - 1)*dx_ave(i) ) + 1
      ubound = int( local_ave_nx( p_upper, i ) *dx_ave(i) )
    else
      lbound = int( (local_ave_nx( p_lower, i ) - 1)*dx_ave(i) ) + 1
      if ( local_ave_nx( p_upper,i ) < ave_nx(i) ) then
        ubound = int( local_ave_nx( p_upper, i )*dx_ave(i) ) + 1
      else
        ubound = global_nx(i)
      endif
    endif
    grange( p_lower, i ) = max( lbound, nx_p( p_lower, i ) )
    grange( p_upper, i ) = min( ubound, nx_p( p_upper, i ) )
  enddo

  Sy = local_ave_nx(p_upper,1) - local_ave_nx(p_lower,1) + 1
  Sz = Sy * ( local_ave_nx(p_upper,2) - local_ave_nx(p_lower,2) + 1 )


  if (  mod_ave(3)  == 0 ) then
     do gi3 = grange( p_lower, 3 ), grange( p_upper, 3 )
        ai3 = int( (gi3-1) * rdx_ave(3) ) + 1
        idx = (ai3 - local_ave_nx(p_lower,3)) * Sz + 1
        call depslab3D( buffer, idx, local_ave_nx, fld, nx_p, grange, gi3, &
                                mod_ave, rdx_ave, rdx_ave(3), envelope )
     enddo
  else
     do gi3 = grange( p_lower, 3 ), grange( p_upper, 3 )
        ai3 = int( (gi3-1) * rdx_ave(3) ) + 1

        if ( gi3 < global_nx(3) ) then
          ai3b = int( gi3 * rdx_ave(3) ) + 1
        else
          ai3b = ai3
        endif

        idx = (ai3 - local_ave_nx(p_lower,3)) * Sz + 1
        if ( ai3b == ai3 ) then
           call depslab3D( buffer, idx, local_ave_nx, fld, nx_p, grange, gi3, &
                                   mod_ave, rdx_ave, rdx_ave(3), envelope )
        else
           ! the local global cell contributes to 2 averaged cells

           ! get fraction of global grid cell in averaged grid cell
           delta = (ai3b - 1) * dx_ave(3) - gi3 + 1

           ! check if they belong to the current chunk
           ! (this is only needed for the last dimension)
           if ( ai3 >= local_ave_nx(p_lower,3) ) then
             call depslab3D( buffer, idx, local_ave_nx, fld, nx_p, grange, gi3, &
                                     mod_ave, rdx_ave, rdx_ave(3) * delta, envelope )
           endif

           if ( ai3b <= local_ave_nx(p_upper,3) ) then
             idx = (ai3b - local_ave_nx(p_lower,3)) * Sz + 1
             call depslab3D( buffer, idx, local_ave_nx, fld, nx_p, grange, gi3, &
                                     mod_ave, rdx_ave, rdx_ave(3) * (1-delta), envelope )
           endif
        endif


     enddo
  endif


! --------

contains

  ! -------------------------------------------------
  subroutine depslab3D( buffer, idx, local_ave_nx, fld, nx_p, grange, gi3, &
                                mod_ave, rdx_ave, norm3, envelope )

     implicit none

     real( p_diag_prec ), dimension(:), intent(inout) :: buffer
     integer, intent(in) :: idx
     integer, dimension(:,:), intent(in) :: local_ave_nx

     real( p_k_fld ), dimension(:,:,:,:), intent(in) :: fld
     integer, dimension(:,:), intent(in) :: nx_p
     integer, dimension(:,:), intent(in) :: grange
     integer, intent(in) :: gi3
     logical, intent(in) :: envelope

     integer, dimension(:), intent(in) :: mod_ave
     real( p_double ), dimension(:), intent(in) :: rdx_ave
     real( p_double ), intent(in) :: norm3

     integer :: Sy, gi2, ai2, ai2b, lidx
     real( p_double ) :: delta

     Sy = local_ave_nx(p_upper,1) - local_ave_nx(p_lower,1) + 1

     if (  mod_ave(2)  == 0 ) then
        do gi2 = grange( p_lower, 2 ), grange( p_upper, 2 )
           ai2 = int( (gi2-1) * rdx_ave(2) ) + 1
           lidx = idx + (ai2 - local_ave_nx(p_lower,2)) * Sy
           call depline3D( buffer, lidx, local_ave_nx, fld, nx_p, grange, gi2, gi3, &
                                        mod_ave, rdx_ave, norm3 * rdx_ave(2), envelope )
        enddo
     else
        do gi2 = grange( p_lower, 2 ), grange( p_upper, 2 )
           ai2 = int( (gi2-1) * rdx_ave(2) ) + 1

           if ( gi2 < global_nx(2) ) then
             ai2b = int( gi2 * rdx_ave(2) ) + 1
           else
             ai2b = ai2
           endif

           lidx = idx + (ai2 - local_ave_nx(p_lower,2)) * Sy
           if ( ai2b == ai2 ) then
              call depline3D( buffer, lidx, local_ave_nx, fld, nx_p, grange, gi2, gi3, &
                                      mod_ave, rdx_ave, norm3 * rdx_ave(2), envelope )
           else
              ! the local global cell contributes to 2 averaged cells

              ! get fraction of global grid cell in averaged grid cell
              delta = (ai2b - 1) * dx_ave(2) - gi2 + 1

              call depline3D( buffer, lidx, local_ave_nx, fld, nx_p, grange, gi2, gi3, &
                                      mod_ave, rdx_ave, norm3 * rdx_ave(2) * delta, envelope )

              lidx = lidx + Sy
              call depline3D( buffer, lidx, local_ave_nx, fld, nx_p, grange, gi2, gi3, &
                                      mod_ave, rdx_ave, norm3 * rdx_ave(2) * (1-delta), envelope )
           endif


        enddo
     endif


  end subroutine depslab3D

  ! -------------------------------------------------
  subroutine depline3D( buffer, idx, local_ave_nx, fld, nx_p, grange, gi2, gi3, &
                                mod_ave, rdx_ave, norm32, envelope )

     implicit none

     real( p_diag_prec ), dimension(:), intent(inout) :: buffer
     integer, intent(in) :: idx
     integer, dimension(:,:), intent(in) :: local_ave_nx

     real( p_k_fld ), dimension(:,:,:,:), intent(in) :: fld
     integer, dimension(:,:), intent(in) :: nx_p
     integer, dimension(:,:), intent(in) :: grange
     integer, intent(in) :: gi2, gi3

     integer, dimension(:), intent(in) :: mod_ave
     real( p_double ), dimension(:), intent(in) :: rdx_ave
     real( p_double ), intent(in) :: norm32
     logical, intent(in) :: envelope

     real( p_double ) :: lnorm, delta
     integer :: gi1, ai1, ai1b, lidx
     real( p_diag_prec ) :: tmp

     if ( envelope ) then
       ! envelope data
       if ( mod_ave(1) == 0 ) then

         do gi1 = grange(p_lower, 1), grange(p_upper, 1)
            ai1 = int( (gi1-1) * rdx_ave(1) ) + 1
            lidx = idx + ai1 - local_ave_nx( p_lower, 1)
            tmp = abs( real ( fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1, &
                                          gi2 - nx_p( p_lower, 2 ) + 1, &
                                          gi3 - nx_p( p_lower, 3 ) + 1), p_diag_prec ) )
            if ( tmp > buffer(lidx) ) then
               buffer(lidx) = tmp
             endif
         enddo

       else

         do gi1 = grange(p_lower, 1), grange(p_upper,1)
            ai1 = int( (gi1-1) * rdx_ave(1) ) + 1

            if ( gi1 < global_nx(1) ) then
              ai1b = int( gi1 * rdx_ave(1) ) + 1
            else
              ai1b = ai1
            endif

            lidx = idx + ai1 - local_ave_nx( p_lower, 1)
            tmp = abs( real ( fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1, &
                                          gi2 - nx_p( p_lower, 2 ) + 1, &
                                          gi3 - nx_p( p_lower, 3 ) + 1), p_diag_prec ) )

            if ( ai1b == ai1 ) then
              if ( tmp > buffer(lidx) ) then
               buffer(lidx) = tmp
              endif
            else
              if ( tmp > buffer(lidx) ) then
                buffer(lidx) = tmp
              endif

              if ( tmp > buffer(lidx+1) ) then
                buffer(lidx+1) = tmp
              endif

            endif
         enddo

       endif

     else

       ! average data
       lnorm = norm32 * rdx_ave(1)
       if ( mod_ave(1) == 0 ) then

         do gi1 = grange(p_lower, 1), grange(p_upper, 1)
            ai1 = int( (gi1-1) * rdx_ave(1) ) + 1
            lidx = idx + ai1 - local_ave_nx( p_lower, 1)

            buffer(lidx) = buffer(lidx) + real ( lnorm * &
                                          fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1, &
                                                      gi2 - nx_p( p_lower, 2 ) + 1, &
                                                      gi3 - nx_p( p_lower, 3 ) + 1), p_diag_prec )
         enddo

       else

         do gi1 = grange(p_lower, 1), grange(p_upper,1)
            ai1 = int( (gi1-1) * rdx_ave(1) ) + 1

            if ( gi1 < global_nx(1) ) then
              ai1b = int( gi1 * rdx_ave(1) ) + 1
            else
              ai1b = ai1
            endif

            lidx = idx + ai1 - local_ave_nx( p_lower, 1)
            if ( ai1b == ai1 ) then
              buffer(lidx) = buffer(lidx) + real ( lnorm * &
                                            fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1, &
                                                        gi2 - nx_p( p_lower, 2 ) + 1, &
                                                        gi3 - nx_p( p_lower, 3 ) + 1 ), p_diag_prec )
            else
              delta = (ai1b - 1) * dx_ave(1) - gi1 + 1
              buffer(lidx) = buffer(lidx) + real ( lnorm * delta * &
                                            fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1, &
                                                        gi2 - nx_p( p_lower, 2 ) + 1, &
                                                        gi3 - nx_p( p_lower, 3 ) + 1 ), p_diag_prec )

              buffer(lidx+1) = buffer(lidx+1) + real ( lnorm * (1-delta) * &
                                            fld( fcomp, gi1 - nx_p( p_lower, 1 ) + 1, &
                                                        gi2 - nx_p( p_lower, 2 ) + 1, &
                                                        gi3 - nx_p( p_lower, 3 ) + 1), p_diag_prec )

            endif
         enddo

       endif
     endif

  end subroutine depline3D
  ! -------------------------------------------------

end subroutine average3D
!--------------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------------
!  save grid average / envelope of the specified field component
!--------------------------------------------------------------------------------------------------
subroutine report_ave_vdf( report, vdf, fcomp, g_space, grid, no_co, &
                           n_ave, envelope )
!--------------------------------------------------------------------------------------------------

  use m_vdf_reportfile

  implicit none

  type( t_vdf ), intent(in)         :: vdf      ! vdf object
  integer, intent(in)               :: fcomp    ! field component to report
  type( t_space ), intent(in)       :: g_space  ! spatial information
  class( t_grid ), intent(in)        :: grid     ! global grid info
  class( t_node_conf ), intent(in)   :: no_co    ! node configuration of the object

  type( t_vdf_report ), intent(in)  :: report   ! report file details

  integer, dimension(:), intent(in) :: n_ave ! # average points
  logical, intent(in)               :: envelope ! do envelope instead of average


  integer, dimension(3) :: ave_nx, global_nx
  integer, dimension(2, 3) :: range, local_ave_nx, remote_ave_nx
  integer :: ave_size, bsize, ldatasize, rdatasize
  integer :: nchunks, chunksize
  integer :: i, i1, i2, i3, node
  integer :: idx_local, idx_remote, idx_total, Sy, Sz
  integer, dimension(2) :: chunk

  integer, dimension(3,3) :: lnx_p

  real(p_diag_prec), dimension(:), pointer :: buffer_full, buffer1, buffer2
  integer :: ierr

  class( t_diag_file ) , allocatable :: diagFile
  type( t_diag_dataset ) :: dataset
  type( t_diag_chunk )   :: dset_chunk

  ! MPI variables
  integer :: ping, request, mpi_rnode
  integer, dimension(mpi_status_size) :: status

  ! get averaged grid global dimensions and total size
  ave_size = 1
  do i = 1, vdf%x_dim_
    global_nx(i) = grid%g_nx( i )
    ave_nx(i) =  global_nx(i)/ n_ave(i)
    ave_size = ave_size * ave_nx(i)
  enddo


  ! get buffer for diagnostics (we actually need 2 buffers so we divide the buffer in 2)
  call get_diag_buffer( buffer_full, bsize )
  bsize = bsize / 2
  buffer1 => buffer_full(1:bsize)
  buffer2 => buffer_full(bsize+1:)

  ! get number of chunks
  nchunks = ceiling( float( ave_size ) / bsize )
  if ( nchunks > ave_nx( vdf%x_dim_ ) ) then
    if ( root( no_co ) ) then
      write(0,*) '(*error*) Averaged grid is too large, n_ave = ', n_ave(1:vdf%x_dim_)
      write(0,*) '(*error*) Try increasing n_ave or increasing diag. buffer size'
    endif
    call abort_program( p_err_alloc )
  endif
  chunksize = ave_nx( vdf%x_dim_ ) / nchunks

  ! create file and dataset
  if (root( no_co )) then

    ! Create the report file
    call init_diag_file_vdf( vdf, report, g_space, grid, diagFile )

    ! Correct the grid size
    do i = 1, vdf%x_dim_
      diagFile%grid%count(i) = ave_nx(i)
    enddo

    ! Open file in single node access
    call diagFile % open( p_diag_create )

    ! Create dataset
    call diagFile % start_cdset( report%name, vdf%x_dim_, ave_nx, diag_type_real( p_diag_prec ), dataset )

  endif

  lnx_p(:,1:vdf%x_dim_) = grid%my_nx(:,1:vdf%x_dim_)

  select case ( vdf%x_dim_ )

    case (1) ! --------------------------------------------------- 1D data ----

      ! loop through all chunks
      chunk(1) = 1
      do
        if ( chunk(1) > ave_nx(1) ) exit

        chunk(2) = min( chunk(1) + chunkSize - 1, ave_nx(1) )
        bsize = ( chunk(2) - chunk(1) + 1 )
        buffer1(1:bsize) = 0
        buffer2(1:bsize) = 0

        range( 1:2 , 1 ) = chunk

        call aveDataRange( global_nx, lnx_p, 1, ave_nx, range, &
                           local_ave_nx, ldatasize )

        ! first node
        if (root( no_co )) then

          ! take care of local data if any
          if ( ldatasize > 0 ) then
            call average1D( vdf%f1(:,1:vdf%nx_(1)), fcomp, buffer2, global_nx, &
                             lnx_p, ave_nx, local_ave_nx, envelope )

             idx_local = 1
             if ( envelope ) then
               ! envelope data
               do i1 = local_ave_nx( p_lower, 1 ), local_ave_nx( p_upper, 1 )
                 idx_total = ( i1 - range( p_lower, 1 ) ) + 1
                 if ( buffer2( idx_local ) > buffer1( idx_total ) ) then
                   buffer1( idx_total ) = buffer2( idx_local )
                 endif
                 idx_local = idx_local + 1
               enddo
             else
               ! average data
               do i1 = local_ave_nx( p_lower, 1 ), local_ave_nx( p_upper, 1 )
                 idx_total = ( i1 - range( p_lower, 1 ) ) + 1
                 buffer1( idx_total ) = buffer1( idx_total ) + buffer2( idx_local )
                 idx_local = idx_local + 1
               enddo
             endif

          endif

          ! receive data from other nodes
          do node = 2, no_num( no_co )

             call aveDataRange( grid%g_nx, nx_p( grid, no_co, node ), 1, ave_nx, range, &
                                remote_ave_nx, rdatasize )

             ! get 0 indexed rank of target node
             mpi_rnode = node - 1

             if (  rdatasize > 0 ) then
                ! post receive
                call mpi_irecv( buffer2, bsize, mpi_real_type( p_diag_prec ), mpi_rnode, 0, comm(no_co), request, ierr )
             endif

             ! send ping to other node

! temporary fix because running out of buffer space on some systems
#if 0
             ! using a buffered send so the communication can be started before the remote
             ! receive is posted
             call mpi_bsend( 0, 1, MPI_INTEGER, mpi_rnode, 0, comm(no_co), ierr)
#else
             ! use a standard send
             call mpi_send( 0, 1, MPI_INTEGER, mpi_rnode, 0, comm(no_co), ierr)
#endif

            if (  rdatasize > 0 ) then
               ! wait for data to arrive
               call mpi_wait( request, status, ierr )
               idx_remote = 1
               if ( envelope ) then
                 ! envelope data
                 do i1 = remote_ave_nx( p_lower, 1 ), remote_ave_nx( p_upper, 1 )
                   idx_total = ( i1 - range( p_lower, 1 ) ) + 1
                   if ( buffer2( idx_remote ) > buffer1( idx_total )) then
                     buffer1( idx_total ) =  buffer2( idx_remote )
                   endif
                   idx_remote = idx_remote + 1
                 enddo
               else
                 ! average data
                 do i1 = remote_ave_nx( p_lower, 1 ), remote_ave_nx( p_upper, 1 )
                   idx_total = ( i1 - range( p_lower, 1 ) ) + 1
                   buffer1( idx_total ) = buffer1( idx_total ) + buffer2( idx_remote )
                   idx_remote = idx_remote + 1
                 enddo
               endif
            endif

          enddo

          ! data is ready in buffer1, write chunk

          ! hyperslab coordinates are 0 indexed
          dset_chunk % start(1)  = chunk(1) - 1
          dset_chunk % count(1)  = chunk(2) - chunk(1) + 1
          dset_chunk % stride(1) = 1

          dset_chunk % data = c_loc( buffer1 )

          call diagFile % write_cdset( dataset, dset_chunk )

        else

          ! other nodes

          ! post receive for ping
          call mpi_irecv( ping, 1, MPI_INTEGER, 0, 0, comm( no_co ), request, ierr )

          ! take care of local data if any
          if ( ldatasize > 0 ) then
             ! get average values
             call average1D( vdf%f1(:,1:vdf%nx_(1)), fcomp, buffer1, global_nx, &
                             lnx_p, ave_nx, local_ave_nx, envelope )

             ! wait for ping to complete (this ensures that node 0 is waiting for the data)
             call mpi_wait( request, status, ierr )

             ! send to node 0
             call mpi_send( buffer1, ldatasize, mpi_real_type( p_diag_prec ), 0, 0, comm( no_co ), ierr )

          else

             ! wait for ping to complete (this keeps nodes in sync)
             call mpi_wait( request, status, ierr )

          endif


        endif

        ! Process next chunk
        chunk(1) = chunk(1) + chunkSize
      enddo

    case (2) ! --------------------------------------------------- 2D data ----

      ! loop through all chunks
      chunk(1) = 1
      do
        if ( chunk(1) > ave_nx(2) ) exit

        chunk(2) = min( chunk(1) + chunkSize - 1, ave_nx(2) )
        bsize = ave_nx(1) * ( chunk(2) - chunk(1) + 1 )
        buffer1(1:bsize) = 0
        buffer2(1:bsize) = 0

        range( 1:2 , 1 ) = (/ 1, ave_nx(1) /)
        range( 1:2 , 2 ) = chunk
        Sy = range(p_upper, 1) - range(p_lower, 1) + 1

        call aveDataRange( global_nx, lnx_p, 2, ave_nx, range, &
                           local_ave_nx, ldatasize )

        ! first node
        if (root( no_co )) then

          ! take care of local data if any
          if ( ldatasize > 0 ) then
            call average2D( vdf%f2(:,1:vdf%nx_(1),1:vdf%nx_(2)), fcomp, buffer2, &
                             global_nx, lnx_p, ave_nx, local_ave_nx, envelope )

             idx_local = 1
             if ( envelope ) then
               ! envelope data
               do i2 = local_ave_nx( p_lower, 2 ), local_ave_nx( p_upper, 2 )
                  do i1 = local_ave_nx( p_lower, 1 ), local_ave_nx( p_upper, 1 )
                    idx_total = ( i1 - range( p_lower, 1 ) ) + &
                                ( i2 - range( p_lower, 2 ) ) * Sy + 1
                    if ( buffer2( idx_local ) > buffer1( idx_total )) then
                      buffer1( idx_total ) =  buffer2( idx_local )
                    endif
                    idx_local = idx_local + 1
                  enddo
               enddo
             else
               ! average data
               do i2 = local_ave_nx( p_lower, 2 ), local_ave_nx( p_upper, 2 )
                  do i1 = local_ave_nx( p_lower, 1 ), local_ave_nx( p_upper, 1 )
                    idx_total = ( i1 - range( p_lower, 1 ) ) + &
                                ( i2 - range( p_lower, 2 ) ) * Sy + 1
                    buffer1( idx_total ) = buffer1( idx_total ) + buffer2( idx_local )
                    idx_local = idx_local + 1
                  enddo
               enddo
             endif

          endif

          ! receive data from other nodes
          do node = 2, no_num( no_co )

             call aveDataRange( grid%g_nx, nx_p( grid, no_co, node ), 2, ave_nx, range, &
                                remote_ave_nx, rdatasize )

             ! get 0 indexed rank of target node
             mpi_rnode = node - 1

             if (  rdatasize > 0 ) then
                ! post receive
                call mpi_irecv( buffer2, bsize, mpi_real_type( p_diag_prec ), mpi_rnode, 0, comm(no_co), request, ierr )
             endif

             ! send ping to other node

! temporary fix because running out of buffer space on some systems
#if 0
             ! using a buffered send so the communication can be started before the remote
             ! receive is posted
             call mpi_bsend( 0, 1, MPI_INTEGER, mpi_rnode, 0, comm(no_co), ierr)
#else
             ! use a standard send
             call mpi_send( 0, 1, MPI_INTEGER, mpi_rnode, 0, comm(no_co), ierr)
#endif

            if (  rdatasize > 0 ) then
               ! wait for data to arrive
               call mpi_wait( request, status, ierr )
               idx_remote = 1
               if (envelope) then
                 ! envelope data
                 do i2 = remote_ave_nx( p_lower, 2 ), remote_ave_nx( p_upper, 2 )
                    do i1 = remote_ave_nx( p_lower, 1 ), remote_ave_nx( p_upper, 1 )
                      idx_total = ( i1 - range( p_lower, 1 ) ) + &
                                  ( i2 - range( p_lower, 2 ) ) * Sy + 1
                      if ( buffer2( idx_remote ) > buffer1( idx_total )) then
                         buffer1( idx_total ) =  buffer2( idx_remote )
                      endif
                      idx_remote = idx_remote + 1
                    enddo
                 enddo
               else
                 ! average data
                 do i2 = remote_ave_nx( p_lower, 2 ), remote_ave_nx( p_upper, 2 )
                    do i1 = remote_ave_nx( p_lower, 1 ), remote_ave_nx( p_upper, 1 )
                      idx_total = ( i1 - range( p_lower, 1 ) ) + &
                                  ( i2 - range( p_lower, 2 ) ) * Sy + 1
                      buffer1( idx_total ) = buffer1( idx_total ) + buffer2( idx_remote )
                      idx_remote = idx_remote + 1
                    enddo
                 enddo
               endif
            endif

          enddo

          ! data is ready in buffer1, write chunk

          ! hyperslab coordinates are 0 indexed
          dset_chunk % start(1:2)  = [ 0, chunk(1) - 1 ]
          dset_chunk % count(1:2)  = [ ave_nx(1), chunk(2) - chunk(1) + 1 ]
          dset_chunk % stride(1:2) = 1

          dset_chunk % data = c_loc( buffer1 )

          call diagFile % write_cdset( dataset, dset_chunk )

        else

          ! other nodes

          ! post receive for ping
          call mpi_irecv( ping, 1, MPI_INTEGER, 0, 0, comm( no_co ), request, ierr )

          ! take care of local data if any
          if ( ldatasize > 0 ) then
             ! get average values
             call average2D( vdf%f2(:,1:vdf%nx_(1),1:vdf%nx_(2)), fcomp, buffer1, &
                             global_nx, lnx_p, ave_nx, local_ave_nx, envelope )

             ! wait for ping to complete (this ensures that node 0 is waiting for the data)
             call mpi_wait( request, status, ierr )

             ! send to node 0
             call mpi_send( buffer1, ldatasize, mpi_real_type( p_diag_prec ), 0, 0, comm( no_co ), ierr )

          else

             ! wait for ping to complete (this keeps nodes in sync)
             call mpi_wait( request, status, ierr )

          endif


        endif

        ! Process next chunk
        chunk(1) = chunk(1) + chunkSize
      enddo

    case (3) ! --------------------------------------------------- 3D data ----


      ! loop through all chunks
      chunk(1) = 1
      do
        if ( chunk(1) > ave_nx(3) ) exit

        chunk(2) = min( chunk(1) + chunkSize - 1, ave_nx(3) )
        bsize = ave_nx(1) * ave_nx(2) * ( chunk(2) - chunk(1) + 1 )
        buffer1(1:bsize) = 0
        buffer2(1:bsize) = 0

        range( 1:2 , 1 ) = (/ 1, ave_nx(1) /)
        range( 1:2 , 2 ) = (/ 1, ave_nx(2) /)
        range( 1:2 , 3 ) = chunk

        Sy = range(p_upper, 1) - range(p_lower, 1) + 1
        Sz = Sy * ( range(p_upper, 2) - range(p_lower, 2) + 1 )

        call aveDataRange( global_nx, lnx_p, 3, ave_nx, range, &
                           local_ave_nx, ldatasize )

        ! first node
        if (root( no_co )) then

          ! take care of local data if any
          if ( ldatasize > 0 ) then
            call average3D( vdf%f3(:,1:vdf%nx_(1),1:vdf%nx_(2),1:vdf%nx_(3)), fcomp, &
                             buffer2, global_nx, lnx_p, ave_nx, local_ave_nx, envelope )

             idx_local = 1
             if ( envelope ) then
               ! envelope data
               do i3 = local_ave_nx( p_lower, 3 ), local_ave_nx( p_upper, 3 )
                 do i2 = local_ave_nx( p_lower, 2 ), local_ave_nx( p_upper, 2 )
                    do i1 = local_ave_nx( p_lower, 1 ), local_ave_nx( p_upper, 1 )
                      idx_total = ( i1 - range( p_lower, 1 ) ) + &
                                  ( i2 - range( p_lower, 2 ) ) * Sy + &
                                  ( i3 - range( p_lower, 3 ) ) * Sz + 1
                      if ( buffer2( idx_local ) > buffer1( idx_total )) then
                         buffer1( idx_total ) =  buffer2( idx_local )
                      endif
                      idx_local = idx_local + 1
                    enddo
                 enddo
               enddo
             else
               ! average data
               do i3 = local_ave_nx( p_lower, 3 ), local_ave_nx( p_upper, 3 )
                 do i2 = local_ave_nx( p_lower, 2 ), local_ave_nx( p_upper, 2 )
                    do i1 = local_ave_nx( p_lower, 1 ), local_ave_nx( p_upper, 1 )
                      idx_total = ( i1 - range( p_lower, 1 ) ) + &
                                  ( i2 - range( p_lower, 2 ) ) * Sy + &
                                  ( i3 - range( p_lower, 3 ) ) * Sz + 1
                      buffer1( idx_total ) = buffer1( idx_total ) + buffer2( idx_local )
                      idx_local = idx_local + 1
                    enddo
                 enddo
               enddo
             endif

          endif

          ! receive data from other nodes
          do node = 2, no_num( no_co )

             call aveDataRange( grid%g_nx, nx_p( grid, no_co, node ), 3, ave_nx, range, &
                                remote_ave_nx, rdatasize )

             ! get 0 indexed rank of target node
             mpi_rnode = node - 1

             if (  rdatasize > 0 ) then
                ! post receive
                call mpi_irecv( buffer2, bsize, mpi_real_type( p_diag_prec ), mpi_rnode, 0, comm(no_co), request, ierr )
             endif

             ! send ping to other node

! temporary fix because running out of buffer space on some systems
#if 0
             ! using a buffered send so the communication can be started before the remote
             ! receive is posted
             call mpi_bsend( 0, 1, MPI_INTEGER, mpi_rnode, 0, comm(no_co), ierr)
#else
             ! use a standard send
             call mpi_send( 0, 1, MPI_INTEGER, mpi_rnode, 0, comm(no_co), ierr)
#endif

            if (  rdatasize > 0 ) then
               ! wait for data to arrive
               call mpi_wait( request, status, ierr )
               idx_remote = 1
               if ( envelope ) then
                 ! envelope data
                 do i3 = remote_ave_nx( p_lower, 3 ), remote_ave_nx( p_upper, 3 )
                   do i2 = remote_ave_nx( p_lower, 2 ), remote_ave_nx( p_upper, 2 )
                      do i1 = remote_ave_nx( p_lower, 1 ), remote_ave_nx( p_upper, 1 )
                        idx_total = ( i1 - range( p_lower, 1 ) ) + &
                                    ( i2 - range( p_lower, 2 ) ) * Sy + &
                                    ( i3 - range( p_lower, 3 ) ) * Sz + 1
                        if ( buffer2( idx_remote ) > buffer1( idx_total )) then
                          buffer1( idx_total ) =  buffer2( idx_remote )
                        endif
                        idx_remote = idx_remote + 1
                      enddo
                   enddo
                 enddo
               else
                 ! average data
                 do i3 = remote_ave_nx( p_lower, 3 ), remote_ave_nx( p_upper, 3 )
                   do i2 = remote_ave_nx( p_lower, 2 ), remote_ave_nx( p_upper, 2 )
                      do i1 = remote_ave_nx( p_lower, 1 ), remote_ave_nx( p_upper, 1 )
                        idx_total = ( i1 - range( p_lower, 1 ) ) + &
                                    ( i2 - range( p_lower, 2 ) ) * Sy + &
                                    ( i3 - range( p_lower, 3 ) ) * Sz + 1
                        buffer1( idx_total ) = buffer1( idx_total ) + buffer2( idx_remote )
                        idx_remote = idx_remote + 1
                      enddo
                   enddo
                 enddo
               endif


            endif

          enddo

          ! data is ready in buffer1, write chunk

          ! hyperslab coordinates are 0 indexed
          dset_chunk % start(1:3)  = [ 0, 0, chunk(1) - 1 ]
          dset_chunk % count(1:3)  = [ ave_nx(1), ave_nx(2), chunk(2) - chunk(1) + 1 ]
          dset_chunk % stride(1:3) = 1

          dset_chunk % data = c_loc( buffer1 )

          call diagFile % write_cdset( dataset, dset_chunk )

        else

          ! other nodes

          ! post receive for ping
          call mpi_irecv( ping, 1, MPI_INTEGER, 0, 0, comm( no_co ), request, ierr )

          ! take care of local data if any
          if ( ldatasize > 0 ) then
             ! get average values
             call average3D( vdf%f3(:,1:vdf%nx_(1),1:vdf%nx_(2),1:vdf%nx_(3)), fcomp, &
                             buffer1, global_nx, lnx_p, ave_nx, local_ave_nx, envelope )

             ! wait for ping to complete (this ensures that node 0 is waiting for the data)
             call mpi_wait( request, status, ierr )

             ! send to node 0
             call mpi_send( buffer1, ldatasize, mpi_real_type( p_diag_prec ), 0, 0, comm( no_co ), ierr )

          else

             ! wait for ping to complete (this keeps nodes in sync)
             call mpi_wait( request, status, ierr )

          endif


        endif

        ! Process next chunk
        chunk(1) = chunk(1) + chunkSize
      enddo


  end select

  if (root( no_co )) then
     ! close dataset on file
     call diagFile % end_cdset( dataset )

     ! close the file
     call diagFile % close( )

     deallocate( diagFile )
  endif


end subroutine report_ave_vdf
!--------------------------------------------------------------------------------------------------


end module m_vdf_average
