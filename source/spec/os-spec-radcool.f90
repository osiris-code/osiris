! m_species_radcool module
! 
! Handles particle advance / deposit for the radiation cooling algorithm
!
! Note u is considered to be time centered with positions (so all quantities are time
! centered)

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_radcool

use m_emf
use m_emf_interpolate

use m_species_define
use m_emf_define
use m_vdf_define

use m_system
use m_parameters

use m_species_current

use m_units

private

interface advance_deposit_radcool
  module procedure advance_deposit_radcool
end interface

public :: advance_deposit_radcool

! -----------------------------------------------------------------------------
contains

!---------------------------------------------------------------------------------------------------
function ntrim(x)
!---------------------------------------------------------------------------------------------------
! Returns the integer shift (-1, 0 or +1) so that the coordinate remains in the [-0.5, 0.5[
! range. This is the fastest implementation (twice as fast as a sequence of ifs) because
! the two if structures compile as conditional moves and can be processed independently.
! This has no precision problem and is only 12% slower than the previous "int(x+1.5)-1"
! routine that would break for x = nearest( 0.5, -1.0 ) 
!---------------------------------------------------------------------------------------------------
  implicit none
  
  real(p_k_part), intent(in) :: x
  integer :: ntrim, a, b

  if ( x < -0.5_p_k_part ) then 
    a = -1
  else 
    a = 0
  endif

  if ( x >= 0.5_p_k_part ) then 
    b = +1
  else 
    b = 0
  endif
  
  ntrim = a+b

end function ntrim
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Push particles and deposit electric current
!---------------------------------------------------------------------------------------------------
subroutine advance_deposit_radcool( this, emf, jay, t, tstep, omega_p0, tid, n_threads )
  
  use m_time_step
  
  implicit none

  class( t_species ), intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  type( t_vdf ), intent(inout) :: jay

  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  
  real(p_double), intent(in) :: omega_p0
  
  integer, intent(in) :: tid        ! local thread id
  integer, intent(in) :: n_threads  ! total number of threads
  
  ! local variables

  real(p_double) :: dtcycle
  integer :: chunk, i0, i1
  
  ! if before push start time return silently
  if ( t < this%push_start_time ) return


   ! Time centered energy diagnostic is not yet available, force recalculation later
   this%energy = p_ene_recalc

   dtcycle = real( dt(tstep), p_k_part ) 

   ! range of particles for each thread
   chunk = ( this%num_par + n_threads - 1 ) / n_threads
   i0    = tid * chunk + 1
   i1    = min( (tid+1) * chunk, this%num_par ) 

   ! Push particles. Boundary crossings will be checked at update_boundary

   select case ( this%coordinates )

   case default
      select case ( p_x_dim )

      !case (1)
      !  call advance_deposit_radcool_1d( this, emf, jay, dtcycle, omega_p0, i0, i1 )

      case (2)
         call advance_deposit_radcool_2d( this, emf, jay, dtcycle, omega_p0, i0, i1 )

      case (3)
         call advance_deposit_radcool_3d( this, emf, jay, dtcycle, omega_p0, i0, i1 )

      case default
         ERROR('Not implemented for x_dim = ',p_x_dim)
         call abort_program(p_err_invalid)

      end select

   case ( p_cylindrical_b )

      ! call advance_deposit_cyl_2d( this, emf, jay, dtcycle, i0, i1 )
      
      ERROR('Not implemented for cylindrical geometry')
      call abort_program(p_err_invalid)

   end select
  
end subroutine advance_deposit_radcool
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine advance_deposit_radcool_2d( this, emf, jay, gdt, omega_p0, i0, i1 )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  ! dummy variables

  class( t_species ),    intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  type( t_vdf ),        intent(inout) :: jay
  real(p_double),   intent(in) :: gdt, omega_p0
  integer, intent(in) :: i0, i1

!       local variables
  real(p_k_part), dimension(p_x_dim) :: rdx 
  integer :: i, pp, np, ptrcur, dxi1, dxi2


  real(p_k_part), dimension(p_x_dim,p_cache_size) :: xbuf, xtemp, xfinal
  integer, dimension(p_x_dim,p_cache_size)        :: dxi, ixtemp, ixfinal
  real(p_k_part), dimension(p_cache_size)      :: rg

  real(p_k_part), dimension(p_cache_size) :: q_weight_temp             !probably not necessary
  
  ! variables from du_dt
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp, ufinal
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: k_p_rung          ! Runge-Kutta coefficients, 3p
  real(p_k_part), dimension(p_x_dim,p_cache_size) :: k_x_rung            ! 3x
  real(p_k_part), dimension(p_cache_size) :: gam_temp

  real(p_k_part) :: Cdamp, tem, tem_x
  

  ! cell size 
  rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
  rdx(2) = real( 1.0_p_double/this%dx(2), p_k_part )
  
  ! get factor that includes timestep and charge-to-mass ratio
  ! note that the charge-to-mass ratio is used since the
  ! momentum p is not the total momentum of a particles but
  ! the momentum per unit restmass (which is the electron mass)

  tem = real( 0.5_p_double * gdt / this%rqm, p_k_part )
  tem_x = real(0.5_p_double * gdt, p_k_part)
  
  Cdamp = real(2.0_p_k_part * cgs_e * cgs_e * omega_p0 / &
                 (3.0_p_k_part * cgs_me * cgs_c**3 ) )

  Cdamp = Cdamp * this%q_real**2 / this%rqm 

  
  ! advance particles
  do ptrcur = i0, i1, p_cache_size

     ! check if last copy of table and set np
     if( ptrcur + p_cache_size > i1 ) then
         np = i1 - ptrcur + 1
     else
         np = p_cache_size
     endif
     
     !---------------------------------------------------------------
     !--------------------First Runge-Kutta step---------------------
     
     ! Get fields at particle position
     call this % get_emf( emf, bp, ep, np, ptrcur )
     
     ! Take the initial momenta         
     pp = ptrcur
     do i=1,np
       utemp(1,i) = this%p(1,pp) 
       utemp(2,i) = this%p(2,pp) 
       utemp(3,i) = this%p(3,pp) 

       gam_temp(i)= 1.0_p_k_part / sqrt(1.0_p_k_part+utemp(1,i)**2+ &
                                    utemp(2,i)**2+ &
                                    utemp(3,i)**2)

       pp = pp + 1
     end do
         
     ! Take initial coordinates
     pp = ptrcur
     do i=1,np
        xtemp(1,i) = this%x(1,pp)
        xtemp(2,i) = this%x(2,pp)

        ixtemp(1,i) = this%ix(1,pp)
        ixtemp(2,i) = this%ix(2,pp)
        
        q_weight_temp(i) = this%q(pp)

        pp = pp + 1
     end do
         
     ! Calculate the first set of Runge-Kutta coefficients
     call rk_get(ep,bp,utemp,gam_temp,np,Cdamp,k_p_rung,k_x_rung,q_weight_temp)  
            
            
     ! Update the final values with first correction
     pp = ptrcur
     do i=1,np
        ufinal(1,i) = this%p(1,pp)+tem*k_p_rung(1,i) * 0.3333333333333333_p_k_part              
        ufinal(2,i) = this%p(2,pp)+tem*k_p_rung(2,i) * 0.3333333333333333_p_k_part
        ufinal(3,i) = this%p(3,pp)+tem*k_p_rung(3,i) * 0.3333333333333333_p_k_part
        
        pp = pp + 1
     end do
            
            
     ! Final value for coordinates update. 
     ! Trimmed only at the end of the whole timestep

     pp = ptrcur
     do i=1,np
        xfinal(1,i) = this%x(1,pp)+tem_x*k_x_rung(1,i) * rdx(1) * 0.3333333333333333_p_k_part
        xfinal(2,i) = this%x(2,pp)+tem_x*k_x_rung(2,i) * rdx(2) * 0.3333333333333333_p_k_part             
     
        ixfinal(1,i) = this%ix(1,pp)
        ixfinal(2,i) = this%ix(2,pp)
       
        pp=pp+1
     end do         
           
     !----------------------------------------------------------------
     !--------------------Second Runge-Kutta step---------------------
     
     ! Update the temporary values of momenta for the next step
     pp = ptrcur
     do i=1,np
        utemp(1,i) = this%p(1,pp)+tem*k_p_rung(1,i)              
        utemp(2,i) = this%p(2,pp)+tem*k_p_rung(2,i)
        utemp(3,i) = this%p(3,pp)+tem*k_p_rung(3,i)
        
        gam_temp(i)= 1.0_p_k_part / sqrt(1.0_p_k_part+utemp(1,i)**2+ &
                                    utemp(2,i)**2+ &
                                    utemp(3,i)**2)
        pp = pp + 1
     end do
           
     ! Updating the coordinates on the way to the second step
         
     pp = ptrcur
     do i=1,np
        xtemp(1,i) = this%x(1,pp)+tem_x*k_x_rung(1,i)*rdx(1)
        xtemp(2,i) = this%x(2,pp)+tem_x*k_x_rung(2,i)*rdx(2)

        dxi1 = ntrim( xtemp(1,i) )    
        dxi2 = ntrim( xtemp(2,i) )
          
        ixtemp(1,i) = this%ix(1,pp) + dxi1
        ixtemp(2,i) = this%ix(2,pp) + dxi2
          
        xtemp(1,i) = xtemp(1,i) - dxi1
        xtemp(2,i) = xtemp(2,i) - dxi2

        pp = pp + 1
     end do

     ! getting the new fields
     !!!!!!!!!!!!!!!!!!!!
     call get_emf( emf, bp, ep, ixtemp, xtemp, np, this%interpolation )
               
            
     ! Calculate the second set of Runge-Kutta coefficients    
     call rk_get(ep,bp,utemp,gam_temp,np,Cdamp,k_p_rung,k_x_rung,q_weight_temp)  
        
     ! Update the final momenta with second correction
     do i=1,np
         ufinal(1,i) = ufinal(1,i)+ tem*k_p_rung(1,i)*0.6666666666666666_p_k_part             
         ufinal(2,i) = ufinal(2,i)+ tem*k_p_rung(2,i)*0.6666666666666666_p_k_part   
         ufinal(3,i) = ufinal(3,i)+ tem*k_p_rung(3,i)*0.6666666666666666_p_k_part   
     end do
            
     ! Update the coordinates with the corrections
     ! Again, we do not trim the coordinates
     pp = ptrcur
     do i=1,np
         xfinal(1,i) = xfinal(1,i) + tem_x * k_x_rung(1,i) * rdx(1) * 0.6666666666666666_p_k_part  
         xfinal(2,i) = xfinal(2,i) + tem_x * k_x_rung(2,i) * rdx(2) * 0.6666666666666666_p_k_part           
       
         pp=pp+1
     end do         

     !----------------------------------------------------------------
     !--------------------Third Runge-Kutta step----------------------
    
     ! Update the temporary values of x, gamma and p for the next step
     
    
     pp = ptrcur
     do i=1,np
        utemp(1,i) = this%p(1,pp)+tem*k_p_rung(1,i)              
        utemp(2,i) = this%p(2,pp)+tem*k_p_rung(2,i)
        utemp(3,i) = this%p(3,pp)+tem*k_p_rung(3,i)
        
        gam_temp(i)= 1.0_p_k_part / sqrt(1.0_p_k_part+utemp(1,i)**2+ &
                                    utemp(2,i)**2+ &
                                    utemp(3,i)**2)
        pp = pp + 1
     end do
       
     ! Updating the values of the coordinates 
     pp = ptrcur
     do i=1,np
        xtemp(1,i) = this%x(1,pp)+tem_x*k_x_rung(1,i)*rdx(1)
        xtemp(2,i) = this%x(2,pp)+tem_x*k_x_rung(2,i)*rdx(2)

        dxi1 = ntrim( xtemp(1,i) )    
        dxi2 = ntrim( xtemp(2,i) )
          
        ixtemp(1,i) = this%ix(1,pp) + dxi1
        ixtemp(2,i) = this%ix(2,pp) + dxi2
          
        xtemp(1,i) = xtemp(1,i) - dxi1
        xtemp(2,i) = xtemp(2,i) - dxi2

        pp = pp + 1
     end do
        
     ! getting the new fields
     !!!!!!!!!!!!!!!!!!!!
     call get_emf( emf, bp, ep, ixtemp, xtemp, np, this%interpolation )
       
     ! Calculate the third set of Runge-Kutta coefficients    
     call rk_get(ep,bp,utemp,gam_temp,np,Cdamp,k_p_rung,k_x_rung,q_weight_temp) 
        
     ! Update the final values with the third correction
     ! Still no trimming 
      do i=1,np
         ufinal(1,i) = ufinal(1,i)+tem*k_p_rung(1,i)*0.6666666666666666_p_k_part                
         ufinal(2,i) = ufinal(2,i)+tem*k_p_rung(2,i)*0.6666666666666666_p_k_part  
         ufinal(3,i) = ufinal(3,i)+tem*k_p_rung(3,i)*0.6666666666666666_p_k_part  
      end do
            
      do i=1,np
         xfinal(1,i) = xfinal(1,i)+tem_x * k_x_rung(1,i) * rdx(1) * 0.6666666666666666_p_k_part  
         xfinal(2,i) = xfinal(2,i)+tem_x * k_x_rung(2,i) * rdx(2) * 0.6666666666666666_p_k_part  
      end do         

      !----------------------------------------------------------------
      !--------------------Fourth Runge-Kutta step---------------------
      
      pp = ptrcur
      do i=1,np
         utemp(1,i) = this%p(1,pp)+2.0_p_k_part*tem*k_p_rung(1,i)              
         utemp(2,i) = this%p(2,pp)+2.0_p_k_part*tem*k_p_rung(2,i)
         utemp(3,i) = this%p(3,pp)+2.0_p_k_part*tem*k_p_rung(3,i)
         
         gam_temp(i)= 1.0_p_k_part / sqrt(1.0_p_k_part+utemp(1,i)**2+ &
                                     utemp(2,i)**2+ &
                                     utemp(3,i)**2)
         pp = pp + 1
      end do
           
            
      pp = ptrcur
      do i=1,np
         xtemp(1,i) = this%x(1,pp)+2.0_p_k_part*tem_x*k_x_rung(1,i)*rdx(1)
         xtemp(2,i) = this%x(2,pp)+2.0_p_k_part*tem_x*k_x_rung(2,i)*rdx(2)

         dxi1 = ntrim( xtemp(1,i) )    
         dxi2 = ntrim( xtemp(2,i) )
           
         ixtemp(1,i) = this%ix(1,pp) + dxi1
         ixtemp(2,i) = this%ix(2,pp) + dxi2
           
         xtemp(1,i) = xtemp(1,i) - dxi1
         xtemp(2,i) = xtemp(2,i) - dxi2

         pp = pp + 1
     end do      
       
     ! getting the new fields
     !!!!!!!!!!!!!!!!!!!!
     call get_emf( emf, bp, ep, ixtemp, xtemp, np, this%interpolation )
       
     ! Calculate the third set of Runge-Kutta coefficients    
     call rk_get(ep,bp,utemp,gam_temp,np,Cdamp,k_p_rung,k_x_rung,q_weight_temp) 
        
     ! Update the final momenta with the fourth correction, also calculating the new position xbuf 
     pp = ptrcur
     do i=1,np
        this%p(1,pp) = ufinal(1,i)+tem*k_p_rung(1,i)*0.3333333333333333_p_k_part            
        this%p(2,pp) = ufinal(2,i)+tem*k_p_rung(2,i)*0.3333333333333333_p_k_part
        this%p(3,pp) = ufinal(3,i)+tem*k_p_rung(3,i)*0.3333333333333333_p_k_part        
        
        pp = pp + 1
     end do
        
     do i=1,np
        xbuf(1,i) = xfinal(1,i)+tem_x*k_x_rung(1,i)*rdx(1)*0.3333333333333333_p_k_part
        xbuf(2,i) = xfinal(2,i)+tem_x*k_x_rung(2,i)*rdx(2)*0.3333333333333333_p_k_part
     end do         
     
     ! Get final gamma - this is required to deposit j3
     pp = ptrcur
     do i=1,np
        rg(i) = 1.0_p_k_part / &
          sqrt( 1.0_p_k_part + &
            this%p(1,pp)**2 + &
            this%p(2,pp)**2 + &
            this%p(3,pp)**2 )                
        
        pp = pp + 1
     end do
           
     ! Get cell crossings - these will also be used in current deposition
     do i=1,np
       dxi(1,i) = ntrim( xbuf(1,i) )
       dxi(2,i) = ntrim( xbuf(2,i) )
     end do
           
     ! Deposit current
     call this % dep_current_2d(  jay, dxi, xbuf, &
                         this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                         this%q(ptrcur:), rg,     &
                         this%p(:,ptrcur:),      &
                         np, gdt )

     ! trim positions and store results
     pp = ptrcur
     do i = 1, np
       this%x(1,pp)  = xbuf(1,i)  - dxi(1,i) 
       this%x(2,pp)  = xbuf(2,i)  - dxi(2,i)
       this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
       this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)
 
       pp = pp + 1
     end do

  enddo

contains

!---------------------------------------------------------------------------------------------------
subroutine rk_get(ep,bp,utemp,gam,np,Cdamp,k_p_rung,k_x_rung,q_weight_temp) 
!---------------------------------------------------------------------------------------------------
! Get the values of coefficients for Runge-Kutta integration
!---------------------------------------------------------------------------------------------------
  implicit none
  
  ! dummy variables
  real(p_k_part), dimension(:,:), intent(in) :: ep, bp, utemp
  real(p_k_part), dimension(:), intent(in) :: gam           !this is 1/gamma
  integer, intent(in) :: np
  real(p_k_part), intent(in) :: Cdamp
  real(p_k_part), dimension(:,:), intent(inout) :: k_p_rung
  real(p_k_part), dimension(:,:), intent(inout) :: k_x_rung
  real(p_k_part), dimension(:) :: q_weight_temp

  !local variables
  integer :: i
  real(p_k_part) :: Fx, Fy, Fz, pE, gam_temp

  do i=1,np
    Fx = ep(1,i)+gam(i)*(utemp(2,i)*bp(3,i)-utemp(3,i)*bp(2,i))
    Fy = ep(2,i)+gam(i)*(utemp(3,i)*bp(1,i)-utemp(1,i)*bp(3,i))
    Fz = ep(3,i)+gam(i)*(utemp(1,i)*bp(2,i)-utemp(2,i)*bp(1,i))
    pE= utemp(1,i)*ep(1,i)+utemp(2,i)*ep(2,i)+utemp(3,i)*ep(3,i)
    gam_temp=1.0_p_k_part/gam(i)

    !gam(i) is 1/gamma       and gam_temp is gamma :) 

    k_p_rung(1,i) = Fx+Cdamp*(ep(2,i)*bp(3,i)-ep(3,i)*bp(2,i)+gam(i)*( &
               bp(2,i)*bp(1,i)*utemp(2,i) &
                -(bp(2,i)**2+bp(3,i)**2)*utemp(1,i)+bp(3,i)*bp(1,i)*utemp(3,i) &
                +ep(1,i)*pE+utemp(1,i)*(pE**2))-utemp(1,i)*(Fx*Fx+Fy*Fy+Fz*Fz)*gam_temp)
    
    k_p_rung(2,i) = Fy + Cdamp*(ep(3,i)*bp(1,i)-ep(1,i)*bp(3,i)+gam(i)*( &
                bp(3,i)*bp(2,i)*utemp(3,i)-(bp(3,i)**2+bp(1,i)**2)*utemp(2,i) &
                +bp(1,i)*bp(2,i)*utemp(1,i)+ep(2,i)*pE+utemp(2,i)*(pE**2))- &
                gam_temp*utemp(2,i)*(Fx*Fx+Fy*Fy+Fz*Fz))
    
    k_p_rung(3,i) = Fz + Cdamp*(ep(1,i)*bp(2,i)-ep(2,i)*bp(1,i)+gam(i)*( &
                bp(1,i)*bp(3,i)*utemp(1,i)-(bp(1,i)**2+bp(2,i)**2)*utemp(3,i) &
                +bp(2,i)*bp(3,i)*utemp(2,i)+ep(3,i)*pE+utemp(3,i)*(pE**2))-&
                gam_temp*utemp(3,i)*(Fx*Fx+Fy*Fy+Fz*Fz))
          
    k_x_rung(1,i) = gam(i)*utemp(1,i)
    k_x_rung(2,i) = gam(i)*utemp(2,i)
  end do
       
!       do i=1,np
!         k_p_rung(1,i) = ep(1,i)+gam(i)*(utemp(2,i)*bp(3,i)-utemp(3,i)*bp(2,i))
!         k_p_rung(2,i) = ep(2,i)+gam(i)*(utemp(3,i)*bp(1,i)-utemp(1,i)*bp(3,i))
!         k_p_rung(3,i) = ep(3,i)+gam(i)*(utemp(1,i)*bp(2,i)-utemp(2,i)*bp(1,i))
!         k_x_rung(1,i) = gam(i)*utemp(1,i)
!         k_x_rung(2,i) = gam(i)*utemp(2,i)
!       end do
    
end subroutine rk_get


end subroutine advance_deposit_radcool_2d
!---------------------------------------------------------------------------------------------------



subroutine rk_get_3d(ep,bp,utemp,gam,np,Cdamp,k_p_rung,k_x_rung,q_weight_temp) 
!---------------------------------------------------------------------------------------------------
! Get the values of coefficients for Runge-Kutta integration 3d
!---------------------------------------------------------------------------------------------------
  implicit none
  
  ! dummy variables
  real(p_k_part), dimension(:,:), intent(in) :: ep, bp, utemp
  real(p_k_part), dimension(:), intent(in) :: gam           !this is 1/gamma
  integer, intent(in) :: np
  real(p_k_part), intent(in) :: Cdamp
  real(p_k_part), dimension(:,:), intent(inout) :: k_p_rung
  real(p_k_part), dimension(:,:), intent(inout) :: k_x_rung
  real(p_k_part), dimension(:), intent(in) :: q_weight_temp
  !local variables
  integer :: i
  real(p_k_part) :: Fx, Fy, Fz, pE, gam_temp

  do i=1,np
    Fx = ep(1,i) +gam(i) * (utemp(2,i)*bp(3,i) - utemp(3,i)*bp(2,i))
    Fy = ep(2,i) +gam(i) * (utemp(3,i)*bp(1,i) - utemp(1,i)*bp(3,i))
    Fz = ep(3,i) +gam(i) * (utemp(1,i)*bp(2,i) - utemp(2,i)*bp(1,i))
    pE= utemp(1,i) * ep(1,i) + utemp(2,i) * ep(2,i) + utemp(3,i) * ep(3,i)
    gam_temp=1.0_p_k_part/gam(i)

    !gam(i) is 1/gamma       and gam_temp is gamma :) 

    k_p_rung(1,i) = Fx+Cdamp*(ep(2,i)*bp(3,i)-ep(3,i)*bp(2,i)+gam(i)*( &
               bp(2,i)*bp(1,i)*utemp(2,i) &
                -(bp(2,i)**2+bp(3,i)**2)*utemp(1,i)+bp(3,i)*bp(1,i)*utemp(3,i) &
                +ep(1,i)*pE+utemp(1,i)*(pE**2))-utemp(1,i)*(Fx*Fx+Fy*Fy+Fz*Fz)*gam_temp)
    
    k_p_rung(2,i) = Fy + Cdamp*(ep(3,i)*bp(1,i)-ep(1,i)*bp(3,i)+gam(i)*( &
                bp(3,i)*bp(2,i)*utemp(3,i)-(bp(3,i)**2+bp(1,i)**2)*utemp(2,i) &
                +bp(1,i)*bp(2,i)*utemp(1,i)+ep(2,i)*pE+utemp(2,i)*(pE**2))- &
                gam_temp*utemp(2,i)*(Fx*Fx+Fy*Fy+Fz*Fz))
    
    k_p_rung(3,i) = Fz + Cdamp*(ep(1,i)*bp(2,i)-ep(2,i)*bp(1,i)+gam(i)*( &
                bp(1,i)*bp(3,i)*utemp(1,i)-(bp(1,i)**2+bp(2,i)**2)*utemp(3,i) &
                +bp(2,i)*bp(3,i)*utemp(2,i)+ep(3,i)*pE+utemp(3,i)*(pE**2))-&
                gam_temp*utemp(3,i)*(Fx*Fx+Fy*Fy+Fz*Fz))
          
    k_x_rung(1,i) = gam(i) * utemp(1,i)
    k_x_rung(2,i) = gam(i) * utemp(2,i)
    k_x_rung(3,i) = gam(i) * utemp(3,i)
    
  end do
       
!       do i=1,np
!         k_p_rung(1,i) = ep(1,i) +gam(i) * (utemp(2,i)*bp(3,i) - utemp(3,i)*bp(2,i))
!         k_p_rung(2,i) = ep(2,i) +gam(i) * (utemp(3,i)*bp(1,i) - utemp(1,i)*bp(3,i))
!         k_p_rung(3,i) = ep(3,i) +gam(i) * (utemp(1,i)*bp(2,i) - utemp(2,i)*bp(1,i))
!         k_x_rung(1,i) = gam(i) * utemp(1,i)
!         k_x_rung(2,i) = gam(i) * utemp(2,i)
!         k_x_rung(3,i) = gam(i) * utemp(3,i)
!       end do
    
end subroutine rk_get_3d
!---------------------------------------------------------------------------------------------------


subroutine push_radcool_3d( this, emf, dt , omega_p0, ptrcur, np, xbuf)

  implicit none
  
  integer, parameter :: rank = 3
  
  ! dummy variables
  class( t_species ),    intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  real(p_double), intent(in) :: dt, omega_p0
  integer, intent(in) :: ptrcur, np
  real(p_k_part), dimension(:,:), intent(inout) :: xbuf 

  ! local variables
  real(p_k_part) :: tem, tem_x

  integer :: i, pp

  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp, ufinal
  real(p_k_part), dimension(rank,p_cache_size) :: xtemp, xfinal  
  real(p_k_part), dimension(p_cache_size) :: q_weight_temp
  integer, dimension(rank,p_cache_size) :: ixtemp, ixfinal
  real(p_k_part), dimension(rank) :: rdx 
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: k_p_rung          ! Runge-Kutta coefficients, 3p
  real(p_k_part), dimension(rank,p_cache_size) :: k_x_rung            ! 3x
  real(p_k_part), dimension(p_cache_size) :: gam_temp

  real(p_k_part) :: Cdamp
  
  integer :: dxi1, dxi2, dxi3
  
  ! executable statements

  rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
  rdx(2) = real( 1.0_p_double/this%dx(2), p_k_part )
  rdx(3) = real( 1.0_p_double/this%dx(3), p_k_part )

     
  ! get factor that includes timestep and charge-to-mass ratio
  ! note that the charge-to-mass ratio is used since the
  ! momentum p is not the total momentum of a particles but
  ! the momentum per unit restmass (which is the electron mass)

  tem = real( 0.5_p_double * dt / this%rqm, p_k_part )
  tem_x = real(0.5_p_double * dt, p_k_part)

  Cdamp = real(2.0_p_k_part * cgs_e * cgs_e * omega_p0 / &
                 (3.0_p_k_part * cgs_me * cgs_c**3 ) )

  Cdamp = Cdamp * this%q_real**2 / this%rqm 


  ! Get fields at particle position
  call this % get_emf( emf, bp, ep, np, ptrcur )
   
  !---------------------------------------------------------------
  !--------------------First Runge-Kutta step---------------------
  
  ! Take the initial momenta         
  pp = ptrcur
  do i=1,np
    utemp(1,i) = this%p(1,pp) 
    utemp(2,i) = this%p(2,pp) 
    utemp(3,i) = this%p(3,pp) 

    gam_temp(i)= 1.0_p_k_part / sqrt(1.0_p_k_part+utemp(1,i)**2+ &
                                 utemp(2,i)**2+ &
                                 utemp(3,i)**2)

    pp = pp + 1
  end do
  
  ! Take initial coordinates
  pp = ptrcur
  do i=1,np
     ! This is overwritten below
     !xtemp(1,i) = this%x(1,pp)
     !xtemp(2,i) = this%x(2,pp)
     !xtemp(3,i) = this%x(3,pp)

     ! This is overwritten below
     !ixtemp(1,i) = this%ix(1,pp)
     !ixtemp(2,i) = this%ix(2,pp)
     !ixtemp(3,i) = this%ix(3,pp)

     q_weight_temp(i) = this%q(pp)

     pp = pp + 1
  end do
         
  ! Calculate the first set of Runge-Kutta coefficients
  call rk_get_3d(ep,bp,utemp,gam_temp,np,Cdamp,k_p_rung,k_x_rung,q_weight_temp)  
            
            
  ! Update the final values with first correction
  pp = ptrcur
  do i=1,np
     ufinal(1,i) = this%p(1,pp) + tem * k_p_rung(1,i) * 0.3333333333333333_p_k_part               
     ufinal(2,i) = this%p(2,pp) + tem * k_p_rung(2,i) * 0.3333333333333333_p_k_part 
     ufinal(3,i) = this%p(3,pp) + tem * k_p_rung(3,i) * 0.3333333333333333_p_k_part 
    
     pp = pp + 1
  end do
            
            
  ! Final value for coordinates update. 
  ! Trimmed only at the end of the whole timestep
  pp = ptrcur
  do i=1,np
     xfinal(1,i) = this%x(1,pp) + tem_x * k_x_rung(1,i) * rdx(1) * 0.3333333333333333_p_k_part 
     xfinal(2,i) = this%x(2,pp) + tem_x * k_x_rung(2,i) * rdx(2) * 0.3333333333333333_p_k_part 
     xfinal(3,i) = this%x(3,pp) + tem_x * k_x_rung(3,i) * rdx(3) * 0.3333333333333333_p_k_part 
  
     ixfinal(1,i) = this%ix(1,pp)
     ixfinal(2,i) = this%ix(2,pp)
     ixfinal(3,i) = this%ix(3,pp)
    
     pp=pp+1
  end do         
           
  !----------------------------------------------------------------
  !--------------------Second Runge-Kutta step---------------------
 
  ! Update the temporary values of momenta for the next step
  pp = ptrcur
  do i=1,np
     utemp(1,i) = this%p(1,pp) + tem * k_p_rung(1,i)              
     utemp(2,i) = this%p(2,pp) + tem * k_p_rung(2,i)
     utemp(3,i) = this%p(3,pp) + tem * k_p_rung(3,i)
    
     gam_temp(i)= 1.0_p_k_part / sqrt(1.0_p_k_part+utemp(1,i)**2+ &
                                 utemp(2,i)**2+ &
                                 utemp(3,i)**2)
     pp = pp + 1
  end do
           
  ! Updating the coordinates on the way to the second step
  ! Important to cover all the cases because of the field interpolation 
      
  pp = ptrcur
  do i=1,np
     xtemp(1,i) = this%x(1,pp) + tem_x * k_x_rung(1,i) * rdx(1)
     xtemp(2,i) = this%x(2,pp) + tem_x * k_x_rung(2,i) * rdx(2)
     xtemp(3,i) = this%x(3,pp) + tem_x * k_x_rung(3,i) * rdx(3)
     
     pp = pp + 1
  end do
  
  pp = ptrcur
  do i=1,np
     dxi1 = ntrim( xtemp(1,i) )    
     dxi2 = ntrim( xtemp(2,i) )
     dxi3 = ntrim( xtemp(3,i) )
       
     ixtemp(1,i) = this%ix(1,pp) + dxi1
     ixtemp(2,i) = this%ix(2,pp) + dxi2
     ixtemp(3,i) = this%ix(3,pp) + dxi3  
       
     xtemp(1,i) = xtemp(1,i) - dxi1
     xtemp(2,i) = xtemp(2,i) - dxi2
     xtemp(3,i) = xtemp(3,i) - dxi3

     pp = pp + 1
  end do      
               

  ! getting the new fields
  !!!!!!!!!!!!!!!!!!!!
  call get_emf( emf, bp, ep, ixtemp, xtemp, np, this%interpolation )
        
  ! Calculate the second set of Runge-Kutta coefficients    
  call rk_get_3d(ep,bp,utemp,gam_temp,np,Cdamp,k_p_rung,k_x_rung,q_weight_temp)  
       
  ! Update the final momenta with second correction
  do i=1,np
     ufinal(1,i) = ufinal(1,i) + tem * k_p_rung(1,i) * 0.6666666666666666_p_k_part             
     ufinal(2,i) = ufinal(2,i) + tem * k_p_rung(2,i) * 0.6666666666666666_p_k_part
     ufinal(3,i) = ufinal(3,i) + tem * k_p_rung(3,i) * 0.6666666666666666_p_k_part
  end do
       
  ! Update the coordinates with the corrections
  ! Again, we do not trim the coordinates

  do i=1,np
     xfinal(1,i) = xfinal(1,i) + tem_x * k_x_rung(1,i) * rdx(1) * 0.6666666666666666_p_k_part
     xfinal(2,i) = xfinal(2,i) + tem_x * k_x_rung(2,i) * rdx(2) * 0.6666666666666666_p_k_part
     xfinal(3,i) = xfinal(3,i) + tem_x * k_x_rung(3,i) * rdx(3) * 0.6666666666666666_p_k_part
  end do         

  !----------------------------------------------------------------
  !--------------------Third Runge-Kutta step----------------------
   
  ! Update the temporary values of x, gamma and p for the next step
  
 
  pp = ptrcur
  do i=1,np
     utemp(1,i) = this%p(1,pp)+tem*k_p_rung(1,i)              
     utemp(2,i) = this%p(2,pp)+tem*k_p_rung(2,i)
     utemp(3,i) = this%p(3,pp)+tem*k_p_rung(3,i)
     
     gam_temp(i)= 1.0_p_k_part / sqrt(1.0_p_k_part+utemp(1,i)**2+ &
                                 utemp(2,i)**2+ &
                                 utemp(3,i)**2)
     pp = pp + 1
  end do
          
  pp = ptrcur
  do i=1,np
     xtemp(1,i) = this%x(1,pp) + tem_x * k_x_rung(1,i) * rdx(1)
     xtemp(2,i) = this%x(2,pp) + tem_x * k_x_rung(2,i) * rdx(2)
     xtemp(3,i) = this%x(3,pp) + tem_x * k_x_rung(3,i) * rdx(3)

     dxi1 = ntrim( xtemp(1,i) )    
     dxi2 = ntrim( xtemp(2,i) )
     dxi3 = ntrim( xtemp(3,i) )
     
     ixtemp(1,i) = this%ix(1,pp) + dxi1
     ixtemp(2,i) = this%ix(2,pp) + dxi2
     ixtemp(3,i) = this%ix(3,pp) + dxi3
       
     xtemp(1,i) = xtemp(1,i) - dxi1
     xtemp(2,i) = xtemp(2,i) - dxi2
     xtemp(3,i) = xtemp(3,i) - dxi3

     pp = pp + 1
  end do      
               
  ! getting the new fields
  !!!!!!!!!!!!!!!!!!!!
  call get_emf( emf, bp, ep, ixtemp, xtemp, np, this%interpolation )
       
  ! Calculate the third set of Runge-Kutta coefficients    
  call rk_get_3d(ep,bp,utemp,gam_temp,np,Cdamp,k_p_rung,k_x_rung,q_weight_temp) 
        
  ! Update the final values with the third correction
  ! Still no trimming 
  do i=1,np
     ufinal(1,i) = ufinal(1,i) + tem * k_p_rung(1,i) * 0.6666666666666666_p_k_part              
     ufinal(2,i) = ufinal(2,i) + tem * k_p_rung(2,i) * 0.6666666666666666_p_k_part
     ufinal(3,i) = ufinal(3,i) + tem * k_p_rung(3,i) * 0.6666666666666666_p_k_part

     xfinal(1,i) = xfinal(1,i) + tem_x * k_x_rung(1,i) * rdx(1) * 0.6666666666666666_p_k_part
     xfinal(2,i) = xfinal(2,i) + tem_x * k_x_rung(2,i) * rdx(2) * 0.6666666666666666_p_k_part
     xfinal(3,i) = xfinal(3,i) + tem_x * k_x_rung(3,i) * rdx(3) * 0.6666666666666666_p_k_part
  end do         

  !----------------------------------------------------------------
  !--------------------Fourth Runge-Kutta step---------------------
  
  pp = ptrcur
  do i=1,np
     utemp(1,i) = this%p(1,pp) + 2.0_p_k_part * tem * k_p_rung(1,i)              
     utemp(2,i) = this%p(2,pp) + 2.0_p_k_part * tem * k_p_rung(2,i)
     utemp(3,i) = this%p(3,pp) + 2.0_p_k_part * tem * k_p_rung(3,i)
     
     gam_temp(i)= 1.0_p_k_part / sqrt(1.0_p_k_part+utemp(1,i)**2+ &
                                 utemp(2,i)**2+ &
                                 utemp(3,i)**2)
     pp = pp + 1
  end do
           
  pp = ptrcur
  do i=1,np
     xtemp(1,i) = this%x(1,pp) + 2.0_p_k_part * tem_x * k_x_rung(1,i) * rdx(1)
     xtemp(2,i) = this%x(2,pp) + 2.0_p_k_part * tem_x * k_x_rung(2,i) * rdx(2)
     xtemp(3,i) = this%x(3,pp) + 2.0_p_k_part * tem_x * k_x_rung(3,i) * rdx(3)

     dxi1 = ntrim( xtemp(1,i) )    
     dxi2 = ntrim( xtemp(2,i) )
     dxi3 = ntrim( xtemp(3,i) )
       
     ixtemp(1,i) = this%ix(1,pp) + dxi1
     ixtemp(2,i) = this%ix(2,pp) + dxi2
     ixtemp(3,i) = this%ix(3,pp) + dxi3
       
     xtemp(1,i) = xtemp(1,i) - dxi1
     xtemp(2,i) = xtemp(2,i) - dxi2
     xtemp(3,i) = xtemp(3,i) - dxi3

     pp = pp + 1
  end do      
         
       
  ! getting the new fields
  !!!!!!!!!!!!!!!!!!!!
  call get_emf( emf, bp, ep, ixtemp, xtemp, np, this%interpolation )
            
  ! Calculate the third set of Runge-Kutta coefficients    
  call rk_get_3d(ep,bp,utemp,gam_temp,np,Cdamp,k_p_rung,k_x_rung,q_weight_temp) 
        
  ! Update the final momenta with the fourth correction, also calculating the new position xbuf 
  pp = ptrcur
  do i=1,np
     this%p(1,pp) = ufinal(1,i) + tem * k_p_rung(1,i) * 0.3333333333333333_p_k_part            
     this%p(2,pp) = ufinal(2,i) + tem * k_p_rung(2,i) * 0.3333333333333333_p_k_part
     this%p(3,pp) = ufinal(3,i) + tem * k_p_rung(3,i) * 0.3333333333333333_p_k_part        
     
     xbuf(1,i) = xfinal(1,i) + tem_x * k_x_rung(1,i) * rdx(1) * 0.3333333333333333_p_k_part        
     xbuf(2,i) = xfinal(2,i) + tem_x * k_x_rung(2,i) * rdx(2) * 0.3333333333333333_p_k_part    
     xbuf(3,i) = xfinal(3,i) + tem_x * k_x_rung(3,i) * rdx(3) * 0.3333333333333333_p_k_part
     
     pp = pp + 1
  end do     
   
end subroutine push_radcool_3d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine advance_deposit_radcool_3d( this, emf, jay, gdt, omega_p0, i0, i1 )
!---------------------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 3

  ! dummy variables
  class( t_species ),    intent(inout) :: this
  class( t_emf ), intent( in )  ::  emf
  type( t_vdf), intent(inout) :: jay
  real(p_double), intent(in) :: gdt, omega_p0
  integer, intent(in) :: i0, i1

  ! local variables
  integer :: np, ptrcur, pp, i


  real(p_k_part), dimension(rank,p_cache_size) :: xbuf
  integer, dimension(rank,p_cache_size) :: dxi


  do ptrcur = i0, i1, p_cache_size

     ! check if last copy of table and set np
     if( ptrcur + p_cache_size > i1 ) then
         np = i1 - ptrcur + 1
     else
         np = p_cache_size
     endif
     
     ! advance momenta and get final positions indexed to original cell
       call push_radcool_3d(this, emf, gdt, omega_p0, ptrcur, np, xbuf)
     
     ! get cell crossings
     do i=1,np
       dxi(1,i) = ntrim( xbuf(1,i) )
       dxi(2,i) = ntrim( xbuf(2,i) )
       dxi(3,i) = ntrim( xbuf(3,i) )
     end do
  
     ! Deposit current
     call this % dep_current_3d( jay, dxi, xbuf, &
                       this%ix(:,ptrcur:), this%x(:,ptrcur:), &
                       this%q(ptrcur:), np, gdt )

     ! trim positions and store results
     pp = ptrcur
     do i = 1, np
       this%x(1,pp)  = xbuf(1,i)  - dxi(1,i) 
       this%x(2,pp)  = xbuf(2,i)  - dxi(2,i)
       this%x(3,pp)  = xbuf(3,i)  - dxi(3,i)
       this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
       this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)
       this%ix(3,pp) = this%ix(3,pp) + dxi(3,i)
 
       pp = pp + 1
     end do

  enddo
  

end subroutine advance_deposit_radcool_3d 
!---------------------------------------------------------------------------------------------------


end module m_species_radcool