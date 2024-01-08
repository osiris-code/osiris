!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     field diagnostics class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

!*******************************************************************************

!-------------------------------------------------------------------------------
subroutine report_energy_emf( this, no_co, tstep, t )
!-------------------------------------------------------------------------------

  use m_emf_define
  use m_time_step
  use m_node_conf
  use m_time_step
  use m_parameters

  use m_emf_diag

  use m_vdf_define
  use m_vdf_math

  use m_logprof

  implicit none

  class( t_emf ), intent(inout)    :: this
  class( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double),     intent(in) :: t

  real(p_double), dimension(2*p_f_dim) :: temp_int

  integer :: ierr

  ! file name and path
  character(len=256) :: full_name, path

  call begin_event( diag_emf_ev )

  ! reports on integrated field energy
  if (test_if_report( tstep, this%diag%ndump_fac_ene_int ) ) then

     call this%fill_data( tstep )

     ! get local e and b field integrals
     call total( this%b, temp_int, pow = 2 )
     call total( this%e, temp_int(p_f_dim+1:), pow = 2 )

     ! sum up results from all nodes
     call reduce_array( no_co, temp_int, operation = p_sum)

     ! save the data
     if ( root(no_co) ) then
        ! normalize energies
        temp_int = temp_int*this%e%dvol()*0.5_p_double

        ! setup path and file names
        path  =  trim(path_hist)
        full_name = trim(path) // 'fld_ene'

        ! Open file and position at the last record
        if ( t == 0.0_p_double ) then

           call mkdir( path, ierr )

           open (unit=file_id_fldene, file=full_name, status = 'REPLACE' , &
                 form='formatted')

           ! Write Header
           write(file_id_fldene, '(A)' ) '! EM field energy per field component'
           write(file_id_fldene,'( A6, 1X,A15,6(1X,A23) )') &
                  'Iter','Time    ','B1','B2','B3',&
                                    'E1','E2','E3'
        else

           open (unit=file_id_fldene, file=full_name, position = 'append', &
                 form='formatted')
        endif

        write(file_id_fldene, '( I6, 1X,g15.8, 6(1X,es23.16) )') n(tstep), t,temp_int
        close(file_id_fldene)

     endif

  endif

  call end_event( diag_emf_ev )

end subroutine report_energy_emf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       report on electro-magnetic field - diagnostic
!-----------------------------------------------------------------------------------------
subroutine report_diag_emf( this, g_space, grid, no_co, tstep, t, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

  use m_system
  use m_emf_define
  use m_space
  use m_grid_define
  use m_node_conf
  use m_time_step

  use m_emf_diag
  use m_logprof

  use m_vdf_define
  use m_vdf_report
  use m_vdf_math
  use m_vdf_comm, only : t_vdf_msg
  use m_emf_poynting
  use m_emf_psi

  ! (*debug*) - Uncomment to debug PML boundaries
  ! use m_vpml

  implicit none

  class( t_emf ),                 intent(inout) :: this

  type( t_space ),     intent(in) :: g_space
  class( t_grid ),      intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double),      intent(in) :: t
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  type( t_vdf_report ), pointer :: rep
  integer :: i

  ! temporary vdf object for diagnostics
  type( t_vdf ) :: vdf_a

  type( t_vdf ), pointer :: psi => null()

  ! executable statements
  call begin_event( diag_emf_ev )

  rep => this%diag%reports
  do
    if ( .not. associated( rep ) ) exit

    if ( if_report( rep, tstep ) ) then
       call this%fill_data( tstep )
       
       select case ( rep%quant )
         case ( p_e1, p_e2, p_e3 )
           call report_vdf( rep, this%e, rep%quant - p_e1 + 1, g_space, grid, no_co, tstep, t )

         case ( p_b1, p_b2, p_b3 )
           call report_vdf( rep, this%b, rep%quant - p_b1 + 1, g_space, grid, no_co, tstep, t )

         case ( p_ext_e1, p_ext_e2, p_ext_e3 )
           call report_vdf( rep, this%ext_e, rep%quant - p_ext_e1 + 1, g_space, grid, no_co, tstep, t )

         case ( p_ext_b1, p_ext_b2, p_ext_b3 )
           call report_vdf( rep, this%ext_b, rep%quant - p_ext_b1 + 1, g_space, grid, no_co, tstep, t )

         case ( p_part_e1, p_part_e2, p_part_e3 )
           call report_vdf( rep, this%e_part, rep%quant - p_part_e1 + 1, g_space, grid, no_co, tstep, t )

         case ( p_part_b1, p_part_b2, p_part_b3 )
           call report_vdf( rep, this%b_part, rep%quant - p_part_b1 + 1, g_space, grid, no_co, tstep, t )

         case ( p_ene_e1, p_ene_e2, p_ene_e3 )
           ! calculate energy in field component
           call vdf_a%new( this%e, f_dim = 1, copy = .true., fc = rep%quant - p_ene_e1 + 1 )
           call pow(vdf_a, 2)

           ! report it

           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )

           ! free temporary memory
           call vdf_a%cleanup()

         case ( p_ene_b1, p_ene_b2, p_ene_b3 )
           ! calculate energy in field component
           call vdf_a % new( this%b, f_dim = 1, copy = .true., fc = rep%quant - p_ene_b1 + 1 )
           call pow(vdf_a, 2)

           ! report it

           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )

           ! free temporary memory
           call vdf_a%cleanup()

         case ( p_ene_e )
           ! calculate energy in electric field
           call vdf_a % new( this%e, f_dim = 1, copy = .true., fc = 1 )
           call pow( vdf_a, 2 )
           do i = 2, 3
             call add_pow( vdf_a, this%e, 2, fc = i )
           enddo

           ! report it

           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )

           ! free temporary memory
           call vdf_a % cleanup()

         case ( p_ene_b )
           ! calculate energy in magnetic field
           call vdf_a % new( this%b, f_dim = 1, copy = .true., fc = 1 )
           call pow( vdf_a, 2 )
           do i = 2, 3
             call add_pow( vdf_a, this%b, 2, fc = i )
           enddo

           ! report it

           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )

           ! free temporary memory
           call vdf_a % cleanup()

         case ( p_ene_emf )
           ! calculate energy in electric field
           call vdf_a % new( this%e, f_dim = 1, copy = .true., fc = 1 )
           call pow( vdf_a, 2 )
           do i = 2, 3
             call add_pow( vdf_a, this%e, 2, fc = i )
           enddo
           do i = 1, 3
             call add_pow( vdf_a, this%b, 2, fc = i )
           enddo

           ! report it

           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )

           ! free temporary memory
           call vdf_a % cleanup()

         case ( p_div_e )
           ! calculate electric field divergence
           call vdf_a % new( this%e, f_dim = 1 )
           call div( this%e, vdf_a )
           ! report it

           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )
           ! free temporary memory
           call vdf_a % cleanup()

         case ( p_div_b )
           ! calculate magnetic field divergence
           call vdf_a % new( this%b, f_dim = 1 )
           call div( this%b, vdf_a, 1 )
           ! report it

           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )
           ! free temporary memory
           call vdf_a%cleanup()

         case ( p_charge_cons )
           ! print *, 'Reporting charge conservation '
            ! call report_vdf( rep, this%f, 1, g_space, grid, no_co, tstep, t )
          ERROR("Charge conservation diagnostic temporarily disabled")
          call abort_program()

         case ( p_psi )
           ! calculate psi diagnostic
           call get_psi( this, n(tstep), no_co, psi, send_msg, recv_msg )
           ! report it

           call report_vdf( rep, psi, 1, g_space, grid, no_co, tstep, t )

           ! psi is just a pointer to the real psi buffer, no need to deallocate the memory
           ! just nullify the pointer
           psi => null()

         case ( p_s1, p_s2, p_s3 )
           ! calculate component of Poynting flux
           call vdf_a % new( this%e, f_dim = 1 )
           call poynting( this, rep%quant - p_s1 + 1, vdf_a )

           ! report it
           call report_vdf( rep, vdf_a, 1, g_space, grid, no_co, tstep, t )

           ! free temporary memory
           call vdf_a % cleanup()

         case default
           ! unknown quantity, must belong to a subclass
              ! print *, 'In report_diag_emf, unknown quantity : ',  trim(rep%name)
           continue

       end select

    endif

    rep => rep%next
  enddo

  ! (*debug*) - Report PML boundary data. Only works in serial runs
  ! call report( this%bnd_con%vpml_all, g_space, grid, no_co, tstep, t, this%e%dx )

  call end_event( diag_emf_ev )

end subroutine report_diag_emf
!-----------------------------------------------------------------------------------------
