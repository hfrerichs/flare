!===============================================================================
! Calculate profile of the safety factor q (inverse rotational transform)
!
! Input (taken from run control file):
!    R_start, R_end	Defines radial domain at the midplane (Z=0)
!    Psi(2)             Alternative definition of radial domain by normalized poloidal flux
!    N_steps            Number of discretization steps used to calculate q
!    Phi_output		Toroidal coordinate [deg] from which field lines are traced
!    N_mult             Number of pol. turns used for averaging
!
!    Trace_Step, Trace_Method, Trace_Coords -> see e.g. connection_length
!    Output_File
!===============================================================================
subroutine safety_factor
  use run_control, only: R_start, R_end, Psi, N_steps, Phi_output, N_mult, &
                         Trace_Step, Trace_Method, Trace_Coords, &
                         Output_File
  use parallel
  use equilibrium
  use fieldline
  implicit none

  integer, parameter :: iu = 99

  type(t_fieldline)  :: F

  real*8, dimension(:,:), allocatable :: Qdata
  real*8  :: dR, r(3), y(3), X(3), Psi_av, Lc
  integer :: i, ierr, n_av


  ! initialize
  if (firstP) then
     write (6, *) 'Calculate safety factor, output in: ', adjustl(trim(Output_File))
     write (6, 1001) N_steps, R_start, R_end
     write (6, 1002) N_mult
     write (6, *)
  endif
  allocate (Qdata(0:N_steps,4))
  Qdata = 0.d0
  if (Psi(2) > 0.d0) then
     y(1)    = 180.d0
     y(2)    = Psi(2)
     y(3)    = Phi_output / pi * 180.d0
     r       = get_cylindrical_coordinates(y, ierr)
     if (ierr > 0) then
        write (6, *) 'error: cannot find real space coordinates for Psi = ', Psi(2)
        stop
     endif
     R_end   = r(1)

     y(1)    = 180.d0
     y(2)    = Psi(1)
     y(3)    = Phi_output / pi * 180.d0
     r       = get_cylindrical_coordinates(y, ierr)
     if (ierr > 0) then
        write (6, *) 'error: cannot find real space coordinates for Psi = ', Psi(1)
        stop
     endif
     R_start = r(1)
  else
     r(2)  = 0.d0
  endif
  dR    = (R_end - R_start) / N_steps
  r(3)  = Phi_output


  ! run main loop
  main_loop: do i=mype,N_steps,nprs
     r(1)       = R_start + i*dR
     Qdata(i,1) = r(1)
     Psi_av     = 0.d0
     n_av       = 0
     Lc         = 0.d0
     write (6, 2000) i, r(1)
     call coord_trans (r, CYLINDRICAL, y, Trace_Coords)
     call F%init(y, Trace_Step, Trace_Method, Trace_Coords)

     trace_loop: do
        Lc = Lc + F%trace_1step()
        if (F%intersect_boundary(X)) exit trace_loop

        Psi_av = Psi_av + F%get_PsiN()
        n_av   = n_av   + 1

        ! stop after N_mult poloidal turns
        if (abs(F%theta_int) > N_mult*pi2) then
           Qdata(i,2) = F%phi_int / F%theta_int
           Qdata(i,3) = Psi_av / n_av
           Qdata(i,4) = Lc / (F%theta_int / pi2)
           exit trace_loop
        endif
     enddo trace_loop
  enddo main_loop


  ! write data
  call wait_pe()
  call sum_real_data (Qdata, (N_steps+1)*4)
  if (firstP) then
     open  (iu, file=Output_File)
     write (iu, '(a)') '# R [cm], q, Psi, Lc'
     do i=0,N_steps
        write (iu, 2001) Qdata(i,:)
     enddo
     close (iu)
  endif

  ! cleanup
  deallocate (Qdata)

 1001 format (3x,'- Using ', i4, ' steps for radial domain: ', f8.3, ' -> ', f8.3, ' cm')
 1002 format (3x,'- Averaging over ', i4, ' toroidal turns')
 2000 format (5x,i5,f10.3)
 2001 format (4e18.10)
end subroutine safety_factor
