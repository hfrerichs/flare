!===============================================================================
! Generete axisymmetric (unperturbed) flux surfaces
!
! Input (taken from run control file):
! Input_Format          = 1: cylindrical coordinates (R[cm], Z[cm], Phi[deg])
!                       = 2: flux coordinates (Theta[deg], PsiN, Phi[deg])
! x_start(1:2)          Reference point on flux surface
!
!
! Trace_Step, Trace_Method -> see e.g. connection length
! N_steps               Max. number of trace steps for flux surface generation (optional)
! Output_File
!===============================================================================
subroutine generate_flux_surface_2D
  use run_control, only: Input_Format, x_start, Trace_Step, N_steps, Trace_Method, Output_File
  use flux_surface_2D
  use equilibrium
  use parallel
  implicit none

  type(t_flux_surface_2D) :: S
  real(real64)            :: r(3), y(3)
  integer                 :: ierr


  if (firstP) then
     write (6, *) 'Generate flux surface contour (2D), output in: ', trim(Output_File)
  else
     return
  endif


  ! set reference point on flux surface
  r = x_start
  if (Input_Format == 2) then
     y = x_start
     r = get_cylindrical_coordinates(y, ierr)
     if (ierr > 0) then
        write (6, *) 'error: cannot get cylindrical coordinates for flux coordinates ', x_start
        stop
     endif
  endif


  if (N_steps > 0) then
     call S%generate(r(1:2), Trace_Step=Trace_Step, N_steps=N_steps, Trace_Method=Trace_Method)
  else
     call S%generate(r(1:2), Trace_Step=Trace_Step, Trace_Method=Trace_Method)
  endif
  call S%plot(filename=Output_File)
  
end subroutine generate_flux_surface_2D
