!===============================================================================
! Generete axisymmetric (unperturbed) flux surfaces
!
! Input (taken from run control file):
! x_start(1:2)          R and Z coordinates [cm] of reference position on flux surface
!
! Trace_Step, Trace_Method -> see e.g. connection length
! N_steps               Max. number of trace steps for flux surface generation (optional)
! Output_File
!===============================================================================
subroutine generate_flux_surface_2D
  use run_control, only: x_start, Trace_Step, N_steps, Trace_Method, Output_File
  use flux_surface_2D
  use parallel
  implicit none

  type(t_flux_surface_2D) :: S


  if (firstP) then
     write (6, *) 'Generate flux surface contour (2D), output in: ', trim(Output_File)
  else
     return
  endif


  if (N_steps > 0) then
     call S%generate(x_start(1:2), Trace_Step=Trace_Step, N_steps=N_steps, Trace_Method=Trace_Method)
  else
     call S%generate(x_start(1:2), Trace_Step=Trace_Step, Trace_Method=Trace_Method)
  endif
  call S%plot(filename=Output_File)
  
end subroutine generate_flux_surface_2D
