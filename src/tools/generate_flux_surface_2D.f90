!===============================================================================
! Generete axisymmetric (unperturbed) flux surfaces
!
! Input (taken from run control file):
! x_start(1:2)          R and Z coordinates [cm] of reference position on flux surface
!
! Trace_Step, Trace_Method -> see e.g. connection length
! Output_File
!===============================================================================
subroutine generate_flux_surface_2D
  use run_control, only: x_start, Trace_Step, Trace_Method, Output_File
  use flux_surface_2D
  implicit none

  type(t_flux_surface_2D) :: S


  call S%generate(x_start(1:2), Trace_Step, Trace_Method)
  call S%write(Output_File=Output_File)
  
end subroutine generate_flux_surface_2D
