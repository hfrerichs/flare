subroutine generate_flux_surface_2D
  use run_control, only: x_start, Trace_Step, Trace_Method, Output_File
  use flux_surface_2D
  implicit none

  type(t_flux_surface_2D) :: S


  call S%generate(x_start, Trace_Step, Trace_Method)
  call S%write(Output_File=Output_File)
  
end subroutine generate_flux_surface_2D
