!===============================================================================
! Generete flux surfaces from field line tracing
!
! Input (taken from run control file):
! x_start               Reference point (R[cm], Z[cm], phi[deg]) on flux surface
! N_points
! ...
!
! Trace_Method -> see e.g. connection length
! Output_File
!===============================================================================
subroutine generate_flux_surface_3D
  use iso_fortran_env
  use run_control, only: x_start, N_points, N_sym, N_mult, N_steps, Trace_Method, Output_File
  use flux_surface_3D
  use math
  implicit none

  type(t_flux_surface_3D) :: S
  real(real64)            :: y0(3)


  y0    = x_start
  y0(3) = y0(3) / 180.d0 * pi
  if (N_steps == 0) N_steps = 10
  call S%generate(y0, N_points, N_sym, N_mult, N_steps, Trace_Method)
  call S%plot(filename=Output_File)
  
end subroutine generate_flux_surface_3D
