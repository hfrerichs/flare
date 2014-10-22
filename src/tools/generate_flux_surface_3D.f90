!===============================================================================
! Generete flux surfaces from field line tracing
!
! Input (taken from run control file):
! x_start               Reference point (R[cm], Z[cm], phi[deg]) on flux surface
! N_points              Number of points for each slice
! N_sym                 Toroidal symmetry number (i.e. 5: 0->72 deg)
! N_mult                Number of slices
! Output_Format         = 1: use poloidal angle to sample points of flux surface slice
!                       = 2: use distance along flux surface slice to sample points
!
! Trace_Method -> see e.g. connection length
! Output_File
!===============================================================================
subroutine generate_flux_surface_3D
  use iso_fortran_env
  use run_control, only: x_start, N_points, N_sym, N_mult, N_steps, Trace_Method, Output_File, Output_Format
  use flux_surface_3D
  use math
  use parallel
  implicit none

  type(t_flux_surface_3D) :: S
  real(real64)            :: y0(3)


  if (firstP) then
     write (6, *) 'Generate flux surface (from field line tracing), output in: ', adjustl(trim(Output_File))
     write (6, 1001) N_sym
     write (6, 1002) N_mult
     write (6, *)
  else
     return
  endif


  y0    = x_start
  y0(3) = y0(3) / 180.d0 * pi
  if (N_steps == 0) N_steps = 10
  call S%generate(y0, N_points, N_sym, N_mult, N_steps, Trace_Method, poloidal_coordinate=Output_Format)
  call S%plot(filename=Output_File)
  
 1001 format (8x,'Toroidal symmetry number:         ',i4)
 1002 format (8x,'Number of slices to be generated: ',i4)
end subroutine generate_flux_surface_3D
