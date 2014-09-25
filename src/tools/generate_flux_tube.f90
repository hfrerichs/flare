!===============================================================================
! Generate finite flux tube
!
! Input (taken from run control file):
!    Grid_File          Initial shape of flux tube
!    N_phi              Number of toroidal slices
!    N_mult             Total length of flux tube will be Delta_Phi = 360 deg / N_mult
!
!    Output_File
!===============================================================================
subroutine generate_flux_tube
  use iso_fortran_env
  use run_control, only: Grid_File, Output_File, N_phi, N_mult
  use parallel
  use grid
  use flux_tube
  use math
  implicit none

  type(t_grid)          :: G
  type(t_cross_section) :: cs0
  type(t_flux_tube)     :: FT


  if (firstP) then
     write (6, *) 'Generate finite flux tube, output in ', trim(Output_File)
     write (6, *)
  else
     return
  endif


  ! load initial cross-section of flux tube
  call cs0%load(Grid_File)


  ! generate flux tube and write to Output_File
  call FT%generate(cs0, N_phi, N_mult)
  call FT%plot(Output_File)

end subroutine generate_flux_tube
