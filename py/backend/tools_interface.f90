subroutine generate_magnetic_axis_interface(x0, nsym, nphi)
  use run_control, only: N_sym, N_phi, x_start
  use types
  implicit none

  real(dp), intent(in) :: x0(3)
  integer,  intent(in) :: nsym, nphi


  N_sym = nsym;   N_phi = nphi;   x_start = x0
  call generate_magnetic_axis()

end subroutine generate_magnetic_axis_interface
