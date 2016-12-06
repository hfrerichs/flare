program benchmark_mesh_spacings
  use iso_fortran_env
  use mesh_spacing
  implicit none

  type(t_spacing) :: M


  ! 1. default spacing
  call M%init('')
  call plot(M, 'spacing_DEFAULT')

  ! 2. equidistant spacing
  call M%init('LINEAR')
  call plot(M, 'spacing_LINEAR')

  ! 3. exponential distribution
  call M%init('EXPONENTIAL 0.5d0')
  call plot(M, 'spacing_EXPONENTIAL_0.5')
  call M%init('EXPONENTIAL 0.125d0')
  call plot(M, 'spacing_EXPONENTIAL_0.125')

  ! 4. Delta-R symmetric distribution
  call M%init('DELTA_R_SYM 0.1 4')
  call plot(M, 'spacing_DELTA_R_SYM_0.1_4')
  call M%init('DELTA_R_SYM 0.2 8')
  call plot(M, 'spacing_DELTA_R_SYM_0.2_8', nsample=19)

  ! 5. Spline-X1 distribution
  call M%init('SPLINE_X1 0.8d0 0.2d0')
  call plot(M, 'spacing_SPLINE_X1_0.8_0.2')
  call M%init('SPLINE_X1 0.9d0 0.7d0')
  call plot(M, 'spacing_SPLINE_X1_0.9_0.7')

  ! 6. X1 (piecewise linear) distribution
  call M%init('X1 0.8 0.2')
  call plot(M, 'spacing_X1_0.8_0.2')
  call M%init('X1 0.9 0.7')
  call plot(M, 'spacing_X1_0.9_0.7')


  contains
  !=====================================================================


  !=====================================================================
  subroutine plot(M, label, nsample)
  type(t_spacing),  intent(in) :: M
  character(len=*), intent(in) :: label
  integer,          intent(in), optional :: nsample

  character(len=len(label)+4)  :: filename
  integer :: ns


  ns = 20
  if (present(nsample)) ns = nsample
  filename = trim(label)//'.plt';  call M%plot(filename=filename, nsample=ns, style='mesh')
  filename = trim(label)//'.dat';  call M%plot(filename=filename, nsample=ns, style='function')

  end subroutine plot
  !=====================================================================

end program benchmark_mesh_spacings
