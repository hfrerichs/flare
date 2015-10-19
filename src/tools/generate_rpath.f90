!===============================================================================
! Generate path across magnetic surfaces
!
! Input (taken from run control file):
!    N_psi              Start at X-point number N_psi
!    Limit              Length of path from X-point
!    N_steps            Number of discretization points used for Grid_File
!    Grid_File
!    Output_File
!===============================================================================
subroutine generate_rpath
  use iso_fortran_env
  use run_control, only: Grid_File, Output_File, N_steps, N_psi, Limit
  use xpaths
  use grid
  implicit none

  type(t_xpath) :: rpath
  type(t_grid)  :: G
  real(real64)  :: t, x(2)
  integer       :: i


  call rpath%generateX(N_psi, ASCENT, LIMIT_LENGTH, Limit, SAMPLE_LENGTH)
  call rpath%plot(filename=Output_File)

  if (N_steps == 0) return
  call G%new(CYLINDRICAL, UNSTRUCTURED, FIXED_COORD3, N_steps+1)
  do i=0,N_steps
     t = 1.d0 * i / N_steps
     call rpath%sample_at(t, x)
     G%x(i+1,1:2) = x
  enddo
  call G%store(filename=Grid_File)

end subroutine generate_rpath
