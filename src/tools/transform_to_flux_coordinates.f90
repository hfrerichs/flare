!===============================================================================
! Perform coordinate transformation (to flux coordinates) for a given set of points
!
! Input (taken from run control file):
!    Grid_File          Set of coordinates to transform
!    Input_Format       = 1: x, y, z [cm]
!                       = 2: R, Z [cm], phi [deg]
!                       = 3: R, Z [cm] at Phi_output [deg]
!    Phi_Output         Default position for input format 3
!    Output_Format      = 1: 0 <= theta <= 360 deg
!                         2: -180 <= theta <= 180 deg
!
!    Output_File
!===============================================================================
subroutine transform_to_flux_coordinates
  use run_control, only: Grid_File, Input_Format, Phi_Output, Output_File, Output_Format
  use parallel
  use equilibrium
  use curve2D
  use grid
  use math
  implicit none

  integer, parameter :: iu = 42

  type(t_grid)  :: G
  real*8, dimension(:,:), allocatable :: X
  character*120 :: str
  real*8  :: r(3), theta, PsiN
  integer :: i, ig, n_grid


  if (firstP) then
     write (6, *) 'Transform to flux coordinates, output in: ', adjustl(trim(Output_File))
     write (6, *)
  endif


  ! read input coordinates
  open  (iu, file=Grid_File)
  read  (iu, '(a120)') str
  close (iu)
  if (str(3:9) .eq. 'grid_id') then
     call G%load(Grid_File)
  else
     call G%load_usr(Grid_File, Input_Format, Phi_Output)
  endif
  n_grid = G%nodes()
  allocate (X(n_grid,3))


  ! transform to flux coordinates
  grid_loop: do ig=mype,n_grid,nprs
     r = G%node(ig+1)

     theta   = get_poloidal_angle(r) * 180.d0 / pi
     if (Output_Format == 1  .and.  theta < 0.d0) theta = theta + 360.d0
     PsiN    = get_PsiN(r)
     X(ig,1) = theta
     X(ig,2) = PsiN
     X(ig,3) = r(3)
  enddo grid_loop


  ! write output data
  call wait_pe()
  call sum_real_data (X, n_grid*3)
  if (firstP) then
     open  (iu, file=Output_File)
     do ig=1,n_grid
        write (iu, 1000) X(ig,:)
     enddo
     close (iu)
  endif


  ! cleanup
  deallocate (X)

 1000 format (3e18.10)
end subroutine transform_to_flux_coordinates
