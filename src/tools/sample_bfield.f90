subroutine sample_bfield
  use bfield
  use run_control, only: Grid_File, Output_File
  use grid
  implicit none


  integer, parameter :: iu = 42
  integer :: iflag

      real*8 :: Bf(3), xvec(3), psi, psi_sepx, psi_axis, R_axis, Z_axis
      real*8 :: tH(3), eH(3), theta, phi, R, bphi



  write (6,1000) Output_File
  call read_grid (Grid_File)

  open  (iu, file=Output_File, err=5010)

  grid_point_loop: do
     call get_next_grid_point (iflag, xvec)
     if (iflag.eq.-1) exit grid_point_loop

     Bf = get_Bf_Cart (xvec)

     Bf = Bf/1.d4	! Gauss -> Tesla
     write (iu,1001) Bf
  enddo grid_point_loop
  close (iu)

  return
 1000 format (3x,'- sample magnetic field, output in: ',a120)
 1001 format (3e16.8)
 5010 write (6,5011) Output_File
 5011 format ('error opening file for output: ', a120)
  stop
end subroutine sample_bfield
