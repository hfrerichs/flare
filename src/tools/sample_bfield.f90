!===============================================================================
! Sample magnetic field
!
! Input (taken from run control file):
!    Grid_File          Sample locations
!    Output_File
!    Output_Format      = 1: Cartesian coordinate system (Bx,By,Bz)
!                       = 2: Cylindrical coordinate system (BR,BZ,Bphi)
!===============================================================================
subroutine sample_bfield
  use run_control, only: Grid_File, Output_File, Output_Format
  use parallel
  use bfield
  use grid
  implicit none

  integer, parameter :: iu = 42

  real*8  :: Bf(3), xvec(3)
  integer :: iflag


  if (firstP) then
    write (6, *) 'Sample magnetic field, output in: ', adjustl(trim(Output_File))
    write (6, *)
  endif
  open  (iu, file=Output_File, err=5010)


  if (Output_Format < 1 .or. Output_Format > 2) then
     write (6, *) 'undefined output format ', Output_Format
     stop
  endif
  write (iu, 1000) COORDINATES(Output_Format)
  call read_grid (Grid_File, use_coordinates=COORDINATES(Output_Format))
  grid_point_loop: do
     call get_next_grid_point (iflag, xvec)
     if (iflag.eq.-1) exit grid_point_loop

     select case (Output_Format)
     case (1)
        Bf = get_Bf_Cart (xvec)
     case (2)
        Bf = get_Bf_Cyl (xvec)
     end select

     Bf = Bf/1.d4	! Gauss -> Tesla
     write (iu,1001) Bf
  enddo grid_point_loop
  close (iu)

  return
 1000 format ('# Magnetic field components [Tesla], coordinate system: ', a12)
 1001 format (3e16.8)
 5010 write (6,5011) Output_File
 5011 format ('error opening file for output: ', a120)
  stop
end subroutine sample_bfield
