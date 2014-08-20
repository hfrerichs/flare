!===============================================================================
! Sample magnetic field
!
! Input (taken from run control file):
!    Grid_File          Sample locations
!
!    Output_File
!    Output_Format      = 1: Cartesian components (Bx,By,Bz), PsiN
!                       = 2: Cylindrical components (BR,BZ,Bphi), PsiN
!                       = 3: (Output_Format 2), Bpol/Btor
!===============================================================================
subroutine sample_bfield
  use run_control, only: Grid_File, Output_File, Output_Format
  use parallel
  use bfield
  use equilibrium
  use grid
  use math
  implicit none

  integer, parameter :: iu = 42

  real*8  :: Bf(3), xvec(3), r(3), PsiN, Bpol, rpol
  integer :: iflag


  if (firstP) then
     write (6, *) 'Sample magnetic field, output in: ', adjustl(trim(Output_File))
     write (6, *)
  endif
  open  (iu, file=Output_File, err=5010)


  if (Output_Format < 1 .or. Output_Format > 3) then
     write (6, *) 'undefined output format ', Output_Format
     stop
  endif
  write (iu, 1000) COORDINATES(min(Output_Format,2))
  call read_grid (Grid_File, use_coordinates=COORDINATES(min(Output_Format,2)))
  grid_point_loop: do
     call get_next_grid_point (iflag, xvec)
     if (iflag.eq.-1) exit grid_point_loop

     select case (Output_Format)
     case (1)
        Bf = get_Bf_Cart (xvec)
     case (2,3)
        Bf = get_Bf_Cyl (xvec)
     end select
     Bf = Bf/1.d4	! Gauss -> Tesla


     ! get normalized poloidal flux
     call coord_trans (xvec, Output_Format, r, CYLINDRICAL)
     PsiN = get_PsiN(r)


     if (Output_Format .le. 2) then
        write (iu,1001) Bf, PsiN
     else
        ! ratio of poloidal to toroidal field
        Bpol = sqrt(Bf(1)**2 + Bf(2)**2)
        rpol = Bpol / Bf(3)
        write (iu,1002) Bf, PsiN, rpol
     endif
  enddo grid_point_loop
  close (iu)

  return
 1000 format ('# Magnetic field components [Tesla], coordinate system: ', a12)
 1001 format (4e22.14)
 1002 format (5e22.14)
 5010 write (6,5011) Output_File
 5011 format ('error opening file for output: ', a120)
  stop
end subroutine sample_bfield
