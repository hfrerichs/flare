!===============================================================================
! Sample magnetic field
!
! Input (taken from run control file):
!    Grid_File          Sample locations
!
!    Output_File
!    Output_Format      = 1: |B|, Cartesian components (Bx,By,Bz), PsiN
!                       = 2: |B|, Cylindrical components (BR,BZ,Bphi), PsiN
!                       = 3: (Output_Format 2), Bpol/Btor
!                       = 4: |grad PsiN|, e_Psi * e_B
!===============================================================================
subroutine sample_bfield
  use iso_fortran_env
  use run_control, only: Grid_File, Output_File, Output_Format
  use parallel
  use math
  use bfield
  use equilibrium
  use grid
  use dataset
  implicit none

  integer, parameter :: iu = 42

  type(t_grid)    :: G
  type(t_dataset) :: D
  real(real64)    :: Bmod, Bf(3), xvec(3), r(3), PsiN, Bpol, rpol
  real(real64)    :: DPsiDR, DPsiDZ
  integer         :: ig, n


  if (firstP) then
     write (6, *) 'Sample magnetic field, output in: ', adjustl(trim(Output_File))
     write (6, *)
     open  (iu, file=Output_File, err=5010)

     if (Output_Format < 1 .or. Output_Format > 4) then
        write (6, *) 'undefined output format ', Output_Format
        stop
     endif
     write (iu, 1000) COORDINATES(min(Output_Format,2))
  endif


  ! load grid
  call G%load(Grid_File)


  ! initialize output data
  n = 5
  if (Output_Format == 3) n = 6
  if (Output_Format == 4) n = 2
  call D%new(G%nodes(), n)


  ! main loop
  grid_loop: do ig=mype+1,G%nodes(),nprs
     xvec = G%node(ig, coordinates=min(Output_Format,2))

     select case (Output_Format)
     case (1)
        Bf = get_Bf_Cart (xvec)
     case (2,3,4)
        Bf = get_Bf_Cyl (xvec)
     end select
     Bf   = Bf/1.d4	! Gauss -> Tesla
     Bmod = sqrt(sum(Bf**2))


     ! get normalized poloidal flux
     call coord_trans (xvec, Output_Format, r, CYLINDRICAL)
     PsiN = get_PsiN(r)


     if (Output_Format <= 3) then
        D%x(ig,1)   = Bmod
        D%x(ig,2:4) = Bf
        D%x(ig,5)   = PsiN
     endif
     if (Output_Format == 3) then
        ! ratio of poloidal to toroidal field
        Bpol = sqrt(Bf(1)**2 + Bf(2)**2)
        rpol = Bpol / Bf(3)
        D%x(ig,6)   = rpol
     endif
     if (Output_Format == 4) then
        DPsiDR      = get_DPsiN(r, 1, 0)
        DPsiDZ      = get_DPsiN(r, 0, 1)
        D%x(ig,1)   = sqrt(DPsiDR**2 + DPsiDZ**2)
        D%x(ig,2)   = (DPsiDR*Bf(1) + DPsiDZ*Bf(2)) / D%x(ig,1) / Bmod
     endif
  enddo grid_loop
  call D%mpi_allreduce()


  ! write data
  if (firstP) then
     call D%plot(iu, formatstr='e22.14')
     close (iu)
  endif

  return
 1000 format ('# Magnetic field components [Tesla], coordinate system: ', a12)
 5010 write (6,5011) Output_File
 5011 format ('error opening file for output: ', a120)
  stop
end subroutine sample_bfield
