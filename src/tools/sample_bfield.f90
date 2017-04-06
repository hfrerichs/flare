!===============================================================================
! Sample magnetic field
!
! Input (taken from run control file):
!    Grid_File          Sample locations
!
!    Output_File
!    Output_Format      = 1: |B|, Cartesian components (Bx,By,Bz), PsiN
!                       = 2: |B|, Cylindrical components (BR,BZ,Bphi), PsiN
!                       = 3: Bpol, Btor, R, Bpol / Btor
!                       = 4: |grad PsiN|, e_Psi * e_B
!                       = 5: divB
!                       = 6: Psi, PsiN, dPsi_dR, dPsi_dZ
!                       = 7: Psi, d2Psi_dR2, d2Psi_dRdZ, d2Psi_dZ2
!                       = 8: toroidal, poloidal and radial components of non-equilibrium field
!                       = 9: Jacobian
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
  real(real64)    :: Bmod, Bf(3), xvec(3), r(3), PsiN, Psi, Bpol, divB, Beq(3), JBf(3,3)
  real(real64)    :: DPsiDR, DPsiDZ, ePsi(2), ePol(2)
  integer         :: ig, n


  if (firstP) then
     write (6, *) 'Sample magnetic field, output in: ', adjustl(trim(Output_File))
     write (6, *)

     if (Output_Format < 1 .or. Output_Format > 9) then
        write (6, *) 'undefined output format ', Output_Format
        stop
     endif

     open  (iu, file=Output_File, err=5010)
     write (iu, 1000) COORDINATES(min(Output_Format,2))
     write (iu, 1001) Output_Format
  endif


  ! load grid
  call G%load(Grid_File)


  ! initialize output data
  n = 5
  if (Output_Format == 3) n = 4
  if (Output_Format == 4) n = 2
  if (Output_Format == 5) n = 4
  if (Output_Format == 6) n = 4
  if (Output_Format == 7) n = 4
  if (Output_Format == 8) n = 4
  if (Output_Format == 9) n = 9
  call D%new(G%nodes(), n)


  ! main loop
  grid_loop: do ig=mype+1,G%nodes(),nprs
     xvec = G%node(ig, coordinates=min(Output_Format,2))

     Bf = 0.d0
     select case (Output_Format)
     case (1)
        Bf = get_Bf_Cart (xvec)
     case (2,3,4,5,6,8)
        Bf = get_Bf_Cyl (xvec)
     end select
     Bf   = Bf/1.d4	! Gauss -> Tesla
     Bmod = sqrt(sum(Bf**2))


     ! get normalized poloidal flux
     call coord_trans (xvec, Output_Format, r, CYLINDRICAL)
     PsiN = get_PsiN(r)

     ! radial direction
     ePsi(1) = get_DPsiN(r, 1, 0)
     ePsi(2) = get_DPsiN(r, 0, 1)
     ePsi    = ePsi / sqrt(sum(ePsi**2))

     ! poloidal field and direction
     Beq     = get_Bf_eq2D(r)/1.d4
     Bpol    = sqrt(Beq(1)**2 + Beq(2)**2)
     ePol    = Beq(1:2) / Bpol


     if (Output_Format < 3) then
        D%x(ig,1)   = Bmod
        D%x(ig,2:4) = Bf
        D%x(ig,5)   = PsiN
     endif
     if (Output_Format == 3) then
        ! ratio of poloidal to toroidal field
        Bpol = sqrt(Bf(1)**2 + Bf(2)**2)
        D%x(ig,1)   = Bpol
        D%x(ig,2)   = Bf(3)
        D%x(ig,3)   = r(1)
        D%x(ig,4)   = Bpol / Bf(3)
     endif
     if (Output_Format == 4) then
        DPsiDR      = get_DPsiN(r, 1, 0)
        DPsiDZ      = get_DPsiN(r, 0, 1)
        D%x(ig,1)   = sqrt(DPsiDR**2 + DPsiDZ**2)
        D%x(ig,2)   = (DPsiDR*Bf(1) + DPsiDZ*Bf(2)) / D%x(ig,1) / Bmod
     endif
     if (Output_Format == 5) then
        divB        = get_divB(r)
        D%x(ig,1)   = r(1)
        D%x(ig,2)   = r(2)
        D%x(ig,3)   = divB / 1.d4 ! Gauss/cm -> Tesla/cm
        D%x(ig,4)   = D%x(ig,3) / Bmod
     endif
     if (Output_Format == 6) then
        Psi         = get_Psi(r)
        DPsiDR      = get_DPsi(r, 1, 0)
        DPsiDZ      = get_DPsi(r, 0, 1)
        D%x(ig,1)   = Psi
        D%x(ig,2)   = PsiN
        D%x(ig,3)   = DPsiDR
        D%x(ig,4)   = DPsiDZ
     endif
     if (Output_Format == 7) then
        D%x(ig,1)   = get_Psi(r)
        D%x(ig,2)   = get_DPsi(r, 2, 0)
        D%x(ig,3)   = get_DPsi(r, 1, 1)
        D%x(ig,4)   = get_DPsi(r, 0, 2)
     endif
     if (Output_Format == 8) then
        ! non-equilibrium field
        Bf          = get_Bf_Cyl_non2D(r)/1.d4

        D%x(ig,1)   = Bf(3)
        D%x(ig,2)   = sum(Bf(1:2)*ePol)
        D%x(ig,3)   = sum(Bf(1:2)*ePsi)
        D%x(ig,4)   = D%x(ig,3) / sqrt(sum(Beq**2))
     endif
     if (Output_Format == 9) then
        JBf         = get_JBf_Cyl(r)
        D%x(ig,1:3) = JBf(1,1:3)
        D%x(ig,4:6) = JBf(2,1:3)
        D%x(ig,7:9) = JBf(3,1:3)
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
 1001 format ('# Output_Format = ', i0)
 5010 write (6,5011) Output_File
 5011 format ('error opening file for output: ', a120)
  stop
end subroutine sample_bfield
