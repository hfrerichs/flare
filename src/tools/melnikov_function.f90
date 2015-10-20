!===============================================================================
! Calculate (subharmonic) Melnikov function
!
! ATTENTION: This subroutine is based on an axisymmetric equilibrium (BF_EQ2D) and
!            a perturbation field from an explicit coil set (BF_COILS)
!
!===============================================================================
subroutine melnikov_function
  use parallel
  use run_control, only: Grid_File, Output_File, Trace_Step
  use grid
  use fieldline
  use dataset
  use math
  use bfield
  use equilibrium
  use polygones
  implicit none

  type(t_fieldline) :: F
  type(t_dataset)   :: D
  type(t_grid) :: G
  real(real64) :: y(3), x(3), Bpert(3), DPsi, GradPsi(2), DPsi_int(-1:1), Bf(3), Bmod, dl, PsiN
  integer      :: ig, idir


  if (firstP) then
     write (6, *) 'Calculate Melnikov function, output in: ', adjustl(trim(Output_File))
     write (6, *)
  endif


  ! load grid
  call G%load(Grid_File)


  ! set equilibrium field
  iconfig = 0
  iconfig(BF_EQ2D) = 1


  call D%new(G%nodes(), 4)
  grid_loop: do ig=mype,G%nodes()-1,nprs
     write (6, *) ig
     y        = G%node(ig+1, coordinates=CYLINDRICAL)
     PsiN     = get_PsiN(y)
     DPsi_int = 0.d0
     do idir=-1,1,2
        call F%init(y, idir*Trace_Step, NM_AdamsBashforth4, CYLINDRICAL)

        trace_loop: do
           dl         = F%trace_1step()
           x          = F%yc

           Bpert      = get_Bcyl_polygones(x)
           GradPsi(1) = get_DPsiN(x, 1, 0)
           GradPsi(2) = get_DPsiN(x, 0, 1)
           Bf         = get_Bf_Cyl(x)
           Bmod       = sqrt(sum(Bf**2))
           DPsi       = Bpert(1)*GradPsi(1) + Bpert(2)*GradPsi(2)
           DPsi       = DPsi / Bmod * dl
           DPsi_int(idir) = DPsi_int(idir) + DPsi

           if (abs(F%theta_int) .ge. pi2) exit trace_loop
           if (PsiN > 1.d0  .and.  F%intersect_boundary()) exit trace_loop
        enddo trace_loop
     enddo

     D%x(ig+1, 1) = DPsi_int(-1)
     D%x(ig+1, 2) = DPsi_int( 1)
     D%x(ig+1, 3) = PsiN + DPsi_int(-1)
     D%x(ig+1, 4) = PsiN + DPsi_int( 1)
  enddo grid_loop
  call D%mpi_allreduce()

  ! write data
  if (firstP) then
     call D%plot(filename=Output_File)
  endif

end subroutine melnikov_function
