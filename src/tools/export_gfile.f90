!===============================================================================
! Export g-file (2D), based on module equilibrium
!    - Evaluate function get_Psi
!    - Use plasma current Ip (only sign is relevant, see CurrentFix in module geqdsk)
!    - Use reference toroidal field Bt at R0
!    - Use 1st axisymmetric boundary
!
! Input (taken from run control file):
!    R_start, R_end     Computational domain, equidistant discretization
!    Z_start, Z_end
!
!    N_R, N_Z           Resolution of computational domain (number of grid nodes)
!
!    Output_File        filename for output g-file
!===============================================================================
subroutine export_gfile
  use iso_fortran_env
  use run_control, only: R_start, R_end, Z_start, Z_end, N_R, N_Z, Output_File
  use equilibrium
  use boundary
  use parallel
  implicit none

  integer, parameter :: iu = 42


  ! internal variables
  real(real64), dimension(:,:), allocatable :: psirz
  real(real64), dimension(:), allocatable :: empty, fpol
  real(real64) :: R, Z, Zmid, DeltaR, DeltaZ, Raxis, Zaxis, r3(3)
  integer      :: i, j, nr, nz


  if (firstP) then
     write (6, 1000) adjustl(Output_File)
     write (6, 1001)
     write (6, 1002) R_start, R_end
     write (6, 1003) Z_start, Z_end
     write (6, 1004) N_R, N_Z
     write (6, *)
  else
     return
  endif


  ! set up dependent geometry parameters
  Zmid   = 0.5d0 * (Z_start + Z_end)
  DeltaR = R_end - R_start
  DeltaZ = Z_end - Z_start
  nr     = N_R
  nz     = N_Z
  ! magnetic axis
  r3     = get_magnetic_axis(0.d0)
  Raxis  = r3(1); Zaxis = r3(2)


  allocate (psirz(nr,nz), empty(nr), fpol(nr))
  empty = 0.d0
  fpol  = Bt * (R0/1.d2)


  ! set up Psi(R,Z) for output
  do j=1,nz
  do i=1,nr
     r3(1)      = R_start + (i-1.d0) / (nr-1.d0) * DeltaR
     r3(2)      = Z_start + (j-1.d0) / (nz-1.d0) * DeltaZ
     psirz(i,j) = get_Psi(r3)
  enddo
  enddo


  ! convert units
  DeltaR  = DeltaR  / 1.d2
  DeltaZ  = DeltaZ  / 1.d2
  R_start = R_start / 1.d2
  R_end   = R_end   / 1.d2
  Raxis   = Raxis   / 1.d2
  Zaxis   = Zaxis   / 1.d2
  R0      = R0      / 1.d2
  Zmid    = Zmid    / 1.d2



  ! write gfile
  open  (iu, file=Output_File)
  write (iu, 2000) 0, nr, nz
  write (iu, 2020) DeltaR, DeltaZ, R0, R_start, Zmid
  write (iu, 2020) Raxis, Zaxis, Psi_axis, Psi_sepx, Bt
  write (iu, 2020) Ip*1.d6, Psi_axis, 0.d0, Raxis, 0.d0
  write (iu, 2020) Zaxis, 0.d0, Psi_sepx, 0.d0, 0.d0

  write (iu, 2020) fpol
  write (iu, 2020) empty
  write (iu, 2020) empty
  write (iu, 2020) empty
  write (iu, 2020) psirz
  write (iu, 2020) empty

  ! plasma boundary (not implemented) and wall/limiter
  if (n_axi == 0) then
     write (iu, 2022) 1, 1
     write (iu, 2020) 0.d0, 0.d0
     write (iu, 2020) 0.d0, 0.d0
  else
     write (iu, 2022) 1, S_axi(1)%n_seg+1
     write (iu, 2020) 0.d0, 0.d0
     write (iu, 2020) (S_axi(1)%x(i,1)/1.d2, S_axi(1)%x(i,2)/1.d2, i=0,S_axi(1)%n_seg)
  endif
  close (iu)
  deallocate (psirz, empty, fpol)

 1000 format(3x,'- Export g-file, output in: ', a)
 1001 format(8x,'Computational domain:')
 1002 format(8x,'R: ',f12.6,' -> ',f12.6)
 1003 format(8x,'Z: ',f12.6,' -> ',f12.6)
 1004 format(8x,'Resolution (number of grid nodes): ',i0,' x ',i0)
 2000 format(48x,3i4)
 2020 format(5e16.9)
 2022 format(2i5)
end subroutine export_gfile
