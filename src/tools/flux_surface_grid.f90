!===============================================================================
! Generate (unperturbed/axisymmetric) flux surface grid (up to PsiN = 1)
!
! Input (taken from run control file):
!
! Theta                 Poloidal range [deg] (Theta(1) -> Theta(2)) in which grid
!                       will be generated
! Psi                   Radial range (Psi(1) -> Psi(2)) in which grid will be generated
! n_theta, n_psi        Poloidal and radial resolution
!
! Phi_output            Toroidal angle [deg]
!
! Grid_File             Output file for computational grid (cylindrical coordinates)
! Output_File           Output file for plotting (flux surface coordinates)
!===============================================================================
subroutine flux_surface_grid
  use iso_fortran_env
  use run_control, only: Theta, Psi, n_theta, n_psi, Phi_output, Grid_File, Output_File
  use flux_surface_2D
!  use magnetic_axis
  use equilibrium
  use parallel
  implicit none

  integer, parameter :: &
     iu1 = 41, &
     iu2 = 42

  type(t_flux_surface_2D) :: S
  real(real64) :: x0(2), xc(3), t, x(3), y(3)
  integer :: i, j, n, ierr


  if (firstP) then
     write (6, *) 'Generate axisymmetric flux surface grid, output in: ', adjustl(trim(Grid_File))
     write (6, *)
  endif


  open  (iu1, file=Grid_File)
  write (iu1, 1000)
  write (iu1, 1001) n_theta * n_psi
  write (iu1, 1002) Phi_Output
 1000 format ('# grid_id = 1       (irregular RZ grid)')
 1001 format ('# resolution:        n_grid  =  ',i10)
 1002 format ('# toroidal position: phi     =  ',f6.2)

  open  (iu2, file=Output_File)
  write (iu2, 2000)
  write (iu2, 2001) n_theta * n_psi
  write (iu2, 2002) Phi_Output
 2000 format ('# grid_id = 80      (irregular theta-psin grid)')
 2001 format ('# resolution:        n_grid  =  ',i10)
 2002 format ('# toroidal position: phi     =  ',f6.2)

  do i=0, n_psi-1
     do j=0, n_theta-1
        y(2) = Psi(1) + 1.d0*i/(n_psi-1) * (Psi(2)-Psi(1))
        y(1) = Theta(1) + 1.d0*j/(n_theta-1) * (Theta(2)-Theta(1))
        y(3) = Phi_output

        x = get_cylindrical_coordinates (y, ierr)
        write (iu1, 3000) x(1:2)
        write (iu2, 3000) y(1:2)
     enddo
  enddo
  close (iu1)
  close (iu2)
 3000 format (2e18.10)


!  x0(1) = 0.1278518850E+03
!  x0(2) = -0.1030043743E+03
!  xc    = get_magnetic_axis(0.d0)
!  call S%generate(x0)
!  call S%sort_loop(xc(1:2))
!  call S%setup_angular_sampling(xc(1:2))
!
!  n = 100
!  do i=0,n
!     t = 1.d0*i/n
!     t = Theta(1) + (Theta(2)-Theta(1)) * t
!     t = t / 360.d0
!     call S%sample_at(t, x)
!  enddo

end subroutine flux_surface_grid
