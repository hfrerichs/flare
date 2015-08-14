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
  use run_control, only: Theta, Psi, n_theta, n_psi, Phi_output, Grid_File, Output_File, Debug
  use flux_surface_2D
!  use magnetic_axis
  use equilibrium
  use parallel
  use dataset
  use grid
  implicit none

  type(t_flux_surface_2D) :: S
  type(t_dataset)         :: D
  type(t_grid)            :: G1, G2, G_debug
  real(real64) :: x0(2), xc(3), t, x(3), y(3)
  integer :: i, j, ig, n, ierr, iterations


  if (firstP) then
     write (6, *) 'Generate axisymmetric flux surface grid, output in: ', adjustl(trim(Grid_File))
     write (6, *)
  endif


  call G1%new(CYLINDRICAL, MESH_2D, TOROIDAL_SLICE, n_psi, n_theta, fixed_coord_value=Phi_Output)
  call G2%new(LOCALL,      MESH_2D, TOROIDAL_SLICE, n_psi, n_theta, fixed_coord_value=Phi_Output)

  call G_debug%new(LOCAL, STRUCTURED, FIXED_COORD3, n_theta, n_psi)
  n  = n_theta * n_psi
  ig = 1
  call D%new(n,2)
  do i=0, n_psi-1
     do j=0, n_theta-1
        y(2) = Psi(1) + 1.d0*i/(n_psi-1) * (Psi(2)-Psi(1))
        y(1) = Theta(1) + 1.d0*j/(n_theta-1) * (Theta(2)-Theta(1))
        y(3) = Phi_output

        x = get_cylindrical_coordinates (y, ierr, iter=iterations)
        !if (ierr > 0) x = get_cylindrical_coordinates (y, ierr, damping=0.01d0, iter=iterations)
        G1%mesh(i,j,1:2) = x(1:2)
        G2%mesh(i,j,1:2) = y(1:2)
        D%x(ig,1) = ierr
        D%x(ig,2) = iterations
        ig        = ig + 1

        G_debug%x1(j+1) = j
        G_debug%x2(i+1) = i
     enddo
  enddo
  call G1%store(filename=Grid_File)
  call G2%store(filename=Output_File)
 3000 format (2e18.10)
  if (Debug) then
     call D%plot(filename='debug.dat')
     call G_debug%store(filename='debug.grid')
  endif


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
