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
  real(real64) :: x0(2), xc(3), t, x(3), y(3), dPsi, dTheta
  integer :: i, j, j2, ig, n, ierr, iterations


  if (firstP) then
     write (6, *) 'Generate axisymmetric flux surface grid, output in: ', adjustl(trim(Grid_File))
     write (6, *)
  endif


  call G1%new(CYLINDRICAL, MESH_2D, TOROIDAL_SLICE, n_psi, n_theta, fixed_coord_value=Phi_Output)
  call G2%new(LOCAL,       MESH_2D, TOROIDAL_SLICE, n_psi, n_theta, fixed_coord_value=Phi_Output)

  call G_debug%new(LOCAL, STRUCTURED, FIXED_COORD3, n_theta, n_psi)
  n  = n_theta * n_psi
  ig = 1
  call D%new(n,2)
  dTheta = 0.d0; if (n_theta > 1) dTheta = (Theta(2)-Theta(1)) / (n_theta - 1)
  dPsi   = 0.d0; if (n_psi   > 1) dPsi   = (Psi(2)  -Psi(1)  ) / (n_psi   - 1)
  radial_loop: do i=0, n_psi-1
     poloidal_loop: do j=0, n_theta-1
        y(2) = Psi(1)   + i * dPsi
        y(1) = Theta(1) + j * dTheta
        y(3) = Phi_output

        x = get_cylindrical_coordinates (y, ierr, iter=iterations)
        !if (ierr > 0) x = get_cylindrical_coordinates (y, ierr, damping=0.01d0, iter=iterations)
        G_debug%x1(j+1) = j
        G_debug%x2(i+1) = i

        ! failed to find cylindrical coordinates for y from function get_cylindrical_coordinates
        if (ierr > 0  .and.  .not.Debug) then
           if (j == 0) then
              write (6, *) 'error: cannot find first point on surface ', i
              stop
           endif
           write (6, 9000) Psi(1)   + i * dPsi, Theta(1) + j * dTheta

           ! generate flux surface from first point
           x(1:2) = G1%mesh(i,0,1:2)
           x(3)   = Phi_output
           write (6, 9001) x(1:2)
           call S%generate_closed(x)
           call S%setup_angular_sampling()
           do j2=0,n_theta-1
              t    = 1.d0 * j2 / n_theta
              y(1) = Theta(1) + j * dTheta
              call S%sample_at(t, x(1:2))
              G1%mesh(i,j2,1:2) = x(1:2)
              G2%mesh(i,j2,1:2) = y(1:2)
           enddo

           exit poloidal_loop
        endif

        G1%mesh(i,j,1:2) = x(1:2)
        G2%mesh(i,j,1:2) = y(1:2)
        D%x(ig,1) = ierr
        D%x(ig,2) = iterations
        ig        = ig + 1
     enddo poloidal_loop
  enddo radial_loop
  call G1%store(filename=Grid_File)
  call G2%store(filename=Output_File)
  if (Debug) then
     call D%plot(filename='debug.dat')
     call G_debug%store(filename='debug.grid')
  endif

 9000 format ('function get_cylindrical_coordinates failed at Psi, Theta = ',f10.5,', ',f10.5)
 9001 format ('flux surface is generated from ',f10.5,', ',f10.5,' by field line tracing')
end subroutine flux_surface_grid
