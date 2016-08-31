subroutine perturbation_spectrum()
  use iso_fortran_env
  use run_control, only: Psi, N_psi, N_theta, Grid_File, Output_File
  use flux_surface_2D
  use grid
  implicit none

  integer, parameter      :: iu = 98

  type(t_flux_surface_2D) :: F
  type(t_grid)            :: G

  !real(real64) :: Bpsi_harm(-50:50)
  real(real64), dimension(:,:), allocatable :: Bpsi_harm
  real(real64), dimension(:),   allocatable :: Bpsi_harm1
  real(real64) :: Psi0
  integer :: i, ig, j, n


  n = 3

  call G%new(LOCAL, UNSTRUCTURED, 3, (2*N_theta+1)*(N_psi+1))

  allocate (Bpsi_harm(0:N_psi,-N_theta:N_theta))
  allocate (Bpsi_harm1(-N_theta:N_theta))

  open  (iu, file=Output_File)
  ig = 0
  do i=0,N_psi
     write (6, *) i
     Psi0 = Psi(1) + (Psi(2)-Psi(1)) / N_psi * i

     call F%setup_theta_map(Psi0)
     write (90, *) -n*F%q, Psi0

     call F%get_spectrum(n, N_theta, Bpsi_harm1)
     Bpsi_harm(i,:) = Bpsi_harm1

     do j=-N_theta,N_theta
        ig = ig + 1
        G%x(ig,1) = j
        G%x(ig,2) = Psi0

        write (iu, *) Bpsi_harm1(j)
     enddo
  enddo
  close (iu)
  call G%store(Grid_File)

  deallocate (Bpsi_harm, Bpsi_harm1)

end subroutine perturbation_spectrum
