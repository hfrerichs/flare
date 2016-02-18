!===============================================================================
! Calculate field line losses (e.g. from former closed magnetic flux surfaces)
!
! Input (taken from run control file):
!
!    N_sym              Toroidal symmetry
!    N_phi              Number of slices between 0 and 360/N_sym deg
!    N_theta            Poloidal resolution of initial surface
!    N_psi              Number of initial surfaces
!    Psi(1:2)           Lower and upper radial boundary
!
!    Limit              Reference length for field line tracing
!    N_steps            Calculate loss fraction at N_steps steps of 'Limit' (default = 1)
!    Trace_Step, Trace_Method, Trace_Coords as ususal
!
!    Output_File
!    Grid_File
!===============================================================================
subroutine field_line_loss
  use iso_fortran_env
  use flux_surface_2D
  use flux_surface_3D
  use run_control, only: N_sym, N_phi, N_theta, Psi, N_psi, N_steps, Limit, &
                         Output_File, Grid_File, Debug
  use equilibrium
  use fieldline
  use parallel
  use grid
  implicit none

  integer, parameter      :: iu = 80

  integer, dimension(:,:), allocatable :: iloss
  type(t_flux_surface_2D), dimension(:), allocatable :: S2D
  type(t_flux_surface_3D), dimension(:), allocatable :: S3D
  type(t_grid) :: G
  real(real64) :: r(3), y(3), DPsi
  integer      :: i, ierr, j, n


  ! check user input, set default values
  if (N_steps <=0) N_steps = 1


  ! reference output to screen
  if (firstP) then
     write (6, *) 'Calculate field line losses, output in: ', adjustl(trim(Output_File))
     write (6, *)
     if (N_Psi == 1) then
        write (6, 1001) Psi(1)
     else
        write (6, 1002) Psi(1), Psi(2)
     endif
     write (6, *)

     write (6, 1003)
     write (6, 1004) N_sym
     write (6, 1005) N_phi, N_theta
     write (6, *)

     if (N_steps == 1) then
        write (6, 1006) Limit/1.d2
     else
        write (6, 1007) N_steps, Limit/1.d2, N_steps*Limit/1.d2
     endif
     write (6, *)
 1001 format(3x,'- Radial position [PsiN]:         ',5x,f0.4)
 1002 format(3x,'- Radial domain [PsiN]:           ',5x,f0.4,' -> ',f0.4)
 1003 format(3x,'- Flux surface discretization:')
 1004 format(8x,'Toroidal symmetry:                ',i0)
 1005 format(8x,'Toroidal x Poloidal resolution:   ',i0,' x ',i0)
 1006 format(3x,'- Reference length [m]:           ',5x,f0.2)
 1007 format(3x,'- Reference length [m]:           ',5x,i0,' x ',f0.2,' = ',f0.2)

     open  (iu, file=Output_File)
     write (iu, 2000)
  endif
  call G%new(CYLINDRICAL, MESH_2D, FIXED_COORD3, N_Psi, N_theta)


  ! generate unperturbed flux surfaces Psi(1) -> Psi(2) for reference points
  allocate (S2D(0:N_psi-1), S3D(0:N_psi-1))
  if (firstP) write (6, *) 'Generating unperturbed flux surfaces for reference points ...'
  y(1) = 0.d0
  y(3) = 0.d0
  DPsi = 0.d0; if (N_psi > 1) DPsi = (Psi(2)-Psi(1)) / (N_psi - 1.d0)
  do i=0,N_psi-1
     y(2) = Psi(1) + i * DPsi
     if (firstP) write (6, *) y(2)

     ! 1. find reference coordinate on flux surface (inboard midplane)
     r    = get_cylindrical_coordinates (y, ierr)
     if (ierr .ne. 0) then
        write (6, *) 'error in subroutine field_line_loss: could not find real space coordinates!'
        write (6, *) r
        stop
     endif


     ! 2. generate flux surface shape
     call S2D(i)%generate_closed(r(1:2))
     if (Debug) call S2D(i)%plot(filename='flux_surfaces.plt', append=.true.)
     call S3D(i)%generate_from_axisymmetric_surface(S2D(i), N_sym, N_phi, N_theta)
     G%mesh(i,:,1:2) = S3D(i)%slice(0)%x(1:N_theta,1:2)
  enddo
  if (firstP) then
     write (6, *) ' -> stored in ', adjustl(trim(Grid_File))
     call G%store(Grid_File)
  endif
  write (6, *)



  ! initialize output variables
  allocate (iloss(-1:1, N_steps))
  iloss = 0

  ! loop over all unperturbed flux surfaces Psi(1) -> Psi(2)
  if (firstP) write (6, *) 'Starting field line loss calculation ...'
  do i=0,N_psi-1
     y(2) = Psi(1) + i * DPsi
     if (firstP) write (6, *) y(2)

     ! 3. calculate field line losses from this unperturbed flux surface
     iloss = S3D(i)%field_line_loss(N_steps, Limit)
     n     = N_phi * N_theta
     if (firstP) then
        do j=1,N_steps
           write (iu, 2001) y(2), j*Limit, 1.d0*iloss(-1,j)/n, 1.d0*iloss(1,j)/n
        enddo
     endif
  enddo


  ! finalize
  if (firstP) then
     write (6, *) 'done'
     close (iu)
  endif
  deallocate (iloss, S2D, S3D)

 2000 format('# PsiN,       L               loss(backward)  loss(forward)')
 2001 format(f12.8,3e16.8)
end subroutine field_line_loss
