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
!    Limit              Max. length of field line tracing
!    N_steps            Calculate loss fraction at N_steps steps of 'Limit'
!    Trace_Step, Trace_Method, Trace_Coords as ususal
!
!    Output_File
!===============================================================================
subroutine field_line_loss
  use iso_fortran_env
  use flux_surface_2D
  use flux_surface_3D
  use run_control, only: N_sym, N_phi, N_theta, Psi, N_psi, N_steps, Limit, Output_File, Debug
  use equilibrium
  use fieldline
  use parallel
  implicit none

  integer, parameter      :: iu = 80

  integer, dimension(:,:), allocatable :: iloss
  type(t_flux_surface_2D) :: S2D
  type(t_flux_surface_3D) :: S3D
  real(real64) :: r(3), y(3), DPsi
  integer      :: i, ierr, j, n


  if (firstP) then
     write (6, *) 'Calculate field line losses within ',Limit/1.d2,' m, output in: ', adjustl(trim(Output_File))
     write (6, *)

     open  (iu, file=Output_File)
     write (iu, 1000)
  endif


  ! check user input, set default values
  if (N_steps <=0) N_steps = 1


  ! initialize output variables
  allocate (iloss(-1:1, N_steps))
  iloss = 0

  ! loop over all unperturbed flux surfaces Psi(1) -> Psi(2)
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
     call S2D%generate_closed(r(1:2))
     if (Debug) call S2D%plot(filename='flux_surfaces.plt', append=.true.)
     call S3D%generate_from_axisymmetric_surface(S2D, N_sym, N_phi, N_theta)

  
     ! 3. calculate field line losses from this unperturbed flux surface
     iloss = S3D%field_line_loss(N_steps, Limit)
     n     = N_phi * N_theta
     if (firstP) then
        do j=1,N_steps
           write (iu, 1001) y(2), j*Limit, 1.d0*iloss(-1,j)/n, 1.d0*iloss(1,j)/n
        enddo
     endif
  enddo


  ! finalize
  if (firstP) then
     close (iu)
  endif
  deallocate (iloss)

 1000 format('# PsiN,       L               loss(backward)  loss(forward)')
 1001 format(f12.8,3e16.8)
end subroutine field_line_loss
