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
!    Trace_Step, Trace_Method, Trace_Coords as ususal
!
!    Output_File
!===============================================================================
subroutine field_line_loss
  use iso_fortran_env
  use flux_surface_2D
  use flux_surface_3D
  use run_control, only: N_sym, N_phi, N_theta, Psi, N_psi, &
                         Trace_Step, Trace_Method, Trace_Coords, Limit, &
                         Output_File
  use fieldline
  use parallel
  implicit none

  integer, parameter      :: iu = 80

  type(t_flux_surface_2D) :: S2D
  type(t_flux_surface_3D) :: S3D
  real(real64) :: r(3), y(3)
  integer      :: i, ierr


  if (firstP) then
     write (6, *) 'Calculate field line losses within ',Limit/1.d2,' m, output in: ', adjustl(trim(Output_File))
     write (6, *)

     open  (iu, file=Output_File)
  endif


  ! loop over all initial surfaces
  y(1) = 0.d0
  y(3) = 0.d0
  do i=0,N_psi-1
     y(2) = Psi(1) + 1.d0*i/(N_psi-1) * (Psi(2)-Psi(1))
     if (firstP) write (6, *) y(2)

     r    = get_cylindrical_coordinates (y, ierr)
     call S2D%generate(r(1:2))
     call S3D%generate_from_axisymmetric_surface(S2D, N_sym, N_phi, N_theta)
  
     call field_line_loss_from_surface(S3D)
  enddo


  ! finalize
  if (firstP) then
     close (iu)
  endif
  contains
!.......................................................................

!.......................................................................
  subroutine field_line_loss_from_surface(S)
  type(t_flux_surface_3D), intent(in) :: S

  type(t_fieldline) :: F
  real(real64)      :: r0(3), y0(3), Lc, PsiNext
  integer :: i, j, idir, n, iloss(-1:1), inext(-1:1), inext1


  iloss = 0
  inext = 0
  n     = 0
  !PsiNext = 0.99158916136325959
  do i=0,S%n_phi-1
  do j=1,S%slice(i)%n_seg
     r0(1:2) = S%slice(i)%x(j,:)
     r0(3)   = S%slice(i)%phi
     n       = n + 1
     if (mod(n,nprs) .ne. mype) cycle
     call coord_trans (r0, CYLINDRICAL, y0, Trace_Coords)

     ! forward and backward tracing
     do idir=-1,1,2
        call F%init(y0, idir*Trace_Step, Trace_Method, Trace_Coords)
        Lc     = 0.d0
        inext1 = 0
        trace_loop: do
           call F%trace_1step()
           Lc = Lc + Trace_Step

           ! stop field line tracing at limit
           if (Lc > Limit) exit trace_loop

           ! check intersection with boundary
           if (F%intersect_boundary()) then
              iloss(idir) = iloss(idir) + 1
              exit trace_loop
           endif

!           ! connection to next (outward) main resonance
!           if (F%get_PsiN() >= PsiNext) then
!              inext1 = inext1 + 1
!           endif
        enddo trace_loop

        if (inext1 > 0) inext(idir) = inext(idir) + 1
     enddo
  enddo
  enddo

  call wait_pe()
  call sum_inte_data (iloss, 3)
  if (firstP) write (iu, *) 1.d0*iloss(-1)/n, 1.d0*iloss(1)/n
  end subroutine field_line_loss_from_surface
!.......................................................................
end subroutine field_line_loss
