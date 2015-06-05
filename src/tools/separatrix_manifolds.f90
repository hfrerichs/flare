!===============================================================================
! Generate invariant manifolds for hyperbolic fixed point
!===============================================================================
subroutine separatrix_manifolds
  use iso_fortran_env
  use run_control, only: N_sym, N_phi, N_psi, Label, Trace_Step, Grid_File
  use equilibrium
  use grid
  use math
  implicit none

  type(t_grid)       :: G
  character(len=120) :: Label0
  real(real64)       :: y(3), yh(3), Dphi, lambda1, lambda2, v1(2), v2(2), v(2)
  integer            :: iPx, orientation, i, idir


  iPx = N_psi
  call Xp(iPx)%analysis(lambda1, lambda2, v1, v2)


  orientation = 1
  if (Xp(iPx)%x(2) > 0.d0) orientation = -1
  Label0 = ''
  if (Label .ne. '') Label0 = '_'//trim(Label)


  ! setup sample point (toroidal discretization)
  call G%new(CYLINDRICAL, SEMI_STRUCTURED, 3, 1, N_phi)


  ! idir = -1: forward stable, backward unstable
  !         1: forward unstable
  v = v2
  Label = 'stable'//trim(Label0)
  do idir=-1,1,2
     do i=1,N_phi
        G%x(i,1:2) = Xp(iPx)%X + 0.1d0 * v
        G%x(i,  3) = -idir * Bt_sign * pi2/N_sym * (i-1) / N_phi
     enddo
     Grid_File = 'grid_'//trim(label)//'.dat'
     call G%store(Grid_File)

     Trace_Step = -1.d0 * Trace_Step
     call separatrix_manifolds_manual()

     v = v1
     Label = 'unstable'//trim(Label0)
  enddo

end subroutine separatrix_manifolds




!===============================================================================
! Generate invariant manifolds for hyperbolic fixed point
!
! Input (taken from run control file):
! Grid_File             toroidal discretization of X-point (+ small offset for tracing)
! Phi_Output            toroidal reference position [deg]
! N_sym                 toroidal symmetry of the configuration
! N_mult                multiplicity of grid point (in terms of 2 pi / N_sym)
!
! Trace_Step, Trace_Method
! Output_File
!===============================================================================
subroutine separatrix_manifolds_manual
  use iso_fortran_env
  use run_control, only: Grid_File, Label, Trace_Step, Trace_Method, Phi_Output, N_sym, N_mult
  use grid
  use fieldline
  use boundary
  use math
  use dataset
  use parallel
  implicit none

  integer, parameter :: iu = 54

  type(t_fieldline)  :: F
  type(t_grid)       :: G
  type(t_dataset)    :: D
  real(real64)       :: y(3), yh(3), Dphi
  integer            :: i, j, k, ind, ierr


  ! initialize
  if (firstP) then
     write (6, *) 'Generate invariant manifolds for hyperbolic fixed point'
     write (6, *)
  endif
  call G%load(Grid_File)
  write (6, *)
  

  call D%new(G%nodes() * N_mult, 3)
  open  (iu, file='strike_point_'//trim(Label)//'.dat')
  grid_loop: do i=1,G%nodes()
     y = G%node(i)
     write (6, *) y(3) / pi * 180.d0

     ! set initial location
     call F%init(y, Trace_Step, Trace_Method, CYLINDRICAL)

     ! initial trace to reference plane
     Dphi = abs(Phi_Output / 180.d0 * pi - y(3))
     call F%trace_Dphi(Dphi, .false., yh, ierr)
     D%x(i,:) = yh


     ! trace to symmetry planes
     Dphi = pi2 / N_sym
     do j=1,N_mult-1
        call F%trace_Dphi(Dphi, .true., yh, ierr)
        if (ierr > 0) then
           write (iu, *) yh
           exit
        endif

        ind = i + j*G%nodes()
        D%x(ind,:) = yh
     enddo
  enddo grid_loop
  close (iu)


  ! write RZ-cut of manifold
  open  (iu, file='manifold_'//trim(Label)//'.dat')
  do j=0,N_mult-1
     do i=1,G%nodes()
        ind = i + j*G%nodes()
        if (D%x(ind,1) > 0.d0) write (iu, *) D%x(ind,1:2)
     enddo
  enddo
  close (iu)

end subroutine separatrix_manifolds_manual
