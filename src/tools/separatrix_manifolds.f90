!===============================================================================
! Generate invariant manifolds for hyperbolic fixed point
!===============================================================================
subroutine separatrix_manifolds
  use iso_fortran_env
  use run_control, only: N_sym, N_phi, N_psi, Label, Trace_Step, Grid_File
  use magnetic_axis
  use equilibrium, only: Xp
  use grid
  use math
  use parallel
  implicit none

  type(t_grid)       :: G
  character(len=120) :: Label0, Label_tmp(-1:1), Grid_File_tmp(-1:1)
  real(real64)       :: y(3), yh(3), Dphi, lambda1, lambda2, v1(2), v2(2), v(2)
  integer            :: iPx, orientation, i, idir, ierr


  ! generate toroidal discretization based on X-point
  if (firstP) then
     iPx = N_psi
     call Xp(iPx)%analysis(lambda1, lambda2, v1, v2, ierr)


     orientation = 1;  if (Xp(iPx)%x(2) > 0.d0) orientation = -1
     Label0 = '';      if (Label .ne. '') Label0 = '_'//trim(Label)


     ! initialize grid
     call G%new(CYLINDRICAL, SEMI_STRUCTURED, 3, 1, N_phi)


     ! idir = -1: forward stable, backward unstable
     !         1: forward unstable
     v = v2
     Label = 'stable'//trim(Label0)
     do idir=-1,1,2
        do i=1,N_phi
           G%x(i,1:2) = Xp(iPx)%X + 0.1d0 * v * orientation
           G%x(i,  3) = -idir * Bt_sign * pi2/N_sym * (i-1) / N_phi
        enddo
        Label_tmp(idir)     = Label
        Grid_File_tmp(idir) = 'grid_'//trim(Label)//'.dat'
        call G%store(Grid_File_tmp(idir))

        v = v1
        Label = 'unstable'//trim(Label0)
     enddo
  endif


  ! run main computation
  do idir=-1,1,2
     Label     = Label_tmp(idir)
     Grid_File = Grid_File_tmp(idir)

     Trace_Step = -1.d0 * Trace_Step
     call separatrix_manifolds_manual()
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
  use run_control, only: Grid_File, Label, Trace_Step, Trace_Method, Phi_Output, N_sym, N_mult, Output_Format
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
  type(t_dataset)    :: D, S
  real(real64)       :: y(3), yh(3), Dphi, x1(2), l, t1, t2, t3
  integer            :: i, j, k, ind, ierr


  ! initialize
  if (firstP) then
     write (6, *) 'Generate invariant manifolds for hyperbolic fixed point'
     write (6, *)
  endif
  call G%load(Grid_File)
  call D%new(G%nodes() * N_mult, 3)
  call S%new(G%nodes(), 4)
  write (6, *)


  ! begin main computation (parallel)
  call wait_pe()
  call cpu_time(t1)
  grid_loop: do i=1+mype,G%nodes(),nprs
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
           S%x(i,1:3) = yh
           S%x(i,  4) = ierr
           exit
        endif

        ind = i + j*G%nodes()
        D%x(ind,:) = yh
     enddo
  enddo grid_loop
  ! finished main computation (parallel)
  call wait_pe()
  call D%mpi_allreduce()
  call S%mpi_allreduce()
  call cpu_time(t2)
  if (firstP .and. nprs>1) write (6, *) 'time for main computation (parallel) [s]: ', t2 - t1


  ! post-processing (serial)
  if (mype > 0) return

  ! write RZ-cut of manifold
  l = huge(1.d0)
  open  (iu, file='manifold_'//trim(Label)//'.dat')
  do j=0,N_mult-1
     do i=1,G%nodes()
        ind = i + j*G%nodes()
        if (D%x(ind,1) > 0.d0) then
           select case(Output_Format)
           case(1,2)
              write (iu, *) D%x(ind,1:2)
           case(3)
              if (l < huge(1.d0)) l = sqrt(sum((D%x(ind,1:2)-x1)**2))
              if (l > 0.1d0) then
                 x1 = D%x(ind,1:2)
                 write (iu, *) x1, ind
                 l  = 0.d0
              endif
           end select
        else
           if (Output_Format == 2) write (iu, *) 'NaN'
        endif
     enddo
  enddo
  close (iu)

  ! write strike points
  open  (iu, file='strike_point_'//trim(Label)//'.dat')
  do i=1,G%nodes()
     if (S%x(i,4) > 0.d0) write (iu, *) S%x(i,1:3)
  enddo
  close (iu)
  call cpu_time(t3)
  if (nprs > 1) write (6, *) 'time for post-processing (serial) [s]: ', t3 - t2

end subroutine separatrix_manifolds_manual
