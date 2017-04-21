!===============================================================================
! Generate magnetic axis file (required for non-axisymmetric equilibrium
!                              configurations)
!
! Input (taken from run control file):
!    Run_Mode           'automatic' or 'manual', see x_start
!
!    N_sym              Toroidal symmetry number, magnetic axis will be generated
!                       for 0 <= Phi <= 360/N_Sym deg
!
!    N_phi              Toroidal resolution (number of segments)
!
!    x_start            Reference position (cylindrical coordinates) from which
!                       the magnetic axis will be approximated
!                       Run_mode = 'automatic': x_start is initial guess and
!                                               will be optimized
!                                  'manual':    x_start remains fixed and must
!                                               be a good enough guess
!
! optional:
!    N_points           Number of sample points for each toroidal slice
!    Trace_Step
!    Trace_Method
!===============================================================================
subroutine generate_magnetic_axis
  use iso_fortran_env
  use run_control, only: Run_mode, x_start, N_sym, N_phi, N_mult, Output_File, &
                         N_points, Trace_Step, Trace_Method, Trace_Coords, &
                         AUTOMATIC, MANUAL
  use math
  use parallel
  use dataset
  use bfield
  use magnetic_axis
  implicit none

  integer, parameter :: iu = 54

  type(t_dataset)    :: D
  character(len=120) :: filename
  real(real64)       :: R, Z, Phi, Bf(3)
  integer            :: i


  Output_File = 'magnetic_axis.dat'
  if (firstP) then
     write (6, *) 'Generate magnetic axis, output in: ', adjustl(trim(Output_File))
     write (6, *)
  else
     return
  endif


  ! This is because poincare_plot appends slice index _i in output file for N_mult > 1 only!
  if (N_phi .le. 1) Then
     write (6,*) 'error: N_phi must be >= 2'
     write (6,*) 'N_phi =', N_phi
     stop
  endif

  ! select reference point for magnetic axis generation
  select case(Run_Mode)
  case(AUTOMATIC)
     call find_magnetic_axis_RZ_at_phi0(x_start)
  case(MANUAL)
  case default
     write (6, *) 'error: invalid Run_mode = ', adjustl(trim(Run_Mode)), '!'
     stop
  end select


  ! set default number of sample points for each toroidal slice
  N_points = 256

  ! Use cylindrical coordinates
  Trace_Coords = CYLINDRICAL

  ! setup toroidal field direction
  Bf      = get_Bf_cyl(x_start)
  Bt_sign = 1
  if (Bf(3) < 0.d0) Bt_sign = -1


  ! 1. Generate Poincare plot at N_mult toroidal slices
  write (6, *) '1st step:'
  N_mult = N_phi
  call poincare_plot


  ! 2. Combine output
  write (6, *)
  write (6, *) '2nd step:'
  write (6, *) 'Combine output'
  open  (iu, file=Output_File)
  write (iu, *) N_mult, N_sym
  do i=0,N_mult-1
     write (filename, '(i8)') i
     filename = 'magnetic_axis_'//trim(adjustl(filename))//'.dat'
     call D%load(filename, 2, output=SILENT)
     R   = sum(D%x(:,1)) / D%nrow
     Z   = sum(D%x(:,2)) / D%nrow
     Phi = 360.d0 / N_sym / N_mult * i
     write (iu, *) R, Z, Phi
  enddo
  close (iu)

end subroutine generate_magnetic_axis
!===============================================================================



!===============================================================================
! Iteratively adjust RZ position of magnetic axis at x0(3)
!===============================================================================
subroutine find_magnetic_axis_RZ_at_phi0(x0)
  use iso_fortran_env
  use run_control, only: N_sym
  use numerics,    only: Trace_Step, Trace_Method, Trace_Coords
  use math
  use fieldline
  implicit none

  real(real64), intent(inout) :: x0(3)

  type(t_fieldline)  :: F
  real(real64)       :: dR, dZ, yout(3), Dphi, dl, RZ0(3), xtol, xtol_want, ts
  integer            :: nr, nz, N_steps, Nrefine_max, is, Ntransit, imin(2), tc


  write (6, *) 'Running iterative approximation of magnetic axis RZ'
  write (6, *)

  ! Field line tracing parameters
  ts = Trace_Step
  ! explicitly set tracing coordinates to FL_ANGLE
  tc = FL_ANGLE;  if (tc .ne. Trace_Coords) ts = ts*pi/180.d0
  Dphi = 2.d0*pi/Real(N_sym,real64)
  N_steps = Nint(Abs(Dphi/ts))
  ts = Dphi/Real(N_steps,real64)
  Ntransit = 10 ! distance is average of these punctures from start

  ! Search parameters. Initial grid size and resolution, number of refinement steps
  nr = 31  ! Best to keep these odd numbers
  nz = 31
  Nrefine_max = 10
  dR = 10.d0
  dZ = 10.d0
  xtol_want = 1.d-5

  write(*,*) 'Initial guess',x0
  Write(*,*) 'Using max', Nrefine_max, ' refinement steps.'
  Write(*,*) '  Using next dr,dz',dr,dz
  Do is = 1,Nrefine_max
    RZ0 = search_grid()
    Write(*,*) '  Step ',is,', RZ axis= ',RZ0(1:2), ' dist = ',RZ0(3)
    xtol = RZ0(3)
    If (abs(xtol) < xtol_want) Exit
    dR = Max(2.d0*Abs(RZ0(1) - x0(1)),xtol_want)
    dZ = Max(2.d0*Abs(RZ0(2) - x0(2)),xtol_want)
    if (imin(1) .eq. nr) dR = dR*2.d0
    if (imin(2) .eq. nz) dZ = dZ*2.d0

    x0(1:2) = RZ0(1:2)
  Enddo
  If (is .eq. Nrefine_max + 1) Write(*,*) 'Reached max number of refinement steps'


  Contains
!.......................................................................
  Function search_grid () Result (RZdist)
  Real(real64) :: RZdist(3)
  integer      :: i, ir, iz, it
  real(real64), Allocatable :: dist(:,:)
  real(real64) :: x1(3)


  x1(3) = x0(3)
  Allocate(dist(nr,nz))
  dist(:,:) = 0.d0
  do ir = 1,nr
    do iz = 1,nz
      x1(1) = x0(1) - dR/2.d0 + Real(ir-1,real64)/Real(nr-1,real64)*dR
      x1(2) = x0(2) - dZ/2.d0 + Real(iz-1,real64)/Real(nz-1,real64)*dZ
      call F%init(x1, ts, Trace_Method, tc)
      Do it = 1,Ntransit
        Do i = 1,N_steps
          dl = F%Trace_1step()
        Enddo
        yout = F%rc
        dist(ir,iz) = dist(ir,iz) + sqrt( (yout(1)-x1(1))**2 + (yout(2)-x1(2))**2)
      Enddo
    enddo
  enddo
  dist = dist/Real(Ntransit,real64)
  imin = Minloc(dist)

  RZdist(1) = x0(1) - dR/2.d0 + Real(imin(1)-1,real64)/Real(nr-1,real64)*dR
  RZdist(2) = x0(2) - dZ/2.d0 + Real(imin(2)-1,real64)/Real(nz-1,real64)*dZ
  RZdist(3) = dist(imin(1),imin(2))
  Deallocate(dist)
  End Function search_grid
!.......................................................................
end subroutine find_magnetic_axis_RZ_at_phi0
!===============================================================================
