!===============================================================================
! Generate magnetic separatrix (axisymmetric version)
!===============================================================================
module separatrix
  use equilibrium
  use curve2D
  use flux_surface_2D
  implicit none

  private
  type, public :: t_separatrix
  !type, extends(t_curve), public :: t_separatrix
     type(t_flux_surface_2D) :: M1, M2, M3, M4
     real(real64) :: Px(2)
     contains
     procedure :: generate, plot
  end type t_separatrix

  type, extends(t_curve), public :: t_gradPsiN_path
     contains
     procedure :: generate => generate_gradPsiN_path
  end type t_gradPsiN_path


  contains
!=======================================================================



!=======================================================================
! Generate separatrix starting from X-point (Px).


! Sections 1 and 2 will be the right
! and left central (i.e. core) parts, respectively, while sections 3 and 4 will
! be the right and left divertor parts.
! Required input:
! Px          =  Coordinates of X-point (R[cm], Z[cm])
! orientation =  1: lower null
!             = -1: upper null
! theta_cut   =  poloidal cut-off angle (-> split core separatrix into left and right segments)
!===============================================================================

  subroutine generate (this, Px, orientation, theta_cut, AltSurf)
  class(t_separatrix)       :: this
  real(real64), intent(in)  :: Px(2)
  integer,      intent(in)  :: orientation
  real(real64), intent(in)  :: theta_cut
  type(t_curve), intent(in), optional :: AltSurf

  real(real64) :: v1(2), v2(2), x0(2), ds, ds0


  call H_eigenvectors(Px, v1, v2)
  ds0     = 0.1d0
  v1      = ds0*v1; v2 = ds0*v2
  this%Px = Px
  ds      = ds0**2 * orientation


  ! right core segment
  x0 = Px + v1 + v2
  call this%M1%generate(x0, -1, ds, AltSurf=AltSurf, theta_cut=theta_cut)

  ! left core segment
  x0 = Px - v1 + v2
  call this%M2%generate(x0,  1, ds, AltSurf=AltSurf, theta_cut=theta_cut)

  ! right divertor leg
  x0 = Px + v1 - v2
  call this%M3%generate(x0,  1, ds, AltSurf=AltSurf)

  ! left divertor leg
  x0 = Px - v1 - v2
  call this%M4%generate(x0, -1, ds, AltSurf=AltSurf)

  end subroutine generate
!-----------------------------------------------------------------------
! calculate eigenvectors v1,v2 of Hessian matrix of pol. magn. flux at x
!-----------------------------------------------------------------------
  subroutine H_eigenvectors (x, v1, v2)
  real(real64), intent(in)  :: x(2)
  real(real64), intent(out) :: v1(2), v2(2)

  real(real64) :: r(3), psi_xx, psi_xy, psi_yy, l1, l2, ac2, ac4, b2


  ! evaluate Hessian at x
  r(1:2) = x
  r(3)   = 0.d0
  psi_xx = get_DPsiN(r, 2, 0)
  psi_xy = get_DPsiN(r, 1, 1)
  psi_yy = get_DPsiN(r, 0, 2)


  ! get eigenvalues l1,l2 of Hessian at X-point
  ac2 = 0.5d0  * (psi_xx + psi_yy)
  ac4 = 0.25d0 * (psi_xx - psi_yy)**2
  b2  = psi_xy**2
  l1  = ac2 + dsqrt(ac4 + b2)
  l2  = ac2 - dsqrt(ac4 + b2)


  ! construct normalized eigenvectors
  ! ISSUE: this might not work if the X-point is straight below the magnetic axis!
  v1(1) = 1.d0
  v1(2) = - (psi_xx - l1) / psi_xy
  v1    = v1 / sqrt(sum(v1**2))

  v2(1) = 1.d0
  v2(2) = - (psi_xx - l2) / psi_xy
  v2    = v2 / sqrt(sum(v2**2))

  end subroutine H_eigenvectors
!-----------------------------------------------------------------------
!=======================================================================



!=======================================================================
  subroutine plot(this, filename_prefix)
  class(t_separatrix)          :: this
  character(len=*), intent(in) :: filename_prefix

  integer, parameter :: iu = 99

  character(len=120) :: filename


  filename = filename_prefix//'_1.txt'
  call this%M1%plot(filename=filename)
  filename = filename_prefix//'_2.txt'
  call this%M2%plot(filename=filename)
  filename = filename_prefix//'_3.txt'
  call this%M3%plot(filename=filename)
  filename = filename_prefix//'_4.txt'
  call this%M4%plot(filename=filename)
  end subroutine plot
!=======================================================================



!=======================================================================
! interface subroutine to normalized grad Psi vector for ODE solver
!=======================================================================
  subroutine ePsi_sub (n, t, y, f)
  integer, intent(in)       :: n
  real(real64), intent(in)  :: t, y(n)
  real(real64), intent(out) :: f(n)

  real(real64) :: r(3), DPsiDR, DPsiDZ

  if (n .ne. 2) then
     write (6, *) 'error in subroutine ePsi_sub: n <> 2!'
     stop
  endif

  r(1:2) = y
  r(3)   = 0.d0
  f(1)   = get_DPsiN(r, 1, 0)
  f(2)   = get_DPsiN(r, 0, 1)
  f      = f / sqrt(sum(f**2))

  end subroutine ePsi_sub
!=======================================================================



!=======================================================================
! Generate path along grad-Psi from Px (X-point)
! orientation = 1: ascent PsiN in left SOL
!             = 2: ascent PsiN in right SOL
!             = 3: descent PsiN to core
!             = 4: descent PsiN to PFR
!=======================================================================
  subroutine generate_gradPsiN_path(this, Px, orientation, L, PsiN)
  use ode_solver
  use run_control, only: Trace_Method, N_steps
  class(t_gradPsiN_path)   :: this
  real(real64), intent(in) :: Px(2)
  integer,      intent(in) :: orientation
  real(real64), intent(in), optional :: L, PsiN

  type(t_ODE)  :: Path
  real(real64) :: v1(2), v2(2), x0(2), y(3), ds, dl, t, Psi0
  integer      :: n_seg, is


  ! 0. initialize
  call H_eigenvectors(Px, v1, v2)
  ! offset from X-point for tracing
  dl    = Px(1) / 1.d2
  ! step size
  ds    = 0.1d0 * dl


  ! 1. select orientation from saddle point (X-point)
  select case(orientation)
  case(1)
     x0    = Px - dl*v1
  case(2)
     x0    = Px + dl*v1
  case(3)
     x0    = Px + dl*v2
     ds    = -ds
  case(4)
     x0    = Px - dl*v2
     ds    = -ds
  case default
     write (6, *) 'error in subroutine t_gradPsiN_path%generate:'
     write (6, *) 'invalid parameter orientation = ', orientation
     stop
  end select


  ! 2.1 determine length/number of segments by final PsiN value
  if (present(PsiN)) then
     n_seg  = 1
     y(1:2) = x0
     y(3)   = 0.d0
     Psi0   = get_PsiN(y)
     call Path%init_ODE(2, x0, ds, ePsi_sub, Trace_Method)
     do
        if (Psi0 < PsiN  .and.  get_PsiN(y) > PsiN) exit
        if (Psi0 > PsiN  .and.  get_PsiN(y) < PsiN) exit
        y(1:2) = Path%next_step()
        n_seg  = n_seg + 1
     enddo

  ! 2.2 expected number of segments
  elseif (present(L)) then
     n_seg = nint((L-dl) / abs(ds)) + 2

  ! either PsiN or L must be given!
  else
     write (6, *) 'error in subroutine t_gradPsiN_Path%generate:'
     write (6, *) 'either parameter L or PsiN must be given!'
     stop
  endif
  call this%new(n_seg)
  this%x(0,:) = Px
  this%x(1,:) = x0


  ! 3. generate grad PsiN path
  call Path%init_ODE(2, x0, ds, ePsi_sub, Trace_Method)
  do is=2,n_seg
     this%x(is,:) = Path%next_step()
     dl           = dl + abs(ds)
  enddo


  ! 4. adjust last node to match L
  t = 1.d0
  if (present(L)) t = (L-dl+abs(ds))/abs(ds)
  this%x(n_seg,:) = (1.d0-t)*this%x(n_seg-1,:) + t*this%x(n_seg,:)

  end subroutine generate_gradPsiN_path
!=======================================================================


end module separatrix
