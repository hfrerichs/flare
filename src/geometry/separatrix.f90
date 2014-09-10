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
  ds      = ds0**2 * Ip_sign * orientation


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
  use equilibrium

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
  subroutine generate_gradPsiN_path(this, Px, L)
  use ode_solver
  use run_control, only: Trace_Method, N_steps
  class(t_gradPsiN_path)   :: this
  real(real64), intent(in) :: Px(2), L

  type(t_ODE)  :: Path
  real(real64) :: v1(2), v2(2), x0(2), ds, dl, t
  integer      :: n_seg, is


  call H_eigenvectors(Px, v1, v2)
  ! offset from X-point for tracing
  dl    = Px(1) / 1.d2
  x0    = Px - dl*v1

  ! step size
  ds    = 0.1d0 * dl

  ! expected number of segments
  n_seg = nint((L-dl) / ds) + 2
  call this%new(n_seg)
  this%x_data(0,:) = Px
  this%x_data(1,:) = x0


  ! generate grad PsiN path
  call Path%init_ODE(2, x0, ds, ePsi_sub, Trace_Method)
  do is=2,n_seg
     this%x_data(is,:) = Path%next_step()
     dl                = dl + ds
  enddo


  ! adjust last node to match L
  t = (L-dl+ds)/ds
  this%x_data(n_seg,:) = (1.d0-t)*this%x_data(n_seg-1,:) + t*this%x_data(n_seg,:)

  end subroutine generate_gradPsiN_path
!=======================================================================


end module separatrix
