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

  public :: ePsi_sub, H_eigenvectors


  contains
!=======================================================================



!=======================================================================
! Generate separatrix starting from X-point (Px).


! Sections 1 and 2 will be the right
! and left central (i.e. core) parts, respectively, while sections 3 and 4 will
! be the right and left divertor parts.
! Required input:
! iPx         =  Number of X-point (must have been set up in module equilibrium)
!                -> Xp(iPx)%X: Coordinates of X-point (R[cm], Z[cm])
! orientation =  1: lower null
!             = -1: upper null
! theta_cut   =  poloidal cut-off angle (-> split core separatrix into left and right segments)
! iconnect    =  number of X-point to which separatrix connects
!                (-> used to set poloidal cut-off angle)
!===============================================================================

  subroutine generate (this, iPx, theta_cut, C_cutl, C_cutr, iconnect)
  use math
  class(t_separatrix)                 :: this
  integer,       intent(in)           :: iPx
  real(real64),  intent(in), optional :: theta_cut
  type(t_curve), intent(in), optional :: C_cutl, C_cutr
  integer,       intent(in), optional :: iconnect

  real(real64) :: Px(2), H(2,2), v1(2), v2(2), x0(2), ds, ds0, theta_cutL, theta_cutR
  integer      :: orientation


  ! set orientation (lower null vs. upper null)
  orientation = 1; if (Xp(iPx)%X(2) > 0.d0) orientation = -1


  ! set cut-off poloidal angles
  theta_cutL = 0.d0; theta_cutR = 0.d0
  if (present(iconnect)) then
     if (Xp(iconnect)%undefined) then
        write (6, *) 'error: cannot connect to undefined X-point!'
        stop
     endif

     if (iconnect == iPx) then
        theta_cutR = Xp(iPx)%theta + pi
        theta_cutL = Xp(iPx)%theta + pi
     else
        theta_cutR = 0.5d0 * (Xp(iconnect)%theta + Xp(iPx)%theta)
        theta_cutL = theta_cutR + pi
     endif
  endif
  if (present(theta_cut)) then
     theta_cutL = theta_cut
     theta_cutR = theta_cut
  endif


  Px = Xp(iPx)%X ! Coordinates of X-point (R[cm], Z[cm])
  H  = Xp(iPx)%H ! Hessian matrix at X-point
  call H_eigenvectors(H, v1, v2)
  ds0     = 0.1d0
  v1      = ds0*v1; v2 = ds0*v2
  this%Px = Px
  ds      = ds0**2


  ! right core segment
  x0 = Px + v1*orientation + v2
  call this%M1%generate(x0,  1, ds, AltSurf=C_cutl, theta_cut=theta_cutR)
  this%M1%x(0,:) = Px

  ! left core segment
  x0 = Px - v1*orientation + v2
  call this%M2%generate(x0, -1, ds, AltSurf=C_cutr, theta_cut=theta_cutL)
  this%M2%x(this%M2%n_seg,:) = Px

  ! right divertor leg
  x0 = Px + v1*orientation - v2
  call this%M3%generate(x0, -1, ds, AltSurf=C_cutr)
  this%M3%x(this%M3%n_seg,:) = Px

  ! left divertor leg
  x0 = Px - v1*orientation - v2
  call this%M4%generate(x0,  1, ds, AltSurf=C_cutl)
  this%M4%x(0,:) = Px

  end subroutine generate
!=======================================================================



!=======================================================================
! calculate eigenvectors v1,v2 of Hessian matrix of pol. magn. flux at x
!=======================================================================
  subroutine H_eigenvectors (H, v1, v2)
  real(real64), intent(in)  :: H(2,2)
  real(real64), intent(out) :: v1(2), v2(2)

  real(real64) :: r(3), psi_xx, psi_xy, psi_yy, l1, l2, ac2, ac4, b2


  psi_xx = H(1,1) / (Psi_sepx-Psi_axis)
  psi_xy = H(1,2) / (Psi_sepx-Psi_axis)
  psi_yy = H(2,2) / (Psi_sepx-Psi_axis)


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
!=======================================================================



!=======================================================================
  subroutine plot(this, filename_prefix, parts)
  class(t_separatrix)          :: this
  character(len=*), intent(in) :: filename_prefix
  logical,          intent(in), optional :: parts

  integer, parameter :: iu = 99

  character(len=120) :: filename


  if (present(parts)  .and.  parts) then
     filename = filename_prefix//'_1.txt'
     call this%M1%plot(filename=filename)
     filename = filename_prefix//'_2.txt'
     call this%M2%plot(filename=filename)
     filename = filename_prefix//'_3.txt'
     call this%M3%plot(filename=filename)
     filename = filename_prefix//'_4.txt'
     call this%M4%plot(filename=filename)
  else
     filename = filename_prefix//'.plt'
     call this%M1%plot(filename=filename)
     call this%M2%plot(filename=filename, append=.true.)
     call this%M3%plot(filename=filename, append=.true.)
     call this%M4%plot(filename=filename, append=.true.)
  endif

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




end module separatrix
