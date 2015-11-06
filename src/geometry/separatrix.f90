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
     procedure :: generate
     procedure :: plot
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
!
! Optional input:
! theta_cut   =  poloidal cut-off angle (-> split core separatrix into left and right segments)
! C_cutl,
! C_cutr      =  cut-off surface for separatrix segments
! iconnect    =  number of X-point to which separatrix connects
!                (-> used to set poloidal cut-off angle)
! offset      =  offset from X-point to start tracing
! trace_step  =  integration step size
!===============================================================================

  subroutine generate (this, iPx, theta_cut, C_cutl, C_cutr, iconnect, offset, trace_step)
  use math
  use run_control, only: Debug
  class(t_separatrix)                 :: this
  integer,       intent(in)           :: iPx
  real(real64),  intent(in), optional :: theta_cut, offset, trace_step
  type(t_curve), intent(in), optional :: C_cutl, C_cutr
  integer,       intent(in), optional :: iconnect

  real(real64) :: Px(2), H(2,2), v1(2), v2(2), x0(2), ds, s1, &
                  theta_cutL, theta_cutR, dtheta
  integer      :: orientation


  ! set orientation (lower null vs. upper null)
  orientation = 1
  if (Xp(iPx)%X(2) > 0.d0) then
     orientation = -1
  endif


  ! set poloidal cut-off angles
  theta_cutL = 0.d0; theta_cutR = 0.d0
  if (present(iconnect)) then
     if (Xp(abs(iconnect))%undefined) then
        write (6, *) 'error: cannot connect to undefined X-point!'
        stop
     endif

     ! connect to itself (cut-off at +/- pi from X-point)
     if (iconnect == iPx) then
        theta_cutR = Xp(iPx)%theta + pi; if (theta_cutR > pi) theta_cutR = theta_cutR - pi2
        theta_cutL = theta_cutR

     ! direct connection to another X-point
     ! cut-off halfway between X-points
     elseif (iconnect > 0) then
        dtheta = Xp(iconnect)%theta - Xp(iPx)%theta; if (dtheta < 0.d0) dtheta = dtheta + pi2
        if (Debug) then
           write (6, *) 'theta1 = ', Xp(iPx)%theta/pi*180.d0
           write (6, *) 'theta2 = ', Xp(iconnect)%theta/pi*180.d0
           write (6, *) 'dtheta = ', dtheta/pi*180.d0
           write (6, *)
        endif

        theta_cutR = Xp(iPx)%theta + 0.5d0*dtheta
        theta_cutL = theta_cutR + pi; if (theta_cutL > pi) theta_cutL = theta_cutL - pi2

     ! use position of another X-point as reference
     elseif (iconnect < 0) then
        theta_cutR = Xp(abs(iconnect))%theta
        theta_cutL = Xp(abs(iconnect))%theta

     ! unset cut-off angle
     else
        theta_cutR = -pi2
        theta_cutL = -pi2
     endif
  endif
  if (present(theta_cut)) then
     theta_cutL = theta_cut
     theta_cutR = theta_cut
  endif
  if (Debug) then
     write (6, *) 'theta_cut = ', theta_cutL/pi*180.d0, theta_cutR/pi*180.d0
  endif


  ! set offset from X-point
  s1 = 0.1d0
  if (present(offset)) then
     s1 = offset
  endif


  ! set trace step
  ds      = 1.d-2
  if (present(trace_step)) then
     ds   = trace_step
  endif


  Px = Xp(iPx)%X ! Coordinates of X-point (R[cm], Z[cm])
  H  = Xp(iPx)%H ! Hessian matrix at X-point
  call H_eigenvectors(H, v1, v2)
  v1      = s1*v1; v2 = s1*v2
  this%Px = Px
  if (Debug) then
     open  (97, file='right.tmp', position='append')
     write (97, *) Px
     write (97, *) Px + v1*orientation + v2
     close (97)
     open  (96, file='left.tmp', position='append')
     write (96, *) Px
     write (96, *) Px - v1*orientation + v2
     close (96)
     open  (95, file='v1.tmp', position='append')
     write (95, *) Px
     write (95, *) Px + v1
     write (95, *)
     close (95)
     open  (94, file='v2.tmp', position='append')
     write (94, *) Px
     write (94, *) Px + v2
     write (94, *)
     close (94)
  endif


  ! right core segment
  x0 = Px + v1*orientation + v2*orientation
  call this%M1%generate(x0,  1, ds, AltSurf=C_cutl, theta_cut=theta_cutR)
  this%M1%x(0,:) = Px

  ! left core segment
  x0 = Px - v1*orientation + v2*orientation
  call this%M2%generate(x0, -1, ds, AltSurf=C_cutr, theta_cut=theta_cutL)
  this%M2%x(this%M2%n_seg,:) = Px

  ! right divertor leg
  x0 = Px + v1*orientation - v2*orientation
  call this%M3%generate(x0, -1, ds, AltSurf=C_cutr)
  this%M3%x(this%M3%n_seg,:) = Px

  ! left divertor leg
  x0 = Px - v1*orientation - v2*orientation
  call this%M4%generate(x0,  1, ds, AltSurf=C_cutl)
  this%M4%x(0,:) = Px

  end subroutine generate
!=======================================================================



!=======================================================================
  subroutine plot(this, filename_prefix, parts)
  class(t_separatrix)          :: this
  character(len=*), intent(in) :: filename_prefix
  logical,          intent(in), optional :: parts

  integer, parameter :: iu = 50

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
