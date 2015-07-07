!===============================================================================
! Generete unperturbed flux surfaces
!===============================================================================
module flux_surface_2D
  use iso_fortran_env
  use curve2D
  implicit none
  private

  integer, parameter, public :: &
     RIGHT_HANDED =  1, &
     LEFT_HANDED  = -1

  type, extends(t_curve) :: t_flux_surface_2D
     real(real64) :: PsiN

     contains
     procedure :: generate => generate_flux_surface_2D
     procedure :: generate_closed
     procedure :: generate_open
     procedure :: setup_sampling
  end type t_flux_surface_2D

  public :: t_flux_surface_2D

  contains
!=======================================================================


!=======================================================================
! Generate axisymmetric flux surface at r = (R,Z [cm]) using a step size
! of Trace_Step and integration method Trace_Method.
!
! Default operation mode is tracing in both directions, but tracing in
! one direction only is set by the optional parameter direction.
!
! An alternate limiting surface can be given by the optional parameter AltSurf.
! An optional cut-off poloidal angle theta_cut can be given.
! Re-tracing of half-open surfaces is optional
!=======================================================================
  recursive subroutine generate_flux_surface_2D(this, r, direction, Trace_Step, N_steps, &
      Trace_Method, AltSurf, theta_cut, retrace, sampling)
  use equilibrium
  use ode_solver
  use boundary
  use math
  class(t_flux_surface_2D) :: this
  real(real64), intent(in) :: r(2)
  
  integer, intent(in), optional       :: direction, N_steps, Trace_Method
  real(real64), intent(in), optional  :: Trace_Step, theta_cut
  type(t_curve), intent(in), optional :: AltSurf
  logical, intent(in), optional       :: retrace
  integer, intent(in), optional       :: sampling

  type(t_ODE) :: F
  real*8, dimension(:,:), allocatable :: tmp
  real*8  :: yl(3), yc(3), thetal, thetac, dtheta, X(3), ds, r3(3)
  integer :: idir, i, nmax, imethod, id, n(-1:1)
  logical :: retrace_from_boundary = .false.


  ! determine trace step
  if (present(Trace_Step)) then
     ds = Trace_Step
  else
     ! TODO: get PsiN from equilibrium and use small trace step close to the separatrix
     ds = length_scale() / 800.d0
  endif


  ! set trace method
  if (present(Trace_Method)) then
     imethod = Trace_Method
  else
     imethod = NM_AdamsBashforth4
  endif


  ! allocate temporary array for output
  if (present(N_steps)) then
     nmax = N_steps
  else
     nmax = N_steps_guess(ds)
  endif
  allocate (tmp(-nmax:nmax,2))


  ! trace direction
  n         = nmax
  if (present(direction)) then
     if (abs(direction) .ne. 1) then
        write (6, *) 'error: invalid parameter direction ', direction
     endif
     n(-direction) = 0
  endif


  ! re-trace from boundary if flux surface is half-open?
  if (present(retrace)) retrace_from_boundary = retrace


  ! initialize variables
  yl(3)     = 0.d0
  yc(3)     = 0.d0
  r3(1:2)   = r
  r3(3)     = 0.d0
  this%PsiN = get_PsiN(r3)


  ! trace in forward and backward direction
  do idir=-1,1,2
     if (n(idir) == 0) cycle
     ! set initial position
     call F%init_ODE(2, r, idir*ds, Bpol_sub, imethod)
     yl(1:2)  = r
     thetal   = get_poloidal_angle(r3)
     tmp(0,:) = F%yc

     do i=1,nmax
        ! progress one step
        yc(1:2)       = F%next_step()
        thetac        = get_poloidal_angle(yc)

        ! save current position
        tmp(idir*i,:) = yc(1:2)

        ! check intersection user defined surface (AltSurf)
        if (present(AltSurf)) then
           if (intersect_curve(yl(1:2), yc(1:2), AltSurf, X(1:2))) then
              tmp(idir*i,:) = X(1:2)
              n(idir)       = i
              exit
           endif

        ! check intersection with boundary
        else
           if (intersect_boundary(yl, yc, X, id)) then
              tmp(idir*i,:) = X(1:2)
              n(idir)       = i
              exit
           endif
        endif


        ! cross cut-off poloidal angle?
        if (present(theta_cut)) then
           dtheta = thetac - thetal
           if (abs(dtheta) > pi) dtheta = dtheta - sign(pi2,dtheta)
           if ((thetal+dtheta-theta_cut)*(thetal-theta_cut) < 0.d0) then
              tmp(idir*i,:) = tmp(idir*(i-1),:) + abs((thetal-theta_cut)/dtheta)*(yc(1:2)-tmp(idir*(i-1),:))
              n(idir)       = i
              exit
           endif
        endif

        ! prepare next step
        yl     = yc
        thetal = thetac
     enddo
  enddo


! save data
  ! either closed flux surface, or flux surface is limited on both side
  if (((n(-1)  < nmax  .and.  n(1)  < nmax)  .or.  &
      (n(-1) == nmax  .and.  n(1) == nmax)) .or. .not.retrace_from_boundary) then
     call this%new(n(-1) + n(1))
     this%x = tmp(-n(-1):n(1),:)
     deallocate (tmp)

  ! flux surface is limited on one side only
  else
     if (n(-1) < nmax) r3(1:2) = tmp(-n(-1),:) + 0.01d0*(tmp(-n(-1),:) - tmp(-n(-1)+1,:))
     if (n( 1) < nmax) r3(1:2) = tmp( n( 1),:) + 0.01d0*(tmp( n( 1),:) - tmp( n( 1)-1,:))
     deallocate (tmp)
     call this%generate (r3(1:2), direction, Trace_Step, N_steps, Trace_Method, AltSurf, theta_cut)
  endif


  ! setup sampling array
  if (present(sampling)) then
     select case(sampling)
     case(ANGLE)
        r3 = get_magnetic_axis(0.d0)
        call this%setup_angular_sampling(r3(1:2))
     case(DISTANCE)
        call this%setup_length_sampling()
     case default
        write (6, *) 'warning: invalid option for sampling in t_flux_surface%generate!'
     end select
  endif

  end subroutine generate_flux_surface_2D
!=======================================================================



!=======================================================================
  subroutine generate_closed(this, r, direction, Trace_Step, N_steps, Trace_Method, AltSurf, retrace)
  use equilibrium
  class(t_flux_surface_2D) :: this
  real(real64), intent(in) :: r(2)
  integer, intent(in)      :: direction

  integer, intent(in), optional       :: N_steps, Trace_Method
  real(real64), intent(in), optional  :: Trace_Step
  type(t_curve), intent(in), optional :: AltSurf
  logical, intent(in), optional       :: retrace

  real(real64) :: theta_cut, r3(3), x1(2), x2(2), dl, dl0


  r3(1:2) = r
  r3(3)   = 0.d0
  theta_cut = get_poloidal_angle(r3)
  call this%generate(r, direction, Trace_Step, N_steps, Trace_Method, AltSurf, theta_cut, retrace)

  ! retrieve step size
  x1  = this%x(0,:)
  x2  = this%x(1,:)
  dl0 = sqrt(sum((x1-x2)**2))

  ! close flux surface
  x1 = this%x(0,:)
  x2 = this%x(this%n_seg,:)
  dl = sqrt(sum((x1-x2)**2))
  this%closed          = .true.
  this%x(this%n_seg,:) = x1
  if (dl > 1.d-2*dl0) then
     write (6, *) 'warning: deviation in closed flux surface > 0.01 * Trace_Step!'
  endif

  end subroutine generate_closed
!=======================================================================



!=======================================================================
  subroutine generate_open(this, r, cut_fw, cut_bw, extend_fw, extend_bw)
  class(t_flux_surface_2D)            :: this
  real(real64), intent(in)            :: r(2)
  type(t_curve), intent(in), optional :: cut_fw, cut_bw
  real(real64), intent(in), optional  :: extend_fw, extend_bw

  type(t_flux_surface_2D) :: f_fw, f_bw


  call f_fw%generate(r,  1, AltSurf=cut_fw)
  call f_bw%generate(r, -1, AltSurf=cut_bw)
  this%t_curve = connect(f_bw%t_curve, f_fw%t_curve)
  call this%setup_length_sampling()

  end subroutine generate_open
!=======================================================================



!=======================================================================
! interface subroutine for poloidal magnetic field
!=======================================================================
  subroutine Bpol_sub (n, t, y, f)
  use equilibrium

  integer, intent(in) :: n
  real*8,  intent(in) :: t, y(n)
  real*8, intent(out) :: f(n)

  real*8 :: Bf(3), Bpol, y3(3)


  y3(1:2) = y
  y3(3)   = 0.d0
  Bf      = get_Bf_eq2D(y3)
  Bpol    = sqrt(Bf(1)**2 + Bf(2)**2)
  f       = - Bf(1:2)/Bpol / Ip_sign

  end subroutine Bpol_sub
!=======================================================================



!===============================================================================
! N_STEPS_GUESS
!
! Guess the required number of trace steps (for flux surface generation),
! based on the step size.
!===============================================================================
  function N_steps_guess (Trace_Step) result(N)
  use equilibrium
  use math
  real*8, intent(in) :: Trace_Step
  integer :: N

  real*8  :: dl, r3(3)


  r3 = get_magnetic_axis(0.d0)
  dl = pi2 * r3(1) / dabs(Trace_Step)
  N  = 2 * int(dl)

  end function N_steps_guess
!===============================================================================



!===============================================================================
! Prepare flux surface for sampling using length weights for the parts near the
! X-points while angular weights are used for the intermediate part.
! The line X-point (x1,x2) to magnetic axis (xc) is used as reference. The transition
! between angle-weighted sampling and length-weighted sampling occurs at an
! angle of DthetaR (lower end) and DthetaL (upper end)
! Dtheta0: reference angular weight (pi2 for full loop)
!===============================================================================
  subroutine setup_sampling(this, x1, x2, xc, DthetaR, DthetaL)
  use math
  class(t_flux_surface_2D) :: this
  real(real64), intent(in) :: x1(2), x2(2), xc(2), DthetaR, DthetaL

  real(real64) :: x(2), phi, dphi, dphi1, dphi2, phi1, phi2, s, f0, g0, w, Dtheta0
  integer      :: i, n, iseg1, iseg2


  ! set reference weight automatically
  phi1  = datan2(x1(2)-xc(2), x1(1)-xc(1))
  phi2  = datan2(x2(2)-xc(2), x2(1)-xc(1))
  Dtheta0 = phi2 - phi1; if (Dtheta0 <= 0.d0) Dtheta0 = Dtheta0 + pi2


  !.....................................................................
  ! 0.a small segments
  if (Dtheta0 < DthetaR+DthetaL) then
     call this%setup_length_sampling()
     return
  endif

  ! 0.b use all angular weights
  if (DthetaR.eq.0.0d0 .and. DthetaL.eq.0.d0) then
     call this%setup_angular_sampling(xc)
     return
  endif

  ! otherwise setup sampling by segment lengths
  n  = this%n_seg
  call this%setup_length_sampling(raw_weights=.true.)
  !.....................................................................


  !.....................................................................
  ! 1. find first and last segment number for angle-weighted domain

  ! find first segment on L
  do i=1,n
     ! upper node of segment
     x     = this%x(i,:)
     phi   = datan2(x(2)-xc(2), x(1)-xc(1))
     dphi1 = phi - phi1
     iseg1 = i
     if (dabs(dphi1) .gt. DthetaR) exit
  enddo

  ! find last segment on L
  do i=n,1,-1
     ! lower node of segment
     x     = this%x(i-1,:)
     phi   = datan2(x(2)-xc(2), x(1)-xc(1))
     dphi2 = phi - phi2
     iseg2 = i
     if (dabs(dphi2) .gt. DthetaL) exit
  enddo
  !.....................................................................


  !.....................................................................
  ! 2. set new (angle) weights in middle part
  x    = this%x(iseg1,:)
  phi1 = datan2(x(2)-xc(2), x(1)-xc(1))
  w    = 0.d0
  do i=iseg1+1,iseg2-1
     x    = this%x(i,:)
     phi2 = datan2(x(2)-xc(2), x(1)-xc(1))

     dphi = phi2 - phi1
     if (dabs(dphi) .gt. pi) dphi = dphi - dsign(pi2,dphi)

     this%w(i) = dphi / Dtheta0
     w         = w + dphi / Dtheta0
     phi1      = phi2
  enddo
  !.....................................................................


  !.....................................................................
  ! 3. update length weights in first and last parts

  ! first part
  s     = sum(this%w(1:iseg1))
  dphi1 = dabs(dphi1) / Dtheta0
  this%w(1:iseg1) = this%w(1:iseg1) / s * dphi1

  !  last part
  s     = sum(this%w(iseg2:n))
  dphi2 = dabs(dphi2) / Dtheta0
  this%w(iseg2:n) = this%w(iseg2:n) / s * dphi2
  !.....................................................................


  !.....................................................................
  ! 4. check normalization
  w = sum(this%w)
  if (w < 1.d0 - 1.d-7) then
     write (6, *) 'error in t_flux_surface_2D%setup_sampling: unexpected sum of weights!'
     write (6, *) 'parameters are:'
     write (6, *) 'x1    = ', x1
     write (6, *) 'x2    = ', x2
     write (6, *) 'xc    = ', xc
     write (6, *) 'DthetaR = ', DthetaR
     write (6, *) 'DthetaL = ', DthetaL
     write (6, *) 'Dtheta0 = ', Dtheta0
     call this%plot(filename='error.plt')
     stop
  else
     do i=1,n
        this%w(i) = this%w(i-1) + this%w(i)
     enddo
     this%w = this%w / this%w(n)
  endif
  !.....................................................................

  end subroutine setup_sampling
!===============================================================================


end module flux_surface_2D
