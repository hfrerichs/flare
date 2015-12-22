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
     LEFT_HANDED  = -1, &
     FORWARD      =  1, &
     BACKWARD     = -1, &
     CCW          =  1001, &
     CW           = -1001

  type, extends(t_curve) :: t_flux_surface_2D
     real(real64) :: PsiN

     contains
     procedure :: generate ! to be replaced by generate_branch
     procedure :: generate_branch
     procedure :: generate_closed
     procedure :: generate_open
     procedure :: setup_sampling
     procedure :: surface
     procedure :: surface_analysis
     procedure :: volume
     procedure :: broadcast
  end type t_flux_surface_2D

  public :: t_flux_surface_2D
  public :: adaptive_step_size

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
  recursive subroutine generate(this, r, direction, Trace_Step, N_steps, &
      Trace_Method, AltSurf, theta_cut, retrace, sampling, AltSurf_bw, AltSurf_fw)
  use magnetic_axis
  use equilibrium, only: length_scale, get_PsiN, get_poloidal_angle
  use ode_solver
  use boundary
  use math
  class(t_flux_surface_2D) :: this
  real(real64), intent(in) :: r(2)
  
  integer, intent(in), optional       :: direction, N_steps, Trace_Method
  real(real64), intent(in), optional  :: Trace_Step, theta_cut
  type(t_curve), target, intent(in), optional :: AltSurf, AltSurf_bw, AltSurf_fw
  logical, intent(in), optional       :: retrace
  integer, intent(in), optional       :: sampling

  type(t_curve), pointer :: C_boundary
  type(t_ODE) :: F
  real*8, dimension(:,:), allocatable :: tmp
  real*8  :: yl(3), yc(3), thetal, thetac, dtheta, X(3), ds, r3(3)
  integer :: idir, i, nmax, imethod, id, n(-1:1)
  logical :: retrace_from_boundary = .false.


  ! determine trace step
  if (present(Trace_Step)) then
     ds = Trace_Step
  else
     ! set small trace step close to the separatrix
     ds = length_scale() * min(1.d-3, abs(1.d0-get_PsiN(r)))
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


  ! set up user defined boundary
  nullify(C_boundary)
  if (present(AltSurf)) C_boundary => AltSurf
  if (present(AltSurf_bw)) C_boundary => AltSurf_bw


  ! re-trace from boundary if flux surface is half-open?
  if (present(retrace)) retrace_from_boundary = retrace


  ! initialize variables
  yl(3)     = 0.d0
  yc(3)     = 0.d0
  r3(1:2)   = r
  r3(3)     = 0.d0
  call this%destroy()
  this%PsiN = get_PsiN(r3)


  ! trace in forward and backward direction
  do idir=-1,1,2
     if (idir == 1  .and.  present(AltSurf_fw)) C_boundary => AltSurf_fw
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

        ! check intersection with user defined surface (AltSurf)
        if (associated(C_boundary)) then
           if (intersect_curve(yl(1:2), yc(1:2), C_boundary, X(1:2))) then
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
     call this%generate (r3(1:2), direction, Trace_Step, N_steps, Trace_Method, AltSurf, &
        theta_cut, sampling=sampling, AltSurf_fw=AltSurf_fw, AltSurf_bw=AltSurf_bw)
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

  end subroutine generate
!=======================================================================



!=======================================================================
! Generate flux surface branch from reference point xinit = (R,Z [cm])
!
! direction =
!    FORWARD         in forward Bpol direction
!    BACKWARD        in backward Bpol direction
!    CCW             in counter-clockwise direction
!    CW              in clockwise direction
!
! x0                 add 0th point (e.g. X-point)
! stop_at_boundary   stop tracing at configuration boundary (defined in
!                    module boundary)
! theta_cutoff       stop tracing at poloidal angle [rad]
! cutoff_boundary    stop tracing at this boundary
! cutoff_X           stop tracing at this distance from an X-point
!
! return ierr =
!    0               successful generation of flux surface branch
!    -1              flux surface branch leaves equilibrium domain
!    -2              flux surface branch connects to X-point
!=======================================================================
  subroutine generate_branch(this, xinit, direction, ierr, &
                             x0, trace_step, &
                             stop_at_boundary, theta_cutoff, cutoff_boundary, cutoff_X, connectX)
  use ode_solver
  use equilibrium, only: get_PsiN, get_poloidal_angle, Ip_sign, &
      leave_equilibrium_domain, nx_max, Xp, get_flux_surface_curvature, correct_PsiN
  use boundary
  use math
  use curve2D
  class(t_flux_surface_2D)   :: this
  real(real64),  intent(in)  :: xinit(2)
  integer,       intent(in)  :: direction
  integer,       intent(out) :: ierr
  real(real64),  intent(in),  optional :: x0(2), trace_step, theta_cutoff, cutoff_X
  logical,       intent(in),  optional :: stop_at_boundary
  type(t_curve), intent(in),  optional :: cutoff_boundary
  integer,       intent(out), optional :: connectX

  integer, parameter       :: chunk_size = 1024

  real(real64), dimension(:,:), allocatable :: xtmp, xtmp_tmp
  type(t_ODE)  :: F
  real(real64) :: L, ds, X(3), dX, t, thetal, thetac, dtheta
  integer      :: i, idir, ix, ierrPsiN, nchunks, ntmp, boundary_id
  logical      :: check_boundary


  ! 0. set up defaults for optional input/output
  ! 0.1. stop tracing at boundary?
  check_boundary = .true.
  if (present(stop_at_boundary)) then
     check_boundary = stop_at_boundary
  endif
  ! 0.2. connecting X-point id
  if (present(connectX)) then
     connectX = 0
  endif


  ! 1. initialize ------------------------------------------------------
  ! 1.1 trace direction
  select case(direction)
  case(FORWARD)
     idir = 1
  case(BACKWARD)
     idir = -1
  case(CCW)
     idir = -Ip_sign
  case(CW)
     idir = Ip_sign
  case default
     write (6, 9000)
     write (6, 9001) direction
     stop
  end select

  ! 1.2 temporary data array
  nchunks     = 1
  ntmp        = chunk_size
  allocate (xtmp(ntmp, 2))

  ! 1.3 set 0th data point (if present)
  i = 0
  L = 0.d0
  if (present(x0)) then
     i           = 1
     xtmp(1,1:2) = x0
     L           = sqrt(sum((xinit-x0)**2))
  endif
  ! set 1st data point
  i           = i + 1
  xtmp(i,1:2) = xinit
  this%PsiN   = get_PsiN(xinit)

  !---------------------------------------------------------------------


  ! 2. generate flux surface branch
  ierr = 0
  call F%init_ODE(2, xinit, 0.d0, epol_sub, NM_RUNGEKUTTA4)
  trace_loop: do
     i = i + 1

     ! A. increase temporary array, if necessary
     if (i > ntmp) then
        allocate (xtmp_tmp(ntmp, 2))
        xtmp_tmp       = xtmp
        deallocate (xtmp)

        nchunks        = nchunks + 1
        allocate (xtmp(nchunks*chunk_size, 2))
        xtmp(1:ntmp,:) = xtmp_tmp
        ntmp           = ntmp + chunk_size
        deallocate (xtmp_tmp)
     endif


     ! B. calculate next trace step size
     if (present(trace_step)) then
        ds       = idir * trace_step
     else
        ! calculate local curvature
        ds       = idir * adaptive_step_size(xtmp(i-1,1:2))
     endif
     ! D.5 check proximity to X-points
     if (present(cutoff_X)) then
        do ix=1,nx_max
           if (Xp(ix)%undefined) cycle

           dX = sqrt(sum((xtmp(i-1,1:2) - Xp(ix)%X)**2)) ! present distance from X-point

           ! stop tracing near X-point
           if (dX < cutoff_X) then
              ! check PsiN at present position vs. Xp(ix)%PsiN()
              xtmp(i,1:2) = Xp(ix)%X
              ierr        = -2
              if (present(connectX)) connectX = ix
              exit trace_loop
           endif

           ! else, adapt step size near X-point
           if (abs(ds) > dX/2.d0) ds = idir * dX/2.d0
        enddo
     endif


     ! C. trace one step
     X(1:2)      = F%step_ds(ds)
     xtmp(i,1:2) = correct_PsiN(X(1:2), this%PsiN, ierrPsiN) ! perform correction step
     if (ierrPsiN > 0) xtmp(i,1:2) = X(1:2)
     F%yc(1:2)   = xtmp(i,1:2) ! update ODE solver
     L           = L + abs(ds)
     thetac      = get_poloidal_angle(xtmp(i,1:2))


     ! D. boundary check
     ! D.1 check configuration boundary
     if (check_boundary) then
     if (intersect_boundary(xtmp(i-1,1:2), xtmp(i,1:2), X, boundary_id)) then
        xtmp(i,1:2) = X(1:2)
        exit
     endif
     endif
     ! D.2 check equilibrium boundary
     if (leave_equilibrium_domain(xtmp(i-1,1:2), xtmp(i,1:2), X(1:2))) then
        xtmp(i,1:2) = X(1:2)
        if (present(connectX)) connectX = -1
        ierr        = -1
        exit
     endif
     ! D.3 check cut-off poloidal angle
     if (present(theta_cutoff)) then
        dtheta = thetac - thetal
        if (abs(dtheta) > pi) dtheta = dtheta - sign(pi2,dtheta)
        if ((thetal+dtheta-theta_cutoff)*(thetal-theta_cutoff) < 0.d0) then
           xtmp(i,:) = xtmp(i-1,:) + abs((thetal-theta_cutoff)/dtheta)*(xtmp(i,:)-xtmp(i-1,:))
           exit
        endif
     endif
     ! D.4 check cut-off boundary
     if (present(cutoff_boundary)) then
     if (intersect_curve(xtmp(i-1,1:2), xtmp(i,1:2), cutoff_boundary, X(1:2))) then
        xtmp(i,1:2) = X(1:2)
        exit
     endif
     endif


     ! Z. prepare for next step
     thetal = thetac
  enddo trace_loop
  call this%new(i-1)
  this%x(:,:) = xtmp(1:i,1:2)
  call this%setup_length_sampling()


  ! 99. cleanup
  deallocate(xtmp)

 9000 format('error in t_flux_surface_2D%generate_branch:')
 9001 format('invalid direction ', i0, '!')
  end subroutine generate_branch
!=======================================================================



!=======================================================================
  function adaptive_step_size(r) result(ds)
  use equilibrium, only: get_flux_surface_curvature
  use math
  real(real64), intent(in) :: r(2)
  real(real64)             :: ds

  real(real64) :: dalpha, ds_min, ds_max


  ! internal parameters
  dalpha = pi/1800.d0 ! resolve 0.1 deg of local curvature
  ds_min = 1.d-6      ! minimum step size
  ds_max = 1.d0       ! maximum step size

  ds     = get_flux_surface_curvature(r) * dalpha
  ds     = max(min(abs(ds), ds_max), ds_min)

  end function adaptive_step_size
!=======================================================================



!=======================================================================
  subroutine generate_closed(this, r, direction, Trace_Step, N_steps, Trace_Method, AltSurf, retrace)
  use equilibrium
  class(t_flux_surface_2D) :: this
  real(real64), intent(in) :: r(2)

  integer, intent(in), optional       :: direction
  integer, intent(in), optional       :: N_steps, Trace_Method
  real(real64), intent(in), optional  :: Trace_Step
  type(t_curve), intent(in), optional :: AltSurf
  logical, intent(in), optional       :: retrace

  real(real64) :: theta_cut, r3(3), x1(2), x2(2), dl, dl0
  integer      :: direction0


  ! set default direction
  direction0 = RIGHT_HANDED
  if (present(direction)) direction0 = direction


  r3(1:2) = r
  r3(3)   = 0.d0
  theta_cut = get_poloidal_angle(r3)
  call this%generate(r, direction0, Trace_Step, N_steps, Trace_Method, AltSurf, theta_cut, retrace)

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
  use magnetic_axis
  use equilibrium, only: get_Bf_eq2D

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
  subroutine epol_sub (n, t, y, f)
  use equilibrium, only: get_Bf_eq2D
  integer,      intent(in)  :: n
  real(real64), intent(in)  :: t, y(n)
  real(real64), intent(out) :: f(n)

  real(real64) :: Bf(3), Bpol, y3(3)


  y3(1:2) = y
  y3(3)   = 0.d0
  Bf      = get_Bf_eq2D(y3)
  Bpol    = sqrt(Bf(1)**2 + Bf(2)**2)
  f       = Bf(1:2)/Bpol

  end subroutine epol_sub
!=======================================================================



!===============================================================================
! N_STEPS_GUESS
!
! Guess the required number of trace steps (for flux surface generation),
! based on the step size.
!===============================================================================
  function N_steps_guess (Trace_Step) result(N)
  use magnetic_axis
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
  subroutine setup_sampling(this, x1, x2, xc, DthetaR, DthetaL, ierr)
  use math
  class(t_flux_surface_2D) :: this
  real(real64), intent(in) :: x1(2), x2(2), xc(2), DthetaR, DthetaL
  integer,      intent(out), optional :: ierr

  real(real64) :: x(2), phi, dphi, dphi1, dphi2, phi1, phi2, s, f0, g0, w, Dtheta0
  integer      :: i, n, iseg1, iseg2


  if (present(ierr)) ierr  = 0
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
     write (6, *) 'flux surface (segment) is written to "error.plt"'
     write (6, *)
     write (6, *) 'geometry parameters are:'
     write (6, *) 'x1 (1st X-point)   =  ', x1
     write (6, *) 'x2 (2nd X-point)   =  ', x2
     write (6, *) 'xc (Magn. Axis)    =  ', xc
     write (6, *) 'DthetaR (transition angle, lower end) [deg] =  ', DthetaR / pi * 180.d0
     write (6, *) 'DthetaL (transition angle, upper end) [deg] =  ', DthetaL / pi * 180.d0
     write (6, *) 'Dtheta0 (reference angular weight) [deg]    =  ', Dtheta0 / pi * 180.d0
     write (6, *)
     call this%plot(filename='error.plt')
     if (present(ierr)) then
        ierr = 1
        return
     else
        stop
     endif
  else
     do i=1,n
        this%w(i) = this%w(i-1) + this%w(i)
     enddo
     this%w = this%w / this%w(n)
  endif
  !.....................................................................

  end subroutine setup_sampling
!===============================================================================



!===============================================================================
  function surface(this) result(area)
  class(t_flux_surface_2D)  :: this
  real(real64)              :: area

  real(real64) :: dPsi


  call this%surface_analysis(area, dPsi)

  end function surface
!===============================================================================



!===============================================================================
  subroutine surface_analysis(this, area, GradPsi)
  use math
  use equilibrium
  class(t_flux_surface_2D)  :: this
  real(real64), intent(out) :: area, GradPsi

  real(real64) :: x(2), dl(2), dA, dPsi(2)
  integer :: i


  area = 0.d0
  GradPsi = 0.d0
  do i=1,this%n_seg
     x    = 0.5d0 * (this%x(i-1,:) + this%x(i,:))
     dl   = this%x(i,:) - this%x(i-1,:)

     dA   = pi2 * sqrt(sum(dl**2)) * x(1)
     area = area + dA

     dPsi(1) = get_dPsiN(x, 1, 0)
     dPsi(2) = get_dPsiN(x, 0, 1)
     GradPsi = GradPsi + dA * sqrt(sum(dPsi**2))
  enddo
  GradPsi = GradPsi / area

  end subroutine surface_analysis
!===============================================================================



!===============================================================================
  function volume(this) result(V)
  use math
  class(t_flux_surface_2D) :: this
  real(real64)             :: V

  real(real64) :: x(2), dl(2), dA, dS(2)
  integer :: i


  V = 0.d0
  do i=1,this%n_seg
     x     = 0.5d0 * (this%x(i-1,:) + this%x(i,:))
     dl    = this%x(i,:) - this%x(i-1,:)
     dA    = pi2 * sqrt(sum(dl**2)) * x(1)

     dS(1) = dl(2);     dS(2) = -dl(1)
     dS    = dS / sqrt(sum(dS**2)) * dA

     V     = V + 0.5d0 * (dS(2) * x(2)  +  dS(1) * x(1) / 2.d0)
  enddo

  end function volume
!===============================================================================



!=======================================================================
  subroutine broadcast(this)
  use parallel
  class(t_flux_surface_2D) :: this

  call this%t_curve%broadcast()
  call broadcast_real_s(this%PsiN)

  end subroutine broadcast
!=======================================================================

end module flux_surface_2D
