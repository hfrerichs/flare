!===============================================================================
! Interface for field line tracing (either by numerical integration or by
! reconstruction)
!===============================================================================
module fieldline
  use iso_fortran_env
  use math
  use ode_solver
  use reconstruct
  implicit none

  integer, parameter :: FL_LINE  = 1
  integer, parameter :: FL_ARC   = 2
  integer, parameter :: FL_ANGLE = 3

  integer, parameter :: &
     FULL_FIELD        = 0, &
     EQUILIBRIUM_FIELD = 1

  integer, parameter :: FL_Reconstruction = 0

  type, extends(t_ODE) :: t_fieldline
     ! integrated toroidal and poloidal angle
     real*8 :: phi0, phi_int, theta0, theta_int, Dphi, lc
     ! last (l) and current (c) state in cylindrical coordinates
     real*8 :: rl(3), rc(3), thetal, thetac, Dtheta
     real(real64) :: PsiNl, PsiNc

     ! toroidal distance to last intersection with symmetry plane
     real*8 :: phit
     real*8 :: phi_sym
     integer :: iplane, sgn

     real(real64) :: Trace_Step
     integer :: Trace_Coords, ierr

     type(t_fluxtube_coords) :: F

     procedure(trace_1step_ODE), pointer :: trace_1step

     contains
     procedure :: init
     procedure :: trace
     procedure :: init_toroidal_tracing
     procedure :: intersect_sym_plane
     procedure :: cross_PsiN
     procedure :: get_PsiN => fieldline_get_PsiN
     procedure :: get_flux_coordinates
     procedure :: trace_Dphi
     procedure :: intersect_boundary => fieldline_intersects_boundary
  end type t_fieldline

  contains
!=======================================================================



!=======================================================================
! trace parameter (Trace_Coords), D = Trace_Step

! 1:     D = sqrt(DX**2 + DY**2 + DZ**2)
! 2:     D = sqrt(DR**2 + DZ**2 + (R*DPHI)**2)
! 3:     D = DPHI
!=======================================================================
  subroutine init (this, y0, Trace_Step, Trace_Method, Trace_Coords, bfield)
  use numerics, only: Trace_Step_default   => Trace_Step, &
                      Trace_Method_default => Trace_Method, &
                      Trace_Coords_default => Trace_Coords
  use exceptions, only: OUT_OF_BOUNDS
  use equilibrium
  class (t_fieldline) :: this
  real*8, intent(in)  :: y0(3)
  real*8, intent(in), optional  :: Trace_Step
  integer, intent(in), optional :: Trace_Method, Trace_Coords, bfield

  real*8  :: y1(3), ds
  integer :: tm, bf, tc


  this%ierr     = 0
  OUT_OF_BOUNDS = .false.


  ! set step size
  ds = Trace_Step_default
  if (present(Trace_Step)) ds = Trace_Step

  ! set numerical method for field line tracing
  tm = Trace_Method_default
  if (present(Trace_Method)) tm = Trace_Method

  ! set discretization coordinates
  this%Trace_Coords = Trace_Coords_default
  if (present(Trace_Coords)) this%Trace_Coords = Trace_Coords

  ! set magnetic field model
  bf = FULL_FIELD
  if (present(bfield)) bf = bfield


  y1 = y0
  ! use field line reconstruction method
  if (tm == FL_Reconstruction) then
     this%trace_1step  => trace_1step_reconstruct

  ! use numerical integration method isolver > 0
  elseif (tm > 0) then
     this%Trace_Step   =  ds

     this%trace_1step  => trace_1step_ODE

     tc = this%Trace_Coords + 10*bf
     select case (tc)
     case (FL_LINE)
        call this%init_ODE (3, y1, ds, Bf_sub_cart,     tm)
     case (FL_ARC)
        call this%init_ODE (3, y1, ds, Bf_sub_cyl,      tm)
     case (FL_ARC   + 10*EQUILIBRIUM_FIELD)
        call this%init_ODE (3, y1, ds, Bfeq_sub_cyl,    tm)
     case (FL_ANGLE)
        call this%init_ODE (3, y1, ds, Bf_sub_cyl_norm, tm)
     case (FL_ANGLE + 10*EQUILIBRIUM_FIELD)
        call this%init_ODE (3, y1, ds, Bfeq_sub_cyl_norm, tm)
     case default
        write (6, *) 'invalid parameter Trace_Coords = ', this%Trace_Coords
        stop
     end select
  endif


  ! initialize internal variables
  call coord_trans (y0, this%Trace_Coords, this%rc, CYLINDRICAL)
  this%phi_int   = 0.d0
  this%phi0      = this%rc(3)
  this%theta_int = 0.d0
  this%thetac    = get_poloidal_angle(this%rc)
  this%theta0    = this%thetac
  this%PsiNc     = get_PsiN(this%rc)
  this%Dphi      = 0.d0
  this%lc        = 0.d0

  if (OUT_OF_BOUNDS) this%ierr = 2

  end subroutine init
!=======================================================================



!=======================================================================
! Trace field line one step size
!=======================================================================
  function trace_1step_ODE(this) result(dl)
  use equilibrium
  use exceptions, only: OUT_OF_BOUNDS
  class(t_fieldline), intent(inout) :: this
  real(real64) :: dl

  real(real64) :: yc(3), Dtheta

  ! save last step
  OUT_OF_BOUNDS  = .false.
  this%rl        = this%rc
  this%thetal    = this%thetac
  this%PsiNl     = this%PsiNc


  ! calculate next step
  yc             = this%next_step()
  call coord_trans (yc, this%Trace_Coords, this%rc, CYLINDRICAL)


  ! update toroidal angle
  this%Dphi      = this%rc(3) - this%rl(3)
  if (abs(this%Dphi) > pi) this%Dphi = this%Dphi - sign(pi2,this%Dphi)
  this%phi_int   = this%phi_int + this%Dphi


  ! update poloidal angle
  this%thetac    = get_poloidal_angle(this%rc)
  Dtheta         = this%thetac - this%thetal
  if (abs(Dtheta) > pi) Dtheta = Dtheta - sign(pi2,Dtheta)
  this%Dtheta    = Dtheta
  this%theta_int = this%theta_int + Dtheta


  ! update radial coordinate
  this%PsiNc     = get_PsiN(this%rc)


  ! update trace step
  dl = this%Trace_Step
  if (this%Trace_Coords == FL_ANGLE) dl = dl * 0.5d0 * (this%rl(1)+this%rc(1))


  if (OUT_OF_BOUNDS) this%ierr = 2

  end function trace_1step_ODE
!=======================================================================
  function trace_1step_reconstruct(this) result(dl)
  class(t_fieldline), intent(inout) :: this
  real(real64)                      :: dl

  end function trace_1step_reconstruct
!=======================================================================



!=======================================================================
  function fieldline_get_PsiN(this) result(PsiN)
  use equilibrium
  class(t_fieldline) :: this
  real*8             :: PsiN

  PsiN = get_PsiN(this%rc)

  end function fieldline_get_PsiN
!=======================================================================



!=======================================================================
  function cross_PsiN(this, PsiN) result(l)
  class(t_fieldline) :: this
  real(real64)       :: PsiN
  logical            :: l


  l = .false.
  if ((this%PsiNl-PsiN)*(this%PsiNc-PsiN) <= 0.d0) l = .true.

  end function cross_PsiN
!=======================================================================



!=======================================================================
  function get_flux_coordinates (this) result(y)
  class(t_fieldline) :: this
  real(real64)       :: y(3)


  y(1) = this%theta0 + this%theta_int
  y(2) = this%PsiNc
  y(3) = this%phi0   + this%phi_int

  end function get_flux_coordinates
!=======================================================================



!=======================================================================
! trace field line for a toroidal distance of Dphi > 0 (direction is not
! checked)
!
! phi_int will be reset after a call to this subroutine
!
! ierr = 1 if stop_at_boundary == .true. and field line intersects boundary
!        2 if OUT_OF_BOUNDS
!=======================================================================
  subroutine trace_Dphi(this, Dphi, stop_at_boundary, yout, ierr)
  use exceptions, only: OUT_OF_BOUNDS
  class (t_fieldline)       :: this
  real(real64), intent(in)  :: Dphi
  logical,      intent(in)  :: stop_at_boundary
  real(real64), intent(out) :: yout(3)
  integer,      intent(out) :: ierr

  real(real64) :: phi_int, f, dl


  ierr         = 0
  yout         = 0.d0
  phi_int      = -this%phi_int ! offset from last call
  this%phi_int = 0.d0          ! reset integration distance
  trace_loop: do
     dl = this%trace_1step()

     ! check intersection with boundary
     if (stop_at_boundary  .and.  this%intersect_boundary(yout)) then
        ierr = 1
        return
     endif

     if (OUT_OF_BOUNDS) then
        yout = this%rl
        ierr = 2
        return
     endif

     ! stop tracing after toroidal distance Dphi
     f = abs(this%phi_int - phi_int) - Dphi
     if (f > 0.d0) then
        f         = f / abs(this%Dphi)
        yout      = f*this%rl + (1.d0-f)*this%rc

        ! store offset for next call
        this%phi_int = f*this%Dphi
        return
     endif
  enddo trace_loop

  end subroutine trace_Dphi
!=======================================================================



!=======================================================================
  subroutine trace(this, Limit, stop_at_boundary)
  class (t_fieldline) :: this
  real*8, intent(in)  :: Limit
  logical, intent(in) :: stop_at_boundary

  real*8 :: yc(3), dl


  this%lc = 0.d0
  trace_loop: do
     dl = this%trace_1step()
     this%lc = this%lc + dl

     if (abs(this%lc) > Limit) exit trace_loop

     if (stop_at_boundary) then
        if (this%intersect_boundary()) then
           exit trace_loop
        endif
     endif
  enddo trace_loop

  end subroutine trace
!=======================================================================



!=======================================================================
  function fieldline_intersects_boundary(this, rcut, id, ielem, tau) result(l)
  use boundary
  class(t_fieldline), intent(inout) :: this
  real*8, intent(out), optional     :: rcut(3)
  integer, intent(out), optional    :: id, ielem
  real(real64), intent(out), optional :: tau
  logical                           :: l

  real*8  :: X(3)
  integer :: id1, ielem1


  l = intersect_boundary (this%rl, this%rc, X, id1, ielem1, tau)
  if (present(rcut)) rcut = X
  if (present(id))   id   = id1
  if (present(ielem)) ielem = ielem1

  end function fieldline_intersects_boundary
!=======================================================================



!=======================================================================
  function intersect_sym_plane(this, icut, rcut, thcut) result(l)
  class(t_fieldline), intent(inout) :: this
  integer, intent(inout) :: icut
  real*8, intent(out)    :: rcut(3), thcut
  logical :: l

  real*8 :: dphi, fff


  dphi      = dabs(this%rc(3) - this%rl(3))
  if (dphi > pi) dphi = pi2 - dphi


  if ((this%phit - this%phi_sym)*(this%phit + dphi - this%phi_sym) .le. 0.d0) then
     fff       = (this%phi_sym - this%phit) / dphi
     this%phit = (1.d0 - fff) * dphi
     l         = .true.
     icut      = icut + this%sgn
     rcut      = this%rl + fff * (this%rc-this%rl)
     thcut     = this%theta0 + this%theta_int-this%Dtheta + fff * this%Dtheta
  else
     this%phit = this%phit + dphi
     l         = .false.
  endif

  end function intersect_sym_plane
!=======================================================================



!=======================================================================
  subroutine init_toroidal_tracing(this, y0, ds, isolver, icoord, nsym, phi_out)
  use magnetic_axis
  use bfield
  use parallel
  !use equilibrium

  class (t_fieldline) :: this
  real*8, intent(in)  :: y0(3), ds, phi_out
  integer, intent(in) :: isolver, icoord, nsym

  real*8 :: Bf0(3)


  call this%init(y0, ds, isolver, icoord)


  ! Bt_sign is undefined (because no magnetic axis has been defined)
  if (Bt_sign == 0) then
     Bt_sign = 1
     Bf0     = get_Bf_Cyl(this%rc)
     if (Bf0(3) < 0.d0) Bt_sign = -1

     if (firstP) write (6, *) 'WARNING: magnetic axis is still undefined, poloidal angles will not be correct!'
  endif

  ! for Poincare plots
  this%iplane = 0
  this%phi_sym = pi2 / nsym
  if (icoord == FL_ANGLE) then
     this%sgn  = Bt_sign
  else
     this%sgn  = nint(sign(1.d0, ds)) * Bt_sign
  endif
  this%phit    = mod(this%rc(3), this%phi_sym)
  this%phit    = this%phit - this%sgn * phi_out

  end subroutine
!=======================================================================



!=======================================================================
  subroutine Bf_sub_cart (n, t, y, f)
  use bfield
  implicit none

  integer, intent(in) :: n
  real*8,  intent(in) :: t, y(n)
  real*8, intent(out) :: f(n)

  f = get_Bf_Cart(y)
  f = f / dsqrt(sum(f**2))

  end subroutine Bf_sub_cart
!=======================================================================



!=======================================================================
  subroutine Bf_sub_cyl (n, t, y, f)
  use bfield
  implicit none

  integer, intent(in) :: n
  real*8,  intent(in) :: t, y(n)
  real*8, intent(out) :: f(n)

  f    = get_Bf_Cyl(y)
  f    = f / dsqrt(sum(f**2))
  f(3) = f(3) / y(1)

  end subroutine Bf_sub_cyl
!=======================================================================



!=======================================================================
  subroutine Bfeq_sub_cyl (n, t, y, f)
  use equilibrium
  implicit none

  integer, intent(in) :: n
  real*8,  intent(in) :: t, y(n)
  real*8, intent(out) :: f(n)

  f    = get_Bf_eq2D(y)
  f    = f / dsqrt(sum(f**2))
  f(3) = f(3) / y(1)

  end subroutine Bfeq_sub_cyl
!=======================================================================



!=======================================================================
  subroutine Bf_sub_cyl_norm (n, t, y, f)
  use bfield
  implicit none

  integer, intent(in) :: n
  real*8,  intent(in) :: t, y(n)
  real*8, intent(out) :: f(n)

  f      = get_Bf_Cyl(y)
  f(1:2) = y(1) * f(1:2) / f(3)
  f(3)   = 1.d0

  end subroutine Bf_sub_cyl_norm
!=======================================================================



!=======================================================================
  subroutine Bfeq_sub_cyl_norm (n, t, y, f)
  use equilibrium
  implicit none

  integer, intent(in) :: n
  real*8,  intent(in) :: t, y(n)
  real*8, intent(out) :: f(n)

  f      = get_Bf_eq2D(y)
  f(1:2) = y(1) * f(1:2) / f(3)
  f(3)   = 1.d0

  end subroutine Bfeq_sub_cyl_norm
!=======================================================================


end module fieldline
