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

  integer, parameter :: FL_Reconstruction = 0

  type, extends(t_ODE) :: t_fieldline
     ! integrated toroidal and poloidal angle
     real*8 :: phi0, phi_int, theta0, theta_int, Dphi
     ! last (l) and current (c) state in cylindrical coordinates
     real*8 :: rl(3), rc(3), thetal, thetac
     real(real64) :: PsiNl, PsiNc

     ! toroidal distance to last intersection with symmetry plane
     real*8 :: phit
     real*8 :: phi_sym
     integer :: iplane, sgn

     integer :: Trace_Coords

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
  subroutine init (this, y0, ds, isolver, icoord)
  use equilibrium
  class (t_fieldline) :: this
  real*8, intent(in)  :: y0(3), ds
  integer, intent(in) :: isolver, icoord

  real*8 :: y1(3)


  y1 = y0
  ! use field line reconstruction method
  if (isolver == FL_Reconstruction) then
     this%trace_1step  => trace_1step_reconstruct

  ! use numerical integration method isolver > 0
  elseif (isolver > 0) then
     this%Trace_Coords = icoord
     this%trace_1step  => trace_1step_ODE

     select case (icoord)
     case (FL_LINE)
        call this%init_ODE (3, y1, ds, Bf_sub_cart, isolver)
     case (FL_ARC)
        call this%init_ODE (3, y1, ds, Bf_sub_cyl, isolver)
     case (FL_ANGLE)
        call this%init_ODE (3, y1, ds, Bf_sub_cyl_norm, isolver)
     case default
        write (6, *) 'invalid parameter icoord = ', icoord
        stop
     end select
  endif


  ! initialize internal variables
  call coord_trans (y0, icoord, this%rc, CYLINDRICAL)
  this%phi_int   = 0.d0
  this%phi0      = this%rc(3)
  this%theta_int = 0.d0
  this%thetac    = get_poloidal_angle(this%rc)
  this%theta0    = this%thetac
  this%PsiNc     = get_PsiN(this%rc)
  this%Dphi      = 0.d0

  end subroutine init
!=======================================================================



!=======================================================================
! Trace field line one step size
!=======================================================================
  subroutine trace_1step_ODE(this)
  use equilibrium
  class(t_fieldline), intent(inout) :: this

  real(real64) :: yc(3), Dtheta

  ! save last step
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
  this%theta_int = this%theta_int + Dtheta


  ! update radial coordinate
  this%PsiNc     = get_PsiN(this%rc)

  end subroutine trace_1step_ODE
!=======================================================================
  subroutine trace_1step_reconstruct(this)
  class(t_fieldline), intent(inout) :: this

  end subroutine trace_1step_reconstruct
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
!=======================================================================
  subroutine trace_Dphi(this, Dphi, stop_at_boundary, yout, ierr)
  class (t_fieldline)       :: this
  real(real64), intent(in)  :: Dphi
  logical,      intent(in)  :: stop_at_boundary
  real(real64), intent(out) :: yout(3)
  integer,      intent(out) :: ierr

  real(real64) :: phi_int, f


  ierr         = 0
  yout         = 0.d0
  phi_int      = -this%phi_int ! offset from last call
  this%phi_int = 0.d0          ! reset integration distance
  trace_loop: do
     call this%trace_1step()

     ! check intersection with boundary
     if (stop_at_boundary  .and.  this%intersect_boundary(yout)) then
        ierr = 1
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

  real*8 :: yc(3), lc


  lc = 0.d0
  trace_loop: do
     call this%trace_1step()
     lc      = lc + this%ds

     if (abs(lc) > Limit) exit trace_loop

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
  function intersect_sym_plane(this, icut, rcut) result(l)
  class(t_fieldline), intent(inout) :: this
  integer, intent(inout) :: icut
  real*8, intent(out)    :: rcut(3)
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
  else
     this%phit = this%phit + dphi
     l         = .false.
  endif

  end function intersect_sym_plane
!=======================================================================



!=======================================================================
  subroutine init_toroidal_tracing(this, y0, ds, isolver, icoord, nsym, phi_out)
  use magnetic_axis
  !use equilibrium

  class (t_fieldline) :: this
  real*8, intent(in)  :: y0(3), ds, phi_out
  integer, intent(in) :: isolver, icoord, nsym


  call this%init(y0, ds, isolver, icoord)


  ! for Poincare plots
  this%iplane = 0
  this%phi_sym = pi2 / nsym
  if (icoord == 2) then
     this%sgn  = Bt_sign
  else
     this%sgn  = nint(sign(1.d0, ds)) * Bt_sign
  endif
  this%phit    = mod(this%rc(3), this%phi_sym)
  this%phit    = this%phit - this%sgn * phi_out * pi2/360.d0

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


end module fieldline
