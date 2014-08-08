!===============================================================================
! Interface for field line tracing (either by numerical integration or by
! reconstruction)
!===============================================================================
module fieldline
  use math
  use ode_solver
  use bfield
  use boundary
  implicit none

  integer, parameter :: FL_LINE  = 1
  integer, parameter :: FL_ARC   = 2
  integer, parameter :: FL_ANGLE = 3

  integer, parameter :: FL_Reconstruction = 0

  type, extends(t_ODE) :: t_fieldline
     ! integrated toroidal and poloidal angle
     real*8 :: phic, thetac
     ! last and current state in cylindrical coordinates
     real*8 :: rl(3), rc(3)

     ! toroidal distance to last intersection with symmetry plane
     real*8 :: phit
     real*8 :: phi_sym
     integer :: iplane

     contains
     procedure :: init, trace, init_toroidal_tracing, intersect_sym_plane

     procedure :: intersect_boundary => fieldline_intersects_boundary
  end type t_fieldline

  contains
!=======================================================================



!=======================================================================
  subroutine init (this, y0, ds, isolver, icoord)
  class (t_fieldline) :: this
  real*8, intent(in)  :: y0(3), ds
  integer, intent(in) :: isolver, icoord

  real*8 :: y1(3)


  y1 = y0
  ! use field line reconstruction method
  if (isolver == FL_Reconstruction) then

  ! use numerical integration method isolver > 0
  elseif (isolver > 0) then
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

  end subroutine init
!=======================================================================



!=======================================================================
  function trace(this) result(yc)
  class(t_fieldline), intent(inout) :: this
  real*8                            :: yc(3)

  this%rl = this%rc
  yc      = this%next_step()
  call coord_trans (yc, Trace_Coords, this%rc, CYLINDRICAL)

  end function trace
!=======================================================================



!=======================================================================
  function fieldline_intersects_boundary(this, rcut) result(l)
  class(t_fieldline), intent(inout) :: this
  real*8, intent(out)               :: rcut(3)
  logical                           :: l

  l = intersect_boundary (this%rl, this%rc, rcut)

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
     icut      = icut + 1
     rcut      = this%rl + fff * (this%rc-this%rl)
  else
     this%phit = this%phit + dphi
     l         = .false.
  endif

  end function intersect_sym_plane
!=======================================================================



!=======================================================================
  subroutine init_toroidal_tracing(this, y0, ds, isolver, icoord, nsym, phi_out)
  use equilibrium

  class (t_fieldline) :: this
  real*8, intent(in)  :: y0(3), ds, phi_out
  integer, intent(in) :: isolver, icoord, nsym

  real*8 :: sgn


  call this%init(y0, ds, isolver, icoord)


  ! for Poincare plots
  this%iplane = 0
  this%phi_sym = pi2 / nsym
  if (icoord == 2) then
     sgn  = 1.d0 * Bt_sign
  else
     sgn  = sign(1.d0, ds) * Bt_sign
  endif
  this%phit    = mod(this%rc(3), this%phi_sym)
  this%phit    = this%phit - sgn * phi_out * pi2/360.d0

  end subroutine
!=======================================================================



!=======================================================================
  subroutine Bf_sub_cart (n, t, y, f)
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
