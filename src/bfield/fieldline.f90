!===============================================================================
! Interface for field line tracing (either by numerical integration or by
! reconstruction)
!===============================================================================
module fieldline
  use ode_solver
  use bfield
  use math
  implicit none

  integer, parameter :: FL_LINE  = 1
  integer, parameter :: FL_ARC   = 2
  integer, parameter :: FL_ANGLE = 3

  integer, parameter :: FL_Reconstruction = 0

  type, extends(t_ODE) :: t_fieldline
     contains
     procedure :: init
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
        y1(3) = y1(3) / 180.d0 * pi	! deg -> rad
        call this%init_ODE (3, y1, ds, Bf_sub_cyl, isolver)
     case (FL_ANGLE)
        y1(3) = y1(3) / 180.d0 * pi	! deg -> rad
        call this%init_ODE (3, y1, ds, Bf_sub_cyl_norm, isolver)
     case default
        write (6, *) 'invalid parameter icoord = ', icoord
        stop
     end select
  endif

  end subroutine init
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
  f(3) = f(3) / y(1) * 180.d0/pi

  end subroutine Bf_sub_cyl
!=======================================================================



!=======================================================================
  subroutine Bf_sub_cyl_norm (n, t, y, f)
  implicit none

  integer, intent(in) :: n
  real*8,  intent(in) :: t, y(n)
  real*8, intent(out) :: f(n)

  f      = get_Bf_Cyl(y)
  f(1:2) = y(1) * f(1:2) / f(3) *pi/180.d0
  f(3)   = 1.d0

  end subroutine Bf_sub_cyl_norm
!=======================================================================


end module fieldline
