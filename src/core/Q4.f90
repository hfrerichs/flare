!===============================================================================
! Basic 4-node quadrilateral (i.e. straight edges)
!
!
!            ^ eta               X(xi,eta) = a0  +  xi*a1  +  eta*a2  +  xi*eta*a3
! 4:(-1, 1)  |     3:( 1, 1)
!       X---------X              a0 = 1/4 * ( x1 + x2 + x3 + x4)
!       |         |              a1 = 1/4 * (-x1 + x2 + x3 - x4)
!       |         |---> xi       a2 = 1/4 * (-x1 - x2 + x3 + x4)
!       |         |              a3 = 1/4 * ( x1 - x2 + x3 - x4)
!       X---------X
! 1:(-1,-1)        2:( 1,-1)
!===============================================================================
module Q4
  use iso_fortran_env
  implicit none
  private

  ! intrinsic coordinates of nodes
  real(real64), parameter :: xi(4) = (/-1,  1,  1,  -1/), eta(4) = (/-1, -1,  1,  1/)

  type, public :: t_Q4
     ! node coordinates
     real(real64) :: x1(2), x2(2), x3(2), x4(2)

     ! shape coefficients
     real(real64) :: a0(2), a1(2), a2(2), a3(2)

     contains
     procedure :: set_nodes
     procedure :: set_shape
     procedure :: set_shape_advanced
     procedure :: area, Rcenter, Zcenter
  end type t_Q4

  contains
!=======================================================================



!=======================================================================
! Set node coordinates and setup shape coefficients
!=======================================================================
  subroutine set_nodes (this, x1, x2, x3, x4)
  class(t_Q4)                            :: this
  real(real64), dimension(2), intent(in) :: x1, x2, x3, x4


  ! set nodes coordinates
  this%x1 = x1
  this%x2 = x2
  this%x3 = x3
  this%x4 = x4

  ! setup shape coefficients
  this%a0 = 0.25d0 * ( x1 + x2 + x3 + x4)
  this%a1 = 0.25d0 * (-x1 + x2 + x3 - x4)
  this%a2 = 0.25d0 * (-x1 - x2 + x3 + x4)
  this%a3 = 0.25d0 * ( x1 - x2 + x3 - x4)

  end subroutine set_nodes
!=======================================================================



!=======================================================================
! Set shape coefficients a0, a1, a2, a3 and setup nodes
!=======================================================================
  subroutine set_shape (this, a0, a1, a2, a3)
  class(t_Q4)                            :: this
  real(real64), dimension(2), intent(in) :: a0, a1, a2, a3


  ! set shape coefficients
  this%a0 = a0
  this%a1 = a1
  this%a2 = a2
  this%a3 = a3

  ! setup nodes
  this%x1 = a0  +  xi(1)*a1  +  eta(1)*a2  +  xi(1)*eta(1)*a3
  this%x2 = a0  +  xi(2)*a1  +  eta(2)*a2  +  xi(2)*eta(2)*a3
  this%x3 = a0  +  xi(3)*a1  +  eta(3)*a2  +  xi(3)*eta(3)*a3
  this%x4 = a0  +  xi(4)*a1  +  eta(4)*a2  +  xi(4)*eta(4)*a3

  end subroutine set_shape
!=======================================================================



!=======================================================================
! Set shape from basic geometry (parallelogram) and distortion measures
!                                            ^ a2 = ey
! alpha1 = alpha2 = 0 provides a rectangle   |
! with positive orientation:                 |-----> a1 = ex
!=======================================================================
  subroutine set_shape_advanced(this, x0, b1, alpha1, b2, alpha2, P, theta)
  use math
  class(t_Q4)              :: this
  real(real64), intent(in) :: x0(2), b1, alpha1, b2, alpha2, P, theta

  real(real64) :: a1(2), a2(2), a3(2), P13, P23


  a1(1) = b1 * cos(alpha1)
  a1(2) = b1 * sin(alpha1)

  a2(1) = b2 * cos(alpha2+pi/2.d0)
  a2(2) = b2 * sin(alpha2+pi/2.d0)

  P13   = P * cos(theta)
  P23   = P * sin(theta)
  a3    = P13*a2 - P23*a1
  call this%set_shape(x0, a1, a2, a3)

  end subroutine set_shape_advanced
!=======================================================================



!=======================================================================
! calculate area of quadrilateral
!=======================================================================
  function area(this)
  class(t_Q4)  :: this
  real(real64) :: area

  real(real64) :: abcd1, abcd2

  abcd1 = (this%x3(1)-this%x2(1))*(this%x2(2)-this%x1(2)) &
        - (this%x2(1)-this%x1(1))*(this%x3(2)-this%x2(2))
  abcd2 = (this%x4(1)-this%x1(1))*(this%x3(2)-this%x4(2)) &
        - (this%x3(1)-this%x4(1))*(this%x4(2)-this%x1(2))
  area  = 0.5d0 * (abcd1 + abcd2)

  end function area
!=======================================================================



!=======================================================================
! calculate center point of quadrilateral (should be this%a0)
!=======================================================================
  function Rcenter(this)
  class(t_Q4)  :: this
  real(real64) :: Rcenter

  Rcenter = (this%x1(1)+this%x2(1)+this%x3(1)+this%x4(1)) / 4.d0
  end function Rcenter
!=======================================================================
  function Zcenter(this)
  class(t_Q4)  :: this
  real(real64) :: Zcenter

  Zcenter = (this%x1(2)+this%x2(2)+this%x3(2)+this%x4(2)) / 4.d0
  end function Zcenter
!=======================================================================

end module Q4
