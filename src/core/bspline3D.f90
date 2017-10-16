!===============================================================================
! B-spline interpolation in 3D (on structured grid)
!===============================================================================
module bspline3D
  use iso_fortran_env
  use bspline
  implicit none
  private


  type, public :: t_bspline3D
     ! coefficients for spline interpolation
     real(real64), dimension(:,:,:), allocatable :: dcoeff
     real(real64), dimension(:),     allocatable :: x1not, x2not, x3not

     integer      :: n1, n2, n3

     ! spline interpolation order
     integer      :: nord

     contains
     procedure :: new
     procedure :: init
     procedure :: broadcast
     procedure :: eval
     procedure :: derivative
     procedure :: cleanup
  end type t_bspline3D

  contains
!=======================================================================



!=======================================================================
! allocate memory for spline coefficients and knots
!
! n1, n2, n3           resolution in x1, x2 and x3 direction
! nord                 spline interpolation order
!=======================================================================
  subroutine new(this, n1, n2, n3, nord)
  class(t_bspline3D)  :: this
  integer, intent(in) :: n1, n2, n3
  integer, intent(in), optional :: nord


  this%n1   = n1
  this%n2   = n2
  this%n3   = n3
  if (present(nord)) then
     this%nord = nord
  else
     this%nord = 5
  endif
  call this%cleanup()


  allocate (this%x1not(n1+this%nord), this%x2not(n2+this%nord), this%x3not(n3+this%nord))
  allocate (this%dcoeff(n1, n2, n3))

  end subroutine new
!=======================================================================



!=======================================================================
! initialize 3D spline from given data (original data is not stored here)
!
! n1, n2, n3           resolution in x1, x2 and x3 direction
! nord                 spline interpolation order
!=======================================================================
  subroutine init(this, n1, n2, n3, x1, x2, x3, D, nord)
  class(t_bspline3D)       :: this
  integer,      intent(in) :: n1, n2, n3
  real(real64), intent(in) :: x1(n1), x2(n2), x3(n3), D(n1,n2,n3)
  integer,      intent(in), optional :: nord


  call this%new(n1, n2, n3, nord)

  call dbsnak(n1, x1, this%nord, this%x1not)
  call dbsnak(n2, x2, this%nord, this%x2not)
  call dbsnak(n3, x3, this%nord, this%x3not)
  call dbs3in(n1, x1, n2, x2, n3, x3, D, &
              n1, n2, this%nord, this%nord, this%nord, this%x1not, this%x2not, this%x3not, this%dcoeff)

  end subroutine init
!=======================================================================



!=======================================================================
  subroutine broadcast(this)
  use parallel
  class(t_bspline3D)   :: this


  call broadcast_inte_s(this%n1)
  call broadcast_inte_s(this%n2)
  call broadcast_inte_s(this%n3)
  call broadcast_inte_s(this%nord)
  if (mype > 0) then
     call this%new(this%n1, this%n2, this%n3, this%nord)
  endif

  ! broadcast spline coefficients
  call broadcast_real  (this%x1not,  this%n1 + this%nord)
  call broadcast_real  (this%x2not,  this%n2 + this%nord)
  call broadcast_real  (this%x3not,  this%n3 + this%nord)
  call broadcast_real  (this%dcoeff, this%n1*this%n2*this%n3)

  end subroutine broadcast
!=======================================================================



!=======================================================================
  function eval(this, x) result(d)
  class(t_bspline3D)       :: this
  real(real64), intent(in) :: x(3)
  real(real64)             :: d


  d    = dbs3dr(0,0,0,x(1),x(2),x(3),this%nord,this%nord,this%nord, &
                this%x1not,this%x2not,this%x3not,this%n1,this%n2,this%n3,this%dcoeff)

  end function eval
!=======================================================================



!=======================================================================
  function derivative(this, x, m1, m2, m3) result(d)
  class(t_bspline3D)       :: this
  real(real64), intent(in) :: x(3)
  integer,      intent(in) :: m1, m2, m3
  real(real64)             :: d


  d    = dbs3dr(m1,m2,m3,x(1),x(2),x(3),this%nord,this%nord,this%nord, &
                this%x1not,this%x2not,this%x3not,this%n1,this%n2,this%n3,this%dcoeff)

  end function derivative
!=======================================================================



!=======================================================================
  subroutine cleanup(this)
  class(t_bspline3D) :: this


  if (allocated(this%x1not))  deallocate(this%x1not)
  if (allocated(this%x2not))  deallocate(this%x2not)
  if (allocated(this%x3not))  deallocate(this%x3not)
  if (allocated(this%dcoeff)) deallocate(this%dcoeff)

  end subroutine cleanup
!=======================================================================

end module bspline3D
