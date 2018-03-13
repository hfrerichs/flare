!===============================================================================
! B-spline interpolation in 1D (on structured grid)
!===============================================================================
module bspline1D
  use iso_fortran_env
  use bspline
  implicit none
  private


  type, public :: t_bspline1D
     ! coefficients for spline interpolation
     real(real64), dimension(:), allocatable :: dcoeff
     real(real64), dimension(:), allocatable :: x1not

     integer      :: n1

     ! spline interpolation order
     integer      :: nord

     contains
     procedure :: new
     procedure :: init
     procedure :: broadcast
     procedure :: eval
     procedure :: derivative
     final     :: cleanup
  end type t_bspline1D

  contains
!=======================================================================



!=======================================================================
! allocate memory for spline coefficients and knots
!
! n1                   resolution in x1 direction
! nord                 spline interpolation order
!=======================================================================
  subroutine new(this, n1, nord)
  class(t_bspline1D)  :: this
  integer, intent(in) :: n1
  integer, intent(in), optional :: nord


  this%n1   = n1
  if (present(nord)) then
     this%nord = nord
  else
     this%nord = 5
  endif
  call cleanup(this)


  allocate (this%x1not(n1+this%nord))
  allocate (this%dcoeff(n1))

  end subroutine new
!=======================================================================



!=======================================================================
! initialize 3D spline from given data (original data is not stored here)
!
! n1                   resolution in x1 direction
! nord                 spline interpolation order
!=======================================================================
  subroutine init(this, n1, x1, D, nord)
  class(t_bspline1D)       :: this
  integer,      intent(in) :: n1
  real(real64), intent(in) :: x1(n1), D(n1)
  integer,      intent(in), optional :: nord


  call this%new(n1, nord)

  call dbsnak(n1, x1, this%nord, this%x1not)
  call dbsint(n1, x1, D, this%nord, this%x1not, this%dcoeff)

  end subroutine init
!=======================================================================



!=======================================================================
  subroutine broadcast(this)
  use parallel
  class(t_bspline1D)   :: this


  call broadcast_inte_s(this%n1)
  call broadcast_inte_s(this%nord)
  if (mype > 0) then
     call this%new(this%n1, this%nord)
  endif

  ! broadcast spline coefficients
  call broadcast_real  (this%x1not,  this%n1 + this%nord)
  call broadcast_real  (this%dcoeff, this%n1)

  end subroutine broadcast
!=======================================================================



!=======================================================================
  function eval(this, x) result(d)
  class(t_bspline1D)       :: this
  real(real64), intent(in) :: x
  real(real64)             :: d


  d    = dbsval(x, this%nord, this%x1not, this%n1, this%dcoeff)

  end function eval
!=======================================================================



!=======================================================================
  function derivative(this, x, m1) result(d)
  class(t_bspline1D)       :: this
  real(real64), intent(in) :: x
  integer,      intent(in) :: m1
  real(real64)             :: d


  d    = dbsder(m1, x, this%nord, this%x1not, this%n1, this%dcoeff)

  end function derivative
!=======================================================================



!=======================================================================
  subroutine cleanup(this)
  type(t_bspline1D) :: this


  if (allocated(this%x1not))  deallocate(this%x1not)
  if (allocated(this%dcoeff)) deallocate(this%dcoeff)

  end subroutine cleanup
!=======================================================================

end module bspline1D
