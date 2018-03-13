!===============================================================================
! B-spline interpolation in 2D (on structured grid)
!===============================================================================
module bspline2D
  use iso_fortran_env
  use bspline
  implicit none
  private


  type, public :: t_bspline2D
     ! coefficients for spline interpolation
     real(real64), dimension(:,:), allocatable :: dcoeff
     real(real64), dimension(:),   allocatable :: xnot, ynot

     integer      :: nx, ny

     ! spline interpolation order
     integer      :: nord

     contains
     procedure :: new
     procedure :: init
     procedure :: broadcast
     procedure :: eval
     procedure :: derivative
     final     :: cleanup
  end type t_bspline2D

  contains
!=======================================================================



!=======================================================================
! allocate memory for spline coefficients and knots
!
! nx, ny               resolution in x and y direction
! nord                 spline interpolation order
!=======================================================================
  subroutine new(this, nx, ny, nord)
  class(t_bspline2D)  :: this
  integer, intent(in) :: nx, ny
  integer, intent(in), optional :: nord


  this%nx   = nx
  this%ny   = ny
  if (present(nord)) then
     this%nord = nord
  else
     this%nord = 5
  endif
  call cleanup(this)


  allocate (this%xnot(nx+this%nord), this%ynot(ny+this%nord))
  allocate (this%dcoeff(nx, ny))

  end subroutine new
!=======================================================================



!=======================================================================
! initialize 2D spline from given data (original data is not stored here)
!
! nx, ny               resolution in x and y direction
! nord                 spline interpolation order
!=======================================================================
  subroutine init(this, nx, ny, x, y, D, nord)
  class(t_bspline2D)       :: this
  integer,      intent(in) :: nx, ny
  real(real64), intent(in) :: x(nx), y(ny), D(nx,ny)
  integer,      intent(in), optional :: nord


  call this%new(nx, ny, nord)

  call dbsnak(nx, x, this%nord, this%xnot)
  call dbsnak(ny, y, this%nord, this%ynot)
  call dbs2in(nx, x, ny, y, D, &
              nx, this%nord, this%nord, this%xnot, this%ynot, this%dcoeff)

  end subroutine init
!=======================================================================



!=======================================================================
  subroutine broadcast(this)
  use parallel
  class(t_bspline2D)   :: this


  call broadcast_inte_s(this%nx)
  call broadcast_inte_s(this%ny)
  call broadcast_inte_s(this%nord)
  if (mype > 0) then
     call this%new(this%nx, this%ny, this%nord)
  endif

  ! broadcast spline coefficients
  call broadcast_real  (this%xnot, this%nx + this%nord)
  call broadcast_real  (this%ynot, this%ny + this%nord)
  call broadcast_real  (this%dcoeff, this%nx*this%ny)

  end subroutine broadcast
!=======================================================================



!=======================================================================
  function eval(this, x) result(d)
  class(t_bspline2D)       :: this
  real(real64), intent(in) :: x(2)
  real(real64)             :: d


  d    = dbs2dr(0,0,x(1),x(2),this%nord,this%nord, &
                this%xnot,this%ynot,this%nx,this%ny,this%dcoeff)

  end function eval
!=======================================================================



!=======================================================================
  function derivative(this, x, mx, my) result(d)
  class(t_bspline2D)       :: this
  real(real64), intent(in) :: x(2)
  integer,      intent(in) :: mx, my
  real(real64)             :: d


  d    = dbs2dr(mx,my,x(1),x(2),this%nord,this%nord, &
                this%xnot,this%ynot,this%nx,this%ny,this%dcoeff)

  end function derivative
!=======================================================================



!=======================================================================
  subroutine cleanup(this)
  type(t_bspline2D) :: this


  if (allocated(this%xnot))   deallocate(this%xnot)
  if (allocated(this%ynot))   deallocate(this%ynot)
  if (allocated(this%dcoeff)) deallocate(this%dcoeff)

  end subroutine cleanup
!=======================================================================

end module bspline2D
