!===============================================================================
! B-spline interpolation in 2D (on structured grid)
!===============================================================================
module bspline2D
  use iso_fortran_env
  use bspline
  implicit none
  private


  type, public :: t_bspline2D
     ! raw data on nodes
     real(real64), dimension(:,:), allocatable :: d
     ! node coordinates
     real(real64), dimension(:),   allocatable :: x, y

     ! coefficients for spline interpolation
     real(real64), dimension(:,:), allocatable :: dcoeff
     real(real64), dimension(:),   allocatable :: xnot, ynot

     ! Rc: center major radius, w: width, h: height
     !real(real64) :: Rc, Zc, w, h
     integer      :: nx, ny

     ! spline interpolation order
     integer      :: nord

     contains
     procedure :: new
     !procedure :: default_geometry
     !procedure :: load
     !procedure :: store
     procedure :: broadcast
     procedure :: setup
     procedure :: eval
  end type t_bspline2D

  contains
!=======================================================================



!=======================================================================
! initialize new dataset for interpolation
!
! nx, ny               resolution in x and y direction
! nord                 spline interpolation order
!=======================================================================
  subroutine new(this, nx, ny, nord, x, y, D)
  class(t_bspline2D)  :: this
  integer, intent(in) :: nx, ny
  integer, intent(in), optional :: nord
  real(real64), intent(in), optional :: x(nx), y(ny), D(nx,ny)

  integer :: i, j, k


  this%nx   = nx
  this%ny   = ny
  if (present(nord)) then
     this%nord = nord
  else
     this%nord = 5
  endif


  if (allocated(this%d)) deallocate(this%d)
  allocate (this%d(nx, ny))

  if (allocated(this%x)) deallocate(this%x, this%y)
  allocate (this%x(nx), this%y(ny))

  if (allocated(this%xnot)) deallocate(this%xnot, this%ynot)
  allocate (this%xnot  (this%nx  +this%nord), &
            this%ynot  (this%ny  +this%nord))

  if (allocated(this%dcoeff)) deallocate(this%dcoeff)
  allocate (this%dcoeff(this%nx, this%ny))

  if (present(x)) this%x = x
  if (present(y)) this%y = y
  if (present(D)) this%d = D

  end subroutine new
!=======================================================================


!=======================================================================
!  subroutine default_geometry

!  do i=1,this%nr
!     this%R(i)   = this%Rc + this%w * (-0.5d0 + 1.d0 * (i-1) / (this%nr-1))
!  enddo

!  do j=1,this%nz
!     this%Z(j)   = this%Zc + this%h * (-0.5d0 + 1.d0 * (j-1) / (this%nz-1))
!  enddo

!  do k=1,this%nphi
!     this%Phi(k) = 2.d0*pi / this%nsym * (k-1) / (this%nphi-1)
!  enddo

!  end subroutine default_geometry
!=======================================================================



!=======================================================================
!  subroutine load(this, filename)
!  class(t_interpolate3D)       :: this
!  character(len=*), intent(in) :: filename
!
!  integer, parameter :: iu = 32
!
!  integer      :: i, nr, nz, nphi, nsym
!  real(real64) :: Rc, Zc, w, h
!
!
!  open  (iu, file=filename)
!  read  (iu, *) nr, nz, nphi, nsym
!  read  (iu, *) Rc, Zc, w,    h
!  call this%new(nr, nz, Rc, Zc, w, h, nphi, nsym)
!  do i=1,this%nphi
!     read  (iu, *) this%d(:,:,i)
!  enddo
!  call this%setup()
!  close (iu)
!
!  end subroutine load
!=======================================================================



!=======================================================================
!  subroutine store(this, filename)
!  class(t_interpolate3D)       :: this
!  character(len=*), intent(in) :: filename
!
!  integer, parameter :: iu = 32
!  integer :: i
!
!
!  open  (iu, file=filename)
!  write (iu, *) this%nr, this%nz, this%nphi, this%nsym
!  write (iu, *) this%Rc, this%Zc, this%w,    this%h
!  do i=1,this%nphi
!     write (iu, *) this%d(:,:,i)
!  enddo
!  close (iu)
!
!  end subroutine store
!=======================================================================



!=======================================================================
  subroutine setup(this)
  class(t_bspline2D)   :: this


  call dbsnak (this%nx, this%x, this%nord, this%xnot)
  call dbsnak (this%ny, this%y, this%nord, this%ynot)

  call dbs2in (this%nx, this%x, this%ny, this%y, this%d, &
               this%nx, this%nord, this%nord, &
               this%xnot, this%ynot, this%dcoeff)

  end subroutine setup
!=======================================================================



!=======================================================================
  subroutine broadcast(this)
  use parallel
  class(t_bspline2D)   :: this


  call broadcast_inte_s(this%nx)
  call broadcast_inte_s(this%ny)
  call broadcast_inte_s(this%nord)

  ! broadcast raw data
  call broadcast_real  (this%x, this%nx)
  call broadcast_real  (this%y, this%ny)
  call broadcast_real  (this%d, this%nx*this%ny)

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

end module bspline2D
