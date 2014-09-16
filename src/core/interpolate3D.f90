!===============================================================================
! Interpolate quantity in a structured cylindrical grid
!===============================================================================
module interpolate3D
  use iso_fortran_env
  use bspline
  implicit none
  private


  type, public :: t_interpolate3D
     ! raw data on nodes
     real(real64), dimension(:,:,:), allocatable :: d
     ! actual node coordinates
     real(real64), dimension(:),     allocatable :: R, Z, Phi

     ! coefficients for spline interpolation
     real(real64), dimension(:,:,:), allocatable :: dcoeff
     real(real64), dimension(:),     allocatable :: Rnot, Znot, Phinot

     ! Rc: center major radius, w: width, h: height
     real(real64) :: Rc, Zc, w, h
     integer      :: nr, nz, nphi, nsym

     ! spline interpolation order
     integer      :: nord

     contains
     procedure :: new, load, store, setup, eval
  end type t_interpolate3D

  contains
!=======================================================================



!=======================================================================
! initialize new dataset for interpolation
!
! nr, nz, nphi         resolution in R, Z and Phi direction
! nsym                 toroidal symmetry number
! Rc, Zc               center of computational box
! w, h                 width and height of computational box
! nord                 spline interpolation order
!=======================================================================
  subroutine new(this, nr, nz, Rc, Zc, w, h, nphi, nsym)
  use math
  class(t_interpolate3D)   :: this
  integer, intent(in)      :: nr, nz, nphi, nsym
  real(real64), intent(in) :: Rc, Zc, w, h

  integer :: i, j, k


  this%nr   = nr
  this%nz   = nz
  this%nphi = nphi
  this%nsym = nsym
  this%nord = 5

  this%Rc   = Rc
  this%Zc   = Zc
  this%w    = w
  this%h    = h


  if (allocated(this%d)) deallocate(this%d)
  allocate (this%d(nr, nz, nphi))

  if (allocated(this%R)) deallocate(this%R, this%Z, this%Phi)
  allocate (this%R(nr), this%Z(nz), this%Phi(nphi))

  if (allocated(this%Rnot)) deallocate(this%Rnot, this%Znot, this%Phinot)
  allocate (this%Rnot  (this%nr  +this%nord), &
            this%Znot  (this%nz  +this%nord), &
            this%Phinot(this%nphi+this%nord))

  do i=1,this%nr
     this%R(i)   = this%Rc + this%w * (-0.5d0 + 1.d0 * (i-1) / (this%nr-1))
  enddo
  call dbsnak (this%nr, this%R, this%nord, this%Rnot)

  do j=1,this%nz
     this%Z(j)   = this%Zc + this%h * (-0.5d0 + 1.d0 * (j-1) / (this%nz-1))
  enddo
  call dbsnak (this%nz, this%Z, this%nord, this%Znot)

  do k=1,this%nphi
     this%Phi(k) = 2.d0*pi / this%nsym * (k-1) / (this%nphi-1)
  enddo
  call dbsnak (this%nphi, this%Phi, this%nord, this%Phinot)

  end subroutine new
!=======================================================================



!=======================================================================
  subroutine load(this, filename)
  class(t_interpolate3D)       :: this
  character(len=*), intent(in) :: filename

  integer, parameter :: iu = 32

  integer      :: nr, nz, nphi, nsym
  real(real64) :: Rc, Zc, w, h


  open  (iu, file=filename)
  read  (iu, *) nr, nz, nphi, nsym
  read  (iu, *) Rc, Zc, w,    h
  call this%new(nr, nz, Rc, Zc, w, h, nphi, nsym)
  read  (iu, *) this%d
  call this%setup()
  close (iu)

  end subroutine load
!=======================================================================



!=======================================================================
  subroutine store(this, filename)
  class(t_interpolate3D)       :: this
  character(len=*), intent(in) :: filename

  integer, parameter :: iu = 32


  open  (iu, file=filename)
  write (iu, *) this%nr, this%nz, this%nphi, this%nsym
  write (iu, *) this%Rc, this%Zc, this%w,    this%h
  write (iu, *) this%d
  close (iu)

  end subroutine store
!=======================================================================



!=======================================================================
  subroutine setup(this)
  class(t_interpolate3D)   :: this


  allocate (this%dcoeff(this%nr, this%nz, this%nphi))
  call dbs3in (this%nr, this%R, this%nz, this%Z, this%nphi, this%Phi, &
               this%d, this%nr, this%nz, this%nord, this%nord, this%nord, &
               this%Rnot, this%Znot, this%Phinot, this%dcoeff)

  end subroutine setup
!=======================================================================



!=======================================================================
  function eval(this, x) result(d)
  use math
  class(t_interpolate3D)   :: this
  real(real64), intent(in) :: x(3)
  real(real64)             :: d

  real(real64) :: phi0


  phi0 = phi_sym(x(3), this%nsym)
  d    = dbs3dr(0,0,0,x(1),x(2),phi0,this%nord,this%nord,this%nord, &
                this%Rnot,this%Znot,this%Phinot,this%nr,this%nz,this%nphi,this%dcoeff)

  end function eval
!=======================================================================

end module interpolate3D
