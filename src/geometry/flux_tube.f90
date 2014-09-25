!===============================================================================
! Finite length (Delta Phi) flux tubes with quadrilateral cross-section.
!
! Such flux tubes are the base for field line reconstruction, and this module
! provides some accuracy analysis tools.
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
module flux_tube
  use iso_fortran_env
  use curve2D
  implicit none
  private


  integer, parameter :: nsub = 10

  ! natural coordinates of nodes
  real(real64), parameter :: xi(4) = (/-1,  1,  1,  -1/), eta(4) = (/-1, -1,  1,  1/)

  type, public :: t_cross_section
     ! node coordinates in real space
     real(real64) :: x(4,2)


     ! Toriodal position
     real(real64) :: phi

     ! shape coefficients
     real(real64) :: a0(2), a1(2), a2(2), a3(2)

     contains
     procedure :: setup_shape
     procedure :: generate => generate_cross_section
     procedure :: load
  end type t_cross_section


  type, public :: t_flux_tube
     type(t_curve) :: F(4)
     real(real64), dimension(:), allocatable     :: phi

     ! initial cross-section of flux tube
     type(t_cross_section) :: cs0

     integer   :: nphi, nlength

     contains
     procedure :: generate
     procedure :: plot
     procedure :: trace_fieldline
     procedure :: analyze_fieldline
!     procedure :: analyze_flux_tube
  end type t_flux_tube


  contains
!=======================================================================



!=======================================================================
! Setup shape coefficients a0, a1, a2, a3 for pre-set nodes
!=======================================================================
  subroutine setup_shape (this)
  class(t_cross_section) :: this


  this%a0 = 0.25d0 * ( this%x(1,:) + this%x(2,:) + this%x(3,:) + this%x(4,:))
  this%a1 = 0.25d0 * (-this%x(1,:) + this%x(2,:) + this%x(3,:) - this%x(4,:))
  this%a2 = 0.25d0 * (-this%x(1,:) - this%x(2,:) + this%x(3,:) + this%x(4,:))
  this%a3 = 0.25d0 * ( this%x(1,:) - this%x(2,:) + this%x(3,:) - this%x(4,:))

  end subroutine setup_shape
!=======================================================================



!=======================================================================
! Generate nodes and shape coefficients from basic geometry
! coefficents and distortion measures.
!                                            ^ a2
! alpha1 = alpha2 = 0 provides a rectangle   |
! with positive orientation:                 |-----> a1
!=======================================================================
  subroutine generate_cross_section (this, x0, phi0, a1, alpha1, a2, alpha2, P, theta)
  use math
  class(t_cross_section)   :: this
  real(real64), intent(in) :: x0(2), phi0, a1, alpha1, a2, alpha2, P, theta

  real(real64) :: P13, P23
  integer      :: i


  ! 1. set shape coefficients
  this%phi   = phi0
  this%a0    = x0

  this%a1(1) = a1 * cos(alpha1)
  this%a1(2) = a1 * sin(alpha1)

  this%a2(1) = a2 * cos(alpha2+pi/2.d0)
  this%a2(2) = a2 * sin(alpha2+pi/2.d0)

  P13        = P * cos(theta)
  P23        = P * sin(theta)
  this%a3    = P13*this%a2 - P23*this%a1


  ! 2. set node coordinates
  do i=1,4
     this%x(i,:) = this%a0  +  xi(i)*this%a1  +  eta(i)*this%a2  +  xi(i)*eta(i)*this%a3
  enddo

  end subroutine generate_cross_section
!=======================================================================



!=======================================================================
  subroutine load (this, filename)
  use grid
  use math
  class(t_cross_section)       :: this
  character(len=*), intent(in) :: filename

  type(t_grid) :: G


  ! load initial cross-section of flux tube
  call G%load(filename)
  if ((G%coordinates .ne. CYLINDRICAL)  .or.  (G%fixed_coord .ne. 3)) then
     write (6, *) 'error: cylindrical grid at fixed toroidal angle (grid_id = 2*3) required!'
     stop
  endif
  if (G%n .ne. 4) then
     write (6, *) 'error: exactly 4 nodes required for initial cross-section of flux tube!'
     stop
  endif
  this%x   = G%x(1:4,1:2)
  this%phi = G%x(1,3)
  call this%setup_shape()

  end subroutine load
!=======================================================================



!=======================================================================
! Generate finite flux tube
!
! Input:
!    cs0       initial cross-section
!    nphi      number of toroidal steps
!    nlength   flux tube will be of length 360deg / nlength
!=======================================================================
  subroutine generate(this, cs0, nphi, nlength)
  use fieldline
  use math
  class(t_flux_tube)                :: this
  type(t_cross_section), intent(in) :: cs0
  integer,               intent(in) :: nphi, nlength

  type(t_fieldline) :: F
  real(real64)      :: dphi
  integer           :: j, k


  ! initialize flux tube
  this%cs0     = cs0
  this%nphi    = nphi
  this%nlength = nlength


  ! initialize local variables
  dphi = pi2 / nlength / nphi / nsub


  ! start new flux tube
  if (allocated(this%phi)) deallocate(this%phi)
  allocate (this%phi(0:nphi))
  do j=0,nphi
     this%phi(j) = cs0%phi + j*dphi*nsub
  enddo


  ! trace field lines from all 4 nodes
  do k=1,4
     this%F(k) = this%trace_fieldline(xi(k), eta(k))
  enddo

  end subroutine generate
!=======================================================================



!=======================================================================
  subroutine plot(this, filename)
  use math
  class(t_flux_tube)           :: this
  character(len=*), intent(in) :: filename

  integer, parameter  :: iu = 32

  integer :: j, k, k4


  open  (iu, file=filename)
  do j=0,this%nphi
     write (iu, 1000) this%phi(j) / pi * 180.d0
     do k=0,4
        k4 = mod(k,4) + 1
        write (iu, *) this%F(k4)%x(j,:)
     enddo
     write (iu, *)
  enddo
  close (iu)

 1000 format('# phi [deg] = ', f10.3)
  end subroutine plot
!=======================================================================



!=======================================================================
  function trace_fieldline(this, xi, eta) result(FR)
  use math
  use fieldline
  class(t_flux_tube)       :: this
  real(real64), intent(in) :: xi, eta
  type(t_curve)            :: FR

  type(t_fieldline)        :: F
  real(real64) :: dphi, x(3)
  integer      :: i, j


  call FR%new(this%nphi)
  dphi = pi2 / this%nlength / this%nphi / nsub

  ! initial coordinates in real space
  x(1:2)    = this%cs0%a0  +  xi*this%cs0%a1  +  eta*this%cs0%a2  +  xi*eta*this%cs0%a3
  x(3)      = this%phi(0)
  FR%x(0,:) = x(1:2)

  call F%init(x, dphi, NM_AdamsBashforth4, FL_ANGLE)
  do j=1,this%nphi
     do i=1,nsub
        call F%trace_1step()
     enddo
     FR%x(j,:) = F%rc(1:2)
  enddo

  end function trace_fieldline
!=======================================================================



!=======================================================================
  function analyze_fieldline(this, xi, eta) result(D)
  use math
  use equilibrium
  use dataset
  class(t_flux_tube)       :: this
  real(real64), intent(in) :: xi, eta
  type(t_curve)            :: F
  type(t_dataset)          :: D

  real(real64) :: x(2), x1(2), x2(2), x3(2), x4(2), a0(2), a1(2), a2(2), a3(2)
  real(real64) :: xref(2), dx(2), y(3), ePsi(2), ePol(2)
  integer :: j


  F = this%trace_fieldline(xi, eta)
  call D%new(this%nphi,3)
  do j=1,this%nphi
     x1 = this%F(1)%x(j,:); x2 = this%F(2)%x(j,:); x3 = this%F(3)%x(j,:); x4 = this%F(4)%x(j,:)
     a0 = 0.25d0 * ( x1 + x2 + x3 + x4)
     a1 = 0.25d0 * (-x1 + x2 + x3 - x4)
     a2 = 0.25d0 * (-x1 - x2 + x3 + x4)
     a3 = 0.25d0 * ( x1 - x2 + x3 - x4)

     ! reconstructed field line position
     x  = a0  +  xi*a1  +  eta*a2  +  xi*eta*a3

     ! reference position
     xref   = F%x(j,:)
     y(1:2) = xref
     y(3)   = this%phi(j)
     ePsi   = get_ePsi(y)
     ePol(1) = ePsi(2); ePol(2) = -ePsi(1)

     dx = x - xref
     D%x(j,1) = sqrt(sum(dx**2))
     D%x(j,2) = sum(ePsi*dx)
     D%x(j,3) = sum(ePol*dx)
  enddo

  end function analyze_fieldline
!=======================================================================



!=======================================================================
!  function analyze_flux_tube(this, cs0, nphi, nlength) result(D)
!  use dataset
!  class(t_flux_tube)                :: this
!  type(t_cross_section), intent(in) :: cs0
!  integer,               intent(in) :: nphi, nlength
!  type(t_dataset)                   :: D
!
!  real(real64) :: xi, eta
!
!
!  call this%generate(cs0, nphi, nlength)
!
!  xi  = 0.d0
!  eta = 0.d0
!  D   = this%analyze_fieldline(xi, eta)
!
!  end function analyze_flux_tube
!=======================================================================

end module flux_tube
