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
     procedure :: area
     procedure :: Rcenter, Zcenter
     procedure :: generate_mesh
  end type t_cross_section


  type, public :: t_flux_tube
     type(t_curve) :: F(4)
     real(real64), dimension(:,:), allocatable   :: Bmod
     real(real64), dimension(:), allocatable     :: phi

     ! initial cross-section of flux tube
     type(t_cross_section) :: cs0

     integer   :: nphi, nlength

     contains
     procedure :: new
     procedure :: initialize
     procedure :: generate
     procedure :: plot
     procedure :: get_cross_section
     procedure :: trace_fieldline
     procedure :: analyze_fieldline
!     procedure :: analyze_flux_tube
     procedure :: get_flux_along_tube
     procedure :: sample_pitch
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
  function area(this)
  class(t_cross_section) :: this
  real(real64)           :: area

  real(real64) :: abcd1, abcd2

  abcd1 = (this%x(3,1)-this%x(2,1))*(this%x(2,2)-this%x(1,2)) &
        - (this%x(2,1)-this%x(1,1))*(this%x(3,2)-this%x(2,2))
  abcd2 = (this%x(4,1)-this%x(1,1))*(this%x(3,2)-this%x(4,2)) &
        - (this%x(3,1)-this%x(4,1))*(this%x(4,2)-this%x(1,2))
  area  = 0.5d0 * (abcd1 + abcd2)

  end function area
!=======================================================================
!=======================================================================
  function Rcenter(this)
  class(t_cross_section) :: this
  real(real64) :: Rcenter

  Rcenter = (this%x(1,1)+this%x(2,1)+this%x(3,1)+this%x(4,1)) / 4.d0
  end function Rcenter
!=======================================================================
  function Zcenter(this)
  class(t_cross_section) :: this
  real(real64) :: Zcenter

  Zcenter = (this%x(1,2)+this%x(2,2)+this%x(3,2)+this%x(4,2)) / 4.d0
  end function Zcenter
!=======================================================================



!=======================================================================
  function generate_mesh(this, nxi, neta) result(M)
  use math
  use grid
  class(t_cross_section) :: this
  integer, intent(in)    :: nxi, neta
  type(t_grid)           :: M

  real(real64) :: xi, eta, x(2)
  integer      :: i, j, ig


  write (6, *) 'generating mesh at ', this%phi
  call M%new(CYLINDRICAL, UNSTRUCTURED, 3, nxi, neta, 1, this%phi)
  ig = 0
  do j=0,neta-1
     eta = -1.d0 + 2.d0*j/(neta-1)
     do i=0,nxi-1
        xi = -1.d0 + 2.d0*i/(nxi-1)

        ig = ig + 1
        x  = this%a0  +  xi*this%a1  +  eta*this%a2  +  xi*eta*this%a3
        M%x(ig,1:2) = x
     enddo
  enddo

  end function generate_mesh
!=======================================================================



!=======================================================================
! create new (empty) flux tube
!
! Input:
!    nphi      number of toriodal steps
!=======================================================================
  subroutine new(this, nphi)
  class(t_flux_tube)  :: this
  integer, intent(in) :: nphi

  integer :: i


  write (6, *) 'creating new flux tube'
  write (6, *) 'number of segments: ', nphi
  if (allocated(this%phi)) deallocate(this%phi)
  allocate (this%phi(0:nphi))
  this%nphi = nphi

  do i=1,4
     call this%F(i)%new(nphi)
  enddo

  if (allocated(this%Bmod)) deallocate(this%Bmod)
  allocate (this%Bmod(0:nphi,4))

  end subroutine new
!=======================================================================



!=======================================================================
  subroutine initialize(this, nphi, phi, R, Z)
  class(t_flux_tube)       :: this
  integer, intent(in)      :: nphi
  real(real64), intent(in) :: phi(0:nphi), R(4,0:nphi), Z(4,0:nphi)


  call this%new(nphi)
  write (6, *) 'initializing flux tube not yet implemented!'
  stop

  end subroutine initialize
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
  use bfield
  use math
  class(t_flux_tube)                :: this
  type(t_cross_section), intent(in) :: cs0
  integer,               intent(in) :: nphi, nlength

  type(t_fieldline) :: F
  real(real64)      :: dphi, x(3), Bf(3)
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
  if (allocated(this%Bmod)) deallocate(this%Bmod)
  allocate (this%Bmod(0:nphi,4))


  ! trace field lines from all 4 nodes
  do k=1,4
     this%F(k) = this%trace_fieldline(xi(k), eta(k))
     do j=0,nphi
        x  = this%F(k)%x(j,1:3)
        Bf = get_Bf_cyl(x)
        this%Bmod(j,k) = sqrt(sum(Bf**2))
     enddo
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
  function get_cross_section(this, iphi) result(cs)
  class(t_flux_tube)    :: this
  integer, intent(in)   :: iphi
  type(t_cross_section) :: cs

  integer :: i


  if (iphi < 0  .or.  iphi > this%nphi) then
     write (6, *) 'error in function t_flux_tube%get_cross_section:'
     write (6, *) 'iphi out of range'
     stop
  endif


  cs%phi = this%phi(iphi)
  do i=1,4
     cs%x(i,1) = this%F(i)%x(iphi,1)
     cs%x(i,2) = this%F(i)%x(iphi,2)
  enddo
  call cs%setup_shape()

  end function get_cross_section
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


!=======================================================================
  function get_flux_along_tube(this) result(PHI)
  use dataset
  use math
  class(t_flux_tube) :: this
  type(t_dataset)    :: PHI

  real(real64) :: Bmod2, Bmod1, R2, R1, Z2, Z1, phi2, phi1, area2, area1, DFL, pitch, flux
  type(t_cross_section) :: cs
  integer :: it


  call PHI%new(this%nphi,5)
  do it=0,this%nphi
     cs    = this%get_cross_section(it)
     area2 = cs%area()
     R2    = cs%Rcenter()
     Z2    = cs%Zcenter()
     Phi2  = cs%phi
     Bmod2 = sum(this%Bmod(it,:)) / 4.d0

     if (it/=0) then
        DFL   = 0.5d0*(R1+R2)*abs(Phi2-Phi1)
        pitch = DFL/sqrt(DFL**2+(R2-R1)**2+(Z2-Z1)**2)
        flux  = (area1+area2)*pitch*(Bmod2+Bmod1)*0.25d0
        PHI%x(it,1) = (Phi1 + Phi2) / 2.d0 / pi * 180.d0
        PHI%x(it,2) = flux
        PHI%x(it,3) = pitch
        PHI%x(it,4) = (area1 + area2) / 2.d0
        PHI%x(it,5) = (Bmod1 + Bmod2) / 2.d0
     endif
     area1 = area2; R1 = R2; Z1 = Z2; phi1 = phi2; Bmod1 = Bmod2
  enddo

  end function get_flux_along_tube
!=======================================================================



!=======================================================================
  function sample_pitch(this, it, nxi, neta) result(pitch)
  use dataset
  use grid
  class(t_flux_tube)  :: this
  integer, intent(in) :: it, nxi, neta
  type(t_dataset)     :: pitch
  
  type(t_grid)          :: M1, M2
  type(t_cross_section) :: cs1, cs2
  real(real64) :: R1, R2, Z1, Z2, phi1, phi2, DFL
  integer :: ig, i, j


  if (it < 0  .or.  it >= this%nphi) then
     write (6, *) 'error: index it out of range!'
     stop
  endif

  call pitch%new(nxi*neta, 3)

  cs1  = this%get_cross_section(it)
  M1   = cs1%generate_mesh(nxi, neta)
  phi1 = cs1%phi
  write (6, *) 'cross section 1 at phi = ', phi1
  
  cs2 = this%get_cross_section(it+1)
  M2  = cs2%generate_mesh(nxi, neta)
  phi2 = cs2%phi
  write (6, *) 'cross section 2 at phi = ', phi2

  ig = 0
  do j=1,neta
     do i=1,nxi
        ig = ig + 1

        R1 = M1%x(ig,1)
        Z1 = M1%x(ig,2)
        R2 = M2%x(ig,1)
        Z2 = M2%x(ig,2)

        DFL   = 0.5d0*(R1+R2)*abs(phi2-phi1)
        pitch%x(ig,1) = i
        pitch%x(ig,2) = j
        pitch%x(ig,3) = DFL/sqrt(DFL**2+(R2-R1)**2+(Z2-Z1)**2)
     enddo
  enddo


  end function sample_pitch
!=======================================================================

end module flux_tube
