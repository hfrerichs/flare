!===============================================================================
! Magnetic Field Coils (defined by closed polygons)
!===============================================================================
module mfc_polygon
  use iso_fortran_env
  implicit none
  private


  real(real64), parameter :: min_dist_to_coil = 1.d-8


  type, public :: t_mfc_polygon
     ! number of polygon segments
     integer :: n_seg

     ! vertices
     real(real64), dimension(:), allocatable :: X, Y, Z

     ! working arrays for segments (direction and length)
     real(real64), dimension(:), pointer     :: dx, dy, dz, rp, rp1

     ! current
     real(real64) :: I0

     contains
     procedure :: load
     procedure :: get_A       ! magnetic vector potential [T m]
     procedure :: get_A_cyl   ! magnetic vector potential [T m]
     procedure :: get_Bf      ! magnetic field vector [Gauss]
     procedure :: get_JBf     ! Jacobian [T/m]
     procedure :: get_JBf_cyl ! Jacobian [T/m]
     procedure :: broadcast
  end type t_mfc_polygon


  contains
!===============================================================================



!===============================================================================
  subroutine load(this, iu, I_scale)
  class(t_mfc_polygon)     :: this
  integer, intent(in)      :: iu
  real(real64), intent(in) :: I_scale

  real(real64) :: I0
  integer      :: i, n


  read  (iu, *) n, I0
  this%n_seg = n
  this%I0    = I0 * I_scale
  allocate (this%X(0:n), this%Y(0:n), this%Z(0:n))
  write (6, 1000) n, this%I0 / 1.d3
  read  (iu ,*)  (this%X(i), this%Y(i), this%Z(i), i=0,n)

  allocate (this%dx(0:n), this%dy(0:n), this%dz(0:n), this%rp(0:n), this%rp1(0:n))

 1000 format (8x,i0,' segments,',5x,'Ic = ',f10.3,' kA')
  end subroutine load
!===============================================================================



!===============================================================================
! return magnetic vector potential [T m] at x (Cartesian coordinates)
!===============================================================================
  function get_A(this, x) result(A)
  class(t_mfc_polygon)     :: this
  real(real64), intent(in) :: x(3)
  real(real64)             :: A(3)

  real(real64) :: F, x1(3), x2(3), a1, b1, c1, lnX
  real(real64) :: L, Ri, Rf, eps
  integer      :: is


  A = 0.d0
  F = this%I0 * 1.d-7 ! mu_0 / 4pi
  do is=1,this%n_seg
     x1(1) = this%X(is-1)
     x1(2) = this%Y(is-1)
     x1(3) = this%Z(is-1)
     x2(1) = this%X(is)
     x2(2) = this%Y(is)
     x2(3) = this%Z(is)

     a1    = sum((x2-x1)**2)
     b1    = sum(-2.d0*(x2-x1)*(x-x1))
     c1    = sum((x-x1)**2)

     lnX   = log(((0.5d0*b1+a1)/sqrt(a1) + sqrt(a1+b1+c1)) / (0.5d0*b1/sqrt(a1) + sqrt(c1)))
     A     = A + F*(x2-x1) / sqrt(a1) * lnX

!     L   = sum((x2-x1)**2)
!     Rf  = sum((x-x2)**2)
!     Ri  = sum((x-x1)**2)
!     eps = L / (Ri + Rf)
!     A   = A + F * (x2-x1)/L * log((1.d0 + eps) / (1.d0 - eps))
  enddo

  end function get_A
!===============================================================================



!===============================================================================
! return magnetic vector potential [T m] at r (Cylindrical coordinates)
!===============================================================================
  function get_A_cyl(this, r) result(A)
  class(t_mfc_polygon)     :: this
  real(real64), intent(in) :: r(3)
  real(real64)             :: A(3)

  real(real64) :: x(3), cosphi, sinphi, Acart(3)


  cosphi = cos(r(3))
  sinphi = sin(r(3))
  x(1)   = r(1)*cosphi
  x(2)   = r(1)*sinphi
  x(3)   = r(2)
  Acart  = this%get_A(x)

  A(1)   =  Acart(1)*cosphi + Acart(2)*sinphi
  A(2)   =  Acart(3)
  A(3)   = -Acart(1)*sinphi + Acart(2)*cosphi

  end function get_A_cyl
!===============================================================================



!===============================================================================
! return magnetic field vector [Gauss] at x (Cartesian coordinates)
!===============================================================================
  function get_Bf(this, x) result(Bf)
  class(t_mfc_polygon)     :: this
  real(real64), intent(in) :: x(3)
  real(real64)             :: Bf(3)

  real(real64), dimension(:), pointer :: dx, dy, dz, rp, rp1

  real(real64) :: rp2, rp12, som, s1, ck, ak, rx, ry, rz, var
  integer :: i, i1, n


  dx => this%dx; dy => this%dy; dz => this%dz
  rp => this%rp; rp1 => this%rp1

  n = this%n_seg
  do i=0,n
     dx(i) = x(1) - this%X(i)
     dy(i) = x(2) - this%Y(i)
     dz(i) = x(3) - this%Z(i)
     rp2   = dx(i)**2 + dy(i)**2 + dz(i)**2
     rp2   = max(rp2,min_dist_to_coil**2)
     rp(i) = sqrt(rp2)
     rp1(i)= 1.d0/rp2
  enddo

  Bf  = 0.d0
  var = this%I0 * 0.1d0 ! mu_0 / 4pi * 100 cm/m * 10000 Gauss/T
  do i=1,n
     i1    = i-1
     rp12  = rp(i1) * rp(i)
     som   = rp(i1) + rp(i)
     s1    = 1.d0 /  (rp12 * (rp12 + dx(i1)*dx(i) &
                                   + dy(i1)*dy(i) &
                                   + dz(i1)*dz(i)))
     ck    = var * s1
     ak    = som * ck
     rx    = dy(i1)*dz(i) - dz(i1)*dy(i)
     ry    = dz(i1)*dx(i) - dx(i1)*dz(i)
     rz    = dx(i1)*dy(i) - dy(i1)*dx(i)
     Bf(1) = Bf(1)  +  ak * rx
     Bf(2) = Bf(2)  +  ak * ry
     Bf(3) = Bf(3)  +  ak * rz
  enddo

  end function get_Bf
!===============================================================================



!===============================================================================
  function get_JBf(this, x) result(J)
  class(t_mfc_polygon)     :: this
  real(real64), intent(in) :: x(3)
  real(real64)             :: J(3,3)

  real(real64), dimension(:), pointer :: dx, dy, dz, rp, rp1

  real(real64) :: rp2, rp12, som, s1, ck, ak, rx, ry, rz, var
  real(real64) :: bk, p(6), dkx, dky, dkz
  integer :: i, i1, n


  dx => this%dx; dy => this%dy; dz => this%dz
  rp => this%rp; rp1 => this%rp1

  n = this%n_seg
  do i=0,n
     dx(i) = x(1) - this%X(i)
     dy(i) = x(2) - this%Y(i)
     dz(i) = x(3) - this%Z(i)
     rp2   = dx(i)**2 + dy(i)**2 + dz(i)**2
     rp2   = max(rp2,min_dist_to_coil**2)
     rp(i) = sqrt(rp2)
     rp1(i)= 1.d0/rp2
  enddo

  J   = 0.d0
  var = this%I0 * 0.1d0 ! mu_0 / 4pi * 100 cm/m * 10000 Gauss/T
  do i=1,n
     i1    = i-1
     rp12  = rp(i1) * rp(i)
     som   = rp(i1) + rp(i)
     s1    = 1.d0 /  (rp12 * (rp12 + dx(i1)*dx(i) &
                                   + dy(i1)*dy(i) &
                                   + dz(i1)*dz(i)))
     ck     = var * s1
     ak     = som * ck
     rx     = dy(i1)*dz(i) - dz(i1)*dy(i)
     ry     = dz(i1)*dx(i) - dx(i1)*dz(i)
     rz     = dx(i1)*dy(i) - dy(i1)*dx(i)

     bk     = -ak *som * s1
     p(1)   = rp(i1) * dx(i)
     p(2)   = rp(i)  * dx(i1)
     dkx    = bk * (p(1)+p(2))  -  ck * (p(2)*rp1(i1) + p(1)*rp1(i))
     p(3)   = rp(i1) * dy(i)
     p(4)   = rp(i)  * dy(i1)
     dky    = bk * (p(3)+p(4))  -  ck * (p(4)*rp1(i1) + p(3)*rp1(i))
     p(5)   = rp(i1) * dz(i)
     p(6)   = rp(i)  * dz(i1)
     dkz    = bk * (p(5)+p(6))  -  ck * (p(6)*rp1(i1) + p(5)*rp1(i))

!     Bf1(4)  = Bf1(4)  +  dkx*rx
!     Bf1(5)  = Bf1(5)  +  dky*rx + ak*(dz(is)  - dz(is1))
!     Bf1(6)  = Bf1(6)  +  dkz*rx + ak*(dy(is1) - dy(is))
!     Bf1(7)  = Bf1(7)  +  dky*ry
!     Bf1(8)  = Bf1(8)  +  dkz*ry + ak*(dx(is)  - dx(is1))

     ! Jacobian is symmetric for current filaments as long as x is not on the filament itself
     J(1,1) = J(1,1)  +  dkx*rx
     J(1,2) = J(1,2)  +  dky*rx + ak*(dz(i)  - dz(i1))
     J(1,3) = J(1,3)  +  dkz*rx + ak*(dy(i1) - dy(i))
     J(2,1) = J(2,1)  +  dky*rx + ak*(dz(i)  - dz(i1))
     J(2,2) = J(2,2)  +  dky*ry
     J(2,3) = J(2,3)  +  dkz*ry + ak*(dx(i)  - dx(i1))
     J(3,1) = J(3,1)  +  dkz*rx + ak*(dy(i1) - dy(i))
     J(3,2) = J(3,2)  +  dkz*ry + ak*(dx(i)  - dx(i1))
  enddo
  J(3,3) = - J(1,1) - J(2,2)
  J      = J / 1.d2 ! Gauss/cm -> T/m

  end function get_JBf
!===============================================================================



!===============================================================================
  function get_JBf_cyl(this, r) result(J)
  class(t_mfc_polygon)     :: this
  real(real64), intent(in) :: r(3)
  real(real64)             :: J(3,3)

  real(real64) :: R1, x(3), Jcart(3,3), B(3), cosphi, sinphi, dBrdx, dBrdy, dBpdx, dBpdy


  R1     = r(1) / 1.e2 ! cm -> m
  cosphi = cos(r(3))
  sinphi = sin(r(3))
  x(1)   = r(1)*cosphi
  x(2)   = r(1)*sinphi
  x(3)   = r(2)
  Jcart  = this%get_JBf(x)       ! [T/m]
  B      = this%get_Bf(x) / 1.d4 ! [T]

  dBrdx  = Jcart(1,1)*cosphi + B(1)/R1*sinphi**2 + Jcart(2,1)*sinphi - B(2)/R1*sinphi*cosphi
  dBrdy  = Jcart(1,2)*cosphi - B(1)/R1*sinphi*cosphi + Jcart(2,2)*sinphi + B(2)/R1*cosphi**2
  dBpdx  = -Jcart(1,1)*sinphi + B(1)/R1*sinphi*cosphi + Jcart(2,1)*cosphi + B(2)/R1*sinphi**2
  dBpdy  = -Jcart(1,2)*sinphi - B(1)/R1*cosphi**2 + Jcart(2,2)*cosphi - B(2)/R1*sinphi*cosphi

  ! Br
  J(1,1) = dBrdx * cosphi + dBrdy * sinphi
  J(1,2) = Jcart(1,3)*cosphi + Jcart(2,3)*sinphi
  J(1,3) = R1*(     -dBrdx*sinphi +      dBrdy*cosphi)

  ! Bz
  J(2,1) =      Jcart(3,1)*cosphi + Jcart(3,2)*sinphi
  J(2,2) =      Jcart(3,3)
  J(2,3) = R1*(-Jcart(3,1)*sinphi + Jcart(3,2)*cosphi)

  ! Bphi
  J(3,1) = dBpdx * cosphi + dBpdy * sinphi
  J(3,2) = -Jcart(1,3)*sinphi + Jcart(2,3)*cosphi
  J(3,3) = R1*(     -dBpdx*sinphi +      dBpdy*cosphi)

  end function get_JBf_cyl
!===============================================================================



!===============================================================================
  subroutine broadcast(this)
  use parallel
  class(t_mfc_polygon) :: this

  integer :: n


  call broadcast_real_s(this%I0)
  call broadcast_inte_s(this%n_seg)
  n = this%n_seg
  if (mype > 0) then
     allocate (this%X(0:n), this%Y(0:n), this%Z(0:n))
  endif
  call broadcast_real  (this%X, n+1)
  call broadcast_real  (this%Y, n+1)
  call broadcast_real  (this%Z, n+1)

  end subroutine broadcast
!===============================================================================

end module mfc_polygon
