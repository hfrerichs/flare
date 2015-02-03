module math
  use iso_fortran_env
  implicit none

  real(real64), parameter :: pi    = 3.14159265358979323846264338328d0, pi2 = 2.d0 * pi

  integer, parameter :: &
     CARTESIAN   = 1, &
     CYLINDRICAL = 2, &
     TORUS       = 4      ! torus coordinates: minor radius, poloidal angle, toroidal angle
  
  integer, parameter :: &
     GEOANGLE    = 1, &
     ARCLENGTH   = 2

  character(len=12), dimension(2), parameter :: &
     COORDINATES = (/ 'cartesian   ', &
                      'cylindrical ' /)

  contains
!=======================================================================
  


!=======================================================================
! Coordinate transformation
! c = 1: Cartesian coordinates [cm]
! c > 1: Cylindrical coordinates [cm,rad]
!=======================================================================
  subroutine coord_trans (y_in, c_in, y_out, c_out)
  real(real64), intent(in)  :: y_in(3)
  integer,      intent(in)  :: c_in
  real(real64), intent(out) :: y_out(3)
  integer,      intent(in)  :: c_out


  if (c_in == 1 .and. c_out == 1) then
     y_out = y_in
     return
  endif

  if (c_in > 1 .and. c_out > 1) then
     y_out = y_in
     return
  endif

  ! input is in Cartesian coordinates, output in cylindrical coordinates
  if (c_in == 1) then
     y_out(1) = sqrt(y_in(1)**2 +  y_in(2)**2)
     y_out(2) = y_in(3)
     y_out(3) = atan2(y_in(2), y_in(1))

  ! input is in cylindrical coordinates, output in Cartesian coordinates
  else
     y_out(1) = y_in(1) * cos(y_in(3))
     y_out(2) = y_in(1) * sin(y_in(3))
     y_out(3) = y_in(2)
  endif


  end subroutine coord_trans
!=======================================================================



!=======================================================================
! Transformation from torus coordinates to cylindrical coordinates
! Input:
!    y(1):     Minor radius [cm]
!    y(2):     Poloidal angle [rad]
!    y(3):     Toroidal angle [rad]
!    R0:       Reference Major radius [cm]
!
! Output:
!    r(1:2):   R [cm], Z [cm] coordinate
!    r(3) = y(3)
!=======================================================================
  subroutine coord_trans_torus (y, R0, r)
  real(real64), intent(in)  :: y(3), R0
  real(real64), intent(out) :: r(3)


  r(1) = R0 + y(1) * cos(y(2))
  r(2) =      y(1) * sin(y(2))
  r(3) = y(3)

  end subroutine coord_trans_torus
!=======================================================================



!=======================================================================
! for angular coordinate phi return phi_sym with 0 <= phi_sym <= 2*pi/n_sym
!=======================================================================
  function phi_sym (phi, n_sym)
  real(real64), intent(in) :: phi
  integer,      intent(in) :: n_sym
  real(real64)             :: phi_sym

  real(real64) :: Dphi
  integer      :: k


  Dphi = pi2 / n_sym
  phi_sym = mod(phi, Dphi)
  if (phi_sym < 0.d0) phi_sym = phi_sym + Dphi

  end function phi_sym
!=======================================================================

end module math
