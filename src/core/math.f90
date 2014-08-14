module math
  implicit none

  real*8, parameter :: pi    = 3.14159265358979323846264338328d0, pi2 = 2.d0 * pi

  integer, parameter :: &
     CARTESIAN   = 1, &
     CYLINDRICAL = 2
  
  contains
!=======================================================================
  


!=======================================================================
! Coordinate transformation
! c = 1: Cartesian coordinates [cm]
! c > 1: Cylindrical coordinates [cm,rad]
!=======================================================================
  subroutine coord_trans (y_in, c_in, y_out, c_out)
  real*8,  intent(in)  :: y_in(3)
  integer, intent(in)  :: c_in
  real*8,  intent(out) :: y_out(3)
  integer, intent(in)  :: c_out


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
! for angular coordinate phi return phi_sym with 0 <= phi_sym <= 2*pi/n_sym
!=======================================================================
  function phi_sym (phi, n_sym)
  real*8, intent(in)  :: phi
  integer, intent(in) :: n_sym
  real*8              :: phi_sym

  real*8  :: Dphi
  integer :: k


  phi_sym = phi
  Dphi    = pi2 / n_sym
  k       = phi_sym / Dphi
  phi_sym = phi_sym - k*Dphi
  if (phi_sym.lt.0) phi_sym = phi_sym + Dphi

  end function phi_sym
!=======================================================================

end module math
