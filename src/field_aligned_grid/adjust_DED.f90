subroutine adjust_DED(C1)
  use iso_fortran_env
  use math
  use curve2D
  implicit none
  type(t_curve), intent(inout) :: C1

  real(real64), parameter :: &
     R_edge  = 153.34465d0, & ! R [cm] of upper DED-edge
     Z_edge  = 42.501011d0, & ! Z [cm] of upper DED-edge
     rDEDlim = 49.1d0,      & ! limiting minor radius [cm] behind DED plates
     R0      = 175.d0,      & ! Major radius [cm]    of Textor vessel center
     Z0      =   0.d0,      & ! Vertical coord. [cm] of Textor vessel center
     Dtheta1 = 7.d0,        & ! Pol. increment [deg] from DED-edge: begin adjustment
     Dtheta2 = 10.d0,       & !                                       end adjustment
     theta0  = pi - atan(Z_edge, R0 - R_edge), &
     thetaS1 = theta0 + Dtheta1/180.d0*pi, &
     thetaS2 = theta0 + Dtheta2/180.d0*pi


  real(real64) :: y0(2), y1(2), theta, r
  integer :: i


  do i=0,C1%n_seg
     y1(1) = C1%x(i,1)
     y1(2) = C1%x(i,2)
     theta = atan2(y1(2) - Z0, y1(1) - R0)
     if (theta < 0) theta = theta + pi2

     ! r: original minor radius
     r = sqrt((y1(1)-R0)**2 + (y1(2)-Z0)**2)


     ! full adjustment
     if (theta > thetaS2  .and.  theta < pi2-thetaS2) then
        r = rDEDlim

     ! partial adjustment on upper edge
     elseif (theta > thetaS1  .and.  theta <= thetaS2) then
        r = r + (rDEDlim-r) * (theta-thetaS1) / (thetaS2-thetaS1)

     ! partial adjustment on lower edge
     elseif (theta < pi2-thetaS1  .and.  theta >= pi2-thetaS2) then
        r = r + (rDEDlim-r) * (pi2-theta-thetaS1) / (thetaS2-thetaS1)

     endif


     ! set new coordinates
     y1(1) = R0 + r*cos(theta)
     y1(2) = Z0 + r*sin(theta)
     C1%x(i,:) = y1
  enddo

end subroutine adjust_DED
