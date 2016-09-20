!===============================================================================
! This file contains a collection of modules provided by Y. Suzuki
!===============================================================================


!===============================================================================
!=module.f90
!
!==Version
!
! $Revision: $
! $Id: $
!
!==Overview
!
!==Reference
!
!==Error Handlings
!
!==Known Bugs
!
!==Note
!
!==TODO
!

MODULE kind_spec
  INTEGER, PARAMETER :: DP =  8
END MODULE kind_spec

MODULE param1
  USE kind_spec
  IMPLICIT NONE
  REAL(DP), PARAMETER :: pi  =  3.141592653589793238462643383279502884197_DP, &
    &                    pi2 =  pi + pi
END MODULE param1

MODULE fline_mod
  USE kind_spec
  IMPLICIT NONE
  CHARACTER(LEN=20) :: mode =  '2D'
  LOGICAL :: lupdown, &
    &        lrtheta, &
    &        lflxqnt, &
    &        lpoin,   &
    &        ltext,   &
    &        lfigout, &
    &        lcontb
  INTEGER :: mr     =  30,  &
    &        nstep  =  20,  &
    &        mcirc  =  100, &
    &        ntheta =  150
  REAL(DP) :: h_in        =  0.01_DP,    &
    &         lc_in       =  1.0e+03_DP, &
    &         drflx       =  0.02_DP,    &
    &         dzflx       =  0.0_DP,     &
    &         dpflx       =  0.0_DP,     &
    &         rstart      =  0.0_DP,     &
    &         zstart      =  0.0_DP,     &
    &         pstart      =  0.0_DP,     &
    &         sstart      =  0.0_DP,     &
    &         sout        =  0.99_DP,    &
    &         rout        =  0.0_DP,     &
    &         zout        =  0.0_DP,     &
    &         pcros_in(8) =  0.0_DP
END MODULE fline_mod

MODULE axis_mod
  USE kind_spec
  IMPLICIT NONE
  REAL(DP) :: rax0,   &
    &         zax0,   &
    &         raxis,  &
    &         zaxis,  &
    &         bpaxis, &
    &         baxis,  &
    &         paxis
END MODULE axis_mod
!===============================================================================



!===============================================================================
!=spline_mod.f90
!
!==Version
!
! $Revision: $
! $Id: $
!
!==Overview
!
!    This subroutine is written by Yasuhiro Suzuki
!        at National Institute for Fusion Science ( NIFS )
!         2010/12/12
!
!    Based program is wrtten by T.Watanabe
!        at NIFS
!         1990/8/29
!
!  Multi-Dimensional Spline Interpolation  ( 4th )
!
!  This module contains following public routines.
!
!   1. splin1 : Initializing 1-dimensional interpolation.
!   2. splin2 : Initializing 2-dimensional interpolation.
!   3. splin3 : Initializing 3-dimensional interpolation.
!   4. sp1df  : 1-dimensional interpolation. 
!   5. sp1dd  : 1-dimensional interpolation and derivation along X. 
!   6. sp2df  : 2-dimensional interpolation. 
!   7. sp2dd  : 2-dimensional interpolation and derivation along X and Y.
!   8. sp3df  : 3-dimensional interpolation.
!   9. sp3dd  : 3-dimensional interpolation and derivation along X, Y and Z.
!
!
!              xsc                          xlc
!               <------------- x ------------>
!      1--+--+--4-------------------------nxxm-3-+--+--nxxm  a(i)
!      +--------+----------------------------+---------+
!               <---------------------------->
!
!
!
!==Reference
!
!  "Multi-Dimension Highly Accurate Spline Interpolation Method"
!  by T.Watanabe  The Japan Society for Industrial and Applied Mathmatics
!
!==Error Handlings
!
!==Known Bugs
!
!==Note
!
!==TODO
!
MODULE spline_mod

  USE kind_spec

  IMPLICIT NONE

  PRIVATE

  INTEGER :: l1d =  1, &
    &        l2d =  1, &
    &        l3d =  3, &
    &        nx1d,     &
    &        nx2d,     &
    &        ny2d,     &
    &        nx3d,     &
    &        ny3d,     &
    &        nz3d
   
  REAL(DP), PARAMETER :: c411 = -5.0_DP      / 2048.0_DP,   &
    &                    c412 =  9611.0_DP   / 737280.0_DP, &
    &                    c413 =  259.0_DP    / 23040.0_DP,  &
    &                    c414 = -6629.0_DP   / 46080.0_DP,  &
    &                    c415 = -7.0_DP      / 1152.0_DP,   &
    &                    c416 =  819.0_DP    / 1024.0_DP,   &
    &                    c417 =  1.0_DP      / 1440.0_DP,   &
    &                    c418 = -1067.0_DP   / 360.0_DP,    &
    &                    c41a =  3941.0_DP   / 576.0_DP,    &
    &                    c41c = -1603.0_DP   / 180.0_DP,    &
    &                    c41e =  901.0_DP    / 180.0_DP,    &
    &                    c421 =  49.0_DP     / 2048.0_DP,   &
    &                    c422 = -70733.0_DP  / 737280.0_DP, &
    &                    c423 = -499.0_DP    / 4608.0_DP,   &
    &                    c424 =  47363.0_DP  / 46080.0_DP,  &
    &                    c425 =  59.0_DP     / 1152.0_DP,   &
    &                    c426 = -86123.0_DP  / 15360.0_DP,  &
    &                    c427 = -1.0_DP      / 288.0_DP,    &
    &                    c431 = -245.0_DP    / 2048.0_DP,   &
    &                    c432 =  27759.0_DP  / 81920.0_DP,  &
    &                    c433 =  1299.0_DP   / 2560.0_DP,   &
    &                    c434 = -50563.0_DP  / 15360.0_DP,  &
    &                    c435 = -15.0_DP     / 128.0_DP,    &
    &                    c436 =  51725.0_DP  / 3072.0_DP,   &
    &                    c437 =  1.0_DP      / 160.0_DP,    &
    &                    c441 =  1225.0_DP   / 2048.0_DP,   &
    &                    c442 = -240077.0_DP / 147456.0_DP, &
    &                    c443 = -1891.0_DP   / 4608.0_DP,   &
    &                    c444 =  52931.0_DP  / 9216.0_DP,   &
    &                    c445 =  83.0_DP     / 1152.0_DP,   &
    &                    c446 = -86251.0_DP  / 3072.0_DP,   &
    &                    c447 = -1.0_DP      / 288.0_DP,    &
    &                    d413 =  c413 * 2.0_DP,                &
    &                    d423 =  c423 * 2.0_DP,                &
    &                    d433 =  c433 * 2.0_DP,                &
    &                    d443 =  c443 * 2.0_DP,                &
    &                    d414 =  c414 * 3.0_DP,                &
    &                    d424 =  c424 * 3.0_DP,                &
    &                    d434 =  c434 * 3.0_DP,                &
    &                    d444 =  c444 * 3.0_DP,                &
    &                    d415 =  c415 * 4.0_DP,                &
    &                    d425 =  c425 * 4.0_DP,                &
    &                    d435 =  c435 * 4.0_DP,                &
    &                    d445 =  c445 * 4.0_DP,                &
    &                    d416 =  c416 * 5.0_DP,                &
    &                    d426 =  c426 * 5.0_DP,                &
    &                    d436 =  c436 * 5.0_DP,                &
    &                    d446 =  c446 * 5.0_DP,                &
    &                    d417 =  c417 * 6.0_DP,                &
    &                    d427 =  c427 * 6.0_DP,                &
    &                    d437 =  c437 * 6.0_DP,                &
    &                    d447 =  c447 * 6.0_DP,                &
    &                    d418 =  c418 * 7.0_DP,                &
    &                    d41a =  c41a * 9.0_DP,                &
    &                    d41c =  c41c * 11.0_DP,               &
    &                    d41e =  c41e * 13.0_DP

  REAL(DP) :: h1x, &
    &         h2x, &
    &         h2y, &
    &         h3x, &
    &         h3y, &
    &         h3z, &
    &         xs1, &
    &         xs2, &
    &         ys2, &
    &         xs3, &
    &         ys3, &
    &         zs3, &
    &         xl1, &
    &         xl2, &
    &         yl2, &
    &         xl3, &
    &         yl3, &
    &         zl3
  REAL(DP), ALLOCATABLE :: f1d(:,:),    &
    &                      f2d(:,:,:),  &
    &                      f3d(:,:,:,:)
  !$acc declare create (f3d)

  PUBLIC :: l1d,    &
    &       l2d,    &
    &       l3d,    &
    &       nx1d,   &
    &       nx2d,   &
    &       ny2d,   &
    &       nx3d,   &
    &       ny3d,   &
    &       nz3d,   &
    &       f1d,    &
    &       f2d,    &
    &       f3d,    &
    &       splin1, &
    &       splin2, &
    &       splin3, &
    &       spl1df, &
    &       spl1dd, &
    &       spl2df, &
    &       spl2dd, &
    &       spl3df, &
    &       spl3dd

CONTAINS

  SUBROUTINE splin1 (xsd, xld & ! (in)
    &               )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: xsd, &
      &                     xld


    xs1 =  xsd
    xl1 =  xld
    h1x = (xld - xsd) / (nx1d - 1)


  END SUBROUTINE splin1

  SUBROUTINE splin2 (xsd, xld, ysd, yld & ! (in)
    &               )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: xsd, &
      &                     xld, &
      &                     ysd, &
      &                     yld


    xs2 =  xsd
    xl2 =  xld
    h2x = (xld - xsd) / (nx2d - 1)
    ys2 =  ysd
    yl2 =  yld
    h2y = (yld - ysd) / (ny2d - 1)


  END SUBROUTINE splin2

  SUBROUTINE splin3 (xsd, xld, ysd, yld, zsd, zld & ! (in)
    &               )
    
    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: xsd, &
      &                     xld, &
      &                     ysd, &
      &                     yld, &
      &                     zsd, &
      &                     zld


    xs3 =  xsd  
    xl3 =  xld  
    h3x = (xld - xsd) / (nx3d - 1)
    ys3 =  ysd  
    yl3 =  yld 
    h3y = (yld - ysd) / (ny3d - 1)
    zs3 =  zsd  
    zl3 =  zld  
    h3z = (zld - zsd) / (nz3d - 1)


  END SUBROUTINE splin3

  SUBROUTINE spl1df (xd, & !(in) 
    &                w0  & !(out)
    &               )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN)  :: xd
    REAL(DP), INTENT(OUT) :: w0(l1d)
!Local variables
    INTEGER, PARAMETER :: m2 =  8
    INTEGER :: l, &
      &        ix
    REAL(DP) :: x,    &
      &         ux,   &
      &         usx,  &
      &         x400, &
      &         x41m, &
      &         x41p, &
      &         x42m, &
      &         x42p, &
      &         x43m, &
      &         x43p, &
      &         x44m, &
      &         x44p


    w0(:) =  0.0_DP


    x =  xd
    IF(x < xs1) RETURN
    IF(x > xl1) RETURN

    ux = (x - xs1) / h1x
    ix =  ux
    IF(ix >= nx1d - 1)THEN
      ix =  nx1d - 2
      ux =  nx1d - 1
    END IF

    ux    =  ux - ix - 0.5_DP
    usx   =  ux * ux

!###########       case for m == 4       #############################

    x400 = (((c41e * usx + c41c) * usx + c41a) * usx + c418) * usx

    x41m = (((x400       + c416) * usx + c414) * usx + c412) * ux
    x41p =  ((c417 * usx + c415) * usx + c413) * usx + c411

    x42m = (((-x400 * 7.0_DP + c426) * usx + c424) * usx + c422) * ux
    x42p =  ((c427  * usx    + c425) * usx + c423) * usx + c421

    x43m = (((x400 * 21.0_DP + c436) * usx + c434) * usx + c432) * ux
    x43p =  ((c437 * usx     + c435) * usx + c433) * usx + c431

    x44m = (((-x400 * 35.0_DP + c446) * usx + c444) * usx + c442) * ux
    x44p =  ((c447  * usx     + c445) * usx + c443) * usx + c441

    loop010 : DO l=1,l1d
      w0(l) = ((x41p + x41m) * f1d(l,ix-2)   &
        &   +  (x42p + x42m) * f1d(l,ix-1)   &
        &   +  (x43p + x43m) * f1d(l,ix  )   &
        &   +  (x44p + x44m) * f1d(l,ix+1)   &
        &   +  (x44p - x44m) * f1d(l,ix+2)   &
        &   +  (x43p - x43m) * f1d(l,ix+3)   &
        &   +  (x42p - x42m) * f1d(l,ix+4)   &
        &   +  (x41p - x41m) * f1d(l,ix+5))
    END DO loop010


  END SUBROUTINE spl1df

  SUBROUTINE spl1dd (xd,    &!(in)
    &                w0, wx &!(out)
    &               )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: xd(2)
    REAL(DP), INTENT(OUT) :: w0(l1d), &
      &                      wx(l1d)
!Local variables
    INTEGER, PARAMETER :: m2 =  8
    INTEGER :: l, &
      &        ix
    REAL(DP) :: x,     &
      &         ux,    &
      &         usx,   &
      &         x400,  &
      &         x41m,  &
      &         x41p,  &
      &         x42m,  &
      &         x42p,  &
      &         x43m,  &
      &         x43p,  &
      &         x44m,  &
      &         x44p,  &
      &         dx400, &
      &         dx41m, &
      &         dx41p, &
      &         dx42m, &
      &         dx42p, &
      &         dx43m, &
      &         dx43p, &
      &         dx44m, &
      &         dx44p


    w0(:) =  0.0_DP
    wx(:) =  0.0_DP


    x =  xd(1)
    IF(x < xs1) RETURN
    IF(x > xl1) RETURN

    ux = (x - xs1) / h1x
    ix =  ux
    IF(ix >= nx1d - 1)THEN
      ix =  nx1d - 2
      ux =  nx1d - 1
    END IF

    ux    =  ux - ix - 0.5_DP
    usx   =  ux * ux

!############       case for m == 4       ###########################

    x400  = (((c41e * usx + c41c) * usx + c41a) * usx + c418) * usx
    dx400 = (((d41e * usx + d41c) * usx + d41a) * usx + d418) * usx

    x41m  = (((x400       + c416) * usx + c414) * usx + c412) * ux
    x41p  =  ((c417 * usx + c415) * usx + c413) * usx + c411
    dx41m =  ((dx400      + d416) * usx + d414) * usx + c412
    dx41p =  ((d417 * usx + d415) * usx + d413) * ux

    x42m  = (((-x400  * 7.0_DP  + c426) * usx + c424) * usx + c422) * ux
    x42p  =  ((c427   * usx     + c425) * usx + c423) * usx + c421
    dx42m =  ((-dx400 * 7.0_DP  + d426) * usx + d424) * usx + c422
    dx42p =  ((d427   * usx     + d425) * usx + d423) * ux

    x43m  = (((x400  * 21.0_DP + c436) * usx + c434) * usx + c432) * ux
    x43p  =  ((c437  * usx     + c435) * usx + c433) * usx + c431
    dx43m =  ((dx400 * 21.0_DP + d436) * usx + d434) * usx + c432
    dx43p =  ((d437  * usx     + d435) * usx + d433) * ux

    x44m  = (((-x400  * 35.0_DP + c446) * usx + c444) * usx + c442) * ux
    x44p  =  ((c447   * usx     + c445) * usx + c443) * usx + c441
    dx44m =  ((-dx400 * 35.0_DP + d446) * usx + d444) * usx + c442
    dx44p =  ((d447   * usx     + d445) * usx + d443) * ux

    loop010 : DO l=1,l1d
      w0(l) = ((x41p + x41m) * f1d(l,ix-2)   &
        & +    (x42p + x42m) * f1d(l,ix-1)   &
        & +    (x43p + x43m) * f1d(l,ix  )   &
        & +    (x44p + x44m) * f1d(l,ix+1)   &
        & +    (x44p - x44m) * f1d(l,ix+2)   &
        & +    (x43p - x43m) * f1d(l,ix+3)   &
        & +    (x42p - x42m) * f1d(l,ix+4)   &
        & +    (x41p - x41m) * f1d(l,ix+5))

      wx(l) = ((dx41p + dx41m) * f1d(l,ix-2)  &
        & +    (dx42p + dx42m) * f1d(l,ix-1)  &
        & +    (dx43p + dx43m) * f1d(l,ix  )  &
        & +    (dx44p + dx44m) * f1d(l,ix+1)  &
        & +    (dx44p - dx44m) * f1d(l,ix+2)  &
        & +    (dx43p - dx43m) * f1d(l,ix+3)  &
        & +    (dx42p - dx42m) * f1d(l,ix+4)  &
        & +    (dx41p - dx41m) * f1d(l,ix+5))

      wx(l) =  wx(l) / h1x

    END DO loop010


  END SUBROUTINE spl1dd

  SUBROUTINE spl2df (xd, & !(in) 
    &                w0  & !(out)
    &               )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN)  :: xd(2)
    REAL(DP), INTENT(OUT) :: w0(l2d)
!Local variables
    INTEGER, PARAMETER :: m2 =  8
    INTEGER :: l,  &
      &        my, &
      &        ly, &
      &        ix, &
      &        iy
    REAL(DP) :: x,     &
      &         y,     &
      &         ux,    &
      &         uy,    &
      &         usx,   &
      &         usy,   &
      &         x400,  &
      &         y400,  &
      &         x41m,  &
      &         y41m,  &
      &         x41p,  &
      &         y41p,  &
      &         x42m,  &
      &         y42m,  &
      &         x42p,  &
      &         y42p,  &
      &         x43m,  &
      &         y43m,  &
      &         x43p,  &
      &         y43p,  &
      &         x44m,  &
      &         y44m,  &
      &         x44p,  &
      &         y44p,  &
      &         cy(8)


    w0(:) =  0.0_DP


    x =  xd(1)
    IF(x < xs2) RETURN
    IF(x > xl2) RETURN

    y =  xd(2)   
    IF(y < ys2) RETURN
    IF(y > yl2) RETURN

    ux = (x - xs2) / h2x
    ix =  ux
    IF(ix >= nx2d - 1)THEN
      ix =  nx2d - 2
      ux =  nx2d - 1
    END IF
 
    uy = (y - ys2) / h2y
    iy =  uy
    IF(iy >= ny2d - 1)THEN
      iy =  ny2d - 2
      uy =  ny2d - 1
    END IF

    ux    =  ux - ix - 0.5_DP
    uy    =  uy - iy - 0.5_DP
    usx   =  ux * ux
    usy   =  uy * uy

!###########       case for m == 4       #############################

    x400 = (((c41e * usx + c41c) * usx + c41a) * usx + c418) * usx
    y400 = (((c41e * usy + c41c) * usy + c41a) * usy + c418) * usy
    
    x41m = (((x400       + c416) * usx + c414) * usx + c412) * ux
    x41p =  ((c417 * usx + c415) * usx + c413) * usx + c411
    y41m = (((y400       + c416) * usy + c414) * usy + c412) * uy
    y41p =  ((c417 * usy + c415) * usy + c413) * usy + c411

    x42m = (((-x400 * 7.0_DP + c426) * usx + c424) * usx + c422) * ux
    x42p =  ((c427  * usx    + c425) * usx + c423) * usx + c421
    y42m = (((-y400 * 7.0_DP + c426) * usy + c424) * usy + c422) * uy
    y42p =  ((c427  * usy    + c425) * usy + c423) * usy + c421

    x43m = (((x400 * 21.0_DP + c436) * usx + c434) * usx + c432) * ux
    x43p =  ((c437 * usx     + c435) * usx + c433) * usx + c431
    y43m = (((y400 * 21.0_DP + c436) * usy + c434) * usy + c432) * uy
    y43p =  ((c437 * usy     + c435) * usy + c433) * usy + c431

    x44m = (((-x400 * 35.0_DP + c446) * usx + c444) * usx + c442) * ux
    x44p =  ((c447  * usx     + c445) * usx + c443) * usx + c441
    y44m = (((-y400 * 35.0_DP + c446) * usy + c444) * usy + c442) * uy
    y44p =  ((c447  * usy     + c445) * usy + c443) * usy + c441

    cy(1) =  y41p + y41m
    cy(2) =  y42p + y42m
    cy(3) =  y43p + y43m
    cy(4) =  y44p + y44m
    cy(5) =  y44p - y44m
    cy(6) =  y43p - y43m
    cy(7) =  y42p - y42m
    cy(8) =  y41p - y41m

    loop010 : DO l=1,l2d
      loop020 : DO my=1,m2
        ly    =  iy + my - 3
        w0(l) = ((x41p + x41m) * f2d(l,ix-2,ly)                    &
          &   +  (x42p + x42m) * f2d(l,ix-1,ly)                    &
          &   +  (x43p + x43m) * f2d(l,ix  ,ly)                    &
          &   +  (x44p + x44m) * f2d(l,ix+1,ly)                    &
          &   +  (x44p - x44m) * f2d(l,ix+2,ly)                    & 
          &   +  (x43p - x43m) * f2d(l,ix+3,ly)                    &
          &   +  (x42p - x42m) * f2d(l,ix+4,ly)                    &
          &   +  (x41p - x41m) * f2d(l,ix+5,ly)) * cy(my) + w0(l)
      END DO loop020
    END DO loop010


  END SUBROUTINE spl2df

  SUBROUTINE spl2dd (xd,        &!(in)
    &                w0, wx, wy &!(out)
    &               )
    
    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: xd(2)
    REAL(DP), INTENT(OUT) :: w0(l2d), &
      &                      wx(l2d), &
      &                      wy(l2d)
!Local variables
    INTEGER, PARAMETER :: m2 =  8
    INTEGER :: l,  &
      &        ix, &
      &        iy, &
      &        my, &
      &        ly
    REAL(DP) :: x,     &
      &         y,     &
      &         ux,    &
      &         uy,    &
      &         usx,   &
      &         usy,   &
      &         x400,  &
      &         y400,  &
      &         x41m,  &
      &         y41m,  &
      &         x41p,  &
      &         y41p,  &
      &         x42m,  &
      &         y42m,  &
      &         x42p,  &
      &         y42p,  &
      &         x43m,  &
      &         y43m,  &
      &         x43p,  &
      &         y43p,  &
      &         x44m,  &
      &         y44m,  &
      &         x44p,  &
      &         y44p,  &
      &         dx400, &
      &         dy400, &
      &         dx41m, &
      &         dy41m, &
      &         dx41p, &
      &         dy41p, &
      &         dx42m, &
      &         dy42m, &
      &         dx42p, &
      &         dy42p, &
      &         dx43m, &
      &         dy43m, &
      &         dx43p, &
      &         dy43p, &
      &         dx44m, &
      &         dy44m, &
      &         dx44p, &
      &         dy44p, &
      &         w0l,   &
      &         wxl,   &
      &         wyl,   &
      &         ww,    &
      &         cy(8), &
      &         dy(8)


    w0(:) =  0.0_DP
    wx(:) =  0.0_DP
    wy(:) =  0.0_DP


    x =  xd(1)
    IF(x < xs2) RETURN
    IF(x > xl2) RETURN
    
    y =  xd(2)
    IF(y < ys2) RETURN
    IF(y > yl2) RETURN

    ux = (x - xs2) / h2x
    ix =  ux
    IF(ix >= nx2d - 1)THEN
      ix =  nx2d - 2
      ux =  nx2d - 1
    END IF
 
    uy = (y - ys2) / h2y
    iy =  uy
    IF(iy >= ny2d - 1)THEN
      iy =  ny2d - 2
      uy =  ny2d - 1
    END IF

    ux  =  ux - ix - 0.5_DP
    uy  =  uy - iy - 0.5_DP

    usx =  ux * ux
    usy =  uy * uy

!############       case for m == 4       ###########################

    x400  = (((c41e * usx + c41c) * usx + c41a) * usx + c418) * usx
    y400  = (((c41e * usy + c41c) * usy + c41a) * usy + c418) * usy
    dx400 = (((d41e * usx + d41c) * usx + d41a) * usx + d418) * usx
    dy400 = (((d41e * usy + d41c) * usy + d41a) * usy + d418) * usy

    x41m  = (((x400       + c416) * usx + c414) * usx + c412) * ux
    y41m  = (((y400       + c416) * usy + c414) * usy + c412) * uy
    x41p  =  ((c417 * usx + c415) * usx + c413) * usx + c411
    y41p  =  ((c417 * usy + c415) * usy + c413) * usy + c411
    dx41m =  ((dx400      + d416) * usx + d414) * usx + c412
    dy41m =  ((dy400      + d416) * usy + d414) * usy + c412
    dx41p =  ((d417 * usx + d415) * usx + d413) * ux
    dy41p =  ((d417 * usy + d415) * usy + d413) * uy

    x42m  = (((-x400  * 7.0_DP + c426) * usx + c424) * usx + c422) * ux
    y42m  = (((-y400  * 7.0_DP + c426) * usy + c424) * usy + c422) * uy
    x42p  =  ((c427   * usx    + c425) * usx + c423) * usx + c421
    y42p  =  ((c427   * usy    + c425) * usy + c423) * usy + c421
    dx42m =  ((-dx400 * 7.0_DP + d426) * usx + d424) * usx + c422
    dy42m =  ((-dy400 * 7.0_DP + d426) * usy + d424) * usy + c422
    dx42p =  ((d427   * usx    + d425) * usx + d423) * ux
    dy42p =  ((d427   * usy    + d425) * usy + d423) * uy

    x43m  = (((x400  * 21.0_DP + c436) * usx + c434) * usx + c432) * ux
    y43m  = (((y400  * 21.0_DP + c436) * usy + c434) * usy + c432) * uy
    x43p  =  ((c437  * usx     + c435) * usx + c433) * usx + c431
    y43p  =  ((c437  * usy     + c435) * usy + c433) * usy + c431
    dx43m =  ((dx400 * 21.0_DP + d436) * usx + d434) * usx + c432
    dy43m =  ((dy400 * 21.0_DP + d436) * usy + d434) * usy + c432
    dx43p =  ((d437  * usx     + d435) * usx + d433) * ux
    dy43p =  ((d437  * usy     + d435) * usy + d433) * uy

    x44m  = (((-x400  * 35.0_DP + c446) * usx + c444) * usx + c442) * ux
    y44m  = (((-y400  * 35.0_DP + c446) * usy + c444) * usy + c442) * uy
    x44p  =  ((c447   * usx     + c445) * usx + c443) * usx + c441
    y44p  =  ((c447   * usy     + c445) * usy + c443) * usy + c441
    dx44m =  ((-dx400 * 35.0_DP + d446) * usx + d444) * usx + c442
    dy44m =  ((-dy400 * 35.0_DP + d446) * usy + d444) * usy + c442
    dx44p =  ((d447   * usx     + d445) * usx + d443) * ux
    dy44p =  (( d447  * usy     + d445) * usy + d443) * uy

    cy(1)   =  y41p +  y41m
    cy(2)   =  y42p +  y42m
    cy(3)   =  y43p +  y43m
    cy(4)   =  y44p +  y44m
    cy(5)   =  y44p -  y44m
    cy(6)   =  y43p -  y43m
    cy(7)   =  y42p -  y42m
    cy(8)   =  y41p -  y41m
    dy(1)   =  dy41p + dy41m
    dy(2)   =  dy42p + dy42m
    dy(3)   =  dy43p + dy43m
    dy(4)   =  dy44p + dy44m
    dy(5)   =  dy44p - dy44m
    dy(6)   =  dy43p - dy43m
    dy(7)   =  dy42p - dy42m
    dy(8)   =  dy41p - dy41m


    loop010 : DO l=1,l2d
      w0l =  0.0_DP
      wxl =  0.0_DP
      wyl =  0.0_DP

      loop020 : DO my=1,m2
        ly  =  iy + my - 3
        ww  = ((x41p + x41m) * f2d(l,ix-2,ly)   &
          & +  (x42p + x42m) * f2d(l,ix-1,ly)   &
          & +  (x43p + x43m) * f2d(l,ix  ,ly)   &
          & +  (x44p + x44m) * f2d(l,ix+1,ly)   &
          & +  (x44p - x44m) * f2d(l,ix+2,ly)   &
          & +  (x43p - x43m) * f2d(l,ix+3,ly)   &
          & +  (x42p - x42m) * f2d(l,ix+4,ly)   &
          & +  (x41p - x41m) * f2d(l,ix+5,ly))

        w0l =  ww * cy(my) + w0l
        wyl =  ww * dy(my) + wyl
        wxl = ((dx41p + dx41m) * f2d(l,ix-2,ly)                  &
          & +  (dx42p + dx42m) * f2d(l,ix-1,ly)                  &
          & +  (dx43p + dx43m) * f2d(l,ix  ,ly)                  &
          & +  (dx44p + dx44m) * f2d(l,ix+1,ly)                  &
          & +  (dx44p - dx44m) * f2d(l,ix+2,ly)                  &
          & +  (dx43p - dx43m) * f2d(l,ix+3,ly)                  &
          & +  (dx42p - dx42m) * f2d(l,ix+4,ly)                  &
          & +  (dx41p - dx41m) * f2d(l,ix+5,ly)) * cy(my) + wxl

      END DO loop020
         
      w0(l) =  w0l
      wx(l) =  wxl / h2x
      wy(l) =  wyl / h2y

    END DO loop010


    RETURN
  END SUBROUTINE spl2dd

  SUBROUTINE spl3df (xd, & !(in) 
    &                w0  & !(out)
    &               )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN)  :: xd(3)
    REAL(DP), INTENT(OUT) :: w0(l3d)
!Local variables
    INTEGER, PARAMETER :: m2 = 8
    INTEGER :: l,  &
      &        mz, &
      &        lz, &
      &        my, &
      &        ly, &
      &        ix, &
      &        iy, &
      &        iz
    REAL(DP) :: x,     &
      &         y,     &
      &         z,     &
      &         ux,    &
      &         uy,    &
      &         uz,    &
      &         usx,   &
      &         usy,   &
      &         usz,   &
      &         x400,  &
      &         y400,  &
      &         z400,  &
      &         x41m,  &
      &         y41m,  &
      &         z41m,  &
      &         x41p,  &
      &         y41p,  &
      &         z41p,  &
      &         x42m,  &
      &         y42m,  &
      &         z42m,  &
      &         x42p,  &
      &         y42p,  &
      &         z42p,  &
      &         x43m,  &
      &         y43m,  &
      &         z43m,  &
      &         x43p,  &
      &         y43p,  &
      &         z43p,  &
      &         x44m,  &
      &         y44m,  &
      &         z44m,  &
      &         x44p,  &
      &         y44p,  &
      &         z44p,  &
      &         cyz,   &
      &         cy(8), &
      &         cz(8)


    w0(:) =  0.0_DP


    x =  xd(1)
    IF(x < xs3) RETURN
    IF(x > xl3) RETURN

    y =  xd(2)   
    IF(y < ys3) RETURN
    IF(y > yl3) RETURN

    z =  xd(3)
    IF(z < zs3) z =  zs3
    IF(z > zl3) z =  zl3

    ux = (x - xs3) / h3x
    ix =  ux
    IF(ix >= nx3d - 1)THEN
      ix =  nx3d - 2
      ux =  nx3d - 1
    END IF
 
    uy = (y - ys3) / h3y
    iy =  uy
    IF(iy >= ny3d - 1)THEN
      iy =  ny3d - 2
      uy =  ny3d - 1
    END IF

    uz = (z - zs3) / h3z
    iz =  uz
    IF(iz >= nz3d - 1)THEN
      iz =  nz3d - 2
      uz =  nz3d - 1
    END IF

    ux  =  ux - ix - 0.5_DP
    uy  =  uy - iy - 0.5_DP
    uz  =  uz - iz - 0.5_DP
    usx =  ux * ux
    usy =  uy * uy
    usz =  uz * uz

!###########       case for m == 4       #############################

    x400 = (((c41e * usx + c41c) * usx + c41a) * usx + c418) * usx
    y400 = (((c41e * usy + c41c) * usy + c41a) * usy + c418) * usy
    z400 = (((c41e * usz + c41c) * usz + c41a) * usz + c418) * usz
    
    x41m = (((x400       + c416) * usx + c414) * usx + c412) * ux
    x41p =  ((c417 * usx + c415) * usx + c413) * usx + c411
    y41m = (((y400       + c416) * usy + c414) * usy + c412) * uy
    y41p =  ((c417 * usy + c415) * usy + c413) * usy + c411
    z41m = (((z400       + c416) * usz + c414) * usz + c412) * uz
    z41p =  ((c417 * usz + c415) * usz + c413) * usz + c411

    x42m = (((-x400 * 7.0_DP + c426) * usx + c424) * usx + c422) * ux
    x42p =  ((c427  * usx    + c425) * usx + c423) * usx + c421
    y42m = (((-y400 * 7.0_DP + c426) * usy + c424) * usy + c422) * uy
    y42p =  ((c427  * usy    + c425) * usy + c423) * usy + c421
    z42m = (((-z400 * 7.0_DP + c426) * usz + c424) * usz + c422) * uz
    z42p =  ((c427  * usz    + c425) * usz + c423) * usz + c421

    x43m = (((x400 * 21.0_DP + c436) * usx + c434) * usx + c432) * ux
    x43p =  ((c437 * usx     + c435) * usx + c433) * usx + c431
    y43m = (((y400 * 21.0_DP + c436) * usy + c434) * usy + c432) * uy
    y43p =  ((c437 * usy     + c435) * usy + c433) * usy + c431
    z43m = (((z400 * 21.0_DP + c436) * usz + c434) * usz + c432) * uz
    z43p =  ((c437 * usz     + c435) * usz + c433) * usz + c431

    x44m = (((-x400 * 35.0_DP + c446) * usx + c444) * usx + c442) * ux
    x44p =  ((c447  * usx     + c445) * usx + c443) * usx + c441
    y44m = (((-y400 * 35.0_DP + c446) * usy + c444) * usy + c442) * uy
    y44p =  ((c447  * usy     + c445) * usy + c443) * usy + c441
    z44m = (((-z400 * 35.0_DP + c446) * usz + c444) * usz + c442) * uz
    z44p =  ((c447  * usz     + c445) * usz + c443) * usz + c441

    cy(1) =  y41p + y41m
    cy(2) =  y42p + y42m
    cy(3) =  y43p + y43m
    cy(4) =  y44p + y44m
    cy(5) =  y44p - y44m
    cy(6) =  y43p - y43m
    cy(7) =  y42p - y42m
    cy(8) =  y41p - y41m
    cz(1) =  z41p + z41m
    cz(2) =  z42p + z42m
    cz(3) =  z43p + z43m
    cz(4) =  z44p + z44m
    cz(5) =  z44p - z44m
    cz(6) =  z43p - z43m
    cz(7) =  z42p - z42m
    cz(8) =  z41p - z41m

    loop010 : DO l=1,l3d
      loop020 : DO mz=1,m2
        lz =  iz + mz - 3
        loop030 : DO my=1,m2
          ly    =  iy + my - 3
          cyz   =  cy(my) * cz(mz)
          w0(l) = ((x41p + x41m) * f3d(l,ix-2,ly,lz)                &
            &   +  (x42p + x42m) * f3d(l,ix-1,ly,lz)                &
            &   +  (x43p + x43m) * f3d(l,ix  ,ly,lz)                &
            &   +  (x44p + x44m) * f3d(l,ix+1,ly,lz)                &
            &   +  (x44p - x44m) * f3d(l,ix+2,ly,lz)                & 
            &   +  (x43p - x43m) * f3d(l,ix+3,ly,lz)                &
            &   +  (x42p - x42m) * f3d(l,ix+4,ly,lz)                &
            &   +  (x41p - x41m) * f3d(l,ix+5,ly,lz)) * cyz + w0(l)
        END DO loop030
      END DO loop020
    END DO loop010


    RETURN
  END SUBROUTINE spl3df

  SUBROUTINE spl3dd (xd,            &!(in)
    &                w0, wx, wy, wz &!(out)
    &               )
    
    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN) :: xd(3)
    REAL(DP), INTENT(OUT) :: w0(l3d), &
      &                      wx(l3d), &
      &                      wy(l3d), &
      &                      wz(l3d)
!Local variables
    INTEGER, PARAMETER :: m2 = 8
    INTEGER :: l,   &
      &        ix,  &
      &        iy,  &
      &        iz,  &
      &        mz,  &
      &        lz,  &
      &        ix1, &
      &        ix2, &
      &        ix4, &
      &        ix5, &
      &        ix6, &
      &        ix7, &
      &        ix8, &
      &        iy1, &
      &        iy2, &
      &        iy4, &
      &        iy5, &
      &        iy6, &
      &        iy7, &
      &        iy8
    REAL(DP) :: x,     &
      &         y,     &
      &         z,     &
      &         ux,    &
      &         uy,    &
      &         uz,    &
      &         usx,   &
      &         usy,   &
      &         usz,   &
      &         x400,  &
      &         y400,  &
      &         z400,  &
      &         x41m,  &
      &         y41m,  &
      &         z41m,  &
      &         x41p,  &
      &         y41p,  &
      &         z41p,  &
      &         x42m,  &
      &         y42m,  &
      &         z42m,  &
      &         x42p,  &
      &         y42p,  &
      &         z42p,  &
      &         x43m,  &
      &         y43m,  &
      &         z43m,  &
      &         x43p,  &
      &         y43p,  &
      &         z43p,  &
      &         x44m,  &
      &         y44m,  &
      &         z44m,  &
      &         x44p,  &
      &         y44p,  &
      &         z44p,  &
      &         dx400, &
      &         dy400, &
      &         dz400, &
      &         dx41m, &
      &         dy41m, &
      &         dz41m, &
      &         dx41p, &
      &         dy41p, &
      &         dz41p, &
      &         dx42m, &
      &         dy42m, &
      &         dz42m, &
      &         dx42p, &
      &         dy42p, &
      &         dz42p, &
      &         dx43m, &
      &         dy43m, &
      &         dz43m, &
      &         dx43p, &
      &         dy43p, &
      &         dz43p, &
      &         dx44m, &
      &         dy44m, &
      &         dz44m, &
      &         dx44p, &
      &         dy44p, &
      &         dz44p, &
      &         cx1,   &
      &         cx2,   &
      &         cx3,   &
      &         cx4,   &
      &         cx5,   &
      &         cx6,   &
      &         cx7,   &
      &         cx8,   &
      &         cy1,   &
      &         cy2,   &
      &         cy3,   &
      &         cy4,   &
      &         cy5,   &
      &         cy6,   &
      &         cy7,   &
      &         cy8,   &
      &         dx1,   &
      &         dx2,   &
      &         dx3,   &
      &         dx4,   &
      &         dx5,   &
      &         dx6,   &
      &         dx7,   &
      &         dx8,   &
      &         dy1,   &
      &         dy2,   &
      &         dy3,   &
      &         dy4,   &
      &         dy5,   &
      &         dy6,   &
      &         dy7,   &
      &         dy8,   &
      &         w0l,   &
      &         wxl,   &
      &         wyl,   &
      &         wzl,   &
      &         w01,   &
      &         w02,   &
      &         w03,   &
      &         w04,   &
      &         w05,   &
      &         w06,   &
      &         w07,   &
      &         w08,   &
      &         w00,   &
      &         cz(8), &
      &         dz(8)


    w0(:) =  0.0_DP
    wx(:) =  0.0_DP
    wy(:) =  0.0_DP
    wz(:) =  0.0_DP

    x =  xd(1)
    IF(x < xs3) RETURN
    IF(x > xl3) RETURN

    y =  xd(2)
    IF(y < ys3) RETURN
    IF(y > yl3) RETURN

    z =  xd(3)
    IF(z < zs3) z =  zs3
    IF(z > zl3) z =  zl3

    ux = (x - xs3) / h3x
    ix =  ux
    IF(ix >= nx3d - 1)THEN
      ix =  nx3d - 2
      ux =  nx3d - 1
    END IF
 
    uy = (y - ys3) / h3y
    iy =  uy
    IF(iy >= ny3d - 1)THEN
      iy =  ny3d - 2
      uy =  ny3d - 1
    END IF

    uz = (z - zs3) / h3z
    iz =  uz
    IF(iz >= nz3d - 1)THEN
      iz =  nz3d - 2
      uz =  nz3d - 1
    END IF

    ux  =  ux - ix - 0.5_DP
    uy  =  uy - iy - 0.5_DP
    uz  =  uz - iz - 0.5_DP

    usx =  ux * ux
    usy =  uy * uy
    usz =  uz * uz

    ix1 =  ix - 2
    ix2 =  ix - 1
    ix4 =  ix + 1
    ix5 =  ix + 2
    ix6 =  ix + 3
    ix7 =  ix + 4
    ix8 =  ix + 5

    iy1 =  iy - 2
    iy2 =  iy - 1
    iy4 =  iy + 1
    iy5 =  iy + 2
    iy6 =  iy + 3
    iy7 =  iy + 4
    iy8 =  iy + 5

!############       case for m == 4       ###########################

    x400  = (((c41e * usx + c41c) * usx + c41a) * usx + c418) * usx
    y400  = (((c41e * usy + c41c) * usy + c41a) * usy + c418) * usy
    z400  = (((c41e * usz + c41c) * usz + c41a) * usz + c418) * usz
    dx400 = (((d41e * usx + d41c) * usx + d41a) * usx + d418) * usx
    dy400 = (((d41e * usy + d41c) * usy + d41a) * usy + d418) * usy
    dz400 = (((d41e * usz + d41c) * usz + d41a) * usz + d418) * usz

    x41m  = (((x400       + c416) * usx + c414) * usx + c412) * ux
    y41m  = (((y400       + c416) * usy + c414) * usy + c412) * uy
    z41m  = (((z400       + c416) * usz + c414) * usz + c412) * uz
    x41p  =  ((c417 * usx + c415) * usx + c413) * usx + c411
    y41p  =  ((c417 * usy + c415) * usy + c413) * usy + c411
    z41p  =  ((c417 * usz + c415) * usz + c413) * usz + c411
    dx41m =  ((dx400      + d416) * usx + d414) * usx + c412
    dy41m =  ((dy400      + d416) * usy + d414) * usy + c412
    dz41m =  ((dz400      + d416) * usz + d414) * usz + c412
    dx41p =  ((d417 * usx + d415) * usx + d413) * ux
    dy41p =  ((d417 * usy + d415) * usy + d413) * uy
    dz41p =  ((d417 * usz + d415) * usz + d413) * uz

    x42m  = (((-x400  * 7.0_DP + c426) * usx + c424) * usx + c422) * ux
    y42m  = (((-y400  * 7.0_DP + c426) * usy + c424) * usy + c422) * uy
    z42m  = (((-z400  * 7.0_DP + c426) * usz + c424) * usz + c422) * uz
    x42p  =  ((c427   * usx    + c425) * usx + c423) * usx + c421
    y42p  =  ((c427   * usy    + c425) * usy + c423) * usy + c421
    z42p  =  ((c427   * usz    + c425) * usz + c423) * usz + c421
    dx42m =  ((-dx400 * 7.0_DP + d426) * usx + d424) * usx + c422
    dy42m =  ((-dy400 * 7.0_DP + d426) * usy + d424) * usy + c422
    dz42m =  ((-dz400 * 7.0_DP + d426) * usz + d424) * usz + c422
    dx42p =  ((d427   * usx    + d425) * usx + d423) * ux
    dy42p =  ((d427   * usy    + d425) * usy + d423) * uy
    dz42p =  ((d427   * usz    + d425) * usz + d423) * uz

    x43m  = (((x400  * 21.0_DP + c436) * usx + c434) * usx + c432) * ux
    y43m  = (((y400  * 21.0_DP + c436) * usy + c434) * usy + c432) * uy
    z43m  = (((z400  * 21.0_DP + c436) * usz + c434) * usz + c432) * uz
    x43p  =  ((c437  * usx     + c435) * usx + c433) * usx + c431
    y43p  =  ((c437  * usy     + c435) * usy + c433) * usy + c431
    z43p  =  ((c437  * usz     + c435) * usz + c433) * usz + c431
    dx43m =  ((dx400 * 21.0_DP + d436) * usx + d434) * usx + c432
    dy43m =  ((dy400 * 21.0_DP + d436) * usy + d434) * usy + c432
    dz43m =  ((dz400 * 21.0_DP + d436) * usz + d434) * usz + c432
    dx43p =  ((d437  * usx     + d435) * usx + d433) * ux
    dy43p =  ((d437  * usy     + d435) * usy + d433) * uy
    dz43p =  ((d437  * usz     + d435) * usz + d433) * uz

    x44m  = (((-x400  * 35.0_DP + c446) * usx + c444) * usx + c442) * ux
    y44m  = (((-y400  * 35.0_DP + c446) * usy + c444) * usy + c442) * uy
    z44m  = (((-z400  * 35.0_DP + c446) * usz + c444) * usz + c442) * uz
    x44p  =  ((c447   * usx     + c445) * usx + c443) * usx + c441
    y44p  =  ((c447   * usy     + c445) * usy + c443) * usy + c441
    z44p  =  ((c447   * usz     + c445) * usz + c443) * usz + c441
    dx44m =  ((-dx400 * 35.0_DP + d446) * usx + d444) * usx + c442
    dy44m =  ((-dy400 * 35.0_DP + d446) * usy + d444) * usy + c442
    dz44m =  ((-dz400 * 35.0_DP + d446) * usz + d444) * usz + c442
    dx44p =  ((d447   * usx     + d445) * usx + d443) * ux
    dy44p =  ((d447   * usy     + d445) * usy + d443) * uy
    dz44p =  ((d447   * usz     + d445) * usz + d443) * uz

    cx1   =  x41p +  x41m
    cx2   =  x42p +  x42m
    cx3   =  x43p +  x43m
    cx4   =  x44p +  x44m
    cx5   =  x44p -  x44m
    cx6   =  x43p -  x43m
    cx7   =  x42p -  x42m
    cx8   =  x41p -  x41m
    cy1   =  y41p +  y41m
    cy2   =  y42p +  y42m
    cy3   =  y43p +  y43m
    cy4   =  y44p +  y44m
    cy5   =  y44p -  y44m
    cy6   =  y43p -  y43m
    cy7   =  y42p -  y42m
    cy8   =  y41p -  y41m
    cz(1) =  z41p +  z41m
    cz(2) =  z42p +  z42m
    cz(3) =  z43p +  z43m
    cz(4) =  z44p +  z44m
    cz(5) =  z44p -  z44m
    cz(6) =  z43p -  z43m
    cz(7) =  z42p -  z42m
    cz(8) =  z41p -  z41m
    dx1   =  dx41p + dx41m
    dx2   =  dx42p + dx42m
    dx3   =  dx43p + dx43m
    dx4   =  dx44p + dx44m
    dx5   =  dx44p - dx44m
    dx6   =  dx43p - dx43m
    dx7   =  dx42p - dx42m
    dx8   =  dx41p - dx41m
    dy1   =  dy41p + dy41m
    dy2   =  dy42p + dy42m
    dy3   =  dy43p + dy43m
    dy4   =  dy44p + dy44m
    dy5   =  dy44p - dy44m
    dy6   =  dy43p - dy43m
    dy7   =  dy42p - dy42m
    dy8   =  dy41p - dy41m
    dz(1) =  dz41p + dz41m
    dz(2) =  dz42p + dz42m
    dz(3) =  dz43p + dz43m
    dz(4) =  dz44p + dz44m
    dz(5) =  dz44p - dz44m
    dz(6) =  dz43p - dz43m
    dz(7) =  dz42p - dz42m
    dz(8) =  dz41p - dz41m


    loop010 : DO l=1,l3d
      w0l =  0.0_DP
      wxl =  0.0_DP
      wyl =  0.0_DP
      wzl =  0.0_DP

      loop020 : DO mz=1,m2
        lz =  iz + mz - 3
        w01 =     cx1 * f3d(l,ix1,iy1,lz) + cx2 * f3d(l,ix2,iy1,lz)       &
          & +     cx3 * f3d(l,ix, iy1,lz) + cx4 * f3d(l,ix4,iy1,lz)       &
          & +     cx5 * f3d(l,ix5,iy1,lz) + cx6 * f3d(l,ix6,iy1,lz)       &
          & +     cx7 * f3d(l,ix7,iy1,lz) + cx8 * f3d(l,ix8,iy1,lz)
        w02 =     cx1 * f3d(l,ix1,iy2,lz) + cx2 * f3d(l,ix2,iy2,lz)       &
          & +     cx3 * f3d(l,ix, iy2,lz) + cx4 * f3d(l,ix4,iy2,lz)       &
          & +     cx5 * f3d(l,ix5,iy2,lz) + cx6 * f3d(l,ix6,iy2,lz)       &
          & +     cx7 * f3d(l,ix7,iy2,lz) + cx8 * f3d(l,ix8,iy2,lz)
        w03 =     cx1 * f3d(l,ix1,iy, lz) + cx2 * f3d(l,ix2,iy, lz)       &
          & +     cx3 * f3d(l,ix, iy, lz) + cx4 * f3d(l,ix4,iy, lz)       &
          & +     cx5 * f3d(l,ix5,iy, lz) + cx6 * f3d(l,ix6,iy, lz)       &
          & +     cx7 * f3d(l,ix7,iy, lz) + cx8 * f3d(l,ix8,iy, lz)
        w04 =     cx1 * f3d(l,ix1,iy4,lz) + cx2 * f3d(l,ix2,iy4,lz)       &
          & +     cx3 * f3d(l,ix, iy4,lz) + cx4 * f3d(l,ix4,iy4,lz)       &
          & +     cx5 * f3d(l,ix5,iy4,lz) + cx6 * f3d(l,ix6,iy4,lz)       &
          & +     cx7 * f3d(l,ix7,iy4,lz) + cx8 * f3d(l,ix8,iy4,lz)
        w05 =     cx1 * f3d(l,ix1,iy5,lz) + cx2 * f3d(l,ix2,iy5,lz)       &
          & +     cx3 * f3d(l,ix, iy5,lz) + cx4 * f3d(l,ix4,iy5,lz)       &
          & +     cx5 * f3d(l,ix5,iy5,lz) + cx6 * f3d(l,ix6,iy5,lz)       &
          & +     cx7 * f3d(l,ix7,iy5,lz) + cx8 * f3d(l,ix8,iy5,lz)
        w06 =     cx1 * f3d(l,ix1,iy6,lz) + cx2 * f3d(l,ix2,iy6,lz)       &
          & +     cx3 * f3d(l,ix, iy6,lz) + cx4 * f3d(l,ix4,iy6,lz)       &
          & +     cx5 * f3d(l,ix5,iy6,lz) + cx6 * f3d(l,ix6,iy6,lz)       &
          & +     cx7 * f3d(l,ix7,iy6,lz) + cx8 * f3d(l,ix8,iy6,lz)
        w07 =     cx1 * f3d(l,ix1,iy7,lz) + cx2 * f3d(l,ix2,iy7,lz)       &
          & +     cx3 * f3d(l,ix, iy7,lz) + cx4 * f3d(l,ix4,iy7,lz)       &
          & +     cx5 * f3d(l,ix5,iy7,lz) + cx6 * f3d(l,ix6,iy7,lz)       &
          & +     cx7 * f3d(l,ix7,iy7,lz) + cx8 * f3d(l,ix8,iy7,lz)
        w08 =     cx1 * f3d(l,ix1,iy8,lz) + cx2 * f3d(l,ix2,iy8,lz)       &
          & +     cx3 * f3d(l,ix, iy8,lz) + cx4 * f3d(l,ix4,iy8,lz)       &
          & +     cx5 * f3d(l,ix5,iy8,lz) + cx6 * f3d(l,ix6,iy8,lz)       &
          & +     cx7 * f3d(l,ix7,iy8,lz) + cx8 * f3d(l,ix8,iy8,lz)

        wxl = ((dx1 * f3d(l,ix1,iy1,lz) + dx2 * f3d(l,ix2,iy1,lz)         &
          & +   dx3 * f3d(l,ix, iy1,lz) + dx4 * f3d(l,ix4,iy1,lz)         &
          & +   dx5 * f3d(l,ix5,iy1,lz) + dx6 * f3d(l,ix6,iy1,lz)         &
          & +   dx7 * f3d(l,ix7,iy1,lz) + dx8 * f3d(l,ix8,iy1,lz)) * cy1  &
          & +  (dx1 * f3d(l,ix1,iy2,lz) + dx2 * f3d(l,ix2,iy2,lz)         &
          & +   dx3 * f3d(l,ix, iy2,lz) + dx4 * f3d(l,ix4,iy2,lz)         &
          & +   dx5 * f3d(l,ix5,iy2,lz) + dx6 * f3d(l,ix6,iy2,lz)         &
          & +   dx7 * f3d(l,ix7,iy2,lz) + dx8 * f3d(l,ix8,iy2,lz)) * cy2  &
          & +  (dx1 * f3d(l,ix1,iy, lz) + dx2 * f3d(l,ix2,iy, lz)         &
          & +   dx3 * f3d(l,ix, iy, lz) + dx4 * f3d(l,ix4,iy, lz)         &
          & +   dx5 * f3d(l,ix5,iy, lz) + dx6 * f3d(l,ix6,iy, lz)         &
          & +   dx7 * f3d(l,ix7,iy, lz) + dx8 * f3d(l,ix8,iy, lz)) * cy3  &
          & +  (dx1 * f3d(l,ix1,iy4,lz) + dx2 * f3d(l,ix2,iy4,lz)         &
          & +   dx3 * f3d(l,ix, iy4,lz) + dx4 * f3d(l,ix4,iy4,lz)         &
          & +   dx5 * f3d(l,ix5,iy4,lz) + dx6 * f3d(l,ix6,iy4,lz)         &
          & +   dx7 * f3d(l,ix7,iy4,lz) + dx8 * f3d(l,ix8,iy4,lz)) * cy4  &
          & +  (dx1 * f3d(l,ix1,iy5,lz) + dx2 * f3d(l,ix2,iy5,lz)         &
          & +   dx3 * f3d(l,ix, iy5,lz) + dx4 * f3d(l,ix4,iy5,lz)         &
          & +   dx5 * f3d(l,ix5,iy5,lz) + dx6 * f3d(l,ix6,iy5,lz)         &
          & +   dx7 * f3d(l,ix7,iy5,lz) + dx8 * f3d(l,ix8,iy5,lz)) * cy5  &
          & +  (dx1 * f3d(l,ix1,iy6,lz) + dx2 * f3d(l,ix2,iy6,lz)         &
          & +   dx3 * f3d(l,ix, iy6,lz) + dx4 * f3d(l,ix4,iy6,lz)         &
          & +   dx5 * f3d(l,ix5,iy6,lz) + dx6 * f3d(l,ix6,iy6,lz)         &
          & +   dx7 * f3d(l,ix7,iy6,lz) + dx8 * f3d(l,ix8,iy6,lz)) * cy6  &
          & +  (dx1 * f3d(l,ix1,iy7,lz) + dx2 * f3d(l,ix2,iy7,lz)         &
          & +   dx3 * f3d(l,ix, iy7,lz) + dx4 * f3d(l,ix4,iy7,lz)         &
          & +   dx5 * f3d(l,ix5,iy7,lz) + dx6 * f3d(l,ix6,iy7,lz)         &
          & +   dx7 * f3d(l,ix7,iy7,lz) + dx8 * f3d(l,ix8,iy7,lz)) * cy7  &
          & +  (dx1 * f3d(l,ix1,iy8,lz) + dx2 * f3d(l,ix2,iy8,lz)         &
          & +   dx3 * f3d(l,ix, iy8,lz) + dx4 * f3d(l,ix4,iy8,lz)         &
          & +   dx5 * f3d(l,ix5,iy8,lz) + dx6 * f3d(l,ix6,iy8,lz)         &
          & +   dx7 * f3d(l,ix7,iy8,lz) + dx8 * f3d(l,ix8,iy8,lz)) * cy8) &
          & *   cz(mz) + wxl

        w00 =  w01 * cy1 + w02 * cy2 + w03 * cy3 + w04 * cy4              &
          & +  w05 * cy5 + w06 * cy6 + w07 * cy7 + w08 * cy8
        w0l =  w00 * cz(mz) + w0l
        wzl =  w00 * dz(mz) + wzl
        wyl = (w01 * dy1 + w02 * dy2 + w03 * dy3 + w04 * dy4              &
          & +  w05 * dy5 + w06 * dy6 + w07 * dy7 + w08 * dy8) * cz(mz)    &
          & +  wyl

      END DO loop020
         
      w0(l) =  w0l
      wx(l) =  wxl / h3x
      wy(l) =  wyl / h3y
      wz(l) =  wzl / h3z

    END DO loop010


    RETURN
  END SUBROUTINE spl3dd

END MODULE spline_mod
!===============================================================================



!===============================================================================
!=cylindrical_coord_mod.f90
!
!==Version
!
! $Revision: $
! $Id: $
!
!==Overview
!
!==Reference
!
!==Error Handlings
!
!==Known Bugs
!
!==Note
!
!==TODO
!
MODULE cylindrical_coord_mod

  USE kind_spec
  USE param1,     ONLY : pi2
  USE spline_mod, ONLY : l3d,    &
    &                    f3d,    &
    &                    nx3d,   &
    &                    ny3d,   &
    &                    nz3d,   &
    &                    splin3, &
    &                    spl3df, &
    &                    spl3dd

  IMPLICIT NONE

  PRIVATE

! File format of the magnetic field
!

  CHARACTER(LEN=20) :: mag_form = '', &
    &                  version  = ''
!----------------------------------------------------------------------------
! A parameter MAG_FORM prescribes the table format of the magnetic field.
! A parameter VERSION prescribes the version for the old legacy format.
!
! 1. mgrid:   file format used in the MAKEGRID code
! 2. mgo:     file format used in KMAG/KMAG2 codes
! 3. mag:     file format for the equilibrium field used in HINT/HINT2 codes
! 3.1. ver2:  legacy HINT2 format
! 4. vac:     file format for the vacuum field used in HINT/HINT2 codes
! 4.1. ver2:  legacy HINT2 format
! 5. movie:   file format used in the MIPS code
! 6. mips_eq: file format for the equilibrium field used in the MIPS code
!----------------------------------------------------------------------------
  LOGICAL :: lflux
!
! A parameter LFLUX switch on/off to read the flux distribution.
!
  CHARACTER(LEN=20) :: flx_form = 'xss'
!
! A parameter FLX_FORM prescribes the table format of the flux distribution
!
! 1. xss:  file format used in the SMAP code
! 2. eqdsk: file format used in the EFIT code
! 3. eqdata: file format used in the TOPICS
!
  LOGICAL :: lsymmetry
!
! 
!
  INTEGER :: mtor,                & ! toroidal field period
    &        nr0b,                & ! grid number along R-direction
    &        nt0b,                & ! grid number along phi-direction
    &        nz0b,                & ! grid number along Z-direction
!
    &        kstep   =  99999,    & ! time steps of MIPS
    &        igrid(4)               ! grid number of MIPS
  REAL(DP) :: bmax,                  &
    &         pi2m,                  & ! one toroidal field period [rad]
    &         rmaxb,                 & ! major radius of the inner side of computational DOmain [m]
    &         rminb,                 & ! major radius of the outer side of computational DOmain [m]
    &         zmaxb,                 & ! height of the lower side of computational DOmain [m]
    &         zminb,                 & ! height of the higher side of the computational DOmain [m]
    &         delrb,                 & ! grid size along R-direction [m]
    &         delzb,                 & ! grid size along Z-direction [m]
    &         badjust   =  1.0_DP,   & ! normalization factor of magnetic field
    &         bnorm     =  3.0_DP,   & ! normalization factor of magnetic field
    &         cj(50)    =  0.0_DP,   & ! work array to store coil current in each coils [A]
!                                        ! CAUTION! array size is 50 (max)
    &         cturn(50) =  1.0_DP,   & ! work array to store coil current in each coils [A]
!                                        ! CAUTION! array size is 50 (max)
    &         cfact(50) =  1.0_DP,   & ! work array to store coil current in each coils [A]
!                                        ! CAUTION! array size is 50 (max)
    &         mbound(4) =  0.0_DP      ! boundary of MIPS
!                                      ! mbound(1) = Rmin [m]
!                                      ! mbound(2) = Rmax [m]
!                                      ! mbound(3) = Zmin [m]
!                                      ! mbound(4) = Zmax [m]
  REAL(DP), ALLOCATABLE :: rg(:), &
    &                      zg(:)


  NAMELIST /nlinp_coil_dat/ mag_form,  &
    &                       flx_form,  &
    &                       version,   &
    &                       lflux,     &
    &                       lsymmetry, &
    &                       cj,        &
    &                       cturn,     &
    &                       cfact,     &
    &                       badjust,   &
    &                       bnorm,     &
    &                       igrid,     &
    &                       mbound,    &
    &                       kstep

  PUBLIC :: mag_form,       &
    &       flx_form,       &
    &       version,        &
    &       lflux,          &
    &       lsymmetry,      &
    &       cj,             &
    &       cturn,          &
    &       cfact,          &
    &       badjust,        &
    &       bnorm,          &
    &       igrid,          &
    &       mbound,         &
    &       kstep,          &
    &       pi2m,           &
    &       mtor,           &
    &       nr0b,           &
    &       nt0b,           &
    &       nz0b,           &
    &       rminb,          &
    &       rmaxb,          &
    &       zminb,          &
    &       zmaxb,          &
    &       delrb,          &
    &       delzb,          &
    &       rg,             &
    &       zg,             &
    &       free_mem_field, &
    &       magset,         &
    &       mgval1,         &
    &       mgval2,         &
    &       mgval3,         &
    &       nlinp_coil_dat

CONTAINS

  SUBROUTINE read_field

    IMPLICIT NONE


    l3d =  4
    IF(lflux)THEN
      l3d =  5
    END IF

    SELECT CASE(TRIM(mag_form))
      CASE('mgrid', "MGRID")
        CALL read_mgrid
      CASE('mgo', "MGO")
        CALL read_mgo
      CASE('mag', "MAG")
        CALL read_hint_eq
      CASE('vac', "VAC")
        CALL read_hint_vac
      CASE('mips', "MIPS")
        CALL read_mips
      CASE('mips_eq', "MIPS_EQ")
        CALL read_mips_eq
      CASE DEFAULT
        CALL read_mgrid
    END SELECT

    IF(lflux)THEN
      SELECT CASE(TRIM(flx_form))
        CASE('xss')
          CALL read_flux
        CASE('eqdsk')
          CALL read_eqdsk
        CASE('eqdata')
          CALL read_eqdata
        CASE DEFAULT
          CALL read_flux
      END SELECT
    END IF


  END SUBROUTINE read_field

  SUBROUTINE read_mgrid

    IMPLICIT NONE

    INTEGER :: nextcur, &
      &        i,       &
      &        j,       &
      &        k,       &
      &        n
    REAL(DP), ALLOCATABLE :: br(:,:,:), &
      &                      bp(:,:,:), &
      &                      bz(:,:,:)
    CHARACTER(LEN=30), ALLOCATABLE :: curlabel(:)
    CHARACTER(LEN=100) :: fmt
    LOGICAL :: lstyle2000


    READ(25) nr0b, nz0b, nt0b, mtor, nextcur
    READ(25) rminb, zminb, rmaxb, zmaxb

    IF(nextcur < 0) lstyle2000 = .true.
    nextcur = ABS(nextcur)

    IF(nextcur > 50) STOP ' nextcur > 50'

    ALLOCATE(curlabel(nextcur))
    ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))

    f3d(:,:,:,:) =  0.0_DP

    READ(25) (curlabel(i), i=1,nextcur)

    ALLOCATE(br(nr0b,nz0b,nt0b), bp(nr0b,nz0b,nt0b), bz(nr0b,nz0b,nt0b))

    PRINT *, ' Reading coil data'
    PRINT *

    PRINT *
    PRINT *, ' coil No    c_I [A/T]      trun  c_factor '
    PRINT *, '--------------------------------------------'
 
    fmt = '(I7, ES17.6, I7, ES12.4)'

    DO n=1,nextcur
      PRINT fmt, n, cj(n), INT(cturn(n)), cfact(n)
    END DO

    PRINT *

    fmt = '(A10, I3, A15, A5, F12.6, A6)'

    loop010 : DO n=1,nextcur 
      PRINT fmt, ' coil NO: ', n, trim(curlabel(n)), ' c_I ', cj(n) * cfact(n) * cturn(n) * badjust * 1.0E-06_DP, ' [MA]'
      IF(lstyle2000)THEN
        READ(25) br, bp, bz
      ELSE
        READ(25) (((br(i,j,k), bz(i,j,k), bp(i,j,k), i=1,nr0b), j=1,nz0b), k=1,nt0b)
      END IF
      f3d(1,1:nr0b,1:nz0b,1:nt0b) =  f3d(1,1:nr0b,1:nz0b,1:nt0b) + cj(n) * cfact(n) * cturn(n) * badjust * br(1:nr0b,1:nz0b,1:nt0b) !< B_R
      f3d(2,1:nr0b,1:nz0b,1:nt0b) =  f3d(2,1:nr0b,1:nz0b,1:nt0b) + cj(n) * cfact(n) * cturn(n) * badjust * bp(1:nr0b,1:nz0b,1:nt0b) !< B_phi
      f3d(3,1:nr0b,1:nz0b,1:nt0b) =  f3d(3,1:nr0b,1:nz0b,1:nt0b) + cj(n) * cfact(n) * cturn(n) * badjust * bz(1:nr0b,1:nz0b,1:nt0b) !< B_Z
    END DO loop010


    DEALLOCATE(br, bp, bz)
    DEALLOCATE(curlabel)


  END SUBROUTINE read_mgrid

  SUBROUTINE read_mgo

    IMPLICIT NONE

    INTEGER :: i, &
      &        j, &
      &        k
    REAL(DP), ALLOCATABLE :: br(:,:,:), &
      &                      bp(:,:,:), &
      &                      bz(:,:,:)
    CHARACTER(LEN=100) :: fmt


    READ(25) nr0b, nz0b, nt0b, mtor
    READ(25) rminb, zminb, rmaxb, zmaxb


    ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))

    f3d(:,:,:,:) =  0.0_DP

    ALLOCATE(br(nr0b,nz0b,nt0b), bp(nr0b,nz0b,nt0b), bz(nr0b,nz0b,nt0b))

    PRINT *, ' Reading magnetic field data as MGO format'
    PRINT *

    loop010 : DO k=1,nt0b
      READ(25) ((br(i,j,k), bz(i,j,k), bp(i,j,k), i=1,nr0b), j=1,nz0b)
    END DO loop010

    f3d(1,1:nr0b,1:nz0b,1:nt0b) =  br(1:nr0b,1:nz0b,1:nt0b) !< B_R
    f3d(2,1:nr0b,1:nz0b,1:nt0b) =  bp(1:nr0b,1:nz0b,1:nt0b) !< B_phi
    f3d(3,1:nr0b,1:nz0b,1:nt0b) =  bz(1:nr0b,1:nz0b,1:nt0b) !< B_Z

    DEALLOCATE(br, bp, bz)


  END SUBROUTINE read_mgo

  SUBROUTINE read_hint_eq

    IMPLICIT NONE

    INTEGER :: nt0bh,   &
      &        nt1b,    &
      &        kstep,   &
      &        i,       &
      &        j,       &
      &        jq,      &
      &        k,       &
      &        itemp(4)
    REAL(DP) :: time,   &
      &         kpitch, &
      &         dt,     &
      &         temp
    REAL(DP), ALLOCATABLE :: br(:,:,:), &
      &                      bp(:,:,:), &
      &                      bz(:,:,:)
    CHARACTER(LEN=100) :: fmt


    SELECT CASE(TRIM(version))

      CASE('ver2')

        READ(25) time
        READ(25) nr0b, nz0b, nt1b, mtor
        READ(25) rminb, zminb, rmaxb, zmaxb, kpitch

        IF(kpitch == 0.5_DP)THEN
          PRINT *
          PRINT *,' !!! kpitch=0.5 !!!: assumed helical symmetry.'
          PRINT *
          nt0bh =  nt1b - 4
          nt0b  =  2 * (nt1b - 5)
        ELSE
          nt0b  =  nt1b - 5
          mtor  =  mtor / kpitch
        END IF

        ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))

        f3d(:,:,:,:) =  0.0_DP

        ALLOCATE(br(nr0b,nz0b,nt1b), bp(nr0b,nz0b,nt1b), bz(nr0b,nz0b,nt1b))

        PRINT *, ' Reading magnetic field data as HINT2 format'
        PRINT *

        REWIND 25

 1      CONTINUE

        READ(25,END=2) time
        READ(25) (itemp(i), i=1,4)
        READ(25) temp, temp, temp, temp, temp
        READ(25) temp

        fmt =  '(A7,ES12.4)' 

        PRINT fmt, ' time= ', time
 
        loop010 : DO k=1,nt1b
          READ(25) ((br(i,j,k), bz(i,j,k), bp(i,j,k), i=1,nr0b), j=1,nz0b)
          READ(25) ((temp, temp, temp, i=1,nr0b), j=1,nz0b)
          READ(25) ((temp, i=1,nr0b), j=1,nz0b)
        END DO loop010

        GOTO 1

 2      CONTINUE

        IF(kpitch == 0.5_DP)THEN
          f3d(1,1:nr0b,1:nz0b,1:nt0bh) =  br(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(2,1:nr0b,1:nz0b,1:nt0bh) =  bp(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(3,1:nr0b,1:nz0b,1:nt0bh) =  bz(1:nr0b,1:nz0b,3:nt1b-2)
          DO k=1,nt0bh-1
            DO j=1,nz0b
              jq =  nz0b - j + 1
              DO i=1,nr0b
                f3d(1,i,j,nt0b+2-k) = -f3d(1,i,jq,k)
                f3d(2,i,j,nt0b+2-k) =  f3d(2,i,jq,k)
                f3d(3,i,j,nt0b+2-k) =  f3d(3,i,jq,k)
              END DO
            END DO
          END DO
        ELSE
          f3d(1,1:nr0b,1:nz0b,1:nt0b+1) =  br(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(2,1:nr0b,1:nz0b,1:nt0b+1) =  bp(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(3,1:nr0b,1:nz0b,1:nt0b+1) =  bz(1:nr0b,1:nz0b,3:nt1b-2)
        END IF

      CASE DEFAULT

        READ(25) kstep
        READ(25) time
        READ(25) nr0b, nz0b, nt0b, mtor
        READ(25) rminb, zminb, rmaxb, zmaxb


        ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))

        f3d(:,:,:,:) =  0.0_DP

        ALLOCATE(br(nr0b,nz0b,nt0b), bp(nr0b,nz0b,nt0b), bz(nr0b,nz0b,nt0b))

        PRINT *, ' Reading magnetic field data as HINT3 format'
        PRINT *

        REWIND 25

 3      CONTINUE

        READ(25,END=4) kstep
        READ(25) time
        READ(25) (itemp(i), i=1,4)
        READ(25) temp, temp, temp, temp

        fmt =  '(A7,ES12.4)'

        PRINT fmt, ' time= ', time

        READ(25) (((br(i,j,k), bp(i,j,k), bz(i,j,k), i=1,nr0b), j=1,nz0b), k=1,nt0b)
        READ(25) (((temp, temp, temp, i=1,nr0b), j=1,nz0b), k=1,nt0b)
        READ(25) (((temp, i=1,nr0b), j=1,nz0b), k=1,nt0b)

        GOTO 3

 4      CONTINUE

        f3d(1,1:nr0b,1:nz0b,1:nt0b) =  br(1:nr0b,1:nz0b,1:nt0b)
        f3d(2,1:nr0b,1:nz0b,1:nt0b) =  bp(1:nr0b,1:nz0b,1:nt0b)
        f3d(3,1:nr0b,1:nz0b,1:nt0b) =  bz(1:nr0b,1:nz0b,1:nt0b)

    END SELECT

    DEALLOCATE(br, bp, bz)


  END SUBROUTINE read_hint_eq

  SUBROUTINE read_hint_vac

    IMPLICIT NONE

    INTEGER :: nt0bh, &
      &        nt1b,  &
      &        i,     &
      &        j,     &
      &        jq,    &
      &        k
    REAL(DP) :: kpitch
    REAL(DP), ALLOCATABLE :: br(:,:,:), &
      &                      bp(:,:,:), &
      &                      bz(:,:,:)
    CHARACTER(LEN=100) :: fmt


    SELECT CASE(TRIM(version))

      CASE('ver2')

        READ(25) nr0b, nz0b, nt1b, mtor
        READ(25) rminb, zminb, rmaxb, zmaxb, kpitch

        IF(kpitch == 0.5_DP)THEN
          PRINT *
          PRINT *,' !!! kpitch=0.5 !!!: assumed helical symmetry.'
          PRINT *
          nt0bh =  nt1b - 4
          nt0b  =  2 * (nt1b - 5)
        ELSE
          nt0b  =  nt1b - 5
          mtor  =  mtor / kpitch
        END IF

        ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))

        f3d(:,:,:,:) =  0.0_DP

        ALLOCATE(br(nr0b,nz0b,nt1b), bp(nr0b,nz0b,nt1b), bz(nr0b,nz0b,nt1b))

        PRINT *, ' Reading magnetic field data as legacy HINT2 format'
        PRINT *

        loop010 : DO k=1,nt1b
          READ(25) ((br(i,j,k), bz(i,j,k), bp(i,j,k), i=1,nr0b), j=1,nz0b)
        END DO loop010

        IF(kpitch == 0.5_DP)THEN
          f3d(1,1:nr0b,1:nz0b,1:nt0bh) =  br(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(2,1:nr0b,1:nz0b,1:nt0bh) =  bp(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(3,1:nr0b,1:nz0b,1:nt0bh) =  bz(1:nr0b,1:nz0b,3:nt1b-2)
          DO k=1,nt0bh-1
            DO j=1,nz0b
              jq =  nz0b - j + 1
              DO i=1,nr0b
                f3d(1,i,j,nt0b+2-k) = -f3d(1,i,jq,k)
                f3d(2,i,j,nt0b+2-k) =  f3d(2,i,jq,k)
                f3d(3,i,j,nt0b+2-k) =  f3d(3,i,jq,k)
              END DO
            END DO
          END DO
        ELSE
          f3d(1,1:nr0b,1:nz0b,1:nt0b+1) =  br(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(2,1:nr0b,1:nz0b,1:nt0b+1) =  bp(1:nr0b,1:nz0b,3:nt1b-2)
          f3d(3,1:nr0b,1:nz0b,1:nt0b+1) =  bz(1:nr0b,1:nz0b,3:nt1b-2)
        END IF

      CASE DEFAULT

        READ(25) nr0b, nz0b, nt0b, mtor
        READ(25) rminb, zminb, rmaxb, zmaxb

        ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))

        f3d(:,:,:,:) =  0.0_DP

        ALLOCATE(br(nr0b,nz0b,nt0b), bp(nr0b,nz0b,nt0b), bz(nr0b,nz0b,nt0b))

        PRINT *, ' Reading vacuum magnetic field data as HINT format'
        PRINT *

        READ(25) (((br(i,j,k), bp(i,j,k), bz(i,j,k), i=1,nr0b), j=1,nz0b), k=1,nt0b)

        f3d(1,1:nr0b,1:nz0b,1:nt0b) =  br(1:nr0b,1:nz0b,1:nt0b)
        f3d(2,1:nr0b,1:nz0b,1:nt0b) =  bp(1:nr0b,1:nz0b,1:nt0b)
        f3d(3,1:nr0b,1:nz0b,1:nt0b) =  bz(1:nr0b,1:nz0b,1:nt0b)

    END SELECT

    DEALLOCATE(br, bp, bz)


  END SUBROUTINE read_hint_vac

  SUBROUTINE read_mips_eq

    IMPLICIT NONE

    INTEGER :: itemp, &
      &        jtemp, &
      &        ktemp
    REAL(DP) :: pminb,  &
      &         pmaxb,  &
      &         dtempr, &
      &         dtempz, &
      &         dtempt
    REAL(DP), ALLOCATABLE :: br(:,:,:), &
      &                      bp(:,:,:), &
      &                      bz(:,:,:), &
      &                      p(:,:,:)


    mtor =  igrid(4)
    IF(mtor <= 0) mtor =  1

    nr0b =  igrid(1)
    nz0b =  igrid(2)
    nt0b =  igrid(3) - 4

    ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))

    f3d(:,:,:,:) =  0.0_DP

    ALLOCATE(br(nr0b,nz0b,nt0b+4), bp(nr0b,nz0b,nt0b+4), bz(nr0b,nz0b,nt0b+4), p(nr0b,nz0b,nt0b+4))

    PRINT *, ' Reading magnetic field data as MIPS format'
    PRINT *

    READ(25) itemp, jtemp, ktemp,    &
      &      rminb, rmaxb,           &
      &      zminb, zmaxb,           &
      &      pminb, pmaxb,           &
      &      dtempr, dtempz, dtempt, &
      &      br, bz, bp, p

    f3d(1,1:nr0b,1:nz0b,1:nt0b) =  br(1:nr0b,1:nz0b,3:nt0b+2)
    f3d(2,1:nr0b,1:nz0b,1:nt0b) =  bp(1:nr0b,1:nz0b,3:nt0b+2)
    f3d(3,1:nr0b,1:nz0b,1:nt0b) =  bz(1:nr0b,1:nz0b,3:nt0b+2)

    DEALLOCATE(br, bp, bz, p)


  END SUBROUTINE read_mips_eq

  SUBROUTINE read_mips

    IMPLICIT NONE

    INTEGER :: i, j, k
    REAL(DP) :: temp, bmax
    REAL(DP), ALLOCATABLE :: br(:,:,:), &
      &                      bp(:,:,:), &
      &                      bz(:,:,:)


    bmax =  2.18529336349835_DP
    !bmax =  2.891266292494260_DP

    READ(25) nr0b, nz0b, nt0b, mtor
    READ(25) rminb, zminb, rmaxb, zmaxb

    ALLOCATE(f3d(l3d,-2:nr0b+3,-2:nz0b+3,-2:nt0b+4))

    f3d(:,:,:,:) =  0.0_DP

    ALLOCATE(br(nr0b,nz0b,nt0b), bp(nr0b,nz0b,nt0b), bz(nr0b,nz0b,nt0b))

    PRINT *, ' Reading vacuum magnetic field data as HINT format'
    PRINT *

    DO k=1,nt0b
      READ(25) ((br(i,j,k), bz(i,j,k), bp(i,j,k), i=1,nr0b), j=1,nz0b)
      READ(25) ((temp, temp, temp, i=1,nr0b), j=1,nz0b)
      READ(25) ((temp, temp, i=1,nr0b), j=1,nz0b)
    END DO

    f3d(1,1:nr0b,1:nz0b,1:nt0b) =  br(1:nr0b,1:nz0b,1:nt0b) * bmax
    f3d(2,1:nr0b,1:nz0b,1:nt0b) =  bp(1:nr0b,1:nz0b,1:nt0b) * bmax
    f3d(3,1:nr0b,1:nz0b,1:nt0b) =  bz(1:nr0b,1:nz0b,1:nt0b) * bmax


    DEALLOCATE(br, bp, bz)


    RETURN
  END SUBROUTINE read_mips

  SUBROUTINE read_flux

    IMPLICIT NONE

    INTEGER :: itemp, &
      &        jtemp, &
      &        ktemp, &
      &        mtemp, &
      &        i,     &
      &        j,     &
      &        k
    REAL(DP) :: r1temp, &
      &         r2temp, &
      &         z1temp, &
      &         z2temp
    REAL(DP), ALLOCATABLE :: ss(:,:,:)
    CHARACTER(LEN=100) :: fmt


    READ(26) itemp, jtemp, ktemp, mtemp
    READ(26) r1temp, z1temp, r2temp, z2temp

    IF(itemp /= nr0b)THEN
      PRINT *, 'read_flux: parameter error  itemp =  ', itemp,  '  nr0b  = ', nr0b
      STOP
    END IF

    IF(jtemp /= nz0b)THEN
      PRINT *, 'read_flux: parameter error  jtemp =  ', jtemp,  '  nz0b  = ', nz0b
      STOP
    END IF

    IF(ktemp /= nt0b)THEN
      PRINT *, 'read_flux: parameter error  ktemp =  ', ktemp,  '  nt0b  = ', nt0b
      STOP
    END IF

    IF(mtemp /= mtor)THEN
      PRINT *, 'read_flux: parameter error  mtemp =  ', mtemp,  '  mtor  = ', mtor
      STOP
    END IF

    IF(r1temp /= rminb)THEN
      PRINT *, 'read_flux: parameter error  r1temp = ', r1temp, '  rminb = ', rminb
      STOP
    END IF

    IF(r2temp /= rmaxb)THEN
      PRINT *, 'read_flux: parameter error  r2temp = ', r2temp, '  rmaxb = ', rmaxb
      STOP
    END IF

    IF(z1temp /= zminb)THEN
      PRINT *, 'read_flux: parameter error  z1temp = ', z1temp, '  zminb = ', zminb
      STOP
    END IF

    IF(z2temp /= zmaxb)THEN
      PRINT *, 'read_flux: parameter error  z2temp = ', z2temp, '  zmaxb = ', zmaxb
      STOP
    END IF

    ALLOCATE(ss(nr0b,nz0b,nt0b))

    ss(:,:,:) =  0.0_DP

    PRINT *, ' Reading normalized toroidal flux as XSS format'
    PRINT *

    loop010 : DO k=1,nt0b
      READ(26) ((ss(i,j,k), i=1,nr0b), j=1,nz0b)
    END DO loop010

    f3d(5,1:nr0b,1:nz0b,1:nt0b) =  ss(1:nr0b,1:nz0b,1:nt0b)


    DEALLOCATE(ss)


  END SUBROUTINE read_flux

  SUBROUTINE read_eqdsk

    IMPLICIT NONE

    INTEGER :: nr, &
      &        nz, &
      &        ns, &
      &        i,  &
      &        j,  &
      &        k
    REAL(DP) :: rdim,  &
      &         zdim,  &
      &         rleft, &
      &         zmid,  &
      &         psi0,  &
      &         psia,  &
      &         temp
    REAL(DP), ALLOCATABLE :: psig(:,:)
    CHARACTER(LEN=10) :: header(6)
    CHARACTER(LEN=100) :: fmt


    PRINT *
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *, '          SUBROUTINE READ_EQDSK                                  '
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *

    READ(26,'( 6A8, 3I4 )') (header(i), i=1,6), j, nr, nz

    ns =  nr

    READ(26,'( 5E16.9 )') rdim, zdim, temp, rleft, zmid
    READ(26,'( 5E16.9 )') temp, temp, psi0, psia, temp
    READ(26,'( 5E16.9 )') temp, temp, temp, temp, temp
    READ(26,'( 5E16.9 )') temp, temp, temp, temp, temp

    PRINT *
    PRINT *,  '   COORDINATES SYSTEMS '
    PRINT *,  '  --------------------------------------------------------'
    PRINT *,  '       RMIN          RMAX          ZMIN          ZMAX '

    fmt = '(2X,4F14.9)'

    PRINT fmt, rleft, rleft + rdim, zmid - 0.5_DP * zdim, zmid + 0.5_DP * zdim

    PRINT *
    PRINT *,  '   RESOLUSIONS'
    PRINT *,  '  ------------------------'

    PRINT *,  '      NR      NZ      NS '

    fmt = '(2X,3I8)'

    PRINT fmt, nr, nz, ns

    PRINT *
    PRINT *,  '   EQUILIBRIUM PARAMETRS '
    PRINT *,  '  ----------------------------'
    PRINT *,  '       PSI0          PSIA '

    fmt = '(2x,2F14.9)'

    PRINT fmt, psi0, psia

    PRINT *
    PRINT *

    ALLOCATE(psig(nr,nz))

    READ(26,'( 5E16.9 )') (temp, i=1,ns)
    READ(26,'( 5E16.9 )') (temp, i=1,ns)
    READ(26,'( 5E16.9 )') (temp, i=1,ns)
    READ(26,'( 5E16.9 )') (temp, i=1,ns)
    READ(26,'( 5E16.9 )') ((psig(i,j), i=1,nr), j=1,nz)

    IF(nr /= nr0b)THEN
      PRINT *, 'read_eqdsk: parameter error  nr          = ', nr,                  '  nr0b  = ', nr0b
      STOP
    END IF

    IF(nz /= nz0b)THEN
      PRINT *, 'read_eqdsk: parameter error  nz          = ', nz,                  '  nz0b  = ', nz0b
      STOP
    END IF

    IF(rleft /= rminb)THEN
      PRINT *, 'read_eqdsk: parameter error  rleft       =', rleft,                '  rminb = ', rminb
      STOP
    END IF

    IF(rleft + rdim /= rmaxb)THEN
      PRINT *, 'read_eqdsk: parameter error  rleft+rdim  =', rleft + rdim,         '  rmaxb = ', rmaxb
      STOP
    END IF

    IF(zmid - 0.5_DP * zdim /= zminb)THEN
      PRINT *, 'read_eqdsk: parameter error  zmid-zdim/2 =', zmid - 0.5_DP * zdim, '  zminb = ', zminb
      STOP
    END IF

    IF(zmid + 0.5_DP * zdim /= zmaxb)THEN
      PRINT *, 'read_eqdsk: parameter error  zmid+zdim/2 =', zmid + 0.5_DP * zdim, '  zminb = ', zminb
      STOP
    END IF

    PRINT *, ' Reading normalized poloidal flux as EQDSK format'
    PRINT *

    !psig(:,:) =  psig(:,:) - psia
    !psig(:,:) =  psig(:,:) / (psia - psi0) + 1.0_DP
    psig(:,:) = 1.0_DP - (psig(:,:) - psia) / (psi0 - psia)

    !DO j=1,nz0b
    !  DO i=1,nr0b
    !    WRITE(88,'(2I5,10E15.7)') i, j, psig(i,j)
    !  END DO
    !  WRITE(88,*)
    !END DO
    !STOP

    DO k=1,nt0b
      f3d(5,1:nr0b,1:nz0b,k) =  psig(1:nr0b,1:nz0b)
    END DO

    DEALLOCATE(psig)


  END SUBROUTINE read_eqdsk

  SUBROUTINE read_eqdata

    IMPLICIT NONE
    INTEGER :: irdm,      &
      &        izdm2,     &
      &        irzdm2,    &
      &        ivdm,      &
      &        i,         &
      &        j,         &
      &        k,         &
      &        itemp(10)
    REAL(DP) :: saxis, &
      &         temp
    REAL(DP), ALLOCATABLE :: rr(:),    &
      &                      zz(:),    &
      &                      psi(:),   &
      &                      psig(:,:)
    CHARACTER(LEN=10) :: header(6)


    READ(26,*) itemp(1), temp
    READ(26,*) irdm, izdm2, irzdm2, itemp(1), itemp(2), ivdm, itemp(3), itemp(4), itemp(5)

    ALLOCATE(rr(irdm), zz(izdm2), psi(irzdm2))

    READ(26,*) (psi(i), temp, i=1,irzdm2)
    READ(26,*) (rr(i), i=1,irdm)
    READ(26,*) (temp, i=1,irdm)
    READ(26,*) (zz(i), i=1,izdm2)

    READ(26,*) (temp, i=1,ivdm)
    READ(26,*) (temp, temp, temp, i=1,ivdm)

    READ(26,*) temp, temp, temp, saxis, temp, temp, temp, temp, temp, temp, temp, temp, temp, temp, temp

    ALLOCATE(psig(irdm,izdm2))

    k =  0
    DO j=1,izdm2
      DO i=1,irdm
        k         =  k + 1
        psig(i,j) =  psi(k)
      END DO
    END DO

    IF(irdm /= nr0b)THEN
      PRINT *, 'read_eqdata: parameter error  irdm      = ', irdm,     '  nr0b  = ', nr0b
      STOP
    END IF

    IF(izdm2 /= nz0b)THEN
      PRINT *, 'read_eqdata: parameter error  izdm2     = ', izdm2,    '  nz0b  = ', nz0b
      STOP
    END IF

    IF(rr(1) /= rminb)THEN
      PRINT *, 'read_eqdata: parameter error  rr(1)     =', rr(1),     '  rminb = ', rminb
      STOP
    END IF

    IF(rr(irdm) /= rmaxb)THEN
      PRINT *, 'read_eqdata: parameter error  rr(irdm)  =', rr(irdm),  '  rmaxb = ', rmaxb
      STOP
    END IF

    IF(zz(1) /= zminb)THEN
      PRINT *, 'read_eqdata: parameter error  zz(1)     =', zz(1),     '  zminb = ', zminb
      STOP
    END IF

    IF(zz(izdm2) /= zmaxb)THEN
      PRINT *, 'read_eqdata: parameter error  zz(izdm2) =', zz(izdm2), '  zminb = ', zminb
      STOP
    END IF

    PRINT *, ' Reading normalized poloidal flux as EQDATA format'
    PRINT *

    psig(:,:) =  psig(:,:) - saxis
    psig(:,:) =  psig(:,:) / ABS(saxis)

    DO k=1,nt0b
      f3d(5,1:nr0b,1:nz0b,k) =  psig(1:nr0b,1:nz0b)
    END DO

    DEALLOCATE(rr, zz, psi, psig)


  END SUBROUTINE read_eqdata

  SUBROUTINE free_mem_field

    IMPLICIT NONE


    DEALLOCATE(f3d)


  END SUBROUTINE free_mem_field

  SUBROUTINE magset

    IMPLICIT NONE

    INTEGER :: i,  &
      &        j,  &
      &        jq, &
      &        k,  &
      &        l
    REAL(DP) :: r,  &
      &         z,  &
      &         bb, &
      &         eps
    CHARACTER(LEN=100) :: fmt


    PRINT *
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *, '          SUBROUTINE magset                                      '
    PRINT *, ' ----------------------------------------------------------------'
    PRINT *


    CALL read_field

    ALLOCATE(rg(nr0b), zg(nz0b))

    pi2m  =  pi2 / mtor

    nx3d =  nr0b
    ny3d =  nz0b
    nz3d =  nt0b + 1

    delrb = (rmaxb - rminb) / (nr0b - 1)
    delzb = (zmaxb - zminb) / (nz0b - 1)

    DO i=1,nr0b
      rg(i) =  rminb + delrb * (i - 1)
    END DO

    DO j=1,nz0b
      zg(j) =  zminb + delzb * (j - 1)
    END DO
 
    CALL splin3(rminb, rmaxb, zminb, zmaxb, 0.0_DP, pi2m)

    f3d(4,:,:,:) =  SQRT(f3d(1,:,:,:)**2 + f3d(2,:,:,:)**2 + f3d(3,:,:,:)**2)

    bmax =  5.0_DP

    DO k=1,nt0b
      DO j=1,nz0b
        DO i=1,nr0b
          bb =  f3d(4,i,j,k)
          IF(bb > bmax)THEN
            f3d(1,i,j,k) =  bmax * f3d(1,i,j,k) / bb
            f3d(2,i,j,k) =  bmax * f3d(2,i,j,k) / bb
            f3d(3,i,j,k) =  bmax * f3d(3,i,j,k) / bb
          END IF
        END DO
      END DO
    END DO

    DO k=1,nt0b
      DO j=1,nz0b
        DO i=1,3
          DO l=1,l3d
            r =  rminb - delrb * (4 - i)
            CALL polint(rg(1:4),         f3d(l,1:4,j,k),         4, r, f3d(l,i-3,j,k),    eps)
            r =  rmaxb + delrb * i
            CALL polint(rg(nr0b-3:nr0b), f3d(l,nr0b-3:nr0b,j,k), 4, r, f3d(l,nr0b+i,j,k), eps)
          END DO
        END DO
      END DO
    END DO

    DO k=1,nt0b
      DO i=-2,nr0b+3
        DO j=1,3
          DO l=1,l3d
            z =  zminb - delzb * (4 - j)
            CALL polint(zg(1:4),         f3d(l,i,1:4,k),         4, z, f3d(l,i,j-3,k),    eps)
            z =  zmaxb + delzb * j
            CALL polint(zg(nz0b-3:nz0b), f3d(l,i,nz0b-3:nz0b,k), 4, z, f3d(l,i,nz0b+j,k), eps)
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE(rg, zg)

    IF(nt0b == 1)THEN
      f3d(:,:,:,nt0b+1) =  f3d(:,:,:,1)
      f3d(:,:,:,nt0b+2) =  f3d(:,:,:,1)
      f3d(:,:,:,nt0b+3) =  f3d(:,:,:,1)
      f3d(:,:,:,nt0b+4) =  f3d(:,:,:,1)

      f3d(:,:,:,-2)     =  f3d(:,:,:,1)
      f3d(:,:,:,-1)     =  f3d(:,:,:,1)
      f3d(:,:,:,0)      =  f3d(:,:,:,1)
    ELSE
      IF(lsymmetry)THEN

        IF(mod(nt0b,2) /= 0)THEN
          PRINT *
          PRINT *, ' parameter nt0b is NOT even number!!!'
          PRINT *
          STOP
        END IF

        f3d(:,:,:,1)        =  0.5_DP * (f3d(:,:,:,2)      + f3d(:,:,:,nt0b))
        f3d(:,:,:,nt0b/2+1) =  0.5_DP * (f3d(:,:,:,nt0b/2) + f3d(:,:,:,nt0b/2+2))
      END IF 

      f3d(:,:,:,nt0b+1) =  f3d(:,:,:,1)
      f3d(:,:,:,nt0b+2) =  f3d(:,:,:,2)
      f3d(:,:,:,nt0b+3) =  f3d(:,:,:,3)
      f3d(:,:,:,nt0b+4) =  f3d(:,:,:,4)

      f3d(:,:,:,-2)     =  f3d(:,:,:,nt0b-2)
      f3d(:,:,:,-1)     =  f3d(:,:,:,nt0b-1)
      f3d(:,:,:,0)      =  f3d(:,:,:,nt0b)
    END IF

    PRINT *
    PRINT *, ' COORDINATES SYSTEMS '
    PRINT *, '-----------------------------------------------------------------'
    PRINT *, '   M    Rmin     Rmax     Zmin     Zmax '

    fmt = '(I5, 4F9.4)'

    PRINT fmt, mtor, rminb, rmaxb, zminb, zmaxb

    PRINT *
    PRINT *, ' RESOLUSIONS'
    PRINT *, '----------------------------------------------------------------------'
    PRINT *, '   nr0b    nt0b    nz0b   delrb    dphi     delzb             '

    fmt = '(3I8, 4F9.4)'
    PRINT fmt, nr0b, nt0b, nz0b, delrb, pi2m / nt0b, delzb
    PRINT *


  END SUBROUTINE magset

  SUBROUTINE mgval1 (r, phi, z,     & ! (in)
    &                br, bp, bz, bb & ! (out)
    &               )

    IMPLICIT NONE

    REAL(DP), INTENT(IN)  :: r,   & ! major radius in computational region [m]
      &                      phi, & ! toroidal angle in computational region [rad]
      &                      z      ! height in computational region [m]
    REAL(DP), INTENT(OUT) :: br, & ! BR component on (R,phi,Z
      &                      bp, & ! Bphi component on (R,phi,Z)
      &                      bz, & ! BZ component on (R,phi,Z)
      &                      bb    ! B strength on (R,phi,Z) [T]
    INTEGER :: iphi
    REAL(DP) :: phi1,   &
      &         xd(3),  &
      &         w0(l3d)


    iphi =  phi / pi2m
    phi1 =  phi - pi2m * iphi
    IF(phi1 < 0.0_DP) phi1 =  phi1 + pi2m
    IF(phi1 >= pi2m)  phi1 =  phi1 - pi2m       

    xd(1) =  r
    xd(2) =  z
    xd(3) =  phi1

    CALL spl3df(xd(:), w0(:))

    br =  w0(1)
    bp =  w0(2)
    bz =  w0(3)
    bb =  w0(4)


    RETURN
  END SUBROUTINE mgval1

  SUBROUTINE mgval2 (r, phi, z,          & ! (in)
    &                b, dbdr, dbdp, dbdz & ! (out)
    &               )

    IMPLICIT NONE
!Arguments
    REAL(DP), INTENT(IN)  :: r,         & ! major radius in computational region [m]
      &                      phi,       & ! toroidal angle in computational region [rad]
      &                      z            ! height in computational region [m]
    REAL(DP), INTENT(OUT) :: b(l3d),    & ! B vector (1:BR,2:BT,3:BZ,4:B)
      &                      dbdr(l3d), & ! R-derivative of B vector (1:BR,2:BZ,3:BT,4:B)
      &                      dbdp(l3d), & ! phi-derivative of B vector (1:BR,2:BZ,3:BT,4:B)
      &                      dbdz(l3d)    ! Z-derivative of B vector (1:BR,2:BZ,3:BT,4:B)
!Local variables
    INTEGER :: iphi
    REAL(DP) :: phi1,    &
      &         xd(3),   &
      &         w0(l3d), &
      &         wx(l3d), &
      &         wy(l3d), &
      &         wz(l3d)


    iphi =  phi / pi2m
    phi1 =  phi - pi2m * iphi
    IF(phi1 < 0.0_DP) phi1 =  phi1 + pi2m
    IF(phi1 >= pi2m)  phi1 =  phi1 - pi2m      

    xd(1) =  r
    xd(2) =  z
    xd(3) =  phi1

    CALL spl3dd(xd(:), w0(:), wx(:), wy(:), wz(:))

    b(:)    =  w0(:)
    dbdr(:) =  wx(:)
    dbdz(:) =  wy(:)
    dbdp(:) =  wz(:)


    RETURN
  END SUBROUTINE mgval2

  SUBROUTINE mgval3 (r, phi, z, & ! (in)
    &                s          & ! (out)
    &               )

    IMPLICIT NONE

    REAL(DP), INTENT(IN)  :: r,   & ! major radius in computational region [m]
      &                      phi, & ! toroidal angle in computational region [rad]
      &                      z      ! height in computational region [m]
    REAL(DP), INTENT(OUT) :: s     ! flux label

    INTEGER :: iphi
    REAL(DP) :: phi1,  &
      &         xd(3),  &
      &         w0(l3d)


    iphi =  phi / pi2m
    phi1 =  phi - pi2m * iphi
    IF(phi1 < 0.0_DP) phi1 =  phi1 + pi2m
    IF(phi1 >= pi2m)  phi1 =  phi1 - pi2m

    xd(1) =  r
    xd(2) =  z
    xd(3) =  phi1

    CALL spl3df(xd(:), w0(:))

    s =  w0(5)


    RETURN
  END SUBROUTINE mgval3

  SUBROUTINE polint (xa, ya, n, x, & ! (in)
    &                y, dy         & ! (out)
    &               )

    IMPLICIT NONE

!Arguments
    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: x,     &
      &                     xa(n), &
      &                     ya(n)
    REAL(DP), INTENT(OUT) :: y,  &
    &                        dy
!Local variables
    INTEGER :: m,  &
      &        ns, &
      &        i
    REAL(DP) :: den,  &
      &         dif,  &
      &         dift, &
      &         ho,   &
      &         hp,   &
      &         w,    &
      &         c(n), &
      &         d(n)
  

    ns  =  1
    dif =  ABS(x - xa(1))

    loop100 : DO i=1,n
      dift =  ABS(x - xa(i))
      IF(dift < dif)THEN
        ns  =  i
        dif =  dift
      END IF
      c(i) =  ya(i)
      d(i) =  ya(i)
    END DO loop100

    y  =  ya(ns)
    ns =  ns - 1

    loop200 : DO m=1,n-1
      loop210 : DO i=1,n-m
        ho  =  xa(i)   - x
        hp  =  xa(i+m) - x
        w   =  c(i+1)  - d(i)
        den =  ho      - hp
        IF(den == 0.0_DP) STOP 'failure in polint'
        den  =  w  / den
        d(i) =  hp * den
        c(i) =  ho * den
      END DO loop210
      IF(2 * ns < n - m)THEN
        dy =  c(ns+1)
      ELSE
        dy =  d(ns)
        ns =  ns - 1
      END IF
      y =  y + dy
    END DO loop200


    RETURN
  END SUBROUTINE polint

END MODULE cylindrical_coord_mod
!===============================================================================
