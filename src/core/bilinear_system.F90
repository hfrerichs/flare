!=======================================================================
! solve xA + s*xB + t*xC + s*t*xD = 0
! with boundary check (bc), i.e. -1<s<1 and 0<t<1
! input:        coefficients xA,xB,xC,xD
! output:	y=(s1,t1) ,y22=(s2,t2)
!		istat = number of valid solutions
!=======================================================================
  subroutine solve_bilinear_system_bc (xA, xB, xC, xD, y, y22, istat)
  real*8, intent(in)  :: xA(2), xB(2), xC(2), xD(2)
  real*8, intent(out) :: y(2), y22(2)
  integer, intent(inout) :: istat

  real*8 :: y1(2), y2(2)
  integer :: is1, is2


  call solve_bilinear_system (xA, xB, xC, xD, y1, y2, istat)
  y = 0.d0
  y22 = 0.d0

  ! no solution found
  if (istat == 0) return


  ! check if solutions are within boundaries
  is1 = 0
  is2 = 0
  ! check first solution
  if (y1(1).ge.-1.d0 .and. y1(1).le.1.d0 .and. &
      y1(2).ge.0.d0 .and. y1(2).le.1.d0) then
     is1 = 1
  endif
  ! check second solution (if it exists)
  if (istat == 2 .and. &
      y2(1).ge.-1.d0 .and. y2(1).le.1.d0 .and. &
      y2(2).ge.0.d0 .and. y2(2).le.1.d0) then
     is2 = 1
  endif


  ! one solution found, and solution is good?
  if (istat == 1) then
     if (is1 == 1) then
        y = y1
        return
     else
        istat = 0
        return
     endif
  endif


  ! two solutions found
  ! both solutions are good?
  if (is1 == 1 .and. is2 == 1) then
     y = y1
     y22 = y2
     !if (y2(2) < y1(2)) y = y2
     return

  ! only 1st solution is good?
  elseif (is1 == 1 .and. is2 == 0) then
     istat = 1
     y = y1
     return

  ! only 2nd solution is good?
  elseif (is1 == 0 .and. is2 == 1) then
     istat = 1
     y = y2
     return

  ! none of the solutions are good?
  else
     istat = 0
     return
  endif

  return
  end subroutine solve_bilinear_system_bc
!=======================================================================



!=======================================================================
! solve xA + s*xB + t*xC + s*t*xD = 0
! input:        coefficients xA,xB,xC,xD
! output:	y1=(s1,t1) ,y2=(s2,t2)
!		istat = number of solutions
!=======================================================================
  subroutine solve_bilinear_system (xA, xB, xC, xD, y1, y2, istat)
  real*8, intent(in)     :: xA(2), xB(2), xC(2), xD(2)
  real*8, intent(out)    :: y1(2), y2(2)
  integer, intent(inout) :: istat

  real*8  :: At, Bt, Ct, Dt, H(2)
  integer :: I


  y1 = 0.d0
  y2 = 0.d0

  At = xD(1)*xC(2) - xC(1)*xD(2)
  Bt = xD(1)*xA(2) - xC(1)*xB(2) + xB(1)*xC(2) - xA(1)*xD(2)
  Ct = xB(1)*xA(2) - xA(1)*xB(2)

  ! Linear elements
  if (At.eq.0.d0) then
     if (Bt.ne.0.d0) then
        istat = 1
        y1(2) = - Ct / Bt

        ! to calculate s = - (xA + t*xC) / (xB + t*xD), choose component I with larger |xB_I|
        I = 2
        if (abs(xB(1)) .gt. abs(xB(2))) I = 1
        y1(1) = - (xA(I) + y1(2)*xC(I)) / (xB(I) + y1(2)*xD(I))
#if defined(DEBUG)
  write (6, *) '1 solution: '
  write (6, *) 's1,t1 = ', y1
#endif
        return
     endif
     istat = 0
     return
  endif


  ! solver for 1/t, which should be more stable for small non-linearities A
  Dt = Bt**2 / Ct**2 / 4.d0 - At / Ct
  if (Dt.eq.0.d0) then
     istat = 1
     y1(2) = - Bt / Ct
  elseif (Dt.gt.0.d0) then
     istat = 2
     y1(2) = - Bt / Ct / 2.d0 + sqrt(Dt)
     y1(2) = 1.d0 / y1(2)
     y2(2) = - Bt / Ct / 2.d0 - sqrt(Dt)
     y2(2) = 1.d0 / y2(2)
     H(1) = xB(1) + y2(2)*xD(1)
     H(2) = xB(2) + y2(2)*xD(2)
     if (abs(H(1)).gt.abs(H(2))) then
        y2(1) = - (xA(1) + y2(2)*xC(1)) / H(1)
     else
        y2(1) = - (xA(2) + y2(2)*xC(2)) / H(2)
     endif
  else
     istat = 0
     return
  endif

  H(1) = xB(1) + y1(2)*xD(1)
  H(2) = xB(2) + y1(2)*xD(2)
  if (abs(H(1)).gt.abs(H(2))) then
     y1(1) = - (xA(1) + y1(2)*xC(1)) / H(1)
  else
     y1(1) = - (xA(2) + y1(2)*xC(2)) / H(2)
  endif

#if defined(DEBUG)
  write (6, *) '2 solutions: '
  write (6, *) 's1,t1 = ', y1
  write (6, *) 's2,t2 = ', y2
#endif
  return
  end subroutine solve_bilinear_system
!=======================================================================
