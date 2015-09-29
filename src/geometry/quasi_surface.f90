!===============================================================================
! Generate homogeneous magnetic surfaces
!===============================================================================
module quasi_surface
  use iso_fortran_env
  use curve2D
  implicit none
  private


  type, extends(t_curve) :: t_quasi_surface
     ! toroidal location
     real(real64) :: Phi

     ! average poloidal flux on surface
     !real(real64) :: PsiN

     ! homogeneity parameter
     real(real64) :: H

     ! reference position (within tolerance)
     real(real64) :: Rout

     contains
     procedure :: generate
     procedure :: smooth_plot
! TODO: spline representation with weights from upt. flux surface
  end type t_quasi_surface


  public :: t_quasi_surface

  contains
!=======================================================================



!=======================================================================
! Function:	generate magnetic surface
!
! Input:
!   r		reference position (R[cm], Z[cm], Phi[rad])
!   tol		allowed tolerance for position of magnetic surface
!   n_sym	toroidal symmetry number
!   Trace_Step
!=======================================================================
  subroutine generate(this, r, tol, N_sym, Trace_Step)
  class(t_quasi_surface)   :: this
  real(real64), intent(in) :: r(3), tol, Trace_Step
  integer,      intent(in) :: n_sym

  integer, parameter :: &
     ntrac_p = 200, &
     np_min  =  50

  real(real64), dimension(0:ntrac_p) :: Rtarg, Ztarg
  real(real64) :: a, d, H
  integer      :: i, it, j, k, np_out


  this%Phi = r(3)
  call homo_dist_P(n_sym, Trace_Step, r(3), r(1), r(2), tol, ntrac_p, &
                   Rtarg, Ztarg, np_min, np_out, H)
  this%H   = H


  call this%new(np_out)
  k  = 0
  it = np_out - 1
  this%x(0,1) = Rtarg(0); this%x(0,2) = Ztarg(0)
  do i=0,np_out-2
     Rtarg(k:it-1) = Rtarg(k+1:it); Ztarg(k:it-1) = Ztarg(k+1:it)
     it = it - 1
     d  = 1.d60
     do j=0,it
        a = (Rtarg(j)-this%x(i,1))**2 + (Ztarg(j)-this%x(i,2))**2
        if (a < d) then
           d = a
           k = j
        endif
     enddo
     this%x(i+1,1) = Rtarg(k); this%x(i+1,2) = Ztarg(k)
  enddo
  this%x(np_out,:) = this%x(0,:)
  this%Rout        = this%x(0,1)

  end subroutine generate
!=======================================================================



!=======================================================================
  subroutine bline_trace(H,PC,RC,ZC,TL,NP,RT,ZT,IERR)
  use fieldline
  real(real64), intent(in)                   :: H,PC,RC,ZC,TL
  integer,      intent(in)                   :: NP
  real(real64), dimension(0:NP), intent(out) :: RT, ZT
  integer,                       intent(out) :: IERR

  type(t_fieldline) :: F
  real(real64)      :: HL, PHI, R3(3)
  integer           :: NFL, I, j


  NFL = TL/ABS(H) 
  HL  = TL/NFL; IF(H<0.) HL=-HL
 
  RT(0) = RC; ZT(0) = ZC
  R3(1) = RC; R3(2) = ZC; R3(3) = PC
  call F%init(R3, HL, NM_AdamsBashforthMoulton, FL_ANGLE)
  DO I=1,NP 
     do j=1,NFL
        call F%trace_1step()
     enddo
     RT(I) = F%yc(1); ZT(I) = F%yc(2)
  ENDDO
  IERR = 0

  end subroutine bline_trace
!=======================================================================



!=======================================================================
! Name:		homo_dist_P
! Function:	Look for a magnetic surface in a range R0 +/- TOL
!		so that the locations of the field-line over NTRAC_P
!		field periods are as homogeinious as possible.
!
! Input:
!   NSYM	toroidal symmetry number
!   HPAS	approximate trace step
!   PHI0, R0, Z0	reference coordinates on magnetic surface
!   TOL		tolerance for position of magnetic surface
!   NTRAC_P	number of field periods to trace field lines
!   NP_MIN	minimum number of points on magnetic surface
!
! Output:
!   NP_OUT	number of points of magnetic surface
!   HOMO	homogeneity parameter
!   RTARG, ZTARG	coordinates of points on magnetic surface
!=======================================================================
  subroutine homo_dist_P(NSYM,HPAS,PHI0,R0,Z0,TOL,NTRAC_P,RTARG,ZTARG,NP_MIN,NP_OUT,HOMO)
  use math
  real(real64), intent(in)  :: HPAS, PHI0, R0, Z0, TOL
  integer,      intent(in)  :: NSYM, NTRAC_P,NP_MIN
  integer,      intent(out) :: NP_OUT
  real(real64), intent(out) :: HOMO
  real(real64), dimension(0:NTRAC_P), intent(out) :: RTARG,ZTARG


  real(real64) :: deviation(-5:5), Rin(-5:5)
  real(real64) :: TRACE_L,R_B,R_TRY,D_B,D_TRY,STEP,R_MIN,D_MIN
  integer      :: IERR,I,I_MIN
  intrinsic    :: MINLOC,MINVAL


  TRACE_L = pi2 / NSYM
  STEP    = TOL / FLOAT(5)
  DO I=-5,5
     RIN(I) = R0+FLOAT(I)*STEP
     CALL BLINE_TRACE(HPAS,PHI0,RIN(I),Z0,TRACE_L,NTRAC_P,RTARG,ZTARG,IERR)
     DEVIATION(I) = HOMOGENEITY(NTRAC_P+1,RTARG,ZTARG,NP_MIN,NP_OUT)
  ENDDO

  I_MIN = minloc(DEVIATION,1)
  I_MIN = I_MIN - 6
  D_MIN = DEVIATION(I_MIN)
  R_MIN = RIN(I_MIN)


  IF (D_MIN >.01  .AND.  ABS(I_MIN) /= 5) THEN
     IF (DEVIATION(I_MIN-1) < DEVIATION(I_MIN+1)) THEN
        D_B = DEVIATION(I_MIN-1)
        R_B = RIN      (I_MIN-1)
     ELSE
        D_B = DEVIATION(I_MIN+1)
        R_B = RIN      (I_MIN+1)
     ENDIF
     LOOK_FOR : DO
        R_TRY = 0.5D0*(R_B+R_MIN)
        CALL BLINE_TRACE(HPAS,PHI0,R_TRY,Z0,TRACE_L,NTRAC_P,RTARG,ZTARG,IERR)
        D_TRY = HOMOGENEITY(NTRAC_P+1,RTARG,ZTARG,NP_MIN,NP_OUT)
        IF (ABS(D_TRY-D_MIN) < 0.001) EXIT LOOK_FOR
        IF (D_TRY < D_MIN) THEN
           R_B   = R_MIN
           D_B   = D_MIN
           R_MIN = R_TRY
           D_MIN = D_TRY
        ELSE
           R_B   = R_TRY
           D_B   = D_TRY
        ENDIF
     ENDDO LOOK_FOR
  ENDIF


  CALL BLINE_TRACE(HPAS,PHI0,R_MIN,Z0,TRACE_L,NTRAC_P,RTARG,ZTARG,IERR)
  HOMO = HOMOGENEITY(NTRAC_P+1,RTARG,ZTARG,NP_MIN,NP_OUT)

  end subroutine homo_dist_P
!=======================================================================



!=======================================================================
! Name:		homogeneity
! Function:	optimize distribution of points on a closed surface
!
! Input:   
!   np		number of points
!   Ri,Zi	R,Z-coordinates of the points
!   np_min	minimum points taken into account
!
! Output:
!   np_out	best homogeneity found for the first np_out points
!=======================================================================
  function homogeneity(np, Ri, Zi, np_min, np_out) result(H)
  integer,      intent(in)  :: np, np_min
  real(real64), intent(in)  :: Ri(np), Zi(np)
  integer,      intent(out) :: np_out
  real(real64)              :: H

  real(real64) :: R(np), Z(np), dist(np)
  real(real64) :: A, D, H_min, Rl, Zl
  integer      :: i, j, il, mp, nnp


  do nnp=np_min,np
     R  = Ri
     Z  = Zi
     il = 1
     mp = nnp

     do i=1,nnp-1
        Rl = R(il); Zl = Z(il)
        do j=il,mp-1
           R(j) = R(j+1); Z(j) = Z(j+1)
        enddo

        mp = mp-1
        D  = 1.d60
        do j=1,mp
           A = (R(j)-Rl)**2 + (Z(j)-Zl)**2
           if (A < D) then
              D  = A
              il = j
           endif
        enddo 
        dist(i) = sqrt(D)
     enddo

     H = 0.d0
     do i=1,nnp-2
        H = H + exp(abs((dist(i+1)-dist(i))/(dist(i+1)+dist(i)))/0.2d0)-1.d0
     enddo
     H = H / float(nnp-2)
     if (nnp == np_min) then 
        H_min  = H
        np_out = nnp
     elseif (H < H_min) then
        H_min  = H
        np_out = nnp
     endif
 
  enddo 
  H = H_min

  end function homogeneity
!=======================================================================



!=======================================================================
  subroutine smooth_plot(this, filename, iu, nsample)
  use equilibrium, only: get_PsiN, get_poloidal_angle
  use math
  use cspline
  class(t_quasi_surface)                 :: this
  character(len=*), intent(in), optional :: filename
  integer,          intent(in), optional :: iu, nsample

  type(t_cspline) :: S
  real(real64)    :: t, x(2), theta, PsiN
  integer         :: i, iu1, n


  ! set unit number for output
  iu1 = 50
  if (present(iu)) iu1 = iu

  ! set number of sample points
  n   = 1000
  if (present(nsample)) n = nsample

  ! open output file
  if (present(filename)) open  (iu1, file=filename)


  call S%setup(this%nodes, periodic=.true.)
  do i=0,n
     t     = 1.d0 * i / n
     x     = S%eval(t)
     PsiN  = get_PsiN(x)
     theta = get_poloidal_angle(x) / pi * 180.d0
     if (theta < 0.d0) theta = theta + 360.d0
     write (iu1, *) x, theta, PsiN
  enddo


  ! close output file
  if (present(filename)) close (iu1)

  end subroutine smooth_plot
!=======================================================================

end module quasi_surface
