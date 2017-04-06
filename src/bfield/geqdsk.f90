!===============================================================================
! Axisymmetric MHD equilibrium, data is provided in G-EQDSK format
!===============================================================================
module geqdsk
  use iso_fortran_env
  use magnetic_axis
  implicit none
  private


  ! spline interpolation order
  integer :: nord = 5


  real(real64), dimension(:,:), allocatable :: psirz, Psicoeff
  real(real64), dimension(:),   allocatable :: fpol, pres, ffprim, pprime, &
     qpsi, fpolcoeff, prescoeff, ffprimcoeff, pprimecoeff, &
     REQD, ZEQD, PsinEQD, rbbbs, zbbbs
  real(real64), dimension(:),   allocatable, target :: rlim, zlim

  real(real64) :: Rdim, Zdim, Rcentr, Rleft, Zmid, Current, Rmaxis, Zmaxis, Simag, Sibry, Bcentr
  real(real64) :: Bt_scale, Ip_scale
  integer      :: nR, nZ, nbbbs, limitr



  public :: &
     geqdsk_load, &
     geqdsk_broadcast, &
     geqdsk_get_Bf, &
     geqdsk_get_JBf,&
     geqdsk_get_Psi, &
     geqdsk_get_DPsi, &
     geqdsk_get_pressure, &
     geqdsk_get_domain, &
     geqdsk_export_boundary, &
     geqdsk_info

  contains
!===============================================================================



!===============================================================================
! Load equilibrum data from file "Data_File"
!
! required input:
! Ip, Bt               Equilibrium can be scaled to user defined plasma current
!                      and toroidal magnetic field (this is sort of optional, no
!                      scaling is done if = 0)
! CurrentFix           Use sign(Current) to determine the plasma current /
!                      poloidal field direction
!
! output:
! Psi_axis, Psi_sepx   Poloidal magnetic flux on axis and separatrix
!
! optional input:
! Header_Format        = STRICT, FREE:
!===============================================================================
  subroutine geqdsk_load(Data_File, Ip, Bt, CurrentFix, &
                         Psi_axis, Psi_sepx, Header_Format)
  use numerics, only: Spline_Order
  use bspline
  use system
  character(len=*), intent(in)  :: Data_File
  real(real64),     intent(inout)  :: Ip, Bt
  logical,          intent(in)  :: CurrentFix
  real(real64),     intent(out) :: Psi_axis, Psi_sepx
  integer,          intent(in), optional :: Header_Format

  integer, parameter :: iu = 17

  real(real64), dimension(:), allocatable :: Rtmp, Ztmp, Psintmp

  character(len=80)  :: tmp
  character(len=8)   :: case_(6)
  real(real64) :: xdum, PsiC, PsiL, PsiR, rc_fix
  integer      :: i, j, idum, iformat
  logical      :: lL, lR, concaveup


  ! 0. set defaults for optional input .................................
  ! 0.1. header format: STRICT or FREE
  iformat = STRICT
  if (present(Header_Format)) iformat = Header_Format
  select case(iformat)
  case(STRICT,FREE)
  case default
     write (6, 9000) iformat
     stop
  end select
 9000 format('error: header format ', i0, ' not defined!')
  ! 0.2. set spline order
  nord    = Spline_Order
  !.....................................................................


  ! 1. read data file ..................................................
  open  (iu, file=Data_File)
  ! 1.1. read header
  select case(iformat)
  case(STRICT)
     read  (iu, 2000) (case_(i),i=1,6), idum, nR, nZ
  case(FREE)
     read  (iu, 2001) (case_(i),i=1,6), tmp
     read  (tmp, *) idum, nR, nZ
  end select
  read  (iu, 2020) Rdim, Zdim, Rcentr, Rleft, Zmid
  read  (iu, 2020) Rmaxis, Zmaxis, Simag, Sibry, Bcentr
  read  (iu, 2020) Current, Simag, xdum, Rmaxis, xdum
  read  (iu, 2020) Zmaxis, xdum, Sibry, xdum, xdum
  call setup_magnetic_axis_2D (Rmaxis*1.d2, Zmaxis*1.d2)

  ! print characteristic values
  write (6, 3000) adjustl(case_(1)), case_(2:6)
  write (6, *)
  write (6, 3001) nR, nZ
  write (6, 3002) Rleft*1.e2, (Rleft+Rdim)*1.d2, (Zmid-Zdim/2.d0)*1.d2, (Zmid+Zdim/2.d0)*1.d2
  write (6, 3003) Rcentr*1.d2
  write (6, 3004) Bcentr

  ! prepare data arrays
  allocate (fpol(nR), pres(nR), ffprim(nR), pprime(nR), qpsi(nR))
  allocate (psirz(nR,nZ))


  ! 1.2. read equilibrium data
  read  (iu, 2020) (fpol(i),i=1,nR)
  read  (iu, 2020) (pres(i),i=1,nR)
  read  (iu, 2020) (ffprim(i),i=1,nR)
  read  (iu, 2020) (pprime(i),i=1,nR)
  read  (iu, 2020) ((psirz(i,j),i=1,nR),j=1,nZ)
  read  (iu, 2020) (qpsi(i),i=1,nR)

  read  (iu, 2022) nbbbs, limitr
  allocate (rbbbs(nbbbs), zbbbs(nbbbs), rlim(limitr), zlim(limitr))
  read  (iu, 2020) (rbbbs(i), zbbbs(i), i=1,nbbbs)
  read  (iu, 2020) (rlim(i), zlim(i), i=1,limitr)
  close (iu)


  ! 1.3. (optional) equilibrium scale factors
  Bt_scale = 1.d0
  if (Bt .ne. 0.d0) then
     Bt_scale = Bt / Bcentr
     Bcentr   = Bcentr * Bt_scale
  else
     Bt       = Bcentr
  endif
  Ip_scale = 1.d0
  if (Ip .ne. 0.d0) then
     Ip_scale = Ip / Current
     Current  = Current * Ip_scale
     Simag    = Simag * Ip_scale
     Sibry    = Sibry * Ip_scale
  else
     Ip       = Current
  endif
  fpol  = fpol  * Bt_scale
  psirz = psirz * Ip_scale
  !.....................................................................


  ! 2. (optional) set plasma current / poloidal field direction ........
  ! implementation based on TRIP3D
  rc_fix    = 1.d0
  ! 2.1. Make psirz concave up, i.e. psi_center > psi_edge
  lL        = .true.
  lR        = .true.
  concaveup = .true.
  PsiC      = psirz(nR/2,nZ/2)
  PsiL      = psirz(nR/10, nZ/2)
  PsiR      = psirz(9*nR/10, nZ/2)
  if (PsiC .gt. PsiL) lL = .false.
  if (PsiC .gt. PsiR) lR = .false.
  if ((lL.eqv. .false.) .and. (lR.eqv. .false.)) concaveup = .false.

  ! flip psirz if not concave up
  if (concaveup.eqv. .false.) then
     if (CurrentFix) then
        psirz  = - psirz
        Simag  = - Simag
        Sibry  = - Sibry
     else
        rc_fix = -1.d0
     endif
  endif
  ! now psirz is concave up, or rc_fix = -1

  ! 2.2 Current fix: (+) -> flip to concave down / convex, (-) -> keep concave up
  if (Current .gt. 0.d0) then
     if (CurrentFix) then
        psirz  = - psirz
        Simag  = - Simag
        Sibry  = - Sibry
     else
        rc_fix = -1.d0 * rc_fix
     endif
  endif

  ! 2.3 update Current to match psirz configuration, if psirz is not modified
  if (.not.CurrentFix) then
     Current = Current * rc_fix
     write (6, 3006)
  endif
  write (6, 3005) Current/1.d6
  !.....................................................................


  ! 3. set up equilibrium ..............................................
  Bt_sign = int(sign(1.d0, Bcentr))
  Ip_sign = int(sign(1.d0, Current))

  ! 3.1. set up 1D arrays for R, Z, Psin values
  allocate (Rtmp(nR), Ztmp(nZ), Psintmp(nR))
  allocate (REQD(nR+nord), ZEQD(nZ+nord), PsinEQD(nR+nord))
  do i=1,nR
     Rtmp(i) = Rleft + (i-1.d0) / (nR-1.d0) * Rdim
  enddo
  call dbsnak (nR, Rtmp, nord, REQD)
  do i=1,nZ
     Ztmp(i) = Zmid - Zdim/2.d0 + (i-1.d0) / (nZ-1.d0) * Zdim
  enddo
  call dbsnak (nZ, Ztmp, nord, ZEQD)
  do i=1,nR
     Psintmp(i) = (i-1.d0) / (nR-1.d0)
  enddo
  call dbsnak (nR, Psintmp, nord, PsinEQD)

  ! 3.2. set up arrays for dependent coefficients
  allocate (fpolcoeff(nR), prescoeff(nR), Psicoeff(nR, nZ))
  allocate (ffprimcoeff(nR), pprimecoeff(nR))

  ! setup spline interpolation for Psi
  call dbs2in(nR,Rtmp,nZ,Ztmp,psirz,nR,nord,nord,REQD,ZEQD,Psicoeff)

  ! setup spline interpolation for Fpol
  call dbsint(nR,Psintmp,fpol,nord,PsinEQD,fpolcoeff)
  call dbsint(nR,Psintmp,pres,nord,PsinEQD,prescoeff)
  call dbsint(nR,Psintmp,ffprim,nord,PsinEQD,ffprimcoeff)
  call dbsint(nR,Psintmp,pprime,nord,PsinEQD,pprimecoeff)

  ! 3.3. cleanup
  deallocate (Rtmp, Ztmp, Psintmp)
  !.....................................................................


  ! 4. set output variables ............................................
  Psi_axis = Simag
  Psi_sepx = Sibry
  !.....................................................................

 2000 format (6a8,3i4)
 2001 format (6a8,a80)
 2020 format (5e16.9)
 2022 format (2i5)
 2024 format (i5,e16.9,i5)
 3000 format (8x,a,5a8)
 3001 format (8x,'Grid resolution: ',10x,i4,' x ',i4, ' nodes')
 3002 format (8x,'Computational Box:          ', &
              'R      = ',f8.3, ' cm  ->  ',f8.3,' cm',/36x, &
              'Z      = ',f8.3, ' cm  ->  ',f8.3,' cm')
 3003 format (8x,'Reference position:         R0     = ',f8.3,' cm')
 3004 format (8x,'Toroidal magnetic field:    Bt(R0) = ',f8.3,' T')
 3005 format (8x,'Total plasma current:       Ip     = ',f8.3,' MA')
 3006 format (8x,'Direction of plasma current is taken from flux', &
              ' distribution!')
  end subroutine geqdsk_load
!===============================================================================



!===============================================================================
! Broadcast data for parallel execution
!===============================================================================
  subroutine geqdsk_broadcast()
  use parallel

  call broadcast_inte_s (nR)
  call broadcast_inte_s (nZ)
  call broadcast_inte_s (nord)
  call broadcast_real_s (Simag)
  call broadcast_real_s (Sibry)
  call broadcast_real_s (Rdim)
  call broadcast_real_s (Zdim)
  call broadcast_real_s (Rcentr)
  call broadcast_real_s (Rleft)
  call broadcast_real_s (Zmid)
  call broadcast_real_s (Rmaxis)
  call broadcast_real_s (Zmaxis)

  if (mype .gt. 0) then
      allocate (REQD(nR+nord), ZEQD(nZ+nord), PsinEQD(nR+nord))
      allocate (fpolcoeff(nR), prescoeff(nR), Psicoeff(nR, nZ))
      allocate (ffprimcoeff(nR), pprimecoeff(nR))
  endif

  call broadcast_real   (REQD, nR+nord)
  call broadcast_real   (ZEQD, nZ+nord)
  call broadcast_real   (PsinEQD, nR+nord)
  call broadcast_real   (fpolcoeff, nR)
  call broadcast_real   (prescoeff, nR)
  call broadcast_real   (ffprimcoeff, nR)
  call broadcast_real   (pprimecoeff, nR)
  call broadcast_real   (Psicoeff, nR*nZ)

  ! limiter/wall configuration
  call broadcast_inte_s (limitr)
  if (mype .gt. 0) then
      allocate (rlim(limitr), zlim(limitr))
  endif
  call broadcast_real   (rlim, limitr)
  call broadcast_real   (zlim, limitr)

  end subroutine geqdsk_broadcast
!===============================================================================



!===============================================================================
! Calculate R,phi,Z components of magnetic field vector [Gauss] at r=(R,Z [cm], phi [rad])
!===============================================================================
  function geqdsk_get_Bf(r) result(Bf)
  use bspline
  real(real64), intent(in) :: r(3)
  real(real64)             :: Bf(3)


  real(real64) :: rr, zz, psi, dpsidr, dpsidz


  ! convert cm -> m
  rr      =  r(1) / 100.d0
  zz      =  r(2) / 100.d0


  ! force limits of computational box
  if (rr.lt.Rleft)          rr = Rleft
  if (rr.gt.Rleft+Rdim)     rr = Rleft+Rdim
  if (zz.lt.Zmid-Zdim/2.d0) zz = Zmid-Zdim/2.d0
  if (zz.gt.Zmid+Zdim/2.d0) zz = Zmid+Zdim/2.d0


  ! poloidal field
  psi     =  dbs2dr(0,0,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
  psi     =  (psi - Simag) / (Sibry - Simag)
  dpsidr  =  dbs2dr(1,0,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
  dpsidz  =  dbs2dr(0,1,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
  Bf(1)   = -dpsidz/rr *1.d4 ! T -> Gauss
  Bf(2)   =  dpsidr/rr *1.d4 ! T -> Gauss
  Bf(3)   =  0.d0


  ! toroidal field
  if (psi.gt.1.d0) psi = 1.d0
  Bf(3)   =  dbsval(psi,nord,PsinEQD,nR,fpolcoeff) / rr
  Bf(3)   =  Bf(3) * 1.d4    ! T -> Gauss

  end function geqdsk_get_Bf
!===============================================================================



!===============================================================================
! Calculate Jacobian(Bf) [T/m] in cylindrical coordinates
!===============================================================================
  function geqdsk_get_JBf(r) result(J)
  use bspline
  real(real64), intent(in) :: r(3)
  real(real64)             :: J(3,3)

  real(real64) :: rr, zz, psi, dpsidr, dpsidz, d2psi(3), F, dFdpsi


  ! convert cm -> m
  rr      =  r(1) / 100.d0
  zz      =  r(2) / 100.d0


  ! force limits of computational box
  if (rr.lt.Rleft)          rr = Rleft
  if (rr.gt.Rleft+Rdim)     rr = Rleft+Rdim
  if (zz.lt.Zmid-Zdim/2.d0) zz = Zmid-Zdim/2.d0
  if (zz.gt.Zmid+Zdim/2.d0) zz = Zmid+Zdim/2.d0


  J       =  0.d0
  psi     =  dbs2dr(0,0,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
  psi     =  (psi - Simag) / (Sibry - Simag)
  if (psi.gt.1.d0) psi = 1.d0
  dpsidr  =  dbs2dr(1,0,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
  dpsidz  =  dbs2dr(0,1,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
  d2psi(1)=  dbs2dr(2,0,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
  d2psi(2)=  dbs2dr(1,1,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
  d2psi(3)=  dbs2dr(0,2,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
  F       =  dbsval(psi,nord,PsinEQD,nR,fpolcoeff)
  dFdpsi  =  dbsder(1,psi,nord,PsinEQD,nR,fpolcoeff)

  ! Br
  J(1,1)  =  -d2psi(2)/rr + dpsidz/rr**2
  J(1,2)  =  -d2psi(3)/rr
  ! Bz
  J(2,1)  =  d2psi(1)/rr - dpsidr/rr**2
  J(2,2)  =  d2psi(2)/rr
  ! Bphi
  J(3,1)  =  dFdpsi * dpsidr / rr - F/rr**2
  J(3,2)  =  dFdpsi * dpsidz / rr

  end function geqdsk_get_JBf
!===============================================================================



!===============================================================================
! Return poloidal magnetic flux at r=(R,Z [cm], phi [rad])
!===============================================================================
  function geqdsk_get_Psi(r) result(Psi)
  use bspline
  real(real64), intent(in) :: r(3)
  real(real64)             :: Psi

  real(real64) :: rr, zz


  ! convert cm -> m
  rr  = r(1) / 100.d0
  zz  = r(2) / 100.d0


  ! force limits of computational box
  if (rr.lt.Rleft)          rr = Rleft
  if (rr.gt.Rleft+Rdim)     rr = Rleft+Rdim
  if (zz.lt.Zmid-Zdim/2.d0) zz = Zmid-Zdim/2.d0
  if (zz.gt.Zmid+Zdim/2.d0) zz = Zmid+Zdim/2.d0


  Psi =  dbs2dr(0,0,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)

  end function geqdsk_get_Psi
!===============================================================================



!===============================================================================
! Return (mR,mZ)-th derivative of poloidal magnetic flux at r=(R,Z [cm])
!===============================================================================
  function geqdsk_get_DPsi(r, mR, mZ) result(DPsi)
  use bspline
  real(real64), intent(in) :: r(2)
  integer,      intent(in) :: mR, mZ
  real(real64)             :: DPsi

  real(real64) :: rr, zz


  ! convert cm -> m
  rr   = r(1) / 100.d0
  zz   = r(2) / 100.d0


  ! force limits of computational box
  if (rr.lt.Rleft)          rr = Rleft
  if (rr.gt.Rleft+Rdim)     rr = Rleft+Rdim
  if (zz.lt.Zmid-Zdim/2.d0) zz = Zmid-Zdim/2.d0
  if (zz.gt.Zmid+Zdim/2.d0) zz = Zmid+Zdim/2.d0


  DPsi = dbs2dr(mR,mZ,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
  DPsi = DPsi / 100.d0**(mR+mZ)

  end function geqdsk_get_DPsi
!===============================================================================



!===============================================================================
! Return Pressure at radial position Psi
! NOT IMPLEMENTED YET
!===============================================================================
  function geqdsk_get_pressure(Psi) result(P)
  real(real64), intent(in) :: Psi
  real(real64)             :: P

  real(real64) :: PsiN


  PsiN =  (psi - Simag) / (Sibry - Simag)
  P    =  dbsval(PsiN,nord,PsinEQD,nR,prescoeff)

  end function geqdsk_get_pressure
!===============================================================================



!===============================================================================
! Return boundaries [cm] of equilibrium domain
!===============================================================================
  subroutine geqdsk_get_domain(Rbox, Zbox)
  real(real64), intent(out) :: Rbox(2), Zbox(2)


  Rbox(1) = Rleft * 1.d2
  Rbox(2) = (Rleft+Rdim) * 1.d2
  Zbox(1) = (Zmid-Zdim/2.d0) * 1.d2
  Zbox(2) = (Zmid+Zdim/2.d0) * 1.d2

  end subroutine geqdsk_get_domain
!===============================================================================



!===============================================================================
! Return the wall configuration provided in g-file
!===============================================================================
  subroutine geqdsk_export_boundary(S)
  use curve2D
  type(t_curve), intent(out) :: S


  call make_2D_curve (limitr, rlim, zlim, S)
  ! m -> cm
  S%x = S%x * 1.d2

  end subroutine geqdsk_export_boundary
!===============================================================================



!===============================================================================
! Write characteristic profiles
!===============================================================================
  subroutine geqdsk_info()
  use math

  real(real64) :: rr, zz, psi, psi0, ffprim0, pprime0, jt
  integer      :: i


  ! write safety factor and pressure profiles
  open  (99, file='q.dat')
  open  (98, file='P.dat')
  do i=1,nR
      write (99, *) (i-1.d0) / (nR-1.d0), qpsi(i)
      write (98, *) (i-1.d0) / (nR-1.d0), pres(i)
  enddo
  close (99)
  close (98)


  ! write current profile
  open  (99, file='jt.txt')
  do i=1,nR
     rr      =  Rleft + (i-1.d0) / (nR-1.d0) * Rdim
     zz      =  0.d0
     psi     =  dbs2dr(0,0,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
     psi     =  (psi - Simag) / (Sibry - Simag)
     psi0    =  psi
     if (psi.gt.1.d0) psi0 = 1.d0

     ffprim0 =  dbsval(psi0,nord,PsinEQD,nR,ffprimcoeff)
     pprime0 =  dbsval(psi0,nord,PsinEQD,nR,pprimecoeff)

     jt      = rr * pprime0 + ffprim0 / rr / pi / 4.d-7

     write (99, 1000) rr, psi, ffprim0, pprime0, jt
  enddo
  close (99)

 1000 format (5e16.8)
  end subroutine geqdsk_info
!===============================================================================

  end module geqdsk
