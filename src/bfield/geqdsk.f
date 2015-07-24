!===============================================================================
! Magnetic field configuration given in G-EQDSK format
!
! THIS MODULE NEED SOME CLEANUP !!!
!===============================================================================
      module geqdsk
      use iso_fortran_env
      use magnetic_axis
      implicit none

      private 

      character*120 ::
     .    G_file = ''
      logical       ::
     .    use_wall = .true.,		! Add provided boundary data to wall configuration.
     .    CurrentFix = .true.		! Use sign(Current) to determine the plasma current direction,
					! and hence the poloidal field direction (based on TRIP3D).
      integer       ::
     .    DiagnosticLevel = 0

      namelist /G_EQDSK_Input/
     .         G_file, use_wall, CurrentFix, DiagnosticLevel

      ! set spline interpolation order
      integer :: nord = 5


      real*8, dimension(:,:), allocatable :: psirz, Psicoeff
      real*8, dimension(:), allocatable :: fpol, pres, ffprim, pprime,
     .    qpsi, fpolcoeff, ffprimcoeff, pprimecoeff,
     .    REQD, ZEQD, PsinEQD, rbbbs, zbbbs
      real*8, dimension(:), allocatable, target :: rlim, zlim

      real*8  :: Rdim, Zdim, Rcentr, Rleft, Zmid, Current,
     .    Rmaxis, Zmaxis, Simag, Sibry, Bcentr

      integer :: nR, nZ, nbbbs, limitr

      logical :: switch_poloidal_field

      public ::
     1    geqdsk_load,
     2    geqdsk_broadcast,
     3    geqdsk_get_Bf,
     4    geqdsk_get_JBf,
     5    geqdsk_get_Psi, 
     6    geqdsk_get_DPsi,
     7    geqdsk_get_domain,
     8    geqdsk_provides_boundary,
     9    geqdsk_export_boundary
!     a    sample_psi1_EQDSK, get_D1Psi_geqdsk,
!     b    sample_psi2_EQDSK, get_D2Psi_geqdsk,
!     4    equi_info_EQDSK, get_equi_domain_EQDSK, magnetic_axis_geqdsk,
!     b    switch_poloidal_field_only_EQDSK, use_wall,
!     5    guess_Xpoint_EQDSK,
!     b    geqdsk_provides_PFC, export_PFC_geqdsk,get_wall_configuration,
!     c    jt_profile,
!     6    Bcentr, Current
      contains
!===============================================================================



!===============================================================================
! Load equilibrum data from file
!===============================================================================
      subroutine geqdsk_load (Data_File, use_PFC_, CurrentFix_, DL_,
     .                        psi_axis, psi_sepx)
      use run_control, only: Spline_Order
      use bspline
      character*120, intent(in) :: Data_File
      logical, intent(in), optional  :: use_PFC_, CurrentFix_
      integer, intent(in), optional  :: DL_
      real*8, intent(out), optional  :: psi_axis,psi_sepx

      integer, parameter :: iu = 17

      real*8, dimension(:), allocatable :: Rtmp, Ztmp, Psintmp

      character*8  :: case_(6)
      integer :: i, j, idum
      real*8  :: xdum, PsiC, PsiL, PsiR, rCurrentFix
      logical :: lL, lR, concaveup


      if (present(use_PFC_)) use_wall = use_PFC_
      if (present(CurrentFix_)) CurrentFix = CurrentFix_
      if (present(DL_)) DiagnosticLevel = DL_
      nord = Spline_Order

      open  (iu, file=Data_File)

      ! read equilibrium configuration
      read  (iu, 2000) (case_(i),i=1,6), idum, nR, nZ
      read  (iu, 2020) Rdim, Zdim, Rcentr, Rleft, Zmid
      read  (iu, 2020) Rmaxis, Zmaxis, Simag, Sibry, Bcentr
      read  (iu, 2020) Current, Simag, xdum, Rmaxis, xdum
      read  (iu, 2020) Zmaxis, xdum, Sibry, xdum, xdum
      call setup_magnetic_axis_2D (Rmaxis*1.d2, Zmaxis*1.d2)


      ! runtime feedback of characteristic values
      write (6, 3000) adjustl(case_(1)), case_(2:6)
      write (6, *)
      write (6, 3001) nR, nZ
      write (6, 3002) Rleft*1.e2, (Rleft+Rdim)*1.d2,
     .                (Zmid-Zdim/2.d0)*1.d2, (Zmid+Zdim/2.d0)*1.d2
      write (6, 3003) Rcentr*1.d2
      write (6, 3004) Bcentr


      allocate (fpol(nR), pres(nR), ffprim(nR), pprime(nR), qpsi(nR))
      allocate (psirz(nR,nZ))

      ! read equilibrium data
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


! diagnostic output ............................................................
      if (DiagnosticLevel.eq.1) then
         open  (99, file='boundary_test.txt')
         do i=1,nbbbs
            write (99, *) rbbbs(i), zbbbs(i)
         enddo
         close (99)
         open  (99, file='wall_test.txt')
         do i=1,limitr
            write (99, *) rlim(i), zlim(i)
         enddo
         close (99)
      endif
      if (DiagnosticLevel.eq.2) then
         open  (99, file='qtest.txt')
         do i=1,nR
             write (99, *) (i-1.d0) / (nR-1.d0), qpsi(i)
         enddo
         close (99)
      endif
!...............................................................................


c-----------------------------------------------------------------------
! Use sign(Current) to determine the plasma current direction, and hence
! the poloidal field direction (based on TRIP3D).
c-----------------------------------------------------------------------
      !if (CurrentFix) then
      rCurrentFix = 1.d0
c- Make psirz concave up, i.e. psi_center > psi_edge
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
            psirz = - psirz
            Simag = - Simag
            Sibry = - Sibry
         else
            rCurrentFix = -1.d0
         endif
      endif
c now psirz is concave up
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
c- Current fix: (+) -> flip to concave down, (-) -> keep concave up
      if (Current .gt. 0.d0) then
         if (CurrentFix) then
            psirz = - psirz
            Simag = - Simag
            Sibry = - Sibry
         else
            rCurrentFix = -1.d0 * rCurrentFix
         endif
      endif
      !endif

      if (.not.CurrentFix) Current = Current * rCurrentFix
      write (6, 3005) Current/1.d6
      if (.not.CurrentFix) then
         write (6, 3006)
      endif
c-----------------------------------------------------------------------
      Bt_sign = int(sign(1.d0, Bcentr))
      Ip_sign = int(sign(1.d0, Current))


      ! setup 1D arrays of R, Z, Psin values
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

      allocate (fpolcoeff(nR), Psicoeff(nR, nZ))
      allocate (ffprimcoeff(nR), pprimecoeff(nR))

      ! setup spline interpolation for Psi
      call dbs2in(nR,Rtmp,nZ,Ztmp,psirz,nR,nord,nord,REQD,ZEQD,Psicoeff)

      ! setup spline interpolation for Fpol
      call dbsint(nR,Psintmp,fpol,nord,PsinEQD,fpolcoeff)
      call dbsint(nR,Psintmp,ffprim,nord,PsinEQD,ffprimcoeff)
      call dbsint(nR,Psintmp,pprime,nord,PsinEQD,pprimecoeff)

      deallocate (Rtmp, Ztmp, Psintmp)

      if (present(psi_axis)) psi_axis = Simag
      if (present(psi_sepx)) psi_sepx = Sibry
      return
 2000 format (6a8,3i4)
 2020 format (5e16.9)
 2022 format (2i5)
 2024 format (i5,e16.9,i5)

 3000 format (8x,a,5a8)
 3001 format (8x,'Grid resolution: ',10x,i4,' x ',i4, ' nodes')
 3002 format (8x,'Computational Box:          ',
     .        'R      = ',f8.3, ' cm  ->  ',f8.3,' cm',/
     .        8x,'                            ',
     .        'Z      = ',f8.3, ' cm  ->  ',f8.3,' cm')
 3003 format (8x,'Reference position:         R0     = ',f8.3,' cm')
 3004 format (8x,'Toroidal magnetic field:    Bt(R0) = ',f8.3,' T')
 3005 format (8x,'Total plasma current:       Ip     = ',f8.3,' MA')
 3006 format (8x,'Direction of plasma current is taken from flux',
     .        ' distribution!')
      end subroutine geqdsk_load
!===============================================================================



!===============================================================================
! Broadcast data for parallel execution
!===============================================================================
      subroutine geqdsk_broadcast
      use parallel

      call broadcast_inte_s (nR)
      call broadcast_inte_s (nZ)
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
          allocate (fpolcoeff(nR), Psicoeff(nR, nZ))
          allocate (ffprimcoeff(nR), pprimecoeff(nR))
      endif

      call broadcast_real   (REQD, nR+nord)
      call broadcast_real   (ZEQD, nZ+nord)
      call broadcast_real   (PsinEQD, nR+nord)
      call broadcast_real   (fpolcoeff, nR)
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

      return
      end subroutine geqdsk_broadcast
!===============================================================================



!===============================================================================
! Calculate R,phi,Z components of magnetic field vector [Gauss] at r=(R,Z [cm], phi [rad])
!===============================================================================
      function geqdsk_get_Bf(r) result(Bf)
      use bspline

      real*8, intent(in)  :: r(3)
      real*8              :: Bf(3)


      real*8 :: rr, zz, psi, dpsidr, dpsidz

      ! convert cm -> m
      rr      =  r(1) / 100.d0
      zz      =  r(2) / 100.d0

      ! force limits of computational box
      if (rr.lt.Rleft)          rr = Rleft
      if (rr.gt.Rleft+Rdim)     rr = Rleft+Rdim
      if (zz.lt.Zmid-Zdim/2.d0) zz = Zmid-Zdim/2.d0
      if (zz.gt.Zmid+Zdim/2.d0) zz = Zmid+Zdim/2.d0


      psi     =  dbs2dr(0,0,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
      psi     =  (psi - Simag) / (Sibry - Simag)
      dpsidr  =  dbs2dr(1,0,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
      dpsidz  =  dbs2dr(0,1,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
      Bf(1)   = -dpsidz/rr *1.d4 ! T -> Gauss
      Bf(2)   =  dpsidr/rr *1.d4 ! T -> Gauss
      Bf(3)   =  0.d0
      if (switch_poloidal_field) return


      ! toroidal field
      if (psi.gt.1.d0) psi = 1.d0
      Bf(3)   =  dbsval(psi,nord,PsinEQD,nR,fpolcoeff) / rr
      Bf(3)   =  Bf(3) * 1.d4 ! T -> Gauss

      end function geqdsk_get_Bf
!===============================================================================



!===============================================================================
      function geqdsk_get_JBf(r) result(J)
      use bspline

      real(real64), intent(in)  :: r(3)
      real(real64)              :: J(3,3)

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


      J       =  J * 1.d4 ! T -> Gauss

      end function geqdsk_get_JBf
!===============================================================================



!===============================================================================
! Sample poloidal magnetic flux at r=(R,Z [cm], phi [rad])
!===============================================================================
      function geqdsk_get_Psi(r) result(Psi)
      use bspline

      real*8, intent(in) :: r(3)
      real*8             :: Psi

      real*8 :: rr, zz


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
! Sample (mR,mZ)-th derivative of poloidal magnetic flux at r=(R,Z [cm])
!===============================================================================
      function geqdsk_get_DPsi(r, mR, mZ) result(DPsi)
      use bspline

      real*8, intent(in)  :: r(2)
      integer, intent(in) :: mR, mZ
      real*8              :: DPsi

      real*8 :: rr, zz


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
! Return boundaries [cm] of equilibrium domain
!===============================================================================
      subroutine geqdsk_get_domain(Rbox, Zbox)
      real*8, intent(out) :: Rbox(2), Zbox(2)


      Rbox(1) = Rleft * 1.d2
      Rbox(2) = (Rleft+Rdim) * 1.d2
      Zbox(1) = (Zmid-Zdim/2.d0) * 1.d2
      Zbox(2) = (Zmid+Zdim/2.d0) * 1.d2

      return
      end subroutine geqdsk_get_domain
!===============================================================================



!===============================================================================
c Return the wall configuration provided in the g-file
!===============================================================================
      function geqdsk_provides_boundary() result(l)
      logical :: l

      l = use_wall
      end function geqdsk_provides_boundary
c-------------------------------------------------------------------------------
      subroutine geqdsk_export_boundary(S)
      use curve2D
      type(t_curve), intent(out) :: S

      call make_2D_curve (limitr, rlim, zlim, S)
      ! m -> cm
      S%x = S%x * 1.d2
      end subroutine geqdsk_export_boundary
!===============================================================================
















c-------------------------------------------------------------------------------
c- sample first derivative of poloidal flux at R, Z
c-------------------------------------------------------------------------------
      subroutine sample_psi1_EQDSK (R, Z, dpsidr, dpsidz)
      use bspline

      real*8, intent(in)  :: R, z
      real*8, intent(out) :: dpsidr, dpsidz

      real*8 :: rr, zz

      ! convert cm -> m
      rr   = R / 100.d0
      zz   = Z / 100.d0

      dpsidr  =  dbs2dr(1,0,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
      dpsidz  =  dbs2dr(0,1,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)


      return
      end subroutine sample_psi1_EQDSK
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c- sample second derivative of poloidal flux at R, Z
c-------------------------------------------------------------------------------
      subroutine sample_psi2_EQDSK (R, Z, dpsi2)
      use bspline

      real*8, intent(in)  :: R, z
      real*8, intent(out) :: dpsi2(3)

      real*8 :: rr, zz

      ! convert cm -> m
      rr   = R / 100.d0
      zz   = Z / 100.d0

      dpsi2(1)  =  dbs2dr(2,0,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
      dpsi2(2)  =  dbs2dr(1,1,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
      dpsi2(3)  =  dbs2dr(0,2,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)


      return
      end subroutine sample_psi2_EQDSK
c-------------------------------------------------------------------------------
      function get_D2Psi_geqdsk(r) result(D2Psi)
      use bspline

      real*8, intent(in) :: r(3)
      real*8             :: D2Psi(3)

      real*8 :: rr, zz


      ! convert cm -> m
      rr  = r(1) / 100.d0
      zz  = r(2) / 100.d0
      D2Psi(1)  =  dbs2dr(2,0,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
      D2Psi(2)  =  dbs2dr(1,1,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)
      D2Psi(3)  =  dbs2dr(0,2,rr,zz,nord,nord,REQD,ZEQD,nR,nZ,Psicoeff)

      end function get_D2Psi_geqdsk
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
      subroutine equi_info_EQDSK (psi_sepx, psi_axis, R_axis, Z_axis)
      real*8, intent(out) :: psi_sepx, psi_axis, R_axis, Z_axis

      ! convert m -> cm
      R_axis = Rmaxis * 1.d2
      Z_axis = Zmaxis * 1.d2
      psi_sepx = Sibry
      psi_axis = Simag

      return
      end subroutine equi_info_EQDSK
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
      function magnetic_axis_geqdsk() result(r)
      real*8 :: r(2)

      ! convert m -> cm
      r(1) = Rmaxis * 1.d2
      r(2) = Zmaxis * 1.d2

      end function magnetic_axis_geqdsk
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
c Switch between poloidal magnetic field component and full configuration
c-------------------------------------------------------------------------------
      subroutine switch_poloidal_field_only_EQDSK (lswitch)
      logical, intent(in) :: lswitch

      switch_poloidal_field = lswitch

      return
      end subroutine switch_poloidal_field_only_EQDSK
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
      subroutine guess_Xpoint_EQDSK (Rx, Zx)
      real*8, intent(out) :: Rx, Zx

      Rx = 0.5d0 * (Rleft + Rcentr) * 1.d2
      Zx = (Zmid - Zdim/3.d0) * 1.d2

      return
      end subroutine guess_Xpoint_EQDSK
c-------------------------------------------------------------------------------


c-------------------------------------------------------------------------------
      subroutine jt_profile
      use bspline
      use math

      integer :: i
      real*8  :: rr, zz, psi, psi0, ffprim0, pprime0, jt

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

      return
 1000 format (5e16.8)
      end subroutine jt_profile
c-------------------------------------------------------------------------------

      end module geqdsk
