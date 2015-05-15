!===============================================================================
! Magnetic field configuration provided by the DIVA MHD code
!===============================================================================
module divamhd
  implicit none

  private

  ! reference toroidal magnetic field [T] and radial position [cm]
  real*8  :: Bt, R0

  ! grid resolution
  integer :: nrark, nzark

  ! dummy grid resolution
  integer, parameter :: arkrres  = 530
  integer, parameter :: arkzres  = 530
  integer, parameter :: arkorder =   5

  real*8 :: &
      arkr    (arkrres), &
      arkz    (arkzres), &
      arkrnot (arkrres+arkorder), &
      arkznot (arkzres+arkorder), &
      arkpsi  (arkrres,arkzres), &
      dummy   (arkrres,arkzres), &
      arkcoef (arkrres,arkzres), &
      arkrl, arkzl, arkru, arkzu


  public :: &
     divamhd_load, &
     divamhd_broadcast, &
     divamhd_get_Bf, &
     divamhd_get_Psi, &
     divamhd_get_DPsi, &
     divamhd_get_domain

  contains
!===============================================================================



!===============================================================================
  subroutine divamhd_load (Data_File, Ip_scale, Bt_, R0_)
  use bspline

  character*120, intent(in) :: Data_File
  real*8, intent(in)        :: Ip_scale, Bt_, R0_

  integer, parameter :: iu = 17

  character*80 :: v80
  character*1  :: dummyChar
  real*8       :: sumipl, arkdr, arkdz
  integer      :: i, j, io


  ! set reference toroidal magnetic field
  Bt = Bt_
  R0 = R0_ / 1.d2

  ! read data from DIVA-equilibrium file
  open  (iu, file=Data_file, iostat=io)
  if (io.ne.0) then
     write (6,2000) Data_file
     stop
  endif
  ! skip header
  read (iu,*)
  read (iu,*)
  read (iu,*)
  read (iu,*)
  read (iu,5050) v80

  ! read grid resolution
  read (v80(32:32),5011) dummyChar
  if (dummyChar.eq.'x') then
     ! 2 digits
     read (v80(29:31),5010) nrark
     read (v80(33:35),5010) nzark
  else
     read (v80(33:33),5011) dummyChar
     if (dummyChar.eq.'x') then
        ! 3 digits
        read (v80(30:32),5010) nrark
        read (v80(34:36),5010) nzark
     else
        write (6, *) 'error: cannot determine grid resolution for DIVA MHD data!'
        stop
     endif
  endif
  write(6,1000) nrark, nzark

  ! read physical position
  read  (iu,        5050) v80
  read  (v80(38:42),5025) arkrl
  read  (v80(44:50),5020) arkzl
  read  (iu,        5050) v80
  read  (v80(38:42),5025) arkru
  read  (v80(44:50),5020) arkzu
  write (6,         1001) arkrl, arkzl, arkru, arkzu
  read  (iu, *)
  read  (iu, *)

  ! read main data
  read  (iu, 5040)((arkpsi(i,j),i=1,nrark+1),j=1,nzark+1)	! psi
  read  (iu, *)
  read  (iu, 5040)((dummy(i,j),i=1,nrark+1),j=1,nzark+1)	! B_R
  read  (iu, *)
  read  (iu, 5040)((dummy(i,j),i=1,nrark+1),j=1,nzark+1)	! B_Z
  read  (iu, *)
  read  (iu, 5040)((dummy(i,j),i=1,nrark+1),j=1,nzark+1)	! j
  close (iu)


  ! calculate total plasma current
  sumipl  = 0.d0
  do 100 i = 1,nrark+1
     do 100 j = 1,nzark+1
        sumipl = sumipl+dummy(i,j)
  100 continue


  ! setup interpolation grid
  arkdr   = (arkru-arkrl)/nrark
  arkdz   = (arkzu-arkzl)/nzark
  do 110 i = 1,nrark+1
     arkr(i)  = arkrl+(i-1)*arkdr
  110 continue
  do 120 i = 1,nzark+1
     arkz(i)  = arkzl+(i-1)*arkdz
  120 continue


  ! scale psi-function for given I_P
  if (Ip_scale.ne.0.d0) then
     sumipl  = sumipl*arkdr*arkdz
     arkpsi = arkpsi / sumipl * Ip_scale
  endif


  call dbsnak(nrark+1,arkr,arkorder,arkrnot)
  call dbsnak(nzark+1,arkz,arkorder,arkznot)
  call dbs2in(nrark+1,arkr,nzark+1,arkz,arkpsi,arkrres,arkorder, &
              arkorder,arkrnot,arkznot,arkcoef)

 1000 format (8x,'grid resolution for psi-function: ', &
              19x,'R : ',i4,' points;    Z : ',i4,' points')
 1001 format (8x,'physical resolution:', &
              18x,' left lower point: (',f7.3,', ',f7.3,')'/ &
              46x,'right upper point: (',f7.3,', ',f7.3,')')
 2000 format ('error reading DIVA-file: ',a120)
 5010 format (i3)
 5011 format (a1)
 5025 format (f5.3)
 5020 format (f7.3)
 5040 format (2x,5e14.6)
 5050 format (a80)
  end subroutine divamhd_load
!===============================================================================



!===============================================================================
  subroutine divamhd_broadcast
  use parallel

  call broadcast_real_s (Bt)
  call broadcast_real_s (R0)

  call broadcast_inte_s (nrark)
  call broadcast_inte_s (nzark)

  call broadcast_real_s (arkrl)
  call broadcast_real_s (arkru)
  call broadcast_real_s (arkzl)
  call broadcast_real_s (arkzu)
  call broadcast_real   (arkr,    arkrres)
  call broadcast_real   (arkz,    arkzres)
  call broadcast_real   (arkrnot, arkrres+arkorder)
  call broadcast_real   (arkznot, arkzres+arkorder)
  call broadcast_real   (arkpsi,  arkrres*arkzres)
  call broadcast_real   (arkcoef, arkrres*arkzres)

  end subroutine divamhd_broadcast
!===============================================================================



!===============================================================================
! Calculate R,phi,Z components of magnetic field vector [Gauss] at r=(R,Z [cm], phi [rad])
!===============================================================================
  function divamhd_get_Bf(r) result(Bf)
  use bspline
  use math

  real*8, intent(in)  :: r(3)
  real*8              :: Bf(3)

  real*8   :: rr, zz, tpr, bbr, dpsidr, dpsidz


  ! convert cm -> m
  rr     = r(1) / 100.d0
  zz     = r(2) / 100.d0
  tpr    = pi2*rr

  ! force limits of computational box
  if (rr.lt.arkrl) rr = arkrl
  if (rr.gt.arkru) rr = arkru
  if (zz.lt.arkzl) rr = arkzl
  if (zz.gt.arkzu) rr = arkzu


  dpsidr = dbs2dr(1,0,rr,zz,arkorder,arkorder,arkrnot,arkznot,nrark+1,nzark+1,arkcoef)
  dpsidz = dbs2dr(0,1,rr,zz,arkorder,arkorder,arkrnot,arkznot,nrark+1,nzark+1,arkcoef)

  Bf(1) = -dpsidz/tpr*1.d4
  Bf(2) =  dpsidr/tpr*1.d4
  Bf(3) =  Bt*R0/r(1)*1.d4

  end function divamhd_get_Bf
!===============================================================================



!===============================================================================
! Sample (derivative of) poloidal magnetic flux at r=(R,Z [cm], phi [rad])
!===============================================================================
  function divamhd_get_Psi(r) result(Psi)
  real*8, intent(in)  :: r(3)
  real*8              :: Psi

  Psi = divamhd_get_DPsi(r(1:2), 0, 0)

  end function divamhd_get_Psi
!===============================================================================
  function divamhd_get_DPsi(r, mR, mZ) result(DPsi)
  use bspline

  real*8, intent(in)  :: r(2)
  integer, intent(in) :: mR, mZ
  real*8              :: DPsi

  real*8 :: rr, zz


  ! convert cm -> m
  rr   = r(1) / 100.d0
  zz   = r(2) / 100.d0
  DPsi = dbs2dr(0,0,rr,zz,arkorder,arkorder,arkrnot,arkznot,nrark+1,nzark+1,arkcoef)
  DPsi = DPsi / 100.d0**(mR+mZ)

  end function divamhd_get_DPsi
!===============================================================================



!===============================================================================
! Return boundaries [cm] of equilibrium domain
!===============================================================================
  subroutine divamhd_get_domain(Rbox, Zbox)
  real*8, intent(out) :: Rbox(2), Zbox(2)


  Rbox(1) = arkrl * 1.d2
  Rbox(2) = arkru * 1.d2
  Zbox(1) = arkzl * 1.d2
  Zbox(2) = arkzu * 1.d2

  end subroutine divamhd_get_domain
!===============================================================================

end module divamhd
