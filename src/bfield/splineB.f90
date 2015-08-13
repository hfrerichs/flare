!===============================================================================
! Use spline interpolation on vector potential or magnetic field vector given on
! a cylindrical grid.
! The spline order can be given by the parameter "spline_order" (default = 5).
!===============================================================================
module splineB
  use iso_fortran_env
  implicit none
  private

  integer, parameter :: n_max = 32

  integer, parameter :: &
     BFIELD = 1, &
     VECTOR_POTENTIAL = 2


  integer            :: n_sets = 0, spline_order = 5, data_type = VECTOR_POTENTIAL
  character(len=120) :: &
     grid_layout     = 'layout.dat', &
     basename(n_max) = '', &
     R_suffix        = '_r.dat', &
     Z_suffix        = '_z.dat', &
     Phi_suffix      = '_phi.dat'

  real(real64)       :: amplitude(n_max) = 1.d0, amplitude0 = 1.d0

  namelist /SplineB_Input/ &
     n_sets, grid_layout, basename, amplitude, amplitude0, data_type, spline_order, &
     R_suffix, Z_suffix, Phi_suffix


  ! internal variables
  !integer, parameter :: nord = 5

  real(real64), dimension(:,:,:), allocatable :: Ar, Az, Aphi, Arcoeff, Azcoeff, Aphicoeff
  real(real64), dimension(:),     allocatable :: Rnot, Znot, Phinot

  real(real64) :: R_min, R_max, Z_min, Z_max
  integer      :: nr, nz, nphi, nsym, nord


  public :: &
     splineB_load, &
     splineB_broadcast, &
     splineB_get_Bf, &
     splineB_get_JBf

  contains
!=======================================================================


!=======================================================================
  subroutine splineB_load (iu, iconfig)
  use run_control, only: Prefix
  use math
  use bspline
  integer, intent(in)  :: iu
  integer, intent(out) :: iconfig

  integer, parameter   :: iu1 = 70

  real(real64), dimension(:,:,:), allocatable :: tmp
  real(real64), dimension(:),     allocatable :: Rtmp, Ztmp, Phitmp
  character(len=120) :: Data_File
  integer :: i, j, k, is, io


  ! read user configuration
  rewind (iu)
  read   (iu, SplineB_Input, end=1000)
  iconfig = 1
  write (6, *)
  write (6, 1001)
  select case(data_type)
  case(VECTOR_POTENTIAL)
     write (6, 1002)
  case(BFIELD)
     write (6, 1003)
  end select
  write (6, 1004) spline_order
  nord = spline_order


  ! read grid layout
  Data_file = trim(Prefix)//grid_layout
  open  (iu1, file=Data_File)
  !read  (iu1, *) nr
  !read  (iu1, *) nz
  !read  (iu1, *) nphi
  !read  (iu1, *) R_min
  !read  (iu1, *) R_max
  !read  (iu1, *) Z_min
  !read  (iu1, *) Z_max
  !read  (iu1, *) nsym
  read  (iu1, *) nr, nz, nphi, nsym, R_min, R_max, Z_min, Z_max
  close (iu1)
  write (6, 1010) nr, nz, nphi
  write (6, 1011) R_min, Z_min
  write (6, 1012) R_max, Z_max
  write (6, 1013) nsym
  write (6, *)


  ! allocate memory for data
  nphi = nphi + 1
  allocate (Ar(nr, nz, nphi), Az(nr, nz, nphi), Aphi(nr, nz, nphi))
  allocate (tmp(nr, nz, nphi))


  ! read data
  Ar   = 0.d0
  Az   = 0.d0
  Aphi = 0.d0
  do is=1,n_sets
     write (6, 1020) basename(is)

     ! Ar
     write (6, 1021)
     Data_File = trim(Prefix)//trim(basename(is))//trim(R_suffix)
     open  (iu1, file=Data_File, iostat=io)
     if (io.ne.0) then
        write (6,2000) Data_File
        stop
     endif
     read  (iu1, *) (((tmp(i,j,k), k=1,nphi-1), j=1,nz), i=1,nr)
     Ar = Ar + amplitude(is)*tmp
     close (iu1)

     ! Az
     write (6, 1022)
     Data_File = trim(Prefix)//trim(basename(is))//trim(Z_suffix)
     open  (iu1, file=Data_File, iostat=io)
     if (io.ne.0) then
        write (6,2000) Data_File
        stop
     endif
     read  (iu1, *) (((tmp(i,j,k), k=1,nphi-1), j=1,nz), i=1,nr)
     Az = Az + amplitude(is)*tmp
     close (iu1)

     ! Aphi
     write (6, 1023)
     Data_File = trim(Prefix)//trim(basename(is))//trim(Phi_suffix)
     open  (iu1, file=Data_File, iostat=io)
     if (io.ne.0) then
        write (6,2000) Data_File
        stop
     endif
     read  (iu1, *) (((tmp(i,j,k), k=1,nphi-1), j=1,nz), i=1,nr)
     Aphi = Aphi + amplitude(is)*tmp
     close (iu1)

     write (6, *)
  enddo
  Ar   = Ar   * amplitude0
  Az   = Az   * amplitude0
  Aphi = Aphi * amplitude0
  deallocate (tmp)


  ! setup data
  ! set A(2*pi/nsym) = A(0)
  Ar(:,:,nphi)   = Ar(:,:,1)
  Az(:,:,nphi)   = Az(:,:,1)
  Aphi(:,:,nphi) = Aphi(:,:,1)

  ! prepare grid nodes
  allocate (Rtmp(nr), Ztmp(nz), Phitmp(nphi))
  allocate (Rnot(nr+nord), Znot(nz+nord), Phinot(nphi+nord))
  do i=1,nr
     Rtmp(i) = R_min + 1.d0*(i-1)/(nr-1) * (R_max - R_min)
  enddo
  call dbsnak (nr,Rtmp,nord, Rnot)
  do j=1,nz
     Ztmp(j) = Z_min + 1.d0*(j-1)/(nz-1) * (Z_max - Z_min)
  enddo
  call dbsnak (nz,Ztmp,nord, Znot)
  do k=1,nphi
     Phitmp(k) = 2.d0*pi/nsym*(k-1)/(nphi-1)
  enddo
  call dbsnak (nphi,Phitmp,nord, Phinot)


  ! calculate coefficients for spline interpolation
  allocate (Arcoeff(nr, nz, nphi), Azcoeff(nr, nz, nphi), Aphicoeff(nr, nz, nphi))

  write (6, *) 'Preparing B-spline coefficients for Ar ...'
  call dbs3in (nr,Rtmp,nz,Ztmp,nphi,Phitmp,Ar,nr,nz,nord,nord,nord,&
               Rnot,Znot,Phinot,Arcoeff)

  write (6, *) 'Preparing B-spline coefficients for Az ...'
  call dbs3in (nr,Rtmp,nz,Ztmp,nphi,Phitmp,Az,nr,nz,nord,nord,nord,&
               Rnot,Znot,Phinot,Azcoeff)

  write (6, *) 'Preparing B-spline coefficients for Aphi ...'
  call dbs3in (nr,Rtmp,nz,Ztmp,nphi,Phitmp,Aphi,nr,nz,nord,nord,&
               nord,Rnot,Znot,Phinot,Aphicoeff)


  ! cleanup
  deallocate (Rtmp, Ztmp, Phitmp)


  return
 1000 iconfig = 0
 1001 format (3x,'- Spline interpolation on cylindrical grid')
 1002 format (8x,'vector potential given on grid nodes')
 1003 format (8x,'magnetic field vector given on grid nodes')
 1004 format (8x,'spline order: ', i4)
 1010 format (8x,'grid resolution:        ',i4,' x ',i4,' x ',i4)
 1011 format (8x,'lower left corner:      ',e11.4,', ',e11.4)
 1012 format (8x,'upper right corner:     ',e11.4,', ',e11.4)
 1013 format (8x,'symmetry: ', i8)
 1020 format (8x,'reading data set: 'a)
 1021 format (8x,'   ... Ar components ...')
 1022 format (8x,'   ... Az components ...')
 1023 format (8x,'   ... Aphi components ...')
 2000 format ('error reading data file: ',a120)
  end subroutine splineB_load
!=======================================================================



!=======================================================================
  subroutine splineB_broadcast
  use parallel


  if (nprs == 1) return


  call broadcast_inte_s (nr)
  call broadcast_inte_s (nz)
  call broadcast_inte_s (nphi)
  call broadcast_inte_s (nsym)
  call broadcast_inte_s (nord)
  call broadcast_inte_s (data_type)

  if (mype == 0) then
     allocate (Rnot(nr+nord), Znot(nz+nord), Phinot(nphi+nord))
     allocate (Arcoeff(nr, nz, nphi), Azcoeff(nr, nz, nphi), Aphicoeff(nr, nz, nphi))
  endif

  call broadcast_real   (Rnot, nr+nord)
  call broadcast_real   (Znot, nz+nord)
  call broadcast_real   (Phinot, nphi+nord)
  call broadcast_real   (Arcoeff, nr*nz*nphi)
  call broadcast_real   (Azcoeff, nr*nz*nphi)
  call broadcast_real   (Aphicoeff, nr*nz*nphi)

  end subroutine splineB_broadcast
!=======================================================================
 


!=======================================================================
  function splineB_get_Bf(r3) result(Bf)
  use bspline
  use math
  real(real64), intent(in)  :: r3(3)
  real(real64)              :: Bf(3)

  real(real64) :: rr, zz, phi
  real(real64) :: dAz_dp, dAp_dz, dAr_dz, dAz_dr, Aphi, dAp_dr, dAr_dp


  ! convert units: cm -> m
  rr  = r3(1) / 100.d0
  zz  = r3(2) / 100.d0
  phi = phi_sym(r3(3),nsym)


  ! vector potential given on grid nodes
  if (data_type == VECTOR_POTENTIAL) then
  dAr_dz = dbs3dr(0,1,0,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Arcoeff)
  dAr_dp = dbs3dr(0,0,1,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Arcoeff)
  dAz_dr = dbs3dr(1,0,0,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Azcoeff)
  dAz_dp = dbs3dr(0,0,1,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Azcoeff)
  dAp_dr = dbs3dr(1,0,0,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Aphicoeff)
  dAp_dz = dbs3dr(0,1,0,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Aphicoeff)
  Aphi   = dbs3dr(0,0,0,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Aphicoeff)

  Bf(1) = 1.d0/rr * dAz_dp - dAp_dz
  Bf(2) = 1.d0/rr * Aphi + dAp_dr - 1.d0/rr * dAr_dp
  Bf(3) = dAr_dz - dAz_dr

  else ! data_type == BFIELD)
  ! interpolate magnetic field itself
  Bf(1) = dbs3dr(0,0,0,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Arcoeff)
  Bf(2) = dbs3dr(0,0,0,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Azcoeff)
  Bf(3) = dbs3dr(0,0,0,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Aphicoeff)
  endif

  ! convert units: T -> Gauss
  Bf    = Bf * 1.d4

  end function splineB_get_Bf
!=======================================================================



!=======================================================================
! Calculate Jacobian(Bf) [T/m] in cylindrical coordinates
!=======================================================================
  function splineB_get_JBf(r3) result(J)
  use bspline
  use math
  real(real64), intent(in)  :: r3(3)
  real(real64)              :: J(3,3)

  real(real64) :: rr, zz, phi
  real(real64) :: dAz_dp, dAp_dz, dAr_dz, dAz_dr, Aphi, dAp_dr, dAr_dp


  ! convert units: cm -> m
  rr  = r3(1) / 100.d0
  zz  = r3(2) / 100.d0
  phi = phi_sym(r3(3),nsym)

  if (data_type == VECTOR_POTENTIAL) then
     write (6, *) 'warning: JBf not yet implemented for vector potential!'
     J = 0.d0

  else
  ! interpolate magnetic field itself
  J(1,1) = dbs3dr(1,0,0,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Arcoeff)
  J(1,2) = dbs3dr(0,1,0,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Arcoeff)
  J(1,3) = dbs3dr(0,0,1,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Arcoeff)
  J(2,1) = dbs3dr(1,0,0,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Azcoeff)
  J(2,2) = dbs3dr(0,1,0,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Azcoeff)
  J(2,3) = dbs3dr(0,0,1,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Azcoeff)
  J(3,1) = dbs3dr(1,0,0,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Aphicoeff)
  J(3,2) = dbs3dr(0,1,0,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Aphicoeff)
  J(3,3) = dbs3dr(0,0,1,rr,zz,phi,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Aphicoeff)
  endif

  end function splineB_get_JBf
!=======================================================================

end module splineB
