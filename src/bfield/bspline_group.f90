module bspline_group
  use iso_fortran_env
  use abstract_bfield
  use bspline3d
  implicit none
  private
  public :: t_bspline_group, t_bspline_potential, t_bspline_bfield


  type, abstract, extends(t_bfield) :: t_bspline_group
     type(t_bspline3d) :: bspline(3)
     integer   :: nfp
     contains
     procedure :: broadcast
  end type t_bspline_group


  ! B-spline interpolation of vector potential (A)
  type, extends(t_bspline_group) :: t_bspline_potential
     contains
     procedure :: get_Bf  => get_Bf_bspline_potential
     procedure :: get_JBf => get_JBf_bspline_potential
  end type t_bspline_potential


  ! B-spline interpolation of magnetic field (B)
  type, extends(t_bspline_group) :: t_bspline_bfield
     contains
     procedure :: get_Bf  => get_Bf_bspline_bfield
     procedure :: get_JBf => get_JBf_bspline_bfield
  end type t_bspline_bfield


  interface t_bspline_group
     procedure constructor
  end interface


  contains
  !=====================================================================


  !=====================================================================
  subroutine broadcast(this)
  use parallel
  class(t_bspline_group)       :: this

  integer :: i


  call this%broadcast_bounds()
  call broadcast_inte_s(this%nfp)
  do i=1,3
     call this%bspline(i)%broadcast()
  enddo

  end subroutine broadcast
  !=====================================================================



  !=====================================================================
  ! Implementation of get_Bf for interpolation of vector potential
  !=====================================================================
  function get_Bf_bspline_potential(this, r) result(Bf)
  use math
  class(t_bspline_potential)   :: this
  real(real64), intent(in)     :: r(3)
  real(real64)                 :: Bf(3)

  real(real64) :: dAz_dp, dAp_dz, dAr_dz, dAz_dr, Aphi, dAp_dr, dAr_dp
  real(real64) :: rsym(3), rx1


  Bf        = 0.d0;   if (this%out_of_bounds(r)) return
  rsym(1:2) = r(1:2)
  rsym(3)   = phi_sym(r(3), this%nfp)
  rx1       = 1.d0 / rsym(1)

  dAr_dz = this%bspline(1)%derivative(rsym, 0, 1, 0)
  dAr_dp = this%bspline(1)%derivative(rsym, 0, 0, 1)
  dAz_dr = this%bspline(2)%derivative(rsym, 1, 0, 0)
  dAz_dp = this%bspline(2)%derivative(rsym, 0, 0, 1)
  dAp_dr = this%bspline(3)%derivative(rsym, 1, 0, 0)
  dAp_dz = this%bspline(3)%derivative(rsym, 0, 1, 0)
  Aphi   = this%bspline(3)%eval(rsym)

  Bf(1) = rx1 * dAz_dp - dAp_dz
  Bf(2) = rx1 * (Aphi - dAr_dp) + dAp_dr
  Bf(3) = dAr_dz - dAz_dr

  end function get_Bf_bspline_potential
  !=====================================================================



  !=====================================================================
  ! Implementation of get_JBf for interpolation of vector potential
  !=====================================================================
  function get_JBf_bspline_potential(this, r) result(JBf)
  use math
  class(t_bspline_potential)   :: this
  real(real64), intent(in)     :: r(3)
  real(real64)                 :: JBf(3,3)

  integer, parameter :: dr(3) = (/1, 0, 0/), dz(3) = (/0, 1, 0/), dp(3) = (/0, 0, 1/)

  real(real64) :: dAr2(3,3), dAz2(3,3), dAp2(3,3)
  real(real64) :: Aphi, dAr_dp, dAz_dp, dAp_dr, dAp_dz, dAp_dp
  real(real64) :: rsym(3), rx1, rx2

  integer :: i, j, ddr, ddz, ddp


  JBf       = 0.d0;   if (this%out_of_bounds(r)) return
  rsym(1:2) = r(1:2)
  rsym(3)   = phi_sym(r(3), this%nfp)
  rx1       = 1.d0 / rsym(1);   rx2 = rx1**2

  dAr_dp = this%bspline(1)%derivative(rsym, 0, 0, 1)
  dAz_dp = this%bspline(2)%derivative(rsym, 0, 0, 1)
  dAp_dr = this%bspline(3)%derivative(rsym, 1, 0, 0)
  dAp_dz = this%bspline(3)%derivative(rsym, 0, 1, 0)
  dAp_dp = this%bspline(3)%derivative(rsym, 0, 0, 1)
  Aphi   = this%bspline(3)%eval(rsym)
  do i=1,3
  do j=i,3
     ddr = dr(i) + dr(j)
     ddz = dz(i) + dz(j)
     ddp = dp(i) + dp(j)
     dAr2(i,j) = this%bspline(1)%derivative(rsym, ddr, ddz, ddp)
     dAz2(i,j) = this%bspline(2)%derivative(rsym, ddr, ddz, ddp)
     dAp2(i,j) = this%bspline(3)%derivative(rsym, ddr, ddz, ddp)
  enddo
  enddo

  JBf(1,1) = -rx2 * dAz_dp + rx1 * dAz2(1,3) - dAp2(1,2)
  JBf(1,2) = rx1 * dAz2(2,3) - dAp2(2,2)
  JBf(1,3) = rx1 * dAz2(3,3) - dAp2(2,3)
  JBf(2,1) = -rx2 * (Aphi - dAr_dp) + rx1 * (dAp_dr - dAr2(1,3)) + dAp2(1,1)
  JBf(2,2) = rx1 * (dAp_dz - dAr2(2,3)) + dAp2(1,2)
  JBf(2,3) = rx1 * (dAp_dp - dAr2(3,3)) + dAp2(1,3)
  JBf(3,1) = dAr2(1,2) - dAz2(1,1)
  JBf(3,2) = dAr2(2,2) - dAz2(1,2)
  JBf(3,3) = dAr2(2,3) - dAz2(1,3)

  end function get_JBf_bspline_potential
  !=====================================================================



  !=====================================================================
  ! Implementation of get_Bf for interpolation of magnetic field vector
  !=====================================================================
  function get_Bf_bspline_bfield(this, r) result(Bf)
  use math
  class(t_bspline_bfield)      :: this
  real(real64), intent(in)     :: r(3)
  real(real64)                 :: Bf(3)

  real(real64) :: rsym(3)


  Bf        = 0.d0;   if (this%out_of_bounds(r)) return
  rsym(1:2) = r(1:2)
  rsym(3)   = phi_sym(r(3), this%nfp)

  Bf(1) = this%bspline(1)%eval(rsym)
  Bf(2) = this%bspline(2)%eval(rsym)
  Bf(3) = this%bspline(3)%eval(rsym)

  end function get_Bf_bspline_bfield
  !=====================================================================



  !=====================================================================
  ! Implementation of get_JBf for interpolation of magnetic field vector
  !=====================================================================
  function get_JBf_bspline_bfield(this, r) result(JBf)
  use math
  class(t_bspline_bfield)      :: this
  real(real64), intent(in)     :: r(3)
  real(real64)                 :: JBf(3,3)

  real(real64) :: rsym(3)


  JBf       = 0.d0;   if (this%out_of_bounds(r)) return
  rsym(1:2) = r(1:2)
  rsym(3)   = phi_sym(r(3), this%nfp)

  JBf(1,1) = this%bspline(1)%derivative(rsym, 1, 0, 0)
  JBf(1,2) = this%bspline(1)%derivative(rsym, 0, 1, 0)
  JBf(1,3) = this%bspline(1)%derivative(rsym, 0, 0, 1)
  JBf(2,1) = this%bspline(2)%derivative(rsym, 1, 0, 0)
  JBf(2,2) = this%bspline(2)%derivative(rsym, 0, 1, 0)
  JBf(2,3) = this%bspline(2)%derivative(rsym, 0, 0, 1)
  JBf(3,1) = this%bspline(3)%derivative(rsym, 1, 0, 0)
  JBf(3,2) = this%bspline(3)%derivative(rsym, 0, 1, 0)
  JBf(3,3) = this%bspline(3)%derivative(rsym, 0, 0, 1)

  end function get_JBf_bspline_bfield
  !=====================================================================



  !=====================================================================
  subroutine read_text_data(filename, &
                     nr, nz, np, nfp, rmin, rmax, zmin, zmax, Qr, Qz, Qp)

  character(len=*), intent(in)  :: filename
  integer,          intent(out) :: nr, nz, np, nfp
  real(real64),     intent(out) :: rmin, rmax, zmin, zmax
  real(real64), dimension(:,:,:), allocatable, intent(out) :: Qr, Qz, Qp

  integer, parameter :: iu = 70

  character(len=256) :: prefix, src(3)
  character(len=72)  :: srcopt
  logical :: ex
  integer :: i, iu_src(3), j, k


  open  (iu, file=filename)
  read  (iu, *) nr, nz, np, nfp, rmin, rmax, zmin, zmax
  srcopt = 'NONE';   src = ''
  read  (iu, *, end=8001) srcopt, src
 8001 continue

  allocate (Qr(nr,nz,np+1), Qz(nr,nz,np+1), Qp(nr,nz,np+1))
  select case(srcopt)
  case('TEXT')
     iu_src = iu
  case('FILE')
     ! basename of data directory
     i = scan(trim(filename), '/', back=.true.)
     prefix = '';   if (i > 0) prefix = filename(1:i)

     ! open source files
     do i=1,3
        iu_src(i) = iu + i
        src(i)    = trim(prefix)//trim(src(i))
        inquire (file=src(i), exist=ex)
        if (.not.ex) then
           write (6, 9002) trim(src(i));   stop
        endif

        open  (iu_src(i), file=src(i))
     enddo
  case default
     write (6, 9001) trim(srcopt);   stop
  end select


  ! read data
  write (6, 1020)
  write (6, 1021)
  read  (iu_src(1), *) (((Qr(i,j,k), k=1,np), j=1,nz), i=1,nr)
  write (6, 1022)
  read  (iu_src(2), *) (((Qz(i,j,k), k=1,np), j=1,nz), i=1,nr)
  write (6, 1023)
  read  (iu_src(3), *) (((Qp(i,j,k), k=1,np), j=1,nz), i=1,nr)
  write (6, *)
  np = np + 1


  ! close data files
  close (iu)
  if (srcopt == 'FILE') then
  do i=1,3
     close (iu_src(i))
  enddo
  endif


  ! set up data at upper boundary
  Qr(:,:,np) = Qr(:,:,1)
  Qz(:,:,np) = Qz(:,:,1)
  Qp(:,:,np) = Qp(:,:,1)


 1020 format (8x,'reading data set: ')
 1021 format (8x,'   ... R components ...')
 1022 format (8x,'   ... Z components ...')
 1023 format (8x,'   ... phi components ...')
 9001 format('error: invalid option ',a,'!')
 9002 format('error: data file ',a,' does not exist!')
  end subroutine read_text_data
  !=====================================================================



  !=====================================================================
  subroutine read_netcdf_data(filename, RDIM_NAME, ZDIM_NAME, PDIM_NAME, NFP_NAME, &
                     RMIN_NAME, RMAX_NAME, ZMIN_NAME, ZMAX_NAME, &
                     QR_NAME, QZ_NAME, QP_NAME, &
                     nr, nz, np, nfp, rmin, rmax, zmin, zmax, Qr, Qz, Qp, &
                     NGROUP_NAME, LABEL_NAME, &
                     scale_subset)
  use iso_fortran_env
  use netcdf
  implicit none

  character(len=*), intent(in)  :: filename
  character(len=*), intent(in)  :: RDIM_NAME, ZDIM_NAME, PDIM_NAME, NFP_NAME, &
                                   RMIN_NAME, RMAX_NAME, ZMIN_NAME, ZMAX_NAME, &
                                   QR_NAME, QZ_NAME, QP_NAME
  integer,          intent(out) :: nr, nz, np, nfp
  real(real64),     intent(out) :: rmin, rmax, zmin, zmax
  real(real64), dimension(:,:,:), allocatable, intent(out) :: Qr, Qz, Qp
  character(len=*),           intent(in), optional :: NGROUP_NAME, LABEL_NAME
  real(real64), dimension(:), intent(in), optional :: scale_subset

  real(real64), dimension(:),     allocatable :: scale_factor
  real(real64), dimension(:,:,:), allocatable :: tmpr, tmpz, tmpp
  character(len=30), dimension(:), allocatable :: curlabel
  character(len=len(QR_NAME)+4) :: vname
  logical :: add_suffix
  real(real64) :: r
  integer :: i, ngroup, ncid, varid, jr, jz


  ! open data file
  call check(nf90_open(filename, NF90_NOWRITE, ncid))

  ! read dimensions
  call check(nf90_inq_dimid(ncid, RDIM_NAME, varid))
  call check(nf90_inquire_dimension(ncid, varid, len=nr))
  call check(nf90_inq_dimid(ncid, ZDIM_NAME, varid))
  call check(nf90_inquire_dimension(ncid, varid, len=nz))
  call check(nf90_inq_dimid(ncid, PDIM_NAME, varid))
  call check(nf90_inquire_dimension(ncid, varid, len=np))

  ! read layout
  call check(nf90_inq_varid(ncid, NFP_NAME, varid))
  call check(nf90_get_var(ncid, varid, nfp))
  call check(nf90_inq_varid(ncid, RMIN_NAME, varid))
  call check(nf90_get_var(ncid, varid, rmin))
  call check(nf90_inq_varid(ncid, RMAX_NAME, varid))
  call check(nf90_get_var(ncid, varid, rmax))
  call check(nf90_inq_varid(ncid, ZMIN_NAME, varid))
  call check(nf90_get_var(ncid, varid, zmin))
  call check(nf90_inq_varid(ncid, ZMAX_NAME, varid))
  call check(nf90_get_var(ncid, varid, zmax))

  ! data groups and scalings
  ngroup = 1
  add_suffix = .false.
  if (present(NGROUP_NAME)) then
     call check(nf90_inq_varid(ncid, NGROUP_NAME, varid))
     call check(nf90_get_var(ncid, varid, ngroup))
     add_suffix = .true.
  endif
  allocate (scale_factor(ngroup));   scale_factor = 1.d0
  if (present(scale_subset)) then
     if (size(scale_subset) < ngroup) then
        write (6, *) 'error: size of variable scale_subset too small!'
        stop
     endif
     scale_factor = scale_subset(1:ngroup)
  endif
  if (present(LABEL_NAME)) then
     allocate (curlabel(ngroup))
     call check(nf90_inq_varid(ncid, LABEL_NAME, varid))
     call check(nf90_get_var(ncid, varid, curlabel))
  endif

  ! read data
  allocate (Qr(nr,nz,np+1), Qz(nr,nz,np+1), Qp(nr,nz,np+1))
  allocate (tmpr(nr,nz,np), tmpz(nr,nz,np), tmpp(nr,nz,np))
  Qr = 0.d0;   Qz = 0.d0;   Qp = 0.d0
  do i=1,ngroup
     if (scale_factor(i) == 0) cycle

     if (present(LABEL_NAME)) write (6, 1010) trim(curlabel(i)), scale_factor(i)

     !if (i>1) scale_factor(i) = 0.d0
     vname = QR_NAME;   if (add_suffix) write (vname, 1000) QR_NAME, i
     call check(nf90_inq_varid(ncid, vname, varid))
     call check(nf90_get_var(ncid, varid, tmpr))
     vname = QZ_NAME;   if (add_suffix) write (vname, 1000) QZ_NAME, i
     call check(nf90_inq_varid(ncid, vname, varid))
     call check(nf90_get_var(ncid, varid, tmpz))
     vname = QP_NAME;   if (add_suffix) write (vname, 1000) QP_NAME, i
     call check(nf90_inq_varid(ncid, vname, varid))
     call check(nf90_get_var(ncid, varid, tmpp))

     Qr(:,:,1:np) = Qr(:,:,1:np) + scale_factor(i) * tmpr
     Qz(:,:,1:np) = Qz(:,:,1:np) + scale_factor(i) * tmpz
     Qp(:,:,1:np) = Qp(:,:,1:np) + scale_factor(i) * tmpp
  enddo
  if (present(LABEL_NAME)) write (6, *)
 1000 format(a,'_',i3.3)
 1010 format(8x,a,', scale factor: ',f0.3)
  np = np+1
  Qr(:,:,np) = Qr(:,:,1)
  Qz(:,:,np) = Qz(:,:,1)
  Qp(:,:,np) = Qp(:,:,1)

  ! close data file and cleanup
  call check( nf90_close(ncid) )
  deallocate (tmpr, tmpz, tmpp, scale_factor)

  contains

  subroutine check(status)
  integer, intent(in) :: status

  if (status /= nf90_noerr) then
     print *, trim(nf90_strerror(status))
     stop "Stopped"
  endif

  end subroutine check

  end subroutine read_netcdf_data
  !=====================================================================



  !=====================================================================
  function constructor(source, source_type, amplify, scale_subset)
  use math
  character(len=*),           intent(in)           :: source, source_type
  real(real64),               intent(in)           :: amplify
  real(real64), dimension(:), intent(in), optional :: scale_subset
  class(t_bspline_group), allocatable              :: constructor

  real(real64), dimension(:,:,:), allocatable :: Qr, Qz, Qp

  real(real64) :: rmin, rmax, zmin, zmax, phi(1)
  integer :: itype, nr, nz, np, nfp, m
  integer :: spline_order


  ! allocate variable and set spline order
  select case(source_type)
  case('bmw', 'mgrid', 'atxt')
     allocate (t_bspline_potential :: constructor);   spline_order = 5

  case('mgrid_bfield', 'btxt')
     allocate (t_bspline_bfield    :: constructor);   spline_order = 4

  case default
     write (6, 9000) trim(source_type)
     stop
  end select
 9000 format('error: unkown format ',a,' for source ',a,'!')


  ! read data
  select case(source_type)
  case('mgrid_bfield')
     write (6, 1000)
     call read_netcdf_data(source, 'rad', 'zee', 'phi', 'nfp', &
                 'rmin', 'rmax', 'zmin', 'zmax', &
                 'br', 'bz', 'bp', &
                 nr, nz, np, nfp, rmin, rmax, zmin, zmax, Qr, Qz, Qp, &
                 'nextcur', 'coil_group', scale_subset)

  case('mgrid')
     write (6, 1001)
     call read_netcdf_data(source, 'rad', 'zee', 'phi', 'nfp', &
                 'rmin', 'rmax', 'zmin', 'zmax', &
                 'ar', 'az', 'ap', &
                 nr, nz, np, nfp, rmin, rmax, zmin, zmax, Qr, Qz, Qp, &
                 'nextcur', 'coil_group', scale_subset)

  case('bmw')
     write (6, 2000)
     call read_netcdf_data(source, 'r', 'z', 'phi', 'nfp', &
                 'rmin', 'rmax', 'zmin', 'zmax', &
                 'ar_grid', 'az_grid', 'ap_grid', &
                 nr, nz, np, nfp, rmin, rmax, zmin, zmax, Qr, Qz, Qp)

  case('atxt')
     write (6, 3000)
     call read_text_data(source, &
                 nr, nz, np, nfp, rmin, rmax, zmin, zmax, Qr, Qz, Qp)


  case('btxt')
     write (6, 4000)
     call read_text_data(source, &
                 nr, nz, np, nfp, rmin, rmax, zmin, zmax, Qr, Qz, Qp)
  end select


  write (6, 8000) amplify
  write (6, 8001) nr, nz, np
  write (6, 8002) rmin*1.d2, rmax*1.d2
  write (6, 8003) zmin*1.d2, zmax*1.d2
  write (6, 8004) nfp
  write (6, *)
  constructor%rmin = rmin;   constructor%rmax = rmax
  constructor%zmin = zmin;   constructor%zmax = zmax
  constructor%nfp  = nfp;    phi = 2.d0 * pi / nfp
  Qr  = amplify * Qr;   Qz  = amplify * Qz;   Qp  = amplify * Qp
  m   = spline_order
  call constructor%bspline(1)%init(nr, nz, np, (/rmin, rmax/), (/zmin, zmax/), phi, Qr, m)
  call constructor%bspline(2)%init(nr, nz, np, (/rmin, rmax/), (/zmin, zmax/), phi, Qz, m)
  call constructor%bspline(3)%init(nr, nz, np, (/rmin, rmax/), (/zmin, zmax/), phi, Qp, m)
  deallocate (Qr, Qz, Qp)

 1000 format(3x,'- Magnetic field from mgrid file')
 1001 format(3x,'- Vector potential from mgrid file')
 2000 format(3x,'- Vector potential from BMW file')
 3000 format(3x,'- Vector potential on cylindrical grid')
 4000 format(3x,'- Magnetic field on cylindrical grid')
 8000 format(8x,'Scale factor: ',f0.3)
 8001 format(8x,'B-Spline interpolation on ',i0,' x ',i0,' x ',i0,' grid')
 8002 format(8x,'Radial domain:   ',f0.3,' -> ',f0.3,' cm')
 8003 format(8x,'Vertical domain: ',f0.3,' -> ',f0.3,' cm')
 8004 format(8x,'Number of field periods: ',i0)
  end function constructor
  !=====================================================================

end module bspline_group
