!===============================================================================
! Field line reconstruction analysis
!===============================================================================
subroutine FLR_analysis
  use iso_fortran_env
  use parallel
  use flux_tube
  use dataset
  use math
  use emc3_grid
  implicit none

  integer, parameter :: iu = 32

  character(len=72)  :: &
     Operation     = 'nothing', &
     Target_File   = 'analysis.conf', &
     Output_File   = '', &
     Grid_File     = '', &
     Cross_Section = ''

  integer            :: &
     iz(2)         = -1, &
     ir(2)         = -1, &
     ip(2)         = -1


  real(real64)       :: x0(2), phi0, a1, alpha1, a2, alpha2, P, theta
  real(real64)       :: xi = 0.d0, eta = 0.d0
  integer            :: iz0, ir12(2), ip12(2), nphi, nlength

  namelist /Analysis_Input/ &
     Operation, Output_File, Grid_File, Cross_Section, &
     x0, phi0, a1, alpha1, a2, alpha2, P, theta, nphi, nlength, &
     xi, eta, &
     iz, ir, ip

  type(t_flux_tube)  :: FT


  ! initialize field line reconstruction analysis
  if (firstP) then
     write (6, *) 'Field line reconstruction analysis, target file: ', trim(Target_File)
     write (6, *)
  else
     return
  endif
  open  (iu, file=Target_File)
  read  (iu, Analysis_Input)
  close (iu)
  call load_emc3_grid()


  ! 1. set default values
  ! 1.1 file names
  !if (Output_File == '') Output_File = trim(Operation)//'.dat'
  !if (Grid_File   == '') Grid_File   = trim(Operation)//'.grid'
  write (6, 1000) Operation
  !write (6, 1001) trim(Output_File), trim(Grid_File)

  ! 1.2 indices (zone range)
  if (iz(1) < 0) then
     iz(1) = 0
     iz(2) = NZONET-1
  else
     if (iz(2) < 0) iz(2) = iz(1)
  endif



  select case(Operation)
  case ('single_fieldline')
     call single_fieldline
  case ('fluxtube')
     call fluxtube
  case ('grid_accuracy')
     do iz0=iz(1),iz(2)
        call setup_indices(iz0, ir12, ip12)
        call grid_accuracy(iz0, ir12, ip12)
     enddo
  case ('flux_conservation')
     do iz0=iz(1),iz(2)
        call setup_indices(iz0, ir12, ip12)
        call flux_conservation(iz0, ir12, ip12)
     enddo
  case ('flux_surface_discretization')
     call flux_surface_discretization()
  case default
     write (6, *) 'error: invalid operation ', trim(Operation), '!'
     stop
  end select

  return
 1000 format(3x,'- Operation: ',a)
 1001 format(8x,'Output files: ',a,' ',a)
  contains
!=======================================================================


!-----------------------------------------------------------------------
  subroutine setup_indices(iz_set, ir12, ip12)

  integer, intent(in)  :: iz_set
  integer, intent(out) :: ir12(2), ip12(2)

  ir12 = ir
  ! 1.2 indices (radial and poloidal range)
  if (ir(1) < 0) then
     ir12(1) = R_SURF_PL_TRANS_RANGE(1,iz_set)
     ir12(2) = R_SURF_PL_TRANS_RANGE(2,iz_set)-1
  else
     if (ir(2) < 0) ir12(2) = ir(1)
  endif

  ip12 = ip
  if (ip(1) < 0) then
     ip12(1) = P_SURF_PL_TRANS_RANGE(1,iz_set)
     ip12(2) = P_SURF_PL_TRANS_RANGE(2,iz_set)-1
  else
     if (ip(2) < 0) ip12(2) = ip(1)
  endif

  end subroutine setup_indices
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
  function load_cross_section() result(cs)
  type (t_cross_section) :: cs


  if (Cross_Section .ne. '') then
     call cs%load(Cross_Section)
     x0 = cs%a0
  else
     ! deg -> rad
     phi0  = phi0 / 180.d0 * Pi
     call cs%generate(x0, phi0, a1, alpha1, a2, alpha2, P, theta)
  endif
  write (6, 1000) cs%a0, cs%a1, cs%a2, cs%a3
  write (6, *)

 1000 format(3x,'- Shape coefficients of flux tube:',/,&
             8x,'a0 = (',f8.3,', ',f8.3,')',/, &
             8x,'a1 = (',f8.3,', ',f8.3,')',/, &
             8x,'a2 = (',f8.3,', ',f8.3,')',/, &
             8x,'a3 = (',f8.3,', ',f8.3,')')
  end function load_cross_section
!-----------------------------------------------------------------------



!=======================================================================
! Analyze field line with intrinsic coordinates (xi, eta) in flux tube
! (cs0, nphi, nlength). The initial cross section cs0 is either given in
! external file Cross_Section, or defined by geometry parameters
! (x0, phi0, a1, alpha1, a2, alpha2, P, theta)
!
! Output: nphi rows of
!     phi[deg]    err[cm]    err_rad[cm]    err_pol[cm]
!
! err = |(R,Z)_reconstructed  -  (R,Z)_traced|
! err_rad: contribution in e_Psi (radial) direction
! err_pol: contribution in poloidal direction (perpendicular to e_Psi)
!=======================================================================
  subroutine single_fieldline

  type(t_cross_section) :: cs0
  type(t_dataset) :: D
  integer         :: j


  ! get initial cross-section and generate flux tube
  cs0 = load_cross_section()
  call FT%generate(cs0, nphi, nlength)
  write (6, 1000) x0, phi0


  ! check intrinsic coordinates
  if ((abs(xi) > 1.d0)  .or.  (abs(eta) > 1.d0)) then
     write (6, *) 'error: intrinsic coordinates xi, eta must be in interval [-1,1]!'
     write (6, *) 'xi  = ', xi
     write (6, *) 'eta = ', eta
     stop
  endif


  ! analyze field line (xi, eta)
  D = FT%analyze_fieldline(xi, eta)
  open  (iu, file=Output_File)
  do j=1,nphi
     write (iu, 2000) FT%phi(j)*180.d0/pi, D%x(j,:)
  enddo
  close (iu)

 1000 format(3x,'- Analyzing single fieldline at (R,Z,phi) = ( ',f8.3,' cm, ',f8.3,' cm, ',f8.3,' deg)')
 2000 format(4e18.10)
  end subroutine single_fieldline
!=======================================================================



!=======================================================================
! Analysis of entire flux tube cross-section at final length 360 deg / nlength
!
! Output: nsample x nsample values of
!     phi[deg]    err[cm]    err_rad[cm]    err_pol[cm]
!=======================================================================
  subroutine fluxtube
  use grid
  use mesh_spacing
  type(t_cross_section) :: cs0
  type(t_dataset) :: Dfl, D
  type(t_grid)    :: G
  integer         :: i, j, ig, nsample


  ! get initial cross-section and generate flux tube
  cs0 = load_cross_section()
  call FT%generate(cs0, 1, nlength)
  write (6, 1000)


  ! sample cross-section of flux tube
  nsample = 64
  call G%new(LOCAL, UNSTRUCTURED, 3, nsample, nsample)
  call D%new(nsample*nsample, 3)
  ig = 0
  do j=0,nsample-1
     eta = -1.d0  +  2.d0 * Equidistant%node(j, nsample-1)
     do i=0,nsample-1
        xi = -1.d0  +  2.d0 * Equidistant%node(i, nsample-1)

        ig = ig + 1
        G%x(ig,1) = xi
        G%x(ig,2) = eta

        Dfl = FT%analyze_fieldline(xi, eta)
        D%x(ig,:) = Dfl%x(1,:)
     enddo
  enddo

  ! write data and grid to output files
  call G%store(Grid_File)
  call D%plot(filename=Output_File)

 1000 format(3x,'- Analyzing flux tube')
  end subroutine fluxtube
!=======================================================================

end subroutine FLR_analysis
!=======================================================================



!=======================================================================
  subroutine flux_conservation(iz, ir, ip)
  use iso_fortran_env
  use emc3_grid
  use dataset
  use grid
  use string
  implicit none

  integer, intent(in) :: iz, ir(2), ip(2)

  real(real64), dimension(:,:),   allocatable :: div, pitch
  real(real64), dimension(:,:,:), allocatable :: NL

  character(len=72) :: Output_File, Grid_File
  type(t_grid)    :: G
  type(t_dataset) :: D
  real(real64)    :: phi0
  integer         :: it0, i, j, ig, n, n1, n2


  Grid_File   = 'zone_'//trim(str(iz))//'.grid'
  Output_File = 'zone_'//trim(str(iz))//'_flux_conservation.dat'

  call load_bfstren()
  write (6, 1000)

  it0 = ZON_TORO(iz) / 2
  allocate (div  (0:ZON_POLO(iz)-1, 0:ZON_RADI(iz)-1), &
            pitch(0:ZON_POLO(iz)-1, 0:ZON_RADI(iz)-1), &
            NL   (0:ZON_POLO(iz)-1, 0:ZON_RADI(iz)-1, 0:ZON_TORO(iz)))
  call check_mesh(iz, div, NL, pitch)

  n1     = ir(2)-ir(1)+1
  n2     = ip(2)-ip(1)+1
  n      = n1 * n2
  ! setup mesh in real space
  phi0   = PHI_PLANE(it0 + PHI_PL_OS(iz))
  call G%new(CYLINDRICAL, MESH_2D, FIXED_COORD3, n1+1, n2+1, fixed_coord_value=phi0)
  do j=ip(1),ip(2)+1
  do i=ir(1),ir(2)+1
     ig = i + (j + it0*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
     G%mesh(i-ir(1),j-ip(1),1) = RG(ig)
     G%mesh(i-ir(1),j-ip(1),2) = ZG(ig)
  enddo
  enddo


  call D%new(n, 5)
  ig = 0
  write (6, 1001) ir(1), ir(2)
  write (6, 1002) ip(1), ip(2)
  do j=ip(1),ip(2)
     do i=ir(1),ir(2)
        ig = ig + 1
        D%x(ig,1) = div(j,i) * 100.d0
        D%x(ig,2) = pitch(j,i)
        D%x(ig,3) = NL(j,i,0)
        D%x(ig,4) = NL(j,i,it0)
        D%x(ig,5) = NL(j,i,ZON_TORO(iz))
     enddo
  enddo
  call D%plot(filename=Output_File)
  call G%store(filename=Grid_File)


  deallocate (div, pitch, NL)
 1000 format(3x,'- Running flux conservation and non-linearity check')
 1001 format(8x,'Radial cell range:   ',i0,' -> ',i0)
 1002 format(8x,'Poloidal cell range: ',i0,' -> ',i0)
  end subroutine flux_conservation
!=======================================================================



!=======================================================================
  subroutine flux_surface_discretization()
  use equilibrium
  use flux_surface_2D
  implicit none

  type(t_flux_surface_2D) :: F1, F2
  real(real64) :: Psi1, Psi2, y(3), r1(3), r2(3)
  integer      :: ierr


  Psi1  = 0.990d0
  Psi2  = 0.995d0

  y     = 0.d0
  y(2)  = Psi1; r1 = get_cylindrical_coordinates(y, ierr)
  if (ierr .ne. 0) then
     write (6, 9000) Psi1, r1
     stop
  endif
  y(2)  = Psi2; r2 = get_cylindrical_coordinates(y, ierr)
  if (ierr .ne. 0) then
     write (6, 9000) Psi2, r2
     stop
  endif


  ! generate flux surfaces Psi1, Psi2
  call F1%generate_closed(r1, RIGHT_HANDED)
  call F2%generate_closed(r2, RIGHT_HANDED)

  return
 9000 format('error in subroutine flux_surface_discretization: could not find real space coordinates for Psi = ', f10.5,//3e10.5)
  end subroutine flux_surface_discretization
!=======================================================================





!=======================================================================
  subroutine grid_accuracy(iz, ir, ip)
  use emc3_grid
  use fieldline
  use equilibrium, only: get_PsiN, get_DPsiN, get_poloidal_angle, get_ePsi, get_rmin
  use dataset
  use grid
  use string
  implicit none

  integer, intent(in) :: iz, ir(2), ip(2)

  type(t_fieldline) :: Fo
  type(t_dataset)   :: D
  type(t_grid)      :: G
  real(real64), dimension(:,:), allocatable :: Fg

  character(len=72) :: Output_File, Grid_File
  real(real64) :: Dphi, xi(3), xo(3), ts, dx(-1:1), dr(-1:1), dt(-1:1), xg(3)
  real(real64) :: dist(2), ePsi(2), ePol(2), rmin
  real(real64) :: dPsidR, dPsidZ, Psio, Psig, dPsi, thetao, thetag, dtheta, phi0
  integer :: i, j, ig, n, n1, n2, idir, it(-1:1), tm, tc, ierr, iscan


  write (6, 1000)
  write (6, *) 'cell range:'
  write (6, *) 'ir = ', ir(1), ' -> ', ir(2)
  write (6, *) 'ip = ', ip(1), ' -> ', ip(2)

  Grid_File   = 'zone_'//trim(str(iz))//'.grid'
  Output_File = 'zone_'//trim(str(iz))//'_accuracy.dat'


  ts = pi2 / 3600.d0
  tm = NM_AdamsBashforth4
  tc = FL_ANGLE


  ! Fo: field line from (o)rdinary differential equation
  ! Fg: field line reconstructed from (g)rid
  n1     = ir(2)-ir(1)+1
  n2     = ip(2)-ip(1)+1
  n      = n1 * n2
  it(-1) = 0
  it( 0) = ZON_TORO(iz) / 2
  it( 1) = ZON_TORO(iz)

  ! setup mesh in real space
  phi0   = PHI_PLANE(it(0) + PHI_PL_OS(iz))
  call G%new(CYLINDRICAL, MESH_2D, FIXED_COORD3, n1+1, n2+1, fixed_coord_value=phi0)
  do j=ip(1),ip(2)+1
  do i=ir(1),ir(2)+1
     ig = i + (j + it(0)*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
     G%mesh(i-ir(1),j-ip(1),1) = RG(ig)
     G%mesh(i-ir(1),j-ip(1),2) = ZG(ig)
  enddo
  enddo

  allocate (Fg(0:ZON_TORO(iz),3))
  call D%new(n, 8)
  ig     = 1
  do j=ip(1),ip(2)
  do i=ir(1),ir(2)
     write (6, *) i, j
     call reconstruct_field_line (iz, i, j, 0.d0, 0.d0, Fg)
     ! initial point for field line tracing
     xi = Fg(it(0),:)

     do idir=-1,1,2
        Dphi = abs(Fg(it(idir), 3) - Fg(it(0), 3))
        call Fo%init(xi, idir*ts, tm, tc)
        call Fo%trace_Dphi(Dphi, .false., xo, ierr)
        if (ierr .ne. 0) then
           write (6, *) 'error in subroutine trace_grid: ', &
                        'trace_Dphi returned error ', ierr
           stop
        endif

        ! evaluate field line tracing vs. reconstruction
        xg = Fg(it(idir),:)

        ! absolute displacement
        dist     = xg(1:2) - xo(1:2)
        dx(idir) = sqrt((xg(1)-xo(1))**2 + (xg(2)-xo(2))**2)

        ! radial and poloidal displacement
        ePsi     = get_ePsi(xi)
        ePol(1) = ePsi(2); ePol(2) = -ePsi(1)

        dPsidR = get_DPsiN(xo, 1, 0)
        dPsidZ = get_DPsiN(xo, 0, 1)
        Psig   = get_PsiN(xg)
        Psio   = get_PsiN(xo)
        dPsi   = Psig - Psio
        dr(idir) = dPsi / sqrt(dPsidR**2 + dPsidZ**2)
!        dr(idir) = sum(ePsi*dist)

        thetag = get_poloidal_angle(xg)
        thetao = get_poloidal_angle(xo)
        dtheta = thetao - thetag; if (abs(dtheta) > pi) dtheta = dtheta - sign(pi2,dtheta)
        rmin   = 0.5d0 * (get_rmin(xo) + get_rmin(xg))
        dt(idir) = dtheta * rmin
!        dt(idir) = sum(ePol*dist)
     enddo
     D%x(ig,1) = get_poloidal_angle(xi)
     D%x(ig,2) = dx(-1)
     D%x(ig,3) = dx( 1)
     D%x(ig,4) = dr(-1)
     D%x(ig,5) = dr( 1)
     D%x(ig,6) = dt(-1)
     D%x(ig,7) = dt( 1)
     D%x(ig,8) = flux_tube_length(iz, i, j)
     ig        = ig + 1
  enddo
  enddo
  call D%plot(filename=Output_File)
  call G%store(filename=Grid_File)

  deallocate (Fg)
 1000 format(3x,'- Running grid accuracy check')
  end subroutine grid_accuracy
!=======================================================================
