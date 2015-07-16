!===============================================================================
! Field line reconstruction analysis
!===============================================================================
subroutine FLR_analysis
  use iso_fortran_env
  use parallel
  use flux_tube
  use dataset
  use math
  implicit none

  integer, parameter :: iu = 32

  character(len=72)  :: &
     Operation     = 'nothing', &
     Target_File   = 'analysis.conf', &
     Output_File   = 'output.dat', &
     Grid_File     = 'output.grid', &
     Cross_Section = ''

  real(real64)       :: x0(2), phi0, a1, alpha1, a2, alpha2, P, theta
  real(real64)       :: xi = 0.d0, eta = 0.d0
  integer            :: nphi, nlength

  namelist /Analysis_Input/ &
     Operation, Output_File, Grid_File, Cross_Section, &
     x0, phi0, a1, alpha1, a2, alpha2, P, theta, nphi, nlength, &
     xi, eta

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


  select case(Operation)
  case ('single_fieldline')
     call single_fieldline
  case ('fluxtube')
     call fluxtube
  case ('grid_accuracy')
     call grid_accuracy()
  case ('flux_conservation')
     call flux_conservation()
  case default
     write (6, *) 'error: invalid operation ', trim(Operation), '!'
     stop
  end select

  return
  contains
!=======================================================================



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



!=======================================================================
  subroutine grid_accuracy()
  use emc3_grid
  use fieldline
  use equilibrium
  use dataset
  use grid
  implicit none

  type(t_fieldline) :: Fo
  type(t_dataset)   :: D
  type(t_grid)      :: G1, G2
  real(real64), dimension(:,:), allocatable :: Fg
  integer, dimension(:), allocatable :: ir_list, ip_list

  real(real64) :: Dphi, xi(3), xo(3), ts, dx(-1:1), dr(-1:1), dt(-1:1), xg(3)
  real(real64) :: dPsidR, dPsidZ, Psio, Psig, dPsi, thetao, thetag, dtheta
  integer :: i, j, ig, n, idir, it(-1:1), iz, tm, tc, ierr, iscan


  write (6, 1000)
  call load_emc3_grid()


  ts = pi2 / 3600.d0
  tm = NM_AdamsBashforth4
  tc = FL_ANGLE


  iz    = 0
  iscan = 3
  select case(iscan)
  ! set up poloidal profile
  case(1)
     n  = ZON_POLO(iz)
     allocate (ip_list(n), ir_list(n))
     ir_list = 6
     !ir_list = ZON_RADI(0) - 4
     do i=1,ZON_POLO(iz)
        ip_list(i) = i-1
     enddo

  ! set up radial-poloidal scan
  case(3)
     n  = ZON_POLO(iz) * ZON_RADI(iz)
     ig = 1
     allocate (ip_list(n), ir_list(n))
     do i=0,ZON_RADI(iz)-1
        do j=0,ZON_POLO(iz)-1
           ip_list(ig) = j
           ir_list(ig) = i
           ig          = ig + 1
        enddo
     enddo

  end select


  ! Fo: field line from (o)rdinary differential equation
  ! Fg: field line reconstructed from (g)rid
  allocate (Fg(0:ZON_TORO(iz),3))
  call D%new(n, 7)
  call G1%new(LOCAL, UNSTRUCTURED, FIXED_COORD3, n)
  call G2%new(LOCAL, UNSTRUCTURED, FIXED_COORD3, n)
  it(-1) = 0
  it( 0) = ZON_TORO(iz) / 2
  it( 1) = ZON_TORO(iz)
  do i=1,n
     write (6, *) ir_list(i), ip_list(i)
     call reconstruct_field_line (iz, ir_list(i), ip_list(i), 0.d0, 0.d0, Fg)
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
        dx(idir) = sqrt((xg(1)-xo(1))**2 + (xg(2)-xo(2))**2)

        ! radial and poloidal displacement
        dPsidR = get_DPsiN(xo, 1, 0)
        dPsidZ = get_DPsiN(xo, 0, 1)
        Psig   = get_PsiN(xg)
        Psio   = get_PsiN(xo)
        dPsi   = Psig - Psio
        dr(idir) = dPsi / sqrt(dPsidR**2 + dPsidZ**2)

        thetag = get_poloidal_angle(xg)
        thetao = get_poloidal_angle(xo)
        dtheta = thetao - thetag; if (abs(dtheta) > pi) dtheta = dtheta - sign(pi2,dtheta)
        dt(idir) = dtheta
     enddo
     D%x(i,1)    = get_poloidal_angle(xi)
     D%x(i,2)    = dx(-1)
     D%x(i,3)    = dx( 1)
     D%x(i,4)    = dr(-1)
     D%x(i,5)    = dr( 1)
     D%x(i,6)    = dt(-1)
     D%x(i,7)    = dt( 1)
     G1%x(i,1) = ip_list(i)
     G1%x(i,2) = ir_list(i)
     G2%x(i,1) = D%x(i,1) / pi * 180.d0
     G2%x(i,2) = get_PsiN(xi)
  enddo
  call D%plot(filename=Output_File)
  call G1%store(filename='cell_indices.dat')
  call G2%store(filename=Grid_File)

  deallocate (ip_list, ir_list)
  deallocate (Fg)
 1000 format(3x,'- Running grid accuracy check')
  end subroutine grid_accuracy
!=======================================================================



!=======================================================================
  subroutine flux_conservation()
  use iso_fortran_env
  use emc3_grid
  use dataset
  use grid
  implicit none

  real(real64), dimension(:,:),   allocatable :: div, pitch
  real(real64), dimension(:,:,:), allocatable :: NL

  type(t_grid)    :: G
  type(t_dataset) :: D
  integer         :: iz, it0, i, j, ig, n


  write (6, 1000)
  call load_emc3_grid()
  call load_bfstren()

  iz  = 0
  it0 = ZON_TORO(iz) / 2
  allocate (div  (0:ZON_POLO(iz)-1, 0:ZON_RADI(iz)-1), &
            pitch(0:ZON_POLO(iz)-1, 0:ZON_RADI(iz)-1), &
            NL   (0:ZON_POLO(iz)-1, 0:ZON_RADI(iz)-1, 0:ZON_TORO(iz)))
  call check_mesh(iz, div, NL, pitch)

  n = (P_SURF_PL_TRANS_RANGE(2,iz) - P_SURF_PL_TRANS_RANGE(1,iz)) * &
      (R_SURF_PL_TRANS_RANGE(2,iz) - R_SURF_PL_TRANS_RANGE(1,iz))
  call D%new(n, 5)
  call G%new(LOCAL, UNSTRUCTURED, FIXED_COORD3, n)
  ig = 0
  do j=P_SURF_PL_TRANS_RANGE(1,iz),P_SURF_PL_TRANS_RANGE(2,iz)-1
     do i=R_SURF_PL_TRANS_RANGE(1,iz),R_SURF_PL_TRANS_RANGE(2,iz)-1
        ig = ig + 1
        G%x(ig,1) = j
        G%x(ig,2) = i
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
  end subroutine flux_conservation
!=======================================================================

end subroutine FLR_analysis
