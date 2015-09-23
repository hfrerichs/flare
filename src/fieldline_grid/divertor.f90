module divertor
  use iso_fortran_env
  use emc3_grid
  use curve2D
  use separatrix
  implicit none

  integer, parameter :: iud     = 72


  ! magnetic separatrix
  type(t_separatrix), dimension(:), allocatable :: S

  ! guiding_surface
  type(t_curve), public      :: C_guide, C_cutL, C_cutR

  ! store magnetic axis
  real(real64) :: Pmag(2)


  contains
!=======================================================================


  !=====================================================================
  ! setup geometry of simulation domain
  ! based on nx X-points Xp(1:nx) from module equilibrium
  ! and how they connect to each other (connectX(ix) = number of X-point to which separatrix from X-point ix connects to)
  !=====================================================================
  subroutine setup_geometry(nx, connectX)
  use math
  use equilibrium
  use boundary
  use fieldline_grid, only: guiding_surface, d_cutL, d_cutR
  use inner_boundary
  integer, intent(in) :: nx, connectX(nx)

  character(len=2) :: Sstr
  real(real64)     :: tmp(3)
  integer          :: ix


  ! 1.a set up guiding surface for divertor legs (C_guide) ---------------------
  if (guiding_surface .ne. '') then
     write (6, 1000)
     call C_guide%load(guiding_surface)
  else if (n_axi > 0) then
     write (6, 1001)
     call C_guide%copy(S_axi(1))
  else
     write (6, *) 'error: cannot determine divertor geometry!'
     write (6, *) 'neither guiding_surface is set, nor is an axisymmetric surface defined.'
     stop
  endif
 1000 format(8x,'User defined guiding surface for divertor strike points')
 1001 format(8x,'First axisymmetric surface used for divertor strike points')

  ! 1.b set up extended guiding surfaces for divertor leg discretization -------
  ! C_cutL, C_cutR
  call C_cutL%copy(C_guide)
  call C_cutL%left_hand_shift(d_cutL(1))
  call C_cutL%plot(filename='C_cutL.plt')
  call C_cutR%copy(C_guide)
  call C_cutR%left_hand_shift(d_cutR(1))
  call C_cutR%plot(filename='C_cutR.plt')


  ! 2.a set up magnetic axis (Pmag) --------------------------------------
  tmp = get_magnetic_axis(0.d0); Pmag = tmp(1:2)
  write (6, 2000) Pmag
 2000 format(8x,'Magnetic axis at: ',2f10.4)

  ! 2.b check X-points
  do ix=1,nx
     tmp(1:2) = Xp(ix)%load()
     write (6, 2001) ix, tmp(1:2)
     write (6, 2002) Xp(ix)%theta/pi*180.d0
  enddo
  write (6, *)
 2001 format(8x,i0,'. X-point at: ',2f10.4)
 2002 format(11x,'-> poloidal angle [deg]: ',f10.4)

  ! 2.c generate separatrix(ces)
  allocate(S(nx))
  do ix=1,nx
     call S(ix)%generate(ix, C_cutL=C_cutL, C_cutR=C_cutR, iconnect=connectX(ix))
     write (Sstr, 2003) ix; call S(ix)%plot(trim(Sstr), parts=.true.)

     call S(ix)%M1%setup_length_sampling()
     call S(ix)%M2%setup_length_sampling()
     call S(ix)%M3%setup_length_sampling()
     call S(ix)%M4%setup_length_sampling()
  enddo
 2003 format('S',i0)


  ! 3. inner boundaries for EMC3 grid
  call load_inner_boundaries(Xp(1)%theta)


  end subroutine setup_geometry
  !=====================================================================



  !=====================================================================
  subroutine check_domain()
  use fieldline_grid, only: layers, rpath, guiding_surface, label

  real(real64) :: x(2), tau, L
  integer      :: i


  do i=0,layers-1
     if (rpath(i)%intersect_curve(C_guide, x, tau)) then
        write (6, 9000) i
        write (6, 9001)

        L = rpath(i)%length() * tau
        write (6, 9002) trim(label(i)), L
     endif
  enddo

 9000 format('WARNING: reference path for radial discretization (rpath_I.plt, I = ',i0,') crosses guiding surface!')
 9001 format('This might cause a problem for adjusting the flux surface discretization to the strike point position!')
 9002 format('chosing d_', a, ' < ', f0.3, ' may be required!')
  end subroutine check_domain
  !=====================================================================



  !=====================================================================
  subroutine divide_SOL2(F, eta, side, alpha, r, C)
  use math
  use flux_surface_2D
  type(t_flux_surface_2D), intent(in)  :: F
  real(real64),            intent(in)  :: eta, alpha, r
  integer,                 intent(in)  :: side
  type(t_curve),           intent(out) :: C(2)

  real(real64) :: l, alpha1, xi(1)


  l = F%length()

  alpha1 = 1.d0 + eta * (alpha - 1.d0)
  select case(side)
  case(1)
     xi  = alpha1 * r / l
  case(-1)
     xi  = 1.d0 - alpha1 * r / l
  case default
     write (6, *) 'error in subroutine divide_SOL2: side = 1 or -1 required!'
     stop
  end select
  !write (6, *) '(xi, eta, alpha, r, l) = ', xi, eta, alpha, r, l

  call F%splitn(2, xi, C)
  !call CR%setup_length_sampling()
  !call C0%setup_sampling(Xp(1)%X, Xp(1)%X, Magnetic_Axis%X, eta, eta, pi2, Dtheta_sampling)
  !call CL%setup_length_sampling()

  end subroutine divide_SOL2
  !=====================================================================



  !=====================================================================
  subroutine divide_SOL3(F, eta, CL, C0, CR, ix1, ix2)
  use math
  use flux_surface_2D
  use equilibrium
  use fieldline_grid, only: Dtheta_sampling, Dtheta_separatrix, alphaL, alphaR
  type(t_flux_surface_2D), intent(in)  :: F
  real(real64),            intent(in)  :: eta
  type(t_curve),           intent(out) :: CL, CR
  type(t_flux_surface_2D), intent(out) :: C0
  integer,                 intent(in)  :: ix1, ix2

  real(real64) :: l, alpha, xiR, xiL, eta1, eta2, DthetaR, DthetaL
  real(real64) :: x(2), thetaR, thetaL
  integer      :: jx1, jx2


  ! set reference length
  l = F%length()
  jx1 = abs(ix1)
  jx2 = abs(ix2)


  ! in outer SOL set eta'=1+eta for divertor legs from inner SOL
  eta1 = eta; eta2 = eta
  if (ix1 < 0) eta1 = eta1 + 1.d0
  if (ix2 < 0) eta2 = eta2 + 1.d0


  ! setup relative coordinates xiL, xiR for divertor legs
  alpha = 1.d0 + eta1 * (alphaR(jx1) - 1.d0)
  xiR   = alpha * S(jx1)%M3%l / l
  DthetaR = eta1 * Dtheta_sampling
  alpha = 1.d0 + eta2 * (alphaL(jx2) - 1.d0)
  xiL   = 1.d0 - alpha * S(jx2)%M4%l / l
  DthetaL = eta2 * Dtheta_sampling


  ! split flux surface in main part and divertor segments
  call F%split3(xiR, xiL, CR, C0%t_curve, CL)


  ! automatic calculation of Dtheta
  if (Dtheta_sampling < 0.d0) then
     x       = C0%x(C0%n_seg,:)
     thetaL  = atan2(x(2)-Pmag(2), x(1)-Pmag(1))
     DthetaL = 2.d0 * (Xp(jx2)%theta - thetaL)

     x       = C0%x(0,:)
     thetaR  = atan2(x(2)-Pmag(2), x(1)-Pmag(1))
     DthetaR = 2.d0 * (thetaR - Xp(jx2)%theta)
  endif


  ! set minimum Dtheta to Dtheta_separatrix
  DthetaR = max(DthetaR, Dtheta_separatrix)
  DthetaL = max(DthetaL, Dtheta_separatrix)


  ! set up segments for sampling
  call CR%setup_length_sampling()
  call C0%setup_sampling(Xp(jx1)%X, Xp(jx2)%X, Pmag, DthetaR, DthetaL)
  call CL%setup_length_sampling()

  !call CL%plot(filename='CL.plt', append=.true.)
  !call C0%plot(filename='C0.plt', append=.true.)
  !call CR%plot(filename='CR.plt', append=.true.)

  end subroutine divide_SOL3
  !=====================================================================




  !=====================================================================
  ! Calculate intersection between divertor leg (of flux surface) and
  ! target surface / guiding surface.
  ! Return relative coordinate (eta) of intersection on divertor leg.
  !=====================================================================
  subroutine divertor_leg_interface(C_leg, C_cut, eta)
  use curve2D
  type(t_curve), intent(in) :: C_leg, C_cut
  real(real64), intent(out) :: eta

  real(real64) :: x(2)


  if (.not.C_leg%intersect_curve(C_cut, x, eta)) then
     write (6, *) 'error: could not find intersection between divertor leg and guiding surface!'
     write (6, *) 'check divertor_leg.plt vs C_cutL.plt/C_cutR.plt!'
     call C_leg%plot(filename='divertor_leg.plt')
     stop
  endif

  end subroutine divertor_leg_interface
  !=====================================================================
  subroutine divertor_leg_discretization(C, eta, np, D)
  use mesh_spacing
  type(t_curve), intent(in)  :: C
  real(real64),  intent(in)  :: eta
  integer,       intent(in)  :: np
  real(real64),  intent(out) :: D(0:np,2)

  type(t_spacing) :: Sguide
  real(real64)    :: xi, xi0, x(2)
  integer         :: j


  call divertor_leg_interface(C, C_guide, xi0)
  call Sguide%init_spline_X1(eta, xi0)
  do j=0,np
     xi     = Sguide%node(j,np)
     call C%sample_at(xi, x)
     D(j,:) = x
  enddo

  end subroutine divertor_leg_discretization
  !=====================================================================


  !=============================================================================
  subroutine make_interface_core(F, ix1, ix2, S, np, D)
  use flux_surface_2D
  use mesh_spacing
  use equilibrium
  use fieldline_grid, only: Dtheta_separatrix
  use math
  type(t_flux_surface_2D), intent(in)  :: F
  integer,                 intent(in)  :: ix1, ix2, np
  type(t_spacing),         intent(in)  :: S
  real(real64),            intent(out) :: D(0:np,2)

  real(real64) :: rho, x(2)
  integer      :: j


  call F%setup_sampling(Xp(ix1)%X, Xp(ix2)%X, Pmag, Dtheta_separatrix, Dtheta_separatrix)

  do j=0,np
     rho = S%node(j,np)

     call F%sample_at(rho, x)
     D(j,:) = x
  enddo

  end subroutine make_interface_core
  !=============================================================================



  !=============================================================================
  ! unperturbed FLUX SURFACES (high pressure region)
  !=============================================================================
  subroutine make_flux_surfaces_HPR(M, nr, np, ir1, ir2, rpath, Sr, Sp)
  use run_control, only: Debug
  use xpaths
  use mesh_spacing
  use flux_surface_2D

  real(real64), dimension(:,:,:), pointer, intent(inout) :: M
  integer, intent(in) :: nr, np, ir1, ir2
  type(t_xpath),   intent(in) :: rpath
  type(t_spacing), intent(in) :: Sr, Sp

  type(t_flux_surface_2D) :: F
  real(real64) :: eta, x(2), xi
  integer      :: i, j


  if (ir2 .ge. ir1) then
     write (6, 1010) ir2, ir1
     write (6, 1011) rpath%length()
  endif
  do i=ir2, ir1, -1
     write (6, *) i
     eta = 1.d0 - Sr%node(i-1,nr-1)
     call rpath%sample_at(eta, x)
     if (Debug) write (iud, *) x
     call F%generate_closed(x, RIGHT_HANDED)
     call F%setup_angular_sampling(Pmag)

     do j=0,np
        xi = Sp%node(j,np)
        call F%sample_at(xi, x)
        M(i,j,:) = x
     enddo
  enddo

 1010 format (8x,'generating high pressure region: ', i0, ' -> ', i0)
 1011 format (8x,'d_width = ',f8.3)
  end subroutine make_flux_surfaces_HPR
  !=============================================================================

  !=============================================================================
  ! inner boundaries and interpolated surfaces (2 -> 1+n_interpolate) (high pressure region)
  !=============================================================================
  subroutine make_interpolated_surfaces(M, nr, np, ir1, ir2, Sr, Sp, C)
  use run_control, only: Debug
  use mesh_spacing
  use flux_surface_2D
  use fieldline_grid, only: nr_perturbed

  real(real64), dimension(:,:,:), pointer, intent(inout) :: M
  integer, intent(in) :: nr, np, ir1, ir2
  type(t_spacing), intent(in) :: Sr, Sp
  type(t_curve), intent(in)   :: C(0:1)

  real(real64) :: eta, xi, x(2), x1(2), x2(2)
  integer      :: i, j


  write (6, 1001) ir1+1, ir2-1
  do j=0,np
     xi = Sp%node(j,np)
     ! innermost surfaces
     do i=0,nr_perturbed-1
        call C(i)%sample_at(xi, x)
        M(i, j, :) = x
     enddo

     ! interpolated surfaces
     x1 = M(ir1, j, :)
     x2 = M(ir2, j, :)
     do i=ir1+1,ir2-1
        eta = Sr%node(i-ir1, nr-1) / Sr%node(ir2-ir1, nr-1)

        M(i,j,:) = x1 + eta * (x2-x1)
     enddo
  enddo

 1001 format (8x,'interpolating from inner boundary to 1st unperturbed flux surface: ', &
              i0, ' -> ', i0)
  end subroutine make_interpolated_surfaces
  !=============================================================================

  !=============================================================================
  ! orthogonal flux surface grid
  !=============================================================================
  subroutine make_ortho_grid(M, nr, np, ir0, ir1, ir2, ip0, ip1, ip2, ipx, ix, rpath, Sr, Sp, periodic)
  use run_control, only: Debug
  use xpaths
  use mesh_spacing
  use equilibrium, only: get_PsiN

  real(real64), dimension(:,:,:), pointer, intent(inout) :: M
  integer, intent(in) :: nr, np, ir0, ir1, ir2, ip0, ip1, ip2, ipx, ix
  type(t_xpath),   intent(in) :: rpath
  type(t_spacing), intent(in) :: Sr, Sp
  logical,         intent(in), optional :: periodic

  type(t_xpath) :: R
  real(real64)  :: x(2), PsiN(0:nr), eta
  integer       :: ir, ip, direction


  ! set direction
  direction = DESCENT
  if (ir1 > ir0) direction = ASCENT


  ! set up nodes along rpath and radial coordinate PsiN
  if (ir2 .ge. ir1) then
     write (6, 1010) ir2, ir1, ip1, ip2
     write (6, 1011) rpath%length()
  endif
  do ir=ir1,ir2
     eta = 1.d0 - Sr%node(ir-1,nr-1)
     call rpath%sample_at(eta, x)
     PsiN(ir) = get_PsiN(x)
     M(ir,ip0,:) = x
  enddo


  ! set up nodes in poloidal range ip1->ip2
  do ip=ip1,ip2
     call R%generate(M(ir0,ip,:), direction, LIMIT_PSIN, rpath%PsiN(2), sampling=SAMPLE_PSIN)

     do ir=ir1,ir2
        call R%sample_at_PsiN(PsiN(ir), x)
        M(ir,ip,:) = x
     enddo
  enddo


  ! set up nodes at poloidal index ix
  if (ipx >= 0) then
     call R%generateX(ix, direction, LIMIT_PSIN, rpath%PsiN(2), sampling=SAMPLE_PSIN)
     do ir=ir1,ir2
        call R%sample_at_PsiN(PsiN(ir), x)
        M(ir,ipx,:) = x
     enddo
  endif


  ! set up periodic boundaries
  if (present(periodic)  .and.  periodic) then
     M(:,np,:) = M(:,0,:)
  endif

 1010 format (8x,'generating high pressure region: (',i0,' -> ',i0,'), (',i0,' -> ',i0,')')
 1011 format (8x,'d_width = ',f8.3)
  end subroutine make_ortho_grid
  !=============================================================================

  !=============================================================================
  ! inner boundaries and interpolated surfaces (2 -> 1+n_interpolate) (high pressure region)
  !=============================================================================
  subroutine make_interpolated_surfaces_ortho(M, nr, np, ir2, Sr, Sp, C, PsiN1)
  use run_control, only: Debug
  use mesh_spacing
  use flux_surface_2D
  use math
  use equilibrium
  use xpaths
  use fieldline_grid, only: nr_perturbed

  real(real64), dimension(:,:,:), pointer, intent(inout) :: M
  integer, intent(in) :: nr, np, ir2
  type(t_spacing), intent(in) :: Sr, Sp
  type(t_curve), intent(in)   :: C(0:1)
  real(real64),  intent(in)   :: PsiN1

  type(t_xpath) :: R
  real(real64)  :: eta, xi, x(2), x1(2), x2(2), theta
  integer       :: i, ir1, j


  ir1 = nr_perturbed-1
  write (6, 1001) ir1+1, ir2-1
  do j=0,np
     write (6, *) j
     x = M(ir2,j,:)
     call R%generate(x, DESCENT_CORE, LIMIT_CURVE, PsiN1, C_limit=C(ir1), sampling=SAMPLE_LENGTH)

     ! interpolated surfaces
     do i=ir1,ir2-1
        eta = 1.d0 - Sr%node(i-1, nr-1) / Sr%node(ir2-1, nr-1)
        call R%sample_at(eta, x)

        M(i,j,:) = x
     enddo

     ! innermost surface
     if (ir1 == 1) then
        theta = atan2(M(ir1,j,2)-Pmag(2), M(ir1,j,1)-Pmag(1)) - Xp(1)%theta
        xi    = theta / pi2; if (xi < 0.d0) xi = xi + 1.d0
        call C(0)%sample_at(xi, x)
        M(0, j, :) = x
     endif
  enddo

 1001 format (8x,'interpolating from inner boundary to 1st unperturbed flux surface: ', &
              i0, ' -> ', i0)
  end subroutine make_interpolated_surfaces_ortho
  !=============================================================================


  !=============================================================================
  ! scrape-off layer
  !=============================================================================
  subroutine make_flux_surfaces_SOL(M, nr, npL, np0, npR, ir1, ir2, rpath, ix1, ix2, Sr, Sp)
  use run_control, only: Debug
  use xpaths
  use mesh_spacing
  use flux_surface_2D
  use fieldline_grid, only: etaR, etaL

  real(real64), dimension(:,:,:), pointer, intent(inout) :: M
  integer, intent(in) :: nr, npL, np0, npR, ir1, ir2, ix1, ix2
  type(t_xpath),   intent(in) :: rpath
  type(t_spacing), intent(in) :: Sr, Sp

  type(t_flux_surface_2D) :: F, C0
  type(t_curve)           :: CL, CR
  type(t_spacing) :: Sdr, Sdl
  real(real64)  :: x0(2), x(2), eta, xi, xiR, xiL, DL(0:npL, 2), DR(0:npR, 2)
  integer       :: i, j


  write (6, 1020) ir1, ir2
  write (6, 1021) rpath%length()
  do i=ir1,ir2
     write (6, *) i
     eta = Sr%node(i,nr)
     call rpath%sample_at(eta, x0)
     if (Debug) write (iud, *) x0
     call F%generate_open(x0, C_cutL, C_cutR)
     call divide_SOL3(F, eta, CL, C0, CR, ix1, ix2)

     ! right divertor leg
     call divertor_leg_discretization(CR, 1.d0-etaR(1), npR, DR)
     M(i,0:npR,:) = DR


     ! main SOL
     do j=0,np0
        xi = Sp%node(j,np0)
        call C0%sample_at(xi, x)
        M(i,npR+j,:) = x
     enddo

     ! left divertor leg
     call divertor_leg_discretization(CL, etaL(1), npL, DL)
     M(i,npR+np0:npR+np0+npL,:) = DL
  enddo

 1020 format (8x,'generating scrape-off layer: ',i0,' -> ', i0)
 1021 format (8x,'d_SOL = ',f8.3)
  end subroutine make_flux_surfaces_SOL
  !=============================================================================

  !=============================================================================
  ! private flux region
  !=============================================================================
  subroutine make_flux_surfaces_PFR(M, nr, npL, npR, ir1, ir2, rpath, Sr, Sp)
  use run_control, only: Debug
  use xpaths
  use mesh_spacing
  use flux_surface_2D
  use fieldline_grid, only: etaR, etaL

  real(real64), dimension(:,:,:), pointer, intent(inout) :: M
  integer, intent(in) :: nr, npL, npR, ir1, ir2
  type(t_xpath),   intent(in) :: rpath
  type(t_spacing), intent(in) :: Sr, Sp

  type(t_flux_surface_2D) :: F
  type(t_spacing) :: Sdr, Sdl
  real(real64) :: eta, xi, xiL, xiR, x0(2), x(2), DL(0:npL, 2), DR(0:npR, 2)
  integer      :: i, j

  write (6, 1030) nr-1
  write (6, 1031) rpath%length()
  do i=0,nr-1
     write (6, *) i
     eta = Sr%node(i,nr)
     call rpath%sample_at(1.d0 - eta, x0)
     if (Debug) write (iud, *) x0

     ! right divertor leg
     call F%generate(x0, -1, AltSurf=C_cutR, sampling=DISTANCE)
     call divertor_leg_discretization(F%t_curve, 1.d0-etaR(1), npR, DR)
     M(i,0:npR,:) = DR


     ! left divertor leg
     call F%generate(x0,  1, AltSurf=C_cutL, sampling=DISTANCE)
     call divertor_leg_discretization(F%t_curve, etaL(1), npL, DL)
     M(i,npR:npR+npL,:) = DL
  enddo

 1030 format (8x,'generating private flux region: 0 -> ', i0)
 1031 format (8x,'d_PFR = ',f8.3)
  end subroutine make_flux_surfaces_PFR
  !=============================================================================



!=======================================================================
! CLOSE_GRID_DOMAIN
!
! Extend the divertor legs to first/last toroidal grid slice in order to
! close the simulation domain (required for module MAPPING in EMC3).
!
! !!! These cells MUST be behind the wall !!!
!=======================================================================
  subroutine close_grid_domain(iz)
  use equilibrium

  integer, intent(in) :: iz

  real(real64), parameter :: d_extend = 2.d0


  integer      :: iside, j, j0(-1:1), k0(-1:1)
  logical      :: debug = .false.


  P_SURF_PL_TRANS_RANGE(1,iz) = 1
  P_SURF_PL_TRANS_RANGE(2,iz) = ZON_POLO(iz)-1


  j0(-1) = 0
  j0( 1) = SRF_POLO(iz)-1
  k0(-1) = SRF_TORO(iz)-1
  k0( 1) = 0
  do iside=-1,1,2
     ! poloidal index of reference surface
     j = j0(int(iside * Ip_sign * Bt_sign))

     ! extend poloidal surface
     call extend_poloidal_surface(d_extend)

     ! distribute nodes in toroidal direction
     call distribute_nodes()
  enddo


  contains
  !.....................................................................

  !.....................................................................
  subroutine distribute_nodes()

  real(real64) :: R1, Z1
  integer      :: i, ig, ig0, k

  write (6, *) 'adjusting nodes to surface (iz, it, ip) = ', iz, k0(iside), j
  do i=0,SRF_RADI(iz)-1
     ig0 = i + (j + k0(iside)*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
     R1 = RG(ig0)
     Z1 = ZG(ig0)

     do k=0,SRF_TORO(iz)-1
        ig = i + (j + k*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
        RG(ig) = R1
        ZG(ig) = Z1
     enddo
     if (debug) write (96+iside, *) R1, Z1
  enddo
  end subroutine distribute_nodes
  !.....................................................................


  !.....................................................................
  subroutine extend_poloidal_surface(d_extend)
  ! => straight connection between first and last radial node
  ! this should fix some problems for module MAPPING in EMC3
  ! TODO: find optimal/minimal value for d_extend

  real(real64), intent(in) :: d_extend


  real(real64), dimension(:), allocatable :: w
  real(real64) :: R1, Z1, R2, Z2
  integer :: i, k, ir1, ir2, ig, ig0


  ir1 = 0
  ir2 = SRF_RADI(iz)-1
  !ir1 = R_SURF_PL_TRANS_RANGE(1,iz)
  !ir2 = R_SURF_PL_TRANS_RANGE(2,iz)


  allocate (w(ir1:ir2))
  w = 0.d0
  k = k0(iside)


  ! set up weight to maintain relative cell widths
  do i=ir1+1,ir2
     ig   = i + (j + k*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
     w(i) = w(i-1) + sqrt((ZG(ig)-ZG(ig-1))**2 + (RG(ig)-RG(ig-1))**2)
     !write (98, *) RG(ig-1), ZG(ig-1)
     !write (98, *) RG(ig), ZG(ig)
  enddo
  w = w / w(ir2)


  ! move first and last node in poloidal direction
  ig = ir1 + (j + k*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
  if (debug) write (96, *) RG(ig), ZG(ig)
  call move_node(iz,ir1,j,k,0,iside,0,d_extend)
  R1 = RG(ig); Z1 = ZG(ig)
  if (debug) then
     write (96, *) RG(ig), ZG(ig)
     write (96, *)
  endif


  ig = ir2 + (j + k*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
  if (debug) write (96, *) RG(ig), ZG(ig)
  call move_node(iz,ir2,j,k,0,iside,0,d_extend)
  R2 = RG(ig); Z2 = ZG(ig)
  if (debug) then
     write (96, *) RG(ig), ZG(ig)
     write (96, *)
  endif


  ! interpolate nodes in between (and maintain relative cell width)
  do i=ir1,ir2
     ig = i + (j + k*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
     RG(ig) = R1 + w(i) * (R2-R1)
     ZG(ig) = Z1 + w(i) * (Z2-Z1)
  enddo


  ! cleanup
  deallocate(w)

  end subroutine extend_poloidal_surface
  !.....................................................................
  end subroutine close_grid_domain
!=======================================================================
!=======================================================================
  subroutine move_node(iz,ir,ip,it,irdir,ipdir,itdir,d)
  use emc3_grid
  integer,      intent(in) :: iz, ir, ip, it, irdir, ipdir, itdir
  real(real64), intent(in) :: d

  real(real64) :: v(2)
  integer      :: ig, ig1


  ig     = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
  ig1    = ir+irdir + (ip+ipdir + (it+itdir)*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)

  v(1)   = RG(ig) - RG(ig1)
  v(2)   = ZG(ig) - ZG(ig1)
  v      = v / sqrt(sum(v**2))

  RG(ig) = RG(ig) + d*v(1)
  ZG(ig) = ZG(ig) + d*v(2)

  end subroutine move_node
!=======================================================================

end module divertor
