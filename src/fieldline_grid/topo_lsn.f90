!===============================================================================
! Lower Single Null configuration: block-structured decomposition with zones for
! high pressure region (HPR), scrape-off layer (SOL) and private flux region (PFR)
!===============================================================================
module topo_lsn
  use iso_fortran_env
  use grid
  use equilibrium
  use separatrix
  use curve2D
  use fieldline_grid, unused => TOPO_LSN
  use inner_boundary
  implicit none
  private

  integer, parameter :: DEFAULT = 0
  integer, parameter :: iud     = 72




  ! coordinates of X-point and magnetic axis
  real(real64) :: Px(2), Pmag(2)

  ! discretization method
  integer :: method = DEFAULT


  ! base grid in high pressure region (HPR), scrape-off layer (SOL) and private flux region (PFR)
  type(t_grid), dimension(:), allocatable :: G_HPR, G_SOL, G_PFR ! (0:blocks-1)

  ! magnetic separatrix
  type(t_separatrix), public :: S
  type(t_curve)      :: S0

  ! guiding_surface
  type(t_curve), public      :: C_guide, C_cutL, C_cutR


  public :: &
     setup_topo_lsn, &
     make_base_grids_lsn, &
     post_process_grid_lsn, &
     make_flux_surfaces_SOL, &
     make_flux_surfaces_PFR, &
     make_interpolated_surfaces

  contains
  !=====================================================================



  !=====================================================================
  subroutine setup_topo_lsn()
  use emc3_grid, only: NZONET
  implicit none

  integer :: iz, iz0, iz1, iz2, ib


  ! 0. setup number of zones for lower single null topology
  layers = 3
  NZONET = blocks * layers


  write (6, 1000)
  write (6, 1001)
  do ib=0,blocks-1
     ! 1. setup resolution for each zone
     ! set derived parameters
     Block(ib)%np(1) = Block(ib)%npR(1) + Block(ib)%np(0) + Block(ib)%npL(1)
     Block(ib)%np(2) = Block(ib)%npR(1)                   + Block(ib)%npL(1)


     ! high pressure region (HPR)
     iz0 = 3*ib
     Zone(iz0)%nr = Block(ib)%nr(0) + nr_EIRENE_core
     Zone(iz0)%np = Block(ib)%np(0)

     ! scrape-off layer (SOL)
     iz1 = iz0 + 1
     Zone(iz1)%nr = Block(ib)%nr(1) + nr_EIRENE_vac
     Zone(iz1)%np = Block(ib)%np(1)

     ! private flux region (PFR)
     iz2 = iz1 + 1
     Zone(iz2)%nr = Block(ib)%nr(2) + nr_EIRENE_vac
     Zone(iz2)%np = Block(ib)%np(2)

     do iz=iz0,iz2
        ! setup toroidal discretization
        call Zone(iz)%setup_default_toroidal_discretization(ib)
     enddo


     ! 2. setup boundaries and connectivity between zones
     do iz=iz0,iz2
        call Zone(iz)%setup_default_boundaries()
     enddo

     ! 2.a high pressure region (HPR)
     Zone(iz0)%isfr(2) = SF_MAPPING
     Zone(iz0)%r_surf_pl_trans_range(2) = Zone(iz0)%nr

     ! 2.b scrape-off layer (SOL)
     Zone(iz1)%isfr(1) = SF_MAPPING
     Zone(iz1)%isfp(1) = SF_VACUUM
     Zone(iz1)%isfp(2) = SF_VACUUM
     Zone(iz1)%r_surf_pl_trans_range(1) = 0
     Zone(iz1)%d_N0    = d_N0(1)

     ! 2.c private flux region (PFR)
     Zone(iz2)%isfr(1) = SF_VACUUM
     Zone(iz2)%isfr(2) = SF_MAPPING
     Zone(iz2)%isfp(1) = SF_VACUUM
     Zone(iz2)%isfp(2) = SF_VACUUM
     Zone(iz2)%r_surf_pl_trans_range(1) = nr_EIRENE_vac
     Zone(iz2)%r_surf_pl_trans_range(2) = Zone(iz2)%nr
     Zone(iz2)%d_N0    = d_N0(2)


     write (6, 1002) ib, Zone(iz0)%nr, Zone(iz0)%np, Zone(iz0)%nt, &
                         Zone(iz1)%nr, Zone(iz1)%np, Zone(iz1)%nt, &
                         Zone(iz2)%nr, Zone(iz2)%np, Zone(iz2)%nt
  enddo
 1000 format(8x,'Grid resolution is (radial x poloidal x toroidal):')
 1001 format(8x,'block #, high pressure region, scrape-off layer, private flux region')
 1002 format(12x,i3,3(3x,'(',i0,' x ',i0,' x ',i0,')'))

  end subroutine setup_topo_lsn
  !=====================================================================



  !=====================================================================
  subroutine setup_domain
  use boundary
  use run_control, only: Debug
  use math
  use flux_surface_2D, only: RIGHT_HANDED
  use equilibrium
  use inner_boundary

  real(real64) :: tmp(3), dx(2)


  ! 1.a setup guiding surface for divertor legs (C_guide) ------------------
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

  ! 1.b setup extended guiding surfaces for divertor leg discretization ----
  ! C_cutL, C_cutR
  call C_cutL%copy(C_guide)
  call C_cutL%left_hand_shift(d_cutL(1))
  call C_cutR%copy(C_guide)
  call C_cutR%left_hand_shift(d_cutR(1))
  if (Debug) then
     call C_cutL%plot(filename='C_cutL.plt')
     call C_cutR%plot(filename='C_cutR.plt')
  endif


  ! 2.a setup magnetic axis (Pmag) --------------------------------------
  tmp = get_magnetic_axis(0.d0); Pmag = tmp(1:2)
  Magnetic_Axis%X = Pmag

  ! 2.b setup X-point (Px, theta0) --------------------------------------
  Px = Xp(1)%load()
  write (6, 2000) Px
  write (6, 2001) Xp(1)%theta

  ! 2.c separatrix (S, S0) ---------------------------------------------
  call S%generate(1, RIGHT_HANDED, pi/2.d0, C_cutL, C_cutR)
  call S%plot('S', parts=.true.)

  ! connect core segments of separatrix
  S0 = connect(S%M1%t_curve, S%M2%t_curve)
  call S0%plot(filename='S0.plt')
  call S0%setup_angular_sampling(Pmag)
  call S%M3%setup_length_sampling()
  call S%M4%setup_length_sampling()
 2000 format(8x,'found magnetic X-point at: ',2f10.4)
 2001 format(11x,'-> poloidal angle [deg]: ',f10.4)


  ! 3. inner boundaries for EMC3 grid
  call load_inner_boundaries(Xp(1)%theta)

  ! 4. setup paths for discretization in radial direction
  ! 4.1 SOL
  dx(1) = Px(2) - Pmag(2)
  dx(2) = Pmag(1) - Px(1)
  dx    = dx / sqrt(sum(dx**2)) * d_SOL(1)
  call rpath(1)%setup_linear(Px, dx)
  ! 4.2 PFR
  dx    = Px - Pmag
  dx    = dx / sqrt(sum(dx**2)) * d_PFR(1)
  call rpath(2)%setup_linear(Px, dx)

  end subroutine setup_domain
  !=====================================================================



  !=====================================================================
  subroutine divide_SOL(F, eta, CL, C0, CR)
  use flux_surface_2D
  use math
  type(t_flux_surface_2D), intent(in)  :: F
  real(real64),            intent(in)  :: eta
  type(t_curve),           intent(out) :: CL, CR
  type(t_flux_surface_2D), intent(out) :: C0

  real(real64) :: l, alpha, xiR, xiL


  l = F%length()

  alpha = 1.d0 + eta * (alphaR(1) - 1.d0)
  xiR   = alpha * S%M3%l / l
  alpha = 1.d0 + eta * (alphaL(1) - 1.d0)
  xiL   = 1.d0 - alpha * S%M4%l / l
  call F%split3(xiR, xiL, CR, C0%t_curve, CL)
  call CR%setup_length_sampling()
  call C0%setup_sampling(Xp(1)%X, Xp(1)%X, Magnetic_Axis%X, eta, eta, pi2, Dtheta_sampling)
  call CL%setup_length_sampling()

  end subroutine divide_SOL
  !=====================================================================



  !=====================================================================
  subroutine make_base_grids_lsn
  use run_control, only: Debug
  use math
  use inner_boundary
  use flux_surface_2D
  use mesh_spacing
  use divertor

  integer, parameter      :: iu = 72

  type(t_flux_surface_2D) :: FS, FSL, FSR, C0
  type(t_curve)           :: CL, CR
  type(t_spacing)         :: Sl, Sr

  real(real64), dimension(:,:,:), pointer :: M_HPR, M_SOL, M_PFR

  real(real64) :: xi, eta, phi, x(2), x0(2), x1(2), x2(2), d_HPR(2), dx(2)
  integer :: i, j, iz, iz0, iz1, iz2, nr0, nr1, nr2, np0, np1, np1l, np1r, np2

  logical :: generate_flux_surfaces_HPR
  logical :: generate_flux_surfaces_SOL
  logical :: generate_flux_surfaces_PFR
  real(real64) :: xiL, xiR
  integer :: iblock


  write (6, 1000)
  if (Debug) then
     open  (iu, file='base_grid_debug.txt')
  endif
  !.....................................................................
  ! 0. initialize geometry
  call setup_domain()
  !.....................................................................


  !.....................................................................
  ! 1. check input
  if (n_interpolate < 0) then
     write (6, *) 'error: n_interpolate must not be negative!'; stop
  endif
  if (n_interpolate > nr(0)-2) then
     write (6, *) 'error: n_interpolate > nr0 - 2!'; stop
  endif
  !.....................................................................


  !.....................................................................
  ! 2. setup working arrays for base grid
  allocate (G_HPR(0:blocks-1), G_SOL(0:blocks-1), G_PFR(0:blocks-1))
  !.....................................................................



  do iblock=0,blocks-1
     write (6, 1001) iblock

     ! set zone indices
     iz0 = iblock*3
     iz1 = iz0 + 1
     iz2 = iz0 + 2

     ! set local variables for resolution
     nr0 = Block(iblock)%nr(0); np0 = Block(iblock)%np(0)
     nr1 = Block(iblock)%nr(1); np1 = Block(iblock)%np(1)
     np1l = Block(iblock)%npL(1); np1r = Block(iblock)%npR(1)
     nr2 = Block(iblock)%nr(2); np2 = Block(iblock)%np(2)
     ! check if radial-poloidal resolution is different from last block
     generate_flux_surfaces_HPR = .true.
     generate_flux_surfaces_SOL = .true.
     generate_flux_surfaces_PFR = .true.
     if (iblock > 0) then
        ! copy unperturbed flux surface discretization if resolution is the same
        if (nr0 == Block(iblock-1)%nr(0)  .and.  np0 == Block(iblock-1)%np(0)) then
           generate_flux_surfaces_HPR = .false.
        endif

        if (np1l == Block(iblock-1)%npL(1) .and. np1r == Block(iblock-1)%npR(1)) then
           if (nr1 == Block(iblock-1)%nr(1)) generate_flux_surfaces_SOL = .false.
           if (nr2 == Block(iblock-1)%nr(2)) generate_flux_surfaces_PFR = .false.
        endif
     endif



     ! cell spacings
     do i=0,2
        iz = iz0 + i
        call Zone(iz)%Sr%init(radial_spacing(i))
        call Zone(iz)%Sp%init(poloidal_spacing(i))
     enddo


     ! initialize base grids in present block
     phi = Block(iblock)%phi_base / 180.d0 * pi
     call G_HPR(iblock)%new(CYLINDRICAL, MESH_2D, 3, nr0+1, np0+1, fixed_coord_value=phi)
     call G_SOL(iblock)%new(CYLINDRICAL, MESH_2D, 3, nr1+1, np1+1, fixed_coord_value=phi)
     call G_PFR(iblock)%new(CYLINDRICAL, MESH_2D, 3, nr2+1, np2+1, fixed_coord_value=phi)
     M_HPR => G_HPR(iblock)%mesh
     M_SOL => G_SOL(iblock)%mesh
     M_PFR => G_PFR(iblock)%mesh


     ! start grid generation
     ! 1. unperturbed separatrix
     call make_separatrix()

     ! 2. unperturbed FLUX SURFACES
     ! 2.a high pressure region (HPR)
     if (generate_flux_surfaces_HPR) then
        call make_flux_surfaces_HPR(M_HPR, nr0, np0, 2+n_interpolate, nr0-1, iz0)
     else
        G_HPR(iblock)%mesh = G_HPR(iblock-1)%mesh
     endif

     ! 2.b scrape-off layer (SOL)
     if (generate_flux_surfaces_SOL) then
        call make_flux_surfaces_SOL(M_SOL, nr1, np1l, np0, np1r, 1, nr1, Zone(iz1)%Sr, Zone(iz0)%Sp)
     else
        G_SOL(iblock)%mesh = G_SOL(iblock-1)%mesh
     endif

     ! 2.c private flux region (PFR)
     if (generate_flux_surfaces_PFR) then
        call make_flux_surfaces_PFR(M_PFR, nr2, np1l, np1r, 1, nr2, rpath(2), Zone(iz2)%Sr, Zone(iz2)%Sp)
     else
        G_PFR(iblock)%mesh = G_PFR(iblock-1)%mesh
     endif

     ! 3. interpolated surfaces
     call make_interpolated_surfaces(M_HPR, nr0, np0, 1, 2+n_interpolate, Zone(iz0)%Sr, Zone(iz0)%Sp, C_in(iblock,:))


     ! output
     call write_base_grid(G_HPR(iblock), iz0)
     call write_base_grid(G_SOL(iblock), iz1)
     call write_base_grid(G_PFR(iblock), iz2)
     write (6, 1002) iblock
  enddo
  if (Debug) close (iu)

 1000 format(//3x,'- Setup for base grids:')
 1001 format(//1x,'Start generation of base grids for block ',i0,' ',32('.'))
 1002 format(1x,'Finished generation of base grids for block ',i0,' ',32('.'),//)
  contains
  !.....................................................................
  subroutine make_separatrix()

  ! 1. discretization of main part of separatrix
  do j=0,np0
     xi = Zone(iz0)%Sp%node(j,np0)

     call S0%sample_at(xi, x)
     M_HPR(nr0,      j, :) = x
     M_SOL(  0, np1r+j, :) = x
  enddo

  ! 2. discretization of right separatrix leg
  call divertor_leg_interface(S%M3%t_curve, C_guide, xiR)
  call Sr%init_spline_X1(etaR(1), 1.d0-xiR)
  do j=0,np1r
     xi = 1.d0 - Sr%node(np1r-j,np1r)
     call S%M3%sample_at(xi, x)
     M_SOL(  0,j,:) = x
     M_PFR(nr2,j,:) = x
  enddo

  ! 3. discretization of left separatrix leg
  call divertor_leg_interface(S%M4%t_curve, C_guide, xiL)
  call Sl%init_spline_X1(etaL(1), xiL)
  do j=1,np1l
     xi = Sl%node(j,np1l)
     call S%M4%sample_at(xi, x)
     M_SOL(  0,np1r + np0 + j,:) = x
     M_PFR(nr2,np1r       + j,:) = x
  enddo
  end subroutine make_separatrix
  !.....................................................................
  end subroutine make_base_grids_lsn
  !=============================================================================


  !=============================================================================
  ! unperturbed FLUX SURFACES (high pressure region)
  !=============================================================================
  subroutine make_flux_surfaces_HPR(M, nr, np, ir1, ir2, iz0)
  use run_control, only: Debug
  use flux_surface_2D
  !use mesh_spacing
  !use divertor

  real(real64), dimension(:,:,:), pointer, intent(inout) :: M
  integer, intent(in) :: nr, np, ir1, ir2, iz0

  type(t_flux_surface_2D) :: F
  real(real64) :: d_HPR(2)
  real(real64) :: eta, x(2), xi
  integer      :: i, j


  ! 1. get radial width at poloidal angle of X-point
  d_HPR = get_d_HPR(Px, Pmag)

  ! 2. generate flux surfaces
  if (ir2 .ge. ir1) then
     write (6, 1010) ir2, ir1
     write (6, 1011) Px, Px + d_HPR
  endif
  do i=ir2, ir1, -1
     write (6, *) i
     eta = 1.d0 - Zone(iz0)%Sr%node(i-1,nr-1)

     x = Px + eta * d_HPR
     if (Debug) write (iud, *) x
     call F%generate_closed(x, RIGHT_HANDED)
     call F%setup_angular_sampling(Pmag)

     do j=0,np
        xi = Zone(iz0)%Sp%node(j,np)
        call F%sample_at(xi, x)
        M(i,j,:) = x
     enddo
  enddo

 1010 format (8x,'generating high pressure region: ', i0, ' -> ', i0)
 1011 format (8x,'from (',f8.3,', ',f8.3,') to (',f8.3,', ',f8.3,')')
  end subroutine make_flux_surfaces_HPR
  !=============================================================================

  !=============================================================================
  ! inner boundaries and interpolated surfaces (2 -> 1+n_interpolate) (high pressure region)
  !=============================================================================
  subroutine make_interpolated_surfaces(M, nr, np, ir1, ir2, Sr, Sp, C)
  use run_control, only: Debug
  use flux_surface_2D

  real(real64), dimension(:,:,:), pointer, intent(inout) :: M
  integer, intent(in) :: nr, np, ir1, ir2
  type(t_spacing), intent(in) :: Sr, Sp
  type(t_curve), intent(in)   :: C(0:1)

  real(real64) :: eta, xi, x(2), x1(2), x2(2)
  integer      :: i, j


  write (6, 1001) 2, 1+n_interpolate
  do j=0,np
     xi = Sp%node(j,np)
     ! innermost surfaces
     do i=0,1
        call C(i)%sample_at(xi, x)
        M(i, j, :) = x
     enddo

     ! interpolated surfaces
     x1 = M(ir1, j, :)
     x2 = M(ir2, j, :)
     do i=ir1+1,ir2-1
        eta = Sr%node(i-1, nr-1) / Sr%node(ir2-1, nr-1)

        M(i,j,:) = x1 + eta * (x2-x1)
     enddo
  enddo

 1001 format (8x,'interpolating from inner boundary to 1st unperturbed flux surface: ', &
              i0, ' -> ', i0)
  end subroutine make_interpolated_surfaces
  !=============================================================================

  !=============================================================================
  ! scrape-off layer
  !=============================================================================
  subroutine make_flux_surfaces_SOL(M, nr, npL, np0, npR, ir1, ir2, Sr, Sp)
  use run_control, only: Debug
  use flux_surface_2D
  use divertor

  real(real64), dimension(:,:,:), pointer, intent(inout) :: M
  integer, intent(in) :: nr, npL, np0, npR, ir1, ir2
  type(t_spacing), intent(in) :: Sr, Sp

  type(t_flux_surface_2D) :: F, C0
  type(t_curve)           :: CL, CR
  type(t_spacing) :: Sdr, Sdl
  real(real64)  :: x0(2), x(2), eta, xi, xiR, xiL
  integer       :: i, j


  write (6, 1020) nr
  write (6, 1021) d_SOL(1)
  do i=ir1,ir2
     write (6, *) i
     eta = Sr%node(i,nr)
     call rpath(1)%sample_at(eta, x0)
     if (Debug) write (iud, *) x0
     call F%generate_open(x0, C_cutL, C_cutR)
     call divide_SOL(F, eta, CL, C0, CR)

     ! right divertor leg
     call divertor_leg_interface(CR, C_guide, xiR)
     call Sdr%init_spline_X1(etaR(1), 1.d0-xiR)
     do j=0,npR
        xi = 1.d0 - Sdr%node(npR-j,npR)
        call CR%sample_at(xi, x)
        M(i,j,:) = x
     enddo

     ! main SOL
     do j=0,np0
        xi = Sp%node(j,np0)
        call C0%sample_at(xi, x)
        M(i,npR+j,:) = x
     enddo

     ! left divertor leg
     call divertor_leg_interface(CL, C_guide, xiL)
     call Sdl%init_spline_X1(etaL(1), xiL)
     do j=1,npL
        xi = Sdl%node(j,npL)
        call CL%sample_at(xi, x)
        M(i,npR+np0+j,:) = x
     enddo
  enddo

 1020 format (8x,'generating scrape-off layer: 1 -> ', i0)
 1021 format (8x,'d_SOL = ',f8.3)
  end subroutine make_flux_surfaces_SOL
  !=============================================================================

  !=============================================================================
  ! private flux region
  !=============================================================================
  subroutine make_flux_surfaces_PFR(M, nr, npL, npR, ir1, ir2, rpath, Sr, Sp)
  use run_control, only: Debug
  use flux_surface_2D
  use divertor

  real(real64), dimension(:,:,:), pointer, intent(inout) :: M
  integer, intent(in) :: nr, npL, npR, ir1, ir2
  type(t_xpath),   intent(in) :: rpath
  type(t_spacing), intent(in) :: Sr, Sp

  type(t_flux_surface_2D) :: F
  type(t_spacing) :: Sdr, Sdl
  real(real64) :: eta, xi, xiL, xiR, x0(2), x(2)
  integer      :: i, j

  write (6, 1030) nr-1
  write (6, 1031) d_PFR(1)
  do i=0,nr-1
     write (6, *) i
     eta = Sr%node(i,nr)
     call rpath%sample_at(1.d0 - eta, x0)
     if (Debug) write (iud, *) x0

     ! right divertor leg
     call F%generate(x0, -1, AltSurf=C_cutR, sampling=DISTANCE)
     call divertor_leg_interface(F%t_curve, C_guide, xiR)
     call Sdr%init_spline_X1(etaR(1), 1.d0-xiR)
     do j=0,npR
        xi = 1.d0 - Sdr%node(npR-j,npR)
        call F%sample_at(xi, x)
        M(i,j,:) = x
     enddo

     ! left divertor leg
     call F%generate(x0,  1, AltSurf=C_cutL, sampling=DISTANCE)
     call divertor_leg_interface(F%t_curve, C_guide, xiL)
     call Sdl%init_spline_X1(etaL(1), xiL)
     do j=1,npL
        xi = Sdl%node(j,npL)
        call F%sample_at(xi, x)
        M(i,npR + j,:) = x
     enddo
  enddo

 1030 format (8x,'generating private flux region: 0 -> ', i0)
 1031 format (8x,'d_PFR = ',f8.3)
  end subroutine make_flux_surfaces_PFR
  !=============================================================================





!===============================================================================
! FIX GRID for M3D-C1 configuration (connect block boundaries)
! This is necessary because a small deviation between field lines starting from
! the exact same location can occur. This effect is related to the order in
! which the FIO library checks the mesh elements and which depends on the
! results of previous searches.
!===============================================================================
  subroutine fix_interfaces_for_m3dc1 ()
  use emc3_grid
  implicit none

  real(real64) :: R, Z, dmax
  integer      :: ib, iz0, iz1, iz2, nt, npL1, np, npR1


  write (6, 1000)
 1000 format(8x,'fixing interfaces between blocks for M3D-C1 configurations...')
  dmax = 0.d0
  do ib=0,blocks-1
     iz0  = 3*ib
     iz1  = iz0  +  1
     iz2  = iz0  +  2
     nt   = ZON_TORO(iz0)

     npL1 = Block(ib)%npL(1)
     np   = Block(ib)%np(0)
     npR1 = Block(ib)%npR(1)

     ! connect right divertor leg
     call fix_interface(iz1, iz2, 0, 0, nt, 0, 0, npR1, 0, ZON_RADI(iz2), dmax)

     ! connect core
     call fix_interface(iz0, iz1, 0, 0, nt, 0, npR1, np, ZON_RADI(iz0), 0, dmax)

     ! connect left divertor leg
     call fix_interface(iz1, iz2, 0, 0, nt, npR1+np, npR1, npL1, 0, ZON_RADI(iz2), dmax)
  enddo

  write (6, *) 'max. deviation: ', dmax
  end subroutine fix_interfaces_for_m3dc1
  !=============================================================================



  !=============================================================================
  subroutine post_process_grid_lsn()
  use divertor

  integer :: iblock, iz, iz1, iz2


  write (6, 1000)
 1000 format(3x,'- Post processing fieldline grid')

  write (6, 1001)
 1001 format(8x,'closing grid at last divertor cells')
  do iblock=0,blocks-1
     iz1 = 3*iblock + 1
     iz2 = 3*iblock + 2
     do iz=iz1,iz2
        call close_grid_domain(iz)
     enddo
  enddo

  call fix_interfaces_for_m3dc1 ()

  end subroutine post_process_grid_lsn
  !=============================================================================

end module topo_lsn
