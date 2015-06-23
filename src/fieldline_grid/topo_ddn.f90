!===============================================================================
! Disconnected Double Null configuration: block-structured decomposition with zones for
! high pressure region (HPR), inner and outer scrape-off layer (SOL) and
! upper and lower private flux regions (PFR)
!===============================================================================
module topo_ddn
  use iso_fortran_env
  use grid
  use separatrix
  use curve2D
  use equilibrium
  use fieldline_grid, unused => TOPO_DDN
  use inner_boundary
  use topo_lsn
  implicit none
  private

  integer, parameter :: DEFAULT = 0
  integer, parameter :: iud = 72


  character(len=*), parameter :: ZONE_LABEL(0:5) = (/ 'HPR  ', 'SOL1 ', 'SOL2a', 'SOL2b', 'PFR1 ', 'PFR2 ' /)


  ! coordinates of X-point and magnetic axis
  real(real64) :: Px1(2), Px2(2), Pmag(2), dtheta



  ! discretization method
  integer :: method = DEFAULT


  ! base grid in high pressure region (HPR), scrape-off layer (SOL) and private flux region (PFR)
  type(t_grid), dimension(:,:), allocatable :: G ! (0:blocks-1,0:layers-1)

  ! magnetic separatrix
  !type(t_separatrix) :: Si, So
  !type(t_separatrix) :: S2
  !type(t_curve)      :: S0

  ! guiding_surface
  !type(t_curve)      :: C_guide, C_cutL, C_cutR


  public :: &
     setup_topo_ddn, &
     make_base_grids_ddn, &
     post_process_grid_ddn

  contains
  !=====================================================================



  !=====================================================================
  subroutine setup_topo_ddn()
  use emc3_grid, only: NZONET
  implicit none

  integer :: iz, iz0, ib


  ! 0. setup number of zones for lower single null topology
  layers = 6
  NZONET = blocks * layers


  write (6, 1000)
  write (6, 1001)
  do ib=0,blocks-1
     ! 1. setup resolution for each zone
     ! set derived parameters
     Block(ib)%np(0) = Block(ib)%npR(0) + Block(ib)%npL(0)
     Block(ib)%np(1) = Block(ib)%npR(1) + Block(ib)%np(0)  + Block(ib)%npL(1)
     Block(ib)%np(2) = Block(ib)%npL(2) + Block(ib)%npL(0) + Block(ib)%npL(1)
     Block(ib)%np(3) = Block(ib)%npR(1) + Block(ib)%npR(0) + Block(ib)%npR(2)
     Block(ib)%np(4) = Block(ib)%npR(1)                    + Block(ib)%npL(1)
     Block(ib)%np(5) = Block(ib)%npR(2)                    + Block(ib)%npL(2)


     ! high pressure region (HPR)
     iz0 = ib * layers
     Zone(iz0)%nr    = Block(ib)%nr(0) + nr_EIRENE_core
     Zone(iz0)%np    = Block(ib)%np(0)

     ! inner scrape-off layer (SOL)
     Zone(iz0+1)%nr  = Block(ib)%nr(1)
     Zone(iz0+1)%np  = Block(ib)%np(1)

     ! outer scrape-off layer (SOL)
     Zone(iz0+2)%nr  = Block(ib)%nr(2) + nr_EIRENE_vac
     Zone(iz0+2)%np  = Block(ib)%np(2)
     Zone(iz0+3)%nr  = Block(ib)%nr(3) + nr_EIRENE_vac
     Zone(iz0+3)%np  = Block(ib)%np(3)

     ! private flux region (PFR)
     Zone(iz0+4)%nr  = Block(ib)%nr(4) + nr_EIRENE_vac
     Zone(iz0+4)%np  = Block(ib)%np(4)
     Zone(iz0+5)%nr  = Block(ib)%nr(5) + nr_EIRENE_vac
     Zone(iz0+5)%np  = Block(ib)%np(5)

     do iz=iz0,iz0+layers-1
        ! setup toroidal discretization
        call Zone(iz)%setup_default_toroidal_discretization(ib)
     enddo


     ! 2. setup boundaries and connectivity between zones
     do iz=iz0,iz0+layers-1
        call Zone(iz)%setup_default_boundaries()
     enddo

     ! 2.a high pressure region (HPR)
     Zone(iz0)%isfr(2) = SF_MAPPING
     Zone(iz0)%r_surf_pl_trans_range(2) = Zone(iz0)%nr

     ! 2.b inner scrape-off layer (SOL)
     Zone(iz0+1)%isfr(1) = SF_MAPPING
     Zone(iz0+1)%isfr(2) = SF_MAPPING
     Zone(iz0+1)%isfp(1) = SF_VACUUM
     Zone(iz0+1)%isfp(2) = SF_VACUUM
     Zone(iz0+1)%r_surf_pl_trans_range(1) = 0
     Zone(iz0+1)%r_surf_pl_trans_range(2) = Zone(iz0+1)%nr
     ! outer scrape-off layer
     do iz=iz0+2,iz0+3
        Zone(iz)%isfr(1) = SF_MAPPING
        Zone(iz)%isfp(1) = SF_VACUUM
        Zone(iz)%isfp(2) = SF_VACUUM
        Zone(iz)%r_surf_pl_trans_range(1) = 0
        Zone(iz)%d_N0    = d_N0(iz-iz0)
     enddo

     ! 2.c private flux region (PFR)
     do iz=iz0+4,iz0+5
        Zone(iz)%isfr(1) = SF_VACUUM
        Zone(iz)%isfr(2) = SF_MAPPING
        Zone(iz)%isfp(1) = SF_VACUUM
        Zone(iz)%isfp(2) = SF_VACUUM
        Zone(iz)%r_surf_pl_trans_range(1) = nr_EIRENE_vac
        Zone(iz)%r_surf_pl_trans_range(2) = Zone(iz)%nr
        Zone(iz)%d_N0    = d_N0(iz-iz0)
     enddo

     write (6, 1002) ib, ZONE_LABEL(0),      Zone(iz0)%nr, Zone(iz0)%np, Zone(iz0)%nt
     do iz=iz0+1,iz0+layers-1
        write (6, 1003)  ZONE_LABEL(iz-iz0), Zone(iz )%nr, Zone(iz )%np, Zone(iz )%nt
     enddo
  enddo

 1000 format(8x,'Grid resolution in')
 1001 format(8x,'block #, zone #, (radial x poloidal x toroidal):')
 1002 format(12x,i3,3x,a6,3x'(',i0,' x ',i0,' x ',i0,')')
 1003 format(12x,   6x,a6,3x'(',i0,' x ',i0,' x ',i0,')')
  end subroutine setup_topo_ddn
  !=====================================================================



  !=====================================================================
  subroutine setup_domain
  use boundary
  use run_control, only: Debug
  use math
  use flux_surface_2D, only: RIGHT_HANDED
  use inner_boundary

  real(real64) :: tmp(3), dx(2)
  integer      :: ix


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
  do ix=1,2
     write (6, 2000) Xp(ix)%load()
     write (6, 2001) Xp(ix)%theta/pi*180.d0
  enddo
  Px1 = Xp(1)%X
  Px2 = Xp(2)%X
  dtheta = Xp(2)%theta - Xp(1)%theta


  ! 2.c separatrix (S, S0) ---------------------------------------------
  call S(1)%generate(1, 1,  Xp(2)%theta, C_cutL, C_cutR)
  call S(1)%plot('Si', parts=.true.)
  call S(2)%generate(2, -1, Xp(1)%theta, C_cutL, C_cutR)
  call S(2)%plot('So', parts=.true.)

  ! connect core segments of separatrix
  !S0 = connect(Si%M1%t_curve, Si%M2%t_curve)
  !call S0%plot(filename='S0.plt')
  !!call S0%setup_angular_sampling(Pmag)
  !call S0%setup_angular_sampling(tmp(1:2))
  call S(1)%M1%setup_angular_sampling(Pmag)
  call S(1)%M2%setup_angular_sampling(Pmag)
  call S(1)%M3%setup_length_sampling()
  call S(1)%M4%setup_length_sampling()
  call S(2)%M1%setup_length_sampling()
  call S(2)%M2%setup_length_sampling()
  call S(2)%M3%setup_length_sampling()
  call S(2)%M4%setup_length_sampling()
 2000 format(8x,'found magnetic X-point at: ',2f10.4)
 2001 format(11x,'-> poloidal angle [deg]: ',f10.4)


  ! 3. inner boundaries for EMC3 grid
  call load_inner_boundaries(Xp(1)%theta)

  ! 4. setup paths for discretization in radial direction
  ! 4.1 SOL
  call rpath(1)%generate(1, ASCENT_LEFT, LIMIT_PSIN, Xp(2)%PsiN())
  call rpath(1)%plot(filename='rpath_1.plt')
  ! 4.2 left outer SOL
  call rpath(2)%generate(2, ASCENT_LEFT, LIMIT_LENGTH, d_SOL(1))
  call rpath(2)%plot(filename='rpath_2.plt')
  ! 4.3 right outer SOL
  call rpath(3)%generate(2, ASCENT_RIGHT, LIMIT_LENGTH, d_SOL(2))
  call rpath(3)%plot(filename='rpath_3.plt')
  ! 4.4 PFR1
  call rpath(4)%generate(1, DESCENT_PFR, LIMIT_LENGTH, d_PFR(1))
  call rpath(4)%plot(filename='rpath_4.plt')
  ! 4.5 PFR2
  call rpath(5)%generate(2, DESCENT_PFR, LIMIT_LENGTH, d_PFR(2))
  call rpath(5)%plot(filename='rpath_5.plt')

  end subroutine setup_domain
  !=====================================================================



  !=====================================================================
  subroutine divide_SOL2(F, eta, side, alpha, r, C)
  use flux_surface_2D
  use math
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
  subroutine make_base_grids_ddn
  use run_control, only: Debug
  use math
  use flux_surface_2D
  use mesh_spacing
  use divertor
  use topo_lsn, only: make_flux_surfaces_PFR, make_interpolated_surfaces

  integer, parameter      :: iu = 72

  type(t_flux_surface_2D) :: FS, FSL, FSR, C0
  type(t_curve)           :: CL, CR
  type(t_spacing)         :: Sl, Sr, Sp12

  real(real64), dimension(:,:,:), pointer :: M_HPR, M_SOL1, M_SOL2a, M_SOL2b, M_PFR1, M_PFR2

  !real(real64) :: xi, eta, phi, x(2), x0(2), x1(2), x2(2), d_HPR(2), dx(2)
  !integer :: i, j, iz, iz0, iz1, iz2, nr0, nr1, nr2, np0, np1, np1l, np1r, np2
  real(real64) :: phi, xi, x(2)
  integer :: i, j, iz, iz0

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
  allocate (G(0:blocks-1, 0:layers-1))
  !.....................................................................



  do iblock=0,blocks-1
     write (6, 1001) iblock

     ! set zone indices
     iz0 = iblock*layers

     ! set local variables for resolution
     call load_local_resolution(iblock)
     !nr0 = Block(iblock)%nr(0); np0 = Block(iblock)%np(0)
     !nr1 = Block(iblock)%nr(1); np1 = Block(iblock)%np(1)
     !np1l = Block(iblock)%npL(1); np1r = Block(iblock)%npR(1)
     !nr2 = Block(iblock)%nr(2); np2 = Block(iblock)%np(2)

     ! check if radial-poloidal resolution is different from last block
     generate_flux_surfaces_HPR = .true.
     generate_flux_surfaces_SOL = .true.
     generate_flux_surfaces_PFR = .true.
!     if (iblock > 0) then
!        ! copy unperturbed flux surface discretization if resolution is the same
!        if (nr0 == Block(iblock-1)%nr(0)  .and.  np0 == Block(iblock-1)%np(0)) then
!           generate_flux_surfaces_HPR = .false.
!        endif
!
!        if (np1l == Block(iblock-1)%npL(1) .and. np1r == Block(iblock-1)%npR(1)) then
!           if (nr1 == Block(iblock-1)%nr(1)) generate_flux_surfaces_SOL = .false.
!           if (nr2 == Block(iblock-1)%nr(2)) generate_flux_surfaces_PFR = .false.
!        endif
!     endif



     ! setup cell spacings
     do i=0,layers-1
        iz = iz0 + i
        call Zone(iz)%Sr%init(radial_spacing(i))
        call Zone(iz)%Sp%init(poloidal_spacing(i))
     enddo


     ! initialize base grids in present block
     phi = Block(iblock)%phi_base / 180.d0 * pi
!     call G_HPR(iblock)%new(CYLINDRICAL, MESH_2D, 3, nr0+1, np0+1, fixed_coord_value=phi)
!     call G_SOL(iblock)%new(CYLINDRICAL, MESH_2D, 3, nr1+1, np1+1, fixed_coord_value=phi)
!     call G_PFR(iblock)%new(CYLINDRICAL, MESH_2D, 3, nr2+1, np2+1, fixed_coord_value=phi)
!     M_HPR => G_HPR(iblock)%mesh
!     M_SOL => G_SOL(iblock)%mesh
!     M_PFR => G_PFR(iblock)%mesh
     do i=0,layers-1
        call G(iblock,i)%new(CYLINDRICAL, MESH_2D, 3, nr(i)+1, np(i)+1, fixed_coord_value=phi)
     enddo
     M_HPR   => G(iblock,0)%mesh
     M_SOL1  => G(iblock,1)%mesh
     M_SOL2a => G(iblock,2)%mesh
     M_SOL2b => G(iblock,3)%mesh
     M_PFR1  => G(iblock,4)%mesh
     M_PFR2  => G(iblock,5)%mesh


     ! start grid generation
     call Sp12%init_X1(1.d0 * npR(0)/np(0), dtheta/pi2)

     ! 1. unperturbed separatrix
     call make_separatrix()
     call make_separatrix2()

     ! 2. unperturbed FLUX SURFACES
     ! 2.a high pressure region (HPR)
     if (generate_flux_surfaces_HPR) then
        call make_flux_surfaces_HPR(M_HPR, nr(0), np(0), 2+n_interpolate, nr(0)-1, Zone(iz0)%Sr, Sp12)
     else
        !G_HPR(iblock)%mesh = G_HPR(iblock-1)%mesh
     endif

     ! 2.b scrape-off layer (SOL)
     if (generate_flux_surfaces_SOL) then
        call make_flux_surfaces_SOL(M_SOL1,nr(1), npL(1), np(0), npR(1), 1, nr(1)-1, rpath(1), 1, 1, Zone(iz0+1)%Sr, Sp12)
     else
        !G_SOL(iblock)%mesh = G_SOL(iblock-1)%mesh
     endif
     call make_flux_surfaces_SOL(M_SOL2a,nr(2), npL(1), npL(0), npL(2), 1, nr(2), rpath(2), 2, 1, Zone(iz0+2)%Sr, Zone(iz0+2)%Sp)
     call make_flux_surfaces_SOL(M_SOL2b,nr(3), npR(2), npR(0), npR(1), 1, nr(3), rpath(3), 1, 2, Zone(iz0+3)%Sr, Zone(iz0+3)%Sp)

     ! 2.c private flux region (PFR)
     if (generate_flux_surfaces_PFR) then
        call make_flux_surfaces_PFR(M_PFR1, nr(4), npL(1), npR(1), 1, nr(4), rpath(4), Zone(iz0+4)%Sr, Zone(iz0+4)%Sp)
     else
        !G_PFR(iblock)%mesh = G_PFR(iblock-1)%mesh
     endif
     call make_flux_surfaces_PFR(M_PFR2, nr(5), npR(2), npL(2), 1, nr(5), rpath(5), Zone(iz0+5)%Sr, Zone(iz0+5)%Sp)

     ! 3. interpolated surfaces
     call make_interpolated_surfaces(M_HPR, nr(0), np(0), 1, 2+n_interpolate, Zone(iz0)%Sr, Sp12, C_in(iblock,:))


     ! output
     do i=0,layers-1
        iz = iz0 + i
        call write_base_grid(G(iblock,i), iz)
     enddo
     write (6, 1002) iblock
  enddo
  if (Debug) close (iu)

 1000 format(//3x,'- Setup for base grids:')
 1001 format(//1x,'Start generation of base grids for block ',i0,' ',32('.'))
 1002 format(1x,'Finished generation of base grids for block ',i0,' ',32('.'),//)
  contains
  !.....................................................................
  subroutine make_separatrix()

  ! 1. discretization of main part of 1st separatrix
  ! 1.a right segment
  do j=0,npR(0)
     xi = Zone(iz0)%Sp%node(j,npR(0))

     call S(1)%M1%sample_at(xi, x)
     M_HPR (nr(0),                   j, :) = x
     M_SOL1(   0 ,          npR(1) + j, :) = x
  enddo
  ! 1.b left segment
  do j=0,npL(0)
     xi = Zone(iz0)%Sp%node(j,npL(0))

     call S(1)%M2%sample_at(xi, x)
     M_HPR (nr(0), npR(0) +          j, :) = x
     M_SOL1(   0 , npR(0) + npR(1) + j, :) = x
  enddo

  ! 2. discretization of right separatrix leg
  call divertor_leg_interface(S(1)%M3%t_curve, C_guide, xiR)
  call Sr%init_spline_X1(etaR(1), 1.d0-xiR)
  do j=0,npR(1)
     xi = 1.d0 - Sr%node(npR(1)-j,npR(1))
     call S(1)%M3%sample_at(xi, x)
     M_SOL1(    0,                   j,:) = x
     M_PFR1(nr(4),                   j,:) = x
  enddo

  ! 3. discretization of left separatrix leg
  call divertor_leg_interface(S(1)%M4%t_curve, C_guide, xiL)
  call Sl%init_spline_X1(etaL(1), xiL)
  do j=1,npL(1)
     xi = Sl%node(j,npL(1))
     call S(1)%M4%sample_at(xi, x)
     M_SOL1(    0, npR(1) + np(0)  + j,:) = x
     M_PFR1(nr(4), npR(1)          + j,:) = x
  enddo

  end subroutine make_separatrix
  !.....................................................................

  !.....................................................................
  subroutine make_separatrix2()

  type(t_flux_surface_2D) :: CR(2), CL(2)
  real(real64) :: x(2), xi, xiR, xiL
  integer :: j


  ! 1. right segments
  call divide_SOL2(S(2)%M2, 1.d0,  1, alphaR(1), S(1)%M3%length(), CR%t_curve)
  !call CR(2)%flip()
  call CR(2)%setup_sampling(Xp(1)%X, Xp(2)%X, Magnetic_Axis%X, 1.d0, 0.d0, dtheta, Dtheta_sampling)
  !call CR(1)%flip()
  call CR(1)%setup_length_sampling()
  !call CR(1)%plot(filename='CR1.plt')
  !call CR(2)%plot(filename='CR2.plt')

  ! 1.1 core segment
  do j=0,npR(0)
     xi = Zone(iz0)%Sp%node(j,npR(0))

     call CR(2)%sample_at(xi, x)
     M_SOL1 (nr(1),          npR(1) + j, :) = x
     M_SOL2b(   0 ,          npR(1) + j, :) = x
  enddo

  ! 1.2 primary divertor segment
  call divertor_leg_interface(CR(1)%t_curve, C_guide, xiR)
  call Sr%init_spline_X1(etaR(1), 1.d0-xiR)
  do j=0,npR(1)
     xi = 1.d0 - Sr%node(npR(1)-j,npR(1))
     call CR(1)%sample_at(xi, x)
     M_SOL1 (nr(1),                   j,:) = x
     M_SOL2b(   0 ,                   j,:) = x
  enddo

  ! 1.3 secondary divertor segment
  call divertor_leg_interface(S(2)%M4%t_curve, C_guide, xiR)
  call Sr%init_spline_X1(etaR(1), xiR)
  do j=0,npR(2)
     xi = Sr%node(j,npR(2))
     call S(2)%M4%sample_at(xi, x)
     M_PFR2 (nr(5), npL(2)          + j,:) = x
     M_SOL2b(   0 , npR(1) + npR(0) + j,:) = x
  enddo


  ! 2. left segments
  call divide_SOL2(S(2)%M1, 1.d0, -1, alphaL(1), S(1)%M4%length(), CL%t_curve)
  !call CL(2)%flip()
  call CL(2)%setup_length_sampling()
  !call CL(1)%flip()
  call CL(1)%setup_sampling(Xp(2)%X, Xp(1)%X, Magnetic_Axis%X, 0.d0, 1.d0, pi2-dtheta, Dtheta_sampling)
  !call CL(1)%plot(filename='CL1.plt')
  !call CL(2)%plot(filename='CL2.plt')

  ! 2.1 core segment
  do j=0,npL(0)
     xi = Zone(iz0)%Sp%node(j,npL(0))

     call CL(1)%sample_at(xi, x)
     M_SOL1 (nr(1), npR(1) + npR(0) + j, :) = x
     M_SOL2a(   0 ,          npL(2) + j, :) = x
  enddo

  ! 2.2 primary divertor segment
  call divertor_leg_interface(CL(2)%t_curve, C_guide, xiL)
  call Sl%init_spline_X1(etaL(1), xiL)
  do j=0,npL(1)
     xi = Sl%node(j,npL(1))
     call CL(2)%sample_at(xi, x)
     M_SOL1 (nr(1), npR(1) + npR(0) + npL(0) + j,:) = x
     M_SOL2a(   0 ,          npL(2) + npL(0) + j,:) = x
  enddo

  ! 2.3 secondary divertor segment
  call divertor_leg_interface(S(2)%M3%t_curve, C_guide, xiL)
  call Sl%init_spline_X1(etaL(1), 1.d0-xiL)
  do j=0,npL(2)
     xi = 1.d0 - Sl%node(npL(2)-j,npL(2))
     call S(2)%M3%sample_at(xi, x)
     M_PFR2 (nr(5),                   j,:) = x
     M_SOL2a(   0 ,                   j,:) = x
  enddo

  end subroutine make_separatrix2
  !.....................................................................

  end subroutine make_base_grids_ddn
  !=============================================================================



  !=============================================================================
  ! unperturbed FLUX SURFACES (high pressure region)
  !=============================================================================
  subroutine make_flux_surfaces_HPR(M, nr, np, ir1, ir2, Sr, Sp)
  use run_control, only: Debug
  use flux_surface_2D
  use mesh_spacing
  !use divertor

  real(real64), dimension(:,:,:), pointer, intent(inout) :: M
  integer, intent(in) :: nr, np, ir1, ir2
  type(t_spacing), intent(in) :: Sr, Sp

  type(t_flux_surface_2D) :: F
  real(real64) :: d_HPR(2)
  real(real64) :: eta, x(2), xi
  integer      :: i, j


  ! 1. get radial width at poloidal angle of X-point
  d_HPR = get_d_HPR(Px1, Pmag)

  ! 2. generate flux surfaces
  if (ir2 .ge. ir1) then
     write (6, 1010) ir2, ir1
     write (6, 1011) Px1, Px1 + d_HPR
  endif
  do i=ir2, ir1, -1
     write (6, *) i
     eta = 1.d0 - Sr%node(i-1,nr-1)

     x = Px1 + eta * d_HPR
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
 1011 format (8x,'from (',f8.3,', ',f8.3,') to (',f8.3,', ',f8.3,')')
  end subroutine make_flux_surfaces_HPR
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
  subroutine post_process_grid_ddn()
  use divertor

  integer :: iblock, iz, iz1, iz2


  write (6, 1000)
 1000 format(3x,'- Post processing fieldline grid')

  write (6, 1001)
 1001 format(8x,'closing grid at last divertor cells')
  do iblock=0,blocks-1
     iz1 = iblock*layers + 1
     iz2 = iblock*layers + layers-1
     do iz=iz1,iz2
        call close_grid_domain(iz)
     enddo
  enddo

  call fix_interfaces_for_m3dc1 ()

  end subroutine post_process_grid_ddn
  !=============================================================================

end module topo_ddn
