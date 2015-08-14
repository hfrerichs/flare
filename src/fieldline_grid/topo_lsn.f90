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

  integer, parameter :: layers_lsn = 3
  integer, parameter :: DEFAULT = 0




  ! coordinates of X-point and magnetic axis
!  real(real64) :: Px(2)
!  real(real64), public :: Pmag(2)

  ! discretization method
  integer :: method = DEFAULT


  ! base grid in high pressure region (HPR), scrape-off layer (SOL) and private flux region (PFR)
  type(t_grid), dimension(:), allocatable :: G_HPR, G_SOL, G_PFR ! (0:blocks-1)

  ! magnetic separatrix
!  type(t_separatrix), public :: S(max_layers)
  type(t_curve)      :: S0

  ! guiding_surface
!  type(t_curve), public      :: C_guide, C_cutL, C_cutR


  public :: &
     setup_topo_lsn, &
     make_base_grids_lsn, &
     post_process_grid_lsn

  contains
  !=====================================================================



  !=====================================================================
  subroutine setup_topo_lsn()
  use emc3_grid, only: NZONET
  implicit none

  integer :: iz, iz0, iz1, iz2, ib


  ! 0. setup number of zones for lower single null topology
  layers = layers_lsn
  NZONET = blocks * layers


  write (6, 1000)
  write (6, 1001)
  do ib=0,blocks-1
     iz0 = ib * layers

     ! 1. set up derived parameters
     Block(ib)%np(1) = Block(ib)%npR(1) + Block(ib)%np(0) + Block(ib)%npL(1)
     Block(ib)%np(2) = Block(ib)%npR(1)                   + Block(ib)%npL(1)


     ! 2. set up zones
     call Zone(iz0+0)%setup(ib, 0, TYPE_HPR, SF_PERIODIC)
     call Zone(iz0+1)%setup(ib, 1, TYPE_SOL, SF_VACUUM)
     call Zone(iz0+2)%setup(ib, 2, TYPE_PFR, SF_VACUUM)


     ! 3. show zone information
     write (6, 1002) ib, Zone(iz0)%nr, Zone(iz0)%np, Zone(iz0)%nt, &
                         Zone(iz0+1)%nr, Zone(iz0+1)%np, Zone(iz0+1)%nt, &
                         Zone(iz0+2)%nr, Zone(iz0+2)%np, Zone(iz0+2)%nt
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
  use divertor

  real(real64) :: dx(2), Px(2)


  ! connect core segments of separatrix
  S0 = connect(S(1)%M1%t_curve, S(1)%M2%t_curve)
  call S0%plot(filename='S0.plt')
  call S0%setup_angular_sampling(Pmag)


  ! 4. setup paths for discretization in radial direction
  ! 4.0 HPR
  Px    = Xp(1)%X
  dx    = get_d_HPR(Px, Pmag)
  call rpath(0)%setup_linear(Px, dx)
  call rpath(0)%plot(filename='rpath_0.plt')
  ! 4.1 SOL
  call rpath(1)%generateX(1, ASCENT_LEFT, LIMIT_LENGTH, d_SOL(1))
  call rpath(1)%plot(filename='rpath_1.plt')
  ! 4.2 PFR
  call rpath(2)%generateX(1, DESCENT_PFR, LIMIT_LENGTH, d_PFR(1))
  call rpath(2)%plot(filename='rpath_2.plt')

  end subroutine setup_domain
  !=====================================================================



  !=====================================================================
  subroutine make_base_grids_lsn
  use run_control, only: Debug
  use math
  use inner_boundary
  use flux_surface_2D
  use mesh_spacing
  use divertor

  integer, parameter      :: nx = 1

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
  integer :: iblock, connectX(nx)


  write (6, 1000)
  if (Debug) then
     open  (iud, file='base_grid_debug.txt')
  endif
  !.....................................................................
  ! 0. initialize geometry
  connectX(1) = 1 ! X-point connects to itself
  call setup_geometry(nx, connectX)
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
        call make_flux_surfaces_HPR(M_HPR, nr0, np0, 2+n_interpolate, nr0-1, rpath(0), Zone(iz0)%Sr, Zone(iz0)%Sp)
     else
        G_HPR(iblock)%mesh = G_HPR(iblock-1)%mesh
     endif

     ! 2.b scrape-off layer (SOL)
     if (generate_flux_surfaces_SOL) then
        call make_flux_surfaces_SOL(M_SOL, nr1, np1l, np0, np1r, 1, nr1, rpath(1), 1, 1, Zone(iz1)%Sr, Zone(iz0)%Sp)
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
  if (Debug) close (iud)

 1000 format(//3x,'- Setup for base grids:')
 1001 format(//1x,'Start generation of base grids for block ',i0,' ',32('.'))
 1002 format(1x,'Finished generation of base grids for block ',i0,' ',32('.'),//)
  contains
  !.....................................................................
  subroutine make_separatrix() ! make_interfaces

  real(real64) :: DL(0:np1l, 2), DR(0:np1r, 2)

  ! 1. discretization of main part of separatrix
  do j=0,np0
     xi = Zone(iz0)%Sp%node(j,np0)

     call S0%sample_at(xi, x)
     M_HPR(nr0,      j, :) = x
     M_SOL(  0, np1r+j, :) = x
  enddo

  ! 2. discretization of right separatrix leg
  call divertor_leg_discretization(S(1)%M3%t_curve, 1.d0-etaR(1), np1r, DR)
  M_SOL(  0, 0:np1r, :) = DR
  M_PFR(nr2, 0:np1r, :) = DR


  ! 3. discretization of left separatrix leg
  call divertor_leg_discretization(S(1)%M4%t_curve, etaL(1), np1l, DL)
  M_SOL(  0, np1r+np0:np1r+np0+np1l, :) = DL
  M_PFR(nr2, np1r    :np1r    +np1l, :) = DL

  end subroutine make_separatrix
  !.....................................................................
  end subroutine make_base_grids_lsn
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
