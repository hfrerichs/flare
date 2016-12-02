!===============================================================================
! Disconnected Double Null configuration: block-structured decomposition with zones for
! high pressure region (HPR), inner and outer scrape-off layer (SOL) and
! upper and lower private flux regions (PFR)
!===============================================================================
module modtopo_ddn
  use iso_fortran_env
  use grid
  use separatrix
  use curve2D
  use equilibrium
  use fieldline_grid, unused => TOPO_DDN, unused2 => TOPO_LSN
  use inner_boundary
  use modtopo_lsn
  implicit none
  private

  integer, parameter :: layers_ddn = 6
  integer, parameter :: DEFAULT = 0
  integer, parameter :: iud = 72


  character(len=*), parameter :: ZONE_LABEL(0:layers_ddn-1) = (/ 'HPR   ', 'SOL(0)', 'SOL(1)', 'SOL(2)', 'PFR(1)', 'PFR(2)' /)


  ! coordinates of X-point and magnetic axis
!  real(real64) :: Px1(2), Px2(2), dtheta
  real(real64) :: dtheta



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

  integer :: iz, iz0, ib, ilayer


  ! 0. setup number of zones for disconnected double null topology
  layers = layers_ddn
  NZONET = blocks * layers
  label(0:layers-1) = ZONE_LABEL


  write (6, 1000)
  write (6, 1001)
  do ib=0,blocks-1
     iz0 = ib * layers

     ! 1. set up derived parameters
     Block(ib)%np(0) = Block(ib)%npR(0) + Block(ib)%npL(0)
     Block(ib)%np(1) = Block(ib)%npR(1) + Block(ib)%np(0)  + Block(ib)%npL(1)
     Block(ib)%np(2) = Block(ib)%npR(1) + Block(ib)%npR(0) + Block(ib)%npR(2)
     Block(ib)%np(3) = Block(ib)%npL(2) + Block(ib)%npL(0) + Block(ib)%npL(1)
     Block(ib)%np(4) = Block(ib)%npR(1)                    + Block(ib)%npL(1)
     Block(ib)%np(5) = Block(ib)%npR(2)                    + Block(ib)%npL(2)


     ! 2. set up zones
     call Zone(iz0+0)%setup(ib, 0, TYPE_HPR,    SF_PERIODIC)
     call Zone(iz0+1)%setup(ib, 1, TYPE_SOLMAP, SF_VACUUM)
     do ilayer=2,3; call Zone(iz0+ilayer)%setup(ib, ilayer, TYPE_SOL, SF_VACUUM); enddo
     do ilayer=4,5; call Zone(iz0+ilayer)%setup(ib, ilayer, TYPE_PFR, SF_VACUUM); enddo


     ! 3. show zone information
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
  use divertor

  real(real64) :: dx(2)


  dtheta = Xp(2)%theta - Xp(1)%theta
  if (dtheta < 0.d0) dtheta = dtheta + pi2


  select case(discretization_method)
  case(POLOIDAL_ANGLE)
     call S(1)%M1%setup_angular_sampling(Pmag)
     call S(1)%M2%setup_angular_sampling(Pmag)
  case(ORTHOGONAL)
     call S(1)%M1%setup_length_sampling()
     call S(1)%M2%setup_length_sampling()
  end select


  ! 4. setup paths for discretization in radial direction
  ! 4.0 HPR
  select case(discretization_method)
  case(POLOIDAL_ANGLE)
     dx    = get_d_HPR(Xp(1)%X, Pmag)
     call rpath(0)%setup_linear(Xp(1)%X, dx)
  case(ORTHOGONAL)
     call rpath(0)%generateX(1, DESCENT_CORE, LIMIT_PSIN, PsiN_in)
  end select
  call rpath(0)%plot(filename='rpath_0.plt')
  ! 4.1 SOL
  call rpath(1)%generateX(1, ASCENT_LEFT, LIMIT_PSIN, Xp(2)%PsiN())
  call rpath(1)%plot(filename='rpath_1.plt')
  ! 4.2 right outer SOL
  call rpath(2)%generateX(2, ASCENT_RIGHT, LIMIT_LENGTH, d_SOL(1))
  call rpath(2)%plot(filename='rpath_2.plt')
  ! 4.3 left outer SOL
  call rpath(3)%generateX(2, ASCENT_LEFT, LIMIT_LENGTH, d_SOL(2))
  call rpath(3)%plot(filename='rpath_3.plt')
  ! 4.4 PFR1
  call rpath(4)%generateX(1, DESCENT_PFR, LIMIT_LENGTH, d_PFR(1))
  call rpath(4)%plot(filename='rpath_4.plt')
  ! 4.5 PFR2
  call rpath(5)%generateX(2, DESCENT_PFR, LIMIT_LENGTH, d_PFR(2))
  call rpath(5)%plot(filename='rpath_5.plt')


  call check_domain()

  end subroutine setup_domain
  !=====================================================================



  !=====================================================================
  subroutine make_base_grids_ddn
  use run_control, only: Debug
  use math
  use flux_surface_2D
  use mesh_spacing
  use divertor

  integer, parameter      :: nx = 2

  real(real64), dimension(:,:,:), pointer :: M_HPR, M_SOL1, M_SOL2a, M_SOL2b, M_PFR1, M_PFR2

  type(t_spacing) :: Sp_HPR, Sp1, Sp2
  real(real64)    :: phi
  integer         :: i, iz, iz0, iblock, connectX(nx), orientationX(nx)


  write (6, 1000)
  if (Debug) then
     open  (iud, file='base_grid_debug.txt')
  endif
  !.....................................................................
  ! 0. initialize geometry
  connectX(1) = -2 ! primary X-point connects to itself, but poloidal angle
  connectX(2) = -2 ! from secondary X-point is used to split segments
  orientationX(1) = DESCENT_CORE
  orientationX(2) = 0
  call setup_geometry(nx, connectX, orientationX)
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


     ! setup cell spacings
     do i=0,layers-1
        iz = iz0 + i
        call Zone(iz)%Sr%init(radial_spacing(i))
     enddo
     ! set up poloidal spacing in core layer
     call Sp1%init(poloidal_spacing(0))
     call Sp2%init(poloidal_spacing(1))
     call Sp_HPR%init_recursive(Sp1, Sp2, 1.d0 * npR(0)/np(0), dtheta/pi2)


     ! initialize base grids in present block
     phi = Block(iblock)%phi_base / 180.d0 * pi
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
     ! 1. unperturbed separatrix
     call make_separatrix()
     call make_separatrix2()

     ! 2. unperturbed FLUX SURFACES
     ! 2.a high pressure region (HPR)
     select case(discretization_method)
     case(POLOIDAL_ANGLE)
     call make_flux_surfaces_HPR(M_HPR, nr(0), np(0), 2+n_interpolate, nr(0)-1, rpath(0), Zone(iz0)%Sr, Sp_HPR)

     case(ORTHOGONAL)
     call make_ortho_grid(M_HPR, nr(0), np(0), nr(0), 2+n_interpolate, nr(0)-1, &
                      0, 1, np(0)-1, NO_X, 0, rpath(0), Zone(iz0)%Sr, Sp_HPR, periodic=.true.)

     end select



     ! 2.b scrape-off layer (SOL)
     select case(discretization_method)
     case(POLOIDAL_ANGLE)
     call make_flux_surfaces_SOL(M_SOL1,nr(1), npL(1), np(0), npR(1), 1, nr(1)-1, rpath(1), 1, 1, Zone(iz0+1)%Sr, Sp_HPR)
     call make_flux_surfaces_SOL(M_SOL2a,nr(2), npR(2), npR(0), npR(1), 1, nr(2), rpath(2), -1, 2, Zone(iz0+2)%Sr, Sp1)
     call make_flux_surfaces_SOL(M_SOL2b,nr(3), npL(1), npL(0), npL(2), 1, nr(3), rpath(3), 2, -1, Zone(iz0+3)%Sr, Sp2)

     case(ORTHOGONAL)

     call make_flux_surfaces_SOLo(M_SOL1,nr(1), npL(1), np(0), npR(1), 1, nr(1), npR(1)+np(0), npR(1), 1, rpath(1), Zone(iz0+1)%Sr)
     call update_separatrix2()

     call make_flux_surfaces_SOLo(M_SOL2a,nr(2), npR(2), npR(0), npR(1),1, nr(2), npR(1)+npR(0), NO_X,1, rpath(2), Zone(iz0+2)%Sr)
     call make_flux_surfaces_SOLo(M_SOL2b,nr(3), npL(1), npL(0), npL(2),1, nr(3), npL(2), NO_X, 1, rpath(3), Zone(iz0+3)%Sr)
     end select


     ! 2.c private flux region (PFR)
     call make_flux_surfaces_PFR(M_PFR1, nr(4), npL(1), npR(1), 1, nr(4), rpath(4), Zone(iz0+4)%Sr, Zone(iz0+4)%Sp)
     call make_flux_surfaces_PFR(M_PFR2, nr(5), npR(2), npL(2), 1, nr(5), rpath(5), Zone(iz0+5)%Sr, Zone(iz0+5)%Sp)

     ! 3. interpolated surfaces
     select case(discretization_method)
     case(POLOIDAL_ANGLE)
        call make_interpolated_surfaces(M_HPR, nr(0), np(0), 1, 2+n_interpolate, Zone(iz0)%Sr, Sp_HPR, C_in(iblock,:))
     case(ORTHOGONAL)
        call make_interpolated_surfaces_ortho(M_HPR, nr(0), np(0), 2+n_interpolate, Zone(iz0)%Sr, Sp_HPR, &
                                           C_in(iblock,:), DPsiN1(iblock,1))
     end select


     ! output
     do i=0,layers-1
        iz = iz0 + i
        call write_base_grid(G(iblock,i), iz)
     enddo
     write (6, 1002) iblock
  enddo
  if (Debug) close (iud)

 1000 format(//3x,'- Setup for base grids:')
 1001 format(//1x,'Start generation of base grids for block ',i0,' ',32('.'))
 1002 format(1x,'Finished generation of base grids for block ',i0,' ',32('.'),//)
  contains
  !.....................................................................
  subroutine make_separatrix()

  real(real64) :: DL(0:npL(1), 2), DR(0:npR(1), 2)
  real(real64) :: xi, x(2)
  integer      :: j


  ! 1. discretization of main part of 1st separatrix
  ! 1.a right segment
  do j=0,npR(0)
     xi = Sp1%node(j,npR(0))

     call S(1)%M1%sample_at(xi, x)
     M_HPR (nr(0),                   j, :) = x
     M_SOL1(   0 ,          npR(1) + j, :) = x
  enddo
  ! 1.b left segment
  do j=0,npL(0)
     xi = Sp2%node(j,npL(0))

     call S(1)%M2%sample_at(xi, x)
     M_HPR (nr(0), npR(0) +          j, :) = x
     M_SOL1(   0 , npR(0) + npR(1) + j, :) = x
  enddo

  ! 2. discretization of right separatrix leg
  call divertor_leg_discretization(S(1)%M3%t_curve, 1.d0-etaR(1), npR(1), DR)
  M_SOL1(    0, 0:npR(1), :) = DR
  M_PFR1(nr(4), 0:npR(1), :) = DR

  ! 3. discretization of left separatrix leg
  call divertor_leg_discretization(S(1)%M4%t_curve, etaL(1), npL(1), DL)
  M_SOL1(    0, npR(1)+np(0):npR(1)+np(0)+npL(1), :) = DL
  M_PFR1(nr(4), npR(1)      :npR(1)      +npL(1), :) = DL

  end subroutine make_separatrix
  !.....................................................................

  !.....................................................................
  subroutine make_separatrix2()

  type(t_flux_surface_2D) :: CR(2), CL(2)
  real(real64) :: DL1(0:npL(1), 2), DR1(0:npR(1), 2), DL2(0:npL(2), 2), DR2(0:npR(2), 2)
  real(real64) :: x(2), xi, xiR, xiL
  integer :: j


  ! 1. right segments
  call divide_SOL2(S(2)%M2, 1.d0,  1, alphaR(1), S(1)%M3%length(), CR%t_curve)
  call CR(2)%setup_sampling(Xp(1)%X, Xp(2)%X, Pmag, Dtheta_sampling, 0.d0)
  call CR(1)%setup_length_sampling()

  ! 1.1 core segment
  do j=0,npR(0)
     xi = Sp1%node(j,npR(0))

     call CR(2)%sample_at(xi, x)
     M_SOL1 (nr(1),          npR(1) + j, :) = x
     M_SOL2a(   0 ,          npR(1) + j, :) = x
  enddo

  ! 1.2 primary divertor segment
  call divertor_leg_discretization(CR(1)%t_curve, 1.d0-etaR(1), npR(1), DR1)
  M_SOL1 (nr(1), 0:npR(1), :) = DR1
  M_SOL2a(   0 , 0:npR(1), :) = DR1

  ! 1.3 secondary divertor segment
  call divertor_leg_discretization(S(2)%M4%t_curve, etaR(1), npR(2), DR2)
  M_PFR2 (nr(5), npL(2)       :npL(2)       +npR(2), :) = DR2
  M_SOL2a(   0 , npR(1)+npR(0):npR(1)+npR(0)+npR(2), :) = DR2


  ! 2. left segments
  call divide_SOL2(S(2)%M1, 1.d0, -1, alphaL(1), S(1)%M4%length(), CL%t_curve)
  call CL(2)%setup_length_sampling()
  call CL(1)%setup_sampling(Xp(2)%X, Xp(1)%X, Pmag, 0.d0, Dtheta_sampling)

  ! 2.1 core segment
  do j=0,npL(0)
     xi = Sp2%node(j,npL(0))

     call CL(1)%sample_at(xi, x)
     M_SOL1 (nr(1), npR(1) + npR(0) + j, :) = x
     M_SOL2b(   0 ,          npL(2) + j, :) = x
  enddo

  ! 2.2 primary divertor segment
  call divertor_leg_discretization(CL(2)%t_curve, etaL(1), npL(1), DL1)
  M_SOL1 (nr(1), npR(1)+npR(0)+npL(0):npR(1)+npR(0)+npL(0)+npL(1), :) = DL1
  M_SOL2b(   0 ,        npL(2)+npL(0):       npL(2)+npL(0)+npL(1), :) = DL1

  ! 2.3 secondary divertor segment
  call divertor_leg_discretization(S(2)%M3%t_curve, 1.d0-etaL(1), npL(2), DL2)
  M_PFR2 (nr(5), 0:npL(2), :) = DL2
  M_SOL2b(   0 , 0:npL(2), :) = DL2

  end subroutine make_separatrix2
  !.....................................................................

  !.....................................................................
  subroutine update_separatrix2
  integer :: i

  ! X-point
  M_SOL1(nr(1), npR(1)+npR(0), :) = S(2)%M1%x(0,:)

  ! right segment
  M_SOL2a(   0 , 0:npR(1)+npR(0), :) = M_SOL1 (nr(1), 0:npR(1)+npR(0), :)
!  do i=0,npR(1)+npR(0)
!     write (80, *) M_SOL2a(   0 , i, :)
!  enddo

  ! left segment
  M_SOL2b(   0 ,        npL(2):       npL(2)+npL(0)+npL(1), :) = &
  M_SOL1 (nr(1), npR(1)+npR(0):npR(1)+npR(0)+npL(0)+npL(1), :)
!  do i=npL(2),npL(2)+npL(0)+npL(1)
!     write (81, *) M_SOL2b(   0 , i, :)
!  enddo
  end subroutine update_separatrix2
  !.....................................................................

  end subroutine make_base_grids_ddn
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
     iz2  = iz0  +  4
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

     ! TODO: connect secondary separatrix
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

end module modtopo_ddn
