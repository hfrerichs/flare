module base_mesh
  use iso_fortran_env
  use grid
  use curve2D
  use separatrix
  use xpaths
  use mesh_interface
  use mod_zone
  use mfs_mesh
  implicit none
  private


  ! magnetic axis
  real(real64)  :: Pmag(2)

  ! X-point(s), separatrix(ces) and radial paths
  type(t_separatrix), dimension(:),   allocatable :: S
  type(t_xpath),      dimension(:,:), allocatable :: R
  integer,            dimension(:),   allocatable :: connectX
  type(t_curve) :: S0, S0L, S0R
  integer       :: nX

  ! guiding surface
  type(t_curve) :: C_guide

  type(t_mfs_mesh), dimension(:), allocatable :: M, Mtmp


  public :: setup_topology
  public :: setup_geometry
  public :: setup_interfaces
  public :: generate_base_mesh

  contains
!=======================================================================



!=======================================================================
  subroutine setup_topology()
  use fieldline_grid

  select case(topology)
  ! lower single null (LSN)
  case(TOPO_LSN, TOPO_LSN1)
     call initialize_zones(6)
     call initialize_interfaces(3) ! radial interfaces
     nX = 1;  allocate(connectX(nX))
     connectX(1) = 1
     layers      = 3

     ! innermost domain
     call Z(1)%setup_boundary(LOWER, POLOIDAL, PERIODIC) ! periodic poloidal boundaries
     call Z(1)%setup_boundary(UPPER, POLOIDAL, PERIODIC) ! periodic poloidal boundaries
     call Z(1)%setup_boundary(LOWER, RADIAL,   CORE)     ! core boundary
     call Z(1)%setup_mapping (UPPER, RADIAL,   Z(2), 1)  ! connect to main SOL

     ! main SOL
     call Z(2)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(2)%setup_mapping (LOWER, POLOIDAL, Z(3), 0)  ! connect to right divertor leg
     call Z(2)%setup_mapping (UPPER, POLOIDAL, Z(4), 0)  ! connect to left divertor leg

     ! right divertor leg (SOL)
     call Z(3)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(3)%setup_boundary(LOWER, POLOIDAL, DIVERTOR) ! divertor target
     call Z(3)%setup_mapping (LOWER, RADIAL,   Z(5), 2)  ! connect to right PFR

     ! left divertor leg (SOL)
     call Z(4)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(4)%setup_boundary(UPPER, POLOIDAL, DIVERTOR) ! divertor target
     call Z(4)%setup_mapping (LOWER, RADIAL,   Z(6), 3)  ! connect to right PFR

     ! right divertor leg (PFR)
     call Z(5)%setup_boundary(LOWER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(5)%setup_boundary(LOWER, POLOIDAL, DIVERTOR) ! divertor target
     call Z(5)%setup_mapping (UPPER, POLOIDAL, Z(6), 0)  ! connect to left PFR

     ! left divertor leg (PFR)
     call Z(6)%setup_boundary(LOWER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(6)%setup_boundary(UPPER, POLOIDAL, DIVERTOR) ! divertor target

     call undefined_zone_boundary_check(.true.)


  ! DDN
  case(TOPO_DDN, TOPO_DDN1)
     call initialize_zones(6)
     call initialize_interfaces(8) ! radial interfaces
     nX = 2;  allocate(connectX(nX))
     connectX(1) = -2
     connectX(2) = -2
     layers      = 6

     ! innermost domain
     call Z(1)%setup_boundary(LOWER, RADIAL,   CORE)     ! core boundary
     call Z(1)%setup_mapping (UPPER, RADIAL,   Z(3), 1)  ! connect to main SOL
     call Z(2)%setup_boundary(LOWER, RADIAL,   CORE)     ! core boundary
     call Z(2)%setup_mapping (UPPER, RADIAL,   Z(4), 2)  ! connect to main SOL
     call Z(1)%setup_mapping (UPPER, POLOIDAL, Z(2), 0)  ! connect left and right segments
     call Z(2)%setup_mapping (UPPER, POLOIDAL, Z(1), 0)  ! connect left and right segments

     !...


  ! CDN
  case(TOPO_CDN, TOPO_CDN1)
     !call initialize_interfaces(3) ! radial interfaces
     nX = 2;  allocate(connectX(nX))
     connectX(1) = 2
     connectX(2) = 1
     !ncell = 6
     !allocate (neighbor(ncell,4,2))

  case default
     write (6, 9000) trim(topology)
     stop
  end select

 9000 format('error: invalid topology ', a, '!')
  end subroutine setup_topology
!=======================================================================



!=======================================================================
! Set up geometry of computational domain:
! Magnetic axis, X-points, separatrix(ces) and radial paths from X-points
! SHOULD THIS BE MOVED TO MODULE fieldline_grid?
!=======================================================================
  subroutine setup_geometry()
  use fieldline_grid, only: guiding_surface, d_SOL, d_PFR
  use boundary,       only: S_axi, n_axi
  use equilibrium,    only: get_magnetic_axis, get_poloidal_angle, get_PsiN, Xp
  use math
  use separatrix
  use inner_boundary
  use string

  character(len=2) :: Sstr
  real(real64)     :: dx, tmp(3), theta_cut, PsiN
  integer          :: ix, ierr, iSOL, iPFR, jx, k, orientation


  ! 1. set up guiding surface for divertor legs (C_guide) --------------
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


  ! 2. CRITICAL POINTS -------------------------------------------------
  ! 2.a set up magnetic axis (Pmag) ------------------------------------
  tmp = get_magnetic_axis(0.d0); Pmag = tmp(1:2)
  write (6, 2000) Pmag
 2000 format(8x,'Magnetic axis at: ',2f10.4)


  ! 2.b check X-points -------------------------------------------------
  do ix=1,nX
     tmp(1:2) = Xp(ix)%load(ierr)
     if (ierr .ne. 0) then
        write (6, 9001) ix
        stop
     endif
     write (6, 2001) ix, tmp(1:2)
     write (6, 2002) Xp(ix)%theta/pi*180.d0
  enddo
  write (6, *)
 2001 format(8x,i0,'. X-point at: ',2f10.4)
 2002 format(11x,'-> poloidal angle [deg]: ',f10.4)
 9001 format('error: ',i0,'. X-point is not defined!')


  ! 3. inner boundaries for plasma transport domain (EMC3) -------------
  call load_inner_boundaries(Xp(1)%theta)
  write (6, *)


  ! 4.1 generate separatrix(ces) ---------------------------------------
  allocate (S(nX), R(nX,4))
  do ix=1,nX
     ! find cut-off poloidal angle for guiding X-point
     jx = connectX(ix)
     if (jx < 0  .and.  abs(jx).ne.ix) then
        orientation = DESCENT_CORE ! TODO: find orientation from PsiN(ix)-PsiN(abs(connectX(ix)))
        call R(abs(jx),orientation)%generateX(abs(jx), orientation, LIMIT_PSIN, Xp(ix)%PsiN())
        theta_cut = get_poloidal_angle(R(abs(jx),orientation)%boundary_node(UPPER))

        call S(ix)%generate(ix, C_cutl=C_guide, C_cutR=C_guide, theta_cut=theta_cut)
     else
        call S(ix)%generate(ix, C_cutl=C_guide, C_cutR=C_guide, iconnect=jx)
     endif
     write (Sstr, 2003) ix; call S(ix)%plot(trim(Sstr), parts=.true.)

     call S(ix)%M1%setup_length_sampling()
     call S(ix)%M2%setup_length_sampling()
     call S(ix)%M3%setup_length_sampling()
     call S(ix)%M4%setup_length_sampling()
  enddo
 2003 format('S',i0)
  if (connectX(1) == 1) then
     ! connect core segments of separatrix
     S0 = connect(S(1)%M1%t_curve, S(1)%M2%t_curve)
     call S0%setup_length_sampling()
  elseif (connectX(1) > 1) then
     ! connect branches for left and right core segment of separatrix
     S0R = connect(S(1)%M1%t_curve, S(2)%M2%t_curve)
     S0L = connect(S(2)%M1%t_curve, S(1)%M2%t_curve)
  else
     ! setup left and right branch of core separatrix if second X-point is used for guidance
     S0R = S(1)%M1%t_curve
     S0L = S(1)%M2%t_curve
  endif


  ! 4.2 generate radial paths from X-points ----------------------------
     ! and set up interfaces between zones
     !call Iface(1)%set_curve(S0)
  write (6, 3000)
  iSOL = 0
  iPFR = 0
  do ix=1,nX
     jx = connectX(ix)


     ! "core-interface"
     if (ix == 1  .or.  jx > 0) then
        write (6, 3010) ix
        call R(ix, DESCENT_CORE)%generateX(ix, DESCENT_CORE, LIMIT_PSIN, PsiN_in)
     else
        if (R(ix, DESCENT_CORE)%n_seg < 0) then
           write (6, *) "WARNING: rpath segment ", ix, DESCENT_CORE, " is not set up!"
           stop
        endif
     endif


     ! scrape-off layer
     if (jx == ix) then
        ! single SOL, left and right branch on same flux surface
        iSOL = iSOL + 1
        write (6, 3020) ix, d_SOL(iSOL)
        call R(ix, ASCENT_LEFT)%generateX(ix, ASCENT_LEFT, LIMIT_LENGTH, d_SOL(iSOL))
        PsiN = get_PsiN(R(ix, ASCENT_LEFT)%boundary_node(boundary=UPPER))
        call R(ix, ASCENT_RIGHT)%generateX(ix, ASCENT_RIGHT, LIMIT_PSIN, PsiN)

     elseif (jx > ix) then
        ! connected X-points, left and right branch on individual flux surfaces
        iSOL = iSOL + 1
        write (6, 3021) ix, jx, d_SOL(iSOL)
        call R(ix, ASCENT_LEFT)%generateX(ix, ASCENT_LEFT, LIMIT_LENGTH, d_SOL(iSOL))
        PsiN = get_PsiN(R(ix, ASCENT_LEFT)%boundary_node(boundary=UPPER))
        call R(jx, ASCENT_RIGHT)%generateX(jx, ASCENT_RIGHT, LIMIT_PSIN, PsiN)

        iSOL = iSOL + 1
        write (6, 3022) ix, jx, d_SOL(iSOL)
        call R(ix, ASCENT_RIGHT)%generateX(ix, ASCENT_RIGHT, LIMIT_LENGTH, d_SOL(iSOL))
        PsiN = get_PsiN(R(ix, ASCENT_RIGHT)%boundary_node(boundary=UPPER))
        call R(jx, ASCENT_LEFT)%generateX(jx, ASCENT_LEFT, LIMIT_PSIN, PsiN)

     elseif (jx == -ix) then
        ! outer SOL with left and right branch on individual flux surfaces
        iSOL = iSOL + 1
        write (6, 3023) ix, d_SOL(iSOL)
        call R(ix, ASCENT_LEFT)%generateX(ix, ASCENT_LEFT, LIMIT_LENGTH, d_SOL(iSOL))

        iSOL = iSOL + 1
        write (6, 3024) ix, d_SOL(iSOL)
        call R(ix, ASCENT_RIGHT)%generateX(ix, ASCENT_RIGHT, LIMIT_LENGTH, d_SOL(iSOL))

     elseif (jx < 0) then
        ! this SOL's boundary is another separatrix
        PsiN = Xp(abs(jx))%PsiN()
        write (6, 3025) ix, abs(jx), PsiN
        call R(ix, ASCENT_LEFT)%generateX(ix, ASCENT_LEFT, LIMIT_PSIN, PsiN)
        call R(ix, ASCENT_RIGHT)%generateX(ix, ASCENT_RIGHT, LIMIT_PSIN, PsiN)
     endif


     ! private flux region
     iPFR = iPFR + 1
     write (6, 3030) ix, d_PFR(iPFR)
     call R(ix, DESCENT_PFR)%generateX(ix, DESCENT_PFR, LIMIT_LENGTH, d_PFR(iPFR))


     ! plot paths
     do k=1,4
        call R(ix, k)%plot(filename='rpath_'//trim(str(ix))//'_'//trim(str(k))//'.plt')
     enddo
  enddo



 3000 format(3x,'- Setting up radial paths for block-structured decomposition')
 3010 format(8x,'generating core segment for X-point ', i0)
 3020 format(8x,'generating SOL segment for X-point ', i0, ' (length = ', f0.3, ' cm)')
 3021 format(8x,'generating left SOL segment for X-points ', i0, ', ', i0, ' (length = ', f0.3, ' cm)')
 3022 format(8x,'generating right SOL segment for X-points ', i0, ', ', i0, ' (length = ', f0.3, ' cm)')
 3023 format(8x,'generating left SOL segment for X-point ', i0, ' (length = ', f0.3, ' cm)')
 3024 format(8x,'generating right SOL segment for X-point ', i0, ' (length = ', f0.3, ' cm)')
 3025 format(8x,'generating SOL segment for X-point ', i0, ' up to separatrix from X-point ', &
                i0, ' at PsiN = ', f0.3)
 3030 format(8x,'generating PFR segment for X-point ', i0, ' (length = ', f0.3, ' cm)')
  end subroutine setup_geometry
!=======================================================================



!=======================================================================
  subroutine setup_interfaces()
  use mesh_interface
  use string


  integer :: ix, jx, i


  i = 0
  do ix=1,nX
     jx = connectX(ix)

     ! connect back to same X-point
     if (jx == ix) then
        if (ix .ne. 1) then
           write (6, 9000) ix
           stop
        endif

        i = i + 1
        call radial_interface(i)%set_curve(S0)
        call radial_interface(i)%setup(ix, ix)

     ! all branches connect to divertor targets
     elseif (jx == -ix) then
        i = i + 1
        call radial_interface(i)%set_curve(S(ix)%M1%t_curve)
        call radial_interface(i)%setup(STRIKE_POINT, ix)
        i = i + 1
        call radial_interface(i)%set_curve(S(ix)%M2%t_curve)
        call radial_interface(i)%setup(ix, STRIKE_POINT)

     ! connect to other X-point OR
     ! main separatrix decomposition is guided by secondary X-point
     elseif (jx > ix  .or.  jx < 0) then
        if (ix .ne. 1) then
           if (jx > ix) then
              write (6, 9001) ix
           else
              write (6, 9002) ix, abs(jx)
           endif
           stop
        endif

        ! right core interface
        i = i + 1
        call radial_interface(i)%set_curve(S0R)
        call radial_interface(i)%setup(ix, jx)

        ! left core interface
        i = i + 1
        call radial_interface(i)%set_curve(S0L)
        call radial_interface(i)%setup(jx, ix)

     ! nothing to be done here anymore
     else

     endif

     ! divertor branches
     i = i + 1
     call radial_interface(i)%set_curve(S(ix)%M3%t_curve)
     call radial_interface(i)%setup(STRIKE_POINT, ix)
     i = i + 1
     call radial_interface(i)%set_curve(S(ix)%M4%t_curve)
     call radial_interface(i)%setup(ix, STRIKE_POINT)


  enddo


  do i=1,radial_interfaces
     call radial_interface(i)%C%plot(filename='I'//trim(str(i))//'.plt')
  enddo

 9000 format('error: seconday X-point ', i0, ' connects back to itself!')
 9001 format('error: seconday X-point ', i0, ' does not connect back to primary one!')
 9002 format('error: seconday X-point ', i0, ' used as guiding point for separatrix ', i0, '!')
  end subroutine setup_interfaces
!=======================================================================



!=======================================================================
  subroutine setup_layers(iblock)
  use fieldline_grid
  integer, intent(in) :: iblock

  real(real64) :: phi
  integer      :: il


  ! initialize block
  call load_local_resolution(iblock)


  il = 0
  layers = 0
!  phi = Block(iblock)%phi_base / 180.d0 * pi
!  do il=0,layers-1
!     !call M(il)%initialize(nr(il), np(il), phi)
!  enddo

  end subroutine setup_layers
!=======================================================================



!=======================================================================
  subroutine generate_base_mesh(iblock)
  use fieldline_grid
  use mesh_spacing
  use inner_boundary
  integer, intent(in) :: iblock

  type(t_mfs_mesh)   :: Mtmp2(2), Mtmp3(3), Mtmp4(4)
  type(t_spacing)    :: Sp, SpL, SpR, Sr

  real(real64) :: phi
  integer      :: il, il0, iz, iz0


  ! initialize block
  il0 = iblock * layers
  allocate (M(0:layers-1))
  !call setup_layers(iblock)


  call load_local_resolution(iblock)

  phi = Block(iblock)%phi_base / 180.d0 * pi
  do il=0,layers-1
     call M(il)%initialize(nr(il), np(il), phi)
  enddo

  ! mesh discretization for poloidal zones
  allocate (Mtmp(nzone))


  ! generate core-interface
  il = 0
  if (connectX(1) == 1) then
     ! single zone
     call Sp%init(poloidal_spacing(il))
     call Mtmp(1)%initialize(nr(il), np(il), phi)
     call Mtmp(1)%setup_boundary_nodes(UPPER, RADIAL, S0, Sp)
  else
     ! left and right sub-zones
     call SpR%init(poloidal_spacing_R(il))
     call Mtmp(1)%initialize(nr(il), npR(il), phi)
     call Mtmp(1)%setup_boundary_nodes(UPPER, RADIAL, S0R, SpR)

     call SpL%init(poloidal_spacing_L(il))
     call Mtmp(2)%initialize(nr(il), npL(il), phi)
     call Mtmp(2)%setup_boundary_nodes(UPPER, RADIAL, S0L, SpL)
  endif


  ! generate "closed" domain
!  call Sr%init(radial_spacing(0))
!  call M(0)%setup_boundary_nodes(POLOIDAL, LOWER, R(1,DESCENT_CORE)%t_curve, Sr, nr(0), 1)
!  if (connectX(1) == 1  .or. connectX(1) < 0) then
!     call M(0)%make_orthogonal_grid(periodic=.true., rrange=(/2+n_interpolate, nr(0)-1/))
!  else
!     call M(0)%make_orthogonal_grid(periodic=.true., rrange=(/2+n_interpolate, nr(0)-1/), addX=(/abs(connectX(1)), npR(0)/))
!  endif
!  call M(0)%make_interpolated_mesh(2+n_interpolate, Sr, C_in(iblock,:), DPsiN1(iblock,1))


!  ! connect core to SOL
!  if (connectX(1) == 1) then
!     call Mtmp3(2)%initialize(nr(1), np(0), phi)
!     call M(0)%connect_to(Mtmp3(2), RADIAL, UPPER_TO_LOWER)
!     call Mtmp3(2)%store(filename='Mtmp3_2.plt')
!  else
!  endif


  !radial_scan: do
  !il     = 0
  iz0    = 1
  !ii     = 1
  !npz    = 0
  !npz(0) = 1
  !call generate_layer(0, 1)
  do il=0,layers-1
     ! set up radial spacings
     call Sr%init(radial_spacing(il))

     !iz0 = ...
     !call M(0)%setup_boundary_nodes(POLOIDAL, LOWER, R(1,DESCENT_CORE)%t_curve, Sr, nr(0), 1)

     call generate_layer(il, iz0, Sr)
     exit
  enddo



  ! write output files
  do il=0,layers-1
     iz = il0 + il
     call write_base_grid(M(il)%t_grid, iz)
  enddo


  ! cleanup
  deallocate (M, Mtmp)

  end subroutine generate_base_mesh
!=======================================================================



!=======================================================================
  subroutine generate_layer(il, iz0, Sr)
  use mesh_spacing
  integer,         intent(in) :: il, iz0
  type(t_spacing), intent(in) :: Sr

  integer :: idir, irside, iri, ipside, ipi, iz, npz(-1:1)


  ! select upper or lower radial boundary for reference nodes
  if (Mtmp(iz0)%ir0 == 0) then
     irside = LOWER
  elseif (Mtmp(iz0)%ir0 > 0) then
     irside = UPPER
  else
     write (6, 9000) il, iz0
     write (6, 9001)
     stop
  endif
  iri = Z(iz0)%rad_bound(irside)
  if (iri == UNDEFINED) then
     write (6, 9000) il, iz0
     write (6, 9002)
     stop
  endif


  ! select upper or lower poloidal boundary -> start from an X-point
  do ipside=-1,1,2
     if (radial_interface(iri)%inode(ipside) > 0) exit
  enddo
  ipi = Z(iz0)%pol_bound(ipside)
  if (ipi == UNDEFINED) then
     write (6, 9000) il, iz0
     write (6, 9003)
     stop
  endif


  ! initialize radial discretization in layer
  call Mtmp(iz0)%setup_boundary_nodes(ipside, POLOIDAL, poloidal_interface(ipi)%C, Sr)


  ! how many poloidal zones in this layer?
  call Z(iz0)%generate_mesh()
  npz    = 0
  npz(0) = 1
  idir_loop: do idir=1,-1,-2
     iz = iz0
     poloidal_scan: do
        if (Z(iz)%map_p(idir) == PERIODIC) exit
        if (Z(iz)%map_p(idir) == DIVERTOR) exit
        if (Z(iz)%map_p(idir) == iz0) exit idir_loop

        iz        = Z(iz)%map_p(idir)
        call Z(iz)%generate_mesh()
        npz(idir) = npz(idir) + 1
     enddo poloidal_scan

!     call Sp%init(poloidal_spacing(0))
!     ! single zone?
!     if (Z(iz)%map_p(UPPER) == PERIODIC) then
!
     !else
     !endif
  enddo idir_loop
  write (6, *) 'poloidal zones in layer 0: ', npz

 9000 format('error in generate_layer for il, iz0 = ', i0, ', ', i0)
 9001 format('undefined reference nodes on radial boundary!')
 9002 format('undefined radial interface!')
 9003 format('undefined poloidal interface!')
  end subroutine generate_layer
!=======================================================================

end module base_mesh
