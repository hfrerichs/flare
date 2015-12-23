module base_mesh
  use iso_fortran_env
  use grid
  use curve2D
  use separatrix
  use xpaths
  use mesh_interface
  use elements
  use mfs_mesh
  use fieldline_grid, only: t_toroidal_discretization
  implicit none
  private


  integer, parameter :: &
     LEFT   =  1, &
     CENTER =  0, &
     RIGHT  = -1

  type t_layer
     ! number of (poloidal) elements in layer
     integer   :: nz

     ! element indices
     integer, dimension(:), allocatable :: iz

     ! base element index (i0), poloidal layer (ipl) and side (ipl_side)
     integer   :: i0

     ! radial and poloidal resolution in layer
     integer   :: nr = UNDEFINED, np = UNDEFINED

     ! toroidal discretization
     type(t_toroidal_discretization) :: T


     contains
     procedure :: initialize
     procedure :: setup_resolution
     procedure :: map_poloidal_resolution
  end type t_layer
  type(t_layer), dimension(:), allocatable :: L


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
  subroutine initialize(this, nz, iz, i0)
  class(t_layer)      :: this
  integer, intent(in) :: nz, iz(nz), i0


  if (nz <= 0) then
     write (6, *) 'error in t_layer%initialize: number of elements must be positive!'
     stop
  endif


  this%nz = nz
  allocate (this%iz(nz))
  this%iz = iz

  ! set base element index
  this%i0 = i0

  end subroutine initialize
!=======================================================================



!=======================================================================
  subroutine setup_resolution(this, il)
  use fieldline_grid, only: nr, np, npR, npL
  class(t_layer)      :: this
  integer, intent(in) :: il

  integer :: i, i1, iz, ipl, ipl0


  ! 1. set radial resolution in this layer..............................
  this%nr = nr(il)
  !.....................................................................


  ! 2. set resolution in elements.......................................
  ! case A: innermost domain
  if (il == 0) then
     if (this%nz == 1) then
        iz = this%iz(1);  Z(iz)%nr = nr(il);  Z(iz)%np = np(il)
        Z(iz)%ipl      = 0
        Z(iz)%ipl_side = CENTER
     elseif (this%nz == 2) then
        iz = this%iz(1);  Z(iz)%nr = nr(il);  Z(iz)%np = npL(il)
        Z(iz)%ipl      = 0
        Z(iz)%ipl_side = LEFT

        iz = this%iz(2);  Z(iz)%nr = nr(il);  Z(iz)%np = npR(il)
        Z(iz)%ipl      = 0
        Z(iz)%ipl_side = RIGHT
     else
        write (6, 9000);  write(6, 9001);  stop
     endif

  ! case B: outer layers -> poloidal resolution is already defined in at least one element
  else
     ! 1. set radial resolution throughout layer
     do i=1,this%nz
        iz = this%iz(i);  Z(iz)%nr = nr(il)
     enddo

     ! 2. set poloidal resolution on lower/right side of layer
     ! 2.1 find index of first element with defined poloidal resolution
     do i1=1,this%nz
        iz = this%iz(i1)
        if (Z(iz)%np /= UNDEFINED) exit
     enddo
     if (i1 > this%nz) then
        write (6, 9000);  write(6, 9002);  stop
     endif
     ipl0 = Z(iz)%ipl; if (Z(iz)%ipl_side == LEFT) ipl0 = ipl0 + 1
     ! 2.2 go backwards and set up poloidal layer and corresponding resolution
     do i=i1-1,1,-1
        iz             = this%iz(i)
        ipl            = ipl0 + i1-i
        Z(iz)%np       = npR(ipl)
        Z(iz)%ipl      = ipl
        Z(iz)%ipl_side = RIGHT
     enddo

     ! 3. set poloidal resolution on upper/left side of layer
     ! 3.1 find index of first element with defined poloidal resolution
     do i1=this%nz,1,-1
        iz = this%iz(i1)
        if (Z(iz)%np /= UNDEFINED) exit
     enddo
     if (i1 < 1) then
        write (6, 9000);  write(6, 9002);  stop
     endif
     ipl0 = Z(iz)%ipl; if (Z(iz)%ipl_side == RIGHT) ipl0 = ipl0 + 1
     ! 3.2 go backwards and set up poloidal layer and corresponding resolution
     do i=i1+1,this%nz
        iz             = this%iz(i)
        ipl            = ipl0 + i-i1
        Z(iz)%np       = npL(ipl)
        Z(iz)%ipl      = ipl
        Z(iz)%ipl_side = LEFT
     enddo
  endif
  !.....................................................................


  ! 3. set poloidal resolution in layer.................................
  this%np = 0
  do i=1,this%nz
     iz      = this%iz(i)
     this%np = this%np + Z(iz)%np
  enddo
  !.....................................................................


 9000 format('error in t_layer%setup_resolution:')
 9001 format('innermost domain with ', i0, ' > 2 elements not supported!')
 9002 format('poloidal resolution undefined in all elements!')
  end subroutine setup_resolution
!=======================================================================



!=======================================================================
  subroutine map_poloidal_resolution(this)
  class(t_layer)      :: this

  integer :: i, iside, iz, iz_map


  do i=1,this%nz
     iz = this%iz(i)

     do iside=-1,1,2
        iz_map = Z(iz)%map_r(iside)
        ! map to another element?
        if (iz_map < 0) cycle

        ! map poloidal resolution
        if (Z(iz_map)%np == UNDEFINED) then
           Z(iz_map)%np       = Z(iz)%np
           Z(iz_map)%ipl      = Z(iz)%ipl
           Z(iz_map)%ipl_side = Z(iz)%ipl_side

        ! poloidal resolution in mapped element already defined?
        elseif (Z(iz_map)%np /= Z(iz)%np) then
           write (6, 9000) iz, iz_map, Z(iz)%np, Z(iz_map)%np;  stop
        endif
     enddo
  enddo

 9000 format('error in t_layer%map_poloidal_resolution:'//, 'element ', i0, ' maps to ', i0, &
             ', but poloidal resolution is ', i0, ' vs. ', i0, '!')
  end subroutine map_poloidal_resolution
!=======================================================================



!=======================================================================
  subroutine setup_topology()
  use fieldline_grid

  ! 1. initialize topology
  layers = -1
  select case(topology)
  ! lower single null (LSN)
  case(TOPO_LSN, TOPO_LSN1)
     call initialize_elements(6)
     call initialize_interfaces(3) ! radial interfaces
     nX = 1;  allocate(connectX(nX))
     connectX(1) = 1

  ! disconnected double null (DDN)
  case(TOPO_DDN, TOPO_DDN1)
     call initialize_elements(16)
     call initialize_interfaces(10) ! radial interfaces
     nX = 2;  allocate(connectX(nX))
     connectX(1) = -2
     connectX(2) = -2

  ! connected double null (CDN)
  case(TOPO_CDN, TOPO_CDN1)
     !call initialize_zones(6)
     !call initialize_interfaces(3) ! radial interfaces
     nX = 2;  allocate(connectX(nX))
     connectX(1) = 2
     connectX(2) = 1

  case default
     write (6, 9000) trim(topology)
     stop
  end select
  call initialize_poloidal_interfaces(4*nX)



  ! 2. setup element topology
  select case(topology)
  ! lower single null (LSN)
  case(TOPO_LSN, TOPO_LSN1)
     ! innermost domain
     call Z(1)%setup_boundary(LOWER, POLOIDAL, PERIODIC, 1) ! periodic poloidal boundaries at interface R1
     call Z(1)%setup_boundary(UPPER, POLOIDAL, PERIODIC, 1) ! periodic poloidal boundaries
     call Z(1)%setup_boundary(LOWER, RADIAL,   CORE)     ! core boundary
     call Z(1)%setup_mapping (UPPER, RADIAL,   Z(2), 1)  ! connect to main SOL at interface I1

     ! main SOL
     call Z(2)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(2)%setup_mapping (LOWER, POLOIDAL, Z(3), 2)     ! connect to right divertor leg at interface R2
     call Z(2)%setup_mapping (UPPER, POLOIDAL, Z(4))  ! connect to left divertor leg

     ! right divertor leg (SOL)
     call Z(3)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(3)%setup_boundary(LOWER, POLOIDAL, DIVERTOR) ! divertor target
     call Z(3)%setup_mapping (LOWER, RADIAL,   Z(5), 2)  ! connect to right PFR at interface I2

     ! left divertor leg (SOL)
     call Z(4)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(4)%setup_boundary(UPPER, POLOIDAL, DIVERTOR) ! divertor target
     call Z(4)%setup_mapping (LOWER, RADIAL,   Z(6), 3)  ! connect to right PFR at interface I3

     ! right divertor leg (PFR)
     call Z(5)%setup_boundary(LOWER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(5)%setup_boundary(LOWER, POLOIDAL, DIVERTOR) ! divertor target
     call Z(5)%setup_mapping (UPPER, POLOIDAL, Z(6), 4)  ! connect to left PFR at interface R4

     ! left divertor leg (PFR)
     call Z(6)%setup_boundary(LOWER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(6)%setup_boundary(UPPER, POLOIDAL, DIVERTOR) ! divertor target


  ! DDN
  case(TOPO_DDN, TOPO_DDN1)
     ! innermost domain
     call Z(1)%setup_boundary (LOWER, RADIAL,   CORE)     ! core boundary
     call Z(1)%setup_mapping  (UPPER, RADIAL,   Z(3), 1)  ! connect to main SOL at interface 1
     call Z(2)%setup_boundary (LOWER, RADIAL,   CORE)     ! core boundary
     call Z(2)%setup_mapping  (UPPER, RADIAL,   Z(4), 2)  ! connect to main SOL at interface 2
     call Z(1)%setup_mapping  (UPPER, POLOIDAL, Z(2))     ! connect left and right segments
     call Z(2)%setup_mapping  (UPPER, POLOIDAL, Z(1), 1)  ! connect left and right segments at interface R1

     ! primary SOL
     call Z(3)%setup_mapping  (UPPER, RADIAL,   Z(5), 8)  ! connect to right secondary SOL at PARTIAL interface I6 (this should only be used to find the poloidal side for the generating radial path!)
     call Z(4)%setup_mapping  (UPPER, RADIAL,   Z(6), 7)  ! connect to left secondary SOL at interface I5 (SAME NOTE AS ABOVE)
     call Z(3)%setup_mapping  (UPPER, POLOIDAL, Z(4))     ! connect left and right segments
     call Z(3)%setup_mapping  (LOWER, POLOIDAL, Z(7), 2)  ! connect to right divertor leg at interface R2
     call Z(7)%setup_boundary (LOWER, POLOIDAL, DIVERTOR) ! right divertor target
     call Z(4)%setup_mapping  (UPPER, POLOIDAL, Z(8))     ! connect to left divertor leg
     call Z(8)%setup_boundary (UPPER, POLOIDAL, DIVERTOR) ! left divertor target
     call Z(7)%setup_mapping  (UPPER, RADIAL,   Z(9), 5)  ! VIRTUAL interface I5
     call Z(8)%setup_mapping  (UPPER, RADIAL,   Z(12),6)  ! VIRTUAL interface I6

     ! secondary SOL
     ! right branch
     call Z(5)%setup_boundary (UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(5)%setup_mapping  (LOWER, POLOIDAL, Z(9))     ! connect to right divertor leg (right branch)
     call Z(9)%setup_boundary (LOWER, POLOIDAL, DIVERTOR) ! right divertor target
     call Z(5)%setup_mapping  (UPPER, POLOIDAL, Z(10), 6)    ! connect to left divertor leg  (right branch) at interface R6
     call Z(10)%setup_boundary(UPPER, POLOIDAL, DIVERTOR) ! left divertor target
     call Z(9)%setup_boundary (UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(10)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain
     ! left branch
     call Z(6)%setup_boundary (UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(6)%setup_mapping  (LOWER, POLOIDAL, Z(11), 7)    ! connect to right divertor leg (left branch) at interface R7
     call Z(11)%setup_boundary(LOWER, POLOIDAL, DIVERTOR) ! right divertor target
     call Z(6)%setup_mapping  (UPPER, POLOIDAL, Z(12))    ! connect to left divertor leg  (left branch)
     call Z(12)%setup_boundary(UPPER, POLOIDAL, DIVERTOR) ! left divertor target
     call Z(11)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain
     call Z(12)%setup_boundary(UPPER, RADIAL,   VACUUM)   ! vacuum domain

     ! primary PFR
     call Z(7)%setup_mapping  (LOWER, RADIAL,   Z(13), 3)  ! connect to right primary PFR at interface I3
     call Z(13)%setup_boundary(LOWER, RADIAL,   VACUUM)    !
     call Z(8)%setup_mapping  (LOWER, RADIAL,   Z(14), 4)  ! connect to left primary PFR at interface I4
     call Z(14)%setup_boundary(LOWER, RADIAL,   VACUUM)    !
     call Z(13)%setup_boundary(LOWER, POLOIDAL, DIVERTOR)  ! right divertor target
     call Z(13)%setup_mapping (UPPER, POLOIDAL, Z(14), 4)  ! connect to left primary PFR at interface R4
     call Z(14)%setup_boundary(UPPER, POLOIDAL, DIVERTOR)  ! left divertor target

     ! secondary PFR
     call Z(10)%setup_mapping (LOWER, RADIAL,   Z(15),10)  ! connect to right secondary PFR at interface I10
     call Z(15)%setup_boundary(LOWER, RADIAL,   VACUUM)    !
     call Z(11)%setup_mapping (LOWER, RADIAL,   Z(16), 9)  ! connect to left secondary PFR at interface I9
     call Z(16)%setup_boundary(LOWER, RADIAL,   VACUUM)    !
     call Z(16)%setup_mapping (UPPER, POLOIDAL, Z(15), 8)  ! connect to left secondary PFR at interface R8
     call Z(16)%setup_boundary(LOWER, POLOIDAL, DIVERTOR)  ! right divertor target
     call Z(15)%setup_boundary(UPPER, POLOIDAL, DIVERTOR)  ! left divertor target


  ! connected double null (CDN)
  case(TOPO_CDN, TOPO_CDN1)
  end select
  call undefined_element_boundary_check(.false.)

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
  integer          :: ix, ierr, iSOL, iPFR, ipi, jx, k, orientation


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
  call initialize_module_mfs_mesh(C_guide)


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
     ! and set up interfaces between elements
     !call Iface(1)%set_curve(S0)
  write (6, 3000)
  iSOL = 0
  iPFR = 0
  ipi  = 0
  do ix=1,nX
     jx = connectX(ix)


     ! "core-interface"
     if (ix == 1  .or.  jx > 0) then
        write (6, 3010) ix
        call R(ix, DESCENT_CORE)%generateX(ix, DESCENT_CORE, LIMIT_PSIN, PsiN_in)
        call R(ix, DESCENT_CORE)%flip()
        call poloidal_interface(ipi+1)%set_curve(R(ix, DESCENT_CORE)%t_curve)
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
        call poloidal_interface(ipi+3)%set_curve(R(ix, ASCENT_LEFT)%t_curve)

        PsiN = get_PsiN(R(ix, ASCENT_LEFT)%boundary_node(boundary=UPPER))
        call R(ix, ASCENT_RIGHT)%generateX(ix, ASCENT_RIGHT, LIMIT_PSIN, PsiN)
        call poloidal_interface(ipi+2)%set_curve(R(ix, ASCENT_RIGHT)%t_curve)

     elseif (jx > ix) then
        ! connected X-points, left and right branch on individual flux surfaces
        iSOL = iSOL + 1
        write (6, 3021) ix, jx, d_SOL(iSOL)
        call R(ix, ASCENT_LEFT)%generateX(ix, ASCENT_LEFT, LIMIT_LENGTH, d_SOL(iSOL))
        call poloidal_interface(ipi+3)%set_curve(R(ix, ASCENT_LEFT)%t_curve)
        PsiN = get_PsiN(R(ix, ASCENT_LEFT)%boundary_node(boundary=UPPER))
        call R(jx, ASCENT_RIGHT)%generateX(jx, ASCENT_RIGHT, LIMIT_PSIN, PsiN)

        iSOL = iSOL + 1
        write (6, 3022) ix, jx, d_SOL(iSOL)
        call R(ix, ASCENT_RIGHT)%generateX(ix, ASCENT_RIGHT, LIMIT_LENGTH, d_SOL(iSOL))
        call poloidal_interface(ipi+2)%set_curve(R(ix, ASCENT_RIGHT)%t_curve)
        PsiN = get_PsiN(R(ix, ASCENT_RIGHT)%boundary_node(boundary=UPPER))
        call R(jx, ASCENT_LEFT)%generateX(jx, ASCENT_LEFT, LIMIT_PSIN, PsiN)

     elseif (jx == -ix) then
        ! outer SOL with left and right branch on individual flux surfaces
        iSOL = iSOL + 1
        write (6, 3023) ix, d_SOL(iSOL)
        call R(ix, ASCENT_LEFT)%generateX(ix, ASCENT_LEFT, LIMIT_LENGTH, d_SOL(iSOL))
        call poloidal_interface(ipi+3)%set_curve(R(ix, ASCENT_LEFT)%t_curve)

        iSOL = iSOL + 1
        write (6, 3024) ix, d_SOL(iSOL)
        call R(ix, ASCENT_RIGHT)%generateX(ix, ASCENT_RIGHT, LIMIT_LENGTH, d_SOL(iSOL))
        call poloidal_interface(ipi+2)%set_curve(R(ix, ASCENT_RIGHT)%t_curve)

     elseif (jx < 0) then
        ! this SOL's boundary is another separatrix
        PsiN = Xp(abs(jx))%PsiN()
        write (6, 3025) ix, abs(jx), PsiN
        call R(ix, ASCENT_LEFT)%generateX(ix, ASCENT_LEFT, LIMIT_PSIN, PsiN)
        call poloidal_interface(ipi+3)%set_curve(R(ix, ASCENT_LEFT)%t_curve)
        call R(ix, ASCENT_RIGHT)%generateX(ix, ASCENT_RIGHT, LIMIT_PSIN, PsiN)
        call poloidal_interface(ipi+2)%set_curve(R(ix, ASCENT_RIGHT)%t_curve)
     endif


     ! private flux region
     iPFR = iPFR + 1
     write (6, 3030) ix, d_PFR(iPFR)
     call R(ix, DESCENT_PFR)%generateX(ix, DESCENT_PFR, LIMIT_LENGTH, d_PFR(iPFR))
     call R(ix, DESCENT_PFR)%flip()
     call poloidal_interface(ipi+4)%set_curve(R(ix, DESCENT_PFR)%t_curve)


     ! plot paths
!     do k=1,4
!        call R(ix, k)%plot(filename='rpath_'//trim(str(ix))//'_'//trim(str(k))//'.plt')
!     enddo

     ipi = ipi + 4
  enddo
  write (6, *)

  do ipi=1,poloidal_interfaces
     call poloidal_interface(ipi)%C%plot(filename='R'//trim(str(ipi))//'.plt')
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


  integer :: ix, ix1, jx, iri, iconnect


  iri = 0
  do ix=1,nX
     jx = connectX(ix)

     ! connect back to same X-point
     if (jx == ix) then
        if (ix .ne. 1) then
           write (6, 9000) ix
           stop
        endif

        iri = iri + 1
        call radial_interface(iri)%set_curve(S0)
        call radial_interface(iri)%setup(ix, ix)

     ! all branches connect to divertor targets
     elseif (jx == -ix) then
        iconnect  = STRIKE_POINT
        ! are these "upstream" branches?
        do ix1=1,ix-1
           if (abs(connectX(ix1)) == ix) then
              iconnect = -ix1
              exit
           endif
        enddo
        iri = iri + 1
        if (iconnect == STRIKE_POINT) call radial_interface(iri)%set_curve(S(ix)%M1%t_curve)
        call radial_interface(iri)%setup(ix, iconnect)
        iri = iri + 1
        if (iconnect == STRIKE_POINT) call radial_interface(iri)%set_curve(S(ix)%M2%t_curve)
        call radial_interface(iri)%setup(iconnect, ix)

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
        iri = iri + 1
        call radial_interface(iri)%set_curve(S0R)
        call radial_interface(iri)%setup(ix, jx)

        ! left core interface
        iri = iri + 1
        call radial_interface(iri)%set_curve(S0L)
        call radial_interface(iri)%setup(jx, ix)

     ! nothing to be done here anymore
     else

     endif

     ! divertor branches
     iri = iri + 1
     call radial_interface(iri)%set_curve(S(ix)%M3%t_curve)
     call radial_interface(iri)%setup(STRIKE_POINT, ix)
     iri = iri + 1
     call radial_interface(iri)%set_curve(S(ix)%M4%t_curve)
     call radial_interface(iri)%setup(ix, STRIKE_POINT)
     ! add divertor branches for outer separatrix
     if (jx < -ix) then
        iri = iri + 1
        call radial_interface(iri)%setup(STRIKE_POINT, -ix)
        iri = iri + 1
        call radial_interface(iri)%setup(-ix, STRIKE_POINT)
     endif


  enddo


  !write (6, *) 'radial interfaces:'
  do iri=1,radial_interfaces
     !write (6, *) iri, radial_interface(iri)%inode(-1), radial_interface(iri)%inode(1)
     call radial_interface(iri)%C%plot(filename='I'//trim(str(iri))//'.plt')
  enddo

 9000 format('error: seconday X-point ', i0, ' connects back to itself!')
 9001 format('error: seconday X-point ', i0, ' does not connect back to primary one!')
 9002 format('error: seconday X-point ', i0, ' used as guiding point for separatrix ', i0, '!')
  end subroutine setup_interfaces
!=======================================================================



!=======================================================================
  subroutine initialize_layers(iblock)
  use fieldline_grid
  integer, intent(in) :: iblock

  integer, parameter :: COUNT_RUN = 1, SETUP_RUN = 2
  integer, dimension(:), allocatable :: markz, izl
  type(t_toroidal_discretization)    :: T
  character(len=1) :: cside(-1:1) = (/'R','C','L'/)
  real(real64) :: phi
  integer      :: i, il, iz, iz0, iz_map, idir, irun, nzl(-1:1)


  ! initialize block
  call load_local_resolution(iblock)
  phi = Block(iblock)%phi_base / 180.d0 * pi


  write (6, 1000)
  allocate (markz(nelement), izl(-nelement:nelement))

  do irun=COUNT_RUN,SETUP_RUN
  markz = 0
  il    = 0
  if (irun == SETUP_RUN) allocate(L(0:layers-1))
  layer_loop: do
     ! exit if no more elements are unmarked
     if (sum(markz) == nelement) exit

     ! start new layer
     il  = il + 1
     izl = 0;  nzl = 0

     ! find first unmarked element
     do iz=1,nelement
        if (markz(iz) == 0) exit
     enddo

     ! set base element in layer
     iz0 = iz;  markz(iz) = 1;  izl(0) = iz0

     ! scan in both poloidal directions
     dir_loop: do idir=-1,1,2
        ! start poloidal scan at base element
        iz = iz0
        poloidal_scan: do
           iz_map = Z(iz)%map_p(idir)
           ! poloidal scan in both directions finished when returning to base element
           if (iz_map == PERIODIC) exit dir_loop
           if (iz_map == iz0) exit dir_loop
           ! poloidal scan in this direction finished at divertor targets
           if (iz_map == DIVERTOR) exit

           iz = iz_map;  markz(iz) = 1
           nzl(idir) = nzl(idir) + 1;  izl(idir*nzl(idir)) = iz
        enddo poloidal_scan
     enddo dir_loop

     if (irun == SETUP_RUN) then
        ! now set up element indices in this layer
        call L(il-1)%initialize(nzl(-1)+1+nzl(1), izl(-nzl(-1):nzl(1)), nzl(-1)+1)

        ! set up resolution in each element
        call L(il-1)%setup_resolution(il-1)

        ! map poloidal resolution in radial direction
        call L(il-1)%map_poloidal_resolution()
     endif
  enddo layer_loop
  layers = il
  enddo


  ! initialize toroidal discretization
  call T%setup(nt, Block(iblock)%it_base, Block(iblock)%phi)
  do il=0,layers-1;  L(il)%T = T;  enddo
  do iz=1,nelement;  Z(iz)%T = T;  enddo


  ! initialize mesh for each domain element
  allocate (Mtmp(nelement))
  allocate (M(0:layers-1))
  do il=0,layers-1
     write (6, 1001) il
     do i=1,L(il)%nz
        iz = L(il)%iz(i)
        call Mtmp(iz)%initialize(Z(iz)%nr, Z(iz)%np, phi)
        write (6, 1002) iz, Z(iz)%np, Z(iz)%nr, Z(iz)%ipl, cside(Z(iz)%ipl_side)
     enddo
     call M(il)%initialize(L(il)%nr, L(il)%np, phi)
  enddo


  ! cleanup
  deallocate (markz, izl)

 1000 format(3x,'- Setting up radial layers:')
 1001 format(8x,'Layer ',i0)
 1002 format(8x,i3,': ',i5,' x ',i3,5x,'(',i0,a1,')')
  end subroutine initialize_layers
!=======================================================================



!=======================================================================
  subroutine generate_base_mesh(iblock)
  use fieldline_grid
  use mesh_spacing
  use inner_boundary
  use string
  integer, intent(in) :: iblock

  integer, dimension(:), allocatable :: markz
  type(t_spacing) :: Sp, SpL, SpR, Sr
  integer         :: i, il, il0, iz, iz0, iz_map, iside


  ! initialize block
  call initialize_layers(iblock)


  ! generate core-interface
  il = 0
  if (connectX(1) == 1) then
     ! single element
     call Sp%init(poloidal_spacing(il))
     call Mtmp(1)%setup_boundary_nodes(UPPER, RADIAL, S0, Sp)
  else
     ! left and right elements
     call SpR%init(poloidal_spacing_R(il))
     call Mtmp(1)%setup_boundary_nodes(UPPER, RADIAL, S0R, SpR)

     call SpL%init(poloidal_spacing_L(il))
     call Mtmp(2)%setup_boundary_nodes(UPPER, RADIAL, S0L, SpL)
  endif


  ! main loop: generate mesh for each layer
  write (6, *)
  allocate (markz(nelement));  markz = 0
  do il=0,layers-1
     ! base element index
     iz0 = L(il)%iz(L(il)%i0)

     ! set up radial spacings for this layer
     call Sr%init(radial_spacing(il))

     ! generate mesh for this layer
     call generate_layer(il, iz0, iblock, Sr)

     ! map radial interface to next element
     do i=1,L(il)%nz
        iz = L(il)%iz(i)
        markz(iz) = 1 ! mark elements in this layer

        do iside=-1,1,2
           iz_map = Z(iz)%map_r(iside)
           ! map to another element?
           if (iz_map < 0) cycle

           ! map poloidal resolution
           if (markz(iz_map) == 0) then
              call Mtmp(iz)%connect_to(Mtmp(iz_map), RADIAL, iside)

              ! status of radial interface
              !write (6, *) 'radial interface in element ', iz, ' side ', iside, ' is ', Z(iz)%rad_bound(iside)
           endif
        enddo
     enddo
  enddo


  ! write output files
  il0 = iblock * layers
  do il=0,layers-1
     iz = il0 + il
     call write_base_grid(M(il)%t_grid, iz)
  enddo


  ! cleanup
  deallocate (M, Mtmp, markz)

  end subroutine generate_base_mesh
!=======================================================================



!=======================================================================
  subroutine generate_layer(il, iz0, iblock, Sr)
  use run_control, only: Debug
  use fieldline_grid, only: poloidal_spacing, poloidal_spacing_L, poloidal_spacing_R
  use mesh_spacing
  use string
  integer,         intent(in) :: il, iz0, iblock
  type(t_spacing), intent(in) :: Sr

  type(t_spacing) :: Sp
  integer :: i, idir, irside, iri, ipside, ipi, ip0, iz, iz_map, npz(-1:1)


  write (6, 1000) il, iz0
  write (6, *) 'reference discretization at:'


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
  write (6, *) 'radial boundary ', irside, ' which is interface ', iri
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
  write (6, *) 'poloidal boundary ', ipside, ' which is interface ', ipi
  if (ipi == UNDEFINED) then
     write (6, 9000) il, iz0
     write (6, 9003)
     stop
  endif


  ! initialize radial discretization in layer
  call Mtmp(iz0)%setup_boundary_nodes(ipside, POLOIDAL, poloidal_interface(ipi)%C, Sr)


  ! generate mesh in base element
  ! no Sp needed in base element
  call Z(iz0)%generate_mesh(Mtmp(iz0), irside, ipside, iblock, Sr, Sp)
  if (Debug) call Mtmp(iz0)%plot_mesh('Mtmp'//trim(str(iz0))//'.plt')


  ! scan through poloidal elements in this layer
  npz    = 0
  npz(0) = 1
  idir_loop: do idir=1,-1,-2
     iz = iz0
     poloidal_scan: do
        iz_map = Z(iz)%map_p(idir)
        ! 1. poloidal boundary of layer?
        ! 1.1 periodic boundaries: connect element back to itself
        if (iz_map == PERIODIC) then
           call Mtmp(iz)%connect_to(Mtmp(iz), POLOIDAL, LOWER_TO_UPPER)
           exit idir_loop
        endif
        ! 1.2 back to initial/base element
        if (iz_map == iz0) then
           call Mtmp(iz)%connect_to(Mtmp(iz0), POLOIDAL, LOWER_TO_UPPER)
           exit idir_loop
        endif
        ! 1.3 divertor targets
        if (iz_map == DIVERTOR) exit

        ! 2. connect mesh to next element
        write (6, *) 'connect element ', iz, ' to ', iz_map
        call Mtmp(iz)%connect_to(Mtmp(iz_map), POLOIDAL, idir)
        iz        = iz_map

        ! 3. generate mesh in next element
        select case(Z(iz)%ipl_side)
        case(LEFT)
           call Sp%init(poloidal_spacing_L(Z(iz)%ipl))
        case(CENTER)
           call Sp%init(poloidal_spacing(Z(iz)%ipl))
        case(RIGHT)
           call Sp%init(poloidal_spacing_R(Z(iz)%ipl))
        end select
        call Z(iz)%generate_mesh(Mtmp(iz), irside, ipside, iblock, Sr, Sp)
        if (Debug) call Mtmp(iz)%plot_mesh('Mtmp'//trim(str(iz))//'.plt')
        npz(idir) = npz(idir) + 1
     enddo poloidal_scan
  enddo idir_loop
  write (6, *)


  ! merge elements
  ip0 = 0
  do i=1,L(il)%nz
     iz  = L(il)%iz(i)
     call M(il)%copy(0, ip0, Mtmp(iz))
     ip0 = ip0 + Z(iz)%np
  enddo

 1000 format('Generate layer ', i0, ' from base element ', i0)
 9000 format('error in generate_layer for il, iz0 = ', i0, ', ', i0)
 9001 format('undefined reference nodes on radial boundary!')
 9002 format('undefined radial interface!')
 9003 format('undefined poloidal interface!')
  end subroutine generate_layer
!=======================================================================

end module base_mesh
