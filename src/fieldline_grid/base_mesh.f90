module base_mesh
  use iso_fortran_env
  use grid
  use curve2D
  use separatrix
  use xpaths
  implicit none
  private

  integer, parameter, public :: &
     RADIAL         = 1, &
     POLOIDAL       = 2, &
     LOWER_TO_UPPER = 1, &
     UPPER_TO_LOWER = 2


  ! orthogonal, flux surface aligned mesh with strike point adjustment
  type, extends(t_grid), public :: t_base_mesh
     integer :: nr, np
     integer :: &
        ir0 = -1, & ! reference radial surface for mesh generation
        ip0 = -1    ! reference poloidal surface for mesh generation
     contains
     procedure :: initialize
     procedure :: connect_to
     !procedure :: connect_partial_to
     procedure :: setup_boundary_nodes
     procedure :: make_orthogonal_grid
     procedure :: make_divertor_grid
     ! procedure :: merge(M1, M2, ...)
  end type t_base_mesh


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


  public :: setup_geometry
  !public :: setup_topology
  public :: generate_base_mesh

  contains
!=======================================================================



!=======================================================================
  subroutine initialize(this, nr, np, phi)
  class(t_base_mesh)       :: this
  integer, intent(in)      :: nr, np
  real(real64), intent(in) :: phi


  this%nr = nr
  this%np = np
  call this%new(CYLINDRICAL, MESH_2D, 3, nr+1, np+1, fixed_coord_value=phi)

  end subroutine initialize
!=======================================================================



!=======================================================================
! Connect this bock to next block (M)
! direction:	boundary direction (radial or poloidal)
! side:		lower-to-upper or upper-to-lower boundary
!=======================================================================
  subroutine connect_to(this, M, direction, side)
  class(t_base_mesh)  :: this, M
  integer, intent(in) :: direction, side

  real(real64), dimension(:,:,:), pointer :: M1, M2
  integer :: ir, ir1, ir2, ip, ip1, ip2


  M1 => this%mesh
  M2 => M%mesh



  select case(direction)
  ! CONNECT IN RADIAL DIRECTION
  case(RADIAL)
     ! 1. check poloidal resolution of interface
     if (this%np .ne. M%np) then
        write (6, 9000)
        write (6, 9010) this%np, M%np
        stop
     endif

     ! 2. setup radial indices
     select case(side)
     case(LOWER_TO_UPPER)
        ir1 = 0
        ir2 = M%nr
     case(UPPER_TO_LOWER)
        ir1 = this%nr
        ir2 = 0
     case default
        write (6, 9000)
        write (6, 9002) side
        stop
     end select

     ! 3. copy boundary nodes
     do ip=0,this%np
        M2(ir2, ip, :) = M1(ir1, ip, :)
     enddo
     M%ir0 = ir2

  ! CONNECT IN POLOIDAL DIRECTION
  case(POLOIDAL)
     ! 1. check radial resolution of interface
     if (this%nr .ne. M%nr) then
        write (6, 9000)
        write (6, 9011) this%nr, M%nr
        stop
     endif

     ! 2. setup poloidal indices
     select case(side)
     case(LOWER_TO_UPPER)
        ip1 = 0
        ip2 = M%np
     case(UPPER_TO_LOWER)
        ip1 = this%np
        ip2 = 0
     case default
        write (6, 9000)
        write (6, 9002) side
        stop
     end select

     ! 3. copy boundary nodes
     do ir=0,this%nr
        M2(ir, ip2, :) = M1(ir, ip2, :)
     enddo
     M%ip0 = ip2

  case default
     write (6, 9000)
     write (6, 9001) direction
     stop
  end select

 9000 format('error in t_base_mesh%connect_to')
 9010 format('invalid poloidal resolution in neighboring blocks: ', i0, 2x, i0)
 9011 format('invalid radial resolution in neighboring blocks: ', i0, 2x, i0)
 9001 format('invalid boundary direction: ', i0)
 9002 format('invalid boundary side id: ', i0)
  end subroutine connect_to
!=======================================================================



!=======================================================================
! Setup boundary nodes in mesh
! boundary_type:	radial or poloidal block boundary
! boundary_side:	lower or upper boundary

! C_boundary:	boundary definition
! spacings:	spacing function for nodes on boundary
!=======================================================================
  subroutine setup_boundary_nodes(this, boundary_type, boundary_side, C_boundary, spacings, i1, i2)
  use mesh_spacing
  use curve2D
  class(t_base_mesh)          :: this
  integer,         intent(in) :: boundary_type, boundary_side
  type(t_curve),   intent(in) :: C_boundary
  type(t_spacing), intent(in) :: spacings
  integer,         intent(in), optional :: i1, i2

  real(real64), dimension(:,:,:), pointer :: M
  real(real64) :: tau, x(2)
  integer      :: ir, ip, i11, i22, idir


  M => this%mesh
  call C_boundary%setup_length_sampling()


  select case(boundary_side)
  case(LOWER)
     ir = 0
     ip = 0
  case(UPPER)
     ir = this%nr
     ip = this%np
  case default
     write (6, *) 'error in t_base_mesh%setup_boundary_nodes:'
     write (6, *) 'invalid boundary side ', boundary_side
     stop
  end select


  idir = 1
  select case(boundary_type)
  case(RADIAL)
     this%ir0 = ir
     i11      = 0
     i22      = this%np
     if (present(i1)) i11 = i1
     if (present(i2)) i22 = i2
     if (i22 < i11) idir = -1
     do ip=i11,i22,idir
        tau = spacings%node(ip-i11, i22-i11)
        call C_boundary%sample_at(tau, x)
        M(ir, ip, :) = x
     enddo
  case(POLOIDAL)
     this%ip0 = ip
     i11      = 0
     i22      = this%nr
     if (present(i1)) i11 = i1
     if (present(i2)) i22 = i2
     if (i22 < i11) idir = -1
     do ir=i11,i22,idir
        tau = spacings%node(ir-i11, i22-i11)
        call C_boundary%sample_at(tau, x)
        M(ir, ip, :) = x
     enddo
  end select

  end subroutine setup_boundary_nodes
!=======================================================================



!=======================================================================
! Make (quasi) orthogonal grid
!=======================================================================
  subroutine make_orthogonal_grid(this, rrange, prange, periodic)
  use equilibrium, only: get_PsiN
  class(t_base_mesh)  :: this
  integer, intent(in), optional :: rrange(2), prange(2)
  logical, intent(in), optional :: periodic

  real(real64), dimension(:,:,:), pointer :: M
  type(t_xpath) :: R
  real(real64)  :: PsiN(0:this%nr), x(2), PsiN_final
  integer       :: ir, ir0, ir1, ir2, ir_final, ip, ip0, ip1, ip2, ipp, direction


  M => this%mesh


  ! set poloidal range
  ip0 = this%ip0
  if (ip0 == 0) then
     ip1 = 1
     ip2 = this%np
     ipp = this%np
  elseif (ip0 == this%np) then
     ip1 = 0
     ip2 = this%np - 1
     ipp = 0
  else
     write (6, 9000)
     write (6, 9001) ip0
     stop
  endif
  if (present(periodic)) then
  if (periodic) then
     ip1 = 1
     ip2 = this%np - 1
  endif
  endif
  if (present(prange)) then
     ip1 = prange(1)
     ip2 = prange(2)
  endif
  write (6, *) 'poloidal range: ', ip0, ip1, ip2


  ! set radial range
  ir0 = this%ir0
  if (ir0 == 0) then
     ir1        = 1
     ir2        = this%nr
     direction  = ASCENT
  elseif (ir0 == this%nr) then
     ir1        = 0
     ir2        = this%nr - 1
     direction  = DESCENT
  else
     write (6, 9000)
     write (6, 9002) ir0
     stop
  endif
  if (present(rrange)) then
     ir1 = rrange(1)
     ir2 = rrange(2)
  endif
  write (6, *) 'radial range: ', ir0, ir1, ir2


  ! setup reference PsiN values
  do ir=ir1,ir2
     x        = M(ir,ip0,:)
     PsiN(ir) = get_PsiN(x)
  enddo

  select case(direction)
  case(ASCENT)
     ir_final = ir2
  case(DESCENT)
     ir_final = ir1
  end select
  PsiN_final = PsiN(ir_final)


  ! set up nodes in poloidal range ip1->ip2 and radial range ir1->ir2
  do ip=ip1,ip2
     call R%generate(M(ir0,ip,:), direction, LIMIT_PSIN, PsiN_final, sampling=SAMPLE_PSIN)

     do ir=ir1,ir2
        call R%sample_at_PsiN(PsiN(ir), x, enforce_boundary=.true.)
        M(ir,ip,:) = x
     enddo
  enddo


  ! set up periodic boundaries
  if (present(periodic)) then
  if (periodic) then
     M(:,ipp,:) = M(:,ip0,:)
  endif
  endif

 9000 format('error in t_base_mesh%make_orthogonal_grid')
 9001 format('invalid poloidal reference index ', i0)
 9002 format('invalid radial reference index ', i0)
  end subroutine make_orthogonal_grid
!=======================================================================



!=======================================================================
! Generate nodes in the base plane so that field lines connect to the
! strike point x0
!=======================================================================
  subroutine align_strike_points(x0, Z, nsub, M)
  use fieldline_grid, only: t_zone
  use fieldline
  real(real64), intent(in)  :: x0(2)
  type(t_zone), intent(in)  :: Z
  integer,      intent(in)  :: nsub
  real(real64), intent(out) :: M(0:Z%nt*nsub, 2)

  type(t_fieldline) :: F
  real(real64)      :: ts, y0(3), y1(3), Dphi
  integer           :: tm, tc, idir, it, it0, it_sub, it_end, ierr


  ! set parameters for field line tracing
  ! (taken from subroutine trace_nodes)
  ts = pi2 / 3600.d0
  tm = NM_AdamsBashforth4
  tc = FL_ANGLE



  ! set initial point
  y0(1:2)    = x0
  y0(3)      = Z%phi(Z%it_base)
  it0        = Z%it_base * nsub
  M(it0,1:2) = y0(1:2)

  do idir=-1,1,2
     call F%init(y0, idir*ts, tm, tc)

     ! trace from base location to zone boundaries
     it_end = 0
     if (idir > 0) it_end = Z%nt
     do it=Z%it_base+idir,it_end, idir
        Dphi = abs(Z%phi(it) - Z%phi(it-idir)) / nsub / 180.d0 * pi

        ! sub-resolution
        do it_sub=1,nsub
           call F%trace_Dphi(Dphi, .false., y1, ierr)
           if (ierr .ne. 0) then
              write (6, 9000) ierr
              stop
           endif

           M((it-idir)*nsub + it_sub*idir,1:2) = y1(1:2)
        enddo
     enddo
  enddo
 9000 format('error in subroutine align_strike_point: ',//, &
             't_fieldline%trace_Dphi returned ierr = ', i0)
  end subroutine align_strike_points
!=======================================================================



!=======================================================================
! Generate grid in divertor legs
!=======================================================================
!  subroutine make_divertor_grid(this, R, Rside, Sr, P, Pside, Sp, Z, ierr)
  subroutine make_divertor_grid(this, Rside, ip0, Z, ierr)
  use fieldline_grid, only: t_zone
  use equilibrium, only: get_PsiN, Ip_sign, Bt_sign
  use curve2D
  use mesh_spacing
  use flux_surface_2D
  class(t_base_mesh)           :: this
!  type(t_curve),   intent(in)  :: R, P
!  integer,         intent(in)  :: Rside, Pside
!  type(t_spacing), intent(in)  :: Sr, Sp
  integer,         intent(in)  :: Rside, ip0
  type(t_zone),    intent(in)  :: Z
  integer,         intent(out) :: ierr

  integer, parameter :: nsub = 2

  type(t_flux_surface_2D) :: F

  real(real64), dimension(:,:,:), pointer :: M
  real(real64), dimension(:,:), allocatable :: MSP
  real(real64)  :: PsiN(0:this%nr), x(2), PsiN_final, tau, L0, L
  integer       :: ir, ir0, ir1, ir2, ip, ips, ip1, ip2, dir
  integer       :: it, its, it_start, it_end, dirT, downstream


  ! check intersection of R with guiding surface
!  if (R%intersect_curve(C_guide, x, tau)) then
!     L    = R%length() * tau
!     write (6, 9000) L
!     ierr = 1
!     return
!  endif


  ! check resolution
  if (Z%nt * nsub > this%np) then
     write (6, 9010) this%np, Z%nt, nsub
     stop
  endif


  ! setup poloidal boundary with reference nodes for flux surfaces
  !call this%setup_boundary_nodes(POLOIDAL, Rside, R, Sr)
  select case(Rside)
  case(LOWER)
     dir        = RIGHT_HANDED
     downstream = 1
     ips        = this%np - Z%nt*nsub
  case(UPPER)
     dir        = LEFT_HANDED
     downstream = -1
     ips        = Z%nt*nsub
     ! ip_ds    = 0
  end select


  ! find downstream direction
  dirT       = dir * Ip_sign * Bt_sign
  it_start   = Z%nt * nsub
  it_end     = 0
  if (dirT > 0) then
     it_start = 0
     it_end   = Z%nt * nsub
  endif


  M => this%mesh
  allocate (MSP(0:Z%nt*nsub, 2))
  do ir=0,this%nr
     x = M(ir,ip0,:)

     ! generate flux surface from "upstream" location x to target
     call F%generate(x, dir, Trace_Step=0.1d0)
     call F%plot(filename='F.plt', append=.true.)

     ! generate nodes from which field lines connect to strike point x
     select case(dir)
     case(RIGHT_HANDED)
        x = F%x(F%n_seg, :)
     case(LEFT_HANDED)
        x = F%x(0,:)
     end select
     !write (99, *) x
     call align_strike_points(x, Z, nsub, MSP)


     ! setup downstream strike point nodes
     do it=it_start,it_end,dirT
        ip = ips + abs(it - it_end)*dir
        M(ir,ip,:) = MSP(it,:)
     enddo


     ! length on flux surface for strike point mesh
     L0 = F%length()
     L  = 0.d0
     do it=Z%it_base*nsub+dirT,it_end,dirT
        L = L + sqrt(sum(  (MSP(it-dirT,:)-MSP(it,:))**2))
        !write (92, *) MSP(it, :)
     enddo

     ! aligned strike point mesh extends beyond upstream reference location?
     if (L > L0) then
        write (6, 9020) ir
        stop
     endif
     ! interpolate nodes on flux surface between upstream orthogonal grid nodes
     ! and downstream strike point nodes
     call F%setup_length_sampling()
     select case(dir)
     case(RIGHT_HANDED)
        F%w = F%w / (1.d0 - L/L0)
     case(LEFT_HANDED)
        F%w = (F%w - L/L0) / (1.d0 - L/L0)
     end select

     do ip=ip0+dir,ips-dir,dir
        tau = 1.d0 * (ip - ip0) / (ips - ip0)
        if (dir == LEFT_HANDED) tau = 1.d0 - tau
        call F%sample_at(tau, x)
        M(ir,ip,:) = x
        !write (93, *) x, tau
     enddo
  enddo

  deallocate (MSP)
 9000 format('ERROR: reference path for radial discretization crosses guiding surface!', &
             'intersection at L = ', f0.3)
 9010 format('ERROR: poloidal grid resolution too small!'//,&
             i0, ' < ', i0, ' x ', i0)
 9020 format('ERROR: aligned strike point mesh extends beyond upstream reference point ', &
             'at radial index ', i0)
  end subroutine make_divertor_grid
!=======================================================================



!=======================================================================
! Set up geometry of computational domain:
! Magnetic axis, X-points, separatrix(ces) and radial paths from X-points
! SHOULD THIS BE MOVED TO MODULE fieldline_grid?
!=======================================================================
  subroutine setup_geometry(nX_, connectX_)
  use fieldline_grid, only: guiding_surface, d_SOL, d_PFR
  use boundary,       only: S_axi, n_axi
  use equilibrium,    only: get_magnetic_axis, get_poloidal_angle, get_PsiN, Xp
  use math
  use separatrix
  use inner_boundary
  use string

  integer, intent(in) :: nX_, connectX_(nX_)

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
  nX = nX_
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
  allocate (connectX(nX), S(nX), R(nX,4))
  connectX = connectX_
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
!=======================================================================
  subroutine generate_base_mesh(iblock)
  use fieldline_grid
  use mesh_spacing
  integer, intent(in) :: iblock

  type(t_base_mesh)   :: M(0:layers-1)
  type(t_spacing)     :: Sp, SpL, SpR, Sr

  real(real64) :: phi
  integer      :: il, iz, iz0


  ! initialize block
  iz0 = iblock * layers
  call load_local_resolution(iblock)

  phi = Block(iblock)%phi_base / 180.d0 * pi
  do il=0,layers-1
     call M(il)%initialize(nr(il), np(il), phi)
  enddo


  ! generate core-interface
  if (connectX(1) == 1) then
     ! single zone
     call Sp%init(poloidal_spacing(0))
     call M(0)%setup_boundary_nodes(RADIAL, UPPER, S0, Sp)
  else
     ! left and right sub-zones
     call SpL%init(poloidal_spacing(0))
     call SpR%init(poloidal_spacing(1))
     call M(0)%setup_boundary_nodes(RADIAL, UPPER, S0R, SpR, 0,      npR(0))
     call M(0)%setup_boundary_nodes(RADIAL, UPPER, S0L, SpL, npR(0), np(0) )
  endif


  ! generate "closed" domain
  call Sr%init(radial_spacing(0))
  if (connectX(1) == 1  .or. connectX(1) < 0) then
     call M(0)%setup_boundary_nodes(POLOIDAL, LOWER, R(1,DESCENT_CORE)%t_curve, Sr, nr(0), 1)
     call M(0)%make_orthogonal_grid(periodic=.true., rrange=(/2+n_interpolate, nr(0)-1/))
  else
  endif



  ! write output files
  do il=0,layers-1
     iz = iz0 + il
     call write_base_grid(M(il)%t_grid, iz)
  enddo

  end subroutine generate_base_mesh
!=======================================================================

end module base_mesh
