! Magnetic Flux Surface mesh
module mfs_mesh
  use iso_fortran_env
  use grid
  use curve2D
  use separatrix
  use xpaths
  use mesh_interface
  implicit none
  private

  integer, parameter, public :: &
     LOWER_TO_UPPER = LOWER, &
     UPPER_TO_LOWER = UPPER

  integer, parameter, public :: &
     RADIAL         = 1, &
     POLOIDAL       = 2, &
     AUTOMATIC      = -1024



  ! orthogonal, flux surface aligned mesh with strike point adjustment
  type, extends(t_grid), public :: t_mfs_mesh
     integer :: nr, np
     integer :: &
        ir0 = -1, & ! reference radial surface for mesh generation
        ip0 = -1    ! reference poloidal surface for mesh generation

     contains
     procedure :: initialize
     procedure :: connect_to
     procedure :: setup_boundary_nodes
     procedure :: make_orthogonal_grid
     procedure :: make_interpolated_mesh
     procedure :: make_interpolated_submesh
     procedure :: make_divertor_grid
     procedure :: copy
  end type t_mfs_mesh


  type(t_curve) :: C_guide

  public :: initialize_module_mfs_mesh

  contains
!=======================================================================



!=======================================================================
  subroutine initialize(this, nr, np, phi)
  class(t_mfs_mesh)        :: this
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
  class(t_mfs_mesh)   :: this, M
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
        M2(ir, ip2, :) = M1(ir, ip1, :)
     enddo
     M%ip0 = ip2

  case default
     write (6, 9000)
     write (6, 9001) direction
     stop
  end select

 9000 format('error in t_mfs_mesh%connect_to')
 9010 format('invalid poloidal resolution in neighboring blocks: ', i0, 2x, i0)
 9011 format('invalid radial resolution in neighboring blocks: ', i0, 2x, i0)
 9001 format('invalid boundary direction: ', i0)
 9002 format('invalid boundary side id: ', i0)
  end subroutine connect_to
!=======================================================================



!=======================================================================
! Setup boundary nodes in mesh
! boundary_side:	lower or upper boundary
! boundary_type:	radial or poloidal block boundary

! C_boundary:	boundary definition
! spacings:	spacing function for nodes on boundary
!=======================================================================
  subroutine setup_boundary_nodes(this, boundary_side, boundary_type, C_boundary, spacings, i1, i2, debug)
  use mesh_spacing
  use curve2D
  class(t_mfs_mesh)           :: this
  integer,         intent(in) :: boundary_side, boundary_type
  type(t_curve),   intent(in) :: C_boundary
  type(t_spacing), intent(in) :: spacings
  integer,         intent(in), optional :: i1, i2
  logical,         intent(in), optional :: debug

  real(real64), dimension(:,:,:), pointer :: M
  real(real64) :: tau, x(2)
  integer      :: ir, ip, i11, i22, idir
  logical      :: screen_output


  screen_output = .false.
  if (present(debug)) then
     screen_output = debug
  endif


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
        if (screen_output) write (6, *) x
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
        if (boundary_side == UPPER) tau = 1.d0 - spacings%node(i22-ir, i22-i11)
        call C_boundary%sample_at(tau, x)
        M(ir, ip, :) = x
        if (screen_output) write (6, *) x
     enddo
  end select

  end subroutine setup_boundary_nodes
!=======================================================================



!=======================================================================
! Make (quasi) orthogonal grid
!=======================================================================
  subroutine make_orthogonal_grid(this, rrange, prange, periodic, addX, debug)
  use equilibrium, only: get_PsiN
  class(t_mfs_mesh)             :: this
  integer, intent(in), optional :: rrange(2), prange(2), addX(2)
  logical, intent(in), optional :: periodic, debug

  real(real64), dimension(:,:,:), pointer :: M
  type(t_xpath) :: R
  real(real64)  :: PsiN(0:this%nr), x(2), PsiN_final
  integer       :: ir, ir0, ir1, ir2, ir_final, ip, ip0, ip1, ip2, ipp, direction, ix, ipx
  logical       :: screen_output


  screen_output = .false.
  if (present(debug)) then
     screen_output = debug
  endif


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


  ! setup reference PsiN values
  do ir=ir1,ir2
     x        = M(ir,ip0,:)
     PsiN(ir) = get_PsiN(x)
     if (screen_output) write (6, *) ir, x, PsiN(ir)
  enddo

  select case(direction)
  case(ASCENT)
     ir_final = ir2
  case(DESCENT)
     ir_final = ir1
  end select
  PsiN_final = PsiN(ir_final)


  ! additional X-point to take into account
  ipx = -1
  if (present(addX)) then
     ix  = addX(1)
     ipx = addX(2)

     ! set poloidal index to opposite boundary from reference boundary ip0
     if (ipx == AUTOMATIC) ipx = (1 - ip0 / this%np) * this%np
  endif


  ! set up nodes in poloidal range ip1->ip2 and radial range ir1->ir2
  write (6, 1000) ir0, ir1, ir2, ip0, ip1, ip2
  do ip=ip1,ip2
     call progress(ip-ip1, ip2-ip1, 0.1d0)
     if (ip == ipx) then
        if (screen_output) write (6, *) ip, ix
        call R%generateX(ix,         direction, LIMIT_PSIN, PsiN_final, sampling=SAMPLE_PSIN)
     else
        if (screen_output) write (6, *) ip, M(ir0,ip,:)
        call R%generate(M(ir0,ip,:), direction, LIMIT_PSIN, PsiN_final, sampling=SAMPLE_PSIN)
     endif

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


 1000 format(8x,'generate orthogonal mesh: (',i0,': ',i0,' -> ',i0,') x (',i0,': ',i0,' -> ',i0,')')
 9000 format('error in t_mfs_mesh%make_orthogonal_grid')
 9001 format('invalid poloidal reference index ', i0)
 9002 format('invalid radial reference index ', i0)
  end subroutine make_orthogonal_grid
!=======================================================================



!=======================================================================
! Generate interpolated mesh to inner simulation boundary (2 -> ir2)
!=======================================================================
  subroutine make_interpolated_mesh(this, ir2, Sr, C, PsiN1_max, prange)
  use equilibrium, only: get_PsiN, Xp
  use mesh_spacing
  class(t_mfs_mesh)           :: this
  integer,         intent(in) :: ir2
  type(t_spacing), intent(in) :: Sr
  type(t_curve),   intent(in) :: C(0:1)
  real(real64),    intent(in) :: PsiN1_max
  integer,         intent(in), optional :: prange(2)

  real(real64), dimension(:,:,:), pointer :: M
  type(t_xpath) :: Rtmp
  real(real64)  :: PsiN2, eta, x(2)
  integer       :: i, ir1, ip1, ip2, j, nr


  ! set defaults for optional input
  ip1 = 0
  ip2 = this%np
  if (present(prange)) then
     ip1 = prange(1)
     ip2 = prange(2)
  endif


  M  => this%mesh
  nr  = this%nr
  ir1 = 1
  write (6, 1001) ir1+1, ir2-1

  ! sanity check
  PsiN2 = get_PsiN(M(ir2,0,:)) ! radial location of innermost unperturbed flux surface
  if (PsiN2 < PsiN1_max) then
     write (6, 9000) ir2, PsiN2
     write (6, 9001) PsiN1_max
     stop
  endif

  do j=ip1,ip2
     call progress(j-ip1, ip2-ip1, 0.1d0)
     x = M(ir2,j,:)
     call Rtmp%generate(x, DESCENT_CORE, LIMIT_CURVE, PsiN1_max, C_limit=C(ir1), sampling=SAMPLE_LENGTH)

     ! interpolated surfaces
     do i=ir1,ir2-1
        eta = 1.d0 - Sr%node(i-1, nr-1) / Sr%node(ir2-1, nr-1)
        call Rtmp%sample_at(eta, x)

        M(i,j,:) = x
     enddo

     ! innermost surface
     x = M(ir1,j,:)
     call Rtmp%generate(x, DESCENT_CORE, LIMIT_CURVE, -1.d0, C_limit=C(0))
     M(0,j,:) = Rtmp%boundary_node(UPPER)
  enddo

 1001 format(8x,'interpolating from inner boundary to 1st unperturbed flux surface: ', &
             i0, ' -> ', i0)
 9000 format('error: last unperturbed flux surface at radial index ', i0, ' is at PsiN = ', &
             f0.3, ' but it must be completely outside of inner simulation boundary!'// &
             'try using a larger n_interpolate!')
 9001 format('outer most point on inner simulation boundary is at PsiN = ', f0.3)
  end subroutine make_interpolated_mesh
!=======================================================================



!=======================================================================
  subroutine make_interpolated_submesh(this, rrange, prange)
  use flux_surface_2D
  class(t_mfs_mesh)           :: this
  integer,         intent(in) :: rrange(2), prange(2)

  real(real64), dimension(:,:,:), pointer :: M
  type(t_curve)           :: C_limit
  type(t_flux_surface_2D) :: F
  real(real64) :: x(2), s, x1(2), x2(2)
  integer      :: ir, ip


  M => this%mesh

  ! set up upper poloidal boundary
  call C_limit%new(rrange(2)-rrange(1))
  do ir=rrange(1),rrange(2)
     C_limit%x(ir-rrange(1), :) = M(ir,prange(2),:)
  enddo

  do ir=rrange(1),rrange(2)
!     x = M(ir,prange(1),:)
!     call F%generate(x, RIGHT_HANDED, AltSurf=C_limit)
!     call F%setup_length_sampling()
!     do ip=prange(1)+1,prange(2)-1
!        s = 1.d0 * (ip-prange(1)) / (prange(2)-prange(1))
!        call F%sample_at(s, x)
!        M(ir,ip,:) = x
!     enddo
     x1 = M(ir,prange(1),:)
     x2 = M(ir,prange(2),:)
     do ip=prange(1)+1,prange(2)-1
        s = 1.d0 * (ip-prange(1)) / (prange(2)-prange(1))
        x = x1 + s * (x2-x1)
        M(ir,ip,:) = x
     enddo
  enddo

  end subroutine make_interpolated_submesh
!=======================================================================



!=======================================================================
! Generate nodes in the base plane so that field lines connect to the
! strike point x0
!
! Z:     toroidal discretization used for strike point adjustment
!=======================================================================
  subroutine align_strike_points(x0, Z, M)
  use fieldline_grid, only: t_toroidal_discretization
  use fieldline
  real(real64), intent(in)  :: x0(2)
  type(t_toroidal_discretization), intent(in)  :: Z
  real(real64), intent(out) :: M(0:Z%nt, 2)

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
  it0        = Z%it_base
  M(it0,1:2) = y0(1:2)

  do idir=-1,1,2
     call F%init(y0, idir*ts, tm, tc)

     ! trace from base location to zone boundaries
     it_end = 0
     if (idir > 0) it_end = Z%nt
     do it=Z%it_base+idir,it_end, idir
        Dphi = abs(Z%phi(it) - Z%phi(it-idir)) / 180.d0 * pi

        call F%trace_Dphi(Dphi, .false., y1, ierr)
        if (ierr .ne. 0) then
           write (6, 9000) ierr
           write (6, *) 'reference point: ', y0
           write (6, *) 'it, idir = ', it, idir
           write (6, *) 'it_end   = ', it_end
           write (6, *) 'it_base  = ', Z%it_base
           write (6, *) 'Dphi[deg]= ', Dphi / pi * 180.d0
           write (6, *) 'present location: ', y1
           stop
        endif

        M(it,1:2) = y1(1:2)
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
  use fieldline_grid, only: t_toroidal_discretization, np_sub_divertor
  use equilibrium, only: get_PsiN, Ip_sign, Bt_sign
  use curve2D
  use mesh_spacing
  use flux_surface_2D
  class(t_mfs_mesh)            :: this
!  type(t_curve),   intent(in)  :: R, P
!  integer,         intent(in)  :: Rside, Pside
!  type(t_spacing), intent(in)  :: Sr, Sp
  integer,         intent(in)  :: Rside, ip0
  type(t_toroidal_discretization),    intent(in)  :: Z
  integer,         intent(out) :: ierr

  integer, parameter :: nskip = 0, nextend = 1, iu_err = 66

  type(t_flux_surface_2D) :: F
  type(t_toroidal_discretization) :: TSP

  real(real64), dimension(:,:,:), pointer :: M
  real(real64), dimension(:,:), allocatable :: MSP
  real(real64)  :: PsiN(0:this%nr), x(2), PsiN_final, tau, L0, L, dphi
  integer       :: ir, ir0, ir1, ir2, ip, ips, ip1, ip2, dir, iextend, np_SP
  integer       :: it, its, it_start, it_end, dirT, downstream, nsub, np_skip


  ! check intersection of R with guiding surface
!  if (R%intersect_curve(C_guide, x, tau)) then
!     L    = R%length() * tau
!     write (6, 9000) L
!     ierr = 1
!     return
!  endif


  ! set up effective resolution for strike point area
  nsub  = np_sub_divertor
  np_SP = Z%nt * nsub / (nskip+1) + nextend


  ! check resolution
  np_skip = ip0;  if (Rside == UPPER) np_skip = this%np - ip0
  if (np_SP > this%np-np_skip) then
     write (6, 9010) this%np, np_skip, Z%nt, nsub, nextend
     stop
  endif


  ! setup poloidal boundary with reference nodes for flux surfaces
  !call this%setup_boundary_nodes(POLOIDAL, Rside, R, Sr)
  select case(Rside)
  case(LOWER)
     dir        = RIGHT_HANDED
     downstream = 1
     ips        = this%np - np_SP
  case(UPPER)
     dir        = LEFT_HANDED
     downstream = -1
     ips        = np_SP
     ! ip_ds    = 0
  end select


  ! find downstream direction
  dirT       = dir * Ip_sign * Bt_sign
  it_start   = np_SP
  it_end     = 0
  if (dirT > 0) then
     it_start = 0
     it_end   = np_SP
  endif


  ! initialize toroidal discretization for strike point adjustment
  call TSP%init(np_SP)
  ! extend poloidal discretization by nextend cells beyond target by extending
  ! the toroidal discretization in the appropriate direction
  iextend     = 0; if (dirT > 0) iextend = nextend
  if (iextend == 0) then
     ! add slice in ccw toroidal direction
     dphi = Z%phi(Z%nt) - Z%phi(Z%nt-1)
     do its=1,nextend
        TSP%phi(np_SP-nextend + its) = Z%phi(Z%nt) + its*dphi
     enddo
  else
     ! add slice in cw toroidal direction
     dphi = Z%phi(1) - Z%phi(0)
     do its=1,nextend
        TSP%phi(nextend - its) = Z%phi(0) - its*dphi
     enddo
  endif

  ! set up toroidal discretization for strike point adjustment
  TSP%it_base = Z%it_base * nsub / (nskip+1)  +  iextend
  do it=0,Z%nt
     TSP%phi(it*nsub + iextend) = Z%phi(it)
  enddo
  ! add sub-resolution
  do it=0,Z%nt-1
     dphi = Z%phi(it+1) - Z%phi(it)

     do its=1,nsub-1
        TSP%phi(it*nsub + its + iextend) = Z%phi(it) + 1.d0 * its/nsub * dphi
     enddo
  enddo


  M => this%mesh
  allocate (MSP(0:TSP%nt, 2))
  do ir=0,this%nr
     x = M(ir,ip0,:)

     ! generate flux surface from "upstream" location x to target
     call F%generate(x, dir, Trace_Step=0.1d0, AltSurf=C_guide)
     !call F%plot(filename='F.plt', append=.true.)

     ! generate nodes from which field lines connect to strike point x
     select case(dir)
     case(RIGHT_HANDED)
        x = F%x(F%n_seg, :)
     case(LEFT_HANDED)
        x = F%x(0,:)
     end select
     !write (99, *) x
     call align_strike_points(x, TSP, MSP)


     ! setup downstream strike point nodes
     do it=it_start,it_end,dirT
        ip = ips + abs(it - it_end)*dir
        M(ir,ip,:) = MSP(it,:)
     enddo


     ! length on flux surface for strike point mesh
     L0 = F%length()
     L  = 0.d0
     do it=TSP%it_base+dirT,it_end,dirT
        L = L + sqrt(sum(  (MSP(it-dirT,:)-MSP(it,:))**2))
        !write (92, *) MSP(it, :)
     enddo

     ! aligned strike point mesh extends beyond upstream reference location?
     if (L > L0) then
        write (6, 9020) ir
        open  (iu_err, file='error_strike_point_mesh1.plt')
        write (iu_err, *) M(ir,ip0,:)
        close (iu_err)
        open  (iu_err, file='error_strike_point_mesh2.plt')
        do it=0,TSP%nt
           write (iu_err, *) MSP(it,:)
        enddo
        close (iu_err)
        open  (iu_err, file='error_upstream_nodes.plt')
        do ir1=0,this%nr
           write (iu_err, *) M(ir1,ip0,:)
        enddo
        close (iu_err)
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

  ierr = 0
  deallocate (MSP)
 9000 format('ERROR: reference path for radial discretization crosses guiding surface!', &
             'intersection at L = ', f0.3)
 9010 format('ERROR: poloidal grid resolution too small!'//,&
             i0, ' - ', i0, ' < ', i0, ' x ', i0, ' + ', i0)
 9020 format('ERROR: aligned strike point mesh extends beyond upstream reference point ', &
             'at radial index ', i0//,&
             'see error_strike_point_mesh.plt')
  end subroutine make_divertor_grid
!=======================================================================



!=======================================================================
! copy mesh M onto this mesh at node index ir0, ip0
!=======================================================================
  subroutine copy(this, ir0, ip0, M)
  class(t_mfs_mesh)            :: this
  integer,          intent(in) :: ir0, ip0
  type(t_mfs_mesh), intent(in) :: M

  integer :: ir, ip


  ! check boundaries
  if (ir0 < 0  .or.  ir0 > this%nr  .or.  ip0 < 0  .or.  ip0 > this%np) then
     write (6, 9000);  write (6, 9001) ir0, ip0, this%nr, this%np;  stop
  endif

  ! check size
  if (ir0+M%nr > this%nr  .or.  ip0+M%np > this%np) then
     write (6, 9000);  write (6, 9002) ir0+M%nr, ip0+M%np, this%nr, this%np;  stop
  endif

  ! copy mesh
  do ir=0,M%nr
  do ip=0,M%np
     this%mesh(ir0+ir, ip0+ip, :) = M%mesh(ir, ip, :)
  enddo
  enddo

 9000 format('error in t_mfs_mesh%copy:')
 9001 format('initial node is outside of mesh!'//'ir0, ip0 = ',i0,', ',i0// &
             'nr,  np  = ',i0,', ',i0)
 9002 format('upper node is outside of mesh!'//'ir,  ip  = ',i0,', ',i0// &
             'nr,  np  = ',i0,', ',i0)
  end subroutine copy
!=======================================================================



!=======================================================================
  subroutine initialize_module_mfs_mesh(C_guide_)
  type(t_curve), intent(in) :: C_guide_


  C_guide = C_guide_

  end subroutine initialize_module_mfs_mesh
!=======================================================================

end module mfs_mesh
