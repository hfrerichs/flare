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

  real(real64), parameter :: &
     COMPRESSION_FACTOR = 0.2d0


  ! orthogonal, flux surface aligned mesh with strike point adjustment
  type, extends(t_grid), public :: t_mfs_mesh
     integer :: nr, np
     integer :: &
        ir0 = -1, & ! reference radial surface for mesh generation
        ip0 = -1    ! reference poloidal surface for mesh generation

     ! "upstream" poloidal neighbor
     class(t_mfs_mesh), pointer  :: upn => null()
     integer                     :: upn_side

     contains
     procedure :: initialize
     procedure :: connect_to
     procedure :: setup_boundary_nodes
     procedure :: plot_boundary
     procedure :: make_orthogonal_grid
     procedure :: make_interpolated_mesh
     procedure :: make_interpolated_submesh
     procedure :: make_divertor_grid
     procedure :: copy
     procedure :: arclength ! calculate arclength on flux surface ir between nodes ip1 and ip2
     procedure :: get_segment
     procedure :: compress  ! compress part of flux surface
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
  class(t_mfs_mesh), target   :: this
  class(t_mfs_mesh)   :: M
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
     M%upn => this
     M%upn_side = side

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
  subroutine plot_boundary(this, direction, side, prefix, iz)
  class(t_mfs_mesh)            :: this
  integer, intent(in)          :: direction, side, iz
  character(len=*), intent(in) :: prefix

  integer, parameter :: iu = 99

  real(real64), dimension(:,:,:), pointer :: M
  character(len=72) :: filename
  integer :: ir, ip


  M => this%mesh

  write (filename, 1000) prefix, iz, side
 1000 format(a,'_Z',i0,'_side',i0)
  open  (iu, file=filename)
  select case(direction)
  case(RADIAL)
     select case(side)
     case(LOWER)
        ir = 0
     case(UPPER)
        ir = this%nr
     end select
     do ip=0,this%np
        write (iu, *) M(ir, ip, :)
     enddo

  case(POLOIDAL)
     select case(side)
     case(LOWER)
        ip = 0
     case(UPPER)
        ip = this%np
     end select
     do ir=0,this%nr
        write (iu, *) M(ir, ip, :)
     enddo
  end select
  close (iu)

  end subroutine plot_boundary
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
!
! TODO: upstream = 0, downstream = np and flip orientation at the end?
!
! Rside: location of seed mesh
!=======================================================================
!  subroutine make_divertor_grid(this, R, Rside, Sr, P, Pside, Sp, Z, ierr)
  subroutine make_divertor_grid(this, Rside, ip0, Sp, Z, ierr)
  use run_control,    only: Debug
  use fieldline_grid, only: t_toroidal_discretization, np_sub_divertor, poloidal_discretization, ORTHOGONAL_AUTO_ADJUST
  use equilibrium, only: get_PsiN, Ip_sign, Bt_sign
  use curve2D
  use mesh_spacing
  use flux_surface_2D
  class(t_mfs_mesh)            :: this
!  type(t_curve),   intent(in)  :: R, P
!  integer,         intent(in)  :: Rside, Pside
!  type(t_spacing), intent(in)  :: Sr, Sp
  integer,         intent(in)  :: Rside, ip0
  type(t_spacing), intent(in)  :: Sp
  type(t_toroidal_discretization),    intent(in)  :: Z
  integer,         intent(out) :: ierr

  integer, parameter :: nskip = 0, nextend = 1, iu_err = 66

  type(t_flux_surface_2D) :: F
  type(t_toroidal_discretization) :: TSP

  real(real64), dimension(:,:,:), pointer :: M
  real(real64), dimension(:,:), allocatable :: MSP
  real(real64)  :: PsiN(0:this%nr), x(2), PsiN_final, tau, xi, L0, L, dphi
  real(real64)  :: Ladjust, Lu1, Lu2, Lcompress
  integer       :: ir, ir0, ir1, ir2, ip, ips, dir, iextend, np_SP
  integer       :: it, its, it_start, it_end, dirT, downstream, nsub, np_skip
  integer       :: ipu


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
     ipu        = 0
  case(UPPER)
     dir        = LEFT_HANDED
     downstream = -1
     ips        = np_SP
     ipu        = this%np
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
     !call F%generate(x, dir, Trace_Step=0.1d0, AltSurf=C_guide)
     call F%generate_branch(x, dir, ierr, reference_direction=CCW_DIRECTION, &
             Trace_Step=0.1d0, cutoff_boundary=C_guide, stop_at_boundary=.false.)
     ! A. successfull trace of flux surface to target
     if (ierr == 0) then
!        select case(dir)
!        case(RIGHT_HANDED)
           x = F%x(F%n_seg, :)
!        case(LEFT_HANDED)
!           x = F%x(0,:)
!        end select
        !write (99, *) x
        L0 = F%length()

     ! B. flux surface out of bounds, trace backwards to target and adjust upstream mesh
     elseif (ierr == -1) then
        write (6, *) 'x0 = ', x
        ! adjust "upstream" location
        ! 1. trace back to boundary
        call F%generate_branch(x, -dir, ierr, reference_direction=CCW_DIRECTION, &
                Trace_Step=0.1d0, cutoff_boundary=C_guide, stop_at_boundary=.false.)

        ! point on boundary
        x  = F%x(F%n_seg, :)
        L0 = -F%length()

     ! C. UNKOWN situation
     else
        write (6, *) 'error: t_flux_surface2D%generate_branch returned ierr = ', ierr
        stop
     endif
     if (Debug) call F%plot(filename='F.plt', append=.true.)


     ! generate nodes from which field lines connect to strike point x
     call align_strike_points(x, TSP, MSP)


     ! setup downstream strike point nodes
     do it=it_start,it_end,dirT
        ip = ips + abs(it - it_end)*dir
        M(ir,ip,:) = MSP(it,:)
     enddo


     ! length on flux surface for strike point mesh
     L  = 0.d0
     do it=TSP%it_base+dirT,it_end,dirT
        L = L + sqrt(sum(  (MSP(it-dirT,:)-MSP(it,:))**2))
        !write (92, *) MSP(it, :)
     enddo

     ! aligned strike point mesh extends beyond upstream reference location?
     ! -> adjust upstream location by L-L0
     Ladjust = L - L0
     if ((Ladjust > 0.d0)  .and.  poloidal_discretization /= ORTHOGONAL_AUTO_ADJUST) then
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
     if (Ladjust > 0.d0) then
        ! 0. how strong should we compress the mesh upstream?
        Lcompress = COMPRESSION_FACTOR * Ladjust
        if (Debug) write (6, 9021) ir, Ladjust, Lcompress
 9021 format('compressing upstream end of flux surface ',i0,' from ',f0.3,' into ',f0.3,' cm further upstream')

        ! 1. inside this mesh
        ! distance from upstream reference to upstream end
        Lu1 = this%arclength(ir,ip0,ipu)
        ! adjustment fits into this mesh
        Lu2 = Ladjust + Lcompress - Lu1
        if (Lu2 < 0.d0) then
           write (6, 9022)
        else
           write (6, 9023) Lu2

           if (Lu2 > Lcompress) then
              write (6, 9024)
              if (.not.associated(this%upn)) then
                 write (6, *) 'error: pointer to upstream mesh not associated!'
                 stop
              endif
              !call this%upn%compress(ir, ip0, ipu, Ladjust - Lu1, Lcompress)
              call adjust_upstream()
           else
              write (6, 9025)
           endif
        endif
 9022 format('adjustment can be done in divertor grid')
 9023 format('adjustment exceeds divertor grid by ',f0.3, ' cm')
 9024 format('complete adjustment located in upstream grid')
 9025 format('partial adjustment located in upstream grid')
     else
     ! interpolate nodes on flux surface between upstream orthogonal grid nodes
     ! and downstream strike point nodes
     call F%setup_length_sampling()
!     select case(dir)
!     case(RIGHT_HANDED)
        F%w = F%w / (1.d0 - L/L0)
!     case(LEFT_HANDED)
!        F%w = (F%w - L/L0) / (1.d0 - L/L0)
!     end select

     do ip=ip0+dir,ips-dir,dir
        tau = 1.d0 * (ip - ip0) / (ips - ip0)
        xi  = Sp%sample(tau)
!        if (dir == LEFT_HANDED) tau = 1.d0 - tau
        call F%sample_at(xi, x)
        M(ir,ip,:) = x
        !write (93, *) x, tau
     enddo
     endif
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
  contains
  !---------------------------------------------------------------------
  subroutine adjust_upstream()
  real(real64), dimension(:), allocatable :: w1, w2
  type(t_curve) :: C
  real(real64) :: dl, x(2), s
  integer :: n1, n2, itmp


  call this%upn%get_segment(ir, dir, Ladjust-Lu1, Lcompress, C, n2, w2)
  call C%plot(filename="segment.plt")
  write (6, *) 'n2  = ', n2
  write (6, *) 'Lu1 = ', Lu1

  write (6, *) 'ips = ', ips ! 41
  write (6, *) 'ipu = ', ipu ! 0
  write (6, *) 'dir = ', dir ! 1
  allocate (w1(abs(ipu-ips+dir)))

  itmp = 1
  write (6, *) 'ip0 = ', ip0
  n1 = abs(ipu-ips+dir)
  write (6, *) 'n1  = ', n1

  ! adjust upstream grid

  ! adjust divertor grid
  call C%setup_length_sampling()
  do ip=ips-dir,ipu,-dir
!     dl = sqrt(sum( (M(ir,ip,:)-M(ir,ip+dir,:))**2 ))
!     write (6, *) itmp, dl
!     w1(itmp) = dl;  itmp = itmp + 1

  !   do ip=ip0+dir,ips-dir,dir
  !      tau = 1.d0 * (ip - ip0) / (ips - ip0)
  !      if (dir == LEFT_HANDED) tau = 1.d0 - tau
     s = 1.d0 * (ip-ips) / (ipu-ips) * n1 / (n1+n2)
     call C%sample_at(s, x)
     M(ir,ip,:) = x
  enddo
!  write (6, *) 'sum(w1) = ', sum(w1)
!  write (6, *) 'sum(w2) = ', sum(w2)

  end subroutine adjust_upstream
  !---------------------------------------------------------------------
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
! calculate arclength on grid/flux surface ir between nodes ip1 and ip2
!=======================================================================
  function arclength(this, ir, ip1, ip2) result(l)
  class(t_mfs_mesh)   :: this
  integer, intent(in) :: ir, ip1, ip2
  real(real64)        :: l

  real(real64), dimension(:,:,:), pointer :: M
  integer :: ip


  M => this%mesh

  ! check inputs
  call irange_check(ir,  0, this%nr, 'ir', 't_mfs_mesh%arclength')
  call irange_check(ip1, 0, this%np, 'ip', 't_mfs_mesh%arclength')
  call irange_check(ip2, 0, this%np, 'ip', 't_mfs_mesh%arclength')


  l = 0.d0
  do ip=min(ip1,ip2)+1,max(ip1,ip2)
     l = l + sqrt(sum( (M(ir,ip,:)-M(ir,ip-1,:))**2 ))
  enddo

  end function arclength
!=======================================================================



!=======================================================================
! return segment of flux surface ir
! length    arclength of segment to be returned
! offset    offset from end of flux surface where segment begins
! side      select from which end of flux surface to start
!=======================================================================
  subroutine get_segment(this, ir, side, offset, length, C, nseg, wseg)
  class(t_mfs_mesh)        :: this
  integer, intent(in)      :: ir, side
  real(real64), intent(in) :: offset, length
  type(t_curve), intent(out) :: C
  integer,       intent(out) :: nseg
  real(real64), dimension(:), allocatable :: wseg


  real(real64), dimension(:,:,:), pointer :: M
  real(real64), dimension(:,:), allocatable :: xtmp
  real(real64), dimension(:), allocatable :: wtmp
  real(real64) :: L, dl, w
  integer :: ip, ip1, ip2, ipdir, itmp, ip_next

  ! check inputs
  call irange_check(ir, 0, this%nr, 'ir', 't_mfs_mesh%get_segment')

  select case(side)
  case(1)
     ip1 = this%np
     ip2 = 0
  case(-1)
     ip1 = 0
     ip2 = this%np
  case default
     write (6, *) 'error: invalid argument side = ', side, '!'
     stop
  end select
  ipdir = -side


  write (6, *) 'offset = ', offset
  write (6, *) 'length = ', length
  write (6, *) 'ipdir  = ', ipdir

  allocate (xtmp(0:this%np,2), wtmp(0:this%np))
  L    = 0.d0
  itmp = 0
  M    => this%mesh
  do ip=ip1+ipdir,ip2,ipdir
     dl = sqrt(sum( (M(ir,ip,:)-M(ir,ip-ipdir,:))**2 ))
     wtmp(ip) = dl
     L  = L + dl

     write (6, *) ip, L
     if (L < offset) cycle

     ! first point
     if (itmp == 0) then
        w = (L - offset) / dl
        xtmp(itmp,:) = w * M(ir,ip-ipdir,:) + (1.d0-w) * M(ir,ip,:)
        write (6, *) 'first point ', ip, xtmp(itmp,:)
        itmp = itmp + 1
     endif

     ! last point
     if (L > offset + length) then
        w = (L - offset - length) / dl
        xtmp(itmp,:) = w * M(ir,ip-ipdir,:) + (1.d0-w) * M(ir,ip,:)
        wtmp(ip)   = wtmp(ip) * (1.d0-w)
        ip_next = ip
        write (6, *) 'last  point ', ip, xtmp(itmp,:)
        exit
     endif

     ! this node
     xtmp(itmp,:) = M(ir,ip,:)
     write (6, *) 'point ', ip, xtmp(itmp,:)
     itmp         = itmp + 1
  enddo
  call make_2D_curve(itmp+1, xtmp(0:itmp,1), xtmp(0:itmp,2), C)


  select case(side)
  case(1)
     nseg = this%np - ip_next
  case(-1)
     nseg = ip_next
  end select

  allocate (wseg(nseg))
  itmp = 1
  do ip=ip1+ipdir,ip_next,ipdir
     wseg(itmp) = wtmp(ip)
     itmp = itmp + 1
  enddo

  ! cleanup
  deallocate (xtmp, wtmp)

  end subroutine get_segment
!=======================================================================



!=======================================================================
! compress flux surface ir between ip1 and ip2
!=======================================================================
  subroutine compress(this, ir, ip1, ip2, offset, length)
  class(t_mfs_mesh)   :: this
  integer, intent(in) :: ir, ip1, ip2
  real(real64)        :: offset, length

  type(t_curve) :: C

  ! check inputs
  call irange_check(ir,  0, this%nr, 'ir', 't_mfs_mesh%compress')
  call irange_check(ip1, 0, this%np, 'ip', 't_mfs_mesh%compress')
  call irange_check(ip2, 0, this%np, 'ip', 't_mfs_mesh%compress')


  !call C = this%get_segment(ir, 1, offset, length)

  end subroutine compress
!=======================================================================



!=======================================================================
  subroutine initialize_module_mfs_mesh(C_guide_)
  type(t_curve), intent(in) :: C_guide_


  C_guide = C_guide_

  end subroutine initialize_module_mfs_mesh
!=======================================================================

end module mfs_mesh
