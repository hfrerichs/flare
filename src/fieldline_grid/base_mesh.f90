module base_mesh
  use iso_fortran_env
  use grid
  implicit none
  private

  integer, parameter, public :: &
     RADIAL         = 1, &
     POLOIDAL       = 2, &
     LOWER          = 1, &
     LOWER_TO_UPPER = 1, &
     UPPER          = 2, &
     UPPER_TO_LOWER = 2


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
  integer      :: ir, ip, i11, i22


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


  select case(boundary_type)
  case(RADIAL)
     this%ir0 = ir
     i11      = 0
     i22      = this%np
     if (present(i1)) i11 = i1
     if (present(i2)) i22 = i2
     do ip=i11,i22
        tau = spacings%node(ip-i11, i22-i11)
        call C_boundary%sample_at(tau, x)
        M(ir, ip, :) = x
        write (6, *) ir, ip, x
     enddo
  case(POLOIDAL)
     this%ip0 = ip
     i11      = 0
     i22      = this%nr
     if (present(i1)) i11 = i1
     if (present(i2)) i22 = i2
     do ir=i11,i22
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
  subroutine make_orthogonal_grid(this, prange)
  use equilibrium, only: get_PsiN
  use xpaths
  class(t_base_mesh)  :: this
  integer, intent(in), optional :: prange(2)

  real(real64), dimension(:,:,:), pointer :: M
  type(t_xpath) :: R
  real(real64)  :: PsiN(0:this%nr), x(2), PsiN_final
  integer       :: ir, ir0, ir1, ir2, ip0, ip, ip1, ip2, direction


  M => this%mesh


  ! set poloidal range
  ip0 = this%ip0
  if (ip0 == 0) then
     ip1 = 1
     ip2 = this%np
  elseif (ip0 == this%np) then
     ip1 = 0
     ip2 = this%np - 1
  else
     write (6, 9000)
     write (6, 9001) ip0
     stop
  endif
  if (present(prange)) then
     ip1 = prange(1)
     ip2 = prange(2)
  endif
  write (6, *) 'poloidal range: ', ip0, ip1, ip2


  ! setup reference PsiN values
  do ir=0,this%nr
     x        = M(ir,ip0,:)
     PsiN(ir) = get_PsiN(x)
  enddo


  ! set radial range
  ir0 = this%ir0
  if (ir0 == 0) then
     ir1        = 1
     ir2        = this%nr
     direction  = ASCENT
     PsiN_final = PsiN(this%nr)
  elseif (ir0 == this%nr) then
     ir1        = 0
     ir2        = this%nr - 1
     direction  = DESCENT
     PsiN_final = PsiN(0)
  else
     write (6, 9000)
     write (6, 9002) ir0
     stop
  endif
  write (6, *) 'radial range: ', ir0, ir1, ir2


  ! set up nodes in poloidal range ip1->ip2
  do ip=ip1,ip2
     call R%generate(M(ir0,ip,:), direction, LIMIT_PSIN, PsiN_final, sampling=SAMPLE_PSIN)

     do ir=ir1,ir2
        call R%sample_at_PsiN(PsiN(ir), x, enforce_boundary=.true.)
        M(ir,ip,:) = x
     enddo
  enddo

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
  use divertor, only: C_guide
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

end module base_mesh
