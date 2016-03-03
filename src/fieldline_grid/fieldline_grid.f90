module fieldline_grid
  use iso_fortran_env
  use system
  use math
  use mesh_spacing
  use xpaths
  implicit none

!.......................................................................
  ! topology definitions
  character(len=*), parameter :: &
     TOPO_SC     = 'simply_connected', &
     TOPO_SC1    = 'sc', &
     TOPO_LSN    = 'lower_single_null', &
     TOPO_LSN1   = 'lsn', &
     TOPO_DDN    = 'disconnected_double_null', &
     TOPO_DDN1   = 'ddn', &
     TOPO_CDN    = 'connected_double_null', &
     TOPO_CDN1   = 'cdn'


  ! Discretization type definitions
  character(len=*), parameter :: &
     POLOIDAL_ANGLE    = 'poloidal_angle_fixed', &
     ORTHOGONAL        = 'quasi_orthogonal'


  ! zone type definitions
  integer, parameter :: &
     TYPE_SINGLE_LAYER = -1, &
     TYPE_HPR          =  0, &
     TYPE_SOL          =  1, &
     TYPE_SOLMAP       = 11, &
     TYPE_PFR          =  2


  ! surface type definitions
  integer, parameter :: &
     SF_PERIODIC =  1, &
     SF_UPDOWN   =  2, &
     SF_MAPPING  =  3, &
     SF_CORE     = -1, &
     SF_VACUUM   = -2


  ! Type of innermost flux surface (exact or quasi flux surface)
  character(len=*), parameter :: &
     SF_EXACT    = 'EXACT', &
     SF_QUASI    = 'QUASI'


  integer, parameter :: &
     max_blocks  = 360, &        ! Maximum number of toroidal blocks
     max_zones   = 360, &        ! Maximum number of zones
     max_layers  = 6             ! Maximum number of layer (zones per block)
!.......................................................................


!.......................................................................
! user defined variables
!.......................................................................
  character(len=80) :: &
     topology                         = TOPO_SC, &
     Innermost_Flux_Surface           = SF_EXACT, &
     discretization_method            = POLOIDAL_ANGLE, &
     radial_spacing(0:max_layers-1)   = '', &
     poloidal_spacing(0:max_layers-1) = '', &
     toroidal_spacing(0:max_layers-1) = '', &
     guiding_surface                  = '', &
     N0_file(0:max_layers-1)          = '', &
     N0_method(0:max_layers-1)        = ''


  integer :: &
     symmetry            =   1, &          ! toroidal symmetry (i.e. 5 => 72 deg)
     blocks              =   1, &          ! number of toroidal blocks
     nt                  =  12, &          ! default toroidal resolution
     np(0:max_layers-1)  = 360, &          ! default poloidal resolution
     npL(0:max_layers-1) =  30, &          ! default poloidal resolution in divertor legs
     npR(0:max_layers-1) =  30, &          !    (L)eft and (R)ight segments
     nr(0:max_layers-1)  =  32, &          ! default radial resolution
     n_interpolate       =   4, &          ! number of interpolated flux surfaces (for the transition between the pair of perturbed flux surfaces at the inner simulation boundary and unperturbed flux surfaces further outside)
     nr_EIRENE_core      =   1, &          ! radial resolution in core (EIRENE only)
     nr_EIRENE_vac       =   1, &          ! radial resolution in vacuum (EIRENE only)
     nr_perturbed        =   2             ! number of perturbed flux surfaces at the inner boundary
! TODO: aligned surfaces for neutrals (nr_EIRENE_core_aligned, ...)

  real(real64) :: &
     phi0                = -360.d0, &      ! lower boundary of simulation domain [deg]
     x_in1(3)            = (/120.d0, 0.d0, 0.d0/), &  ! reference points (R[cm], Z[cm], phi[deg]) ...
     x_in2(3)            = (/119.d0, 0.d0, 0.d0/), &  ! ... on 1st and 2nd innermost flux surfaces
     d_SOL(2)            =   24.d0, &      ! radial width of scrape-off layer
     d_PFR(2)            =   15.d0, &      ! radial width of private flux region
     d_N0(0:max_layers-1)=   10.d0, &      ! radial width of vacuum region
     d_extend(0:max_layers-1,-1:1) = -1.d0, &      ! poloidal extension of divertor leg (used in close_grid_domain)
     d_cutL(2)           =    6.d0, &      ! cut-off length for flux surfaces behind the wall
     d_cutR(2)           =    8.d0, &      !    (L)eft and (R)ight segments
     alphaL(2)           =    0.9d0, &     ! Relative length of divertor legs at outermost boundary
     alphaR(2)           =    1.0d0, &     !    (L)eft and (R)ight segments
     etaL(2)             =    0.8d0, &     ! fraction of cells in front of the target
     etaR(2)             =    0.8d0, &     !    (L)eft and (R)ight segments
     Dtheta_sampling     =    20.d0, &     ! Transition between angular and length weighted sampling of flux surfaces
     Dtheta_separatrix   =     0.d0        ! ... same on separatrix

  logical :: &
     extend_alpha_SOL2   =  .true.



  ! user defined input for individual blocks
  type t_block_input
     integer :: &
        nr(0:max_layers-1)  = -1, & ! radial resolution
        np(0:max_layers-1)  = -1, & ! poloidal resolution
        npL(0:max_layers-1) = -1, & ! poloidal resolution in divertor legs
        npR(0:max_layers-1) = -1, &
        nt      = -1, &             ! number of cells in toroidal direction
        it_base = -1                ! 0 <= index of base grid position <= nt

     real(real64) :: &
        width   = -360.d0          ! toroidal width of block [deg]
  end type t_block_input
!.......................................................................


!.......................................................................
! internal variables
!.......................................................................
  character(len=*), parameter :: config_file = 'grid.conf'


  ! extended block data
  type, extends(t_block_input) :: t_block
     ! toroidal cell spacings
     type(t_spacing) :: S = Equidistant

     real(real64) :: &
        phi_base    = -360.d0, &          ! position of base grid [deg]
        phi_left    = -360.d0, &
        phi_right   = -360.d0

     real(real64), dimension(:), allocatable :: phi
  end type t_block
  type(t_block) :: Block(0:max_blocks-1)


  ! zone data
  type :: t_zone
     ! resolution in radial, poloidal and toroidal direction
     integer :: nr, np, nt

     ! connectivity between zones (surface types = periodic, mapping, ...)
     integer :: isfr(2), isfp(2), isft(2)

     ! zone type information
     integer :: itypeR, itypeP

     ! surface indices for plasma transport range
     integer :: r_surf_pl_trans_range(2), p_surf_pl_trans_range(2)

     ! toroidal discretization
     integer :: it_base ! index of base grid position
     real(real64), dimension(:), allocatable :: phi

     ! spacing for radial and poloidal discretization
     type(t_spacing) :: Sr = Equidistant, Sp = Equidistant

     ! additional domain for neutral particles
     real(real64) :: d_N0 = 0.d0
     real(real64) :: d_extend(-1:1) = 0.d0
     character(len=80) :: N0_file, N0_method

     contains
     procedure :: setup
  end type t_zone
  type(t_zone) :: Zone(0:max_zones-1)


  character(len=32) :: label(0:max_layers-1) = ''
  type(t_xpath)     :: rpath(0:max_layers-1)
  logical      :: default_decomposition
  real(real64) :: Delta_phi_sim
  integer      :: layers
!.......................................................................

  contains
!=======================================================================



!=======================================================================
! Set up zone for block "iblock" and layer "ilayer".
! Zone layout itypeR = TYPE_SINGLE_LAYER, TYPE_HPR, TYPE_SOL, TYPE_PFR
!             itypeP = SF_PERIODIC, SF_VACUUM
!=======================================================================
  subroutine setup(this, iblock, ilayer, itypeR, itypeP)
  class(t_zone)       :: this
  integer, intent(in) :: iblock, ilayer, itypeR, itypeP

  integer :: add1 = 0, add2 = 0


  ! 1. grid resolution
  call setup_resolution()


  ! 2. default toroidal discretization
  allocate (this%phi(0:this%nt))
  this%phi     = Block(iblock)%phi
  this%it_base = Block(iblock)%it_base


  ! 3.1 plamsa transport domain
  this%r_surf_pl_trans_range(1) = add1
  this%r_surf_pl_trans_range(2) = this%nr - add2
  this%p_surf_pl_trans_range(1) = 0
  this%p_surf_pl_trans_range(2) = this%np

  ! 3.2 set parameters for additional neutral domain
  this%d_N0      = d_N0(ilayer)
  this%d_extend  = d_extend(ilayer,-1:1)
  this%N0_file   = N0_file(ilayer)
  this%N0_method = N0_method(ilayer)


  ! 4. boundaries and zone type
  this%itypeR = itypeR
  this%itypeP = itypeP
  call setup_boundaries()

  contains
  !---------------------------------------------------------------------
  subroutine setup_resolution()


  ! toroidal resolution
  this%nt      = Block(iblock)%nt


  ! radial and poloidal resolution (depends on type of zone)
  add1 = 0; add2 = 0
  select case(itypeR)
  ! single layer: add core domain on lower and vaccum domain on upper radial boundary
  case(TYPE_SINGLE_LAYER)
     add1 = nr_EIRENE_core
     add2 = nr_EIRENE_vac

  ! high pressure region (HPR): add core domain on lower radial boundary
  case(TYPE_HPR)
     add1 = nr_EIRENE_core

  ! scrape-off layer (SOL): add vacuum domain on upper radial boundary
  case(TYPE_SOL)
     add2 = nr_EIRENE_vac

  ! private flux region (PFR): add vacuum domain on lower radial boundary
  case(TYPE_PFR)
     add1 = nr_EIRENE_vac
  end select
  this%nr = Block(iblock)%nr(ilayer) + add1 + add2
  this%np = Block(iblock)%np(ilayer)

  end subroutine setup_resolution
  !---------------------------------------------------------------------
  subroutine setup_boundaries()


  ! 1. radial boundaries
  select case(itypeR)
  case(TYPE_SINGLE_LAYER)
     this%isfr(1) = SF_CORE
     this%isfr(2) = SF_VACUUM

  case(TYPE_HPR)
     this%isfr(1) = SF_CORE
     this%isfr(2) = SF_MAPPING

  case(TYPE_SOL)
     this%isfr(1) = SF_MAPPING
     this%isfr(2) = SF_VACUUM

  case(TYPE_SOLMAP)
     this%isfr(1) = SF_MAPPING
     this%isfr(2) = SF_MAPPING

  case(TYPE_PFR)
     this%isfr(1) = SF_VACUUM
     this%isfr(2) = SF_MAPPING
  end select


  ! 2. poloidal boundaries
  this%isfp = itypeP


  ! 3. toroidal boundaries
  this%isft(1) = SF_MAPPING
  this%isft(2) = SF_MAPPING

  end subroutine setup_boundaries
  !---------------------------------------------------------------------
  end subroutine setup
!=======================================================================



!=======================================================================
  subroutine setup_grid_configuration

  integer, parameter :: iu = 12

  type(t_block_input) :: Block(0:max_blocks-1)

  namelist /FieldlineGrid_Input/ &
     topology, symmetry, blocks, Block, &
     phi0, x_in1, x_in2, d_SOL, d_PFR, d_N0, N0_file, N0_method, d_extend, &
     nt, np, npL, npR, nr, nr_EIRENE_core, nr_EIRENE_vac, &
     n_interpolate, nr_perturbed, &
     radial_spacing, poloidal_spacing, toroidal_spacing, &
     d_cutL, d_cutR, etaL, etaR, alphaL, alphaR, extend_alpha_SOL2, &
     Dtheta_sampling, Dtheta_separatrix, &
     discretization_method, guiding_surface


  ! 1. read user configuration from input file
  open  (iu, file='run_input', err=9000)
  read  (iu, FieldlineGrid_Input, err=9000, end=9100)
  close (iu)
  if (blocks > max_blocks) then
     write (6, *) 'error: number of blocks exceeds maximum'
     write (6, *) blocks, ' > ', max_blocks
     stop
  endif
  if (nr_perturbed < 1  .or. nr_perturbed > 2) then
     write (6, *) 'error: nr_perturbed = 1 or 2 required!'
     stop
  endif
  Dtheta_sampling   = Dtheta_sampling / 180.d0 * pi
  Dtheta_separatrix = Dtheta_separatrix / 180.d0 * pi


  ! 2. setup size and position of toroidal blocks
  call setup_toroidal_blocks(Block)


  ! 3. setup mesh spacing functions
  call setup_mesh_spacing()


  return
 1000 format(3x,'- Topology of configuration: ',a)
 9000 write (6, *) 'error while reading input file ', trim(config_file), '!'
  stop
 9100 write (6, *) 'error: cannot find FieldlineGrid_Input namelist in run control file!'
  stop
  end subroutine setup_grid_configuration
!=======================================================================



!=======================================================================
  subroutine setup_toroidal_blocks(Block_input)
  type(t_block_input), intent(in) :: Block_input(0:max_blocks-1)

  real(real64) :: tmp, xi, Dphi
  integer      :: ib, it, i


  ! 0. Initialize variable Block
  default_decomposition = .true.
  Block%t_block_input = Block_input
  do ib=0,blocks-1
     ! set toroidal resolution
     call set_value(Block(ib)%nt, nt)

     ! set index of base plane
     if (Block(ib)%it_base == -1) then
        Block(ib)%it_base = Block(ib)%nt / 2
     else
        default_decomposition = .false.
     endif
     ! check input
     if (Block(ib)%it_base < 0  .or.  Block(ib)%it_base > Block(ib)%nt) then
        write (6, 9002) ib, Block(ib)%it_base, Block(ib)%nt
        stop
     endif
  enddo
 9001 format('error: invalid toroidal resolution in block ',i0,': ',i0)
 9002 format('error: invalid base index i0 = ',i0, '(0 <= i0 <= ',i0,' required!)')



  ! 1. set total size of simulation domain
  Delta_phi_sim         = real(360, real64) / symmetry


  ! 2. set size of toroidal blocks
  ! DEFAULT: Delta_phi_sim / (number of blocks), split equally in forward and backward direction
  tmp = 0.d0
  do ib=0,blocks-1
     if (Block(ib)%width < 0.d0) then
        Block(ib)%width = Delta_phi_sim / blocks
     else
        default_decomposition = .false.
     endif
     tmp             = tmp + Block(ib)%width
  enddo
  if (abs(Delta_phi_sim - tmp) > epsilon_r64) then
     write (6, 9003)
     stop
  endif
 9003 format('error: block sizes do not add up to size of simulation domain!')


  ! 3. set lower boundary of simulation domain
  ! DEFAULT: neg. half of first block
  if (phi0 == -360.d0) then
     phi0 = -Block(0)%width / 2.d0
  else
     default_decomposition = .false.
  endif


  ! 4a. set absolute position of blocks
  Block(0)%phi_left  = phi0
  Block(0)%phi_right = phi0 + Block(0)%width
  do ib=1,blocks-1
     Block(ib)%phi_left  = Block(ib-1)%phi_right
     Block(ib)%phi_right = Block(ib)%phi_left + Block(ib)%width
  enddo
  ! 4b. set position of base planes
  ! TODO: non-default toroidal spacing
  do ib=0,blocks-1
     xi   = Block(ib)%S%node(Block(ib)%it_base, Block(ib)%nt)
     Dphi = Block(ib)%phi_right - Block(ib)%phi_left
     Block(ib)%phi_base  = Block(ib)%phi_left + xi * Dphi
  enddo


  ! 5. setup toroidal discretization
  do ib=0,blocks-1
     nt = Block(ib)%nt
     allocate (Block(ib)%phi(0:nt))

     do it=0,Block(ib)%nt
        xi   = Block(ib)%S%node(it, Block(ib)%nt)
        Dphi = Block(ib)%phi_right - Block(ib)%phi_left
        Block(ib)%phi(it) = Block(ib)%phi_left + xi * Dphi
     enddo
  enddo


  ! 6. radial and poloidal resolution
  do ib=0,blocks-1
     do i=0,max_layers-1; call set_value(Block(ib)%nr(i),  nr(i)); enddo
     do i=0,max_layers-1; call set_value(Block(ib)%np(i),  np(i)); enddo
     do i=0,max_layers-1; call set_value(Block(ib)%npL(i), npL(i)); enddo
     do i=0,max_layers-1; call set_value(Block(ib)%npR(i), npR(i)); enddo
  enddo


  ! 7. output to screen
  write (6, *)
  write (6, 1000) phi0, phi0 + Delta_phi_sim
  write (6, 1001)
  if (default_decomposition) write (6, 1004)
  do ib=0,blocks-1
     write (6, 1002) ib, Block(ib)%phi_base, Block(ib)%phi_left, Block(ib)%phi_right
  enddo
 1000 format (3x,'- Decomposition of simulation domain (',f7.3,' -> ',f7.3,' deg):')
 1001 format (8x,'block #, base location [deg], domain [deg]')
 1002 format (8x,      i7,5x,f7.3,':',5x,f7.3,' -> ',f7.3)
 1004 format (8x,'using default decomposition')

  contains
!.......................................................................
  subroutine set_value(n, n_default)
  integer, intent(inout) :: n
  integer, intent(in)    :: n_default

  if (n == -1) n = n_default
  if (n <= 0) then
     write (6, 9000) n
     stop
  endif
 9000 format('error: non-positive value ',i0,' not allowed!')

  end subroutine set_value
!.......................................................................
  end subroutine setup_toroidal_blocks
!=======================================================================



!=======================================================================
  subroutine load_local_resolution(iblock)
  integer, intent(in) :: iblock

  integer :: ilayer


  do ilayer=0,max_layers-1
     nr(ilayer)  = Block(iblock)%nr(ilayer)
     np(ilayer)  = Block(iblock)%np(ilayer)
     npL(ilayer) = Block(iblock)%npL(ilayer)
     npR(ilayer) = Block(iblock)%npR(ilayer)
  enddo

  end subroutine load_local_resolution
!=======================================================================



!=======================================================================
  subroutine setup_mesh_spacing()
  end subroutine setup_mesh_spacing
!=======================================================================


!=======================================================================
  subroutine initialize_emc3_grid()
  use emc3_grid
  use grid
  use string
  use system

  type(t_grid) :: G
  real(real64) :: phi, delta
  integer :: iz, it, itz, ip, ip0, ir, ir0, ig, nr1, nr2, np1, np2


  ! 1a. allocate grid resolution arrays
  allocate (SRF_RADI(0:NZONET-1),SRF_POLO(0:NZONET-1),SRF_TORO(0:NZONET-1), &
            ZON_RADI(0:NZONET-1),ZON_POLO(0:NZONET-1),ZON_TORO(0:NZONET-1) )
  allocate (R_SURF_PL_TRANS_RANGE(2,0:NZONET-1), P_SURF_PL_TRANS_RANGE(2,0:NZONET-1))
  ! 1b. set grid resolution and plasma transport range
  do iz=0,NZONET-1
     ZON_RADI(iz) = Zone(iz)%nr
     ZON_POLO(iz) = Zone(iz)%np
     ZON_TORO(iz) = Zone(iz)%nt

     R_SURF_PL_TRANS_RANGE(1:2,iz) = Zone(iz)%r_surf_pl_trans_range(1:2)
     P_SURF_PL_TRANS_RANGE(1:2,iz) = Zone(iz)%p_surf_pl_trans_range(1:2)
  enddo
  SRF_RADI = ZON_RADI + 1
  SRF_POLO = ZON_POLO + 1
  SRF_TORO = ZON_TORO + 1


  ! 2. setup zone offset arrays
  allocate (PHI_PL_OS(0:NZONET), GRID_P_OS(0:NZONET), MESH_P_OS(0:NZONET))
  PHI_PL_OS = 0
  GRID_P_OS = 0
  MESH_P_OS = 0
  do iz=1,NZONET
     PHI_PL_OS(iz) = PHI_PL_OS(iz-1) + SRF_TORO(iz-1)
     GRID_P_OS(iz) = GRID_P_OS(iz-1) + SRF_RADI(iz-1)*SRF_POLO(iz-1)*SRF_TORO(iz-1)
     MESH_P_OS(iz) = MESH_P_OS(iz-1) + ZON_RADI(iz-1)*ZON_POLO(iz-1)*ZON_TORO(iz-1)
  enddo


  ! 3. allocate main arrays
  allocate (PHI_PLANE(0:PHI_PL_OS(NZONET)-1), &
                   RG(0:GRID_P_OS(NZONET)-1), &
                   ZG(0:GRID_P_OS(NZONET)-1))


  ! 4. load field lines
  do iz=0,NZONET-1
     call G%load('fieldlines_'//trim(str(iz))//'.dat', silent=.true.)


     ! 4.1 check input
     ! 4.1.1 radial resolution
     nr1 = R_SURF_PL_TRANS_RANGE(1,iz)
     nr2 = R_SURF_PL_TRANS_RANGE(2,iz)
     if (G%n1-1 /= nr2-nr1) then
        write (6, *) 'error: mismatching radial resolution: ', G%n1
        write (6, *) 'expected index range for aligned grid: ', nr1, '->', nr2
        stop
     endif
     ! 4.1.2 poloidal resolution
     np1 = P_SURF_PL_TRANS_RANGE(1,iz)
     np2 = P_SURF_PL_TRANS_RANGE(2,iz)
     if (G%n2-1 /= np2-np1) then
        write (6, *) 'error: mismatching poloidal resolution: ', G%n2
        write (6, *) 'expected index range for aligned grid: ', np1, '->', np2
        stop
     endif
     ! 4.1.3 toroidal resolution
     nt = Zone(iz)%nt
     if (G%n3-1 /= nt) then
        write (6, *) 'error: mismatching toroidal resolution: ', G%n3
        write (6, *) 'expected resolution is: ', nt
        stop
     endif


     ! 4.2 set position of slices
     delta = machine_precision*1000.d0
     do it=0,nt
        phi = Zone(iz)%phi(it)
        if (abs(phi-G%x3(it)) > delta) then
           write (6, *) 'error: mismatching toroidal positions: ', phi, G%x3(it)
           write (6, *) 'at slice ', it, ' in zone ', iz
           write (6, *) 'delta phi = ', abs(phi-G%x3(it)), ' > ', delta
           stop
        endif

        itz = it + PHI_PL_OS(iz)
        PHI_PLANE(itz) = phi / 180.d0 * pi
     enddo


     ! 4.3 setup grid from field lines
     ir0 = R_SURF_PL_TRANS_RANGE(1,iz)
     ip0 = P_SURF_PL_TRANS_RANGE(1,iz)
     do it=0,nt
     do ip=0,G%n2-1
     do ir=0,G%n1-1
        ig = ir0 + ir + (ip0 + ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        RG(ig) = G%mesh3D(ir,ip,it,1)
        ZG(ig) = G%mesh3D(ir,ip,it,2)
     enddo
     enddo
     enddo
  enddo

  end subroutine initialize_emc3_grid
!=======================================================================



!=======================================================================
! WRITE INPUT FILES (input.geo, input.n0g, input.par)
!=======================================================================
  subroutine write_emc3_input_files
  use emc3_grid
  implicit none

  integer, parameter :: iu = 72


  call write_input_geo()
  call write_input_n0g()
  call write_input_par()

  contains
  !---------------------------------------------------------------------
  subroutine write_input_geo
  integer :: ir, ip, it, iz, irun, n


  open  (iu, file='input.geo')
  write (iu, 1000)
 1000 format ('* geometry information for EMC3')

  ! 1. geometry, mesh resolution
  write (iu, 9999)
  write (iu, 1001)
  write (iu, 9999)
  write (iu, 1002)
  write (iu, 1003) NZONET
  write (iu, 1004)
  do iz=0,NZONET-1
     write (iu, *) SRF_RADI(iz), SRF_POLO(iz), SRF_TORO(iz)
  enddo
 1001 format ('*** 1. grid resolution')
 1002 format ('* number of zones/blocks')
 1003 format (i0)
 1004 format ('* number of radial, poloidal and toroidal grid points')


  ! 2. surface definitions
  write (iu, 9999)
  write (iu, 2000)
  write (iu, 9999)
  ! 2.1 non default surfaces (periodic, mapping, ...)
  write (iu, 2001)
  ! 2.1.a - radial
  write (iu, 2002)
  n = 0
  do irun=0,1
     ! write number of non default radial surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,NZONET-1
        if (Zone(iz)%isfr(1) > 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) 0, iz, Zone(iz)%isfr(1)
              write (iu, *) 0, ZON_POLO(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
        if (Zone(iz)%isfr(2) > 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) SRF_RADI(iz)-1, iz, Zone(iz)%isfr(2)
              write (iu, *) 0, ZON_POLO(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
     enddo
  enddo
  ! 2.1.b - poloidal
  write (iu, 2003)
  n = 0
  do irun=0,1
     ! write number of non default poloidal surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,NZONET-1
        if (Zone(iz)%isfp(1) > 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) 0, iz, Zone(iz)%isfp(1)
              write (iu, *) 0, ZON_RADI(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
        if (Zone(iz)%isfp(2) > 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) SRF_POLO(iz)-1, iz, Zone(iz)%isfp(2)
              write (iu, *) 0, ZON_RADI(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
     enddo
  enddo
  ! 2.1.c - toroidal
  write (iu, 2004)
  n = 0
  do irun=0,1
     ! write number of non default toroidal surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,NZONET-1
        if (Zone(iz)%isft(1) > 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) 0, iz, Zone(iz)%isft(1)
              write (iu, *) 0, ZON_RADI(iz)-1, 0, ZON_POLO(iz)-1
           endif
        endif
        if (Zone(iz)%isft(2) > 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) SRF_TORO(iz)-1, iz, Zone(iz)%isft(2)
              write (iu, *) 0, ZON_RADI(iz)-1, 0, ZON_POLO(iz)-1
           endif
        endif
     enddo
  enddo

  ! 2.2 non transparent surfaces (boundary conditions)
  write (iu, 2005)
  ! 2.2.a - radial
  write (iu, 2002)
  n = 0
  do irun=0,1
     ! write number of non transparent radial surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,NZONET-1
        if (Zone(iz)%isfr(1) < 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) R_SURF_PL_TRANS_RANGE(1,iz), iz, 1
              write (iu, *) 0, ZON_POLO(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
        if (Zone(iz)%isfr(2) < 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) R_SURF_PL_TRANS_RANGE(2,iz), iz, -1
              write (iu, *) 0, ZON_POLO(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
     enddo
  enddo
  ! 2.2.b - poloidal
  write (iu, 2003)
  n = 0
  do irun=0,1
     ! write number of non transparent poloidal surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,NZONET-1
        if (Zone(iz)%isfp(1) < 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) P_SURF_PL_TRANS_RANGE(1,iz), iz, 1
              write (iu, *) 0, ZON_RADI(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
        if (Zone(iz)%isfp(2) < 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) P_SURF_PL_TRANS_RANGE(2,iz), iz, -1
              write (iu, *) 0, ZON_RADI(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
     enddo
  enddo
  ! 2.2.c - toroidal
  write (iu, 2004)
  write (iu, *) 0

  ! 2.3 plate surfaces
  write (iu, 2006)
  write (iu, 2002)
  write (iu, *) -3 ! user defined
  write (iu, 2003)
  write (iu, *) -3 ! user defined
  write (iu, 2004)
  write (iu, *) -3 ! user defined
 2000 format ('*** 2. surface definitions')
 2001 format ('*** 2.1 non default surface')
 2002 format ('* radial')
 2003 format ('* poloidal')
 2004 format ('* toroidal')
 2005 format ('*** 2.2 non transparent surface (Boundary condition must be defined)')
 2006 format ('*** 2.3 plate surface (Bohm Boundary condition)')


  ! 3. physical cell definition
  write (iu, 9999)
  write (iu, 3000)
  write (iu, 9999)
  write (iu, *) -1
  write (iu, 3001)
  write (iu, 3002) .true.
 3000 format ('*** 3. physical cell definition')
 3001 format ('* run cell check?')
 3002 format (L1)
  close (iu)
 9999 format ('*',32('-'))

  end subroutine write_input_geo
  !---------------------------------------------------------------------
  subroutine write_input_n0g

  integer :: ir, iz, irun, n


  open  (iu, file='input.N0G')
  write (iu, 1000)
  write (iu, 9999)
  write (iu, 1001)
  write (iu, 9999)
  write (iu, 1002)
  write (iu, 1003)
  write (iu, 1004)
 1000 format ('******** additional geometry and parameters for EIRENE ****')
 1001 format ('*** 1. non-transparent surfaces for neutral particles')
 1002 format ('*  non-transparent surfaces with informations about')
 1003 format ('*  this surface being defined in EIRENE. The surface')
 1004 format ('*  number must be indicated here.')

  ! 1. non-transparent radial surfaces
  write (iu, 1012)
  n = 0
  do irun=0,1
     ! write number of non default radial surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,NZONET-1
        if (Zone(iz)%isfr(1) < 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, 1015) 0, iz, -Zone(iz)%isfr(1)
              write (iu, *) 0, ZON_POLO(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
        if (Zone(iz)%isfr(2) < 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, 1015) SRF_RADI(iz)-1, iz, -Zone(iz)%isfr(2)
              write (iu, *) 0, ZON_POLO(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
     enddo
  enddo
  ! non-transparent poloidal and toroidal surfaces
  write (iu, 1013)
  write (iu, *) 0
  write (iu, 1014)
  write (iu, *) 0
 1012 format ('* radial')
 1013 format ('* poloidal')
 1014 format ('* toroidal')
 1015 format (2i8,4x,'EIRENE_SF',i0)


  ! 2. additional physical cells for neutrals
  !ne0 = 1.d14;	Te0 = 4.d3;	Ti0 = 4.d3;	M0  = 0.d0
  !ne1 = 1.d7;	Te1 = 1.d-1;	Ti1 = 1.d-1;	M1  = 0.d0
  write (iu, 9999)
  write (iu, 2000)
  write (iu, 9999)
  write (iu, 2001)
  write (iu, 2002)
  n = 0
  do irun=0,1
     ! write number of additional cell blocks
     if (irun == 1) write (iu, *) n, 70

     ! confined region
     do iz=0,NZONET-1
        if (Zone(iz)%isfr(1) == SF_CORE) then
           if (irun == 0) then
              n = n + R_SURF_PL_TRANS_RANGE(1,iz)
           else
              do ir=0,R_SURF_PL_TRANS_RANGE(1,iz)-1
                 write (iu, *) 2, 1
                 write (iu, 2003) iz, ir, ir+1, 1, &
                           0, ZON_POLO(iz), ZON_POLO(iz), &
                           0, ZON_TORO(iz), ZON_TORO(iz)
                 write (iu, 2004) iz, ir
              enddo
           endif
        endif
     enddo

     ! vacuum region (collect what's left)
     do iz=0,NZONET-1
        if (irun == 0) then
           n = n + 1
        else
           write (iu, *) 2, 0
           write (iu, 2003) iz, 0, SRF_RADI(iz)-1, 1, &
                     0, ZON_POLO(iz), ZON_POLO(iz), &
                     0, ZON_TORO(iz), ZON_TORO(iz)
           write (iu, 2005) iz
        endif
     enddo
  enddo
 2000 format ('*** 2. DEFINE ADDITIONAL PHYSICAL CELLS FOR NEUTRALS')
 2001 format ('*   ZONE  R1    R2    DR    P1    P2    DP    T1    T2    DT')
 2002 format ('* ne       Te      Ti        M')
 2003 format (10i6)
 2004 format ('EIRENE_CORE_',i0,'_',i0)
 2005 format ('EIRENE_VACUUM_',i0)


  ! 3. Neutral sources
  write (iu, 9999)
  write (iu, 3000)
  write (iu, 9999)
  write (iu, 3001)
  write (iu, *) 0, 0, 1
 3000 format ('*** 3. Neutral Source distribution')
 3001 format ('* N0S NS_PLACE  NSSIDE')


  ! 4. Additional surfaces
  write (iu, 9999)
  write (iu, 4000)
  write (iu, 9999)
  write (iu, 4001)
  close (iu)
 4000 format ('*** 4 Additional surfaces')
 4001 format ('./../../geometry/ADD_SF_N0')
  close (iu)
 9999 format ('*',32('-'))

  end subroutine write_input_n0g
  !---------------------------------------------------------------------
  subroutine write_input_par

  integer :: ir, iz, irun, n


  open  (iu, file='input.PAR.6')
  write (iu, 9999)
  write (iu, 1000)
  write (iu, 9999)
 1000 format ('*** 6. boundary contitions')

  ! 1. particle transport
  write (iu, 1001)
 1001 format ('*** 6.1 particle transport for main ions + impurities')
  write (iu, 1002)
  write (iu, 1003)
  do irun=0,1
     ! write number of non default radial surfaces
     if (irun == 1) write (iu, *) n
     n = 0

     do iz=0,NZONET-1
        if (Zone(iz)%isfr(1) < 0) then
           n = n + 1
           if (irun == 1) then
              write (iu, 1004) 1, n, -Zone(iz)%isfr(1), -Zone(iz)%isfr(1)
           endif
        endif
        if (Zone(iz)%isfr(2) < 0) then
           n = n + 1
           if (irun == 1) then
              write (iu, 1004) 1, n, -Zone(iz)%isfr(2), -Zone(iz)%isfr(2)
           endif
        endif
     enddo
  enddo
 1002 format ('* Main plasma ions')
 1003 format ('NBUND_TYE    CBUND_COE')
 1004 format (2i4,4x,'EMC3_SF',i0,4x,'EMC3_SF',i0,'_PVAL')

  ! 2. energy transport
  write (iu, 2001)
 2001 format ('*** 6.2 energy transport for el. + ions')
  do irun=0,1
     ! write number of non default radial surfaces
     if (irun == 1) write (iu, *) n
     n = 0

     do iz=0,NZONET-1
        if (Zone(iz)%isfr(1) < 0) then
           n = n + 1
           if (irun == 1) then
              write (iu, 2002) 1, n, -Zone(iz)%isfr(1), -Zone(iz)%isfr(1), -Zone(iz)%isfr(1)
           endif
        endif
        if (Zone(iz)%isfr(2) < 0) then
           n = n + 1
           if (irun == 1) then
              write (iu, 2002) 1, n, -Zone(iz)%isfr(2), -Zone(iz)%isfr(2), -Zone(iz)%isfr(2)
           endif
        endif
     enddo
  enddo
 2002 format (2i4,4x,'EMC3_SF',i0,4x,'EMC3_SF',i0,'_EEVAL',4x,'EMC3_SF',i0,'_EIVAL')

  ! 3. momentum transport
  write (iu, 3001)
 3001 format ('*** 6.3 momentum transport')
  do irun=0,1
     ! write number of non default radial surfaces
     if (irun == 1) write (iu, *) n
     n = 0

     do iz=0,NZONET-1
        if (Zone(iz)%isfr(1) < 0) then
           n = n + 1
           if (irun == 1) then
              write (iu, 3002) 1, n, -Zone(iz)%isfr(1), -Zone(iz)%isfr(1)
           endif
        endif
        if (Zone(iz)%isfr(2) < 0) then
           n = n + 1
           if (irun == 1) then
              write (iu, 3002) 1, n, -Zone(iz)%isfr(2), -Zone(iz)%isfr(2)
           endif
        endif
     enddo
  enddo
 3002 format (2i4,4x,'EMC3_SF',i0,4x,'EMC3_SF',i0,'_MVAL')

  close (iu)
 9999 format ('*',32('-'))

  end subroutine write_input_par
  !---------------------------------------------------------------------
  end subroutine write_emc3_input_files
!=======================================================================



!=======================================================================
! WRITE BASE GRID
!=======================================================================
  subroutine write_base_grid(G, iz)
  use grid
  type(t_grid), intent(in) :: G
  integer, intent(in)      :: iz

  character(len=72) :: filename


  ! write grid file for field line tracing
  write (filename, 9000) iz
  call G%store(filename=filename)

  ! write grid file for plotting
  write (filename, 9001) iz
  call G%plot_mesh(filename)

 9000 format ('base_grid_',i0,'.dat')
 9001 format ('base_grid_',i0,'.plt')
  end subroutine write_base_grid
!=======================================================================



!=======================================================================
! GENERATE PLATE DEFINITIONS
!=======================================================================
  subroutine generate_plates()
  use iso_fortran_env
  use emc3_grid
  use boundary
  use curve2D
  use math
  use run_control, only: Debug
  use string
  use dataset
  implicit none

  integer, parameter :: iu = 78
  real(real64), parameter :: l0 = 10.d-10

  type(t_curve), dimension(:,:), allocatable :: C

  real(real64), dimension(:), allocatable   :: RC_TEM, ZC_TEM
  integer, dimension(:), allocatable   :: ID_TEM
  integer, dimension(:,:), allocatable :: iindex ! number of cells behind a plate
  integer, dimension(:), allocatable   :: knumb  ! cell index in flux tube for plate cells

  logical      :: plate_cell
  real(real64) :: x(3), phi
  integer      :: nr, np, nt, iz, i, j, k, l, l1, l2, irun, icut, ig(8), ic


  write (6, 9000)
 9000 format(3x,'- Set up plate surfaces')


  allocate (ID_TEM(0:MESH_P_OS(NZONET)-1))
  allocate (RC_TEM(0:MESH_P_OS(NZONET)-1))
  allocate (ZC_TEM(0:MESH_P_OS(NZONET)-1))
  ID_TEM = 0
  open  (iu, file='plates.dat')
  do iz=0,NZONET-1
     nr = ZON_RADI(iz)
     np = ZON_POLO(iz)
     nt = ZON_TORO(iz)
     allocate (iindex(0:nr-1, 0:np-1))
     iindex = 0

     write (6, 1000) iz, nr, np, nt
 1000 format (8x,'zone ',i0,' (',i0,' x ',i0,' x ',i0,')')

     ! setup slices for Q4-type surfaces
     if (n_quad > 0) then
     allocate (C(0:nt-1, n_quad))
     do k=0,nt-1
     do l=1,n_quad
         phi  = (PHI_PLANE(k+1+PHI_PL_OS(iz)) + PHI_PLANE(k+PHI_PL_OS(iz))) / 2.d0
         phi  = phi_sym(phi, symmetry)
         C(k,l) = S_quad(l)%slice(phi)

         if (Debug) then
            call C(k,l)%plot(filename='debug/Q4surf_'//trim(str(l))//'_zone'//trim(str(iz))// &
                                      '_t'//trim(str(k)))
         endif
     enddo
     enddo
     endif


     ! (0) count plate cells, (1) setup plate cells
     do irun=0,1
        icut = 0
        ! loop over all flux tubes
        do i=R_SURF_PL_TRANS_RANGE(1,iz),R_SURF_PL_TRANS_RANGE(2,iz)-1
        do j=P_SURF_PL_TRANS_RANGE(1,iz),P_SURF_PL_TRANS_RANGE(2,iz)-1
        do k=0,nt-1
           ig(1)   = i + (j + k*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
           ig(2)   = ig(1) + 1
           ig(3)   = ig(2) + SRF_RADI(iz)
           ig(4)   = ig(3) - 1
           ig(5:8) = ig(1:4) + SRF_POLO(iz)*SRF_RADI(iz)

           x(1)    = sum(RG(ig))/8.d0
           x(2)    = sum(ZG(ig))/8.d0
           x(3)    = (PHI_PLANE(k+1+PHI_PL_OS(iz)) + PHI_PLANE(k+PHI_PL_OS(iz))) / 2.d0
           x(3)    = phi_sym(x(3), symmetry)

           ic      = i + (j + k*ZON_POLO(iz))*ZON_RADI(iz) + MESH_P_OS(iz)
           RC_TEM(ic) = x(1)
           ZC_TEM(ic) = x(2)

           plate_cell = outside_boundary()

           if (plate_cell) then
              ID_TEM(ic) = 1
              icut = icut + 1
              if (irun == 1) then
                 iindex(i,j) = iindex(i,j) + 1
                 knumb(icut) = k
              endif
           endif
        enddo
        enddo
        enddo

        ! allocate knumb array after 1st run
        if (irun == 0) then
           allocate (knumb(icut))
           knumb = 0

        ! write plate cells after 2nd run
        else
           icut = 0
           do i=0,nr-1
           do j=0,np-1
              if (iindex(i,j) .ne. 0) then
                 l1 = icut + 1
                 l2 = icut + iindex(i,j)
                 write (iu, *) iz, i, j, iindex(i,j),(knumb(l),l=l1,l2)
                 icut = icut + iindex(i,j)
              endif
           enddo
           enddo

           write (6, 2000) icut
           2000 format (8x,i0,' plate cells')
           deallocate (knumb, iindex)
        endif
     enddo

     ! cleanup slices of Q4-type surfaces
     if (n_quad > 0) then
     do k=0,nt-1
     do l=1,n_quad
        call C(k,l)%destroy()
     enddo
     enddo
     deallocate (C)
     endif
  enddo
  close (iu)


  ! additional output for debugging
  if (Debug) then
     do iz=0,NZONET-1
     do k=0,ZON_TORO(iz)-1
        open  (99, file='debug/SlicePlasma_zone'//trim(str(iz))//'_t'//trim(str(k))//'.plt')
        open  (98, file='debug/SliceBoundary_zone'//trim(str(iz))//'_t'//trim(str(k))//'.plt')

        do i=R_SURF_PL_TRANS_RANGE(1,iz),R_SURF_PL_TRANS_RANGE(2,iz)-1
        do j=P_SURF_PL_TRANS_RANGE(1,iz),P_SURF_PL_TRANS_RANGE(2,iz)-1
           ic      = i + (j + k*ZON_POLO(iz))*ZON_RADI(iz) + MESH_P_OS(iz)
           if (ID_TEM(ic) == 1) then
              write (98, *) RC_TEM(ic), ZC_TEM(ic)
           else
              write (99, *) RC_TEM(ic), ZC_TEM(ic)
           endif
        enddo
        enddo
        close (99)
        close (98)
     enddo
     enddo

     call plate_check_all()
  endif
  deallocate (ID_TEM, RC_TEM, ZC_TEM)
  contains
  !-------------------------------------------------------------------
  function outside_boundary()
  ! input:
  !    k: toroidal index
  !    C: slices of Q4-type surfaces at phi(k)
  !    x: reference point
  logical :: outside_boundary, side


  ! set default
  outside_boundary = .false.


  ! check axisymmetric (L2-type) surfaces
  do l=1,n_axi
     side = boundary_side(l) == 1
     if (S_axi(l)%outside(x(1:2)) .eqv. side) then
        outside_boundary = .true.
        return
     endif
  enddo

  ! check Q4-type surfaces
  do l=1,n_quad
     if (C(k,l)%n_seg <= 0) cycle
     side = boundary_side(n_axi + n_block + l) == 1
     if (C(k,l)%outside(x(1:2)) .eqv. side) then
        outside_boundary = .true.
        return
     endif
  enddo

  ! check block limiters (CSG-type)
  do l=1,n_block
     if (bl_outside(l, x)) then
        outside_boundary = .true.
        return
     endif
  enddo

  end function outside_boundary
  !-------------------------------------------------------------------


  !-------------------------------------------------------------------
  subroutine outside_boundary_check

  integer, parameter :: iu = 99, iu2 = 98, n = 101

  character(len=72) :: s
  real(real64) :: x(2), Phi
  integer :: i

  allocate (C(0:n-1, n_quad))
  open  (iu, file='fl+.dat')
  open  (iu2, file='fl+B.dat')
  read  (iu, *) s
  do i=0,n-1
     read (iu, *) x, Phi

     phi  = phi_sym(phi, symmetry)
     C(i,1) = S_quad(1)%slice(phi)

     if (C(i,1)%outside(x)) then
        !write (iu2, *) Phi, 1
     else
        !write (iu2, *) Phi, 0
        write (iu2, *) x, Phi
     endif
  enddo
  close (iu)
  close (iu2)

  stop
  end subroutine outside_boundary_check
  !-------------------------------------------------------------------


  !-------------------------------------------------------------------
  subroutine plate_check_all
  integer :: iz, k, j, i, ic1, ic2, iplate


  write (6, *) 'running plate checks ...'

  ID_TEM = ID_TEM*2 - 1
  iplate = 0

  do iz=0,NZONET-1
  do i=0,ZON_RADI(iz)-1
  do j=0,ZON_POLO(iz)-1
     do k=0,ZON_TORO(iz)-2
        ic1 = i + (j +  k   *ZON_POLO(iz))*ZON_RADI(iz) + MESH_P_OS(iz)
        ic2 = i + (j + (k+1)*ZON_POLO(iz))*ZON_RADI(iz) + MESH_P_OS(iz)

        if (ID_TEM(ic1)*ID_TEM(ic2) < 0) then
           iplate = iplate + 1
           write (6, *) iplate
           call plate_check(iz,i,j,k+1)
        endif
     enddo
  enddo
  enddo
  enddo

  end subroutine plate_check_all
  !-------------------------------------------------------------------
  subroutine plate_check(iz, ir, ip, jt)
  use Q4
  use math
  use grid
  use fieldline
  use dataset
  integer, intent(in) :: iz, ir, ip, jt

  real(real64), parameter :: Limit = 360.d0

  type(t_Q4)   :: Q
  type(t_grid) :: G
  type(t_fieldline) :: F
  type(t_dataset)   :: D
  real(real64) :: x1(2), x2(2), x3(2), x4(2), phi, y(3), ts, Lc, Lcsav, dl
  integer      :: ig(4), i, idir, n, nsuccess


  ! calculate node indices
  ig(1)   = ir + (ip + jt*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
  ig(2)   = ig(1) + 1
  ig(3)   = ig(2) + SRF_RADI(iz)
  ig(4)   = ig(3) - 1

  ! get node coordinates
  x1(1)   = RG(ig(1)); x1(2)   = ZG(ig(1))
  x2(1)   = RG(ig(2)); x2(2)   = ZG(ig(2))
  x3(1)   = RG(ig(3)); x3(2)   = ZG(ig(3))
  x4(1)   = RG(ig(4)); x4(2)   = ZG(ig(4))
  phi     = PHI_PLANE(jt+PHI_PL_OS(iz))

  ! setup quadrilateral
  call Q%set_nodes(x1, x2, x3, x4)

  ! generate mesh
  n = 10
  G = Q%generate_mesh(n, n, phi)
  G%coordinates       = CYLINDRICAL
  !G%fixed_coord       = 3
  !G%fixed_coord_value = phi / 180.d0 * pi

  ! sample connection length on grid
  ts       = 1.d0
  Lcsav    = 0.d0
  nsuccess = 0
  call D%new(G%nodes(),5)
  do i=1,G%nodes()
     y = G%node(i)

     D%x(i,1) = y(1)
     D%x(i,2) = y(2)

     ! trace field line in both directions
     do idir=-1,1,2
        call F%init(y, idir*ts, NM_AdamsBashforth4, FL_ARC)

        trace_loop: do
           dl = F%trace_1step()
           Lc = F%phi_int * 180.d0 / pi

           if (abs(Lc) > Limit) exit trace_loop

           if (F%intersect_boundary()) exit trace_loop
        enddo trace_loop

        D%x(i,3 + (idir+1)/2) = Lc
     enddo
     D%x(i,5) = min(abs(D%x(i,3)), abs(D%x(i,4)))

     ! success frequency, average shortest toroidal distance
     if (D%x(i,5) < Limit) then
        nsuccess = nsuccess + 1
        Lcsav    = Lcsav    + D%x(i,5)
     endif
  enddo
  Lcsav = Lcsav / nsuccess
  write (6, *) 'success frequency = ', 100.d0 * nsuccess / G%nodes(), ' %'
  write (6, *) 'average toroidal distance to plates [deg] = ', Lcsav
  call G%store('plate1.grid')
  call D%plot(filename='plate1.dat')

  call D%destroy()
  call G%destroy()
  stop

  end subroutine plate_check
  !-------------------------------------------------------------------
  end subroutine generate_plates
!=======================================================================



!=======================================================================
! SAMPLE MAGNETIC FIELD STRENGTH ON GRID NODES (bfield.dat)
!=======================================================================
  subroutine sample_bfield_on_emc3_grid
  use iso_fortran_env
  use emc3_grid
  use bfield
  use equilibrium, only: get_PsiN
  use math, only: pi
  implicit none

  integer, parameter :: iu = 72

  real(real64) :: x(3), Bf(3)
  integer :: ir, ip, it, iz, ig


  write (6, *) 'sampling magnetic field strength on grid ...'
  if (.not.allocated(BFSTREN)) allocate(BFSTREN(0:GRID_P_OS(NZONET)-1))
  BFSTREN = 0.d0
  if (.not.allocated(PSI_N))   allocate(PSI_N(0:GRID_P_OS(NZONET)-1))
  PSI_N   = 0.d0


  do iz=0,NZONET-1
  do it=0,SRF_TORO(iz)-1
     write (6, *) iz, it
     x(3) = PHI_PLANE(it + PHI_PL_OS(iz))
     do ip=0,SRF_POLO(iz)-1
     do ir=R_SURF_PL_TRANS_RANGE(1,iz),R_SURF_PL_TRANS_RANGE(2,iz)
        ig = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)

        x(1) = RG(ig)
        x(2) = ZG(ig)
        Bf   = get_Bf_cyl(x)

        BFSTREN(ig) = sqrt(sum(Bf**2))
        PSI_N(ig)   = get_PsiN(x)
     enddo
     enddo
  enddo
  enddo
  ! convert units: Gauss -> T
  BFSTREN = BFSTREN / 1.d4


  open  (iu, file='bfield.dat')
  write (iu, *) BFSTREN
  close (iu)

  open  (iu, file='psiN.dat')
  write (iu, *) PSI_N
  close (iu)

  end subroutine sample_bfield_on_emc3_grid
!=======================================================================

end module fieldline_grid
