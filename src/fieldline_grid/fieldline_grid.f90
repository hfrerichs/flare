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
     TOPO_STEL   = 'stellarator', &
     TOPO_STEL1  = 'scstel', &
     TOPO_LSN    = 'lower_single_null', &
     TOPO_LSN1   = 'lsn', &
     TOPO_DDN    = 'disconnected_double_null', &
     TOPO_DDN1   = 'ddn', &
     TOPO_CDN    = 'connected_double_null', &
     TOPO_CDN1   = 'cdn', &
     TOPO_DSFP   = 'disconnected_snowflake_plus', &
     TOPO_DSFP1  = 'dsf+'


  ! Discretization type definitions
  character(len=*), parameter :: &
     POLOIDAL_ANGLE_FIXED = 'poloidal_angle_fixed', &
     ORTHOGONAL           = 'quasi_orthogonal', &
     ORTHOGONAL_STRICT    = 'quasi_orthogonal_strict', &
     POLOIDAL_ANGLE       = 'poloidal_angle', &
     ARC_LENGTH           = 'arc_length'


  ! zone type definitions
  integer, parameter :: &
     TYPE_SINGLE_LAYER = -1, &
     TYPE_HPR          =  0, &
     TYPE_IMD          =  0, &
     TYPE_SOL          =  1, &
     TYPE_SOLMAP       = 11, &
     TYPE_PFR          =  2, &
     TYPE_PFRMAP       = 22


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

  ! Method for core domain
  character(len=*), parameter :: &
     CORE_EXTRAPOLATE   = 'EXTRAPOLATE', &
     CORE_FLUX_SURFACES = 'FLUX_SURFACES'

  ! Plate definition methods
  character(len=*), parameter :: &
     PLATES_DEFAULT          = 'CELL_CENTER_OUTSIDE', &
     PLATES_RADIAL_INTERSECT = 'RADIAL_INTERSECTION'


  integer, parameter :: &
     max_blocks  = 360, &        ! Maximum number of toroidal blocks
     max_zones   = 360, &        ! Maximum number of zones
     max_layers  = 9             ! Maximum number of layer (zones per block)
!.......................................................................


!.......................................................................
! user defined variables
!.......................................................................
  character(len=256) :: &
     topology                         = TOPO_SC, &
     Innermost_Flux_Surface           = SF_EXACT, &
     x_in_coordinates                 = COORDINATES(2),&
     discretization_method            = POLOIDAL_ANGLE, &
     poloidal_discretization          = POLOIDAL_ANGLE, &
     radial_discretization            = POLOIDAL_ANGLE_FIXED, &
     radial_spacing(-1:max_layers-1)   = '', &
     poloidal_spacing(0:max_layers-1) = '', &
     poloidal_spacing_L(0:max_layers-1) = '', &
     poloidal_spacing_R(0:max_layers-1) = '', &
     toroidal_spacing(0:max_layers-1) = '', &
     guiding_surface                  = '', &
     N0_file(0:max_layers-1)          = '', &
     vacuum_domain(0:max_layers-1)    = '', &
     core_domain                      = CORE_EXTRAPOLATE, &
     plate_generator                  = PLATES_DEFAULT


  integer :: &
     symmetry            =   1, &          ! toroidal symmetry (i.e. 5 => 72 deg)
     blocks              =   1, &          ! number of toroidal blocks
     nt                  =  12, &          ! default toroidal resolution
     np(0:max_layers-1)  = 360, &          ! default poloidal resolution
     npL(0:max_layers-1) =  30, &          ! default poloidal resolution in divertor legs
     npR(0:max_layers-1) =  30, &          !    (L)eft and (R)ight segments
     nr(0:max_layers-1)  =  32, &          ! default radial resolution
     n_interpolate       =   4, &          ! number of interpolated flux surfaces (for the transition between the pair of perturbed flux surfaces at the inner simulation boundary and unperturbed flux surfaces further outside)
     np_ortho_divertor   =  10, &          ! number of orthogonal surfaces in divertor legs
     np_sub_divertor     =   2, &          ! sub-resolution in target aligned domain
     nr_EIRENE_core      =   1, &          ! radial resolution in core (EIRENE only)
     nr_EIRENE_vac(0:max_layers-1)       =   1, &          ! radial resolution in vacuum (EIRENE only)
     nr_perturbed        =   2, &          ! number of perturbed flux surfaces at the inner boundary
     plate_format        =   1             ! format for plate definition file

  real(real64) :: &
     phi0                = -360.d0, &      ! lower boundary of simulation domain [deg]
     x_in1(3)            = (/120.d0, 0.d0, 0.d0/), &  ! reference points (R[cm], Z[cm], phi[deg]) ...
     x_in2(3)            = (/119.d0, 0.d0, 0.d0/), &  ! ... on 1st and 2nd innermost flux surfaces
     d_SOL(2)            =   24.d0, &      ! radial width of scrape-off layer
     d_PFR(4)            =   15.d0, &      ! radial width of private flux region
     d_N0(0:max_layers-1)=   10.d0, &      ! radial width of vacuum region
     d_extend(0:max_layers-1,-1:1) = -1.d0, &      ! poloidal extension of divertor leg (used in close_grid_domain)
     d_cutL(2)           =    6.d0, &      ! cut-off length for flux surfaces behind the wall
     d_cutR(2)           =    8.d0, &      !    (L)eft and (R)ight segments
     alphaL(3)           =    0.9d0, &     ! Relative length of divertor legs at outermost boundary
     alphaR(3)           =    1.0d0, &     !    (L)eft and (R)ight segments
     etaL(3)             =    0.8d0, &     ! fraction of cells in front of the target
     etaR(3)             =    0.8d0, &     !    (L)eft and (R)ight segments
     Dtheta_sampling     =    20.d0, &     ! Transition between angular and length weighted sampling of flux surfaces
     Dtheta_separatrix   =     0.d0        ! ... same on separatrix

  logical :: &
     extend_alpha_SOL2   =  .true., &
     stellarator_symmetry = .false.



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


  ! toroidal discretization
  type t_toroidal_discretization
     integer :: nt, it_base
     real(real64), dimension(:), pointer :: phi => null()

     contains
     procedure :: init
     procedure :: setup => setup_toroidal_discretization
  end type t_toroidal_discretization


  ! zone data
  type :: t_zone
     ! resolution in radial, poloidal and toroidal direction
     integer :: nr, np, nt, nr_vac

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
     character(len=80) :: N0_file, vacuum_domain

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
  subroutine init(this, nt)
  class(t_toroidal_discretization) :: this
  integer,      intent(in)         :: nt


  this%nt      = nt
  if (associated(this%phi)) deallocate(this%phi)
  allocate(this%phi(0:nt))

  end subroutine init
!=======================================================================



!=======================================================================
  subroutine setup_toroidal_discretization(this, nt, it_base, phi)
  class(t_toroidal_discretization) :: this
  integer,      intent(in)         :: nt, it_base
  real(real64), intent(in)         :: phi(0:nt)


  call this%init(nt)
  this%it_base = it_base
  this%phi     = phi

  end subroutine setup_toroidal_discretization
!=======================================================================



!=======================================================================
! Set up zone for block "iblock" and layer "ilayer".
! Zone layout itypeR = TYPE_SINGLE_LAYER, TYPE_HPR, TYPE_SOL, TYPE_PFR
!             itypeP = SF_PERIODIC, SF_VACUUM
!=======================================================================
  subroutine setup(this, iblock, ilayer, itypeR, itypeP)
  class(t_zone)       :: this
  integer, intent(in) :: iblock, ilayer, itypeR, itypeP

  integer :: nr_add1 = 0, nr_add2 = 0, np_add = 0


  ! 1. grid resolution
  call setup_resolution()


  ! 2. default toroidal discretization
  allocate (this%phi(0:this%nt))
  this%phi     = Block(iblock)%phi
  this%it_base = Block(iblock)%it_base


  ! 3.1 plamsa transport domain
  this%r_surf_pl_trans_range(1) = nr_add1
  this%r_surf_pl_trans_range(2) = this%nr - nr_add2
  this%p_surf_pl_trans_range(1) = np_add
  this%p_surf_pl_trans_range(2) = this%np - np_add

  ! 3.2 set parameters for additional neutral domain
  this%d_N0      = d_N0(ilayer)
  this%d_extend  = d_extend(ilayer,-1:1)
  this%N0_file   = N0_file(ilayer)
  this%vacuum_domain = vacuum_domain(ilayer)
  this%nr_vac    = nr_EIRENE_vac(ilayer)


  ! 4. boundaries and zone type
  this%itypeR = itypeR
  this%itypeP = itypeP
  call setup_boundaries()

  contains
  !---------------------------------------------------------------------
  subroutine setup_resolution()


  ! toroidal resolution
  this%nt      = Block(iblock)%nt


  ! radial resolution (depends on type of zone)
  nr_add1 = 0; nr_add2 = 0
  select case(itypeR)
  ! single layer: add core domain on lower and vaccum domain on upper radial boundary
  case(TYPE_SINGLE_LAYER)
     nr_add1 = nr_EIRENE_core
     nr_add2 = nr_EIRENE_vac(ilayer)

  ! high pressure region (HPR): add core domain on lower radial boundary
  case(TYPE_HPR)
     nr_add1 = nr_EIRENE_core

  ! scrape-off layer (SOL): add vacuum domain on upper radial boundary
  case(TYPE_SOL)
     nr_add2 = nr_EIRENE_vac(ilayer)

  ! private flux region (PFR): add vacuum domain on lower radial boundary
  case(TYPE_PFR)
     nr_add1 = nr_EIRENE_vac(ilayer)
  end select
  this%nr = Block(iblock)%nr(ilayer) + nr_add1 + nr_add2


  ! poloidal resolution (depends on type of zone)
  np_add = 0
  !if (itypeP == SF_VACUUM) np_add = 1
  this%np = Block(iblock)%np(ilayer) + 2*np_add

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

  case(TYPE_SOLMAP, TYPE_PFRMAP)
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
  use run_control, only: run_control_file
  use equilibrium, only: get_cylindrical_coordinates

  integer, parameter :: iu = 12

  type(t_block_input) :: Block(0:max_blocks-1)
  real(real64) :: tmp(3)
  integer      :: ierr

  namelist /FieldlineGrid_Input/ &
     topology, symmetry, stellarator_symmetry, blocks, Block, &
     phi0, x_in1, x_in2, x_in_coordinates, d_SOL, d_PFR, d_N0, N0_file, vacuum_domain, d_extend, &
     nt, np, npL, npR, nr, nr_EIRENE_core, nr_EIRENE_vac, core_domain, &
     n_interpolate, nr_perturbed, plate_generator, plate_format, &
     np_ortho_divertor, np_sub_divertor, &
     radial_spacing, poloidal_spacing, poloidal_spacing_L, poloidal_spacing_R, toroidal_spacing, &
     d_cutL, d_cutR, etaL, etaR, alphaL, alphaR, extend_alpha_SOL2, &
     Dtheta_sampling, Dtheta_separatrix, &
     discretization_method, poloidal_discretization, radial_discretization, guiding_surface


  ! 1. read user configuration from input file
  open  (iu, file=run_control_file, err=9000)
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
  select case(x_in_coordinates)
  ! cylindrical coordinates
  case(COORDINATES(2))
     ! nothing to be done here

  ! toroidal / flux coordinates
  case(COORDINATES(3))
     write (6, 1000)
     tmp = x_in1;  x_in1 = get_cylindrical_coordinates(tmp, ierr)
     if (ierr == 0) then
        write (6, 1001) tmp, x_in1
     else
        write (6, 9001) tmp;  stop
     endif

     tmp = x_in2;  x_in2 = get_cylindrical_coordinates(tmp, ierr)
     if (ierr == 0) then
        write (6, 1001) tmp, x_in2
     else
        write (6, 9001) tmp;  stop
     endif

  case default
     write (6, *) 'error: invalid x_in_coordinates = ', trim(x_in_coordinates)
     stop
  end select


  ! 2. setup size and position of toroidal blocks
  call setup_toroidal_blocks(Block)


  ! 3. setup mesh spacing functions
  call setup_mesh_spacing()


  return
 1000 format(3x,'- Calculating cylindrical coordinates for reference points:')
 1001 format(8x,'( ',f0.3,', ',f0.3,', ',f0.3,')  ->  ( ',f0.3,', ',f0.3,', ',f0.3,')')
 9001 format('error: cannot convert ( ',f0.3,', ',f0.3,', ',f0.3,') to cylindrical coordinates!')
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

        ! default base plane in 1st block in stellarator symmetric configurations is 0
        if (stellarator_symmetry  .and.  ib==0) Block(ib)%it_base = 0
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
  if (stellarator_symmetry) Delta_phi_sim = Delta_phi_sim / 2


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
  ! DEFAULT: neg. half of simulation domain
  if (phi0 == -360.d0) then
     phi0 = -Delta_phi_sim / 2.d0

     ! default lower boundary is 0.0 deg for stellarator symmetric configurations
     if (stellarator_symmetry) phi0 = 0.d0
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
  write (iu, *) -plate_format ! user defined
  write (iu, 2003)
  write (iu, *) -plate_format ! user defined
  write (iu, 2004)
  write (iu, *) -plate_format ! user defined
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
  subroutine check_emc3_grid()
  use emc3_grid

  real(real64) :: R, R0, Z, B, B0
  integer      :: iz, ir, ip, it(2), jt, ig, igc(2)


  R0 = sum(RG) / GRID_P_OS(NZONET)
  if (allocated(BFSTREN)) B0 = sum(BFSTREN) / GRID_P_OS(NZONET)
  ! check up/down symmetry
  do iz=0,NZONET-1
     it(1) = 0
     it(2) = SRF_TORO(iz)-1
     do jt=1,2
     if (Zone(iz)%isft(jt) == SF_UPDOWN) then
        do ir=0,SRF_RADI(iz)-1
           ig = ir + it(jt)*SRF_POLO(iz)*SRF_RADI(iz) + GRID_P_OS(iz)
           do ip=0,ZON_POLO(iz)/2
              igc(1) = ig +               ip  * SRF_RADI(iz)
              igc(2) = ig + (ZON_POLO(iz)-ip) * SRF_RADI(iz)

              ! check geometry
              if (abs(RG(igc(1))-RG(igc(2)))/R0 > 1.d8  .or. &
                  abs(ZG(igc(1))-RG(igc(2)))/R0 > 1.d8) then
                 write (6, *) 'error: deviation from up/down symmetry too large!'
                 write (6, *) 'at grid node iz,ir,ip,it = ', iz, ir, ip, it(jt)
                 write (6, *) 'R,Z = ', RG(igc(1)), ZG(igc(1))
                 write (6, *) 'R,Z = ', RG(igc(2)), ZG(igc(2))
                 stop
              else
                 R     = 0.5d0 * (RG(igc(1)) + RG(igc(2)))
                 Z     = 0.5d0 * (ZG(igc(1)) - ZG(igc(2)))
                 RG(igc(1)) = R;  ZG(igc(1)) = Z
                 RG(igc(2)) = R;  ZG(igc(2)) =-Z
              endif

              ! check field strength (only if it has been set up at this point)
              if (allocated(BFSTREN)) then
              if (abs(BFSTREN(igc(2))-BFSTREN(igc(1)))/B0 > 1.d8) then
                 write (6, *) 'error: deviation from up/down symmetry too large!'
                 write (6, *) 'at grid node iz,ir,ip,it = ', iz, ir, ip, it(jt)
                 write (6, *) 'R,Z,B = ', RG(igc(1)), ZG(igc(1)), BFSTREN(igc(1))
                 write (6, *) 'R,Z,B = ', RG(igc(2)), ZG(igc(2)), BFSTREN(igc(2))
                 stop
              else
                 B      = 0.5d0 * (BFSTREN(igc(1)) + BFSTREN(igc(2)))
                 BFSTREN(igc(1)) = B
                 BFSTREN(igc(2)) = B
              endif
              endif
           enddo
        enddo
     endif
     enddo
  enddo

  end subroutine check_emc3_grid
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


  call check_emc3_grid()
  open  (iu, file='bfield.dat')
  write (iu, *) BFSTREN
  close (iu)

  open  (iu, file='psiN.dat')
  write (iu, *) PSI_N
  close (iu)

  end subroutine sample_bfield_on_emc3_grid
!=======================================================================

end module fieldline_grid
