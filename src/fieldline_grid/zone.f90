module mod_zone
  use iso_fortran_env
  use curve2D, only: UPPER, LOWER
  use mesh_interface
  use mfs_mesh
  implicit none
  private


  integer, parameter, public :: &
     PERIODIC       = -1, &
     CORE           = -2, &
     VACUUM         = -3, &
     DIVERTOR       = -4, &
     UNDEFINED      = -100


  ! grid layout: zone = structured domain of cells
  ! This is not the same as a zone in EMC3! EMC3 zones are referred to as layer here, which can be made out of several zones connected in poloidal direction
  type, public :: t_zone
     ! zone number
     integer :: id

     ! zone neighbors
     integer :: map_r(-1:1) = UNDEFINED, map_p(-1:1) = UNDEFINED
     ! zone boundaries
     integer :: rad_bound(-1:1) = UNDEFINED, pol_bound(-1:1) = UNDEFINED
     !integer           :: neighbor_id(4), neighbor_surf(4)
     !integer           :: mapping(4)

     !type(t_mesh_interface) :: generating_element
     contains
     procedure :: setup_mapping
     procedure :: setup_boundary
     procedure :: generate_mesh
     procedure :: define
  end type t_zone
  type(t_zone), dimension(:), allocatable, public :: Z
  integer, public :: nzone


  public :: initialize_zones
  public :: setup_zones
  public :: undefined_zone_boundary_check

  contains
!=======================================================================



!=======================================================================
! define surface mapping
!=======================================================================
  subroutine define(this, surface, mapping)
  class(t_zone)          :: this
  integer, intent(inout) :: surface
  integer, intent(in)    :: mapping


  if (surface == UNDEFINED) then
     surface = mapping
  else
     write (6, *) 'error: mapping/boundary surface is alread defined!'
     write (6, *) 'zone id = ', this%id
     write (6, *) 'surface = ', surface
     write (6, *) 'value   = ', mapping
     stop
  endif

  end subroutine define
!=======================================================================



!=======================================================================
! set up mapping between zones
!=======================================================================
  subroutine setup_mapping(this, side, boundary, zone, iinterface)
  class(t_zone)               :: this
  integer,      intent(in)    :: side, boundary
  type(t_zone), intent(inout) :: zone
  integer,      intent(in)    :: iinterface


  ! check input for side
  if (side.ne.UPPER .and. side.ne.LOWER) then
     write (6, *) 'error in t_zone%setup_mapping: invalid side = ', side
     stop
  endif


  ! set up mapping between zones
  select case(boundary)
  case(RADIAL)
     call this%define(this%map_r(side),  zone%id) ! map this zone to neighbor zone
     call this%define(zone%map_r(-side), this%id) ! set up return map
     call this%define(this%rad_bound(side),  iinterface)
     call this%define(zone%rad_bound(-side), iinterface)
  case(POLOIDAL)
     call this%define(this%map_p(side),  zone%id) ! map this zone to neighbor zone
     call this%define(zone%map_p(-side), this%id) ! set up return map
     call this%define(this%pol_bound(side),  iinterface)
     call this%define(zone%pol_bound(-side), iinterface)
  case default
     write (6, *) 'error in t_zone%setup_mapping: invalid boundary = ', boundary
     stop
  end select


  ! set up interface
  ! ...

  end subroutine setup_mapping
!=======================================================================



!=======================================================================
! set up zone boundary (boundary_type = CORE, VACUUM, PERIODIC, DIVERTOR)
!=======================================================================
  subroutine setup_boundary(this, side, boundary, boundary_type, iinterface)
  class(t_zone)               :: this
  integer,      intent(in)    :: side, boundary, boundary_type
  integer,      intent(in), optional :: iinterface


  ! check input for side
  if (side.ne.UPPER .and. side.ne.LOWER) then
     write (6, *) 'error in t_zone%setup_boundary: invalid side = ', side
     stop
  endif


  ! check input for boundary_type
  select case(boundary_type)
  case(PERIODIC,CORE,VACUUM,DIVERTOR)
  case default
     write (6, *) 'error in t_zone%setup_boundary: invalid boundary_type = ', boundary_type
     stop
  end select


  ! set up domain boundary
  select case(boundary)
  case(RADIAL)
     call this%define(this%map_r(side), boundary_type)
  case(POLOIDAL)
     call this%define(this%map_p(side), boundary_type)
     if (present(iinterface)) call this%define(this%pol_bound(side),  iinterface)
  case default
     write (6, *) 'error in t_zone%setup_boundary: invalid boundary = ', boundary
     stop
  end select

  end subroutine setup_boundary
!=======================================================================



!=======================================================================
! Generate mesh from reference discretization on radial and poloidal
! boundaries. M%ir0 and M%ip0 must be set!
! irside = side of radial reference surface
! ipside = side of poloidal reference surface
!=======================================================================
  subroutine generate_mesh(this, M, irside, ipside, iblock, Sr)
  use fieldline_grid, only: n_interpolate
  use inner_boundary, only: C_in, DPsiN1
  use mesh_spacing
  class(t_zone)                   :: this
  type(t_mfs_mesh), intent(inout) :: M
  integer,          intent(in)    :: irside, ipside, iblock
  type(t_spacing),  intent(in)    :: Sr

  integer :: iri, ix(-1:1), addX(2)


  iri = this%rad_bound(irside)
  ix  = radial_interface(iri)%inode


  ! 1. strike point on lower poloidal side
  if (ix(LOWER) == STRIKE_POINT) then

  ! 2. strike point on upper poloidal side
  elseif (ix(UPPER) == STRIKE_POINT) then

  ! 3. intermediate domain
  else
     ! connect to X-point?
     addX = UNDEFINED
     if (ix(-ipside) > 0) then
        addX(1) = ix(-ipside)
        addX(2) = AUTOMATIC
     endif


     if (this%map_r(LOWER) == CORE) then
        call M%make_orthogonal_grid(rrange=(/2+n_interpolate, M%nr-1/), addX=addX)
        call M%make_interpolated_mesh(2+n_interpolate, Sr, C_in(iblock,:), DPsiN1(iblock,1))

     else
        call M%make_orthogonal_grid(addX=addX)
     endif
  endif

 9000 format('error in t_zone%generate_mesh:')
 9001 format('mesh generation for map_p = ', i0, ', ', i0, ' not implemented!')
  end subroutine generate_mesh
!=======================================================================



!=======================================================================
! initialize zone array, and set up zone IDs
!=======================================================================
  subroutine initialize_zones(n)
  integer, intent(in) :: n

  integer :: i


  nzone = n
  allocate (Z(n))
  do i=1,n
     Z(i)%id = i
  enddo

  end subroutine initialize_zones
!=======================================================================



!=======================================================================
! automatic setup of zone topology
!=======================================================================
  subroutine setup_zones(nX, connectX)
  integer, intent(in) :: nX, connectX(nX)


  integer :: ix, iz, jx, iri, ipi0


  iz   = 0
  ipi0 = 0
  iri  = 0
  do ix=1,nX
     jx = connectX(ix)

     ! max 8 zones per X-point, 2 for each GradPsiN direction

     ! 1. towards core
     ! 1.1 connect back to same X-point
     if (ix == jx) then
        iz = iz + 1
        !call Z(iz)%setup_boundary(LOWER, POLOIDAL, PERIODIC, ipi0+1) ! periodic pol. boundary
        !call Z(iz)%setup_boundary(UPPER, POLOIDAL, PERIODIC, ipi0+1) ! periodic pol. boundary
        call Z(iz)%setup_mapping (UPPER, POLOIDAL, Z(iz), 0)       ! connect to next segment
        call Z(iz)%setup_boundary(LOWER, RADIAL,   CORE)             ! core boundary
        call Z(iz)%setup_mapping (UPPER, RADIAL,   Z(iz+1),  iri)    ! connect to main SOL

     ! 1.2 connect to second X-point, or use second X-point as guiding point
     elseif (abs(jx) > ix) then
        iz = iz + 1
        call Z(iz)%setup_mapping (UPPER, POLOIDAL, Z(iz+1), 0)       ! connect to next segment
        call Z(iz)%setup_boundary(LOWER, RADIAL,   CORE)             ! core boundary
        call Z(iz)%setup_mapping (UPPER, RADIAL,   Z(iz+2),  iri)    ! connect to main SOL
     endif
  enddo

  end subroutine setup_zones
!=======================================================================



!=======================================================================
! check left over zone boundaries
!=======================================================================
  subroutine undefined_zone_boundary_check(debug)
  logical, intent(in) :: debug

  integer :: iz, iside


  do iz=1,nzone
     do iside=-1,1,2
        ! radial boundaries
        if (debug) write (6, 1001) iz, iside, Z(iz)%map_r(iside)
        if (Z(iz)%map_r(iside) == UNDEFINED) then
           write (6, 9001)
           write (6, 9003) iz, iside
           stop
        endif

        ! poloidal boundaries
        if (debug) write (6, 1002) iz, iside, Z(iz)%map_p(iside)
        if (Z(iz)%map_p(iside) == UNDEFINED) then
           write (6, 9002)
           write (6, 9003) iz, iside
           stop
        endif
     enddo
  enddo

 1001 format('zone ', i0, ',   radial boundary ', i2, ' connects to zone ', i0)
 1002 format('zone ', i0, ', poloidal boundary ', i2, ' connects to zone ', i0)
 9001 format('error: undefined radial boundary!')
 9002 format('error: undefined poloidal boundary!')
 9003 format('zone = ', i0, ', iside = ', i0)
  end subroutine undefined_zone_boundary_check
!=======================================================================

end module mod_zone
