module elements
  use iso_fortran_env
  use curve2D, only: t_curve, UPPER, LOWER
  use mesh_interface
  use mfs_mesh
  use fieldline_grid, only: t_toroidal_discretization
  implicit none
  private


  integer, parameter, public :: &
     PERIODIC       = -1, &
     CORE           = -2, &
     VACUUM         = -3, &
     DIVERTOR       = -4


  ! grid layout: element = structured domain of cells
  type, public :: t_element
     ! element id
     integer :: id

     ! mapping to neighbor elements
     integer :: map_r(-1:1) = UNDEFINED, map_p(-1:1) = UNDEFINED
     ! boundary geometry
     integer :: rad_bound(-1:1) = UNDEFINED, pol_bound(-1:1) = UNDEFINED
     !integer           :: neighbor_id(4), neighbor_surf(4)
     !integer           :: mapping(4)

     ! radial and poloidal resolution in element
     integer :: nr = UNDEFINED, np = UNDEFINED
     ! radial layer index
     integer :: irl = UNDEFINED
     ! poloidal layer/block index
     integer :: ipl = UNDEFINED, ipl_side = 0

     ! toroidal discretization
     type(t_toroidal_discretization) :: T

     !type(t_mesh_interface) :: generating_element
     contains
     procedure :: setup_mapping
     procedure :: setup_boundary
     procedure :: generate_mesh
     procedure :: define
  end type t_element
  type(t_element), dimension(:), allocatable, public :: Z
  integer, public :: nelement


  public :: initialize_elements
  public :: setup_elements
  public :: undefined_element_boundary_check

  contains
!=======================================================================



!=======================================================================
! define surface mapping
!=======================================================================
  subroutine define(this, surface, mapping)
  class(t_element)       :: this
  integer, intent(inout) :: surface
  integer, intent(in)    :: mapping


  if (surface == UNDEFINED) then
     surface = mapping
  else
     write (6, *) 'error: mapping/boundary surface is alread defined!'
     write (6, *) 'element id = ', this%id
     write (6, *) 'surface = ', surface
     write (6, *) 'value   = ', mapping
     stop
  endif

  end subroutine define
!=======================================================================



!=======================================================================
! set up mapping between elements
!=======================================================================
  subroutine setup_mapping(this, side, boundary, element, interface_id)
  class(t_element)               :: this
  integer,         intent(in)    :: side, boundary
  type(t_element), intent(inout) :: element
  integer,         intent(in), optional :: interface_id

  integer :: id


  ! set default values for optional input
  id = UNDEFINED
  if (present(interface_id)) id = interface_id


  ! check input for side
  if (side.ne.UPPER .and. side.ne.LOWER) then
     write (6, *) 'error in t_element%setup_mapping: invalid side = ', side
     stop
  endif


  ! set up mapping between elements
  select case(boundary)
  case(RADIAL)
     if (id /= UNDEFINED) then
     if (id < 1  .or.  id > radial_interfaces) then
        write (6, *) 'error in t_element%setup_mapping: invalid rad. interface id = ', id
        stop
     endif
     endif
     call this%define(this%map_r(side),  element%id) ! map this element to neighbor element
     call this%define(element%map_r(-side), this%id) ! set up return map
     call this%define(this%rad_bound(side),  id)
     call this%define(element%rad_bound(-side), id)
  case(POLOIDAL)
     if (id /= UNDEFINED) then
     if (id < 1  .or.  id > poloidal_interfaces) then
        write (6, *) 'error in t_element%setup_mapping: invalid pol. interface id = ', id
        stop
     endif
     endif
     call this%define(this%map_p(side),  element%id) ! map this element to neighbor element
     call this%define(element%map_p(-side), this%id) ! set up return map
     call this%define(this%pol_bound(side),  id)
     call this%define(element%pol_bound(-side), id)
  case default
     write (6, *) 'error in t_element%setup_mapping: invalid boundary = ', boundary
     stop
  end select


  end subroutine setup_mapping
!=======================================================================



!=======================================================================
! set up element boundary (boundary_type = CORE, VACUUM, PERIODIC, DIVERTOR)
!=======================================================================
  subroutine setup_boundary(this, side, boundary, boundary_type, iinterface)
  class(t_element)            :: this
  integer,      intent(in)    :: side, boundary, boundary_type
  integer,      intent(in), optional :: iinterface


  ! check input for side
  if (side.ne.UPPER .and. side.ne.LOWER) then
     write (6, *) 'error in t_element%setup_boundary: invalid side = ', side
     stop
  endif


  ! check input for boundary_type
  select case(boundary_type)
  case(PERIODIC,CORE,VACUUM,DIVERTOR)
  case default
     write (6, *) 'error in t_element%setup_boundary: invalid boundary_type = ', boundary_type
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
     write (6, *) 'error in t_element%setup_boundary: invalid boundary = ', boundary
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
  subroutine generate_mesh(this, M, irside, ipside, iblock, Sr, Sp, U, debug)
  use fieldline_grid, only: n_interpolate, np_ortho_divertor, Zone
  use inner_boundary, only: C_in, DPsiN1
  use mesh_spacing
  class(t_element)                :: this
  type(t_mfs_mesh), intent(inout) :: M
  integer,          intent(in)    :: irside, ipside, iblock
  type(t_spacing),  intent(in)    :: Sr, Sp
  type(t_curve),    intent(in)    :: U
  logical,          intent(in), optional :: debug

  integer :: ir0, iri, iri2, ix(-1:1), ix2(-1:1), iside, addX(2), np0, np1, ierr
  integer :: npA_range(2)


  iri = this%rad_bound(irside)
  if (iri == UNDEFINED) then
     write (6, 9000);  write (6, 9003);  stop
  endif
  ix  = radial_interface(iri)%inode


  ! 0. check if generating boundary nodes are set
  ir0 = M%ir0
  if (M%ir0 < 0) then
     ! check if interface geometry is defined
     if (radial_interface(iri)%geometry_undefined()) then
        write (6, 9000);  write (6, 9002) ierr;  stop
     endif
     write (6, *) 'setting boundary nodes from radial interface ', iri
     call M%setup_boundary_nodes(irside, RADIAL, radial_interface(iri)%C, Equidistant, debug=debug)
  endif
  ! quasi-orhogonal mesh in upstream divertor legs
  np0  = np_ortho_divertor
  np1  = 0
  npA_range = 0
  ! transition PFR between two X-points
  iri2 = this%rad_bound(-irside)
  if (iri2 /= UNDEFINED) then
     iside = 0
     if (ix(LOWER) == STRIKE_POINT) iside = LOWER
     if (ix(UPPER) == STRIKE_POINT) iside = UPPER
     ix2 = radial_interface(iri2)%inode

     if (iside /= 0) then
        if (ix2(-iside) > 0) then
           np1 = 1
           np0 = max(np0,np1)
           write (6, *) 'TEST'
        endif
     endif
  endif


  ! 1. strike point on lower poloidal side
  if (ix(LOWER) == STRIKE_POINT) then
     write (6, *) 'strike point on lower poloidal side'
     if (np0 > 0) then
        npA_range = (/this%np-np0, this%np-1-np1/)
!        call M%make_orthogonal_grid(prange=(/this%np-np0, this%np-1-np1/), debug=debug)
     endif
!     if (np1 > 0) then
!        call M%make_interpolated_submesh((/0,this%nr-1/), (/this%np-1-np1, this%np/))
!     endif
     call M%make_divertor_grid(UPPER, this%np-np0, npA_range, Sp, U, this%T, ir0, ierr)
     if (ierr /= 0) then
        write (6, 9000);  write (6, 9001) ierr;  stop
     endif

  ! 2. strike point on upper poloidal side
  elseif (ix(UPPER) == STRIKE_POINT) then
     write (6, *) 'strike point on upper poloidal side'
     if (np0 > 0) then
        npA_range = (/1+np1,np0/)
!        call M%make_orthogonal_grid(prange=(/1+np1,np0/), debug=debug)
     endif
!     if (np1 > 0) then
!        call M%make_interpolated_submesh((/0,this%nr-1/), (/0, 1+np1/))
!     endif
     call M%make_divertor_grid(LOWER, np0, npA_range, Sp, U, this%T, ir0, ierr)
     if (ierr /= 0) then
        write (6, 9000);  write (6, 9001) ierr;  stop
     endif


  ! 3. intermediate domain
  else
     ! connect to X-point?
     addX = UNDEFINED
     iside = 0
     if (ix(-ipside) > 0) then
        addX(1) = ix(-ipside)
        addX(2) = AUTOMATIC
        write (6, *) 'upstream segment connecting X-point ', ix(ipside), ' to ', ix(-ipside)
        if (ix(-ipside)*ipside > ix(ipside)*ipside) iside = 1

     ! guide by X-point?
     elseif (iri2 > 0) then
        addX(1) = ix(-ipside)
        addX(2) = AUTOMATIC
        write (6, *) 'upstream segment guided by X-point ', -ix(-ipside)

     else
        write (6, *) 'upstream segment connecting to X-point ', ix(ipside)
     endif


     if (this%map_r(LOWER) == CORE) then
        call M%make_orthogonal_grid(rrange=(/2+n_interpolate, M%nr-1/), addX=addX)
        call M%make_interpolated_mesh(2+n_interpolate, Sr, C_in(iblock,:), DPsiN1(iblock,1))

     else
        call M%make_orthogonal_grid(addX=addX, debug=debug, side=iside)
     endif
  endif

 9000 format('error in t_element%generate_mesh:')
 9001 format('make_divertor_grid returned ierr = ', i0)
 9002 format('geometry of radial interface is undefined! ierr = ',i0)
 9003 format('radial interface is undefined!')
  end subroutine generate_mesh
!=======================================================================



!=======================================================================
! initialize element array, and set up element IDs
!=======================================================================
  subroutine initialize_elements(n)
  integer, intent(in) :: n

  integer :: i


  nelement = n
  allocate (Z(n))
  do i=1,n
     Z(i)%id = i
  enddo

  end subroutine initialize_elements
!=======================================================================



!=======================================================================
! automatic setup of element topology
!=======================================================================
  subroutine setup_elements(nX, connectX)
  integer, intent(in) :: nX, connectX(nX)


  integer :: ix, iz, jx, iri, ipi0


  iz   = 0
  ipi0 = 0
  iri  = 0
  do ix=1,nX
     jx = connectX(ix)

     ! max 8 elements per X-point, 2 for each GradPsiN direction

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

  end subroutine setup_elements
!=======================================================================



!=======================================================================
! check left over element boundaries
!=======================================================================
  subroutine undefined_element_boundary_check(debug)
  logical, intent(in) :: debug

  integer :: iz, iside


  do iz=1,nelement
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

 1001 format('element ', i0, ',   radial boundary ', i2, ' connects to element ', i0)
 1002 format('element ', i0, ', poloidal boundary ', i2, ' connects to element ', i0)
 9001 format('error: undefined radial boundary!')
 9002 format('error: undefined poloidal boundary!')
 9003 format('element = ', i0, ', iside = ', i0)
  end subroutine undefined_element_boundary_check
!=======================================================================

end module elements
