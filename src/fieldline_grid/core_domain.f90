!===============================================================================
! Set up core domain for EIRENE
!===============================================================================
  subroutine setup_core_domain()
  use emc3_grid, only: NZONET
  use fieldline_grid
  implicit none

  integer :: iz


  write (6, 1000)
  do iz=0,NZONET-1
  if (Zone(iz)%isfr(1) == SF_CORE) then
     select case(core_domain)
     case(CORE_EXTRAPOLATE)
        call setup_core_domain_extrapolate(iz, nr_EIRENE_core)

     case(CORE_FLUX_SURFACES)
        call setup_core_domain_from_flux_surfaces(iz, nr_EIRENE_core)

     case default
        write (6, 9000) trim(core_domain)
        stop
     end select
  endif
  enddo

 1000 format(3x,' - Setting up core domain for EIRENE')
 1001 format(8x,'by extrapolating from inner EMC3 boundary')
 1002 format(8x,'by generating flux surfaces')
 9000 format('error: invalid method "',a,'" for core domain!')
  end subroutine setup_core_domain
  !=====================================================================



  !=====================================================================
  subroutine setup_core_domain_extrapolate(iz, nr_core)
  use iso_fortran_env
  use emc3_grid
  use fieldline_grid, only: alpha_core
  implicit none

  integer, intent(in) :: iz, nr_core


  integer, parameter :: &
     MAGNETIC_AXIS0   = 1, &
     GEOMETRIC_CENTER = 2


  real(real64) :: phi, r0(2), x(2), dx(2)
  integer :: i, j, k, ig, method = GEOMETRIC_CENTER


  write (6, 1000) iz, nr_core
 1000 format(8x,'zone ',i0,': ',i0,' core cell(s)')

  do k=0,SRF_TORO(iz)-1
     phi = PHI_PLANE(k+PHI_PL_OS(iz))
     r0  = get_r0(phi)

     do j=0,SRF_POLO(iz)-1
        ig = (j + k*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
        dx(1)  = RG(ig+nr_core) - r0(1)
        dx(2)  = ZG(ig+nr_core) - r0(2)

        do i=0,nr_core-1
           x        = r0 + alpha_core * dx * (1.d0 + 1.d0 * i / nr_core)
           RG(ig+i) = x(1)
           ZG(ig+i) = x(2)
        enddo
     enddo
  enddo

  return
  contains
  !---------------------------------------------------------------------
  function get_r0(phi)
  use magnetic_axis, only: get_magnetic_axis
  real(real64), intent(in) :: phi
  real(real64)             :: get_r0(2)

  real(real64) :: r3(3), w, w_tot
  integer      :: j, ig2


  get_r0 = 0.d0

  select case(method)
  case(MAGNETIC_AXIS0)
     r3     = get_magnetic_axis(phi)
     get_r0 = r3(1:2)

  case(GEOMETRIC_CENTER)
     get_r0 = 0.d0
     w_tot  = 0.d0
     do j=0,ZON_POLO(iz)-1
        ig  = nr_core + (j + k*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
        ig2 = ig + SRF_RADI(iz)
        w   = sqrt((RG(ig)-RG(ig2))**2 + (ZG(ig)-ZG(ig2))**2)

        get_r0(1) = get_r0(1) + w*RG(ig)
        get_r0(2) = get_r0(2) + w*ZG(ig)
        w_tot     = w_tot     + w
     enddo
     get_r0 = get_r0 / w_tot
  end select

  end function get_r0
  !---------------------------------------------------------------------
  end subroutine setup_core_domain_extrapolate
  !=====================================================================



  !=====================================================================
  subroutine setup_core_domain_from_flux_surfaces(iz, nr_core)
  use iso_fortran_env
  use emc3_grid
  use fieldline_grid
  use magnetic_axis
  use flux_surface_3D
  use run_control, only: N_points, Trace_Method
  use mesh_spacing
  use grid
  use curve2D
  implicit none

  integer, intent(in) :: iz, nr_core

  type(t_flux_surface_3D) :: F
  type(t_spacing)         :: S
  type(t_grid)            :: B, G
  real(real64) :: x0(3), xmag(3), x(3), r
  integer      :: ig, ir, ip, ipc, it


  ! set default number of points for flux surfaces
  if (N_points == 0) N_points = 1000

  ! set poloidal discetization parameter
  select case(discretization_method)
  case(POLOIDAL_ANGLE)
     ipc = ANGLE
  case(ARC_LENGTH)
     ipc = DISTANCE
  end select


  ! set up reference point on inner EMC3 boundary (x0) and magnetic axis (xmag)
  ig = nr_core + Zone(iz)%it_base * SRF_POLO(iz) * SRF_RADI(iz)   +  GRID_P_OS(iz)
  x0(1) = RG(ig)
  x0(2) = ZG(ig)
  x0(3) = Zone(iz)%phi(Zone(iz)%it_base) / 180.d0 * pi
  xmag  = get_magnetic_axis(x0(3))
  call S%init(radial_spacing(-1))


  ! initialize base grid for core region
  call B%new(CYLINDRICAL, MESH_2D, 3, nr_core, SRF_POLO(iz), fixed_coord_value=x0(3))


  ! generate flux surfaces for core region -> setup base grid
  do ir=0,nr_core-1
     r = S%node(ir,nr_core)
     r = alpha_core + (1.d0 - alpha_core) * r
     x = xmag + r * (x0-xmag)

     call F%generate(x, N_points, symmetry, 1, 360, Trace_Method, poloidal_coordinate=ipc, resample=SRF_POLO(iz))
     do ip=0,SRF_POLO(iz)-1
        B%mesh(ir,ip,:) = F%slice(0)%x(ip,:)
     enddo
  enddo
  call B%store('core_flux_surfaces.plt')


  ! trace nodes from base grid and set up 3D grid
  call trace_nodes_1zone(iz, B, G)
  do it=0,SRF_TORO(iz)-1
  do ip=0,SRF_POLO(iz)-1
  do ir=0,nr_core-1
     ig = ir + (ip + it * SRF_POLO(iz)) * SRF_RADI(iz)   +  GRID_P_OS(iz)
     RG(ig) = G%mesh3D(ir,ip,it,1)
     ZG(ig) = G%mesh3D(ir,ip,it,2)
  enddo
  enddo
  enddo

  end subroutine setup_core_domain_from_flux_surfaces
  !=====================================================================
