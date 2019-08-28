!===============================================================================
! Set up core domain for EIRENE
!===============================================================================
  subroutine setup_core_domain()
  use fieldline, only: FULL_FIELD, EQUILIBRIUM_FIELD
  use emc3_grid, only: NZONET
  use fieldline_grid
  use string
  implicit none

  ! relative radial position of core surface
  real(real64), parameter :: alpha_core_default = 0.5d0

  character(len=len(core_domain)) :: command, argument, stmp
  real(real64) :: alpha_core
  integer      :: iz, alpha_type


  write (6, 1000)
  call read_command(core_domain, 1, command, argument)
  if (command == "USER_DEFINED") then
     call user_defined_core_domain(argument)
     return
  endif


  if (argument == 'undefined') then
     alpha_core = alpha_core_default
     alpha_type = 0
  else
     read  (argument, *) stmp, alpha_core
     select case(stmp)
     case('relative')
        alpha_type = 0
     case('PsiN')
        alpha_type = 1
     case default
        write (6, 9001) trim(command);   stop
     end select
  endif
 9001 format("error: keyword 'PsiN' or 'relative' expected in argument in ",a)


  do iz=0,NZONET-1
  if (Zone(iz)%isfr(1) == SF_CORE) then
     select case(command)
     case(CORE_EXTRAPOLATE)
        write (6, 1001) iz, alpha_core
        call setup_core_domain_extrapolate(iz, nr_EIRENE_core, alpha_core)

     case(CORE_FLUX_SURFACES)
        write (6, 1002) iz, alpha_core
        call setup_core_domain_from_flux_surfaces(iz, nr_EIRENE_core, alpha_core, alpha_type, FULL_FIELD)

     case(CORE_EQ_FLUX_SURFACES)
        write (6, 1002) iz, alpha_core
        call setup_core_domain_from_flux_surfaces(iz, nr_EIRENE_core, alpha_core, alpha_type, EQUILIBRIUM_FIELD)

     case default
        write (6, 9000) trim(core_domain)
        stop
     end select
  endif
  enddo
  write (6, *)

 1000 format(3x,' - Setting up core domain for EIRENE')
 1001 format(8x,'Zone ',i0,': extrapolating from inner EMC3 boundary to r = ',f0.3)
 1002 format(8x,'Zone ',i0,': generating flux surfaces to r = ',f0.3)
 9000 format('error: invalid method "',a,'" for core domain!')
  end subroutine setup_core_domain
  !=====================================================================



  !=====================================================================
  subroutine setup_core_domain_extrapolate(iz, nr_core, alpha_core)
  use iso_fortran_env
  use emc3_grid
  implicit none

  integer,      intent(in) :: iz, nr_core
  real(real64), intent(in) :: alpha_core


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
  subroutine setup_core_domain_from_flux_surfaces(iz, nr_core, alpha_core, alpha_type, bfield)
  use iso_fortran_env
  use emc3_grid
  use fieldline_grid
  use magnetic_axis
  use flux_surface_3D
  use run_control, only: N_points, Trace_Method
  use mesh_spacing
  use grid
  use curve2D
  use equilibrium, only: get_PsiN, get_cylindrical_coordinates
  implicit none

  integer,      intent(in) :: iz, nr_core, bfield, alpha_type
  real(real64), intent(in) :: alpha_core

  type(t_flux_surface_3D) :: F
  type(t_spacing)         :: Sr, Sp
  type(t_curve)           :: C
  type(t_grid)            :: B, G
  character(len=80)       :: filename
  real(real64) :: x0(3), xmag(3), x(3), y(3), r, w, PsiN1, theta
  integer      :: ierr, ig, ir, ip, ipc, it, updown_symmetry


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
  theta = atan2(x0(2)-xmag(2), x0(1)-xmag(1))
  call Sr%init(radial_spacing(-1))
  ! set up poloidal spacing
  call C%new(ZON_POLO(iz))
  C%x(0,:) = 0.d0
  PsiN1    = 0
  do ip=1,SRF_POLO(iz)-1
     ig = ig + SRF_RADI(iz)
     w  = sqrt((RG(ig)-RG(ig-SRF_RADI(iz)))**2 + (ZG(ig)-ZG(ig-SRF_RADI(iz)))**2)
     C%x(ip,1) = 1.d0 * ip / (SRF_POLO(iz)-1)
     C%x(ip,2) = C%x(ip-1,2) + w
     if (alpha_type == 1) PsiN1 = PsiN1 + get_PsiN((/ RG(ig), ZG(ig)/))
  enddo
  C%x(:,2) = C%x(:,2) / C%x(SRF_POLO(iz)-1,2)
  PsiN1    = PsiN1 / (SRF_POLO(iz)-1)
  call Sp%init_manual_function(C)


  ! initialize base grid for core region
  call B%new(CYLINDRICAL, MESH_2D, 3, nr_core, SRF_POLO(iz), fixed_coord_value=x0(3))


  ! up/down symmetric flux surfaces ?
  updown_symmetry = -1
  if (Zone(iz)%it_base == 0  .and.  Zone(iz)%isft(1) == SF_UPDOWN) then
     updown_symmetry = 0
     write (6, *) 'applying up/down symmetry...'
  endif


  ! generate flux surfaces for core region -> setup base grid
  do ir=0,nr_core-1
     r = Sr%node(ir,nr_core)
     select case(alpha_type)
     case(0)
        r = alpha_core + (1.d0 - alpha_core) * r
        x = xmag + r * (x0-xmag)
     case(1)
        y(1) = theta*180.d0/pi;   y(3) = x0(3)*180.d0/pi
        y(2) = alpha_core + (PsiN1 - alpha_core) * r
        x    = get_cylindrical_coordinates(y, ierr)
        if (ierr > 0) then
           write (6, *) 'error: cannot find coordinates for flux surface PsiN = ', y(2)
           stop
        endif
     end select

     call F%generate(x, N_points, symmetry, 1, 360, Trace_Method, poloidal_coordinate=ipc, &
        resample=SRF_POLO(iz), spacings=Sp, updown_symmetry=updown_symmetry, bfield=bfield)
     do ip=0,SRF_POLO(iz)-1
        B%mesh(ir,ip,:) = F%slice(0)%x(ip,:)
     enddo
  enddo
  write (filename, 2000) iz
 2000 format('core_flux_surfaces_',i0,'.plt')
  call B%store(filename)


  ! trace nodes from base grid and set up 3D grid
  call trace_nodes_1zone(iz, B, G, bfield)
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



  !=====================================================================
  subroutine user_defined_core_domain(arguments)
  use emc3_grid, only: NZONET, SRF_RADI, SRF_POLO, SRF_TORO, GRID_P_OS, RG, ZG
  use fieldline, only: FULL_FIELD
  use fieldline_grid
  use grid
  implicit none
  character(len=*), intent(in) :: arguments

  integer, parameter :: iu = 99

  character(len=128) :: filename(0:NZONET-1)
  type(t_grid)       :: B(0:NZONET-1), G
  real(real64)       :: phi, x(2)
  integer            :: iz, ir, ip, it, ig, ios, nr_core


  ! read filenames for each zone from argument list
  read  (arguments, *, iostat=ios) filename
  if (ios /= 0) then
     write (6, 9001) NZONET, trim(arguments);  stop
  endif
 9001 format("error: cannot read filenames for ",i0," zones from '",a,"'")


  ! zone loop
  nr_core = nr_EIRENE_core
  do iz=0,NZONET-1
     ! initialze base mesh for core region
     write (6, 1001) iz
     phi = Zone(iz)%phi(Zone(iz)%it_base) / 180.d0 * pi
     call B(iz)%new(CYLINDRICAL, MESH_2D, 3, nr_core, SRF_POLO(iz), fixed_coord_value=phi)


     ! read base mesh from user defined location
     open  (iu, file=filename(iz))
     do ir=0,nr_core-1
     do ip=0,SRF_POLO(iz)-1
        read  (iu, *) x
        B(iz)%mesh(ir,ip,:) = x
     enddo
     enddo
     close (iu)


     ! trace nodes from base grid and set up 3D grid
     call trace_nodes_1zone(iz, B(iz), G, FULL_FIELD)
     do it=0,SRF_TORO(iz)-1
     do ip=0,SRF_POLO(iz)-1
     do ir=0,nr_core-1
        ig = ir + (ip + it * SRF_POLO(iz)) * SRF_RADI(iz)   +  GRID_P_OS(iz)
        RG(ig) = G%mesh3D(ir,ip,it,1)
        ZG(ig) = G%mesh3D(ir,ip,it,2)
     enddo
     enddo
     enddo
  enddo
 1001 format(8x,'Zone ',i0,' ...')

  end subroutine user_defined_core_domain
  !=====================================================================
