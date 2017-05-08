!===============================================================================
! Simply connected grid layout
!===============================================================================
module modtopo_stel
  use iso_fortran_env
  use fieldline_grid, unused => TOPO_STEL
  use inner_boundary
  implicit none
  private


  ! outer boundary
  type(t_curve) :: B


  public :: &
     setup_topo_stel, &
     make_base_grids_stel, &
     post_process_grid_stel

  contains
  !=====================================================================


  
  !=====================================================================
  subroutine setup_topo_stel()
  use emc3_grid, only: NZONET
  implicit none

  integer :: iz, ib


  ! 0. setup number of zones
  NZONET = blocks


  write (6, 1000)
  do ib=0,blocks-1
     ! 1. set up derived parameters
     iz = ib
     Zone(iz)%nr = Block(ib)%nr(0) + nr_EIRENE_core + nr_EIRENE_vac
     Zone(iz)%np = Block(ib)%np(0)


     ! 2. set up zones
     call Zone(iz)%setup(ib, 0, TYPE_SINGLE_LAYER, SF_PERIODIC)


     ! 3. show zone information
     write (6, 1002) iz, Zone(iz)%nr, Zone(iz)%np, Zone(iz)%nt
  enddo
  Zone(0)%isft(1) = SF_UPDOWN
  write (6, *)

 1000 format(8x,'Grid resolution is (radial x poloidal x toroidal):')
 1002 format(12x,i3,3x,'(',i0,' x ',i0,' x ',i0,')')
  end subroutine setup_topo_stel
  !=====================================================================



  !=====================================================================
  subroutine setup_domain(iblock)
  use run_control, only: Debug, N_points, Trace_Method
  use equilibrium
  use boundary
  use poincare_set
  use divertor, only: Pmag
  use string

  integer, intent(in)  :: iblock

  character(len=256)   :: filename, command, argument
  type(t_poincare_set) :: P
  real(real64)         :: phi0, x1(2), tmp(3), theta0, dl
  integer              :: i, iboundary, n, nsample


  ! 0. initialize magnetic axis
  phi0   = Block(iblock)%phi_base / 180.d0 * pi
  tmp    = get_magnetic_axis(phi0); Pmag = tmp(1:2)
  x1     = x_in2(1:2)
  theta0 = get_poloidal_angle(x_in2)


  ! 1. set up outer simulation boundary
  write (6, 1000)
  n = get_commands(guiding_surface(iblock))
  do i=1,n
     call read_command(guiding_surface(iblock), i, command, argument)
     select case(command)
     case('LOAD')
        write (6, 1001) trim(argument)
        call B%load(argument, output=SILENT)


     case('EXPAND')
        read  (argument, *, err=9000) dl
        write (6, 1002) dl
        ! negative sign for expanding surface with nodes in counter-clockwise direction
        call B%left_hand_shift(-dl)


     case('SORT')
        write (6, 1003)
        call B%sort_loop(Pmag)


     case('FLUX_SURFACE')
        tmp = read_vector(argument, 3)
        write (6, 1004) tmp
        ! set default number of points
        if (N_points == 0) N_points = 1000
        call P%generate(tmp, N_points, symmetry, 1, 3600/symmetry, Trace_Method, .false.)
        call B%new(P%slice(0)%nrow-1)
        B%x = P%slice(0)%x(:,1:2)


     case('RESAMPLE')
        read  (argument, *, err=9010) nsample
        write (6, 1006) nsample
        call B%resample(nsample)


     case('PLOT')
        write (6, 1005) trim(argument)
        call B%plot(filename=argument)


     case('BOUNDARY')
        read  (argument, *, err=9010) iboundary
        B = boundary_slice(iboundary, phi0)
        write (6, 1007)


     case default
        write (6, *) 'error: invalid command ', trim(command), ' for guiding surface!'
        stop
     end select

     ! DEBUGGING OUTPUT
     if (Debug) then
        write (filename, 8000) i
        call B%plot(filename=filename)
     endif
  enddo
  write (6, *)
  if (B%n_seg < 0) then
     write (6, *) 'error: outer boundary undefined!'
     stop
  endif


  ! 2. set up inner simulation boundary and sampling
  if (iblock > 0) return
  select case(poloidal_discretization)
  case(POLOIDAL_ANGLE)
     ! use geometric poloidal angle as reference coordinate
     theta0 = 0.d0
     call load_inner_boundaries(theta0)

     ! setup outer boundary
     call B%setup_angular_sampling(Pmag)

  case(ARC_LENGTH)
     theta0 = 0.d0
     call load_inner_boundaries(theta0, DISTANCE)

     ! setup outer boundary
     call B%setup_length_sampling()

  end select


  return
 1000 format(3x,'- Outer simulation boundary:')
 1001 format(8x,'loading from file "',a,'"')
 1002 format(8x,'expanding surface by dl=',f0.3)
 1003 format(8x,'sorting points with respect to geometric poloidal angle')
 1004 format(8x,'generating flux surace from reference point at (',f0.3,', ',f0.3,', ',f0.3,')')
 1005 format(8x,'writing boundary surface to file "',a,'"')
 1006 format(8x,'resampling surface with ',i0,' points')
 1007 format(8x,'loading shape from boundary ',i0)
 8000 format('DEBUG_OUTER_BOUNDARY_STEP',i0,'.PLT')
 9000 write (6, 9001) trim(argument);  stop
 9001 format('error: cannot obtain floating point value from argument ', a)
 9010 write (6, 9011) trim(argument);  stop
 9011 format('error: cannot obtain integer value from argument ', a)
  end subroutine setup_domain
  !=====================================================================


  !=====================================================================
  subroutine make_base_grids_stel()
  use grid
  use divertor, only: make_flux_surfaces_HPR, make_interpolated_surfaces, &
                      make_ortho_grid

  real(real64), dimension(:,:,:), pointer :: M
  type(t_grid) :: G(0:blocks-1)
  real(real64) :: phi
  integer      :: iblock, iz


  do iblock=0,blocks-1
     write (6, *) iblock
     call setup_domain(iblock)


     ! set zone indix
     iz  = iblock

     ! set local variables for resolution
     call load_local_resolution(iblock)


     ! cell spacings
     call Zone(iz)%Sp%init(poloidal_spacing(0))
     call Zone(iz)%Sr%init(radial_spacing(0))


     ! initialize base grid in present block
     phi = Block(iblock)%phi_base / 180.d0 * pi
     call G(iblock)%new(CYLINDRICAL, MESH_2D, 3, nr(0)+1, np(0)+1, fixed_coord_value=phi)
     M => G(iblock)%mesh


     ! set up discretization
     call setup_inner_boundary(iblock, 1, Zone(iz)%Sp, G(iblock))
     call sample_flux_surface(M, B, nr(0), np(0), Zone(iz)%Sp)
     call interpolate_flux_surfaces(M, nr(0), np(0), 1, nr(0), Zone(iz)%Sr)
     call setup_inner_boundary0(iblock, G(iblock))
     if (iblock == 0) call force_up_down_symmetry(M, nr(0), np(0))


     ! output
     call write_base_grid(G(iblock), iblock)
  enddo

  end subroutine make_base_grids_stel
  !=====================================================================



  !=====================================================================
  subroutine sample_flux_surface(M, fs, ir, np, Sp)
  use mesh_spacing
  real(real64), dimension(:,:,:), pointer, intent(inout) :: M
  type(t_curve),                           intent(in)    :: fs
  integer,                                 intent(in)    :: ir, np
  type(t_spacing),                         intent(in)    :: Sp

  real(real64) :: xi, x(2)
  integer :: j


  do j=0,np
     xi = Sp%node(j,np)
     call fs%sample_at(xi, x)
     M(ir,j,:) = x
  enddo

  end subroutine sample_flux_surface
  !=====================================================================



  !=====================================================================
  subroutine interpolate_flux_surfaces(M, nr, np, ir1, ir2, Sr)
  real(real64), dimension(:,:,:), pointer, intent(inout) :: M
  integer,                                 intent(in)    :: nr, np, ir1, ir2
  type(t_spacing),                         intent(in)    :: Sr

  real(real64) :: x1(2), x2(2), s
  integer      :: i, j


  do j=0,np
     x1 = M(ir1, j, :)
     x2 = M(ir2, j, :)
     do i=ir1+1,ir2-1
        s = Sr%node(i-ir1, ir2-ir1)
        M(i, j, :) = x1 + s * (x2-x1)
     enddo
  enddo

  end subroutine interpolate_flux_surfaces
  !=====================================================================



  !=====================================================================
  subroutine force_up_down_symmetry(M, nr, np)
  real(real64), dimension(:,:,:), pointer, intent(inout) :: M
  integer,                                 intent(in)    :: nr, np

  real(real64) :: x(2)
  integer :: i, j, j2



  do i=0,nr
     do j=0,np/2
        j2 = np-j
        x(1) = 0.5d0 * (M(i, j, 1) + M(i, j2, 1))
        x(2) = 0.5d0 * (M(i, j, 2) - M(i, j2, 2))

        M(i,  j, :) = x
        M(i, j2, 1) = x(1)
        M(i, j2, 2) = -x(2)
     enddo
  enddo

  end subroutine force_up_down_symmetry
  !=====================================================================



  !=============================================================================
  subroutine post_process_grid_stel()

  return
  end subroutine post_process_grid_stel
  !=============================================================================



end module modtopo_stel
