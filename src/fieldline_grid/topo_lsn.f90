!===============================================================================
! Lower Single Null configuration: block-structured decomposition with zones for
! high pressure region (HPR), scrape-off layer (SOL) and private flux region (PFR)
!===============================================================================
module topo_lsn
  use iso_fortran_env
  use grid
  use separatrix
  use curve2D
  use fieldline_grid, unused => TOPO_LSN
  implicit none
  private

  integer, parameter :: DEFAULT = 0



! more parameters related to the grid generation process
      integer :: n_int  = 4			! use n_int interpolated surfaces for the transition between perturbed flux surfaces at the inner simulation boundary and unperturbed flux surfaces in the main part.




!...............................................................................
! derived parameters                                                           .
  integer :: np1, np2
!...............................................................................




  ! coordinates of X-point and magnetic axis
  real(real64) :: Px(2), Pmag(2)

  ! discretization method
  integer :: method = DEFAULT


  ! base grid in high pressure region (HPR), scrape-off layer (SOL) and private flux region (PFR)
  type(t_grid), dimension(:), allocatable :: G_HPR, G_SOL, G_PFR ! (0:blocks-1)

  ! magnetic separatrix
  type(t_separatrix) :: S
  type(t_curve)      :: S0

  ! guiding_surface
  type(t_curve)      :: C_guide, C_cutL, C_cutR


  public :: &
     setup_topo_lsn, &
     make_base_grids_lsn

  contains
  !=====================================================================



  !=====================================================================
  subroutine setup_topo_lsn()
  use emc3_grid, only: NZONET
  implicit none

  integer :: iz, iz0, iz1, iz2, ib


  ! 0. setup number of zones for lower single null topology
  NZONET = blocks * 3


  ! 1. setup resolution for each zone
  np1 = npR(1) + np(0) + npL(1)
  np2 = npR(1) +         npL(1)
  write (6, 1000)
  write (6, 1001)
  do ib=0,blocks-1
     ! high pressure region (HPR)
     iz0 = 3*ib
     if (Zone(iz0)%nr == -1) Zone(iz0)%nr = nr(0)
     if (Zone(iz0)%np == -1) Zone(iz0)%np = np(0)
     !if (Zone(iz0)%nt == -1) Zone(iz0)%nt = nt
     Zone(iz0)%nt = Block(ib)%nt

     ! scrape-off layer (SOL)
     iz1 = iz0 + 1
     if (Zone(iz1)%nr == -1) Zone(iz1)%nr = nr(1)
     if (Zone(iz1)%np == -1) Zone(iz1)%np = np1
     !if (Zone(iz1)%nt == -1) Zone(iz1)%nt = nt
     Zone(iz1)%nt = Block(ib)%nt

     ! private flux region (PFR)
     iz2 = iz1 + 1
     if (Zone(iz2)%nr == -1) Zone(iz2)%nr = nr(2)
     if (Zone(iz2)%np == -1) Zone(iz2)%np = np2
     !if (Zone(iz2)%nt == -1) Zone(iz2)%nt = nt
     Zone(iz2)%nt = Block(ib)%nt

     ! setup toroidal discretization
     do iz=iz0,iz2
        allocate (Zone(iz)%phi(0:Zone(iz)%nt))
        Zone(iz)%phi = Block(ib)%phi
        Zone(iz)%it_base = Block(ib)%it_base
     enddo


     !write (6, 1002) ib, nr0, np0, nt, nr1, np1, nt, nr2, np2, nt
     write (6, 1002) ib, Zone(iz0)%nr, Zone(iz0)%np, Zone(iz0)%nt, &
                         Zone(iz1)%nr, Zone(iz1)%np, Zone(iz1)%nt, &
                         Zone(iz2)%nr, Zone(iz2)%np, Zone(iz2)%nt
  enddo
 1000 format(8x,'Grid resolution is (radial x poloidal x toroidal):')
 1001 format(8x,'block #, high pressure region, scrape-off layer, private flux region')
 1002 format(12x,i3,3(3x,'(',i0,' x ',i0,' x ',i0,')'))


  ! 2. setup connectivity between zones
  do ib=0,blocks-1
     ! high pressure region (HPR)
     iz = 3 * ib
     Zone(iz)%isfr(1) = SF_CORE
     Zone(iz)%isfr(2) = SF_MAPPING
     Zone(iz)%isfp(1) = SF_PERIODIC
     Zone(iz)%isfp(2) = SF_PERIODIC
     Zone(iz)%isft(1) = SF_MAPPING
     Zone(iz)%isft(2) = SF_MAPPING
     Zone(iz)%r_surf_pl_trans_range(1) = nr_EIRENE_core
     Zone(iz)%r_surf_pl_trans_range(2) = Zone(iz)%nr
     Zone(iz)%p_surf_pl_trans_range(1) = 0
     Zone(iz)%p_surf_pl_trans_range(2) = Zone(iz)%np

     ! scrape-off layer (SOL)
     iz = 3 * ib + 1
     Zone(iz)%isfr(1) = SF_MAPPING
     Zone(iz)%isfr(2) = SF_VACUUM
     Zone(iz)%isfp(1) = SF_VACUUM
     Zone(iz)%isfp(2) = SF_VACUUM
     Zone(iz)%isft(1) = SF_MAPPING
     Zone(iz)%isft(2) = SF_MAPPING
     Zone(iz)%r_surf_pl_trans_range(1) = 0
     Zone(iz)%r_surf_pl_trans_range(2) = Zone(iz)%nr - nr_EIRENE_vac
     Zone(iz)%p_surf_pl_trans_range(1) = 0
     Zone(iz)%p_surf_pl_trans_range(2) = Zone(iz)%np

     ! private flux region (PFR)
     iz = 3 * ib + 2
     Zone(iz)%isfr(1) = SF_VACUUM
     Zone(iz)%isfr(2) = SF_MAPPING
     Zone(iz)%isfp(1) = SF_VACUUM
     Zone(iz)%isfp(2) = SF_VACUUM
     Zone(iz)%isft(1) = SF_MAPPING
     Zone(iz)%isft(2) = SF_MAPPING
     Zone(iz)%r_surf_pl_trans_range(1) = nr_EIRENE_vac
     Zone(iz)%r_surf_pl_trans_range(2) = Zone(iz)%nr
     Zone(iz)%p_surf_pl_trans_range(1) = 0
     Zone(iz)%p_surf_pl_trans_range(2) = Zone(iz)%np
  enddo

  end subroutine setup_topo_lsn
  !=====================================================================



  !=====================================================================
  subroutine setup_domain
  use boundary
  use run_control, only: Debug
  use math
  use flux_surface_2D, only: RIGHT_HANDED
  use equilibrium
  use inner_boundary

  real(real64) :: tmp(3), Rbox(2), Zbox(2), Px0(2), theta0


  ! 1.a setup guiding surface for divertor legs (C_guide) ------------------
  if (guiding_surface .ne. '') then
     write (6, 1000)
     call C_guide%load(guiding_surface)
  else if (n_axi > 1) then
     write (6, 1001)
     call C_guide%copy(S_axi(1))
  else
     write (6, *) 'error: cannot determine divertor geometry!'
     write (6, *) 'neither guiding_surface is set, nor is an axisymmetric surface defined.'
     stop
  endif
 1000 format(8x,'User defined guiding surface for divertor strike points')
 1001 format(8x,'First axisymmetric surface used for divertor strike points')

  ! 1.b setup extended guiding surfaces for divertor leg discretization ----
  ! C_cutL, C_cutR
  call C_cutL%copy(C_guide)
  call C_cutL%left_hand_shift(d_cutL(1))
  call C_cutR%copy(C_guide)
  call C_cutR%left_hand_shift(d_cutR(1))
  if (Debug) then
     call C_cutL%plot(filename='C_cutL.plt')
     call C_cutR%plot(filename='C_cutR.plt')
  endif


  ! 2.a setup magnetic axis (Pmag) --------------------------------------
  tmp = get_magnetic_axis(0.d0); Pmag = tmp(1:2)

  ! 2.b setup X-point (Px, theta0) --------------------------------------
  call get_domain(Rbox, Zbox)
  Px0(1) = Rbox(1) + 1.d0/3.d0 * (Rbox(2)-Rbox(1))
  Px0(2) = Zbox(1) + 1.d0/6.d0 * (Zbox(2)-Zbox(1))
  Px     = find_X(Px0)
  write (6, 2000) Px

  theta0 = atan2(Px(2) - Pmag(2), Px(1) - Pmag(1))
  write (6, 2001) theta0/pi*180.d0

  ! 2.c separatrix (S, S0) ---------------------------------------------
  call S%generate(Px, RIGHT_HANDED, pi/2.d0, C_cutL, C_cutR)
  call S%plot('S', parts=.true.)

  ! connect core segments of separatrix
  S0 = connect(S%M1%t_curve, S%M2%t_curve)
  call S0%plot(filename='S0.plt')
  call S0%setup_angular_sampling(Pmag)
  call S%M3%setup_length_sampling()
  call S%M4%setup_length_sampling()
 2000 format(8x,'found magnetic X-point at: ',2f10.4)
 2001 format(11x,'-> poloidal angle [deg]: ',f10.4)


  ! 3. inner boundaries for EMC3 grid
  call load_inner_boundaries(Pmag, theta0)

  end subroutine setup_domain
  !=====================================================================



  !=====================================================================
  subroutine divertor_leg_interface(C_leg, C_cut, eta)
  type(t_curve), intent(in) :: C_leg, C_cut
  real(real64), intent(out) :: eta

  real(real64) :: x(2)


  if (.not.C_leg%intersect_curve(C_cut, x, eta)) then
     write (6, *) 'error: could not find intersection between divertor leg and guiding surface!'
     call C_leg%plot(filename='divertor_leg.plt')
     stop
  endif

  end subroutine divertor_leg_interface
  !=====================================================================



  !=====================================================================
  subroutine divide_SOL(F, eta, CL, C0, CR)
  use flux_surface_2D
  use math
  type(t_flux_surface_2D), intent(in)  :: F
  real(real64),            intent(in)  :: eta
  type(t_curve),           intent(out) :: CL, CR
  type(t_flux_surface_2D), intent(out) :: C0

  real(real64) :: l, alpha, xiR, xiL


  l = F%length()

  alpha = 1.d0 + eta * (alphaR(1) - 1.d0)
  xiR   = alpha * S%M3%l / l
  alpha = 1.d0 + eta * (alphaL(1) - 1.d0)
  xiL   = 1.d0 - alpha * S%M4%l / l
  call F%split3(xiR, xiL, CR, C0%t_curve, CL)
  call CR%setup_length_sampling()
  !call C0%setup_angular_sampling(Pmag)
  call C0%setup_sampling(Px, Px, Pmag, eta, eta, pi2)
  call CL%setup_length_sampling()

  end subroutine divide_SOL
  !=====================================================================



  !=====================================================================
  subroutine make_base_grids_lsn
  use math
  use inner_boundary
  use flux_surface_2D
  use mesh_spacing

  type(t_flux_surface_2D) :: FS, FSL, FSR, C0
  type(t_curve)           :: CL, CR
  type(t_spacing)         :: Sl, Sr

  real(real64), dimension(:,:,:), pointer :: M_HPR, M_SOL, M_PFR

  character(len=72)   :: filename
  real(real64) :: xi, eta, phi, x(2), x0(2), x1(2), x2(2), d_HPR(2), dx(2)
  integer :: i, j, iz, iz1, iz2, nr0, nr1, nr2, np0, np1l, np1r

  real(real64) :: xiL, xiR
  integer :: iblock


  !.....................................................................
  ! 0. initialize geometry
  call setup_domain()
  !.....................................................................


  !.....................................................................
  ! 1. check input
  if (n_int < 0) then
     write (6, *) 'error: n_int must not be negative!'; stop
  endif
  if (n_int > nr(0)-2) then
     write (6, *) 'error: n_int > nr0 - 2!'; stop
  endif
  !.....................................................................


  !.....................................................................
  ! 2. setup working arrays for base grid
  allocate (G_HPR(0:blocks-1), G_SOL(0:blocks-1), G_PFR(0:blocks-1))
  !.....................................................................



  do iblock=0,blocks-1
     phi = Block(iblock)%phi_base

     ! 0. initialize grids
     iz = iblock*3
     iz1 = iz + 1
     iz2 = iz + 2

     !nr0 = Zone(iz)%nr;  np0 = Zone(iz)%np
     !nr1 = Zone(iz1)%nr; np1 = Zone(iz1)%np
     !nr2 = Zone(iz2)%nr; np2 = Zone(iz2)%np
     nr0 = nr(0); np0 = np(0)
     nr1 = nr(1)
     nr2 = nr(2)
     np1l = npL(1)
     np1r = npR(1)

     ! cell spacings
     call Zone(iz1)%Sr%init(radial_spacing(1))


     call G_HPR(iblock)%new(CYLINDRICAL, MESH_2D, 3, nr0+1, np0+1, fixed_coord_value=phi)
     call G_SOL(iblock)%new(CYLINDRICAL, MESH_2D, 3, nr1+1, np1+1, fixed_coord_value=phi)
     call G_PFR(iblock)%new(CYLINDRICAL, MESH_2D, 3, nr2+1, np2+1, fixed_coord_value=phi)
     M_HPR => G_HPR(iblock)%mesh
     M_SOL => G_SOL(iblock)%mesh
     M_PFR => G_PFR(iblock)%mesh

     !call make_HPR_grid()

     call make_separatrix()
     call make_flux_surfaces_HPR()
     call make_interpolated_surfaces()
     call make_flux_surfaces_SOL()
     call make_flux_surfaces_PFR()

     ! output
     write (filename, 9000) iz
     call G_HPR(iblock)%store(filename=filename)
     write (filename, 9001) iz
     call G_HPR(iblock)%plot_mesh(filename)

     write (filename, 9000) iz1
     call G_SOL(iblock)%store(filename=filename)
     write (filename, 9001) iz1
     call G_SOL(iblock)%plot_mesh(filename)

     write (filename, 9000) iz2
     call G_PFR(iblock)%store(filename=filename)
     write (filename, 9001) iz2
     call G_PFR(iblock)%plot_mesh(filename)
  enddo

 9000 format ('base_grid_',i0,'.dat')
 9001 format ('base_grid_',i0,'.plt')
  contains
  !=====================================================================



  !=====================================================================
  !subroutine make_HPR_grid (iblock)
  !subroutine make_HPR_grid



  !.....................................................................
  subroutine make_separatrix()

  ! 1. discretization of main part of separatrix
  do j=0,np0
     xi = Zone(iz)%Sp%node(j,np0)

     call S0%sample_at(xi, x)
     M_HPR(nr0,      j, :) = x
     M_SOL(  0, np1r+j, :) = x
  enddo

  ! 2. discretization of right separatrix leg
  call divertor_leg_interface(S%M3%t_curve, C_guide, xiR)
  call Sr%init_spline_X1(etaR(1), 1.d0-xiR)
  do j=0,np1r
     xi = 1.d0 - Sr%node(np1r-j,np1r)
     call S%M3%sample_at(xi, x)
     M_SOL(  0,j,:) = x
     M_PFR(nr2,j,:) = x
  enddo

  ! 3. discretization of left separatrix leg
  call divertor_leg_interface(S%M4%t_curve, C_guide, xiL)
  call Sl%init_spline_X1(etaL(1), xiL)
  do j=1,np1l
     xi = Sl%node(j,np1l)
     call S%M4%sample_at(xi, x)
     M_SOL(  0,np1r + np0 + j,:) = x
     M_PFR(nr2,np1r       + j,:) = x
  enddo
  end subroutine make_separatrix
  !.....................................................................


  !.....................................................................
  subroutine make_flux_surfaces_HPR()
  ! unperturbed FLUX SURFACES (high pressure region)

  ! 1. get radial width at poloidal angle of X-point
!  d_HPR = 0.d0
!  do i=0,blocks-1
!     call C_in(i,1)%sample_at(0.d0, x)
!     d_HPR = d_HPR + x / blocks
!  enddo
!  d_HPR = d_HPR - Px
  d_HPR = get_d_HPR(Px, Pmag)

  ! 2. generate flux surfaces
  if (nr0-1 .ge. 2+n_int) write (6, 1000) nr0-1, 2+n_int
  do i=nr0-1, 2+n_int, -1
     write (6, *) i
     eta = 1.d0 - Zone(iz)%Sr%node(i-1,nr0-1)

     x = Px + eta * d_HPR
     call FS%generate_closed(x, RIGHT_HANDED)
     call FS%setup_angular_sampling(Pmag)

     do j=0,np0
        xi = Zone(iz)%Sp%node(j,np0)
        call FS%sample_at(xi, x)
        M_HPR(i,j,:) = x
     enddo
  enddo

 1000 format (8x,'generating unperturbed flux surfaces: ', i0, ' -> ', i0)
  end subroutine make_flux_surfaces_HPR
  !.....................................................................


  !.....................................................................
  subroutine make_interpolated_surfaces()
  ! inner boundaries and interpolated surfaces (2 -> 1+n_int) (high pressure region)

  write (6, 1001) 2, 1+n_int
  do j=0,np0
     xi = Zone(iz)%Sp%node(j,np0)
     ! innermost surfaces
     do i=0,1
        call C_in(iblock,i)%sample_at(xi, x)
        M_HPR(i, j, :) = x
     enddo

     ! interpolated surfaces
     x1 = M_HPR(1,      j,:)
     x2 = M_HPR(2+n_int,j,:)
     do i=2,1+n_int
        eta = Zone(iz)%Sr%node(    i-1,nr0-1) / Zone(iz)%Sr%node(1+n_int,nr0-1)

        M_HPR(i,j,:) = x1 + eta * (x2-x1)
     enddo
  enddo

 1001 format (8x,'interpolating from inner boundary to 1st unperturbed flux surface:, ', &
              i0, ' -> ', i0)
  end subroutine make_interpolated_surfaces
  !.....................................................................


  !.....................................................................
  subroutine make_flux_surfaces_SOL
  ! scrape-off layer

  dx(1) = Px(2) - Pmag(2)
  dx(2) = Pmag(1) - Px(1)
  dx    = dx / sqrt(sum(dx**2)) * d_SOL(1)
  write (6, 1002) nr1
  do i=1,nr1
     write (6, *) i
     eta = Zone(iz1)%Sr%node(i,nr1)
     x0  = Px + eta * dx
     call FS%generate_open(x0, C_cutL, C_cutR)
     call divide_SOL(FS, eta, CL, C0, CR)

     ! right divertor leg
     call divertor_leg_interface(CR, C_guide, xiR)
     call Sr%init_spline_X1(etaR(1), 1.d0-xiR)
     do j=0,np1r
        xi = 1.d0 - Sr%node(np1r-j,np1r)
        call CR%sample_at(xi, x)
        M_SOL(i,j,:) = x
     enddo

     ! main SOL
     do j=0,np0
        xi = Zone(iz)%Sp%node(j,np0)
        call C0%sample_at(xi, x)
        M_SOL(i,np1r+j,:) = x
     enddo

     ! left divertor leg
     call divertor_leg_interface(CL, C_guide, xiL)
     call Sl%init_spline_X1(etaL(1), xiL)
     do j=1,np1l
        xi = Sl%node(j,np1l)
        call CL%sample_at(xi, x)
        M_SOL(i,np1r+np0+j,:) = x
     enddo
  enddo

 1002 format (8x,'generating scrape-off layer: 1 -> ', i0)
  end subroutine make_flux_surfaces_SOL
  !.....................................................................


  !.....................................................................
  subroutine make_flux_surfaces_PFR
  ! private flux region

  dx = Px - Pmag
  dx = dx / sqrt(sum(dx**2)) * d_PFR(1)
  write (6, 1003) nr2-1
  do i=0,nr2-1
     write (6, *) i
     eta = Zone(iz2)%Sr%node(i,nr2)

     x0 = Px + (1.d0-eta) * dx
     ! right divertor leg
     call FSR%generate(x0, -1, AltSurf=C_cutR, sampling=DISTANCE)
     call divertor_leg_interface(FSR%t_curve, C_guide, xiR)
     call Sr%init_spline_X1(etaR(1), 1.d0-xiR)
     do j=0,np1r
        xi = 1.d0 - Sr%node(np1r-j,np1r)
        call FSR%sample_at(xi, x)
        M_PFR(i,j,:) = x
     enddo

     ! left divertor leg
     call FSL%generate(x0,  1, AltSurf=C_cutL, sampling=DISTANCE)
     call divertor_leg_interface(FSL%t_curve, C_guide, xiL)
     call Sl%init_spline_X1(etaL(1), xiL)
     do j=1,np1l
        xi = Sl%node(j,np1l)
        call FSL%sample_at(xi, x)
        M_PFR(i,np1r + j,:) = x
     enddo
  enddo

 1003 format (8x,'generating private flux region: 0 -> ', i0)
  end subroutine make_flux_surfaces_PFR
  !.....................................................................

  end subroutine make_base_grids_lsn
  !=============================================================================

end module topo_lsn
