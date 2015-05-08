!===============================================================================
! Lower Single Null configuration: block-structured decomposition with zones for
! high pressure region (HPR), scrape-off layer (SOL) and private flux region (PFR)
!===============================================================================
module topo_lsn
  use iso_fortran_env
  use grid
  use separatrix
  use curve2D
  use fieldline_grid, only: blocks, Block, Zone
  implicit none

  private

  integer, parameter :: DEFAULT = 0

!...............................................................................
! user defined parameters (to be set via configuration file)                   .

! np: poloidal resolution ...
  integer :: &
      np0     =  60, &                  ! ... in high pressure region
      np1l    =  70, &                  ! ... in left divertor leg
      np1r    =  70                     ! ... in right divertor leg


! nr: radial resolution ...
  integer :: &
      nr0     =  10, &                  ! ... in high pressure region
      nr1     =  20, &                  ! ... in SOL region
      nr2     =  20                     ! ... in private flux region

!...............................................................................


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

  public :: &
     make_base_grids

  contains
  !=====================================================================



  !=====================================================================
  subroutine make_base_grids
  use math
  use equilibrium
  use inner_boundary


  real(real64) :: Rbox(2), Zbox(2), Px0(2), tmp(3), theta0
  integer :: iblock


  ! set derived parameters
  np1 = np1r + np0 + np1l
  np2 = np1r +       np1l


!.......................................................................
! 0. initialize geometry
!.......................................................................
  ! check input
  if (n_int < 0) then
     write (6, *) 'error: n_int must not be negative!'; stop
  endif
  if (n_int > nr0-2) then
     write (6, *) 'error: n_int > nr0 - 2!'; stop
  endif
  
!- setup magnetic axis -------------------------------------------------
  tmp = get_magnetic_axis(0.d0); Pmag = tmp(1:2)


!- setup X-point -------------------------------------------------------
  call get_domain(Rbox, Zbox)
  Px0(1) = Rbox(1) + 1.d0/3.d0 * (Rbox(2)-Rbox(1))
  Px0(2) = Zbox(1) + 1.d0/6.d0 * (Zbox(2)-Zbox(1))
  Px     = find_X(Px0)
  write (6, *) 'found magnetic X-point at: ', Px

  theta0 = atan2(Px(2) - Pmag(2), Px(1) - Pmag(1))
  write (6, *) '   -> poloidal angle: ', theta0
!      eXm = Pmag - Px
!      eXm = eXm / dsqrt(sum(eXm**2))
!-----------------------------------------------------------------------

  call S%generate(Px, 1, pi/2.d0)
  !call S%plot('S', parts=.true.)

  ! connect core segments of separatrix
  !call S%M1%flip()
  !call S%M2%flip()
  S0 = connect(S%M1%t_curve, S%M2%t_curve)
  call S0%plot(filename='S0.plt')
  call S0%setup_angular_sampling(Pmag)

  call load_inner_boundaries(Pmag, theta0)
!.......................................................................



!.......................................................................
! setup working arrays for base grid
  allocate (G_HPR(0:blocks-1), G_SOL(0:blocks-1), G_PFR(0:blocks-1))
!.......................................................................



  do iblock=0,blocks-1
     !call make_HPR_grid(iblock)
     call make_HPR_grid()
  enddo


  contains
  !end subroutine make_base_grids
  !=====================================================================



  !=====================================================================
  !subroutine make_HPR_grid (iblock)
  subroutine make_HPR_grid
  use math
  use inner_boundary
  use flux_surface_2D

  !integer, intent(in) :: iblock

  type(t_flux_surface_2D) :: FS

  real(real64), dimension(:,:,:), pointer :: M_HPR, M_SOL, M_PFR

  character(len=72)   :: filename
  real(real64) :: xi, eta, phi, x(2), x1(2), x2(2), d_HPR(2)
  integer :: i, j, iz, iz1, iz2


  phi = Block(iblock)%phi_base

  ! 0. initialize grids
  iz = iblock*3
  iz1 = iz + 1
  iz2 = iz + 2

  nr0 = Zone(iz)%nr;  np0 = Zone(iz)%np
  nr1 = Zone(iz1)%nr; np1 = Zone(iz1)%np
  nr2 = Zone(iz2)%nr; np2 = Zone(iz2)%np


  call G_HPR(iblock)%new(CYLINDRICAL, MESH_2D, 3, nr0+1, np0+1, fixed_coord_value=phi)
  call G_SOL(iblock)%new(CYLINDRICAL, MESH_2D, 3, nr1+1, np1+1, fixed_coord_value=phi)
  call G_PFR(iblock)%new(CYLINDRICAL, MESH_2D, 3, nr2+1, np2+1, fixed_coord_value=phi)
  M_HPR => G_HPR(iblock)%mesh
  M_SOL => G_SOL(iblock)%mesh
  M_PFR => G_PFR(iblock)%mesh


  ! 1. generate discretization on separatrix and innermost boundaries
  do j=0,np0
     xi = Zone(iz)%Sp%node(j,np0)

     call S0%sample_at(xi, x)
     M_HPR(nr0,      j, :) = x
     M_SOL(  0, np1r+j, :) = x

     do i=0,1
        call C_in(iblock,i)%sample_at(xi, x)
        M_HPR(i, j, :) = x
     enddo
  enddo


  ! 2. FLUX SURFACES
  d_HPR = 0.d0
  do i=0,blocks-1
     call C_in(0,1)%sample_at(0.d0, x)
     d_HPR = d_HPR + x / blocks
  enddo
  d_HPR = d_HPR - Px

  if (nr0-1 .ge. 2+n_int) write (6, 1000) nr0-1, 2+n_int
  do i=nr0-1, 2+n_int, -1
     eta = 1.d0 - Zone(iz)%Sr%node(i-1,nr0-1)

     x = Px + eta * d_HPR
     write (6, *) i, eta, x
     call FS%generate(x, 1, theta_cut=theta0)
     FS%x(FS%n_seg,:) = FS%x(0,:)
     call FS%plot(filename='test_fs.plt')
     !stop
     call FS%setup_angular_sampling(Pmag)

     do j=0,np0
        xi = Zone(iz)%Sp%node(j,np0)
        call FS%sample_at(xi, x)
        M_HPR(i,j,:) = x
     enddo
  enddo




  ! 3. interpolate (2 -> 1+n_int)
  do j=0,np0
     !xi = Zone(iz)%Sp%node(j,np0)

     x1 = M_HPR(1,      j,:)
     x2 = M_HPR(2+n_int,j,:)
     do i=2,1+n_int
        eta = (Zone(iz)%Sr%node(    i-1,nr0-1) - Zone(iz)%Sr%node(0,nr0-1)) &
            / (Zone(iz)%Sr%node(1+n_int,nr0-1) - Zone(iz)%Sr%node(0,nr0-1))

        M_HPR(i,j,:) = x1 + eta * (x2-x1)
     enddo
  enddo


  ! 4. output
  write (filename, 9000) iz
  call G_HPR(iblock)%store(filename=filename)
  write (filename, 9001) iz
  call G_HPR(iblock)%plot_mesh(filename)

 1000 format (8x,'generating unperturbed flux surfaces: ', i0, ' -> ', i0)
 9000 format ('base_grid_',i0,'.dat')
 9001 format ('base_grid_',i0,'.plt')
  end subroutine make_HPR_grid
  !=====================================================================
  end subroutine make_base_grids

end module topo_lsn













!===============================================================================
! Lower Single Null configuration: block-structured decomposition with zones for
! confined region (CNF), scrape-off layer (SOL) and private flux region (PFR)
!===============================================================================
subroutine setup_topology_lsn()
  use fieldline_grid, Zone_ => Zone
  use emc3_grid
  implicit none

  ! unit number for grid configuration file
  integer, parameter :: iu = 12

!...............................................................................
! user defined parameters (to be set via configuration file)                   .

  ! radial resolution ...
  integer :: &
     nr0     =   10, &                    ! ... of confined region
     nr1     =   20, &                    ! ... of SOL region
     nr2     =   20                       ! ... of private flux region


  ! poloidal resolution ...
  integer :: &
     np0     =  180, &                    ! ... of confined region
     np1l    =   30, &                    ! ... of left divertor leg
     np1r    =   30                       ! ... of right divertor leg


  real(real64) :: &
     DSOL = 24.0, &            ! width of scrape-off layer (SOL)
     DPFR = 12.0               ! width of private flux reagion (PFR)


  type(t_zone_input) :: Zone(0:max_zones-1)
!...............................................................................

  namelist /Block_Resolution/ &
      Zone, nr0, nr1, nr2, np0, np1l, np1r, nt, &
      DSOL, DPFR

!...............................................................................
! internal variables

  integer :: np1, np2, iz, iz0, iz1, iz2, ib
!...............................................................................


  ! 0. setup number of zones for lower single null topology
  NZONET = blocks * 3


  ! 1. read user configuration from input file
  open  (iu, file=config_file)
  read  (iu, Block_Resolution)
  close (iu)
  np1 = np1r + np0 + np1l
  np2 = np1r + np1l


  ! 3. setup resolution for each zone
  write (6, 1000)
  write (6, 1001)
  do ib=0,blocks-1
     ! confined region
     iz0 = 3*ib
     if (Zone(iz0)%nr == -1) Zone(iz0)%nr = nr0
     if (Zone(iz0)%np == -1) Zone(iz0)%np = np0
     !if (Zone(iz0)%nt == -1) Zone(iz0)%nt = nt
     Zone(iz0)%nt = Block(ib)%nt

     ! scrape-off layer
     iz1 = iz0 + 1
     if (Zone(iz1)%nr == -1) Zone(iz1)%nr = nr1
     if (Zone(iz1)%np == -1) Zone(iz1)%np = np1
     !if (Zone(iz1)%nt == -1) Zone(iz1)%nt = nt
     Zone(iz1)%nt = Block(ib)%nt

     ! private flux region
     iz2 = iz1 + 1
     if (Zone(iz2)%nr == -1) Zone(iz2)%nr = nr2
     if (Zone(iz2)%np == -1) Zone(iz2)%np = np2
     !if (Zone(iz2)%nt == -1) Zone(iz2)%nt = nt
     Zone(iz2)%nt = Block(ib)%nt

     ! setup toroidal discretization
     do iz=iz0,iz2
        allocate (Zone_(iz)%phi(0:Zone(iz)%nt))
        Zone_(iz)%phi = Block(ib)%phi
        Zone_(iz)%it_base = Block(ib)%it_base
     enddo


     !write (6, 1002) ib, nr0, np0, nt, nr1, np1, nt, nr2, np2, nt
     write (6, 1002) ib, Zone(iz0)%nr, Zone(iz0)%np, Zone(iz0)%nt, &
                         Zone(iz1)%nr, Zone(iz1)%np, Zone(iz1)%nt, &
                         Zone(iz2)%nr, Zone(iz2)%np, Zone(iz2)%nt
  enddo

 1000 format(8x,'Grid resolution is (radial x poloidal x toroidal):')
 1001 format(8x,'block #, confined region, scrape-off layer, private flux region')
 1002 format(12x,i3,3(3x,'(',i0,' x ',i0,' x ',i0,')'))


  ! setup connectivity between zones
  do ib=0,blocks-1
     ! upstream region
     iz = 3 * ib
     Zone_(iz)%isfr(1) = SF_CORE
     Zone_(iz)%isfr(2) = SF_MAPPING
     Zone_(iz)%isfp(1) = SF_PERIODIC
     Zone_(iz)%isfp(2) = SF_PERIODIC
     Zone_(iz)%isft(1) = SF_MAPPING
     Zone_(iz)%isft(2) = SF_MAPPING
     Zone_(iz)%r_surf_pl_trans_range(1) = nr_EIRENE_core
     Zone_(iz)%r_surf_pl_trans_range(2) = Zone(iz)%nr
     Zone_(iz)%p_surf_pl_trans_range(1) = 0
     Zone_(iz)%p_surf_pl_trans_range(2) = Zone(iz)%np

     ! scrape-off layer
     iz = 3 * ib + 1
     Zone_(iz)%isfr(1) = SF_MAPPING
     Zone_(iz)%isfr(2) = SF_VACUUM
     Zone_(iz)%isfp(1) = SF_VACUUM
     Zone_(iz)%isfp(2) = SF_VACUUM
     Zone_(iz)%isft(1) = SF_MAPPING
     Zone_(iz)%isft(2) = SF_MAPPING
     Zone_(iz)%r_surf_pl_trans_range(1) = 0
     Zone_(iz)%r_surf_pl_trans_range(2) = Zone(iz)%nr - nr_EIRENE_vac
     Zone_(iz)%p_surf_pl_trans_range(1) = 0
     Zone_(iz)%p_surf_pl_trans_range(2) = Zone(iz)%np

     ! private flux region
     iz = 3 * ib + 2
     Zone_(iz)%isfr(1) = SF_VACUUM
     Zone_(iz)%isfr(2) = SF_MAPPING
     Zone_(iz)%isfp(1) = SF_VACUUM
     Zone_(iz)%isfp(2) = SF_VACUUM
     Zone_(iz)%isft(1) = SF_MAPPING
     Zone_(iz)%isft(2) = SF_MAPPING
     Zone_(iz)%r_surf_pl_trans_range(1) = nr_EIRENE_vac
     Zone_(iz)%r_surf_pl_trans_range(2) = Zone(iz)%nr
     Zone_(iz)%p_surf_pl_trans_range(1) = 0
     Zone_(iz)%p_surf_pl_trans_range(2) = Zone(iz)%np
  enddo

  Zone_%t_zone_input = Zone
end subroutine setup_topology_lsn
