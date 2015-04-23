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
