!===============================================================================
! Simply connected grid layout
!===============================================================================
subroutine setup_topology_sc
  use fieldline_grid, Zone_ => Zone
  use emc3_grid
  implicit none

  ! unit number for grid configuration file
  integer, parameter :: iu = 12

!...............................................................................
! user defined parameters (to be set via configuration file)                   .

  ! default resolution
  integer :: &
     nr     =   16, &                   ! radial resolution
     np     =  180                      ! poloidal resolution

  type(t_zone_input) :: Zone(0:max_zones-1)
!...............................................................................

  namelist /Block_Resolution/ &
      Zone, nr, np, nt

!...............................................................................
! internal variables

  integer :: iz, ib
!...............................................................................


  ! 0. setup number of zones
  NZONET = blocks


  ! 1. read user configuration from input file
  open  (iu, file=config_file)
  read  (iu, Block_Resolution)
  close (iu)


  ! 3. setup resolution for each zone
  write (6, 1000)
  do ib=0,blocks-1
     iz = ib
     if (Zone(iz)%nr == -1) Zone(iz)%nr = nr
     if (Zone(iz)%np == -1) Zone(iz)%np = np
     !if (Zone(iz)%nt == -1) Zone(iz)%nt = nt
     Zone(iz)%nt = Block(ib)%nt

     ! setup toroidal discretization
     allocate (Zone_(iz)%phi(0:Zone(iz)%nt))
     Zone_(iz)%phi     = Block(ib)%phi
     Zone_(iz)%it_base = Block(ib)%it_base

     write (6, 1002) iz, Zone(iz)%nr, Zone(iz)%np, Zone(iz)%nt
  enddo

 1000 format(8x,'Grid resolution is (radial x poloidal x toroidal):')
 1002 format(12x,i3,3x,'(',i0,' x ',i0,' x ',i0,')')


  ! setup connectivity between zones
  do iz=0,blocks-1
     Zone_(iz)%isfr(1) = SF_CORE
     Zone_(iz)%isfr(2) = SF_VACUUM
     Zone_(iz)%isfp(1) = SF_PERIODIC
     Zone_(iz)%isfp(2) = SF_PERIODIC
     Zone_(iz)%isft(1) = SF_MAPPING
     Zone_(iz)%isft(2) = SF_MAPPING

     Zone_(iz)%r_surf_pl_trans_range(1) = nr_EIRENE_core
     Zone_(iz)%r_surf_pl_trans_range(2) = Zone(iz)%nr - nr_EIRENE_vac
     Zone_(iz)%p_surf_pl_trans_range(1) = 0
     Zone_(iz)%p_surf_pl_trans_range(2) = Zone(iz)%np
  enddo

  Zone_%t_zone_input = Zone
end subroutine setup_topology_sc
