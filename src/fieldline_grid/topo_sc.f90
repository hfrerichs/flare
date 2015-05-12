!===============================================================================
! Simply connected grid layout
!===============================================================================
module topo_sc
  use iso_fortran_env
  implicit none
  private

  public :: &
     setup_topo_sc, &
     make_base_grids_sc

  contains
  !=====================================================================


  
  !=====================================================================
  subroutine setup_topo_sc()
  use fieldline_grid
  use emc3_grid, only: NZONET
  implicit none

  integer :: iz, ib


  ! 0. setup number of zones
  NZONET = blocks


  ! 1. setup resolution for each zone
  write (6, 1000)
  do ib=0,blocks-1
     iz = ib
     if (Zone(iz)%nr == -1) Zone(iz)%nr = nr(0)
     if (Zone(iz)%np == -1) Zone(iz)%np = np(0)
     !if (Zone(iz)%nt == -1) Zone(iz)%nt = nt
     Zone(iz)%nt = Block(ib)%nt

     ! setup toroidal discretization
     allocate (Zone(iz)%phi(0:Zone(iz)%nt))
     Zone(iz)%phi     = Block(ib)%phi
     Zone(iz)%it_base = Block(ib)%it_base

     write (6, 1002) iz, Zone(iz)%nr, Zone(iz)%np, Zone(iz)%nt
  enddo
 1000 format(8x,'Grid resolution is (radial x poloidal x toroidal):')
 1002 format(12x,i3,3x,'(',i0,' x ',i0,' x ',i0,')')


  ! 2. setup connectivity between zones
  do iz=0,blocks-1
     Zone(iz)%isfr(1) = SF_CORE
     Zone(iz)%isfr(2) = SF_VACUUM
     Zone(iz)%isfp(1) = SF_PERIODIC
     Zone(iz)%isfp(2) = SF_PERIODIC
     Zone(iz)%isft(1) = SF_MAPPING
     Zone(iz)%isft(2) = SF_MAPPING

     Zone(iz)%r_surf_pl_trans_range(1) = nr_EIRENE_core
     Zone(iz)%r_surf_pl_trans_range(2) = Zone(iz)%nr - nr_EIRENE_vac
     Zone(iz)%p_surf_pl_trans_range(1) = 0
     Zone(iz)%p_surf_pl_trans_range(2) = Zone(iz)%np
  enddo

  end subroutine setup_topo_sc
  !=====================================================================



  !=====================================================================
  subroutine make_base_grids_sc()
  write (6, *) 'make_base_grids_sc to be implmented!'
  stop
  end subroutine make_base_grids_sc
  !=====================================================================

end module topo_sc
