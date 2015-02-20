!module layout_sc
!  use field_aligned_grid
!  implicit none
!
!  !private
!
!  !type(t_layout) :: L
!
!end module layout_sc


subroutine setup_layout_sc (tZ)
  use field_aligned_grid
  use emc3_grid
  implicit none

  type(t_zone), dimension(0:N_block-1), intent(inout) :: tZ

  integer :: iz


  do iz=0,N_block-1
     tZ(iz)%irsfa = CORE_BOUNDARY_SF
     tZ(iz)%irsfb = VACUUM_BOUNDARY_SF
     tZ(iz)%ipsfa = PERIODIC_SF
     tZ(iz)%ipsfb = PERIODIC_SF
     tZ(iz)%itsfa = MAPPING_SF
     tZ(iz)%itsfb = MAPPING_SF

     R_SURF_PL_TRANS_RANGE(1,iz) = nr_EIRENE_core
     R_SURF_PL_TRANS_RANGE(2,iz) = SRF_RADI(iz)-1 - nr_EIRENE_SOL
     P_SURF_PL_TRANS_RANGE(1,iz) = 0
     P_SURF_PL_TRANS_RANGE(2,iz) = SRF_POLO(iz)-1
  enddo

end subroutine setup_layout_sc
