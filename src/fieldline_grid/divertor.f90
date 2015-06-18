module divertor
  use iso_fortran_env
  use emc3_grid
  implicit none

  contains
!=======================================================================


  !=====================================================================
  ! Calculate intersection between divertor leg (of flux surface) and
  ! target surface.
  ! Return relative coordinate (eta) of intersection on divertor leg.
  !=====================================================================
  subroutine divertor_leg_interface(C_leg, C_cut, eta)
  use curve2D
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



!=======================================================================
! CLOSE_GRID_DOMAIN
!
! Extend the divertor legs to first/last toroidal grid slice in order to
! close the simulation domain (required for module MAPPING in EMC3).
!
! !!! These cells MUST be behind the wall !!!
!=======================================================================
  subroutine close_grid_domain(iz)
  use equilibrium

  integer, intent(in) :: iz

  integer :: j0(-1:1), i, j, k, ig, ig0


  j0(-1) = 0
  j0( 1) = SRF_POLO(iz)-1

  ! Adjust divertor leg to the first toroidal slice
  j = j0(int(Ip_sign * Bt_sign))
  k = 0
  write (6, *) 'adjusting to surface iz it ip = ', iz, k, j
  do i=0,SRF_RADI(iz)-1
     k = 0
     ig0 = i + (j + k*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
     do k=1,SRF_TORO(iz)-1
        ig = i + (j + k*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
        RG(ig) = RG(ig0)
        ZG(ig) = ZG(ig0)
     enddo
  enddo

  ! Adjust divertor leg to the last toroidal surface
  j = j0(int(-1 * Ip_sign * Bt_sign))
  k = SRF_TORO(iz)-1
  write (6, *) 'adjusting to surface iz it ip = ', iz, k, j
  do i=0,SRF_RADI(iz)-1
     k = SRF_TORO(iz)-1
     ig0 = i + (j + k*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
     do k=0,SRF_TORO(iz)-2
        ig = i + (j + k*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
        RG(ig) = RG(ig0)
        ZG(ig) = ZG(ig0)
     enddo
  enddo

  end subroutine close_grid_domain
!=======================================================================

end module divertor
