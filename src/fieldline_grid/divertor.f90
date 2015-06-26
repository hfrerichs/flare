module divertor
  use iso_fortran_env
  use emc3_grid
  use fieldline_grid, only: guiding_surface, d_cutL, d_cutR
  use curve2D
  use separatrix
  implicit none

  ! magnetic separatrix
  type(t_separatrix), dimension(:), allocatable :: S

  ! guiding_surface
  type(t_curve), public      :: C_guide, C_cutL, C_cutR

  ! store magnetic axis
  real(real64) :: Pmag(2)

  contains
!=======================================================================


  !=====================================================================
  ! setup geometry of simulation domain
  ! based on nx X-points Xp(1:nx) from module equilibrium
  ! and how they connect to each other (connectX(ix) = number of X-point to which separatrix from X-point ix connects to)
  !=====================================================================
  subroutine setup_geometry(nx, connectX)
  use math
  use equilibrium
  use boundary
  use inner_boundary
  integer, intent(in) :: nx, connectX(nx)

  character(len=2) :: Sstr
  real(real64)     :: tmp(3)
  integer          :: ix


  ! 1.a set up guiding surface for divertor legs (C_guide) ---------------------
  if (guiding_surface .ne. '') then
     write (6, 1000)
     call C_guide%load(guiding_surface)
  else if (n_axi > 0) then
     write (6, 1001)
     call C_guide%copy(S_axi(1))
  else
     write (6, *) 'error: cannot determine divertor geometry!'
     write (6, *) 'neither guiding_surface is set, nor is an axisymmetric surface defined.'
     stop
  endif
 1000 format(8x,'User defined guiding surface for divertor strike points')
 1001 format(8x,'First axisymmetric surface used for divertor strike points')

  ! 1.b set up extended guiding surfaces for divertor leg discretization -------
  ! C_cutL, C_cutR
  call C_cutL%copy(C_guide)
  call C_cutL%left_hand_shift(d_cutL(1))
  call C_cutL%plot(filename='C_cutL.plt')
  call C_cutR%copy(C_guide)
  call C_cutR%left_hand_shift(d_cutR(1))
  call C_cutR%plot(filename='C_cutR.plt')


  ! 2.a set up magnetic axis (Pmag) --------------------------------------
  tmp = get_magnetic_axis(0.d0); Pmag = tmp(1:2)
  write (6, 2000) Pmag
 2000 format(8x,'Magnetic axis at: ',2f10.4)

  ! 2.b check X-points
  do ix=1,nx
     tmp(1:2) = Xp(ix)%load()
     write (6, 2001) ix, tmp(1:2)
     write (6, 2002) Xp(ix)%theta/pi*180.d0
  enddo
 2001 format(8x,i0,'. X-point at: ',2f10.4)
 2002 format(11x,'-> poloidal angle [deg]: ',f10.4)

  ! 2.c generate separatrix(ces)
  allocate(S(nx))
  do ix=1,nx
     call S(ix)%generate(ix, C_cutL=C_cutL, C_cutR=C_cutR, iconnect=connectX(ix))
     write (Sstr, 2003) ix; call S(ix)%plot(trim(Sstr), parts=.true.)

     call S(ix)%M1%setup_length_sampling()
     call S(ix)%M2%setup_length_sampling()
     call S(ix)%M3%setup_length_sampling()
     call S(ix)%M4%setup_length_sampling()
  enddo
 2003 format('S',i0)


  ! 3. inner boundaries for EMC3 grid
  call load_inner_boundaries(Xp(1)%theta)


  end subroutine setup_geometry
  !=====================================================================



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
