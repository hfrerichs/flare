!===============================================================================
! Interface to variables in EMC3 for fieldline grid
!===============================================================================
module emc3_grid
  implicit none

  integer, save :: NZONET

  integer, dimension(:), allocatable, save :: &
     SRF_RADI,            SRF_POLO,              SRF_TORO, &
     ZON_RADI,            ZON_POLO,              ZON_TORO, &
     NRS_OFF,             NPS_OFF,               NTS_OFF, &
     GRID_P_OS,           MESH_P_OS,             PHI_PL_OS

  integer, dimension(:,:), allocatable, save :: &
     R_SURF_PL_TRANS_RANGE, &
     P_SURF_PL_TRANS_RANGE

  real*8, dimension(:), allocatable, save :: &
     RG,                  ZG,                PHI_PLANE, &
     BFSTREN,        PSI_POLOIDAL


  contains
!=======================================================================



!=======================================================================
  subroutine scrape(nu,readfi)
  implicit none
  integer, intent(in) :: nu
  character(len=72)   :: readfi


  do
     read  (nu, '(a72)') readfi
     if (readfi(1:1) == '*') then
        if (readfi(2:3) == '**') then
           write (6, *) readfi
        endif
     else
        exit
     endif
  enddo

  end subroutine scrape
!===============================================================================



!===============================================================================
! LOAD_GRID_LAYOUT (from input file "input.geo")
!===============================================================================
  subroutine load_grid_layout
  implicit none

  integer, parameter :: iu = 24

  character(len=72)  :: readfi
  integer :: i, iz, nz


  open  (iu, file='input.geo')
  call scrape(iu, readfi)
  read (readfi, *) NZONET


  allocate (SRF_RADI(0:NZONET-1),SRF_POLO(0:NZONET-1),SRF_TORO(0:NZONET-1), &
            ZON_RADI(0:NZONET-1),ZON_POLO(0:NZONET-1),ZON_TORO(0:NZONET-1))
  do iz=0,NZONET-1
     call scrape(iu, readfi)
     read  (readfi, *) SRF_RADI(iz), SRF_POLO(iz), SRF_TORO(iz)
     write (6 , *)     SRF_RADI(iz), SRF_POLO(iz), SRF_TORO(iz)
  enddo
  ZON_RADI(iz) = SRF_RADI(iz) - 1
  ZON_POLO(iz) = SRF_POLO(iz) - 1
  ZON_TORO(iz) = SRF_TORO(iz) - 1
  close (iu)

  allocate (GRID_P_OS(0:NZONET), PHI_PL_OS(0:NZONET))
      
  PHI_PL_OS(0) = 0
  GRID_P_OS(0) = 0
  do i=1,NZONET
     PHI_PL_OS(i) = PHI_PL_OS(i-1) + SRF_TORO(i-1)
     GRID_P_OS(i) = GRID_P_OS(i-1) + SRF_RADI(i-1)*SRF_POLO(i-1)*SRF_TORO(i-1)
  enddo

  ! 3. allocate main arrays
  allocate (PHI_PLANE(0:PHI_PL_OS(NZONET)-1), &
                   RG(0:GRID_P_OS(NZONET)-1), &
                   ZG(0:GRID_P_OS(NZONET)-1))

  end subroutine load_grid_layout
!===============================================================================



!===============================================================================
! LOAD_EMC3_GRID (from file grid3D.dat)
!===============================================================================
  subroutine load_emc3_grid
  use math
  implicit none

  integer, parameter :: iu = 24

  integer :: iz, i, i1, i2, k, nr, np, nt


  open  (iu, file='grid3D.dat')
  do iz=0,NZONET-1
     read  (iu, *) nr, np, nt
     do k=0,nt-1
        read (iu, *) PHI_PLANE(k+PHI_PL_OS(iz))
        ! convert deg to rad
        PHI_PLANE(k+PHI_PL_OS(iz)) = PHI_PLANE(k+PHI_PL_OS(iz)) * pi / 180.d0

        i1 = GRID_P_OS(iz) + k*SRF_RADI(iz)*SRF_POLO(iz)
        i2 = i1 + SRF_RADI(iz)*SRF_POLO(iz) - 1
        read (iu, *) (RG(i), i=i1,i2)
        read (iu, *) (ZG(i), i=i1,i2)
     enddo
  enddo
  close (iu)

  end subroutine load_emc3_grid
!=======================================================================



!=======================================================================
#ifdef FLARE
  function export_flux_tube(iz, ir, ip) result(flux_tube)
  use flux_tube
  use math
  integer, intent(in) :: iz, ir, ip
  type(t_flux_tube)   :: flux_tube

  integer :: i, it, ig(4)


  call flux_tube%new(ZON_TORO(iz))
  flux_tube%phi = PHI_PLANE(PHI_PL_OS(iz):PHI_PL_OS(iz+1)-1) / 180.d0 * pi
  do it=0,ZON_TORO(iz)
     ig(1) = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
     ig(2) = ig(1) + SRF_RADI(iz)
     ig(3) = ig(2) + 1
     ig(4) = ig(1) + 1

     do i=1,4
        flux_tube%F(i)%x(it,1) = RG(ig(i))
        flux_tube%F(i)%x(it,2) = ZG(ig(i))

        flux_tube%Bmod(it,i)   = BFSTREN(ig(i))
     enddo
  enddo

  end function export_flux_tube
#endif
!=======================================================================

end module emc3_grid
