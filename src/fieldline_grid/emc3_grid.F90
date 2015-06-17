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
  subroutine load_grid_layout()
  implicit none

  integer, parameter :: iu = 24

  character(len=72)  :: readfi
  integer :: i, iz, ir, itmp, n


  write (6, 1000)
 1000 format(3x,'- Loading EMC3-EIRENE grid ...')
  open  (iu, file='input.geo')
  call scrape(iu, readfi)
  read (readfi, *) NZONET


  allocate (SRF_RADI(0:NZONET-1),SRF_POLO(0:NZONET-1),SRF_TORO(0:NZONET-1), &
            ZON_RADI(0:NZONET-1),ZON_POLO(0:NZONET-1),ZON_TORO(0:NZONET-1))
  do iz=0,NZONET-1
     call scrape(iu, readfi)
     read  (readfi, *) SRF_RADI(iz), SRF_POLO(iz), SRF_TORO(iz)
     write (6 , *)     SRF_RADI(iz), SRF_POLO(iz), SRF_TORO(iz)
     ZON_RADI(iz) = SRF_RADI(iz) - 1
     ZON_POLO(iz) = SRF_POLO(iz) - 1
     ZON_TORO(iz) = SRF_TORO(iz) - 1
  enddo

  allocate (GRID_P_OS(0:NZONET), MESH_P_OS(0:NZONET), PHI_PL_OS(0:NZONET))
      
  PHI_PL_OS(0) = 0
  GRID_P_OS(0) = 0
  MESH_P_OS(0) = 0
  do i=1,NZONET
     PHI_PL_OS(i) = PHI_PL_OS(i-1) + SRF_TORO(i-1)
     GRID_P_OS(i) = GRID_P_OS(i-1) + SRF_RADI(i-1)*SRF_POLO(i-1)*SRF_TORO(i-1)
     MESH_P_OS(i) = GRID_P_OS(i-1) + ZON_RADI(i-1)*ZON_POLO(i-1)*ZON_TORO(i-1)
  enddo

  ! 3. allocate main arrays
  allocate (PHI_PLANE(0:PHI_PL_OS(NZONET)-1), &
                   RG(0:GRID_P_OS(NZONET)-1), &
                   ZG(0:GRID_P_OS(NZONET)-1))




  ! initialize radial plasma range (r_surf_pl_trans_range)
  allocate (R_SURF_PL_TRANS_RANGE(2,0:NZONET-1), P_SURF_PL_TRANS_RANGE(2,0:NZONET-1))
  do iz=0,NZONET-1
     R_SURF_PL_TRANS_RANGE(1,iz) = SRF_RADI(iz)-1
     R_SURF_PL_TRANS_RANGE(2,iz) = 0
     P_SURF_PL_TRANS_RANGE(1,iz) = 0
     P_SURF_PL_TRANS_RANGE(2,iz) = SRF_POLO(iz)-1
  enddo


  ! non default surfaces
  ! 1. radial
  call scrape(iu, readfi)
  read (readfi, *) n
  do i=1,n
     read  (iu, *) ir, iz, itmp
     write (6,  *) ir, iz, itmp
     read  (iu, *) readfi
     R_SURF_PL_TRANS_RANGE(1,iz) = min(R_SURF_PL_TRANS_RANGE(1,iz), ir)
     R_SURF_PL_TRANS_RANGE(2,iz) = max(R_SURF_PL_TRANS_RANGE(2,iz), ir)
  enddo

  ! 2. poloidal
  call scrape(iu, readfi)
  read (readfi, *) n
  do i=1,n
     read (iu, *) readfi
     read (iu, *) readfi
  enddo

  ! 3. toroidal
  call scrape(iu, readfi)
  read (readfi, *) n
  do i=1,n
     read (iu, *) readfi
     read (iu, *) readfi
  enddo


  ! non transparent surfaces
  ! 1. radial
  call scrape(iu, readfi)
  read (readfi, *) n
  do i=1,n
     read  (iu, *) ir, iz, itmp
     write (6,  *) ir, iz, itmp
     read  (iu, *) readfi
     R_SURF_PL_TRANS_RANGE(1,iz) = min(R_SURF_PL_TRANS_RANGE(1,iz), ir)
     R_SURF_PL_TRANS_RANGE(2,iz) = max(R_SURF_PL_TRANS_RANGE(2,iz), ir)
  enddo

  close (iu)

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


  call load_grid_layout()
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
! WRITE_EMC3_GRID (write 3D field aligned grid to file "grid3D.dat")
!=======================================================================
  subroutine write_emc3_grid
  use math
  implicit none

  integer, parameter :: iu = 24

  integer :: iz, nt, i, j, k, l


  ! write data to file
  open  (iu, file='grid3D.dat')
  do iz=0,NZONET-1
     nt = SRF_TORO(iz)
     write (iu, *) SRF_RADI(iz), SRF_POLO(iz), SRF_TORO(iz)
     do k=0,nt-1
        ! write toroidal position of slices in deg
        write (iu, *) PHI_PLANE(k+PHI_PL_OS(iz)) / pi * 180.d0

        i = k*SRF_POLO(iz)*SRF_RADI(iz) + GRID_P_OS(iz)
        j = i + SRF_POLO(iz)*SRF_RADI(iz) - 1
        write (iu, '(6f12.6)') (RG(l), l=i,j)
        write (iu, '(6f12.6)') (ZG(l), l=i,j)
     enddo
  enddo
  close (iu)

  end subroutine write_emc3_grid
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



!=======================================================================
! fix nodes at interface between zone izA and izB
! toroidal cell range:  itA0 -> itA0 + nt (in zone izA)
!                       itB0 -> itB0 + nt (in zone izB)
! poloidal cell range:  ipA0 -> ipA0 + np (in zone izA)
!                       ipB0 -> ipB0 + np (in zone izB)
! radial surface index: irA (in zone izA)
!                       irB (in zone izB)
! OUTPUT: dmax = updated maximum displacement between nodes
!=======================================================================
  subroutine fix_interface(izA, izB, itA0, itB0, nt, ipA0, ipB0, np, irA, irB, dmax)
  use iso_fortran_env
  integer,      intent(in)    :: izA, izB, itA0, itB0, nt, ipA0, ipB0, np, irA, irB
  real(real64), intent(inout) :: dmax

  real(real64) :: d, R, Z
  integer :: it, itA, itB, ip, ipA, ipB, igA, igB


  do ip=0,np
     ipA = ipA0 + ip
     ipB = ipB0 + ip
     do it=0,nt
        itA = itA0 + it
        itB = itB0 + it

        igA = irA + (ipA + itA*SRF_POLO(izA))*SRF_RADI(izA)  +  GRID_P_OS(izA)
        igB = irB + (ipB + itB*SRF_POLO(izB))*SRF_RADI(izB)  +  GRID_P_OS(izB)


        ! distance between nodes, save max. distance
        d = sqrt((RG(igA)-RG(igB))**2 - (ZG(igA)-ZG(igB))**2)
        if (d.gt.dmax) dmax = d

        ! replace nodes by mean value
        R = 0.5d0 * (RG(igA) + RG(igB))
        RG(igA) = R; RG(igB) = R
        Z = 0.5d0 * (ZG(igA) + ZG(igB))
        ZG(igA) = Z; ZG(igB) = Z
     enddo
  enddo

  end subroutine fix_interface
!=======================================================================

end module emc3_grid
