module emc3_grid
  implicit none

  integer, save :: NZONET, NDOMAIN, NSYM

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
! LOAD_EMC3_GRID
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
!=======================================================================


#ifdef FLARE
!=======================================================================
  subroutine setup_emc3_grid_layout
  use field_aligned_grid
  use emc3_grid


  integer :: iz


  write (6, *) 'running subroutine setup_emc3_grid_layout'

  NZONET = N_block

  ! 1. setup resolution for each zone
  allocate ( &
     SRF_RADI(0:NZONET-1),SRF_POLO(0:NZONET-1),SRF_TORO(0:NZONET-1), &
     ZON_RADI(0:NZONET-1),ZON_POLO(0:NZONET-1),ZON_TORO(0:NZONET-1) )
  allocate (R_SURF_PL_TRANS_RANGE(2,0:NZONET-1))
  allocate (P_SURF_PL_TRANS_RANGE(2,0:NZONET-1))
  R_SURF_PL_TRANS_RANGE = 0
  P_SURF_PL_TRANS_RANGE = 0

  do iz=0,NZONET-1
     ZON_RADI(IZ) = TD(iz+1)%nr
     ZON_POLO(IZ) = TD(iz+1)%np
     ZON_TORO(IZ) = TD(iz+1)%nt(0)
  enddo
  SRF_RADI = ZON_RADI + 1
  SRF_POLO = ZON_POLO + 1
  SRF_TORO = ZON_TORO + 1


  ! 2. setup offset arrays
  allocate (PHI_PL_OS(0:NZONET), GRID_P_OS(0:NZONET), MESH_P_OS(0:NZONET))
  PHI_PL_OS(0) = 0
  GRID_P_OS(0) = 0
  MESH_P_OS(0) = 0
  do IZ=1,NZONET
    PHI_PL_OS(IZ) = PHI_PL_OS(IZ-1) + SRF_TORO(IZ-1)
    GRID_P_OS(IZ) = GRID_P_OS(IZ-1) + SRF_RADI(IZ-1)*SRF_POLO(IZ-1)*SRF_TORO(IZ-1)
    MESH_P_OS(IZ) = MESH_P_OS(IZ-1) + ZON_RADI(IZ-1)*ZON_POLO(IZ-1)*ZON_TORO(IZ-1)
  enddo


  ! 3. allocate main arrays
  allocate (PHI_PLANE(0:PHI_PL_OS(NZONET)-1), &
                   RG(0:GRID_P_OS(NZONET)-1), &
                   ZG(0:GRID_P_OS(NZONET)-1))

  end subroutine setup_emc3_grid_layout
!=======================================================================




!=======================================================================
! WRITE_EMC3_GRID (write 3D field aligned grid to file "grid.data"
!===============================================================================
  subroutine write_emc3_grid
  use emc3_grid
  use math
  implicit none

  integer, parameter :: iu = 24

  integer :: iz, it, i, j, k, l


  ! write data to file
  open  (iu, file='grid3D.dat')
  do iz=0,NZONET-1
     it = ZON_TORO(iz)
     write (iu, *) SRF_RADI(iz), SRF_POLO(iz), SRF_TORO(iz)
     do k=0,it
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
  subroutine sample_bfield_on_emc3_grid
  use iso_fortran_env
  use emc3_grid
  use bfield
  use math, only: pi
  implicit none

  integer, parameter :: iu = 72

  real(real64) :: x(3), Bf(3)
  integer :: ir, ip, it, iz, ig


  write (6, *) 'sampling magnetic field strength on grid ...'
  if (.not.allocated(BFSTREN)) allocate(BFSTREN(0:GRID_P_OS(NZONET)-1))
  BFSTREN = 0.d0


  do iz=0,NZONET-1
  do it=0,SRF_TORO(iz)-1
     write (6, *) iz, it
     x(3) = PHI_PLANE(it + PHI_PL_OS(iz))
     do ip=0,SRF_POLO(iz)-1
     !do ir=0,SRF_RADI(iz)-1
     do ir=R_SURF_PL_TRANS_RANGE(1,iz),R_SURF_PL_TRANS_RANGE(2,iz)
        ig = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)

        x(1) = RG(ig)
        x(2) = ZG(ig)
        Bf   = get_Bf_cyl(x)

        BFSTREN(ig) = sqrt(sum(Bf**2))
     enddo
     enddo
  enddo
  enddo
  ! convert units: Gauss -> T
  !BFSTREN = BFSTREN / 1.d4 ! TODO: check why it's already in T


  open  (iu, file='bfield.dat')
  write (iu, *) BFSTREN
  close (iu)

  end subroutine sample_bfield_on_emc3_grid
!=======================================================================



!=======================================================================
  subroutine write_emc3_input_files
  use field_aligned_grid, only: TD
  use emc3_grid
  implicit none

  integer, parameter :: iu = 72


  call write_input_geo()
  call write_input_n0g()
  call write_input_par()

  contains
  !---------------------------------------------------------------------
  subroutine write_input_geo
  integer :: ir, ip, it, iz, irun, n


  open  (iu, file='input.geo')
  write (iu, 1000)
 1000 format ('* geometry information for EMC3')

  ! 1. geometry, mesh resolution
  write (iu, 9999)
  write (iu, 1001)
  write (iu, 9999)
  write (iu, 1002)
  write (iu, 1003) NZONET
  write (iu, 1004)
  do iz=0,NZONET-1
     write (iu, *) SRF_RADI(iz), SRF_POLO(iz), SRF_TORO(iz)
  enddo
 1001 format ('*** 1. grid resolution')
 1002 format ('* number of zones/blocks')
 1003 format (i0)
 1004 format ('* number of radial, poloidal and toroidal grid points')


  ! 2. surface definitions
  write (iu, 9999)
  write (iu, 2000)
  write (iu, 9999)
  ! 2.1 non default surfaces (periodic, mapping, ...)
  write (iu, 2001)
  ! 2.1.a - radial
  write (iu, 2002)
  n = 0
  do irun=0,1
     ! write number of non default radial surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,NZONET-1
        if (TD(iz+1)%irsfa > 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) 0, iz, TD(iz+1)%irsfa
              write (iu, *) 0, ZON_POLO(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
        if (TD(iz+1)%irsfb > 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) SRF_RADI(iz)-1, iz, TD(iz+1)%irsfb
              write (iu, *) 0, ZON_POLO(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
     enddo
  enddo
  ! 2.1.b - poloidal
  write (iu, 2003)
  n = 0
  do irun=0,1
     ! write number of non default poloidal surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,NZONET-1
        if (TD(iz+1)%ipsfa > 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) 0, iz, TD(iz+1)%ipsfa
              write (iu, *) 0, ZON_RADI(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
        if (TD(iz+1)%ipsfb > 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) SRF_POLO(iz)-1, iz, TD(iz+1)%ipsfb
              write (iu, *) 0, ZON_RADI(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
     enddo
  enddo
  ! 2.1.c - toroidal
  write (iu, 2004)
  n = 0
  do irun=0,1
     ! write number of non default toroidal surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,NZONET-1
        if (TD(iz+1)%itsfa > 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) 0, iz, TD(iz+1)%itsfa
              write (iu, *) 0, ZON_RADI(iz)-1, 0, ZON_POLO(iz)-1
           endif
        endif
        if (TD(iz+1)%itsfb > 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) SRF_TORO(iz)-1, iz, TD(iz+1)%itsfb
              write (iu, *) 0, ZON_RADI(iz)-1, 0, ZON_POLO(iz)-1
           endif
        endif
     enddo
  enddo

  ! 2.2 non transparent surfaces (boundary conditions)
  write (iu, 2005)
  ! 2.2.a - radial
  write (iu, 2002)
  n = 0
  do irun=0,1
     ! write number of non transparent radial surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,NZONET-1
        if (TD(iz+1)%irsfa < 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) R_SURF_PL_TRANS_RANGE(1,iz), iz, 1
              write (iu, *) 0, ZON_POLO(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
        if (TD(iz+1)%irsfb < 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) R_SURF_PL_TRANS_RANGE(2,iz), iz, -1
              write (iu, *) 0, ZON_POLO(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
     enddo
  enddo
  ! 2.2.b - poloidal
  write (iu, 2003)
  write (iu, *) 0
  ! 2.2.c - toroidal
  write (iu, 2004)
  write (iu, *) 0

  ! 2.3 plate surfaces
  write (iu, 2006)
  write (iu, 2002)
  write (iu, *) -3 ! user defined
  write (iu, 2003)
  write (iu, *) -3 ! user defined
  write (iu, 2004)
  write (iu, *) -3 ! user defined
 2000 format ('*** 2. surface definitions')
 2001 format ('*** 2.1 non default surface')
 2002 format ('* radial')
 2003 format ('* poloidal')
 2004 format ('* toroidal')
 2005 format ('*** 2.2 non transparent surface (Boundary condition must be defined)')
 2006 format ('*** 2.3 plate surface (Bohm Boundary condition)')


  ! 3. physical cell definition
  write (iu, 9999)
  write (iu, 3000)
  write (iu, 9999)
  write (iu, *) -1
  write (iu, 3001)
  write (iu, 3002) .true.
 3000 format ('*** 3. physical cell definition')
 3001 format ('* run cell check?')
 3002 format (L1)
  close (iu)
 9999 format ('*',32('-'))

  end subroutine write_input_geo
  !---------------------------------------------------------------------
  subroutine write_input_n0g

  integer :: ir, iz, irun, n


  open  (iu, file='input.n0g')
  write (iu, 1000)
  write (iu, 9999)
  write (iu, 1001)
  write (iu, 9999)
  write (iu, 1002)
  write (iu, 1003)
  write (iu, 1004)
 1000 format ('******** additional geometry and parameters for EIRENE ****')
 1001 format ('*** 1. non-transparent surfaces for neutral particles')
 1002 format ('*  non-transparent surfaces with informations about')
 1003 format ('*  this surface being defined in EIRENE. The surface')
 1004 format ('*  number must be indicated here.')

  ! 1. non-transparent radial surfaces
  write (iu, 1012)
  n = 0
  do irun=0,1
     ! write number of non default radial surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,NZONET-1
        if (TD(iz+1)%irsfa < 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, 1015) 0, iz, -TD(iz+1)%irsfa
              write (iu, *) 0, ZON_POLO(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
        if (TD(iz+1)%irsfb < 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, 1015) SRF_RADI(iz)-1, iz, -TD(iz+1)%irsfb
              write (iu, *) 0, ZON_POLO(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
     enddo
  enddo
  ! non-transparent poloidal and toroidal surfaces
  write (iu, 1013)
  write (iu, *) 0
  write (iu, 1014)
  write (iu, *) 0
 1012 format ('* radial')
 1013 format ('* poloidal')
 1014 format ('* toroidal')
 1015 format (2i8,4x,'EIRENE_SF',i0)


  ! 2. additional physical cells for neutrals
  !ne0 = 1.d14;	Te0 = 4.d3;	Ti0 = 4.d3;	M0  = 0.d0
  !ne1 = 1.d7;	Te1 = 1.d-1;	Ti1 = 1.d-1;	M1  = 0.d0
  write (iu, 9999)
  write (iu, 2000)
  write (iu, 9999)
  write (iu, 2001)
  write (iu, 2002)
  n = 0
  do irun=0,1
     ! write number of additional cell blocks
     if (irun == 1) write (iu, *) n, 70

     ! confined region
     do iz=0,NZONET-1
        if (TD(iz+1)%irsfa < 0  .and.  R_SURF_PL_TRANS_RANGE(1,iz) > 0) then
           if (irun == 0) then
              n = n + R_SURF_PL_TRANS_RANGE(1,iz)
           else
              do ir=0,R_SURF_PL_TRANS_RANGE(1,iz)-1
                 write (iu, *) 2, 1
                 write (iu, 2003) iz, ir, ir+1, 1, &
                           0, ZON_POLO(iz), ZON_POLO(iz), &
                           0, ZON_TORO(iz), ZON_TORO(iz)
                 write (iu, 2004) iz, ir
              enddo
           endif
        endif
     enddo

     ! vacuum region (collect what's left)
     do iz=0,NZONET-1
        if (irun == 0) then
           n = n + 1
        else
           write (iu, *) 2, 0
           write (iu, 2003) iz, 0, SRF_RADI(iz)-1, 1, &
                     0, ZON_POLO(iz), ZON_POLO(iz), &
                     0, ZON_TORO(iz), ZON_TORO(iz)
           write (iu, 2005) iz
        endif
     enddo
  enddo
 2000 format ('*** 2. DEFINE ADDITIONAL PHYSICAL CELLS FOR NEUTRALS')
 2001 format ('*   ZONE  R1    R2    DR    P1    P2    DP    T1    T2    DT')
 2002 format ('* ne       Te      Ti        M')
 2003 format (10i6)
 2004 format ('EIRENE_CORE_',i0,'_',i0)
 2005 format ('EIRENE_VACUUM_',i0)


  ! 3. Neutral sources
  write (iu, 9999)
  write (iu, 3000)
  write (iu, 9999)
  write (iu, 3001)
  write (iu, *) 0, 0, -1
 3000 format ('*** 3. Neutral Source distribution')
 3001 format ('* N0S NS_PLACE  NSSIDE')


  ! 4. Additional surfaces
  write (iu, 9999)
  write (iu, 4000)
  write (iu, 9999)
  write (iu, 4001)
  close (iu)
 4000 format ('*** 4 Additional surfaces')
 4001 format ('./../../geometry/ADD_SF_N0')
  close (iu)
 9999 format ('*',32('-'))

  end subroutine write_input_n0g
  !---------------------------------------------------------------------
  subroutine write_input_par

  integer :: ir, iz, irun, n


  open  (iu, file='input.par.6')
  write (iu, 9999)
  write (iu, 1000)
  write (iu, 9999)
 1000 format ('*** 6. boundary contitions')

  ! 1. particle transport
  write (iu, 1001)
 1001 format ('*** 6.1 particle transport for main ions + impurities')
  write (iu, 1002)
  write (iu, 1003)
  do irun=0,1
     ! write number of non default radial surfaces
     if (irun == 1) write (iu, *) n
     n = 0

     do iz=0,NZONET-1
        if (TD(iz+1)%irsfa < 0) then
           n = n + 1
           if (irun == 1) then
              write (iu, 1004) 1, n, -TD(iz+1)%irsfa, -TD(iz+1)%irsfa
           endif
        endif
        if (TD(iz+1)%irsfb < 0) then
           n = n + 1
           if (irun == 1) then
              write (iu, 1004) 1, n, -TD(iz+1)%irsfb, -TD(iz+1)%irsfb
           endif
        endif
     enddo
  enddo
 1002 format ('* Main plasma ions')
 1003 format ('NBUND_TYE    CBUND_COE')
 1004 format (2i4,4x,'EMC3_SF',i0,4x,'EMC3_SF',i0,'_PVAL')

  ! 2. energy transport
  write (iu, 2001)
 2001 format ('*** 6.2 energy transport for el. + ions')
  do irun=0,1
     ! write number of non default radial surfaces
     if (irun == 1) write (iu, *) n
     n = 0

     do iz=0,NZONET-1
        if (TD(iz+1)%irsfa < 0) then
           n = n + 1
           if (irun == 1) then
              write (iu, 2002) 1, n, -TD(iz+1)%irsfa, -TD(iz+1)%irsfa, -TD(iz+1)%irsfa
           endif
        endif
        if (TD(iz+1)%irsfb < 0) then
           n = n + 1
           if (irun == 1) then
              write (iu, 2002) 1, n, -TD(iz+1)%irsfb, -TD(iz+1)%irsfb, -TD(iz+1)%irsfb
           endif
        endif
     enddo
  enddo
 2002 format (2i4,4x,'EMC3_SF',i0,4x,'EMC3_SF',i0,'_EEVAL',4x,'EMC3_SF',i0,'_EIVAL')

  ! 3. momentum transport
  write (iu, 3001)
 3001 format ('*** 6.3 momentum transport')
  do irun=0,1
     ! write number of non default radial surfaces
     if (irun == 1) write (iu, *) n
     n = 0

     do iz=0,NZONET-1
        if (TD(iz+1)%irsfa < 0) then
           n = n + 1
           if (irun == 1) then
              write (iu, 3002) 1, n, -TD(iz+1)%irsfa, -TD(iz+1)%irsfa
           endif
        endif
        if (TD(iz+1)%irsfb < 0) then
           n = n + 1
           if (irun == 1) then
              write (iu, 3002) 1, n, -TD(iz+1)%irsfb, -TD(iz+1)%irsfb
           endif
        endif
     enddo
  enddo
 3002 format (2i4,4x,'EMC3_SF',i0,4x,'EMC3_SF',i0,'_MVAL')

  close (iu)
 9999 format ('*',32('-'))

  end subroutine write_input_par
  !---------------------------------------------------------------------
  end subroutine write_emc3_input_files
!=======================================================================



!=======================================================================
  subroutine generate_plates()
  use iso_fortran_env
  use emc3_grid
  use field_aligned_grid, only: N_sym
  use boundary
  use curve2D
  use math
  use run_control, only: Debug
  use string
  use dataset
  implicit none

  integer, parameter :: iu = 78
  real(real64), parameter :: l0 = 10.d-10

  type(t_curve), dimension(:,:), allocatable :: C

  real(real64), dimension(:), allocatable   :: RC_TEM, ZC_TEM
  integer, dimension(:), allocatable   :: ID_TEM
  integer, dimension(:,:), allocatable :: iindex ! number of cells behind a plate
  integer, dimension(:), allocatable   :: knumb  ! cell index in flux tube for plate cells

  logical      :: plate_cell
  real(real64) :: x(3), phi
  integer      :: nr, np, nt, iz, i, j, k, l, l1, l2, irun, icut, ig(8), ic


  !call outside_boundary_check()
  allocate (ID_TEM(0:MESH_P_OS(NZONET)-1))
  allocate (RC_TEM(0:MESH_P_OS(NZONET)-1))
  allocate (ZC_TEM(0:MESH_P_OS(NZONET)-1))
  ID_TEM = 0
  open  (iu, file='plates.dat')
  do iz=0,NZONET-1
  !do iz=0,0
     nr = ZON_RADI(iz)
     np = ZON_POLO(iz)
     nt = ZON_TORO(iz)
     allocate (iindex(0:nr-1, 0:np-1))
     iindex = 0

     write (6, 1000) iz, nr, np, nt
     1000 format (8x,'zone ',i0,' (',i0,' x ',i0,' x ',i0,')')

     ! setup slices for Q4-type surfaces
     allocate (C(0:nt-1, n_quad))
     do k=0,nt-1
     do l=1,n_quad
         phi  = (PHI_PLANE(k+1+PHI_PL_OS(iz)) + PHI_PLANE(k+PHI_PL_OS(iz))) / 2.d0
         phi  = phi_sym(phi, N_sym)
         C(k,l) = S_quad(l)%slice(phi)

         if (Debug) then
            call C(k,l)%plot(filename='debug/Q4surf_'//trim(str(l))//'_zone'//trim(str(iz))// &
                                      '_t'//trim(str(k)))
         endif
     enddo
     enddo


     ! (0) count plate cells, (1) setup plate cells
     do irun=0,1
        icut = 0
        ! loop over all flux tubes
        do i=R_SURF_PL_TRANS_RANGE(1,iz),R_SURF_PL_TRANS_RANGE(2,iz)-1
        do j=P_SURF_PL_TRANS_RANGE(1,iz),P_SURF_PL_TRANS_RANGE(2,iz)-1
        !do i=0,nr-1
        !do j=0,np-1
        !do i=23,23
        !do j=120,120
        do k=0,nt-1
           ig(1)   = i + (j + k*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
           ig(2)   = ig(1) + 1
           ig(3)   = ig(2) + SRF_RADI(iz)
           ig(4)   = ig(3) - 1
           ig(5:8) = ig(1:4) + SRF_POLO(iz)*SRF_RADI(iz)

           x(1)    = sum(RG(ig))/8.d0
           x(2)    = sum(ZG(ig))/8.d0
           x(3)    = (PHI_PLANE(k+1+PHI_PL_OS(iz)) + PHI_PLANE(k+PHI_PL_OS(iz))) / 2.d0
           !write (6, *) x
           x(3)    = phi_sym(x(3), N_sym)

           ic      = i + (j + k*ZON_POLO(iz))*ZON_RADI(iz) + MESH_P_OS(iz)
           RC_TEM(ic) = x(1)
           ZC_TEM(ic) = x(2)

           !plate_cell = .false.
           plate_cell = outside_boundary()

           if (plate_cell) then
              ID_TEM(ic) = 1
              icut = icut + 1
              if (irun == 1) then
                 iindex(i,j) = iindex(i,j) + 1
                 knumb(icut) = k
              endif
           endif
        enddo
        enddo
        enddo

        ! allocate knumb array after 1st run
        if (irun == 0) then
           allocate (knumb(icut))
           knumb = 0

        ! write plate cells after 2nd run
        else
           icut = 0
           do i=0,nr-1
           do j=0,np-1
              if (iindex(i,j) .ne. 0) then
                 l1 = icut + 1
                 l2 = icut + iindex(i,j)
                 write (iu, *) iz, i, j, iindex(i,j),(knumb(l),l=l1,l2)
                 icut = icut + iindex(i,j)
              endif
           enddo
           enddo

           write (6, 2000) icut
           2000 format (8x,i0,' plate cells')
           deallocate (knumb, iindex)
        endif
     enddo

     ! cleanup slices of Q4-type surfaces
     do k=0,nt-1
     do l=1,n_quad
        call C(k,l)%destroy()
     enddo
     enddo
     deallocate (C)
  enddo
  close (iu)


  ! additional output for debugging
  if (Debug) then
     do iz=0,NZONET-1
     do k=0,ZON_TORO(iz)-1
        open  (99, file='debug/SlicePlasma_zone'//trim(str(iz))//'_t'//trim(str(k))//'.plt')
        open  (98, file='debug/SliceBoundary_zone'//trim(str(iz))//'_t'//trim(str(k))//'.plt')

        do i=R_SURF_PL_TRANS_RANGE(1,iz),R_SURF_PL_TRANS_RANGE(2,iz)-1
        do j=P_SURF_PL_TRANS_RANGE(1,iz),P_SURF_PL_TRANS_RANGE(2,iz)-1
           ic      = i + (j + k*ZON_POLO(iz))*ZON_RADI(iz) + MESH_P_OS(iz)
           if (ID_TEM(ic) == 1) then
              write (98, *) RC_TEM(ic), ZC_TEM(ic)
           else
              write (99, *) RC_TEM(ic), ZC_TEM(ic)
           endif
        enddo
        enddo
        close (99)
        close (98)
     enddo
     enddo

     call plate_check_all()
  endif
  deallocate (ID_TEM, RC_TEM, ZC_TEM)
  contains
  !-------------------------------------------------------------------
  function outside_boundary()
  ! input:
  !    k: toroidal index
  !    C: slices of Q4-type surfaces at phi(k)
  !    x: reference point
  logical :: outside_boundary


  ! set default
  outside_boundary = .false.


  ! check axisymmetric (L2-type) surfaces
  do l=1,n_axi
     if (S_axi(l)%outside(x(1:2))) then
        outside_boundary = .true.
        return
     endif
  enddo

  ! check Q4-type surfaces
  do l=1,n_quad
     if (C(k,l)%outside(x(1:2))) then
        outside_boundary = .true.
        return
     endif
  enddo

  ! check block limiters (CSG-type)
  do l=1,n_block
     if (bl_outside(l, x)) then
        outside_boundary = .true.
        return
     endif
  enddo

  end function outside_boundary
  !-------------------------------------------------------------------


  !-------------------------------------------------------------------
  subroutine outside_boundary_check
  use field_aligned_grid, only: N_sym

  integer, parameter :: iu = 99, iu2 = 98, n = 101

  character(len=72) :: s
  real(real64) :: x(2), Phi
  integer :: i

  allocate (C(0:n-1, n_quad))
  open  (iu, file='fl+.dat')
  open  (iu2, file='fl+B.dat')
  read  (iu, *) s
  do i=0,n-1
     read (iu, *) x, Phi

     phi  = phi_sym(phi, N_sym)
     C(i,1) = S_quad(1)%slice(phi)

     if (C(i,1)%outside(x)) then
        !write (iu2, *) Phi, 1
     else
        !write (iu2, *) Phi, 0
        write (iu2, *) x, Phi
     endif
  enddo
  close (iu)
  close (iu2)

  stop
  end subroutine outside_boundary_check
  !-------------------------------------------------------------------


  !-------------------------------------------------------------------
  subroutine plate_check_all
  integer :: iz, k, j, i, ic1, ic2, iplate


  write (6, *) 'running plate checks ...'

  ID_TEM = ID_TEM*2 - 1
  iplate = 0

  do iz=0,NZONET-1
  do i=0,ZON_RADI(iz)-1
  do j=0,ZON_POLO(iz)-1
     do k=0,ZON_TORO(iz)-2
        ic1 = i + (j +  k   *ZON_POLO(iz))*ZON_RADI(iz) + MESH_P_OS(iz)
        ic2 = i + (j + (k+1)*ZON_POLO(iz))*ZON_RADI(iz) + MESH_P_OS(iz)

        if (ID_TEM(ic1)*ID_TEM(ic2) < 0) then
           iplate = iplate + 1
           write (6, *) iplate
           call plate_check(iz,i,j,k+1)
        endif
     enddo
  enddo
  enddo
  enddo

  end subroutine plate_check_all
  !-------------------------------------------------------------------
  subroutine plate_check(iz, ir, ip, jt)
  use Q4
  use math
  use grid
  use fieldline
  use dataset
  integer, intent(in) :: iz, ir, ip, jt

  real(real64), parameter :: Limit = 360.d0

  type(t_Q4)   :: Q
  type(t_grid) :: G
  type(t_fieldline) :: F
  type(t_dataset)   :: D
  real(real64) :: x1(2), x2(2), x3(2), x4(2), phi, y(3), ts, Lc, Lcsav
  integer      :: ig(4), i, idir, n, nsuccess


  ! calculate node indices
  ig(1)   = ir + (ip + jt*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
  ig(2)   = ig(1) + 1
  ig(3)   = ig(2) + SRF_RADI(iz)
  ig(4)   = ig(3) - 1

  ! get node coordinates
  x1(1)   = RG(ig(1)); x1(2)   = ZG(ig(1))
  x2(1)   = RG(ig(2)); x2(2)   = ZG(ig(2))
  x3(1)   = RG(ig(3)); x3(2)   = ZG(ig(3))
  x4(1)   = RG(ig(4)); x4(2)   = ZG(ig(4))
  phi     = PHI_PLANE(jt+PHI_PL_OS(iz))

  ! setup quadrilateral
  call Q%set_nodes(x1, x2, x3, x4)

  ! generate mesh
  n = 10
  G = Q%generate_mesh(n, n, phi)
  G%coordinates       = CYLINDRICAL
  !G%fixed_coord       = 3
  !G%fixed_coord_value = phi / 180.d0 * pi

  ! sample connection length on grid
  ts       = 1.d0
  Lcsav    = 0.d0
  nsuccess = 0
  call D%new(G%nodes(),5)
  do i=1,G%nodes()
     y = G%node(i)

     D%x(i,1) = y(1)
     D%x(i,2) = y(2)

     ! trace field line in both directions
     do idir=-1,1,2
        call F%init(y, idir*ts, NM_AdamsBashforth4, FL_ARC)

        trace_loop: do
           call F%trace_1step()
           Lc = F%phi_int * 180.d0 / pi

           if (abs(Lc) > Limit) exit trace_loop

           if (F%intersect_boundary()) exit trace_loop
        enddo trace_loop

        D%x(i,3 + (idir+1)/2) = Lc
     enddo
     D%x(i,5) = min(abs(D%x(i,3)), abs(D%x(i,4)))

     ! success frequency, average shortest toroidal distance
     if (D%x(i,5) < Limit) then
        nsuccess = nsuccess + 1
        Lcsav    = Lcsav    + D%x(i,5)
     endif
  enddo
  Lcsav = Lcsav / nsuccess
  write (6, *) 'success frequency = ', 100.d0 * nsuccess / G%nodes(), ' %'
  write (6, *) 'average toroidal distance to plates [deg] = ', Lcsav
  call G%store('plate1.grid')
  call D%plot(filename='plate1.dat')

  call D%destroy()
  call G%destroy()
  stop

  end subroutine plate_check
  !-------------------------------------------------------------------
  end subroutine generate_plates
!=======================================================================

#endif
