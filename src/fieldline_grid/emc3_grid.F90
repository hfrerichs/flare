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
     RG,                  ZG,                    PHI_PLANE, &
     BFSTREN,             PSI_N

  logical, private, save :: already_loaded = .false.


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
  integer :: i, iz, ir, ip, itmp, n


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
     P_SURF_PL_TRANS_RANGE(1,iz) = SRF_POLO(iz)-1
     P_SURF_PL_TRANS_RANGE(2,iz) = 0
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
     read  (iu, *) ip, iz, itmp
     write (6,  *) ip, iz, itmp
     read  (iu, *) readfi
     P_SURF_PL_TRANS_RANGE(1,iz) = min(P_SURF_PL_TRANS_RANGE(1,iz), ip)
     P_SURF_PL_TRANS_RANGE(2,iz) = max(P_SURF_PL_TRANS_RANGE(2,iz), ip)
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

  ! 2. poloidal
  call scrape(iu, readfi)
  read (readfi, *) n
  do i=1,n
     read  (iu, *) ip, iz, itmp
     write (6,  *) ip, iz, itmp
     read  (iu, *) readfi
     P_SURF_PL_TRANS_RANGE(1,iz) = min(P_SURF_PL_TRANS_RANGE(1,iz), ip)
     P_SURF_PL_TRANS_RANGE(2,iz) = max(P_SURF_PL_TRANS_RANGE(2,iz), ip)
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


  if (already_loaded) return

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
  already_loaded = .true.

  end subroutine load_emc3_grid
!=======================================================================



!=======================================================================
  subroutine load_bfstren()
  implicit none

  integer, parameter :: iu = 24


  if (allocated(BFSTREN)) deallocate(BFSTREN)
  allocate (BFSTREN(0:GRID_P_OS(NZONET)-1))

  open  (iu, file='bfield.dat')
  read  (iu, *) BFSTREN
  close (iu)

  end subroutine load_bfstren
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
  already_loaded = .true.

  end subroutine write_emc3_grid
!=======================================================================



#ifdef FLARE
!=======================================================================
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
!=======================================================================



!=======================================================================
  subroutine export_slice(iz, it, ir1, ir2, G)
  use grid

  integer,      intent(in)  :: iz, it, ir1, ir2
  type(t_grid), intent(out) :: G

  integer :: nr, np, ir, ip, ig


  nr = ir2 - ir1 + 1
  np = SRF_POLO(iz)
  call G%new(2, MESH_2D, 3, nr, np, fixed_coord_value=PHI_PLANE(it+PHI_PL_OS(iz)))
  do ir=ir1,ir2
  do ip=0,SRF_POLO(iz)-1
     ig = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
     G%mesh(ir-ir1,ip,1) = RG(ig)
     G%mesh(ir-ir1,ip,2) = ZG(ig)
  enddo
  enddo

  end subroutine export_slice
!=======================================================================



!=======================================================================
  function flux_tube_length(iz, ir, ip) result(L)
  use iso_fortran_env
  integer, intent(in) :: iz, ir, ip
  real(real64)        :: L

  real(real64) :: R1, Z1, phi1, R2, Z2, phi2, dl
  integer      :: it, ig(4)

            ! get length of finite, interpolated field line
  L = 0.d0
  do it=0,ZON_TORO(iz)-1
     ig(1) = ir+(ip+it*SRF_POLO(iz))*SRF_RADI(iz)+GRID_P_OS(iz)
     ig(2) = ig(1) + SRF_RADI(iz)
     ig(3) = ig(2) + 1
     ig(4) = ig(1) + 1

     R2    = 0.25d0 * sum(RG(ig))
     Z2    = 0.25d0 * sum(ZG(ig))
     phi2  = PHI_PLANE(PHI_PL_OS(iz)+it+1)

     if (it > 0) then
        dl = dsqrt((0.5d0*(R1+R2)*(phi2-phi1))**2  +  (R2-R1)**2  +  (Z2-Z1)**2)
        L  = L + dl
     endif

     R1    = R2
     Z1    = Z2
     phi1  = phi2
  enddo

  end function flux_tube_length
!=======================================================================



!=======================================================================
      SUBROUTINE RZ_REAL_COORDINATES(NZ0,JR0,JP0,RJ0,PJ0,TJ0,R,Z)

      INTEGER,INTENT(IN) :: NZ0,JR0,JP0
      REAL*8, INTENT(IN) :: RJ0,PJ0,TJ0
      REAL*8, INTENT(OUT):: R,Z

      INTEGER :: JT0, I1, I2, I3, I4
      REAL*8  :: AR, BR, CR, DR, AZ, BZ, CZ, DZ, R1, Z1, R2, Z2, XT
      JT0 = TJ0
      XT  = TJ0 - JT0

      I1= JR0+(JP0+JT0*SRF_POLO(NZ0))*SRF_RADI(NZ0)+GRID_P_OS(NZ0)
      I2= I1 + SRF_RADI(NZ0)
      I3= I2 + 1
      I4= I1 + 1

      AR  = 0.25*(RG(I1)+RG(I2)+RG(I3)+RG(I4))
      BR  = 0.25*(RG(I3)+RG(I4)-RG(I1)-RG(I2))
      CR  = 0.25*(RG(I3)+RG(I2)-RG(I1)-RG(I4))
      DR  = 0.25*(RG(I1)+RG(I3)-RG(I2)-RG(I4))

      AZ  = 0.25*(ZG(I1)+ZG(I2)+ZG(I3)+ZG(I4))
      BZ  = 0.25*(ZG(I3)+ZG(I4)-ZG(I1)-ZG(I2))
      CZ  = 0.25*(ZG(I3)+ZG(I2)-ZG(I1)-ZG(I4))
      DZ  = 0.25*(ZG(I1)+ZG(I3)-ZG(I2)-ZG(I4))
      R1  = AR + BR*RJ0 + CR*PJ0 + DR*RJ0*PJ0
      Z1  = AZ + BZ*RJ0 + CZ*PJ0 + DZ*RJ0*PJ0

      IF ( XT == 0. ) THEN
        R   = R1
        Z   = Z1
        RETURN
      ENDIF

      I1= I1 + SRF_POLO(NZ0)*SRF_RADI(NZ0)
      I2= I1 + SRF_RADI(NZ0)
      I3= I2 + 1
      I4= I1 + 1

      AR  = 0.25*(RG(I1)+RG(I2)+RG(I3)+RG(I4))
      BR  = 0.25*(RG(I3)+RG(I4)-RG(I1)-RG(I2))
      CR  = 0.25*(RG(I3)+RG(I2)-RG(I1)-RG(I4))
      DR  = 0.25*(RG(I1)+RG(I3)-RG(I2)-RG(I4))

      AZ  = 0.25*(ZG(I1)+ZG(I2)+ZG(I3)+ZG(I4))
      BZ  = 0.25*(ZG(I3)+ZG(I4)-ZG(I1)-ZG(I2))
      CZ  = 0.25*(ZG(I3)+ZG(I2)-ZG(I1)-ZG(I4))
      DZ  = 0.25*(ZG(I1)+ZG(I3)-ZG(I2)-ZG(I4))

      R2  = AR + BR*RJ0 + CR*PJ0 + DR*RJ0*PJ0
      Z2  = AZ + BZ*RJ0 + CZ*PJ0 + DZ*RJ0*PJ0

      R   = R1 + (R2-R1)*XT
      Z   = Z1 + (Z2-Z1)*XT
      RETURN
      END SUBROUTINE RZ_REAL_COORDINATES
!=======================================================================



!=======================================================================
  subroutine reconstruct_field_line(iz, ir, ip, rj, pj, x)
  use iso_fortran_env
  integer,      intent(in)  :: iz, ir, ip
  real(real64), intent(in)  :: rj, pj
  real(real64), intent(out) :: x(0:ZON_TORO(iz),3)

  real(real64) :: R, Z, tj
  integer :: it


  do it=0,ZON_TORO(iz)
     tj      = 1.d0 * it
     call RZ_REAL_COORDINATES(iz, ir, ip, rj, pj, tj, R, Z)
     x(it,1) = R
     x(it,2) = Z
     x(it,3) = PHI_PLANE(it + PHI_PL_OS(iz))
  enddo

  end subroutine reconstruct_field_line
!=======================================================================



!=======================================================================
  subroutine check_mesh(iz, div, NL, pitch)
  use iso_fortran_env
  implicit none

  integer,      intent(in)  :: iz
  real(real64), intent(out) :: div(0:ZON_POLO(iz)-1, 0:ZON_RADI(iz)-1)
  real(real64), intent(out) :: pitch(0:ZON_POLO(iz)-1, 0:ZON_RADI(iz)-1)
  real(real64), intent(out) :: NL(0:ZON_POLO(iz)-1, 0:ZON_RADI(iz)-1, 0:ZON_TORO(iz))

      INTEGER :: I,J,K,I1,I2,I3,I4,IS,IC,IG
      REAL*8  :: A1,A2,B1,B2,C1,C2,D1,D2,D,ABCD1,ABCD2, &
                 NON_LINEAR,NON_LINEAR_AVER, &
                 AREA1,AREA2,R1,R2,Z1,Z2,BF1,BF2,FI1,FI2,DFL, &
                 AVER_F,DIVMAX,DIVAVE,F_SUM    
      REAL*8,DIMENSION(:  ),ALLOCATABLE :: FLUX,DIV_CELL
      LOGICAL :: S_S,PROB                      

      DIVMAX = 0.
      DIVAVE = 0.
      F_SUM  = 0.
      div    = 0.d0
      NL     = 0.d0
      NON_LINEAR_AVER = 0.

      ALLOCATE(FLUX(ZON_TORO(IZ)))

      DO 10 I =0,ZON_RADI(IZ)-1
      DO 10 J =0,ZON_POLO(IZ)-1
      DO 20 K =0,ZON_TORO(IZ)
         I1 = I + (J+K*SRF_POLO(IZ))*SRF_RADI(IZ) + GRID_P_OS(IZ)
         I2 = I1+ SRF_RADI(IZ)
         I3 = I2 + 1
         I4 = I1 + 1

         ABCD1 = (RG(I3)-RG(I2))*(ZG(I2)-ZG(I1)) - (RG(I2)-RG(I1))*(ZG(I3)-ZG(I2))
         ABCD2 = (RG(I4)-RG(I1))*(ZG(I3)-ZG(I4)) - (RG(I3)-RG(I4))*(ZG(I4)-ZG(I1))

         A1 = 0.25*(RG(I3)+RG(I4)-RG(I1)-RG(I2))
         B1 = 0.25*(RG(I3)+RG(I2)-RG(I1)-RG(I4))
         C1 = 0.25*(RG(I1)+RG(I3)-RG(I2)-RG(I4))

         A2 = 0.25*(ZG(I3)+ZG(I4)-ZG(I1)-ZG(I2))
         B2 = 0.25*(ZG(I3)+ZG(I2)-ZG(I1)-ZG(I4))
         C2 = 0.25*(ZG(I1)+ZG(I3)-ZG(I2)-ZG(I4))

         D  = ABS(A1*B2 - A2*B1)

         D1 = -(C1*B2-C2*B1)
         D2 = -(A1*C2-A2*C1)
         NON_LINEAR = ABS(D1)+ABS(D2)
         NL(J,I,K) = NON_LINEAR/D
         NON_LINEAR_AVER = NON_LINEAR_AVER + NON_LINEAR/D

! check flux conservation in plasma transport range

         IF(I>=R_SURF_PL_TRANS_RANGE(1,IZ) .AND. I< R_SURF_PL_TRANS_RANGE(2,IZ) .AND. &
            J>=P_SURF_PL_TRANS_RANGE(1,IZ) .AND. J< P_SURF_PL_TRANS_RANGE(2,IZ) ) THEN 

            AREA2=0.5*ABS(ABCD1+ABCD2)
            R2   = 0.25*(RG(I1)+RG(I2)+RG(I3)+RG(I4))
            Z2   = 0.25*(ZG(I1)+ZG(I2)+ZG(I3)+ZG(I4))
            FI2  = PHI_PLANE(PHI_PL_OS(IZ)+K)
            BF2  = 0.25*( BFSTREN(I1)+BFSTREN(I2)+BFSTREN(I3)+BFSTREN(I4) )

            IF(K/=0) THEN
             DFL    = 0.5*(R1+R2)*ABS(FI2-FI1)
             PITCH(J,I)  = DFL/SQRT(DFL**2+(R2-R1)**2+(Z2-Z1)**2)
             FLUX(K)= (AREA1+AREA2)*PITCH(J,I)*(BF2+BF1)*0.25
            ENDIF
            AREA1=AREA2; R1=R2; Z1=Z2; FI1=FI2; BF1=BF2
         ENDIF
20    CONTINUE
      IF(I>=R_SURF_PL_TRANS_RANGE(1,IZ) .AND. I< R_SURF_PL_TRANS_RANGE(2,IZ) .AND. &
         J>=P_SURF_PL_TRANS_RANGE(1,IZ) .AND. J< P_SURF_PL_TRANS_RANGE(2,IZ) ) THEN
         AVER_F = SUM(FLUX)/FLOAT(ZON_TORO(IZ))
         DIV(J,I) = (MAXVAL(FLUX)-MINVAL(FLUX))/AVER_F
         DIVMAX = MAX(DIVMAX,DIV(J,I))
         DIVAVE = DIVAVE+MAXVAL(FLUX)-MINVAL(FLUX)   
         F_SUM  = F_SUM + AVER_F
      ENDIF
10    CONTINUE
      DEALLOCATE(FLUX)


      NON_LINEAR_AVER = NON_LINEAR_AVER/FLOAT(MESH_P_OS(NZONET))
      WRITE(6,*)' 1.1 mesh shape '
      WRITE(6,'(A18,F8.3,A2)')'   Non-linearity:', NON_LINEAR_AVER*100.,' %'

      DIVAVE = DIVAVE/(F_SUM+1.E-10)*100.
      WRITE(6,*) ' 1.2 Check accuracy in flux conservation'
      WRITE(6,'(A18,F8.3,A2)')'   Flux conserved:',100.-DIVAVE,' %'
      WRITE(6,'(A18,F8.3,A2)')'   max. diviation:',DIVMAX*100.,' %'
      IF(DIVAVE>1. .OR. DIVMAX>0.1) THEN
        WRITE(6,*)
        WRITE(6,*)' Mag. flux not conserved inside flux tubes'
      ENDIF

  end subroutine check_mesh
!=======================================================================
#endif



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
