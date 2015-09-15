C Name    : GEOMETRY_PL
C Function: Fine geometric grids 
C Contains:
C   Subroutines  : READ_GEOMETRY_PL
C   Functions    : None
C Use     :
C   Subroutines  : SCRAPE             STOP_ALL 
C   Functions    : None
C   Common blocks: None
C   Modules      : PHYS_CONST  
C

      MODULE GEOMETRY_PL
      IMPLICIT NONE
C 1. Global mesh
C  Name         type  dimension   remark
C NZONET         I(S)     0     total zone number
C SRF_RADI       I(D)     1     radial surface number
C SRF_POLO       I(D)     1     poloi.                
C SRF_TORO       I(D)     1     toroi.               
C ZON_RADI       I(D)     1     radial intervals         
C ZON_POLO       I(D)     1     poloi.                
C ZON_TORO       I(D)     1     toroi.               
C GRID_P_OS      I(D)     1     offset of index of grid points 
C MESH_P_OS      I(D)     1     offset of index of meshes
C PHI_PL_OS      I(D)     1     offset of index of phi-cross-sections
C PHI_PLANE      R(D)     1     PHI angle of each cross-section
C RG             R(D)     1     R-coordinates of grid points
C ZG             R(D)     1     Z-coordinates     
C VOL3D          R(D)     1     volume of geometric cell
C BFSTREN        R(D)     1     B-field strength 
      INTEGER,SAVE                           :: NZONET
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE  ::
     .        SRF_RADI,            SRF_POLO,              SRF_TORO 
     .,       ZON_RADI,            ZON_POLO,              ZON_TORO
     .,      GRID_P_OS,           MESH_P_OS,             PHI_PL_OS
      REAL*8, DIMENSION(:),ALLOCATABLE,SAVE  ::
     .       PHI_PLANE,                  RG,                    ZG    
     .,      VOL3D,                 BFSTREN,                 LCELL

      CONTAINS                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE READ_GEOMETRY_PL(IU1,IU2,IU3)
      USE PHYS_CONST
      IMPLICIT NONE
      INTEGER,INTENT(IN)       :: IU1,IU2,IU3
C IU1   : uint number of file containing general information

      INTEGER :: IZ,    IO,     I1,     I2,   I,    K
     .,          N1,    N2,     N3,    IERR
      CHARACTER*72 :: READFI,MESSAG 
C------------------------------ ---------------------------------------
C *** 1. fein grid mesh
CV: NZONET
      CALL SCRAPE(IU1,READFI,*992) 
      READ(READFI,*) NZONET              

      IF    (NZONET.LE.0 ) THEN
        MESSAG = 'NZONET.LE.0'
        GOTO 991
      ENDIF

      WRITE(6,'(A7,I6)')'NZONET:',NZONET
      WRITE(6,*        )'   IZ   SRF_RADI  SRF_POLO  SRF_TORO'

CV: SRF_RADI
CV: SRF_POLO
CV: SRF_TORO
CV: ZON_RADI
CV: ZON_POLO
CV: ZON_TORO
      ALLOCATE (
     . SRF_RADI(0:NZONET-1),SRF_POLO(0:NZONET-1),SRF_TORO(0:NZONET-1)
     .,ZON_RADI(0:NZONET-1),ZON_POLO(0:NZONET-1),ZON_TORO(0:NZONET-1)
     .         )
      DO I=0,NZONET-1
          CALL SCRAPE(IU1,READFI,*992) 
          READ(READFI,*) SRF_RADI(I),SRF_POLO(I),SRF_TORO(I)
          IF( SRF_RADI(I) < 2 .OR. SRF_POLO(I) < 2
     .                        .OR. SRF_TORO(I) < 2 )THEN 
            MESSAG = 'SRF_RADI,_POLO or SRF_TORO < 2' 
            GOTO 991
          ENDIF 
          ZON_RADI(I) = SRF_RADI(I) - 1
          ZON_POLO(I) = SRF_POLO(I) - 1
          ZON_TORO(I) = SRF_TORO(I) - 1
          WRITE(6,'(I5,3I10)')I,SRF_RADI(I),SRF_POLO(I),SRF_TORO(I)
      ENDDO

CV: GRID_P_OS
CV: MESH_P_OS
CV: PHI_PL_OS
      ALLOCATE (
     .  GRID_P_OS(0:NZONET),MESH_P_OS(0:NZONET),PHI_PL_OS(0:NZONET)
     .         )

      GRID_P_OS(0) = 0
      MESH_P_OS(0) = 0
      PHI_PL_OS(0) = 0
      DO IZ=1,NZONET
        PHI_PL_OS(IZ) = PHI_PL_OS(IZ-1)+SRF_TORO(IZ-1)
        MESH_P_OS(IZ) = MESH_P_OS(IZ-1) 
     .+                 ZON_RADI(IZ-1)*ZON_POLO(IZ-1)*ZON_TORO(IZ-1)
        GRID_P_OS(IZ) = GRID_P_OS(IZ-1) 
     .+                 SRF_RADI(IZ-1)*SRF_POLO(IZ-1)*SRF_TORO(IZ-1)
      ENDDO

CV: PHI_PLANE
CV: RG
CV: ZG
      ALLOCATE (PHI_PLANE(0:PHI_PL_OS(NZONET)-1) 
     .,                RG(0:GRID_P_OS(NZONET)-1)
     .,                ZG(0:GRID_P_OS(NZONET)-1)
     .,           BFSTREN(0:GRID_P_OS(NZONET)-1)
     .         )

        DO IZ=0,NZONET-1
        READ(IU2,*,IOSTAT=IO)N1,N2,N3 
        IF(IO /= 0 ) GOTO 990

        IF( N1 /= SRF_RADI(IZ) .OR. N2 /= SRF_POLO(IZ)
     .   .OR. N3 /= SRF_TORO(IZ) ) THEN
          MESSAG = 'Sf. number not meet that in GRID_3D_DATA'
          GOTO 991 
        ENDIF

        DO K =0,SRF_TORO(IZ)-1
          READ(IU2,*,IOSTAT=IO)PHI_PLANE(PHI_PL_OS(IZ)+K)
          IF(IO /= 0 ) GOTO 990
        
          I1 = GRID_P_OS(IZ) + K*SRF_RADI(IZ)*SRF_POLO(IZ)
          I2 = I1  + SRF_RADI(IZ)*SRF_POLO(IZ) - 1
             
          READ(IU2,*,IOSTAT=IO) RG(I1:I2)
          IF(IO /= 0 ) GOTO 990
          READ(IU2,*,IOSTAT=IO) ZG(I1:I2)
          IF(IO /= 0 ) GOTO 990
        ENDDO 
        ENDDO
      PHI_PLANE = PHI_PLANE*PI/180.D0

C check grid points
      CALL CHECK_GRID_POINTS(IERR)
      IF(IERR/=0) GOTO 993

C magnetic field strength
      READ(IU3,*) BFSTREN

      RETURN
C Error 
990   WRITE(6,*)' ERROR IN READ_GRID IN MODULE GEOMETRY_PL'
      MESSAG = 'END OF FILE GRID_3D_DATA is reached'
      CALL STOP_ALL(MESSAG)
      RETURN
991   WRITE(6,*)' ERROR IN READ_GRID IN MODULE GEOMETRY_PL'
      CALL STOP_ALL(MESSAG)
      RETURN
992   WRITE(6,*)' ERROR IN READ_GRID IN MODULE GEOMETRY_PL'
      MESSAG = 'END OF FILE IU1 is reached'
      CALL STOP_ALL(MESSAG)
      RETURN
993   WRITE(6,*)' Program terminated !!!'                       
      MESSAG = '        '          
      CALL STOP_ALL(MESSAG)
      RETURN
      END SUBROUTINE READ_GEOMETRY_PL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SETUP_GEOMETRY_PL()
      INTEGER :: IZ,  I,   J,   K,   I1,   I2,   I3,   I4
     .,          IG
      REAL*8  :: R1, R2, Z1, Z2, FI1, FI2, AREA1, AREA2, DFL

      ALLOCATE ( VOL3D(0:GRID_P_OS(NZONET)-1) )
      ALLOCATE ( LCELL(0:GRID_P_OS(NZONET)-1) )
      DO 10 IZ=0,NZONET-1
      DO 10 I =0,ZON_RADI(IZ) - 1
      DO 10 J =0,ZON_POLO(IZ) - 1
         K  = 0         
         I1 = I + (J+K*SRF_POLO(IZ))*SRF_RADI(IZ)+GRID_P_OS(IZ)
         I2 = I1+ SRF_RADI(IZ)
         I3 = I2 + 1
         I4 = I1 + 1
         AREA1 = 0.5*ABS( (RG(I3)-RG(I2))*(ZG(I2)-ZG(I1))
     .-                   (RG(I2)-RG(I1))*(ZG(I3)-ZG(I2)) )
     .+          0.5*ABS( (RG(I4)-RG(I1))*(ZG(I3)-ZG(I4))
     .-                   (RG(I3)-RG(I4))*(ZG(I4)-ZG(I1)) )
         FI1 = PHI_PLANE(PHI_PL_OS(IZ)+K)
         R1  = 0.25*(RG(I1)+RG(I2)+RG(I3)+RG(I4))
         Z1  = 0.25*(ZG(I1)+ZG(I2)+ZG(I3)+ZG(I4))
      DO 10 K =0,ZON_TORO(IZ) - 1
         IG=I+(J+K*ZON_POLO(IZ))*ZON_RADI(IZ)+MESH_P_OS(IZ)

         I1 = I + (J+(K+1)*SRF_POLO(IZ))*SRF_RADI(IZ)+GRID_P_OS(IZ)
         I2 = I1+ SRF_RADI(IZ)
         I3 = I2 + 1
         I4 = I1 + 1
         AREA2 = 0.5*ABS( (RG(I3)-RG(I2))*(ZG(I2)-ZG(I1))
     .-                   (RG(I2)-RG(I1))*(ZG(I3)-ZG(I2)) )
     .+          0.5*ABS( (RG(I4)-RG(I1))*(ZG(I3)-ZG(I4))
     .-                   (RG(I3)-RG(I4))*(ZG(I4)-ZG(I1)) )

         FI2 = PHI_PLANE(PHI_PL_OS(IZ)+K+1)
         R2  = 0.25*(RG(I1)+RG(I2)+RG(I3)+RG(I4))
         VOL3D(IG) = 0.5*(R1+R2)*ABS(FI2-FI1)*0.5*(AREA1+AREA2)
         DFL       = 0.5*(R1+R2)*(FI2-FI1)
         LCELL(IG) = SQRT(DFL**2 + (R2-R1)**2 + (Z2-Z1)**2)

         AREA1 = AREA2
         FI1   = FI2
         R1    = R2
         Z1    = Z2
10    CONTINUE 
 
      RETURN 
      END SUBROUTINE SETUP_GEOMETRY_PL

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE RZ_REAL_COORDINATES(NZ0,JR0,JP0,RJ0,PJ0,TJ0,R,Z)

      INTEGER,INTENT(IN) :: NZ0,JR0,JP0
      REAL*8, INTENT(IN) :: RJ0,PJ0,TJ0 
      REAL*8, INTENT(OUT):: R,Z

      INTEGER :: JT0
     I,          I1,    I2,    I3,    I4
      REAL*8  :: AR, BR, CR, DR, AZ, BZ, CZ, DZ, R1, Z1, R2, Z2, XT
C---------------------------------------------------------------------
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

      FUNCTION LENGTH_FINE_CELL(IZ,I,J,K)
      IMPLICIT NONE
      REAL*8             :: LENGTH_FINE_CELL
      INTEGER,INTENT(IN) :: IZ,I,J,K
      INTEGER:: I1,I2,I3,I4
      REAL*8 :: R1,R2,FI1,FI2,Z1,Z2

      I1 = I + (J+K*SRF_POLO(IZ))*SRF_RADI(IZ)+GRID_P_OS(IZ)
      I2 = I1+ SRF_RADI(IZ)
      I3 = I2+ 1
      I4 = I1+ 1
      FI1= PHI_PLANE(PHI_PL_OS(IZ)+K)
      R1 = 0.25*(RG(I1)+RG(I2)+RG(I3)+RG(I4))
      Z1 = 0.25*(ZG(I1)+ZG(I2)+ZG(I3)+ZG(I4))

      I1 = I + (J+(K+1)*SRF_POLO(IZ))*SRF_RADI(IZ)+GRID_P_OS(IZ)
      I2 = I1+ SRF_RADI(IZ)
      I3 = I2 + 1
      I4 = I1 + 1
      FI2= PHI_PLANE(PHI_PL_OS(IZ)+K+1)
      R2 = 0.25*(RG(I1)+RG(I2)+RG(I3)+RG(I4))
      Z2 = 0.25*(ZG(I1)+ZG(I2)+ZG(I3)+ZG(I4))

      LENGTH_FINE_CELL=
     .SQRT((0.5*(R1+R2)*(FI2-FI1))**2+(R2-R1)**2+(Z2-Z1)**2)
      RETURN

      END FUNCTION LENGTH_FINE_CELL

      END MODULE GEOMETRY_PL
