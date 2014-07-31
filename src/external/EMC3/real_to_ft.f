CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      NAME: REAL_TO_FT                                                C
C  FUNCTION: Coordinate transformation from real to flux-tube          C
C  Y.Feng 22.03.2010                                                   C
C Input: R,Z,PHI real coordinates                                      C
C     SITUATION=0: Nothing more known than (R,Z,PHI)                   C
C                  but,all NZ0,JR0,JP0,JT0,RJ0,PJ0,TJ0 are required    C
C               1: Nothing more known than (R,Z,PHI)                   C
C                  but,only NZ0,JR0,JP0,JT0 need to be estimated       C
C               2: initial NZ0,JR0,JP0,JT0,TJ0 roughly known           C
C                  exact NZ0,JR0,JP0,JT0,RJ0,PJ0,TJ0   are required    C
C OUTPUT: NZ0,JR0,JP0,JT0,RJ0,PJ0,TJ0                                  C
C   IERR=1: the point(R,Z,PHI) outside of the roroidal domain range    C
C        2: outside of the radial or poloidal range                    C
C        3: error id getting (RJ0,PJ0)                                 C
C        4: deviation > 1.E-6 cm                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE REAL_TO_FT
     .(R,Z,PHI,NZ0,JR0,JP0,JT0,RJ0,PJ0,TJ0,SITUATION,IERR)
      USE GEOMETRY_PL
      USE SURFACE_PL  

      IMPLICIT NONE
      REAL*8, INTENT(IN ):: R,Z,PHI
      INTEGER,INTENT(INOUT):: NZ0,JR0,JP0,JT0
      REAL*8, INTENT(INOUT):: RJ0,PJ0,TJ0
      INTEGER,INTENT(IN   ):: SITUATION
      INTEGER,INTENT(  OUT):: IERR

      INTEGER :: IZ,K1,K2,K_OFF,IG_OFF,WSP,IDR,IDP,ISTATU,I,ISURF,ITR
     .,          ISEARCH,STEPS,JT0_SAVE
      REAL*8  :: F,DIST_MIN,D,RP,ZP,JACOB,A,B,C,AR,BR,CR,DR,AZ,BZ,CZ,DZ
     .,          F1,F2,RS,ZS,TJ0_SAVE,F1_SAVE,F2_SAVE
      REAL*8,DIMENSION(:), ALLOCATABLE :: DL
      REAL*8,DIMENSION(5) :: R4,Z4 
      LOGICAL :: FOUND
      INTEGER,PARAMETER :: TRACE_MAX=1000

      STEPS = 2
      IF(SITUATION==2) STEPS=1
 
      QUICK_SLOW_SEARCH : DO ISEARCH=1,STEPS
      IERR= 0       
C Only necessary if NZ0,JR0,JP0,JT0 are unkown 
      IF(SITUATION==2) GOTO 15 

      JT0 = -1
      DIST_MIN = 1.E10
      DO 10 IZ=0,NZONET-1
        K_OFF = PHI_PL_OS(IZ)
        K1    = 0
        K2    = ZON_TORO(IZ)

        IF(    PHI< PHI_PLANE(K_OFF+K1) 
     .   .OR.  PHI> PHI_PLANE(K_OFF+K2) ) GOTO 10

        SEARCH_FOR_JT0 : DO
          JT0  = (K1+K2)/2
C In this way,JT0 always <  ZON_TORO(IZ)
          IF(K1+1==k2) THEN  
            EXIT SEARCH_FOR_JT0
          ELSEIF( PHI .LE. PHI_PLANE(K_OFF+JT0) ) THEN 
             K2 = JT0 
           ELSE
             K1 = JT0
           ENDIF
        ENDDO SEARCH_FOR_JT0 
     
        F   = (PHI-PHI_PLANE(K_OFF+JT0))
     .  /     (PHI_PLANE(K_OFF+JT0+1)-PHI_PLANE(K_OFF+JT0))
        TJ0       = JT0  +  F 
C JT0,TJ0 found
C--------------------------------------------------------------------
        F1        = 1.D0 - F
        F2        = F
       IF(ISEARCH==1) THEN 
       WSP = SRF_RADI(IZ)*SRF_POLO(IZ)
       IDP = (SRF_POLO(IZ)-2)/10
       ALLOCATE( DL(0:IDP) )
       IG_OFF=JT0*WSP+GRID_P_OS(IZ)
       ISURF = (SRF_RADI(IZ)-1)/2
       DO I=0,IDP         
        K1 = ISURF + I*10*SRF_RADI(IZ) + IG_OFF
        K2 = K1  + WSP
        DL(I)  = (F1*RG(K1) + F2*RG(K2) - R )**2   
     .+          (F1*ZG(K1) + F2*ZG(K2) - Z )**2
       ENDDO
 
       D = MINVAL(DL)
       IF(D < DIST_MIN) THEN 
         DIST_MIN = D
         NZ0 =  IZ
         JT0_SAVE = JT0
         JP0 =  (MINLOC(DL,1)-1)*10
         JR0 =  MIN(ISURF,SRF_RADI(IZ)-2)
         TJ0_SAVE = TJ0
         F1_SAVE  = F1
         F2_SAVE  = F2
       ENDIF
       DEALLOCATE(DL)
       ELSE ! ISEARCH==2
       WSP= SRF_RADI(IZ)*SRF_POLO(IZ)
       ALLOCATE( DL(0:WSP-1) )
       IG_OFF=JT0*WSP+GRID_P_OS(IZ)
       DL(0:WSP-1) = (R-   F2*RG(IG_OFF+WSP:IG_OFF+2*WSP-1)
     .-                    F1*RG(IG_OFF    :IG_OFF+WSP  -1))**2
     .+              (Z-   F2*ZG(IG_OFF+WSP:IG_OFF+2*WSP-1)
     .-                    F1*ZG(IG_OFF    :IG_OFF+WSP  -1))**2
 
       D = MINVAL(DL)
       IF(D < DIST_MIN) THEN 
         DIST_MIN = D
         NZ0 = IZ
         JT0_SAVE = JT0
         K1  = MINLOC(DL,1)-1
         JP0 = MIN(K1/SRF_RADI(IZ),ZON_POLO(IZ)-1)
         JR0 = MIN(K1-JP0*SRF_RADI(IZ),ZON_RADI(IZ)-1)
         TJ0_SAVE = TJ0
         F1_SAVE  = F1
         F2_SAVE  = F2
       ENDIF
       DEALLOCATE(DL)
       ENDIF ! ISEARCH  
10    CONTINUE 
C point outside of the toroidal range: IERR = 1
      IF(JT0 == -1) THEN 
        IERR = 1
        EXIT QUICK_SLOW_SEARCH
      ENDIF
C JT0,TJ0, NZ0 is determined
         JT0 =  JT0_SAVE
         TJ0 =  TJ0_SAVE
         F1  =  F1_SAVE
         F2  =  F2_SAVE
C JR0,JP0 estimated
15    CONTINUE 
      IF(SITUATION==2) THEN   
        F         = TJ0 - JT0
        F1        = 1.D0 - F
        F2        = F
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Find the termination cell: JR0,JP0
C Trace from cell center to (R,Z)
      IERR = 2 ! tracing stopped by a non-transparent surface

      TRACE_P : DO ITR=1,TRACE_MAX
      WSP   = SRF_RADI(NZ0)*SRF_POLO(NZ0)
      IG_OFF= JT0*WSP+GRID_P_OS(NZ0)

      K1    = IG_OFF+JP0*SRF_RADI(NZ0) + JR0
      K2    = K1+WSP
      R4(1) = F1*RG(K1  )+F2*RG(K2) 
      Z4(1) = F1*ZG(K1  )+F2*ZG(K2) 
      R4(4) = F1*RG(K1+1)+F2*RG(K2+1) 
      Z4(4) = F1*ZG(K1+1)+F2*ZG(K2+1) 

      K1    = K1 + SRF_RADI(NZ0)
      K2    = K1 + WSP
      R4(2) = F1*RG(K1  )+F2*RG(K2  ) 
      Z4(2) = F1*ZG(K1  )+F2*ZG(K2  ) 
      R4(3) = F1*RG(K1+1)+F2*RG(K2+1) 
      Z4(3) = F1*ZG(K1+1)+F2*ZG(K2+1) 

      RP   =SUM(R4(1:4))*0.25D0 - R
      ZP   =SUM(Z4(1:4))*0.25D0 - Z

      R4(5) = R4(1)
      Z4(5) = Z4(1)

      FOUND = .TRUE.
      DO 20 I=1,4
      JACOB = RP*(Z4(I+1)-Z4(I))-ZP*(R4(I+1)-R4(I))
      IF(JACOB==0.) GOTO 20

      A=( (R4(I)-R)*(Z4(I+1)-Z4(I))-(Z4(I)-Z)*(R4(I+1)-R4(I)) )/JACOB
      IF(A<0. .OR. A>1.) GOTO 20

      B=( (R4(I)-R)*ZP -(Z4(I)-Z)*RP )/JACOB
      IF(B<0. .OR. B>1.) GOTO 20

      FOUND = .FALSE.
      RJ0   = 0.5
      PJ0   = 0.5
      IF    (I==1) THEN
C Lower radial
        IF(JR0 == 0) THEN 
         ISURF =JR0+(JP0+JT0*ZON_POLO(NZ0))*SRF_RADI(NZ0)+NRS_OFF(NZ0)
         ISTATU=IDSURR(ISURF); IF(ISTATU<10) EXIT TRACE_P
         IDR = -1
         CALL R_SF_JUMP(ISTATU,IDR,NZ0,JR0,JP0,JT0,RJ0,PJ0,TJ0)
        ELSE
         JR0 = JR0 - 1   
        ENDIF
      ELSEIF(I==3) THEN
C Upper radial
        IF(JR0 == ZON_RADI(NZ0)-1) THEN 
         ISURF =JR0+1+(JP0+JT0*ZON_POLO(NZ0))*SRF_RADI(NZ0)+NRS_OFF(NZ0)
         ISTATU=IDSURR(ISURF); IF(ISTATU<10) EXIT TRACE_P
         IDR =  1
         CALL R_SF_JUMP(ISTATU,IDR,NZ0,JR0,JP0,JT0,RJ0,PJ0,TJ0)
        ELSE
         JR0 = JR0 + 1   
        ENDIF
      ELSEIF(I==4) THEN
C Lower poloidal 
        IF(JP0 == 0) THEN 
         ISURF =JR0+(JP0+JT0*SRF_POLO(NZ0))*ZON_RADI(NZ0)+NPS_OFF(NZ0)
         ISTATU=IDSURP(ISURF)
         IDP   = -1
         IF(ISTATU==1 .OR. ISTATU>=10) THEN
         CALL P_SF_JUMP(ISTATU,IDP,NZ0,JR0,JP0,JT0,RJ0,PJ0,TJ0)
         ELSE
         EXIT TRACE_P
         ENDIF
        ELSE
         JP0 = JP0 - 1
        ENDIF
      ELSEIF(I==2) THEN
C Upper poloidal 
        IF(JP0 == ZON_POLO(NZ0)-1) THEN 
         ISURF =JR0+(JP0+1+JT0*SRF_POLO(NZ0))*ZON_RADI(NZ0)+NPS_OFF(NZ0)
         ISTATU=IDSURP(ISURF)
         IDP   =  1
         IF(ISTATU==1 .OR. ISTATU>=10) THEN
         CALL P_SF_JUMP(ISTATU,IDP,NZ0,JR0,JP0,JT0,RJ0,PJ0,TJ0)
         ELSE
         EXIT TRACE_P
         ENDIF
        ELSE
         JP0 = JP0 + 1
        ENDIF
      ENDIF          
20    CONTINUE 
      IF(FOUND) THEN
        IERR = 0
        EXIT TRACE_P           
      ENDIF
      ENDDO TRACE_P
CC calculate RJ0,PJ0 if SITUATION/=1

      IF(IERR==0 .AND. SITUATION/=1 ) THEN
       AR  = RP
       BR  = 0.25*(R4(3)+R4(4)-R4(1)-R4(2))
       CR  = 0.25*(R4(3)+R4(2)-R4(1)-R4(4))
       DR  = 0.25*(R4(1)+R4(3)-R4(2)-R4(4))

       AZ  = ZP 
       BZ  = 0.25*(Z4(3)+Z4(4)-Z4(1)-Z4(2))
       CZ  = 0.25*(Z4(3)+Z4(2)-Z4(1)-Z4(4))
       DZ  = 0.25*(Z4(1)+Z4(3)-Z4(2)-Z4(4)) 

       a = BZ*DR-DZ*BR
       b = BZ*CR-CZ*BR+AZ*DR-DZ*AR
       c = AZ*CR-CZ*AR
       RJ0=-2.*C/B/(1.+SQRT(1.-4.*A*C/B**2))

       a = CZ*DR-DZ*CR
       b = CZ*BR-BZ*CR+AZ*DR-DZ*AR
       c = AZ*BR-BZ*AR
       PJ0=-2.*C/B/(1.+SQRT(1.-4.*A*C/B**2))

       IF(ABS(RJ0)>1. .AND. ABS(PJ0)>1. ) IERR = 3 
       call RZ_REAL_COORDINATES(NZ0,JR0,JP0,RJ0,PJ0,TJ0,RS,ZS)
       IF(SQRT( (R-RS)**2+(Z-ZS)**2 ) > 1.E-6 ) IERR = 4
         
      ENDIF
 
      IF(IERR==0) EXIT QUICK_SLOW_SEARCH
      ENDDO QUICK_SLOW_SEARCH

      RETURN

      END SUBROUTINE REAL_TO_FT 
