!===============================================================================
! Interpolate magnetic field
!===============================================================================
module interpolateB
  implicit none

  private

  integer, parameter :: n_max = 32

  integer       :: n_sets = 0
  character*120 :: filename(n_max)
  real*8        :: amplitude(n_max) = 1.d0, amplitude0 = 1.d0

  namelist /InterpolateB_Input/ &
     n_sets, filename, amplitude, amplitude0


  ! internal variables
  real*8, dimension(:,:), allocatable :: BA
  real*8  :: DELTA_PHI, DELTA_R, DELTA_Z, RBOX(3), ZBOX(3), BO(20,3), FIPER
  integer :: NF, NR, NZ, NS

  integer :: &
     LFALT = -1, &
     LRALT = -1, &
     LZALT = -1

  public :: &
     interpolateB_load, &
     interpolateB_broadcast, &
     interpolateB_get_Bf, &
     generate_field_from_coils, &
     write_ascii_data, &
     write_binary_data

  contains
!=======================================================================


!=======================================================================
  subroutine interpolateB_load (iu, iconfig)
  use run_control, only: Prefix
  use math
  integer, intent(in)  :: iu
  integer, intent(out) :: iconfig

  integer, parameter   :: iu1 = 70

  real*8, dimension(:,:), allocatable :: rtmp
  character*120 :: Data_File
  character*72  :: title, dummy
  integer :: i, io, j, ntmp(4), m


  ! read user configuration
  rewind (iu)
  read   (iu, InterpolateB_Input, end=1000)
  iconfig = 1
  write (6, *)
  write (6, 1001)


  ! read data
  do i=1,n_sets
     Data_File = trim(Prefix)//filename(i)
     open  (iu1, file=Data_File, iostat=io)
     !open  (iu1, file=Data_File, iostat=io, form='unformatted')
     if (io.ne.0) then
        write (6,2000) Data_File
        stop
     endif

     read  (iu1, '(a)') title
     write (6, 1002) title
     read  (iu1, *) DELTA_PHI, DELTA_R, DELTA_Z, RBOX(1), RBOX(3), RBOX(2), ZBOX(1), ZBOX(3), &
                 NF, NR, NS, NZ, ntmp
     ZBOX(2) = 0.5d0 * (ZBOX(1) + ZBOX(3))

     do j=1,4
        if (ntmp(j).gt.0) read (iu1, *) dummy
     enddo

     if (i == 1) then
        m     = NF * NR * NZ
        FIPER = pi2 / NS
        allocate (BA(8, 0:m-1), rtmp(8, 0:m-1))
        BA = 0.d0
        !write (6, 1003) RBOX(1), RBOX(3)
        !write (6, 1004) ZBOX(1), ZBOX(3)
        !write (6, 1005) NS
     else
        if (m .ne. NF * NR * NZ) then
           write (6, *) 'error: all data files must have the same resolution!'
           stop
        endif
     endif

     read  (iu1, *) rtmp
     BA = BA + amplitude(i)*rtmp
     close (iu1)
  enddo
  BA = BA * amplitude0


  deallocate (rtmp)
  return
 1000 iconfig = 0
 1001 format (3x,'- Pre-calculated field for interpolation:')
 1002 format (8x,a)
 1003 format (8x,'Rbox: ',2f10.3)
 1004 format (8x,'Zbox: ',2f10.3)
 1005 format (8x,'Symmetry: ', i8)
 2000 format ('error reading data file: ',a120)
  end subroutine interpolateB_load
!=======================================================================



!=======================================================================
  subroutine interpolateB_broadcast
  use parallel

  integer :: m


  if (nprs == 1) return

  call broadcast_real_s (DELTA_PHI)
  call broadcast_real_s (DELTA_R)
  call broadcast_real_s (DELTA_Z)
  call broadcast_real   (RBOX, 3)
  call broadcast_real   (ZBOX, 3)
  call broadcast_real_s (FIPER)

  call broadcast_inte_s (NF)
  call broadcast_inte_s (NR)
  call broadcast_inte_s (NZ)
  call broadcast_inte_s (NS)

  call wait_pe()
  m = NF * NR * NZ
  if (mype > 0) allocate (BA(8, 0:m-1))
  call broadcast_real   (BA, 8*m)

  end subroutine interpolateB_broadcast
!=======================================================================
 


!=======================================================================
  function interpolateB_get_Bf(r3) result(Bf)
  use math
  real*8, intent(in)  :: r3(3)
  real*8              :: Bf(3)

  REAL*8  :: R,Z,PHI
  LOGICAL :: OUTSIDE

  REAL*8  :: BI (8,8), FQ (2), RQ (2), ZQ (2), FR (4), FZ (4), FRZ (20)
  REAL*8  :: R1,R2,RQ1,ZQ1,FQ1,ALZ,ALR,ALF,F
  INTEGER :: L,I,J,K    
  INTEGER,DIMENSION(8) :: CORNE_P


  R   = r3(1)
  Z   = r3(2)
  PHI = r3(3)
  if (PHI.lt.0.d0) then
     K   = -PHI/pi2 + 1
     PHI = PHI + K*pi2
  endif

  Bf  = 0.d0
  OUTSIDE=R<RBOX(1) .OR. R>RBOX(3) .OR. Z<ZBOX(1) .OR. Z>ZBOX(3) 
  IF(OUTSIDE) RETURN                                

  K   = PHI/FIPER
  F   = PHI -K*FIPER
  ALR=(R-RBOX(1))/DELTA_R; ALZ=(Z-ZBOX(1))/DELTA_Z; ALF=F/DELTA_PHI      
  I  = ALR               ; J  = ALZ               ; K  = ALF
  RQ(2)=ALR-I            ; ZQ(2)=ALZ-J            ; FQ(2)=ALF-K  
  RQ(1)=1.-RQ(2)         ; ZQ(1)=1.-ZQ(2)         ; FQ(1)=1.-FQ(2)

  RQ1    = RQ(1)*RQ(2)*DELTA_R*0.5
  ZQ1    = ZQ(1)*ZQ(2)*DELTA_Z*0.5
  FQ1    = FQ(1)*FQ(2)*DELTA_PHI*0.5

  IF (K.NE.LFALT .OR. I.NE.LRALT .OR. J.NE.LZALT) THEN
     LFALT  = K    
     LRALT  = I    
     LZALT  = J

     R1     = I*DELTA_R + RBOX(1)
     R2     = R1   + DELTA_R
     CORNE_P(1) = K  + (I     + J   *NR)*NF
     CORNE_P(3) = K  + (I+1   + J   *NR)*NF
     CORNE_P(5) = K  + (I     +(J+1)*NR)*NF
     CORNE_P(7) = K  + (I+1   +(J+1)*NR)*NF

     K  = mod(K+1,NF)
     CORNE_P(2) = K  + (I     + J   *NR)*NF
     CORNE_P(4) = K  + (I+1   + J   *NR)*NF
     CORNE_P(6) = K  + (I     +(J+1)*NR)*NF
     CORNE_P(8) = K  + (I+1   +(J+1)*NR)*NF

     DO L = 1,8
        BI(1:8,L) = BA(L,CORNE_P(1:8))
     ENDDO                                  

     BO(1:8,1:3)=  BI(1:8,1:3)
     BO( 9,1)=-R1*(BI(1,6)+BI(1,8)-BI(2,6)-BI(2,8))-BI(1,2)+BI(2,2)
     BO( 9,2)= R1*(BI(1,4)-BI(2,4))                +BI(1,1)-BI(2,1)
     BO( 9,3)= R1*(BI(1,5)-BI(2,5))
     BO(10,1)=     BI(1,4)-BI(3,4)
     BO(10,2)=     BI(1,6)-BI(3,6)
     BO(10,3)=     BI(1,7)-BI(3,7)
     BO(11,1)=     BI(1,5)-BI(5,5)
     BO(11,2)=     BI(1,7)-BI(5,7)
     BO(11,3)=     BI(1,8)-BI(5,8)
     BO(12,1)=-R2*(BI(3,6)+BI(3,8)-BI(4,6)-BI(4,8))-BI(3,2)+BI(4,2)
     BO(12,2)= R2*(BI(3,4)-BI(4,4))                +BI(3,1)-BI(4,1)
     BO(12,3)= R2*(BI(3,5)-BI(4,5))
     BO(13,1)=     BI(2,4)-BI(4,4)
     BO(13,2)=     BI(2,6)-BI(4,6)
     BO(13,3)=     BI(2,7)-BI(4,7)
     BO(14,1)=     BI(2,5)-BI(6,5)
     BO(14,2)=     BI(2,7)-BI(6,7)
     BO(14,3)=     BI(2,8)-BI(6,8)
     BO(15,1)=-R1*(BI(5,6)+BI(5,8)-BI(6,6)-BI(6,8))-BI(5,2)+BI(6,2)
     BO(15,2)= R1*(BI(5,4)-BI(6,4))                +BI(5,1)-BI(6,1)
     BO(15,3)= R1*(BI(5,5)-BI(6,5))
     BO(16,1)=     BI(5,4)-BI(7,4)
     BO(16,2)=     BI(5,6)-BI(7,6)
     BO(16,3)=     BI(5,7)-BI(7,7)
     BO(17,1)=     BI(3,5)-BI(7,5)
     BO(17,2)=     BI(3,7)-BI(7,7)
     BO(17,3)=     BI(3,8)-BI(7,8)
     BO(18,1)=-R2*(BI(7,6)+BI(7,8)-BI(8,6)-BI(8,8))-BI(7,2)+BI(8,2)
     BO(18,2)= R2*(BI(7,4)-BI(8,4))                +BI(7,1)-BI(8,1)
     BO(18,3)= R2*(BI(7,5)-BI(8,5))
     BO(19,1)=     BI(6,4)-BI(8,4)
     BO(19,2)=     BI(6,6)-BI(8,6)
     BO(19,3)=     BI(6,7)-BI(8,7)
     BO(20,1)=     BI(4,5)-BI(8,5)
     BO(20,2)=     BI(4,7)-BI(8,7)
     BO(20,3)=     BI(4,8)-BI(8,8)
  ENDIF

  FR (1:2) = FQ (1:2) * RQ (1)
  FR (3:4) = FQ (1:2) * RQ (2)
  FZ (1:2) = FQ (1:2) * ZQ (1)
  FZ (3:4) = FQ (1:2) * ZQ (2)

  FRZ ( 1:4) = FR (1:4) * ZQ (1)
  FRZ ( 5:8) = FR (1:4) * ZQ (2)
  FRZ ( 9  ) = RQ (1) * ZQ (1) * FQ1
  FRZ (12  ) = RQ (2) * ZQ (1) * FQ1
  FRZ (15  ) = RQ (1) * ZQ (2) * FQ1
  FRZ (18  ) = RQ (2) * ZQ (2) * FQ1
  FRZ (10  ) = FZ (1) * RQ1
  FRZ (13  ) = FZ (2) * RQ1
  FRZ (16  ) = FZ (3) * RQ1
  FRZ (19  ) = FZ (4) * RQ1
  FRZ (11  ) = FR (1) * ZQ1
  FRZ (14  ) = FR (2) * ZQ1
  FRZ (17  ) = FR (3) * ZQ1
  FRZ (20  ) = FR (4) * ZQ1

  Bf(3)  =  sum( FRZ(1:20)*BO(1:20,1) )
  Bf(1)  =  sum( FRZ(1:20)*BO(1:20,2) )
  Bf(2)  =  sum( FRZ(1:20)*BO(1:20,3) )

  end function interpolateB_get_Bf
!=======================================================================



!=======================================================================
  subroutine generate_field_from_coils (N_sym, N_R, N_Z, N_phi, R_center, width, height)
  use polygones
  use math
  integer, intent(in) :: N_sym, N_R, N_Z, N_phi
  real*8, intent(in)  :: R_center, width, height

  real*8  :: Bf(3), r(3), CosPhi, SinPhi
  integer :: i, j, k, m, iba


  DELTA_PHI = pi2    / (N_sym * N_phi)
  DELTA_R   = width  / (N_R - 1)
  DELTA_Z   = height / (N_Z - 1)

  RBOX(1)   = R_center - 0.5d0*width
  RBOX(2)   = R_center
  RBOX(3)   = R_center + 0.5d0*width

  ZBOX(1)   = -0.5d0*height
  ZBOX(2)   = 0.d0
  ZBOX(3)   =  0.5d0*height

  NF        = N_phi
  NR        = N_R
  NS        = N_sym
  NZ        = N_Z

  m         = N_R * N_Z * N_phi
  allocate (BA(8, 0:m-1))
  BA = 0.d0
  do k=1,N_phi
     r(3)   = 1.d0 * (k-1) * DELTA_PHI
     CosPhi = cos(r(3))
     SinPhi = sin(r(3))
     write (6, *) r(3) / pi * 180.d0

     do j=1,N_Z
        r(2) = ZBOX(1) + (j-1) * DELTA_Z
        do i=1,N_R
           r(1) = RBOX(1) + (i-1) * DELTA_R
           Bf   = get_Bcyl_polygones(r)

           iba       = (k-1)  +  ((i-1) + (j-1)*N_R)*N_phi
           BA(1,iba) = Bf(3)
           BA(2,iba) = Bf(1)
           BA(3,iba) = Bf(2)
        enddo
     enddo
  enddo

  end subroutine generate_field_from_coils
!=======================================================================



!=======================================================================
  subroutine write_ascii_data (title, Output_File)
  character*72,  intent(in) :: title
  character*120, intent(in) :: Output_File

  integer, parameter :: iu = 42


  open  (iu, file=Output_File)
  write (iu, '(a72)') title
  write (iu, *) DELTA_PHI, DELTA_R, DELTA_Z, RBOX(1), RBOX(3), RBOX(2), ZBOX(1), ZBOX(3), &
                NF, NR, NS, NZ, 0, 0, 0, 0
  write (iu, *) BA
  close (iu)

  end subroutine write_ascii_data
!=======================================================================
  subroutine write_binary_data (title, Output_File)
  character*72,  intent(in) :: title
  character*120, intent(in) :: Output_File

  integer, parameter :: iu = 42


  open  (iu, file=Output_File, form='unformatted')
  write (iu) title
  write (iu) DELTA_PHI, DELTA_R, DELTA_Z, RBOX(1), RBOX(3), RBOX(2), ZBOX(1), ZBOX(3), &
             NF, NR, NS, NZ, 0, 0, 0, 0
  write (iu) BA
  close (iu)

  end subroutine write_binary_data
!=======================================================================

end module interpolateB
