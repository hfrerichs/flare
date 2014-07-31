!===============================================================================
! module:	pfc
!
! description:	Provide setup of plasma facing components (PFC)
!
! provides:
!	setup_pfc
!
!===============================================================================
module pfc
  use curve2D
  implicit none

  integer, parameter :: N_PFC_MAX        = 256, &
                        N_PFC_AXISYM_ELE = 1, &
                        N_PFC_BLOCK_LIM  = 2, &
                        N_PFC_TRI_ELE    = 3, &
                        N_PFC_QUAD_ELE   = 4

  character*120 :: pfc_file(N_PFC_MAX) = ''
  integer       :: pfc_type(N_PFC_MAX) = 0
  integer       :: n_pfc = 0

  namelist /PFC_Input/ n_pfc, pfc_file, pfc_type

  type(t_curve), dimension(:), allocatable :: S_axi

  integer :: n_axi, n_block, n_tri, n_quad

  contains
!=======================================================================


!=======================================================================
  subroutine setup_pfc()
  use parallel

  ! load configuration on first processor
  if (mype == 0) call load_pfc()

  call wait_pe()
  !call broadcast_inte (i_config, BF_MAX_CONFIG)

  end subroutine setup_pfc
!=======================================================================


!=======================================================================
  subroutine load_pfc()
  use run_control
!  use pfc_axisym
  use equilibrium

  integer, parameter :: iu = 24

!  type (t_pfc_axisym) :: T
  character*120 :: pfc_dir(3)
  character*80  :: header
  integer :: io, i, j, irun

  write (6,1000)
  write (6, *) 'Plasma facing components (PFC) input: '


  pfc_dir(1) = trim(Prefix)
  pfc_dir(2) = trim(Prefix)//trim(PFC_sub_dir)//'/'
  pfc_dir(3) = './'

  ! 1: get number of pfc components of each type
  ! 2: read data
  do irun=1,2
     n_axi   = 0
     n_block = 0
     n_tri   = 0
     n_quad  = 0


     ! 0. check if pfc is provided by equilibrium
     if (equilibrium_provides_PFC()) then
        n_axi = n_axi + 1
        if (irun == 2) then
           call export_PFC(S_axi(n_axi))
           write (6, *)
           write (6, 3000) 'Provided by equilibrium'
           write (6, 3001) S_axi(n_axi)%n_seg
        endif
     endif


     ! 1. check base directory for input
     ! 2. check PFC_sub_dir for input
     ! 3. check local directory for input (if not equivalent to base directory)
     do i=1,3
        if (i == 3  .and.  pfc_dir(1) == './') exit

        ! read namelist input
        open  (iu, file=trim(pfc_dir(i))//PFC_input_file, status='old', iostat=io)
        if (io.ne.0) cycle
        read  (iu, PFC_Input, end=2000)
        close (iu)
        pfc_file = trim(pfc_dir(i))//pfc_file

        ! check max. input number
        if (n_pfc > N_PFC_MAX) then
           write (6, *) 'error: n_pfc = ', n_pfc, ' > ', N_PFC_MAX
           stop
        endif

        ! sort input according to pfc type
        do j=1,n_pfc
           if (irun == 2) write (6, *)

           select case (pfc_type(j))
           ! axisymmetric surfaces
           case (N_PFC_AXISYM_ELE)
              n_axi   = n_axi   + 1
              if (irun == 2) then
                 call S_axi(n_axi)%read(pfc_file(j), report=.false., header=header)
                 if (header .ne. '') then
                    write (6, 3000) trim(header)
                 else
                    write (6, 3000) trim(pfc_file(j))
                 endif
                 write (6, 3001) S_axi(n_axi)%n_seg
              endif

           ! local block-limiters
           case (N_PFC_BLOCK_LIM)
              n_block = n_block + 1
              if (irun == 2) then
                    write (6, 3000) trim(pfc_file(j))
              endif

           ! mesh of triangular elements
           case (N_PFC_TRI_ELE)
              n_tri   = n_tri   + 1
              if (irun == 2) then
                    write (6, 3000) trim(pfc_file(j))
              endif

           ! mesh of quadrilateral elements
           case (N_PFC_QUAD_ELE)
              n_quad  = n_quad  + 1
              if (irun == 2) then
                    write (6, 3000) trim(pfc_file(j))
              endif

           case default
              write (6, *) 'pfc type ', pfc_type(j), ' not supported!'
              stop
           end select
        enddo
     enddo


     ! allocate memory after 1st run
     if (irun == 1) then
        if (n_axi > 0) allocate (S_axi(n_axi))
        !if (n_block > 0)
     endif
  enddo







!  ! 1a. check base directory for input
!  pfc_dir = trim(Prefix)
!  open  (iu, file=trim(pfc_dir)//PFC_input_file, status='old', iostat=io)
!
!  ! 1b. check PFC_sub_dir for input
!  pfc_dir = trim(Prefix)//trim(PFC_sub_dir)//'/'
!  if (io.ne.0) then
!     open  (iu, file=trim(pfc_dir)//PFC_input_file, status='old', iostat=io)
!  endif
!
!  ! 1c. check local directory for input (if not equivalent to base directory)
!  do i=1,3
!     if (i == 3  .and.  dir(1) == './') exit
!  enddo
!
!  ! 2. read PFC configuration
!  if (io == 0) then
!     read  (iu, PFC_Input, end=2000)
!     close (iu)
!     pfc_file = trim(pfc_dir)//pfc_file
!  else
!     write (6, *)
!     write (6, *) 'no PFC input provided'
!  endif
!  if (n_pfc > N_PFC_MAX) then
!     write (6, *) 'error: n_pfc = ', n_pfc, ' > ', N_PFC_MAX
!     stop
!  endif


  ! 3. sort input according to pfc type
  !call T%init()

  return
 1000 format (/ '========================================================================')
 2000 write (6, *) 'error while reading PFC input file!'
 3000 format ('   - ',a)
 3001 format (8x,'Number of segments: ',i8)
  stop
  end subroutine load_pfc
!=======================================================================



!=======================================================================
! check intersection (X) of trajectory r1->r2 with plasma facing components
!=======================================================================
  function intersect_pfc(r1, r2, X) result(l)
  real*8, intent(in)  :: r1(3), r2(3)
  real*8, intent(out) :: X(3)
  logical :: l

  real*8  :: phih, rz(2), t
  integer :: i


  l = .false.

  ! check intersection with axisymmetric surfaces
  do i=1,n_axi
     if (intersect_curve (r1, r2, S_axi(i), rz, t)) then
        l      = .true.
        X(1:2) = rz
        X(3)   = r1(3) + t*(r2(3)-r1(3))
        return
     endif
  enddo

  end function intersect_pfc
!=======================================================================


end module pfc
