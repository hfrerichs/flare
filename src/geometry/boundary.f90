!===============================================================================
! module:	boundary
!
! description:	Interface for boundaries (i.e. divertor targets, limiters, ...)
!               
!
! provides:
!	setup_boundary
!
!===============================================================================
module boundary
  use curve2D
  implicit none

  integer, parameter :: N_BNDRY_MAX        = 256, &
                        BNDRY_AXISYM_ELE = 1, &
                        BNDRY_BLOCK_LIM  = 2, &
                        BNDRY_TRI_ELE    = 3, &
                        BNDRY_QUAD_ELE   = 4

  character*120 :: boundary_file(N_BNDRY_MAX) = ''
  integer       :: boundary_type(N_BNDRY_MAX) = 0
  integer       :: n_boundary = 0

  namelist /Boundary_Input/ n_boundary, boundary_file, boundary_type

  type(t_curve), dimension(:), allocatable :: S_axi

  integer :: n_axi, n_block, n_tri, n_quad

  contains
!=======================================================================


!=======================================================================
  subroutine setup_boundary()
  use parallel

  ! load configuration on first processor
  if (mype == 0) call load_boundary()

  call wait_pe()
  !call broadcast_inte (i_config, BF_MAX_CONFIG)

  end subroutine setup_boundary
!=======================================================================


!=======================================================================
  subroutine load_boundary()
  use run_control
  use equilibrium

  integer, parameter :: iu = 24

  character*120 :: boundary_dir(3)
  character*80  :: header
  integer :: io, i, j, irun

  write (6,1000)
  write (6, *) 'Boundaries (divertor targets, limiters, vessel): '


  boundary_dir(1) = trim(Prefix)
  boundary_dir(2) = trim(Prefix)//trim(Boundary_sub_dir)//'/'
  boundary_dir(3) = './'

  ! 1: get number of components for each boundary type
  ! 2: read data
  do irun=1,2
     n_axi   = 0
     n_block = 0
     n_tri   = 0
     n_quad  = 0


     ! 0. check if boundary is provided by equilibrium
     if (equilibrium_provides_PFC()) then
        n_axi = n_axi + 1
        if (irun == 2) then
           call export_PFC(S_axi(n_axi))
           write (6, *)
           write (6, 3000) 'Axisymmetric surface (provided by equilibrium)'
           write (6, 3001) S_axi(n_axi)%n_seg
        endif
     endif


     ! 1. check base directory for input
     ! 2. check Boundary_sub_dir for input
     ! 3. check local directory for input (if not equivalent to base directory)
     do i=1,3
        if (i == 3  .and.  boundary_dir(1) == './') exit

        ! read namelist input
        open  (iu, file=trim(boundary_dir(i))//Boundary_input_file, status='old', iostat=io)
        if (io.ne.0) cycle
        read  (iu, Boundary_Input, end=2000)
        close (iu)
        boundary_file = trim(boundary_dir(i))//boundary_file

        ! check max. input number
        if (n_boundary > N_BNDRY_MAX) then
           write (6, *) 'error: n_boundary = ', n_boundary, ' > ', N_BNDRY_MAX
           stop
        endif

        ! sort input according to boundary type
        do j=1,n_boundary
           if (irun == 2) write (6, *)

           select case (boundary_type(j))
           ! axisymmetric surfaces
           case (BNDRY_AXISYM_ELE)
              n_axi   = n_axi   + 1
              if (irun == 2) then
                 call S_axi(n_axi)%read(boundary_file(j), report=.false., header=header)
                 if (header .ne. '') then
                    write (6, 3000) trim(header)
                 else
                    write (6, 3000) trim(boundary_file(j))
                 endif
                 write (6, 3001) S_axi(n_axi)%n_seg
              endif

           ! local block-limiters
           case (BNDRY_BLOCK_LIM)
              n_block = n_block + 1
              if (irun == 2) then
                    write (6, 3000) trim(boundary_file(j))
              endif

           ! mesh of triangular elements
           case (BNDRY_TRI_ELE)
              n_tri   = n_tri   + 1
              if (irun == 2) then
                    write (6, 3000) trim(boundary_file(j))
              endif

           ! mesh of quadrilateral elements
           case (BNDRY_QUAD_ELE)
              n_quad  = n_quad  + 1
              if (irun == 2) then
                    write (6, 3000) trim(boundary_file(j))
              endif

           case default
              write (6, *) 'boundary type ', boundary_type(j), ' not supported!'
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


  return
 1000 format (/ '========================================================================')
 2000 write (6, *) 'error while reading boundary input file!'
 3000 format ('   - ',a)
 3001 format (8x,'Number of segments: ',i8)
  stop
  end subroutine load_boundary
!=======================================================================



!=======================================================================
! check intersection (X) of trajectory r1->r2 with boundaries
!=======================================================================
  function intersect_boundary(r1, r2, X) result(l)
  real*8, intent(in)  :: r1(3), r2(3)
  real*8, intent(out) :: X(3)
  logical :: l

  real*8  :: phih, rz(2), t
  integer :: i


  l = .false.

  ! check intersection with axisymmetric surfaces
  do i=1,n_axi
     if (intersect_curve (r1(1:2), r2(1:2), S_axi(i), rz, t)) then
     !if (intersect_axisym_surf (r1, r2, S_axi(i), X)) then
        l      = .true.
        X(1:2) = rz
        X(3)   = r1(3) + t*(r2(3)-r1(3))
        return
     endif
  enddo

  end function intersect_boundary
!=======================================================================


end module boundary
