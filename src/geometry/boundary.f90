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
  use block_limiter, only: setup_block_limiter, center_bl, center_s, dist_p_oripocs, &
                           oripo_p, radius_s, normv_p, num_p, num_s, radius_bl, bl_filename, &
                           max_num_p, max_num_s, check_intersection, allocate_bl_arrays
  use quad_ele
  implicit none

  integer, parameter :: N_BNDRY_MAX        = 256, &
                        BNDRY_AXISYM_ELE = 1, &
                        BNDRY_BLOCK_LIM  = 2, &
                        BNDRY_TRI_ELE    = 3, &
                        BNDRY_QUAD_ELE   = 4, &
                        BNDRY_QUAD_ELE_STELLARATOR_SYM   = 40

  character*120 :: boundary_file(N_BNDRY_MAX) = ''
  integer       :: boundary_type(N_BNDRY_MAX) = 0
  integer       :: n_boundary = 0

  namelist /Boundary_Input/ n_boundary, boundary_file, boundary_type

  type(t_curve), dimension(:), allocatable :: S_axi
  type(t_quad_ele), dimension(:), allocatable :: S_quad

  integer :: n_axi, n_block, n_tri, n_quad

  contains
!=======================================================================


!=======================================================================
  subroutine setup_boundary()
  use parallel

  ! load configuration on first processor
  if (mype == 0) call load_boundary()
  call broadcast_boundary()

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

  ! irun = 1: get number of components for each boundary type
  !        2: read data
  do irun=1,2
     n_axi   = 0
     n_block = 0
     n_tri   = 0
     n_quad  = 0


     ! 0. check if boundary is provided by equilibrium
     if (equilibrium_provides_boundary()) then
        n_axi = n_axi + 1
        if (irun == 2) then
           call export_boundary(S_axi(n_axi))
           write (6, *)
           write (6, 3000) 'Axisymmetric surface (provided by equilibrium)'
           write (6, 3001) S_axi(n_axi)%n_seg
        endif
     endif


     ! i = 1. check base directory for input
     !     2. check Boundary_sub_dir for input
     !     3. check local directory for input (if not equivalent to base directory)
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
                 call S_axi(n_axi)%load(boundary_file(j), report=.false., header=header)
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
                 bl_filename = boundary_file(j)
                 call setup_block_limiter (bl_filename, n_block)
                    !write (6, 3000) trim(boundary_file(j))
              endif

           ! mesh of triangular elements
           case (BNDRY_TRI_ELE)
              n_tri   = n_tri   + 1
              if (irun == 2) then
                 if (header .ne. '') then
                    write (6, 3000) trim(header)
                 else
                    write (6, 3000) trim(boundary_file(j))
                 endif
              endif

           ! mesh of quadrilateral elements
           case (BNDRY_QUAD_ELE,BNDRY_QUAD_ELE_STELLARATOR_SYM)
              n_quad  = n_quad  + 1
              if (irun == 2) then
                 call S_quad(n_quad)%load(boundary_file(j), title=header)
                 if (header .ne. '') then
                    write (6, 3000) trim(header)
                 else
                    write (6, 3000) trim(boundary_file(j))
                 endif
              endif

           case default
              write (6, *) 'boundary type ', boundary_type(j), ' not supported!'
              stop
           end select

           ! add stellarator symmetric elements
           if (boundary_type(j) == BNDRY_QUAD_ELE_STELLARATOR_SYM) then
              n_quad  = n_quad  + 1
              if (irun == 2) then
                 S_quad(n_quad) = S_quad(n_quad-1)%get_stellarator_symmetric_element()
              endif
           endif
        enddo
     enddo


     ! allocate memory after 1st run
     if (irun == 1) then
        if (n_axi > 0)   allocate (S_axi(n_axi))
        if (n_quad > 0)  allocate (S_quad(n_quad))
        if (n_block > 0) call allocate_bl_arrays (n_block)
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
  subroutine broadcast_boundary()
  use parallel

  if (nprs == 1) return

  call wait_pe()
  call broadcast_axisym_surf()
  call broadcast_block_limiters()
  call broadcast_quad_ele()

  contains
!-----------------------------------------------------------------------
! axisymmetric surfaces
  subroutine broadcast_axisym_surf()

  integer :: i, n, m


  call broadcast_inte_s (n_axi)
  if (n_axi == 0) return

  if (mype > 0) allocate (S_axi(n_axi))
  do i=1,n_axi
     call broadcast_inte_s (S_axi(i)%n_seg)
     call broadcast_inte_s (S_axi(i)%n_dim)

     n = S_axi(i)%n_seg
     m = S_axi(i)%n_dim
     if (mype > 0) then
        allocate (S_axi(i)%x_data(0:n,m))
        !allocate (S_axi(i)%w_seg (1:n))
     endif
     call broadcast_real  (S_axi(i)%x_data, (n+1)*m)
     !call broadcast_real  (S_axi(i)%w_seg,   n)
  enddo

  end subroutine broadcast_axisym_surf
!-----------------------------------------------------------------------
! block limiters
  subroutine broadcast_block_limiters

  call broadcast_inte_s (n_block)
  if (n_block == 0) return

  if (mype > 0) call allocate_bl_arrays (n_block)
  call wait_pe()
  call broadcast_inte (num_p,          n_block            )
  call broadcast_real (oripo_p,        n_block*max_num_p*3)
  call broadcast_real (normv_p,        n_block*max_num_p*3)
  call broadcast_real (dist_p_oripocs, n_block*max_num_p  )
  call broadcast_real (center_bl,      n_block          *3)
  call broadcast_real (radius_bl,      n_block            )
  call broadcast_inte (num_s,          n_block            )
  call broadcast_real (radius_s,       n_block*max_num_s  )
  call broadcast_real (center_s,       n_block*max_num_s*3)

  end subroutine broadcast_block_limiters
!-----------------------------------------------------------------------
! quadrilateral elements
  subroutine broadcast_quad_ele

  integer :: i, n, m


  call broadcast_inte_s (n_quad)
  if (n_quad == 0) return

  if (mype > 0) allocate (S_quad(n_quad))
  do i=1,n_quad
     call broadcast_inte_s (S_quad(i)%n_phi)
     call broadcast_inte_s (S_quad(i)%n_RZ)
     call broadcast_inte_s (S_quad(i)%n_sym)

     n = S_quad(i)%n_phi
     m = S_quad(i)%n_RZ
     if (mype > 0) then
        allocate (S_quad(i)%phi(0:n))
        allocate (S_quad(i)%R(0:n,0:m))
        allocate (S_quad(i)%Z(0:n,0:m))
        allocate (S_quad(i)%cA(n,m,2))
        allocate (S_quad(i)%cB(n,m,2))
        allocate (S_quad(i)%cC(n,m,2))
        allocate (S_quad(i)%cD(n,m,2))
     endif
     call broadcast_real  (S_quad(i)%phi, n+1)
     call broadcast_real  (S_quad(i)%R,  (n+1)*(m+1))
     call broadcast_real  (S_quad(i)%Z,  (n+1)*(m+1))
     call broadcast_real  (S_quad(i)%cA,  n   * m   * 2)
     call broadcast_real  (S_quad(i)%cB,  n   * m   * 2)
     call broadcast_real  (S_quad(i)%cC,  n   * m   * 2)
     call broadcast_real  (S_quad(i)%cD,  n   * m   * 2)
  enddo

  end subroutine broadcast_quad_ele
!-----------------------------------------------------------------------
  end subroutine broadcast_boundary
!=======================================================================



!=======================================================================
! check intersection (X) of trajectory r1->r2 with boundaries
!=======================================================================
  function intersect_boundary(r1, r2, X) result(l)
  use math
  real*8, intent(in)  :: r1(3), r2(3)
  real*8, intent(out) :: X(3)
  logical :: l

  real*8  :: phih, rz(2), t, x1(3), x2(3), lambda_min
  integer :: i, i_ip


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


  ! check intersections with block limiter(s)
  call coord_trans (r1, CYLINDRICAL, x1, CARTESIAN)
  call coord_trans (r2, CYLINDRICAL, x2, CARTESIAN)
  do i=1,n_block
     if (check_intersection(x1, x2, X, lambda_min, i_ip, i)) then
        x1 = X
        call coord_trans (x1, CARTESIAN, X, CYLINDRICAL)
        l = .true.
        return
     endif
  enddo


  ! check intersection with mesh of quadrilateral elements
  do i=1,n_quad
     if (S_quad(i)%intersect(r1, r2, X)) then
        l = .true.
        return
     endif
  enddo


  end function intersect_boundary
!=======================================================================


end module boundary
