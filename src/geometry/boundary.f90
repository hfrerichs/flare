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
  use iso_fortran_env
  use dataset
  use curve2D
  use block_limiter, only: setup_block_limiter, center_bl, center_s, dist_p_oripocs, &
                           oripo_p, radius_s, normv_p, num_p, num_s, radius_bl, bl_filename, &
                           max_num_p, max_num_s, check_intersection, allocate_bl_arrays, bl_outside
  use quad_ele
  implicit none

  integer, parameter :: N_BNDRY_MAX        = 256, &
                        BNDRY_AXISYM_ELE = 1, &
                        BNDRY_BLOCK_LIM  = 2, &
                        BNDRY_TRI_ELE    = 3, &
                        BNDRY_QUAD_ELE   = 4, &
                        BNDRY_QUAD_ELE_STELLARATOR_SYM   = 40

!  type t_boundary
!     ! ... data pointer ...
!
!     ! boundary type (AXISYM_ELE, BLOCK_ELE, TRI_ELE, QUAD_ELE)
!     integer :: itype
!     character(len=80) :: label
!  end type t_boundary


  character*120 :: boundary_file(N_BNDRY_MAX) = ''
  integer       :: boundary_type(N_BNDRY_MAX) = 0

  ! type 1 and 4 surfaces can be used to model a limiter. Then boundary_side = -1 defines
  ! the outside of this volume (MUST BE A CLOSED CONTOUR IN R-Z PLANE) as inside the
  ! plasma domain
  integer       :: boundary_side(N_BNDRY_MAX) = 1
  integer       :: n_boundary = 0

  namelist /Boundary_Input/ n_boundary, boundary_file, boundary_type, boundary_side

  type(t_curve), dimension(:), allocatable :: S_axi
  type(t_quad_ele), dimension(:), allocatable :: S_quad
  character(len=80), dimension(:), allocatable :: L_axi, L_block, L_tri, L_quad

  integer, dimension(:), allocatable :: elem_os

  integer :: n_axi, n_block, n_tri, n_quad, n_boundary1, n_elem

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
  use run_control, only: Prefix, Boundary_Prefix
  use equilibrium
  include '../config.h'

  integer, parameter :: iu = 24

  character*120 :: boundary_dir(4)
  character*80  :: header
  integer :: io, i, j, irun, n

  write (6,1000)
  write (6, *) 'Boundaries (divertor targets, limiters, vessel): '


  boundary_dir(1) = trim(Boundary_Prefix)
  boundary_dir(2) = trim(Prefix)
  boundary_dir(3) = trim(Prefix)//trim(Boundary_sub_dir)//'/'
  boundary_dir(4) = './'

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
           call S_axi(n_axi)%closed_check()
           write (6, *)
           write (6, 3000) 'Axisymmetric surface (provided by equilibrium)'
           write (6, 3001) S_axi(n_axi)%n_seg
        endif
     endif


     ! i = 1. check user defined boundary directory for input
     !     2. check base directory (magnetic configuration) for input
     !     3. check Boundary_sub_dir for input
     !     4. check local directory for input (if not equivalent to base directory)
     do i=1,4
        if (boundary_dir(i) == '') cycle
        if (i == 4  .and.  boundary_dir(2) == './') exit

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
                 call S_axi(n_axi)%load(boundary_file(j), output=SILENT, header=header)
                 if (header .ne. '') then
                    write (6, 3000) trim(header)
                    L_axi(n_axi) = trim(header)
                 else
                    write (6, 3000) trim(boundary_file(j))
                    L_axi(n_axi) = trim(boundary_file(j))
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
                 L_block(n_block) = trim(boundary_file(j))
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
                    L_quad(n_quad) = trim(header)
                 else
                    write (6, 3000) trim(boundary_file(j))
                    L_quad(n_quad) = trim(boundary_file(j))
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
        if (n_axi > 0)   allocate (S_axi(n_axi), L_axi(n_axi))
        if (n_quad > 0)  allocate (S_quad(n_quad), L_quad(n_quad))
        if (n_block > 0) then
           call allocate_bl_arrays (n_block)
           allocate (L_block(n_block));  L_block = ''
        endif
     endif
  enddo


  ! setup offset for boundary elements
  ! actual number of boundaries
  ! (might be larger than n_boundary for stellarator symmetric elements)
  n_boundary1 = n_axi + n_quad + n_tri + n_block
  allocate (elem_os(0:n_boundary1))
  elem_os = 0
  do i=1,n_axi
     elem_os(i) = elem_os(i-1) + S_axi(i)%n_seg
  enddo
  do i=1,n_block
     elem_os(n_axi+i) = elem_os(n_axi+i-1) + 0
  enddo
  do i=1,n_quad
     n = S_quad(i)%n_phi * S_quad(i)%n_RZ
     elem_os(n_axi+n_block+i) = elem_os(n_axi+n_block+i-1) + n
  enddo
  n_elem = elem_os(n_boundary1)

  write (6, *) 'elem_os = '
  do i=0,n_boundary1
     write (6, *) elem_os(i)
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
!     call broadcast_inte_s (S_axi(i)%n_seg)
!     call broadcast_inte_s (S_axi(i)%n_dim)
!
!     n = S_axi(i)%n_seg
!     m = S_axi(i)%n_dim
!     if (mype > 0) then
!        allocate (S_axi(i)%x(0:n,m))
!        !allocate (S_axi(i)%w_seg (1:n))
!     endif
!     call broadcast_real  (S_axi(i)%x_data, (n+1)*m)
!     !call broadcast_real  (S_axi(i)%w_seg,   n)
     call S_axi(i)%broadcast()
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
! check intersection (X) of trajectory r1->r2 with axisymmetric boundaries
!=======================================================================
! optional output:
!    X:      coordinates of intersection with boundary
!    id:     boundary number
!    ielem:  element number on boundary
!    tau:    relative coordinate along trajectory r1->r2
!=======================================================================
  function intersect_axisymmetric_boundary(r1, r2, X, id, ielem, tau) result(l)
  use math
  real*8, intent(in)   :: r1(2), r2(2)
  real*8, intent(out)  :: X(2)
  integer, intent(out) :: id
  integer, intent(out), optional :: ielem
  real(real64), intent(out), optional :: tau
  logical :: l

  real*8  :: phih, rz(2), t, x1(3), x2(3), lambda_min
  integer :: i, ish


  l  = .false.
  id = 0
  if (present(ielem)) ielem = -1

  ! check intersection with axisymmetric surfaces
  do i=1,n_axi
     if (intersect_curve (r1, r2, S_axi(i), rz, t, ish=ish)) then
     !if (intersect_axisym_surf (r1, r2, S_axi(i), X)) then
        l      = .true.
        X(1:2) = rz
        id     = i
        if (present(ielem)) ielem  = ish
        if (present(tau))   tau    = t
        return
     endif
  enddo

  end function intersect_axisymmetric_boundary
!=======================================================================



!=======================================================================
! check intersection (X) of trajectory r1->r2 with boundaries
!=======================================================================
! optional output:
!    X:      coordinates of intersection with boundary
!    id:     boundary number
!    ielem:  element number on boundary
!    tau:    relative coordinate along trajectory r1->r2
!=======================================================================
  function intersect_boundary(r1, r2, X, id, ielem, tau) result(l)
  use math
  real*8, intent(in)   :: r1(3), r2(3)
  real*8, intent(out)  :: X(3)
  integer, intent(out) :: id
  integer, intent(out), optional :: ielem
  real(real64), intent(out), optional :: tau
  logical :: l

  real*8  :: phih, rz(2), t, x1(3), x2(3), lambda_min
  integer :: i, i_ip, ish


  l  = .false.
  id = 0
  if (present(ielem)) ielem = -1

  ! check intersection with axisymmetric surfaces
  do i=1,n_axi
     if (intersect_curve (r1(1:2), r2(1:2), S_axi(i), rz, t, ish=ish)) then
     !if (intersect_axisym_surf (r1, r2, S_axi(i), X)) then
        l      = .true.
        X(1:2) = rz
        X(3)   = r1(3) + t*(r2(3)-r1(3))
        id     = i
        if (present(ielem)) ielem  = ish
        if (present(tau))   tau    = t
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
        l  = .true.
        id = n_axi + i
        if (present(tau))   tau    = lambda_min
        return
     endif
  enddo


  ! check intersection with mesh of quadrilateral elements
  do i=1,n_quad
     if (S_quad(i)%intersect(r1, r2, X, ish, tau)) then
        l  = .true.
        id = n_axi + n_block + i
        if (ish < 0) then
           write (6, *) i, id
           write (6, *) r1
           write (6, *) r2
           write (6, *) X
           write (6, *) ish
           stop
        endif
        if (present(ielem)) ielem  = ish
        return
     endif
  enddo


  end function intersect_boundary
!=======================================================================



!=======================================================================
  subroutine surface_plot (n, W)
  integer, intent(in) :: n
  real(real64), dimension(0:n_elem-1,n), intent(in) :: W

  integer, parameter :: iu = 42


  real(real64) :: Wsurf(n)
  integer :: is, i1, i2, j, ielem, ielem0, ix, iy


  do is=1,n_boundary1
     i1    = elem_os(is-1)
     i2    = elem_os(is)-1

     write (6, *) 'surface ', is, ': element ', i1, '->', i2
     do j=1,n
        Wsurf(j) = sum(W(i1:i2,j))
     enddo
     write (6, *) 'total = ', Wsurf

     if (sum(Wsurf) > 0.d0) then
        open (iu, file="surface_plot.txt")
        ! check surface type
        do ielem=i1,i2
           ielem0 = ielem - i1

           ix = ielem0 / S_quad(1)%n_RZ
           iy = ielem0 - ix * S_quad(1)%n_RZ

           if (iy == 0) write(iu, *)
           write (iu, *) ix, iy, W(ielem,:)
        enddo
        close (iu)
     endif
  enddo

  end subroutine surface_plot
!=======================================================================



!=======================================================================
! check if point r is outside of configuration boundary
!=======================================================================
  function outside_boundary(r)
  use iso_fortran_env
  real(real64), intent(in) :: r(3)
  logical                  :: outside_boundary

  type(t_curve) :: C
  logical       :: side
  real(real64)  :: phi
  integer       :: l


  ! 0. initialize
  outside_boundary = .false.


  ! 1. check axisymmetric (L2-type) surfaces
  do l=1,n_axi
     side = boundary_side(l) == 1
     if (S_axi(l)%outside(r(1:2)) .eqv. side) then
        outside_boundary = .true.
        return
     endif
  enddo


  ! 2. check block limiters (CSG-type)
  do l=1,n_block
     if (bl_outside(l, r)) then
        outside_boundary = .true.
        return
     endif
  enddo


  ! 3.c heck Q4-type boundaries
  phi = r(3)
  do l=1,n_quad
     side = boundary_side(n_axi + n_block + l) == 1
     C    = S_quad(l)%slice(phi)
     if (C%outside(r(1:2)) .eqv. side) then
        outside_boundary = .true.
        return
     endif
  enddo

  end function outside_boundary
!=======================================================================



!=======================================================================
! BOUNDARY_SLICE: return outline of boundary iboundary in phi-plane
!=======================================================================
  function boundary_slice(iboundary, phi) result(S)
  integer,      intent(in) :: iboundary
  real(real64), intent(in) :: phi
  type(t_curve)            :: S

  integer :: i


  if (iboundary <= n_axi) then
     S = S_axi(iboundary)

  elseif (iboundary <= n_axi + n_block) then
     write (6, *) 'error: slicing of block-limiter not yet implemented!'

  elseif (iboundary <= n_axi + n_block + n_quad) then
     i = iboundary - n_axi - n_block
     S = S_quad(i)%slice(phi)

  else
     write (6, *) 'error: boundary id ', iboundary, ' out of range!'
     stop
  endif

  end function boundary_slice
!=======================================================================



!=======================================================================
  function boundary_in_zone(iboundary, phi1, phi2)
  integer,      intent(in) :: iboundary
  real(real64), intent(in) :: phi1, phi2
  logical                  :: boundary_in_zone

  integer :: i


  if (iboundary <= n_axi) then
     boundary_in_zone = .true.

  elseif (iboundary <= n_axi + n_block) then
     ! TODO: check location of block limiter
     boundary_in_zone = .true.

  elseif (iboundary <= n_axi + n_block + n_quad) then
     i = iboundary - n_axi - n_block
     if (S_quad(i)%phi(0) >= phi1 .and. S_quad(i)%phi(S_quad(i)%n_phi) <= phi2) then
        boundary_in_zone = .true.
     else
        boundary_in_zone = .false.
     endif

  else
     write (6, *) 'error: boundary id ', iboundary, ' out of range!'
     stop
  endif
  end function boundary_in_zone
!=======================================================================



!=======================================================================
  function boundary_label(iboundary) result(L)
  integer, intent(in) :: iboundary
  character(len=80)   :: L

  integer :: i


  if (iboundary <= n_axi) then
     L = L_axi(iboundary)

  elseif (iboundary <= n_axi + n_block) then
     i = iboundary - n_axi
     L = L_block(i)

  elseif (iboundary <= n_axi + n_block + n_quad) then
     i = iboundary - n_axi - n_block
     L = L_quad(i)

  else
     write (6, *) 'error: boundary id ', iboundary, ' out of range!'
     stop
  endif

  end function boundary_label
!=======================================================================

end module boundary
