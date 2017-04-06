!===============================================================================
module slicer
  use iso_fortran_env
  use emc3_grid
  implicit none
  private

  integer, parameter :: &
     UNDEFINED = 0, &
     ON_NODE   = 1, &
     ON_EDGE   = 2


  type, public :: t_polygon
     integer :: n = 0, i = 0
     real(real64), allocatable :: x(:,:)
!     integer,      allocatable :: itype(:)

     contains
     procedure :: new
     procedure :: add_node
     procedure :: insert_node
     procedure :: current_node
     procedure :: current_type
     procedure :: write
     procedure :: flip
  end type t_polygon

  public :: Zslice_cell

  contains
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine new(this, n)
  class(t_polygon)    :: this
  integer, intent(in) :: n

  this%n = n
  if (allocated(this%x)) deallocate(this%x)
!  if (allocated(this%itype)) deallocate(this%itype)
  if (n <= 0) return
  !allocate (this%x(0:n-1,2), this%itype(0:n-1))
  allocate (this%x(0:n-1,2))
  this%x = 1.2345d-6

  end subroutine new
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  !subroutine add_node(this, x, itype)
  subroutine add_node(this, x)
  class(t_polygon)         :: this
  real(real64), intent(in) :: x(2)
  !integer,      intent(in) :: itype

  integer :: i


  i = this%i

  if (i+1 > this%n) then
     write (6, *) 'error: cannot add node to polygon, number of nodes exceeded!'
     stop
  endif
  this%x(i,:)   = x
  !this%itype(i) = itype
  this%i = i + 1

  end subroutine add_node
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  subroutine insert_node(this, x)
  class(t_polygon)         :: this
  real(real64), intent(in) :: x(2)

  real(real64) :: x1(2), x2(2), p1(2), p2(2)
  integer :: i, i2, j, j2, n
  logical :: XX1, XX2


  ! last existing node
  n = this%i - 1


  ! scan through all possible new triangles (x1,x,x2)
  scan_loop: do i=0,n
     i2 = mod(i+1,n+1)
     x1 = this%x(i ,:)
     x2 = this%x(i2,:)

     ! check if v1 (x->x1) or v2 (x->x2) intersect any segments
     do j=0,n
        j2 = mod(j+1,n+1)
        p1 = this%x(j ,:)
        p2 = this%x(j2,:)

        !write (6, *) 'x  = ', x
        !write (6, *) 'x1 = ', x1
        !write (6, *) 'x2 = ', x2
        !write (6, *) 'p1 = ', p1
        !write (6, *) 'p2 = ', p2
        XX1 = .false.;  XX2 = .false.
        if (i  .ne. j  .and.  i  .ne. j2) XX1 = intersect(x, x1, p1, p2)
        if (i2 .ne. j  .and.  i2 .ne. j2) XX2 = intersect(x, x2, p1, p2)
        !write (6, *) i, j, XX1, XX2
        if (XX1 .or. XX2) cycle scan_loop
     enddo

     ! here is a good position to insert this point
!     write (6, *) 'position = ', i+1
!     do j=0,n
!        write (6, *) j, this%x(j,:)
!     enddo


     ! move all following points
     call this%add_node(x)
     do j=n,i+1,-1
        this%x(j+1,:) = this%x(j,:)
     enddo
     this%x(i+1,:) = x
     return
  enddo scan_loop

  write (6, *) 'error: cannot insert point to polygon!'
  write (6, *) 'x = ', x
  do j=0,n
     write (6, *) j, this%x(j,:)
  enddo
  stop
  contains
  !.....................................................................
  function intersect(x1, x2, p1, p2)
  real(real64), intent(in) :: x1(2), x2(2), p1(2), p2(2)
  logical                  :: intersect

  real(real64) :: det, dp(2), dx(2), b(2), xi(2)


  intersect = .false.
  dp = p2 - p1;  dx = x1 - x2
  det = dp(1)*dx(2) - dp(2)*dx(1)
  if (det == 0.d0) return

  b = x1 - p1
  xi(1) = 1.d0/det * ( dx(2)*b(1) - dx(1)*b(2))
  xi(2) = 1.d0/det * (-dp(2)*b(1) + dp(1)*b(2))
  !write (6, *) 'xi = ', xi
  if (xi(1) > 0.d0  .and.  xi(1) < 1.d0  .and.  xi(2) > 0.d0  .and.  xi(2) < 1.d0) intersect = .true.

  end function intersect
  !.....................................................................
  end subroutine insert_node
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  function current_node(this)
  class(t_polygon) :: this
  integer          :: current_node

  current_node = this%i + 1

  end function
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  function current_type(this)
  class(t_polygon) :: this
  integer          :: current_type

  integer :: i

  current_type = UNDEFINED
  i = this%i
  !if (i < this%n) current_type = this%itype(i)

  end function current_type
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  subroutine write(this, iu)
  class(t_polygon) :: this
  integer, intent(in) :: iu

  integer :: i


  do i=0,this%n
     write (iu, *) this%x(mod(i,this%n),:)
  enddo
  write (iu, *)

  end subroutine write
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  subroutine flip(this)
  class(t_polygon) :: this

  real(real64) :: x(0:this%n-1,2)
  integer :: i


  if (this%i .ne. this%n) then
     write (6, *) 'error: polygon has not been set up properly yet!'
     write (6, *) 'i = ', this%i
     write (6, *) 'n = ', this%n
     stop
  endif

  do i=0,this%n-1
     x(i,:) = this%x(this%n-1-i,:)
  enddo
  this%x = x

  end subroutine flip
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  subroutine Zslice_cell(ir, ip, it, iz, Z0, P)
  integer,         intent(in)  :: ir, ip, it, iz
  real(real64),    intent(in)  :: Z0
  type(t_polygon), intent(out) :: P

  integer, parameter :: &
     IMAP(0:7) = (/-1,  0,  1,  1,  1,  0, -1, -1/), &
     JMAP(0:7) = (/-1, -1, -1,  0,  1,  1,  1,  0/)

  real(real64) :: x(2), xi(3), xi1(3), xi2(3), v1(2), v2(2)
  real(real64) :: xcube(-1:1, -1:1, -1:1, 3)
  real(real64) :: phiA, phiB, R(8), R1, R2, Z(8), Z1, Z2, E(12,3), t
  !integer :: N(0:1, 0:1, 0:1), last_face, next_face
  integer :: Ncube(-1:1, -1:1, -1:1), Nc, nnode, nedge, jnode, jedge
  integer :: i, i1, i2, ii, j, j1, jj, k, k1, k2, ic, ig(8), m, mm


  phiA    = PHI_PLANE(it + PHI_PL_OS(iz))
  phiB    = PHI_PLANE(it+1 + PHI_PL_OS(iz))
  !ic      = ir + (ip + it*ZON_POLO(iz))*ZON_RADI(iz)  +  MESH_P_OS(iz)

  ig(1)   = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
  ig(2)   = ig(1) + 1
  ig(3)   = ig(2) + SRF_RADI(iz)
  ig(4)   = ig(1) + SRF_RADI(iz)
  ig(5:8) = ig(1:4) + SRF_POLO(iz)*SRF_RADI(iz)
  R       = RG(ig)
  Z       = ZG(ig)


  P%n = 0
  ! no intersection with Z0-plane
  if (minval(Z) > Z0  .or.  maxval(Z) < Z0) then
     return
  endif


  ! setup cell nodes and check if nodes are in Z0-plane ................
  Ncube = 0
  xcube = 0.d0
  do i=-1,1,2
  do j=-1,1,2
  do k=-1,1,2
     ic = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
     if (i == 1) ic = ic + 1
     if (j == 1) ic = ic + SRF_RADI(iz)
     if (k == 1) ic = ic + SRF_RADI(iz)*SRF_POLO(iz)

     xcube(i,j,k,1) = RG(ic)
     xcube(i,j,k,2) = ZG(ic)
     k2 = 0;  if (k == 1) k2 = 1
     xcube(i,j,k,3) = PHI_PLANE(it+k2 + PHI_PL_OS(iz))

     if (Z1 == ZG(ic)) Ncube(i,j,k) = 1
  enddo
  enddo
  enddo
  Nc = sum(Ncube)
  ! ....................................................................


  ! find intersections of edges with Z0-plane ..........................
  ! toroidal edges
  k = 0
  do i=-1,1,2
  do j=-1,1,2
     ! cell nodes at the end of this edge are already counted
     if (sum(Ncube(i,j,:)) == 2) cycle

     Z1 = xcube(i,j,-1,2)
     Z2 = xcube(i,j, 1,2)
     if (Z1 .ne. Z2) then
        t = (Z0 - Z1) / (Z2 - Z1)
        if (t > 0.d0  .and.  t < 1.d0) then
           xcube(i,j,0,:) = xcube(i,j,-1,:) + t * (xcube(i,j,1,:) - xcube(i,j,-1,:))
           Ncube(i,j,0)   = 1
        endif
     endif
  enddo
  enddo

  ! poloidal edges
  j = 0
  do i=-1,1,2
  do k=-1,1,2
     ! cell nodes at the end of this edge are already counted
     if (sum(Ncube(i,:,k)) == 2) cycle

     Z1 = xcube(i,-1,k,2)
     Z2 = xcube(i, 1,k,2)
     if (Z1 .ne. Z2) then
        t = (Z0 - Z1) / (Z2 - Z1)
        if (t > 0.d0  .and.  t < 1.d0) then
           xcube(i,0,k,:) = xcube(i,-1,k,:) + t * (xcube(i,1,k,:) - xcube(i,-1,k,:))
           Ncube(i,0,k)   = 1
        endif
     endif
  enddo
  enddo

  ! radial edges
  i = 0
  do j=-1,1,2
  do k=-1,1,2
     ! cell nodes at the end of this edge are already counted
     if (sum(Ncube(:,j,k)) == 2) cycle

     Z1 = xcube(-1,j,k,2)
     Z2 = xcube( 1,j,k,2)
     if (Z1 .ne. Z2) then
        t = (Z0 - Z1) / (Z2 - Z1)
        if (t > 0.d0  .and.  t < 1.d0) then
           xcube(0,j,k,:) = xcube(-1,j,k,:) + t * (xcube(1,j,k,:) - xcube(-1,j,k,:))
           Ncube(0,j,k)   = 1
        endif
     endif
  enddo
  enddo
  ! ....................................................................


  ! initialize polygon for intersection with Z0-plane
  m = sum(Ncube)
!  if (m > 0) then
!     write (89, *) R(1), phiA, Z(1), m, Nc
!  endif
!  return
!  if (m == 6) then
!     call plot_cell(ir,ip,it,iz, 91)
!     do i=-1,1;  do j=-1,1;  do k=-1,1
!     !if ((is_edge(i,j,k)  .or.  is_node(i,j,k))  .and.  Ncube(i,j,k) == 1) then
!     if (Ncube(i,j,k) == 1) then
!        write (92, *) xcube(i,j,k,:)
!     endif
!     enddo;  enddo;  enddo
!     stop
!  endif
!  return
  call P%new(m)


  ! only one cell node is touching Z0-plane
  if (m == 1) then
     return
  endif


  ! only one edge
  if (m == 2) then
!     do i=0,1
!     do j=0,1
!     do k=0,1
!        x(1) = xn(i,j,k,3);  x(2) = xn(i,j,k,1)
!        if (N(i,j,k) == 1) call P%add_node(x, ON_NODE)
!     enddo
!     enddo
!     enddo
!     if (P%current_node() .ne. 2) then
!     if (sum(N) .ne. m) then
     if (Nc .ne. m) then
        write (6, *) 'error: weird geometry N1E1!'
        stop
     endif

     return
  endif


  ! DIAGNOSTIC OUTPUT
  do i=-1,1;  do j=-1,1;  do k=-1,1
     if (Ncube(i,j,k) == 1) then
        x(1) = xcube(i,j,k,3);  x(2) = xcube(i,j,k,1)
        write (98, *) x
     endif
  enddo;  enddo;  enddo


  ! exactly 3 intersections with Z0-plane
  if (m == 3) then
  do i=-1,1;  do j=-1,1;  do k=-1,1
     if (Ncube(i,j,k) == 1) then
        x(1) = xcube(i,j,k,3);  x(2) = xcube(i,j,k,1)
        call P%add_node(x)
     endif
  enddo;  enddo;  enddo

  ! define orientation of polygon to be counter-clockwise (not strictly necessary)
  v1 = P%x(1,:) - P%x(0,:);  v2 = P%x(2,:) - P%x(0,:)
  if (v1(1)*v2(2) - v1(2)*v2(1) < 0.d0) call P%flip()

  ! done
  return
  endif


  ! more than 3 intersections with Z0-plane
  mm = 0
  ! first triangle
  loop3: do i=-1,1;  do j=-1,1;  do k=-1,1
  if ((is_node(i,j,k)  .or.  is_edge(i,j,k))  .and.  Ncube(i,j,k) == 1) then
     Ncube(i,j,k) = Ncube(i,j,k) + 1
     x(1) = xcube(i,j,k,3);  x(2) = xcube(i,j,k,1)
     call P%add_node(x)
     mm   = mm + 1
     if (mm == 3) exit loop3
  endif
  enddo;  enddo;  enddo loop3

  ! add mode triangles
  loopm: do i=-1,1;  do j=-1,1;  do k=-1,1
  if ((is_node(i,j,k)  .or.  is_edge(i,j,k))  .and.  Ncube(i,j,k) == 1) then
     x(1) = xcube(i,j,k,3);  x(2) = xcube(i,j,k,1)
     call P%insert_node(x)
     !call P%add_node(x)
     ! find segment where this point needs to be inserted


     !Ncube(i,j,k) = Ncube(i,j,k) + 1
     !x(1) = xcube(i,j,k,3);  x(2) = xcube(i,j,k,1)
     !call P%add_node(x)
     !mm   = mm + 1
     !if (mm == 3) exit loopm
  endif
  enddo;  enddo;  enddo loopm



  return









  ! set up faces from cell nodes and edges
  do i=-1,1;  do j=-1,1;  do k=-1,1
  if (is_node(i,j,k)  .or.  is_edge(i,j,k)) then
     Ncube(i,0,0) = Ncube(i,0,0) + Ncube(i,j,k)
     Ncube(0,j,0) = Ncube(0,j,0) + Ncube(i,j,k)
     Ncube(0,0,k) = Ncube(0,0,k) + Ncube(i,j,k)
     Ncube(0,0,0) = Ncube(0,0,0) + Ncube(i,j,k)
  endif
  enddo;  enddo;  enddo
  if (Ncube(0,0,0) .ne. m) then
     write (6, *) 'error: m \= Ncube(0,0,0), this should not happen!'
     stop
  endif

  ! select primary face
  i1 = 0;  j1 = 0;  k1 = 0;  m = 0
  do i=-1,1;  do j=-1,1;  do k=-1,1
  if (is_face(i,j,k)) then
     if (Ncube(i,j,k) > m) then
        m  = Ncube(i,j,k)
        i1 = i;  j1 = j;  k1 = k
     endif
  endif
  enddo;  enddo;  enddo


  ! this should not happen
  if (m > 4) then
     write (6, *) 'error: weird geometry!'
     write (6, *) Ncube
     stop
  endif


  ! all 4 cell nodes on this face are in Z0-plane
  if (m == 4) then
     ! double check that geometry is sane
     if (Ncube(0,0,0) > 4) then
        write (6, *) 'error: weird geometry!'
        write (6, *) Ncube
        stop
     endif
     Nc = 0
     do ii=-1,1,2;  do jj=-1,1,2
        call map_to_face(ii,jj, i1,j1,k1, i,j,k)
        Nc = Nc + Ncube(i,j,k)
     enddo;  enddo
     if (Nc .ne. 4) then
        write (6, *) 'error: weird geometry on face (',i,',',j,',',k,')!'
        stop
     endif

     do m=1,4
        ii = (-1)**((m+2)/2)
        jj = (-1)**((m+1)/2)
        call map_to_face(ii,jj, i1,j1,k1, i,j,k)
        x(1) = xcube(i,j,k,3);  x(2) = xcube(i,j,k,1)
        call P%add_node(x)
     enddo
     return
  endif


!  ! all intersections with Z0-plane are on this face
  if (m == 3  .and.  Ncube(0,0,0) == 3) then
     write (6, *) 'error: this case should have already been coverd!'
     stop
!     do m=1,8
!        ii = IMAP(m);  jj = JMAP(m)
!        call map_to_face(ii,jj, i1,j1,k1, i,j,k)
!        if (Ncube(i,j,k) == 1) then
!           x(1) = xcube(i,j,k,3);  x(2) = xcube(i,j,k,1)
!           call P%add_node(x)
!        endif
!     enddo
!     return
  endif


  ! 3 intersections with Z0-Plane on this face + additional intersection on other face
  ! orientation matters!
  if (m == 3  .and.  Ncube(0,0,0) > 3) then
     ! count cell nodes and edges intersecting with Z0-plane
     nnode = 0;  nedge = 0
     do mm=0,7
        ii = IMAP(mm);  jj = JMAP(mm)
        call map_to_face(ii,jj, i1,j1,k1, i,j,k)

        if (is_edge(i,j,k)) then
           if (Ncube(i,j,k) == 1) then
              nedge = nedge + 1
           else
              jedge = mm
           endif
        endif
        if (is_node(i,j,k)) then
           if (Ncube(i,j,k) == 1) then
              nnode = nnode + 1
           else
              jnode = mm
           endif
        endif
     enddo

     ! connect 3 nodes (jnode marks the missing node)
     if (nnode == 3) then
        do mm=jnode+1,jnode+7
           ii = IMAP(mod(mm,8));  jj = JMAP(mod(mm,8))
           call map_to_face(ii,jj, i1,j1,k1, i,j,k)
           if (is_node(i,j,k)) then
              x(1) = xcube(i,j,k,3);  x(2) = xcube(i,j,k,1)
              call P%add_node(x)
              Ncube(i,j,k) = 0

              ! move to next face...
              if (P%i == 3) then
                 call connect_face_at_edge(i1,j1,k1, IMAP(mod(mm+1,8)), JMAP(mod(mm+1,8)))
              endif
           endif
        enddo


     ! connect 3 intersection points on edges (jedge marks the missing edge)
     elseif (nedge == 3) then
        do mm=jedge+1,jedge+7
           ii = IMAP(mod(mm,8));  jj = JMAP(mod(mm,8))
           call map_to_face(ii,jj, i1,j1,k1, i,j,k)
           if (is_edge(i,j,k)) then
              x(1) = xcube(i,j,k,3);  x(2) = xcube(i,j,k,1)
              call P%add_node(x)
              Ncube(i,j,k) = 0

              ! move to next face...
              if (P%i == 3) then
                 call connect_face_at_edge(i1,j1,k1, ii, jj)
              endif
           endif
        enddo


     ! edgeX->node->edgeX
     elseif (nedge == 2  .and.  nnode == 1) then


     ! ???node->edgeX->node

     endif

  ! set up first segment of polygon from the two intersection points
  elseif (m == 2) then

  else
     write (6, *) 'error: this should not happen!'
     write (6, *) 'm = ', m
     stop
  endif



  ! search face and set up first node of polygon
  loop1: do ii=-1,1;  do jj=-1,1
     call map_to_face(ii,jj, i1,j1,k1, i,j,k)
     if ((is_node(i,j,k) .or. is_edge(i,j,k)) .and.  Ncube(i,j,k) == 1) then
        exit loop1
     endif
  enddo;  enddo loop1





!  ! start with cell node if possible
!  if (Nc > 0) then
!     outer_loop: do i=0,1
!     do j=0,1
!     do k=0,1
!        !x(1) = xn(i,j,k,3);  x(2) = xn(i,j,k,1)
!        if (N(i,j,k) == 1) then
!           xi1(1) = 1.d0 * i;  xi1(2) = 1.d0 * j;  xi1(3) = 1.d0 * k
!           !call P%add_node(x, ON_NODE)
!           call P%add_node(x)
!           N(i,j,k) = 0
!           exit outer_loop
!        endif
!     enddo
!     enddo
!     enddo outer_loop
!
!     ! set next_face
!
!
!  ! otherwise start on edge
!  else
!     do i=1,12
!        if (E(i,1) >= 0.d0) then
!           !call P%add_node(E(i,2:3), ON_EDGE)
!           call P%add_node(E(i,2:3))
!           E(i,1) = -1.d0
!           exit
!        endif
!     enddo
!  endif
!  last_face = -1
!
!
!  ! set up following nodes
!  do j=2,m
!     ! 1st segment
!     if (last_face == -1) then
!     if (P%current_node() == ON_NODE) then
!     else
!     endif
!     endif
!
!
!     ! find_next
!     !call P%add_node(next_node)
!  enddo


  contains
  !.....................................................................
  function is_node(i,j,k)
  integer :: i,j,k
  logical :: is_node

  is_node = .false.
  if (i*j*k /= 0) is_node = .true.

  end function is_node
  !.....................................................................
  function is_edge(i,j,k)
  integer :: i,j,k
  logical :: is_edge

  is_edge = .false.
  if (abs(i*j + j*k + k*i) == 1) is_edge = .true.

  end function is_edge
  !.....................................................................
  function is_face(i,j,k)
  integer :: i,j,k
  logical :: is_face

  is_face = .false.
  if (abs(i + j + k) == 1) is_face = .true.

  end function is_face
  !.....................................................................
  function face_side(i,j,k)
  integer :: i,j,k,face_side

  if (.not.is_face(i,j,k)) then
     write (6, *) 'error: face_side called with (',i,',',j,',',k,') which is not a face!'
     stop
  endif
  if (i .ne. 0) face_side = sign(1, i)
  if (j .ne. 0) face_side = sign(2, j)
  if (k .ne. 0) face_side = sign(3, k)

  end function face_side
  !.....................................................................
  subroutine map_to_face(ii,jj, i1,j1,k1, i,j,k)
  integer, intent(in)  :: ii,jj, i1,j1,k1
  integer, intent(out) :: i,j,k

  integer :: m

  if (.not.is_face(i1,j1,k1)) then
     write (6, *) 'error: face_side called with (',i1,',',j1,',',k1,') which is not a face!'
     stop
  endif
  m = face_side(i1,j1,k1)
  select case(abs(m))
  case(1)
     i = sign(1,m);  j = ii;  k = jj
  case(2)
     i = ii;  j = sign(1,m);  k = jj
  case(3)
     i = ii;  j = jj;  k = sign(1,m)
  end select

  end subroutine map_to_face
  !.....................................................................
  subroutine connect_face_at_edge(i1,j1,k1, ii,jj)
  integer, intent(inout) :: i1,j1,k1
  integer, intent(in)    :: ii,jj

  integer :: face


  if (.not.is_face(i1,j1,k1)) then
     write (6, *) 'error: connect_face_at_edge called with invalid (i1,j1,k1) = (',i1,',',j1,',',k1,')!'
     stop
  endif


  face = face_side(i1,j1,k1)
  i1 = 0;  j1 = 0;  k1 = 0
  select case(abs(face))
  case(1)
     j1 = ii
     k1 = jj

  case(2)
     i1 = ii
     k1 = jj

  case(3)
     i1 = ii
     j1 = jj
  end select

  end subroutine connect_face_at_edge
  !.....................................................................
  end subroutine Zslice_cell
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine plot_cell(ir, ip, it, iz, iu)
  integer, intent(in) :: ir, ip, it, iz, iu

  real(real64) :: phi(0:1)
  integer :: i1, i2, ig(0:7), j


  ig(0)   = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
  ig(1)   = ig(0) + 1
  ig(2)   = ig(1) + SRF_RADI(iz)
  ig(3)   = ig(0) + SRF_RADI(iz)
  ig(4:7) = ig(0:3) + SRF_POLO(iz)*SRF_RADI(iz)
  phi(0)  = PHI_PLANE(it   + PHI_PL_OS(iz))
  phi(1)  = PHI_PLANE(it+1 + PHI_PL_OS(iz))

  do j=0,1
     do i1=0,4
        i2 = mod(i1,4)
        write (iu, *) RG(ig(i2+j*4)), ZG(ig(i2+j*4)), phi(j)
     enddo
     write (iu, *)
  enddo
  do i1=0,3
     write (iu, *) RG(ig(i1)),   ZG(ig(i1)),   phi(0)
     write (iu, *) RG(ig(i1+4)), ZG(ig(i1+4)), phi(1)
     write (iu, *)
  enddo

  end subroutine plot_cell
  !---------------------------------------------------------------------

end module slicer
!===============================================================================




!===============================================================================
program slice_grid
  use iso_fortran_env
  use emc3_grid
  use slicer
  implicit none

  type(t_polygon) :: P
  real(real64)    :: Z
  integer         :: iz, ir, ip, it, ic


  call load_emc3_grid()


  Z = 0.d0
  do iz=0,NZONET-1
  do it=0,ZON_TORO(iz)-1
  do ip=0,ZON_POLO(iz)-1
  do ir=0,ZON_RADI(iz)-1
  !do ir=4,4
     !write (6, *)
     !write (6, *)
     !write (6, *)
     !write (6, *) ir, ip, it, iz
     call Zslice_cell(ir, ip, it, iz, Z, P)

     if (P%n > 2) call P%write(99)
  enddo
  enddo
  enddo
  enddo


end program slice_grid
!===============================================================================
