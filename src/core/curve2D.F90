!===============================================================================
! module:	curve2D
!
! description:	Operations on 2D curves (line segments)
!===============================================================================
module curve2D
  use iso_fortran_env
  use dataset

#if defined(FLARE)
  use math
  use parallel
  implicit none
#else
  implicit none

  real(real64), parameter :: &
     pi   = 3.14159265358979323846264338328d0, &
     pi2  = 2.d0 * pi
#endif
  private


  integer, parameter, public :: &
     ANGLE     = 1, &
     DISTANCE  = 2


  type, public :: t_curve
     ! number of line segments (n_seg), number of coordinates (n_dim)
     integer  :: n_seg = -1, n_dim = 0

     ! length of curve
     real(real64) :: l = 0.d0

     ! nodes along the curve, dimension(0:n_seg, 1:n_dim)
     type (t_dataset) :: nodes
     real(real64), dimension(:,:), pointer :: x => null()

     ! relative weight for each line segment, dimension(1:n_seg)
     !real(real64), dimension(:),   pointer :: w_seg => null()
     ! relative coordinate along curve (0:n_seg): 0->1
     real(real64), dimension(:),   pointer :: w => null()

     ! treat curve as closed loop?
     logical :: closed = .false.

     contains

     procedure :: load
     procedure :: closed_check
     procedure :: duplicate_node_check
     procedure :: new
     procedure :: destroy
     procedure :: copy
     procedure :: broadcast
     procedure :: plot
     procedure :: sort_loop
     procedure :: sort_by_angle
     procedure :: sort_by_distance
     procedure :: expand
     procedure :: flip
     procedure :: left_hand_shift
     procedure :: get_distance_to
     procedure :: setup_angular_sampling
     procedure :: setup_length_sampling
     procedure :: setup_length_sampling_curvature_weighted
     procedure :: setup_segment_sampling
     procedure :: setup_coordinate_sampling
     procedure :: sample_at
     procedure :: split3, splitn
     procedure :: split3seg, splitnseg
     procedure :: length
     procedure :: area
     procedure :: outside
     procedure :: intersect_curve => t_curve_intersect_curve
  end type t_curve

  type(t_curve), public, parameter :: Empty_curve = t_curve(0,0,0.d0,Empty_dataset,null(),null())


  public :: intersect_curve, intersect_lines, make_2D_curve, connect
  public :: SILENT, VERBOSE

  contains
!=======================================================================



!=======================================================================
  subroutine load(this, data_file, output, header)
  class (t_curve),  intent(inout)         :: this
  character(len=*), intent(in)            :: data_file
  integer,          intent(in),  optional :: output
  character(len=*), intent(out), optional :: header

  real(real64) :: dl, x1(2), x2(2)


  call this%nodes%load(data_file, 2, output, header, -1)
  this%n_seg =  this%nodes%nrow-1
  this%n_dim =  2
  this%x     => this%nodes%x
  call this%closed_check()
  call this%duplicate_node_check()

  end subroutine load
!=======================================================================



!=======================================================================
  subroutine closed_check(this)
  class (t_curve),  intent(inout)         :: this

  real(real64) :: dl, x1(2), x2(2)


  if (this%n_seg < 0) return
  ! check if curve is closed
  x1 = this%x(0,:)
  x2 = this%x(this%n_seg,:)
  dl = sqrt(sum((x1-x2)**2))
  if (dl < epsilon(real(1.0,real64))) then
     this%closed          = .true.
     this%x(0,:)          = 0.5d0 * (x1+x2)
     this%x(this%n_seg,:) = 0.5d0 * (x1+x2)
  endif

  end subroutine closed_check
!=======================================================================



!=======================================================================
  subroutine duplicate_node_check(this)
  class (t_curve),  intent(inout)         :: this

  real(real64), dimension(:,:), allocatable :: xtmp
  real(real64) :: dl, x1(2), x2(2)
  integer      :: idrop(this%n_seg), ndrop, i, j


  ndrop = 0
  ! find duplicate nodes
  do i=1,this%n_seg
     x1 = this%x(i-1,:)
     x2 = this%x(i  ,:)
     dl = sqrt(sum((x1-x2)**2))
     if (dl < epsilon(real(1.0,real64))) then
        ndrop        = ndrop + 1
        idrop(ndrop) = i
     endif
  enddo


  ! remove duplicate nodes
  if (ndrop == 0) return
  allocate (xtmp(0:this%n_seg-ndrop,2))
  j     = 0
  ndrop = 1
  do i=0,this%n_seg
     if (i == idrop(ndrop)) then
        write (6, *) 'removing duplicate node at ', this%x(i,:)
        ndrop = ndrop + 1
        cycle
     endif

     xtmp(j,:) = this%x(i,:)
     j         = j+1
  enddo

  ! copy xtmp
  call this%new(j-1)
  this%x = xtmp


  ! cleanup
  deallocate (xtmp)

  end subroutine duplicate_node_check
!=======================================================================



!=======================================================================
  subroutine new(this, n_seg, sampling)
  class(t_curve)      ::  this
  integer, intent(in) :: n_seg
  logical, intent(in), optional :: sampling


  this%n_seg = n_seg
  this%n_dim = 2
  call this%nodes%new(n_seg+1, 2, -1)
  this%x     => this%nodes%x

  if (present(sampling) .and. sampling) then
     if (associated(this%w)) deallocate(this%w)
     allocate (this%w(0:n_seg))
  endif

  end subroutine new
!=======================================================================



!=======================================================================
  subroutine destroy(this)
  class(t_curve)      ::  this

  call this%nodes%destroy()
  nullify(this%x)
  !if (associated(this%w_seg)) deallocate(this%w_seg)
  if (associated(this%w)) deallocate(this%w)

  end subroutine destroy
!=======================================================================



!=======================================================================
  subroutine copy(this, C)
  class(t_curve) ::  this
  type(t_curve)  :: C

  integer :: n_seg


  n_seg = C%n_seg
  call this%new(n_seg)
  this%x = C%x
  this%closed = C%closed

  end subroutine copy
!=======================================================================



!=======================================================================
  subroutine broadcast(this)
  class(t_curve) :: this

  integer :: n, m


#if defined(FLARE)
  call broadcast_inte_s (this%n_seg)
  call broadcast_inte_s (this%n_dim)
  n = this%n_seg
  m = this%n_dim

  if (mype > 0) then
     allocate (this%nodes%x(0:n,m))
     this%x => this%nodes%x
  endif
  call broadcast_real  (this%nodes%x, (n+1)*m)

#endif

  end subroutine broadcast
!=======================================================================



!=======================================================================
! Generate 2D curve from two 1D arrays (xarr, yarr) of size n
!=======================================================================
  subroutine make_2D_curve (n, xarr, yarr, C)
  integer,       intent(in)  :: n
  real(real64),  intent(in)  :: xarr(n), yarr(n)
  type(t_curve), intent(out) :: C


  call C%new(n-1)
  C%x(:,1) = xarr
  C%x(:,2) = yarr

  end subroutine make_2D_curve
!=======================================================================



!=======================================================================
! INTERSECT_CURVE
!
! calculate intersection of segment x1->x2 with curve C.
!
! output (optional):
! xh:	intersection point
! th:   position on segment x1->x2
! sh:   position on segment ish of curve C
! ish:  segment number of C where intersection point is located
!
! optional input:
! intersect_mode = 0: check for intersection with segment x1->x2 only (default)
!                  1: check for intersection also beyond x2
!                 -1: check for intersection along the full line defined by x1->x2
! is_start:           start search for intersection at this segment
!===============================================================================
  function intersect_curve (x1, x2, C, xh, th, sh, ish, intersect_mode, is_start)
  real(real64),  intent(in)            :: x1(2), x2(2)
  type(t_curve), intent(in)            :: C
  real(real64),  intent(out), optional :: xh(2), th, sh
  integer,       intent(out), optional :: ish
  integer,       intent(in ), optional :: intersect_mode, is_start
  logical                              :: intersect_curve

  real(real64) :: t, s, xl1(2), xl2(2), xh0(2), th0
  integer      :: i, is, is0, mode


  intersect_curve = .false.
  t   = 0.d0
  s   = 0.d0
  th0 = 1.d99

  mode = 0
  if (present(intersect_mode)) mode = intersect_mode

  is0  = 0
  if (present(is_start)) is0 = is_start

  do i=0,C%n_seg-1
     is  = mod(i + is0, C%n_seg)
     xl1 = C%x(is  ,:)
     xl2 = C%x(is+1,:)
     if (intersect_lines(x1,x2,xl1,xl2,t,s,xh0)) then
        ! intersection with actual segment on curve
        if (s.ge.0.d0 .and. s.le.1.d0) then

           if ((mode == -1)  .or.  &
               (mode ==  0  .and.  t.ge.0.d0 .and. t.le.1.d0)  .or.  &
               (mode ==  1  .and.  t.ge.0.d0)) then

              intersect_curve = .true.
              if (abs(t) < abs(th0)) then
              th0 = t
              if (present(xh)) then
                 xh = xh0
              endif
              if (present(th)) then
                 th = t
              endif
              if (present(sh)) then
                 sh = s
              endif
              if (present(ish)) then
                 ish = is+1
              endif
              endif
           endif
        endif
     endif
  enddo

  end function intersect_curve
!=======================================================================



!=======================================================================
! INTERSECT_LINES
!
! calculate intersection of two lines, each defined by two points
! L:	line through L1 and L2
! M:	line through M1 and M2
!
! output:
! l, m:	relative position on L and M, respectively
! X:    intersection point (optional)
!=======================================================================
  function intersect_lines (L1, L2, M1, M2, l, m, X)
  real(real64), intent(in)  :: L1(2), L2(2), M1(2), M2(2)
  real(real64), intent(out) :: l, m
  real(real64), optional    :: X(2)
  logical                   :: intersect_lines

  real(real64) :: e1(2), e2(2), d1, d2, edet, dx(2), n0(2), b


  intersect_lines = .false.
  l   = 0.d0
  m   = 0.d0

  e1   = L2 - L1
  d1   = dsqrt(sum(e1**2))
  e1   = e1 / d1
  e2   = M2 - M1
  d2   = dsqrt(sum(e2**2))
  if (d2 == 0.d0) return
  e2   = e2 / d2
  edet = - e1(1)*e2(2) + e1(2)*e2(1)
  dx   = M1 - L1


  ! check if lines are parallel
  if (edet.eq.0.d0) then
     n0(1) = -e1(2)
     n0(2) =  e1(1)
     b = sum(n0*dx)

     ! check if lines are identical
     if (b.eq.0.d0) then
        intersect_lines = .true.
        l = sum(dx*e1) / d1
        m = sum(-dx*e2) / d2
        if (present(X)) then
           X = L + l * e1 * d1
        endif
     endif

     return
  endif


  ! lines are not parallel
  intersect_lines = .true.
  l    = - e2(2)*dx(1) + e2(1)*dx(2)
  l    = l / (edet * d1)
  m    = - e1(2)*dx(1) + e1(1)*dx(2)
  m    = m / (edet * d2)

  if (present(X)) then
     X = L1 + l * e1 * d1
  endif

  end function intersect_lines
!=======================================================================



!=======================================================================
! Calculate (first) intersection with curve C
! output:
!     x      intersection point
!     tau    relative coordinate from first node
!=======================================================================
  function t_curve_intersect_curve(this, C, x, tau)
  class(t_curve)            :: this
  type(t_curve), intent(in) :: C
  real(real64), intent(out) :: x(2), tau
  logical                   :: t_curve_intersect_curve

  real(real64) :: x1(2), x2(2)
  integer :: is


  t_curve_intersect_curve = .false.
  x   =  0.d0
  tau =  0.d0
  do is=1,this%n_seg
     x1 = this%x(is-1, :)
     x2 = this%x(is  , :)
     if (intersect_curve(x1, x2, C, x)) then
        tau = tau + dsqrt(sum((x-x1)**2))
        tau = tau / this%length()
        t_curve_intersect_curve = .true.
        exit
     endif

     tau = tau + sqrt(sum((x2-x1)**2))
  enddo


  end function t_curve_intersect_curve
!=======================================================================



!=======================================================================
  subroutine plot(this, iu, filename, append, length)
  class(t_curve) :: this
  integer,          intent(in), optional :: iu
  character(len=*), intent(in), optional :: filename
  logical,          intent(in), optional :: append, length

  real(real64) :: l
  integer      :: i, iu0


  ! set default unit number for output
  iu0 = 50

  ! Unit number given for output?
  if (present(iu)) iu0 = iu

  ! Output_File given?
  if (present(filename)) then
     if (present(append)  .and.  append) then
        open  (iu0, file=filename, position='append')
        write (iu0, *)
     else
        open  (iu0, file=filename)
     endif
  endif


  ! write data
  l = 0.d0
  do i=0,this%n_seg
     if (present(length)  .and.  length) then
        if (i>0) l = l + sqrt(sum((this%x(i,:)-this%x(i-1,:))**2))
        write (iu0, *) this%x(i,:), l
     elseif (associated(this%w)) then
        write (iu0, *) this%x(i,:), this%w(i)
     else
        write (iu0, *) this%x(i,:)
     endif
  enddo


  ! Output_File given?
  if (present(filename)) then
     close (iu0)
  endif

  end subroutine plot
!=======================================================================



!=======================================================================
! sort points on a closed curve 'C' with regard to center/reference point 'xc'
! and direction 'd'
!=======================================================================
  subroutine sort_loop (C, xc, d, method)
  implicit none

  class(t_curve), intent(inout)        :: C
  real(real64),   intent(in), optional :: xc(2), d(2)
  integer,        intent(in), optional :: method

  integer :: sort_method


  sort_method = ANGLE
  if (present(method)) sort_method = method


  select case (sort_method)
  case(ANGLE)
     call C%sort_by_angle(xc, d)
  case(DISTANCE)
     call C%sort_by_distance(xc)
  case default
     write (6, *) 'error in t_curve%sort_loop: method ', sort_method, ' undefined!'
     stop
  end select

  end subroutine sort_loop
!=======================================================================
  subroutine sort_by_angle (C, xc, d)
  implicit none

  class(t_curve), intent(inout)        :: C
  real(real64),   intent(in), optional :: xc(2), d(2)

  real(real64), dimension(:,:), allocatable :: tmp_arr
  real(real64) :: xr(2), dr(2), xmax, xmin, ymax, ymin
  real(real64) :: theta_0, theta, dxl(2), x_0(2), tmp(3), t
  integer      :: i, j, n


  n = C%n_seg
  if (n <= 0) return


  ! check if reference point x_c is given
  if (present(xc)) then
     xr = xc
  ! set default values otherwise
  else
     ! calculate characteristic parameters
     xmin = minval(C%x(:,1))
     xmax = maxval(C%x(:,1))
     ymin = minval(C%x(:,2))
     ymax = maxval(C%x(:,2))

     xr(1) = 0.5d0 * (xmin + xmax)
     xr(2) = 0.5d0 * (ymin + ymax)
  endif


  ! check if reference direction is given
  if (present(d)) then
     dr    = d
  ! set default values otherwise
  else
     dr(1) = 1.d0
     dr(2) = 0.d0
  endif


  ! allocate memory for temporary array
  allocate (tmp_arr(0:n, 3))
  tmp_arr(:,1:2) = C%x

  ! calculate angles w.r.t. xr
  theta_0  = datan2(dr(2), dr(1))
  do i=0,n
     dxl   = C%x(i,:) - xr
     theta = datan2(dxl(2), dxl(1))
     if (theta .lt. theta_0) theta = theta + pi2
     tmp_arr(i,3) = theta
  enddo

  ! sort x,y data with respect to theta values
  do 100 j=1,n
  do 100 i=0,j
     if (tmp_arr(i,3) .gt. tmp_arr(j,3)) then
        tmp          = tmp_arr(i,:)
        tmp_arr(i,:) = tmp_arr(j,:)
        tmp_arr(j,:) = tmp
     endif
  100 continue


  ! close loop
  t    = (tmp_arr(0,3) - tmp_arr(n,3) + pi2)
  if (t.eq.0.d0) then
     ! loop already closed
     C%x(0:n,:) = tmp_arr(0:n,1:2)
  else
     t    = (theta_0      - tmp_arr(n,3) + pi2) / t
     x_0  = tmp_arr(n,1:2)
     dxl  = tmp_arr(0,1:2) - x_0
     x_0  = x_0 + t * dxl

     ! copy data to output variable L
     if (x_0(1).eq.tmp_arr(0,1) .and. x_0(2).eq.tmp_arr(0,2)) then
        call C%new(n+1)
        C%x(0:  n,:) = tmp_arr(0:n,1:2)
        C%x(  n+1,:) = x_0
     else
        call C%new(n+2)
        C%x(0    ,:) = x_0
        C%x(1:n+1,:) = tmp_arr(0:n,1:2)
        C%x(  n+2,:) = x_0
     endif
  endif
  C%closed = .true.
  call C%setup_angular_sampling(xc)

  end subroutine sort_by_angle
!=======================================================================
  subroutine sort_by_distance (C, xc)
  implicit none

  class(t_curve), intent(inout)        :: C
  real(real64),   intent(in), optional :: xc(2)


  real(real64), dimension(:,:), allocatable :: x
  integer,      dimension(:),   allocatable :: imark, jmark
  real(real64) :: xr(2), x0(2), t, v(2)
  integer      :: i, inb, n


  if (present(xc)) then
     xr = xc
  else
     xr(1) = 0.5d0 * (minval(C%x(:,1)) + maxval(C%x(:,1)))
     xr(2) = 0.5d0 * (minval(C%x(:,2)) + maxval(C%x(:,2)))
  endif


  n = C%n_seg
  allocate (imark(0:n), jmark(0:n), x(0:n,2))
  imark = -1
  jmark = -1
  x     = 0.d0

  inb = find_first(0.d0)
  v(1) = 0.d0
  v(2) = 1.d0
  imark(0)   = inb
  jmark(inb) = 0
  !write (97, *) C%x(inb,:)
  x(0,:) = C%x(inb,:)
  do i=0,n-1
     !write (6, *) i
     call find_closest_neighbor(imark(i), v, inb)
     imark(i+1) = inb
     jmark(inb) = i

     !write (97, *) C%x(inb,:)
     x(i+1,:) = C%x(inb,:)
  enddo

  ! close loop
  t  = (xr(2) - x(0,2)) / (x(n,2) - x(0,2))
  x0 = x(0,:) + t * (x(n,:) - x(0,:))
  x(0,:) = x0
  x(n,:) = x0
  C%closed = .true.

  C%x = x



  ! check orientation
  if (x(1,2)-x(0,2) > 0) then
!     write (6, *) 1
  else
!     write (6, *) -1
     call C%flip()
  endif

  deallocate (imark, jmark, x)
  call C%setup_length_sampling()

  contains
  !.....................................................................
  subroutine find_closest_neighbor(i0, v, inb)
  integer, intent(in)  :: i0
  real(real64), intent(inout) :: v(2)
  integer, intent(out) :: inb

  real(real64) :: x0(2), d, dmin, dmin_fallback, vnew(2)
  integer      :: j, inb_fallback


  x0   = C%x(i0,:)
  dmin = 1.d99
  dmin_fallback = 1.d99
  inb  = -1
  do j=0,n
     ! only check unmarked elements
     if (jmark(j) == -1) then
        vnew = x0-C%x(j,:)
        d = sqrt(sum(vnew**2))
        if (d < dmin  .and.  sum(vnew*v) > 0.d0) then
           dmin = d
           inb  = j
        elseif (d < dmin_fallback) then
           dmin_fallback = d
           inb_fallback  = j
        endif
     endif
  enddo
  v = x0-C%x(inb,:)

  if (inb == -1) inb = inb_fallback

  if (inb == -1) then
     write (6, *) 'error: inb = -1'
     stop
  endif
  !write (6, *) i0, inb, dmin

  end subroutine find_closest_neighbor
  !.....................................................................
  function find_first(theta0) result(i0)

  real(real64) :: x(2), theta, theta0, dtheta
  integer :: j, i0

  dtheta = 1.d99
  do j=0,n
     x = C%x(j,:)
     theta = atan2(x(2)-xr(2), x(1)-xr(1))
     if (abs(theta-theta0) < dtheta  .and.  theta > 0) then
        dtheta = abs(theta-theta0)
        i0     = j
     endif
  enddo

  end function find_first
  end subroutine sort_by_distance
!=======================================================================



!=======================================================================
! Expand curve/loop by scaling the distance to reference point x0
!=======================================================================
  subroutine expand(this, x0, dl)
  class(t_curve)           :: this
  real(real64), intent(in) :: x0(2), dl

  real(real64) :: d(2), d0
  integer :: i


  do i=0,this%n_seg
     d           = this%x(i,:)-x0
     d0          = sqrt(sum(d**2))
     this%x(i,:) = x0 + (d0+dl)/d0 * d
  enddo

  end subroutine expand
!=======================================================================



!=======================================================================
! flip orientation of curve
!=======================================================================
  subroutine flip(this)
  class(t_curve)           :: this

  real(real64) :: xtmp(this%n_dim)
  integer      :: i, j


  do i=0,this%n_seg/2
     j = this%n_seg - i

     xtmp        = this%x(i,:)
     this%x(i,:) = this%x(j,:)
     this%x(j,:) = xtmp
  enddo

  end subroutine flip
!=======================================================================



!=======================================================================
! Shift each line segment in direction of its left-hand side normal vector and
! re-connect nodes.
!
!                     ^ normal vector
!                     |
! segment: ====================>
!
!=======================================================================
  recursive subroutine left_hand_shift(this, dl)
#ifdef DEBUG
  use string
#endif
  class(t_curve)           :: this
  real(real64), intent(in) :: dl

  real(real64), parameter  :: emin = 1.d-8

  integer, save :: IDEBUG = 0

  real(real64), dimension(:,:,:), allocatable :: x_new
  real(real64), dimension(:,:),   allocatable :: ts_new, x_tmp
  type(t_curve) :: Ctmp
  real(real64)  :: el(2), en(2), x11(2), x12(2), x21(2), x22(2), xh(2), t, s, d
  integer       :: k, n, k2, kmax, n2, i_remove


! 1. initialize local variables
  n = this%n_seg
  if (this%n_dim .ne. 2) then
     write (6, *) 'error in subroutine t_curve%left_hand_shift: n_dim = 2 expected!'
     stop
  endif
  call Ctmp%copy(this)
#ifdef DEBUG
  write (6, *) 'left_hand_shift', IDEBUG
#endif


! 2. working array x_new(i,j,k): raw nodes positions
  ! i:	coordinate
  ! j:	lower (1) and upper (2) node
  ! k:	segment number
  allocate (x_new(2,2,0:n-1))
  do k=0,n-1
     ! direction of line segment
     el = Ctmp%x(k+1,:) - Ctmp%x(k,:)
     d  = dsqrt(sum(el**2))
     if (d == 0.d0) then
        write (6, *) 'error in t_curve%left_hand_shift: duplicate nodes!'
        stop
     endif
     el = el / d

     ! left-hand normal vector
     en(1) = - el(2)
     en(2) =   el(1)

     ! shift each segment by dl in left-hand normal direction
     x_new(:,1,k) = Ctmp%x(k  ,:) + en*dl
     x_new(:,2,k) = Ctmp%x(k+1,:) + en*dl
  enddo
  ! set 1st and last node
  Ctmp%x(0,:) = x_new(:,1,0)
  Ctmp%x(n,:) = x_new(:,2,n-1)


! 3. working array ts_new: connect shifted line segments
  allocate (ts_new(2,0:n-1))
  ts_new(1,0)   = 0.d0
  ts_new(2,n-1) = 1.d0
  kmax          = n-1
  if (this%closed) kmax = n
  ! calculate new nodes (at the intersection of the 2 lines defined by the shiftet adjacent segments)
  do k=1,kmax
     k2  = mod(k,n)
     x11 = x_new(:,1,k-1)
     x12 = x_new(:,2,k-1)
     x21 = x_new(:,1,  k2)
     x22 = x_new(:,2,  k2)

     s   = (x12(2)-x11(2))*(x22(1)-x21(1)) - (x12(1)-x11(1))*(x22(2)-x21(2))
     s   = s / (abs(x12(2)-x11(2))+abs(x22(1)-x21(1))) &
             / (abs(x12(1)-x11(1))+abs(x22(2)-x21(2)))

     ! check if segments are parallel (then use x12 as new node)
     if (abs(s).lt.emin) then
        Ctmp%x(k,:)   = x12
        ts_new(2,k-1) = 1.d0
        ts_new(1,k2)   = 0.d0

     ! calculate intersection point
     elseif (intersect_lines (x11, x12, x21, x22, t, s, xh)) then
        Ctmp%x(k,:)   = xh
        ts_new(2,k-1) = t
        ts_new(1,k2)   = s

     else
        write (6, *) 'error in t_curve%left_hand_shift: no intersection!'
        write (6, *) 'L_1: ', x11, x12
        write (6, *) 'L_2: ', x21, x22
        stop
     endif
     if (k == n) Ctmp%x(0,:) = Ctmp%x(k,:)
  enddo
#ifdef DEBUG
  call Ctmp%plot(filename='tmp1_'//trim(str(IDEBUG))//'.plt')
#endif


! 4. check misaligned nodes
  k2       = 0
  i_remove = 0
  allocate (x_tmp(0:n,2))
  x_tmp(k2,:) = this%x(0,:)
  do k=0,n-1
     if (ts_new(1,k) .gt. ts_new(2,k)) then
        ! remove original segment
#ifdef DEBUG
        write (6, *) 'removing segment ', k
#endif
        i_remove = i_remove + 1

        ! first segment (-> skip first node)
        if (k.eq.0) then
           x_tmp(k2,:) = this%x(k+1,:)

        ! last segment (-> skip last node)
        elseif (k.eq.n-1) then
           x_tmp(k2,:) = this%x(k,:)

        ! internal segment (remove segment by extending the left and right neighbour segments)
        else
           x11 = this%x(k-1,:)
           x12 = this%x(k  ,:)
           x21 = this%x(k+1,:)
           x22 = this%x(k+2,:)
           if (intersect_lines (x11, x12, x21, x22, t, s, xh))then
              x_tmp(k2,:) = xh
           else
              x_tmp(k2,:) = 0.5d0 * (x12 + x21)
           endif
        endif
     ! just copy node data if this segment's alignment is ok
     else
        k2 = k2 + 1
        x_tmp(k2,:) = this%x(k+1,:)
     endif
  enddo

! 4.2 recursively call left_hand_shift without misaligned nodes
  n2 = n - i_remove
  if (i_remove .ne. 0) then
     call Ctmp%new(n2)
     Ctmp%x = x_tmp(0:n2,:)
#ifdef DEBUG
     call Ctmp%plot(filename='tmp2_'//trim(str(IDEBUG))//'.plt')
     IDEBUG = IDEBUG + 1
#endif
     call Ctmp%left_hand_shift(dl)
  endif


! 5. cleanup
  call this%copy(Ctmp)
  deallocate (x_new, ts_new, x_tmp)
  call Ctmp%destroy()

! 6. close curve (optional)
  if (this%closed) then
     n = this%n_seg
     x11 = 0.5d0 * (this%x(0,:) + this%x(n,:))
     this%x(0,:) = x11
     this%x(n,:) = x11
  endif

  end subroutine left_hand_shift
!=======================================================================



!=======================================================================
  function get_distance_to(this, p) result(d)
  class(t_curve)           :: this
  real(real64), intent(in) :: p(2)
  real(real64)             :: d

  real(real64) :: x1(2), x2(2), d1, dist, en(2), en0(2), en1(2), ex(2), dx, t
  integer      :: i, n


  ! 0. initialize
  d  = 1.d99
  d1 = 1.d99
  n  = this%n_seg


  ! 1. minimum distance to nodes of L (required for convex segments)
  do i=1,n-1
     x1    = this%x(i,:)
     dist  = dsqrt(sum((x1-p)**2))

     ex    = this%x(i+1,:) - this%x(i-1,:)
     ex    = ex / dsqrt(sum(ex**2))
     en(1) = ex(2)
     en(2) = -ex(1)
     if (sum(en * (p - x1)).lt.0.d0) dist = -dist

     if (dabs(dist).lt.dabs(d1)) then
        d1 = dist
        en1 = en
     endif
  enddo


  ! 2. minimum distance to segments
  do i=1,n
     x1    = this%x(i-1,:)
     x2    = this%x(i  ,:)
     ex    = x2 - x1
     dx    = dsqrt(sum(ex**2))
     ex    = ex / dx
     t     = sum(ex * (x2-p)) / dx

     ! no perpendicular line through p intersecting segment i
     if (t.lt.0.d0 .or. t.gt.1.d0) cycle


     en(1) = ex(2)
     en(2) = -ex(1)
     dist  = sum(en * (p - x1))

     if (dabs(dist).lt.dabs(d)) then
        d   = dist
        en0 = en
     endif
  enddo


  ! 3. result
  if (dabs(d1).lt.dabs(d)) then
     d=d1
     en0 = en1
  endif


  end function get_distance_to
!=======================================================================



!=======================================================================
! prepare sampling along L using the angle with regard to a reference
! point 'x_c' as weight factor
!=======================================================================
  subroutine setup_angular_sampling (L, x_c_)
  implicit none

  class(t_curve), intent(inout) :: L
  real(real64),   intent(in), optional    :: x_c_(2)

  real(real64) :: w_tot, x(2), phi1, phi2, dphi, x_c(2)
  integer      :: i, n


  n = L%n_seg

  ! check if reference point x_c is given
  if (present(x_c_)) then
     x_c = x_c_
  ! set default values otherwise
  else
     x_c(1) = sum(L%x(:,1)) / (n+1)
     x_c(2) = sum(L%x(:,2)) / (n+1)
  endif


  ! allocate memory for weight array
  !if (associated(L%w_seg)) deallocate(L%w_seg)
  !allocate (L%w_seg(n))
  if (associated(L%w)) deallocate(L%w)
  allocate (L%w(0:n))

  ! setup weight array
  w_tot = 0.d0
  x     = L%x(0,:)
  phi1  = datan2(x(2)-x_c(2), x(1)-x_c(1))
  L%w     = 0.d0
  do i=1,n
     x     = L%x(i,:)
     phi2  = datan2(x(2)-x_c(2), x(1)-x_c(1))
     dphi  = phi2 - phi1
     if (dabs(dphi) .gt. pi) dphi = dphi - dsign(pi2,dphi)

     !L%w_seg(i) = dphi
     L%w(i)     = L%w(i-1) + dphi
     !w_tot      = w_tot + dphi
     phi1  = phi2
  enddo
  !L%w_seg  = L%w_seg / w_tot
  L%w  = L%w / L%w(n)

  end subroutine setup_angular_sampling
!=======================================================================



!=======================================================================
! prepare sampling along L using the segment lengths as weight factor
!=======================================================================
  subroutine setup_length_sampling(L, raw_weights)
  class(t_curve) :: L

  logical, intent(in), optional :: raw_weights

  logical      :: integrated_weights
  real(real64) :: s, dx(L%n_dim)
  integer      :: i, n


  integrated_weights = .true.
  if (present(raw_weights)) then
     integrated_weights = .not.raw_weights
  endif


  ! allocate memory for weight array
  n = L%n_seg
  !if (associated(L%w_seg)) deallocate(L%w_seg)
  !allocate (L%w_seg(n))
  if (associated(L%w)) deallocate(L%w)
  allocate (L%w(0:n))

  ! setup weight array
  !w_tot = 0.d0
  L%l = 0.d0
  L%w = 0.d0
  do i=1,n
     dx    = L%x(i,:) - L%x(i-1,:)
     s     = dsqrt(sum(dx**2))

     L%l   = L%l + s
     if (integrated_weights) then
        L%w(i) = L%l
     else
        L%w(i) = s
     endif


     !L%w_seg(i) = s
     !L%w(i) = L%w(i-1) + s
     !w_tot      = w_tot + s
  enddo
  !L%w_seg  = L%w_seg / w_tot
  !L%l      = w_tot
  !L%l      = L%w(n)
  ! normalize weights
  if (integrated_weights) then
     L%w  = L%w / L%w(n)
  endif

  end subroutine setup_length_sampling
!=======================================================================



!=======================================================================
  subroutine setup_segment_sampling(L)
  class(t_curve) :: L

  integer      :: i, n


  ! allocate memory for weight array
  n = L%n_seg
  if (associated(L%w)) deallocate(L%w)
  allocate (L%w(0:n))

  do i=0,n
     L%w(i) = 1.d0 * i / n
  enddo


  end subroutine setup_segment_sampling
!=======================================================================




!=======================================================================
! prepare sampling along L using the ic-th coordinate as weight factor
!=======================================================================
  subroutine setup_coordinate_sampling(L, ic)
  class(t_curve)      :: L
  integer, intent(in) :: ic

  real(real64) :: w_tot, s, dx(L%n_dim)
  integer      :: i, n


  if (ic > L%n_dim) then
     write (6, *) 'error: cannot use ', ic, '-th coordinate when only ', &
                  L%n_dim, ' are defined!'
     stop
  endif

  ! allocate memory for weight array
  n = L%n_seg
  !if (associated(L%w_seg)) deallocate(L%w_seg)
  !allocate (L%w_seg(n))
  if (associated(L%w)) deallocate(L%w)
  allocate (L%w(0:n))

  ! setup weight array
  w_tot = 0.d0
  L%w = 0.d0
  do i=1,n
     s          = L%x(i,ic) - L%x(i-1,ic)
     !L%w_seg(i) = s
     L%w(i) = L%w(i-1) + s
     !w_tot      = w_tot + s
  enddo
  !L%w_seg  = L%w_seg / w_tot
  L%w  = L%w / L%w(n)

  end subroutine setup_coordinate_sampling
!=======================================================================



!=======================================================================
! return point 'x' on curve 'L' at position 'xi'
!
! optional output:
!    et		tangent vector
!=======================================================================
  subroutine sample_at (L, xi, x, et)
  use search
  implicit none

  class(t_curve), intent(in) :: L
  real(real64),  intent(in)  :: xi
  real(real64),  intent(out) :: x(L%n_dim)
  real(real64),  intent(out), optional :: et(L%n_dim)

  real(real64), parameter :: delta = 0.01d0
  real(real64) :: w
  integer      :: i, n, ierr
  integer      :: i2


  ! check if curve is set up for sampling
  if (.not.associated(L%w)) then
     write (6, *) 'error: curve is not set up for sampling!'
     stop
  endif


  ! find index i of line segment associated with t
  n = L%n_seg
  i = binary_interval_search(0, n, L%w, xi, ierr)
  if (ierr .ne. 0) then
     x = 0.d0
     return
  endif
!  n    = L%n_seg
!  wt   = t
!  wint = 0.d0
!  do i=1,n
!     wint = wint + L%w_seg(i)
!     if (wint .ge. wt) exit
!  enddo
!  if (i.eq.n+1) i=n


  ! calculate relative position on line segment
  !wt   = (wint - wt) / L%w_seg(i)
  w   = (L%w(i+1) - xi) / (L%w(i+1) - L%w(i))

  ! return point x
  !x    = L%x(i-1,:) * wt + L%x(i,:) * (1.d0 - wt)
  x    = L%x(i,:) * w + L%x(i+1,:) * (1.d0 - w)


  ! optional output: tangent vector
  if (present(et)) then
     !x1 = (L%x(i,:)-L%x(i-1,:))/L%w_seg(i)

     ! tangent vector parallel to line segment
     et = L%x(i+1,:)-L%x(i,:)


     ! transition in the vicinity of the nodes
     ! ... near lower node
     if (xi < delta) then
        i2 = i-1
        if (i2 < 0) then
           i2 = 0
           if (L%closed) i2 = n - 1
        endif

        et = L%x(i+1,:)-L%x(i2,:)
     endif
     ! ... near upper node
     if (xi > 1.d0-delta) then
        i2 = i+2
        if (i2 > n) then
           i2 = n
           if (L%closed) i2 = 1
        endif

        et = L%x(i2,:)-L%x(i,:)
     endif


     ! normalize tangent vector
     et = et / sqrt(sum(et**2))
  endif

  end subroutine sample_at
!=======================================================================



!=======================================================================
! split curve into 3 elements              |---C1---|---C2---|---C3---|
! at intrinsic coordinates xiA,xiB    xi = 0        xiA      xiB      1
!=======================================================================
  subroutine split3(this, xiA, xiB, C1, C2, C3)
  use search
  class(t_curve)             :: this
  real(real64),  intent(in)  :: xiA, xiB
  type(t_curve), intent(out) :: C1, C2, C3

  real(real64) :: tA, tB
  integer :: iA, iB, ierr, n


  ! 1. check input
  if (xiA < 0.d0  .or.  xiA > 1.d0) then
     write (6, *) 'error in t_curve%split3: xiA must be in [0,1]!'
     stop
  endif
  if (xiB < 0.d0  .or.  xiB > 1.d0) then
     write (6, *) 'error in t_curve%split3: xiB must be in [0,1]!'
     stop
  endif
  if (xiA > xiB) then
     write (6, *) 'error in t_curve%split3: xiA > xiB not allowed!'
     stop
  endif


  ! 2. find segment to split
  n  = this%n_seg
  iA = binary_interval_search(0, n, this%w, xiA, ierr)
  iB = binary_interval_search(0, n, this%w, xiB, ierr)
  tA   = (this%w(iA+1) - xiA) / (this%w(iA+1) - this%w(iA))
  tB   = (this%w(iB+1) - xiB) / (this%w(iB+1) - this%w(iB))


  ! 3. generate new curves
  call this%split3seg(iA, iB, tA, tB, C1, C2, C3)

  end subroutine split3
!=======================================================================


!=======================================================================
! split at position t within segment i
  !  |     |     |  x   |     |     |  x   |     |     |
  !  0           iA tA  iA+1        iB tB  iB+1        n
!=======================================================================
  subroutine split3seg(this, iA, iB, tA, tB, C1, C2, C3)
  class(t_curve)             :: this
  integer,       intent(in)  :: iA, iB
  real(real64),  intent(in)  :: tA, tB
  type(t_curve), intent(out) :: C1, C2, C3

  real(real64) :: x(this%n_dim)


  call C1%new(iA+1)
  C1%x            = this%x(0:iA+1,:)
  x               = this%x(iA,:) * tA + this%x(iA+1,:) * (1.d0 - tA)
  C1%x(iA+1,:)    = x

  call C2%new(iB-iA+1)
  C2%x            = this%x(iA:iB+1,:)
  C2%x(0,:)       = x
  x               = this%x(iB,:) * tB + this%x(iB+1,:) * (1.d0 - tB)
  C2%x(iB-iA+1,:) = x

  call C3%new(this%n_seg-iB)
  C3%x            = this%x(iB:this%n_seg,:)
  C3%x(0,:)       = x

  end subroutine split3seg
!=======================================================================



!=======================================================================
! split curve into n elements              |---C1---|-......-|---Cn---|
! at intrinsic coordinates xiA,xiB    xi = 0        xi1      xin-1    1
!=======================================================================
  subroutine splitn(this, n, xi, C)
  use search
  class(t_curve)             :: this
  integer,       intent(in)  :: n
  real(real64),  intent(in)  :: xi(n-1)
  type(t_curve), intent(out) :: C(n)

  real(real64) :: xxi(0:n), tsplit(n-1)
  integer :: i, ierr, n_seg, isplit(n-1)


  ! 0. initialize
  xxi(0)     = 0.d0
  xxi(1:n-1) = xi
  xxi(n)     = 1.d0


  ! 1. check input
  do i=1,n
     if (xxi(i) < xxi(i-1)) then
        write (6, *) 'error in t_curve%splitn: elements of xi must be in [0,1] and monotonically increasing!'
        write (6, *) 'xi = ', xi
        stop
     endif
  enddo


  ! 2. find segment to split
  n_seg  = this%n_seg
  do i=1,n-1
     isplit(i) = binary_interval_search(0, n_seg, this%w, xi(i), ierr)
     tsplit(i) = (xi(i) - this%w(isplit(i))) / (this%w(isplit(i)+1) - this%w(isplit(i)))
  enddo


  ! 3. generate new curves
  call this%splitnseg(n, isplit, tsplit, C)

  end subroutine splitn
!=======================================================================


!=======================================================================
! split at position tsplit within segment isplit
  !  |     |     |  x   |     |     |  x   |     |     |
  !  0           iA tA  iA+1        iB tB  iB+1        n
!=======================================================================
  subroutine splitnseg(this, n, isplit, tsplit, C)
  class(t_curve)             :: this
  integer,       intent(in)  :: n, isplit(n-1)
  real(real64),  intent(in)  :: tsplit(n-1)
  type(t_curve), intent(out) :: C(n)

  real(real64) :: x(this%n_dim), tsplit1(0:n), tA, tB
  integer      :: i, isplit1(0:n), iA, iB, m


  isplit1(0)     = 0
  isplit1(1:n-1) = isplit
  isplit1(n)     = this%n_seg-1
  tsplit1(0)     = 0.d0
  tsplit1(1:n-1) = tsplit
  tsplit1(n)     = 1.d0

  do i=1,n
     iA     = isplit1(i-1)
     iB     = isplit1(i)

     ! check input
     if (iA < 0  .or.  iA >= this%n_seg) then
        write (6, *) 'error in t_curve2D%splitnseg: invalid segment number!'
        write (6, *) 'isplit(', i-1, ') = ', iA
     endif
     if (iB < 0  .or.  iB >= this%n_seg) then
        write (6, *) 'error in t_curve2D%splitnseg: invalid segment number!'
        write (6, *) 'isplit(', i, ') = ', iB
     endif

     m      = iB - iA + 1
     call C(i)%new(m)
     C(i)%x = this%x(iA:iB+1,:)

     ! adapt first and last node
     tA     = tsplit1(i-1)
     tB     = tsplit1(i)
     C(i)%x(0,:) = this%x(iA,:) * (1.d0-tA) + this%x(iA+1,:) * tA
     C(i)%x(m,:) = this%x(iB,:) * (1.d0-tB) + this%x(iB+1,:) * tB
  enddo

  end subroutine splitnseg
!=======================================================================




!=======================================================================
  function length(this) result(L)
  class(t_curve) :: this
  real(real64)   :: L

  integer :: i


  L = 0.d0
  do i=1,this%n_seg
     L = L + sqrt(sum((this%x(i,:)-this%x(i-1,:))**2))
  enddo

  end function length
!=======================================================================



!=======================================================================
  function area(this) result(A)
  class(t_curve) :: this
  real(real64)   :: A

  real(real64)   :: dA, dl(2), x(2)
  integer :: i


  A = 0.d0
  if (.not. this%closed) return
  do i=1,this%n_seg
     dl = this%x(i,:) - this%x(i-1,:)
     x  = 0.5d0*(this%x(i,:) + this%x(i-1,:))
     dA = 0.5d0 * (dl(2)*x(1) - dl(1)*x(2))

     A = A + dA
  enddo

  end function area
!=======================================================================



!=======================================================================
  function outside(this, x)
  use math
  class(t_curve)                         :: this
  real(real64), intent(in), dimension(2) :: x
  logical                                :: outside

  real(real64), parameter :: eps = 1.d-8

  real(real64) :: Iphi, phi, phi0, dphi, dx(2)
  integer      :: i, n


  ! set default
  outside = .false.
  if (.not.this%closed) then
     write (6, *) 'warning: "outside" of open curve not defined!'
     return
  endif


  ! integrate
  n    = this%n_seg
  dx   = x - this%x(0,1:2)
  Iphi = 0.d0
  phi0 = atan2(dx(2), dx(1))
  do i=1,n
     dx = x - this%x(i,1:2)
     phi  = atan2(dx(2), dx(1))
     dphi = phi - phi0
     if (dphi.gt.pi)  dphi = dphi - 2*pi
     if (dphi.lt.-pi) dphi = dphi + 2*pi

     Iphi = Iphi + dphi
     phi0 = phi
  enddo

  if (abs(Iphi) < eps) then
     outside = .true.
  endif

  end function outside
!=======================================================================



!=======================================================================
  function connect(C1, C2) result(C)
  type(t_curve), intent(in) :: C1, C2
  type(t_curve)             :: C

  real(real64) :: x1(C1%n_dim), x2(C2%n_dim), dl
  integer :: i2, n, n1, n2


  if (C1%n_dim .ne. C2%n_dim) then
     write (6, *) 'error in module curve2D, connect:'
     write (6, *) 'cannot connect curves with unequal dimension!'
     stop
  endif


  ! check if end node of C1 matches first node of C2
  n1 = C1%n_seg
  n2 = C2%n_seg
  x1 = C1%x(n1,:)
  x2 = C1%x(0,:)
  dl = sqrt(sum((x1-x2)**2))
  i2 = 0
  n  = n1 + n2 + 1
  if (dl < epsilon(real(1.0,real64))) then
     i2 = 1
     n  = n - 1
  endif


  call C%new(n)
  C%x(0:n1,:)   = C1%x(0:n1,:)
  C%x(n1+1:n,:) = C2%x(i2:n2,:)

  end function connect
!=======================================================================



!=======================================================================
  subroutine setup_length_sampling_curvature_weighted(this)
  class(t_curve) :: this

  integer, parameter :: iu = 31

  type(t_dataset) :: kappa
  real(real64) :: dx, dy, ddx, ddy, x(2), x1(2), x2(2), L, t
  integer :: i, i1, i2, n


  n = this%n_seg
  if (associated(this%w)) deallocate(this%w)
  allocate (this%w(0:n))
  call kappa%new(n+1, 2, -1)


!  open  (iu, file=filename)
  L = 0.d0
  do i=0,n
     x = this%x(i,:)

     i1 = i+1; if (i == n) i1 = 1
     i2 = i-1; if (i == 0) i2 = n-1
     x1 = this%x(i1,:)
     x2 = this%x(i2,:)
     t  = sqrt(sum((x2-x1)**2))
     dx = x1(1) - x2(1)
     dy = x1(2) - x2(2)
     ddx = x1(1) - 2.d0*x(1) + x2(1)
     ddy = x1(2) - 2.d0*x(2) + x2(2)

     kappa%x(i,1) = (dx*ddy - dy*ddx) / (dx**2 + dy**2)**1.5d0
     kappa%x(i,2) = t
!     write (iu, *) L, kappa
!!     L = L + t

!!     this%w(i) = t * (1.d-2 + kappa)
  enddo
!  close (iu)
  call kappa%plot(filename='kappa0.plt')
  call kappa%smooth(2,1)
  call kappa%plot(filename='kappa2.plt')


  ! integrate and normalize weights
  this%w(0) = 0.d0
  do i=1,n
     t = kappa%x(i,2)
     this%w(i) = this%w(i-1) + t * (1.d-2 + kappa%x(i,1))
  enddo
  this%w = this%w / this%w(n)

  end subroutine setup_length_sampling_curvature_weighted
!=======================================================================

end module curve2D
