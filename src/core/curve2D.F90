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

     contains

     procedure :: load
     procedure :: new
     procedure :: destroy
     procedure :: copy
     procedure :: broadcast
     procedure :: plot
     procedure :: sort_loop
     procedure :: expand
     procedure :: left_hand_shift
     procedure :: get_distance_to
     procedure :: setup_angular_sampling
     procedure :: setup_length_sampling
     procedure :: setup_coordinate_sampling
     procedure :: sample_at
     procedure :: split3
     procedure :: length
  end type t_curve

  type(t_curve), public, parameter :: Empty_curve = t_curve(0,0,0.d0,Empty_dataset,null(),null())


  public :: intersect_curve, make_2D_curve

  contains
!=======================================================================



!=======================================================================
  subroutine load(this, data_file, output, header)
  class (t_curve),  intent(inout)         :: this
  character(len=*), intent(in)            :: data_file
  integer,          intent(in),  optional :: output
  character(len=*), intent(out), optional :: header


  call this%nodes%load(data_file, 2, output, header, -1)
  this%n_seg =  this%nodes%nrow-1
  this%n_dim =  2
  this%x     => this%nodes%x

  end subroutine load
!=======================================================================



!=======================================================================
  subroutine new(this, n_seg)
  class(t_curve)      ::  this
  integer, intent(in) :: n_seg


  this%n_seg = n_seg
  this%n_dim = 2
  call this%nodes%new(n_seg+1, 2, -1)
  this%x     => this%nodes%x

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

  end subroutine copy
!=======================================================================



!=======================================================================
  subroutine broadcast(this)
  class(t_curve) :: this

  integer :: n, m


#if defined(FLARE)
  n = this%n_seg
  m = this%n_dim

  call broadcast_inte_s (this%n_seg)
  call broadcast_inte_s (this%n_dim)

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
!===============================================================================
  function intersect_curve (x1, x2, C, xh, th, sh, ish)
  real(real64),  intent(in)            :: x1(2), x2(2)
  type(t_curve), intent(in)            :: C
  real(real64),  intent(out), optional :: xh(2), th, sh
  integer,       intent(out), optional :: ish
  logical                              :: intersect_curve

  real(real64) :: t, s, xl1(2), xl2(2), xh0(2)
  integer      :: is


  intersect_curve = .false.
  t = 0.d0
  s = 0.d0

  do is=1,C%n_seg
     xl1 = C%x(is-1,:)
     xl2 = C%x(is  ,:)
     if (intersect_lines(x1,x2,xl1,xl2,t,s,xh0)) then
        if (t.ge.0.d0 .and. t.le.1.d0 .and. s.ge.0.d0 .and. s.le.1.d0) then
           intersect_curve = .true.
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
              ish = is
           endif
           return
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
  subroutine plot(this, iu, filename, append)
  class(t_curve) :: this
  integer,          intent(in), optional :: iu
  character(len=*), intent(in), optional :: filename
  logical,          intent(in), optional :: append

  integer :: i, iu0


  ! set default unit number for output
  iu0 = 99

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
  do i=0,this%n_seg
     write (iu0, *) this%x(i,:)
  enddo


  ! Output_File given?
  if (present(filename)) then
     close (iu0)
  endif

  end subroutine plot
!=======================================================================



!=======================================================================
! sort points on a closed contour line 'L' with regard to center 'x_c'
! and direction 'dx'
!=======================================================================
  subroutine sort_loop (L, x_c_, dx_)
  implicit none

  class(t_curve), intent(inout)          :: L
  real(real64),     intent(in), optional :: x_c_(2), dx_(2)

  real(real64), dimension(:,:), allocatable :: tmp_arr
  real(real64) :: x_c(2), dx(2)
  real(real64) :: theta_0, theta, dxl(2), x_0(2), tmp(3), t
  integer      :: i, j, n


  n = L%n_seg


  ! check if reference point x_c is given
  if (present(x_c_)) then
     x_c = x_c_
  ! set default values otherwise
  else
     x_c(1) = sum(L%x(:,1)) / (n+1)
     x_c(2) = sum(L%x(:,2)) / (n+1)
  endif


  ! check if reference direction is given
  if (present(dx_)) then
     dx    = dx_
  ! set default values otherwise
  else
     dx(1) = 1.d0
     dx(2) = 0.d0
  endif


  ! allocate memory for temporary array
  allocate (tmp_arr(0:n, 3))
  tmp_arr(:,1:2) = L%x

  ! calculate angles w.r.t. x_c
  theta_0  = datan2(dx(2), dx(1))
  do i=0,n
     dxl   = L%x(i,:) - x_c
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
     L%x(0:n,:) = tmp_arr(0:n,1:2)
  else
     t    = (theta_0      - tmp_arr(n,3) + pi2) / t
     x_0  = tmp_arr(n,1:2)
     dxl  = tmp_arr(0,1:2) - x_0
     x_0  = x_0 + t * dxl

     ! copy data to output variable L
     if (x_0(1).eq.tmp_arr(0,1) .and. x_0(2).eq.tmp_arr(0,2)) then
        call L%new(n+1)
        L%x(0:  n,:) = tmp_arr(0:n,1:2)
        L%x(  n+1,:) = x_0
     else
        call L%new(n+2)
        L%x(0    ,:) = x_0
        L%x(1:n+1,:) = tmp_arr(0:n,1:2)
        L%x(  n+2,:) = x_0
     endif
  endif

  end subroutine sort_loop
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
! Shift each line segment in direction of its left-hand side normal vector and
! re-connect nodes.
!
!                     ^ normal vector
!                     |
! segment: ====================>
!
!=======================================================================
  recursive subroutine left_hand_shift(this, dl)
  class(t_curve)           :: this
  real(real64), intent(in) :: dl

  real(real64), parameter  :: emin = 1.d-8

  real(real64), dimension(:,:,:), allocatable :: x_new
  real(real64), dimension(:,:),   allocatable :: ts_new, x_tmp
  type(t_curve) :: Ctmp
  real(real64)  :: el(2), en(2), x11(2), x12(2), x21(2), x22(2), xh(2), t, s
  integer       :: k, n, k2, n2, i_remove


! 1. initialize local variables
  n = this%n_seg
  if (this%n_dim .ne. 2) then
     write (6, *) 'error in subroutine t_curve%left_hand_shift: n_dim = 2 expected!'
     stop
  endif
  call Ctmp%copy(this)


! 2. working array x_new(i,j,k): raw nodes positions
  ! i:	coordinate
  ! j:	lower (1) and upper (2) node
  ! k:	segment number
  allocate (x_new(2,2,0:n-1))
  do k=0,n-1
     ! direction of line segment
     el = Ctmp%x(k+1,:) - Ctmp%x(k,:)
     el = el / dsqrt(sum(el**2))

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
  ! calculate new nodes (at the intersection of the 2 lines defined by the shiftet adjacent segments)
  do k=1,n-1
     x11 = x_new(:,1,k-1)
     x12 = x_new(:,2,k-1)
     x21 = x_new(:,1,  k)
     x22 = x_new(:,2,  k)

     s   = (x12(2)-x11(2))*(x22(1)-x21(1)) - (x12(1)-x11(1))*(x22(2)-x21(2))
     s   = s / (abs(x12(2)-x11(2))+abs(x22(1)-x21(1))) &
             / (abs(x12(1)-x11(1))+abs(x22(2)-x21(2)))

     ! check if segments are parallel (then use x12 as new node)
     if (abs(s).lt.emin) then
        Ctmp%x(k,:)   = x12
        ts_new(2,k-1) = 1.d0
        ts_new(1,k)   = 0.d0

     ! calculate intersection point
     elseif (intersect_lines (x11, x12, x21, x22, t, s, xh)) then
        Ctmp%x(k,:)   = xh
        ts_new(2,k-1) = t
        ts_new(1,k)   = s

     else
        write (6, *) 'error in t_curve%left_hand_shift: no intersection!'
        write (6, *) 'L_1: ', x11, x12
        write (6, *) 'L_2: ', x21, x22
        stop
     endif
  enddo


! 4. check misaligned nodes
  k2       = 0
  i_remove = 0
  allocate (x_tmp(0:n,2))
  x_tmp(k2,:) = this%x(0,:)
  do k=0,n-1
     if (ts_new(1,k) .gt. ts_new(2,k)) then
        ! remove original segment
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
     call Ctmp%left_hand_shift(dl)
  endif


! 5. cleanup
  call this%copy(Ctmp)
  deallocate (x_new, ts_new, x_tmp)
  call Ctmp%destroy()

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
  subroutine setup_length_sampling(L)
  class(t_curve) :: L

  real(real64) :: w_tot, s, dx(L%n_dim)
  integer      :: i, n


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
     dx    = L%x(i,:) - L%x(i-1,:)
     s     = dsqrt(sum(dx**2))

     !L%w_seg(i) = s
     L%w(i) = L%w(i-1) + s
     w_tot      = w_tot + s
  enddo
  !L%w_seg  = L%w_seg / w_tot
  !L%l      = w_tot
  L%l      = L%w(n)
  L%w  = L%w / L%w(n)

  end subroutine setup_length_sampling
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
! sample data point 'x' from curve 'L' at position 't'
!=======================================================================
  subroutine sample_at (L, t, x, x1)
  use search
  implicit none

  class(t_curve), intent(in) :: L
  real(real64),  intent(in)  :: t
  real(real64),  intent(out) :: x(L%n_dim)
  real(real64),  intent(out), optional :: x1(L%n_dim)

  real(real64) :: wt, wint
  integer      :: i, n, ierr


  ! find index i of line segment associated with t
  n = L%n_seg
  i = binary_interval_search(0, n, L%w, t, ierr)
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
  wt   = (L%w(i+1) - t) / (L%w(i+1) - L%w(i))

  ! return data point
  !x    = L%x(i-1,:) * wt + L%x(i,:) * (1.d0 - wt)
  x    = L%x(i,:) * wt + L%x(i+1,:) * (1.d0 - wt)

  ! optional output: tangent vector
  if (present(x1)) then
     !x1 = (L%x(i,:)-L%x(i-1,:))/L%w_seg(i)
     x1 = (L%x(i+1,:)-L%x(i,:))/L%w(i)
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

  real(real64) :: tA, tB, x(this%n_dim)
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
  !  |     |     |  x   |     |     |  x   |     |     |
  !  0           iA xiA iA+1        iB xiB iB+1        n
  call C1%new(iA+1)
  C1%x            = this%x(0:iA+1,:)
  x               = this%x(iA,:) * tA + this%x(iA+1,:) * (1.d0 - tA)
  C1%x(iA+1,:)    = x

  call C2%new(iB-iA+1)
  C2%x            = this%x(iA:iB+1,:)
  C2%x(0,:)       = x
  x               = this%x(iB,:) * tB + this%x(iB+1,:) * (1.d0 - tB)
  C2%x(iB-iA+1,:) = x

  call C3%new(n-iB)
  C3%x            = this%x(iB:n,:)
  C3%x(0,:)       = x

  end subroutine split3
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

end module curve2D
