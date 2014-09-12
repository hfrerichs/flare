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

     ! nodes along the curve, dimension(0:n_seg, 1:n_dim)
     type (t_dataset) :: nodes
     real(real64), dimension(:,:), pointer :: x => null()

     ! relative weight for each line segment, dimension(1:n_seg)
     real(real64), dimension(:),   pointer :: w_seg => null()

     contains

     procedure :: load
     procedure :: new
     procedure :: broadcast
     procedure :: plot => curve2D_plot
     procedure :: sort_loop
     procedure :: setup_angular_sampling
     procedure :: setup_length_sampling
     procedure :: sample_at
     procedure :: length
  end type t_curve


  public :: intersect_curve, make_2D_curve

  contains
!=======================================================================



!=======================================================================
  subroutine load(this, data_file, output, header)
  class (t_curve),  intent(inout)         :: this
  character(len=*), intent(in)            :: data_file
  integer,          intent(in),  optional :: output
  character(len=*), intent(out), optional :: header


  this%n_seg =  this%nodes%nrow-1
  this%n_dim =  2
  call this%nodes%load(data_file, 2, output, header, -1)
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
  subroutine curve2D_plot(this, iu, filename)
  class(t_curve)                      :: this
  integer, intent(in), optional       :: iu
  character(len=*), intent(in), optional :: filename

  integer :: i, iu0


  ! set default unit number for output
  iu0 = 99

  ! Unit number given for output?
  if (present(iu)) iu0 = iu

  ! Output_File given?
  if (present(filename)) then
     open  (iu0, file=filename)
  endif


  ! write data
  do i=0,this%n_seg
     write (iu0, *) this%x(i,:)
  enddo


  ! Output_File given?
  if (present(filename)) then
     close (iu0)
  endif

  end subroutine curve2D_plot
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
  if (associated(L%w_seg)) deallocate(L%w_seg)
  allocate (L%w_seg(n))

  ! setup weight array
  w_tot = 0.d0
  x     = L%x(0,:)
  phi1  = datan2(x(2)-x_c(2), x(1)-x_c(1))
  do i=1,n
     x     = L%x(i,:)
     phi2  = datan2(x(2)-x_c(2), x(1)-x_c(1))
     dphi  = phi2 - phi1
     if (dabs(dphi) .gt. pi) dphi = dphi - dsign(pi2,dphi)

     L%w_seg(i) = dphi
     w_tot      = w_tot + dphi
     phi1  = phi2
  enddo
  L%w_seg  = L%w_seg / w_tot

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
  if (associated(L%w_seg)) deallocate(L%w_seg)
  allocate (L%w_seg(n))

  ! setup weight array
  w_tot = 0.d0
  do i=1,n
     dx    = L%x(i,:) - L%x(i-1,:)
     s     = dsqrt(sum(dx**2))

     L%w_seg(i) = s
     w_tot      = w_tot + s
  enddo
  L%w_seg  = L%w_seg / w_tot

  end subroutine setup_length_sampling
!=======================================================================



!=======================================================================
! sample data point 'x' from curve 'L' at position 't'
!=======================================================================
  subroutine sample_at (L, t, x, x1)
  implicit none

  class(t_curve), intent(in) :: L
  real(real64),  intent(in)  :: t
  real(real64),  intent(out) :: x(L%n_dim)
  real(real64),  intent(out), optional :: x1(L%n_dim)

  real(real64) :: wt, wint
  integer      :: i, n, ierr


  ! find index i of line segment associated with t
!  n = L%n_seg
!  i = binary_interval_search(1, n, L%w_seg, t, ierr)
!  if (ierr .ne. 0) then
!     x = 0.d0
!     return
!  endif
  n    = L%n_seg
  wt   = t
  wint = 0.d0
  do i=1,n
     wint = wint + L%w_seg(i)
     if (wint .ge. wt) exit
  enddo
  if (i.eq.n+1) i=n


  ! calculate relative position on line segment
  wt   = (wint - wt) / L%w_seg(i)

  ! return data point
  x    = L%x(i-1,:) * wt + L%x(i,:) * (1.d0 - wt)

  ! optional output: tangent vector
  if (present(x1)) then
     x1 = (L%x(i,:)-L%x(i-1,:))/L%w_seg(i)
  endif

  end subroutine sample_at
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
