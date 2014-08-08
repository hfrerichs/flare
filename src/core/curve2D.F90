!===============================================================================
! module:	curve2D
!
! description:	Operations on 2D curves (line segments)
!===============================================================================
module curve2D

#if defined(FLARE)
  use math
  implicit none
#else
  implicit none

  real*8, parameter :: &
     pi   = 3.14159265358979323846264338328d0, &
     pi2  = 2.d0 * pi
#endif


  !type, extends(t_nodes/t_curve) :: t_spline
  type t_curve
     ! number of line segments (n_seg), number of coordinates (n_dim)
     integer  :: n_seg, n_dim

     ! nodes along the curve, dimension(0:n_seg, 1:n_dim)
     real*8, dimension(:,:), pointer :: x_data => null()

     ! relative weight for each line segment, dimension(1:n_seg)
     real*8, dimension(:),   pointer :: w_seg => null()

     contains

     procedure :: read
  end type t_curve

  contains
!=======================================================================



!=======================================================================
  subroutine read(this, data_file, columns, report, header)
  class (t_curve), intent(inout) :: this
  character*120,   intent(in)    :: data_file
  integer,         intent(in), optional  :: columns
  logical,         intent(in), optional  :: report
  character*80,    intent(out), optional :: header

  integer, parameter :: iu = 42

  real*8, dimension(:), allocatable :: tmp
  character*120 :: str
  integer       :: i, j, ncount, ncol, icom
  logical       :: lreport


  ! display messages
  lreport = .true.
  if (present(report)) then
     if (report .eqv. .false.) lreport = .false.
  endif


  ! number of data columns
  ncol = 2
  if (present(columns)) then
     ncol = columns
  endif


  ! header
  if (present(header)) then
     header = ''
  endif


  ! parse date file to get number of data lines
  open  (iu, file=data_file)
  ncount = 0
  icom   = 0
  parse_loop: do
     read (iu, 4000, end=2000) str
     do while (str(1:1) .eq. '#')
        if (icom == 0  .and. present(header)) header = str(3:82)
        read (iu, *, end=2000) str
        icom = icom + 1
     enddo
     ncount = ncount + 1
  enddo parse_loop
 2000 rewind(iu)
  if (lreport) write (6,1000) ncount, data_file(1:len_trim(data_file))


  ! allocate memory
  this%n_seg = ncount-1
  if (associated(this%x_data)) deallocate (this%x_data)
  allocate (this%x_data(0:ncount-1,ncol))
  allocate (tmp(ncol))


  ! read actual data
  j = 0
  read_loop: do i=1,ncount+icom
     read (iu, 4000) str
     if (str(1:1) .ne. '#') then
        read (str, *) tmp
        this%x_data(j,:) = tmp
        j = j + 1
     endif
  enddo read_loop
  close (iu)
  this%n_dim = ncol
  deallocate (tmp)


  return
 1000 format ('found ',i6,' data lines in file: ',a)
 4000 format (a120)
  end subroutine read
!=======================================================================



!=======================================================================
! Generate 2D curve from two 1D arrays (xarr, yarr) of size n
!=======================================================================
  subroutine make_2D_curve (n, xarr, yarr, C)
  integer, intent(in)        :: n
  real*8,  intent(in)        :: xarr(n), yarr(n)
  type(t_curve), intent(out) :: C

  C%n_seg = n-1
  C%n_dim = 2
  allocate (C%x_data(0:n-1,2))
  C%x_data(:,1) = xarr
  C%x_data(:,2) = yarr

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
  real*8,        intent(in)            :: x1(2), x2(2)
  type(t_curve), intent(in)            :: C
  real*8,        intent(out), optional :: xh(2), th, sh
  integer,       intent(out), optional :: ish
  logical                              :: intersect_curve

  real*8  :: t, s, xl1(2), xl2(2), xh0(2)
  integer :: is


  intersect_curve = .false.
  t = 0.d0
  s = 0.d0

  do is=1,C%n_seg
     xl1 = C%x_data(is-1,:)
     xl2 = C%x_data(is  ,:)
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
  real*8, intent(in)  :: L1(2), L2(2), M1(2), M2(2)
  real*8, intent(out) :: l, m
  real*8, optional    :: X(2)
  logical             :: intersect_lines

  real*8 :: e1(2), e2(2), d1, d2, edet, dx(2), n0(2), b


  intersect_lines = .false.
  l   = 0.d0
  m   = 0.d0

  e1   = L2 - L1
  d1   = dsqrt(sum(e1**2))
  e1   = e1 / d1
  e2   = M2 - M1
  d2   = dsqrt(sum(e2**2))
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


end module curve2D
