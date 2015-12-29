!===============================================================================
! This "mesh spacing" module provides the necessary tools to control the
! distance between grid nodes.
!===============================================================================
module mesh_spacing
  use iso_fortran_env
  use curve2D
  implicit none
  private


  integer, parameter :: &
     LINEAR      = 0, &
     EXPONENTIAL = 1, &
     SPLINE_X1   = 2, &
     DELTA_R_SYM = 3, &
     X1          = 4, &
     D_CURVE     = 5, &
     USER_DEF    = -1, &
     S_RECURSIVE = -2


  type, public :: t_spacing
     integer :: mode = 0

     integer :: nc ! number of coefficients
     real(real64), dimension(:), allocatable :: c ! internal coefficients
     type(t_curve) :: D
     type(t_spacing), dimension(:), pointer  :: S

     contains
     procedure init
     procedure node
     procedure sample
     procedure plot
     procedure init_spline_x1
     procedure init_X1
     procedure init_recursive
     procedure init_Dcurve
  end type t_spacing
  
  type(t_spacing), public, parameter :: Equidistant = t_spacing(0,0,null(),Empty_curve,null())

  contains
!=======================================================================
  


!=======================================================================
  subroutine init(this, mode)
  use string
  class(t_spacing) :: this
  character(len=*) :: mode

  character(len=256) :: s1, s2, s3
  real(real64) :: eps, kap, del, phi0, dphi, Delta, R
  integer :: iB


  if (allocated(this%c)) deallocate(this%c)
  iB = len_trim(mode)
  if (iB > 256) iB = 256


  ! equidistant spacing
  if (mode == '') then
     this%mode = LINEAR

  ! exponential spacing function
  elseif (mode(1:4) == 'exp:') then
     this%mode = EXPONENTIAL
     this%nc   = 1
     allocate(this%c(this%nc))
     read (mode(5:iB), *, err=5000) this%c(1)

  ! spline with reference node
!  elseif (mode(1:3) == 'X1:') then
     !this%mode = SPLINE_X1
     !call this%init_spline_x1()

  ! Delta-R type spacing function (increase resolution in Delta domain by factor R)
  ! Symmetric version
  elseif (mode(1:8) == 'Delta-R:') then
     this%mode = DELTA_R_SYM
     this%nc   = 2
     allocate (this%c(this%nc))
     s1        = parse_string(mode(9:iB),1)
     s2        = parse_string(mode(9:iB),2)
     read (s1, *, err=5000) this%c(1)
     read (s2, *, err=5000) this%c(2)
     write (6, *) 'Delta-R: ', this%c


  ! Delta-R type spacing function (increase resolution in Delta domain by factor R)
  ! at lower boundary
  elseif (mode(1:9) == 'Delta-R0:') then
     this%mode = X1
     this%nc   = 2
     allocate (this%c(this%nc))
     s1        = parse_string(mode(10:iB),1)
     s2        = parse_string(mode(10:iB),2)
     read (s1, *, err=5000) Delta
     read (s2, *, err=5000) R
     this%c(1) = R * Delta / (1.d0 + Delta * (R-1.d0))
     this%c(2) = Delta
     write (6, *) 'Delta-R0: ', Delta, R


  ! Delta-R type spacing function (increase resolution in Delta domain by factor R)
  ! at upper boundary
  elseif (mode(1:9) == 'Delta-R1:') then
     this%mode = X1
     this%nc   = 2
     allocate (this%c(this%nc))
     s1        = parse_string(mode(10:iB),1)
     s2        = parse_string(mode(10:iB),2)
     read (s1, *, err=5000) Delta
     read (s2, *, err=5000) R
     this%c(1) = (1.d0 - Delta) / (1.d0 + Delta * (R-1.d0))
     this%c(2) = 1.d0 - Delta
     write (6, *) 'Delta-R1: ', Delta, R


  ! 1 additional node
  elseif (mode(1:3) == 'X1:') then
     this%mode = X1
     this%nc   = 2
     allocate (this%c(this%nc))
     s1        = parse_string(mode(4:iB),1)
     s2        = parse_string(mode(4:iB),2)
     read (s1, *, err=5000) this%c(1)
     read (s2, *, err=5000) this%c(2)
     write (6, *) 'X1: ', this%c


  ! D-curve
  elseif (mode(1:8) == 'D-curve:') then
     s1        = parse_string(mode(9:iB),1)
     read (s1, *, err=5000) eps
     s1        = parse_string(mode(9:iB),2)
     read (s1, *, err=5000) kap
     s1        = parse_string(mode(9:iB),3)
     read (s1, *, err=5000) del
     s1        = parse_string(mode(9:iB),4)
     read (s1, *, err=5000) phi0
     s1        = parse_string(mode(9:iB),5)
     read (s1, *, err=5000) dphi
     write (6, *) 'D-curve: ', eps, kap, del, phi0, dphi
     call this%init_Dcurve(eps, kap, del, phi0, dphi)


  ! user defined (external) spacing function
  elseif (mode(1:5) == 'file:') then
     this%mode = USER_DEF
     call this%D%load(mode(6:iB))
     call check_manual_stretching_function(this%D)
     call this%D%setup_coordinate_sampling(1)

  ! undefined mode string
  else
     write (6, *) 'error: undefined mode string "', trim(mode), '"!'
     stop
  endif

  return
 5000 write (6, *) 'error in t_spacing%init: could not read value from string ', this%mode
  stop
  end subroutine init
!=======================================================================



!=======================================================================
  subroutine init_spline_x1(this, eta1, xi1, ierr)
  class(t_spacing)         :: this
  real(real64), intent(in) :: eta1, xi1
  integer,      intent(out), optional :: ierr

  real(real64), parameter  :: beta = 2.d0
  real(real64) :: a1, b1, a2, b2, c2


  if (present(ierr)) ierr = 0

  this%mode = SPLINE_X1
  this%nc   = 2
  if (allocated(this%c)) deallocate(this%c)
  allocate(this%c(this%nc))
  this%c(1) = eta1
  this%c(2) = xi1

  a2 = (1.d0 - xi1) / (1.d0 - eta1)**2 - 1.d0 / beta / (1.d0 - eta1)
  b2 = 1.d0 / beta - 2.d0 * eta1 * a2
  c2 = 1.d0 - a2 - b2

  b1 = 2.d0 * xi1 / eta1 - 1.d0 / beta
  a1 = (xi1 - eta1 * b1) / eta1**2
  ! check if slope at eta=0 is positive
  if (xi1 < eta1 / 2.d0 / beta) then
     write (6, *) 'error: spacing function is ill-defined!'
     write (6, *) 'eta1, xi1, b1 = ', eta1, xi1, b1
     if (present(ierr)) then
        ierr = 1
     else
        stop
     endif
  endif
  ! check if slope at eta=1 is positive
  if (2.d0 * (1.d0-xi1) / (1.d0-eta1) - 1.d0/beta  .le.  0.d0) then
     write (6, *) 'error: spacing function is ill-defined!'
     write (6, *) 'eta1, xi1, a2, b2, c2 = ', eta1, xi1, a2, b2, c2
     if (present(ierr)) then
        ierr = 2
     else
        stop
     endif
  endif
  end subroutine init_spline_x1
!=======================================================================



!=======================================================================
  subroutine init_X1(this, eta1, xi1)
  class(t_spacing)         :: this
  real(real64), intent(in) :: eta1, xi1


  this%mode = X1
  this%nc   = 2
  if (allocated(this%c)) deallocate(this%c)
  allocate(this%c(this%nc))
  this%c(1) = eta1
  this%c(2) = xi1

  end subroutine init_X1
!=======================================================================



!=======================================================================
  subroutine init_recursive(this, S_left, S_right, eta1, rho1)
  class(t_spacing)            :: this
  type(t_spacing), intent(in) :: S_left, S_right
  real(real64),    intent(in) :: eta1, rho1


  ! 0. check input
  if (eta1 < 0.d0  .or.  eta1 > 1.d0) then
     write (6, *) 'error in t_spacing%init_recursive: invalid parameter!'
     write (6, *) 'eta1 = ', eta1
     stop
  endif
  if (rho1 < 0.d0  .or.  rho1 > 1.d0) then
     write (6, *) 'error in t_spacing%init_recursive: invalid parameter!'
     write (6, *) 'rho1 = ', rho1
     stop
  endif


  ! 1. recursive definition of spacing function
  this%mode = S_RECURSIVE
  allocate (this%S(2))
  this%S(1) = S_left
  this%S(2) = S_right

  if (allocated(this%c)) deallocate(this%c)
  allocate (this%c(2))
  this%c(1) = eta1
  this%c(2) = rho1

  end subroutine init_recursive
!=======================================================================



!=======================================================================
  subroutine init_Dcurve(this, eps, kap, del, phi0, dphi)
  use math
  class(t_spacing)            :: this
  real(real64), intent(in)    :: eps, kap, del, phi0, dphi

  integer, parameter :: n = 720
  real(real64)       :: alp, t, K, x1, x2, y1, y2, curv, tau
  integer            :: i


  this%mode = USER_DEF
  alp       = asin(del)

  call this%D%new(n)
  do i=1,n
     this%D%x(i,2) = 1.d0 * i / n

     t    = (i-0.5d0) / n
     tau  = (dphi * t  +  phi0) / 180.d0 * pi
     x1   = -eps*sin(tau + alp*sin(tau)) * (1.d0 + alp*cos(tau))
     x2   = -eps*cos(tau + alp*sin(tau)) * (1.d0 + alp*cos(tau))**2 + eps*sin(tau + alp*sin(tau)) * alp * sin(tau)
     y1   = eps*kap*cos(tau)
     y2   = -eps*kap*sin(tau)
     curv = (x1*y2 - y1*x2) / (x1**2 + y1**2)**(1.5d0)


     !K    = 1.d0/curv
     K    = curv
     this%D%x(i,1) = this%D%x(i-1,1) + K
  enddo
  this%D%x(:,1) = this%D%x(:,1) / this%D%x(n,1)

  !call this%D%plot(filename='Dcurve.tmp')
  call check_manual_stretching_function(this%D)
  call this%D%setup_coordinate_sampling(1)

  end subroutine init_Dcurve
!=======================================================================



!=======================================================================
  subroutine check_manual_stretching_function (D)
  type(t_curve), intent(in) :: D

  real(real64) :: dx(2)
  integer      :: i, n


  n = D%n_seg

  ! check x(0) = 0
  if (D%x(0,1).ne.0.d0 .or. D%x(0,2).ne.0.d0) then
     write (6, 9000)
     write (6, *) 'first data point must be (0.0, 0.0)!'
     stop
  endif

  ! check x(1) = 1
  if (D%x(n,1).ne.1.d0 .or. D%x(n,2).ne.1.d0) then
     write (6, 9000)
     write (6, *) 'last data point must be (1.0, 1.0)!'
     stop
  endif

  ! check if x(xi) is strictly increasing
  do i=0,n-1
     dx = D%x(i+1,:) - D%x(i,:)
     if (dx(1).le.0.d0) then
        write (6, 9000)
        write (6, *) 'x-coordinate not strictly increasing at node ', i+1
        stop
     endif
     if (dx(2).le.0.d0) then
        write (6, 9000)
        write (6, *) 'xi-coordinate not strictly increasing at node ', i+1
        stop
     endif
     if (D%x(i,1).ge.1.d0) then
        write (6, 9000)
        write (6, *) 'interior xi-coordinate is equal or exceeds 1 at node ', i
        stop
     endif
     if (D%x(i,2).ge.1.d0) then
        write (6, 9000)
        write (6, *) 'interior x-coordinate is equal or exceeds 1 at node ', i
        stop
     endif
  enddo


  9000 format ('error: invalid stretching function!')
  end subroutine check_manual_stretching_function
!=======================================================================



!=======================================================================
  function node(this, i, n) result(xi)
  class(t_spacing)    :: this
  integer, intent(in) :: i, n
  real(real64)        :: xi

  real(real64) :: t, x(2)


  t  = 1.d0 * i / n
  xi = this%sample(t)

  end function node
!=======================================================================



!=======================================================================
  function sample(this, t) result(xi)
  class(t_spacing)         :: this
  real(real64), intent(in) :: t
  real(real64)             :: xi

  real(real64) :: x(2)


  select case(this%mode)
  ! linear spacing
  case(LINEAR)
     xi = t

  ! exponential spacings
  case(EXPONENTIAL)
     xi = sample_exp(t, this%c(1))

  ! spline with reference node
  case(SPLINE_X1)
     xi = sample_spline_X1(t, this%c(1), this%c(2))

  ! Delta-R type spacing function (increase resolution in Delta domain by factor R)
  case(DELTA_R_SYM)
     xi = sample_Delta_R_sym(t, this%c(1), this%c(2))

  ! 1 additional node
  case(X1)
     xi = sample_X1(t, this%c(1), this%c(2))

  ! user defined spacings
  case(USER_DEF)
     call this%D%sample_at(t, x)
     xi = x(2)

  ! recursive definition
  case(S_RECURSIVE)
     xi = sample_recursive(t, this%S(1), this%S(2), this%c(1), this%c(2))
  end select

  end function sample
!=======================================================================



!=======================================================================
  function sample_exp(t, lambda) result(xi)
  real(real64), intent(in) :: t, lambda

  real(real64) :: xi, S


  S  = exp(1.d0/lambda) - 1.d0
  xi = (exp(t/lambda) - 1.d0)/S

  end function sample_exp
!=======================================================================



!=======================================================================
  function sample_spline_X1(eta, eta1, xi1) result(xi)
  real(real64), intent(in) :: eta, eta1, xi1
  real(real64)             :: xi

  real(real64), parameter  :: beta = 2.d0 ! dxi/deta (eta1) = 0.5

  real(real64) :: a1, a2, b1, b2, c2


  a2 = (1.d0 - xi1) / (1.d0 - eta1)**2 - 1.d0 / beta / (1.d0 - eta1)
  b2 = 1.d0 / beta - 2.d0 * eta1 * a2
  c2 = 1.d0 - a2 - b2

  b1 = 2.d0 * xi1 / eta1 - 1.d0 / beta
  a1 = (xi1 - eta1 * b1) / eta1**2
  ! TODO: this is already checked in init_spline_X1 and could be removed
  ! check if slope at eta=0 is positive
  if (xi1 < eta1 / 2.d0 / beta) then
     write (6, *) 'error: spacing function is ill-defined!'
     write (6, *) 'eta1, xi1, b1 = ', eta1, xi1, b1
     stop
  endif
  ! check if slope at eta=1 is positive
  if (2.d0 * (1.d0-xi1) / (1.d0-eta1) - 1.d0/beta  .le.  0.d0) then
     write (6, *) 'error: spacing function is ill-defined!'
     write (6, *) 'eta1, xi1, a2, b2, c2 = ', eta1, xi1, a2, b2, c2
     stop
  endif

  if (eta .lt. eta1) then
     xi = a1 * eta**2.d0  +  b1 * eta
  else
     xi = a2 * eta**2.d0  +  b2 * eta  +  c2
  endif

  end function sample_spline_X1
!=======================================================================



!=======================================================================
  function sample_Delta_R_sym(eta, Delta, R) result(xi)
  real(real64), intent(in) :: eta, Delta, R
  real(real64)             :: xi

  real(real64) :: eta1


  eta1 = R * Delta / (1.d0 + (R-1.d0) * 2.d0 * Delta)

  if (eta < eta1) then
     xi = eta * Delta / eta1
  elseif (eta > 1.d0-eta1) then
     xi = eta * Delta / eta1 + 1.d0 - Delta / eta1
  else
     xi = Delta + (eta-eta1) * (1.d0 - 2.d0*Delta) / (1.d0 - 2.d0*eta1)
  endif

  end function sample_Delta_R_sym
!=======================================================================



!=======================================================================
  function sample_X1(eta, eta1, r1) result(xi)
  real(real64), intent(in) :: eta, eta1, r1
  real(real64)             :: xi


  if (eta < eta1) then
     xi = eta * r1 / eta1
  else
     xi = r1 + (eta-eta1) * (1.d0-r1) / (1.d0-eta1)
  endif

  end function sample_X1
!=======================================================================



!=======================================================================
  function sample_recursive(eta, S_left, S_right, eta1, rho1) result(xi)
  real(real64),    intent(in) :: eta, eta1, rho1
  type(t_spacing), intent(in) :: S_left, S_right
  real(real64)                :: xi

  real(real64) :: t

  if (eta < eta1) then
     t  = eta / eta1
     xi = S_left%sample(t) * rho1
  else
     t  = (eta-eta1) / (1.d0-eta1)
     xi = rho1 + S_right%sample(t) * (1.d0-rho1)
  endif

  end function sample_recursive
!=======================================================================



!=======================================================================
  subroutine plot(this, iu, filename, nsample)
  class(t_spacing) :: this
  integer, intent(in), optional          :: iu, nsample
  character(len=*), intent(in), optional :: filename

  real(real64) :: t, xi
  integer :: i, n, iu0 = 99


  if (present(iu)) iu0 = iu
  if (present(filename)) then
     open  (iu0, file=filename)
  endif

  ! set number of sampling points
  n = 20
  if (present(nsample)) n = nsample


  ! write nodes
  do i=0,n
     t  = 1.d0 * i / n
     xi = this%node(i, n)
     write (iu0, *) t, xi
  enddo
  write (iu0, *)


  ! write horizontal and vertical bars
  do i=1,n-1
     t  = 1.d0 * i / n
     xi = this%node(i, n)

     write (iu0, *) t, 0.d0
     write (iu0, *) t, xi
     write (iu0, *) 0.d0, xi
     write (iu0, *)
  enddo


  if (present(filename)) then
     close (iu0)
  endif

  end subroutine plot
!=======================================================================

end module mesh_spacing
