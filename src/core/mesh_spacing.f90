!===============================================================================
! This "mesh spacing" module provides the necessary tools to control the
! distance between grid nodes.
!===============================================================================
module mesh_spacing
  use iso_fortran_env
  use curve2D
  implicit none
  private


  character(len=*), parameter :: &
     LINEAR      = 'LINEAR', &
     EXPONENTIAL = 'EXPONENTIAL', &
     SPLINE_X1   = 'SPLINE_X1', &
     DELTA_R_SYM = 'DELTA_R_SYM', &
     DELTA_R_LB  = 'DELTA_R_LB', &
     DELTA_R_UB  = 'DELTA_R_UB', &
     X1          = 'X1', &
     X1_SMOOTH   = 'X1_SMOOTH', &
     D_CURVE     = 'D_CURVE', &
     USER_DEF    = 'LOAD', &
     S_RECURSIVE = 'RECURSIVE'


  type, public :: t_spacing
     character(len=256) :: dist = LINEAR

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
     procedure init_manual
  end type t_spacing
  
  type(t_spacing), public, parameter :: Equidistant = t_spacing(LINEAR,0,null(),Empty_curve,null())

  contains
!=======================================================================
  


!=======================================================================
  subroutine init(this, distribution)
  use string
  class(t_spacing) :: this
  character(len=*) :: distribution

  character(len=max(len(distribution), len(LINEAR))) :: distribution_type
  character(len=len(distribution))                   :: arguments

  real(real64) :: eps, kap, del, phi0, dphi, Delta, R, eta1, xi1, beta


  if (allocated(this%c)) deallocate(this%c)
  ! default distribution is linear i.e. equidistant
  arguments = ''
  if (distribution == '') then
     distribution_type = LINEAR
  else
     read (distribution, *) distribution_type
     if (len_trim(distribution_type)+2  <=  len_trim(distribution)) then
        arguments         = distribution(len_trim(distribution_type)+2:len_trim(distribution))
     endif
  endif


  this%dist = distribution_type
  select case(distribution_type)
  ! equidistant spacing
  case(LINEAR)
     ! nothing to be done here


  ! exponential spacing function
  case(EXPONENTIAL)
     this%nc   = 1
     allocate(this%c(this%nc))
     read (arguments, *, err=5000, end=5000) this%c


  ! spline with reference node
  case(SPLINE_X1)
     read (arguments, *, err=5000, end=5000) eta1, xi1, beta
     call this%init_spline_X1(eta1, xi1, beta)


  ! Delta-R type spacing function (increase resolution in Delta domain by factor R)
  ! Symmetric version
  case(DELTA_R_SYM)
     this%nc   = 2
     allocate (this%c(this%nc))
     read (arguments, *, err=5000, end=5000) this%c


  ! Delta-R type spacing function (increase resolution in Delta domain by factor R)
  ! at lower boundary
  case(DELTA_R_LB)
     this%dist = X1
     this%nc   = 2
     allocate (this%c(this%nc))
     read (arguments, *, err=5000, end=5000) Delta, R
     this%c(1) = R * Delta / (1.d0 + Delta * (R-1.d0))
     this%c(2) = Delta


  ! Delta-R type spacing function (increase resolution in Delta domain by factor R)
  ! at upper boundary
  case(DELTA_R_UB)
     this%dist = X1
     this%nc   = 2
     allocate (this%c(this%nc))
     read (arguments, *, err=5000, end=5000) Delta, R
     this%c(1) = (1.d0 - Delta) / (1.d0 + Delta * (R-1.d0))
     this%c(2) = 1.d0 - Delta


  ! 1 additional node (piecewise linear)
  case(X1)
     this%nc   = 2
     allocate (this%c(this%nc))
     read (arguments, *, err=5000, end=5000) this%c


  ! 1 additional node (linear segment + smooth extension)
  case(X1_SMOOTH)
     this%nc   = 2
     allocate (this%c(this%nc))
     read (arguments, *, err=5000, end=5000) this%c


  ! D-curve
  case(D_CURVE)
     read (arguments, *, err=5000, end=5000) eps, kap, del, phi0, dphi
     call this%init_Dcurve(eps, kap, del, phi0, dphi)


  ! user defined (external) spacing function
  case(USER_DEF)
     call this%D%load(arguments)
     call check_manual_stretching_function(this%D)
     call this%D%setup_coordinate_sampling(1)


  ! undefined mode string
  case default
     write (6, *) 'error: undefined distribution "', trim(distribution), '"!'
     stop
  end select

  return
 5000 write (6, *) 'error in t_spacing%init: could not read value from string "', trim(arguments), '"'
  stop
  end subroutine init
!=======================================================================



!=======================================================================
  subroutine init_spline_x1(this, eta1, xi1, beta, ierr)
  class(t_spacing)         :: this
  real(real64), intent(in) :: eta1, xi1, beta
  integer,      intent(out), optional :: ierr

  real(real64) :: a1, b1, a2, b2, c2


  if (present(ierr)) ierr = 0

  this%dist = SPLINE_X1
  this%nc   = 3
  if (allocated(this%c)) deallocate(this%c)
  allocate(this%c(this%nc))
  this%c(1) = eta1
  this%c(2) = xi1
  this%c(3) = beta

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


  this%dist = X1
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
  this%dist = S_RECURSIVE
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


  this%dist = USER_DEF
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
  subroutine init_manual(this, C)
  class(t_spacing)            :: this
  type(t_curve), intent(in)   :: C


  this%dist = USER_DEF

  this%D = C
  call check_manual_stretching_function(this%D)
  call this%D%setup_coordinate_sampling(1)

  end subroutine init_manual
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


  if (i == 0) then
     xi = 0.d0
     return
  elseif (i == n) then
     xi = 1.d0
     return
  endif

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


  select case(this%dist)
  ! linear spacing
  case(LINEAR)
     xi = t

  ! exponential spacings
  case(EXPONENTIAL)
     xi = sample_exp(t, this%c(1))

  ! spline with reference node
  case(SPLINE_X1)
     xi = sample_spline_X1(t, this%c(1), this%c(2), this%c(3))

  ! Delta-R type spacing function (increase resolution in Delta domain by factor R)
  case(DELTA_R_SYM)
     xi = sample_Delta_R_sym(t, this%c(1), this%c(2))

  ! 1 additional node
  case(X1)
     xi = sample_X1(t, this%c(1), this%c(2))

  ! 1 additional node (linear segment + smooth extension)
  case(X1_SMOOTH)
     xi = sample_X1_SMOOTH(t, this%c(1), this%c(2))

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
  function sample_spline_X1(eta, eta1, xi1, beta) result(xi)
  real(real64), intent(in) :: eta, eta1, xi1, beta
  real(real64)             :: xi

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
  function sample_X1_smooth(eta, eta1, xi1) result(xi)
  real(real64), intent(in) :: eta, eta1, xi1
  real(real64)             :: xi

  real(real64) :: a, b, c


  if (eta < eta1) then
     xi = eta * xi1 / eta1
  else
     a  = (eta1 - xi1) / (1.d0 - eta1)**2 / eta1
     b  = xi1 / eta1 - 2.d0 * a * eta1
     c  = 1.d0 - a - b
     xi = a*eta**2 + b*eta + c
  endif

  end function sample_X1_smooth
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
  subroutine plot(this, iu, filename, nsample, style)
  class(t_spacing) :: this
  integer, intent(in), optional          :: iu, nsample
  character(len=*), intent(in), optional :: filename, style

  integer, parameter :: &
     STYLE_FUNCTION = 1, &
     STYLE_MESH     = 2

  real(real64) :: t, xi
  integer :: i, n, istyle, iu0 = 99


  if (present(iu)) iu0 = iu
  if (present(filename)) then
     open  (iu0, file=filename)
  endif
  istyle = STYLE_MESH
  if (present(style)) then
     select case(style)
     case('function')
        istyle = STYLE_FUNCTION
     case('mesh')
        istyle = STYLE_MESH
     case default
        write (6, *) 'error in t_spacing%plot: invalid style ', trim(style)
        stop
     end select
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


  if (istyle == STYLE_MESH) then
  ! write horizontal and vertical bars
  do i=1,n-1
     t  = 1.d0 * i / n
     xi = this%node(i, n)

     write (iu0, *) t, 0.d0
     write (iu0, *) t, xi
     write (iu0, *) 0.d0, xi
     write (iu0, *)
  enddo
  endif


  if (present(filename)) then
     close (iu0)
  endif

  end subroutine plot
!=======================================================================

end module mesh_spacing
