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
     EXPONENTIAL = 1, &
     SPLINE_X1   = 2, &
     DELTA_R_SYM = 3, &
     X1          = 4, &
     USER_DEF    = -1


  type, public :: t_spacing
     integer :: mode = 0

     integer :: nc ! number of coefficients
     real(real64), dimension(:), allocatable :: c ! internal coefficients
     type(t_curve) :: D

     contains
     procedure init, node, plot
     procedure init_spline_x1
     procedure init_X1
  end type t_spacing
  
  type(t_spacing), public, parameter :: Equidistant = t_spacing(0,0,null(),Empty_curve)

  contains
!=======================================================================
  


!=======================================================================
  subroutine init(this, mode)
  use string
  class(t_spacing) :: this
  character(len=*) :: mode

  character(len=256) :: s1, s2
  integer :: iB


  if (allocated(this%c)) deallocate(this%c)
  iB = len_trim(mode)
  if (iB > 256) iB = 256


  ! equidistant spacing
  if (mode == '') then
     this%mode = 0

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
  elseif (mode(1:8) == 'Delta-R:') then
     this%mode = DELTA_R_SYM
     this%nc   = 2
     allocate (this%c(this%nc))
     s1        = parse_string(mode(9:iB),1)
     s2        = parse_string(mode(9:iB),2)
     read (s1, *, err=5000) this%c(1)
     read (s2, *, err=5000) this%c(2)
     write (6, *) 'Delta-R: ', this%c


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
  subroutine init_spline_x1(this, eta1, xi1)
  class(t_spacing)         :: this
  real(real64), intent(in) :: eta1, xi1


  this%mode = SPLINE_X1
  this%nc   = 2
  if (allocated(this%c)) deallocate(this%c)
  allocate(this%c(this%nc))
  this%c(1) = eta1
  this%c(2) = xi1

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


  ! equidistant spacings
  t  = 1.d0 * i / n
  xi = t


  select case(this%mode)
  ! exponential spacings
  case(EXPONENTIAL)
     xi = sample_exp(i, n, this%c(1))

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

  end select

  end function node
!=======================================================================



!=======================================================================
  function sample_exp(i, n, lambda) result(xi)
  integer,      intent(in) :: i, n
  real(real64), intent(in) :: lambda

  real(real64) :: xi, t, S


  S  = exp(1.d0/lambda) - 1.d0
  t  = 1.d0 * i / n
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
  subroutine plot(this, iu, filename, nsample)
  class(t_spacing) :: this
  integer, intent(in), optional          :: iu, nsample
  character(len=*), intent(in), optional :: filename

  real(real64) :: t, xi
  integer :: i, iu0 = 99


  if (present(iu)) iu0 = iu
  if (present(filename)) then
     open  (iu0, file=filename)
  endif


  ! write nodes
  do i=0,nsample
     t  = 1.d0 * i / nsample
     xi = this%node(i, nsample)
     write (iu0, *) t, xi
  enddo
  write (iu0, *)


  ! write horizontal and vertical bars
  do i=1,nsample-1
     t  = 1.d0 * i / nsample
     xi = this%node(i, nsample)

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
