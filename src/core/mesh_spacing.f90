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
     USER_DEF    = -1


  type, public :: t_spacing
     integer :: mode = 0

     integer :: nc ! number of coefficients
     real(real64), dimension(:), allocatable :: c ! internal coefficients
     type(t_curve) :: D

     contains
     procedure init, node, plot
  end type t_spacing
  
  type(t_spacing), public, parameter :: Equidistant = t_spacing(0,0,null(),Empty_curve)

  contains
!=======================================================================
  


!=======================================================================
  subroutine init(this, mode)
  class(t_spacing) :: this
  character(len=*) :: mode

  integer :: iB


  if (allocated(this%c)) deallocate(this%c)
  iB = len_trim(mode)


  ! equidistant spacing
  if (mode == '') then
     this%mode = 0

  ! exponential spacing function
  elseif (mode(1:3) == 'exp:') then
     this%mode = EXPONENTIAL
     this%nc   = 1
     allocate(this%c(this%nc))
     read (mode(5:iB), *, err=5000) this%c(1)

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
