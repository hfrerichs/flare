!===============================================================================
! This "mesh spacing" module provides the necessary tools to control the
! distance between grid nodes.
!===============================================================================
module mesh_spacing
  use iso_fortran_env
  use curve2D
  implicit none
  private

  type, public :: t_spacing
     integer :: mode = 0
     type(t_curve) :: D

     contains
     procedure init, node, plot
  end type t_spacing
  
  type(t_spacing), public, parameter :: Equidistant = t_spacing(0,Empty_curve)

  contains
!=======================================================================
  


!=======================================================================
  subroutine init(this, mode)
  class(t_spacing) :: this
  character(len=*) :: mode

  integer :: iB


  ! equidistant spacing
  if (mode == '') then
     this%mode = 0

  ! user defined (external) spacing function
  elseif (mode(1:5) == 'file:') then
     this%mode = -1
     iB = len_trim(mode)
     call this%D%load(mode(6:iB))
     call check_manual_stretching_function(this%D)
     call this%D%setup_coordinate_sampling(1)

  ! undefined mode string
  else
     write (6, *) 'error: undefined mode string "', trim(mode), '"!'
     stop
  endif

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


  ! equidistant sampling
  t  = 1.d0 * i / n
  xi = t


  ! user defined stretching function
  if (this%mode == -1) then
     call this%D%sample_at(t, x)
     xi = x(2)
  endif

  end function node
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
