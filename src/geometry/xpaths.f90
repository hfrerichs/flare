!===============================================================================
! Generate paths from the X-point in gradPsi direction
!===============================================================================
module xpaths
  use iso_fortran_env
  use curve2D
  use separatrix
  implicit none
  private


  integer, parameter, public :: &
     ASCENT_LEFT  = 1, &
     ASCENT_RIGHT = 2, &
     DESCENT_CORE = 3, &
     DESCENT_PFR  = 4

  integer, parameter, public :: &
     LIMIT_PSIN   = 1, &
     LIMIT_LENGTH = 2


  type, public, extends(t_curve) :: t_xpath
     contains
     procedure :: generate
     procedure :: setup_linear
     procedure :: setup_PsiN_sampling
  end type t_xpath

  !type(t_radial_path), dimension(:), allocatable :: radial_path


!  public :: &
!     setup_radial_paths

  contains
!=======================================================================



!=======================================================================
! Generate path along grad-Psi from Px (X-point)
! orientation = 1: ascent PsiN in left SOL direction
!             = 2: ascent PsiN in right SOL direction
!             = 3: descent PsiN to core
!             = 4: descent PsiN to PFR
!=======================================================================
  subroutine generate(this, iPx, orientation, limit_type, limit_val)
  use ode_solver
  use equilibrium
  use run_control, only: Trace_Method, N_steps
  class(t_xpath)           :: this
  integer,      intent(in) :: iPx, orientation, limit_type
  real(real64), intent(in) :: limit_val

  type(t_ODE)  :: Path
  real(real64) :: Px(2), H(2,2), v1(2), v2(2), x0(2), y(3), ds, dl, t, Psi0, PsiN, L
  integer      :: n_seg, is


  ! 0. initialize
  Px = Xp(iPx)%X
  H  = Xp(iPx)%H
  call H_eigenvectors(H, v1, v2)
  ! offset from X-point for tracing
  dl    = Px(1) / 1.d2
  ! step size
  ds    = 0.1d0 * dl
  ! position of X-point with respect to midplane
  if (Px(2) > 0.d0) v2 = - v2


  ! 1. select orientation from saddle point (X-point)
  select case(orientation)
  case(ASCENT_LEFT)
     x0    = Px - dl*v1
  case(ASCENT_RIGHT)
     x0    = Px + dl*v1
  case(DESCENT_CORE)
     x0    = Px + dl*v2
     ds    = -ds
  case(DESCENT_PFR)
     x0    = Px - dl*v2
     ds    = -ds
  case default
     write (6, *) 'error in subroutine t_gradPsiN_path%generate:'
     write (6, *) 'invalid parameter orientation = ', orientation
     stop
  end select


  ! 2.1 determine length/number of segments by final PsiN value
  select case(limit_type)
  case(LIMIT_PSIN)
     PsiN = limit_val

     n_seg  = 1
     y(1:2) = x0
     y(3)   = 0.d0
     Psi0   = get_PsiN(y)
     call Path%init_ODE(2, x0, ds, ePsi_sub, Trace_Method)
     do
        if (Psi0 < PsiN  .and.  get_PsiN(y) > PsiN) exit
        if (Psi0 > PsiN  .and.  get_PsiN(y) < PsiN) exit
        y(1:2) = Path%next_step()
        n_seg  = n_seg + 1
     enddo

  ! 2.2 expected number of segments from curve length
  case(LIMIT_LENGTH)
     L     = limit_val
     n_seg = nint((L-dl) / abs(ds)) + 1

  ! either PsiN or L must be given!
  case default
     write (6, *) 'error in subroutine t_gradPsiN_Path%generate:'
     write (6, *) 'limit must be either LIMIT_PSIN or LIMIT_LENGTH!'
     stop
  end select
  call this%new(n_seg)
  this%x(0,:) = Px
  this%x(1,:) = x0


  ! 3. generate grad PsiN path
  call Path%init_ODE(2, x0, ds, ePsi_sub, Trace_Method)
  do is=2,n_seg
     this%x(is,:) = Path%next_step()
     dl           = dl + abs(ds)
  enddo


  ! 4. adjust last node to match L
  t = 1.d0
  if (limit_type == LIMIT_LENGTH) t = (L-dl+abs(ds))/abs(ds)
  this%x(n_seg,:) = (1.d0-t)*this%x(n_seg-1,:) + t*this%x(n_seg,:)


  ! 5. setup sampling
  call this%setup_length_sampling()

  end subroutine generate
!=======================================================================



!=======================================================================
  subroutine setup_linear(this, x0, d)
  class(t_xpath)           :: this
  real(real64), intent(in) :: x0(2), d(2)


  call this%new(1)
  this%x(0,:) = x0
  this%x(1,:) = x0 + d
  call this%setup_length_sampling()

  end subroutine setup_linear
!=======================================================================



!=======================================================================
!  subroutine setup(this, mode, x0, d)
!  class(t_radial_path)         :: this
!  character(len=*), intent(in) :: mode
!  real(real64),     intent(in), optional :: x0(2), d
!
!
!  select case(mode)
!!  case('default')
!!     if (.not.present(x0)) then
!!        write (6, 9000)
!!        write (6, *) 'parameter x0 required for default operation mode!'
!!        stop
!!     endif
!  case('GradPsi')
!  case default
!     write (6, 9000)
!     write (6, *) ' operation mode "', trim(mode), '" not defined!'
!     stop
!  end select
!
! 9000 format('error in t_radial_path%setup:')
!  end subroutine setup
!=======================================================================


!=======================================================================
!  subroutine setup_radial_paths (n, mode, default_mode)
!  end subroutine setup_radial_paths
!=======================================================================



!=======================================================================
  subroutine setup_PsiN_sampling(this)
  use equilibrium
  class(t_xpath) :: this

  real(real64) :: x(3), PsiN, PsiN1, w
  integer :: i, n


  n = this%n_seg
  if (associated(this%w)) deallocate(this%w)
  allocate (this%w(0:n))
  this%w = 0

  x(3) = 0.d0
  call this%setup_length_sampling()
  open  (99, file='test_PsiN.plt')
  do i=0,this%n_seg
     x(1:2) = this%x(i,1:2)
     PsiN   = get_PsiN(x)

     write (99, *) this%w(i), PsiN, x(1:2)
     if (i>0) this%w(i) = this%w(i-1) + PsiN - PsiN1

     PsiN1 = PsiN
  enddo
  this%w = this%w / this%w(n)
  close (99)


  open  (99, file='test_PsiN_sample.plt')
  n = 10
  do i=0,n
     w = 1.d0 * i / n
     call this%sample_at(w, x(1:2))
     PsiN = get_PsiN(x)
     write (99, *) x(1:2), PsiN
  enddo
  close (99)

  end subroutine setup_PsiN_sampling
!=======================================================================

end module xpaths
