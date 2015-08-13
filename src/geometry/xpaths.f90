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
     ASCENT        = 1, &
     ASCENT_LEFT   = 1, &
     ASCENT_RIGHT  = 2, &
     DESCENT       = 3, &
     DESCENT_CORE  = 3, &
     DESCENT_PFR   = 4

  integer, parameter, public :: &
     LIMIT_PSIN    = 1, &
     LIMIT_LENGTH  = 2, &
     LIMIT_CURVE   = 3, &
     SAMPLE_PSIN   = 1, &
     SAMPLE_LENGTH = 2


  type, public, extends(t_curve) :: t_xpath
     real(real64) :: PsiN(2)
     contains
     procedure :: generateX
     procedure :: generate
     procedure :: setup_linear
     procedure :: setup_PsiN_sampling
     procedure :: sample_at_PsiN
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
  subroutine generateX(this, iPx, orientation, limit_type, limit_val, sampling)
  use ode_solver
  use equilibrium
  use run_control, only: Trace_Method, N_steps
  class(t_xpath)           :: this
  integer,      intent(in) :: iPx, orientation, limit_type
  real(real64), intent(in) :: limit_val
  integer,      intent(in), optional :: sampling

  type(t_ODE)  :: Path
  real(real64) :: Px(2), H(2,2), v1(2), v2(2), x0(2), dl


  ! 0. initialize
  Px = Xp(iPx)%X
  H  = Xp(iPx)%H
  call H_eigenvectors(H, v1, v2)
  ! offset from X-point for tracing
  dl    = Px(1) / 1.d3
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
  case(DESCENT_PFR)
     x0    = Px - dl*v2
  case default
     write (6, *) 'error in subroutine t_gradPsiN_path%generate:'
     write (6, *) 'invalid parameter orientation = ', orientation
     stop
  end select


  ! 2. generate path from x0
  call this%generate(x0, orientation, limit_type, limit_val, Px, sampling)

  end subroutine generateX
!=======================================================================



!=======================================================================
! Generate path along grad-Psi from x0 (non-hyperbolic point)
! direction = 1,2: ascent PsiN
!           = 3,4: descent PsiN
!=======================================================================
  subroutine generate(this, xinit, direction, limit_type, limit_val, x0, sampling, C_limit)
  use ode_solver
  use equilibrium
  use run_control, only: Trace_Method, N_steps
  class(t_xpath)           :: this
  real(real64),  intent(in) :: xinit(2), limit_val
  integer,       intent(in) :: direction, limit_type
  real(real64),  intent(in), optional :: x0(2)
  integer,       intent(in), optional :: sampling
  type(t_curve), intent(in), optional :: C_limit

  integer, parameter :: n_tmp0 = 1000

  real(real64), dimension(:,:), allocatable :: tmp, tmp_tmp
  type(t_ODE)  :: Path
  real(real64) :: y(3), ds, t, Psi0, PsiN, L
  integer      :: i, n_seg, is, n_tmp, sampling_method


  ! 0. check input
  select case(limit_type)
  case(LIMIT_PSIN,LIMIT_LENGTH)

  case(LIMIT_CURVE)
     if (.not.present(C_limit)) then
        write (6, *) 'error in subroutine t_xpath%generate:'
        write (6, *) 'C_limit must be set for limit_type = LIMIT_CURVE!'
        stop
     endif

  case default
     write (6, *) 'error in subroutine t_xpath%generate:'
     write (6, *) 'limit must be either LIMIT_PSIN, LIMIT_LENGTH or LIMIT_CURVE!'
     write (6, *) 'limit_type = ', limit_type
     stop
  end select


  ! 1. initialize
  ! step size
  ds    = sqrt(sum(xinit**2)) * 0.5d-4
  if (direction == DESCENT_CORE  .or.  direction == DESCENT_PFR) then
     ds = - ds
  endif

  ! temporary data array
  n_tmp = n_tmp0
  allocate (tmp(n_tmp, 3))

  i     = 0
  L     = 0.d0
  ! set 0th data point (if present)
  if (present(x0)) then
     i     = 1
     tmp(1,1:2) = x0
     tmp(1,  3) = get_PsiN(x0)
     L          = sqrt(sum((xinit-x0)**2))
  endif

  ! set 1st data point
  i = i + 1
  tmp(i,1:2) = xinit
  tmp(i,  3) = get_PsiN(xinit)


  ! 2. generate grad PsiN path
  call Path%init_ODE(2, xinit, ds, ePsi_sub, Trace_Method)
  do
     i = i + 1

     ! increase temporary array, if necessary
     if (i > n_tmp) then
        allocate (tmp_tmp(n_tmp, 3))
        tmp_tmp = tmp
        deallocate (tmp)

        allocate (tmp(n_tmp+n_tmp0, 3))
        tmp(1:n_tmp,:) = tmp_tmp
        n_tmp          = n_tmp + n_tmp0
        deallocate (tmp_tmp)
     endif

     ! calculate next step
     tmp(i,1:2)  = Path%next_step()
     tmp(i,  3)  = get_PsiN(tmp(i,1:2))
     L           = L  + abs(ds)

     ! last step?
     select case(limit_type)
     ! 1. final PsiN value given
     case(LIMIT_PSIN)
        t = (limit_val - tmp(i-1,3)) / (tmp(i,3) - tmp(i-1,3))
        ! PsiN rises above limit
        if (tmp(1,3) < limit_val  .and.  tmp(i,3) > limit_val) exit

        ! PsiN falls below limit
        if (tmp(1,3) > limit_val  .and.  tmp(i,3) < limit_val) exit

     ! 2. curve length given
     case(LIMIT_LENGTH)
        t = (limit_val - L + abs(ds)) / abs(ds)
        if (L >= limit_val) exit

     ! 3. limiting contour given, but at least to PsiN=limit_val
     case(LIMIT_CURVE)
        if ((limit_val - tmp(i,3))*(limit_val - tmp(1,3)) < 0.d0) then
           if (intersect_curve(tmp(i-1,1:2), tmp(i,1:2), C_limit, th=t)) exit
        endif

     end select
  enddo
  n_seg = i-1
  call this%new(n_seg)
  this%x(:,:) = tmp(1:i,1:2)


  ! 3. adjust last node to match boundary condition
  this%x(n_seg,:) = (1.d0-t)*this%x(n_seg-1,:) + t*this%x(n_seg,:)
  this%PsiN(1) = tmp(1,3)
  this%PsiN(2) = tmp(i,3)


  ! 4. setup sampling
  sampling_method = SAMPLE_LENGTH
  if (present(sampling)) sampling_method = sampling
  select case(sampling_method)
  case(SAMPLE_PSIN)
     allocate (this%w(0:n_seg))
     this%w = (tmp(1:i,3) - this%PsiN(1)) / (this%PsiN(2) - this%PsiN(1))

  case(SAMPLE_LENGTH)
     call this%setup_length_sampling()

  case default
     write (6, *) 'error in subroutine t_xpath%generate:'
     write (6, *) 'sampling must be either SAMPLE_PSIN or SAMPLE_LENGTH!'
     stop
  end select


  ! 5. cleanup
  deallocate (tmp)

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



!=======================================================================
  subroutine sample_at_PsiN(this, PsiN, x)
  class(t_xpath)            :: this
  real(real64), intent(in)  :: PsiN
  real(real64), intent(out) :: x(2)

  real(real64) :: eta


  eta = (PsiN - this%PsiN(1)) / (this%PsiN(2) - this%PsiN(1))
  call this%sample_at(eta, x)

  end subroutine sample_at_PsiN
!=======================================================================

end module xpaths
