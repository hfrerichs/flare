module cspline
  use iso_fortran_env
  use dataset
  use fgsl
  implicit none
  private


  integer, parameter, public :: &
     ABSOLUTE      = 0, &
     NORMALIZED    = 1


  type, public :: t_cspline
     real(real64) :: L1, L2

     integer(fgsl_size_t)    :: n, ndim
     type(fgsl_spline), dimension(:), allocatable :: S
     type(fgsl_interp_accel), dimension(:), allocatable :: acc

     contains
     procedure :: setup
     procedure :: destroy
     procedure :: eval
     procedure :: plot
  end type t_cspline

  contains
!=======================================================================



!=======================================================================
  subroutine setup(this, D, base_variable, periodic)
  class(t_cspline)    :: this
  type(t_dataset)     :: D
  integer, intent(in), optional :: base_variable
  logical, intent(in), optional :: periodic

  real(real64), dimension(:), allocatable :: t
  integer(fgsl_int)    :: i, i0, m0, stat


  this%n    = D%nrow
  this%ndim = D%ncol
  if (allocated(this%S))   deallocate(this%S)
  if (allocated(this%acc)) deallocate(this%acc)
  allocate (this%S(D%ncol), this%acc(D%ncol))
  if (present(periodic) .and. periodic) then
     do i=1,this%ndim
        this%S(i) = fgsl_spline_alloc(fgsl_interp_cspline_periodic, this%n)
     enddo
  else
     do i=1,this%ndim
        this%S(i) = fgsl_spline_alloc(fgsl_interp_cspline, this%n)
     enddo
  endif


  ! check input for base variable
  m0 = 0
  if (present(base_variable)) m0 = base_variable
  if (m0 < 0  .or.  m0 > D%ncol) then
     write (6, *) 'error in t_cspline%setup: 0 <= base_variable <= D%ncol required'
     write (6, *) 'm0 = ', m0
     stop
  endif


  ! setup reference coordinate
  allocate (t(this%n))
  if (m0 == 0) then
     t  = 0.d0
     i0 = D%nrow_offset
     do i=2,this%n
        t(i) = t(i-1) + sqrt(sum((D%x(i+i0,:)-D%x(i+i0-1,:))**2))
     enddo
  else
     t = D%x(:,m0)
  endif
  this%L1 = t(1)
  this%L2 = t(this%n)


  ! initialize spline
  do i=1,this%ndim
     stat   = fgsl_spline_init(this%S(i), t, D%x(:,i), this%n)
  enddo


  ! cleanup
  deallocate (t)

  end subroutine setup
!=======================================================================



!=======================================================================
  subroutine destroy(this)
  class(t_cspline) :: this

  integer :: i


  do i=1,this%ndim
     call fgsl_spline_free(this%S(i))
     call fgsl_interp_accel_free(this%acc(i))
  enddo

  end subroutine destroy
!=======================================================================



!=======================================================================
  function eval(this, t, base, derivative) result(y)
  class(t_cspline)         :: this
  real(real64), intent(in) :: t
  real(real64)             :: y(this%ndim)
  integer,      intent(in), optional :: base, derivative

  real(real64) :: l
  integer      :: i, m


  ! evaluate location is absolute or normalized
  l = this%L1 + t *  (this%L2 - this%L1)
  if (present(base)) then
     if (base == ABSOLUTE) l = t
  endif


  ! evaluate derivate?
  m = 0
  if (present(derivative)) m = derivative


  do i=1,this%ndim
  select case(m)
  case(0)
     y(i) = fgsl_spline_eval(this%S(i), l, this%acc(i))
  case(1)
     y(i) = fgsl_spline_eval_deriv(this%S(i), l, this%acc(i))
  case(2)
     y(i) = fgsl_spline_eval_deriv2(this%S(i), l, this%acc(i))

  case default
     write (6, *) 'error: derivative = ', m, ' not defined!'
     stop
  end select
  enddo

  end function eval
!=======================================================================



!=======================================================================
  subroutine plot(this, filename, nsample)
  class(t_cspline)             :: this
  character(len=*), intent(in) :: filename
  integer,          intent(in) :: nsample

  integer, parameter :: iu = 50

  real(real64) :: t
  integer :: i

  open  (iu, file=filename)
  do i=0,nsample
     t = 1.d0 * i / nsample
     write (iu, *) this%eval(t)
  enddo
  close (iu)

  end subroutine plot
!=======================================================================

end module cspline
