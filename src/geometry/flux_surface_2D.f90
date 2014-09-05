!===============================================================================
! Generete unperturbed flux surfaces
!===============================================================================
module flux_surface_2D
  use iso_fortran_env
  use curve2D
  implicit none

  private
  type, extends(t_curve) :: t_flux_surface_2D
     real(real64) :: PsiN

     contains
     procedure :: generate => generate_flux_surface_2D
  end type t_flux_surface_2D

  public :: t_flux_surface_2D

  contains
!=======================================================================


!=======================================================================
! Generate axisymmetric flux surface at r = (R,Z [cm]) using a step size
! of Trace_Step and integration method Trace_Method.
! Default operation mode is tracing in both directions, but tracing in
! one direction only is set by the optional parameter direction.
! An alternate limiting surface can be given by the optional parameter AltSurf
!=======================================================================
  subroutine generate_flux_surface_2D(this, r, direction, Trace_Step, Trace_Method, AltSurf)
  use equilibrium
  use ode_solver
  use boundary
  class(t_flux_surface_2D) :: this
  real(real64), intent(in) :: r(2)
  
  integer, intent(in), optional       :: direction, Trace_Method
  real(real64), intent(in), optional  :: Trace_Step
  type(t_curve), intent(in), optional :: AltSurf

  type(t_ODE) :: F
  real*8, dimension(:,:), allocatable :: tmp
  real*8  :: yl(3), yc(3), X(3), ds, r3(3)
  integer :: idir, i, nmax, imethod, id, n(-1:1)


  ! determine trace step
  if (present(Trace_Step)) then
     ds = Trace_Step
  else
     ds = length_scale() / 200.d0
  endif


  ! set trace method
  if (present(Trace_Method)) then
     imethod = Trace_Method
  else
     imethod = NM_AdamsBashforth4
  endif


  ! allocate temporary array for output
  nmax = N_steps_guess(ds) / 2
  allocate (tmp(-nmax:nmax,2))


  ! trace direction
  n         = nmax
  if (present(direction)) then
     if (abs(direction) .ne. 1) then
        write (6, *) 'error: invalid parameter direction ', direction
     endif
     n(-direction) = 0
  endif


  ! initialize variables
  yl(3)     = 0.d0
  yc(3)     = 0.d0
  r3(1:2)   = r
  r3(3)     = 0.d0
  this%PsiN = get_PsiN(r3)


  ! trace in forward and backward direction
  do idir=-1,1,2
     ! set initial position
     call F%init_ODE(2, r, idir*ds, Bpol_sub, imethod)
     yl(1:2)  = r
     tmp(0,:) = F%yc

     do i=1,nmax
        ! progress one step
        yc(1:2)       = F%next_step()

        ! save current position
        tmp(idir*i,:) = yc(1:2)

        ! check intersection with boundary
        if (intersect_boundary(yl, yc, X, id)) then
           tmp(idir*i,:) = X(1:2)
           n(idir)       = i
           exit
        endif

        ! prepare next step
        yl = yc
     enddo
  enddo


! save data
  allocate (this%x_data(0:n(-1)+n(1),2))
  this%x_data = tmp(-n(-1):n(1),:)
  this%n_seg  = n(-1)+n(1)
  this%n_dim  = 2
  deallocate (tmp)

  end subroutine generate_flux_surface_2D
!=======================================================================



!=======================================================================
! interface subroutine for poloidal magnetic field
!=======================================================================
  subroutine Bpol_sub (n, t, y, f)
  use equilibrium

  integer, intent(in) :: n
  real*8,  intent(in) :: t, y(n)
  real*8, intent(out) :: f(n)

  real*8 :: Bf(3), Bpol


  Bf = get_Bf_eq2D(y)
  Bpol = sqrt(Bf(1)**2 + Bf(2)**2)
  f    = Bf(1:2)/Bpol

  end subroutine Bpol_sub
!=======================================================================



!===============================================================================
! N_STEPS_GUESS
!
! Guess the required number of trace steps (for flux surface generation),
! based on the step size.
!===============================================================================
  function N_steps_guess (Trace_Step) result(N)
  use equilibrium
  real*8, intent(in) :: Trace_Step
  integer :: N

  real*8  :: dl, Rbox(2), Zbox(2)


  call get_domain (Rbox, Zbox)
  dl = 2.0d0 * ((Zbox(2)-Zbox(1)) + (Rbox(2)-Rbox(1)))
  dl = dl / dabs(Trace_Step)
  N  = int(dl)

  end function N_steps_guess
!===============================================================================


end module flux_surface_2D
