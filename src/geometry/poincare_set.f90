!===============================================================================
! Generate and evaluate Poincare sets
!===============================================================================
module poincare_set
  use iso_fortran_env
  use dataset
  implicit none

  private

  type, public :: t_poincare_set
     type(t_dataset), dimension(:), allocatable :: slice

     integer :: nslice, npoints, nsym
     contains
     procedure generate, plot
  end type t_poincare_set


  contains
!=======================================================================



!=======================================================================
! y0          initial/reference point
! npoints     number of points per slice
! nsym        toroidal symmetry number
! nslice      number of slices within 0..360/nsym deg
! nsteps      number of step between slices
! solver      id of ODE solver
!=======================================================================
  subroutine generate(this, y0, npoints, nsym, nslice, nsteps, solver)
  use fieldline
  class(t_poincare_set)    :: this
  real(real64), intent(in) :: y0(3)
  integer,      intent(in) :: npoints, nsym, nslice, nsteps, solver

  type (t_fieldline) :: F
  real(real64)       :: ds, L, theta
  integer            :: islice, ipoint, istep


  ! prepare memory
  if (allocated(this%slice)) deallocate(this%slice)
  allocate (this%slice(0:nslice-1))
  do islice=0,nslice-1
     call this%slice(islice)%new(npoints, 4)
  enddo
  this%nslice  = nslice
  this%npoints = npoints
  this%nsym    = nsym


  ds = pi2 / nsym / nslice / nsteps
  L  = 0.d0
  call F%init(y0, ds, solver, FL_ANGLE)

  main_loop: do ipoint=1,npoints
     slice_loop: do islice=0,nslice-1
        steps: do istep=1,nsteps
           call F%trace_1step()
           if (F%intersect_boundary()) exit main_loop
        enddo steps

        this%slice(islice)%x(ipoint,1) = F%rc(1)
        this%slice(islice)%x(ipoint,2) = F%rc(2)
        theta = F%thetac / pi * 180.d0
        if (theta < 0.d0) theta = theta + 360.d0
        this%slice(islice)%x(ipoint,3) = theta
        this%slice(islice)%x(ipoint,4) = F%PsiNc
     enddo slice_loop
  enddo main_loop

  end subroutine generate
!=======================================================================


!=======================================================================
  subroutine plot(this, iu)
  class(t_poincare_set)         :: this
  integer, intent(in), optional :: iu

  integer :: i


  do i=0,this%nslice-1
     call this%slice(i)%plot(iu=iu)
     write (iu, *)
  enddo

  end subroutine plot
!=======================================================================

end module poincare_set
