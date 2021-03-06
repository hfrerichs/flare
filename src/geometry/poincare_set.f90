!===============================================================================
! Generate and evaluate Poincare sets
!===============================================================================
module poincare_set
  use iso_fortran_env
  use dataset
  implicit none

  private

  type, extends(t_dataset) :: t_slice
     ! toroidal position of slice [rad]
     real(real64) :: phi

     ! actual number of points on slice
     integer      :: npoints
  end type t_slice


  type, public :: t_poincare_set
     type(t_slice), dimension(:), allocatable :: slice

     integer :: nslice, npoints_max, nsym

     contains
     procedure :: generate
     procedure :: plot
  end type t_poincare_set


  contains
!=======================================================================



!=======================================================================
! y0          initial/reference point
! npoints     number of points per slice
! nsym        toroidal symmetry number
! nslice      number of slices within 0..360/nsym deg
! nsteps      number of trace steps between slices
! solver      id of ODE solver
! offset      toroidal offset between reference point and first location for Poincare set
!=======================================================================
  subroutine generate(this, y0, npoints, nsym, nslice, nsteps, solver, stop_at_boundary, offset, bfield)
  use fieldline
  use equilibrium, only: get_poloidal_angle
  class(t_poincare_set)    :: this
  real(real64), intent(in) :: y0(3)
  integer,      intent(in) :: npoints, nsym, nslice, nsteps, solver
  logical,      intent(in) :: stop_at_boundary
  real(real64), intent(in), optional :: offset
  integer,      intent(in), optional :: bfield

  type (t_fieldline) :: F
  real(real64)       :: y1(3), ds, ds1, dl, theta, theta0
  integer            :: islice, islice0, ipoint, istep, n


  ! prepare memory
  if (allocated(this%slice)) deallocate(this%slice)
  allocate (this%slice(0:nslice-1))
  do islice=0,nslice-1
     call this%slice(islice)%new(npoints, 4)
     this%slice(islice)%npoints = 0
     this%slice(islice)%phi     = y0(3) + islice * pi2 / nsym / nslice
  enddo
  this%nslice      = nslice
  this%npoints_max = npoints
  this%nsym        = nsym


  ! optional offset between reference point and location of first Poincare set
  if (present(offset)) then
     dl = offset
     do
        if (dl >= 0.d0) exit
        dl = dl + pi2 / nsym
     enddo
  else
     dl = 0.d0
  endif


  ds = pi2 / nsym / nslice / nsteps
  ! trace initial/reference point to location of first Poincare set
  if (dl > 0.d0) then
     n   = ceiling(dl / ds)
     ds1 = dl / n
     call F%init(y0, ds, solver, FL_ANGLE, bfield=bfield)
     do istep=1,n
        dl = F%trace_1step()
        if (stop_at_boundary  .and.  F%intersect_boundary()) return
     enddo
     y1  = F%rc
  else
     y1  = y0
  endif


  call F%init(y1, ds, solver, FL_ANGLE, bfield=bfield)

  theta0 = get_poloidal_angle(y1) / pi * 180.d0
  main_loop: do ipoint=1,npoints
     slice_loop: do islice=1,nslice
        steps: do istep=1,nsteps
           dl = F%trace_1step()
           if (stop_at_boundary  .and.  F%intersect_boundary()) exit main_loop
        enddo steps

        islice0 = mod(islice,nslice)
        this%slice(islice0)%npoints     = ipoint
        this%slice(islice0)%x(ipoint,1) = F%rc(1)
        this%slice(islice0)%x(ipoint,2) = F%rc(2)
        theta = F%thetac / pi * 180.d0
        if (theta < theta0) theta = theta + 360.d0
        this%slice(islice0)%x(ipoint,3) = theta
        this%slice(islice0)%x(ipoint,4) = F%PsiNc
     enddo slice_loop
  enddo main_loop


  ! re-arrange points (sort points using poloidal angle)
  do islice=0,nslice-1
     call this%slice(islice)%sort_rows(3)
  enddo

  end subroutine generate
!=======================================================================


!=======================================================================
  subroutine plot(this, iu)
  class(t_poincare_set)         :: this
  integer, intent(in), optional :: iu

  integer :: i, iu0

  iu0 = 97
  do i=0,this%nslice-1
     call this%slice(i)%plot(iu=iu0, nelem=this%slice(i)%npoints)
     write (iu0, *)
  enddo

  end subroutine plot
!=======================================================================

end module poincare_set
