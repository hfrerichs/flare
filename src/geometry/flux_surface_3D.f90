!===============================================================================
! Generate flux surfaces
!===============================================================================
module flux_surface_3D
  use iso_fortran_env
  use flux_surface_2D
  use dataset
  use curve2D
  use poincare_set
  implicit none

  private


  type, extends(t_curve) :: t_surface_cut
     real(real64) :: phi
  end type t_surface_cut


  type, public :: t_flux_surface_3D
     type(t_surface_cut), dimension(:), allocatable :: slice

     real(real64) :: PsiN
     integer      :: n_sym, n_phi

     contains
     procedure :: new, generate, plot
     procedure :: generate_from_axisymmetric_surface
     !procedure :: expand
  end type t_flux_surface_3D
  !end type t_poincare_set

  contains
!=======================================================================



!=======================================================================
  subroutine new(this, n_sym, n_phi)
  use math
  class(t_flux_surface_3D) :: this
  integer, intent(in)      :: n_sym, n_phi

  integer :: i


  this%n_sym = n_sym
  this%n_phi = n_phi
  allocate (this%slice(0:n_phi-1))
  do i=0,n_phi-1
     this%slice(i)%phi = pi2 / n_sym / n_phi * i
  enddo

  end subroutine new
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
  use equilibrium
  class(t_flux_surface_3D) :: this
  real(real64), intent(in) :: y0(3)
  integer,      intent(in) :: npoints, nsym, nslice, nsteps, solver

  type(t_poincare_set)     :: P
  real(real64) :: Maxis(3)
  integer      :: i, n


  call this%new(nsym, nslice)
  n     = nsteps
  if (n == 0) n = 16
  call P%generate(y0, npoints, nsym, nslice, n, solver, .true.)


  ! set flux surface label (given by normalized poloidal flux)
  this%PsiN = get_PsiN(y0)


  ! check if each slice is complete
  do i=0,nslice-1
     if (P%slice(i)%npoints .ne. npoints) then
        write (6, *) 'error: incomplete flux surface!'
        write (6, *) 'number of elements: ', P%slice(i)%npoints, '/', npoints
        stop
     endif
     call this%slice(i)%new (npoints-1)
     this%slice(i)%phi = P%slice(i)%phi
     this%slice(i)%x   = P%slice(i)%x(:,1:2)
     Maxis             = get_magnetic_axis(this%slice(i)%phi)
     call this%slice(i)%sort_loop(Maxis(1:2))
     call this%slice(i)%expand(Maxis(1:2), 2.d0)
  enddo

  end subroutine generate
!=======================================================================



!=======================================================================
  subroutine plot(this, iu, filename)
  class(t_flux_surface_3D)               :: this
  integer,          intent(in), optional :: iu
  character(len=*), intent(in), optional :: filename

  integer :: i, iu0


  ! set default unit number for output
  iu0 = 90

  ! Unit number given for output?
  if (present(iu)) iu0 = iu

  ! Output_File given?
  if (present(filename)) then
     open  (iu0, file=filename)
  endif

  ! write data
  do i=0,this%n_phi-1
     call this%slice(i)%plot(iu=iu0)
     write (iu0, *)
  enddo

  ! Output_File given?
  if (present(filename)) close (iu0)

  end subroutine plot
!=======================================================================



!=======================================================================
  subroutine generate_from_axisymmetric_surface(this, fs2D, n_sym, n_phi, n_theta)
  class(t_flux_surface_3D) :: this
  type(t_flux_surface_2D)  :: fs2D
  integer, intent(in)      :: n_sym, n_phi, n_theta

  real(real64) :: t, x(2)
  integer :: i, j


  call this%new(n_sym, n_phi)
  this%PsiN = fs2D%PsiN
  call fs2D%sort_loop()
  call fs2D%setup_angular_sampling()

  do i=0,n_phi-1
     call this%slice(i)%new(n_theta)
     do j=0,n_theta
        t = 1.d0 * j / n_theta
        call fs2D%sample_at(t, x)
        this%slice(i)%x(j,:) = x
     enddo
  enddo

  end subroutine generate_from_axisymmetric_surface
!=======================================================================

end module flux_surface_3D
