!===============================================================================
! Generate flux surfaces
!===============================================================================
module flux_surface_3D
  use iso_fortran_env
  use flux_surface_2D
  use dataset
  implicit none

  private
  type, extends(t_dataset), public :: t_poincare_set
     real(real64) :: phi
  end type t_poincare_set



  type, public :: t_flux_surface_3D
!  type, extends(t_dataset), public :: t_flux_surface_3D
!  type, extends(t_dataset), public :: t_poincare_set
!  type, public :: t_poincare_set
!     type(t_dataset), dimension(:), allocatable :: slice
     integer :: n_sym, n_phi
     real(real64) :: PsiN
     type(t_poincare_set), dimension(:), allocatable :: slice
     contains
     procedure :: new
     procedure :: generate, write
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


  ! cleanup old data
  if (allocated(this%slice)) then
     do i=0,n_phi-1
        call this%slice(i)%destroy()
     enddo
  endif


  this%n_sym = n_sym
  this%n_phi = n_phi
  allocate (this%slice(0:n_phi-1))
  do i=0,n_phi-1
     this%slice(i)%phi = pi2 / n_sym / n_phi * i
  enddo

  end subroutine new
!=======================================================================



!=======================================================================
  subroutine generate(this, x0, Trace_Method, Trace_Step, N_sym, N_mult, N_points)
  !class(t_poincare_set) :: this
  class(t_flux_surface_3D) :: this
  real(real64), intent(in) :: x0(3), Trace_Step
  integer,      intent(in) :: Trace_Method, N_sym, N_mult, N_points

  integer :: i


  !if (allocated(this%x)) deallocate(this%x)
!  if (allocated(this%slice)) deallocate(this%slice)
!  allocate (this%slice(N_mult))
!  do i=1,N_mult
!     this%slice%n_row = 4
!  enddo


  end subroutine generate
!=======================================================================



!=======================================================================
  subroutine write(this, iu, Output_File)
  class(t_flux_surface_3D) :: this
  integer, intent(in), optional        :: iu
  character*120, intent(in), optional  :: Output_File

  integer :: i, j, iu0


  ! set default unit number for output
  iu0 = 90

  ! Unit number given for output?
  if (present(iu)) iu0 = iu

  ! Output_File given?
  if (present(Output_File)) then
     open  (iu0, file=Output_File)
  endif

  ! write data
  do i=0,this%n_phi-1
     do j=1,this%slice(i)%nrow
        write (iu0, *) this%slice(i)%x(j,:), this%slice(i)%phi
     enddo
     write (iu0, *)
  enddo

  ! Output_File given?
  if (present(Output_File)) close (iu0)

  end subroutine write
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
     call this%slice(i)%new(n_theta,2)
     do j=1,n_theta
        t = 1.d0 * j / n_theta
        call fs2D%sample_at(t, x)
        this%slice(i)%x(j,:) = x
     enddo
  enddo

  end subroutine generate_from_axisymmetric_surface
!=======================================================================

end module flux_surface_3D
