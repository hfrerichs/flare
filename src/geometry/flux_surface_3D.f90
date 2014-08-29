!===============================================================================
! Generate flux surfaces
!===============================================================================
module flux_surface_3D
  use iso_fortran_env
  use dataset
  implicit none

  private
  type, extends(t_dataset), public :: t_poincare_set
!  type, public :: t_poincare_set
!     type(t_dataset), dimension(:), allocatable :: slice
     contains
     procedure :: generate
  end type t_poincare_set

  contains
!=======================================================================


!=======================================================================
  subroutine generate(this, x0, Trace_Method, Trace_Step, N_sym, N_mult, N_points)
  class(t_poincare_set) :: this
  real(real64), intent(in) :: x0(3), Trace_Step
  integer,      intent(in) :: Trace_Method, N_sym, N_mult, N_points

  integer :: i


  if (allocated(this%x)) deallocate(this%x)
!  if (allocated(this%slice)) deallocate(this%slice)
!  allocate (this%slice(N_mult))
!  do i=1,N_mult
!     this%slice%n_row = 4
!  enddo


  end subroutine generate
!=======================================================================

end module flux_surface_3D
