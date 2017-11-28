module boundary_interface
  implicit none

  integer, parameter :: dp = selected_real_kind(15,307) ! double precision


  ! boundary slice data
  real(dp), dimension(:,:), allocatable :: x
  integer :: n


  contains
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function num_boundaries()
  use boundary
  integer :: num_boundaries


  num_boundaries = n_boundary1

  end function num_boundaries
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine get_slice(iboundary, phi, ierr)
  use curve2D
  use boundary
  integer,  intent(in)  :: iboundary
  real(dp), intent(in)  :: phi
  integer,  intent(out) :: ierr

  type(t_curve) :: S


  if (iboundary < 1  .or.  iboundary > n_boundary1) then
     n    = 0
     ierr = 1
     return
  endif

  ierr = 0
  if (allocated(x)) deallocate(x)
  S = boundary_slice(iboundary, phi)
  n = S%n_seg+1
  allocate (x(n,2))
  x = S%x

  end subroutine get_slice
  !---------------------------------------------------------------------
  
end module boundary_interface
