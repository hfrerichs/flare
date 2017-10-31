module grid_interface
  implicit none

  integer :: layout
  integer :: n, n1, n2, n3
  real*8, dimension(:,:), allocatable :: x
  real*8, dimension(:),   allocatable :: x1, x2, x3

  contains


  subroutine load(filename)
  use grid
  character(len=*), intent(in)  :: filename

  type(t_grid) :: G


  ! cleanup
  if (allocated(x))  deallocate(x)
  if (allocated(x1)) deallocate(x1, x2, x3)


  ! load grid layout and nodes
  call G%load(filename)
  layout = G%layout
  n = G%n
  allocate (x(n,3))
  x = G%x


  ! set up structured grid
  if (layout == STRUCTURED) then
     n1 = G%n1
     n2 = G%n2
     n3 = G%n3
     allocate (x1(n1), x2(n2), x3(n3))
     x1 = G%x1
     x2 = G%x2
     if (G%fixed_coord == 0) then
        x3 = G%x3
     else
        x3 = G%fixed_coord_value
     endif
  endif

  end subroutine load

end module grid_interface
