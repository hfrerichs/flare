module grid_interface
  implicit none

  integer, parameter :: dp = selected_real_kind(15,307)


  integer :: layout, coordinates
  integer :: n, n1, n2, n3, coord1, coord2, fixed_coord
  real(dp), dimension(:,:), allocatable :: x
  real(dp), dimension(:),   allocatable :: x1, x2, x3
  character(len=120) :: coord_label1, coord_label2

  contains
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine load(filename)
  use grid
  character(len=*), intent(in)  :: filename

  type(t_grid) :: G


  ! cleanup
  if (allocated(x))  deallocate(x)
  if (allocated(x1)) deallocate(x1)
  if (allocated(x2)) deallocate(x2)
  if (allocated(x3)) deallocate(x3)


  ! load grid layout and nodes
  call G%load(filename)
  layout = G%layout
  coordinates = G%coordinates
  n = G%n
  allocate (x(n,3))
  x = G%x
  coord1 = G%coord1
  coord2 = G%coord2
  fixed_coord = G%fixed_coord

  n1 = G%n1
  n2 = G%n2
  n3 = G%n3
  coord_label1 = adjustl(G%coord_label(coord1))
  coord_label2 = adjustl(G%coord_label(coord2))

  ! set up structured grid
  if (layout == STRUCTURED) then
     allocate (x1(n1), x2(n2), x3(n3))
     x1 = G%x1
     x2 = G%x2
     if (G%fixed_coord == 0) then
        x3 = G%x3
     else
        x3 = G%fixed_coord_value
     endif

  ! set up semi-structured grid (structured in 1 coordinate)
  elseif (layout == SEMI_STRUCTURED) then
     allocate (x1(n1), x2(n2))
     x1 = G%x1
     x2 = G%x2

  ! set up 2D mesh
  elseif (layout == MESH_2D) then
     allocate (x1(n), x2(n))
     x1 = G%x(:,1)
     x2 = G%x(:,2)

  endif

  end subroutine load
  !---------------------------------------------------------------------

end module grid_interface
