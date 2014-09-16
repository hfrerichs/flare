!===============================================================================
! This module provides an interface between computations on grid nodes and the
! grid topology and coordinate system. Consequently, no information about the
! grid layout and coordinate system is required outside of this module.
!
! Grid nodes can be defined in Cartesian, cylindrical and toroidal/poloidal
! coordinates. The grid layout can be irregular/unstructured (1D array)
! or regular/structured (2D/3D array).
!
! Heinke Frerichs (hfrerichs at wisc.edu)
!===============================================================================
module grid
  use iso_fortran_env
  use math
  use parallel
  implicit none
  private


  ! grid layout
  integer, parameter :: &
     IRREGULAR  = 1, &         ! (x1_i, x2_i, x3_i), i=1..n
     REGULAR_2D = 2, &
     REGULAR_3D = 3            ! (x1_i, x2_j, x3_k), i=1..n1, j=1..n2, k=1..n3


  type t_grid
     ! grid nodes
     real(real64), dimension(:,:), allocatable :: x

     ! internal grid nodes
     real(real64), dimension(:),   allocatable :: x1, x2, x3

     ! number of nodes
     integer :: n, n1, n2, n3

     ! define internal layout and coordinate system
     integer :: layout, coordinates

     contains
     procedure :: new
     procedure :: load
  end type t_grid

  contains
!=======================================================================



!=======================================================================
  subroutine new(this, layout, n1, n2, n3)
  class(t_grid)                 :: this
  integer, intent(in)           :: layout, n1
  integer, intent(in), optional :: n2, n3


  this%layout = layout

  end subroutine new
!=======================================================================



!=======================================================================
  subroutine load(this, filename)
  class(t_grid)                :: this
  character(len=*), intent(in) :: filename
  end subroutine load
!=======================================================================

end module grid
