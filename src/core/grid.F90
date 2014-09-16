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


  ! local coordinates
  integer, public, parameter :: &
     LOCAL      = 0            ! local coordinates, for plotting only

  ! grid layout
  integer, public, parameter :: &
     STRUCTURED_2D   = 2, &
     STRUCTURED_3D   = 3, &
     UNSTRUCTURED_3D = 30, &
     UNSTRUCTURED_2D = 20
     !IRREGULAR  = 1, &         ! (x1_i, x2_i, x3_i), i=1..n
           ! (x1_i, x2_j, x3_k), i=1..n1, j=1..n2, k=1..n3


  type, public :: t_grid
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
     procedure :: new_unstructured
     procedure :: load
     procedure :: store
  end type t_grid

  contains
!=======================================================================



!=======================================================================
  subroutine new(this, coordinates, layout, n1, n2, n3)
  class(t_grid)                 :: this
  integer, intent(in)           :: coordinates, layout, n1
  integer, intent(in), optional :: n2, n3


  this%coordinates = coordinates
  this%n1          = n1
  this%n           = n1
  !this%layout      = UNSTRUCTURED_3D
  this%layout      = layout
  if (present(n2)) then
     this%n2     = n2
     this%n      = this%n * n2
  endif
  if (present(n3)) then
     this%n3     = n3
     this%n      = this%n * n3
  endif


  select case (layout)
  case(UNSTRUCTURED_3D)
     ! unstructured grid nodes (3D)
     if (allocated(this%x)) deallocate(this%x)
     allocate (this%x(this%n,3))
  case(UNSTRUCTURED_2D)
     ! unstructured grid nodes (2D)
     if (allocated(this%x)) deallocate(this%x)
     allocate (this%x(this%n,2))

     ! hyperplane reference value
     if (allocated(this%x3)) deallocate(this%x3)
     allocate (this%x3(1))
     this%n3 = 1
  case(STRUCTURED_2D)
     if (allocated(this%x1)) deallocate(this%x1)
     if (allocated(this%x2)) deallocate(this%x2)
     if (allocated(this%x3)) deallocate(this%x3)
     allocate (this%x1(this%n1))

     if (.not.present(n2)) then
        write (6, *) 'error: resolution parameter n2 missing for regular 2D grid!'
        stop
     endif
     allocate (this%x2(this%n2))
     allocate (this%x3(      1))
  end select


  end subroutine new
!=======================================================================



!=======================================================================
  subroutine new_unstructured(this, coordinates, n)
  class(t_grid)                 :: this
  integer, intent(in)           :: coordinates, n


  this%coordinates = coordinates
  !this%layout      = UNSTRUCTURED_3D
  this%layout      = 100 * coordinates
  if (allocated(this%x)) deallocate(this%x)
  allocate (this%x(n,3))

  end subroutine new_unstructured
!=======================================================================



!=======================================================================
  subroutine load(this, filename)
  class(t_grid)                 :: this
  character(len=*), intent(in)  :: filename
  end subroutine load
!=======================================================================



!=======================================================================
  subroutine store(this, filename, header)
  class(t_grid)                :: this
  character(len=*), intent(in) :: filename
  character(len=*), intent(in), optional :: header

  integer, parameter :: iu = 32

  integer :: grid_id, i


! open output file .............................................
  open  (iu, file=filename)


! write header .................................................
  grid_id = this%coordinates * 100  +  this%layout
  select case(this%coordinates)
  case(CARTESIAN)
     write (iu, 1000) grid_id, 'Cartesian coordinates: x[cm], y[cm], z[cm]'
  case(CYLINDRICAL)
     write (iu, 1000) grid_id, 'cylindrical coordinates: R[cm], Z[cm], Phi[deg]'
  case(LOCAL)
     if (present(header)) then
        write (iu, 1000) grid_id, trim(header)
     else
        write (iu, 1000) grid_id, ''
     endif
     !write (iu, *)
     !write (iu, *)
  end select
 1000 format ('# grid_id = ', i4, 4x, '(', a, ')')


!  select case(this%layout)
!  ! unstructured 3D grids
!  case(UNSTRUCTURED_3D)
!     select case(this%coordinates)
!     case(CARTESIAN)
!        write (iu, 1001)
!     case(CYLINDRICAL)
!     end select
!  end select

!  select case(this%layout)
!  ! 1. RZ-slice (irregular) at constant phi
!  case(1)
!  ! 1000. Cylindrical coordinates (irregular)
!  case(100)
!     write (iu, 1100)
!     write (iu, 0001) this%n
!  end select


! write grid nodes .............................................
  select case(this%layout)
  ! unstructured 3D grids
  case(UNSTRUCTURED_3D)
     write (iu, 2001) this%n
     do i=1,this%n
        write (iu, 3003) this%x(i,:)
     enddo
  case(UNSTRUCTURED_2D)
     write (iu, 2001) this%n
     do i=1,this%n
        write (iu, 3002) this%x(i,:)
     enddo
  end select


! close output file ............................................
  close (iu)

 1100 format ('# grid_id = 1000    (cylindrical coordinates: R[cm], Z[cm], Phi[rad])')
 2001 format ('# grid resolution:   n_grid  =  ',i10)
 3001 format (1e18.10)
 3002 format (2e18.10)
 3003 format (3e18.10)
  end subroutine store
!=======================================================================

end module grid
