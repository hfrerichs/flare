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
     UNSTRUCTURED_3D = 6, &
     UNSTRUCTURED_2D = 5
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
     procedure :: node                   ! return node coordinates
     procedure :: nodes                  ! return number of grid nodes
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


  ! allocate memory for grid nodes
  if (allocated(this%x)) deallocate(this%x)
  allocate (this%x(this%n,3))

!  select case (layout)
!  case(UNSTRUCTURED_3D)
!     ! unstructured grid nodes (3D)
!     if (allocated(this%x)) deallocate(this%x)
!     allocate (this%x(this%n,3))
!  case(UNSTRUCTURED_2D)
!     ! unstructured grid nodes (2D)
!     if (allocated(this%x)) deallocate(this%x)
!     allocate (this%x(this%n,2))
!
!     ! hyperplane reference value
!     if (allocated(this%x3)) deallocate(this%x3)
!     allocate (this%x3(1))
!     this%n3 = 1
!  case(STRUCTURED_2D)
!     if (allocated(this%x1)) deallocate(this%x1)
!     if (allocated(this%x2)) deallocate(this%x2)
!     if (allocated(this%x3)) deallocate(this%x3)
!     allocate (this%x1(this%n1))
!
!     if (.not.present(n2)) then
!        write (6, *) 'error: resolution parameter n2 missing for regular 2D grid!'
!        stop
!     endif
!     allocate (this%x2(this%n2))
!     allocate (this%x3(      1))
!  end select


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

  integer, parameter :: iu = 32


  if (firstP) then
     call read_grid()
  endif
  call broadcast_grid()

  contains
!-----------------------------------------------------------------------
  subroutine read_grid

  character(len=120) :: str
  real(real64)       :: y3(3), y2(2), y1(1), x0
  integer :: grid_id, coordinates, layout, fixed_coord, coord1, coord2, i, n


  ! 1. open grid file
  write (6,1000) adjustl(trim(filename))
  open  (iu, file=filename, err=5000)


  ! 2. read grid header and determine grid layout and coordinates
  read  (iu, 2000) str
  write (6,2001) str(3:74)
  if (str(3:9).ne.'grid_id') then
     write (6,*) 'error: grid type not defined!'
     stop
  endif
  ! 2.1 collect information from grid id
  read  (str(13:16), '(i4)') grid_id
  coordinates = grid_id / 100
  grid_id     = grid_id - 100*coordinates
  layout      = grid_id / 10
  fixed_coord = grid_id - 10*layout
  ! 2.2 setup coordinate indices
  select case(fixed_coord)
  case(0)
  case(1)
     coord1 = 2; coord2 = 3
  case(2)
     coord1 = 1; coord2 = 3
  case(3)
     coord1 = 1; coord2 = 2
  case default
     write (6, *) 'error: invalid id = ', fixed_coord, ' for fixed coordinate!'
     stop
  end select


  ! 3. read grid resolution and grid nodes
  select case(layout)
  !.....................................................................
  case (UNSTRUCTURED_3D)
     ! unstructured 3D grid, n: total number of grid nodes
     call iscrape (iu, n)
     call this%new(coordinates, layout, n)

     ! read all grid nodes and convert to cylindrical coordinates
     do i=1,n
        read  (iu, *) y3

        call coord_trans (y3, coordinates, y3, CYLINDRICAL)
        this%x(i,:) = y3
     enddo

  !.....................................................................
  case (UNSTRUCTURED_2D)
     ! unstructured 2D grid, one coordinate is fixed, n: total number of grid nodes
     call iscrape (iu, n)
     call this%new(coordinates, layout, n)

     ! read fixed coordinate
     call rscrape (iu, x0)

     ! read all grid nodes
     do i=1,n
        read  (iu, *) y2
        this%x(i,coord1)      = y2(1)
        this%x(i,coord2)      = y2(2)
        this%x(i,fixed_coord) = x0
     enddo
  !.....................................................................
  case default
     write (6, *) 'error: invalid grid layout ', this%layout, '!'
     stop
  end select


  ! 4. close grid file
  close (iu)


  ! 5. convert to cylindrical coordinates




  return
 1000 format (3x,'- Using grid file: ',a)
! 1001 format ('progress:             ',i4,' %')
 2000 format (a120)
 2001 format (8x,a72)
 5000 write (6,5001) filename
 5001 format ('error reading grid file: ', a120)
      stop
  end subroutine read_grid
!-----------------------------------------------------------------------
  subroutine broadcast_grid
  end subroutine broadcast_grid
!-----------------------------------------------------------------------
  end subroutine load
!-----------------------------------------------------------------------
! internal routine to read integer values from header
!-----------------------------------------------------------------------
  subroutine iscrape (iu, iout)
  integer, intent(in)  :: iu
  integer, intent(out) :: iout
  character(len=82) :: str
  read  (iu, 2000) str
  write (6,2001) str(3:82)
  read  (str(33:42),*) iout
 2000 format (a82)
 2001 format (8x,a80)
  end subroutine iscrape
!-----------------------------------------------------------------------
! internal routine to read real values from header
!-----------------------------------------------------------------------
  subroutine rscrape (iu, rout)
  integer,      intent(in)  :: iu
  real(real64), intent(out) :: rout
  character(len=82) :: str
  read  (iu, 2000) str
  write (6,2001) str(3:82)
  read  (str(33:42),*) rout
 2000 format (a82)
 2001 format (8x,a80)
  end subroutine rscrape
!-----------------------------------------------------------------------
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
  case default
     write (6, *) 'to be implemented ...'
     stop
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



!=======================================================================
! return coordinates of node i
! cylindrical coordinates R[cm], Z[cm], Phi[rad] are used by default
!=======================================================================
  function node(this, i, coordinates) result(x)
  class(t_grid)       :: this
  integer, intent(in) :: i
  integer, intent(in), optional :: coordinates
  real(real64)        :: x(3)


  x = this%x(i,:)

  if (present(coordinates)) then
     call coord_trans (x, CYLINDRICAl, x, coordinates)
  endif

  end function node
!=======================================================================



!=======================================================================
  function nodes(this) result(n)
  class(t_grid)       :: this
  integer             :: n

  n = this%n
  end function nodes
!=======================================================================

end module grid
