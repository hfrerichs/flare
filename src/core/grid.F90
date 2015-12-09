!===============================================================================
! This module provides an interface between computations on grid nodes and the
! grid topology and coordinate system. Consequently, no information about the
! grid layout and coordinate system is required outside of this module.
!
! Grid nodes can be defined in Cartesian, cylindrical and toroidal/poloidal
! coordinates. The grid layout can be irregular/unstructured (1D array)
! or regular/structured (2D/3D array).
!
  ! grid id = IJK
  ! I: select internal coordinate system
  ! J: unstructured (1), semi-structured (2), structured (3), 2D mesh (4)
  ! K: select fixed or dominant coordinate

! Heinke Frerichs (hfrerichs at wisc.edu)
!===============================================================================
module grid
  use iso_fortran_env
  use math
#ifdef MPI
  use parallel
  implicit none
#else
  implicit none

  logical, parameter :: firstP = .true.
  integer            :: nprs   = 1
#endif
  private


  ! local coordinates
  integer, public, parameter :: &
     LOCAL      = 0
  public :: CYLINDRICAL, CARTESIAN

  ! grid layout
  integer, public, parameter :: &
     UNSTRUCTURED    = 1, &
     SEMI_STRUCTURED = 2, &
     STRUCTURED      = 3, &
     MESH_2D         = 4, &
     MESH_3D         = 5

  ! 2D vs. 3D grids
  integer, public, parameter :: &
     FIXED_COORD1    = 1, &
     FIXED_COORD2    = 2, &
     FIXED_COORD3    = 3, &
     TOROIDAL_SLICE  = FIXED_COORD3, &
     DEFAULT_GRID    = 0


  type, public :: t_grid
     ! grid nodes
     real(real64), dimension(:,:), allocatable :: x

     ! internal grid nodes
     real(real64), dimension(:),   allocatable :: x1, x2, x3
     real(real64), dimension(:,:,:), pointer :: mesh => null()
     real(real64), dimension(:,:,:,:), pointer :: mesh3D => null()

     ! number of nodes
     integer :: n, n1, n2, n3

     ! define internal layout and coordinate system
     integer :: layout, coordinates, fixed_coord, coord1, coord2
     real(real64) :: fixed_coord_value

     contains
     procedure :: new
     procedure :: load
     procedure :: load_usr
     procedure, private :: setup_mesh
     procedure, private :: setup_mesh3D
     procedure :: setup_structured_grid
     procedure :: store
     procedure :: destroy
     procedure :: plot_mesh
     procedure :: node                   ! return node coordinates
     procedure :: nodes                  ! return number of grid nodes
  end type t_grid

  contains
!=======================================================================



!=======================================================================
  subroutine new(this, coordinates, layout, fixed_coord, n1, n2, n3, fixed_coord_value, mesh)
  class(t_grid)                 :: this
  integer, intent(in)           :: coordinates, layout, fixed_coord, n1
  integer, intent(in), optional :: n2, n3
  real(real64), intent(in), optional :: fixed_coord_value
  logical, intent(in), optional :: mesh


  this%coordinates = coordinates
  this%n1          = n1
  this%n2          = 1
  this%n3          = 1
  this%n           = n1
  this%layout      = layout
  this%fixed_coord = fixed_coord
  if (present(n2)) then
     this%n2     = n2
     this%n      = this%n * n2
     !if (present(mesh) .and. mesh) then
     if (layout == MESH_2D) then
        if (associated(this%mesh)) deallocate(this%mesh)
        if (fixed_coord > 0) then
           allocate (this%mesh(0:n1-1, 0:n2-1, 2))
        else
           allocate (this%mesh(0:n1-1, 0:n2-1, 3))
        endif
        this%mesh = 0.d0
     endif
  endif
  if (present(n3)) then
     this%n3     = n3
     this%n      = this%n * n3
  endif


  ! setup coordinate indices
  select case(fixed_coord)
  case(1)
     this%coord1 = 2; this%coord2 = 3
  case(2)
     this%coord1 = 3; this%coord2 = 1
  case(0,3)
     this%coord1 = 1; this%coord2 = 2
  case default
     write (6, *) 'error: invalid id = ', fixed_coord, ' for fixed coordinate!'
     stop
  end select
  if (present(fixed_coord_value)) then
     this%fixed_coord_value = fixed_coord_value
  else
     this%fixed_coord_value = 0.d0
  endif


  ! allocate memory for grid nodes
  if (allocated(this%x)) deallocate(this%x)
  allocate (this%x(this%n,3))
  this%x = 0.d0
  if (fixed_coord > 0) this%x(:,fixed_coord) = this%fixed_coord_value

  select case (layout)
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
!  if (layout == STRUCTURED) then
  case(STRUCTURED)
     if (allocated(this%x1)) deallocate(this%x1)
     if (allocated(this%x2)) deallocate(this%x2)
     if (allocated(this%x3)) deallocate(this%x3)
     allocate (this%x1(this%n1))
     this%x1 = 0.d0
     allocate (this%x2(this%n2))
     this%x2 = 0.d0
     allocate (this%x3(this%n3))
     this%x3 = 0.d0
!  endif
!  if (layout == SEMI_STRUCTURED) then
  case(SEMI_STRUCTURED)
     if (allocated(this%x2)) deallocate(this%x2)
     allocate (this%x2(this%n2))
!  endif

  case(MESH_3D)
     if (associated(this%mesh3D)) deallocate(this%mesh3D)
     if (allocated(this%x3)) deallocate(this%x3)
     if (fixed_coord > 0) then
        allocate (this%mesh3D(0:this%n1-1, 0:this%n2-1, 0:this%n3-1, 2))
        allocate (this%x3(0:this%n3-1))
        this%x3 = 0.d0
     else
        allocate (this%mesh3D(0:this%n1-1, 0:this%n2-1, 0:this%n3-1, 3))
     endif
     this%mesh3D = 0.d0

  end select
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
  subroutine load(this, filename, silent)
  class(t_grid)                 :: this
  character(len=*), intent(in)  :: filename
  logical,          intent(in), optional :: silent

  integer, parameter :: iu = 32

  logical            :: screen_output


  ! set up internal variables
  screen_output = .true.
  if (present(silent) .and. silent) screen_output = .false.


  if (firstP) then
     call read_grid()
  endif
  if (nprs > 1) call broadcast_grid()

  contains
!-----------------------------------------------------------------------
  subroutine read_grid()

  character(len=120) :: str
  real(real64)       :: y3(3), r(3), y2(2), y1(1), x0, R0
  integer :: grid_id, coordinates, layout, fixed_coord, &
             i, j, k, ig, n, n1, n2, n3


  ! 1. open grid file
  if (screen_output) write (6,1000) adjustl(trim(filename))
  open  (iu, file=filename, err=5000)


  ! 2. read grid header and determine grid layout and coordinates
  read  (iu, 2000) str
  if (screen_output) write (6,2001) str(3:74)
  if (str(3:9).ne.'grid_id') then
     write (6,*) 'error: grid type not defined!'
     stop
  endif

  ! 2.1 collect information from grid id = IJK
  read  (str(13:16), '(i4)') grid_id
  coordinates = grid_id / 100                  ! I: select internal coordinate system
  grid_id     = grid_id - 100*coordinates
  layout      = grid_id / 10                   ! J: unstructured, semi-structured, structured
  fixed_coord = grid_id - 10*layout            ! K: select fixed or dominant coordinate

  ! 2.3 read reference parameters
  select case(coordinates)
  case(TORUS)
     call rscrape (iu, R0)
  end select


  ! 3. read grid resolution and grid nodes
  !.....................................................................
  ! 3.1a. unstructured 3D grid, n: total number of grid nodes
  if (layout == UNSTRUCTURED  .and.  fixed_coord == 0) then
     call iscrape (iu, n1)
     call this%new(coordinates, layout, fixed_coord, n1)

     ! read all grid nodes
     do i=1,n1
        read  (iu, *) y3
        this%x(i,:) = y3
     enddo

  !.....................................................................
  ! 3.1b. unstructured 2D grid, one coordinate is fixed, n: total number of grid nodes
  elseif (layout == UNSTRUCTURED  .and.  fixed_coord > 0) then
     call iscrape (iu, n)
     call this%new(coordinates, layout, fixed_coord, n)

     ! read fixed coordinate
     call rscrape (iu, x0)
     this%fixed_coord_value = x0

     ! read all grid nodes
     do i=1,n
        read  (iu, *) y2
        this%x(i, this%coord1) = y2(1)
        this%x(i, this%coord2) = y2(2)
        this%x(i, fixed_coord) = x0
     enddo

  !.....................................................................
  ! 3.2. semi-structured grid, list of (x(coord1), x(coord2)), then list of x(fixed_coord)
  elseif (layout == SEMI_STRUCTURED) then
     if (fixed_coord == 0) then
        write (6, *) 'error: fixed_coord > 0 required for semi-structured grids!'
        stop
     endif
     call iscrape (iu, n1)
     call iscrape (iu, n2)
     call this%new(coordinates, layout, fixed_coord, n1, n2)

     ! read list of (x(coord1), x(coord2))
     do i=1,n1
        read  (iu, *) y2

        ! set y2 for all n2 values of x(fixed_coord)
        do j=1,n2
           this%x((j-1)*n1 + i, this%coord1) = y2(1)
           this%x((j-1)*n1 + i, this%coord2) = y2(2)
        enddo
     enddo

     ! read list of x(fixed_coord)
     do j=1,n2
        read  (iu, *) y1
        this%x2(j) = y1(1)

        ! set y1 for all n1 values of x(coord1), x(coord2)
        do i=1,n1
           this%x((j-1)*n1 + i, fixed_coord) = y1(1)
        enddo
     enddo

  !.....................................................................
  ! 3.3. structured grid, separate list of coordinates
  elseif (layout == STRUCTURED) then
     call iscrape (iu, n1)
     call iscrape (iu, n2)
     ! all coordinates are structured
     if (fixed_coord == 0) then
        call iscrape (iu, n3)
     ! one coordinate is fixed
     else
        n3 = 1
        call rscrape (iu, x0)
        this%fixed_coord_value = x0
     endif

     call this%new(coordinates, layout, fixed_coord, n1, n2, n3)

     ! read 1st coordinates
     do i=1,n1
        read  (iu, *) this%x1(i)
     enddo

     ! read 2nd coordinates
     do j=1,n2
        read  (iu, *) this%x2(j)
     enddo

     if (fixed_coord == 0) then
        ! read 3rd coordinates
        do k=1,n3
           read  (iu, *) this%x3(k)
        enddo
     endif

     ! distribute to main grid
     call this%setup_structured_grid()

  !.....................................................................
  ! 3.4. unstructured mesh, n1*n2: total number of grid nodes
  elseif (layout == MESH_2D  .and.  fixed_coord > 0) then
     call iscrape (iu, n1)
     call iscrape (iu, n2)
     !call this%new(coordinates, layout, fixed_coord, n1, n2, mesh=.true.)
     call this%new(coordinates, layout, fixed_coord, n1, n2)

     ! read fixed coordinate
     call rscrape (iu, x0)
     this%fixed_coord_value = x0

     ! read all grid nodes
     do j=0,n2-1
        do i=0,n1-1
           read  (iu, *) y2
           this%mesh(i,j,1) = y2(1)
           this%mesh(i,j,2) = y2(2)
           !this%mesh(i,j,this%coord1) = y2(1)
           !this%mesh(i,j,this%coord2) = y2(2)
           !this%mesh(i,j,fixed_coord) = x0
           !this%x(i, this%coord1) = y2(1)
           !this%x(i, this%coord2) = y2(2)
           !this%x(i, fixed_coord) = x0
        enddo
     enddo
     call this%setup_mesh()

  !.....................................................................
  ! 3.5. unstructured 3D mesh, n1*n2*n3: total number of grid nodes
  elseif (layout == MESH_3D  .and.  fixed_coord == 3) then
     call iscrape (iu, n1)
     call iscrape (iu, n2)
     call iscrape (iu, n3)
     call this%new(coordinates, layout, fixed_coord, n1, n2, n3)

     ! read all grid nodes
     do k=0,n3-1
        read  (iu, *) this%x3(k)
        do j=0,n2-1
        do i=0,n1-1
           read  (iu, *) y2
           this%mesh3D(i,j,k,1:2) = y2
        enddo
        enddo
     enddo
     call this%setup_mesh3D()

  !.....................................................................
  else
     write (6, *) 'error: invalid grid layout ', layout, '!'
     stop
  endif


  ! 4. close grid file
  close (iu)


  ! 5. convert to cylindrical coordinates
  select case(coordinates)
  case(CARTESIAN)
     do i=1,this%n
        y3          = this%x(i,:)
        call coord_trans(y3, CARTESIAN, r, CYLINDRICAL)
        this%x(i,:) = r
     enddo

  case(CYLINDRICAL)
     ! deg -> rad
     this%x(:,3) = this%x(:,3) / 180.d0 * pi
     if (fixed_coord == 3) this%fixed_coord_value = this%fixed_coord_value / 180.d0 * pi

  case(TORUS)
     ! poloidal, toroidal angle: deg -> rad
     this%x(:,2:3) = this%x(:,2:3) / 180.d0 * pi

     do i=1,this%n
        y3          = this%x(i,:)
        call coord_trans_torus (y3, R0, r)
        this%x(i,:) = r
     enddo
  case default
  end select

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

#ifdef MPI
  call wait_pe()
  call broadcast_inte_s(this%n)
  call broadcast_inte_s(this%coordinates)
  call broadcast_inte_s(this%layout)
  call broadcast_inte_s(this%fixed_coord)
  if (mype > 0) then
     allocate (this%x(this%n, 3))
  endif
  call broadcast_real  (this%x, this%n*3)
#endif

  end subroutine broadcast_grid
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! internal routine to read integer values from header
!-----------------------------------------------------------------------
  subroutine iscrape (iu, iout)
  integer, intent(in)  :: iu
  integer, intent(out) :: iout
  character(len=82) :: str
  read  (iu, 2000) str
  if (screen_output) write (6,2001) str(3:82)
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
  if (screen_output) write (6,2001) str(3:82)
  read  (str(33:42),*) rout
 2000 format (a82)
 2001 format (8x,a80)
  end subroutine rscrape
!-----------------------------------------------------------------------
  end subroutine load
!=======================================================================



!=======================================================================
! iformat = 1:	X, Y, Z [cm]
!           2:  R, Z [cm], phi [deg]
!           3:  R, Z [cm] at phi_default
!=======================================================================
  subroutine load_usr(this, filename, iformat, phi0)
#ifdef FLARE
  use dataset
#endif
  class(t_grid)                :: this
  character(len=*), intent(in) :: filename
  integer,          intent(in) :: iformat
  real(real64),     intent(in), optional :: phi0

  integer, parameter :: iu = 10

#ifdef FLARE
  type(t_dataset) :: D
  real(real64)    :: phi, yi(3), yo(3)
  integer         :: columns, i, n


  phi = 0.d0
  if (present(phi0)) phi = phi0


  ! read data
  columns = 3
  if (iformat == 3) columns = 2
  call D%load(filename, columns)


  ! setup grid
  n = D%nrow
  call this%new(CYLINDRICAL, UNSTRUCTURED, DEFAULT_GRID, n)
  do i=1,n
     yi(1) = D%x(i,1)
     yi(2) = D%x(i,2)
     select case (iformat)
     case(1,2)
        yi(3) = D%x(i-1,3) / 180.d0 * pi
     case(3)
        yi(3) = phi / 180.d0 * pi
     end select
     call coord_trans (yi, iformat, yo, CYLINDRICAL)
     this%x(i,:) = yo
  enddo


  ! cleanup
  call D%destroy()

#endif
  end subroutine load_usr
!=======================================================================



!=======================================================================
  subroutine setup_structured_grid(this)
  class(t_grid)                 :: this

  integer :: i, j, k, ig


  ! distribute to main grid
  do i=1,this%n1
  do j=1,this%n2
  do k=1,this%n3
     ig = (k-1)*this%n1*this%n2  +  (j-1)*this%n1  +  i
     this%x(ig, this%coord1)         = this%x1(i)
     this%x(ig, this%coord2)         = this%x2(j)
     if (this%fixed_coord == 0) then
        this%x(ig, 3)                = this%x3(k)
     else
        this%x(ig, this%fixed_coord) = this%fixed_coord_value
     endif
  enddo
  enddo
  enddo

  end subroutine setup_structured_grid
!=======================================================================



!=======================================================================
  subroutine setup_mesh(this)
  class(t_grid)                 :: this

  integer :: i, j, ig


  if (.not.associated(this%mesh)) then
     write (6, *) 'error in subroutine t_grid%setup_mesh: mesh undefined!'
     stop
  endif


  ig = 0
  do j=0,this%n2-1
     do i=0,this%n1-1
        ig = ig + 1
        this%x(ig,this%coord1)      = this%mesh(i,j,1)
        this%x(ig,this%coord2)      = this%mesh(i,j,2)
        if (this%fixed_coord > 0) then
           this%x(ig,this%fixed_coord) = this%fixed_coord_value
        else
           this%x(ig,3)                = this%mesh(i,j,3)
        endif
     enddo
  enddo

  end subroutine setup_mesh
!=======================================================================



!=======================================================================
  subroutine setup_mesh3D(this)
  class(t_grid)                 :: this

  integer :: i, j, k, ig


  if (.not.associated(this%mesh3D)) then
     write (6, *) 'error in subroutine t_grid%setup_mesh3D: mesh3D undefined!'
     stop
  endif


  ig = 0
  do k=0,this%n3-1
  do j=0,this%n2-1
  do i=0,this%n1-1
     ig = ig + 1
     this%x(ig,this%coord1)      = this%mesh3D(i,j,k,1)
     this%x(ig,this%coord2)      = this%mesh3D(i,j,k,2)
     if (this%fixed_coord > 0) then
        this%x(ig,this%fixed_coord) = this%x3(k)
     else
        this%x(ig,3)                = this%mesh3D(i,j,k,3)
     endif
  enddo
  enddo
  enddo

  end subroutine setup_mesh3D
!=======================================================================



!=======================================================================
  subroutine store(this, filename, header)
  class(t_grid)                :: this
  character(len=*), intent(in) :: filename
  character(len=*), intent(in), optional :: header

  integer, parameter :: iu = 32

  real(real64), dimension(:,:), allocatable :: xout
  real(real64) :: phi, fixed_coord_value_out
  integer :: grid_id, i, j, k, layout, coord1, coord2


! open output file .............................................
  open  (iu, file=filename)


! write header .................................................
  grid_id = this%coordinates * 100  +  this%layout * 10  +  this%fixed_coord
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
  end select
 1000 format ('# grid_id = ', i4, 4x, '(', a, ')')


! convert output units
  fixed_coord_value_out = this%fixed_coord_value
  allocate (xout(this%n, 3))
  xout = this%x
  if (this%coordinates == CYLINDRICAL) then
     ! rad -> deg
     if (this%fixed_coord == 3) then
        fixed_coord_value_out = fixed_coord_value_out / pi * 180.d0 
     endif
     xout(:,3) = xout(:,3) * 180.d0 / pi
  endif


! write grid nodes .............................................
  layout = this%layout
  ! 1.a unstructured 3D grids
  if (layout == UNSTRUCTURED  .and.  this%fixed_coord == 0) then
     write (iu, 2001) this%n
     do i=1,this%n
        write (iu, 3003) xout(i,:)
     enddo

  ! 1.b unstructured 2D grids, one coordinate fixed
  elseif (layout == UNSTRUCTURED  .and.  this%fixed_coord > 0) then
     write (iu, 2001) this%n
     write (iu, 2002) COORD_STR(this%fixed_coord, this%coordinates), fixed_coord_value_out
     do i=1,this%n
        write (iu, 3002) xout(i, this%coord1), xout(i, this%coord2)
     enddo

  ! 2. semi-structured grid, list of (x(coord1), x(coord2)), then list of x(fixed_coord)
  elseif (layout == SEMI_STRUCTURED) then
     write (iu, 2003) this%n1
     write (iu, 2004) this%n2
     do i=1,this%n1
        write (iu, 3002) xout(i, this%coord1), xout(i, this%coord2)
     enddo
     do j=1,this%n2
        write (iu, 3002) xout((j-1)*this%n1+1, this%fixed_coord)
     enddo

  ! 3. structured grid, list of (x(coord1), x(coord2)), then list of x(fixed_coord)
  elseif (layout == STRUCTURED) then
     ! header
     write (iu, 2003) this%n1
     write (iu, 2004) this%n2
     if (this%fixed_coord == 0) then
        write (iu, 2005) this%n3
     else
        write (iu, 2002) COORD_STR(this%fixed_coord, this%coordinates), fixed_coord_value_out
     endif

     ! body
     do i=1,this%n1
        write (iu, 3002) this%x1(i)
     enddo
     do i=1,this%n2
        write (iu, 3002) this%x2(i)
     enddo
     if (this%fixed_coord == 0) then
        do i=1,this%n3
           if (this%coordinates == CYLINDRICAL) then
              write (iu, 3002) this%x3(i) / pi * 180.d0
           else
              write (iu, 3002) this%x3(i)
           endif
        enddo
     endif

  ! 4. unstructured meshs, one coordinate fixed
  elseif (layout == MESH_2D  .and.  this%fixed_coord > 0) then
     write (iu, 2003) this%n1
     write (iu, 2004) this%n2
     write (iu, 2002) COORD_STR(this%fixed_coord, this%coordinates), fixed_coord_value_out
     do j=0,this%n2-1
     do i=0,this%n1-1
        write (iu, 3002) this%mesh(i,j,1), this%mesh(i,j,2)
     enddo
     enddo

  ! 5. unstructured 3D meshs (set of 2D slices)
  elseif (layout == MESH_3D  .and.  this%fixed_coord == 3) then
     write (iu, 2003) this%n1
     write (iu, 2004) this%n2
     write (iu, 2005) this%n3
     do k=0,this%n3-1
        write (iu, 3001) this%x3(k) / pi * 180.d0
        do j=0,this%n2-1
        do i=0,this%n1-1
           write (iu, 3002) this%mesh3D(i,j,k,1), this%mesh3D(i,j,k,2)
        enddo
        enddo
     enddo

  else
     write (6, *) 'to be implemented ...'
     stop
  endif


! close output file ............................................
  close (iu)
  deallocate (xout)

 1100 format ('# grid_id = 1000    (cylindrical coordinates: R[cm], Z[cm], Phi[rad])')
 2001 format ('# grid resolution:   n_grid  =  ',i10)
 2002 format ('# ',a12,' =                ',f10.5)
 2003 format ('# grid resolution:   n1      =  ',i10)
 2004 format ('#                    n2      =  ',i10)
 2005 format ('#                    n3      =  ',i10)
 3001 format (1e18.10)
 3002 format (2e18.10)
 3003 format (3e18.10)
  end subroutine store
!=======================================================================



!=======================================================================
  subroutine destroy(this)
  class(t_grid)                :: this


  if (allocated(this%x))  deallocate(this%x)
  if (allocated(this%x1)) deallocate(this%x1)
  if (allocated(this%x2)) deallocate(this%x2)
  if (allocated(this%x3)) deallocate(this%x3)
  if (associated(this%mesh)) deallocate(this%mesh)
  if (associated(this%mesh3D)) deallocate(this%mesh3D)

  end subroutine destroy
!=======================================================================



!=======================================================================
  subroutine plot_mesh(this, filename)
  class(t_grid)                :: this
  character(len=*), intent(in) :: filename

  integer, parameter :: iu = 42

  integer :: i, j


  open  (iu, file=filename)
  ! write rows
  do i=0,this%n1-1
     do j=0,this%n2-1
        write (iu, *) this%mesh(i,j,:)
     enddo
     write (iu, *)
  enddo

  ! write columns
  do j=0,this%n2-1
     do i=0,this%n1-1
        write (iu, *) this%mesh(i,j,:)
     enddo
     write (iu, *)
  enddo
  close (iu)

  end subroutine plot_mesh
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
     call coord_trans (this%x(i,:), CYLINDRICAl, x, coordinates)
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
