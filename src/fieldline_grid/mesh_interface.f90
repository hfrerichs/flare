module mesh_interface
!module zone_interface
  use iso_fortran_env
  use curve2D
  implicit none
  private


  integer, parameter, public :: &
     STRIKE_POINT = 0


  type, public :: t_mesh_interface
     integer       :: id

     ! geometric definition of interface
     type(t_curve) :: C
     integer       :: inode(-1:1) ! lower and upper end node type:  > 0 X-point
                                  !                                 = 0 strike point
                                  !                                 < guiding point (other X-point)

     ! adjacent zones:
     ! Z(-1) = zone on lower side of interface, i.e. interface is upper zone boundary
     ! Z( 1) = zone on upper side of interface, i.e. interface is lower zone boundary
     integer :: Z(-1,1) = -1


     ! discretization of interface
     real(real64), dimension(:,:), allocatable :: x
     ! number of nodes
     integer :: n

     contains
     procedure :: setup
     procedure :: set_curve
     procedure :: setup_discretization
  end type t_mesh_interface
  !@public :: setup_interfaces

  ! radial interfaces between zones
  type(t_mesh_interface), dimension(:), allocatable, public :: radial_interface
  integer, public :: radial_interfaces
  ! poloidal interfaces between zones
  type(t_mesh_interface), dimension(:), allocatable, public :: poloidal_interface
  integer, public :: poloidal_interfaces


  ! number of X-points in computation domain, and their connectivity
  integer, dimension(:), allocatable :: connectX
  integer :: nX




  public :: initialize_interfaces

  contains
!=======================================================================



!=======================================================================
  subroutine setup(this, lower_boundary, upper_boundary)
  use equilibrium, only: Xp, nx_max
  class(t_mesh_interface) :: this
  integer, intent(in)     :: lower_boundary, upper_boundary

  real(real64) :: X(2)
  integer      :: ix, iside, inode(-1:1), ierr


  inode(-1) = lower_boundary
  inode( 1) = upper_boundary

  do iside=-1,1,2
     ix = inode(iside)
     ! check input
     if (abs(ix) > nx_max) then
        write (6, 9001) lower_boundary, upper_boundary
        stop
     endif

     ! boundary is X-point or guiding point from another X-point
     if (abs(ix) > 0) then
        X  = Xp(abs(ix))%load(ierr)
        if (ierr .ne. 0) then
           write (6, 9002) abs(ix)
           stop
        endif
        this%inode(iside) = ix
   
     ! lower boundary is divertor target
     else
        this%inode(iside) = ix
     endif
  enddo
   

  ! 2. generate interface geometry
  ! call F%generate_branch
  ! for stability: connect 2 X-points by joining curves!


 9001 format('error in t_mesh_interface%setup: invalid X-point IDs = ', 2i0)
 9002 format('error in t_mesh_interface%setup: X-point ', i0, 'is not defined!')
  end subroutine setup
!=======================================================================



!=======================================================================
  subroutine set_curve(this, C)
  class(t_mesh_interface)   :: this
  type(t_curve), intent(in) :: C


  call C%copy_to(this%C)

  end subroutine set_curve
!=======================================================================



!=======================================================================
  subroutine setup_discretization(this, n, spacings)
  use mesh_spacing
  class(t_mesh_interface)   :: this
  integer,       intent(in) :: n
  type(t_spacing)           :: spacings

  real(real64) :: tau, x(2)
  integer      :: i


  this%n = n
  do i=0,n
     tau = spacings%node(i,n)
     call this%C%sample_at(tau, x)
     this%x(i,1:2) = x
  enddo

  end subroutine setup_discretization
!=======================================================================



!=======================================================================
  subroutine initialize_interfaces(n)
  integer, intent(in) :: n


  allocate (radial_interface(n))
  radial_interfaces = n

  end subroutine initialize_interfaces
!=======================================================================



!=======================================================================
  subroutine setup_topology()
  use fieldline_grid


  select case(topology)
  ! lower single null (LSN)
  case(TOPO_LSN, TOPO_LSN1)
     call initialize_interfaces(3) ! radial interfaces
     nX = 1;  allocate(connectX(nX))
     connectX(1) = 1


  ! DDN
  case(TOPO_DDN, TOPO_DDN1)
     call initialize_interfaces(8) ! radial interfaces
     nX = 2;  allocate(connectX(nX))
     connectX(1) = -2
     connectX(2) = -2


  ! CDN
  case(TOPO_CDN, TOPO_CDN1)
     !call initialize_interfaces(3) ! radial interfaces
     nX = 2;  allocate(connectX(nX))
     connectX(1) = 2
     connectX(2) = 1

  case default
     write (6, 9000) trim(topology)
     stop
  end select


  call setup_zones()

 9000 format('error: invalid topology ', a, '!')
  end subroutine setup_topology
!=======================================================================



!=======================================================================
  subroutine setup_zones

  integer :: i, iz

  ! set up zone numbers for interfaces
  iz = 1
  do i=1,radial_interfaces
  enddo

  end subroutine setup_zones
!=======================================================================



!=======================================================================
!  subroutine setup_interfaces(nx, connectX)
!  use equilibrium, only: Xp
!  use math
!  use separatrix
!  !use fieldline_grid, only: C_guide
!  integer, intent(in) :: nx, connectX(nx)
!
!  type(t_separatrix)  :: S(nx)
!  real(real64)        :: X(2)
!  integer             :: ix, ierr, jx
!
!
!  ! 1. check definition of X-points
!  do ix=1,nX
!     X = Xp(ix)%load(ierr)
!     if (ierr .ne. 0) then
!        write (6, 9001) ix
!        stop
!     endif
!     write (6, 2001) ix, X
!     write (6, 2002) Xp(ix)%theta/pi*180.d0
!  enddo
!  write (6, *)
! 2001 format(8x,i0,'. X-point at: ',2f10.4)
! 2002 format(11x,'-> poloidal angle [deg]: ',f10.4)
! 9001 format('error: ',i0,'. X-point is not defined!')
!
!
!  ! 2. generate separatrix for each X-point
!  do ix=1,nX
!     jx = connectX(ix)
!     !call S(ix)%generate(ix, C_cutl=C_guide, C_cutr=C_guid, iconnect=jx)
!     call S(ix)%generate(ix, iconnect=jx)
!  enddo
!
!  end subroutine setup_interfaces
!=======================================================================

end module mesh_interface
