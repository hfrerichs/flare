module mesh_interface
  use iso_fortran_env
  use curve2D
  implicit none
  private


  type, public :: t_mesh_interface
     ! geometric definition of interface
     type(t_curve) :: C
     integer       :: inode(-1:1) ! lower and upper end node type:  > 0 X-point
                                  !                                <= 0 strike point


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

  ! interfaces between zones
  type(t_mesh_interface), dimension(:), allocatable, public :: Iface
  integer, public :: interfaces

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
     if (ix > nx_max) then
        write (6, 9001) lower_boundary, upper_boundary
        stop
     endif

     ! boundary is X-point
     if (ix > 0) then
        X  = Xp(ix)%load(ierr)
        if (ierr .ne. 0) then
           write (6, 9002) ix
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


  !call this%C%copy(C)

  call C%copy_to(this%C)

  end subroutine set_curve
!=======================================================================



!=======================================================================
  subroutine setup_discretization(this, C, n)
  class(t_mesh_interface)   :: this
  type(t_curve), intent(in) :: C
  integer,       intent(in) :: n
  end subroutine setup_discretization
!=======================================================================



!=======================================================================
  subroutine initialize_interfaces(n)
  integer, intent(in) :: n


  allocate (Iface(n))
  interfaces = n

  end subroutine initialize_interfaces
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
