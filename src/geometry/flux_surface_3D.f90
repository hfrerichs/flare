!===============================================================================
! Generate flux surfaces
!===============================================================================
module flux_surface_3D
  use iso_fortran_env
  use flux_surface_2D
  use dataset
  use curve2D
  use poincare_set
  implicit none

  private


  ! pre-sampled distance to flux surface
  type t_distance
     ! raw distance on nodes
     real(real64), dimension(:,:,:), allocatable :: d

     ! coefficients for spline interpolation
     real(real64), dimension(:,:,:), allocatable :: dcoeff
     real(real64), dimension(:),     allocatable :: Rnot, Znot, Phinot

     ! Rc: center major radius, w: width, h: height
     real(real64) :: Rc, w, h
     integer      :: nr, nz, nphi, nsym

     ! interpolation order
     integer      :: nord

     contains
     procedure    :: new      => new_distance
     procedure    :: setup    => setup_distance
     procedure    :: evaluate => evaluate_distance
  end type t_distance


  ! flux surface contour at toroidal position phi
  type, extends(t_curve) :: t_surface_cut
     real(real64) :: phi
  end type t_surface_cut


  ! full flux surface
  type, public :: t_flux_surface_3D
     type(t_surface_cut), dimension(:), allocatable :: slice

     real(real64) :: PsiN
     integer      :: n_sym, n_phi

     ! pre-sampled distance to flux surface
     type (t_distance) :: distance

     contains
     procedure :: new, generate, plot
     procedure :: generate_from_axisymmetric_surface
     procedure :: get_distance_to
     procedure :: load_distance_to
     procedure :: sample_distance_to
     !procedure :: expand
  end type t_flux_surface_3D
  !end type t_poincare_set

  contains
!=======================================================================



!=======================================================================
  subroutine new(this, n_sym, n_phi)
  use math
  class(t_flux_surface_3D) :: this
  integer, intent(in)      :: n_sym, n_phi

  integer :: i


  this%n_sym = n_sym
  this%n_phi = n_phi
  allocate (this%slice(0:n_phi-1))
  do i=0,n_phi-1
     this%slice(i)%phi = pi2 / n_sym / n_phi * i
  enddo

  end subroutine new
!=======================================================================



!=======================================================================
! y0          initial/reference point
! npoints     number of points per slice
! nsym        toroidal symmetry number
! nslice      number of slices within 0..360/nsym deg
! nsteps      number of step between slices
! solver      id of ODE solver
!=======================================================================
  subroutine generate(this, y0, npoints, nsym, nslice, nsteps, solver)
  use equilibrium
  class(t_flux_surface_3D) :: this
  real(real64), intent(in) :: y0(3)
  integer,      intent(in) :: npoints, nsym, nslice, nsteps, solver

  type(t_poincare_set)     :: P
  real(real64) :: Maxis(3)
  integer      :: i, n


  call this%new(nsym, nslice)
  n     = nsteps
  if (n == 0) n = 16
  call P%generate(y0, npoints, nsym, nslice, n, solver, .true.)


  ! set flux surface label (given by normalized poloidal flux)
  this%PsiN = get_PsiN(y0)


  ! check if each slice is complete
  do i=0,nslice-1
     if (P%slice(i)%npoints .ne. npoints) then
        write (6, *) 'error: incomplete flux surface!'
        write (6, *) 'number of elements: ', P%slice(i)%npoints, '/', npoints
        stop
     endif
     call this%slice(i)%new (npoints-1)
     this%slice(i)%phi = P%slice(i)%phi
     this%slice(i)%x   = P%slice(i)%x(:,1:2)
     Maxis             = get_magnetic_axis(this%slice(i)%phi)
     call this%slice(i)%sort_loop(Maxis(1:2))
     !call this%slice(i)%expand(Maxis(1:2), 2.d0)
  enddo

  end subroutine generate
!=======================================================================



!=======================================================================
  subroutine plot(this, iu, filename)
  class(t_flux_surface_3D)               :: this
  integer,          intent(in), optional :: iu
  character(len=*), intent(in), optional :: filename

  integer :: i, iu0


  ! set default unit number for output
  iu0 = 90

  ! Unit number given for output?
  if (present(iu)) iu0 = iu

  ! Output_File given?
  if (present(filename)) then
     open  (iu0, file=filename)
  endif

  ! write data
  do i=0,this%n_phi-1
     call this%slice(i)%plot(iu=iu0)
     write (iu0, *)
  enddo

  ! Output_File given?
  if (present(filename)) close (iu0)

  end subroutine plot
!=======================================================================



!=======================================================================
  subroutine generate_from_axisymmetric_surface(this, fs2D, n_sym, n_phi, n_theta)
  class(t_flux_surface_3D) :: this
  type(t_flux_surface_2D)  :: fs2D
  integer, intent(in)      :: n_sym, n_phi, n_theta

  real(real64) :: t, x(2)
  integer :: i, j


  call this%new(n_sym, n_phi)
  this%PsiN = fs2D%PsiN
  call fs2D%sort_loop()
  call fs2D%setup_angular_sampling()

  do i=0,n_phi-1
     call this%slice(i)%new(n_theta)
     do j=0,n_theta
        t = 1.d0 * j / n_theta
        call fs2D%sample_at(t, x)
        this%slice(i)%x(j,:) = x
     enddo
  enddo

  end subroutine generate_from_axisymmetric_surface
!=======================================================================



!=======================================================================
  function get_distance_to(this, p) result(d)
  class(t_flux_surface_3D) :: this
  real(real64), intent(in) :: p(3)
  real(real64)             :: d


  d = this%distance%evaluate(p)

  end function get_distance_to
!=======================================================================



!=======================================================================
  subroutine load_distance_to(this)
  class(t_flux_surface_3D) :: this
  end subroutine load_distance_to
!=======================================================================



!=======================================================================
  subroutine sample_distance_to(this, grid)
  use grid
  use math
  class(t_flux_surface_3D) :: this
  character(len=*), intent(in) :: grid

  integer, parameter :: iu = 32

  real(real64) :: &
     R_center = 1.d0, &        ! center of computational box [cm]
     width    = 1.d0, &        ! width and height of computational box [cm]
     height   = 1.d0
  integer      :: N_R, N_Z

  namelist /Grid_Layout/ &
     N_R, N_Z, R_center, width, height

  real(real64) :: x(2), R0, Z0
  integer      :: iflag, i, j, k


  write (6, *) 'Sample distance to flux surface'

  open  (iu, file=grid)
  read  (iu, Grid_Layout)
  close (iu)
  call this%distance%new(N_R, N_Z, R_center, width, height, this%n_phi, this%n_sym)


  ! calculate distances
  R0 = R_center - 0.5d0 * width
  Z0 =          - 0.5d0 * height
  do k=0,this%n_phi-1
     do i=1,N_R
        do j=1,N_Z
           x(1) = R0 + 1.d0*(i-1)/(N_R-1) * width
           x(2) = Z0 + 1.d0*(j-1)/(N_Z-1) * height

           this%distance%d(i,j,k+1) = this%slice(k)%get_distance_to(x)
        enddo
     enddo
  enddo



!  call read_grid (grid, log_progress=.false., use_coordinates=COORDINATES(CYLINDRICAL))
!  grid_loop: do
!     call get_next_grid_point (iflag, x)
!     if (iflag .ne. 0) exit grid_loop
!
!     write (99, *) this%slice(0)%get_distance_to(x(1:2))
!  enddo grid_loop
  end subroutine sample_distance_to
!=======================================================================



!=======================================================================
  subroutine new_distance(this, nr, nz, Rc, w, h, nphi, nsym)
  class(t_distance)        :: this
  integer, intent(in)      :: nr, nz, nphi, nsym
  real(real64), intent(in) :: Rc, w, h


  this%nr   = nr
  this%nz   = nz
  this%nphi = nphi+1
  this%nsym = nsym

  this%Rc   = Rc
  this%w    = w
  this%h    = h

  this%nord = 5
  if (allocated(this%d)) deallocate (this%d)
  allocate (this%d(nr, nz, nphi+1))

  end subroutine new_distance
!=======================================================================



!=======================================================================
  subroutine setup_distance(this)
  use bspline
  use math
  class(t_distance) :: this

  real(real64), dimension(:), allocatable :: Rtmp, Ztmp, Phitmp

  integer :: i, j, k


  allocate (Rtmp(this%nr), Ztmp(this%nz), Phitmp(this%nphi))
  allocate (this%Rnot  (this%nr  +this%nord), &
            this%Znot  (this%nz  +this%nord), &
            this%Phinot(this%nphi+this%nord))

  do i=1,this%nr
     Rtmp(i)   = this%Rc + this%w * (-0.5d0 + (i-1) / (this%nr-1))
  enddo
  call dbsnak (this%nr, Rtmp, this%nord, this%Rnot)

  do j=1,this%nz
     Ztmp(j)   =           this%h * (-0.5d0 + (j-1) / (this%nz-1))
  enddo
  call dbsnak (this%nz, Ztmp, this%nord, this%Znot)

  do k=1,this%nphi
     Phitmp(k) = 2.d0*pi / this%nsym * (k-1) / (this%nphi-1)
  enddo
  call dbsnak (this%nphi, Phitmp, this%nord, this%Phinot)


  allocate (this%dcoeff(this%nr, this%nz, this%nphi))
  call dbs3in (this%nr, Rtmp, this%nz, Ztmp, this%nphi, Phitmp, &
               this%d, this%nr, this%nz, this%nord, this%nord, this%nord, &
               this%Rnot, this%Znot, this%Phinot, this%dcoeff)
  deallocate (Rtmp, Ztmp, Phitmp)

  end subroutine setup_distance
!=======================================================================



!=======================================================================
  function evaluate_distance(this, x) result(d)
  use bspline
  use math
  class(t_distance)        :: this
  real(real64), intent(in) :: x(3)
  real(real64)             :: d

  real(real64) :: phi0


  phi0 = mod(x(3),pi2/this%nsym)
  d    = dbs3dr(0,0,0,x(1),x(2),phi0,this%nord,this%nord,this%nord, &
                this%Rnot,this%Znot,this%Phinot,this%nr,this%nz,this%nphi,this%dcoeff)

  end function evaluate_distance
!=======================================================================

end module flux_surface_3D
