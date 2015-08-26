!===============================================================================
! Generate (non-axisymmetric) flux surfaces, based on module poincare_set
!===============================================================================
module flux_surface_3D
  use iso_fortran_env
  use dataset
  use curve2D
  use poincare_set
  use interpolate3D
  implicit none
  private


  ! flux surface contour at toroidal position phi
  type, extends(t_curve) :: t_surface_cut
     real(real64) :: phi
  end type t_surface_cut


  ! full flux surface
  type, public :: t_flux_surface_3D
     type(t_surface_cut), dimension(:), allocatable :: slice
     !type(t_slice), dimension(:), allocatable :: slice

     real(real64) :: PsiN
     integer      :: n_sym, n_phi, n_theta

     ! pre-sampled distance to flux surface
     type (t_interpolate3D) :: distance

     contains
     procedure :: new
     procedure :: load
     procedure :: generate
     procedure :: plot
     procedure :: destroy
     procedure :: generate_from_axisymmetric_surface
     procedure :: get_distance_to
     procedure :: load_distance_to
     procedure :: field_line_loss !(in the presence of magnetic perturbations)
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


  call this%destroy()
  this%n_sym = n_sym
  this%n_phi = n_phi
  allocate (this%slice(0:n_phi-1))
  do i=0,n_phi-1
     this%slice(i)%phi = pi2 / n_sym / n_phi * i
  enddo

  end subroutine new
!=======================================================================



!=======================================================================
  subroutine load(this, filename)
  class(t_flux_surface_3D)     :: this
  character(len=*), intent(in) :: filename

  integer, parameter :: iu = 32

  character(len=80)  :: s
  integer :: i, j, n_phi, n_sym, n_theta


  open  (iu, file=filename)
  ! read flux surface resolution
  read  (iu, 1000) s
  read  (s(3:80), *) n_phi, n_theta, n_sym
  this%n_theta = n_theta
  call this%new(n_sym, n_phi)

  ! read flux surface label (i.e. normlized poloidal flux)
  read  (iu, 1000) s
  read  (s(3:80), *) this%PsiN

  ! read slices
  do i=0,n_phi-1
     read (iu, 1000) s
     call this%slice(i)%new(n_theta-1)
     do j=0,n_theta-1
        read (iu, *) this%slice(i)%x(j,:)
     enddo
     read (iu, 1000) s
  enddo
  close (iu)

 1000 format(a80)
  end subroutine load
!=======================================================================



!=======================================================================
! y0          initial/reference point
! npoints     number of points per slice
! nsym        toroidal symmetry number
! nslice      number of slices within 0..360/nsym deg
! nsteps      number of trace steps between slices
! solver      id of ODE solver
!=======================================================================
  subroutine generate(this, y0, npoints, nsym, nslice, nsteps, solver, poloidal_coordinate)
  use magnetic_axis
  use equilibrium, only: get_PsiN
  use mesh_spacing
  class(t_flux_surface_3D) :: this
  real(real64), intent(in) :: y0(3)
  integer,      intent(in) :: npoints, nsym, nslice, nsteps, solver
  integer,      intent(in), optional :: poloidal_coordinate

  type(t_poincare_set)     :: P
  type(t_curve)            :: C
  real(real64) :: Maxis(3), phi, t, x(2)
  integer      :: i, j, n, ipc


  ! set poloidal coordinate
  ipc = ANGLE
  if (present(poloidal_coordinate)) ipc = poloidal_coordinate
  if (ipc < ANGLE  .or.  ipc > DISTANCE) then
     write (6, *) 'error in t_flux_surface_3D%generate: ', &
                  'invalid poloidal coordinate id ', ipc, '!'
     stop
  endif


  ! initialize and generate Poincare plot for y0
  call this%new(nsym, nslice)
  n     = nsteps
  if (n == 0) n = 16
  call P%generate(y0, npoints, nsym, nslice, n, solver, .true.)


  ! set flux surface label (given by normalized poloidal flux)
  this%PsiN = get_PsiN(y0)


  ! check if each slice is complete, sort points and sample equidistant points on surface
  this%n_theta = npoints
  do i=0,nslice-1
     if (P%slice(i)%npoints .ne. npoints) then
        write (6, *) 'error: incomplete flux surface!'
        write (6, *) 'number of elements: ', P%slice(i)%npoints, '/', npoints
        stop
     endif
     phi   = P%slice(i)%phi
     Maxis = get_magnetic_axis(phi)
     call C%new(npoints-1)
     C%x   = P%slice(i)%x(:,1:2)
     call C%sort_loop(Maxis(1:2), method=ipc)


     call this%slice(i)%new(npoints-1)
     this%slice(i)%phi = phi
     do j=0,npoints-1
        t = Equidistant%node(j, npoints-1)
        call C%sample_at(t, x)
        this%slice(i)%x(j,:) = x
     enddo
  enddo

  end subroutine generate
!=======================================================================



!=======================================================================
  subroutine plot(this, iu, filename, output_format)
  use math
  class(t_flux_surface_3D)               :: this
  integer,          intent(in), optional :: iu, output_format
  character(len=*), intent(in), optional :: filename

  integer :: i, j, iu0, oformat


  ! set default unit number for output, and default format
  iu0 = 90
  oformat = 1
  if (present(output_format)) oformat = output_format

  ! Unit number given for output?
  if (present(iu)) iu0 = iu

  ! Output_File given?
  if (present(filename)) then
     open  (iu0, file=filename)
  endif


  ! write flux surface resolution and label
  write (iu0, 1000) this%n_phi, this%n_theta, this%n_sym
  write (iu0, 1001) this%PsiN


  ! write data
  do i=0,this%n_phi-1
     select case (oformat)
     case(1)
        write (iu0, 1002) this%slice(i)%phi / pi * 180.d0
        call this%slice(i)%plot(iu=iu0)
        write (iu0, *)
     case(2)
        do j=0,this%slice(i)%n_seg-1
           write (iu0, *) this%slice(i)%x(j,:), this%slice(i)%phi / pi * 180.d0
        enddo
     case default
     end select
  enddo

  ! Output_File given?
  if (present(filename)) close (iu0)

 1000 format('# ',3i8)
 1001 format('# ',f10.5)
 1002 format('# phi = ',f10.5)
  end subroutine plot
!=======================================================================



!=======================================================================
  subroutine destroy(this)
  class(t_flux_surface_3D) :: this

  integer :: i, n1, n2


  if (allocated(this%slice)) then
     n1 = lbound(this%slice, 1)
     n2 = ubound(this%slice, 1)
     do i=n1,n2
        call this%slice(i)%destroy()
     enddo
     deallocate(this%slice)
  endif

  !call this%distance%destroy()

  end subroutine destroy
!=======================================================================



!=======================================================================
  subroutine generate_from_axisymmetric_surface(this, fs2D, n_sym, n_phi, n_theta)
  use flux_surface_2D
  class(t_flux_surface_3D) :: this
  type(t_flux_surface_2D)  :: fs2D
  integer, intent(in)      :: n_sym, n_phi, n_theta

  real(real64) :: t, x(2)
  integer :: i, j


  call this%new(n_sym, n_phi)
  this%PsiN    = fs2D%PsiN
  this%n_theta = n_theta
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


  d = this%distance%eval(p)

  end function get_distance_to
!=======================================================================



!=======================================================================
  subroutine load_distance_to(this)
  class(t_flux_surface_3D) :: this


  call this%distance%load(filename='distance.dat')

  end subroutine load_distance_to
!=======================================================================



!=======================================================================
  function field_line_loss(S, nlimit, Limit) result(iloss)
  use run_control, only: Trace_Step, Trace_Method, Trace_Coords
  use fieldline
  use parallel
  implicit none
  class(t_flux_surface_3D) :: S
  integer,      intent(in) :: nlimit
  real(real64), intent(in) :: Limit
  integer                  :: iloss(-1:1, nlimit)

  type(t_fieldline) :: F
  real(real64)      :: r0(3), y0(3), Lc, PsiNext
  integer :: i, j, idir, n, m


  iloss = 0
  n     = 0
  do i=0,S%n_phi-1
  do j=1,S%slice(i)%n_seg
     r0(1:2) = S%slice(i)%x(j,:)
     r0(3)   = S%slice(i)%phi
     n       = n + 1
     if (mod(n,nprs) .ne. mype) cycle
     call coord_trans (r0, CYLINDRICAL, y0, Trace_Coords)

     ! forward and backward tracing
     do idir=-1,1,2
        call F%init(y0, idir*Trace_Step, Trace_Method, Trace_Coords)
        Lc = 0.d0
        m  = 1
        trace_loop: do
           call F%trace_1step()
           Lc = Lc + Trace_Step

           ! stop field line tracing at limit nlimit*Limit
           if (Lc > m*Limit) m = m+1
           if (m > nlimit) exit trace_loop

           ! check intersection with boundary
           if (F%intersect_boundary()) then
              iloss(idir,m:nlimit) = iloss(idir,m:nlimit) + 1
              exit trace_loop
           endif
        enddo trace_loop
     enddo
  enddo
  enddo

  call wait_pe()
  call sum_inte_data (iloss, 3*nlimit)
  end function field_line_loss
!=======================================================================

end module flux_surface_3D
