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

     contains
     procedure :: setup
     procedure :: setup_updown_symmetric
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
  subroutine setup(this, phi, nsample, D, S)
  use search
  use mesh_spacing
  class(t_surface_cut)         :: this
  real(real64),     intent(in) :: phi
  integer,          intent(in) :: nsample
  class(t_dataset), intent(in) :: D
  type(t_spacing),  intent(in), optional :: S

  type(t_curve) :: C
  real(real64)  :: t, x(2)
  integer       :: j


  ! initialize flux surface cut
  call this%new(nsample-1)
  this%phi = phi

  ! set up curve for re-sampling
  call C%new(D%nrow)
  C%x = D%x(:,1:2);  C%x(D%nrow,1:2) = D%x(1,1:2)
  call C%setup_length_sampling()

  ! re-sample points on flux surface
  do j=0,nsample-1
     if (present(S)) then
        t = S%node(j, nsample-1)
     else
        t = Equidistant%node(j, nsample-1)
     endif
     call C%sample_at(t, x)
     this%x(j,:) = x
  enddo

  end subroutine setup
!=======================================================================



!=======================================================================
  subroutine setup_updown_symmetric(this, phi, nsample, D, S)
  use search
  use mesh_spacing
  class(t_surface_cut)         :: this
  real(real64),     intent(in) :: phi
  integer,          intent(in) :: nsample
  class(t_dataset), intent(in) :: D
  type(t_spacing),  intent(in), optional :: S

  type(t_curve) :: C
  real(real64)  :: t, x(2)
  integer       :: i180, ierr, j


  ! initialize flux surface cut
  call this%new(nsample-1)
  this%phi = phi

  ! find index for theta = 180 deg
  i180 = binary_interval_search(1, D%nrow, D%x(:,3), 180.d0, ierr)
  if (ierr .ne. 0) then
     call D%plot(filename='D_sort.err')
     write (6, *) 'error: cannot find theta = 180 deg on flux surface!'
     stop
  endif

  ! set up curve for re-sampling
  call C%new(i180+1)
  C%x(1:i180,1:2) = D%x(1:i180,1:2)

  ! set 1st node
  D%x(D%nrow,3) = D%x(D%nrow,3) - 360.d0
  t = D%x(1,2) / (D%x(1,2) - D%x(D%nrow,2))
  C%x(0,1:2) = D%x(1,1:2) + t * (D%x(D%nrow,1:2) - D%x(1,1:2))

  ! set last node
  t =  - D%x(i180,2) / (D%x(i180+1,2) - D%x(i180,2))
  C%x(i180+1,1:2) = D%x(i180,1:2) + t * (D%x(i180+1,1:2) - D%x(i180,1:2))

  call C%setup_length_sampling()
  ! re-sample points on flux surface
  do j=0,nsample-1
     if (j < nsample / 2) then
        if (present(S)) then
           t = S%node(j, nsample-1) * 2.d0
        else
           t = Equidistant%node(j, nsample-1) * 2.d0
        endif
        call C%sample_at(t, x)
        this%x(j,:) = x
     else
        if (present(S)) then
           t = (1.d0 - S%node(j, nsample-1)) * 2.d0
        else
           t = (1.d0 - Equidistant%node(j, nsample-1)) * 2.d0
        endif
        call C%sample_at(t, x)
        this%x(j,1) = x(1)
        this%x(j,2) = -x(2)
     endif
  enddo

  end subroutine setup_updown_symmetric
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
! required input:
! y0          initial/reference point
! npoints     number of points per slice
! nsym        toroidal symmetry number
! nslice      number of slices within 0..360/nsym deg
! nsteps      number of trace steps between slices
! solver      id of ODE solver
!
! optional input:
! poloidal_coordinate
! resample              resolution for re-sampled flux surface
! updown_symmetry       select slice which requires up/down symmetric representation
! stop_at_boundary	(.false. allows to generate a full flux surface beyond the boundary)
!=======================================================================
  subroutine generate(this, y0, npoints, nsym, nslice, nsteps, solver, &
      poloidal_coordinate, resample, spacings, updown_symmetry, stop_at_boundary)
  use equilibrium, only: get_PsiN
  use mesh_spacing
  class(t_flux_surface_3D) :: this
  real(real64), intent(in) :: y0(3)
  integer,      intent(in) :: npoints, nsym, nslice, nsteps, solver
  integer,      intent(in), optional :: poloidal_coordinate
  integer,      intent(in), optional :: resample, updown_symmetry
  logical,      intent(in), optional :: stop_at_boundary
  type(t_spacing), intent(in), optional :: spacings

  type(t_poincare_set)     :: P
  type(t_curve)            :: C
  real(real64) :: phi, t, x(2)
  integer      :: i, j, n, ipc, nsample, nupdown


  ! set up resolution for re-sampling of flux surface
  nsample = npoints
  if (present(resample)) nsample = resample


  ! select slice which requires up/down symmetric representation
  nupdown = -1
  if (present(updown_symmetry)) nupdown = updown_symmetry


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
  if (present(stop_at_boundary)) then
     call P%generate(y0, npoints, nsym, nslice, n, solver, stop_at_boundary)
  else
     call P%generate(y0, npoints, nsym, nslice, n, solver, .true.)
  endif


  ! set flux surface label (given by normalized poloidal flux)
  this%PsiN = get_PsiN(y0)


  ! check if each slice is complete, sort points and sample equidistant points on surface
  this%n_theta = nsample
  do i=0,nslice-1
     if (P%slice(i)%npoints .ne. npoints) then
        write (6, *) 'error: incomplete flux surface!'
        write (6, *) 'number of elements: ', P%slice(i)%npoints, '/', npoints
        stop
     endif
     phi   = P%slice(i)%phi

     ! set up flux surface slice from Poincare plot P%slice(i)
     if (nupdown .ne. i) then
        call this%slice(i)%setup(phi, nsample, P%slice(i), spacings)

     ! generate up/down symmetric representation for this slice
     else
        call this%slice(i)%setup_updown_symmetric(phi, nsample, P%slice(i), spacings)
     endif
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
           Lc = Lc + F%trace_1step()

           ! stop field line tracing at limit nlimit*Limit
           if (abs(Lc) > m*Limit) m = m+1
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
