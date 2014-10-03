!===============================================================================
! Anything related to the generation of magnetic field aligned grids
!===============================================================================
module field_aligned_grid
  use iso_fortran_env
  use equilibrium
  use separatrix
  use curve2D
  use mesh_spacing
  use fieldline
  use flux_surface_2D
  implicit none

  real(real64), parameter :: &
     epsilon_r64 = epsilon(real(1.0,real64))


  character(len=*), parameter :: &
     LSN    = 'LSN', &
     DDN    = 'DDN', &
     SIMPLE = 'SIMPLE'
  
  character(len=*), parameter :: &
     EXACT = 'EXACT', &
     QUASI = 'QUASI'

 
  integer, parameter :: &
     N_block_max = 360         ! Maximum number of toroidal blocks


  ! Type of innermost flux surface (exact or quasi)
  character(len=12) :: &
     Layout = SIMPLE, &
     Innermost_Flux_Surface = EXACT

  character(len=120) :: &
     Base_Grid_Alignment = '', &
     Radial_Discretization = 'auto', &
     Radial_Spacing        = '', &
     Poloidal_Spacing      = '', &
     Zones                 = 'automatic'


  integer :: &
     N_sym   = 1, &              ! Toroidal symmetry (i.e. 5 => 72 deg)
     N_block = 1, &              ! Number of toroidal blocks
     ir1     = 0, &              ! Radial index of inner(1) and outer(2) ...
     ir2     = 0, &              ! ... guiding surfaces
     nt      = 16, &             ! default toroidal resolution
     np      = 360, &            ! default poloidal resolution
     nr      = 32, &             ! default radial resolution
     ntx(N_block_max,-1:1) = -1  ! user defined toroidal resolution for individual blocks

  real(real64) :: &
     x_in1(3)                    = 0.d0, &  ! reference points (R[cm], Z[cm], phi[deg]) ...
     x_in2(3)                    = 0.d0, &  ! ... on 1st and 2nd innermost flux surfaces
     x_g1(3)                     = 0.d0, &  ! reference points (R[cm], Z[cm], phi[deg]) ...
     x_g2(3)                     = 0.d0, &  ! ... for inner(1) and outer(2) guiding surfaces
     x_out(3)                    = 0.d0, &  ! reference point on outer boundary
     Phi0                        = 0.d0, &  ! lower boundary of simulation domain
     Delta_Phi(N_block_max,-1:1) = 0.d0     ! user defined (non-default) toroidal block width [deg]


  ! Stellarator symmetry for 1st base grid
  logical :: &
     stellarator_symmetry = .false.


  ! Grid_Setup / Grid_Configuration
  namelist /Basic_Input/ &
     Layout, N_sym, N_block, x_in1, x_in2, Innermost_Flux_Surface, &
     nt, ntx, np, nr, ir1, ir2, x_g1, x_g2, x_out, &
     Phi0, Delta_Phi, Base_Grid_Alignment, Radial_Discretization, Radial_Spacing, Poloidal_Spacing

! internal variables
  logical :: default_decomposition


  ! toroidal discretization (within zones)
  type t_zone
     integer :: nt(-1:1), np, nr
     real(real64), dimension(:), allocatable :: phi
  end type t_zone
  type(t_zone), dimension(:), allocatable :: TD


  real(real64) :: Delta_Phi_Sim, &
     Phi_base(N_block_max)       = 0.d0  ! user defined (non-default) base grid locations [deg]

  ! alignment for base grid
  type(t_curve) :: cAlign

  ! radial discretization setup
  type(t_gradPsiN_path) :: Radial_Path

  contains
!=======================================================================



!=======================================================================
  subroutine load_usr_conf

  integer, parameter :: iu = 12

  integer :: i


  ! 1. read from input file
  open  (iu, file='grid.conf', err=9000)
  read  (iu, Basic_Input, err=9000)
  close (iu)


  ! 2. setup internal variables
  ! default decomposition (if 1st Delta_Phi = 0)
  default_decomposition = (Delta_Phi(1,-1) + Delta_phi(1,1) == 0.d0)
  Delta_Phi_Sim         = real(360, real64) / N_sym

  ! set block size and position of base grids
  write (6, *)
  if (default_decomposition) then
     Phi0 = -0.5d0 * Delta_Phi_Sim / N_block
     write (6, 1000) Phi0, Phi0+Delta_Phi_Sim

     Delta_Phi = Delta_Phi_Sim / N_block / 2.d0
     do i=1,N_block
        Phi_base(i) = Delta_Phi_Sim / N_block * (i-1)
     enddo
  else
     write (6, 1001) Phi0, Phi0+Delta_Phi_Sim

     Phi_base(1) = Phi0 + Delta_Phi(1,-1)
     do i=2,N_block
        Phi_base(i) = Phi_base(i-1) + Delta_Phi(i-1,1) + Delta_Phi(i,-1)
     enddo
  endif

  ! output to screen
  write (6, 1002)
  do i=1,N_block
     write (6, 1003) i, Phi_base(i), Phi_base(i)-Delta_Phi(i,-1), Phi_base(i)+Delta_Phi(i,1)
  enddo
 1000 format (3x,'- Default decomposition of simulation domain (',f7.3,' -> ',f7.3,' deg):')
 1001 format (3x,'- User defined decomposition of simulation domain (',f7.3,' -> ',f7.3,' deg):')
 1002 format (8x,'block #, base location [deg], domain [deg]')
 1003 format (8x,      i7,5x,f7.3,':',5x,f7.3,' -> ',f7.3)


  ! 2.2. toroidal resolution
  if (mod(nt,2) .ne. 0) write (6, *) 'warning: odd toroidal resolution!'
  do i=1,N_block
     if (ntx(i,-1) == -1) ntx(i,-1) = nt/2
     if (ntx(i, 1) == -1) ntx(i, 1) = nt/2
     ntx(i,0) = ntx(i,-1) + ntx(i,1)
  enddo
  call setup_toroidal_discretization



  ! 2.3. base grid alignment
  if (Base_Grid_Alignment .ne. '') then
     call cAlign%load(Base_Grid_Alignment)
  else
     call cAlign%new(0)
  endif

  ! 3. sanity check of user defined input
  call check_usr_conf()

  return
 9000 write (6, *) 'error while reading input file grid.conf!'
  stop
  end subroutine load_usr_conf
!=======================================================================



!=======================================================================
! Setup toroidal discretization
!=======================================================================
  subroutine setup_toroidal_discretization

  integer :: i, i1, i2, j, jdir


!  write (6, *) 'Toroidal discretization:'
!  do i=1,N_block
!     write (6, *) i, ntx(i,-1), ntx(i,1)
!  enddo


  allocate (TD(N_block))
  do i=1,N_block
     TD(i)%nt = ntx(i,:)
     i1       = -ntx(i,-1)
     i2       =  ntx(i, 1)
     allocate (TD(i)%phi(i1:i2))

     TD(i)%phi(0) = Phi_base(i)
     do jdir=-1,1,2
     do j=1,TD(i)%nt(jdir)
        TD(i)%phi(jdir*j) = TD(i)%phi(0) + jdir*j*Delta_Phi(i,jdir) / TD(i)%nt(jdir)
     enddo
     enddo
  enddo


!  write (6, *) 'discretization details:'
!  do i=1,N_block
!     write (6, *) 'block ', i
!     do j=-TD(i)%nt(-1),TD(i)%nt(1)
!        write (6, *) j, TD(i)%phi(j)
!     enddo
!  enddo

  end subroutine setup_toroidal_discretization
!=======================================================================



!=======================================================================
! sanity check for user defined input
!=======================================================================
  subroutine check_usr_conf

  real(real64) :: dPhi, f


  if (.not. default_decomposition) then
     ! Simulation domain must add up to 360 deg / N_sym
     dPhi = sum(Delta_Phi(1:N_block,-1)) + sum(Delta_Phi(1:N_block,1))
     f    = abs((dPhi - Delta_Phi_Sim) / Delta_Phi_sim)
     if (f > epsilon_r64) then
        write (6, *) "error: block sizes don't add up to size of simulation domain", &
                     ' (', Delta_Phi_Sim, ' deg)'
        stop
     endif
  endif

  end subroutine check_usr_conf
!=======================================================================



!=======================================================================
! Generate pair of innermost boundaries
!=======================================================================
  subroutine generate_innermost_boundaries

  select case (Innermost_Flux_Surface)
  case (EXACT)
     call exact_surface (x_in1, 'lcfs0')
     call exact_surface (x_in2, 'lcfs1')
  case (QUASI)
     call quasi_surface (x_in1, 'lcfs0')
     call quasi_surface (x_in2, 'lcfs1')
  case default
     write (6, *) 'error: flux surface type ', trim(Innermost_Flux_Surface), ' not defined!'
     stop
  end select

  contains
!-----------------------------------------------------------------------
  subroutine exact_surface (x, s5)
  use run_control, only: N_mult, N_sym_RC => N_sym, N_points, x_start, Phi_output, Output_File

  real(real64),     intent(in) :: x(3)
  character(len=5), intent(in) :: s5

  character(len=3) :: sblock
  integer :: iblock


  ! set default number of points for Poincare plot
  if (N_points == 0) N_points = 1000

  N_sym_RC = N_sym
  x_start  = x
  if (default_decomposition) then
     N_mult   = N_block
     Output_File = s5//'.txt'
     if (N_block == 1) Output_File = s5//'_0.txt'
     call poincare_plot
  else
     N_mult   = 1
     do iblock = 1,N_block
        write (sblock, '(i3)') iblock-1
        Output_File = s5//'_'//trim(adjustl(sblock))//'.txt'
        Phi_output  = Phi_base(iblock)
        call poincare_plot
     enddo
  endif


  end subroutine exact_surface
!-----------------------------------------------------------------------
  subroutine quasi_surface (x, s5)
  real(real64),     intent(in) :: x(3)
  character(len=5), intent(in) :: s5

  write (6, *) 'generation of quasi flux surfaces not yet implemented!'
  stop
  end subroutine quasi_surface
!-----------------------------------------------------------------------
  end subroutine generate_innermost_boundaries
!=======================================================================



!=======================================================================
  subroutine generate_layout


  write (6, *)
  write (6, 1000) Layout
 1000 format (3x,'- Layout for base grid: ',a)

  select case (Layout)
  case ('1')
     call base_grid_1
  case ('DDN')
     call base_grid_ddn
  case default
     write (6, *) 'error: base grid layout ', trim(Layout), ' not supported!'
     stop
  end select

  ! base grid post-processing
  contains
!-----------------------------------------------------------------------
  subroutine base_grid_1
  use iso_fortran_env
  use mesh_spacing
  use grid
  use string
  implicit none

  type(t_flux_surface_2D) :: B_out, B_in1, B_in2, F
  type(t_spacing)         :: radg, rad1, rad2, pol
  type(t_grid)            :: G
  character(len=4)        :: zstr
  real(real64) :: r(3), xi, eta, y1(2), y2(2), MagAxis(2)
  integer :: i, j, iz, i1, i2


  ! 0. initialize base grid and spacing function
  r       = get_magnetic_axis(x_in1(3)/180.d0*pi)
  MagAxis = r(1:2)
  call G%new(CYLINDRICAL, UNSTRUCTURED, 3, nr+1, np+1, mesh=.true.)
  call rad1%init(parse_string(Radial_Spacing, 1))
  call radg%init(parse_string(Radial_Spacing, 2))
  call rad2%init(parse_string(Radial_Spacing, 3))
  call pol%init(parse_string(Poloidal_Spacing, 1))


  ! 1. generate outer boundary
  write (6, 1000) x_out(1:2)
  call B_out%generate(x_out(1:2), AltSurf=Empty_curve, theta_cut=0.d0)
  call adjust_DED(B_out%t_curve)
  call B_out%sort_loop(MagAxis)
  call B_out%setup_angular_sampling(MagAxis)
  call B_out%plot(filename='outer_bounary.plt')


  ! 2.A generate mesh in inner domain ir1->ir2 from flux surfaces
  !!! GuidingSurfaces='AFS/NAFS'
  if (ir1 > 0) then
  if (ir2 >= nr) then
     write (6, *) 'error in subroutine base_grid_1: ir2 >= nr not allowed!'
     stop
  endif
  write (6, 1001) x_g1(1:2)
  write (6, 1002) x_g2(1:2)
  do i=ir1,ir2
     xi = radg%node(i-ir1, ir2-ir1)
     r  = x_g1 + xi*(x_g2-x_g1)

     ! axisymmetric flux surfaces
     call F%generate(r(1:2), AltSurf=Empty_curve, theta_cut=0.d0)
     call F%sort_loop(MagAxis)
     call F%setup_angular_sampling(MagAxis)
     if (i == ir1) call F%plot(filename='guiding_surface_1.plt')
     if (i == ir2) call F%plot(filename='guiding_surface_2.plt')

     do j=0,np
        eta = pol%node(j,np)
        call F%sample_at(eta, y2)
        G%mesh(i,j,:) = y2
     enddo
  enddo

  ! 2.1. generate mesh between outer boundary and 2nd guiding surface
  do j=0,np
     eta = pol%node(j,np)

     ! 2nd guiding surface
     y1 = G%mesh(ir2,j,:)

     ! outer boundary
     call B_out%sample_at(eta, y2)

     do i=ir2+1,nr
        xi = rad2%node(i-ir2, nr-ir2)
        G%mesh(i,j,:) = y1 + xi*(y2-y1)
     enddo
  enddo

  ! 2.2. generate mesh between inner boundaries and 1st guiding surface
  do iz=0,N_block-1
     write (zstr, '(i4)') iz
     ! load inner boundaries
     call B_in1%load('lcfs0_'//trim(adjustl(zstr))//'.txt')
     call B_in1%sort_loop(MagAxis)
     call B_in1%setup_angular_sampling(MagAxis)
     call B_in2%load('lcfs1_'//trim(adjustl(zstr))//'.txt')
     call B_in2%sort_loop(MagAxis)
     call B_in2%setup_angular_sampling(MagAxis)

     ! interpolate between inner boundary and 1st guiding surface
     do j=0,np
        eta = pol%node(j,np)

        ! 1st inner boundary
        call B_in1%sample_at(eta, y1)
        G%mesh(0,j,:) = y1

        ! 2nd inner boundary
        call B_in2%sample_at(eta, y1)

        ! 1st guiding surface
        y2 = G%mesh(ir1,j,:)

        do i=1,ir1-1
           xi = rad1%node(i-1, ir1-1)
           G%mesh(i,j,:) = y1 + xi*(y2-y1)
        enddo
     enddo

     ! write base grids
     G%fixed_coord_value = Phi_base(iz+1)
     call G%setup_mesh()
     call G%store('base_grid_'//trim(adjustl(zstr))//'.dat')
     call G%plot_mesh('base_grid_'//trim(adjustl(zstr))//'.plt')
  enddo

  ! 2.B interpolate between inner and outer boundaries
  else
  do iz=0,N_block-1
     write (zstr, '(i4)') iz
     ! load inner boundaries
     call B_in1%load('lcfs0_'//trim(adjustl(zstr))//'.txt')
     call B_in1%sort_loop(MagAxis)
     call B_in1%setup_angular_sampling(MagAxis)
     call B_in2%load('lcfs1_'//trim(adjustl(zstr))//'.txt')
     call B_in2%sort_loop(MagAxis)
     call B_in2%setup_angular_sampling(MagAxis)

     ! interpolate between inner boundary and 1st guiding surface
     do j=0,np
        eta = pol%node(j,np)

        ! 1st inner boundary
        call B_in1%sample_at(eta, y1)
        G%mesh(0,j,:) = y1

        ! 2nd inner boundary
        call B_in2%sample_at(eta, y1)

        ! 1st guiding surface
        call B_out%sample_at(eta, y2)

        do i=1,nr
           xi = rad1%node(i-1, nr-1)
           G%mesh(i,j,:) = y1 + xi*(y2-y1)
        enddo
     enddo

     ! write base grids
     G%fixed_coord_value = Phi_base(iz+1) / 180.d0 * pi
     call G%setup_mesh()
     call G%store('base_grid_'//trim(adjustl(zstr))//'.dat')
     call G%plot_mesh('base_grid_'//trim(adjustl(zstr))//'.plt')
  enddo
  endif




 1000 format(3x,'- Generating outer boundary at (R,Z) = (',f8.3,', ',f8.3,')')
 1001 format(3x,'- Generating 1st guiding surface at (R,Z) = (',f8.3,', ',f8.3,')')
 1002 format(3x,'- Generating 2nd guiding surface at (R,Z) = (',f8.3,', ',f8.3,')')
  end subroutine base_grid_1
!-----------------------------------------------------------------------
  subroutine base_grid_ddn
  use iso_fortran_env
  implicit none

  type t_X
     real(real64) :: X(2), theta, PsiN
     integer      :: orientation
  end type t_X


  real(real64), dimension(:,:), allocatable :: xsp
  real(real64), dimension(:,:,:), allocatable :: xsp_mesh
  type(t_separatrix) :: Si, So
  type(t_flux_surface_2D) :: F
  type(t_X)    :: Xp(2)
  real(real64) :: X(3), t
  integer      :: iXpi, iXpo, i

  
  ! upper X-point
  X(1:2) = find_uX()
  X(3)   = 0.d0
  Xp(1)%X           = X(1:2)
  Xp(1)%theta       = get_poloidal_angle(X)
  Xp(1)%orientation = -1
  Xp(1)%PsiN        = get_PsiN(X)

  ! lower X-point
  X(1:2) = find_lX()
  X(3)   = 0.d0
  Xp(2)%X           = X(1:2)
  Xp(2)%theta       = get_poloidal_angle(X)
  Xp(2)%orientation = 1
  Xp(2)%PsiN        = get_PsiN(X)

  ! set inner and outer X-point id
  if (Xp(1)%PsiN < Xp(2)%PsiN) then
     iXpi = 1; iXpo = 2
  else
     iXpi = 2; iXpo = 1
  endif


  ! generate inner separatrix
  write (6, 1000) Xp(iXpi)%X
  call Si%generate(Xp(iXpi)%X, Xp(iXpi)%orientation, Xp(iXpo)%theta, cAlign)
  call Si%plot('Si')


  ! generate outer separatrix
  write (6, 1001) Xp(iXpo)%X
  call So%generate(Xp(iXpo)%X, Xp(iXpo)%orientation, Xp(iXpi)%theta, cAlign)
  call So%plot('So')


  ! generate paths for discretization in radial direction
  call Radial_Path%generate(Xp(iXpi)%X, 1, PsiN=Xp(iXpo)%PsiN)
  call Radial_Path%plot(filename='radial_path_1.txt')

  call split_str(Radial_Spacing)
  call Radial_Path%setup_length_sampling()
  do i=1,9
     t = Equidistant%node(i,10)
     call Radial_Path%sample_at(t, X(1:2))
     call F%generate(X(1:2), Trace_Step=0.1d0)
     !call F%plot(98)
     !write (98, *)

     write (90, *) F%x(0,:)
     write (91, *) F%x(F%n_seg,:)
  enddo

  ! generate 1st SOL strike points
  allocate (xsp(1, 2))
  allocate (xsp_mesh(1, -TD(3)%nt(-1):TD(3)%nt(1), 2))
  xsp(1,:) = F%x(0,:)
  call generate_strike_point_mesh(1, TD(3), xsp, xsp_mesh)
  do i=-TD(3)%nt(-1),TD(3)%nt(1)
     write (80, *) xsp_mesh(1,i,:)
  enddo
  deallocate (xsp, xsp_mesh)


  ! generate outer boundary
  call Radial_Path%sample_at(1.d0, X(1:2))
  call Radial_Path%generate(X(1:2), 1, L=20.d0)
  call Radial_Path%setup_length_sampling()
  call Radial_Path%plot(filename='radial_path_2a.txt')
  call Radial_Path%sample_at(1.d0, X(1:2))
  call F%generate(X(1:2), Trace_Step=0.1d0, AltSurf=cAlign)
  call F%plot(filename='outer_boundary.plt')



  call Radial_Path%generate(Xp(iXpi)%X, 3, PsiN=0.8d0)
  call Radial_Path%plot(filename='radial_path_0.txt')
  call Radial_Path%generate(Xp(iXpi)%X, 4, L=20.d0)
  call Radial_Path%plot(filename='radial_path_3a.txt')


 1000 format(8x,'Inner X-point: (',f8.3,', ',f8.3,')')
 1001 format(8x,'Outer X-point: (',f8.3,', ',f8.3,')')
  end subroutine base_grid_ddn
!-----------------------------------------------------------------------
  end subroutine generate_layout
!=======================================================================



!=======================================================================
! Generate discretization of strike point for base grid
! Input:
!    n		radial resolution
!    Z		zone information (i.e. toroidal resolution and positions)
!    xsp(n,2)	initial (R,Z) coordinates of strike points (at Z%phi(0))
! Output:
!    xsp_mesh	mesh around strike points from field line tracing
!=======================================================================
  subroutine generate_strike_point_mesh (n, Z, xsp, xsp_mesh, Trace_Step, Trace_Method, Trace_Coords)
  integer, intent(in)       :: n
  type(t_zone), intent(in)  :: Z
  real(real64), intent(in)  :: xsp(n, 2)
  real(real64), intent(out) :: xsp_mesh(n, -Z%nt(-1):Z%nt(1), 2)
  real(real64), intent(in), optional :: Trace_Step
  integer, intent(in), optional :: Trace_Method, Trace_Coords

  type(t_fieldline) :: F
  real(real64) :: ts, y0(3), Dphi, y1(3)
  integer :: is, jdir, j, tm, tc, ierr


  ! set trace method
  tm = NM_AdamsBashforth4
  if (present(Trace_Method)) tm = Trace_Method


  ! set trace coordinates
  tc = FL_ANGLE
  if (present(Trace_Coords)) tc = Trace_Coords
  if (tc == CARTESIAN) then
     write (6, *) 'error: tracing in Cartesian coordinates not supported!'
  endif


  ! set trace step
  ts = pi2 / 3600.d0
  if (present(Trace_Step)) ts = abs(Trace_Step)


  ! update trace direction according to Bt direction
  ! -> i.e. ts steps are in positive toroidal direction
  if (tc .ne. FL_ANGLE) then
     ts = ts * Bt_sign
  endif


  y0(3) = Z%phi(0)/180.d0*pi
  do is=1,n
     y0(1:2) = xsp(is,:)
     xsp_mesh(is,0,:) = y0(1:2)

     ! jdir = -1: negative toroidal direction
     !         1: positive toroidal direction
     do jdir=-1,1,2
        call F%init(y0, jdir*ts, tm, tc)
        do j=1,Z%nt(jdir)
           Dphi = abs(Z%phi(jdir*j) - Z%phi(jdir*(j-1))) / 180.d0 * pi
           call F%trace_Dphi(Dphi, .false., y1, ierr)
           if (ierr .ne. 0) then
              write (6, *) 'error in subroutine generate_strike_point: ', &
                           'trace_Dphi returned error ', ierr
              stop
           endif
           write (6, *) jdir*Dphi/pi*180.d0, y1(3)/pi*180.d0, F%rc(3)/pi*180.d0, F%phi_int/pi*180.d0

           xsp_mesh(is,jdir*j,:) = y1(1:2)
        enddo
     enddo
  enddo

  end subroutine generate_strike_point_mesh
!=======================================================================



!=======================================================================
  subroutine split_str (str)
  character(len=*) :: str

  ! keyword=value
  end subroutine split_str
!=======================================================================


end module field_aligned_grid
