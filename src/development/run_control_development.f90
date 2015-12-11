subroutine run_control_development (Run_Type)
  use parallel
  implicit none

  character(len=120), intent(in) :: Run_Type


  if (firstP) write (6, 1000)
 1000 format ('WARNING: RUNNING SUB-PROGRAM IN DEVELOPMENT STAGE!')


  select case (Run_Type)
!  case ('hyperbolic_fixed_point')
!     call hyperbolic_fixed_point
  case ('TEST_setup_domain')
     call TEST_setup_domain()
  case ('TEST_base_grid')
     call TEST_base_grid()
  case ('TEST_flux_surface_2D')
     call TEST_flux_surface_2D()
  case ('TEST_correct_PsiN')
     call TEST_correct_PsiN()
  case ('TEST_adaptive_step_size')
     call TEST_adaptive_step_size()
  case default
     write (6, *) 'run type "', trim(Run_Type), '" not defined!'
     stop
  end select

end subroutine run_control_development




subroutine TEST_setup_domain()
  use fieldline_grid
  use divertor
  use xpaths
  use equilibrium
  implicit none

  type(t_xpath) :: rpath_test
  integer, dimension(:), allocatable :: connectX
  real(real64) :: PsiN, x(2)
  integer :: ix, nx, i, n


  nx = 2
  allocate (connectX(nx))
  connectX(1) = 2
  connectX(2) = 1

  call setup_grid_configuration()
  call setup_geometry(nx, connectX)

  call rpath_test%generateX(1, DESCENT_CORE, LIMIT_PSIN, 0.01d0, 1)
  call rpath_test%plot(filename='rpath1_to_axis.plt')
!  call rpath_test%generateX(2, DESCENT_CORE, LIMIT_PSIN, 0.01d0)
!  call rpath_test%plot(filename='rpath2_to_axis.plt')

  deallocate (connectX)




  ! test sample_at_PsiN
  n = 10
  do i=0,n-1
     PsiN = 1.d0 - 1.d0*i/n
     call rpath_test%sample_at_PsiN(PsiN, x)
     write (99, *) x, PsiN, get_PsiN(x)
  enddo
end subroutine TEST_setup_domain




!===============================================================================
  subroutine TEST_base_grid
  use iso_fortran_env
  use base_grid
  implicit none


  nX = 2
  allocate (connectX(nX))
  ! lsn
  nX = 1
  connectX(1) = 1

  ! ddn
  !connectX(1) = -2
  !connectX(2) = -2

  ! ddn
  !connectX(1) = 2
  !connectX(2) = 1

  call setup_grid_configuration()
  call make_base_grids_auto()

  end subroutine TEST_base_grid
!===============================================================================



!===============================================================================
  subroutine TEST_flux_surface_2D()
  use iso_fortran_env
  use flux_surface_2D
  use equilibrium
  use separatrix
  implicit none

  type(t_flux_surface_2D) :: F
  type(t_separatrix)      :: S
  real(real64)            :: r(3), y(3)
  integer                 :: ierr


  y(1) = 0.d0
  y(2) = 1.1d0
  y(3) = 0.d0
  r    = get_cylindrical_coordinates(y, ierr)
  if (ierr > 0) then
     write (6, *) 'error: ierr > 0!'
     stop
  endif


  !call F%generate(r(1:2), direction=RIGHT_HANDED)
  !call F%plot(filename='F_RH.plt')
  !call F%generate(r(1:2), direction=LEFT_HANDED)
  !call F%plot(filename='F_LH.plt')

  !call F%generate_branch(r(1:2), FORWARD, ierr, stop_at_boundary=.false.)
  !call F%plot(filename='test_forward.plt')
  !call F%generate_branch(r(1:2), BACKWARD, ierr, stop_at_boundary=.false.)
  !call F%plot(filename='test_backward.plt')
  !call F%generate_branch(r(1:2), CCW, ierr, stop_at_boundary=.false.)
  !call F%plot(filename='test_ccw.plt')
  !call F%generate_branch(r(1:2), CW, ierr, stop_at_boundary=.false.)
  !call F%plot(filename='test_cw.plt')

  call S%generate_new(1, debug=.true.)
  call S%plot(filename_prefix='S', parts=.true.)

  end subroutine TEST_flux_surface_2D
!===============================================================================



!===============================================================================
  subroutine TEST_correct_PsiN()
  use iso_fortran_env
  use equilibrium, only: get_PsiN, correct_PsiN
  use grid
  use fgsl
  implicit none

  type(fgsl_rng)      :: r
  type(fgsl_rng_type) :: t
  real(fgsl_double)   :: dx, dy

  type(t_grid) :: G
  real(real64) :: x(3), xp(2), xc(2), PsiN, PsiNc, dPsiN, ds, ds_min, ds_max
  integer      :: i, ierr, iterations


  call G%load('grid.dat')


  ! set up rng
  t = fgsl_rng_env_setup()
  t = fgsl_rng_default
  r = fgsl_rng_alloc (t)


  dPsiN = 0.001d0
  ds_min = HUGE(1.d0)
  ds_max = 0.d0
  ds     = 0.1d0
  do i=1,G%nodes()
     write (6, *) i
     x    = G%node(i)
     PsiN = get_PsiN(x)

     ! perturb node
     call fgsl_ran_dir_2D(r, dx, dy)

     xp(1) = x(1) + ds*dx;     xp(2) = x(2) + ds*dy
     write (98, *) xp(1:2)

     xc    = correct_PsiN(xp(1:2), PsiN, ierr, iterations=iterations)
     if (ierr > 0) then
        xc    = correct_PsiN(xp(1:2), PsiN, ierr, iterations=iterations, debug=80)
        stop
     endif
     PsiNc = get_PsiN(xc)
     write (99, *) x(1:2), xc, abs(PsiN-PsiNc), PsiNc, iterations, ierr
!
!     if (ierr == 0) then
!        ds = sqrt(sum((x(1:2)-xc)**2))
!        if (ds > ds_max) ds_max = ds
!        if (ds < ds_min) ds_min = ds
!        write (99, *) xc(1:2), PsiN+dPsiN, iterations
!     else
!        write (98, *) x(1:2), ierr
!     endif
  enddo

!  write (6, *) 'min. correction step: ', ds_min
!  write (6, *) 'max. correction step: ', ds_max

  call fgsl_rng_free(r)
  end subroutine TEST_correct_PsiN
!===============================================================================




!===============================================================================
  subroutine TEST_adaptive_step_size()
  use iso_fortran_env
  use run_control, only: Grid_File, Output_File
  use grid
  use dataset
  use flux_surface_2D
  implicit none

  type(t_grid)    :: G
  type(t_dataset) :: D

  real(real64)    :: r(3)
  integer :: i, n


  call G%load(Grid_File)
  n = G%nodes()

  call D%new(n,1)
  do i=1,n
     r        = G%node(i)
     D%x(i,1) = adaptive_step_size(r(1:2))
  enddo
  call D%store(filename=Output_File)

  end subroutine TEST_adaptive_step_size
!===============================================================================
