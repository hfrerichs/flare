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
