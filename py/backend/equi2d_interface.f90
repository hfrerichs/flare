module equi2d_interface
  use types
  implicit none

  ! critical points of equilibrium
  !real(dp), dimension(:,:), allocatable :: xc
  !integer :: nc
  real(dp), dimension(:,:), allocatable :: xx, xo
  integer :: nx, no

  contains
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  ! Load equilibrium data file
  ! Error codes:
  !   1:  invalid format string given
  !   2:  no format string given and guess failed
  !   3:  loading data file failed
  !---------------------------------------------------------------------
  subroutine load(data_file, data_format, ierr)
  use equilibrium_format
  use equilibrium, only: load_equilibrium_data, setup_equilibrium, i_equi, setup_magnetic_axis
  use boundary
  character(len=*), intent(in)  :: data_file, data_format
  integer,          intent(out) :: ierr


  i_equi = get_equilibrium_format_from_string(data_format)
  ! error 1: invalid data format
  if (i_equi == EQ_UNDEFINED) then
     ierr = 1
     return
  endif

  ! guess equilibrium data format, if necessary
  if (i_equi == EQ_GUESS) then
     i_equi = get_equilibrium_format(data_file)

     ! error 2: guess failed
     if (i_equi == EQ_UNDEFINED) then
        ierr = 2
        return
     endif
  endif


  ! load equilibrium data
  ierr = 0
  call load_equilibrium_data(data_file, ierr)
  if (ierr > 0) return

  ! set up backend
  call setup_equilibrium()
  call setup_magnetic_axis()
  call setup_boundary()

  end subroutine load
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  subroutine get_psi(r, z, nr, nz, psi)
  use equilibrium, only: get_Psi_backend => get_Psi
  real(dp),         intent(in)  :: r(nr), z(nz)
  real(dp),         intent(out) :: psi(nz,nr)
  integer                       :: nr, nz
!f2py intent(in) :: r
!f2py intent(in) :: z
!f2py integer intent(hide),depend(r) :: nr=shape(r,0)
!f2py integer intent(hide),depend(z) :: nz=shape(z,0)

  real(dp) :: r3(3)
  integer  :: i, j


  r3(3) = 0.d0
  do i=1,nr
  do j=1,nz
     r3(1)    = r(i)
     r3(2)    = z(j)
     psi(j,i) = get_Psi_backend(r3)
  enddo
  enddo

  end subroutine get_psi
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  subroutine auto_setup_xpoints(nR, nZ)
  use equilibrium, only: find_hyperbolic_points
  integer, intent(in) :: nR, nZ


  call find_hyperbolic_points(nR, nZ, .true., .false.)

  end subroutine auto_setup_xpoints
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  ! create configuration file bfield.conf
  !---------------------------------------------------------------------
  subroutine write_config_file(data_file)
  use equilibrium_format, only: wcf => write_config_file
  use equilibrium, only: i_equi
  character(len=*), intent(in)  :: data_file


  call wcf(data_file, i_equi)

  end subroutine write_config_file
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  ! initialize magnetic configuration for new equilibrium
  !---------------------------------------------------------------------
  subroutine init(ierr)
  use equilibrium_format
  integer,          intent(out) :: ierr


  call initialize_equilibrium()
  ierr = 0

  end subroutine init
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  subroutine get_domain(rmin, rmax, zmin, zmax)
  use equilibrium, only: EQBox
  real(dp), intent(out) :: rmin, rmax, zmin, zmax

  rmin = EQBox(1,1);   rmax = EQBox(1,2)
  zmin = EQBox(2,1);   zmax = EQBox(2,2)
  end subroutine get_domain
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine get_nx(nx)
  use equilibrium, only: Xp, nX_max
  integer, intent(out) :: nx

  integer :: i


  do i=1,nX_max
     if (Xp(i)%undefined) return
     nx = i
  enddo

  end subroutine get_nx
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine cp_analysis(nr, nz, r, z, i)
  use equilibrium, only: EQBox
  use grid
  use dataset
  integer,  intent(in)  :: nr, nz
  real(dp), intent(out) :: r(nr), z(nz), i(nr*nz)

  type(t_grid)    :: G
  type(t_dataset) :: D
  integer :: j, nc


  call G%create_rlinspace(nr, EQBox(1,1), EQBox(1,2), nz, EQBox(2,1), EQBox(2,2))
  r = G%x1
  z = G%x2
  call critical_point_analysis(G, D, nc)
  i = D%x(:,1)

  ! save list of critical points
  !if (allocated(xc)) deallocate(xc)
  !allocate (xc(nc,2))
  !xc = D%x(1:nc, 2:3)
  nx = 0;   no = 0
  do j=1,nc
     if (D%x(j,4) == 0) then
         nx = nx + 1
     else
         no = no + 1
     endif
  enddo

  if (allocated(xx)) deallocate(xx, xo)
  allocate (xx(nx,2), xo(no,2))
  nx = 0;   no = 0
  do j=1,nc
     if (D%x(j,4) == 0) then
         nx = nx + 1
         xx(nx,:) = D%x(j,2:3)
     else
         no = no + 1
         xo(no,:) = D%x(j,2:3)
     endif
  enddo

  end subroutine cp_analysis
  !---------------------------------------------------------------------

end module equi2d_interface
