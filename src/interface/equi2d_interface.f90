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
  subroutine get_equilibrium_format_from_string(sformat, i_equi)
  use equilibrium_format
  character(len=*), intent(in)  :: sformat
  integer,          intent(out) :: i_equi


  select case(sformat)
  case (S_GEQDSK, S_GEQDSK_FREE)
     i_equi   = EQ_GEQDSK
  case (S_DIVAMHD)
     i_equi   = EQ_DIVAMHD
  case (S_SONNET)
     i_equi   = EQ_SONNET
  case (S_M3DC1)
     i_equi   = EQ_M3DC1
  case (S_AMHD)
     i_equi   = EQ_AMHD
  case default
     i_equi   = EQ_UNDEFINED
  end select

  end subroutine get_equilibrium_format_from_string
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  subroutine guess_equilibrium_format(filename, i_equi)
  use equilibrium_format
  character(len=*), intent(in)  :: filename
  integer,          intent(out) :: i_equi


  i_equi = get_equilibrium_format(filename)

  end subroutine guess_equilibrium_format
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  subroutine load(filename, iformat, ierr)
  use equilibrium_format
  use equilibrium, only: load_equilibrium_data, setup_equilibrium, i_equi
  character(len=*), intent(in)  :: filename
  integer,          intent(in)  :: iformat
  integer,          intent(out) :: ierr


  ierr   = 0
  i_equi = iformat
  call load_equilibrium_data(filename, ierr)
  call setup_equilibrium()

  end subroutine load
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