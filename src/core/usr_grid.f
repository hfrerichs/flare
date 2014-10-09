!-----------------------------------------------------------------------
!     module grid
!-------------------------------------------------------------------------------
!>		This module provides an interface between computations on grid
!>		nodes and the grid topology, i.e. so that no information about
!>		the grid layout and coordinate system is required outside this
!>		module.
!
!>		Grid nodes can be defined in Cartesian, cylindrical and toroidal/poloidal
!>		coordinates. The grid layout can be irregular/unstructured (1D array)
!>		or regular/structured (2D array).
!> \author 	Heinke Frerichs (h.frerichs at fz-juelich.de)
!-------------------------------------------------------------------------------
      module usr_grid
      use math
      implicit none


!> Array with grid nodes
      real*8, dimension(:,:), pointer :: grid_data
!> Total number of grid nodes
      integer :: n_grid


      private
! grid id and resolutions
      integer :: grid_id, n_xyz, n_RZ, n_rt, n_R, n_Z, n_t, n_p, n_RZphi

! unit number for grid file
      integer, parameter :: giun = 10

! internal variables
      integer :: icount, i10
      logical :: log_progress_ = .true.


      public ::
     1    read_grid, COORDINATES,
     2    get_next_grid_point,
     3    new_grid, n_grid, grid_data,
     4    read_grid_usr

      contains
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> internal routine to read integer values from header
!-----------------------------------------------------------------------
      subroutine iscrape (iun, iout)
      integer, intent(in)  :: iun
      integer, intent(out) :: iout
      character*82 :: str
      read  (iun, 2000) str
      write (6,2001) str(3:82)
      read  (str(33:42),*) iout
 2000 format (a82)
 2001 format (8x,a80)
      end subroutine iscrape
!-----------------------------------------------------------------------
!> internal routine to read real*8 values from header
!-----------------------------------------------------------------------
      subroutine rscrape (iun, rout)
      integer, intent(in)  :: iun
      real*8, intent(out) :: rout
      character*82 :: str
      read  (iun, 2000) str
      write (6,2001) str(3:82)
      read  (str(33:42),*) rout
 2000 format (a82)
 2001 format (8x,a80)
      end subroutine rscrape
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Read grid points from data file "Grid_File" and store in "grid_data"
!-----------------------------------------------------------------------
      subroutine read_grid (Grid_File, log_progress, use_coordinates,
     .                      n_nodes)
      use parallel
      character*120, intent(in)     :: Grid_File
      logical, intent(in), optional :: log_progress
      character*12, intent(in), optional :: use_coordinates
      integer, intent(out), optional     :: n_nodes

      character*120 :: str
      integer :: i, j, my_coordinates
      real*8  :: R, Z, phi, cos_phi, sin_phi, rmin, theta,
     .           r_center, r_min, L_tmp
      real*8, dimension(:), allocatable :: R_tmp, Z_tmp, t_tmp, p_tmp
! unit number for grid file
      integer, parameter :: iun = 10
      real*8,  parameter :: pi  = 3.14159265358979323846264338328d0

      integer, parameter ::
     1    i_cartesian = 0,
     2    i_cylindrical = 1

      if (present(log_progress))
     .   log_progress_ = log_progress .and. (mype.eq.0)
      if (associated(grid_data)) deallocate(grid_data)
      

! select output coordinate system for grid_data
      my_coordinates = i_cartesian
      if (present(use_coordinates)) then
         select case (use_coordinates)
         case(COORDINATES(1))
            my_coordinates = i_cartesian
         case(COORDINATES(2))
            my_coordinates = i_cylindrical
         case default
            write (6,*) 'error in read_grid. parameter coordinates = ',
     .                  use_coordinates
            stop
         end select
      endif



      if (mype.eq.0) then
! open grid file
      write (6,1000) adjustl(trim(Grid_File))
      open  (iun, file=Grid_File, err=5000)

! read grid header
      read  (iun, 2000) str
      write (6,2001) str(3:74)
      if (str(3:9).ne.'grid_id') then
         write (6,*) 'no grid type defined'
         stop
      endif
      read  (str(13:14), '(i2)') grid_id

! read actual grid data
      select case (grid_id)
c-----------------------------------------------
c 0: irregular (1D) xyz-grid
      case(0)
          call iscrape(iun, n_xyz)
          n_grid = n_xyz
          allocate (grid_data(n_grid, 3))

          do icount=1,n_xyz
             read (iun, *) grid_data(icount,:)
          enddo

          if (my_coordinates.eq.i_cylindrical) then
             do icount=1,n_xyz
                R   = (grid_data(icount,1)**2 + grid_data(icount,2)**2)
                Z   =  grid_data(icount,3)
                phi = atan2(grid_data(icount,2), grid_data(icount,1))
                grid_data(icount,1) = R
                grid_data(icount,2) = Z
                grid_data(icount,3) = phi
             enddo
          endif
c-----------------------------------------------
c 1: irregular (1D) RZ-grid
      case(1)
          call iscrape(iun, n_RZ)
          n_grid = n_RZ
          allocate (grid_data(n_grid, 3))

          call rscrape(iun, phi)
          cos_phi = cos(phi/180.d0*pi)
          sin_phi = sin(phi/180.d0*pi)

          select case (my_coordinates)
          case (i_cartesian)
             do icount=1,n_RZ
                read (iun,*) R, Z
                grid_data(icount,1) = R * cos_phi
                grid_data(icount,2) = R * sin_phi
                grid_data(icount,3) = Z
             enddo
          case (i_cylindrical)
             do icount=1,n_RZ
                read (iun,*) R, Z
                grid_data(icount,1) = R
                grid_data(icount,2) = Z
                grid_data(icount,3) = phi/180.d0*pi
             enddo
          end select
c-----------------------------------------------
c 2: irregular (1D) rmin-theta-grid
      case(2)
          call iscrape(iun, n_rt)
          n_grid = n_rt
          allocate (grid_data(n_grid, 3))

          call rscrape(iun, r_center)
          call rscrape(iun, phi)
          cos_phi = cos(phi/180.d0*pi)
          sin_phi = sin(phi/180.d0*pi)

          select case (my_coordinates)
          case (i_cartesian)
             do icount=1,n_rt
                !read (iun,3001) rmin, theta
                read (iun,*) rmin, theta
                R = r_center + rmin * cos(theta/180.d0*pi)
                Z =            rmin * sin(theta/180.d0*pi)
                grid_data(icount,1) = R * cos_phi
                grid_data(icount,2) = R * sin_phi
                grid_data(icount,3) = Z
             enddo
          case (i_cylindrical)
             do icount=1,n_rt
                read (iun,*) rmin, theta
                R = r_center + rmin * cos(theta/180.d0*pi)
                Z =            rmin * sin(theta/180.d0*pi)
                grid_data(icount,1) = R
                grid_data(icount,2) = Z
                grid_data(icount,3) = phi/180.d0*pi
             enddo
          end select
c-----------------------------------------------
c 3: regular (2D) RZ-grid
      case (3)
          call iscrape(iun, n_R)
          call iscrape(iun, n_Z)
          n_grid = n_R * n_Z
          allocate (grid_data(n_grid, 3),
     .              R_tmp(n_R), Z_tmp(n_Z))

          call rscrape(iun, phi)
          cos_phi = cos(phi/180.d0*pi)
          sin_phi = sin(phi/180.d0*pi)

          do icount=1,n_R
             read (iun,3002) R_tmp(icount)
          enddo
          do icount=1,n_Z
             read (iun,3002) Z_tmp(icount)
          enddo

          select case (my_coordinates)
          case (i_cartesian)
             do i=0,n_R-1
             do j=0,n_Z-1
                grid_data(i*n_Z + j + 1,1) = R_tmp(i+1) * cos_phi
                grid_data(i*n_Z + j + 1,2) = R_tmp(i+1) * sin_phi
                grid_data(i*n_Z + j + 1,3) = Z_tmp(j+1)
             enddo
             enddo
          case (i_cylindrical)
             do i=0,n_R-1
             do j=0,n_Z-1
                grid_data(i*n_Z + j + 1,1) = R_tmp(i+1)
                grid_data(i*n_Z + j + 1,2) = Z_tmp(j+1)
                grid_data(i*n_Z + j + 1,3) = phi/180.d0*pi
             enddo
             enddo
          end select

          deallocate (R_tmp, Z_tmp)
c-----------------------------------------------
c 4: regular (2D) rmin-theta-grid
      case (4)
          call iscrape(iun, n_t)
          call iscrape(iun, n_r)
          n_grid = n_t * n_r
          allocate (grid_data(n_grid, 3),
     .              t_tmp(n_t), r_tmp(n_r))

          call rscrape(iun, r_center)
          call rscrape(iun, phi)
          cos_phi = cos(phi/180.d0*pi)
          sin_phi = sin(phi/180.d0*pi)

          do icount=1,n_r
             read (iun,3002) r_tmp(icount)
          enddo
          do icount=1,n_t
             read (iun,3002) t_tmp(icount)
          enddo

          select case (my_coordinates)
          case (i_cartesian)
             do i=0,n_t-1
             do j=0,n_r-1
                R = r_center + r_tmp(j+1) * cos(t_tmp(i+1)/180.d0*pi)
                Z =            r_tmp(j+1) * sin(t_tmp(i+1)/180.d0*pi)
                grid_data(i*n_r + j + 1,1) = R * cos_phi
                grid_data(i*n_r + j + 1,2) = R * sin_phi
                grid_data(i*n_r + j + 1,3) = Z
             enddo
             enddo
          case (i_cylindrical)
             do i=0,n_t-1
             do j=0,n_r-1
                R = r_center + r_tmp(j+1) * cos(t_tmp(i+1)/180.d0*pi)
                Z =            r_tmp(j+1) * sin(t_tmp(i+1)/180.d0*pi)
                grid_data(i*n_r + j + 1,1) = R
                grid_data(i*n_r + j + 1,2) = Z
                grid_data(i*n_r + j + 1,3) = phi/180.d0*pi
             enddo
             enddo
          end select

          deallocate (t_tmp, r_tmp)
c-----------------------------------------------
c 5: regular (2D) phi-theta-grid
      case (5)
          call iscrape(iun, n_p)
          call iscrape(iun, n_t)
          n_grid = n_p * n_t
          allocate (grid_data(n_grid, 3),
     .              p_tmp(n_p), t_tmp(n_t))

          call rscrape(iun, r_center)
          call rscrape(iun, r_min)

          do icount=1,n_p
             read (iun,3002) p_tmp(icount)
          enddo
          do icount=1,n_t
             read (iun,3002) t_tmp(icount)
          enddo

          select case (my_coordinates)
          case (i_cartesian)
             do i=0,n_p-1
             do j=0,n_t-1
                cos_phi = cos(p_tmp(i+1)/180.d0*pi)
                sin_phi = sin(p_tmp(i+1)/180.d0*pi)

                R = r_center + r_min * cos(t_tmp(j+1)/180.d0*pi)
                Z =          - r_min * sin(t_tmp(j+1)/180.d0*pi)
                grid_data(i*n_t + j + 1,1) = R * cos_phi
                grid_data(i*n_t + j + 1,2) = R * sin_phi
                grid_data(i*n_t + j + 1,3) = Z
             enddo
             enddo
          case (i_cylindrical)
             do i=0,n_p-1
             do j=0,n_t-1
                R = r_center + r_min * cos(t_tmp(j+1)/180.d0*pi)
                Z =          - r_min * sin(t_tmp(j+1)/180.d0*pi)

                grid_data(i*n_t + j + 1,1) = R
                grid_data(i*n_t + j + 1,2) = Z
                grid_data(i*n_t + j + 1,3) = p_tmp(i+1)/180.d0*pi
             enddo
             enddo
          end select

          deallocate (p_tmp, t_tmp)
c-----------------------------------------------
c 6: regular (2D) toroidal RZ-grid
      case (6)
          call iscrape(iun, n_RZ)
          call iscrape(iun, n_p)
          n_grid = n_RZ * n_p

          allocate (grid_data(n_grid, 3),
     .              R_tmp(n_RZ), Z_tmp(n_RZ), p_tmp(n_p))

          do icount=1,n_RZ
             read (iun,3003) R_tmp(icount), Z_tmp(icount), L_tmp
          enddo
          do icount=1,n_p
             read (iun,3002) p_tmp(icount)
          enddo

          select case (my_coordinates)
          case (i_cartesian)
             do i=0,n_RZ-1
             do j=0,n_p-1
                cos_phi = cos(p_tmp(j+1)/180.d0*pi)
                sin_phi = sin(p_tmp(j+1)/180.d0*pi)

                grid_data(j*n_RZ + i + 1,1) = R_tmp(i+1) * cos_phi
                grid_data(j*n_RZ + i + 1,2) = R_tmp(i+1) * sin_phi
                grid_data(j*n_RZ + i + 1,3) = Z_tmp(i+1)
             enddo
             enddo
          case (i_cylindrical)
             do i=0,n_RZ-1
             do j=0,n_p-1
                grid_data(j*n_RZ + i + 1,1) = R_tmp(i+1)
                grid_data(j*n_RZ + i + 1,2) = Z_tmp(i+1)
                grid_data(j*n_RZ + i + 1,3) = p_tmp(j+1)/180.d0*pi
             enddo
             enddo
          end select

          deallocate (R_tmp, Z_tmp, p_tmp)
c-----------------------------------------------
c 9: irregular (1D) RZphi-grid
      case(9)
          call iscrape(iun, n_RZphi)
          n_grid = n_RZphi
          allocate (grid_data(n_grid, 3))

          do icount=1,n_RZphi
             read (iun, *) grid_data(icount,:)
          enddo

          if (my_coordinates.eq.i_cartesian) then
             do icount=1,n_RZphi
                R   = grid_data(icount,1)
                Z   = grid_data(icount,2)
                phi = grid_data(icount,3)
                grid_data(icount,1) = R * cos(phi)
                grid_data(icount,2) = R * sin(phi)
                grid_data(icount,3) = Z
             enddo
          endif
      end select
c-----------------------------------------------
    
      close (iun)
      endif	! mype = 0

      call wait_pe
      call broadcast_inte_s (n_grid)
      if (mype .gt. 0) then
          allocate (grid_data(n_grid, 3))
      endif
      call broadcast_real   (grid_data, n_grid*3)

      icount = mype + 1 - nprs
      i10    = 0
      if (mype .eq. 0) write (6, *)
      if (log_progress_) write (6,1001) i10


      ! optional return number of grid points / nodes
      if (present(n_nodes)) then
         n_nodes = n_grid
      endif

      return
 1000 format (3x,'- Using grid file: ',a)
 1001 format ('progress:             ',i4,' %')
 2000 format (a120)
 2001 format (8x,a72)
 2008 format (i8)
 3000 format (3e18.10)
 3001 format (2e18.10)
 3002 format (1e18.10)
 3003 format (2e18.10,2x,f8.3)
 5000 write (6,5001) Grid_File
 5001 format ('error reading grid file: ', a120)
      stop
      end subroutine read_grid
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Returns the Cartesian coordinates of the next grid point in "x" and 
!> set "iflag = 0". If no more grid points are left then "iflag = -1".
!
!> The current grid point number is stored in the internal variable "icount".
!-----------------------------------------------------------------------
      subroutine get_next_grid_point (iflag, x)
      use parallel
      integer, intent(out)  :: iflag
      real*8, intent(out) :: x(3)
      real*8 :: r10

      x = 0.d0

      icount = icount + nprs
      iflag  = 0
      if (icount.gt.n_grid) then
         iflag = -1
         if (log_progress_) write (6,1001)
         return
      endif

      x = grid_data(icount,:)

      r10 = 10.d0 * icount/n_grid
      if (r10.ge.i10+1) then
         i10 = i10 + 1
         if (log_progress_) write (6,1000) i10*10
      endif

      return
 1000 format ('                      ',i4,' %')
 1001 format ('finished')
      end subroutine get_next_grid_point
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Create a new grid (allocate memory, reset counter and zero grid nodes.
!
!> The optional parameter "log_progress" controls if progress during
!> computations on grid nodes (i.e. when calling "get_next_grid_point")
!> is logged to standard output.
!-----------------------------------------------------------------------
      function new_grid (n_grid_, log_progress)
      use parallel
      integer, intent(in) :: n_grid_
      logical, intent(in), optional :: log_progress
      real*8, dimension(:,:), pointer :: new_grid

      if (present(log_progress)) log_progress_ = log_progress

      if (associated(grid_data)) deallocate (grid_data)
      new_grid => null()
      n_grid   =  n_grid_

      if (n_grid.gt.0) then
         allocate (grid_data(n_grid,3))
         grid_data = 0.d0
         new_grid  => grid_data
      endif
      icount = mype + 1 - nprs

      return
      end function new_grid
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! iformat = 1:	X, Y, Z [cm]
!           2:  R, Z [cm], phi [deg]
!           3:  R, Z [cm] at phi_default
! oformat = CARTESIAN, CYLINDRICAL
!-----------------------------------------------------------------------
      subroutine read_grid_usr (filename, iformat, oformat, phi_default)
      use dataset
      use math
      character*120, intent(in)    :: filename
      integer, intent(in)          :: iformat, oformat
      real*8, intent(in), optional :: phi_default

      integer, parameter :: iu = 10
      real*8, dimension(:,:), pointer :: G
      type(t_dataset)      :: D
      real*8  :: phi, yi(3), yo(3)
      integer :: icol, i, n


      phi = 0.d0
      if (present(phi_default)) phi = phi_default

      ! read data
      icol = 3
      if (iformat == 3) icol = 2
      call D%load(filename, icol)

      ! setup grid
      n = D%nrow
      G => new_grid(n)
      do i=1,n
         yi(1) = D%x(i,1)
         yi(2) = D%x(i,2)
         select case (iformat)
         case(1,2)
            yi(3) = D%x(i-1,3) / 180.d0 * pi
         case(3)
            yi(3) = phi / 180.d0 * pi
         end select
         call coord_trans (yi, iformat, yo, oformat)
         G(i,:) = yo
      enddo

      end subroutine read_grid_usr
!-----------------------------------------------------------------------

      end module usr_grid
