!===============================================================================
! Module for block limiter support
!
! Initial version (2012-07-30) by Felix Hasenbeck
! Extended (2012-08-28) for spherical block limiter by Heinke Frerichs
!===============================================================================
      MODULE block_limiter
!       USE const
      IMPLICIT NONE
      !SAVE

! cs1: coordinate system of the reactor
! cs2: coordinate system of the block limiter

      integer :: form_type,				! type of the block limiter
     &   num_sides = 5					! number of sides of the block limiter (default value: 5, used to make footprint grids)
      real*8 :: side_a, side_b, side_c			! side lengths of the cuboid
      real*8 :: length_d				! additionally needed length
      real*8 :: phi, theta_bl, r_min			! reference point of the cuboid in the torus
      real*8 :: alpha, beta				! orientation of the cuboid concerning the reference point

      real*8, dimension(7) :: area_side = 0.			! area of the block limiter sides
      real*8, parameter :: R_max = 175.				! main radius of the tokamak
      								! (centre point: centre point of the block limiter)
      integer :: ncount						! to count how many lines with data are in the data file
!       logical :: check_intersection				! function which calculates the intersection point. ...
      								! if there is an ip, check_intersection = true
      logical :: close_enough = .true.				! to check whether the line segment is in the nearer ...
      								! circumcircle of the block limiter 
      character*255 :: bl_file					! name of the file where the block data is stored
      character(len=120), allocatable, dimension(:), save :: data_lines	! lines which contain data from the file ...
      									! where the block limiter data is stored

! ADDITIONAL BLOCK LIMITER!!!  
      real*8, dimension(:,:,:), allocatable, save :: oripo_p, normv_p	! data to define the planes of the geometric form; ...
      								! THE NORMAL VECTORS HAVE TO BE ALREADY NORMALIZED!!!
      								! SHOULD BE DONE SO IN THE SUBROUTINES WHICH
      								! BUILD UP THE BLOCK LIMITER 
      integer, dimension(:), allocatable, save :: num_p					! number of planes of which the block limiter consists
      integer, dimension(:), allocatable, save :: num_s					! number of spheres of which the block limiter consists
      integer, parameter :: max_num_p = 10		! maximum number of planes per block limiter
      integer, parameter :: max_num_s = 10		! maximum number of spheres per block limiter
      real*8, dimension(:,:), allocatable, save :: dist_p_oripocs	! distance between planes and origin point of the cs
      real*8, dimension(:,:), allocatable, save :: center_bl				! center of the block limiter
      real*8, dimension(:), allocatable, save :: radius_bl					! radius of the sphere containing the block limiter
      real*8, dimension(:,:), allocatable, save :: radius_s		! radii of the sub-spheres
      real*8, dimension(:,:,:), allocatable, save :: center_s		! center position of the sub-spheres


      namelist /BlockLimiterData/
     &   form_type, side_a, side_b, side_c, length_d, phi, theta_bl,
     &   r_min, alpha, beta

! variables for bl_triangles_eirene
      integer, dimension(:,:), allocatable, save :: ori_imp_vtx	! origin planes of the important vertices
      real*8, dimension(:,:), allocatable, save :: imp_vtx		! important vertices (no double entries, lie on the surface of the geometric form)
      integer :: num_imp_vtx					! number of relevant vertices
      real*8, external :: sort_sin				! function to sort the values of sin (used by SR bl_triangles_eirene)

! variables for footprint_grids
      character*255 :: 
!     &   fp_machine,						! machine directory where the block limiter data is stored
!     &   fp_shot,						! shot directory where the block limiter data is stored
!     &   fp_bl_filename,						! name of the block limiter filename
     &   bl_filename,						! name of the block limiter filename
     &   fp_grid,						! name of the footprint grid file
     &   fp_grid_swop						! name of the swung-open footprint grid file
     
      integer, dimension(6) ::
     &   x_res, y_res						! arrays to store x- and y-resolution in
     
      real*8 :: alpha_w, alpha_r				! angle of wedge and roof

      integer :: first_side					! index of first side for making the footprint grids (for numeration consult sketch #2)	
     
      real*8, dimension(3), parameter ::
     &   e_x = (/ 1., 0., 0. /),				! unit vectors in x- and y-direction
     &   e_y = (/ 0., 1., 0. /),
     &   e_z = (/ 0., 0., 1. /),
     &   e_0 = (/ 0., 0., 0. /)

      namelist /FootprintGridData/
!     &   fp_machine, fp_shot, fp_bl_filename, fp_grid, fp_grid_swop,
     &   bl_filename, fp_grid, fp_grid_swop,
     &   x_res, y_res, first_side

      real*8, parameter ::
     &   pi   = 3.14159265358979323846264338328d0,
     &   pi2  = 2.d0 * pi
!------------------------------------------------------------------
!------------------------------------------------------------------
! DECLARATION OF  FUNCTIONS AND SUBROUTINES
!------------------------------------------------------------------
!------------------------------------------------------------------

      CONTAINS

!------------------------------------------------------------------
! SUBROUTINE TO MAKE GRIDS FOR GENERATING THE MAGNETIC FOOTPRINT
! OF THE BLOCK LIMITER SURFACE
!------------------------------------------------------------------
      SUBROUTINE make_footprint_grids(dist_edge, i_bl)
      IMPLICIT NONE

      real*8, intent(in) ::
     &   dist_edge				! security distance from the edge of block limiter sides
      
      integer ::
     &   iun = 70, iun2 = 71,
     &   i, j, k,				! counting variables
     &   n = 0,					! total number of grid points of the footprint grid
     &   num_total_po = 0,			! maximum number of grid points of the footprint grid (some might be cut away) 
     &   i_start, i_end				! indices to determine the cornerpoints of a block limiter plane
    
      integer, intent(in) ::
     &   i_bl					! can be used to make work with more than one bl

      real*8 ::
     &   w, r,					! side length of wedge and roof
     &   x_step, y_step,			! step length in x- and y-direction
     &   dist_plot = 5.				! distance from each side of another for the swopped grid

      real*8, parameter :: limit = 10.d-5	! value considered as zero

      real*8, dimension(3) ::
     &   disloc_v_z,				! dislocation of grid points from fp_grid_swop to fp_grid in z-direction
     &   oripo_cs,				! origin point of cs2 in cs1
     &   grid_po, grid_po_swop,			! grid point and grid point of swopped grid
     &   p_1, p_2, p_3, p_4			! corner points of a rectangular zone of the swopped footprint grid

      real*8, dimension(6) ::
     &   x_side, y_side				! side lengths for making the grid

      real*8, dimension(3,3) ::
     &   side_rot_m,				! rotation matrix for rotating one side
     &   cs_rot_m				! rotation matrix for transformation from cs2 to cs1

      real*8, dimension(6,3) ::
     &   start_po,
     &   disloc_v, disloc_v_swop,		! vectors to dislocate the grid points
     &   side_rot_v

      character*255 :: 
     &   fp_datafile,				! name of the file where the data for making the footprint grids is stored
     &   bl_datafile,				! name of the file where the block limiter data is stored
     &   fp_swop_zones,				! name of the file where the data of the rectangular zones of the swopped footprint grid should be stored to
     &   fp_grid_cornerpoints_planes,		! name of the file where the cornerpoints of the planes of the block limiter should be stored to
     &   homedir,
     &	 str					! temporarily used string

      real*8, dimension(:,:), allocatable ::
     &   grid_po_l, grid_po_swop_l		! all grid points and all grid points of swopped grid
     
      logical :: result

      write (6, 1000) i_bl
      p_1 = 0.
      p_2 = 0.
      p_3 = 0.
      p_4 = 0.
      
1000  format (3x, '- generating footprint grids for block limiter ', i4)
0815  format(t9, a/, t9, a/, t9, a)				! standard output format
3000  format('# grid_id = 0  (irregular (1D) xyz-grid) / ',  ! file header
     &       'x-resolution = (', 6i5, ') / y-resolution = (', 6i5, ')'/,
     &       '# number of grid points = ', i12)
3001  format(3e18.10)
3002  format(4e18.10)

! set names where to get the data from and where to write the data
!      bl_datafile = 'block_limiter.dat'
      fp_datafile = 'footprint_grid.conf'
      fp_swop_zones = 'fp_swop_zones.dat'
      fp_grid_cornerpoints_planes = "fp_grid_cornerpoints_planes.dat"

! read footprint grid data from data file
!      call read_footprint_grid_data(fp_datafile)

! read block limiter data and plot it
!      call getenv("HOME", homedir)
!      bl_datafile = bl_filename
!      call setup_block_limiter(bl_datafile, i_bl)

      IF (form_type == 4) then
         first_side = 5
!         write (6, *) 'first_side is set to 5 for spherical limiters!'
         call make_footprint_grid_sphereX (x_res(5),y_res(5),dist_edge)
         return
      ENDIF

! set num_sides (default value: 5)
      IF (form_type == 3) num_sides = 6

      IF (first_side > num_sides) THEN
         write(*, 0815) 'SUBROUTINE make_footprint_grids: ',
     &      'first_side > num_sides!',
     &      'Please modify >> footprint_grid.conf <<!'
      END IF

! check if resolution is greater than 0 for every side 
      DO i = first_side, num_sides
         IF (x_res(i) <= 0 .or. y_res(i) <= 0) THEN
            write(*, 0815) 'SUBROUTINE make_footprint_grids: ',
     &         'The grid resolution for one side is 0!',
     &          'Please modify >> footprint_grid.conf <<!'
            STOP
         END IF
      END DO

! calculate maximum number of grid points and allocate related arrays
      DO i = first_side, num_sides
         num_total_po = num_total_po+x_res(i)*y_res(i)
      END DO
      allocate(grid_po_l(num_total_po,3))
      allocate(grid_po_swop_l(num_total_po,3))
      grid_po_l = 0.
      grid_po_swop_l = 0.

! calculation of side length for wedge and roof
      w = sqrt(side_a**2 + length_d**2)			! side length wedge
      r = sqrt((side_a/2.)**2 + length_d**2)		! side length roof
!      print*, 'r = ', r

! define starting points
      start_po(1,:) = (/ -side_c/2., side_a/2., 0.d0 /) 
!       start_po(1,:) = (/ 1, 1, 1 /)
      start_po(2,:) = (/ -side_b/2., side_c/2., 0.d0 /) 
      start_po(3,:) = start_po(1,:)
      start_po(4,:) = start_po(2,:)
      start_po(5,:) = (/ -side_b/2., side_a/2., 0.d0 /)

! define the side lengths
      x_side(1) = side_c
      x_side(2) = side_b
      x_side(3) = x_side(1)
      x_side(4) = x_side(2)
      x_side(5) = x_side(2)
      y_side(1) = side_a
      y_side(2) = side_c
      y_side(3) = y_side(1)
      y_side(4) = y_side(2)
      y_side(5) = y_side(1)

! define side rotation vectors
      side_rot_v(1,:) = -e_y*pi/2.
      side_rot_v(2,:) = e_x*pi/2.
      side_rot_v(3,:) = -side_rot_v(1,:)
      side_rot_v(4,:) = -side_rot_v(2,:)
      side_rot_v(5,:) = e_0

! define dislocation vectors
      disloc_v_z = side_c/2.*(-e_z)
      disloc_v(1,:) = side_b/2.*e_x + disloc_v_z
      disloc_v(2,:) = side_a/2.*e_y + disloc_v_z
      disloc_v(3,:) = -side_b/2.*e_x + disloc_v_z
      disloc_v(4,:) = -side_a/2.*e_y + disloc_v_z
      disloc_v(5,:) = e_0

! define swopped dislocation vectors
      disloc_v_swop(1,:) = ((side_b+side_c)/2.+dist_plot)*e_x
      disloc_v_swop(2,:) = ((side_a+side_c)/2.+dist_plot)*e_y
      disloc_v_swop(3,:) = -disloc_v_swop(1,:)
      disloc_v_swop(4,:) = -disloc_v_swop(2,:)
      disloc_v_swop(5,:) = e_0

! define the area of the block limiter sides
      area_side(5) = x_side(5)*y_side(5)

! make modifications for wedge and roof block limiter
      SELECT CASE (form_type)
      CASE (2)
         start_po(5,:) = (/ -side_b/2., w/2., 0.d0 /)
         disloc_v(5,:) = length_d/2.*(-e_z)
         disloc_v_swop(2,:) = ((w+side_c)/2.+dist_plot)*e_y
         disloc_v_swop(4,:) = -disloc_v_swop(2,:)
         side_rot_v(5,:) = -alpha_w*e_x
         y_side(5) = w
         area_side(5) = w*x_side(5)
      CASE (3)
         start_po(5,:) = (/ -side_b/2., r/2., 0.d0 /)
         start_po(6,:) = start_po(5,:)
         disloc_v(5,:) = side_a/4.*(-e_y) + length_d/2.*(-e_z)
         disloc_v(6,:) = side_a/4.*e_y + length_d/2.*(-e_z)
         disloc_v_swop(2,:) = ((side_c/2.+r)+dist_plot)*e_y
         disloc_v_swop(4,:) = -disloc_v_swop(2,:)
         disloc_v_swop(5,:) = (r/2.)*(-e_y)		! distance of roof surfaces for the swopped grid can be defined
         disloc_v_swop(6,:) = -disloc_v_swop(5,:)
         side_rot_v(5,:) = -alpha_r*e_x
         side_rot_v(6,:) = alpha_r*e_x
         x_side(6) = x_side(5)
         y_side(5) = r
         y_side(6) = y_side(5)
         area_side(5) = r*x_side(5)
         area_side(6) = area_side(5)
      END SELECT

! get data for transformation from cs2 to cs1
      CALL get_cs_transformation_data(oripo_cs, cs_rot_m)

! write grid_po and grid_po_swop to an array
      DO i = first_side, num_sides						! choose starting point for i to define the sides for the footprint grid
         x_step = (x_side(i) - 2.*dist_edge)/(1.*x_res(i)-1.)			! define step size
         y_step = (y_side(i) - 2.*dist_edge)/(1.*y_res(i)-1.)
         CALL get_rotation_matrix(side_rot_v(i,:), side_rot_m)			! get side rotation matrix
         start_po(i,:) = start_po(i,:)+dist_edge*e_x-dist_edge*e_y
         DO j = 0, x_res(i)-1							! calculate grid points
            DO k = 0, y_res(i)-1
               grid_po_swop = start_po(i,:)+j*x_step*e_x-k*y_step*e_y		! get first version of grid point of the swopped grid ...
               grid_po = matmul(side_rot_m,grid_po_swop)+disloc_v(i,:)		! ... rotate und dislocate it to get grid point in cs2 ...
               grid_po = matmul(cs_rot_m, grid_po) + oripo_cs			! ... calculate grid point position in cs1
               CALL check_point_behind_bl(grid_po, limit, i_bl, result)		! check if grid point is on the block limiter ...
               IF (result) THEN					! ... if yes, write grid point to a file
                  n = n + 1
                  grid_po_l(n,:) = grid_po+dist_edge*normv_p(i_bl,i,:)
                  grid_po_swop_l(n,:) = grid_po_swop+disloc_v_swop(i,:)
               END IF
            END DO
         END DO
      END DO

! write grid points and swopped grid points to a file
      open(iun, file=fp_grid)
      write(iun, fmt=3000) x_res, y_res, n
      DO i = 1, n
         write(iun, fmt=3001) grid_po_l(i,:)
      END DO
!      print*, 'Footprint grid written!'
      close(iun)

      open(iun, file=fp_grid_swop)
      write(iun, fmt=3000) x_res, y_res, n
      DO i = 1, n
         write(iun, fmt=3001) grid_po_swop_l(i,:)
      END DO
!      print*, 'Swopped footprint grid written!'
      close(iun)

! get cornerpoints of the block limiter sides
!      open(iun, file=fp_grid_cornerpoints_planes)
!      k = first_side - 1
!      DO j = first_side, num_sides
!         i_start = 1
!         DO i = first_side, k
!            i_start = i_start + x_res(i)*y_res(i)
!         END DO
!         i_end = i_start - 1 + x_res(i)*y_res(i)
!         write(iun, fmt=3001) grid_po_l(i_start+y_res(j)-1, :)		! P_1
!         write(iun, fmt=3001) grid_po_l(i_start, :)			! P_2
!         write(iun, fmt=3001) grid_po_l(i_end-y_res(j)+1, :)		! P_3
!         write(iun, fmt=3001) grid_po_l(i_end, :)			! P_4
!         k = k + 1
!      END DO
!      close(iun)

! write rectangular zones of the swopped grid to a file (needed for emc3 post processing,
! flux footprint of the block limiter)
!      open(iun, file=fp_swop_zones)
!!      open(iun2, file='fp_zones_plot.dat')				! test
!      print*, "fp_zones_plot.dat written!"
!      i_start = 0
!      DO i = first_side, num_sides
!      print*, 'x_res = ', x_res(i)
!      print*, 'y_res = ', y_res(i)
!         DO j = 0, x_res(i)-2
!            DO k = 1, y_res(i)-1
!               p_1 = grid_po_swop_l(i_start + j*y_res(i) + k, :)
!               p_2 = grid_po_swop_l(i_start + j*y_res(i) + k + 1, :)
!               p_3 = grid_po_swop_l(i_start + (j+1)*y_res(i) + k, :)
!               p_4 = grid_po_swop_l(i_start + (j+1)*y_res(i) + k + 1, :)
!               write(iun, fmt=3002) p_1(1), p_2(1), p_3(1), p_4(1)
!               write(iun, fmt=3002) p_1(2), p_2(2), p_3(2), p_4(2)
!!               write(iun2, fmt=3001) p_1				! test
!!               write(iun2, fmt=3001) p_2				! test
!!               write(iun2, fmt=3001) p_4				! test
!!               write(iun2, fmt=3001) p_3				! test
!!               write(iun2, fmt=3001) p_1				! test
!            END DO
!         END DO
!         i_start = i_start + (x_res(i)*y_res(i))
!      END DO   
!      close(iun)
!!      close(iun2)							!test

      RETURN
      END SUBROUTINE make_footprint_grids
!------------------------------------------------------------------

!------------------------------------------------------------------
!------------------------------------------------------------------
      subroutine make_footprint_grid_sphereX (nx, ny, dr)
      implicit none

      !integer, parameter :: nx = 255, ny = 255
      integer, intent(in) :: nx, ny
      real*8, intent(in)  :: dr

      real*8, dimension(3) ::
     &   oripo_cs				! origin point of cs2 in cs1

      real*8, dimension(3,3) ::
     &   cs_rot_m				! rotation matrix for transformation from cs2 to cs1


      integer :: i, j
      real*8  :: x, y, z, f, r, x0(3)

      CALL get_cs_transformation_data(oripo_cs, cs_rot_m)

      f = (side_a**2 + side_b**2) / length_d / 8.d0 - length_d / 2.d0
      r = length_d + f
      !r = r * 1.00000001d0

      !open  (99, file='grid.data')
      !open  (98, file='grid_3D.data')
      open  (99, file=fp_grid_swop)
      open  (98, file=fp_grid)
      write (99, 3000) nx, ny, nx*ny
      write (98, 3000) nx, ny, nx*ny
      do j=0,ny-1
      y = (-0.5d0 + 1.d0 * j / (ny-1)) * side_a
      do i=0,nx-1
         x = (-0.5d0 + 1.d0 * i / (nx-1)) * side_b
         z = sqrt(r**2 - x**2 - y**2) - r
         write (99, 3001) x, y, z

         x0(1) = x
         x0(2) = y
         x0(3) = z + dr
         x0    = matmul(cs_rot_m, x0) + oripo_cs
         write (98, 3001) x0
      enddo
      enddo
      close (99)
      close (98)

      return
3000  format('# grid_id = 0  (irregular (1D) xyz-grid) / ',  ! file header
     &       'x-resolution = (', i5, ') / y-resolution = (', i5, ')'/,
     &       '# number of grid points = ', i12)
3001  format(3e18.10)
      end subroutine make_footprint_grid_sphereX
!------------------------------------------------------------------

!------------------------------------------------------------------
! SUBROUTINE TO READ THE BLOCK LIMITER DATA
!------------------------------------------------------------------
      SUBROUTINE read_footprint_grid_data(fp_datafile)	
      IMPLICIT NONE

      character*255, intent(in) :: fp_datafile
      integer :: iun = 33, iun2 = 32, ios

0815  format(t9, a, a, a, a, a, a, a, a, a, a, a)	! standard output format

! get the data from the block limiter file
      open(iun, file=fp_datafile, iostat=ios)
      IF (ios /= 0) THEN
         write(*, 0815) 'SUBROUTINE read_footprint_grid_data: ',
     &     'Could not read the footprint grid data file!'
         STOP
      END IF
      read(iun, FootprintGridData)
      close(iun)

      END SUBROUTINE read_footprint_grid_data
!------------------------------------------------------------------

!------------------------------------------------------------------
!------------------------------------------------------------------
      subroutine make_triangles_sphereX (nx, ny)
      implicit none

      !integer, parameter :: nx = 255, ny = 255
      integer, intent(in) :: nx, ny

      real*8, dimension(3) ::
     &   oripo_cs				! origin point of cs2 in cs1

      real*8, dimension(3,3) ::
     &   cs_rot_m				! rotation matrix for transformation from cs2 to cs1


      integer :: i, j
      real*8  :: x, y, z, f, r, x0(3)
      real*8  :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4

      CALL get_cs_transformation_data(oripo_cs, cs_rot_m)

      f = (side_a**2 + side_b**2) / length_d / 8.d0 - length_d / 2.d0
      r = length_d + f

      open  (99, file='triangles.dat')
!      open  (99, file=fp_grid_swop)
!      open  (98, file=fp_grid)
!      write (99, 3000) nx, ny, nx*ny
!      write (98, 3000) nx, ny, nx*ny
      do j=0,ny-1
         y1 = (-0.5d0 + 1.d0 *  j    / ny) * side_a
         y2 = y1
         y3 = (-0.5d0 + 1.d0 * (j+1) / ny) * side_a
         y4 = y3
         do i=0,nx-1
            x1 = (-0.5d0 + 1.d0 *  i    / nx) * side_b
            x2 = (-0.5d0 + 1.d0 * (i+1) / nx) * side_b
            x3 = x2
            x4 = x1

            z1 = sqrt(r**2 - x1**2 - y1**2) - r
            z2 = sqrt(r**2 - x2**2 - y2**2) - r
            z3 = sqrt(r**2 - x3**2 - y3**2) - r
            z4 = sqrt(r**2 - x4**2 - y4**2) - r
            write (99, 3001) x1, y1, z1
            write (99, 3001) x2, y2, z2
            write (99, 3001) x3, y3, z3
            write (99, 3001) x1, y1, z1
            write (99, *)
            write (99, *)
            write (99, 3001) x1, y1, z1
            write (99, 3001) x3, y3, z3
            write (99, 3001) x4, y4, z4
            write (99, 3001) x1, y1, z1
            write (99, *)
            write (99, *)

            x0(1) = x
            x0(2) = y
            x0(3) = z
            x0    = matmul(cs_rot_m, x0) + oripo_cs
!         write (98, 3001) x0
         enddo
      enddo


      ! side part 1 + 3
      j = 0
      y1 = - 0.5d0 * side_a
      y3 =   0.5d0 * side_a
      do i=0,nx-1
         x1 = (-0.5d0 + 1.d0 *  i    / nx) * side_b
         x2 = (-0.5d0 + 1.d0 * (i+1) / nx) * side_b
         x3 = x2
         x4 = x1

         z3 = sqrt(r**2 - x3**2 - y1**2) - r
         z4 = sqrt(r**2 - x4**2 - y1**2) - r
         z1 = -side_c
         z2 = -side_c

         write (99, 3001) x1, y1, z1
         write (99, 3001) x2, y1, z2
         write (99, 3001) x3, y1, z3
         write (99, 3001) x1, y1, z1
         write (99, *)
         write (99, *)
         write (99, 3001) x1, y1, z1
         write (99, 3001) x3, y1, z3
         write (99, 3001) x4, y1, z4
         write (99, 3001) x1, y1, z1
         write (99, *)
         write (99, *)
         write (99, 3001) x1, y3, z1
         write (99, 3001) x2, y3, z2
         write (99, 3001) x3, y3, z3
         write (99, 3001) x1, y3, z1
         write (99, *)
         write (99, *)
         write (99, 3001) x1, y3, z1
         write (99, 3001) x3, y3, z3
         write (99, 3001) x4, y3, z4
         write (99, 3001) x1, y3, z1
         write (99, *)
         write (99, *)
      enddo

      ! side part 2 + 4
      i = 0
      x2 =   0.5d0 * side_b
      x4 = - 0.5d0 * side_b
      do j=0,ny-1
         y1 = (-0.5d0 + 1.d0 *  j    / ny) * side_a
         y2 = (-0.5d0 + 1.d0 * (j+1) / ny) * side_a
         y3 = y2
         y4 = y1

         z3 = sqrt(r**2 - x2**2 - y3**2) - r
         z4 = sqrt(r**2 - x2**2 - y4**2) - r
         z1 = -side_c
         z2 = -side_c

         write (99, 3001) x2, y1, z1
         write (99, 3001) x2, y2, z2
         write (99, 3001) x2, y3, z3
         write (99, 3001) x2, y1, z1
         write (99, *)
         write (99, *)
         write (99, 3001) x2, y1, z1
         write (99, 3001) x2, y3, z3
         write (99, 3001) x2, y4, z4
         write (99, 3001) x2, y1, z1
         write (99, *)
         write (99, *)
         write (99, 3001) x4, y1, z1
         write (99, 3001) x4, y2, z2
         write (99, 3001) x4, y3, z3
         write (99, 3001) x4, y1, z1
         write (99, *)
         write (99, *)
         write (99, 3001) x4, y1, z1
         write (99, 3001) x4, y3, z3
         write (99, 3001) x4, y4, z4
         write (99, 3001) x4, y1, z1
         write (99, *)
         write (99, *)
      enddo
      close (99)
!      close (98)

      return
3000  format('# grid_id = 0  (irregular (1D) xyz-grid) / ',  ! file header
     &       'x-resolution = (', i5, ') / y-resolution = (', i5, ')'/,
     &       '# number of grid points = ', i12)
3001  format(3e18.10)
      end subroutine make_triangles_sphereX
!------------------------------------------------------------------

!----------------------------------------------------------------
! SUBROUTINE TO CALCULATE THE TRIANGLES THE BLOCK LIMITER IS MADE
! OF AND WRITE THEM TO A FILE FOR EIRENE
!----------------------------------------------------------------
      SUBROUTINE bl_triangles_eirene(i_bl)
      IMPLICIT NONE
      
      integer, dimension(:), allocatable, save :: num_vtx_p	! number of vertices per plane
      integer ::
     &   i, j, k, count_vtx,					! counting variables
     &   pos_vtx,						! position of the vertex in the ordered list
     &   iun = 80, ios
      real*8, dimension(:,:,:), allocatable, save ::
     &   vtx_p,							! vertices with respect to their origin planes (origin plane = first index)
     &   vtx_p_dl,						! dislocated vertices
     &   vtx_p_o						! vertices with respect to their origin planes, sorted in mathematically...
     								! ... negative sense (left-hand rule)
      real*8 ::
     &   abs1, abs2, abs3,					! temporarily used variables
     &   sin_theta, cos_theta,					! sinus and cosinus of the angle between two vectors
     &   sp							! scalar product of two vectors

      real*8, dimension(3) ::
     &   cp,							! cross product of two vectors
     &   v1, v2

      real*8, dimension(:,:), allocatable, save ::
     &   dl_v,
     &   ord_vtx_p_dl						! order value for each vertex on the respective plane

      integer ::
     &   i_bl							! index of the block limiter

      character*120 ::
     &   file_eirene,						! name of the file where the triangles data for EIRENE (ADD_SF_N0) is stored in
     &   file_gnuplot						! name of the file where the triangles data for Gnuplot is stored in

2001  format(6f17.10)
2002  format(3f17.10)
2003  format(t3, a, a)

      IF (i_bl == 1) THEN
         file_eirene = 'bl_triangles.dat'
         file_gnuplot = 'bl_triangles_plot.dat'
          
      ELSE IF (i_bl == 2) THEN
         file_eirene = 'bl_triangles_2.dat'
         file_gnuplot = 'bl_triangles_plot_2.dat'
      ELSE
         print*, 'No correct i_bl!'
         STOP
      END IF

      allocate(dl_v(num_p(i_bl),3))
      allocate(num_vtx_p(num_p(i_bl)))
      num_vtx_p = 0
      dl_v = 0.

! get the number of vertices belonging to each plane
      DO i = 0, num_imp_vtx-1
         DO j = 1, 3
            num_vtx_p(ori_imp_vtx(i,j)) = num_vtx_p(ori_imp_vtx(i,j))+1
         END DO
      END DO

      allocate(vtx_p(num_p(i_bl), maxval(num_vtx_p),3))
      allocate(vtx_p_o(num_p(i_bl), maxval(num_vtx_p),3))
      allocate(vtx_p_dl(num_p(i_bl), maxval(num_vtx_p),3))
      allocate(ord_vtx_p_dl(num_p(i_bl), maxval(num_vtx_p)))
      ord_vtx_p_dl = 0.

! save the vertices with respect to their origin planes
      DO i = 1, num_p(i_bl)
         count_vtx = 1
         DO j = 0, num_imp_vtx-1
            DO k = 1, 3
               IF (ori_imp_vtx(j,k) == i) THEN
                  vtx_p(i,count_vtx,:) = imp_vtx(j,:)
                  count_vtx = count_vtx + 1
               END IF
            END DO
         END DO
      END DO

! calculate the dislocation vectors
      DO i = 1, num_p(i_bl)
         DO j = 1, num_vtx_p(i)
            dl_v(i,:) = dl_v(i,:)+1./num_vtx_p(i)*vtx_p(i,j,:)
         END DO
      END DO

! calculate the dislocated vectors
      DO i = 1, num_p(i_bl)
         DO j = 1, num_vtx_p(i)
            vtx_p_dl(i,j,:) = vtx_p(i,j,:) - dl_v(i,:)
         END DO
      END DO

! find out the order of the vertices of each plane
      DO i = 1, num_p(i_bl)
         abs1 = sqrt(dot_product(vtx_p_dl(i,1,:),vtx_p_dl(i,1,:)))		! absolute value of starting vector
         DO j = 2, num_vtx_p(i)
            v1 = vtx_p_dl(i,1,:)
            v2 = vtx_p_dl(i,j,:)
            CALL cross_product(vtx_p_dl(i,1,:), vtx_p_dl(i,j,:), cp)	! cross product of starting vector and vector whose 'order' shall be found out
            abs2 = sqrt(dot_product(vtx_p_dl(i,j,:),vtx_p_dl(i,j,:)))
            sp = dot_product(vtx_p_dl(i,1,:),vtx_p_dl(i,j,:))
            cos_theta = sp/(abs1*abs2)
            sin_theta = dot_product(cp, normv_p(i_bl,i,:))/(abs1*abs2)
            ord_vtx_p_dl(i,j) = sort_sin_val(sin_theta, cos_theta)			! get the sort value of sin_theta
         END DO
      END DO

! put the vertices into order (mathematically negative sense referring to the normal vector of the plane)
      DO i = 1, num_p(i_bl)
         vtx_p_o(i,1,:) = vtx_p(i,1,:)
         DO j = 2, num_vtx_p(i)
            pos_vtx = 2
            DO k = 2, num_vtx_p(i)
               IF(ord_vtx_p_dl(i,j)<ord_vtx_p_dl(i,k) .and. j.ne.k)THEN
                  pos_vtx = pos_vtx + 1
               END IF
            END DO
            vtx_p_o(i,pos_vtx,:) = vtx_p(i,j,:)
         END DO
      END DO

! write the vertices in the right order into a file for EIRENE
! The points of the triangle have to be written in clockwise order as seen looking
! along the normal vector of the plane from the origin (which points to the outside of the
! block limiter) (mathematically negative sense, right-hand-rule). Otherwise EIRENE will
! start particles inside the block limiter!!!
      open(iun, file=file_eirene, iostat=ios)
      IF (ios /= 0) THEN
         write(*,*) 'SUBROUTINE bl_triangles_eirene:',
     &     'Could not read the footprint grid data file!'
         STOP
      END IF
      DO i = 1, num_p(i_bl)
         DO j = 3, num_vtx_p(i)
            write(iun,fmt=2001) vtx_p_o(i,1,:), vtx_p_o(i,j-1,:)	! clockwise order!
            write(iun,fmt=2002) vtx_p_o(i,j,:) 
         END DO
      END DO

! write the vertices in the right order into a file for plotting with gnuplot
      open(iun, file=file_gnuplot, iostat=ios)
      IF (ios /= 0) THEN
         write(*,*) 'SUBROUTINE bl_triangles_eirene:',
     &     'Could not read the footprint grid data file!'
         STOP
      END IF
      DO i = 1, num_p(i_bl)
         DO j = 3, num_vtx_p(i)
            write(iun,*) vtx_p_o(i,1,:)
            write(iun,*) vtx_p_o(i,j-1,:)
            write(iun,*) vtx_p_o(i,j,:)
            write(iun,*) vtx_p_o(i,1,:)
            write(iun,*)
            write(iun,*)
         END DO
      END DO
      close(iun)

      write(6,fmt=2003) 'Output file for EIRENE (ADD_SF_N0): ', 
     &                   file_eirene
      write(6,fmt=2003) 'Output file for Gnuplot:            ',
     &                   file_gnuplot
      write(6,*)

      RETURN
      END SUBROUTINE bl_triangles_eirene
!--------------------------------------------------------------------------

!----------------------------------------------------------------
! FUNCTION TO CALCULATE THE POSITION OF A VERTEX WITHIN AN
! ORDERED LIST OF VERTICES USING THE SINUS
!----------------------------------------------------------------
      real*8 FUNCTION sort_sin_val(sin_ang, cos_ang)

      real*8, intent(in) :: sin_ang, cos_ang
      real*8, parameter ::
     &   new_0 = -1.d-10,
     &   new_1 = 1. + 1.d-10

!       print*
!       print*, 'sin_ang = ', sin_ang
!       print*, 'cos_ang = ', cos_ang

      IF (new_0<sin_ang .and. sin_ang<new_1 .and. cos_ang>new_0)THEN
         sort_sin_val = sin_ang
      ELSE IF (new_0<sin_ang .and. sin_ang<new_1.and.cos_ang<new_0)THEN
         sort_sin_val = 2.-sin_ang
      ELSE IF (-new_1<sin_ang.and.sin_ang<new_0.and.cos_ang<new_0)THEN
         sort_sin_val = 2.-sin_ang
      ELSE IF (-new_1<sin_ang.and.sin_ang<new_0.and.cos_ang>new_0) THEN
         sort_sin_val = 4.+sin_ang
      ELSE
         print*, 'Something wrong occurred in SR sort_sin_val!!!'
      END IF

!       IF (0..le.sin_ang .and. sin_ang.le.1. .and. cos_ang.ge.0.)THEN
!          sort_sin_val = sin_ang
!       ELSE IF (0..le.sin_ang.and.sin_ang.le.1. .and.cos_ang.lt.0.)THEN
!          sort_sin_val = 2.-sin_ang
!       ELSE IF (-1..le.sin_ang.and.sin_ang.lt.0. .and.cos_ang.le.0.)THEN
!          sort_sin_val = 2.-sin_ang
!       ELSE IF (-1..le.sin_ang.and.sin_ang.lt.0..and.cos_ang.gt.0.)THEN
!          sort_sin_val = 4.+sin_ang
!       ELSE
!          print*, 'Something wrong occurred in SR sort_sin_val!!!'
!       END IF
      
      RETURN
      END FUNCTION sort_sin_val 
!-----------------------------------------------------------------

!------------------------------------------------------------------
! SUBROUTINE TO CALCULATE IF A GIVEN POINT IS ON THE INSIDE OF THE 
! BLOCK LIMITER
!------------------------------------------------------------------
      SUBROUTINE check_point_behind_bl(x, limit, i_bl, result)
      IMPLICIT NONE

      real*8, dimension(3), intent(in) :: x	! given point in cartesian coordinates
      real*8, intent(in) :: limit		! value considered as zero
      real*8 :: u				! temporary storing variable
      integer, intent(in) :: i_bl		! index of block limiter
      integer :: i				! counting variable
      integer :: iun				! file number
      logical, intent(out) :: result		! if the given point is on the inside then result is true
      logical :: check_position			! if the given point is on the outside then check_position is true

      check_position = .false.

! check if point x is on the outside of the block limiter
!   1. check planar elements
      DO i=1, num_p(i_bl)
         u = dot_product(x, normv_p(i_bl,i,:))-dist_p_oripocs(i_bl,i)	! u = distance between plane i and point x
         IF (u > limit) THEN						! is the point x on the outside of the plane?
            check_position = .true.
            EXIT
         END IF
      END DO

!   2. check spherical elements
      DO i=1, num_s(i_bl)
         u = dot_product(x-center_s(i_bl,i,:), x-center_s(i_bl,i,:))    ! u = distance between point x and center of
         u = sqrt(u)							! spherical element
         IF (u > radius_s(i_bl,i)) THEN       ! is the point x on the outside of the sphere?
            check_position = .true.
            EXIT
         ENDIF
      ENDDO


! set result
      IF (.not.check_position) THEN
         result = .true.
!         open(iun, file="cutaway_bl.dat", access="append")
!         write(iun, *) x
!         close(iun)
      ELSE
         result = .false.
      END IF

      RETURN
      END SUBROUTINE check_point_behind_bl
!------------------------------------------------------------------

!------------------------------------------------------------------
! FUNCTION TO CALCULATE THE INTERSECTION POINTS OF A LINE SEGMENT 
! GEOMETRIC FORM CONSISTING OF PLANES
!------------------------------------------------------------------
      logical FUNCTION check_intersection(xi, xe, ip, lambda_min, 
     &                                    i_ip, i_bl)

      real*8, dimension(3), intent(in) :: xi, xe	! initial and ending point of the line segment
      real*8, dimension(3), intent(inout) :: ip		! intersection point
      real*8, intent(inout) :: lambda_min		! value of the smallest lambda (to chose the right ip)
      real*8, dimension(max_num_p, 3) :: ip_list	! list of intersection points, temporary saving variable
      real*8, dimension(max_num_p) :: lambda		! parameter of the intersection point
      real*8, dimension(3) :: s				! slope of the line segment
      real*8, parameter :: limit = 10.d-10		! value considered as zero
      real*8 :: u, v, p, q, A, s0, dx(3)		! temporary used variables
!       integer, intent(inout) :: num_reach		! number of line segments out of reach
!       integer, intent(inout) :: num_ip		! number of intersection points
!       integer, intent(inout) :: num_call_ci		! number of calls of subroutine check_intersection
! if these variables should be used they have to be added to the list of arguments of the functions
      integer :: num_rel				! number of relevant planes
      integer, intent(out) :: i_ip			! indice of the intersection point with the smallest lambda
      integer, intent(in) :: i_bl			! indice of the block limiter
      integer :: i, j, k				! counting variables
      logical, dimension(max_num_p) :: rel		! relevance of the plane
      logical :: ls_is_part_of_plane = .false.		! is the line segment part of the plane?

      real*8, dimension(max_num_s)  :: lambda_s	! relative intersection points with the spheres
      logical, dimension(max_num_s) :: rel_s            ! relevance of the sphere
      real*8, dimension(max_num_s, 3) :: ip_list_s	! list of intersection points with spheres, temporary saving variable
      real*8 :: l1, l2


! initialization of some variables
      ip = (/ 666., 666., 666. /)
      lambda_min = 1.
      s = xe - xi
      check_intersection = .false.
      rel = .true.
      ip_list = 0.
      lambda = 0.
      i_ip = 0
      num_rel = 0
      k = 1
      ! for intersection with spheres
      lambda_s  = 0.d0
      rel_s     = .true.
      ip_list_s = 0.d0

! setting control variable
!       num_call_ci = num_call_ci + 1

! check whether the line segment is in the nearer circumcircle of the block limiter
      u = sqrt(dot_product(xi-center_bl(i_bl,:), xi-center_bl(i_bl,:)))
      IF (u > sqrt(dot_product(s, s)) + radius_bl(i_bl) + limit) THEN
!          print*, 'The line segment is >>out of reach<<!'
         RETURN
      END IF

! calculate lambda, check if the slope of the line segment is parallel to one of the planes
      DO i = 1, num_p(i_bl)
         u = dot_product(s, normv_p(i_bl, i,:))
         v = dot_product((xi - oripo_p(i_bl, i,:)), normv_p(i_bl, i,:))
         IF (u .ne. 0.) THEN	
            lambda(i) = -v/u
         ELSE IF (u .eq. 0. .and. v == 0.) THEN		! these points have to be considered, but they will appear ...
            ls_is_part_of_plane = .true.		! as an intersection point with another plane
            rel(i) = .false.
         ELSE
            rel(i) = .false.
         END IF
      END DO

! calculate lambda_s
      do i = 1, num_s(i_bl)
         q  = dot_product(xi-center_s(i_bl,i,:), xi-center_s(i_bl,i,:))
         s0 = sqrt(dot_product(s, s))
         p  = 2.d0 * dot_product(xi-center_s(i_bl,i,:), s) / s0
         A  = p**2 / 4.d0 - q + radius_s(i_bl,i)**2
         if (A.ge.0.d0) then
            A = sqrt(A)
            !lambda_s(i_bl,1) = (- p / 2.d0 + A) / s0
            !lambda_s(i_bl,2) = (- p / 2.d0 - A) / s0

            ! l2 < l1
            l1 = (- p / 2.d0 + A) / s0
            l2 = (- p / 2.d0 - A) / s0


            if (l2.ge.0.d0 .and. l2.le.1.d0) then
               lambda_s(i) = l2
            else if (l1.ge.0.d0 .and. l1.le.1.d0) then
               lambda_s(i) = l1
            else
               rel_s(i) = .false.
            endif
         else
            rel_s(i) = .false.
         endif
      enddo


! check if 0 <= lambda <= 1, eliminate multiply appearing lambdas
      DO i = 1, num_p(i_bl)
         IF (rel(i)) THEN

            IF (0. <= lambda(i) .and. lambda(i) <= 1.) THEN
               DO j=1, i-1
                  IF (lambda(i)==lambda(j) .and. rel(j)) THEN
                     rel(i) = .false.
                     EXIT
                  END IF
               END DO
            ELSE IF (lambda(i) < 0. .or. 1. < lambda(i)) THEN
               rel(i) = .false.
            END IF

         END IF
      END DO
! check if 0 <= lambda_s <= 1, eliminate multiply appearing lambdas
!      DO i = 1, num_s(i_bl)
!         IF (rel_s(i) == .true.) THEN
!
!            IF (0. <= lambda_s(i) .and. lambda_s(i) <= 1.) THEN
!               DO j=1, i-1
!                  IF (lambda_s(i)==lambda_s(j).and.rel_s(j)==.true.)THEN
!                     rel_s(i) = .false.
!                     EXIT
!                  END IF
!               END DO
!            ELSE IF (lambda_s(i) < 0. .or. 1. < lambda_s(i)) THEN
!               rel_s(i) = .false.
!            END IF
!
!         END IF
!      END DO

! calculate coordinates of the intersection point
      DO i = 1, num_p(i_bl)
         IF (rel(i)) THEN
           ip_list(i, :) = xi + lambda(i)*s
           num_rel = num_rel+1
         END IF
      END DO
      DO i = 1, num_s(i_bl)
         IF (rel_s(i)) THEN
           ip_list_s(i, :) = xi + lambda_s(i)*s
           num_rel = num_rel+1
         END IF
      END DO
      !if (num_rel.gt.0) write (6, *) num_rel

! check whether the intersection points are relevant or not (lie outside of the cuboid or not)
      ! 1st: check planes
      DO i=1, num_p(i_bl)
         ! check distance to all planes
         IF (rel(i)) THEN
            DO j=1, num_p(i_bl)
               v = dot_product(ip_list(i, :), normv_p(i_bl, j,:))
               u = v - dist_p_oripocs(i_bl, j)	! u = distance between plane j and intersection point i
               IF (u > limit) THEN				! is the intersection point on the outside of the cuboid?
                  rel(i) = .false.
                  num_rel = num_rel-1
                  EXIT
               END IF
            END DO
         END IF
         ! check distance to to all spheres
         if (rel(i)) then
            do j=1, num_s(i_bl)
               dx = ip_list(i, :) - center_s(i_bl, j, :)	! distance between intersection point i and center of sphere j
               u  = dot_product (dx, dx) - radius_s(i_bl, j)**2	! distance to sphere surface
               if (u > limit) then
                  rel(i) = .false.
                  num_rel = num_rel-1
                  exit
               endif
            enddo
         endif
         IF (rel(i)) THEN				! get the smallest lambda and the associated index of ...
            IF (lambda(i) <= lambda_min) THEN			! the intersection point
               lambda_min = lambda(i)
!               i_ip = i
            END IF
         END IF
      END DO
      ! 2nd: check spheres
      do i=1, num_s(i_bl)
         ! check distance to all planes
         IF (rel_s(i)) THEN
            DO j=1, num_p(i_bl)
               v = dot_product(ip_list_s(i, :), normv_p(i_bl, j,:))
               u = v - dist_p_oripocs(i_bl, j)	! u = distance between plane j and intersection point i
               IF (u > limit) THEN				! is the intersection point on the outside of the cuboid?
                  rel(i) = .false.
                  num_rel = num_rel-1
                  EXIT
               END IF
            END DO
         END IF
         ! check distance to to all spheres
         if (rel_s(i)) then
            do j=1, num_s(i_bl)
               dx = ip_list_s(i, :) - center_s(i_bl, j, :)	! distance between intersection point i and center of sphere j
               u  = dot_product (dx, dx) - radius_s(i_bl, j)**2	! distance to sphere surface
               if (u > limit) then
                  rel(i) = .false.
                  num_rel = num_rel-1
                  !write (6, *) 'new num_rel = ', num_rel
                  exit
               endif
            enddo
         endif
         IF (rel_s(i)) THEN				! get the smallest lambda and the associated index of ...
            IF (lambda_s(i) <= lambda_min) THEN			! the intersection point
               lambda_min = lambda_s(i)
!               i_ip_s = i
            END IF
         END IF
      enddo

! output
      IF (num_rel == 0) THEN	
         RETURN
      ELSE
         check_intersection = .true.
!         ip = ip_list(i_ip, :)
         ip = xi + lambda_min * s
         !write (6, *) lambda_min
      END IF

      RETURN
      END FUNCTION check_intersection
!-----------------------------------------------------------------------

!------------------------------------------------------------------
! SUBROUTINE TO READ THE BLOCK LIMITER DATA
!------------------------------------------------------------------
      SUBROUTINE read_block_limiter_data(data_file)	
      IMPLICIT NONE	

      character*255, intent(in) :: data_file
      integer :: iun = 33, iun2 = 32, ios

0815  format(t9, a, a, a, a, a, a, a, a, a, a, a)	! standard output format
1000  format(t9, '|' a10, '   =   ', f10.2, a4, ' |')
1001  format(t9, '|' a10, '   =   ', i10, '     |')
1002  format(t9, 'BLOCK LIMITER DATA') 
1003  format(t9, '|',tr32, '|')
2000  format(t9, '==================================')


! get the data from the block limiter file
      open(iun, file=data_file, iostat=ios)
      IF (ios /= 0) THEN
         write(*, 0815) 'SUBROUTINE read_block_limiter_data:',
     &     'Could not read the block limiter data file!'
         write(6, *) data_file
         STOP
      END IF
      read(iun, BlockLimiterData)
      close(iun)

      write(*,*)
      write(*, fmt=1002)
      write(*, fmt=2000)
      write(*, fmt=1001) 'form type', form_type
      write(*, fmt=1003)
      write(*, fmt=1000) 'side a', side_a, ' cm'
      write(*, fmt=1000) 'side b', side_b, ' cm'
      write(*, fmt=1000) 'side c', side_c, ' cm'
      write(*, fmt=1000) 'length d', length_d, ' cm'
      write(*, fmt=1003)
      write(*, fmt=1000) 'phi', phi, ' deg'
      write(*, fmt=1000) 'theta', theta_bl, ' deg'
      write(*, fmt=1000) 'r_min', r_min, ' cm'
      write(*, fmt=1003)
      write(*, fmt=1000) 'alpha', alpha, ' deg'
      write(*, fmt=1000) 'beta', beta, ' deg'
      write(*, fmt=2000)
      write(*,*)


! convert the angles from degree to radian
      phi = phi*pi2/360.
      theta_bl = theta_bl*pi2/360.
      alpha = alpha*pi2/360.
      beta = beta*pi2/360.

      END SUBROUTINE read_block_limiter_data
!------------------------------------------------------------------

!-----------------------------------------------------------------------
! SUBROUTINE TO ALLOCATE THE ARRAYS TO STORE THE BLOCK LIMITER DATA IN
!-----------------------------------------------------------------------
      SUBROUTINE allocate_bl_arrays(num_bl)			
      IMPLICIT NONE
      
      integer, intent(in) :: num_bl

! allocation of the needed arrays to store the block limiter data in
      allocate(num_p(num_bl))
      allocate(oripo_p(num_bl, max_num_p, 3))
      allocate(normv_p(num_bl, max_num_p, 3))
      allocate(dist_p_oripocs(num_bl, max_num_p))
      allocate(center_bl(num_bl, 3))				! center or the block limiter´
      allocate(radius_bl(num_bl))					! radius of the sphere containing the block limiter

      allocate(num_s(num_bl))
      allocate(radius_s(num_bl, max_num_s))
      allocate(center_s(num_bl, max_num_s, 3))
      num_s = 0
      
      END SUBROUTINE allocate_bl_arrays
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! SUBROUTINE TO BUILD UP THE BLOCK LIMITER WITH THE GIVEN DATA
!-----------------------------------------------------------------------
      SUBROUTINE setup_block_limiter(data_file, i_bl)			
      IMPLICIT NONE						

      character*255, intent(in) :: data_file		! name of the file where the block limiter data is stored
      integer, intent(in) :: i_bl			! index of the block limiter which should be build up
      integer :: i					! integer for testing 

1000  format(t8, 'Error while reading ', a, a, a, '!')	! error prompts if input file has been
0815  format(t9, a, a, a, a, a, a, a, a, a, a, a)	! standard output format
1001  format(t8, 'You have manipulated the input file in such a way'/, ! manipulated in a bad way
     &       t8, 'that it could not be read properly.'/, 
     &       t8, 'Please refer to the original form and try again.')

! get the block limiter data
      CALL read_block_limiter_data(data_file)

      SELECT CASE (form_type)
!-----------
         CASE (1)
         write(*, fmt=0815) 'The block limiter is of type CUBOID.'
         CALL build_cuboid(i_bl)						! get the coordinates of the cuboid
!-----------
         CASE (2)
         write(*, fmt=0815) 'The block limiter is of type WEDGE.'	! get the coordinates of the wedge
         CALL build_wedge(i_bl)
!-----------
         CASE (3)
         write(*, fmt=0815) 'The block limiter is of type ROOF.'	! get the coordinates of the roof
         CALL build_roof(i_bl)
!-----------
         CASE (4)
         write(*, fmt=0815) 'The block limiter is of type SPHERE.'	! get the coordinates of the roof
         CALL build_sphereX(i_bl)
!-----------
         CASE DEFAULT							! error output
         write(*, fmt=0815) 'SUBROUTINE setup_block_limiter:',
     &      'You have put in a wrong number for ',
     &      'choosing the geometric form of the block limiter.'
         STOP
!-----------
      END SELECT

! HAS TO BE CHANGED FOR ADDITIONAL BLOCK LIMITER
      if (form_type .ne. 4) call gnuplot_data_block_limiter(i_bl)

      RETURN
      END SUBROUTINE setup_block_limiter
!------------------------------------------------------------------

!-----------------------------------------------------------------------
! SUBROUTINE FOR BUILDING UP A CUBOID
!-----------------------------------------------------------------------
      SUBROUTINE build_cuboid(i_bl)
      IMPLICIT NONE

      integer, intent(in) :: i_bl			! index of block limiter
      integer :: i
      
      num_p(i_bl) = 6

! defining the normal vectors and origin points of the planes in the local coordinate system
! NORMAL VECTORS HAVE TO BE DEFINED AS NORMALIZED VECTORS !!!
      normv_p(i_bl, 1,:) = e_x
      normv_p(i_bl, 3,:) = -normv_p(i_bl, 1,:)
      normv_p(i_bl, 2,:) = e_y
      normv_p(i_bl, 4,:) = -normv_p(i_bl, 2,:)
      normv_p(i_bl, 5,:) = e_z
      normv_p(i_bl, 6,:) = -normv_p(i_bl, 5,:)

      oripo_p(i_bl, 1, :) = side_b/2 * normv_p(i_bl, 1,:)
      oripo_p(i_bl, 3, :) = side_b/2 * normv_p(i_bl, 3,:)
      oripo_p(i_bl, 2, :) = side_a/2 * normv_p(i_bl, 2,:)
      oripo_p(i_bl, 4, :) = side_a/2 * normv_p(i_bl, 4,:)
      oripo_p(i_bl, 5, :) = (/ 0., 0., 0. /)
      oripo_p(i_bl, 6, :) = side_c * normv_p(i_bl, 6,:)

! get data to check if a point is in the nearer circumcircle of the block limiter
      center_bl(i_bl, :) = side_c/2 * normv_p(i_bl, 6, :)
      radius_bl(i_bl) = sqrt(side_a**2 + side_b**2 + side_c**2)/2.

! get the coordinates of the cuboid in the cs of the reactor
      CALL cs_transformation(i_bl)

! calculation of the distance between the plane and the origin point of the cs (Hesse normal form) ...
! (needed to check if a vertex lies on the surface of the block limiter or not)
      DO i=1, num_p(i_bl)
         dist_p_oripocs(i_bl, i)
     &      = dot_product(normv_p(i_bl, i,:), oripo_p(i_bl, i,:))
      END DO

      RETURN
      END SUBROUTINE build_cuboid
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! SUBROUTINE FOR BUILDING UP A WEDGE
!-----------------------------------------------------------------------
      SUBROUTINE build_wedge(i_bl)
      IMPLICIT NONE

      real*8, dimension(3, 3) :: rot_mat
      integer, intent(in) :: i_bl			! index of block limiter
      integer :: i

! allocation of the needed arrays to store the block limiter data in
      num_p(i_bl) = 6							

! defining the normal vectors and origin points of the planes in the local coordinate system
! NORMAL VECTORS HAVE TO BE DEFINED AS NORMALIZED VECTORS !!!
      normv_p(i_bl, 1,:) = e_x
      normv_p(i_bl, 3,:) = -normv_p(i_bl, 1,:)
      normv_p(i_bl, 2,:) = e_y
      normv_p(i_bl, 4,:) = -normv_p(i_bl, 2,:)
      normv_p(i_bl, 5,:) = e_z
      normv_p(i_bl, 6,:) = -normv_p(i_bl, 5,:)

      oripo_p(i_bl, 1,:) = side_b/2. * normv_p(i_bl, 1,:)
      oripo_p(i_bl, 3,:) = side_b/2. * normv_p(i_bl, 3,:)
      oripo_p(i_bl, 2,:) = side_a/2. * normv_p(i_bl, 2,:)
      oripo_p(i_bl, 4,:) = side_a/2. * normv_p(i_bl, 4,:)
      oripo_p(i_bl, 5,:) = (/ 0.d0, (-side_a/2), (-length_d) /)
      oripo_p(i_bl, 6, :) = side_c*normv_p(i_bl, 6,:)

! get data to check if a point is in the nearer circumcircle of the block limiter
      center_bl(i_bl, :) = (side_c+length_d)/2.*normv_p(i_bl, 6,:)
      radius_bl(i_bl)=sqrt(side_a**2+side_b**2+(side_c+length_d)**2)/2.

! build the inclined plane of the wedge
      alpha_w = atan(length_d/side_a)
      CALL get_rotation_matrix(-alpha_w*(/ 1., 0., 0. /), rot_mat)
      normv_p(i_bl, 5,:) = matmul(rot_mat, normv_p(i_bl, 5,:))

! get the coordinates of the cuboid in the cs of the reactor
      CALL cs_transformation(i_bl)

! calculation of the distance between the plane and the origin point of the cs (Hesse normal form) ...
! (needed to check if a vertex lies on the surface of the block limiter or not)
      DO i=1, num_p(i_bl)
         dist_p_oripocs(i_bl, i)
     &      = dot_product(normv_p(i_bl, i,:), oripo_p(i_bl, i,:))
      END DO

      RETURN
      END SUBROUTINE build_wedge
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! SUBROUTINE FOR BUILDING UP A ROOF
!-----------------------------------------------------------------------
      SUBROUTINE build_roof(i_bl)
      IMPLICIT NONE

      real*8, dimension(3, 3) :: rot_mat
      real*8, dimension(3) :: test_v
      integer :: i
      integer, intent(in) :: i_bl			! index of block limiter

! allocation of the needed arrays to store the block limiter data in
      num_p(i_bl) = 7							

! defining the normal vectors and origin points of the planes in the local coordinate system
! NORMAL VECTORS HAVE TO BE DEFINED AS NORMALIZED VECTORS !!!
      normv_p(i_bl, 1,:) = e_x
      normv_p(i_bl, 3,:) = -normv_p(i_bl, 1,:)
      normv_p(i_bl, 2,:) = e_y
      normv_p(i_bl, 4,:) = -normv_p(i_bl, 2,:)
      normv_p(i_bl, 5,:) = e_z
      normv_p(i_bl, 6,:) = e_z
      normv_p(i_bl, 7,:) = -normv_p(i_bl, 5,:)

      oripo_p(i_bl, 1,:) = side_b/2. * normv_p(i_bl, 1,:)
      oripo_p(i_bl, 3,:) = side_b/2. * normv_p(i_bl, 3,:)
      oripo_p(i_bl, 2,:) = side_a/2. * normv_p(i_bl, 2,:)
      oripo_p(i_bl, 4,:) = side_a/2. * normv_p(i_bl, 4,:)
      oripo_p(i_bl, 5,:) = (/ 0.d0, -side_a/2., -length_d /)
      oripo_p(i_bl, 6,:) = (/ 0.d0, side_a/2., -length_d /)
      oripo_p(i_bl, 7,:) = side_c * normv_p(i_bl, 7,:)

! get data to check if a point is in the nearer circumcircle of the block limiter
      center_bl(i_bl, :) = (side_c+length_d)/2.*normv_p(i_bl, 7,:)
      radius_bl(i_bl) = sqrt(side_a**2 + side_b**2 
     &                    + (side_c+length_d)**2)/2.

! build the inclined planes of the roof
      alpha_r = atan(length_d/(side_a/2))
      CALL get_rotation_matrix(-alpha_r*(/ 1., 0., 0. /), rot_mat)
      normv_p(i_bl, 5,:) = matmul(rot_mat, normv_p(i_bl, 5,:))
      CALL get_rotation_matrix(alpha_r*(/ 1., 0., 0. /), rot_mat)
      normv_p(i_bl, 6,:) = matmul(rot_mat, normv_p(i_bl, 6,:))

! get the coordinates of the cuboid in the cs of the reactor
      CALL cs_transformation(i_bl)

! calculation of the distance between the plane and the origin point of the cs (Hesse normal form) ...
! (needed to check if a vertex lies on the surface of the block limiter or not)
      DO i=1, num_p(i_bl)
         dist_p_oripocs(i_bl, i)
     &      = dot_product(normv_p(i_bl, i,:), oripo_p(i_bl, i,:))
      END DO

      RETURN
      END SUBROUTINE build_roof
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! SUBROUTINE FOR BUILDING UP A SPHERE-LIMITER
!-----------------------------------------------------------------------
      SUBROUTINE build_sphereX(i_bl)
      IMPLICIT NONE

      integer, intent(in) :: i_bl			! index of block limiter
      integer :: i
      real*8  :: f, df
      
      num_p(i_bl) = 5

! defining the normal vectors and origin points of the planes in the local coordinate system
! NORMAL VECTORS HAVE TO BE DEFINED AS NORMALIZED VECTORS !!!
      normv_p(i_bl, 1,:) = e_x
      normv_p(i_bl, 3,:) = -normv_p(i_bl, 1,:)
      normv_p(i_bl, 2,:) = e_y
      normv_p(i_bl, 4,:) = -normv_p(i_bl, 2,:)
      normv_p(i_bl, 5,:) = -e_z

      oripo_p(i_bl, 1, :) = side_b/2 * normv_p(i_bl, 1,:)
      oripo_p(i_bl, 3, :) = side_b/2 * normv_p(i_bl, 3,:)
      oripo_p(i_bl, 2, :) = side_a/2 * normv_p(i_bl, 2,:)
      oripo_p(i_bl, 4, :) = side_a/2 * normv_p(i_bl, 4,:)
      oripo_p(i_bl, 5, :) = side_c * normv_p(i_bl, 5,:)

! get data to check if a point is in the nearer circumcircle of the block limiter
      center_bl(i_bl, :) = side_c/2 * normv_p(i_bl, 5, :)
      radius_bl(i_bl) = sqrt(side_a**2 + side_b**2 + side_c**2)/2.


! add the spherical surface
      f = (side_a**2 + side_b**2) / length_d / 8.d0 - length_d / 2.d0

      num_s(i_bl) = 1
      radius_s(i_bl, 1)    = length_d + f
      center_s(i_bl, 1, :) = (/ 0.d0, 0.d0, -radius_s(i_bl, 1) /)

      !write (6, *) 'radius = ', radius_s(i_bl, 1)

      ! adapt the radius of the cirumsphere of the block limiter
      df = abs(f - side_c)
      do i=1,num_s(i_bl)
         if (radius_bl(i_bl) .lt. df+radius_s(i_bl,i)) then
            radius_bl(i_bl) = df+radius_s(i_bl,i)
         endif
      enddo



! get the coordinates of the cuboid in the cs of the reactor
      CALL cs_transformation(i_bl)

! calculation of the distance between the plane and the origin point of the cs (Hesse normal form) ...
! (needed to check if a vertex lies on the surface of the block limiter or not)
      DO i=1, num_p(i_bl)
         dist_p_oripocs(i_bl, i)
     &      = dot_product(normv_p(i_bl, i,:), oripo_p(i_bl, i,:))
      END DO
      !write (6, *) 'center_s = ', center_s(i_bl, 1, :)

      !call make_footprint_grid_sphereX

      RETURN
      END SUBROUTINE build_sphereX
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! SUBROUTINE TO TRANSFORM COORDINATES FROM THE INTERNAL BLOCK LIMITER CS
! CS2 TO THE GLOBAL CS OF THE REACTOR
!-----------------------------------------------------------------------
      SUBROUTINE cs_transformation(i_bl)
      IMPLICIT NONE

      integer, intent(in) :: i_bl				! index of block limiter
      real*8, dimension(3) :: oripo_cs					! origin point of cs2 in cs1
      real*8, dimension(3, 3) :: rot_mat				! final rotation matrix between cs1 and cs2
      integer :: i							! counting variable
      integer :: iun = 40						!
      integer :: ios

! get the coordinates of the reference point of the block limiter
      CALL get_cs_transformation_data(oripo_cs, rot_mat)
      
! transform all coordinates
      DO i = 1, num_p(i_bl)
         normv_p(i_bl, i,:) = matmul(rot_mat, normv_p(i_bl, i,:))
      END DO

      DO i = 1, num_p(i_bl)
         oripo_p(i_bl,i,:)= matmul(rot_mat,oripo_p(i_bl,i,:))+oripo_cs
      END DO

      center_bl(i_bl,:) = matmul(rot_mat, center_bl(i_bl,:))+oripo_cs

      do i = 1, num_s(i_bl)
        center_s(i_bl,i,:)= matmul(rot_mat, center_s(i_bl,i,:))+oripo_cs
      enddo

      RETURN
      END SUBROUTINE cs_transformation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! SUBROUTINE TO CALCULATE THE POSITION OF THE ORIGIN POINT OF CS2 IN CS1
! AND TO GET THE TRANSFORMATION MATRIX FROM CS2 TO CS1
!-----------------------------------------------------------------------
      SUBROUTINE get_cs_transformation_data(oripo, rot_mat)
      IMPLICIT NONE 

      real*8, dimension(3), intent(out) :: oripo			! origin point of cs2
      real*8, dimension(3) :: v1 = 0., v2 = 0., v3 = 0.			! temporarily used vectors
      real*8 :: theta2 = 0.							! angle needed to calculate v2
      real*8, dimension(3) :: y_axis_cs2		! used to get rot_mat2
      real*8, dimension(3, 3) :: 
     &   rot_mat1 = 0., rot_mat2 = 0., rot_mat3 = 0.		! temporary used matrix for rotation
      real*8, dimension(3, 3), intent(out) :: rot_mat				! final rotation matrix between cs1 and cs2
      integer :: iun = 28, ios
      
      y_axis_cs2 = (/ 0., 1., 0. /)
      
      CALL spheri_to_cart(R_max, pi/2., phi, v1)			! calculating v1

      theta2 = pi/2. - theta_bl
      IF (theta2 <= pi/2. .and. theta2 >= 0.) THEN			! calculating v2
         CALL spheri_to_cart(r_min, theta2, phi, v2)
      ELSE IF (theta2 < 0. .and. theta2 >= -pi) THEN
         CALL spheri_to_cart(r_min, -theta2, phi+pi, v2)
      ELSE IF (theta2 < -pi .and. theta2 >= -3.*pi/2.) THEN
         CALL spheri_to_cart(r_min, theta2+pi2, phi, v2)
      ELSE
         print*, 'SR get_oripo_cs:'
         print*, 'The angle theta exceeds 360 degrees!'
         print*, 'Please control the block limiter data file!'
      END IF
      oripo = v1 + v2						! calculating the origin point of cs2      
      
! write v1 and v2 and the joint between the two into a file for plotting
!      open(iun, file='v1_plot.dat', iostat=ios)
!      IF (ios /= 0) THEN
!         print*, 'v1_plot.dat could not be created!'
!         STOP
!      END IF
!      write(iun, *) (/ 0., 0., 0. /)
!      write(iun, *) v1
!      close(iun)
!
!      open(iun, file='v2_plot.dat')
!      write(iun, *) v1
!      write(iun, *) (v1+v2)
!      close(iun)
!
!      open(iun, file='joint_v1v2_plot.dat')
!      write(iun, *) v1
!      close(iun)

      v3 = (/ 0.d0, 0.d0, -phi /)
      CALL get_rotation_matrix(v3, rot_mat1)				! get matrix to rotate x' and y' around z by phi
      y_axis_cs2 = matmul(rot_mat1, y_axis_cs2)				! get matrix to rotate x' and z' around y'-axis by theta
      CALL get_rotation_matrix(y_axis_cs2*(theta_bl+pi/2.+alpha),		
     &   rot_mat2)
      v3 = v2*(beta/r_min)
      CALL get_rotation_matrix(v3, rot_mat3)				! get matrix  to rotate x', y', z' around v2 by beta
      rot_mat = matmul(rot_mat2, rot_mat1)				! get the complete rotation matrix
      rot_mat = matmul(rot_mat3, rot_mat)
     
      RETURN
      END SUBROUTINE get_cs_transformation_data
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! SUBROUTINE TO CONVERT A VECTOR FROM CYLINDRICAL TO CARTESIAN COORDINATES
!-----------------------------------------------------------------------
      SUBROUTINE cyndri_to_cart(r, phi2, z, converted_v)
      IMPLICIT NONE

      real*8, intent(in) :: r, phi2, z
      real*8, dimension(3), intent(out) :: converted_v

      converted_v(1) = r*cos(phi2)
      converted_v(2) = r*sin(phi2)
      converted_v(3) = z

      RETURN
      END SUBROUTINE cyndri_to_cart
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! SUBROUTINE TO CONVERT A VECTOR FROM SPHERICAL TO CARTESIAN COORDINATES
!-----------------------------------------------------------------------
      SUBROUTINE spheri_to_cart(r, theta2, phi2, converted_v)
      IMPLICIT NONE

      real*8, intent(in) :: r, theta2, phi2
      real*8, dimension(3), intent(out) :: converted_v

      IF (theta2>pi) THEN
         print*, 'SUBROUTINE spheri_to_cart:'
         print*, 'The angle theta exceeds 180 degrees!'
         print*, 'The vector could not be converted!'
         STOP
      END IF

      converted_v(1) = r*sin(theta2)*cos(phi2)
      converted_v(2) = r*sin(theta2)*sin(phi2)
      converted_v(3) = r*cos(theta2)

      RETURN
      END SUBROUTINE spheri_to_cart
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! SUBROUTINE TO CALCULATE A ROTATION MATRIX 
!-----------------------------------------------------------------------
      SUBROUTINE get_rotation_matrix(W, R)

*  Given:
*     W        d(3)      rotation vector (Note 1)
*  Returned:
*     R        d(3,3)    rotation matrix
*
*  Notes:
*
*  1) A rotation matrix describes a rotation through some angle about
*     some arbitrary axis called the Euler axis.  The "rotation vector"
*     supplied to this routine has the same direction as the Euler axis,
*     and its magnitude is the angle in radians.
*
*  2) If W is null, the unit matrix is returned.
*
*  3) The vector rotates counterclockwise as seen looking along the
*     rotation vector from the origin (mathematically negative sense,
*     left-hand-rule)
*
*  This revision:  2008 May 10
*
*  Copyright (C) 2008 IAU SOFA Review Board.  See notes at end.

      IMPLICIT NONE

      real*8, dimension(3) :: W
      real*8, dimension(3, 3) :: R
      real*8 :: X, Y, Z, PHI, S, C, F

*  Euler angle (magnitude of rotation vector) and functions.
      X = W(1)
      Y = W(2)
      Z = W(3)
      PHI = SQRT(X*X + Y*Y + Z*Z)
      S = SIN(PHI)
      C = COS(PHI)
      F = 1D0 - C

*  Euler axis (direction of rotation vector), perhaps null.
      IF ( PHI .NE. 0D0 ) THEN
         X = X / PHI
         Y = Y / PHI
         Z = Z / PHI
      END IF

*  Form the rotation matrix.
      R(1,1) = X*X*F + C
      R(1,2) = X*Y*F + Z*S
      R(1,3) = X*Z*F - Y*S
      R(2,1) = Y*X*F - Z*S
      R(2,2) = Y*Y*F + C
      R(2,3) = Y*Z*F + X*S
      R(3,1) = Z*X*F + Y*S
      R(3,2) = Z*Y*F - X*S
      R(3,3) = Z*Z*F + C

*  Finished.

      RETURN
      END SUBROUTINE get_rotation_matrix
!-----------------------------------------------------------------------

!--------------------------------------------------------------------------
! SUBROUTINE TO GET THE DATA OF THE BLOCK LIMITER FOR PLOTTING WITH GNUPLOT
!--------------------------------------------------------------------------
      SUBROUTINE gnuplot_data_block_limiter(i_bl)
      IMPLICIT NONE

      real*8, dimension(:,:), allocatable, save :: vtx		! intersection points of the straight lines
      real*8, dimension(:,:), allocatable, save :: v1, v2			! vectors of the parameter form of the plane
      real*8, dimension(:,:), allocatable, save :: oripo_sl, slo_sl	! origin point and slope of the straight lines
      
      real*8, dimension(3) :: p1, p2				! temporarily used vectors
      real*8, parameter :: limit = 10.d-10			! value considered as zero
      real*8 :: lambda						! parameter of the straight line
      real*8 :: b, c, d, u, v					! temporary storage variable
      integer, dimension(max_num_p**3, 2) :: origin_sl		! origin planes of the straight lines
      integer, dimension(max_num_p**3, 3) :: origin_vtx		! origin planes of the vertices
      integer, parameter :: iun_2 = 41, iun_3 = 40		! unit number of file where the plotting data is stored
      integer, intent(in) :: i_bl				! index of the block limiter which should be plotted
      integer :: num_sl						! number of straight lines
      integer :: num_vtx					! number of vertices of the straight lines
      integer :: initial_vtx_index				! the index of the vertex where we started following ...
      								! along the side-limiting traverse of the form
      integer :: x						! index of the new adjacent vertex
      integer :: i_1						! indices of the origin plane of the vertex to which ...
      								! the traverse belongs
      integer :: i_2						! first: second index of the former vertex, ...
                                                        	! then: index of an origin plane of the adjacent vertex 
      integer :: i, j, k, l, e, f				! counting variable
      integer :: alloc_error = 666				! to check if the allocation has been successful
      integer :: ios = 666					! needed for checking if the file could be openend
      logical, dimension(:), allocatable, save :: rel_vtx	! relevance of the vertex
      logical :: no_intersection_line = .false.			! are the two planes intersecting each other?
      character*120 :: char_i_bl				! convert i_bl to a character to store the data in a file with an appropriate name
      
      if (allocated (v1)) deallocate (v1, v2, oripo_sl, slo_sl)
      allocate(
     &         v1(num_p(i_bl)-1, 3),
     &         v2(num_p(i_bl)-1, 3),
     &	       oripo_sl(num_p(i_bl)**3, 3),
     &	       slo_sl(num_p(i_bl)**3, 3)
     &        )
     
      num_vtx = 1
      num_sl = 1

1000  format(/, t9,'Number of vertices for block limiter plotting = ',
     &       i5, /)

! get the vectors of the parameter form of the plane
      DO i = 1, num_p(i_bl)-1
         CALL normal_to_para(normv_p(i_bl,i,:), v1(i,:), v2(i,:))
      END DO

! calculate the intersection lines of the planes with each other
      DO i = 1, num_p(i_bl)-1
         DO j = i+1, num_p(i_bl)
            b = dot_product(oripo_p(i_bl,j,:), normv_p(i_bl,j,:)) -
     &          dot_product(oripo_p(i_bl,i,:), normv_p(i_bl,j,:))
            c = dot_product(v2(i,:), normv_p(i_bl,j,:))
            d = dot_product(v1(i,:), normv_p(i_bl,j,:))

            IF (abs(b) < limit) THEN				! to small values are set to zero, otherwise ...
               b = 0.						! computing problems with too large/too small numbers ...
            END IF						! occur
            IF (abs(c) < limit) THEN
               c = 0.
            END IF
            IF (abs(d) < limit) THEN
               d = 0.
            END IF

            IF (c .eq. 0. .and. d .ne. 0. .and. b .ne. 0.) THEN
               p1 = oripo_p(i_bl,i,:) + (b/d)*v1(i,:)
               p2 = oripo_p(i_bl,i,:) + (b/d)*v1(i,:) + v2(i,:)
            ELSE IF (c .ne. 0. .and. d .eq. 0. .and. b .ne. 0.) THEN
               p1 = oripo_p(i_bl,i,:) + (b/c)*v2(i,:)
               p2 = oripo_p(i_bl,i,:) + (b/c)*v2(i,:) + v1(i,:)
            ELSE IF (c .ne. 0. .and. d .ne. 0.) THEN
               p1 = oripo_p(i_bl,i,:) + (b/d)*v1(i,:)
               p2 = oripo_p(i_bl,i,:) + (1./d)*(b-c)*v1(i,:) + v2(i,:)
            ELSE
               no_intersection_line = .true.
            END IF

            IF (dot_product((p1-p2),(p1-p2)) == 0.) THEN
               no_intersection_line = .true.
            END IF

            IF (.not.no_intersection_line) THEN
               oripo_sl(num_sl,:) = p1
               u = sqrt(dot_product((p2-p1),(p2-p1)))
               slo_sl(num_sl,:) = (p2-p1)/u
               origin_sl(num_sl,1) = i
               origin_sl(num_sl,2) = j
               num_sl = num_sl+1
            ELSE
               no_intersection_line = .false.
            END IF

         END DO
      END DO
      num_sl = num_sl - 1					! has to be diminished by 1 because starting value was 1

! allocation of the vertex variables
      if (allocated(vtx)) deallocate(vtx)
      allocate(vtx(num_sl**3, 3), stat=alloc_error)
      IF (alloc_error /= 0) THEN
         print*, 'SUBROUTINE gnuplot_data_block_limiter:'
         print*, 'The allocation of >>vtx<< was not successful!'
      END IF

      if (allocated(rel_vtx)) deallocate(rel_vtx)
      allocate(rel_vtx(num_sl**3), stat=alloc_error)
      IF (alloc_error /= 0) THEN
         print*, 'SUBROUTINE gnuplot_data_block_limiter:'
         print*, 'The allocation of >>rel_vtx<< was not successful!'
      END IF
      rel_vtx = .true.						! all vertices are relevant in the beginning

! calculate the intersection points of the straight lines with the planes
      DO i = 1, num_sl
         DO j = 1, num_p(i_bl)
            IF (origin_sl(i,1) /= j .and. origin_sl(i,2) /= j) THEN
               u=dot_product(slo_sl(i,:), normv_p(i_bl,j,:))
               v=dot_product((oripo_sl(i,:)-oripo_p(i_bl,j,:)),
     &                        normv_p(i_bl,j,:))
               IF (u .ne. 0) THEN	
                  lambda = -v/u
                  vtx(num_vtx, :) = oripo_sl(i,:) + lambda*slo_sl(i,:)
                  origin_vtx(num_vtx, 1) = origin_sl(i,1)
                  origin_vtx(num_vtx, 2) = origin_sl(i,2)
                  origin_vtx(num_vtx, 3) = j
                  num_vtx = num_vtx + 1
               END IF
            END IF
         END DO
      END DO
      num_vtx = num_vtx - 1					! has to be diminished by 1 because starting value was 1
      num_imp_vtx = num_vtx

! eliminate multiply appearing vertices
      DO i = 1, num_vtx
         IF (rel_vtx(i)) THEN
            DO j = i+1, num_vtx
               u = dot_product((vtx(i,:)-vtx(j,:)),(vtx(i,:)-vtx(j,:)))
               IF (sqrt(u) < 1.d-7) THEN
                  rel_vtx(j) = .false.
                  num_imp_vtx = num_imp_vtx - 1
               END IF
            END DO
         END IF
      END DO

! ! calculation of the distance between the plane and the origin point of the cs (Hesse normal form) ...
! ! (needed to check if a vertex lies on the surface of the block limiter or not)
!       DO i=1, num_p
!          dist_p_oripocs(i)=dot_product(normv_p(i,:),oripo_p(i,:))
!       END DO

! check whether the intersection points are relevant or not (lie on the surface of the block limiter or not)
      DO i=1, num_vtx
         IF (rel_vtx(i)) THEN
            DO j=1, num_p(i_bl)
               v = dot_product(vtx(i,:), normv_p(i_bl,j,:))
               u = v - dist_p_oripocs(i_bl,j)			! u = distance between plane i and intersection point j
               IF (u > limit) THEN				! is the intersection point on the outside of the cuboid?
                  rel_vtx(i) = .false.
                  num_imp_vtx = num_imp_vtx - 1
                  EXIT
               END IF
            END DO
         END IF
      END DO

! allocation of the vertex variables
      if (allocated(imp_vtx)) deallocate(imp_vtx)
      allocate(imp_vtx(0:(num_imp_vtx-1), 3), stat=alloc_error)
      IF (alloc_error /= 0) THEN
         print*, 'SUBROUTINE gnuplot_data_block_limiter:'
         print*, 'The allocation of >>imp_vtx<< was not successful!'
      END IF

      if (allocated(ori_imp_vtx)) deallocate(ori_imp_vtx)
      allocate(ori_imp_vtx(0:(num_imp_vtx-1), 3), stat=alloc_error)
      IF (alloc_error /= 0) THEN
         print*, 'SUBROUTINE gnuplot_data_block_limiter:'
         print*, 'Allocation of >>ori_imp_vtx<< was not successful!'
      END IF

! write the relevant vertices into a new array imp_vtx
      j = 0
      DO i = 1, num_vtx
         IF (rel_vtx(i)) THEN
            imp_vtx(j,:) = vtx(i,:)
            ori_imp_vtx(j,:) = origin_vtx(i,:)
            j = j+1
         END IF
      END DO

! write the data into a file so that gnuplot can plot it appropriately
      IF (i_bl == 1) THEN
         open(iun_2, file='bl_plot.dat')
      ELSE
         write(char_i_bl, *) i_bl
         char_i_bl = adjustl(char_i_bl)
         open(iun_2, file='bl_plot_'//trim(char_i_bl)//'.dat')
      END IF

      DO i_1 = 1, num_p(i_bl)						! chose one plane ...

         first_vtx: DO j = 0, num_imp_vtx-1			! search all the indices ...
            DO e = 1, 3						! to find a vertex ...
               IF (ori_imp_vtx(j, e) == i_1) THEN		! who has the above chosen plane as an origin plane ...
                  initial_vtx_index = j
                  write(iun_2, *) imp_vtx(j, :)			! write the vertex to a file ...
                  f = mod(e+1, 3)+1				! get another index of an origin plane of chosen vertex
                  i_2 = ori_imp_vtx(j, f)
                  EXIT first_vtx
               END IF
            END DO
                    END DO first_vtx

         DO e = 1, num_imp_vtx + 1				! could have been an only 'DO' loop, but it is better ...
      								! to limit the maximum times of execution
            CALL find_adjacent_vertex(j, i_1, i_2, ori_imp_vtx,
     &                                num_imp_vtx, x)
            write(iun_2, *) imp_vtx(x,:)
            j = x
            IF (j == initial_vtx_index) THEN
               write(iun_2, *)
               write(iun_2, *)
               EXIT
            END IF
         END DO

      END DO

      close(iun_2)

      write(*, fmt=1000) num_imp_vtx

      IF ((form_type == 1 .and. num_imp_vtx /= 8) .or.
     &    (form_type == 2 .and. num_imp_vtx /= 8) .or.
     &    (form_type == 3 .and. num_imp_vtx /= 10)) THEN
         print*, 'SUBROUTINE gnuplot_data_block_limiter:'
         print*, 'The number of block limiter vertices is wrong!'
         print*, 'Please check >>block_limiter.dat<<!'
         STOP
      END IF

      RETURN
      END SUBROUTINE gnuplot_data_block_limiter 
!--------------------------------------------------------------------------

!----------------------------------------------------------------
! SUBROUTINE TO GET TWO VECTORS FOR THE PARAMETER FORM OF A PLANE
! OUT OF A NORMAL VECTOR
!----------------------------------------------------------------
      SUBROUTINE normal_to_para(normv, v1, v2)
      IMPLICIT NONE

      real*8, dimension(0:2), intent(in) :: normv		! normal vector of the plane
      real*8, dimension(0:2), intent(out) :: v1, v2		! vectors of the parameter form of the plane
      real*8, parameter :: limit = 10.d-10		! value considered as zero 
      integer :: i

! get the vectors of the parameter form of the plane
      i = maxloc(abs(normv), dim = 1) - 1
      v1(mod(i+1, 3)) = normv(i)
      v1(i) = -normv(mod(i+1,3))
      v1(mod(i+2, 3)) = 0
      CALL cross_product(normv, v1, v2)

! error output
      IF (abs(dot_product(v1, v2)) > limit) THEN
         print*, 'SUBROUTINE normal_to_para'
         print*, 'The two vectors of the parameter form '
         print*, 'are not perpendicular!'
      END IF

      RETURN
      END SUBROUTINE normal_to_para
!--------------------------------------------------------------------------

!----------------------------------------------------------
! FUNCTION WHICH FINDS A VERTEX ADJACENT TO A GIVEN VERTEX 
! [USED FOR PRODUCING THE PLOT DATA FOR GNUPLOT]
!----------------------------------------------------------
      SUBROUTINE find_adjacent_vertex(x_old, i_1, i_2, 
     &                                ori_vtx, num_vtx, x_new)
      IMPLICIT NONE

      integer, intent(inout) :: x_old			! index of the origin vertex (to which an adjacent vertex ...
      							! should be found)
      integer, intent(in) :: i_1			! indices of the origin plane of the vertex to which the ...
      							! traverse belongs
      integer, intent(inout) :: i_2			! first: second index of the former vertex, ...
                                                        ! then: index of an origin plane of the adjacent vertex ...
      							! which should be used in the next iteration of the subroutine
      integer, dimension(0:(num_vtx-1), 3), intent(in) :: ori_vtx	! origin planes of the vertices
      integer, intent(in) :: num_vtx			! number of vertices
      integer, intent(out) :: x_new			! index of the new adjacent vertex
      integer :: m, n, y				! counting variables

! find the index of an adjacent vertex
      vtx_loop: DO y = 1, num_vtx
         DO m = 1, 3
            x_new = mod((y + x_old), num_vtx)		! special counting variable so that the adjacent vertex ...
            IF (ori_vtx(x_new, m) == i_1) THEN		! will not be the initial point
               DO n = 1, 3
                  IF (ori_vtx(x_new, n) == i_2) THEN
                      EXIT vtx_loop
                  END IF
               END DO
            END IF
         END DO
                END DO vtx_loop

! find the new index i_2 which should be used for the next iteration
      DO m = 1, 3
         IF (ori_vtx(x_new,m)/=i_1 .and. ori_vtx(x_new, m)/=i_2) THEN
            i_2 = ori_vtx(x_new, m)
            EXIT
         END IF
      END DO

      RETURN
      END SUBROUTINE find_adjacent_vertex
!----------------------------------------------------------

!----------------------------------------------------------------
! SUBROUTINE TO CALCULATE THE CROSS PRODUCT OF TWO VECTORS
!----------------------------------------------------------------
      SUBROUTINE cross_product(A, B, C)

      IMPLICIT NONE						! no default typing

      REAL*8, DIMENSION(3), INTENT (IN)    :: A			! multiplicand 3-vector
      REAL*8, DIMENSION(3), INTENT (IN)    :: B			! multiplier 3-vector
      REAL*8, DIMENSION(3), INTENT (OUT)   :: C			! result: 3-vector cross product

! compute cross product components
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)

      RETURN
      END SUBROUTINE cross_product
!----------------------------------------------------------------

      END MODULE block_limiter





      subroutine check_point_bl
      use block_limiter
      implicit none

      real*8, parameter :: l0 = 10.d-10
      real*8  :: x(3), phi0, R, Z, R1, R2, Z1, Z2, cphi, sphi
      integer :: i, j, n, m
      logical :: lcheck

      m = 100
      n = 100
      R1 = 165.d0
      R2 = 185.d0
      Z1 = -55.d0
      Z2 = -45.d0
      phi0 = 45.d0 / 180.d0 * pi
      cphi = cos(phi0)
      sphi = sin(phi0)

      open  (99, file='test_bl.txt')
      do j=0,m
         Z    = Z1 + j * (Z2-Z1) / m
         x(3) = Z
         do i=0,n
            R = R1 + i * (R2-R1) / n
            x(1) = R * cphi
            x(2) = R * sphi
            call check_point_behind_bl(x, l0, 1, lcheck)
            if (.not.lcheck) write (99, *) R, Z
            !write (99, *) R, Z, lcheck
         enddo
      enddo
      close (99)

      !SUBROUTINE check_point_behind_bl(x, limit, i_bl, result)
      return
      end subroutine check_point_bl

      subroutine test_bl_triangles_eirene
      use block_limiter
      !call bl_triangles_eirene (1)
      call make_triangles_sphereX (10, 10)
      end subroutine test_bl_triangles_eirene
