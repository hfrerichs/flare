!===============================================================================
! Precalculate distance to surface on grid nodes
!
! Input (taken from run control file):
!    Grid_File          Reference surface (t_flux_surface_3D)
!    R_start, R_end
!    Z_start, Z_end     Definition of computational box
!    N_R, N_Z           Resolution for computational box
!===============================================================================
subroutine setup_distance_to_surface
  use iso_fortran_env
  use run_control, only: R_start, R_end, Z_start, Z_end, N_R, N_Z, N_phi, N_sym, &
                         Grid_File, Output_File
  use flux_surface_3D
  use interpolate3D
  use parallel
  implicit none

  type(t_interpolate3D)   :: distance
  type(t_flux_surface_3D) :: S
  real(real64) :: Rc, Zc, w, h, r(2)
  integer      :: i, j, k


  ! initialize
  if (firstP) then
     write (6, *) 'Calculating distance to surface "', trim(Grid_File), '"'
  else
     return
  endif


  ! load (flux) surface
  call S%load(filename=Grid_File)


  ! check user input
  if (N_phi .ne. S%n_phi) then
     write (6, *) 'error: toroidal resolution does not match!', N_phi, S%n_phi
     stop
  endif
  if (N_sym .ne. S%n_sym) then
     write (6, *) 'error: toroidal symmetry does not match!', N_sym, S%n_sym
     stop
  endif


  ! initialize interpolation scheme
  Rc = 0.5d0 * (R_start + R_end)
  Zc = 0.5d0 * (Z_start + Z_end)
  w  = R_end - R_start
  h  = Z_end - Z_start
  call distance%new(N_R+1, N_Z+1, Rc, Zc, w, h, N_phi+1, N_sym)


  ! calculate distances from grid nodes to flux surface
  write (6, 1000) N_phi
  do k=1,N_phi
     write (6, 1001) k-1
     do j=1,N_Z+1
        do i=1,N_R+1
           r(1) = distance%R(i)
           r(2) = distance%Z(j)

           distance%d(i,j,k) = S%slice(k-1)%get_distance_to(r)
        enddo
     enddo
  enddo
 1000 format (3x,'- Number of slices: ', i4)
 1001 format (8x,i4)
  distance%d(:,:,N_phi+1) = distance%d(:,:,1)


  ! save output
  call distance%store(filename='distance.dat')

end subroutine setup_distance_to_surface
!===============================================================================





!===============================================================================
!===============================================================================
subroutine evaluate_distance_to_surface
  use iso_fortran_env
  use run_control, only: Grid_File, Output_File
  use usr_grid
  use parallel
  use interpolate3D
  use math
  implicit none

  integer, parameter    :: iu = 32

  type(t_interpolate3D) :: distance
  real(real64) :: y(3)
  integer      :: iflag


  ! initialize
  if (firstP) then
     write (6, *) 'Evaluating distance to flux surface "', trim(Grid_File), '"'
  endif


  ! load grid
  call read_grid (Grid_File, use_coordinates=COORDINATES(CYLINDRICAL))


  ! load precalculated distances
  call distance%load(filename='distance.dat')


  ! evaluate on grid nodes
  open  (iu, file=Output_File)
  grid_loop: do
     call get_next_grid_point (iflag, y)
     if (iflag .ne. 0) exit grid_loop

     write (iu, *) distance%eval(y)
  enddo grid_loop
  close (iu)

end subroutine evaluate_distance_to_surface
