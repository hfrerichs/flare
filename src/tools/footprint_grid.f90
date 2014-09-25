!===============================================================================
subroutine footprint_grid
  use iso_fortran_env
  use run_control, only: Grid_File, Output_File, n_theta, n_phi, offset
  use parallel
  implicit none

  integer, parameter :: iu = 32, iu2 = 33


  if (firstP) then
     write (6, *) 'Generate footprint grid, output in: ', adjustl(trim(Output_File)), &
                  ', ', adjustl(trim(Grid_File))
     write (6, *)
  else
     return
  endif

  call footprint_grid_Q(1)
  contains
!=======================================================================
  subroutine footprint_grid_Q (iele)
  use boundary
  use mesh_spacing
  use math
  use equilibrium
  use quad_ele
  use grid
  use curve2D
  integer, intent(in) :: iele

  type(t_grid)     :: G_sample, G_plot
  type(t_quad_ele) :: S
  type(t_curve)    :: C
  real(real64) :: tau, xi, phi, theta, r(3), d(2)
  integer      :: i, j, n, ig


  if (iele < 1  .or.  iele > n_quad) then
     write (6, *) 'error: surface element ', iele, ' does not exist!'
     stop
  endif


  ! copy surface (to be modified locally)
  call S%new(S_quad(iele)%n_phi, S_quad(iele)%n_RZ, S_quad(iele)%n_sym)
  S%phi = S_quad(iele)%phi
  S%R   = S_quad(iele)%R
  S%Z   = S_quad(iele)%Z


  ! shrink surface (because starting points for field line tracing must not lay on the boundary surface)
!  do i=0,S%n_phi
!     r(3) = S%phi(i)
!     r    = get_magnetic_axis(r(3))
!     do j=0,S%n_RZ
!        d(1)     = S%R(i,j) - r(1)
!        d(2)     = S%Z(i,j) - r(2)
!
!        S%R(i,j) = r(1) + (1.d0 - 1.d-3) * d(1)
!        S%Z(i,j) = r(2) + (1.d0 - 1.d-3) * d(2)
!     enddo
!  enddo
  !call S%left_hand_shift(1.d-3)


  n = (n_phi+1) * (n_theta+1)
  ! output file for cylindrical coordinates
  !open  (iu, file=Grid_File)
  !write (iu, 1000)
  !write (iu, 1001) (n_phi+1) * (n_theta+1)
  call G_sample%new(CYLINDRICAL, UNSTRUCTURED, 0, n)

  ! output file for local coordinates
  !open  (iu2, file=Output_File)
  !write (iu2, 2000)
  !write (iu2, 2001) n_phi+1
  !write (iu2, 2002) n_theta+1
  call G_plot%new(LOCAL, UNSTRUCTURED, 3, n)

  call S%setup_coefficients()
  ig = 0
  do i=0,n_phi
     tau = Equidistant%node(i, n_phi)
     phi = S_quad(iele)%sample_phi(tau)
     C   = S_quad(iele)%slice(phi)

     call C%left_hand_shift(offset)
     call C%setup_length_sampling()

     do j=0,n_theta
        ig = ig + 1
        xi = Equidistant%node(j, n_theta)

        r(3) = phi
        call C%sample_at(xi, r(1:2))
!        r  = S%sample(tau, xi)
        G_sample%x(ig,:) = r
        !write (iu, *) r

!        phi   = r(3) / pi * 180.d0
        theta = get_poloidal_angle(r) / pi * 180.d0
        if (theta < 0) theta = theta + 360.d0
        !write (iu2, *) phi, theta
        G_plot%x(ig,1) = phi / pi * 180.d0
        G_plot%x(ig,2) = theta
     enddo

     call C%destroy()
  enddo
  !close (iu)
  !close (iu2)

  call G_sample%store(filename=Grid_File)
  call G_plot%store(filename=Output_File, header='local coordinates: Phi[deg], Theta[deg]')

! 1000 format ('# grid_id = 9       (cylindrical coordinates: R[cm], Z[cm], Phi[rad])')
! 1001 format ('# grid resolution:   n_RZphi =  ',i10)
! 2000 format ('# grid_id = 20      (local coordinates: Phi[deg], Theta[deg]')
! 2001 format ('# Toroidal angle     n_phi   =  ',i10)
! 2002 format ('# Poloidal angle     n_theta =  ',i10)
  end subroutine footprint_grid_Q
!=======================================================================
end subroutine footprint_grid
