!===============================================================================
! Simply connected grid layout
!===============================================================================
module modtopo_sc
  use iso_fortran_env
  use fieldline_grid, unused => TOPO_SC
  use inner_boundary
  implicit none
  private


  character(len=*), parameter :: &
     TEST_METHOD_LENGTH = 'test_method_length'


  public :: &
     setup_topo_sc, &
     make_base_grids_sc, &
     post_process_grid_sc

  contains
  !=====================================================================


  
  !=====================================================================
  subroutine setup_topo_sc()
  use emc3_grid, only: NZONET
  implicit none

  integer :: iz, ib


  ! 0. setup number of zones
  NZONET = blocks


  write (6, 1000)
  do ib=0,blocks-1
     ! 1. set up derived parameters
     iz = ib
     Zone(iz)%nr = Block(ib)%nr(0) + nr_EIRENE_core + nr_EIRENE_vac(0)
     Zone(iz)%np = Block(ib)%np(0)


     ! 2. set up zones
     call Zone(iz)%setup(ib, 0, TYPE_SINGLE_LAYER, SF_PERIODIC)


     ! 3. show zone information
     write (6, 1002) iz, Zone(iz)%nr, Zone(iz)%np, Zone(iz)%nt
  enddo
 1000 format(8x,'Grid resolution is (radial x poloidal x toroidal):')
 1002 format(12x,i3,3x,'(',i0,' x ',i0,' x ',i0,')')

  end subroutine setup_topo_sc
  !=====================================================================



  !=====================================================================
  subroutine setup_domain()
  use equilibrium
  use divertor, only: Pmag

  real(real64) :: x1(2), x2(2), d(2), tmp(3), theta0


  tmp    = get_magnetic_axis(0.d0); Pmag = tmp(1:2)
  x1     = x_in2(1:2)
  theta0 = get_poloidal_angle(x_in2)
  select case(discretization_method)
  case(POLOIDAL_ANGLE)
     call load_inner_boundaries(theta0)
     x2    = x1
     x2(1) = x2(1) - d_SOL(1)
     d     = x1-x2
     call rpath(0)%setup_linear(x2, d)

  case(ORTHOGONAL)
     call load_inner_boundaries(theta0)
     call rpath(0)%generate(x1, ASCENT_LEFT, LIMIT_LENGTH, d_SOL(1))
     call rpath(0)%flip()

  case(TEST_METHOD_LENGTH)
     call load_inner_boundaries()

  end select
  call rpath(0)%plot(filename='rpath_0.plt')

  end subroutine setup_domain
  !=====================================================================


  !=====================================================================
  subroutine make_base_grids_sc()
  use grid
  use divertor, only: make_flux_surfaces_HPR, make_interpolated_surfaces, &
                      make_ortho_grid

  real(real64), dimension(:,:,:), pointer :: M
  type(t_grid) :: G(0:blocks-1)
  real(real64) :: phi
  integer      :: iblock, iz


  call setup_domain()

  do iblock=0,blocks-1
     write (6, *) iblock

     ! set zone indix
     iz  = iblock

     ! set local variables for resolution
     call load_local_resolution(iblock)


     ! cell spacings
     call Zone(iz)%Sp%init(poloidal_spacing(0))
     call Zone(iz)%Sr%init(radial_spacing(0))


     ! initialize base grid in present block
     phi = Block(iblock)%phi_base / 180.d0 * pi
     call G(iblock)%new(CYLINDRICAL, MESH_2D, 3, nr(0)+1, np(0)+1, fixed_coord_value=phi)
     M => G(iblock)%mesh


     select case(discretization_method)
     case (POLOIDAL_ANGLE)
        call make_flux_surfaces_HPR(M, nr(0), np(0), 2+n_interpolate, nr(0), rpath(0), Zone(iz)%Sr, Zone(iz)%Sp)
        call make_interpolated_surfaces(M, nr(0), np(0), 1, 2+n_interpolate, Zone(iz)%Sr, Zone(iz)%Sp, C_in(iblock,:))

     case (ORTHOGONAL)
        call setup_inner_boundaries(G(iblock), iblock, 0, Zone(iz)%Sp)
!        call make_ortho_grid(M, nr(0), np(0), nr(0), 2+n_interpolate, nr(0), &
!                             -1, 0, np(0), -1, 0, rpath(0), Zone(iz)%Sr, Zone(iz)%Sp)
        call make_ortho_grid(M, nr(0), np(0), 1, 2, nr(0), &
                             0, 1, np(0)-1, -1, 0, rpath(0), Zone(iz)%Sr, Zone(iz)%Sp, periodic=.true.)
     case (TEST_METHOD_LENGTH)
        call make_grid_TEST_METHOD_LENGTH(M, nr(0), np(0))
     end select

     ! output
     call write_base_grid(G(iblock), iblock)
  enddo

  end subroutine make_base_grids_sc
  !=====================================================================



  !=====================================================================
  subroutine make_grid_TEST_METHOD_LENGTH(M, nr, np)
  use curve2D
  use inner_boundary
  use grid
  use math
  use dataset
  use run_control, only: Debug

  real(real64), dimension(:,:,:), pointer, intent(inout) :: M
  integer, intent(in) :: nr, np


  type(t_dataset) :: w
  type(t_grid)  :: G(0:blocks-1)
  type(t_curve) :: C, C0
  character(len=72) :: filename
  real(real64) :: eta, xi, phi, dr, x(2), v1(2), v2(2), cosa, x0(2), x1(2), x2(2)
  real(real64) :: et(2), en(2), xh(2), th, sh, xi1
  integer :: iblock, iz, i, j, n, j3(-1:1), ish
  logical :: l

     !call C_in(iblock,1)%plot(filename='fsin1_0.plt')

     call C_in(iblock,0)%setup_length_sampling_curvature_weighted()
     !call C_in(iblock,1)%setup_length_sampling_curvature_weighted()




     call w%new(np+1, 2, -1)
     call C0%new(np)
     C0%closed = .true.
     do j=0,np
        xi = Zone(iz)%Sp%node(j,np)
        call C_in(iblock,0)%sample_at(xi, x, et)
        M(0,j,:) = x

        !call C_in(iblock,1)%sample_at(xi, x, et)
        !call C0%curvature('kappa1.plt')
        !G(iblock)%mesh(1,j,:) = x
        !C0%x(j,:) = x

        en(1) =  et(2)
        en(2) = -et(1)
        l = intersect_curve(x, x+en, C_in(iblock,1), xh, th, sh, ish, 1)
        if (l) then
           M(1,j,:) = xh
           C0%x(j,:) = xh
           !xi1 = C_in(iblock,0)%w(ish-1) + sh*(C_in(iblock,0)%w(ish) - C_in(iblock,0)%w(ish-1))
           !write (99, *) xi, th, xh, xi1
           !w%x(j,1) = xi
           !w%x(j,2) = xi1

           !write (98, *) x
           !write (98, *) xh
           !write (98, *)
        else
           write (6, *) 'error: no intersection found for node ', j, '!'
        endif
     enddo

!     w%x(0,:) = 0.d0
!     call w%sort_rows(2)
!     call w%plot(filename='w.plt')
!     do j=0,np
!        xi = w%x(j,2)
!        call C_in(iblock,0)%sample_at(xi, x)
!        G(iblock)%mesh(0,j,:) = x
!     enddo


     do i=2,nr
        eta = Zone(iz)%Sr%node(i-1,nr-1)
        dr  = eta * D_SOL(1)

        !call C%copy(C_in(iblock, 1))
        call C%copy(C0)
        call C%left_hand_shift(-dr)
!        write (filename, 1000) i
! 1000   format('expanded_sf_',i0,'.plt')
!        call C%plot(filename=filename)

        if (np .ne. C%n_seg) then
           write (6, *) 'error: nodes were dropped in subroutine left_hand_shift!'
           stop
        endif
        do j=0,np
           M(i,j,:) = C%x(j,:)
        enddo
     enddo

  end subroutine make_grid_TEST_METHOD_LENGTH
  !=============================================================================



  !=============================================================================
  subroutine post_process_grid_sc()

  return
  end subroutine post_process_grid_sc
  !=============================================================================



end module modtopo_sc
