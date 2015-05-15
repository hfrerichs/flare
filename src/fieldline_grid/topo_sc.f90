!===============================================================================
! Simply connected grid layout
!===============================================================================
module topo_sc
  use iso_fortran_env
  use fieldline_grid, unused => TOPO_SC
  implicit none
  private

  public :: &
     setup_topo_sc, &
     make_base_grids_sc

  contains
  !=====================================================================


  
  !=====================================================================
  subroutine setup_topo_sc()
  use emc3_grid, only: NZONET
  implicit none

  integer :: iz, ib


  ! 0. setup number of zones
  NZONET = blocks


  ! 1. setup resolution for each zone
  write (6, 1000)
  do ib=0,blocks-1
     iz = ib
     if (Zone(iz)%nr == -1) Zone(iz)%nr = nr(0)
     if (Zone(iz)%np == -1) Zone(iz)%np = np(0)
     !if (Zone(iz)%nt == -1) Zone(iz)%nt = nt
     Zone(iz)%nt = Block(ib)%nt

     ! setup toroidal discretization
     allocate (Zone(iz)%phi(0:Zone(iz)%nt))
     Zone(iz)%phi     = Block(ib)%phi
     Zone(iz)%it_base = Block(ib)%it_base

     write (6, 1002) iz, Zone(iz)%nr, Zone(iz)%np, Zone(iz)%nt
  enddo
 1000 format(8x,'Grid resolution is (radial x poloidal x toroidal):')
 1002 format(12x,i3,3x,'(',i0,' x ',i0,' x ',i0,')')


  ! 2. setup connectivity between zones
  do iz=0,blocks-1
     Zone(iz)%isfr(1) = SF_CORE
     Zone(iz)%isfr(2) = SF_VACUUM
     Zone(iz)%isfp(1) = SF_PERIODIC
     Zone(iz)%isfp(2) = SF_PERIODIC
     Zone(iz)%isft(1) = SF_MAPPING
     Zone(iz)%isft(2) = SF_MAPPING

     Zone(iz)%r_surf_pl_trans_range(1) = nr_EIRENE_core
     Zone(iz)%r_surf_pl_trans_range(2) = Zone(iz)%nr - nr_EIRENE_vac
     Zone(iz)%p_surf_pl_trans_range(1) = 0
     Zone(iz)%p_surf_pl_trans_range(2) = Zone(iz)%np
  enddo

  end subroutine setup_topo_sc
  !=====================================================================



  !=====================================================================
  subroutine setup_domain()
  use inner_boundary

  ! discretization by length
  call load_inner_boundaries()
  ! discretization by poloidal angle
  !call load_inner_boundaries(0.d0)

  end subroutine setup_domain
  !=====================================================================


  !=====================================================================
  subroutine make_base_grids_sc()
  use curve2D
  use inner_boundary
  use grid
  use math
  use dataset
  use run_control, only: Debug

  type(t_dataset) :: w
  type(t_grid)  :: G(0:blocks-1)
  type(t_curve) :: C, C0
  character(len=72) :: filename
  real(real64) :: eta, xi, phi, dr, x(2), v1(2), v2(2), cosa, x0(2), x1(2), x2(2)
  real(real64) :: et(2), en(2), xh(2), th, sh, xi1
  integer :: iblock, iz, i, j, n, nr, np, j3(-1:1), ish
  logical :: l


  call setup_domain()

  do iblock=0,blocks-1
  !do iblock=0,0
     write (6, *) iblock
     iz  = iblock
     phi = Block(iblock)%phi_base

     nr = Block(iblock)%nr(0)
     np = Block(iblock)%np(0)
     call G(iblock)%new(CYLINDRICAL, MESH_2D, 3, nr+1, np+1, fixed_coord_value=phi)

     !call C_in(iblock,1)%plot(filename='fsin1_0.plt')
     call C_in(iblock,0)%setup_length_sampling_curvature_weighted()
     !call C_in(iblock,1)%setup_length_sampling_curvature_weighted()


     ! cell spacings
!     call Zone(iz)%Sp%init(poloidal_spacing(0))
     call Zone(iz)%Sr%init(radial_spacing(0))



     call w%new(np+1, 2, -1)
     call C0%new(np)
     C0%closed = .true.
     do j=0,np
        xi = Zone(iz)%Sp%node(j,np)
        call C_in(iblock,0)%sample_at(xi, x, et)
        G(iblock)%mesh(0,j,:) = x

        !call C_in(iblock,1)%sample_at(xi, x, et)
        !call C0%curvature('kappa1.plt')
        !G(iblock)%mesh(1,j,:) = x
        !C0%x(j,:) = x

        en(1) =  et(2)
        en(2) = -et(1)
        l = intersect_curve(x, x+en, C_in(iblock,1), xh, th, sh, ish, 1)
        if (l) then
           G(iblock)%mesh(1,j,:) = xh
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
           G(iblock)%mesh(i,j,:) = C%x(j,:)
        enddo
     enddo
     !call G(iblock)%plot_mesh(filename='test.plt')
     call write_base_grid(G(iblock), iblock)
  enddo

  end subroutine make_base_grids_sc
  !=====================================================================

  !.....................................................................
  subroutine write_base_grid(G, iz)
  use grid
  type(t_grid), intent(in) :: G
  integer, intent(in)      :: iz

  character(len=72) :: filename


  ! write grid file for field line tracing
  write (filename, 9000) iz
  call G%store(filename=filename)

  ! write grid file for plotting
  write (filename, 9001) iz
  call G%plot_mesh(filename)

 9000 format ('base_grid_',i0,'.dat')
 9001 format ('base_grid_',i0,'.plt')
  end subroutine write_base_grid
  !.....................................................................



end module topo_sc
