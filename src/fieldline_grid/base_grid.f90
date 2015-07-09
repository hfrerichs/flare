!===============================================================================
! Generate base grid for selected zone
!===============================================================================
module base_grid
  use iso_fortran_env
  use fieldline_grid
  use flux_surface_2D
  use divertor
  use grid
!  use curve2D
  implicit none

!  type, public :: t_zone_base
!     type(t_curve)   :: radial_boundary(2), poloidal_boundary(2)
!     type(t_spacing) :: radial_spacing, poloidal_spacing
!
!     type(t_grid)    :: G
!  end type t_zone_base


  type :: t_interface_radial
     ! id of connected zones/layers, radial node index
     integer :: iz(2), ir(2)

     ! poloidal and toroidal range of interface (ip0->ip0+np, it0->it0+nt)
     integer :: ip0(2), np, it0(2), nt

     ! interface type (0: strike point, >0: X-point)
     integer :: itype(2)

     ! real space geometry of interface
     type(t_flux_surface_2D), pointer :: F
  end type t_interface_radial
  type(t_interface_radial), dimension(:), allocatable :: interface_radial

  ! numer of interfaces between zones
  integer :: interfaces


  ! connection of separatrix from X-points
  integer, dimension(:), allocatable :: connectX

  ! number of X-points
  integer :: nX


  procedure(), pointer :: setup_topology

  contains
  !=====================================================================



  !=====================================================================
  ! allocate base grid with nr x np cells at toroidal position phi
  !=====================================================================
  subroutine new_base_grid(B, nr, np, phi)
  type(t_grid), intent(inout) :: B
  integer,      intent(in)    :: nr, np
  real(real64), intent(in)    :: phi

  call B%new(CYLINDRICAL, MESH_2D, FIXED_COORD3, nr+1, np+1, fixed_coord_value=phi)

  end subroutine new_base_grid
  !=====================================================================



  !=====================================================================
  subroutine make_base_grids_auto
  use grid
  use mesh_spacing

  ! base grids
  type(t_grid), dimension(:,:), allocatable ::G !(0:blocks-1,0:layers-1)

  ! mesh pointer
  real(real64), dimension(:,:,:), pointer   :: M1, M2

  type(t_spacing) :: Sp
  real(real64) :: phi
  integer      :: iblock, iz0, iz, il, ix1, ix2, i


  write (6, 1000)

  !.....................................................................
  ! 1. check input
  if (n_interpolate < 0) then
     write (6, *) 'error: n_interpolate must not be negative!'; stop
  endif
! nr needs to be checked for each zone
!  if (n_interpolate > nr(0)-2) then
!     write (6, *) 'error: n_interpolate > nr0 - 2!'; stop
!  endif
  !.....................................................................


  !.....................................................................
  ! 2. initialize geometry
  call setup_geometry(nX, connectX)
  !tmpcall setup_domain()
  !.....................................................................


  !.....................................................................
  ! 3. setup working arrays for base grid
  allocate (G(0:blocks-1, 0:layers-1))
  !.....................................................................


  !.....................................................................
  ! 4. main loop
  do iblock=0,blocks-1
     write (6, 1001) iblock

     ! set index for 1st zone in block
     iz0 = iblock*layers

     ! set local variables for resolution (nr, np, npL, npR)
     call load_local_resolution(iblock)


     ! setup cell spacings (radial direction)
     do il=0,layers-1
        iz = iz0 + il
        call Zone(iz)%Sr%init(radial_spacing(il))
     enddo
     !TODO: setup cell spacings (poloidal direction)


     ! initialize base grids in present block
     phi = Block(iblock)%phi_base / 180.d0 * pi
     do il=0,layers-1
        call new_base_grid(G(iblock,il), nr(il), np(il), phi)
     enddo

     ! set up mesh pointers (or not)
     ! ...


     ! start grid generation
     ! 1. make interfaces between zones
     write (6, 2001)
     do i=0,interfaces-1
        M1 => G(iblock, interface_radial(i)%iz(1))%mesh
        M2 => G(iblock, interface_radial(i)%iz(2))%mesh
        call make_interface(interface_radial(i), Sp, M1, M2)
     enddo
     write (6, *)


     ! 2. generate radial/poloidal discretization
     ! 2.1 confined region / high pressure region / poloidally periodic region
     write (6, 2002)
     write (6, 2003)
     iz = iz0
     ! start at X-point 1
     ix1 = 1
     ix2 = connectX(ix1)
     do
        ! connect directly to another X-point
        if (ix2 > 0  .and.  ix1 > 0) then
           write (6, 2004) abs(ix1), abs(ix2)

        ! use X-point ix1 as guidance
        elseif (ix1 < 0) then
           write (6, 2005) abs(ix1), abs(ix2)

        ! use X-point ix2 as guidance
        else
           write (6, 2006) abs(ix1), abs(ix2)
        endif



        ! do stuff...
        !write (6, *) 'do stuff ', ix1, ix2


        ! are we back at 1st X-point?
        if (ix2 == 1) exit
        ! otherwise go to next segment
        ix1 = ix2
        if (ix2 > 0) then
           ix2 = connectX(abs(ix1))
        else
           ix2 = 1
        endif
     enddo
     write (6, *)


     ! 2.2 scrape-off layer and private flux region
     do il=1,layers-1
        iz = iz0 + il
        select case(Zone(iz)%itypeR)
        case(TYPE_SOL)
        case(TYPE_SOLMAP)
        case(TYPE_PFR)
        end select
     enddo


     ! 99. output
     do i=0,layers-1
        iz = iz0 + i
        call write_base_grid(G(iblock,i), iz)
     enddo
     write (6, 1002) iblock
  enddo


 1000 format(//3x,'- Setup for base grids:')
 1001 format(//1x,'Start generation of base grids for block ',i0,' ',32('.'))
 1002 format(1x,'Finished generation of base grids for block ',i0,' ',32('.'),//)
 2001 format(3x,'- Generating interfaces between zones')
 2002 format(3x,'- Generating discretization of high pressure region')
 2003 format(8x,'connecting X-points')
 2004 format(10x,i0,' -> ',i0)
 2005 format(10x,'(',i0,') -> ',i0)
 2006 format(10x,i0,' -> (',i0,')')
  end subroutine make_base_grids_auto
  !=====================================================================



  !=====================================================================
  subroutine make_interface(I, Sp, M1, M2)
  use equilibrium
  type(t_interface_radial), intent(in) :: I
  type(t_spacing),          intent(in) :: Sp
  real(real64), dimension(:,:,:), pointer, intent(in) :: M1, M2

  type(t_flux_surface_2D), pointer :: F
  real(real64) :: D(0:I%np,2)
  integer      :: ir(2), ip0(2), ix1, ix2, np


  ! set up indices and resolution parameters
  ir  = I%ir
  ip0 = I%ip0
  ix1 = abs(I%itype(1));  ix2 = abs(I%itype(2))
  np  = I%np
  F   => I%F


  ! 1. connection between two X-points
  if (ix1 > 0  .and.  ix2 > 0) then
     call make_interface_core(F, ix1, ix2, Sp, np, D)

  ! 2. connect from strike point to X-point
  elseif (ix1 == 0  .and.  ix2 > 0) then
     call divertor_leg_discretization(F%t_curve, 1.d0 - etaR(ix2), np, D)

  ! 3. connect from X-point to strike point
  elseif (ix1 > 0  .and.  ix2 == 0) then
     call divertor_leg_discretization(F%t_curve, etaL(ix1), np, D)

  ! 4. connect strike point to another strike point
  else
     write (6, *) 'error in subroutine make_interface:'
     write (6, *) 'strike point to strike point connection not implemented!'
     stop
  endif


  M1(ir(1), ip0(1):ip0(1)+np, :) = D
  M2(ir(2), ip0(2):ip0(1)+np, :) = D

  end subroutine make_interface
  !=====================================================================



  !=====================================================================
  subroutine make_HPR_grid(G)
  type(t_grid), intent(inout) :: G

  real(real64), dimension(:,:,:), pointer :: M


  M => G%mesh
  select case(discretization_method)
  ! 1. use fixed poloidal angle for poloidal discretization
  case (POLOIDAL_ANGLE)
     if (Dtheta_separatrix > 0.d0) then
        write (6, *) 'error: Dtheta_separatrix > 0 not implemented for discretization by poloidal angle!'
        stop
     endif
     ! 1.1 unperturbed flux surfaces
!     call make_flux_surfaces_HPR(M_HPR, nr(0), np(0), 2+n_interpolate, nr(0)-1, rpath(0), Zone(iz0)%Sr, Sp_HPR)

     ! 1.2 interpolated surfaces
!        call make_interpolated_surfaces(M_HPR, nr(0), np(0), 1, 2+n_interpolate, Zone(iz0)%Sr, Sp_HPR, C_in(iblock,:))



  ! 2. use gradient(Psi) descent for poloidal discretization
  case (ORTHOGONAL)
     ! 2.1 unperturbed flux surfaces (nodes at fixed PsiN)
!     call make_ortho_grid(M_HPR, nr(0), np(0), nr(0), 2+n_interpolate, nr(0)-1, &
!                          0, 1, npR(0)-1, npR(0), 2, rpath(0), Zone(iz0)%Sr, Sp_HPR)
!     call make_ortho_grid(M_HPR, nr(0), np(0), nr(0), 2+n_interpolate, nr(0)-1, &
!                          np(0), npR(0)+1, np(0)-1, npR(0), 2, rpath(0), Zone(iz0)%Sr, Sp_HPR)

     ! 2.2 interpolated surfaces (along grad PsiN path)
!     call make_interpolated_surfaces_ortho(M_HPR, nr(0), np(0), 2+n_interpolate, Zone(iz0)%Sr, Sp_HPR, &
!                                           C_in(iblock,:), DPsiN1(iblock,1))

  end select

  end subroutine make_HPR_grid
  !=====================================================================

end module base_grid
