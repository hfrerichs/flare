!===============================================================================
! Interface to magnetic field data
!
! The magnetic field strength is given in T (tesla). Spatial coordinates are
! in cm (centi meter), the azimuth angle is given in rad (radians).
!===============================================================================
module bfield
  use run_control
  use parallel
  use reconstruct
  use geqdsk
  use polygones
  implicit none

  integer, parameter :: &
     BF_RECONSTRUCT = 0, &
     BF_GEQDSK      = 1, &
     BF_DIVA        = 2, &
     BF_PARAMETRIC  = 3, &
     BF_COILS       = 4, &
     BF_GRID_A      = 5

  integer, parameter :: BF_MAX_CONFIG = 5


  integer :: i_config(0:BF_MAX_CONFIG)

  contains
!=======================================================================



!=======================================================================
  subroutine setup_bfield_configuration

  integer, parameter :: iu = 24

  ! load configuration on first processor
  if (mype == 0) then
     open  (iu, file=Bfield_input_file)
     write (6,1000)
     write (6, *) 'Magnetic field input: '
     call load_reconstruct_config (iu, i_config(BF_RECONSTRUCT))
     call read_G_EQDSK_config     (iu, i_config(BF_GEQDSK),     Prefix)
     call read_polygones_config   (iu, i_config(BF_COILS),      Prefix)
     close (iu)
  endif

  call wait_pe()
  call broadcast_inte (i_config, BF_MAX_CONFIG)

 1000 format (/ '========================================================================')
  end subroutine setup_bfield_configuration
!=======================================================================



!=======================================================================
! return magnetic field components using Cylindrical coordinates
!=======================================================================
  function get_Bf_Cyl(r) result(Bf)
  real*8, intent(in) :: r(3)
  real*8             :: Bf(3)


  Bf = 0.d0

  if (i_config(BF_GEQDSK) == 1) Bf = Bf + get_Bcyl_geqdsk(r)
  if (i_config(BF_COILS)  == 1) Bf = Bf + get_Bcyl_polygones(r)


  end function get_Bf_Cyl
!=======================================================================



!=======================================================================
! return magnetic field components using Cartesian coordinates
!=======================================================================
  function get_Bf_Cart(x) result(Bf)
  real*8, intent(in) :: x(3)
  real*8             :: Bf(3)


  Bf = 0.d0

  if (i_config(BF_GEQDSK) == 1) Bf = Bf + get_Bcart_geqdsk(x)
  if (i_config(BF_COILS)  == 1) Bf = Bf + get_Bcart_polygones(x)

  end function get_Bf_Cart
!=======================================================================


end module bfield
