!===============================================================================
! Interface to magnetic field data
!
! The magnetic field strength is given in Gauss. Spatial coordinates are
! in cm (centi meter), the azimuth angle is given in rad (radians).
!===============================================================================
module bfield
  implicit none

  integer, parameter :: &
     BF_RECONSTRUCT = 0, &
     BF_EQ2D        = 1, &
     BF_PARAMETRIC  = 3, &
     BF_COILS       = 4, &
     BF_GRID_A      = 5

  integer, parameter :: BF_MAX_CONFIG = 5


  integer :: iconfig(0:BF_MAX_CONFIG)

  contains
!=======================================================================



!=======================================================================
  subroutine setup_bfield_configuration
  use run_control
  use parallel
  use reconstruct
  use equilibrium
  use polygones

  integer, parameter :: iu = 24

  ! load configuration on first processor
  if (firstP) then
     open  (iu, file=Bfield_input_file)
     write (6,1000)
     write (6, *) 'Magnetic field input: '
     call load_reconstruct_config (iu, iconfig(BF_RECONSTRUCT))
     call load_equilibrium_config (iu, iconfig(BF_EQ2D))
     call read_polygones_config   (iu, iconfig(BF_COILS),      Prefix)
     close (iu)
  endif

  call wait_pe()
  call broadcast_inte (iconfig, BF_MAX_CONFIG)

  if (iconfig(BF_RECONSTRUCT) == 1) call broadcast_mod_reconstruct()
  if (iconfig(BF_EQ2D       ) == 1) call broadcast_mod_equilibrium()
  if (iconfig(BF_COILS      ) == 1) call broadcast_mod_polygones()



 1000 format (/ '========================================================================')
  end subroutine setup_bfield_configuration
!=======================================================================



!=======================================================================
! return magnetic field components using Cylindrical coordinates
!=======================================================================
  function get_Bf_Cyl(r) result(Bf)
  use equilibrium
  use polygones
  real*8, intent(in) :: r(3)
  real*8             :: Bf(3)


  Bf = 0.d0

  if (iconfig(BF_EQ2D)   == 1) Bf = Bf + get_Bcyl_eq2D(r)
  if (iconfig(BF_COILS)  == 1) Bf = Bf + get_Bcyl_polygones(r)


  end function get_Bf_Cyl
!=======================================================================



!=======================================================================
! return magnetic field components using Cartesian coordinates
!=======================================================================
  function get_Bf_Cart(x) result(Bf)
  use equilibrium
  use polygones
  real*8, intent(in) :: x(3)
  real*8             :: Bf(3)


  Bf = 0.d0

  if (iconfig(BF_EQ2D)   == 1) Bf = Bf + get_Bcart_eq2D(x)
  if (iconfig(BF_COILS)  == 1) Bf = Bf + get_Bcart_polygones(x)

  end function get_Bf_Cart
!=======================================================================


end module bfield
