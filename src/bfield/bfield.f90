!===============================================================================
! Interface to magnetic field data
!
! The magnetic field strength is given in Gauss. Spatial coordinates are
! in cm (centi meter), the azimuth angle is given in rad (radians).
!===============================================================================
module bfield
  use iso_fortran_env
  implicit none

  integer, parameter :: &
     BF_RECONSTRUCT  = 0, &
     BF_EQ2D         = 1, &
     BF_COILS        = 2, &
     BF_M3DC1        = 3, &
     BF_INTERPOLATEB = 4, &
     BF_SPLINEB      = 5, &
     BF_HINT         = 6, &
     BF_TOR_HARMONICS= 7, &
     BF_MARSF        = 8

  integer, parameter :: BF_MAX_CONFIG = 8


  integer :: iconfig(0:BF_MAX_CONFIG), icall(2)

  contains
!=======================================================================



!=======================================================================
  subroutine setup_bfield_configuration
  use run_control
  use parallel
  use reconstruct
  use magnetic_axis
  use equilibrium, only: load_equilibrium_config, broadcast_mod_equilibrium
  use polygones
  use m3dc1
  use interpolateB
  use splineB
  use HINT
  use toroidal_harmonics
  use marsf

  integer, parameter :: iu = 24

  logical      :: ex
  real(real64) :: r(3), Bf(3)


  icall = 0
  ! load configuration on first processor
  if (firstP) then
     inquire (file=Bfield_input_file, exist=ex)
     if (.not.ex) then
        write (6, *) 'configuration file: ', trim(Bfield_input_file)
        write (6, *) 'error: magnetic configuration does not exist!'
        stop
     endif

     open  (iu, file=Bfield_input_file)
     write (6,1000)
     write (6, *) 'Magnetic field input: '
     call load_reconstruct_config (iu, iconfig(BF_RECONSTRUCT ))
     call        m3dc1_load       (iu, iconfig(BF_M3DC1       ))
     call load_equilibrium_config (iu, iconfig(BF_EQ2D        ))
     call read_polygones_config   (iu, iconfig(BF_COILS       ),      Prefix)
     call interpolateB_load       (iu, iconfig(BF_INTERPOLATEB))
     call      splineB_load       (iu, iconfig(BF_SPLINEB     ))
     call         HINT_load       (iu, iconfig(BF_HINT        ))
     call tor_harmonics_load      (iu, iconfig(BF_TOR_HARMONICS))
     call        marsf_load       (iu, iconfig(BF_MARSF       ))
     close (iu)


     ! supplemental equilibrium information
     ! direction of toroidal magnetic field and plasma current
     ! if Bt_sign (and Ip_sign) is not set by equilibrium
     if (Bt_sign == 0) then
        r(3) = 0.d0
        r    = get_magnetic_axis(r(3))
        if (r(1) > 0.d0) then
           Bf   = get_Bf_cyl(r)

           Bt_sign = 1
           if (Bf(3) < 0.d0) Bt_sign = -1
           write (6, 1001) Bt_sign

           ! guess poloidal field direction:
           ! get BZ left of magnetic axis
           r(1) = 0.9d0 * r(1)
           Bf   = get_Bf_cyl(r)

           Ip_sign = 1
           ! switch sign if BZ is negative
           if (Bf(2) < 0.d0) Ip_sign = -1
           write (6, 1002) Ip_sign

        ! magnetic axis undefined, toroidal and poloidal field directions remain undefined!
        else
           write (6, 1009)
        endif
     endif
  endif


  call wait_pe()
  call broadcast_inte (iconfig, BF_MAX_CONFIG+1)

  if (iconfig(BF_RECONSTRUCT ) == 1) call broadcast_mod_reconstruct()
  if (iconfig(BF_EQ2D        ) == 1) call broadcast_mod_equilibrium()
  if (iconfig(BF_COILS       ) == 1) call broadcast_mod_polygones()
  if (iconfig(BF_M3DC1       ) == 1) call        m3dc1_broadcast()
  if (iconfig(BF_INTERPOLATEB) == 1) call interpolateB_broadcast()
  if (iconfig(BF_SPLINEB     ) == 1) call      splineB_broadcast()
  if (iconfig(BF_HINT        ) == 1) call         HINT_broadcast()
  if (iconfig(BF_TOR_HARMONICS) == 1) call tor_harmonics_broadcast()
  if (iconfig(BF_MARSF       ) == 1) call        marsf_broadcast()



 1000 format (/ '========================================================================')
 1001 format (8x,'Automatic setup of toroidal field direction: ',i2)
 1002 format (8x,'Automatic setup of poloidal field direction: ',i2)
 1009 format (8x,'WARNING: toroidal and poloidal field direction undefined!')
  end subroutine setup_bfield_configuration
!=======================================================================



!=======================================================================
! return magnetic field components using Cylindrical coordinates
!=======================================================================
  function get_Bf_Cyl(r) result(Bf)
  use equilibrium
  use polygones
  use m3dc1
  use interpolateB
  use splineB
  use HINT
  use toroidal_harmonics
  use marsf
  real*8, intent(in) :: r(3)
  real*8             :: Bf(3)


  Bf = 0.d0

  if (iconfig(BF_EQ2D)          == 1) Bf = Bf + get_Bf_eq2D(r)
  if (iconfig(BF_COILS)         == 1) Bf = Bf + get_Bcyl_polygones(r)
  if (iconfig(BF_M3DC1)         == 1) Bf = Bf +        m3dc1_get_Bf(r)
  if (iconfig(BF_INTERPOLATEB)  == 1) Bf = Bf + interpolateB_get_Bf(r)
  if (iconfig(BF_SPLINEB)       == 1) Bf = Bf +      splineB_get_Bf(r)
  if (iconfig(BF_HINT)          == 1) Bf = Bf +         HINT_get_Bf(r)
  if (iconfig(BF_TOR_HARMONICS) == 1) Bf = Bf + tor_harmonics_get_Bf(r)
  if (iconfig(BF_MARSF)         == 1) Bf = Bf +        marsf_get_Bf(r)
  icall(1) = icall(1) + 1


  end function get_Bf_Cyl
!=======================================================================
  function get_Bf_Cyl_non2D(r) result(Bf)
  use equilibrium
  use polygones
  use m3dc1
  use interpolateB
  use splineB
  use HINT
  use toroidal_harmonics
  use marsf
  real*8, intent(in) :: r(3)
  real*8             :: Bf(3)


  Bf = 0.d0

  if (iconfig(BF_COILS)         == 1) Bf = Bf + get_Bcyl_polygones(r)
  if (iconfig(BF_M3DC1)         == 1) Bf = Bf +        m3dc1_get_Bf(r)
  if (iconfig(BF_INTERPOLATEB)  == 1) Bf = Bf + interpolateB_get_Bf(r)
  if (iconfig(BF_SPLINEB)       == 1) Bf = Bf +      splineB_get_Bf(r)
  if (iconfig(BF_HINT)          == 1) Bf = Bf +         HINT_get_Bf(r)
  if (iconfig(BF_TOR_HARMONICS) == 1) Bf = Bf + tor_harmonics_get_Bf(r)
  if (iconfig(BF_MARSF)         == 1) Bf = Bf +        marsf_get_Bf(r)
  icall(1) = icall(1) + 1


  end function get_Bf_Cyl_non2D
!=======================================================================



!=======================================================================
! calculate Jacobian
! WARNING: presently only for interpolateA - components
!=======================================================================
  function get_JBf_Cyl(r) result(J)
  use equilibrium
  use splineB
  use interpolateB
  use toroidal_harmonics
  real(real64), intent(in) :: r(3)
  real(real64)             :: J(3,3)


  J = 0.d0

  if (iconfig(BF_EQ2D   )       == 1) J = J + get_JBf_eq2D(r)
  if (iconfig(BF_SPLINEB)       == 1) J = J + splineB_get_JBf(r)
  if (iconfig(BF_INTERPOLATEB)  == 1) J = J + interpolateB_get_JBf(r) / 100.d0
  if (iconfig(BF_TOR_HARMONICS) == 1) J = J + tor_harmonics_get_JBf(r)

  end function get_JBf_Cyl
!========================================================================



!========================================================================
! calculate divergence B (due to numerical representation of magnetic
! field data)
!========================================================================
  function get_divB(r) result(divB)
  real(real64), intent(in) :: r(3)
  real(real64)             :: divB

  real(real64) :: J(3,3), Bf(3)


  divB = 0.d0
  Bf   = get_Bf_Cyl(r)         ! in Gauss
  J    = get_JBf_Cyl(r) * 1.d2 ! in Gauss/cm
  divB = J(1,1) + Bf(1) / r(1) + J(3,3) / r(1) + J(2,2)

  end function get_divB
!========================================================================



!=======================================================================
! return magnetic field components using Cartesian coordinates
!=======================================================================
  function get_Bf_Cart(x) result(Bf)
  use equilibrium
  use polygones
  use m3dc1
  use interpolateB
  use splineB
  use HINT
  use toroidal_harmonics
  use marsf
  real*8, intent(in) :: x(3)
  real*8             :: Bf(3)

  real*8 :: r(3), Bcyl(3), sin_phi, cos_phi


  ! coordinate transformation (Cartesian to cylindrical coordinates)
  r(1)    = dsqrt(x(1)**2 + x(2)**2)
  r(2)    = x(3)
  r(3)    = datan2(x(2), x(1))
  cos_phi = x(1) / r(1)
  sin_phi = x(2) / r(1)


  ! collect all magnetic field components
  Bf      = 0.d0
  Bcyl    = 0.d0
  if (iconfig(BF_EQ2D)          == 1) Bcyl = Bcyl + get_Bf_eq2D(r)
  if (iconfig(BF_COILS)         == 1) Bf   = Bf   + get_Bcart_polygones(x)
  if (iconfig(BF_M3DC1)         == 1) Bcyl = Bcyl +        m3dc1_get_Bf(r)
  if (iconfig(BF_INTERPOLATEB)  == 1) Bcyl = Bcyl + interpolateB_get_Bf(r)
  if (iconfig(BF_SPLINEB)       == 1) Bcyl = Bcyl +      splineB_get_Bf(r)
  if (iconfig(BF_HINT)          == 1) Bcyl = Bcyl +         HINT_get_Bf(r)
  if (iconfig(BF_TOR_HARMONICS) == 1) Bcyl = Bcyl + tor_harmonics_get_Bf(r)
  if (iconfig(BF_MARSF)         == 1) Bcyl = Bcyl +        marsf_get_Bf(r)


  ! combine Cartesian and cylindrical components
  Bf(1) = Bf(1) + Bcyl(1) * cos_phi - Bcyl(3) * sin_phi
  Bf(2) = Bf(2) + Bcyl(1) * sin_phi + Bcyl(3) * cos_phi
  Bf(3) = Bf(3) + Bcyl(2)
  icall(2) = icall(2) + 1

  end function get_Bf_Cart
!=======================================================================



!=======================================================================
  subroutine reset_counter
  icall = 0
  end subroutine reset_counter
!=======================================================================



!=======================================================================
  subroutine finished_bfield
  use m3dc1
  use parallel
  use marsf

  if (iconfig(BF_M3DC1)         == 1) call m3dc1_close()
  if (iconfig(BF_MARSF)         == 1) call marsf_close()

  call wait_pe()
  call sum_inte_data(icall,2)
  if (firstP) then
     if (icall(1) > 0) write (6, *) icall(1), ' calls to Bf_cyl'
     if (icall(2) > 0) write (6, *) icall(2), ' calls to Bf_cart'
  endif


  end subroutine finished_bfield
!=======================================================================

end module bfield
