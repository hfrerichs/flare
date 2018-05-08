!===============================================================================
! Interface to magnetic field data
!
! The magnetic field strength is given in Gauss. Spatial coordinates are
! in cm (centi meter), the azimuth angle is given in rad (radians).
!===============================================================================
module bfield
  use iso_fortran_env
  use abstract_bfield
  implicit none
  private

  integer, public, parameter :: &
     BF_RECONSTRUCT  = 0, &
     BF_EQ2D         = 1, &
     BF_COILS        = 2, &
     BF_M3DC1        = 3, &
     BF_INTERPOLATEB = 4, &
     BF_SPLINEB      = 5, &
     BF_HINT         = 6, &
     BF_TOR_HARMONICS= 7, &
     BF_MARSF        = 8

  integer, public, parameter :: BF_MAX_CONFIG = 8


  integer, public :: iconfig(0:BF_MAX_CONFIG), icall(2)


  type t_bfield_wrapper
     class(t_bfield), allocatable :: content
     contains
     procedure :: broadcast => broadcast_wrapper
     procedure :: get_Bf    => get_Bf_wrapper
     procedure :: get_JBf   => get_JBf_wrapper
  end type t_bfield_wrapper
  type(t_bfield_wrapper), dimension(:), allocatable :: bfield_group
  integer :: ngroup


  public :: &
     setup_bfield_configuration, &
     get_Bf_Cyl, &
     get_JBf_Cyl, &
     get_Bf_Cyl_non2D, &
     get_Bf_Cart, &
     get_divB, &
     finished_bfield

  contains
!=======================================================================


!=======================================================================
  subroutine broadcast_wrapper(this)
  use parallel
  use bspline_group
  class(t_bfield_wrapper) :: this

  integer, parameter :: &
     TYPE_BSPLINE_POTENTIAL = 1, &
     TYPE_BSPLINE_BFIELD = 2

  integer :: itype
  associate (bf => this%content)


  ! get dynamic type of bfield group from process 0
  if (mype == 0) then
     select type (bf)
     type is(t_bspline_potential)
        itype = TYPE_BSPLINE_POTENTIAL

     type is(t_bspline_bfield)
        itype = TYPE_BSPLINE_BFIELD

     class default
        write (6, 9000)
        stop

     end select
  endif
 9000 format("error: unsupported dynamic type in t_bfield_wrapper%broadcast!")
  call broadcast_inte_s(itype)


  ! set dynamic type on process > 0
  if (mype > 0) then
     select case(itype)
     case(TYPE_BSPLINE_POTENTIAL)
        allocate (t_bspline_potential :: this%content)

     case(TYPE_BSPLINE_BFIELD)
        allocate (t_bspline_bfield    :: this%content)

     end select
  endif


  ! broadcast data
  call bf%broadcast()

  end associate
  end subroutine broadcast_wrapper
!=======================================================================



!=======================================================================
  function get_Bf_wrapper(this, r) result(Bf)
  class(t_bfield_wrapper)  :: this
  real(real64), intent(in) :: r(3)
  real(real64)             :: Bf(3)


  Bf = this%content%get_Bf(r)

  end function get_Bf_wrapper
!=======================================================================



!=======================================================================
  function get_JBf_wrapper(this, r) result(JBf)
  class(t_bfield_wrapper)  :: this
  real(real64), intent(in) :: r(3)
  real(real64)             :: JBf(3,3)


  JBf = this%content%get_JBf(r)

  end function get_JBf_wrapper
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
  use HINT
  use toroidal_harmonics
  use marsf

  integer, parameter :: iu = 24

  logical      :: ex
  real(real64) :: r(3), Bf(3)
  integer      :: i, irun


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
     call         HINT_load       (iu, iconfig(BF_HINT        ))
     call tor_harmonics_load      (iu, iconfig(BF_TOR_HARMONICS))
     call        marsf_load       (iu, iconfig(BF_MARSF       ))

     do irun=1,2
        ngroup = 0
        call read_bspline_input(iu, irun)

        if (irun == 1) allocate(bfield_group(ngroup))
     enddo
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


  if (nprs == 1) return
  call wait_pe()
  call broadcast_inte (iconfig, BF_MAX_CONFIG+1)

  if (iconfig(BF_RECONSTRUCT ) == 1) call broadcast_mod_reconstruct()
  if (iconfig(BF_EQ2D        ) == 1) call broadcast_mod_equilibrium()
  if (iconfig(BF_COILS       ) == 1) call broadcast_mod_polygones()
  if (iconfig(BF_M3DC1       ) == 1) call        m3dc1_broadcast()
  if (iconfig(BF_INTERPOLATEB) == 1) call interpolateB_broadcast()
  if (iconfig(BF_HINT        ) == 1) call         HINT_broadcast()
  if (iconfig(BF_TOR_HARMONICS) == 1) call tor_harmonics_broadcast()
  if (iconfig(BF_MARSF       ) == 1) call        marsf_broadcast()

  call broadcast_inte_s(ngroup)
  if (mype > 0) allocate(bfield_group(ngroup))
  do i=1,ngroup
     call bfield_group(i)%broadcast()
  enddo


 1000 format (/ '========================================================================')
 1001 format (8x,'Automatic setup of toroidal field direction: ',i2)
 1002 format (8x,'Automatic setup of poloidal field direction: ',i2)
 1009 format (8x,'WARNING: toroidal and poloidal field direction undefined!')
  end subroutine setup_bfield_configuration
!=======================================================================



!=======================================================================
  subroutine read_bspline_input(iu, irun)
  use run_control, only: Prefix
  use bspline_group
  use string
  integer, intent(in) :: iu, irun

  integer, parameter :: nsource_max = 128, ngroup_max = 256

  character(len=256), dimension(nsource_max) :: source = ''
  real(real64),       dimension(nsource_max) :: amplify = 1.d0
  real(real64),       dimension(nsource_max,ngroup_max) :: scale_subset = 1.d0
  character(len=256) :: filename, src_type
  integer :: i

  namelist /Bspline_Input/ source, amplify, scale_subset


  rewind (iu);   read (iu, Bspline_Input, end=9000)
  do i=1,nsource_max
  if (source(i) /= '') then
     ngroup = ngroup + 1
     if (irun == 1) cycle

     src_type = parse_string(source(i), 1, ':')
     filename = trim(Prefix)//parse_string(source(i), 2, ':')
     bfield_group(ngroup)%content = t_bspline_group(filename, src_type, amplify(i), scale_subset(i,1:ngroup_max))
  endif
  enddo

 9000 return
  end subroutine read_bspline_input
!=======================================================================



!=======================================================================
! return magnetic field components using Cylindrical coordinates
!=======================================================================
  function get_Bf_Cyl(r) result(Bf)
  use equilibrium
  use polygones
  use m3dc1
  use interpolateB
  use HINT
  use toroidal_harmonics
  use marsf
  real*8, intent(in) :: r(3)
  real*8             :: Bf(3), Bf_T(3), r_m(3)
  integer :: i


  Bf = 0.d0

  if (iconfig(BF_EQ2D)          == 1) Bf = Bf + get_Bf_eq2D(r)
  if (iconfig(BF_COILS)         == 1) Bf = Bf + get_Bcyl_polygones(r)
  if (iconfig(BF_M3DC1)         == 1) Bf = Bf +        m3dc1_get_Bf(r)
  if (iconfig(BF_INTERPOLATEB)  == 1) Bf = Bf + interpolateB_get_Bf(r)
  if (iconfig(BF_HINT)          == 1) Bf = Bf +         HINT_get_Bf(r)
  if (iconfig(BF_TOR_HARMONICS) == 1) Bf = Bf + tor_harmonics_get_Bf(r)
  if (iconfig(BF_MARSF)         == 1) Bf = Bf +        marsf_get_Bf(r)

  r_m(1:2) = r(1:2) / 1.d2;   r_m(3) = r(3)
  Bf_T = 0.d0
  do i=1,ngroup
     Bf_T = Bf_T + bfield_group(i)%get_Bf(r_m)
  enddo
  Bf = Bf + Bf_T * 1.d4
  icall(1) = icall(1) + 1


  end function get_Bf_Cyl
!=======================================================================
  function get_Bf_Cyl_non2D(r) result(Bf)
  use equilibrium
  use polygones
  use m3dc1
  use interpolateB
  use HINT
  use toroidal_harmonics
  use marsf
  real*8, intent(in) :: r(3)
  real*8             :: Bf(3)

  real(real64) :: r_m(3), Bf_T(3)
  integer :: i


  Bf = 0.d0

  if (iconfig(BF_COILS)         == 1) Bf = Bf + get_Bcyl_polygones(r)
  if (iconfig(BF_M3DC1)         == 1) Bf = Bf +        m3dc1_get_Bf(r)
  if (iconfig(BF_INTERPOLATEB)  == 1) Bf = Bf + interpolateB_get_Bf(r)
  if (iconfig(BF_HINT)          == 1) Bf = Bf +         HINT_get_Bf(r)
  if (iconfig(BF_TOR_HARMONICS) == 1) Bf = Bf + tor_harmonics_get_Bf(r)
  if (iconfig(BF_MARSF)         == 1) Bf = Bf +        marsf_get_Bf(r)

  r_m(1:2) = r(1:2) / 1.d2;   r_m(3) = r(3)
  Bf_T = 0.d0
  do i=1,ngroup
     associate (bf => bfield_group(i)%content)
     select type(bf)
     class is(t_equi2d)
        cycle

     class default
        Bf_T = Bf_T + bfield_group(i)%get_Bf(r_m)
     end select
     end associate
  enddo
  Bf = Bf + Bf_T * 1.d4
  icall(1) = icall(1) + 1


  end function get_Bf_Cyl_non2D
!=======================================================================



!=======================================================================
! calculate Jacobian
! WARNING: presently only for interpolateA - components
!=======================================================================
  function get_JBf_Cyl(r) result(J)
  use equilibrium
  use polygones
  use interpolateB
  use toroidal_harmonics
  real(real64), intent(in) :: r(3)
  real(real64)             :: J(3,3)

  real(real64) :: r_m(3)
  integer :: i


  J = 0.d0

  if (iconfig(BF_EQ2D   )       == 1) J = J + get_JBf_eq2D(r)
  if (iconfig(BF_COILS)         == 1) J = J + get_JBf_cyl_polygones(r)
  if (iconfig(BF_INTERPOLATEB)  == 1) J = J + interpolateB_get_JBf(r) / 100.d0
  if (iconfig(BF_TOR_HARMONICS) == 1) J = J + tor_harmonics_get_JBf(r)

  r_m(1:2) = r(1:2) / 1.d2;   r_m(3) = r(3)
  do i=1,ngroup
     J = J + bfield_group(i)%get_JBf(r_m)
  enddo

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
  real(real64), intent(in) :: x(3)
  real(real64)             :: Bf(3)

  real(real64) :: r(3), Bcyl(3), sin_phi, cos_phi


  ! coordinate transformation (Cartesian to cylindrical coordinates)
  r(1)    = dsqrt(x(1)**2 + x(2)**2)
  r(2)    = x(3)
  r(3)    = datan2(x(2), x(1))
  cos_phi = x(1) / r(1)
  sin_phi = x(2) / r(1)


  Bcyl    = get_Bf_Cyl(r)
  Bf(1)   = Bcyl(1) * cos_phi - Bcyl(3) * sin_phi
  Bf(2)   = Bcyl(1) * sin_phi + Bcyl(3) * cos_phi
  Bf(3)   = Bcyl(2)
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
