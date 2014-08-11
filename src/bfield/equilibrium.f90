!===============================================================================
!===============================================================================
module equilibrium
  use curve2D
  implicit none


!...............................................................................
! user defined parameters (to be set via configuration file)                   .

  character*120 :: &
     Data_File        = ''
  character*12  :: &
     Data_Format      = ''

  real*8 :: &
     R_axis_usr       = 0.d0, &        ! user defined position of magnetic axis
     Z_axis_usr       = 0.d0

  logical :: &
     use_PFC          = .true., &
     Current_Fix      = .true.

  integer :: &
     Diagnostic_Level = 0

  namelist /Equilibrium_Input/ &
     Data_File, Data_Format, use_PFC, Current_Fix, Diagnostic_Level, &
     R_axis_usr, Z_axis_usr
!...............................................................................



  ! Direction of toroidal magnetic field and plasma current
  ! +1: positive direction, i.e. counter-clockwise
  ! -1: negative direction, i.e. clockwise
  !  0: no equilibrium defined
  integer :: &
     Bt_sign  = 0, &
     Ip_sign  = 0


  ! equilibrium type
  integer :: i_equi = -1


  real*8 :: R_axis, Z_axis, Psi_sepx, Psi_axis


  ! Interface for specific functions to be set externally
  interface
     function get_Psi_interface(r) result(Psi)
     real*8, intent(in) :: r(3)
     real*8             :: Psi
     end function get_Psi_interface
     function get_D1Psi_interface(r) result(D1Psi)
     real*8, intent(in) :: r(3)
     real*8             :: D1Psi(2)
     end function get_D1Psi_interface
     function get_D2Psi_interface(r) result(D2Psi)
     real*8, intent(in) :: r(3)
     real*8             :: D2Psi(3)
     end function get_D2Psi_interface

     function Psi_axis_interface(phi)
     real*8, intent(in), optional :: phi
     real*8                       :: Psi_axis_interface
     end function Psi_axis_interface

     function logical_inquiry() result(l)
     logical :: l
     end function logical_inquiry

     subroutine export_curve(S)
     import :: t_curve
     type(t_curve), intent(out) :: S
     end subroutine export_curve

     function get_Bf_interface(y) result(Bf)
     real*8, intent(in) :: y(3)
     real*8             :: Bf(3)
     end function
  end interface


  ! get equilibrium magnetic field
  procedure(get_Bf_interface), pointer :: &
     ! ... in Cartesian coordinates
     get_BCart_eq2D, &
     ! ... in cylindrical coordinates
     get_BCyl_eq2D



  ! return poloidal magnetic flux
  procedure(get_Psi_interface), pointer :: get_Psi
  procedure(get_D1Psi_interface), pointer :: get_D1Psi
  procedure(get_D2Psi_interface), pointer :: get_D2Psi

!  ! return poloidal magnetic flux at magnetic axis
!  procedure(Psi_axis_interface), pointer :: Psi_axis

  ! inquire if equilibrium provides PFC setup
  procedure(logical_inquiry), pointer :: &
     equilibrium_provides_PFC => equilibrium_provides_PFC_generic
  procedure(export_curve), pointer    :: export_PFC


  logical, save :: initialized = .false.
!, get_equi_domain
!
  contains
!=======================================================================



!=======================================================================
  subroutine load_equilibrium_config (iu, iconfig)
  use run_control, only: Prefix
  use geqdsk
  integer, intent(in)  :: iu
  integer, intent(out) :: iconfig


! read user configuration
  rewind (iu)
  read   (iu, Equilibrium_Input, end=1000)
  iconfig = 1
  write (6, *)
  write (6, 1001)
  Data_File = trim(Prefix)//Data_File


! set default values
  export_PFC => null()


! determine equilibrium type
   i_equi = 1


!...
!  select case(Data_Format)
!  case ('geqdsk')
!  case default
!     write (6, *) 'error: ', Data_Format, ' is not a valid equilibrium type!'
!     stop
!  end select


! load equilibrium data
  call load_mod_geqdsk (Data_File, use_PFC, Current_Fix, Diagnostic_Level, R_axis, Z_axis, Psi_axis, Psi_sepx)


! set user defined values, if present
  if (R_axis_usr.gt.0.d0) then
     R_axis = R_axis_usr
     Z_axis = Z_axis_usr
  endif


! set dependent variables
  !psi_axis = pol_flux(magnetic_axis())
  !Psi_axis = get_Psi(magnetic_axis())
  !Psi_sepx = 

  call setup_equilibrium()

  return
 1000 iconfig = 0
 1001 format ('   - Axisymmetric (2D) MHD equilibrium:')
  end subroutine load_equilibrium_config
!=======================================================================



!=======================================================================
  subroutine setup_equilibrium()
  use geqdsk

  ! select case equilibrium
  get_BCart_eq2D => get_BCart_geqdsk
  get_BCyl_eq2D  => get_BCyl_geqdsk
  get_Psi        => get_Psi_geqdsk
  get_D1Psi      => get_D1Psi_geqdsk
  get_D2Psi      => get_D2Psi_geqdsk
  equilibrium_provides_PFC => geqdsk_provides_PFC
  export_PFC               => export_PFC_geqdsk

  end subroutine setup_equilibrium
!=======================================================================



!=======================================================================
  subroutine broadcast_mod_equilibrium()
  use parallel
  use geqdsk


  if (nprs == 1) return

  call broadcast_real_s (R_axis)
  call broadcast_real_s (Z_axis)
  call broadcast_real_s (Psi_axis)
  call broadcast_real_s (Psi_sepx)
  call broadcast_inte_s (i_equi)
  call wait_pe()


  ! select case equilibrium
  select case(i_equi)
  case(1)
     call broadcast_mod_geqdsk()
  end select
  if (mype > 0) call setup_equilibrium()

  end subroutine broadcast_mod_equilibrium
!=======================================================================



!=======================================================================
! PFC provided by equilibrium
!=======================================================================
  function equilibrium_provides_PFC_generic() result(l)
  logical :: l

  l = .false.
  end function equilibrium_provides_PFC_generic
!=======================================================================
! export axisymmetric boundary provided by equilibrium
!=======================================================================
!  subroutine export_PFC (S)
!  type(t_curve), intent(out) :: S
!  end subroutine export_PFC
!=======================================================================



!=======================================================================
  function pol_flux(r) result(psi)
  real*8, intent(in) :: r(3)
  real*8             :: psi

  end function pol_flux
!=======================================================================


!=======================================================================
! Provide position r=(R,Z,phi) [cm,cm,rad] of magnetic axis
! Optional input: phi (for non-axisymmetric equilibria)
! Default: phi = 0.0
!=======================================================================
  function magnetic_axis(phi) result(r)
  real*8, intent(in), optional :: phi
  real*8                       :: r(3)

  real*8 :: phi1


  ! scan optional arguments
  phi1 = 0.d0
  if (present(phi)) phi1 = phi
  r(3) = phi1

  ! set magnetic axis
  r(1) = R_axis
  r(2) = Z_axis


  end function magnetic_axis
!=======================================================================


!=======================================================================
  subroutine Bf_pol_sub (n, s, y, f)
  integer, intent(in) :: n
  real*8, intent(in)  :: s, y(n)
  real*8, intent(out) :: f(n)

  ! n = 2, y(1) = R, y(2) = Z
  end subroutine Bf_pol_sub
!=======================================================================

end module equilibrium
