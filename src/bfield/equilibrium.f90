!===============================================================================
!===============================================================================
module equilibrium
  use bfield
  use curve2D
  implicit none


  ! Direction of toroidal magnetic field and plasma current
  ! +1: positive direction, i.e. counter-clockwise
  ! -1: negative direction, i.e. clockwise
  !  0: no equilibrium defined
  integer :: &
     Bt_sign  = 0, &
     Ip_sign  = 0


  ! equilibrium id (in module bfield)
  integer :: i_equi = -1


  ! Interface for specific functions to be set externally
  interface
     function get_Psi_interface(r)
     real*8, intent(in) :: r(3)
     real*8             :: get_Psi_interface
     end function get_Psi_interface

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
  end interface

  ! return poloidal magnetic flux
  procedure(get_Psi_interface), pointer :: get_Psi

  ! return poloidal magnetic flux at magnetic axis
  procedure(Psi_axis_interface), pointer :: Psi_axis

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
  subroutine setup_equilibrium (report)
  use geqdsk
  logical, optional :: report

  integer :: i


  export_PFC => null()
  i = i_config(BF_GEQDSK) + i_config(BF_DIVA)

  ! no equilibrium defined?
  if (i == 0) return

  ! multiple equilibrium definitions?
  if (i > 1) then
     write (6, *) 'error: multiple equilibrium fields defined!'
     stop
  endif
  
  ! find the correct equilibrium id and setup some pointers
  if (i_config(BF_GEQDSK) == 1) then
     i_equi = BF_GEQDSK
     equilibrium_provides_PFC => geqdsk_provides_PFC
     export_PFC               => export_PFC_geqdsk
  elseif (i_config(BF_DIVA) == 1) then
     i_equi = BF_DIVA
  endif

  end subroutine setup_equilibrium
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
  subroutine Bf_pol_sub (n, s, y, f)
  integer, intent(in) :: n
  real*8, intent(in)  :: s, y(n)
  real*8, intent(out) :: f(n)

  ! n = 2, y(1) = R, y(2) = Z
  end subroutine Bf_pol_sub
!=======================================================================


end module equilibrium
