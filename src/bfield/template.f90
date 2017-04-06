module TEMPLATE
  use iso_fortran_env
  implicit none

  private

  public :: TEMPLATE_load, &
            TEMPLATE_broadcast, &
            TEMPLATE_get_Bf, &
            TEMPLATE_get_JBf

  contains
  !=====================================================================



  !=====================================================================
  ! load configuration and setup related variables
  !=====================================================================
  subroutine TEMPLATE_load (iu, iconfig)
  integer, intent(in)  :: iu
  integer, intent(out) :: iconfig


  iconfig = 0

  end subroutine TEMPLATE_load
  !=====================================================================



  !=====================================================================
  ! broadcast data for parallel execution
  !=====================================================================
  subroutine TEMPLATE_broadcast()
  end subroutine TEMPLATE_broadcast
  !=====================================================================



  !=====================================================================
  ! return magnetic field components using cylindrical coordinates (R, Z, phi)
  !=====================================================================
  function TEMPLATE_get_Bf(r) result(Bf)
  real(real64), intent(in) :: r(3)
  real(real64)             :: Bf(3)

  end function TEMPLATE_get_Bf
  !=====================================================================



  !=====================================================================
  ! return Jacobian of magnetic field [T/m] in cylindrical coordinates
  !=====================================================================
  function TEMPLATE_get_JBf(r) result(JBf)
  real(real64), intent(in) :: r(3)
  real(real64)             :: JBf(3,3)

  end function TEMPLATE_get_JBf
  !=====================================================================



  !=====================================================================
  function get_Psi_TEMPLATE (r)
  real*8, intent(in) :: r(3)
  real*8             :: get_Psi_TEMPLATE

  end function get_Psi_TEMPLATE
  !=====================================================================

  !=====================================================================
  function get_DPsi_TEMPLATE (r, n1, n2)
  real*8, intent(in)  :: r(3)
  integer, intent(in) :: n1, n2
  real*8              :: get_DPsi_TEMPLATE

  end function get_DPsi_TEMPLATE
  !=====================================================================

end module TEMPLATE
