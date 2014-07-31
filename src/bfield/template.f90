module TEMPLATE
  implicit none

  private

  public :: load_TEMPLATE_configuration, &
            get_Bf_r3_TEMPLATE, get_Bf_x3_TEMPLATE

  contains
  !=====================================================================



  !=====================================================================
  ! load configuration and setup related variables
  !=====================================================================
  subroutine load_TEMPLATE_configuration (iu)
  integer, intent(in) :: iu
  end subroutine load_TEMPLATE_configuration
  !=====================================================================


  !=====================================================================
  ! return magnetic field components using cylindrical coordinates (R, Z, phi)
  !=====================================================================
  function get_Bf_r3_TEMPLATE(r)
  real*8             :: get_Bf_r3_TEMPLATE(3)
  real*8, intent(in) :: r(3)

  end function get_Bf_r3_TEMPLATE
  !=====================================================================


  !=====================================================================
  ! return magnetic field components using Cartesian coordinates
  !=====================================================================
  function get_Bf_x3_TEMPLATE(x)
  real*8             :: get_Bf_x3_TEMPLATE(3)
  real*8, intent(in) :: x(3)

  end function get_Bf_x3_TEMPLATE
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
