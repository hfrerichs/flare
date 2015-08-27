module system
  use iso_fortran_env
  implicit none

  real(real64), parameter :: &
     epsilon_r64 = epsilon(real(1.0,real64)), &
     machine_precision = epsilon_r64

end module system
