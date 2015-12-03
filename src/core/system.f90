!=======================================================================
module system
  use iso_fortran_env
  implicit none

  real(real64), parameter :: &
     epsilon_r64 = epsilon(real(1.0,real64)), &
     machine_precision = epsilon_r64
end module system
!=======================================================================





!=======================================================================
subroutine progress(i, n, delta)
  use iso_fortran_env
  implicit none
  integer,      intent(in) :: i, n
  real(real64), intent(in) :: delta
  !integer,      intent(in), optional :: incr

  integer, save :: D
  real(real64)  :: F


  if (i == 0) D = 0

  F = 1.d0 * i / n
  if (F >= D*delta  .or.  i==n) then
     !write (6, *) int(F*100), i, n
     !write (6, 1000) int(F*100)
     write (6, 2000) F*100
     D = D + 1
  endif

 1000 format(8x,i3,' %')
 2000 format(8x,f0.2,' %')
end subroutine progress
!=======================================================================
