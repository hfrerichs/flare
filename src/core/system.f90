!=======================================================================
module system
  use iso_fortran_env
  implicit none

  real(real64), parameter :: &
     epsilon_r64 = epsilon(real(1.0,real64)), &
     machine_precision = epsilon_r64

  integer, parameter :: &
     STRICT      = 0, &
     FREE        = 1


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



!=======================================================================
subroutine irange_check(i, i1, i2, label, procedure_name)
  integer, intent(in) :: i, i1, i2
  character(len=*), intent(in) :: label, procedure_name


  if (i < i1  .or.  i > i2) then
     write (6, 9000) trim(label), i, i1, trim(label), i2, procedure_name
     stop
  endif

 9000 format('error: ',a,' = ',i0,' when ',i0,' <= ',a,' <= ',i0,' required in ',a,'!')

end subroutine irange_check
!=======================================================================
