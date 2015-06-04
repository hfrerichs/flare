!===============================================================================
! Generate (axisymmetric) separatrix
!===============================================================================
subroutine generate_separatrix
  use iso_fortran_env
  use run_control, only: x_start, Output_File
  use equilibrium
  use separatrix
  use parallel
  implicit none

  character(len=1)   :: c
  type(t_separatrix) :: S
  real(real64) :: X(2)
  integer      :: i


  if (firstP) then
     write (6, *) 'Generate (axisymmetric) separatrix'
     write (6, *)
  else
     return
  endif


  do i=1,nx_max
     if (Xp(i)%R_estimate <= 0.d0) cycle

     write (6, *) i, Xp(i)%X
     write (c, '(i0)') i
     call S%generate(i, 1, 2.d0)
     call S%plot('separatrix_X'//trim(c), parts=.true.)
  enddo

end subroutine generate_separatrix
