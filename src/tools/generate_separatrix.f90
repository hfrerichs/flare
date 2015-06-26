!===============================================================================
! Generate (axisymmetric) separatrix
!===============================================================================
subroutine generate_separatrix
  use iso_fortran_env
  use run_control, only: x_start, Output_File, Label
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


  if (Label .ne. '') Label = trim(Label)//'_'
  do i=1,nx_max
     if (Xp(i)%undefined) cycle

     write (6, *) i, Xp(i)%X
     write (c, '(i0)') i
     call S%generate(i, 2.d0)
     call S%plot('separatrix_'//trim(Label)//'X'//trim(c), parts=.true.)
  enddo

end subroutine generate_separatrix
