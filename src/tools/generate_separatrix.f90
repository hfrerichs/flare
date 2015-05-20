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

  type(t_separatrix) :: S
  real(real64) :: X(2)


  if (firstP) then
     write (6, *) 'Generate (axisymmetric) separatrix'
     write (6, *)
  else
     return
  endif


  ! find X-point
  X = find_lX()
  write (6, 1000) X
  ! generate separatrix
  call S%generate(X, 1, 2.d0)
  call S%plot('separatrix_lX', parts=.true.)


  ! now find upper X-point
  X = find_uX()
  write (6, 1000) X
  ! generate separatrix
  call S%generate(X, -1, -2.d0)
  call S%plot('separatrix_uX', parts=.true.)

 1000 format('Found X-point at: (',f8.3,', ',f8.3,')')
end subroutine generate_separatrix
