!===============================================================================
! Generate (axisymmetric) separatrix
!===============================================================================
subroutine generate_separatrix
  use iso_fortran_env
  use run_control, only: x_start, Output_File, Output_Format, Label, offset, Trace_Step
  use equilibrium
  use separatrix
  use parallel
  implicit none

  character(len=1)   :: c
  type(t_separatrix) :: S
  real(real64) :: X(2), ts
  integer      :: i


  if (firstP) then
     write (6, *) 'Generate (axisymmetric) separatrix'
     write (6, *)
  else
     return
  endif


  ! set trace step
  ts = Trace_Step
  if (ts > offset/10.d0) ts = offset / 10.d0


  if (Label .ne. '') Label = trim(Label)//'_'
  do i=1,nx_max
     if (Xp(i)%undefined) cycle

     write (6, *) i, Xp(i)%X
     write (c, '(i0)') i
     call S%generate(i, 2.d0, offset=offset, trace_step=ts)
     if (Output_Format == 1) then
        call S%plot('separatrix_'//trim(Label)//'X'//trim(c), parts=.true.)
     else
        call S%plot('separatrix_'//trim(Label)//'X'//trim(c), parts=.false.)
     endif
  enddo

end subroutine generate_separatrix
