!===============================================================================
! Generate (axisymmetric) separatrix
!===============================================================================
subroutine generate_separatrix
  use iso_fortran_env
  use run_control, only: N_psi, N_theta, Output_File, Output_Format, Label, offset, Trace_Step
  use equilibrium
  use separatrix
  use parallel
  use curve2D
  use grid
  implicit none

  character(len=1)   :: c
  type(t_separatrix) :: S
  type(t_curve)      :: S0
  type(t_grid)       :: G
  real(real64)       :: X(2), ts, t
  integer            :: i, n(2), j


  if (firstP) then
     write (6, *) 'Generate (axisymmetric) separatrix'
     write (6, *)
  else
     return
  endif


  ! set trace step
  ts = Trace_Step
  if (ts > offset/10.d0) ts = offset / 10.d0


  ! set base label
  if (Label .ne. '') Label = trim(Label)//'_'


  n = N_psi
  if (N_psi <= 0) then
     n(1) = 1
     n(2) = nx_max
  endif
  if (N_psi > nx_max) then
     write (6, 9001) N_psi
     stop
  endif

  do i=n(1),n(2)
     if (Xp(i)%undefined) cycle
     write (6, *) i, Xp(i)%X
     write (c, '(i0)') i
     call S%generate(i, 2.d0, offset=offset, trace_step=ts)
     select case(Output_Format)
     ! plot all branches in separate files
     case(1)
        call S%plot('separatrix_'//trim(Label)//'X'//trim(c), parts=.true.)

     ! plot all branches in one file
     case(2)
        call S%plot('separatrix_'//trim(Label)//'X'//trim(c), parts=.false.)

     ! generate grid on main part of separatrix
     case(3)
        S0 = connect(S%M1%t_curve, S%M2%t_curve)
        write (6, 1000) S0%length()
        call S0%setup_length_sampling()
        call G%new(CYLINDRICAL, UNSTRUCTURED, FIXED_COORD3, N_theta+1)
        do j=0,N_theta
           t = 1.d0 * j / N_theta
           call S0%sample_at(t, x)
           G%x(j+1,1:2) = x
        enddo
        call G%store(filename='separatrix_'//trim(Label)//'X'//trim(c)//'.grid')

     case default
        write (6, 9000) Output_Format
        stop
     end select
  enddo

 1000 format(8x,'Separatrix length [cm]: ', f12.5)
 9000 format('error in subroutine generate_separatrix: invalid Output_Format = ',i0)
 9001 format('error in subroutine generate_separatrix: invalid X-point N_psi = ',i0)
end subroutine generate_separatrix
