!===============================================================================
! Generate magnetic axis file
!
! Input (taken from run control file):
!    N_sym              Toroidal symmetry number, magnetic axis will be generated
!                       for 0 <= Phi <= 360/N_Sym deg
!    N_mult             Toroidal resolution (number of segments)
!    x_start            Reference position (cylindrical coordinates) from which
!                       the magnetic axis will be approximated
!
! optional:
!    N_points           Number of sample points for each toroidal slice
!    Trace_Step
!    Trace_Method
!===============================================================================
subroutine generate_magnetic_axis
  use iso_fortran_env
  use run_control, only: x_start, N_sym, N_mult, Output_File, &
                         N_points, Trace_Step, Trace_Method, Trace_Coords
  use math
  use parallel
  use dataset
  implicit none

  integer, parameter :: iu = 54

  type(t_dataset)    :: D
  character(len=120) :: filename
  real(real64)       :: R, Z, Phi
  integer            :: i


  Output_File = 'magnetic_axis.dat'
  if (firstP) then
     write (6, *) 'Generate magnetic axis, output in: ', adjustl(trim(Output_File))
     write (6, *)
  else
     return
  endif


  ! set default number of sample points for each toroidal slice
  N_points = 256

  ! Use cylindrical coordinates
  Trace_Coords = CYLINDRICAL


  ! 1. Generate Poincare plot at N_mult toroidal slices
  write (6, *) '1st step:'
  call poincare_plot


  ! 2. Combine output
  write (6, *)
  write (6, *) '2nd step:'
  write (6, *) 'Combine output'
  open  (iu, file=Output_File)
  write (iu, *) N_mult, N_sym
  do i=0,N_mult-1
     write (filename, '(i8)') i
     filename = 'magnetic_axis_'//trim(adjustl(filename))//'.dat'
     call D%load(filename, 2, report=.false.)
     R   = sum(D%x(:,1)) / D%nrow
     Z   = sum(D%x(:,2)) / D%nrow
     Phi = 360.d0 / N_sym / N_mult * i
     write (iu, *) R, Z, Phi
  enddo
  close (iu)

end subroutine generate_magnetic_axis
