!===============================================================================
! Generate homoclinic tangle for hyperbolic fixed point
!
! Input (taken from run control file):
! Grid_File             toroidal discretization of X-point (+ small offset for tracing)
! Phi_Output            toroidal reference position [deg]
!
! Trace_Step, Trace_Method
! Output_File
!===============================================================================
subroutine homoclinic_tangle
  use iso_fortran_env
  use run_control, only: Grid_File, Output_File, Trace_Step, Trace_Method, Phi_Output
  use grid
  use fieldline
  use boundary
  use parallel
  implicit none

  integer, parameter :: iu = 54

  type(t_fieldline)  :: F
  type(t_grid)       :: G
  real(real64)       :: y(3), yh(3), Dphi
  integer            :: i, ierr


  ! initialize
  if (firstP) then
     write (6, *) 'Generate homoclinic tangle for hyperbolic fixed point'
     write (6, *)
  endif
  call G%load(Grid_File)
  write (6, *)
  

  open  (iu, file=Output_File)
  grid_loop: do i=1,G%nodes()
     y = G%node(i)
     write (6, *) y(3) / pi * 180.d0

     ! set initial location
     call F%init(y, Trace_Step, Trace_Method, CYLINDRICAL)

     Dphi = abs(Phi_Output / 180.d0 * pi - y(3))
     call F%trace_Dphi(Dphi, .true., yh, ierr)

     if (ierr == 0) then
        write (iu, *) yh
     else
        write (iu, *)
     endif
  enddo grid_loop
  close (iu)

  return
end subroutine homoclinic_tangle
