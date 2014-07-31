!===============================================================================
! Field line tracing
!
! Input (taken from run control file):
!    x_start(3)         Initial position for tracing (see Trace_Coords)
!
!    Trace_Step         Size of trace steps (see Trace_Coords)
!
!    N_steps            Number of trace steps
!
!    Trace_Method	> 0: Integration (see module ODE_solver)
!                       = 0: Reconstruction from field aligned grid
!
!    Trace_Coords	= 1: solve d(x,y,z)/dl   = (Bx,By,Bx)/B
!                            x_start: Cartesian coordinates [cm]
!                            Trace_step: straight segment [cm]
!
!                       = 2: solve d(R,Z,phi)/ds = (BR,BZ,Bphi/R)/B
!                            x_start: Cylindrical coordinates [cm,deg]
!                            Trace_step: arc segment [cm]
!
!                       = 3: solve d(R,Z)/dphi   = (BR,BZ)*R/Bphi
!                            x_start: Cylindrical coordinates [cm,deg]
!                            Trace_step: arc segment [deg]
!    Output_File
!    Output_Format      = 1: Cartesian coordinates
!                       = 2: Cylindrical coordinates
!===============================================================================
subroutine trace_bline
  use run_control, only: x_start, Trace_Step, N_steps, Trace_Method, Trace_Coords, &
                         Output_File, Output_Format
  use fieldline
  use parallel
  use math
  implicit none

  integer, parameter :: iu = 42

  type(t_fieldline) :: F
  real*8  :: x(3), y(3)
  integer :: i


  ! Field line tracing is performed on the 1st processor only!
  if (firstP) then
     write (6, 1000), Output_File
  else
     return
  endif


  open  (iu, file=Output_File)
  call F%init(x_start, Trace_Step, Trace_Method, Trace_Coords)
  call coord_trans (x_start, Trace_Coords, y, Output_Format)
  write (iu, 1002) y

  do i=1,N_steps
     x = F%next_step()

     call coord_trans (x, Trace_Coords, y, Output_Format)
     write (iu, 1002) y
  enddo
  close (iu)

 1000 format (3x,'- tracing field line, output in: ',a120)
 1002 format (3(e25.18,1x))
end subroutine trace_bline
