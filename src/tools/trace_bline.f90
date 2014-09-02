!===============================================================================
! Field line tracing
!
! Input (taken from run control file):
!    x_start(3)         Initial position for tracing (see Trace_Coords)
!    Grid_File          Provide initial positions for several field lines
!			Either "x_start" or "Grid_File" is required
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
!
!    Output_File
!    Output_Format      = 1: Cartesian coordinates
!                       = 2: Cylindrical coordinates
!===============================================================================
subroutine trace_bline
  use run_control, only: x_start, Trace_Step, N_steps, Trace_Method, Trace_Coords, &
                         Grid_File, Output_File, Output_Format
  use fieldline
  use parallel
  use math
  use grid
  implicit none

  integer, parameter :: iu = 42

  type(t_fieldline) :: F
  real*8, dimension(:,:), pointer :: grid_ptr
  real*8  :: x(3), y(3)
  integer :: i, iflag


  ! Field line tracing is performed on the 1st processor only!
  if (firstP) then
     write (6, *),'Tracing field line(s), output in: ', adjustl(trim(Output_File))
     write (6, *)
     write (6,  1001) N_steps
     write (6, *)
  else
     return
  endif
  open  (iu, file=Output_File)


  ! select initial position(s) for tracing
  if (sum(x_start) .ne. 0.d0) then
     grid_ptr => new_grid(1, log_progress=.false.)
     grid_ptr(1,:) = x_start
  else
     call read_grid (Grid_file, log_progress=.false., use_coordinates=COORDINATES(min(Trace_Coords,2)))
  endif

 
  ! main loop
  write (6,  1002)
  field_line_loop: do
     call get_next_grid_point (iflag, x_start)
     if (iflag < 0) exit field_line_loop

     ! trace one field line
     call F%init(x_start, Trace_Step, Trace_Method, Trace_Coords)
     call coord_trans (x_start, Trace_Coords, y, Output_Format)
     write (iu, 1003) y
     write (6,  1004) y

     do i=1,N_steps
        x = F%next_step()

        call coord_trans (x, Trace_Coords, y, Output_Format)
        write (iu, 1003) y
     enddo
     write (iu, *)
  enddo field_line_loop
  close (iu)

 1001 format (8x,'Number of trace steps: ',i8)
 1002 format (8x,'Initial position:')
 1003 format (3(e25.18,1x))
 1004 format (8x,3(e18.10,1x))
end subroutine trace_bline
