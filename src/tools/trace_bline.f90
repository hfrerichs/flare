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
!    N_steps            Number of trace steps, if N_steps < 0 then tracing is
!                       stopped at the boundary while |N_steps| is the upper limit.
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
!    Output_Format      = 0: Magnetic flux coordinates (Theta, PsiN, Phi)
!                       = 1: Cartesian coordinates
!                       = 2: Cylindrical coordinates
!===============================================================================
subroutine trace_bline
  use iso_fortran_env
  use run_control, only: x_start, Trace_Step, N_steps, Trace_Method, Trace_Coords, &
                         Grid_File, Output_File, Output_Format
  use fieldline
  use parallel
  use math
  use usr_grid
  implicit none

  integer, parameter :: iu = 42

  type(t_fieldline) :: F
  real(real64), dimension(:,:), pointer :: grid_ptr
  integer       :: i, iflag
  logical       :: Stop_at_Boundary


  ! Initialize
  Stop_at_Boundary = .false.
  if (N_steps < 0) then
     Stop_at_Boundary = .true.
     N_steps = abs(N_steps)
  endif


  ! Field line tracing is performed on the 1st processor only!
  if (firstP) then
     write (6, *),'Tracing field line(s), output in: ', adjustl(trim(Output_File))
     write (6, *)
     write (6,  1001) N_steps
     write (6, *)
  else
     return
  endif

  ! Open output file and write information about coordinate system
  open  (iu, file=Output_File)
  select case(Output_Format)
  case(0)
     write(iu, 2000)
  case(1)
     write(iu, 2001)
  case(2)
     write(iu, 2002)
  case(3)
     write(iu, 2002)
  end select


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
     write (iu, 1003) output_coordinates()
     write (6,  1004) output_coordinates()

     do i=1,N_steps
        call F%trace_1step()

        ! intersect boundary?
        if (Stop_at_Boundary  .and.  F%intersect_boundary()) then
           write (iu, 1003) output_coordinates()
           write (6, *) 'Field line tracing is stopped at the boundary'
           exit
        endif

        write (iu, 1003) output_coordinates()
     enddo
     write (iu, *)
  enddo field_line_loop
  close (iu)

 1001 format (8x,'Number of trace steps: ',i8)
 1002 format (8x,'Initial position:')
 1003 format (3(e25.18,1x))
 1004 format (8x,3(e18.10,1x))
 2000 format ('# Flux coordinates: Theta[rad], PsiN, Phi[rad]')
 2001 format ('# Cartesian coordinates: x[cm], y[cm], z[cm]')
 2002 format ('# Cylindrical coordinates: R[cm], Z[cm], Phi[rad]')
  contains
!.......................................................................
  function output_coordinates () result (y)
  real(real64) :: y(3)


  if (Output_Format == 0) then
     y    = F%get_flux_coordinates()
  else
     call coord_trans (F%yc, Trace_Coords, y, Output_Format)
  endif

  end function output_coordinates
!.......................................................................
end subroutine trace_bline
