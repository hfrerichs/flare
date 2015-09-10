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
  use grid
  implicit none

  integer, parameter :: iu = 42

  type(t_fieldline) :: F
  type(t_grid)      :: G
  real(real64), dimension(:,:), pointer :: grid_ptr
  real(real64)  :: rb(3), tau, y(3), l, dl
  integer       :: i, ig
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
     call G%new(Trace_Coords, UNSTRUCTURED, 0, 1)
     G%x(1,:) = x_start
  else
     call G%load(Grid_file)
  endif

 
  ! main loop
  write (6,  1002)
  field_line_loop: do ig=1,G%nodes()
     x_start = G%node(ig, coordinates=min(Trace_Coords,2))

     ! trace one field line
     call F%init(x_start, Trace_Step, Trace_Method, Trace_Coords)
     y = output_coordinates()
     write (iu, 1003) y
     write (6,  1004) y

     l = 0.d0
     do i=1,N_steps
        dl = F%trace_1step()

        ! intersect boundary?
        if (Stop_at_Boundary  .and.  F%intersect_boundary(tau=tau)) then
           rb = output_coordinates()
           rb = y + tau * (rb-y)
           l  = l + tau * dl
           write (iu, 1003) rb
           write (6, *) 'Field line tracing is stopped at the boundary'
           exit
        endif

        l = l + dl
        y = output_coordinates()
        write (iu, 1003) y
     enddo
     write (iu, *)
     write (6, 3000) l
  enddo field_line_loop
  close (iu)

 1001 format (8x,'Number of trace steps: ',i8)
 1002 format (8x,'Initial position:')
 1003 format (3(e25.18,1x))
 1004 format (8x,3(e18.10,1x))
 2000 format ('# Flux coordinates: Theta[rad], PsiN, Phi[rad]')
 2001 format ('# Cartesian coordinates: x[cm], y[cm], z[cm]')
 2002 format ('# Cylindrical coordinates: R[cm], Z[cm], Phi[rad]')
 3000 format (8x,'trace length: ',f12.4)
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
