module run_control
  use parallel
  implicit none
  include '../config.h'

  integer, parameter :: &
     INERVOUS   =  1, &
     IMODERATE  = 10, &
     IDONTPANIC = 42

  ! user defined variables
  character*120 :: &
     Machine        = ' ', &        ! select input directory (1st part)
     Configuration  = ' ', &        ! select input directory (2nd part)
     Run_Type       = ' ', &        ! select sub-program to execute
     Output_File    = '', &
     Grid_File      = ''

  real*8 :: &
     x_start(3)     = 0.d0, &       ! initial position for field line tracing
     Trace_Step     = 1.d0, &       ! step size for field line tracing
     Limit          = 2.d4, &       ! maximum distance for field line tracing (in one direction)
     R_start        = 0.d0, &       ! radial start- and
     R_end          = 0.d0, &       ! end position for Poincare plots
     Phi_output     = 0.d0          ! Reference plane for Poincare plots


  integer :: &
     N_steps        = 1000, &       ! Number of discrete steps
     N_points       = 0, &          ! Max. number of points for Poincare plots
     N_sym          = 1, &          ! Toroidal symmetry factor (for Poincare plots)
     N_mult         = 1, &          !
     Trace_Method   = 3, &          ! Method for field line tracing (see module fieldline)
     Trace_Coords   = 2, &          ! Coordinate system for field line tracing (see module fieldline)
     Output_Format  = 1, &          ! See individual tools
     Panic_Level    = IMODERATE



  ! internal variables
  character*120 :: Prefix, &
                   Bfield_input_file


  namelist /RunControl/ &
     Machine, Configuration, &
     Run_Type, Output_File, Grid_File, Output_Format, Panic_Level, &
     x_start, Trace_Step, Trace_Method, Trace_Coords, N_steps, Limit, &
     R_start, R_end, Phi_output, N_points, N_sym, N_mult

  contains
!=======================================================================


!=======================================================================
  subroutine load_run_control()

  integer, parameter :: iu = 23
  character*255      :: homedir


  ! load run control on first processor
  if (firstP) then
     open  (iu, file='run_input', err=5000)
     read  (iu, RunControl, end=5000)
     close (iu)

     if (Machine .ne. ' ') then
        write (6, *) 'Machine:                ', trim(Machine)
        write (6, *) 'Configuration:          ', trim(Configuration)
        call getenv("HOME", homedir)
        Prefix = trim(homedir)//'/'//base_dir//'/'//trim(Machine)//'/'// &
                 trim(Configuration)//'/'
     else
        Prefix = './'
     endif

     Bfield_input_file = trim(Prefix)//'bfield.conf'
  endif


  ! broadcase data to other processors
  call wait_pe()
  call broadcast_char   (Run_Type   , 120)
  call broadcast_char   (Grid_File  , 120)
  call broadcast_char   (Output_File, 120)
  call broadcast_real   (x_start    ,   3)
  call broadcast_real_s (Trace_Step      )
  call broadcast_real_s (Limit           )
  call broadcast_real_s (R_start         )
  call broadcast_real_s (R_end           )
  call broadcast_real_s (Phi_output      )
  call broadcast_inte_s (N_steps         )
  call broadcast_inte_s (N_points        )
  call broadcast_inte_s (N_sym           )
  call broadcast_inte_s (N_mult          )
  call broadcast_inte_s (Trace_Method    )
  call broadcast_inte_s (Trace_Coords    )
  call broadcast_inte_s (Output_Format   )

  return
 5000 write  (6,5001)
 5001 format ('error reading control table from input file')
  stop
  end subroutine load_run_control
!=======================================================================



!=======================================================================
  subroutine run_control_main()

  if (firstP) then
     write (6, 1000)
  endif


  select case (Run_Type)
  case ('sample_bfield')
     call sample_bfield
  case ('trace_bline')
     call trace_bline
  case ('poincare_plot')
     call poincare_plot
  case ('connection_length')
     call connection_length
  case ('get_equi_info_2D')
     call get_equi_info_2D
  case default
     write (6, *) 'run type "', trim(Run_Type), '" not defined!'
     stop
  end select

 1000 format (/ '========================================================================')
  end subroutine run_control_main
!=======================================================================

end module run_control
