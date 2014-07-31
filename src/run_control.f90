module run_control
  use parallel
  implicit none
  include '../config.h'

  ! user defined variables
  character*120 :: &
     Machine        = ' ', &        ! select input directory (1st part)
     Magnetic_Setup = ' ', &        ! select input directory (2nd part)
     Run_Type       = ' ', &        ! select sub-program to execute
     Output_File    = '', &
     Grid_File      = ''

  real*8 :: &
     x_start(3)     = 0.d0, &       ! initial position for field line tracing
     Trace_Step     = 1.d0          ! step size for field line tracing


  integer :: &
     N_steps        = 1000, &       ! Number of discrete steps
     Trace_Method   = 3, &          ! Method for field line tracing (see module fieldline)
     Trace_Coords   = 2, &          ! Coordinate system for field line tracing (see module fieldline)
     Output_Format  = 1             ! See individual tools



  ! internal variables
  character*120 :: Prefix, &
                   Bfield_input_file, &
                   PFC_input_file, &
                   PFC_sub_dir


  namelist /RunControl/ &
     Machine, Magnetic_Setup, &
     Run_Type, Output_File, Grid_File, Output_Format, &
     x_start, Trace_Step, Trace_Method, Trace_Coords, N_steps

  contains
!=======================================================================


!=======================================================================
  subroutine load_run_control()

  integer, parameter :: iu = 23
  character*255      :: homedir


  ! load run control on first processor
  if (mype == 0) then
     open  (iu, file='run_input', err=5000)
     read  (iu, RunControl, end=5000)
     close (iu)

     if (Machine .ne. ' ') then
        write (6, *) 'Machine:                ', trim(Machine)
        write (6, *) 'Magnetic configuration: ', trim(Magnetic_Setup)
        call getenv("HOME", homedir)
        Prefix = trim(homedir)//'/'//base_dir//'/'//trim(Machine)//'/'// &
                 trim(Magnetic_Setup)//'/'
     else
        Prefix = './'
     endif

     Bfield_input_file = trim(Prefix)//'bfield.conf'
     PFC_input_file    = 'pfc.conf'
     PFC_sub_dir       = 'pfc'
  endif


  ! broadcase data to other processors
  call wait_pe()
  !call broadcast_char   (Run_Type, 120)

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
  case default
     write (6, *) 'run type "', trim(Run_Type), '" not defined!'
     stop
  end select

 1000 format (/ '========================================================================')
  end subroutine run_control_main
!=======================================================================

end module run_control
