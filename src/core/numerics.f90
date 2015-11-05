module numerics
  use iso_fortran_env
  use math
  implicit none


  integer, parameter :: &
      EULER                   = 1, &
      RUNGE_KUTTA_4           = 2, &
      ADAMS_BASHFORTH_4       = 3, &
      ADAMS_BASHFORTH_MOULTON = 4, &
      DLSODE                  = 5

!  integer, parameter :: FL_LINE  = 1
!  integer, parameter :: FL_ARC   = 2
!  integer, parameter :: FL_ANGLE = 3

  real(real64) :: &
     Trace_Step     = 1.d0       ! step size for field line tracing [cm or deg]


  integer :: &
     Trace_Method   = ADAMS_BASHFORTH_4, &          ! Integrator used for field line tracing (see module fieldline)
     Trace_Coords   = CYLINDRICAL          ! Coordinate system for field line tracing (see module fieldline)


  contains
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine load_numerics()
  use parallel

  integer, parameter :: iu = 23

  namelist /NumericsControl/ Trace_Step, Trace_Method, Trace_Coords


  ! load numerical parameters on first processor
  if (firstP) then
     open  (iu, file='run_input')
     read  (iu, NumericsControl, end=1000)
 1000 continue
     close (iu)


     ! convert units: deg -> rad
     if (Trace_Coords == 3) then
        Trace_Step = Trace_Step / 180.d0 * pi
     endif
  endif


  ! broadcast data to other processors
  call wait_pe()
  call broadcast_real_s (Trace_Step      )
  call broadcast_inte_s (Trace_Coords    )
  call broadcast_inte_s (Trace_Method    )

  end subroutine load_numerics
  !---------------------------------------------------------------------

end module numerics
