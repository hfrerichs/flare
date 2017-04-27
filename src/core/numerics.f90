module numerics
  use iso_fortran_env
  use math
  implicit none


  real(real64), parameter :: &
     m_to_cm = 1.d2,   cm_to_m = 1.d0 / m_to_cm


  integer, parameter :: &
      EULER                   = 1, &
      RUNGE_KUTTA_4           = 2, &
      ADAMS_BASHFORTH_4       = 3, &
      ADAMS_BASHFORTH_MOULTON = 4, &
      DLSODE                  = 5

  character(len=*), dimension(5), parameter :: &
      INTEGRATION_METHOD = (/'Euler                         ', &
                             '4th order Runge-Kutta         ', &
                             '4 step Adams-Bashforth        ', &
                             '4 step Adams-Bashforth-Moulton', &
                             'DLSODE                        '/)

!  integer, parameter :: FL_LINE  = 1
!  integer, parameter :: FL_ARC   = 2
!  integer, parameter :: FL_ANGLE = 3

  real(real64) :: &
     Trace_Step     = 1.d0       ! step size for field line tracing [cm or deg]


  integer :: &
     Trace_Method   = ADAMS_BASHFORTH_4, &   ! Integrator used for field line tracing (see module fieldline)
     Trace_Coords   = CYLINDRICAL, &         ! Coordinate system for field line tracing (see module fieldline)
     Spline_Order   = 5


  logical :: &
     OUT_OF_BOUNDS  = .false.


  contains
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine load_numerics()
  use parallel

  integer, parameter :: iu = 23

  namelist /NumericsControl/ &
     Trace_Step, Trace_Method, Trace_Coords, &
     Spline_Order


  ! load numerical parameters on first processor
  if (firstP) then
     write (6, 1000)
     open  (iu, file='run_input')
     read  (iu, NumericsControl, end=1001)
 1001 continue
     close (iu)
     write (6, 1002) Trace_Step, trim(UNITS(Trace_Coords)), trim(COORDINATES(Trace_Coords))
     write (6, 1003) trim(INTEGRATION_METHOD(Trace_Method))


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
  call broadcast_inte_s (Spline_Order    )

 1000 format(3x,'- Numerics:')
 1002 format(8x,'Field line integration: step size = ', g0.4, 1x, a, ' in ', a, ' coordinates')
 1003 format(8x,'Integration method:     ',a)
  end subroutine load_numerics
  !---------------------------------------------------------------------

end module numerics
