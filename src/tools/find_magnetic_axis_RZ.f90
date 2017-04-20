!===============================================================================
! Find magnetic axis RZ in phi_start = x_start(3) plane
!
! Input (taken from run control file):
!    N_sym              Toroidal symmetry number
!    x_start            Initial guess position (cylindrical coordinates) from 
!                       which the magnetic axis will be approximated
!
! optional:
!    Trace_Step
!    Trace_Method
!===============================================================================
subroutine find_magnetic_axis_RZ
  use iso_fortran_env
  use run_control, only: x_start, N_sym, Output_File, N_steps, &
                         Trace_Step, Trace_Method, Trace_Coords
  use math
  use parallel
  use fieldline
  implicit none

  integer, parameter :: iu = 54

  type(t_fieldline)  :: F
  real(real64)       :: x0(3), dR, dZ, yout(3), Dphi, dl, RZ0(3), xtol, xtol_want
  integer            :: nr, nz, Nrefine_max, is, Ntransit, imin(2)

  Output_File = 'magnetic_axis_RZ.dat'
  if (firstP) then
    write (6, *) 'Finding magnetic axis RZ, output in: ', adjustl(trim(Output_File))
    write (6, *)
  else
    return
  endif  

  ! Field line tracing parameters
  Trace_Coords = FL_ANGLE
  Trace_Step = Trace_Step*pi/180.d0
  Dphi = 2.d0*pi/Real(N_sym,real64)
  N_steps = Nint(Abs(Dphi/Trace_Step))
  Trace_step = Dphi/Real(N_steps,real64)
  Ntransit = 10 ! distance is average of these punctures from start

  ! Search parameters. Initial grid size and resolution, number of refinement steps
  nr = 31  ! Best to keep these odd numbers
  nz = 31
  Nrefine_max = 10
  dR = 10.d0
  dZ = 10.d0
  xtol_want = 1.d-5
  
  x0 = x_start
  write(*,*) 'Initial guess',x0
  Write(*,*) 'Using max', Nrefine_max, ' refinement steps.'
  Write(*,*) '  Using next dr,dz',dr,dz
  Do is = 1,Nrefine_max
    RZ0 = search_grid()    
    Write(*,*) '  Step ',is,', RZ axis= ',RZ0(1:2), ' dist = ',RZ0(3)
    xtol = RZ0(3)
    If (abs(xtol) < xtol_want) Exit
    dR = Max(2.d0*Abs(RZ0(1) - x0(1)),xtol_want)
    dZ = Max(2.d0*Abs(RZ0(2) - x0(2)),xtol_want)
    if (imin(1) .eq. nr) dR = dR*2.d0
    if (imin(2) .eq. nz) dZ = dZ*2.d0
    
    x0(1:2) = RZ0(1:2)
  Enddo
  If (is .eq. Nrefine_max + 1) Write(*,*) 'Reached max number of refinement steps'
  
  open  (iu, file=Output_File)
  write (iu, *) RZ0(1),RZ0(2),x_start(3)
  close (iu)
  
  Contains
!.......................................................................
  Function search_grid () Result (RZdist)
    Real(real64) :: RZdist(3)
    integer      :: i, ir, iz, it
    real(real64), Allocatable :: dist(:,:)
    

    Allocate(dist(nr,nz))
    dist(:,:) = 0.d0
    do ir = 1,nr
      do iz = 1,nz
        x_start(1) = x0(1) - dR/2.d0 + Real(ir-1,real64)/Real(nr-1,real64)*dR
        x_start(2) = x0(2) - dZ/2.d0 + Real(iz-1,real64)/Real(nz-1,real64)*dZ
        call F%init(x_start, Trace_Step, Trace_Method, Trace_Coords)
        Do it = 1,Ntransit
          Do i = 1,N_steps
            dl = F%Trace_1step()
          Enddo
          yout = F%rc
          dist(ir,iz) = dist(ir,iz) + sqrt( (yout(1)-x_start(1))**2 + (yout(2)-x_start(2))**2)
        Enddo
      enddo
    enddo
    dist = dist/Real(Ntransit,real64)
    imin = Minloc(dist)

    RZdist(1) = x0(1) - dR/2.d0 + Real(imin(1)-1,real64)/Real(nr-1,real64)*dR
    RZdist(2) = x0(2) - dZ/2.d0 + Real(imin(2)-1,real64)/Real(nz-1,real64)*dZ
    RZdist(3) = dist(imin(1),imin(2))
    Deallocate(dist)
  End Function search_grid
!.......................................................................  
end subroutine find_magnetic_axis_RZ

