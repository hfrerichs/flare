!===============================================================================
! Generate homogeneous magnetic surfaces
!
! Input (taken from run control file):
!    N_sym		Toroidal symmetry number
!    R_start, R_end	Radial range in which magnetic surfaces are generated
!    N_steps            Steps in radial range (-> N_steps+1 surfaces will be generated)
!    Phi_output		Toroidal location [deg] of surface cut
!    tolerance		Tolerance for radial location of magnetic surface
!    N_points           Number of sample points for magnetic surface
!
!    Trace_Step		Size of trace step [deg]
!    Output_File
!===============================================================================
subroutine quasi_surfaces()
  use iso_fortran_env
  use run_control, only: N_sym, R_start, R_end, N_steps, Phi_output, tolerance, N_points, Trace_Step, Output_File
  use parallel
  use math
  use quasi_surface
  use cspline
  use equilibrium, only: get_PsiN, get_poloidal_angle
  implicit none

  integer, parameter    :: iu = 42

  type(t_quasi_surface) :: Q
  type(t_cspline)       :: S
  real(real64)          :: r(3), t, theta, PsiN
  integer               :: i, istep


  if (firstP) then
     write (6, *) 'Generate homogeneous magnetic surfaces, output in: ', adjustl(trim(Output_File))
     write (6, *)
  else
     return
  endif


  write (6, 1001) N_sym
  write (6, 2001) R_start, R_end, N_steps
  write (6, 3001)
  open  (iu, file=Output_File)
  write (iu, 2002)
  do istep=0,N_steps
     ! reference point on magnetic surface
     r(1) = R_start + (R_end - R_start) * istep / N_steps
     r(2) = 0.d0
     r(3) = Phi_output
     write (6, 3002) istep, r(1)

     ! generate magnetic surface at r
     call Q%generate(r, tolerance, N_sym, Trace_Step)
     write (6, 3003) Q%Rout, Q%H

     ! smooth sampling of points on magnetic surface
     call S%setup(Q%nodes, periodic=.true.)
     do i=0,N_points
        t       = 1.d0 * i / N_points
        r(1:2)  = S%eval(t)
        PsiN    = get_PsiN(r)
        theta   = get_poloidal_angle(r) / pi * 180.d0
        if (theta < 0.d0) theta = theta + 360.d0
        write (iu, *) r(1:2), theta, PsiN
     enddo
     write (iu, *)
  enddo
  close (iu)

 1001 format (8x,'Toroidal symmetry number:         ',i4)
 2001 format (8x,'Radial domain:',5x,'R_start = ',f6.2,5x, &
              'R_end = ',f6.2,5x,'with ',i4,' steps'/)
 2002 format ('# R [cm], Z [cm], Theta [deg], PsiN')
 3001 format (8x,'#       R_initial [cm]     R_surface [cm]     Homogeneity')
 3002 format (8x,i0,7x,f0.4)
 3003 format (35x,f0.4,8x,f0.4)

end subroutine quasi_surfaces
