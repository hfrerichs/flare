!===============================================================================
! Field line reconstruction analysis
!
! Input (taken from run control file):
!    Output_File        Configuration file for field line analysis
!===============================================================================
subroutine FLR_analysis
  use iso_fortran_env
  use run_control, only: Target_File => Output_File
  use parallel
  use flux_tube
  use dataset
  use math
  implicit none

  integer, parameter :: iu = 32

  character(len=72)  :: &
     Operation     = 'nothing', &
     Output_File   = 'output.dat', &
     Grid_File     = 'output.grid', &
     Cross_Section = ''

  real(real64)       :: x0(2), phi0, a1, alpha1, a2, alpha2, P, theta
  real(real64)       :: xi = 0.d0, eta = 0.d0
  integer            :: nphi, nlength

  namelist /Analysis_Input/ &
     Operation, Output_File, Grid_File, Cross_Section, &
     x0, phi0, a1, alpha1, a2, alpha2, P, theta, nphi, nlength, &
     xi, eta

  type(t_flux_tube)  :: FT


  ! initialize field line reconstruction analysis
  if (firstP) then
     write (6, *) 'Field line reconstruction analysis, target file: ', trim(Target_File)
     write (6, *)
  else
     return
  endif
  open  (iu, file=Target_File)
  read  (iu, Analysis_Input)
  close (iu)


  select case(Operation)
  case ('single_fieldline')
     call single_fieldline
  case ('fluxtube')
     call fluxtube
  case default
     write (6, *) 'error: invalid operation ', trim(Operation), '!'
     stop
  end select

  return
  contains
!=======================================================================



!-----------------------------------------------------------------------
  function load_cross_section() result(cs)
  type (t_cross_section) :: cs


  if (Cross_Section .ne. '') then
     call cs%load(Cross_Section)
     x0 = cs%a0
  else
     ! deg -> rad
     phi0  = phi0 / 180.d0 * Pi
     call cs%generate(x0, phi0, a1, alpha1, a2, alpha2, P, theta)
  endif
  write (6, 1000) cs%a0, cs%a1, cs%a2, cs%a3
  write (6, *)

 1000 format(3x,'- Shape coefficients of flux tube:',/,&
             8x,'a0 = (',f8.3,', ',f8.3,')',/, &
             8x,'a1 = (',f8.3,', ',f8.3,')',/, &
             8x,'a2 = (',f8.3,', ',f8.3,')',/, &
             8x,'a3 = (',f8.3,', ',f8.3,')')
  end function load_cross_section
!-----------------------------------------------------------------------



!=======================================================================
! Analyze field line with intrinsic coordinates (xi, eta) in flux tube
! (cs0, nphi, nlength). The initial cross section cs0 is either given in
! external file Cross_Section, or defined by geometry parameters
! (x0, phi0, a1, alpha1, a2, alpha2, P, theta)
!
! Output: nphi rows of
!     phi[deg]    err[cm]    err_rad[cm]    err_pol[cm]
!
! err = |(R,Z)_reconstructed  -  (R,Z)_traced|
! err_rad: contribution in e_Psi (radial) direction
! err_pol: contribution in poloidal direction (perpendicular to e_Psi)
!=======================================================================
  subroutine single_fieldline

  type(t_cross_section) :: cs0
  type(t_dataset) :: D
  integer         :: j


  ! get initial cross-section and generate flux tube
  cs0 = load_cross_section()
  call FT%generate(cs0, nphi, nlength)
  write (6, 1000) x0, phi0


  ! check intrinsic coordinates
  if ((abs(xi) > 1.d0)  .or.  (abs(eta) > 1.d0)) then
     write (6, *) 'error: intrinsic coordinates xi, eta must be in interval [-1,1]!'
     write (6, *) 'xi  = ', xi
     write (6, *) 'eta = ', eta
     stop
  endif


  ! analyze field line (xi, eta)
  D = FT%analyze_fieldline(xi, eta)
  open  (iu, file=Output_File)
  do j=1,nphi
     write (iu, 2000) FT%phi(j)*180.d0/pi, D%x(j,:)
  enddo
  close (iu)

 1000 format(3x,'- Analyzing single fieldline at (R,Z,phi) = ( ',f8.3,' cm, ',f8.3,' cm, ',f8.3,' deg)')
 2000 format(4e18.10)
  end subroutine single_fieldline
!=======================================================================



!=======================================================================
! Analysis of entire flux tube cross-section at final length 360 deg / nlength
!
! Output: nsample x nsample values of
!     phi[deg]    err[cm]    err_rad[cm]    err_pol[cm]
!=======================================================================
  subroutine fluxtube
  use grid
  use mesh_spacing
  type(t_cross_section) :: cs0
  type(t_dataset) :: Dfl, D
  type(t_grid)    :: G
  integer         :: i, j, ig, nsample


  ! get initial cross-section and generate flux tube
  cs0 = load_cross_section()
  call FT%generate(cs0, 1, nlength)
  write (6, 1000)


  ! sample cross-section of flux tube
  nsample = 64
  call G%new(LOCAL, UNSTRUCTURED, 3, nsample, nsample)
  call D%new(nsample*nsample, 3)
  ig = 0
  do j=0,nsample-1
     eta = -1.d0  +  2.d0 * Equidistant%node(j, nsample-1)
     do i=0,nsample-1
        xi = -1.d0  +  2.d0 * Equidistant%node(i, nsample-1)

        ig = ig + 1
        G%x(ig,1) = xi
        G%x(ig,2) = eta

        Dfl = FT%analyze_fieldline(xi, eta)
        D%x(ig,:) = Dfl%x(1,:)
     enddo
  enddo

  ! write data and grid to output files
  call G%store(Grid_File)
  call D%plot(filename=Output_File)

 1000 format(3x,'- Analyzing flux tube')
  end subroutine fluxtube
!=======================================================================

end subroutine FLR_analysis
