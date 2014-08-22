!===============================================================================
! User defined magnetic axis
!===============================================================================
module magnetic_axis
  use iso_fortran_env
  use bspline
  implicit none

  private

  ! for axisymmetric (2D) configuration
  real(real64) :: &
     R_axis = -1.d0, &
     Z_axis =  0.d0

  ! for non-axisymmetric (3D) configuration
  integer, parameter :: nord = 5
  real(real64), dimension(:), allocatable :: Phinot, Rcoeff, Zcoeff
  integer :: N_phi, N_sym


  procedure(magnetic_axis_default), pointer :: get_magnetic_axis => magnetic_axis_default

  public :: &
     setup_magnetic_axis_2D, &
     load_magnetic_axis_3D, &
     get_magnetic_axis

  contains
!=======================================================================



!=======================================================================
! Default function for magnetic axis
!=======================================================================
  function magnetic_axis_default(phi) result(r)
  use run_control
  real(real64), intent(in) :: phi
  real(real64)             :: r(3)


  r = 0.d0
  if (Panic_Level < IMODERATE) then
     write (6, *) 'error: magnetic axis undefined!'
     stop
  endif

  end function magnetic_axis_default
!=======================================================================



!=======================================================================
! Set position (R [cm], Z [cm]) of magnetic axis
!=======================================================================
  subroutine setup_magnetic_axis_2D (R, Z)
  real(real64), intent(in) :: R, Z

  R_axis = R
  Z_axis = Z
  get_magnetic_axis => magnetic_axis_2D

  return
  end subroutine setup_magnetic_axis_2D
!=======================================================================



!=======================================================================
! Provide position r=(R,Z,phi) [cm,cm,rad] of magnetic axis
!=======================================================================
  function magnetic_axis_2D(phi) result(r)
  real(real64), intent(in) :: phi
  real(real64)             :: r(3)


  r(1) = R_axis
  r(2) = Z_axis
  r(3) = phi

  end function magnetic_axis_2D
!=======================================================================



!=======================================================================
  subroutine load_magnetic_axis_3D
  use math

  integer, parameter :: iu = 54

  real(real64), dimension(:), allocatable :: Phitmp, Rtmp, Ztmp
  integer :: i


  open  (iu, file='axis.dat')
  read  (iu, *) N_phi, N_sym
  N_phi = N_phi + 1
  allocate (Phitmp(N_phi), Rtmp(N_phi), Ztmp(N_phi))
  allocate (Phinot(N_phi+nord))

  do i=1,N_phi-1
     read (iu, *) Rtmp(i), Ztmp(i)
  enddo
  Rtmp(N_phi) = Rtmp(1)
  Ztmp(N_phi) = Ztmp(1)
  close (iu)

  do i=1,N_phi
     Phitmp(i) = 1.d0 * (i-1) / (N_phi-1) / N_sym * pi2
  enddo
  call dbsnak (N_phi, Phitmp, nord, Phinot)

  allocate (Rcoeff(N_phi), Zcoeff(N_phi))
  call dbsint (N_phi, Phitmp, Rtmp, nord, Phinot, Rcoeff)
  call dbsint (N_phi, Phitmp, Ztmp, nord, Phinot, Zcoeff)
  deallocate (Phitmp, Rtmp, Ztmp)

  get_magnetic_axis => magnetic_axis_3D

  end subroutine load_magnetic_axis_3D
!=======================================================================



!=======================================================================
! Provide position r=(R,Z,phi) [cm,cm,rad] of magnetic axis
!=======================================================================
  function magnetic_axis_3D(phi) result(r)
  use math
  real(real64), intent(in) :: phi
  real(real64)             :: r(3)


  r(3) = phi_sym(phi, N_sym)
  r(1) = dbsval (r(3), nord, Phinot, N_phi, Rcoeff)
  r(2) = dbsval (r(3), nord, Phinot, N_phi, Zcoeff)

  end function magnetic_axis_3D
!=======================================================================

end module magnetic_axis
