module sonnet
  use iso_fortran_env
  use abstract_bfield
  use bspline2D
  implicit none
  private


  type, extends(t_equi2d), public :: t_sonnet
     ! no. of grid points in radial and vertical direction
     integer :: jm, km

     ! radial and vertical coordinates of grid points [m]
     real(real64), dimension(:), allocatable :: R, Z

     ! flux per radiant at grid points and at plasma boundary [Wb/rad]
     real(real64), dimension(:,:), allocatable :: psi
     real(real64) :: psib

     ! toroidal magnetic field [T] at reference location
     real(real64) :: btf
     ! major radius of reference location for btf
     real(real64) :: rtf

     ! 2D spline interpolation of flux
     type(t_bspline2D) :: S

     contains
     procedure :: load
     procedure :: broadcast
     procedure :: get_Bf
     procedure :: get_JBf
     procedure :: get_Psi
     procedure :: get_DPsi
     procedure :: cleanup
     procedure :: info
  end type t_sonnet

  contains
  !=====================================================================


  !=====================================================================
  ! load configuration and setup related variables
  !=====================================================================
  subroutine load(this, source, ierr)
  use numerics, only: m_to_cm
  class(t_sonnet)               :: this
  character(len=*), intent(in)  :: source
  integer,          intent(out) :: ierr

  integer, parameter :: iu = 17

  character(len=120) :: tmp
  logical :: ex
  integer :: i, jm, km


  ierr = 0
  inquire (file=source, exist=ex)
  if (.not.ex) then
     write (6, 9000) trim(source)
     ierr = 1
     return
  endif
 9000 format('error: file ',a,' does not exist!')


  open  (iu, file=source)

  ! skip past header lines
  do i=1,10
     read (iu, 1000) tmp
  enddo

  ! read grid resolution
  read  (iu, 1000) tmp;  read (tmp(18:26), *) jm;  this%jm = jm
  read  (iu, 1000) tmp;  read (tmp(18:26), *) km;  this%km = km
  write (6, 3001) jm, km


  ! allocate memory for data
  allocate (this%R(jm), this%Z(km), this%psi(jm,km))


  ! read characteristic parameters
  read  (iu, 1000) tmp;  read (tmp(14:40), *) this%psib
  read  (iu, 1000) tmp;  read (tmp(14:40), *) this%btf
  read  (iu, 1000) tmp;  read (tmp(14:40), *) this%rtf
  write (6, 3003) this%btf
  write (6, 3004) this%rtf * m_to_cm


  ! read grid nodes
  do i=1,2;  read (iu, 1000) tmp;  enddo
  read  (iu, *) this%R
  do i=1,2;  read (iu, 1000) tmp;  enddo
  read  (iu, *) this%Z
  ! set up boundaries
  this%Rmin = minval(this%R)*m_to_cm;  this%Rmax = maxval(this%R)*m_to_cm
  this%Zmin = minval(this%Z)*m_to_cm;  this%Zmax = maxval(this%Z)*m_to_cm


  ! read flux
  do i=1,2;  read (iu, 1000) tmp;  enddo
  read  (iu, *) this%psi;  this%psi = this%psi + this%psib
  close (iu)


  ! set up 2D spline interpolation
  call this%S%init(jm, km, this%R, this%Z, this%psi)


1000 format (a120)
3001 format (8x,'Grid resolution: ',19x,'R : ',i4,' points;   Z : ',i4,' points')
3003 format (8x,'Vacuum toroidal magnetic field in Tesla at R0:      ',2x,f8.3)
3004 format (8x,'Reference position R0 in cm:                        ',2x,f8.3)
  end subroutine load
  !=====================================================================



  !=====================================================================
  ! broadcast data for parallel execution
  !=====================================================================
  subroutine broadcast(this)
  use parallel
  class(t_sonnet) :: this


  call broadcast_inte_s(this%jm)
  call broadcast_inte_s(this%km)
  call broadcast_real_s(this%psib)
  call broadcast_real_s(this%btf)
  call broadcast_real_s(this%rtf)

  if (mype > 0) then
      allocate (this%R(this%jm), this%Z(this%km), this%psi(this%jm,this%km))
  endif
  call broadcast_real  (this%R, this%jm)
  call broadcast_real  (this%Z, this%km)
  call broadcast_real  (this%psi, this%jm*this%km)
  call this%S%broadcast()

  end subroutine broadcast
  !=====================================================================



  !=====================================================================
  ! return magnetic field components (BR, BZ, Bphi) [T] at location
  ! (R[cm], Z[cm], phi[rad])
  !=====================================================================
  function get_Bf(this, r) result(Bf)
  use numerics, only: cm_to_m
  class(t_sonnet)          :: this
  real(real64), intent(in) :: r(3)
  real(real64)             :: Bf(3)

  real(real64) :: dpsidr, dpsidz


  ! poloidal field
  dpsidr = this%get_DPsi(r, 1, 0)
  dpsidz = this%get_DPsi(r, 0, 1)
  Bf(1)  = -dpsidz / r(1) / cm_to_m**2
  Bf(2)  =  dpsidr / r(1) / cm_to_m**2


  ! toroidal field
  Bf(3)   =  this%btf * this%rtf / (r(1) * cm_to_m)

  end function get_Bf
  !=====================================================================



  !=====================================================================
  ! return Jacobian of magnetic field [T/m] in cylindrical coordinates
  !=====================================================================
  function get_JBf(this, r) result(JBf)
  class(t_sonnet)          :: this
  real(real64), intent(in) :: r(3)
  real(real64)             :: JBf(3,3)


  ! NOT IMPLEMENTED YET
  JBf = 0.d0

  end function get_JBf
  !=====================================================================



  !=====================================================================
  ! return poloidal flux [Wb]
  !=====================================================================
  function get_Psi(this, r)
  class(t_sonnet)          :: this
  real(real64), intent(in) :: r(3)
  real(real64)             :: get_Psi


  get_Psi = this%get_DPsi(r, 0, 0)

  end function get_Psi
  !=====================================================================



  !=====================================================================
  ! return (mR,mZ)-th derivative of poloidal flux  [Wb / cm**(mR+mZ)]
  !=====================================================================
  function get_DPsi(this, r, mR, mZ)
  use numerics, only: m_to_cm, cm_to_m
  class(t_sonnet)          :: this
  real(real64), intent(in) :: r(3)
  integer,      intent(in) :: mR, mZ
  real(real64)             :: get_DPsi

  real(real64) :: x(2)


  get_DPsi = -1.d0
  if (this%out_of_bounds(r)) return

  x        = r(1:2) * cm_to_m
  get_DPsi = this%S%derivative(x, mR, mZ) / m_to_cm**(mR+mZ)

  end function get_DPsi
  !=====================================================================



  !=====================================================================
  subroutine cleanup(this)
  class(t_sonnet) :: this


  if (allocated(this%R)) deallocate(this%R)
  if (allocated(this%Z)) deallocate(this%Z)
  if (allocated(this%psi)) deallocate(this%psi)

  end subroutine cleanup
  !=====================================================================



  !=====================================================================
  subroutine info(this, R0)
  use numerics, only: m_to_cm
  class(t_sonnet)           :: this
  real(real64), intent(out) :: R0


  R0 = this%rtf * m_to_cm

  end subroutine info
  !=====================================================================

end module sonnet
