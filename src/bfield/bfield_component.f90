!===============================================================================
! Template definition for magnetic field components
!===============================================================================
module bfield_component
  use iso_fortran_env
  implicit none
  private


  type, public :: t_bfield
     ! boundaries for domain of definition
     real(real64) :: Rmin = 0.d0, Rmax = huge(1.d0), Zmin = -huge(1.d0), Zmax = - huge(1.d0)

     contains
     procedure :: load
     procedure :: broadcast
     procedure :: out_of_bounds
     procedure :: get_Bf
     procedure :: get_JBf
     procedure :: get_Psi
     procedure :: get_DPsi
     procedure :: cleanup
  end type t_bfield

  contains
  !=====================================================================


  !=====================================================================
  ! load configuration and setup related variables
  !=====================================================================
  subroutine load(this, filename, ierr)
  class(t_bfield)               :: this
  character(len=*), intent(in)  :: filename
  integer,          intent(out) :: ierr

  end subroutine load
  !=====================================================================



  !=====================================================================
  ! broadcast data for parallel execution
  !=====================================================================
  subroutine broadcast(this)
  class(t_bfield) :: this

  end subroutine broadcast
  !=====================================================================



  !=====================================================================
  function out_of_bounds(this, r)
  use exceptions, only: OUT_OF_BOUNDS_FLAG => OUT_OF_BOUNDS
  class(t_bfield)                        :: this
  real(real64), dimension(:), intent(in) :: r
  logical                                :: out_of_bounds


  if (size(r) < 2) then
     write (6, 9000) size(r)
     stop
  endif
 9000 format('error: invalid call to out_of_bounds with size(r) = ',i0)


  out_of_bounds = .false.
  if (r(1) < this%Rmin  .or.  r(1) > this%Rmax  .or. &
      r(2) < this%Zmin  .or.  r(2) > this%Zmax) then
     out_of_bounds = .true.
     OUT_OF_BOUNDS_FLAG = .true.
  endif

  end function out_of_bounds
  !=====================================================================



  !=====================================================================
  ! return magnetic field components (BR, BZ, Bphi) [T] at location
  ! (R[cm], Z[cm], phi[rad])
  !=====================================================================
  function get_Bf(this, r) result(Bf)
  class(t_bfield)          :: this
  real(real64), intent(in) :: r(3)
  real(real64)             :: Bf(3)


  Bf = 0.d0

  end function get_Bf
  !=====================================================================



  !=====================================================================
  ! return Jacobian of magnetic field [T/m] in cylindrical coordinates
  !=====================================================================
  function get_JBf(this, r) result(JBf)
  class(t_bfield)          :: this
  real(real64), intent(in) :: r(3)
  real(real64)             :: JBf(3,3)


  JBf = 0.d0

  end function get_JBf
  !=====================================================================



  !=====================================================================
  ! return poloidal flux [Wb]
  !=====================================================================
  function get_Psi(this, r)
  class(t_bfield)          :: this
  real(real64), intent(in) :: r(3)
  real(real64)             :: get_Psi


  get_Psi = 0.d0

  end function get_Psi
  !=====================================================================



  !=====================================================================
  ! return (mR,mZ)-th derivative of poloidal flux [Wb / cm**(mR+mZ)]
  !=====================================================================
  function get_DPsi(this, r, mR, mZ)
  class(t_bfield)          :: this
  real(real64), intent(in) :: r(3)
  integer,      intent(in) :: mR, mZ
  real(real64)             :: get_DPsi


  get_DPsi = 0.d0

  end function get_DPsi
  !=====================================================================



  !=====================================================================
  subroutine cleanup(this)
  class(t_bfield) :: this

  end subroutine cleanup
  !=====================================================================

end module bfield_component
