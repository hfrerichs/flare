!===============================================================================
! Template definition for magnetic field components
!===============================================================================
module abstract_bfield
  use iso_fortran_env
  implicit none
  private


  type, abstract, public :: t_bfield
     ! boundaries for domain of definition
     real(real64) :: Rmin = 0.d0, Rmax = huge(1.d0), Zmin = -huge(1.d0), Zmax = huge(1.d0)

     contains
     procedure :: broadcast_bounds
     procedure :: out_of_bounds
     procedure(broadcast), deferred :: broadcast
     procedure(get_Bf),    deferred :: get_Bf
     procedure(get_JBf),   deferred :: get_JBf
  end type t_bfield


  ! axisymmetric equilibrium data
  type, abstract, extends(t_bfield), public :: t_equi2d
     contains
     procedure(get_Psi),   deferred :: get_Psi
     procedure(get_DPsi),  deferred :: get_DPsi
  end type t_equi2d


  abstract interface
     ! broadcast data for parallel execution
     subroutine broadcast(this)
     import t_bfield
     class(t_bfield) :: this
     end subroutine broadcast


     ! return magnetic field components (BR, BZ, Bphi) [T] at (R[m], Z[m], phi[rad])
     function get_Bf(this, r) result(Bf)
     use iso_fortran_env
     import t_bfield
     class(t_bfield)               :: this
     real(real64), intent(in)      :: r(3)
     real(real64)                  :: Bf(3)
     end function get_Bf


     ! return Jacobian of magnetic field [T/m] in cylindrical coordinates
     function get_JBf(this, r) result(JBf)
     use iso_fortran_env
     import t_bfield
     class(t_bfield)               :: this
     real(real64), intent(in)      :: r(3)
     real(real64)                  :: JBf(3,3)
     end function get_JBf


     ! return poloidal flux [Wb]
     function get_Psi(this, r)
     use iso_fortran_env
     import t_equi2d
     class(t_equi2d)               :: this
     real(real64), intent(in)      :: r(3)
     real(real64)                  :: get_Psi
     end function get_Psi


     ! return (mR,mZ)-th derivative of poloidal flux [Wb / cm**(mR+mZ)]
     function get_DPsi(this, r, mR, mZ)
     use iso_fortran_env
     import t_equi2d
     class(t_equi2d)               :: this
     real(real64), intent(in)      :: r(3)
     integer,      intent(in)      :: mR, mZ
     real(real64)                  :: get_DPsi
     end function get_DPsi
  end interface
  contains
  !=====================================================================


  !=====================================================================
  ! broadcast data for parallel execution
  !=====================================================================
  subroutine broadcast_bounds(this)
  use parallel
  class(t_bfield) :: this


  call broadcast_real_s(this%Rmin)
  call broadcast_real_s(this%Rmax)
  call broadcast_real_s(this%Zmin)
  call broadcast_real_s(this%Zmax)

  end subroutine broadcast_bounds
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

end module abstract_bfield
