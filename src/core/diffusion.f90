module diffusion
  use iso_fortran_env
  implicit none
  private

  integer, parameter, public :: &
     DIFFUSION_RZ_CIRCLE = 1, &
     DIFFUSION_RZ_GAUSS  = 2


  public :: diffusion_step, field_line_diffusion, setup_diffusion

  logical :: field_line_diffusion = .false.
  !logical :: field_line_diffusion = .true.


  real(real64) :: D = 0.1d0
  !integer      :: itype = DIFFUSION_RZ_GAUSS
  integer      :: itype = DIFFUSION_RZ_CIRCLE


  contains
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine setup_diffusion(diffusion_coefficient, diffusion_type)
  use parallel
  real(real64), intent(in) :: diffusion_coefficient
  integer,      intent(in) :: diffusion_type


  D     = diffusion_coefficient
  itype = diffusion_type
  call broadcast_real_s(D)
  call broadcast_inte_s(itype)
  if (D <= 0.d0) return


  select case(itype)
  case(DIFFUSION_RZ_CIRCLE, DIFFUSION_RZ_GAUSS)
     write (6, 1000) itype, D
  case default
     write (6, 9000) itype;   stop
  end select
 1000 format('field line diffusion ',i0,' with Dfl = ',f0.3)
 9000 format("error: diffusion type ",a," not supported!")
  field_line_diffusion = .true.

  end subroutine setup_diffusion
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  ! r:  current position of field line (R[cm], Z[cm], Phi[rad])
  ! dl: trace step (arclength) along field line
  !---------------------------------------------------------------------
  subroutine diffusion_step(r, dl)
  use math
  real(real64), intent(inout) :: r(3)
  real(real64), intent(in)    :: dl

  real(real64) :: ds0, ds(3), phi, rgauss


  ds = 0.d0
  select case(itype)
  case(DIFFUSION_RZ_CIRCLE)
     ds0    = sqrt(4.d0 * D * dl)
     phi    = rand() * pi2
     ds(1)  = ds0 * cos(phi)
     ds(2)  = ds0 * sin(phi)

  case(DIFFUSION_RZ_GAUSS)
     ds0    = sqrt(2.d0 * D * dl)
     phi    = rand() * pi2
     rgauss = min(sqrt(-2.d0 * log(rand())), 10.d0)
     ds(1)  = ds0 * rgauss * cos(phi)
     ds(2)  = ds0 * rgauss * sin(phi)

  end select
  r = r + ds

  end subroutine diffusion_step
  !---------------------------------------------------------------------

end module diffusion
