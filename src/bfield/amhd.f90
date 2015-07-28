! "One size fits all" analytical solutions to the Grad-Shafranov equiation, A.J. Cerfon et al., Physics of Plasmas 17, 032502 (2010)
module amhd
  use iso_fortran_env
  use fgsl
  implicit none

  private


  ! toroidal magnetic field and plasma current
  real(real64) :: Ip, Bt, R0

  real(real64) :: PSI0_scale


  ! equilibrium shape coefficients (input)
  real(real64) :: &
     eps          = 0.4d0, & ! inverse aspect ratio = normlized midplane minor radius
     del          = 0.4d0, & ! triangularity
     kap          = 1.8d0, & ! elongation
     A            = 0.d0, &
     xsep         = 0.d0, &  ! position of lower X-point
     ysep         = 0.d0, &
     scale_manual = 1.d0
  logical      :: up_down_symmetry = .true.

  namelist /AMHD_Input/ &
     eps, del, kap, A, xsep, ysep, scale_manual, up_down_symmetry


  ! derived parameters (from equilibrium shape coefficients)
  real(real64), dimension(:), allocatable :: ceq
  integer :: n


  public :: &
     amhd_load, &
     amhd_get_Bf, &
     amhd_get_Psi, &
     amhd_get_DPsi, &
     amhd_get_pressure, &
     amhd_get_domain, &
     amhd_broadcast

  contains
!===============================================================================



!===============================================================================
  subroutine amhd_load(iu, iconfig, Ip_, Bt_, R0_)
  use magnetic_axis
  integer,      intent(in)  :: iu
  integer,      intent(out) :: iconfig
  real(real64), intent(in)  :: Ip_, Bt_, R0_

  integer(fgsl_size_t) :: npass

  ! set up toroidal magnetic field and plasma current
  Ip = Ip_
  Bt = Bt_
  R0 = R0_
  write (6, 1000)
  write (6, 1001) Bt
  write (6, 1002) R0
  write (6, 1003) Ip
  write (6, *)
  Bt_sign = int(sign(1.d0, Bt))
  Ip_sign = int(sign(1.d0, Ip))


  ! load configuration file
  rewind(iu)
  read  (iu, AMHD_Input, end=9000)
  iconfig = 1
  write (6, 2001) kap
  write (6, 2002) del
  write (6, 2003) eps
  write (6, *)


  ! set up derived coefficients
  n     = 12
  if (up_down_symmetry) n = 7
  npass = n
  call setup_amhd(npass)


  return
 1000 format(3x,'- Analytic MHD equilibrium')
 1001 format(8x,'Toroidal magnetic field (T): ',f8.3)
 1002 format(8x,'Reference Major Radius (cm): ',f8.3)
 1003 format(8x,'Plasma current (MA):         ',f8.3)
 2001 format(8x,'Elongation (kappa):          ',f8.3)
 2002 format(8x,'Triangularity (delta):       ',f8.3)
 2003 format(8x,'1 / Aspect ratio (epsilon):  ',f8.3)
 9000 iconfig = 0
  end subroutine amhd_load
!===============================================================================



!===============================================================================
  subroutine setup_amhd(n)
  integer(fgsl_size_t), intent(in) :: n

  real(fgsl_double), target :: M(n, n), b(n), x(n), vpsi(0:12)
  type(fgsl_matrix)         :: Mfgsl
  type(fgsl_vector)         :: bfgsl, xfgsl
  integer(fgsl_int)         :: stat, sig
  type(fgsl_permutation)    :: p

  real(real64)              :: alp, N1, N2, N3, r(3), Bf(3), Bpol


  ! set default position of X-point
  if (xsep == 0.d0) xsep = 1.d0 - 1.1d0*del*eps
  if (ysep == 0.d0) ysep = -1.1d0*kap*eps


  ! derived parameters
  alp = asin(del)
  N1  = - (1.d0 + alp)**2 / eps / kap**2
  N2  =   (1.d0 - alp)**2 / eps / kap**2
  N3  = - kap / eps / cos(alp)**2


  ! setup linear system
  ! 1. outer equatorial point
  vpsi   = Psi_osfa(1.d0 + eps, 0.d0)
  M(:,1) = vpsi(1:n);  b(1)   = -vpsi(0)
  ! 2. inner equatorial point
  vpsi   = Psi_osfa(1.d0 - eps, 0.d0)
  M(:,2) = vpsi(1:n);  b(2)   = -vpsi(0)
  ! 3. X-point
  vpsi   = Psi_osfa(xsep, ysep)
  M(:,3) = vpsi(1:n);  b(3)   = -vpsi(0)
  ! 4. BZ=0 at X-point
  vpsi   = Psix_osfa(xsep, ysep)
  M(:,4) = vpsi(1:n);  b(4)   = -vpsi(0)
  ! 5. outer equatorial point curvature
  vpsi   = Psiyy_osfa(1.d0 + eps, 0.d0)
  M(:,5) = vpsi(1:n);  b(5)   = -vpsi(0)
  vpsi   = Psix_osfa(1.d0 + eps, 0.d0)
  M(:,5) = M(:,5) + N1*vpsi(1:n)
  b(5)   = b(5)   - N1*vpsi(0)
  ! 6. inner equatorial point curvature
  vpsi   = Psiyy_osfa(1.d0 - eps, 0.d0)
  M(:,6) = vpsi(1:n);  b(6)   = -vpsi(0)
  vpsi   = Psix_osfa(1.d0 - eps, 0.d0)
  M(:,6) = M(:,6) + N2*vpsi(1:n)
  b(6)   = b(6)   - N2*vpsi(0)
  ! 7. BR=0 at X-point
  vpsi   = Psiy_osfa(xsep, ysep)
  M(:,7) = vpsi(1:n);  b(7)   = -vpsi(0)
  if (.not.up_down_symmetry) then
     ! 8. high point
     vpsi = Psi_osfa(1.d0 - del*eps, kap*eps)
     M(:,8) = vpsi(1:n);  b(8)   = -vpsi(0)
     ! 9. high point maximum
     vpsi = Psix_osfa(1.d0 - del*eps, kap*eps)
     M(:,9) = vpsi(1:n);  b(9)   = -vpsi(0)
     ! 10. high point curvature
     vpsi = Psixx_osfa(1.d0 - del*eps, kap*eps)
     M(:,10) = vpsi(1:n);  b(10)   = -vpsi(0)
     vpsi = Psiy_osfa(1.d0 - del*eps, kap*eps)
     M(:,10) = M(:,10) + N3*vpsi(1:n)
     b(10)   = b(10)   - N3*vpsi(0)
     ! 11. outer equatorial point: up-down symmetry
     vpsi = Psiy_osfa(1.d0 + eps, 0.d0)
     M(:,11) = vpsi(1:n);  b(11)   = -vpsi(0)
     ! 12. inner equatorial point: up-down symmetry
     vpsi = Psiy_osfa(1.d0 - eps, 0.d0)
     M(:,12) = vpsi(1:n);  b(12)   = -vpsi(0)
  endif


  ! solve linear system
  Mfgsl  = fgsl_matrix_init(type=1.0_fgsl_double)
  bfgsl  = fgsl_vector_init(type=1.0_fgsl_double)
  xfgsl  = fgsl_vector_init(type=1.0_fgsl_double)
  p      = fgsl_permutation_alloc(n)
  stat   = fgsl_matrix_align(M, n, n, n, Mfgsl)
  stat   = fgsl_vector_align(b, n, bfgsl, n, 0_fgsl_size_t, 1_fgsl_size_t)
  stat   = fgsl_vector_align(x, n, xfgsl, n, 0_fgsl_size_t, 1_fgsl_size_t)

  stat   = fgsl_linalg_LU_decomp(Mfgsl, p, sig)
  stat   = fgsl_linalg_LU_solve(Mfgsl, p, bfgsl, xfgsl)

  ! gather results
  allocate(ceq(n))
  ceq    = x


  ! cleanup
  call fgsl_matrix_free(Mfgsl)
  call fgsl_vector_free(bfgsl)
  call fgsl_vector_free(xfgsl)
  call fgsl_permutation_free(p)


  ! calculate scaling factor
  PSI0_scale = 1.d0
  r(1) = R0 * (1.d0-eps)
  r(2) = 0.d0
  r(3) = 0.d0
  Bf   = amhd_get_Bf(r)
  Bpol = Ip * 2.d5 / R0 / eps
  PSI0_scale = Bpol / Bf(2) * scale_manual

  end subroutine setup_amhd
!===============================================================================



!===============================================================================
  function Psi_osfa(x, y)
  real(real64), intent(in) :: x, y
  real(real64)             :: Psi_osfa(0:12)


  Psi_osfa(0) = x**4 / 8.d0 + A * (0.5d0*x**2*log(x) - x**4/8.d0)
  Psi_osfa(1) = 1.d0
  Psi_osfa(2) = x**2
  Psi_osfa(3) = y**2 - x**2 * log(x)
  Psi_osfa(4) = x**4 - 4.d0 * x**2 * y**2
  Psi_osfa(5) = 2.d0*y**4 - 9.d0*y**2 * x**2 + 3.d0*x**4 *log(x) - 12.d0*x**2 * y**2 *log(x)
  Psi_osfa(6) = x**6 - 12.d0*x**4 * y**2 + 8.d0*x**2 * y**4
  Psi_osfa(7) = 8.d0*y**6 - 140.d0*y**4 * x**2 + 75.d0*y**2 * x**4 - 15.d0*x**6 *log(x) &
              + 180.d0*x**4 * y**2 *log(x) - 120.d0*x**2 * y**4 *log(x)
  Psi_osfa(8) = y
  Psi_osfa(9) = y * x**2
  Psi_osfa(10)= y**3 - 3.d0*y * x**2 *log(x)
  Psi_osfa(11)= 3.d0*y * x**4 - 4.d0* y**3 * x**2
  Psi_osfa(12)= 8.d0*y**5 - 45.d0*y *x**4 - 80.d0*y**3 * x**2 *log(x) + 60.d0*y * x**4 *log(x)

  end function Psi_osfa
!===============================================================================
  function Psix_osfa(x, y)
  real(real64), intent(in) :: x, y
  real(real64)             :: Psix_osfa(0:12)


  Psix_osfa(0) = x**3 / 2.d0 + A * (x*log(x) + x/2.d0 - x**3/2.d0)
  Psix_osfa(1) = 0.d0
  Psix_osfa(2) = 2.d0*x
  Psix_osfa(3) = -2.d0*x *log(x) - x
  Psix_osfa(4) = 4.d0*x**3 - 8.d0 * x * y**2
  Psix_osfa(5) = -30.d0*y**2 * x + 12.d0*x**3 *log(x) + 3.d0*x**3 - 24.d0*x * y**2 *log(x)
  Psix_osfa(6) = 6.d0*x**5 - 48.d0*x**3 * y**2 + 16.d0*x * y**4
  Psix_osfa(7) = -400.d0*y**4 * x + 480.d0*y**2 * x**3 - 90.d0*x**5 *log(x) &
               - 15.d0*x**5 + 720.d0*x**3 * y**2 *log(x) - 240.d0*x * y**4 *log(x)
  Psix_osfa(8) = 0.d0
  Psix_osfa(9) = 2.d0 * y * x
  Psix_osfa(10)= -6.d0*y * x * log(x) - 3.d0* y * x
  Psix_osfa(11)= 12.d0*y* x**3 - 8.d0*y**3 * x
  Psix_osfa(12)= -120.d0*y * x**3 - 160.d0*y**3 * x * log(x) - 80.d0*y**3 * x + 240.d0*y * x**3 *log(x)

  end function Psix_osfa
!===============================================================================
  function Psixx_osfa(x, y)
  real(real64), intent(in) :: x, y
  real(real64)             :: Psixx_osfa(0:12)


  Psixx_osfa(0) = 3.d0 * x**2 / 2.d0 + A * (log(x) + 3.d0/2.d0 - 3.d0*x**2/2.d0)
  Psixx_osfa(1) = 0.d0
  Psixx_osfa(2) = 2.d0
  Psixx_osfa(3) = -2.d0 *log(x) - 3.d0
  Psixx_osfa(4) = 12.d0*x**2 - 8.d0 * y**2
  Psixx_osfa(5) = -54.d0*y**2 + 36.d0*x**2 *log(x) + 21.d0*x**2 - 24.d0 * y**2 *log(x)
  Psixx_osfa(6) = 30.d0*x**4 - 144.d0*x**2 * y**2 + 16.d0 * y**4
  Psixx_osfa(7) = -640.d0*y**4 + 2160.d0*y**2 * x**2 - 450.d0*x**4 *log(x) &
                - 165.d0*x**4 + 2160.d0*x**2 * y**2 *log(x) - 240.d0 * y**4 *log(x)
  Psixx_osfa(8) = 0.d0
  Psixx_osfa(9) = 2.d0*y
  Psixx_osfa(10)= -6.d0*y *log(x) - 9.d0*y
  Psixx_osfa(11)= 36.d0*y *x**2 - 8.d0*y**3
  Psixx_osfa(12)= -120.d0*y *x**2 - 160.d0*y**3 *log(x) - 240.d0*y**3 + 720.d0*y *x**2 *log(x)

  end function Psixx_osfa
!===============================================================================
  function Psiy_osfa(x, y)
  real(real64), intent(in) :: x, y
  real(real64)             :: Psiy_osfa(0:12)


  Psiy_osfa(0) = 0.d0
  Psiy_osfa(1) = 0.d0
  Psiy_osfa(2) = 0.d0
  Psiy_osfa(3) = 2.d0*y
  Psiy_osfa(4) = -8.d0 * x**2 * y
  Psiy_osfa(5) = 8.d0*y**3 - 18.d0*y * x**2 - 24.d0*x**2 * y *log(x)
  Psiy_osfa(6) = -24.d0*x**4 * y + 32.d0*x**2 * y**3
  Psiy_osfa(7) = 48.d0*y**5 - 560.d0*y**3 * x**2 + 150.d0*y * x**4 &
               + 360.d0*x**4 * y *log(x) - 480.d0*x**2 * y**3 *log(x)
  Psiy_osfa(8) = 1.d0
  Psiy_osfa(9) = x**2
  Psiy_osfa(10)= 3.d0*y**2 - 3.d0*x**2 *log(x)
  Psiy_osfa(11)= 3.d0*x**4 - 12.d0*y**2 * x**2
  Psiy_osfa(12)= 40.d0*y**4 - 45.d0*x**4 - 240.d0*y**2 * x**2 *log(x) + 60.d0*x**4 *log(x)

  end function Psiy_osfa
!===============================================================================
  function Psiyy_osfa(x, y)
  real(real64), intent(in) :: x, y
  real(real64)             :: Psiyy_osfa(0:12)


  Psiyy_osfa(0) = 0.d0
  Psiyy_osfa(1) = 0.d0
  Psiyy_osfa(2) = 0.d0
  Psiyy_osfa(3) = 2.d0
  Psiyy_osfa(4) = -8.d0 * x**2
  Psiyy_osfa(5) = 24.d0*y**2 - 18.d0 * x**2 - 24.d0*x**2 * log(x)
  Psiyy_osfa(6) = -24.d0*x**4 + 96.d0*x**2 * y**2
  Psiyy_osfa(7) = 240.d0*y**4 - 1680.d0*y**2 * x**2 + 150.d0 * x**4 &
                + 360.d0*x**4 *log(x) - 1440.d0*x**2 * y**2 *log(x)
  Psiyy_osfa(8) = 0.d0
  Psiyy_osfa(9) = 0.d0
  Psiyy_osfa(10)= 6.d0*y
  Psiyy_osfa(11)= -24.d0*y * x**2
  Psiyy_osfa(12)= 160.d0*y**3 - 480.d0*y *x**2 *log(x)

  end function Psiyy_osfa
!===============================================================================
  function Psixy_osfa(x, y)
  real(real64), intent(in) :: x, y
  real(real64)             :: Psixy_osfa(0:12)


  Psixy_osfa(0) = 0.d0
  Psixy_osfa(1) = 0.d0
  Psixy_osfa(2) = 0.d0
  Psixy_osfa(3) = 0.d0
  Psixy_osfa(4) = -16.d0 * x * y
  Psixy_osfa(5) = -60.d0 * y * x - 48.d0*x * y * log(x)
  Psixy_osfa(6) = -96.d0*x**3 * y + 64.d0*x * y**3
  Psixy_osfa(7) = -1600.d0*x * y**3 + 960.d0*x**3 * y + 1440.d0*x**3 * y *log(x) - 960.d0 * x * y**3 *log(x)
  Psixy_osfa(8) = 0.d0
  Psixy_osfa(9) = 2.d0*x
  Psixy_osfa(10)= -6.d0*x *log(x) - 3.d0*x
  Psixy_osfa(11)= 12.d0* x**3 - 24.d0*x * y**2
  Psixy_osfa(12)= -120.d0 * x**3 - 480.d0*y**2 *x *log(x) - 240.d0*y**2 * x + 240.d0*x**3 *log(x)

  end function Psixy_osfa
!===============================================================================



!===============================================================================
! Calculate R,phi,Z components of magnetic field vector [Gauss] at r=(R,Z [cm], phi [rad])
!===============================================================================
  function amhd_get_Bf(r) result(Bf)

  real(real64), intent(in)  :: r(3)
  real(real64)              :: Bf(3)

  real(real64), parameter   :: m_to_cm = 1.d2
  real(real64) :: dPsi_dR, dPsi_dZ, R1


  ! toroidal magnetic field
  Bf      = 0.d0
  Bf(3)   =  Bt*R0/r(1)


  ! poloidal magnetic field
  R1      = r(1) / m_to_cm
  dPsi_dR = amhd_get_DPsi(r, 1, 0) * m_to_cm
  dPsi_dZ = amhd_get_DPsi(r, 0, 1) * m_to_cm
  Bf(1)   = Bf(1) - dPsi_dZ / R1
  Bf(2)   = Bf(2) + dPsi_dR / R1


  ! Tesla -> Gauss
  Bf      = Bf * 1.d4

  end function amhd_get_Bf
!===============================================================================



!===============================================================================
! Sample (derivative of) poloidal magnetic flux at r=(R,Z [cm], phi [rad])
!===============================================================================
  function amhd_get_Psi(r) result(Psi)
  use fgsl
  real(real64), intent(in)  :: r(3)
  real(real64)              :: Psi


  Psi = amhd_get_DPsi(r(1:2), 0, 0)

  end function amhd_get_Psi
!===============================================================================
  function amhd_get_DPsi(r, mR, mZ) result(DPsi)
  use fgsl
  real(real64), intent(in) :: r(2)
  integer,      intent(in) :: mR, mZ
  real(real64)             :: DPsi

  real(real64) :: x, y, vPsi(0:12)
  integer :: i


  x = r(1) / R0
  y = r(2) / R0
  if (mR == 0  .and.  mZ == 0) then
     vPsi   = Psi_osfa(x, y)
  elseif (mR == 1  .and.  mZ == 0) then
     vPsi   = Psix_osfa(x, y) / R0
  elseif (mR == 0  .and.  mZ == 1) then
     vPsi   = Psiy_osfa(x, y) / R0
  elseif (mR == 2  .and.  mZ == 0) then
     vPsi   = Psixx_osfa(x, y) / R0**2
  elseif (mR == 1  .and.  mZ == 1) then
     vPsi   = Psixy_osfa(x, y) / R0**2
  elseif (mR == 0  .and.  mZ == 2) then
     vPsi   = Psiyy_osfa(x, y) / R0**2
  else
     vPsi   = 0.d0
  endif
  DPsi = PSI0_scale*(vPsi(0) + sum(ceq(1:n)*vPsi(1:n)))

  end function amhd_get_DPsi
!===============================================================================



!===============================================================================
! Return pressure on flux surface Psi
!===============================================================================
  function amhd_get_pressure(Psi) result(P)
  use math
  real(real64), intent(in) :: Psi
  real(real64)             :: P


  P = PSI0_scale**2 / (4.d-7*pi) / (R0/1.d2)**4 * (1.d0-A) * Psi

  end function amhd_get_pressure
!===============================================================================



!===============================================================================
! Return boundaries [cm] of equilibrium domain
!===============================================================================
  subroutine amhd_get_domain(Rbox, Zbox)
  real(real64), intent(out) :: Rbox(2), Zbox(2)


  Rbox(1) = R0*max(1.d0 - 2.d0*eps, 1.d-1)
  Rbox(2) = R0*(1.d0 + 2.d0*eps)
  Zbox(1) = -R0 * 2.d0*kap*eps
  Zbox(2) =  R0 * 2.d0*kap*eps

  end subroutine amhd_get_domain
!===============================================================================


  
!===============================================================================
  subroutine amhd_broadcast
  use parallel


  call broadcast_real_s (Ip)
  call broadcast_real_s (Bt)
  call broadcast_real_s (R0)
  call broadcast_real_s (PSI0_scale)
  call broadcast_real_s (scale_manual)
  call broadcast_real_s (A)
  call broadcast_inte_s (n)
  call broadcast_real   (ceq,int(n))
  call broadcast_logi   (up_down_symmetry)

  end subroutine amhd_broadcast
!===============================================================================

end module amhd
