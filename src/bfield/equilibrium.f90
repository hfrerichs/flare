!===============================================================================
! Equilibrium related functions and subroutines
!===============================================================================
module equilibrium
  use iso_fortran_env
  use curve2D
  use magnetic_axis
  implicit none


  character(len=*), parameter :: &
     S_GEQDSK     = 'geqdsk', &
     S_DIVAMHD    = 'divamhd', &
     S_JETEQ      = 'jeteq', &
     S_M3DC1      = 'm3dc1'

  integer, parameter :: &
     EQ_GUESS     = -1, &
     EQ_UNDEFINED = 0, &
     EQ_GEQDSK    = 1, &
     EQ_DIVAMHD   = 2, &
     EQ_JET       = 3, &
     EQ_M3DC1     = 4


  integer, parameter :: nX_max = 10


  type t_Xpoint
     real(real64) :: R_estimate = 0.d0, Z_estimate = 0.d0
     real(real64) :: X(2) = -1.d0, Psi, H(2,2), theta
     logical      :: undefined = .true.

     contains
     procedure :: analysis
     procedure :: load
     procedure :: PsiN => t_Xpoint__PsiN
  end type t_Xpoint


!...............................................................................
! user defined parameters (to be set via configuration file)                   .

  character*120 :: &
     Data_File        = ''
  character*12  :: &
     Data_Format      = ''
  character(len=256) :: Magnetic_Axis_File = ''

  real(real64) :: &
     R_axis       = 0.d0, &        ! user defined position of magnetic axis
     Z_axis       = 0.d0, &
     R_sepx       = 0.d0, &        ! user defined position of separatrix
     Z_sepx       = 0.d0, &
     Bt, R0, &   ! reference toroidal magnetic field [T] and radial position [cm]
     Ip           = 0.d0           ! plasma current [A] (equilibrium will be re-scaled)


  type(t_Xpoint) :: Xp(nX_max), Magnetic_Axis


  logical :: &
     use_boundary     = .true., &
     Current_Fix      = .true.

  integer :: &
     Diagnostic_Level = 0

  namelist /Equilibrium_Input/ &
     Data_File, Data_Format, use_boundary, Current_Fix, Diagnostic_Level, &
     R_axis, Z_axis, R_sepx, Z_sepx, Bt, R0, Ip, Xp, Magnetic_Axis, &
     Magnetic_Axis_File
!...............................................................................



!...............................................................................
! public variables                                                             .

  ! Magnetic axis (in axisymmetric configuration)
  type t_Maxis
     real*8 :: R, Z
  end type t_Maxis
  type (t_Maxis) :: Maxis2D


  ! Direction of toroidal magnetic field (Bt) and plasma current (Ip)
  ! +1: positive direction, i.e. counter-clockwise
  ! -1: negative direction, i.e. clockwise
  !  0: no equilibrium defined
! temporarily moved to module magnetic axis
!  integer :: &
!     Bt_sign  = 0, &
!     Ip_sign  = 0


  ! equilibrium type
  integer :: i_equi = EQ_UNDEFINED


  ! Position of magnetic axis, poloidal magnetic flux at separatrix and magnetic axis
  !real*8 :: R_axis, Z_axis, Psi_sepx, Psi_axis
  real*8 :: &
     Psi_sepx = 1.d0, &
     Psi_axis = 0.d0








!...............................................................................
! Interfaces for functions/subroutines from specific equilibrium types         .

  ! get equilibrium magnetic field in cylindrical coordinates
  procedure(default_get_Bf), pointer :: get_Bf_eq2D  => default_get_Bf

  ! get poloidal magnetic flux
  procedure(default_get_Psi), pointer :: get_Psi => default_get_Psi

  ! get derivative of poloidal magnetic flux
  procedure(default_get_DPsi), pointer :: get_DPsi => default_get_DPsi

!  ! return poloidal magnetic flux at magnetic axis
!  procedure(Psi_axis_interface), pointer :: Psi_axis

  ! Return boundaries [cm] of equilibrium domain
  procedure(default_get_domain), pointer :: get_domain => default_get_domain

  ! inquire boundary setup from equilibrium data
  interface
     function logical_inquiry() result(l)
     logical :: l
     end function logical_inquiry

     subroutine export_curve(S)
     import :: t_curve
     type(t_curve), intent(out) :: S
     end subroutine export_curve
  end interface
  procedure(logical_inquiry), pointer :: &
     equilibrium_provides_boundary => default_equilibrium_provides_boundary
  procedure(export_curve), pointer    :: export_boundary

  ! Broadcast data for parallel execution
  procedure(), pointer :: broadcast_equilibrium
!...............................................................................


  logical, save :: initialized = .false.
!, get_equi_domain
!


  integer :: ipanic = 0

  contains
!=======================================================================



!=======================================================================
! Load equilibrium configuration
!=======================================================================
  subroutine load_equilibrium_config (iu, iconfig)
  use run_control, only: Prefix, Debug
  integer, intent(in)  :: iu
  integer, intent(out) :: iconfig

  integer              :: ix
  real(real64)         :: r(3), x0(2)
  real(real64)         :: Rbox(2), Zbox(2)


! 1. read user configuration
  rewind (iu)
  read   (iu, Equilibrium_Input, end=1000)
  iconfig = 1
  write (6, *)
  write (6, 1001)

! 1.b find equilibrium type
  select case(Data_Format)
  case (S_GEQDSK)
     i_equi = EQ_GEQDSK
  case (S_DIVAMHD)
     i_equi = EQ_DIVAMHD
  case (S_M3DC1)
     i_equi = EQ_M3DC1
  case ('')
     i_equi = EQ_GUESS
  case default
     write (6, *) 'error: ', Data_Format, ' is not a valid equilibrium type!'
     stop
  end select


! set default values
  export_boundary => null()


! 2. load equilibrium data (if provided) ...............................
  if (Data_File .ne. '') call load_equilibrium_data()
  call setup_equilibrium()


! 3. setup magnetic axis ...............................................
  ! 3.1 find magnetic axis from estimated position
  if (Magnetic_Axis%R_estimate > 0.d0) then
     x0(1)             = Magnetic_Axis%R_estimate
     x0(2)             = Magnetic_Axis%Z_estimate
     Magnetic_Axis%X   = find_x(x0)

     r(1:2)            = Magnetic_Axis%X; r(3) = 0.d0
     Magnetic_Axis%Psi = get_Psi(r)
     Psi_Axis          = Magnetic_Axis%Psi
     call setup_magnetic_axis_2D (r(1), r(2))
  endif

  ! setup user defined magnetic axis (if provided) .....................
  ! 3.2. axisymmetric configuration
  if (R_axis > 0.d0) then
     call setup_magnetic_axis_2D (R_axis, Z_axis)

     ! set pol. magn. flux on axis
     r(1) = R_axis
     r(2) = Z_axis
     r(3) = 0.d0
     Psi_axis = get_Psi(r)
  endif

  ! 3.3. non-axisymmetric configuration
  if (Magnetic_Axis_File .ne. '') then
     Data_File = trim(Prefix)//Magnetic_Axis_File
     call load_magnetic_axis_3D (Data_File)

     ! pol. magn. flux not supported for non-axisymmetric configurations
     Psi_axis = 0.d0
  endif
  write (6, 3000) Psi_axis


! 4. set up X-points ...................................................
  write (6, 4000)
  do ix=1,nx_max
     if (Xp(ix)%R_estimate <= 0.d0) then
        if (ix == 1) then
           call get_domain (Rbox, Zbox)
           Xp(ix)%R_estimate = Rbox(1) + 1.d0/3.d0 * (Rbox(2)-Rbox(1))
           Xp(ix)%Z_estimate = Zbox(1) + 1.d0/6.d0 * (Zbox(2)-Zbox(1))
        else
           cycle
        endif
     endif

     x0(1)        = Xp(ix)%R_estimate
     x0(2)        = Xp(ix)%Z_estimate
     Xp(ix)%X     = find_x(x0, Hout=Xp(ix)%H)
     if (Xp(ix)%X(1) > 0.d0) Xp(ix)%undefined = .false.

     r(1:2)       = Xp(ix)%X; r(3) = 0.d0
     Xp(ix)%Psi   = get_Psi(r)
     Xp(ix)%theta = get_poloidal_angle(r)
     if (ix == 1) Psi_sepx = Xp(ix)%Psi

     write (6, 4001) Xp(ix)%X, Xp(ix)%Psi
     if (Debug) then
        write (6, 4002) Xp(ix)%H(1,1), Xp(ix)%H(1,2)
        write (6, 4002) Xp(ix)%H(2,1), Xp(ix)%H(2,2)
     endif
  enddo


! 5. set dependent variables ...........................................
  ! pol. magn. flux at separatrix
  if (R_sepx > 0.d0) then
     r(1) = R_sepx
     r(2) = Z_sepx
     r(3) = 0.d0
     Psi_sepx = get_Psi(r)
  endif
  write (6, 5000) Psi_sepx


  return
 1000 iconfig = 0
 1001 format ('   - Equilibrium configuration:')
 3000 format (8x,'Psi_axis = ', e12.4)
 4000 format (3x,'- Configuring X-point(s): (R, Z, Psi)')
 4001 format (8x,2f12.4,2x,e12.4)
 4002 format (8x,2e12.4)
 5000 format (8x,'Psi_sepx = ', e12.4)
  end subroutine load_equilibrium_config
!=======================================================================



!=======================================================================
  subroutine load_equilibrium_data
  use run_control, only: Prefix
  use geqdsk
  use divamhd

  integer, parameter :: iu_scan = 17

  character*80 :: s


  Data_File = trim(Prefix)//Data_File

! determine equilibrium type (if not provided) .........................
  if (i_equi == EQ_GUESS) then
     open  (iu_scan, file=Data_file)
     read  (iu_scan, '(a80)') s
     if (s(3:5) == 'TEQ'  .or.  s(3:6) == 'EFIT') then
        i_equi = EQ_GEQDSK
     elseif (s(5:11) == 'jm   :=') then
        i_equi = EQ_JET
     else
        read  (iu_scan, '(a80)') s
        if (s(4:9) == 'File: ') then
           i_equi = EQ_DIVAMHD
        else
           i_equi = EQ_UNDEFINED
        endif
     endif
     close (iu_scan)
  endif
! ... determine equilibrium type (done) ................................



! load equilibrium data
  select case (i_equi)
  case (EQ_GEQDSK)
     call geqdsk_load (Data_File, use_boundary, Current_Fix, Diagnostic_Level, Psi_axis, Psi_sepx)
  case (EQ_DIVAMHD)
     call divamhd_load (Data_File, Ip, Bt, R0)

  case (EQ_M3DC1)
     ! nothing to be done here

  case default
     write (6, *) 'error: cannot determine equilibrium type!'
     stop
  end select

  end subroutine load_equilibrium_data
!=======================================================================



!=======================================================================
! Setup procedure pointers
!=======================================================================
  subroutine setup_equilibrium()
  use geqdsk
  use divamhd
  use m3dc1

  ! select case equilibrium
  select case (i_equi)
  case (EQ_GEQDSK)
     get_Bf_eq2D                   => geqdsk_get_Bf
     get_Psi                       => geqdsk_get_Psi
     get_DPsi                      => geqdsk_get_DPsi
     get_domain                    => geqdsk_get_domain
     equilibrium_provides_boundary => geqdsk_provides_boundary
     export_boundary               => geqdsk_export_boundary
     broadcast_equilibrium         => geqdsk_broadcast
  case (EQ_DIVAMHD)
     get_Bf_eq2D                   => divamhd_get_Bf
     get_Psi                       => divamhd_get_Psi
     get_DPsi                      => divamhd_get_DPsi
     get_domain                    => divamhd_get_domain
     broadcast_equilibrium         => divamhd_broadcast
  case (EQ_M3DC1)
     get_Bf_eq2D                   => m3dc1_get_Bf_eq2D
     get_Psi                       => m3dc1_get_Psi
     get_DPsi                      => m3dc1_get_DPsi
  end select

  end subroutine setup_equilibrium
!=======================================================================



!=======================================================================
! Broadcast equilibrium data
!=======================================================================
  subroutine broadcast_mod_equilibrium()
  use parallel

  integer :: i


  if (nprs == 1) return

  call broadcast_real_s (Psi_axis)
  call broadcast_real_s (Psi_sepx)
  call broadcast_inte_s (i_equi)
  do i=1,nX_max
     call broadcast_real_s (Xp(i)%R_estimate)
     call broadcast_real_s (Xp(i)%Z_estimate)
     call broadcast_real_s (Xp(i)%Psi)
     call broadcast_real   (Xp(i)%X, 2)
     call broadcast_real   (Xp(i)%H, 4)
  enddo
  call broadcast_magnetic_axis()
  call wait_pe()



  if (mype > 0) call setup_equilibrium()
  call broadcast_equilibrium()

  end subroutine broadcast_mod_equilibrium
!=======================================================================



!=======================================================================
! Sample magnetic field vector
!=======================================================================
  function default_get_Bf(r) result(Bf)
  real*8, intent(in)  :: r(3)
  real*8              :: Bf(3)

  Bf = 0.d0
  if (ipanic > 0) then
     write (6, *) 'error: magnetic field function not defined!'
     stop
  endif

  end function default_get_Bf
!=======================================================================



!=======================================================================
! Sample poloidal magnetic flux at r=(R,Z [cm], phi [rad])
!===============================================================================
  function default_get_Psi(r) result(Psi)
  real*8, intent(in)  :: r(3)
  real*8              :: Psi

  Psi = 0.d0
  if (ipanic > 0) then
     write (6, *) 'error: poloidal magnetic flux function not defined!'
     stop
  endif

  end function default_get_Psi
!=======================================================================



!=======================================================================
! Sample normalized poloidal magnetic flux at r=(R,Z [cm], phi [rad])
!===============================================================================
  function get_PsiN(r0) result(PsiN)
  real(real64), dimension(:), intent(in) :: r0
  real(real64)                           :: PsiN

  real(real64) :: r(3)


  select case(size(r0))
  case(2)
     r(1:2) = r0
     r(  3) = 0.d0
  case(3)
     r      = r0
  case default
     write (6, *) 'error in get_PsiN: invalid size of argument r0!'
     write (6, *) 'size(r0) = ', size(r0)
  end select

  PsiN = (get_Psi(r) - Psi_axis) / (Psi_sepx - Psi_axis)

  end function get_PsiN
!=======================================================================



!=======================================================================
! Sample (nR,nZ)-th derivative of poloidal magnetic flux at r=(R,Z [cm])
!=======================================================================
  function default_get_DPsi (r, nR, nZ) result(DPsi)
  real*8, intent(in)  :: r(2)
  integer, intent(in) :: nR, nZ
  real*8              :: DPsi

  DPsi = 0.d0
  if (ipanic > 0) then
     write (6, *) 'error: derivative of poloidal magnetic flux function not defined!'
     stop
  endif

  end function default_get_DPsi
!=======================================================================
  function get_DPsiN(r, nR, nZ) result(DPsiN)
  real*8, intent(in)  :: r(3)
  integer, intent(in) :: nR, nZ
  real*8              :: DPsiN

  DPsiN = get_DPsi(r(1:2), nR, nZ) / (Psi_sepx - Psi_axis)

  end function get_DPsiN
!=======================================================================
  function get_ePsi(r) result(ePsi)
  real(real64), intent(in) :: r(3)
  real(real64)             :: ePsi(2)

  real(real64) :: D


  ePsi(1) = get_DPsiN(r, 1, 0)
  ePsi(2) = get_DPsiN(r, 0, 1)
  D       = sqrt(sum(ePsi**2))

  if (D > 0.d0) ePsi    = ePsi / D

  end function get_ePsi
!=======================================================================
  function get_H(x) result(H)
  real(real64), intent(in) :: x(2)
  real(real64)             :: H(2,2)


  H(1,1) = get_DPsi(x, 2, 0)
  H(1,2) = get_DPsi(x, 1, 1)
  H(2,2) = get_DPsi(x, 0, 2)
  H(2,1) = H(1,2)
  if (H(1,1) == 0.d0  .and.  H(1,2) == 0.d0  .and.  H(2,2) == 0.d0) then
     H = approximate_H(x)
  endif

  end function get_H
!=======================================================================
  function approximate_H(x) result(H)
  use run_control, only: Debug
  real(real64), intent(in) :: x(2)
  real(real64) :: H(2,2)

  real(real64), parameter :: gamma_1 = 1.d-1
  real(real64) :: xi(2), dfdxi(2), dfdyi(2)

  integer :: i, j


  do j=1,2
     do i=1,2
        xi       = x
        xi(3-j)  = x(3-j) + (-1.d0)**i * gamma_1

        dfdxi(i) = get_DPsi(xi, 1, 0)
        dfdyi(i) = get_DPsi(xi, 0, 1)
     enddo
     H(1,3-j) = (dfdxi(2) - dfdxi(1)) / 2.d0 / gamma_1
     H(2,3-j) = (dfdyi(2) - dfdyi(1)) / 2.d0 / gamma_1
  enddo
  H(1,2) = 0.5d0 * (H(1,2) + H(2,1))
  H(2,1) = H(1,2)

  if (Debug) then
     write (6, *) 'H_approximate = '
     write (6, *) H(1,1), H(1,2)
     write (6, *) H(2,1), H(2,2)
  endif

  end function approximate_H
!=======================================================================



!=======================================================================
! Return poloidal angle [rad] at r=(R,Z [cm], phi [rad])
!=======================================================================
  function get_poloidal_angle(r) result(theta)
  real*8, intent(in) :: r(3)
  real*8             :: theta, Maxis(3)

  Maxis = get_magnetic_axis(r(3))
  theta = atan2(r(2) - Maxis(2), r(1) - Maxis(1))

  end function get_poloidal_angle
!=======================================================================



!=======================================================================
! Get cylindrical coordinates (R[cm], Z[cm], Phi[rad]) for flux
! coordinates (Theta[deg], PsiN, Phi[deg])
!=======================================================================
  function get_cylindrical_coordinates(y, ierr, r0) result(r)
  use iso_fortran_env
  use math
  implicit none

  real(real64), intent(inout)        :: y(3)
  integer,      intent(out)          :: ierr
  real(real64), intent(in), optional :: r0(3)
  real(real64)                       :: r(3)

  integer, parameter :: imax = 160
  real(real64), parameter :: tolerance = 1.d-10
  real(real64), parameter :: damping   = 0.9d0

  real(real64) :: dl, dpsi_dR, dpsi_dZ, dpsi_dl
  real(real64) :: M(3), Theta, PsiN, dr(2), beta

  integer :: i


  ierr  = 0

  ! set start point for approximation
  if (present(r0)) then
     r = r0
  else
     ! start near magnetic axis
     r(3)  = y(3) / 180.d0*pi
     M     = get_magnetic_axis(r(3))
     dr(1) = cos(y(1)/180.d0*pi)
     dr(2) = sin(y(1)/180.d0*pi)
     dl    = 0.2d0 * length_scale()

     r(1:2)= M(1:2) + dl*dr
     PsiN  = get_PsiN(r)
  endif


  do i=1,imax
     dpsi_dR = get_DPsiN(r, 1, 0)
     dpsi_dZ = get_DPsiN(r, 0, 1)
     dr(1) = cos(y(1)/180.d0*pi)
     dr(2) = sin(y(1)/180.d0*pi)
     dpsi_dl = dpsi_dR*dr(1) + dpsi_dZ*dr(2)

     beta    = y(2) - PsiN
     dr      = dr * beta / dpsi_dl * damping

     r(1:2)  = r(1:2) + dr
     PsiN    = get_PsiN(r)
     Theta   = get_poloidal_angle(r) / pi*180.d0
     if (Theta < 0) Theta = Theta + 360.d0

     if (abs(beta) <= tolerance) exit
  enddo
  ! update input parameters to match output coordinates
  y(1) = Theta
  y(2) = PsiN

  if (abs(beta) > tolerance) then
     ierr = 1
  endif

  end function get_cylindrical_coordinates
!=======================================================================



!=======================================================================
! Return boundaries [cm] of equilibrium domain
!=======================================================================
  subroutine default_get_domain (Rbox, Zbox)
  real(real64), intent(out) :: Rbox(2), Zbox(2)


  Rbox(1) = 0.d0
  Rbox(2) = huge(0.d0)
  Zbox(1) = -huge(0.d0)/2.d0
  Zbox(2) = huge(0.d0)
  if (ipanic > 0) then
     write (6, *) 'error: equilibrium domain not defined!'
     stop
  endif

  end subroutine default_get_domain
!=======================================================================



!=======================================================================
! boundary provided by equilibrium
!=======================================================================
  function default_equilibrium_provides_boundary() result(l)
  logical :: l

  l = .false.
  end function default_equilibrium_provides_boundary
!=======================================================================
! export axisymmetric boundary provided by equilibrium
!=======================================================================
!  subroutine export_PFC (S)
!  type(t_curve), intent(out) :: S
!  end subroutine export_PFC
!=======================================================================







!=======================================================================
! Find the X-point of a magnetic configuration, provided an initial guess X
!
! This subroutine applies the Newton-method for the iterative approximation of
! a critical point
!=======================================================================
  function find_X (X0, verbose, Hout) result(X)
  use bspline
  implicit none

  real(real64), intent(in)            :: X0(2)
  real(real64)                        :: X(2)
  logical,      intent(in),  optional :: verbose
  real(real64), intent(out), optional :: Hout(2,2)


  real(real64), parameter :: &
      gamma_0 = 1.d0, &
      delta   = 1.d-8
  integer, parameter :: nmax = 2000


  real(real64) :: xn(2), dx(2), dfdx, dfdy, H(2,2), Hdisc, dxmod, Rbox(2), Zbox(2)
  logical :: screen_output
  integer :: n


  ! setup screen output
  screen_output = .false.
  if (present(verbose)) screen_output = .true.


  ! initialize
  call get_domain (Rbox, Zbox)
  xn = X0
  if (screen_output) write (6, *) 'Initial guess for X-point: ', xn
  n  = 0

  approximation_loop: do
     ! check boundaries
     if (xn(1).lt.Rbox(1) .or. xn(1).gt.Rbox(2) .or. &
         xn(2).lt.Zbox(1) .or. xn(2).gt.Zbox(2)) then
        X = -1.d0
        return
     endif

     ! calculate the gradient
     dfdx = get_DPsi(xn, 1, 0)
     dfdy = get_DPsi(xn, 0, 1)

     ! calculate elements of the Hessian matrix
     H(1,1) = get_DPsi(xn, 2, 0)
     H(1,2) = get_DPsi(xn, 1, 1)
     H(2,2) = get_DPsi(xn, 0, 2)
     H(2,1) = H(1,2)
     Hdisc  = H(1,1) * H(2,2) - H(1,2)*H(2,1)
     if (Hdisc.eq.0) then
        !X = find_X_BFGS(xn, Hout)
        X = find_X_SR1(xn, Hout)
        return
     endif

     ! calculate increment
     dx(1)  =   H(2,2)*dfdx - H(1,2)*dfdy
     dx(2)  = - H(2,1)*dfdx + H(1,1)*dfdy

     xn     = xn - gamma_0*dx/Hdisc
     dxmod  = sqrt(sum(dx**2))/dabs(Hdisc)
     n      = n + 1

     if (dxmod .lt. delta) exit approximation_loop

     if (n.gt.nmax) exit approximation_loop
  enddo approximation_loop

  X = xn
  if (present(Hout)) Hout = H
  end function find_X
!=======================================================================
  function find_X_BFGS(X0, Hout) result(X)
  use run_control, only: Debug
  real(real64), intent(in)            :: X0(2)
  real(real64)                        :: X(2)
  real(real64), intent(out), optional :: Hout(2,2)

  real(real64), parameter :: &
      gamma_1 = 1.d0, &
      delta   = 1.d-13
  integer, parameter :: &
      nmax = 1000000

  real(real64) :: dxmod, P(2), S(2), B(2,2), Bdisc, dfdx, dfdy, dfdx1, dfdy1
  real(real64) :: DB(2,2), Y(2), ys, sBs
  integer :: i, j, k


  if (Debug) open  (90, file='X_BFGS.tmp')
  X = X0
  B = approximate_H(X) ! get an initial approximate
  Bdisc  = B(1,1) * B(2,2) - B(1,2)*B(2,1)
  if (Bdisc == 0.d0) then
     write (6, *) 'error: initial estimate of Hessian matrix has determinant = 0!'
     stop
  endif
  ! calculate the gradient
  dfdx = get_DPsi(X, 1, 0)
  dfdy = get_DPsi(X, 0, 1)
  if (Debug) write (90, *) 0, X, dfdx, dfdy


  approximation_loop: do k=1,nmax
     ! 1. obtain direction (solve B * p = - grad f)
     P(1)  =   B(2,2)*dfdx - B(1,2)*dfdy
     P(2)  = - B(2,1)*dfdx + B(1,1)*dfdy
     P     = - P / Bdisc

     ! 2. use constant gamma_1

     ! 3. calculate increment
     S     = gamma_1 * P
     X     = X + S
     dxmod = sqrt(sum(S**2))
     if (Debug) write (90, *) k, X, S, dxmod
     if (dxmod < delta) exit

     ! 4. calculate local gradient
     dfdx1 = get_DPsi(X, 1, 0)
     dfdy1 = get_DPsi(X, 0, 1)
     Y(1)  = dfdx1 - dfdx
     Y(2)  = dfdy1 - dfdy
     dfdx  = dfdx1
     dfdy  = dfdy1

     ! 5. update Hessian approximation
     ys = sum(Y*S)
     sBs = B(1,1)*S(1)**2 + (B(1,2)+B(2,1))*S(1)*S(2) + B(2,2)*S(2)**2
     do i=1,2
     do j=1,2
        DB(i,j) = Y(i)*Y(j) / ys - (B(i,1)*S(1)+B(i,2)*S(2)) * (S(1)*B(1,j)+S(2)*B(2,j)) / sBs
     enddo
     enddo
     B = B + DB
     Bdisc  = B(1,1) * B(2,2) - B(1,2)*B(2,1)
     if (present(Hout)) Hout = B
  enddo approximation_loop
  if (Debug) close (90)

  end function find_X_BFGS
!=======================================================================
  function find_X_SR1(X0, Bout) result(X)
  use run_control, only: Debug
  real(real64), intent(in)            :: X0(2)
  real(real64)                        :: X(2)
  real(real64), intent(out), optional :: Bout(2,2)

  real(real64), parameter :: &
      gamma_1 = 1.d0, &
      delta   = 1.d-13
  integer, parameter :: &
      nmax = 1000000

  real(real64) :: dxmod, P(2), S(2), B(2,2), Bdisc, dfdx, dfdy, dfdx1, dfdy1
  real(real64) :: DB(2,2), Y(2), v(2), vs, smod, vmod
  integer :: i, j, k


  if (Debug) open  (90, file='X_SR1.tmp')
  X = X0

  ! calculate initial approximation of Hessian (B)
  B = approximate_H(X) ! get an initial approximate
  Bdisc  = B(1,1) * B(2,2) - B(1,2)*B(2,1)
  if (Bdisc == 0.d0) then
     write (6, *) 'error: initial estimate of Hessian matrix has determinant = 0!'
     stop
  endif
  ! calculate the gradient
  dfdx = get_DPsi(X, 1, 0)
  dfdy = get_DPsi(X, 0, 1)
  if (Debug) write (90, *) 0, X, dfdx, dfdy

  approximation_loop: do k=1,nmax
     ! 1. obtain direction (solve B * p = - grad f)
     P(1)  =   B(2,2)*dfdx - B(1,2)*dfdy
     P(2)  = - B(2,1)*dfdx + B(1,1)*dfdy
     P     = - P / Bdisc

     ! 2. use constant gamma_1

     ! 3. calculate increment
     S     = gamma_1 * P
     X     = X + S
     dxmod = sqrt(sum(S**2))
     if (Debug) write (90, *) k, X, S, dxmod
     if (dxmod < delta) exit

     ! 4. calculate local gradient
     dfdx1 = get_DPsi(X, 1, 0)
     dfdy1 = get_DPsi(X, 0, 1)
     Y(1)  = dfdx1 - dfdx
     Y(2)  = dfdy1 - dfdy
     dfdx  = dfdx1
     dfdy  = dfdy1

     ! 5. update Hessian approximation
     ! v = Y - B*S
     do i=1,2
        v(i)  = Y(i) - (B(i,1)*S(1) + B(i,2)*S(2))
     enddo
     vs = sum(v*S)

     smod = sqrt(sum(s**2))
     vmod = sqrt(sum(v**2))
     if (abs(vs) > 1.d-8*smod*vmod) then
     do i=1,2
     do j=1,2
        DB(i,j) = v(i)*v(j) / vs
     enddo
     enddo
     B = B + DB
     Bdisc  = B(1,1) * B(2,2) - B(1,2)*B(2,1)
     endif
     if (present(Bout)) Bout = B
  enddo approximation_loop
  if (Debug) close (90)

  end function find_X_SR1
!=======================================================================
  function find_lX() result(X)
  real(real64) :: X(2)

  real(real64) :: Rbox(2), Zbox(2), X0(2)


  call get_domain (Rbox, Zbox)
  ! try to find lower X at relative coordinate (1/3, 1/6) from lower left corner
  X0(1) = Rbox(1) + 1.d0/3.d0 * (Rbox(2)-Rbox(1))
  X0(2) = Zbox(1) + 1.d0/6.d0 * (Zbox(2)-Zbox(1))
  X     = find_X(X0)

  end function find_lX
!=======================================================================
  function find_uX() result(X)
  real(real64) :: X(2)

  real(real64) :: Rbox(2), Zbox(2), X0(2)


  call get_domain (Rbox, Zbox)
  ! try to find upper X at relative coordinate (1/3, 5/6) from lower left corner
  X0(1) = Rbox(1) + 1.d0/3.d0 * (Rbox(2)-Rbox(1))
  X0(2) = Zbox(1) + 5.d0/6.d0 * (Zbox(2)-Zbox(1))
  X     = find_X(X0)

  end function find_uX
!=======================================================================



!=======================================================================
  subroutine find_hyperbolic_points()
  real(real64)         :: Rbox(2), Zbox(2)

  integer, parameter :: nR = 20, nZ = 20, iu = 54

  type(t_Xpoint) :: Xp
  real(real64)   :: x(2), H(2,2), xk(nR*nZ, 2), r, lambda1, lambda2, v1(2), v2(2)
  real(real64)   :: DPsi, DPsi1, r3(3)
  integer        :: i, j, k, ind, ierr, iPsi


  call get_domain (Rbox, Zbox)
  write (6, 1000)
  write (6, 1001) Rbox
  write (6, 1001) Zbox

  ind = 0
  open  (iu, file='hyperbolic_points.dat')
  loop2: do i=0, nR
  loop1: do j=0, nZ
     x(1) = Rbox(1) + (Rbox(2)-Rbox(1)) * i / nR
     x(2) = Zbox(1) + (Zbox(2)-Zbox(1)) * j / nZ
     x = find_X(x, Hout=H)

     if (x(1) < 0.d0) cycle ! not a valid critical point

     ! check if present critical point is identical to previous ones
     do k=1,ind
        r = sqrt(sum((xk(k,:)-x)**2))
        if (r < 1.d-8) cycle loop1
     enddo

     ! so this is a new critical point, run analysis
     Xp%X = x; Xp%H = H
     call Xp%analysis(lambda1, lambda2, v1, v2, ierr)
     if (ierr .ne. 0) cycle ! this is not a hyperbolic point

     ! add present point to list
     ind = ind + 1
     xk(ind,:) = x
     write (6, 1002) ind, x, lambda1, lambda2
     write (iu, *) x, lambda1, lambda2
  enddo loop1
  enddo loop2
  close (iu)


!  ! find primary X-point
!  DPsi1 = huge(1.d0)
!  iPsi  = 0
!  r3(3) = 0.d0
!  do k=1,ind
!     r3(1:2) = xk(k,:)
!     DPsi = abs(get_Psi(r3)-Psi_axis)
!     if (DPsi < DPsi1) then
!        DPsi1 = DPsi
!        iPsi  = k
!     endif
!  enddo
!  write (6, *) 'primary X-point is: ', xk(iPsi,:)
  write (6, *)

 1000 format(3x,'- Running search for hyperbolic points in domain:')
 1001 format(8x,2f12.4)
 1002 format(8x,i0,4x,'(',f10.4,', ',f10.4,')',4x,'l1 = ',e12.4,',',4x,'l2 = ',e12.4)
  end subroutine find_hyperbolic_points
!=======================================================================



!=======================================================================
  function length_scale () result(L)
  real(real64) :: L

  real(real64) :: M(3)


  L = 1.d0
  M = get_magnetic_axis(0.d0)
  if (M(1) > 0.d0) L = M(1)

  end function length_scale
!=======================================================================




!=======================================================================
  subroutine analysis(this, lambda1, lambda2, v1, v2, ierr)
  use math
  use run_control, only: Debug
  class(t_Xpoint) :: this
  real(real64), intent(out) :: lambda1, lambda2, v1(2), v2(2)
  integer,      intent(out) :: ierr

  real(real64) :: A(2,2), P, Q, phi(2)


  ierr = 0
  A(1,1) = -this%H(1,2); A(2,1) = this%H(1,1)
  A(1,2) = -this%H(2,2); A(2,2) = this%H(2,1)

  P = 0.5d0 * (A(1,1) + A(2,2))
  Q = P**2 - (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  if (Q < 0.d0) then
     if (Debug) then
        write (6, *) 'error: no real eigenvalues!'
        write (6, *) 'failed to analyze X-point!'
        write (6, *) 'A ='
        write (6, *) A(1,1), A(1,2)
        write (6, *) A(2,1), A(2,2)
        write (6, *) 'P, Q = ', P, Q
     endif
     ierr = 1
     return
  endif
  lambda1 = P + sqrt(Q)
  lambda2 = P - sqrt(Q)

  phi(1) = atan2(-A(1,1) + lambda1, A(1,2)); if (phi(1) < 0.d0) phi(1) = phi(1) + pi
  phi(2) = atan2(-A(1,1) + lambda2, A(1,2)); if (phi(2) < 0.d0) phi(2) = phi(2) + pi

  v1(1) = cos(phi(1))
  v1(2) = sin(phi(1))
  v2(1) = cos(phi(2))
  v2(2) = sin(phi(2))


  if (Debug) then
     open  (98, file='A.tmp')
     write (98, *) A(1,1), A(1,2)
     write (98, *) A(2,1), A(2,2)
     close (98)
     open  (97, file='unstable.tmp')
     write (97, *) this%X
     write (97, *) this%X + v1
     close (97)
     open  (96, file='stable.tmp')
     write (96, *) this%X
     write (96, *) this%X + v2
     close (96)
  endif

  end subroutine analysis
!=======================================================================



!=======================================================================
  function load(this) result(X)
  class(t_Xpoint) :: this
  real(real64)    :: X(2)


  if (this%undefined) then
     X = -1.d0
     write (6, *) 'error: X-point is not defined!'
     stop
  else
     X = this%X
  endif

  end function load
!=======================================================================



!=======================================================================
  function t_Xpoint__PsiN(this) result(PsiN)
  class(t_Xpoint) :: this
  real(real64)    :: PsiN


  PsiN = (this%Psi - Psi_axis) / (Psi_sepx - Psi_axis)

  end function t_Xpoint__PsiN
!=======================================================================







!=======================================================================
  function pol_flux(r) result(psi)
  real*8, intent(in) :: r(3)
  real*8             :: psi

  end function pol_flux
!=======================================================================


!=======================================================================
  subroutine Bf_pol_sub (n, s, y, f)
  integer, intent(in) :: n
  real*8, intent(in)  :: s, y(n)
  real*8, intent(out) :: f(n)

  ! n = 2, y(1) = R, y(2) = Z
  end subroutine Bf_pol_sub
!=======================================================================

end module equilibrium
