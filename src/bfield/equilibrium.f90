!===============================================================================
! Equilibrium related functions and subroutines
!===============================================================================
module equilibrium
  use iso_fortran_env
  use system
  use equilibrium_format
  use magnetic_axis
  use curve2D
  use abstract_bfield
  use sonnet
  implicit none


  integer, parameter :: nX_max = 20


  type t_Xpoint ! critical point
     real(real64) :: R_estimate = 0.d0, Z_estimate = 0.d0
     real(real64) :: X(2) = -1.d0, Psi, H(2,2), theta
     logical      :: undefined = .true.

     contains
     procedure :: analysis
     procedure :: load
     procedure :: setup
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
     Bt           = 0.d0, &        ! reference toroidal magnetic field [T]
     R0           = 0.d0, &        ! and radial position [cm]
     Ip           = 0.d0           ! plasma current [A] (equilibrium will be re-scaled)


  type(t_Xpoint) :: Xp(nX_max), M

  class(t_equi2d), pointer :: Bequi => null()

  logical :: &
     use_boundary     = .true., &
     Current_Fix      = .true.

  integer :: &
     Diagnostic_Level = 0

  namelist /Equilibrium_Input/ &
     Data_File, Data_Format, use_boundary, Current_Fix, Diagnostic_Level, &
     R_axis, Z_axis, R_sepx, Z_sepx, Bt, R0, Ip, Xp, M, &
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

  real(real64) :: EQBox(2,2)







!...............................................................................
! Interfaces for functions/subroutines from specific equilibrium types         .

  ! get equilibrium magnetic field in cylindrical coordinates
  procedure(default_get_Bf), pointer  :: get_Bf_eq2D
  procedure(default_get_JBf), pointer :: get_JBf_eq2D

  ! get poloidal magnetic flux
  procedure(default_get_Psi), pointer :: get_Psi

  ! get derivative of poloidal magnetic flux
  procedure(default_get_DPsi), pointer :: get_DPsi

!  ! return poloidal magnetic flux at magnetic axis
!  procedure(Psi_axis_interface), pointer :: Psi_axis
  ! pressure profile
  procedure(default_pressure), pointer :: get_pressure

  ! Return boundaries [cm] of equilibrium domain
  procedure(default_get_domain), pointer :: get_domain

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
  procedure(export_curve), pointer    :: export_boundary

  ! Broadcast data for parallel execution
  procedure(), pointer :: broadcast_equilibrium
  procedure(), pointer :: equilibrium_info
  procedure(), pointer :: post_setup_equilibrium
!...............................................................................


  logical, save :: initialized = .false.
!, get_equi_domain
!


  integer :: ipanic = 0

  integer, private :: i_format = FREE

  contains
!=======================================================================



!=======================================================================
! Load equilibrium configuration
!=======================================================================
  subroutine load_equilibrium_config (iu, iconfig)
  use run_control, only: Prefix, use_boundary_from_equilibrium
  use system
  use m3dc1
  use amhd
  integer, intent(in)  :: iu
  integer, intent(out) :: iconfig

  character(len=256)   :: filename
  integer              :: ix, ierr
  real(real64)         :: r(3)


! 0. initialize
  i_format         = STRICT
  get_Bf_eq2D      => default_get_Bf
  get_JBf_eq2D     => default_get_JBf
  get_Psi          => default_get_Psi
  get_DPsi         => default_get_DPsi
  get_domain       => default_get_domain
  get_pressure     => default_pressure
  export_boundary  => null()
  equilibrium_info => null()
  post_setup_equilibrium        => null()
  call initialize_magnetic_axis()


! 1. read user configuration
  rewind (iu)
  read   (iu, Equilibrium_Input, end=1000)
  iconfig = 1
  write (6, *)
  write (6, 1001)

! 1.b find equilibrium type
  select case(Data_Format)
  case (S_GEQDSK)
     i_equi   = EQ_GEQDSK
  case (S_GEQDSK_FREE)
     i_equi   = EQ_GEQDSK
     i_format = FREE
  case (S_DIVAMHD)
     i_equi   = EQ_DIVAMHD
  case (S_SONNET)
     i_equi   = EQ_SONNET
  case (S_M3DC1)
     i_equi   = EQ_M3DC1
  case (S_AMHD)
     i_equi   = EQ_AMHD
  case ('')
     i_equi   = EQ_GUESS
  case default
     write (6, *) 'error: ', Data_Format, ' is not a valid equilibrium type!'
     stop
  end select


! set default values
  use_boundary = (use_boundary .and. use_boundary_from_equilibrium)


! 2. load equilibrium data (if provided) ...............................
! determine equilibrium type (if not provided) .........................
  filename = trim(Prefix)//Data_File
  if (i_equi == EQ_GUESS .and. Data_File .ne. '') i_equi = get_equilibrium_format(filename)
! ... determine equilibrium type (done) ................................

  select case(i_equi)
  case (EQ_M3DC1)
     if (.not.m3dc1_loaded()) then
        write (6, *) 'error: M3D-C1 data not loaded, cannot set up equilbrium!'
        stop
     endif

  case (EQ_AMHD)
     call amhd_load (iu, iconfig, Ip, Bt, R0)
     if (M%R_estimate <= 0.d0) M%R_estimate = R0
     if (Xp(1)%R_estimate <= 0.d0) call amhd_get_Xp1(Xp(1)%R_estimate, Xp(1)%Z_estimate)

  case default
     if (Data_File == '') then
        select case(i_equi)
        case(EQ_GEQDSK, EQ_SONNET, EQ_DIVAMHD)
           write (6, *) 'error: data file required for data format ', Data_Format
           stop
        case default
           use_boundary = .false.
        end select
     else
        call load_equilibrium_data(filename, ierr)
        if (ierr > 0) stop
     endif
  end select
  call setup_equilibrium()


! 3. setup magnetic axis ...............................................
  call setup_magnetic_axis()


! 4. set up X-points ...................................................
  call setup_xpoints()


! 5. set dependent variables ...........................................
  ! pol. magn. flux at separatrix
  if (R_sepx > 0.d0) then
     r(1) = R_sepx
     r(2) = Z_sepx
     r(3) = 0.d0
     Psi_sepx = get_Psi(r)
  endif
  write (6, 5000) Psi_sepx


  ! 6. post-setup
  if (associated(post_setup_equilibrium)) call post_setup_equilibrium(Psi_axis, Psi_sepx)


  return
 1000 iconfig  = 0
  use_boundary = .false.
 1001 format (3x,'- Equilibrium configuration:')
 5000 format (8x,'Psi_sepx = ', e12.4)
  end subroutine load_equilibrium_config
!=======================================================================



!=======================================================================
  subroutine load_equilibrium_data(filename, ierr)
  use run_control, only: Prefix
  use geqdsk
  use divamhd

  character(len=*), intent(in)  :: filename
  integer,          intent(out) :: ierr


! load equilibrium data
  ierr = 0
  select case (i_equi)
  case (EQ_GEQDSK)
     call geqdsk_load (filename, Ip, Bt, Current_Fix, Psi_axis, Psi_sepx, Header_Format=i_format)
  case (EQ_DIVAMHD)
     call divamhd_load (filename, Ip, Bt, R0)

  case (EQ_SONNET)
     allocate (t_sonnet :: Bequi)
     select type(Bequi)
     class is (t_sonnet)
        call Bequi%load(filename, ierr)
        if (ierr > 0) return

        call Bequi%info(R0)
     end select
     if (M%R_estimate <= 0.d0) M%R_estimate = R0

  case default
     write (6, *) 'error: equilibrium format undefined!'
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
  use amhd

  real(real64)         :: Rbox(2), Zbox(2)


  ! select case equilibrium
  select case (i_equi)
  case (EQ_GEQDSK)
     get_Bf_eq2D                   => geqdsk_get_Bf
     get_JBf_eq2D                  => geqdsk_get_JBf
     get_Psi                       => geqdsk_get_Psi
     get_DPsi                      => geqdsk_get_DPsi
     get_pressure                  => geqdsk_get_pressure
     get_domain                    => geqdsk_get_domain
     export_boundary               => geqdsk_export_boundary
     broadcast_equilibrium         => geqdsk_broadcast
     equilibrium_info              => geqdsk_info
  case (EQ_DIVAMHD)
     get_Bf_eq2D                   => divamhd_get_Bf
     get_Psi                       => divamhd_get_Psi
     get_DPsi                      => divamhd_get_DPsi
     get_domain                    => divamhd_get_domain
     use_boundary                  = .false.
     broadcast_equilibrium         => divamhd_broadcast
  case (EQ_SONNET)
     get_Bf_eq2D                   => get_Bf_WRAPPER
     get_JBf_eq2D                  => get_JBf_WRAPPER
     get_Psi                       => get_Psi_WRAPPER
     get_DPsi                      => get_DPsi_WRAPPER
     get_domain                    => get_domain_WRAPPER
     use_boundary                  = .false.
  case (EQ_M3DC1)
     get_Bf_eq2D                   => m3dc1_get_Bf_eq2D
     get_Psi                       => m3dc1_get_Psi
     get_DPsi                      => m3dc1_get_DPsi
     use_boundary                  = .false.
     broadcast_equilibrium         => m3dc1_broadcast
  case (EQ_AMHD)
     get_Bf_eq2D                   => amhd_get_Bf
     get_Psi                       => amhd_get_Psi
     get_DPsi                      => amhd_get_DPsi
     get_pressure                  => amhd_get_pressure
     get_domain                    => amhd_get_domain
     use_boundary                  = .false.
     broadcast_equilibrium         => amhd_broadcast
     post_setup_equilibrium        => amhd_post_setup_equilibrium
  end select


  call get_domain(Rbox, Zbox)
  EQBox(1,:) = Rbox
  EQBox(2,:) = Zbox

  end subroutine setup_equilibrium
!=======================================================================



!=======================================================================
  subroutine setup_magnetic_axis()
  use run_control, only: Prefix

  real(real64)         :: r(3), x0(2)


  ! 3.1 find magnetic axis from estimated position
  if (M%R_estimate > 0.d0) then
     x0(1)    = M%R_estimate
     x0(2)    = M%Z_estimate
     M%X      = find_x(x0)

     r(1:2)   = M%X; r(3) = 0.d0
     M%Psi    = get_Psi(r)
     Psi_Axis = M%Psi
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
  write (6, *)

 3000 format (8x,'Psi_axis = ', e12.4)
  end subroutine setup_magnetic_axis
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
  if (i_equi == EQ_SONNET) then
     call Bequi%broadcast()
  else
     call broadcast_equilibrium()
  endif

  end subroutine broadcast_mod_equilibrium
!=======================================================================



!=======================================================================
  subroutine setup_xpoints()
  use run_control, only: Prefix, Debug

  integer, parameter :: iu = 32

  character(len=120) :: xpoint_file
  logical      :: ex
  real(real64) :: x(2)
  integer      :: ix, nx


  write (6, 4000)

  ! 1. pre-defined X-points
  xpoint_file = trim(Prefix)//'xpoints.dat'
  inquire (file=xpoint_file, exist=ex)

  ! 1.1. run subroutine initialize_equilibrium if file doesn't exist
  ! only for axisymmetric equilibrium with poloidal flux function
  select case(i_equi)
  case(EQ_GEQDSK, EQ_SONNET, EQ_DIVAMHD)
     if (.not.ex) call initialize_equilibrium()
     inquire (file=xpoint_file, exist=ex)
  case default
  end select

  ! 1.2. read pre-defined X-points
  if (ex) then
  open  (iu, file=xpoint_file)
  read  (iu, *) nx
  do ix=1,nx
     read (iu, *) x
     call Xp(ix)%setup(x)
     if (Xp(ix)%undefined) then
        write (6, 9001) ix;  stop
     endif

     ! update Psi_sepx
     if (ix == 1) Psi_sepx = Xp(1)%Psi

     ! screen output
     write (6, 4001) ix, Xp(ix)%X, Xp(ix)%PsiN()
  enddo
  close (iu)
  endif


  ! 2. user defined X-points
  do ix=1,nx_max
     if (Xp(ix)%R_estimate <= 0.d0) cycle

     ! set up X-point from estimate
     call Xp(ix)%setup()
     if (Xp(ix)%undefined) then
        write (6, 9002) ix;  stop
     endif

     ! update Psi_sepx
     if (ix == 1) Psi_sepx = Xp(1)%Psi

     ! screen output
     write (6, 4001) ix, Xp(ix)%X, Xp(ix)%PsiN()
     if (Debug) then
        write (6, 4002) Xp(ix)%H(1,1), Xp(ix)%H(1,2)
        write (6, 4002) Xp(ix)%H(2,1), Xp(ix)%H(2,2)
     endif
  enddo

 4000 format (3x,'- Configuring X-point(s): (R, Z, Psi)')
 4001 format(8x,i0,'. ',2f12.4,2x,e12.4)
 4002 format(8x,2e12.4)
 9001 format('error in subroutine setup_xpoints: cannot set up prefined X-point ',i0,'!')
 9002 format('error in subroutine setup_xpoints: cannot set up user defined X-point ',i0,'!')
  end subroutine setup_xpoints
!=======================================================================



!=======================================================================
! WRAPPER functions for Bequi (this is a temporary solution until all
! equilibrium configurations are implemented as derived type extended
! from t_bfield)
!=======================================================================
  function get_Bf_WRAPPER(r) result(Bf)
  real(real64), intent(in)  :: r(3)
  real(real64)              :: Bf(3)

  Bf = Bequi%get_Bf(r) * 1.d4 ! T -> Gauss

  end function get_Bf_WRAPPER
!=======================================================================
  function get_JBf_WRAPPER(r) result(JBf)
  real(real64), intent(in)  :: r(3)
  real(real64)              :: JBf(3,3)

  JBf = Bequi%get_JBf(r)

  end function get_JBf_WRAPPER
!=======================================================================
  function get_Psi_WRAPPER(r) result(Psi)
  real(real64), intent(in)  :: r(3)
  real(real64)              :: Psi

  Psi = Bequi%get_Psi(r)

  end function get_Psi_WRAPPER
!=======================================================================
  function get_DPsi_WRAPPER(r, mR, mZ) result(DPsi)
  real(real64), intent(in)  :: r(2)
  integer,      intent(in)  :: mR, mZ
  real(real64)              :: DPsi

  DPsi = Bequi%get_DPsi(get_r3(r), mR, mZ)

  end function get_DPsi_WRAPPER
!=======================================================================
  subroutine get_domain_WRAPPER (Rbox, Zbox)
  real(real64), intent(out) :: Rbox(2), Zbox(2)

  Rbox(1) = Bequi%Rmin;  Rbox(2) = Bequi%Rmax
  Zbox(1) = Bequi%Zmin;  Zbox(2) = Bequi%Zmax

  end subroutine get_domain_WRAPPER
!=======================================================================



!=======================================================================
  function get_r3(r)
  real(real64), dimension(:), intent(in) :: r
  real(real64)                           :: get_r3(3)


  select case(size(r))
  case(2)
     get_r3(1:2) = r
     get_r3(  3) = 0.d0
  case(3)
     get_r3      = r
  case default
     write (6, *) 'error in function get_r3: invalid size of argument r!'
     write (6, *) 'size(r) = ', size(r)
     stop
  end select

  end function get_r3
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
! Jacobian of magnetic field
!=======================================================================
  function default_get_JBf(r) result(J)
  real(real64), intent(in)  :: r(3)
  real(real64)              :: J(3,3)

  J = 0.d0

  end function default_get_JBf
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
  function get_PsiN(r) result(PsiN)
  real(real64), dimension(:), intent(in) :: r
  real(real64)                           :: PsiN

  real(real64) :: r3(3)


  select case(size(r))
  case(2)
     r3(1:2) = r
     r3(  3) = 0.d0
  case(3)
     r3      = r
  case default
     write (6, *) 'error in get_PsiN: invalid size of argument r0!'
     write (6, *) 'size(r) = ', size(r)
     stop
  end select

  PsiN = (get_Psi(r3) - Psi_axis) / (Psi_sepx - Psi_axis)

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
  real(real64), dimension(:), intent(in)  :: r
  integer, intent(in)                     :: nR, nZ
  real(real64)                            :: DPsiN

  real(real64) :: r2(2)


  select case(size(r))
  case(2,3)
     r2 = r(1:2)
  case default
     write (6, *) 'error in get_DPsiN: invalid size of argument r!'
     write (6, *) 'size(r) = ', size(r)
     stop
  end select

  DPsiN = get_DPsi(r2, nR, nZ) / (Psi_sepx - Psi_axis)

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
  function get_HN(x) result(HN)
  real(real64), intent(in) :: x(2)
  real(real64)             :: HN(2,2)


  HN = get_H(x) / (Psi_sepx - Psi_axis)

  end function get_HN
!=======================================================================
  function approximate_H(x) result(H)
  use run_control, only: Debug
  real(real64), intent(in) :: x(2)
  real(real64) :: H(2,2)

  real(real64), parameter :: gamma_1 = 1.d-1
  real(real64) :: xi(2), dfdxi(2), dfdyi(2)

  integer :: i, j


  H = 0.d0
  do j=1,2
     do i=1,2
        xi       = x
        xi(3-j)  = x(3-j) + (-1.d0)**i * gamma_1

        if (outside_domain(xi)) return
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
! default toroidal position: phi = 0
!=======================================================================
  function get_poloidal_angle(r) result(theta)
  real(real64), dimension(:), intent(in) :: r
  real(real64)                           :: theta

  real(real64) :: phi, Maxis(3)


  select case(size(r))
  case(2)
     phi = 0.d0
  case(3)
     phi = r(3)
  case default
     write (6, *) 'error in get_poloidal_angle: invalid size of argument r0!'
     write (6, *) 'size(r) = ', size(r)
     stop
  end select


  theta = 0.d0;   if (.not.associated(get_magnetic_axis)) return
  Maxis = get_magnetic_axis(phi)
  theta = atan2(r(2) - Maxis(2), r(1) - Maxis(1))

  end function get_poloidal_angle
!=======================================================================



!=======================================================================
! Return minor radius [cm] at r=(R,Z [cm], phi [rad])
!=======================================================================
  function get_rmin(r) result(rmin)
  real(real64), intent(in) :: r(3)
  real(real64)             :: rmin, Maxis(3)

  Maxis = get_magnetic_axis(r(3))
  rmin  = sqrt((r(1)-Maxis(1))**2 + (r(2)-Maxis(2))**2)

  end function get_rmin
!=======================================================================



!=======================================================================
! return local curvature of flux surface (in R-Z plane)
!=======================================================================
  function get_flux_surface_curvature(r) result(curv)
  real(real64), dimension(:), intent(in) :: r
  real(real64)                           :: curv

  real(real64) :: v(2), x1, y1, gradx1(2), grady1(2), x2, y2, kappa

  if (size(r) < 2  .or. size(r) > 3) then
     write (6, 9000)
     write (6, 9001) size(r)
     stop
  endif


  v(1)      = - get_DPsiN(r, 0, 1)
  v(2)      =   get_DPsiN(r, 1, 0)
  x1        = v(1)
  y1        = v(2)
  gradx1(1) = - get_DPsiN(r, 1, 1)
  gradx1(2) = - get_DPsiN(r, 0, 2)
  grady1(2) =   get_DPsiN(r, 2, 0)
  grady1(2) = - gradx1(1)
  x2        = sum(v * gradx1)
  y2        = sum(v * grady1)
  kappa     = (x1*y2 - y1*x2) / (x1**2 + y1**2)**1.5d0
  curv      = 1.d0/kappa

 9000 format('error in get_flux_surface_curvature:')
 9001 format('invalid size(r) = ', i0, '!')
  end function get_flux_surface_curvature
!=======================================================================



!=======================================================================
! Get cylindrical coordinates (R[cm], Z[cm], Phi[rad]) for flux
! coordinates (Theta[deg], PsiN, Phi[deg])
!
! ierr <= 0:	successfull operation
!	< 0:	at least one step with 1st order approximation
!	> 0:	exceed required accuracy
!
! optional output:
! iter:		number of iterations performed
!=======================================================================
  function get_cylindrical_coordinates(y, ierr, r0, order, damping, iter) result(r)
  use iso_fortran_env
  use math
  implicit none

  real(real64), intent(inout)         :: y(3)
  integer,      intent(out)           :: ierr
  real(real64), intent(in),  optional :: r0(3)
  real(real64)                        :: r(3)
  integer,      intent(in),  optional :: order
  real(real64), intent(in),  optional :: damping
  integer,      intent(out), optional :: iter

  integer, parameter :: imax = 100
  real(real64), parameter :: tolerance = 1.d-10

  real(real64) :: dpsi_dR, dpsi_dZ, dpsi_dl, dpsi_dl2, H(2,2), P, Q
  real(real64) :: M(3), Theta, PsiN, er(2), beta, delta, d

  integer :: i, o


  ierr  = 0


  ! set order of approximation
  o = 2
  if (present(order)) o = order
  if (o < 1  .or.  o > 2) then
     write (6, *) 'error in subroutine get_cylindrical_coordinates:'
     write (6, *) 'invalid parameter order = ', o
     stop
  endif


  ! set damping factor
  d = 1.d0
  if (present(damping)) d = damping


  ! set start point for approximation
  if (present(r0)) then
     r = r0
     write (6, *) 'warning: start point not supported in present implementation!'
  endif

  ! start near magnetic axis
  r(3)  = y(3) / 180.d0*pi
  M     = get_magnetic_axis(r(3))
  er(1) = cos(y(1)/180.d0*pi)
  er(2) = sin(y(1)/180.d0*pi)
  delta = 0.2d0 * length_scale()

  r(1:2)= M(1:2) + delta * er
  PsiN  = get_PsiN(r)


  do i=1,imax
     dpsi_dR = get_DPsiN(r, 1, 0)
     dpsi_dZ = get_DPsiN(r, 0, 1)

     ! Hessian
     H       = get_HN(r(1:2))

     beta    = y(2) - PsiN
     dpsi_dl = dpsi_dR*er(1) + dpsi_dZ*er(2)
     dpsi_dl2= H(1,1)*er(1)**2  +  2.d0*H(1,2)*er(1)*er(2)  +  H(2,2)*er(2)**2

     ! 1st order approximation
     delta   = beta / dpsi_dl * d
     Q       = 0.d0
     ! 2nd order approximation
     if (o == 2) then
        P       = - dpsi_dl / dpsi_dl2
        Q       = P**2 + 2.d0 * beta / dpsi_dl2
     endif
     if (Q > 0.d0) then
        delta = P - sign(sqrt(Q),P)
     else
        ierr = -1
        ! if (present(ierr)) then
        !    ierr = -1
        ! else
        !    write (6, *) 'error: ...
        !    stop
        ! endif
     endif

     r(1:2)  = r(1:2) + delta * er
     PsiN    = get_PsiN(r)
     if (PsiN > max(1.d0,y(2))) then
        r(1:2)  = r(1:2) - delta * er * (PsiN - max(1.d0,y(2))) / beta
     endif
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

  if (present(iter)) iter = i

  end function get_cylindrical_coordinates
!=======================================================================
  subroutine get_cylindrical_coordinates_error(ierr)
  integer, intent(in) :: ierr

  write (6, *) 'error in subroutine get_cylindrical_coordinates!'
  stop

  end subroutine get_cylindrical_coordinates_error
!=======================================================================



!=======================================================================
! from r0, go in gradPsiN direction to match PsiN_target
!
! return ierr <= 0:   correction step is successfull
!               -1:   at least one iteration with 1st order approximation
!                1:   required accuracy exceeded aver max. number of iterations
!                2:   required accuracy (ds) exceeded although iterations may be successful
!                3:   out of bounds
!=======================================================================
  function correct_PsiN(r0, PsiN_target, ierr, ds, delta_PsiN, iterations, debug) result(rc)
  use math
  use exceptions, only: OUT_OF_BOUNDS
  real(real64), intent(in)  :: r0(2)
  real(real64), intent(in)  :: PsiN_target
  integer,      intent(out) :: ierr
  real(real64)              :: rc(2)
  real(real64), intent(in),  optional :: ds, delta_PsiN
  integer,      intent(out), optional :: iterations
  integer,      intent(in),  optional :: debug

  integer, parameter                  :: imax = 100
  real(real64), parameter             :: damping = 0.2d0


  real(real64) :: PsiN, dPsiN, dPsi_dR, dPsi_dZ, ePsi(2), H(2,2), dPsi_dl, dPsi_dl2, P, Q, delta, tolerance
  integer      :: iter
  logical      :: debug_output


  ! 0. set up optional input
  ! 0.1 accuracy
  tolerance = 1.d-10
  if (present(delta_PsiN)) tolerance = delta_PsiN

  ! 0.2 debugging
  debug_output = .false.
  if (present(debug)) then
     debug_output = .true.
     if (debug < 10  .or.  debug > 99) then
        write (6, *) 'error: unit number for debugging must be in [10,99]!'
        stop
     endif
  endif


  ! 1. initialize
  OUT_OF_BOUNDS = .false.
  rc    = r0
  PsiN  = get_PsiN(rc)
  ierr  = 0
  if (debug_output) then
    write (debug, *) rc
  endif


  ! 2. iterative search for PsiN_target
  do iter=1,imax
     dPsiN    = PsiN_target - PsiN

     ! are we there yet?
     if (abs(dPsiN) < tolerance) exit

     ! if not, go in this direction
     dPsi_dR  = get_DPsiN(rc, 1, 0)
     dPsi_dZ  = get_DPsiN(rc, 0, 1)
     ePsi(1)  = dPsi_dR / sqrt(dPsi_dR**2 + dPsi_dZ**2)
     ePsi(2)  = dPsi_dZ / sqrt(dPsi_dR**2 + dPsi_dZ**2)

     ! Hessian
     H        = get_HN(rc)
     dPsi_dl  = dPsi_dR*ePsi(1) + dPsi_dZ*ePsi(2)
     dPsi_dl2 = H(1,1)*ePsi(1)**2  +  2.d0*H(1,2)*ePsi(1)*ePsi(2)  +  H(2,2)*ePsi(2)**2
     if (OUT_OF_BOUNDS) then
        ierr = 3
        return
     endif

     ! 1st order approximation
     delta    = dPsiN / dPsi_dl * damping
     ! 2nd order approximation
     P        = - dPsi_dl / dPsi_dl2
     Q        = P**2 + 2.d0 * dPsiN / dPsi_dl2
     if (Q > 0.d0) then
        delta = P - sign(sqrt(Q),P)
     else
        ierr  = -1
     endif

     if (debug_output) then
        write (debug, *) rc, delta, dPsi_dl, dPsi_dl2
     endif
     rc       = rc + delta * ePsi
     PsiN     = get_PsiN(rc)
  enddo


  ! 3. update output variables
  ! 3.1 exceed required accuracy (iteration limit)?
  if (iter > imax) ierr = 1
  ! output number of iterations
  if (present(iterations)) iterations = iter
  ! 3.2 exceed required accuracy (step size)
  if (present(ds)) then
     if (sqrt(sum((rc-r0)**2)) > ds) ierr = 2
  endif

  end function correct_PsiN
!=======================================================================



!=======================================================================
  function default_pressure(Psi) result(P)
  real(real64), intent(in) :: Psi
  real(real64)             :: P

  P = 0.d0

  end function default_pressure
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
  function outside_domain(x)
  real(real64), intent(in) :: x(2)
  logical                  :: outside_domain


  if (i_equi == EQ_SONNET) then
     outside_domain = Bequi%out_of_bounds(get_r3(x))
  endif

  outside_domain = .false.
  if (x(1) < EQBox(1,1)  .or.  x(1) > EQBox(1,2)  .or. &
      x(2) < EQBox(2,1)  .or.  x(2) > EQBox(2,2)) then
     outside_domain = .true.
  endif

  end function
!=======================================================================
  function leave_equilibrium_domain(r1, r2, X) result(leave)
  real(real64), dimension(:),        intent(in)  :: r1, r2
  real(real64), dimension(size(r1)), intent(out) :: X
  logical                                        :: leave

  real(real64) :: phi1, phi2, tv, th, t

  if (size(r1) .ne. size(r2)) then
     write (6, 9000)
     write (6, 9001)
     stop
  endif

  select case(size(r1))
  case(2)
     phi1 = 0.d0;  phi2 = 0.d0
  case(3)
     phi1 = r1(3); phi2 = r2(3)
  case default
     write (6, 9000)
     write (6, 9002) size(r1)
     stop
  end select

  leave = .false.
  tv    = 2.d0
  th    = 2.d0
  ! check left vertical boundary
  if (r1(1) > EQBox(1,1)  .and.  r2(1) < EQBox(1,1)) then
     tv    = (EQBox(1,1)-r1(1)) / (r2(1) - r1(1))
     leave = .true.
  endif
  ! check right vertical boundary
  if (r1(1) < EQBox(1,2)  .and.  r2(1) > EQBox(1,2)) then
     tv    = (EQBox(1,2)-r1(1)) / (r2(1) - r1(1))
     leave = .true.
  endif
  ! check lower horizontal boundary
  if (r1(2) > EQBox(2,1)  .and.  r2(2) < EQBox(2,1)) then
     th    = (EQBox(2,1)-r1(2)) / (r2(2) - r1(2))
     leave = .true.
  endif
  ! check upper horizontal boundary
  if (r1(2) < EQBox(2,2)  .and.  r2(2) > EQBox(2,2)) then
     th    = (EQBox(2,2)-r1(2)) / (r2(2) - r1(2))
     leave = .true.
  endif

  X = 0.d0
  if (leave) then
     t = min(tv, th)
     X(1:2) = r1(1:2) + t * (r2(1:2) - r1(1:2))
     X(3)   = phi1    + t * (phi2    - phi1)
  endif

 9000 format('error in leave_equilibrium_domain:')
 9001 format('arguments r1 and r2 must have the same size!')
 9002 format('invalide argument size ', i0)
  end function leave_equilibrium_domain
!=======================================================================



!=======================================================================
! boundary provided by equilibrium
!=======================================================================
  function equilibrium_provides_boundary() result(l)
  logical :: l


  l = use_boundary

  end function equilibrium_provides_boundary
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
     X = -1.d0
     return
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
     ! check definition domain boundary
     if (outside_domain(X)) then
        X = -1.d0
        return
     endif

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
  subroutine find_hyperbolic_points(nR, nZ, setup_Xpoints, file_output)
  use exceptions
  integer, intent(in)  :: nR, nZ
  logical, intent(in)  :: setup_Xpoints, file_output

  integer, parameter   :: iu = 54

  real(real64)         :: Rbox(2), Zbox(2)
  type(t_Xpoint) :: Xp0
  real(real64)   :: x(2), H(2,2), xh(nR*nZ, 2), r, lambda1, lambda2, v1(2), v2(2)
  real(real64)   :: DPsi, DPsi1, r3(3), xm(nR*nZ, 2)
  integer        :: i, j, k, indh, indm, ierr, iPsi


  call get_domain (Rbox, Zbox)
  write (6, 1000) Rbox, Zbox

  indh   = 0
  indm   = 0
  r3(3) = 0.d0
  if (file_output) then
     open  (iu, file='hyperbolic_points.dat')
     write (iu, 1001)
  endif
  write (6,  1002)
  loop2: do i=0, nR
  loop1: do j=0, nZ
     x(1) = Rbox(1) + (Rbox(2)-Rbox(1)) * i / nR
     x(2) = Zbox(1) + (Zbox(2)-Zbox(1)) * j / nZ
     x = find_X(x, Hout=H)

     if (x(1) < 0.d0) cycle ! not a valid critical point

     ! check if present critical point is identical to previous ones
     ! check all hyperbolic points
     do k=1,indh
        r = sqrt(sum((xh(k,:)-x)**2))
        if (r < 1.d-5) cycle loop1
     enddo
     ! check all minima/maxima
     do k=1,indm
        r = sqrt(sum((xm(k,:)-x)**2))
        if (r < 1.d-5) cycle loop1
     enddo

     ! so this is a new critical point, run analysis
     Xp0%X   = x;  Xp0%H = H
     r3(1:2) = x;  Xp0%Psi = get_Psi(r3);  Xp0%theta = get_poloidal_angle(r3)
     call Xp0%analysis(lambda1, lambda2, v1, v2, ierr)
     if (ierr .ne. 0) then
        indm       = indm + 1
        xm(indm,:) = x
        cycle ! this is not a hyperbolic point
     endif

     ! add present point to list
     indh = indh + 1
     xh(indh,:) = x
     write (6, 1003) indh, x, Xp0%PsiN(), lambda1, lambda2
     if (file_output) write (iu, *) x, Xp0%PsiN(), lambda1, lambda2

     if (setup_Xpoints) then
        if (indh > nx_max) then
           write (6, *) 'error: number of hyperbolic points exceeds limit!'
           stop
        endif

        Xp(indh)           = Xp0
        Xp(indh)%undefined = .false.
     endif
  enddo loop1
  enddo loop2
  if (file_output) close (iu)


  ! store minima/maxima
  if (file_output) then
  open  (iu, file='minima_and_maxima.dat')
  do k=1,indm
     x = xm(k,:)
     write (iu, *) x, get_PsiN(x)
  enddo
  close (iu)
  endif
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

  ! reset out of bounds flag
  call reset(OUT_OF_BOUNDS)

 1000 format(3x,'- Running search for hyperbolic points in domain: ',&
             '(',f0.2,' -> ',f0.2,') x (',f0.2,' -> ',f0.2,')')
 1001 format('# Hyperbolic points: R, Z, PsiN, Eigenvalues(l1, l2)')
 1002 format(13x,'(  R       ,  Z        )     PsiN        Eigenvalues')
 1003 format(8x,i0,4x,'(',f10.4,', ',f10.4,')',1x,f12.6,4x,'l1 = ',e12.4,',',4x,'l2 = ',e12.4)
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
! run analysis of Jacobian at X-point
!
! lambda1, v1	eigenvalue and eigenvector in unstable direction
! lambda2, v2   eigenvalue and eigenvector in stable direction
! eigenvectors are facing towards the magnetic axis
!
! ierr = 0: successfull
!        1: eigenvalues are non-real
!=======================================================================
  subroutine analysis(this, lambda1, lambda2, v1, v2, ierr)
  use math
  use run_control, only: Debug
  class(t_Xpoint) :: this
  real(real64), intent(out) :: lambda1, lambda2, v1(2), v2(2)
  integer,      intent(out) :: ierr

  real(real64) :: A(2,2), P, Q, phi(2), l(2), delta
  integer      :: i


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
     lambda1 = P
     lambda2 = Q
     ierr    = 1
     return
  endif
  ! calculate eigenvalues
  lambda1 = P + sqrt(Q);  l(1) = lambda1
  lambda2 = P - sqrt(Q);  l(2) = lambda2

  ! calculate eigenvectors
  do i=1,2
     phi(i) = atan2(-A(1,1) + l(i), A(1,2))
     !write (80, *) phi(i)

     ! find orientation with respect to magnetic axis (poloidal angle of X-point)
     delta  = phi(i) - this%theta
     ! make eigenvector facing towards magnetic axis
     if (abs(delta) < pi/2.d0  .or.  abs(delta) > 3.d0*pi/2.d0) then
        phi(i) = phi(i) - sign(pi, delta)
     endif
     !write (80, *) phi(i), delta, this%theta
  enddo

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
     open  (95, file='lambda.tmp')
     write (95, *) lambda1, lambda2
     close (95)
  endif

  end subroutine analysis
!=======================================================================



!=======================================================================
  function load(this, ierr) result(X)
  class(t_Xpoint) :: this
  real(real64)    :: X(2)
  integer, intent(out) :: ierr


  ierr = 0
  if (this%undefined) then
     X = -1.d0
     ierr = 1
  else
     X = this%X
  endif

  end function load
!=======================================================================



!=======================================================================
  subroutine setup(this, x_guess)
  class(t_Xpoint)          :: this
  real(real64), intent(in), optional :: x_guess(2)

  real(real64) :: x(2), r(3)


  x(1) = this%R_estimate;  x(2) = this%Z_estimate
  if (present(x_guess)) x = x_guess
  this%X = find_x(x, Hout=this%H)
  if (this%X(1) > 0.d0) then
     this%undefined = .false.

     r(1:2)         = this%X;  r(3) = 0.d0
     this%Psi       = get_Psi(r)
     this%theta     = get_poloidal_angle(this%X)
  endif

  end subroutine setup
!=======================================================================



!=======================================================================
  function t_Xpoint__PsiN(this) result(PsiN)
  class(t_Xpoint) :: this
  real(real64)    :: PsiN


  PsiN = (this%Psi - Psi_axis) / (Psi_sepx - Psi_axis)

  end function t_Xpoint__PsiN
!=======================================================================




!=======================================================================
! calculate eigenvectors v1,v2 of Hessian matrix of pol. magn. flux at x
! v1 points in ascending PsiN direction to the right
! v2 points in descending PsiN direction upwards
!=======================================================================
  subroutine H_eigenvectors (H, v1, v2)
  real(real64), intent(in)  :: H(2,2)
  real(real64), intent(out) :: v1(2), v2(2)

  real(real64) :: r(3), psi_xx, psi_xy, psi_yy, l1, l2, ac2, ac4, b2


  psi_xx = H(1,1) / (Psi_sepx-Psi_axis)
  psi_xy = H(1,2) / (Psi_sepx-Psi_axis)
  psi_yy = H(2,2) / (Psi_sepx-Psi_axis)


  ! get eigenvalues l1,l2 of Hessian at X-point
  ac2 = 0.5d0  * (psi_xx + psi_yy)
  ac4 = 0.25d0 * (psi_xx - psi_yy)**2
  b2  = psi_xy**2
  l1  = ac2 + dsqrt(ac4 + b2)
  l2  = ac2 - dsqrt(ac4 + b2)


  ! construct normalized eigenvectors
  ! ISSUE: this might not work if the X-point is straight below the magnetic axis!
  v1(1) = 1.d0
  v1(2) = - (psi_xx - l1) / psi_xy
  v1    = v1 / sqrt(sum(v1**2))

! construct v2 so that it is pointing upwards
  v2(2) = 1.d0
  v2(1) = - psi_xy / (psi_xx - l2)
  v2    = v2 / sqrt(sum(v2**2))

  end subroutine H_eigenvectors
!=======================================================================






!=======================================================================
!  function pol_flux(r) result(psi)
!  real*8, intent(in) :: r(3)
!  real*8             :: psi
!
!  end function pol_flux
!=======================================================================


!=======================================================================
!  subroutine Bf_pol_sub (n, s, y, f)
!  integer, intent(in) :: n
!  real*8, intent(in)  :: s, y(n)
!  real*8, intent(out) :: f(n)
!
!  ! n = 2, y(1) = R, y(2) = Z
!  end subroutine Bf_pol_sub
!=======================================================================

end module equilibrium
