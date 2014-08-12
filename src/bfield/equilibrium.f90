!===============================================================================
! Equilibrium related functions and subroutines
!===============================================================================
module equilibrium
  use curve2D
  implicit none


!...............................................................................
! user defined parameters (to be set via configuration file)                   .

  character*120 :: &
     Data_File        = ''
  character*12  :: &
     Data_Format      = ''

  real*8 :: &
     R_axis_usr       = 0.d0, &        ! user defined position of magnetic axis
     Z_axis_usr       = 0.d0

  logical :: &
     use_boundary     = .true., &
     Current_Fix      = .true.

  integer :: &
     Diagnostic_Level = 0

  namelist /Equilibrium_Input/ &
     Data_File, Data_Format, use_boundary, Current_Fix, Diagnostic_Level, &
     R_axis_usr, Z_axis_usr
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
  integer :: &
     Bt_sign  = 0, &
     Ip_sign  = 0


  ! equilibrium type
  integer :: i_equi = -1


  ! Position of magnetic axis, poloidal magnetic flux at separatrix and magnetic axis
  real*8 :: R_axis, Z_axis, Psi_sepx, Psi_axis








!...............................................................................
! Interfaces for functions/subroutines from specific equilibrium types         .

  ! get equilibrium magnetic field in Cartesian coordinates
  procedure(default_get_Bf), pointer :: get_BCart_eq2D => default_get_Bf

  ! get equilibrium magnetic field in cylindrical coordinates
  procedure(default_get_Bf), pointer :: get_BCyl_eq2D  => default_get_Bf

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


  integer, parameter :: &
     EQ_GEQDSK = 1, &
     EQ_DIVA   = 2

  integer :: ipanic = 0

  contains
!=======================================================================



!=======================================================================
! Load equilibrium configuration
!=======================================================================
  subroutine load_equilibrium_config (iu, iconfig)
  use run_control, only: Prefix
  use geqdsk
  integer, intent(in)  :: iu
  integer, intent(out) :: iconfig


! read user configuration
  rewind (iu)
  read   (iu, Equilibrium_Input, end=1000)
  iconfig = 1
  write (6, *)
  write (6, 1001)
  Data_File = trim(Prefix)//Data_File


! set default values
  export_boundary => null()


! determine equilibrium type
   i_equi = EQ_GEQDSK


!...
!  select case(Data_Format)
!  case ('geqdsk')
!  case default
!     write (6, *) 'error: ', Data_Format, ' is not a valid equilibrium type!'
!     stop
!  end select


! load equilibrium data
  call geqdsk_load (Data_File, use_boundary, Current_Fix, Diagnostic_Level, R_axis, Z_axis, Psi_axis, Psi_sepx)
  call setup_equilibrium()


! set user defined values, if present
  if (R_axis_usr.gt.0.d0) then
     R_axis = R_axis_usr
     Z_axis = Z_axis_usr
  endif


! set dependent variables
  !psi_axis = pol_flux(magnetic_axis())
  !Psi_axis = get_Psi(magnetic_axis())
  !Psi_sepx = 


  return
 1000 iconfig = 0
 1001 format ('   - Axisymmetric (2D) MHD equilibrium:')
  end subroutine load_equilibrium_config
!=======================================================================



!=======================================================================
! Setup procedure pointers
!=======================================================================
  subroutine setup_equilibrium()
  use geqdsk

  ! select case equilibrium
  select case (i_equi)
  case (EQ_GEQDSK)
     get_BCart_eq2D                => geqdsk_get_BCart
     get_BCyl_eq2D                 => geqdsk_get_BCyl
     get_Psi                       => geqdsk_get_Psi
     get_DPsi                      => geqdsk_get_DPsi
     get_domain                    => geqdsk_get_domain
     equilibrium_provides_boundary => geqdsk_provides_boundary
     export_boundary               => geqdsk_export_boundary
     broadcast_equilibrium         => geqdsk_broadcast
  end select

  end subroutine setup_equilibrium
!=======================================================================



!=======================================================================
! Broadcast equilibrium data
!=======================================================================
  subroutine broadcast_mod_equilibrium()
  use parallel
  use geqdsk


  if (nprs == 1) return

  call broadcast_real_s (R_axis)
  call broadcast_real_s (Z_axis)
  call broadcast_real_s (Psi_axis)
  call broadcast_real_s (Psi_sepx)
  call broadcast_inte_s (i_equi)
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
! Sample (nR,nZ)-th derivative of poloidal magnetic flux at r=(R,Z [cm], phi [rad])
!=======================================================================
  function default_get_DPsi (r, nR, nZ) result(DPsi)
  real*8, intent(in)  :: r(3)
  integer, intent(in) :: nR, nZ
  real*8              :: DPsi

  DPsi = 0.d0
  if (ipanic > 0) then
     write (6, *) 'error: derivative of poloidal magnetic flux function not defined!'
     stop
  endif

  end function default_get_DPsi
!=======================================================================



!=======================================================================
! Return boundaries [cm] of equilibrium domain
!=======================================================================
  subroutine default_get_domain (Rbox, Zbox)
  real*8, intent(out) :: Rbox(2), Zbox(2)

  Rbox = 0.d0
  Zbox = 0.d0
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
  function pol_flux(r) result(psi)
  real*8, intent(in) :: r(3)
  real*8             :: psi

  end function pol_flux
!=======================================================================


!=======================================================================
! Provide position r=(R,Z,phi) [cm,cm,rad] of magnetic axis
! Optional input: phi (for non-axisymmetric equilibria)
! Default: phi = 0.0
!=======================================================================
  function magnetic_axis(phi) result(r)
  real*8, intent(in), optional :: phi
  real*8                       :: r(3)

  real*8 :: phi1


  ! scan optional arguments
  phi1 = 0.d0
  if (present(phi)) phi1 = phi
  r(3) = phi1

  ! set magnetic axis
  r(1) = R_axis
  r(2) = Z_axis


  end function magnetic_axis
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
