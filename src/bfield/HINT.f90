!===============================================================================
! Interface for magnetic fields calculated by the HINT code
!===============================================================================
module HINT
  use iso_fortran_env
  implicit none

  private
  character(len=120) :: data_file = ''
  character(len=20)  :: data_format = 'mag'

  public :: HINT_load, &
            HINT_broadcast, &
            HINT_get_Bf, &
            HINT_get_JBf

  contains
  !=====================================================================



  !=====================================================================
  ! load configuration and setup related variables
  !=====================================================================
  subroutine HINT_load (iu, iconfig)
  use cylindrical_coord_mod, only: magset, mag_form
  use run_control, only: Prefix
  integer, intent(in)  :: iu
  integer, intent(out) :: iconfig

  namelist /Hint_Input/ &
     data_file, data_format


  rewind (iu)
  read   (iu, HINT_Input, end=1000)
  iconfig  = 1
  mag_form = data_format
  write (6, *)
  write (6, 1001)
  open  (25, file=trim(Prefix)//data_file, form='unformatted', status='old')
  call magset()
  close (25)

  return
 1000 iconfig = 0
 1001 format(3x,'- Magnetic field from HINT code')
  end subroutine HINT_load
  !=====================================================================



  !=====================================================================
  !=====================================================================
  subroutine HINT_broadcast()
  use parallel


  if (nprs ==1) return

  if (mype == 0) then
     write (6, *) 'error: parallel execution not supported for magnetic field from HINT!'
  endif
  stop

  end subroutine HINT_broadcast
  !=====================================================================



  !=====================================================================
  ! return magnetic field components using cylindrical coordinates (R, Z, phi)
  !=====================================================================
  function HINT_get_Bf(r) result(Bf)
  use kind_spec
  use cylindrical_coord_mod, only: mgval1, rminb, rmaxb, zminb, zmaxb
  use exceptions,            only: OUT_OF_BOUNDS
  real(real64), intent(in) :: r(3)
  real(real64)             :: Bf(3)

  real(dp) :: r1, phi, z, br, bphi, bz, bb


  r1    = r(1) / 100.d0
  z     = r(2) / 100.d0
  phi   = r(3)
  Bf    = 0.d0
  if (r1 < rminb  .or.  r1 > rmaxb  .or. &
      z  < zminb  .or.  z  > zmaxb) then
     OUT_OF_BOUNDS = .true.
     return
  endif


  call mgval1(r1, phi, z, br, bphi, bz, bb)
  Bf(1) = br
  Bf(2) = bz
  Bf(3) = bphi
  Bf    = Bf * 1.d4 ! convert T -> Gauss

  end function HINT_get_Bf
  !=====================================================================



  !=====================================================================
  ! return Jacobian of magnetic field [T/m] in cylindrical coordinates
  !=====================================================================
  function HINT_get_JBf(r) result(JBf)
  use kind_spec
  use spline_mod,            only: l3d
  use cylindrical_coord_mod, only: mgval2, rminb, rmaxb, zminb, zmaxb
  use exceptions,            only: OUT_OF_BOUNDS
  real(real64), intent(in) :: r(3)
  real(real64)             :: JBf(3,3)

  real(dp) :: r1, phi, z, b(l3d), dbdr(l3d), dbdp(l3d), dbdz(l3d)


  r1    = r(1) / 100.d0
  z     = r(2) / 100.d0
  phi   = r(3)
  JBf   = 0.d0
  if (r1 < rminb  .or.  r1 > rmaxb  .or. &
      z  < zminb  .or.  z  > zmaxb) then
     OUT_OF_BOUNDS = .true.
     return
  endif

  call mgval2(r1, phi, z, b, dbdr, dbdp, dbdz)
  JBf(1:3,1) = dbdr(1:3)
  JBf(1:3,2) = dbdz(1:3)
  JBf(1:3,3) = dbdp(1:3)

  end function HINT_get_JBf
  !=====================================================================

end module HINT
