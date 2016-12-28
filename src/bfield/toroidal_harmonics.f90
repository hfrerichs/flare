module toroidal_harmonics
  use iso_fortran_env
  use bspline2D
  implicit none
  private

  type, public :: t_harmonic
     ! toroidal mode number
     integer :: mode_number

     ! resolution in R and Z direction
     integer :: nr, nz

     real(real64) :: R_range(0:1), Z_range(0:1), scale_factor
     type(t_bspline2D) :: Bn(3,2)

     contains
     procedure :: load
     procedure :: get_Bf
     procedure :: broadcast
  end type t_harmonic

  type(t_harmonic), dimension(:), allocatable :: H
  integer :: n_harmonics = 0


  public :: &
     tor_harmonics_load, &
     tor_harmonics_broadcast, &
     tor_harmonics_get_Bf, &
     tor_harmonics_get_JBf

  contains
!===============================================================================



!===============================================================================
  subroutine load(this, filename, mode_number, scale_factor)
  class(t_harmonic)            :: this
  character(len=*), intent(in) :: filename
  integer,          intent(in) :: mode_number
  real(real64),     intent(in) :: scale_factor

  integer, parameter :: iu = 99

  character(len=80)  :: str

  real(real64), dimension(:,:), allocatable :: B_Rc, B_Rs, B_Zc, B_Zs, B_pc, B_ps
  real(real64), dimension(:,:), allocatable :: R, Z

  real(real64), dimension(:),   allocatable :: R1D, Z1D
  real(real64) :: A, A2, Atmp, dA, dR, dZ, v1(2), v2(2), b, b2, dB, btmp
  integer :: i, j, l, n, nr, nz


  open  (iu, file=filename)
  ! 1. read header
  read  (iu, 1000) str
  read  (iu, 1000) str
  read  (iu, 1000) str
  read  (str( 9:16), *) nr
  read  (str(21:28), *) nz
  this%nr = nr;  this%nz = nz
  write (6, *) this%nr, this%nz
  read  (iu, 1000) str
  read  (iu, 1000) str
 1000 format(a80)
  this%mode_number  = mode_number
  this%scale_factor = scale_factor


  ! 2. read data
  allocate (R(nr,nz), Z(nr,nz))
  allocate (B_Rc(nr,nz), B_Rs(nr,nz))
  allocate (B_Zc(nr,nz), B_Zs(nr,nz))
  allocate (B_Pc(nr,nz), B_Ps(nr,nz))
  do i=1,nr
  do j=1,nz
     read (iu, *) l, R(i,j), Z(i,j), B_Rc(i,j), B_Rs(i,j), &
                  B_Zc(i,j), B_Zs(i,j), B_Pc(i,j), B_Ps(i,j)
  enddo
  enddo
  close (iu)


  ! 3. setup bspline
  A  = 0.d0
  A2 = 0.d0
  n  = 0
  do i=1,nr-1
  do j=1,nz-1
     dR = 0.5d0*(R(i+1,j)-R(i,j) + R(i+1,j+1)-R(i,j+1))
     dZ = 0.5d0*(Z(i,j+1)-Z(i,j) + Z(i+1,j+1)-Z(i+1,j))
     dA = dR * dZ

     n    = n  + 1
     Atmp = A
     A    = A  + (dA - A) / n
     A2   = A2 + (dA - Atmp) * (dA - A)
  enddo
  enddo
  A2 = A2 / n
  write (6, *) "A = ", A, " +/- ", sqrt(A2)

  ! check orthogonality
  A  = 0.d0
  A2 = 0.d0
  b  = 0.d0
  b2 = 0.d0
  n  = 0
  do i=1,nr-1
  do j=1,nz-1
     v1(1) = 0.5d0 * (R(i+1,j)+R(i+1,j+1)-R(i,j)-R(i,j+1))
     v1(2) = 0.5d0 * (Z(i+1,j)+Z(i+1,j+1)-Z(i,j)-Z(i,j+1))
     v2(1) = 0.5d0 * (R(i,j+1)+R(i+1,j+1)-R(i,j)-R(i+1,j))
     v2(2) = 0.5d0 * (Z(i,j+1)+Z(i+1,j+1)-Z(i,j)-Z(i+1,j))

     n     = n + 1
     dA    = abs(v1(1)*v2(2) - v1(2)*v2(1))
     Atmp  = A
     A     = A  + (dA - A) / n
     A2    = A2 + (dA - Atmp) * (dA - A)

     dB    = sum(v1*v2)
     btmp  = b
     b     = b + (dB-B) / n
     b2    = b2 + (db - btmp) * (db - b)
  enddo
  enddo
  A2 = A2 / n
  b2 = b2 / n
  write (6, *) 'A_check = ', A, ' +/- ', sqrt(A2)
  write (6, *) 'b_check = ', b, ' +/- ', sqrt(b2)


  allocate (R1D(nr), Z1D(nz))
  do i=1,nr
     R1D(i) = R(i,1)
  enddo
  do j=1,nz
     Z1D(j) = Z(1,j)
  enddo

  call this%Bn(1,1)%new(nr, nz, 4, R1D, Z1D, B_Rc)
  call this%Bn(1,2)%new(nr, nz, 4, R1D, Z1D, B_Rs)
  call this%Bn(2,1)%new(nr, nz, 4, R1D, Z1D, B_Zc)
  call this%Bn(2,2)%new(nr, nz, 4, R1D, Z1D, B_Zs)
  call this%Bn(3,1)%new(nr, nz, 4, R1D, Z1D, B_Pc)
  call this%Bn(3,2)%new(nr, nz, 4, R1D, Z1D, B_Ps)
  do i=1,3
  do j=1,2
     call this%Bn(i,j)%setup()
  enddo
  enddo

  this%R_range(0) = R( 1, 1)
  this%R_range(1) = R(nr, 1)
  this%Z_range(0) = Z( 1, 1)
  this%Z_range(1) = Z( 1,nz)

  deallocate (R1D, Z1D)
  deallocate (R, Z, B_Rc, B_Rs, B_Zc, B_Zs, B_Pc, B_Ps)

  end subroutine load
!===============================================================================



!===============================================================================
  function get_Bf(this, x) result(Bf)
  class(t_harmonic)        :: this
  real(real64), intent(in) :: x(3)
  real(real64)             :: Bf(3)

  real(real64) :: B_Rc, B_Rs, B_Zc, B_Zs, B_pc, B_ps
  real(real64) :: x2(2), nphi, cosnphi, sinnphi


  x2   = x(1:2)
  nphi = this%mode_number * x(3)
  cosnphi = cos(nphi)
  sinnphi = sin(nphi)

  B_Rc = this%Bn(1,1)%eval(x2)
  B_Rs = this%Bn(1,2)%eval(x2)
  B_Zc = this%Bn(2,1)%eval(x2)
  B_Zs = this%Bn(2,2)%eval(x2)
  B_Pc = this%Bn(3,1)%eval(x2)
  B_Ps = this%Bn(3,2)%eval(x2)

  Bf(1) = B_Rc * cosnphi  +  B_Rs * sinnphi
  Bf(2) = B_Zc * cosnphi  +  B_Zs * sinnphi
  Bf(3) = B_Pc * cosnphi  +  B_Ps * sinnphi
  Bf    = Bf * this%scale_factor

  end function get_Bf
!===============================================================================



!===============================================================================
!  function get_JBf(r3) result(J)
!  class(t_harmonic)        :: this
!  real(real64), intent(in) :: x(3)
!  real(real64)             :: J(3,3)
!
!  real(real64) :: x2(2)
!
!
!  x2   = x(1:2)
!  nphi = this%n * x(3)
!  cosnphi = cos(nphi)
!  sinnphi = sin(nphi)
!  J    = 0.d0
!
!  end function get_JBf
!===============================================================================



!===============================================================================
  subroutine broadcast(this)
  use parallel
  class(t_harmonic) :: this

  integer :: i, j


  call broadcast_inte_s(this%mode_number)
  call broadcast_inte_s(this%nr)
  call broadcast_inte_s(this%nz)
  call broadcast_real  (this%R_range, 2)
  call broadcast_real  (this%Z_range, 2)
  call broadcast_real_s(this%scale_factor)
  do i=1,3;  do j=1,2
     call this%Bn(i,j)%broadcast()
  enddo;  enddo

  end subroutine broadcast
!===============================================================================



!===============================================================================
  subroutine tor_harmonics_load(iu, iconfig)
  use run_control, only: Prefix
  integer, intent(in)  :: iu
  integer, intent(out) :: iconfig

  integer, parameter :: MAX_HARMONICS = 128

  character(len=256) :: filename(MAX_HARMONICS)
  real(real64)       :: scale_factor(MAX_HARMONICS) = 1.d0
  integer            :: mode_number(MAX_HARMONICS)

  integer :: i


  namelist /Toroidal_Harmonics_Input/ &
     n_harmonics, mode_number, filename, scale_factor


  rewind (iu)
  read   (iu, Toroidal_Harmonics_Input, end=9000)
  if (n_harmonics > MAX_HARMONICS) then
     write (6, *) 'error: max. number of toroidal harmonics exceeded!'
     write (6, *) 'n_harmonics = ', n_harmonics, ' <= ', MAX_HARMONICS, ' required!'
     stop
  endif
  iconfig = 1
  if (n_harmonics <= 0) return
  write (6, *)
  write (6, 1000)
  write (6, 1001)


  allocate (H(n_harmonics))
  do i=1,n_harmonics
     write (6, 2000) mode_number(i), trim(filename(i))

     filename(i) = trim(Prefix)//filename(i)
     call H(i)%load(filename(i), mode_number(i), scale_factor(i))
  enddo

  return
 1000 format(3x,'- Toroidal harmonics:')
 1001 format(8x,'mode number')
 2000 format(8x,i0,4x,a)
 9000 iconfig = 0
  end subroutine tor_harmonics_load
!===============================================================================



!===============================================================================
  subroutine tor_harmonics_broadcast()
  use parallel

  integer :: i


  call broadcast_inte_s(n_harmonics)
  do i=1,n_harmonics
     call H(i)%broadcast()
  enddo

  end subroutine tor_harmonics_broadcast
!===============================================================================



!===============================================================================
  function tor_harmonics_get_Bf(x) result(Bf)
  real(real64), intent(in) :: x(3)
  real(real64)             :: Bf(3)

  real(real64) :: x_SI(3)
  integer :: i


  x_SI(1:2) = x(1:2) / 100.d0
  x_SI(  3) = x(  3)
  Bf = 0.d0
  do i=1,n_harmonics
     Bf = Bf + H(i)%get_Bf(x_SI) * 1.d4
  enddo

  end function tor_harmonics_get_Bf
!===============================================================================



!===============================================================================
  function tor_harmonics_get_JBf(x) result(JBf)
  real(real64), intent(in) :: x(3)
  real(real64)             :: JBf(3,3)

  integer :: i


  JBf = 0.d0
  do i=1,n_harmonics
     !JBf = JBf + H(i)%get_JBf(x)
  enddo

  end function tor_harmonics_get_JBf
!===============================================================================

end module toroidal_harmonics
