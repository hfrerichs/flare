!===============================================================================
! Interface for magnetic field data calculated with the M3D-C1 code
!
! Module fusion_io is required
!===============================================================================
module m3dc1
  use iso_fortran_env
#ifdef FIO
  use fusion_io
#endif
  implicit none

  private 

  ! time slice to read (0 = vacuum only, 1 = with response)
  integer :: timeslice = 1

  ! factor by which to multiply perturbed (linear) part of solution
  real*8  :: factor = 1.d0



  ! list of input files
  integer, parameter :: n_sets_max = 16

  integer :: n_sets = 0
  character(len=256) :: filename(n_sets_max)
  real*8       :: amplitude(n_sets_max) = 1.d0
  real*8       :: phase(n_sets_max)     = 0.d0	! offset in deg

  namelist /M3DC1_Input/ &
           timeslice, factor, &
           n_sets, filename, amplitude, phase

#ifdef FIO
  type(fio_search_hint), allocatable :: hint(:)
#endif
  integer :: isrcA(n_sets_max), imagA(n_sets_max), iEQ, iPsi

  logical, save :: loaded = .false.

  public :: &
     m3dc1_load, &
     m3dc1_loaded, &
     m3dc1_broadcast, &
     m3dc1_get_Bf, &
     m3dc1_get_Bf_eq2D, &
     m3dc1_get_Psi, &
     m3dc1_get_DPsi, &
     m3dc1_close

  contains
!===============================================================================



!===============================================================================
  subroutine m3dc1_load (iu, iconfig)
  use run_control, only: Prefix

  integer, intent(in)  :: iu
  integer, intent(out) :: iconfig

  character(len=256) :: Input_File
  integer :: ierr, i
  real*8  :: f

  rewind (iu)
  read   (iu, M3DC1_Input, end=1000)
  iconfig = 1
  write (6, *)
  write (6, 1001)
  write (6, 1002) timeslice
  write (6, 1003) factor
  write (6, 1004) n_sets

#ifdef FIO
  allocate (hint(n_sets))
  do i=1,n_sets
     ! Open M3D-C1 source
     Input_File = Prefix(1:len_trim(Prefix))//filename(i)
     call fio_open_source_f(FIO_M3DC1_SOURCE, trim(Input_file), isrcA(i), ierr)
     if(ierr.ne.0) stop

     ! Set source options
     call fio_get_options_f(isrcA(i), ierr)
     call fio_set_int_option_f(FIO_TIMESLICE, timeslice, ierr)


     ! Read equilibrium field
     if (i == 1) then
        call fio_set_int_option_f(FIO_PART, FIO_EQUILIBRIUM_ONLY, ierr)
        call fio_get_field_f(isrcA(i),FIO_VECTOR_POTENTIAL,iPsi,ierr)

        call fio_set_int_option_f(FIO_PART, FIO_EQUILIBRIUM_ONLY, ierr)
        call fio_get_field_f(isrcA(i),FIO_MAGNETIC_FIELD,iEQ,ierr)
     endif


     ! Read perturbation field
     call fio_set_int_option_f(FIO_PART, FIO_PERTURBED_ONLY, ierr)
     f = factor * amplitude(i)
     call fio_set_real_option_f(FIO_LINEAR_SCALE, f, ierr)
     call fio_get_field_f(isrcA(i),FIO_MAGNETIC_FIELD,imagA(i),ierr)


     ! initialize search hints
     call fio_allocate_search_hint_f(isrcA(i), hint(i), ierr)

     write (6, 1005) i, amplitude(i), phase(i)
  enddo
  loaded = .true.
#else
  write (6, *) 'error: FLARE compiled without M3D-C1 support!'
  stop
#endif

  return
 1000 iconfig = 0
 1001 format ('   - Magnetic field from M3D-C1:')
 1002 format (8x,'using time slice:             ',i4)
 1003 format (8x,'overall scale factor:         ',e11.4)
 1004 format (8x,'number of sub-sets:           ',i4)
 1005 format (8x,i4,' amplitude factor = ',f7.3,', phase [deg] = ',f7.3)
end subroutine m3dc1_load
!===============================================================================



!===============================================================================
  function m3dc1_loaded()
  logical :: m3dc1_loaded

  m3dc1_loaded = loaded

  end function m3dc1_loaded
!===============================================================================



!===============================================================================
  subroutine m3dc1_broadcast
  use parallel

  integer :: i


  if (nprs .gt. 1) then
     if (firstP) then
        write (6, *) 'parallel execution not supported for magnetic field data from M3D-C1!'
     endif
     stop
  endif

  call broadcast_inte_s (n_sets)
  call broadcast_inte_s (timeslice)
  call broadcast_real_s (factor)
  call broadcast_real   (amplitude, n_sets_max)
  call broadcast_real   (phase,     n_sets_max)
  do i=1,n_sets_max
     call broadcast_char(filename(i), 256)
  enddo

  return
  end subroutine m3dc1_broadcast
!===============================================================================



!===============================================================================
! Calculate R,phi,Z components of magnetic field vector [Gauss] at r=(R,Z [cm], phi [rad])
!===============================================================================
  function m3dc1_get_Bf(r) result(Bf)
  use math
  real*8, intent(in)  :: r(3)
  real*8              :: Bf(3)

  real*8  :: R3(3), B3(3)
  integer :: i, ierr


  Bf      = 0.d0
#ifdef FIO
  ! convert cm -> m
  R3(1)   = r(1) /1.d2
  R3(3)   = r(2) /1.d2

  do i=1,n_sets
     R3(2) = r(3) + phase(i) / 180.d0 * pi
     call fio_eval_field_f(imagA(i), R3, B3, ierr, hint=hint(i))
     if (ierr.gt.0) cycle

     !          (R,Z,phi)   (R,phi,Z)
     Bf(1)   =  Bf(1)     + B3(1)
     Bf(2)   =  Bf(2)     + B3(3)
     Bf(3)   =  Bf(3)     + B3(2)
  enddo
  ! convert T -> Gauss
  Bf = Bf * 1.d4
#endif

  end function m3dc1_get_Bf
!===============================================================================



!===============================================================================
  function m3dc1_get_Bf_eq2D(r) result(Bf)
  real(real64), intent(in)  :: r(3)
  real(real64)              :: Bf(3)

  real(real64) :: R3(3), B3(3)
  integer      :: ierr


  Bf      = 0.d0
#ifdef FIO
  ! convert cm -> m
  R3(1)   = r(1) /1.d2
  R3(3)   = r(2) /1.d2
  R3(2)   = 0.d0

  call fio_eval_field_f(iEQ, R3, B3, ierr, hint=hint(1))
  ! (R,Z,phi)   (R,phi,Z)
  ! convert T -> Gauss
  Bf(1)   =     B3(1) * 1.d4
  Bf(2)   =     B3(3) * 1.d4
  Bf(3)   =     B3(2) * 1.d4
#endif

  end function m3dc1_get_Bf_eq2D
!===============================================================================
  function m3dc1_get_Psi(r) result(Psi)
  real(real64), intent(in) :: r(3)
  real(real64) :: Psi

  real(real64) :: r3(3), B3(3)
  integer      :: ierr


  Psi     = 0.d0
#ifdef FIO
  ! convert cm -> m
  R3(1)   = r(1) /1.d2
  R3(3)   = r(2) /1.d2
  R3(2)   = r(3)

  call fio_eval_field_f(iPsi, R3, B3, ierr, hint=hint(1))
  Psi     = B3(2) * R3(1)
#endif

  end function m3dc1_get_Psi
!===============================================================================
  function m3dc1_get_DPsi(r, nR, nZ) result(DPsi)
  real(real64), intent(in) :: r(2)
  integer,      intent(in) :: nR, nZ
  real(real64) :: DPsi

  real(real64) :: r3(3), B3(3)
  integer      :: ierr


  DPsi    = 0.d0
  r3(1:2) = r
  r3(3)   = 0.d0
  B3      = m3dc1_get_Bf_eq2D(r3) / 1.d6 ! -> back to T*m

  if (nR == 1  .and.  nZ == 0) then
     DPsi =  B3(2) * (r3(1) / 1.d2)
  elseif (nR == 0  .and.  nZ == 1) then
     DPsi = -B3(1) * (r3(1) / 1.d2)
  endif

  end function m3dc1_get_DPsi
!===============================================================================



!===============================================================================
  subroutine m3dc1_close

  integer :: ierr, i

#ifdef FIO
  call fio_close_field_f(iEQ, ierr)
  do i=1,n_sets
     call fio_close_field_f(imagA(i), ierr)
     call fio_close_source_f(isrcA(i), ierr)
  enddo
  deallocate (hint)
#endif

  return
  end subroutine m3dc1_close
!===============================================================================

end module m3dc1
