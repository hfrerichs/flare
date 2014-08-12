!===============================================================================
! Interface for magnetic field data calculated with the M3D-C1 code
!
! Module fusion_io is required
!===============================================================================
module m3dc1
#ifdef M3DC1
  use fusion_io
#endif
  implicit none

  private 

  ! time slice to read (0 = vacuum only, 1 = with response)
  integer :: itime = 1

  ! factor by which to multiply perturbed (linear) part of solution
  real*8  :: factor = 1.d0



  ! list of input files
  integer, parameter :: n_sets_max = 16

  integer :: n_sets = 0
  character*80 :: filename(n_sets_max)
  real*8       :: amplitude(n_sets_max) = 1.d0
  real*8       :: phase(n_sets_max)     = 0.d0	! offset in deg

  namelist /M3DC1_Input/ &
           itime, factor, &
           n_sets, filename, amplitude, phase

  integer :: isrcA(n_sets_max), imagA(n_sets_max)

  public :: &
     m3dc1_load, &
     m3dc1_broadcast, &
     m3dc1_get_Bf, &
     m3dc1_close

  contains
!===============================================================================



!===============================================================================
  subroutine m3dc1_load (iu, iconfig)
  use run_control, only: Prefix

  integer, intent(in)  :: iu
  integer, intent(out) :: iconfig

  character*80 :: Input_File
  integer :: ierr, i
  real*8  :: f

  rewind (iu)
  read   (iu, M3DC1_Input, end=1000)
  iconfig = 1
  write (6, *)
  write (6, 1001)
  write (6, 1002) itime
  write (6, 1003) factor

#ifdef M3DC1
  do i=1,n_sets
     ! Open M3D-C1 source
     Input_File = Prefix(1:len_trim(Prefix))//filename(i)
     call fio_open_source_f(FIO_M3DC1_SOURCE, trim(Input_file), isrcA(i), ierr)
     if(ierr.ne.0) stop

     ! Set options
     call fio_get_options_f(isrcA(i), ierr)
     call fio_set_int_option_f(FIO_TIMESLICE, itime, ierr)
     call fio_set_int_option_f(FIO_PART, FIO_PERTURBED_ONLY, ierr)
     f = factor * amplitude(i)
     call fio_set_real_option_f(FIO_LINEAR_SCALE, f, ierr)

     ! Read magnetic field
     call fio_get_field_f(isrcA(i),FIO_MAGNETIC_FIELD,imagA(i),ierr)
  enddo
#else
  write (6, *) 'error: FLARE compiled without M3D-C1 support!'
  stop
#endif

  return
 1000 iconfig = 0
 1001 format ('   - Magnetic field from M3D-C1:')
 1002 format ('        using time slice:             ',i4)
 1003 format ('        scale factor:                 ',e11.4)
end subroutine m3dc1_load
!===============================================================================



!===============================================================================
  subroutine m3dc1_broadcast
  use parallel

  if (nprs .gt. 1) then
     if (firstP) then
        write (6, *) 'parallel execution not supported for magnetic field data from M3D-C1!'
     endif
     stop
  endif

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
#ifdef M3DC1
  ! convert cm -> m
  R3(1)   = r(1) /1.d2
  R3(3)   = r(2) /1.d2

  do i=1,n_sets
     R3(2) = r(3) + phase(i) / 180.d0 * pi
     call fio_eval_field_f(imagA(i), R3, B3, ierr)
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
  subroutine m3dc1_close

  integer :: ierr, i

#ifdef M3DC1
  do i=1,n_sets
     call fio_close_field_f(imagA(i), ierr)
     call fio_close_source_f(isrcA(i), ierr)
  enddo
#endif

  return
  end subroutine m3dc1_close
!===============================================================================

end module m3dc1
