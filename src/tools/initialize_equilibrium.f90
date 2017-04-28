subroutine initialize_equilibrium()
  use iso_fortran_env
  use run_control, only: N_R, N_Z, Debug, Prefix
  use equilibrium
  use boundary
  use separatrix
  use string
  implicit none

  integer, parameter :: iu = 54

  type(t_separatrix) :: S(nx_max), Stmp
  type(t_Xpoint)     :: Xtmp
  real(real64)       :: r(3)
  integer :: nR, nZ, ix, ib, icx(nx_max), nx


  ! 1. find hyperbolic points / X-points
  nR = 20;    if (N_R > 1) nR = N_R
  nZ = 20;    if (N_Z > 1) nZ = N_Z
  call find_hyperbolic_points(nR, nZ, .true.)


  ! 2. generate separatrix for each X-point
  write (6, 1020)
  do ix=1,nx_max
     if (Xp(ix)%undefined) cycle
     write (6, 1021) ix

     call S(ix)%generate_iX(ix)
     call S(ix)%plot(filename_prefix='S'//trim(str(ix)))
     !if (Debug) call S(ix)%plot(filename_prefix='S'//trim(str(ix)))
     !write (6, *) ix, S(ix)%connectB
  enddo


  ! 3. sort out irrelevant X-points
  icx = -1
  nx  = 0
  do ix=1,nx_max
     if (Xp(ix)%undefined) cycle

     ! X-point outside configuration boundary?
     r(1:2) = Xp(ix)%X;  r(3) = 0.d0
     if (outside_boundary(r)) cycle

     icx(ix) = 0
     nx      = nx + 1
  enddo
  ! check
  if (Debug) then
  write (6, *) 'relevant X-points:'
  do ix=1,nx_max
     if (icx(ix) == 0) write (6, *) ix, S(ix)%connectB
  enddo
  write (6, *)
  endif


  ! 4. find main X-point(x)
  open  (iu, file=trim(Prefix)//'xpoints.dat')
  write (iu, *) nx
  do ix=1,nx_max
     if (Xp(ix)%undefined) cycle

     do ib=1,4
        if (S(ix)%connectB(ib) > 0  .and.  icx(ix) == 0) then
           icx(ix) = 1
           write (iu, *) Xp(ix)%X
        endif
     enddo
  enddo
  ! check
  if (Debug) then
  write (6, *) 'main X-point(s):'
  do ix=1,nx_max
     if (icx(ix) == 1) write (6, *) ix, S(ix)%connectB
  enddo
  write (6, *)
  endif


  ! 5. output secondary X-points
  do ix=1,nx_max
     if (icx(ix) == 0) write (iu, *) Xp(ix)%X
  enddo
  close (iu)
  write (6, *)


  ! 6. set up X-points
  !call setup_xpoints()

 1020 format(3x,'- Generate separatrix for X-point')
 1021 format(8x,i0)
end subroutine initialize_equilibrium





subroutine import_equilibrium()
  use equilibrium_format
  use equilibrium

  integer, parameter :: iu = 99
  character(len=128) :: filename
  character(len=32)  :: s, fileformat


  call get_command_argument(2, filename)
  write (6, *) "Importing equilibrium from file ", trim(filename)


  call get_command_argument(3, s)
  select case(s)
  case('')
     i_equi = EQ_GUESS

!  case('format')
!     call get_command_argument(4, fileformat)

  case default
     write (6, *) 'error: invalid argument ', trim(s), '!'
     stop
  end select
  call load_equilibrium_data(filename)
  call setup_equilibrium()
  call setup_magnetic_axis()


  call initialize_equilibrium()
  open  (iu, file='bfield.conf')
  write (iu, 1000)
  write (iu, 1001) trim(filename)
  write (iu, 1002)
  close (iu)

 1000 format("&Equilibrium_Input")
 1001 format(2x,"Data_file    = '",a,"'")
 1002 format("/")
end subroutine import_equilibrium
