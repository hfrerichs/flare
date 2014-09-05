!===============================================================================
!===============================================================================
module field_aligned_grid
  use iso_fortran_env
  implicit none

  real(real64), parameter :: &
     epsilon_r64 = epsilon(real(1.0,real64))


  character(len=*), parameter :: &
     LSN    = 'LSN', &
     DDN    = 'DDN', &
     SIMPLE = 'SIMPLE'
  
  character(len=*), parameter :: &
     EXACT = 'EXACT', &
     QUASI = 'QUASI'

 
  integer, parameter :: &
     N_block_max = 360         ! Maximum number of toroidal blocks


  ! Type of innermost flux surface (exact or quasi)
  character(len=72) :: &
     Layout = SIMPLE, &
     Innermost_Flux_Surface = EXACT


  integer :: &
     N_sym = 1, &              ! Toroidal symmetry (i.e. 5 => 72 deg)
     N_block = 1               ! Number of toroidal blocks

  real(real64) :: &
     x_in1(3)                    = 0.d0, &  ! reference points (R[cm], Z[cm], phi[deg]) ...
     x_in2(3)                    = 0.d0, &  ! ... on 1st and 2nd innermost flux surfaces
     Phi0                        = 0.d0, &  ! lower boundary of simulation domain
     Delta_Phi(N_block_max,-1:1) = 0.d0     ! user defined (non-default) toroidal block width [deg]


  ! Stellarator symmetry for 1st base grid
  logical :: &
     stellarator_symmetry = .false.


  namelist /Basic_Input/ &
     Layout, N_sym, N_block, x_in1, x_in2, Innermost_Flux_Surface, &
     Phi0, Delta_Phi

! internal variables
  logical :: default_decomposition

  real(real64) :: Delta_Phi_Sim, &
     Phi_base(N_block_max)       = 0.d0  ! user defined (non-default) base grid locations [deg]

  contains
!=======================================================================


!=======================================================================
  subroutine load_usr_conf

  integer, parameter :: iu = 12

  integer :: i


  ! 1. read from input file
  open  (iu, file='grid.conf', err=9000)
  read  (iu, Basic_Input, err=9000)
  close (iu)


  ! 2. setup internal variables
  ! default decomposition (if 1st Delta_Phi = 0)
  default_decomposition = (Delta_Phi(1,-1) + Delta_phi(1,1) == 0.d0)
  Delta_Phi_Sim         = real(360, real64) / N_sym

  ! set block size and position of base grids
  write (6, *)
  if (default_decomposition) then
     Phi0 = -0.5d0 * Delta_Phi_Sim / N_block
     write (6, 1000) Phi0, Phi0+Delta_Phi_Sim

     Delta_Phi = Delta_Phi_Sim / N_block / 2.d0
     do i=1,N_block
        Phi_base(i) = Delta_Phi_Sim / N_block * (i-1)
     enddo
  else
     write (6, 1001) Phi0, Phi0+Delta_Phi_Sim

     Phi_base(1) = Phi0 + Delta_Phi(1,-1)
     do i=2,N_block
        Phi_base(i) = Phi_base(i-1) + Delta_Phi(i-1,1) + Delta_Phi(i,-1)
     enddo
  endif

  ! output to screen
  write (6, 1002)
  do i=1,N_block
     write (6, 1003) i, Phi_base(i), Phi_base(i)-Delta_Phi(i,-1), Phi_base(i)+Delta_Phi(i,1)
  enddo
 1000 format (3x,'- Default decomposition of simulation domain (',f7.3,' -> ',f7.3,' deg):')
 1001 format (3x,'- User defined decomposition of simulation domain (',f7.3,' -> ',f7.3,' deg):')
 1002 format (8x,'block #, base location [deg], domain [deg]')
 1003 format (8x,      i7,5x,f7.3,':',5x,f7.3,' -> ',f7.3)


  ! 3. sanity check of user defined input
  call check_usr_conf()

  return
 9000 write (6, *) 'error while reading input file grid.conf!'
  stop
  end subroutine load_usr_conf
!=======================================================================



!=======================================================================
! sanity check for user defined input
!=======================================================================
  subroutine check_usr_conf

  real(real64) :: dPhi, f


  if (.not. default_decomposition) then
     ! Simulation domain must add up to 360 deg / N_sym
     dPhi = sum(Delta_Phi(1:N_block,-1)) + sum(Delta_Phi(1:N_block,1))
     f    = abs((dPhi - Delta_Phi_Sim) / Delta_Phi_sim)
     if (f > epsilon_r64) then
        write (6, *) "error: block sizes don't add up to size of simulation domain", &
                     ' (', Delta_Phi_Sim, ' deg)'
        stop
     endif
  endif

  end subroutine check_usr_conf
!=======================================================================



!=======================================================================
  subroutine generate_innermost_boundaries

  select case (Innermost_Flux_Surface)
  case (EXACT)
     call exact_surface (x_in1, 'lcfs0')
     call exact_surface (x_in2, 'lcfs1')
  case (QUASI)
     call quasi_surface (x_in1, 'lcfs0')
     call quasi_surface (x_in2, 'lcfs1')
  case default
     write (6, *) 'error: flux surface type ', trim(Innermost_Flux_Surface), ' not defined!'
     stop
  end select

  contains
!-----------------------------------------------------------------------
  subroutine exact_surface (x, s5)
  use run_control, only: N_mult, N_sym_RC => N_sym, N_points, x_start, Phi_output, Output_File

  real(real64),     intent(in) :: x(3)
  character(len=5), intent(in) :: s5

  character(len=3) :: sblock
  integer :: iblock


  ! set default number of points for Poincare plot
  if (N_points == 0) N_points = 1000

  N_sym_RC = N_sym
  x_start  = x
  if (default_decomposition) then
     N_mult   = N_block
     Output_File = s5//'.txt'
     if (N_block == 1) Output_File = s5//'_0.txt'
     call poincare_plot
  else
     N_mult   = 1
     do iblock = 1,N_block
        write (sblock, '(i3)') iblock-1
        Output_File = s5//'_'//trim(adjustl(sblock))//'.txt'
        Phi_output  = Phi_base(iblock)
        call poincare_plot
     enddo
  endif


  end subroutine exact_surface
!-----------------------------------------------------------------------
  subroutine quasi_surface (x, s5)
  real(real64),     intent(in) :: x(3)
  character(len=5), intent(in) :: s5

  write (6, *) 'generation of quasi flux surfaces not yet implemented!'
  stop
  end subroutine quasi_surface
!-----------------------------------------------------------------------
  end subroutine generate_innermost_boundaries
!=======================================================================

end module field_aligned_grid
