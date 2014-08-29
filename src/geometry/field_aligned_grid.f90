!===============================================================================
!===============================================================================
module field_aligned_grid
  use iso_fortran_env
  implicit none

  character(len=*), parameter :: &
     LSN    = 'LSN', &
     DDN    = 'DDN', &
     SIMPLE = 'SIMPLE'
  
  character(len=*), parameter :: &
     EXACT = 'EXACT', &
     QUASI = 'QUASI'

 
  ! Type of innermost flux surface (exact or quasi)
  character(len=72) :: &
     Layout = SIMPLE, &
     Innermost_Flux_Surface = EXACT


  integer :: &
     N_sym = 1, &
     N_block = 1

  real(real64) :: &
     x_in1(3), &               ! reference points (R[cm], Z[cm], phi[deg])
     x_in2(3)                  ! on 1st and 2nd innermost flux surfaces


  ! Stellarator symmetry for 1st base grid
  logical :: &
     stellarator_symmetry = .false.


  namelist /Basic_Input/ &
     Layout, N_sym, N_block, x_in1, x_in2, Innermost_Flux_Surface

  contains
!=======================================================================


!=======================================================================
  subroutine load_usr_conf

  integer, parameter :: iu = 12


  open  (iu, file='grid.conf', err=9000)
  read  (iu, Basic_Input, err=9000)
  close (iu)

  return
 9000 write (6, *) 'error while reading input file grid.conf!'
  stop
  end subroutine load_usr_conf
!=======================================================================



!=======================================================================
  subroutine generate_innermost_boundaries

  select case (Innermost_Flux_Surface)
  case (EXACT)
     call exact_surface
  case (QUASI)
  case default
     write (6, *) 'error: flux surface type ', trim(Innermost_Flux_Surface), ' not defined!'
     stop
  end select

  contains
!-----------------------------------------------------------------------
  subroutine exact_surface
  use run_control, only: x_start

  integer :: iblock


  do iblock = 1,N_block
  enddo

  end subroutine exact_surface
!-----------------------------------------------------------------------
  end subroutine generate_innermost_boundaries
!=======================================================================

end module field_aligned_grid
