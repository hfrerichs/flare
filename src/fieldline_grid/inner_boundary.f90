!module inner_boundary
!  use iso_fortran_env
!  implicit none
!
!  real(real64) :: &
!     Theta0, &                 ! Reference poloidal angle [deg]
!     x_guide1(3), &            ! Reference point (R,Z,phi) on inner guiding surface
!     x_guide2(3)               !                       ... on outer guiding surface
!
!
!  namelist /Simple_Layout/ &
!     Theta0, x_guide1, x_guide2
!
!
!end module inner_boundary


!=======================================================================
! Generate pair of innermost boundaries
!=======================================================================
  subroutine generate_innermost_boundaries
  use iso_fortran_env
  use fieldline_grid
  implicit none

  integer, parameter :: iu = 12

!  character(len=*), parameter :: &
!     EXACT = 'EXACT', &
!     QUASI = 'QUASI'


  character(len=72)  :: position, type

  real(real64) :: &
     x_in1(3)                    = (/120.d0, 0.d0, 0.d0/), &  ! reference points (R[cm], Z[cm], phi[deg]) ...
     x_in2(3)                    = (/119.d0, 0.d0, 0.d0/)     ! ... on 1st and 2nd innermost flux surfaces

  namelist /Inner_Boundary/ &
     x_in1, x_in2, position, type


  write (6, *)
  write (6, 1000)
  write (6, *)

  open  (iu, file=config_file, err=9000)
  read  (iu, Inner_boundary, end=9000)
  close (iu)


  select case (Innermost_Flux_Surface)
  case (SF_EXACT)
     call exact_surface (x_in1, 'fsin0')
     call exact_surface (x_in2, 'fsin1')
  case (SF_QUASI)
     call quasi_surface (x_in1, 'fsin0')
     call quasi_surface (x_in2, 'fsin1')
  case default
     write (6, *) 'error: flux surface type ', trim(Innermost_Flux_Surface), ' not defined!'
     stop
  end select
  write (6, *) 'finished generating innermost boundaries'

  return
 1000 format(3x,'- Generate inner simulation boundaries (based on Poincare plots)')
 9000 write (6, *) 'error while reading input file ', trim(config_file), '!'
  stop
  contains
!-----------------------------------------------------------------------
  subroutine exact_surface (x, s5)
  use run_control, only: N_mult, N_sym, N_points, x_start, Phi_output, Output_File

  real(real64),     intent(in) :: x(3)
  character(len=5), intent(in) :: s5

  character(len=3) :: sblock
  integer :: iblock


  ! set default number of points for Poincare plot
  if (N_points == 0) N_points = 1000

  N_sym    = symmetry
  x_start  = x
  if (default_decomposition) then
     N_mult   = blocks
     Output_File = s5//'.txt'
     if (blocks == 1) Output_File = s5//'_0.txt'
     call poincare_plot()
  else
     N_mult   = 1
     do iblock = 0,blocks-1
        write (sblock, '(i3)') iblock
        Output_File = s5//'_'//trim(adjustl(sblock))//'.txt'
        Phi_output  = Block(iblock)%phi_base
        call poincare_plot()
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
