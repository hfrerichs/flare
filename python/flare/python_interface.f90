module types
  implicit none
  integer, parameter :: dp = selected_real_kind(15,307) ! double precision
end module types



subroutine init()
  use parallel
  implicit none

!  call initial_parallel()
  call print_version()

end subroutine init



subroutine init_bfield()
  use run_control, only: Bfield_input_file, Prefix
  use bfield
  use boundary
  implicit none

  !Bfield_input_file = '/home/heinke/Database/Mag/DIII-D/mockup/bfield.conf'
  !Bfield_input_file = '/home/heinke/Database/Mag/DIII-D/148712/axi/bfield.conf'
  Prefix = '/home/heinke/Database/Mag/DIII-D/148712/axi/'
  Bfield_input_file = trim(Prefix)//'bfield.conf'
  call setup_bfield_configuration()

  call setup_boundary()

end subroutine init_bfield



subroutine trace_bline_py(x, n)
  use types
  use run_control, only: x_start, N_steps
  implicit none

  real(dp), intent(in) :: x(3)
  integer,      intent(in) :: n


  x_start = x
  N_steps = n
  write (6, *) 'x, n = ', x, n
  call trace_bline()

end subroutine trace_bline_py



subroutine connection_length_py(n, x0, d)
  use types
  use fieldline
  use run_control, only: Trace_Step
  implicit none

  integer,  intent(in)  :: n
  real(dp), intent(in)  :: x0(n, 3)
  real(dp), intent(out) :: d(n, 3)

  type(t_fieldline) :: F
  real(dp) :: lc(-1:1)
  integer  :: i, idir


  do i=1,n
     do idir=-1,1,2
        call F%init(x0(i,:), idir*Trace_Step)
        call F%trace(400.d2, .true.)
        lc(idir) = F%lc
     enddo
     d(i,1) = abs(lc(-1)) + abs(lc(1))
     !d(i,1) = x0(i,1)
     d(i,2) = lc(-1)
     d(i,3) = lc( 1)
     write (6, *) x0(i,:), lc(-1), lc(1)
  enddo

end subroutine connection_length_py
