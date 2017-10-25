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
  use run_control, only: Bfield_input_file
  use bfield
  implicit none

  Bfield_input_file = '/home/heinke/Database/Mag/DIII-D/mockup/bfield.conf'
  call setup_bfield_configuration()

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
