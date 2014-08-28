subroutine run_control_development (Run_Type)
  use parallel
  implicit none

  character(len=120), intent(in) :: Run_Type


  select case (Run_Type)
  case ('hyperbolic_fixed_point')
     if (firstP) write (6, 1000)
     call hyperbolic_fixed_point
  case default
     write (6, *) 'run type "', trim(Run_Type), '" not defined!'
     stop
  end select

 1000 format ('WARNING: RUNNING SUB-PROGRAM IN DEVELOPMENT STAGE!')
end subroutine run_control_development
