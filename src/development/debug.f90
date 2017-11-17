subroutine debug_subroutine(subroutine_name)
  use iso_fortran_env
  implicit none

  character(len=*), intent(in) :: subroutine_name


  write (6, 1000) subroutine_name
  select case(subroutine_name)
  case('left_hand_shift')
     call DEBUG_left_hand_shift()

  case default
     write (6, *) 'error: debugging of this subroutine is not implemented yet!'
     stop
  end select


 1000 format("Debugging subroutine '",a,"'")
end subroutine debug_subroutine





subroutine DEBUG_left_hand_shift()
  use iso_fortran_env
  use curve2D
  implicit none

  type(t_curve)      :: C
  character(len=256) :: filename
  real(real64)       :: dl


  write (6, *) 'enter name of data file: '
  read  (5, *) filename
  call C%load(filename)

  write (6, *) 'enter distance: '
  read  (5, *) dl

  write (6, *) 'enter name of output file: '
  read  (5, *) filename

  call C%left_hand_shift(dl)
  call C%plot(filename=filename)

end subroutine DEBUG_left_hand_shift
