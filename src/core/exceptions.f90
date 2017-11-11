module exceptions
  implicit none

  logical, dimension(1), target :: EXCEPTION = .false.
  logical, pointer :: OUT_OF_BOUNDS => EXCEPTION(1)

  contains
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine raise_exception(flag)
  logical, intent(inout) :: flag

  flag = .true.
  end subroutine raise_exception
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  subroutine reset(flag)
  logical, intent(inout) :: flag

  flag = .false.
  end subroutine reset
  !---------------------------------------------------------------------

end module exceptions
