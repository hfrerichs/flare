module exceptions
  implicit none

  logical, dimension(1), target :: EXCEPTION = .false.
  logical, pointer :: OUT_OF_BOUNDS => EXCEPTION(1)

end module exceptions
