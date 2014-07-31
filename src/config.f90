!===============================================================================
! Provide paths to magnetic field data and plasma facing components
!===============================================================================
module paths
  implicit none

  character*120 :: Prefix, &
                   Bfield_input_file

  contains
  !=====================================================================


  !=====================================================================
  subroutine initialize

  Bfield_input_file = trim(Prefix)
  end subroutine initialize
  !=====================================================================

end module paths
