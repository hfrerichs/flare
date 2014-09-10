!===============================================================================
!===============================================================================
module mesh_spacing
  use iso_fortran_env
  implicit none
  private

  type, public :: t_spacing
     integer :: mode = 0

     contains
     procedure init, node
  end type t_spacing
  
  type(t_spacing), public, parameter :: Equidistant = t_spacing(0)

  contains
!=======================================================================
  


!=======================================================================
  subroutine init(this, mode)
  class(t_spacing) :: this
  character(len=*) :: mode


  ! equidistant spacing
  if (mode == '') then
     this%mode = 0
  endif

  end subroutine init
!=======================================================================



!=======================================================================
  function node(this, i, n) result(xi)
  class(t_spacing)    :: this
  integer, intent(in) :: i, n
  real(real64)        :: xi

  real(real64) :: t


  t  = 1.d0 * i / n
  xi = t
  end function node
!=======================================================================

end module mesh_spacing
