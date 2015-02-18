!===============================================================================
! string related functions
!===============================================================================
module string
  implicit none

  contains
!=======================================================================



!=======================================================================
! parse string str for n-th sub-string delimited by ","
!=======================================================================
  function parse_string(str, n)
  character(len=*), intent(in)  :: str
  integer,          intent(in)  :: n
  character(len=len_trim(str))  :: parse_string

  integer :: i, i1_tmp(n+1), nstr, i1, i2


  ! initialize upper and lower indices
  nstr      = len_trim(str)
  i1_tmp(1) = 1
  i2        = nstr
  do i=1,n
     ! get upper index of i-th sub-string
     i2 = i1_tmp(i) + scan(str(i1_tmp(i):nstr), ',') - 2
     ! adjust upper index if no further delimiter is found
     if (i2 == i1_tmp(i) - 2) i2 = nstr

     ! set new lower index for next sub-string
     i1_tmp(i+1) = i2 + 2
     ! adjust new lower index if upper index is already at the end
     if (i1_tmp(i+1) > nstr) i1_tmp(i+1) = i1_tmp(i)
  enddo
  i1 = i1_tmp(n)

  parse_string = str(i1:i2)

  end function parse_string
!=======================================================================



!=======================================================================
! Convert an integer to string
!=======================================================================
  character(len=20) function str(k)
  integer, intent(in) :: k
  write (str, *) k
  str = adjustl(str)
  end function str
!=======================================================================

end module string
