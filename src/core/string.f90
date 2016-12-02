!===============================================================================
! string related functions
!===============================================================================
module string
  implicit none

  contains
!=======================================================================



!=======================================================================
! parse string str for n-th sub-string delimited by "," or delimiter
!=======================================================================
  function parse_string(str, n, delimiter)
  character(len=*), intent(in)  :: str
  integer,          intent(in)  :: n
  character(len=1), intent(in), optional :: delimiter
  character(len=len_trim(str))  :: parse_string

  character(len=1) :: X
  integer :: i, i1_tmp(n+1), nstr, i1, i2


  X = ','
  if (present(delimiter)) X = delimiter


  parse_string = ''
  ! initialize upper and lower indices
  nstr      = len_trim(str)
  i1_tmp(1) = 1
  i2        = nstr
  do i=1,n
     ! get upper index of i-th sub-string
     i2 = i1_tmp(i) + scan(str(i1_tmp(i):nstr), X) - 2
     ! adjust upper index if no further delimiter is found
     if (i2 == i1_tmp(i) - 2) i2 = nstr

     ! set new lower index for next sub-string
     i1_tmp(i+1) = i2 + 2
     ! adjust new lower index if upper index is already at the end
     if (i1_tmp(i+1) > nstr) exit
  enddo
  i1 = i1_tmp(n)

  parse_string = str(i1:i2)

  end function parse_string
!=======================================================================



!=======================================================================
! return number of commands (sub-strings delimited by ";")
!=======================================================================
  function get_commands(string) result(n)
  character(len=*), intent(in) :: string
  integer                      :: n

  character(len=len(string))   :: stmp

  n = 0
  do
     stmp = parse_string(string, n+1, ';')
     if (stmp == '') exit

     n = n + 1
  enddo

  end function get_commands
!=======================================================================



!=======================================================================
  subroutine read_command(string, n, command, argument)
  character(len=*),           intent(in)  :: string
  integer,                    intent(in)  :: n
  character(len=len(string)), intent(out) :: command, argument

  character(len=len(string)) :: stmp


  stmp     = parse_string(string, n, ';')
  command  = 'undefined'
  argument = 'undefined'
  read (stmp, *) command
  if (len_trim(command)+2 <= len_trim(stmp)) then
     argument = stmp(len_trim(command)+2:len_trim(stmp))
  endif

  end subroutine read_command
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
