module equilibrium_format
  use iso_fortran_env
  implicit none

  character(len=*), parameter :: EQ_FORMAT(-1:5) = &
     (/'guess    ', &
       'undefined', &
       'geqdsk   ', &
       'divamhd  ', &
       'sonnet   ', &
       'm3dc1    ', &
       'amhd     '/)

  character(len=*), parameter :: &
     S_GEQDSK       = 'geqdsk', &
     S_GEQDSK_FREE  = 'geqdsk*', &
     S_DIVAMHD      = 'divamhd', &
     S_SONNET       = 'sonnet', &
     S_M3DC1        = 'm3dc1', &
     S_AMHD         = 'amhd'

  integer, parameter :: &
     EQ_GUESS       = -1, &
     EQ_UNDEFINED   = 0, &
     EQ_GEQDSK      = 1, &
     EQ_DIVAMHD     = 2, &
     EQ_SONNET      = 3, &
     EQ_M3DC1       = 4, &
     EQ_AMHD        = 5


  contains
!=======================================================================


!=======================================================================
  function get_equilibrium_format(filename) result(i_equi)
  character(len=*), intent(in) :: filename
  integer                      :: i_equi

  integer, parameter :: iu = 50

  character(len=80)  :: s, sformat
  logical :: ex


  i_equi = EQ_UNDEFINED
  inquire(file=filename, exist=ex)
  if (.not.ex) then
     write (6, *) 'error: equilibrium file ', trim(filename), ' does not exist!'
     stop
  endif


  sformat = ''
  open  (iu, file=filename)
  read  (iu, 1000) s;  if (s.ne.'') read  (s, *) sformat
  if (s(3:5) == 'TEQ'  .or.  sformat(1:4) == 'EFIT') then
     i_equi = EQ_GEQDSK
  elseif (s(5:11) == 'jm   :=') then
     i_equi = EQ_SONNET
  else
     read  (iu, 1000) s
     if (s(4:9) == 'File: ') then
        i_equi = EQ_DIVAMHD
     else
        i_equi = EQ_UNDEFINED
     endif
  endif
  close (iu)

 1000 format(a80)
  end function get_equilibrium_format
!=======================================================================

end module equilibrium_format
