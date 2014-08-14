!===============================================================================
!===============================================================================
module search
  implicit none

  contains
!=======================================================================



!=======================================================================
! search array A(0:n) for index "ind" with A(ind) <= key < A(ind+1)
!=======================================================================
  recursive function binary_interval_search(imin, imax, A, key, ierr) result(ind)
  integer, intent(in)    :: imin, imax
  real*8, intent(in)     :: A(imin:imax)
  real*8, intent(in)     :: key
  integer, intent(inout) :: ierr
  integer                :: ind

  integer :: imid


  ierr = 0
  ! test if array domain is empty
  if (imax <= imin) then
     ind  = -1
     ierr = 1
     return
  endif

  ! check boundaries
  if (key < A(imin)  .or.  key > A(imax)) then
     ind  = -1
     ierr = 1
     return
  endif


  ! only one interval left
  if (imax == imin+1) then
     ind  = imin
     return
  endif


#if defined(DEBUG)
  write (6, *) 'imin, imax = ', imin, imax
#endif
  ! calculate midpoint to cut array in half
  imid = imin + ((imax - imin) / 2)
  if (key < A(imid)) then
     ind = binary_interval_search (imin, imid, A(imin:imid), key, ierr)
  else
     ind = binary_interval_search (imid, imax, A(imid:imax), key, ierr)
  endif

  return
  end function binary_interval_search
!=======================================================================

end module search
