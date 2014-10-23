  recursive subroutine quicksort(A, n, m, i1, i2, jpivot)
  use iso_fortran_env
  implicit none

  real(real64), dimension(n,m) :: A
  integer, intent(in)          :: n, m, i1, i2, jpivot

  integer :: ip


  if (i1 < i2) then
     ip = partition(A, n, m, i1, i2, jpivot)
     call quicksort(A, n, m, i1, ip - 1, jpivot)
     call quicksort(A, n, m, ip, i2,     jpivot)
  endif

  contains

  function partition(A, n, m, i1, i2, jpivot) result(ip)
  use iso_fortran_env
  implicit none

  real(real64), dimension(n,m) :: A
  integer, intent(in)          :: n, m, i1, i2, jpivot
  integer                      :: ip

  real(real64) :: tmp(m), x
  integer :: i, j


  x = A(i1, jpivot)
  i = i1 - 1
  j = i2 + 1

  do
     j = j-1
     do
        if (A(j,jpivot) <= x) exit
        j = j-1
     enddo
     i = i+1
     do
        if (A(i,jpivot) >= x) exit
        i = i+1
     enddo
     if (i < j) then
        ! exchange A(i) and A(j)
        tmp    = A(i,:)
        A(i,:) = A(j,:)
        A(j,:) = tmp
     elseif (i == j) then
        ip = i + 1
        return
     else
        ip = i
        return
     endif
  enddo

  end function partition

  end subroutine quicksort
