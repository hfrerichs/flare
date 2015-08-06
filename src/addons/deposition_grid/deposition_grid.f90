program deposition_grid
  use curve2D
  use grid
  use iso_fortran_env
  implicit none

  integer, parameter :: iu = 42

  character(len=120) :: grid_file
  type(t_grid) :: G
  real(real64) :: phi, phi_start, phi_end, L
  integer      :: i, j, nt, np


  ! enter filename for initial footprint grid
  write (6, *) 'enter filename for initial footprint grid:'
  read  (5, *) grid_file

  ! read footprint grid
  call G%load(grid_file)
  if (G%coordinates .ne. CYLINDRICAL) then
     write (6, *) 'error: '//trim(grid_file)//' is not a cylindrical grid!'
     stop
  endif
  if (G%layout .ne. SEMI_STRUCTURED) then
     write (6, *) 'error: unepxected grid layout ', G%layout
     stop
  endif
  if (G%fixed_coord .ne. 3) then
     write (6, *) 'error: unexpected toroidal coordinate id ', G%fixed_coord
     stop
  endif



  ! setup toroidal boundary of simulation domain
  nt        = G%n2
  np        = G%n1
  phi_start = G%x2(1)
  phi_end   = G%x2(nt)

  ! write output
  open  (iu, file='output.dat')
  write (iu, FMT="(a,f6.2,a,f6.2)") "deposition: phi =", phi_start," -> ",phi_end
  write (iu, *)
  write (iu, FMT='(3(I5,1x),3(F10.5,1x),a)') nt, np, 1, &
           0.0,0.0,0.0,"  : NPoints_tor  NPoints_pol  NPeriod  DR  DZ  DL"

  do i=1,nt
     phi = G%x2(i)
     L   = 0.d0
     write (iu, *) phi
     do j=1,np
        if (j>1) L = L + sqrt(sum((G%x(j,1:2)-G%x(j-1,1:2))**2))
        write (iu, *) G%x(j,1:2), L
     enddo
  enddo
  close (iu)

end program deposition_grid
