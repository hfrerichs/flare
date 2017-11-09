subroutine equi2d_geometry(n)
  use iso_fortran_env
  use equilibrium, only: get_cylindrical_coordinates
  use flux_surface_2D
  use parallel
  implicit none
  integer, intent(in) :: n

  type(t_flux_surface_2D)  :: F


  integer, parameter  :: iu = 97

  real(real64) :: y(3), r(3), PsiN(n), area(n), dPsi(n)
  integer      :: i, ierr


  if (firstP) then
     write (6, *) 'Generating flux surface contours inside separatrix'
  else
     return
  endif


  open  (iu, file='flux_surfaces.geo')
  do i=1,n
     write (6, *) i, ' / ', n
     PsiN(i) = (i-0.5d0) / n

     ! find coordinates in real space for point i
     y(1) = 180.d0
     y(2) = PsiN(i)
     y(3) = 0.d0
     r    = get_cylindrical_coordinates(y, ierr)
     if (ierr > 0) write (6, *) 'warning: reference point at PsiN = ', y(2), ' exceeds required accuracy!'

     ! generate flux surface
     call F%generate_closed(r(1:2), RIGHT_HANDED)
     call F%setup_length_sampling()
     call F%surface_analysis(area(i), dPsi(i))
     call F%plot(iu);   write (iu, *)
  enddo
  close (iu)


  open  (iu, file='flux_surfaces.dat')
  do i=1,n
     write (iu, *) PsiN(i), area(i), dPsi(i)
  enddo
  close (iu)

end subroutine equi2d_geometry
