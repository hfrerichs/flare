subroutine benchmark_mesh_spacing(ierr)
  use iso_fortran_env
  use mesh_spacing
  implicit none

  integer, intent(out) :: ierr

  integer, parameter   :: iu = 42

  type(t_spacing)      :: S
  character(len=256)   :: mode
  integer :: i, n, nsample


  open  (iu, file='mesh_spacing.ctrl')
  read  (iu, *) n, nsample
  do i=1,n
     read (iu, *) mode
     call S%init(mode)
     call S%plot(filename=trim(mode)//'.plt', nsample=nsample)
  enddo
  close (iu)

  ierr = 0

end subroutine benchmark_mesh_spacing
