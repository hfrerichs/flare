subroutine plot_boundary
  use boundary
  use run_control, only: Output_File, Phi_Output
  implicit none

  character*120 :: filename
  character*12  :: istr
  integer :: i


  ! plot axisymmetric surfaces
  do i=1,n_axi
     write (istr, '(i12)') i
     filename = 'boundary_axi_'//trim(adjustl(istr))//'.plt'
     call S_axi(i)%plot(filename=filename)
  enddo


  ! plot profile of quadrilateral mesh
  do i=1,n_quad
     write (istr, '(i12)') i
     filename = 'boundary_quad_'//trim(adjustl(istr))//'.plt'
     call S_quad(i)%plot_at(Phi_Output, filename)
  enddo

end subroutine plot_boundary
