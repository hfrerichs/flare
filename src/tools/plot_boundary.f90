subroutine plot_boundary
  use boundary
  use run_control, only: Output_File, Phi_Output, Output_Format
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
  select case(Output_Format)
  case(1)
     do i=1,n_quad
        write (istr, '(i12)') i
        filename = 'boundary_quad_'//trim(adjustl(istr))//'.plt'
        call S_quad(i)%plot_at(Phi_Output, filename)
     enddo
  case(2)
     do i=1,n_quad
        write (istr, '(i12)') i
        filename = 'boundary_quad3D_'//trim(adjustl(istr))//'.plt'
        call S_quad(i)%plot(filename, 1)
     enddo
  case(3)
     do i=1,n_quad
        write (istr, '(i12)') i
        filename = 'boundary_quad3DX_'//trim(adjustl(istr))//'.plt'
        call S_quad(i)%plot(filename, 2)
     enddo
  case default 
     write (6, *) 'Output_Format = ', Output_Format, ' undefined!'
     stop
  end select

end subroutine plot_boundary
