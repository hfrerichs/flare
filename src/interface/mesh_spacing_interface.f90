module mesh_spacing_interface
  use types
  implicit none

  contains
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine generate(cmd, n, output_file, output_format)
  use mesh_spacing
  character(len=*), intent(in) :: cmd, output_file, output_format
  integer,          intent(in) :: n

  type(t_spacing) :: M


  call M%init(cmd)
  call M%plot(filename=output_file, nsample=n, style=output_format)

  end subroutine generate
  !---------------------------------------------------------------------

end module mesh_spacing_interface
