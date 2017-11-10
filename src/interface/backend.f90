subroutine init(machine, configuration)
  use run_control, only: load_run_control, data_dir, Prefix, Bfield_input_file
  use bfield
  use boundary
  implicit none
  character(len=*), intent(in), optional :: machine, configuration


  if (present(machine) .and. present(configuration)) then
     Prefix = data_dir//'/'//trim(machine)//'/'//trim(configuration)//'/'
     Bfield_input_file = trim(Prefix)//'bfield.conf'
  else
     call load_run_control('run.conf')
  endif
  call setup_bfield_configuration()
  call setup_boundary()

end subroutine init
