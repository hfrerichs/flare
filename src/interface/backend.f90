subroutine init()
  use run_control, only: load_run_control
  use bfield
  use boundary
  implicit none


  call load_run_control('run.conf')
  call setup_bfield_configuration()
  call setup_boundary()

end subroutine init
