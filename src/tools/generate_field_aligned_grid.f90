subroutine generate_field_aligend_grid (run_level)
  use parallel
  use fieldline_grid
  use emc3_grid
  use topo_sc
  use topo_lsn
  implicit none

  integer, intent(in) :: run_level

  procedure(), pointer :: make_base_grids
  logical :: level(8)


  if (firstP) then
     write (6, *) 'Generating fieldline grid ...'
     write (6, *) '... executing run level ', run_level
  endif


  call setup_grid_configuration()
  select case(topology)
  case(TOPO_SC, TOPO_SC1)
     call setup_topo_sc()
     make_base_grids => make_base_grids_sc
  case(TOPO_LSN, TOPO_LSN1)
     call setup_topo_lsn()
     make_base_grids => make_base_grids_lsn
  case default
     write (6, *) 'error: grid topology ', trim(topology), ' not supported!'
     stop
  end select


  ! select run level
  level = .false.
  if (run_level == 0) then
     level = .true.
  elseif (run_level > 0  .and.  run_level <=8) then
     level(run_level) = .true.
  else
     write (6, *) 'error: run level ', run_level, ' not implemented!'
     stop
  endif


  ! Level 1: generate pair of innermost boundaries
  if (level(1)) then
     call generate_innermost_boundaries()
  endif


  ! Level 2: generate base grids
  if (level(2)) then
     call make_base_grids()
  endif


  ! Level 3: generate 3D grid from field line tracing
  if (level(3)) then
     call setup_emc3_grid_layout()
     call trace_nodes()
     call write_emc3_grid()
     call write_emc3_input_files()
  endif


  ! Level 4: generate vacuum domain (used by EIRENE only)
  if (level(4)) then
     call load_emc3_grid()
     call vacuum_domain_for_EIRENE()
     call write_emc3_grid()
  endif


  ! Level 5: sample magnetic field strength on grid
  if (level(5)) then
     call load_emc3_grid()
     call sample_bfield_on_emc3_grid()
  endif


  ! Level 6: generate plate surfaces for EMC3
  if (level(6)) then
     call load_emc3_grid()
     call generate_plates()
  endif


  !---------------------------------------------

  ! Level 7: re-generate input files
  if (level(7)) then
     call write_emc3_input_files()
  endif


  ! Level 8: run grid checks
  if (level(8)) then
     call load_emc3_grid()
     call check_emc3_grid()
  endif


end subroutine generate_field_aligend_grid
