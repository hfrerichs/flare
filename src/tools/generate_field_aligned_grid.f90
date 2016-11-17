subroutine generate_field_aligend_grid (run_level, run_level_end)
  use parallel
  use fieldline_grid
  use emc3_grid
  use modtopo_sc
  use modtopo_lsn
  use modtopo_ddn
  use modtopo_cdn
  use modtopo_dsfp
  use modtopo_stel
  implicit none

  integer, intent(inout) :: run_level, run_level_end

  integer, parameter     :: max_level = 9

  procedure(), pointer :: make_base_grids, post_process_grid
  integer :: ilevel
  logical :: level(max_level)


  ! select run level
  level = .false.
  if (run_level == 0) then
     level(1:6) = .true.
  elseif (run_level > 0  .and.  run_level <=max_level) then
     run_level_end = max(run_level,run_level_end)
     do ilevel=run_level,run_level_end
        level(ilevel) = .true.
     enddo
  else
     write (6, *) 'error: run level ', run_level, ' not implemented!'
     stop
  endif
  if (firstP) then
     write (6, *) 'Generating fieldline grid ...'
     write (6, *) '... executing run level ', run_level, ' -> ', run_level_end
  endif


  ! initialize grid configuration
  call setup_grid_configuration()
  select case(topology)
  case(TOPO_SC, TOPO_SC1)
     call setup_topo_sc()
     make_base_grids   => make_base_grids_sc
     post_process_grid => post_process_grid_sc
  case(TOPO_LSN, TOPO_LSN1)
     call setup_topo_lsn()
     make_base_grids   => make_base_grids_lsn
     post_process_grid => post_process_grid_lsn
  case(TOPO_DDN, TOPO_DDN1)
     call setup_topo_ddn()
     make_base_grids   => make_base_grids_ddn
     post_process_grid => post_process_grid_ddn
  case(TOPO_CDN, TOPO_CDN1)
     call setup_topo_cdn()
     make_base_grids   => make_base_grids_cdn
     post_process_grid => post_process_grid_cdn
  case(TOPO_DSFP, TOPO_DSFP1)
     call setup_topo_dsfp()
     make_base_grids   => make_base_grids_dsfp
     post_process_grid => post_process_grid_dsfp
  case(TOPO_STEL, TOPO_STEL1)
     call setup_topo_stel()
     make_base_grids   => make_base_grids_stel
     post_process_grid => post_process_grid_stel
  case default
     write (6, *) 'error: grid topology ', trim(topology), ' not supported!'
     stop
  end select


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
     call trace_nodes(post_process_grid)
  endif


  ! Level 4: generate vacuum domain (used by EIRENE only)
  if (level(4)) then
     call initialize_emc3_grid()
     call vacuum_and_core_domain_for_EIRENE()
     call write_emc3_grid()
     call write_emc3_input_files()
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

  ! Level 7: generate vacuum domain from existing grid (grid3D.dat, input.geo)
  if (level(7)) then
     call load_emc3_grid()
     call vacuum_and_core_domain_for_EIRENE()
     call write_emc3_grid()
     call write_emc3_input_files()
  endif


  ! Level 8: run grid checks
  if (level(8)) then
     call load_emc3_grid()
     call check_emc3_grid()
  endif


  ! Level 9: re-run post processing
  if (level(9)) then
     call load_emc3_grid()
     call post_process_grid()
     call write_emc3_grid()
  endif


end subroutine generate_field_aligend_grid
