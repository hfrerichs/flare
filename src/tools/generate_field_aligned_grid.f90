subroutine generate_field_aligend_grid (run_level)
  use parallel
  use field_aligned_grid
  implicit none

  integer, intent(in) :: run_level

  logical :: level(4)


  if (firstP) then
     write (6, *) 'Generating field aligned grid ...'
     write (6, *) '... executing run level ', run_level
  endif


  call load_usr_conf


  ! select run level
  level = .false.
  if (run_level == 0) then
     level = .true.
  elseif (run_level > 0  .and.  run_level <=4) then
     level(run_level) = .true.
  else
     write (6, *) 'error: run level ', run_level, ' not implemented!'
     stop
  endif


  ! Level 1: generate innermost boundaries (for field line reconstruction)
  if (level(1)) then
     call generate_innermost_boundaries
  endif


  ! Level 2: generate base layout (outer boundary + separatrix for block-structure)
  if (level(2)) then
  endif

!
!      call initialize_grid
!      call setup_domain
!
!      if (istage.eq.-1 .or. istage.eq.0) then
!         call prepare_inner_boundaries
!      endif
!
!      ! step 1: 2D base grid
!      if (istage.eq.-1 .or. istage.eq.1) then
!         call make_base_grids
!      endif
!
!      ! step 2: 3D field aligned grid
!      if (istage.eq.-1  .or. istage.eq.2) then
!         call write_surface_mappings
!         call write_grid_info
!         call trace_grids
!         call close_grid_domain
!         call extend_grid
!         call write_emc3_grid
!         call sample_bfield_on_emc3_grid (PlX, elXm)
!         call make_plate_info
!      endif
!
!      ! step 3: add additional region for neutral particles
!      if (istage.eq.3) then
!         call load_emc3_grid
!         call extend_grid
!         call write_emc3_grid
!      endif
!
!      if (istage.eq.4) then
!         call load_emc3_grid
!         call sample_bfield_on_emc3_grid (PlX, elXm)
!         !call make_plate_info
!      endif
!

end subroutine generate_field_aligend_grid
