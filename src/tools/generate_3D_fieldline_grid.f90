subroutine generate_3D_fieldline_grid (run_level, run_level_end)
  use parallel
  use fieldline_grid
  use emc3_grid
  use modtopo_sc
  use modtopo_lsn
  use modtopo_ddn
  use modtopo_cdn
  use modtopo_dsfp
  use modtopo_stel
  use base_mesh, only: make_base_mesh_generic
  implicit none

  integer, intent(inout) :: run_level, run_level_end

  integer, parameter     :: max_level = 10

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
     if (mesh_generator == LEGACY) then
        call make_base_grids()
     elseif (mesh_generator == ORTHOGONAL) then
        call make_base_mesh_generic()
     else
        write (6, *) 'error: invalid mesh generator ', mesh_generator
        stop
     endif
  endif


  ! Level 3: generate 3D grid from field line tracing
  if (level(3)) then
     call trace_nodes(post_process_grid)
  endif


  ! Level 4: generate vacuum domain (used by EIRENE only)
  if (level(4)) then
     call initialize_emc3_grid()
     call setup_core_domain()
     call vacuum_domain_for_EIRENE()
     call check_emc3_grid()
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
     call setup_core_domain()
     call vacuum_domain_for_EIRENE()
     call check_emc3_grid()
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
     call check_emc3_grid()
     call write_emc3_grid()
  endif


  ! Level 10: define divertor reservoirs
  if (level(10)) then
     call load_emc3_grid()
     call define_divertor_reservoirs()
  endif

end subroutine generate_3D_fieldline_grid





subroutine define_divertor_reservoirs()
  use fieldline_grid, only: Zone, npR, npL, topology, TOPO_LSN, TOPO_LSN1
  use emc3_grid
  use curve2D
  use run_control, only: Debug
  implicit none

  type(t_curve) :: CL, CR
  integer, dimension(:), allocatable :: res
  real*8  :: x1(2), x2(2), x(2)
  integer :: ir, ip, it, iz, im, ig, nr


  write (6, *) "setup divertor reservoirs..."
  if (topology /= TOPO_LSN  .and.  topology /= TOPO_LSN1) then
     write (6, 9001);   stop
  endif
 9001 format("error: divertor reservoirs implemented only for single null geometry!")


  allocate (res(0:MESH_P_OS(NZONET)-1), source=0)
  do iz=0,NZONET-1
     if (mod(iz,3) == 0) cycle


  it = Zone(iz)%it_base
  ip = npR(1)
  nr = Zone(iz)%nr
  call CR%new(nr)
  do ir=0,nr
     ig = ir +(ip+it*SRF_POLO(iz))*SRF_RADI(iz)+GRID_P_OS(iz)
     CR%x(ir,1) = RG(ig)
     CR%x(ir,2) = ZG(ig)
  enddo

  ip = ZON_POLO(iz)-npL(1)
  call CL%new(nr)
  do ir=0,nr
     ig = ir +(ip+it*SRF_POLO(iz))*SRF_RADI(iz)+GRID_P_OS(iz)
     CL%x(ir,1) = RG(ig)
     CL%x(ir,2) = ZG(ig)
  enddo


  do ir=0,nr-1
  do it=0,ZON_TORO(iz)-1
     ! right divertor leg
     call RZ_REAL_COORDINATES(iz, ir, 0, 0.d0,0.d0, it+0.5d0, x1(1), x1(2))
     do ip=0,ZON_POLO(iz)-2
        im = ir + (ip+it*ZON_POLO(iz))*ZON_RADI(iz) + MESH_P_OS(iz)
        res(im) = 1
        call RZ_REAL_COORDINATES(iz, ir, ip+1, 0.d0,0.d0, it+0.5d0, x2(1), x2(2))
        if (intersect_curve(x1, x2, CR, x)) then
           if (Debug) write (97, *) x
           exit
        endif

        x1 = x2
     enddo
     !write (6, *) ir, it, ip

     ! left divertor leg
     call RZ_REAL_COORDINATES(iz, ir, ZON_POLO(iz)-1, 0.d0,0.d0, it+0.5d0, x1(1), x1(2))
     do ip=ZON_POLO(iz)-1,1,-1
        im = ir + (ip+it*ZON_POLO(iz))*ZON_RADI(iz) + MESH_P_OS(iz)
        res(im) = 2
        call RZ_REAL_COORDINATES(iz, ir, ip-1, 0.d0,0.d0, it+0.5d0, x2(1), x2(2))
        if (intersect_curve(x1, x2, CL, x)) then
           if (Debug) write (97, *) x
           exit
        endif

        x1 = x2
     enddo
     !write (6, *) ir, it, ip
  enddo
  enddo
  enddo


  open  (99, file="DIVERTOR_RESERVOIRS")
  write (99, *) res
  close (99)

end subroutine define_divertor_reservoirs
