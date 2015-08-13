  program plot_grid
  use iso_fortran_env
  use emc3_grid
  use grid
  implicit none

  integer, parameter :: iu = 99

  character(len=72) :: grid_file, plot_file
  character(len=4)  :: s_iz, s_it
  real(real64), dimension(:), allocatable :: R_grid, Z_grid, phi_grid
  type(t_grid) :: G
  real(real64) :: c
  integer      :: iz, ir, ir1, ir2, ip, it, ig, idomain, nr, np


  call load_emc3_grid()


! select output
  write (6, *) 'r_surf_pl_trans_range:'
  do iz=0,nzonet-1
     write (6, *) iz, R_SURF_PL_TRANS_RANGE(:,iz)
  enddo

  write (6, *) 'select zone:'
  read  (5, *) iz

  write (6, *) 'select slice:'
  read  (5, *) it

  write (s_iz, '(i4)') iz
  write (s_it,  '(i4)') it

  write (6, *) 'select domain (0: plasma, 1: plasma+neutral gas)'
  read  (5, *) idomain
  if (idomain == 0) then
     ir1 = R_SURF_PL_TRANS_RANGE(1,iz)
     ir2 = R_SURF_PL_TRANS_RANGE(2,iz)
  elseif (idomain == 1) then
     ir1 = 0
     ir2 = SRF_RADI(iz)-1
  else
     stop
  endif
  write (6, *) 'radial domain: ', ir1, ' -> ', ir2
  nr  = ir2 - ir1 + 1
  np  = SRF_POLO(iz)


  c = 3.14159265358979323846264338328d0 / 180.d0
  write (6, *) 'phi = ', PHI_PLANE(it+PHI_PL_OS(iz)) / c


  call G%new(2, MESH_2D, 3, nr, np, fixed_coord_value=PHI_PLANE(it+PHI_PL_OS(iz)))
  do ir=ir1,ir2
  do ip=0,SRF_POLO(iz)-1
     ig = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
     G%mesh(ir-ir1,ip,1) = RG(ig)
     G%mesh(ir-ir1,ip,2) = ZG(ig)
  enddo
  enddo

  !plot_file = 'grid_'//trim(adjustl(s_iz))//'_'//trim(adjustl(s_it))//'.dat'
  !call G%store(filename=plot_file)
  plot_file = 'grid_'//trim(adjustl(s_iz))//'_'//trim(adjustl(s_it))//'.plt'
  call G%plot_mesh(plot_file)

  end program
