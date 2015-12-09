  program plot_grid
  use iso_fortran_env
  use emc3_grid
  implicit none

  integer      :: iz, iz1, iz2, it, it0, it1, it2, dit, itmax, idomain


  call load_emc3_grid()


  ! 1. select output zone(s)
  write (6, 1000) NZONET-1
  read  (5, *) iz
  if (iz < -1  .or.  iz > NZONET-1) then
     write (6, *) 'error: invalid zone number!'
     stop
  endif
  iz1 = iz
  iz2 = iz
  if (iz == -1) then
     iz1 = 0
     iz2 = NZONET-1
  endif
 1000 format(3x,'select zone  (0 -> ',i0,', -1: all zones):')


  ! 2. select output slice(s)
  if (iz == -1) then
     itmax = minval(SRF_TORO)
  else
     itmax = SRF_TORO(iz)
  endif
  write (6, 1001) itmax-1
  read  (5, *) it0
  if (it0 < -2  .or.  it0 > itmax) then
     write (6, *) 'error: invalid slice id!'
     stop
  endif
 1001 format(3x,'select slice (0 -> ',i0,', -1: all slices, -2: first and last slice):')


  ! 3. select output mode
  write (6, *) 'select domain (0: plasma, 1: plasma+neutral gas)'
  read  (5, *) idomain
  if (idomain < 0  .or.  idomain > 1) then
     write (6, *) 'error: invalid domain id!'
     stop
  endif


  ! 4. plot loop
  do iz=iz1,iz2
     it1 = it0
     it2 = it0
     dit = 1
     if (it0 < 0) then
        it1 = 0
        it2 = SRF_TORO(iz)-1
     endif
     if (it0 == -2) dit = it2

     do it=it1,it2,dit
        call plot_slice(iz, it, idomain)
     enddo
  enddo

  contains
  !---------------------------------------------------------------------
  subroutine plot_slice(iz, it, idomain)
  use grid
  integer, intent(in) :: iz, it, idomain

  character(len=72) :: grid_file, plot_file
  character(len=4)  :: s_iz, s_it

  type(t_grid) :: G
  real(real64) :: c
  integer      :: ir, ir1, ir2, ip, ig, nr, np


  write (s_iz, '(i4)') iz
  write (s_it,  '(i4)') it


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

  end subroutine plot_slice
  !---------------------------------------------------------------------
  end program
