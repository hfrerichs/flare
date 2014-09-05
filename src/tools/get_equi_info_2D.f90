!===============================================================================
! Write equilibrium information
!===============================================================================
subroutine get_equi_info_2D
  use equilibrium
  use boundary
  use parallel
  implicit none

  integer, parameter :: iu = 42


  character*120 :: fout, sboundary
  character*8   :: cid
  real*8  :: Rbox(2), Zbox(2), Px(2), r(3), Psi
  integer :: i, j, nR, nZ


  if (firstP) then
     write (6, *) 'Write equilibrium information'
     write (6, *)
  else
     return
  endif


  call get_domain (Rbox, Zbox)
  ! get X-point
  Px = find_X()
  write (6, 9000) Px


  ! write grid with sample locations
  nR = 128
  nZ = 128
  fout = 'equi_info.grid'
  open  (iu, file=fout)
  write (iu, 1000)
  write (iu, 1001) nR
  write (iu, 1002) nZ
  write (iu, 1003) 0.d0
  do i=0,nR-1
     r(1) = Rbox(1) + (Rbox(2)-Rbox(1)) * i / (nR-1)
     write (iu, 1004) r(1)
  enddo
  do j=0,nZ-1
     r(2) = Zbox(1) + (Zbox(2)-Zbox(1)) * j / (nZ-1)
     write (iu, 1004) r(2)
  enddo
  close (iu)


  ! write normalized poloidal magnetic flux
  fout = 'equi_info.data'
  r(3) = 0.d0
  open  (iu, file=fout)
  do i=0,nR-1
  do j=0,nZ-1
     r(1) = Rbox(1) + (Rbox(2)-Rbox(1)) * i / (nR-1)
     r(2) = Zbox(1) + (Zbox(2)-Zbox(1)) * j / (nZ-1)
     Psi = get_Psi(r)
     Psi = (Psi - Psi_axis) / (Psi_sepx - Psi_axis)
     write (iu, 1004) Psi
  enddo
  enddo
  close (iu)


  ! write boundary data
  sboundary = 'boundary=['
  do i=1,n_axi
     write (cid, '(i8)') i
     fout      = 'boundary_'//trim(adjustl(cid))//'.txt'
     sboundary = trim(sboundary)//"'"//trim(fout)//"'"
     if (i < n_axi) sboundary = trim(sboundary)//', '

     open  (iu, file=fout)
     do j=0,S_axi(i)%n_seg
        write (iu, *) S_axi(i)%x_data(j,1:2)
     enddo
     close (iu)
  enddo
  sboundary = trim(sboundary)//'], $'


  ! write plot script
  open  (iu, file='plot_equi_info.sh')
  write (iu, 2000)
  write (iu, 2001)
  write (iu, 2002)
  write (iu, 2003)
  write (iu, 2004)
  write (iu, 2005)
  write (iu, 2006)
  write (iu, 2007)
  write (iu, 2008) sboundary
  write (iu, 2009)
  write (iu, 2010)
  close (iu)


 9000 format (3x,'- found X-point at: ',2f10.3)
 1000 format ('# grid_id = 3       (regular RZ grid)')
 1001 format ('# R resolution:      n_R     =  ',i8)
 1002 format ('# Z resolution:      n_Z     =  ',i8)
 1003 format ('# phi position:      phi     =  ',f7.2)
 1004 format (e18.10)

 2000 format ('#!/bin/bash')
 2001 format ('if [ ! -e plot_data.pro ]; then')
 2002 format ('echo "Please make (symbolic) link to plot_data.pro!"')
 2003 format ('echo "This file can be found in the FLARE directory under templates/plot"')
 2004 format ('exit')
 2005 format ('fi')
 2006 format ('idl << EOF')
 2007 format ("plot_data, 'equi_info.grid', 'equi_info.data', 0, zbound=[0,2], clevels=[0.9, 1.0, 1.1], $")
 2008 format (a120)
 2009 format ("ps_plot='equi_info.eps'")
 2010 format ('EOF')
end subroutine get_equi_info_2D
