subroutine check_emc3_grid
  use iso_fortran_env
  use emc3_grid
  use flux_tube
  use math
  use grid
  use dataset
  implicit none

  integer, parameter :: iu = 72
  type(t_flux_tube) :: FT, FT2
  type(t_cross_section) :: cs0
  type(t_grid) :: G
  type(t_dataset) :: D
  type(t_cross_section) :: cs
  integer :: iz
  integer :: ir, ip, it, ig(4), i


  allocate (BFSTREN(0:GRID_P_OS(NZONET)-1))
  open  (iu, file='bfield.dat')
  read  (iu, *) BFSTREN
  close (iu)



!  select case(check_type)
!  case(flux_tube)
      FT = export_flux_tube(5, 18, 338)
      call FT%plot('flux_tube.plt')
      D = FT%get_flux_along_tube()
      call D%plot(filename='post_proc/flux_along_tube.dat')

      cs = FT%get_cross_section(ZON_TORO(5)/2)
      cs%a1 = cs%a1/10.d0
      cs%a2 = cs%a2/10.d0
      cs%a3 = cs%a3/10.d0
      call FT2%generate(cs, 8, 72)
      call FT2%plot('flux_tube_small.plt')
      D = FT2%get_flux_along_tube()
      call D%plot(filename='post_proc/flux_along_small_tube.dat')

      D = FT%sample_pitch(8, 100, 100)
      call D%plot(filename='post_proc/pitch.dat')
      return
!  case(all_cells)
!  end select







  do iz=0,NZONET-1
     call flux_conservation(iz)
  enddo

  iz = 5
  call FT%new(ZON_TORO(iz))
  FT%phi = PHI_PLANE(PHI_PL_OS(iz):PHI_PL_OS(iz+1)-1) / 180.d0 * pi
  ir = 18
  ip = 338
  do it=0,ZON_TORO(iz)
     ig(1) = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
     ig(2) = ig(1) + 1
     ig(4) = ig(1) + SRF_RADI(iz)
     ig(3) = ig(4) + 1

     do i=1,4
        FT%F(i)%x(it,1) = RG(ig(i))
        FT%F(i)%x(it,2) = ZG(ig(i))
     enddo
  enddo

  call FT%plot('flux_tube.plt')
  cs0 = FT%get_cross_section(ZON_TORO(iz)/2)
  call cs0%setup_shape()
  G = cs0%generate_mesh(10,10)
  call G%store('mesh0.dat')


  contains
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  ! this subroutine is taken from the file check.f (EMC3)
  !---------------------------------------------------------------------
  subroutine flux_conservation(iz)
  use grid
  use dataset
  use math
  integer, intent(in) :: iz

  real(real64), dimension(:,:,:), allocatable :: nl_cell
  real(real64), dimension(:,:),   allocatable :: div
  real(real64), dimension(:),     allocatable :: flux, flux_save
  real(real64) :: divmax, divave, f_sum, non_linear_ave
  real(real64) :: abcd1, abcd2, a1, a2, b1, b2, c1, c2, d1, d2, d, &
                  non_linear, area1, area2, R1, R2, Z1, Z2, Phi1, Phi2, &
                  BF1, BF2, DFL, flux_ave, pitch
  integer      :: i, j, k, i1, i2, i3, i4
  type(t_grid) :: G
  type(t_dataset) :: data_out
  character(len=120) :: filename
  integer :: imax, jmax


  divmax = 0.d0
  divave = 0.d0
  f_sum  = 0.d0
  non_linear_ave = 0.d0
  allocate (flux(ZON_TORO(iz)))
  allocate (flux_save(ZON_TORO(iz)))
  allocate (div(P_SURF_PL_TRANS_RANGE(1,iz):P_SURF_PL_TRANS_RANGE(2,iz)-1, &
                R_SURF_PL_TRANS_RANGE(1,iz):R_SURF_PL_TRANS_RANGE(2,iz)-1))
  allocate (nl_cell(0:ZON_TORO(iz), 0:ZON_POLO(iz)-1, 0:ZON_RADI(iz)-1))

  do i=0,ZON_RADI(iz)-1
  do j=0,ZON_POLO(iz)-1
     do k=0,ZON_TORO(iz)
        i1 = i + (j + k*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        i2 = i1 + SRF_RADI(iz)
        i3 = i2 + 1
        i4 = i1 + 1

        abcd1 = (RG(I3)-RG(I2))*(ZG(I2)-ZG(I1)) &
              - (RG(I2)-RG(I1))*(ZG(I3)-ZG(I2))
        abcd2 = (RG(I4)-RG(I1))*(ZG(I3)-ZG(I4)) &
              - (RG(I3)-RG(I4))*(ZG(I4)-ZG(I1))

        a1 = 0.25d0*(RG(I3)+RG(I4)-RG(I1)-RG(I2))
        b1 = 0.25d0*(RG(I3)+RG(I2)-RG(I1)-RG(I4))
        c1 = 0.25d0*(RG(I1)+RG(I3)-RG(I2)-RG(I4))

        a2 = 0.25d0*(ZG(I3)+ZG(I4)-ZG(I1)-ZG(I2))
        b2 = 0.25d0*(ZG(I3)+ZG(I2)-ZG(I1)-ZG(I4))
        c2 = 0.25d0*(ZG(I1)+ZG(I3)-ZG(I2)-ZG(I4))

        d  = abs(a1*b2 - a2*b1)

        d1 = -(c1*b2-c2*b1)
        d2 = -(a1*c2-a2*c1)
        non_linear = abs(d1)+abs(d2)
        nl_cell(k,j,i) = non_linear/d
        non_linear_ave = non_linear_ave + non_linear/d

        ! check flux conservation in plasma transport range
        if (i>=R_SURF_PL_TRANS_RANGE(1,iz) .and. &
            i< R_SURF_PL_TRANS_RANGE(2,iz) .and. &
            j>=P_SURF_PL_TRANS_RANGE(1,iz) .and. &
            j< P_SURF_PL_TRANS_RANGE(2,iz)) then 

           area2 = 0.5d0*abs(abcd1 + abcd2)
           R2    = 0.25d0*(RG(i1)+RG(i2)+RG(i3)+RG(i4))
           Z2    = 0.25d0*(ZG(i1)+ZG(i2)+ZG(i3)+ZG(i4))
           Phi2  = PHI_PLANE(PHI_PL_OS(iz)+k)
           Phi2  = Phi2 / 180.d0 * pi	! TODO: consistent units for toroidal angle
           BF2   = 0.25d0*(BFSTREN(i1) + BFSTREN(i2) &
                         + BFSTREN(i3) + BFSTREN(i4))

           if (k/=0) then
              DFL    = 0.5d0*(R1+R2)*abs(Phi2-Phi1)
              pitch  = DFL/sqrt(DFL**2+(R2-R1)**2+(Z2-Z1)**2)
              flux(k)= (area1+area2)*pitch*(BF2+BF1)*0.25d0
           endif
           area1 = area2; R1 = R2; Z1 = Z2; Phi1 = Phi2; BF1 = BF2
        endif
     enddo

     if (i>=R_SURF_PL_TRANS_RANGE(1,iz) .and. &
         i< R_SURF_PL_TRANS_RANGE(2,iz) .and. &
         j>=P_SURF_PL_TRANS_RANGE(1,iz) .and. &
         j< P_SURF_PL_TRANS_RANGE(2,iz)) then

        flux_ave = sum(flux)/float(ZON_TORO(iz))
        div(j,i) = (maxval(flux)-minval(flux))/flux_ave
        if (div(j,i) > divmax) then
           imax = i
           jmax = j
           flux_save = flux
        endif
        divmax   = max(divmax, div(j,i))
        divave   = divave + maxval(flux)-minval(flux)   
        f_sum    = f_sum + flux_ave
     endif
  enddo
  enddo


  non_linear_ave = non_linear_ave/float(MESH_P_OS(NZONET))
  write (6,*) ' 1.1 mesh shape '
  write (6,'(a18,f8.3,a2)') '   Non-linearity:', non_linear_ave*100.d0,' %'

  divave = divave/(f_sum+1.d-10)*100.d0
  write (6,*) ' 1.2 Check accuracy in flux conservation'
  write (6,'(a18,f8.3,a2)') '   Flux conserved:',100.d0-divave,' %'
  write (6,'(a18,f8.3,a2)') '   max. diviation:',divmax*100.d0,' %'
  write (6,'(a18,i0,2x,i0)')'       at (i,j) = ',imax,jmax
  write (6, *) 'test ', div(jmax,imax)
  write (filename, 2000) iz
  2000 format('post_proc/max_flux_deviation_',i0,'.dat')
  open  (99, file=filename)
  write (99, *) '# ', iz, jmax, imax
  do k=1,ZON_TORO(iz)
     write (99, *) flux_save(k)
  enddo
  close (99)


  ! TODO: check up/down symmetric surfaces

  ! write data
!  write (filename, 1000) iz
!  1000 format('post_proc/grid_',i0,'.dat')
!  call G%new(LOCAL, STRUCTURED, 3, ZON_POLO(iz), ZON_RADI(iz), 1)

  write (filename, 1000) iz
  1000 format('post_proc/div_',i0,'.dat')
  open  (iu, file=filename)
  do j=P_SURF_PL_TRANS_RANGE(1,iz),P_SURF_PL_TRANS_RANGE(2,iz)-1
     do i=R_SURF_PL_TRANS_RANGE(1,iz),R_SURF_PL_TRANS_RANGE(2,iz)-1
!  do j=0,ZON_POLO(iz)-1
!     do i=0,ZON_RADI(iz)-1
        write (iu, *) j, i, div(j,i)*100.d0, nl_cell(:,j,i)
     enddo
     write (iu, *)
  enddo
  close (iu)


  ! cleanup
  deallocate (flux, div, nl_cell)
  deallocate (flux_save)

  end subroutine flux_conservation
  !---------------------------------------------------------------------

end subroutine check_emc3_grid
