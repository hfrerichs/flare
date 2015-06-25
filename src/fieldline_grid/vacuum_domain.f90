!===============================================================================
! Set up additional domain for EIRENE (vacuum region in far SOL and PFR, and
! core region)
!===============================================================================
subroutine vacuum_and_core_domain_for_EIRENE
  use emc3_grid
  use fieldline_grid
  implicit none

  integer :: iz


  write (6, 1000)
 1000 format(3x,' - Set up additional domain for EIRENE')
  do iz=0,NZONET-1
     if (Zone(iz)%isfr(1) == SF_CORE) then
        call setup_core_domain(iz, nr_EIRENE_core)
     endif

     if (Zone(iz)%isfr(1) == SF_VACUUM) then
        call setup_vacuum_domain(iz, nr_EIRENE_vac, 1)
     endif
     if (Zone(iz)%isfr(2) == SF_VACUUM) then
        call setup_vacuum_domain(iz, nr_EIRENE_vac, 2)
     endif
  enddo

end subroutine vacuum_and_core_domain_for_EIRENE
!===============================================================================



!===============================================================================
subroutine setup_core_domain(iz, nr_core)
  use iso_fortran_env
  use emc3_grid
  use fieldline_grid
  implicit none

  integer, intent(in) :: iz, nr_core


  integer, parameter :: &
     MAGNETIC_AXIS    = 1, &
     GEOMETRIC_CENTER = 2


  real(real64) :: phi, r0(2), x(2), dx(2)
  integer :: i, j, k, ig, method = GEOMETRIC_CENTER


  write (6, 1000) iz, nr_core
 1000 format(8x,'zone ',i0,': ',i0,' core cell(s)')

  do k=0,SRF_TORO(iz)-1
     phi = PHI_PLANE(k+PHI_PL_OS(iz))
     r0  = get_r0(phi)

     do j=0,SRF_POLO(iz)-1
        ig = (j + k*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
        dx(1)  = RG(ig+nr_core) - r0(1)
        dx(2)  = ZG(ig+nr_core) - r0(2)

        do i=0,nr_core-1
           x        = r0 + 0.5d0 * dx * (1.d0 + 1.d0 * i / nr_core)
           RG(ig+i) = x(1)
           ZG(ig+i) = x(2)
        enddo
     enddo
  enddo

  return
  contains
  !---------------------------------------------------------------------
  function get_r0(phi)
  use equilibrium, only: get_magnetic_axis
  real(real64), intent(in) :: phi
  real(real64)             :: get_r0(2)

  real(real64) :: r3(3), w, w_tot
  integer      :: j, ig2


  get_r0 = 0.d0

  select case(method)
  case(MAGNETIC_AXIS)
     r3     = get_magnetic_axis(phi)
     get_r0 = r3(1:2)

  case(GEOMETRIC_CENTER)
     get_r0 = 0.d0
     w_tot  = 0.d0
     do j=0,ZON_POLO(iz)-1
        ig  = nr_core + (j + k*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
        ig2 = ig + SRF_RADI(iz)
        w   = sqrt((RG(ig)-RG(ig2))**2 + (ZG(ig)-ZG(ig2))**2)

        get_r0(1) = get_r0(1) + w*RG(ig)
        get_r0(2) = get_r0(2) + w*ZG(ig)
        w_tot     = w_tot     + w
     enddo
     get_r0 = get_r0 / w_tot
  end select

  end function get_r0
  !---------------------------------------------------------------------
end subroutine setup_core_domain
!===============================================================================



!===============================================================================
subroutine setup_vacuum_domain(iz, nr_vac, boundary)
  use iso_fortran_env
  use fieldline_grid
  implicit none

  character(len=*), parameter :: s_boundary(2) = (/ 'lower', 'upper' /)
  integer, intent(in) :: iz, nr_vac, boundary

  integer, parameter :: &
     UPSCALE = 1, &
     MANUAL  = -1

  real(real64) :: dl
  integer      :: Method, ir0, idir, ir2


  Method = UPSCALE
  dl     = Zone(iz)%d_N0
  if (Zone(iz)%N0_file .ne. '') then
     Method = MANUAL
  endif


  ! set surface indices and increment
  ! ir0:  surface index for EMC3 boundary
  ! idir: index direction for EIRENE-only domain
  ! ir2:  surface index for EIRENE boundary
  select case(boundary)
  ! lower boundary
  case(1)
     ir0  = nr_vac
     idir = -1
     ir2  = 0

  ! upper boundary
  case(2)
     ir0  = Zone(iz)%nr - nr_vac
     idir = 1
     ir2  = Zone(iz)%nr

  case default
     write (6, *) 'error: invalid argument boundary = ', boundary
  end select



  write (6, 1000) iz, nr_vac, s_boundary(boundary), dl
 1000 format(8x,'zone ',i0,': ',i0,' vacuum cell(s) at ',a,' boundary, D = ',f8.3)

  select case (Method)
  case (UPSCALE)
      call vacuum_domain_by_upscale(iz, ir0, idir, ir2, dl)
  case (MANUAL)
      call vacuum_domain_manual(iz, ir0, idir, ir2, Zone(iz)%N0_file)
  end select


end subroutine setup_vacuum_domain
!===============================================================================


!===============================================================================
subroutine vacuum_domain_by_upscale(iz, ir0, idir, ir2, dl)
  use iso_fortran_env
  use emc3_grid
  use curve2D
  use run_control, only: Debug
  use string
  implicit none

  integer, intent(in)      :: iz, ir0, idir, ir2
  real(real64), intent(in) :: dl

  logical :: resample = .true.

  real(real64), dimension(:), allocatable :: xi
  real(real64), dimension(:,:), allocatable :: en
  type(t_curve) :: C
  real(real64)  :: DR, DZ, w, x(2), rho
  integer       :: it, ip, ir, ir1, ig, ig0


  allocate (xi(0:SRF_POLO(iz)-1))
  ir1 = ir0 + idir


  ! loop over all toroidal slices
  allocate (en(0:SRF_POLO(iz)-1,2))
  do it=0,SRF_TORO(iz)-1
     call C%new(ZON_POLO(iz))
     ! poloidal loop (setup nodes for curve blowup)
     do ip=0,SRF_POLO(iz)-1
        ig = ir0 + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        C%x(ip,1) = RG(ig)
        C%x(ip,2) = ZG(ig)
        en(ip,1)  = RG(ig) - RG(ig-idir)
        en(ip,2)  = ZG(ig) - ZG(ig-idir)
     enddo
     call C%closed_check()
     if (Debug) then
        call C%plot(filename='debug/VacuumBoundary_'//trim(str(iz))//'_'//trim(str(it))//'.raw')
     endif

     ! poloidal loop (setup segment weights)
     xi(0) = 0.d0
     do ip=1,SRF_POLO(iz)-1
        ig = ir0 + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        DR = RG(ig) - RG(ig - SRF_RADI(iz))
        DZ = ZG(ig) - ZG(ig - SRF_RADI(iz))
        xi(ip) = xi(ip-1) + sqrt(DR**2 + DZ**2)
     enddo
     xi = xi / xi(SRF_POLO(iz)-1)

     call C%left_hand_shift(-idir*dl)
     if (Debug) then
        call C%plot(filename='debug/VacuumBoundary_'//trim(str(iz))//'_'//trim(str(it))//'.plt')
     endif
     if (resample.eqv..false.  .and.  ZON_POLO(iz).ne.C%n_seg) then
        write (6, *) 'error: nodes were dropped in subroutine left_hand_shift!'
        write (6, *) 'iz, it = ', iz, it
        stop
     endif
     call C%setup_length_sampling()

     ! poloidal loop (set new grid nodes)
     do ip=0,SRF_POLO(iz)-1
        ig0 = ir0 + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        if (resample) then
           call C%sample_at(xi(ip), x)
        else
           x = C%x(ip,1:2)
        endif

        do ir=ir1,ir2,idir
           rho = 1.d0 * (ir-ir0) / (ir2-ir1+idir)
           ig  = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)

           RG(ig) = RG(ig0) + rho * (x(1) - RG(ig0))
           ZG(ig) = ZG(ig0) + rho * (x(2) - ZG(ig0))
        enddo
     enddo

     ! cleanup
     call C%destroy()
  enddo
  deallocate (en)

end subroutine vacuum_domain_by_upscale
!===============================================================================



!===============================================================================
subroutine vacuum_domain_manual(iz, ir0, idir, ir2, boundary_file)
  use iso_fortran_env
  use emc3_grid
  use curve2D
  use run_control, only: Debug
  use string
  implicit none

  integer, intent(in)      :: iz, ir0, idir, ir2
  character(len=*), intent(in) :: boundary_file

  type(t_curve) :: C
  real(real64), dimension(:), allocatable :: eta
  real(real64)  :: DR, DZ, x(2), rho
  integer       :: ig, ig0, it, ip, ir, ir1


  call C%load(boundary_file)
  call C%setup_length_sampling()
  allocate (eta(0:SRF_POLO(iz)-1))
  ir1 = ir0 + idir


  do it=0,SRF_TORO(iz)-1
     ! poloidal loop (setup segment weights)
     eta(0) = 0.d0
     do ip=1,SRF_POLO(iz)-1
        ig = ir0 + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        DR = RG(ig) - RG(ig - SRF_RADI(iz))
        DZ = ZG(ig) - ZG(ig - SRF_RADI(iz))
        eta(ip) = eta(ip-1) + sqrt(DR**2 + DZ**2)
     enddo
     eta = eta / eta(SRF_POLO(iz)-1)


     ! poloidal loop (set new grid nodes)
     do ip=0,SRF_POLO(iz)-1
        ig0 = ir0 + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        call C%sample_at(eta(ip), x)

        do ir=ir1,ir2,idir
           rho = 1.d0 * (ir-ir0) / (ir2-ir1+idir)
           ig  = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)

           RG(ig) = RG(ig0) + rho * (x(1) - RG(ig0))
           ZG(ig) = ZG(ig0) + rho * (x(2) - ZG(ig0))
        enddo
     enddo
  enddo

  deallocate (eta)
end subroutine vacuum_domain_manual
