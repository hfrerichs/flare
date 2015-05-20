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
        call setup_vacuum_domain(iz, nr_EIRENE_vac, 1, Zone(iz)%d_N0, 1)
     endif
     if (Zone(iz)%isfr(2) == SF_VACUUM) then
        call setup_vacuum_domain(iz, nr_EIRENE_vac, 2, Zone(iz)%d_N0, 1)
     endif
  enddo

end subroutine vacuum_and_core_domain_for_EIRENE
!===============================================================================



!===============================================================================
subroutine setup_core_domain(iz, nr_core)
  use iso_fortran_env
  use emc3_grid
  use fieldline_grid
  use equilibrium
  implicit none

  integer, intent(in) :: iz, nr_core

  real(real64) :: phi, r3(3), x(2), dx(2)
  integer :: i, j, k, ig


  write (6, 1000) iz, nr_core
 1000 format(8x,'zone ',i0,': ',i0,' core cell(s)')

  do k=0,SRF_TORO(iz)-1
     phi = PHI_PLANE(k+PHI_PL_OS(iz))
     r3  = get_magnetic_axis(phi)

     do j=0,SRF_POLO(iz)-1
        ig = i + (j + k*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
        dx(1) = RG(ig+1) - r3(1)
        dx(2) = ZG(ig+1) - r3(2)

        x     = r3(1:2) + 0.5d0 * dx
        RG(ig) = x(1)
        ZG(ig) = x(2)
     enddo
  enddo

end subroutine setup_core_domain
!===============================================================================



!===============================================================================
subroutine setup_vacuum_domain(iz, nr_vac, boundary, dl, Method)
  use iso_fortran_env
  implicit none

  integer, intent(in) :: iz, nr_vac, boundary, Method
  real(real64), intent(in) :: dl

  integer, parameter :: &
     UPSCALE = 1


  write (6, 1000) iz, nr_vac
 1000 format(8x,'zone ',i0,': ',i0,' vacuum cell(s)')

  select case (Method)
  case (UPSCALE)
      call vacuum_domain_by_upscale(iz, nr_vac, boundary, dl)
  end select


end subroutine setup_vacuum_domain
!===============================================================================


!===============================================================================
subroutine vacuum_domain_by_upscale(iz, nrvac, boundary, dl)
  use iso_fortran_env
  use emc3_grid
  use curve2D
  use run_control, only: Debug
  use string
  implicit none

  integer, intent(in)      :: iz, nrvac, boundary
  real(real64), intent(in) :: dl

  real(real64), dimension(:), allocatable :: xi
  type(t_curve) :: C
  real(real64)  :: DR, DZ, w, x(2)
  integer :: it, ip, ir0, ir, ig, idir, irend


  ! initialize internal parameters
  select case(boundary)
  ! lower boundary
  case(1)
    ir0  = nrvac
    idir = -1
    irend = 0

  ! upper boundary
  case(2)
     ir0 = SRF_RADI(iz)-1 - nrvac
    idir = 1
    irend = SRF_RADI(iz)-1

  case default
     write (6, *) 'error: invalid argument boundary = ', boundary
  end select
  allocate (xi(0:SRF_POLO(iz)-1))


  ! loop over all toroidal slices
  do it=0,SRF_TORO(iz)-1
     call C%new(ZON_POLO(iz))
     !TODO: !C%closed = .true.
     ! poloidal loop (setup nodes for curve blowup)
     do ip=0,SRF_POLO(iz)-1
        ig = ir0 + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        C%x(ip,1) = RG(ig)
        C%x(ip,2) = ZG(ig)
     enddo
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
     call C%setup_length_sampling()

     ! poloidal loop (set new grid nodes)
     do ip=0,SRF_POLO(iz)-1
     do ir=ir0+idir,irend,idir
        ig = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)

        call C%sample_at(xi(ip), x)
        RG(ig) = x(1)
        ZG(ig) = x(2)
     enddo
     enddo

     ! cleanup
     call C%destroy()
  enddo

end subroutine vacuum_domain_by_upscale
!===============================================================================
