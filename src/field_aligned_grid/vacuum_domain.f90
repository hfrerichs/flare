!===============================================================================
! Generate vacuum domain (used by EIRENE only)
!===============================================================================
subroutine vacuum_domain_for_EIRENE
  use emc3_grid
  use field_aligned_grid, only: TD, EDGE_BOUNDARY_SF
  implicit none

  integer :: iz


  do iz=0,NZONET-1
     if (TD(iz+1)%irsfb == EDGE_BOUNDARY_SF) then
        call setup_vacuum_domain(iz, 1, 1)
     endif
  enddo

end subroutine vacuum_domain_for_EIRENE
!===============================================================================



!===============================================================================
subroutine setup_vacuum_domain(iz, nrvac, Method)
  use iso_fortran_env
  implicit none

  integer, intent(in) :: iz, nrvac, Method

  integer, parameter :: &
     UPSCALE = 1


  write (6, *) 'setting up vacuum domain for zone ', iz
  select case (Method)
  case (UPSCALE)
      call vacuum_domain_by_upscale(iz, nrvac, -30.d0)
  end select


end subroutine setup_vacuum_domain
!===============================================================================


!===============================================================================
subroutine vacuum_domain_by_upscale(iz, nrvac, dl)
  use iso_fortran_env
  use emc3_grid
  use curve2D
  implicit none

  integer, intent(in)      :: iz, nrvac
  real(real64), intent(in) :: dl

  real(real64), dimension(:), allocatable :: xi
  type(t_curve) :: C
  real(real64)  :: DR, DZ, w, x(2)
  integer :: it, ip, ir0, ir, ig


  ! initialize internal parameters
  ir0 = SRF_RADI(iz)-1 - nrvac
  allocate (xi(0:SRF_POLO(iz)-1))


  ! loop over all toroidal slices
  do it=0,SRF_TORO(iz)-1
     call C%new(ZON_POLO(iz))
     C%closed = .true.
     ! poloidal loop (setup nodes for curve blowup)
     do ip=0,SRF_POLO(iz)-1
        ig = ir0 + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        C%x(ip,1) = RG(ig)
        C%x(ip,2) = ZG(ig)
     enddo

     ! poloidal loop (setup segment weights)
     xi(0) = 0.d0
     do ip=1,SRF_POLO(iz)-1
        ig = ir0 + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        DR = RG(ig) - RG(ig - SRF_RADI(iz))
        DZ = ZG(ig) - ZG(ig - SRF_RADI(iz))
        xi(ip) = xi(ip-1) + sqrt(DR**2 + DZ**2)
     enddo
     xi = xi / xi(SRF_POLO(iz)-1)

     call C%left_hand_shift(dl)
     call C%setup_length_sampling()

     ! poloidal loop (set new grid nodes)
     do ip=0,SRF_POLO(iz)-1
     do ir=ir0+1,SRF_RADI(iz)-1
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
