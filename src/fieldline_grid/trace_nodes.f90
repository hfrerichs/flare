!===============================================================================
! generate 3D magnetic field aligned grid from base grids (slices) by field line
! tracing
!===============================================================================
subroutine trace_slices()
  use emc3_grid
  use field_aligned_grid
  implicit none

  integer :: iz


  do iz=0,NZONET-1
     call trace_slice(iz, TD(iz+1), base_grid(iz))
  enddo

  contains
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! input:
!    iz		zone number
!    tZ		zone information structure
!    B          base grid nodes
!
! output:
!    PHI_PLANE(PHI_PL_OS(iz):PHI_PL_OS(iz+1)-1)
!    RG       (GRID_PL_OS(iz):GRID_PL_OS(iz+1)-1)
!    ZG       (GRID_PL_OS(iz):GRID_PL_OS(iz+1)-1)
!-----------------------------------------------------------------------
  subroutine trace_slice (iz, tZ, B)
  use iso_fortran_env
  use field_aligned_grid, only: t_zone
  use emc3_grid
  use grid
  use fieldline
  implicit none

  integer, intent(in) :: iz
  type(t_zone), intent(in) :: tZ
  type(t_grid), intent(in) :: B


  integer, parameter :: &
     RADI_POLO = 1, &
     POLO_RADI = 2
  integer, parameter :: order = RADI_POLO


  type(t_fieldline) :: F
  real(real64) :: ts, y0(3), y1(3), Dphi
  integer :: i, j, jdir, tm, tc, ierr, ir, ip, it, ig, nr, np


  ! sanity check
  nr = SRF_RADI(iz)
  np = SRF_POLO(iz)
  if (B%nodes() .ne. nr*np) then
     write (6, *) 'error in subroutine trace_grid: inconsistent grid resolution!'
     write (6, *) B%nodes(), ' <> ', nr, ' * ', np
     stop
  endif



  ! setup toroidal segments
  do jdir=-1,1,2
     do j=0,tZ%nt(jdir)
        it = jdir*j + tZ%nt(-1)
        PHI_PLANE(it + PHI_PL_OS(iz)) = tZ%phi(jdir*j) / 180.d0 * pi
     enddo
  enddo


  ts = pi2 / 3600.d0
  tm = NM_AdamsBashforth4
  tc = FL_ANGLE

  write (6, *) 'progress:'
  ! loop over all grid nodes
  do i=1,B%nodes()
     !write (6, *) i, '/', B%nodes()
     y0 = B%node(i)

     ! get radial and poloidal grid indices
     select case (order)
     case(RADI_POLO)
        ip = (i-1) / nr
        ir = (i-1) - ip*nr
     case(POLO_RADI)
        ir = (i-1) / np
        ip = (i-1) - ir*np
     end select
     write (6, *) ir, ip, iz


     ! set base nodes
     it = tZ%nt(-1)
     ig = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
     RG(ig) = y0(1)
     ZG(ig) = y0(2)


     ! jdir = -1: negative (clockwise) toroidal direction
     !         1: positive (counter-clockwise) toroidal direction
     do jdir=-1,1,2
        ! initialize field line at grid node
        call F%init(y0, jdir*ts, tm, tc)

        do j=1,tZ%nt(jdir)
           it = jdir*j + tZ%nt(-1)
           Dphi = abs(tZ%phi(jdir*j) - tZ%phi(jdir*(j-1))) / 180.d0 * pi
           call F%trace_Dphi(Dphi, .false., y1, ierr)
           if (ierr .ne. 0) then
              write (6, *) 'error in subroutine trace_grid: ', &
                           'trace_Dphi returned error ', ierr
              stop
           endif

           ig = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
           RG(ig) = y1(1)
           ZG(ig) = y1(2)
        enddo
     enddo
  enddo
  write (6, *) 'done'

  end subroutine trace_slice
!-----------------------------------------------------------------------

end subroutine trace_slices
