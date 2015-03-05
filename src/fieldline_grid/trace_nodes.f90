!===============================================================================
! generate 3D grid by field line tracing from 2D base grids
!===============================================================================
subroutine trace_nodes()
  use emc3_grid
  use fieldline_grid
  use fieldline
  use grid
  use string
  use math
  implicit none

  type(t_grid) :: B
  real(real64) :: phi1, phi2, ts
  integer :: iz, nr1, nr2, np1, np2, tm, tc


  write (6, *)
  write (6, 1000)
  write (6, *)

  ! 1. setup toroidal discretization
  ! -> this is done in setup_toroidal_discretization in module fieldline_grid
  ! check setup
!  do ib=0,blocks-1
!     do it=0,Block(ib)%nt
!        write (6, *) ib, it, Block(ib)%phi(it)
!     enddo
!  enddo


  ! set parameters for field line tracing
  ts = pi2 / 3600.d0
  tm = NM_AdamsBashforth4
  tc = FL_ANGLE


  do iz=0,NZONET-1
     call B%load('base_grid_'//trim(str(iz))//'.dat')

     ! check input
     ! 1. layout
     if (B%layout /= MESH_2D) then
        write (6, *) 'error: unexpected grid layout ', B%layout
        stop
     endif
     if (B%fixed_coord /= 3) then
        write (6, *) 'error: grid is not located at fixed toroidal location'
        stop
     endif

     ! 2. toroidal position
     phi1 = Zone(iz)%phi(Zone(iz)%it_base)
     phi2 = B%fixed_coord_value / pi * 180.d0
     if (abs((phi1-phi2)/phi2) > 1.d-6) then
        write (6, *) 'error: unexpected toroidal position of base grid: ', phi2
        write (6, *) 'expected position: ', phi1
        stop
     endif

     ! 3.1 radial resolution
     nr1 = R_SURF_PL_TRANS_RANGE(1,iz)
     nr2 = R_SURF_PL_TRANS_RANGE(2,iz)
     if (B%n1-1 /= nr2-nr1) then
        write (6, *) 'error: mismatching radial resolution: ', B%n1
        write (6, *) 'expected index range for aligned grid: ', nr1, '->', nr2
        stop
     endif
     ! 3.2 poloidal resolution
     np1 = P_SURF_PL_TRANS_RANGE(1,iz)
     np2 = P_SURF_PL_TRANS_RANGE(2,iz)
     if (B%n2-1 /= np2-np1) then
        write (6, *) 'error: mismatching poloidal resolution: ', B%n2
        write (6, *) 'expected index range for aligned grid: ', np1, '->', np2
        stop
     endif


     call trace_nodes_1zone()
  enddo

 1000 format(3x,'- Tracing fieldlines from base grids ...')
  contains
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
! input:
!    iz		zone number
!    Zone(iz)   zone information structure
!    B          base grid nodes
!
! output:
!    PHI_PLANE(PHI_PL_OS(iz):PHI_PL_OS(iz+1)-1)
!    RG       (GRID_PL_OS(iz):GRID_PL_OS(iz+1)-1)
!    ZG       (GRID_PL_OS(iz):GRID_PL_OS(iz+1)-1)
!-----------------------------------------------------------------------
  subroutine trace_nodes_1zone()

  type(t_fieldline) :: F
  real(real64) :: y0(3), y1(3), Dphi
  integer :: ir0, ir, ip0, ip, it, idir, it_end, ig, ierr


  ! 1. setup toroidal segments
  do it=0,Zone(iz)%nt
     PHI_PLANE(it + PHI_PL_OS(iz)) = Zone(iz)%phi(it) / 180.d0 * pi
  enddo


  ! 2. field line tracing from base grid
  ir0 = R_SURF_PL_TRANS_RANGE(1,iz)
  ip0 = P_SURF_PL_TRANS_RANGE(1,iz)
  write (6, *) 'progress:'
  ! loop over all grid nodes
  do ir=0,B%n1-1
     write (6, *) ir, ' / ', B%n1-1
     do ip=0,B%n2-1
        y0(1:2) = B%mesh(ir,ip,:)
        y0(3)   = B%fixed_coord_value

        ! set base nodes
        it = Zone(iz)%it_base
        ig = ir0 + ir + (ip0 + ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        RG(ig) = y0(1)
        ZG(ig) = y0(2)

        ! idir = -1: negative (clockwise) toroidal direction
        !         1: positive (counter-clockwise) toroidal direction
        do idir=-1,1,2
           ! initialize field line at grid node
           call F%init(y0, idir*ts, tm, tc)

           it_end = 0
           if (idir > 0) it_end = Zone(iz)%nt
           do it=Zone(iz)%it_base+idir,it_end,idir
              Dphi = abs(Zone(iz)%phi(it) - Zone(iz)%phi(it-idir)) / 180.d0 * pi
              call F%trace_Dphi(Dphi, .false., y1, ierr)
              if (ierr .ne. 0) then
                 write (6, *) 'error in subroutine trace_grid: ', &
                              'trace_Dphi returned error ', ierr
                 stop
              endif

              ig = ir0 + ir + (ip0 + ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
              RG(ig) = y1(1)
              ZG(ig) = y1(2)
           enddo
        enddo
     enddo
  enddo

  end subroutine trace_nodes_1zone
!-----------------------------------------------------------------------
end subroutine trace_nodes
