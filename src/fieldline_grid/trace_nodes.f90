!===============================================================================
! generate 3D grid by field line tracing from 2D base grids
!===============================================================================
subroutine trace_nodes(post_process_grid)
  use emc3_grid, only: NZONET
  use fieldline_grid
  use fieldline
  use grid
  use string
  use math
  implicit none

  interface
    subroutine post_process_grid()
    end subroutine post_process_grid
  end interface

  type(t_grid), dimension(:), allocatable :: G3D
  type(t_grid) :: B
  real(real64) :: phi1, phi2, ts
  integer :: iz, nr1, nr2, np1, np2, tm, tc


  write (6, *)
  write (6, 1000)
  write (6, *)


  ! set parameters for field line tracing
  ts = pi2 / 3600.d0
  tm = NM_AdamsBashforth4
  tc = FL_ANGLE


  allocate (G3D(0:NZONET-1))
  do iz=0,NZONET-1
     ! 1. load base grid
     call B%load('base_grid_'//trim(str(iz))//'.dat')


     ! 2. check input --------------------------------------------------
     ! 2.1 layout
     if (B%layout /= MESH_2D) then
        write (6, *) 'error: unexpected grid layout ', B%layout
        stop
     endif
     if (B%fixed_coord /= 3) then
        write (6, *) 'error: grid is not located at fixed toroidal location'
        stop
     endif

     ! 2.2 toroidal position of base grid
     phi1 = Zone(iz)%phi(Zone(iz)%it_base)
     phi2 = B%fixed_coord_value / pi * 180.d0
     if (abs((phi1-phi2)/phi2) > 1.d-6) then
        write (6, *) 'error: unexpected toroidal position of base grid: ', phi2, ' deg'
        write (6, *) 'expected position: ', phi1, ' deg'
        stop
     endif
     ! -----------------------------------------------------------------


     ! 3. generate 3D grid
     call trace_nodes_1zone(G3D(iz))
  enddo


  !call post_process_grid()


  ! store fieldlines
  do iz=0,NZONET-1
     call G3D(iz)%store('fieldlines_'//trim(str(iz))//'.dat')
  enddo


  ! cleanup
  deallocate (G3D)
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
!    G          slices with nodes on field lines
!-----------------------------------------------------------------------
  subroutine trace_nodes_1zone(G)
  type(t_grid), intent(out) :: G

  type(t_fieldline) :: F
  real(real64) :: y0(3), y1(3), Dphi
  integer :: ir, ip, it, idir, it_end, ig, ierr


  nt = Zone(iz)%nt
  call G%new(CYLINDRICAL, MESH_3D, 3, B%n1, B%n2, nt+1)


  ! 1. setup toroidal segments
  do it=0,nt
     G%x3(it) = Zone(iz)%phi(it) / 180.d0 * pi
  enddo


  ! 2. field line tracing from base grid
  write (6, *) 'progress:'
  ! loop over all grid nodes
  do ir=0,B%n1-1
     write (6, *) ir, ' / ', B%n1-1
     do ip=0,B%n2-1
        y0(1:2) = B%mesh(ir,ip,:)
        y0(3)   = B%fixed_coord_value

        ! set base nodes
        it = Zone(iz)%it_base
        G%mesh3D(ir,ip,it,1:2) = y0(1:2)

        ! idir = -1: negative (clockwise) toroidal direction
        !         1: positive (counter-clockwise) toroidal direction
        do idir=-1,1,2
           ! initialize field line at grid node
           call F%init(y0, idir*ts, tm, tc)

           it_end = 0
           if (idir > 0) it_end = Zone(iz)%nt
           ! trace field line throughout zone
           do it=Zone(iz)%it_base+idir,it_end,idir
              Dphi = abs(Zone(iz)%phi(it) - Zone(iz)%phi(it-idir)) / 180.d0 * pi
              call F%trace_Dphi(Dphi, .false., y1, ierr)
              if (ierr .ne. 0) then
                 write (6, *) 'error in subroutine trace_grid: ', &
                              'trace_Dphi returned error ', ierr
                 stop
              endif

              G%mesh3D(ir,ip,it,1:2) = y1(1:2)
           enddo
        enddo
     enddo
  enddo

  end subroutine trace_nodes_1zone
!-----------------------------------------------------------------------
end subroutine trace_nodes
