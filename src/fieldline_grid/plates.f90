!===============================================================================
! Generate plate definitions for finite flux tube grids (for EMC3-EIRENE)
!===============================================================================
module plates
  use iso_fortran_env
  use emc3_grid
  use boundary
  implicit none

  character(len=*), parameter :: plates_file = 'plates.dat'

  !integer, dimension(:), allocatable :: ID_TEM


  contains
!=======================================================================



!=======================================================================
  subroutine initialize_plates()

  if (allocated(ID_TEM)) deallocate(ID_TEM)
  allocate (ID_TEM(0:MESH_P_OS(NZONET)-1))
  ID_TEM = 0

  end subroutine initialize_plates
!=======================================================================



!=======================================================================
  subroutine map_slice_to_target(iz, it0)
  use boundary
  integer, intent(in) :: iz, it0

  real(real64) :: r1(3), r2(3), X(3), tau
  integer :: ir, ip, it, idir, it1(-1:1), id, ig(4), ic


  write (6, 1000) it0, iz
 1000 format(2x,'- mapping toroidal slice ',i0,' in zone ',i0,' to target plates')
  write (6, *) R_SURF_PL_TRANS_RANGE(1:2,iz)
  write (6, *) P_SURF_PL_TRANS_RANGE(1:2,iz)

  it1(-1) = 0
  it1( 1) = ZON_TORO(iz)

  do ir=R_SURF_PL_TRANS_RANGE(1,iz),R_SURF_PL_TRANS_RANGE(2,iz)-1
  do ip=P_SURF_PL_TRANS_RANGE(1,iz),P_SURF_PL_TRANS_RANGE(2,iz)-1
     ig(1)   = ir + (ip + it0*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
     ig(2)   = ig(1) + 1
     ig(3)   = ig(2) + SRF_RADI(iz)
     ig(4)   = ig(3) - 1
     r1(1)   = sum(RG(ig))/4.d0
     r1(2)   = sum(ZG(ig))/4.d0
     r1(3)   = PHI_PLANE(it0 + PHI_PL_OS(iz))

     do idir=-1,1,2
     ! 1. scan through flux tube until central field line hits a target plate
     it1(0) = -1
     do it=it0,it1(idir)-idir,idir
        ig(1)   = ir + (ip + (it+idir)*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        ig(2)   = ig(1) + 1
        ig(3)   = ig(2) + SRF_RADI(iz)
        ig(4)   = ig(3) - 1
        r2(1)   = sum(RG(ig))/4.d0
        r2(2)   = sum(ZG(ig))/4.d0
        r2(3)   = PHI_PLANE(it+idir + PHI_PL_OS(iz))

        ! intersect boundary at slice it
        if (intersect_boundary(r1, r2, X, id, tau=tau)) then
!           ! mark toroidal surface closest to target plate
!           if (tau < 0.5d0) then
!              it1(0) = it
!           else
!              it1(0) = it+idir
!           endif
           ! mark toroidal surface beyond intersection with plate (except for 1st and last slice)
           it1(0) = it+idir
           if (it1(0) == it1(idir)) it1(0) = it
           exit
        endif

        r1 = r2
     enddo

     ! 2. mark remaining cells in flux tube
     if (it1(0) > 0  .and.  it1(0) < ZON_TORO(iz)) then
        do it=it1(0),it1(idir)-idir,idir
           ic = ir+(ip+it*ZON_POLO(iz))*ZON_RADI(iz)+MESH_P_OS(iz)
           if (idir < 0) ic = ic - ZON_POLO(iz)*ZON_RADI(iz)
           ID_TEM(ic) = 1
        enddo
     endif
     enddo
  enddo
  enddo

  end subroutine map_slice_to_target
!=======================================================================



!=======================================================================
  subroutine minimum_radial_coordinate(iz)
  use curve2D
  use run_control, only: Debug
  integer, intent(in) :: iz

  real(real64), dimension(:,:), allocatable :: rtouch

  integer, parameter :: nts = 10

  character(len=80)  :: filename
  type(t_curve) :: C
  real(real64)  :: t, phi, rmin, x1(2), x2(2), th
  integer       :: ir, ip, it, its, ib, ig(4)


  write (6, 1000) iz

  allocate (rtouch(0:ZON_POLO(iz)-1, 0:ZON_TORO(iz)-1))
  rtouch = SRF_RADI(iz)-1

  do ib=2,n_boundary1
     if (boundary_in_zone(ib, PHI_PLANE(PHI_PL_OS(iz)), PHI_PLANE(PHI_PL_OS(iz+1)-1))) then
        write (6, *) boundary_label(ib)
     else
        cycle
     endif


  ! toroidal loop over cells
  do it=0,ZON_TORO(iz)-1
  ! sub-resolution in cells
  do its=1,nts
     t   = (float(its)-0.5d0) / nts
     phi = (1.d0-t) * PHI_PLANE(it+PHI_PL_OS(iz))  +  t * PHI_PLANE(it+1+PHI_PL_OS(iz))

     !do ib=1,n_boundary1
     !do ib=2,n_boundary1
        C = boundary_slice(ib, phi)
        if (C%n_seg < 0) cycle
        if (Debug) then
           write (filename, 9000) iz, it, its, ib
           call C%plot(filename=filename)
        endif

        ! set up grid coordinates for boundary segment ib
        do ip=0,ZON_POLO(iz)-1
           ir    = 0
           ig(1) = (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
           ig(2) = ig(1) + SRF_RADI(iz)
           ig(3) = ig(1) + SRF_RADI(iz) * SRF_POLO(iz)
           ig(4) = ig(3) + SRF_RADI(iz)

           ! middle point on lower radial cell boundary (ir, ip, it)
           x1(1) = (1.d0-t) * (RG(ig(1)) + RG(ig(2))) / 2.d0 &
                 +       t  * (RG(ig(3)) + RG(ig(4))) / 2.d0
           x1(2) = (1.d0-t) * (ZG(ig(1)) + ZG(ig(2))) / 2.d0 &
                 +       t  * (ZG(ig(3)) + ZG(ig(4))) / 2.d0

           do ir=1,SRF_RADI(iz)-1
              ig(1) = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
              ig(2) = ig(1) + SRF_RADI(iz)
              ig(3) = ig(1) + SRF_RADI(iz) * SRF_POLO(iz)
              ig(4) = ig(3) + SRF_RADI(iz)

              ! middle point on upper radial cell boundary (ir, ip, it)
              x2(1) = (1.d0-t) * (RG(ig(1)) + RG(ig(2))) / 2.d0 &
                    +       t  * (RG(ig(3)) + RG(ig(4))) / 2.d0
              x2(2) = (1.d0-t) * (ZG(ig(1)) + ZG(ig(2))) / 2.d0 &
                    +       t  * (ZG(ig(3)) + ZG(ig(4))) / 2.d0

              if (intersect_curve(x1, x2, C, th=th)) then
                 rtouch(ip,it) = min(rtouch(ip,it), th+ir-1)
                 rmin          = min(rmin, rtouch(ip,it))
              endif
              x1 = x2
           enddo
        enddo
     !enddo
  enddo
  enddo
  enddo


  if (Debug) then
     open  (99, file='DEBUG_RTOUCH_FLARE')
     do ip=0,ZON_POLO(iz)-1
     do it=0,ZON_TORO(iz)-1
     write (99, *) rtouch(ip,it)+1
     enddo
     enddo
     close (99)
  endif

  ! "upload" rtouch to ID_TEM
  do it=0,ZON_TORO(iz)-1
     write (6, 4000) it, minval(rtouch(:,it))
     do ip=0,ZON_POLO(iz)-1
        th = rtouch(ip,it)
        ir = th
        if (th-ir > 0.5d0) ir = ir + 1

        ib = (ip + it*ZON_POLO(iz)) * ZON_RADI(iz)  +  MESH_P_OS(iz)
        ID_TEM(ib+ir:ib+ZON_RADI(iz)-1) = 1
     enddo
  enddo


  ! cleanup
  deallocate (rtouch)

 1000 format(8x,'Zone ',i0,': intersect radial mesh coordinate with plates')
 4000 format(8x,i0,4x,f12.5)
 9000 format('DEBUG_PLATE_Z',i0,'_T',i0,'-',i0,'_B',i0,'.PLT')
  end subroutine minimum_radial_coordinate
!=======================================================================



!=======================================================================
! kmin: minimum number of plasma cells between target cells
!=======================================================================
  subroutine write_plates_format1(kmin)
  integer, intent(in) :: kmin

  integer, parameter :: iu = 78

  integer, dimension(:), allocatable :: KBEG, KEND
  integer :: ir, ip, it, iz, ns, ic, k, k1, k2


  open  (iu, file=plates_file)
  do iz=0,NZONET-1
  allocate (KBEG(ZON_TORO(iz)), KEND(ZON_TORO(iz)))

  do ir=0,ZON_RADI(iz)-1
  do ip=0,ZON_POLO(iz)-1
     ! 1. set up KBEG and KEND from ID_TEM
     ns = 0
     do it=0,ZON_TORO(iz)-2
        ic = ir + (ip + it * ZON_POLO(iz)) * ZON_RADI(iz)  +  MESH_P_OS(iz)

        ! first cell is already plate cell
        if (it == 0  .and.  ID_TEM(ic) == 1) then
           ns       = ns + 1
           KBEG(ns) = it
        endif

        ! transition between plasma and plate cells
        if (ID_TEM(ic) /= ID_TEM(ic + ZON_POLO(iz)*ZON_RADI(iz))) then
           ! plasma -> target
           if (ID_TEM(ic) == 0) then
              ns       = ns + 1
              KBEG(ns) = it + 1

           ! target -> plasma
           else
              KEND(ns) = it
           endif
        endif
     enddo
     ! last cell is plate cell
     it = ZON_TORO(iz)-1
     ic = ir + (ip + it * ZON_POLO(iz)) * ZON_RADI(iz)  +  MESH_P_OS(iz)
     if (ID_TEM(ic) == 1) then
        KEND(ns) = it
     endif


     ! 2. adjust definition for kmin
     if (ns > 0) then
        ! update KBEG, KEND
        k2 = ns
        k1 = 1
        do
           if (k1 >= k2) exit

           if (KBEG(k1+1) - KEND(k1) < kmin) then
              KEND(k1) = KEND(k1+1)
              k2       = k2 - 1
              do k=k1+2,k2
                 KBEG(k) = KBEG(k+1)
                 KEND(k) = KEND(k+1)
              enddo
           else
              k1 = k1 + 1
           endif
        enddo
        ns = k2

        ! update ID_TEM
        do k=1,ns
           do it=KBEG(k),KEND(k)
              ic = ir + (ip + it * ZON_POLO(iz)) * ZON_RADI(iz)  +  MESH_P_OS(iz)
              ID_TEM(ic) = 1
           enddo
        enddo
     endif


     ! 3. write plate definitions
     if (ns > 0) then
        write (iu, 1000) iz, ir, ip, 2*ns, (KBEG(k), KEND(k), k=1,ns)
     endif
  enddo
  enddo

  deallocate (KBEG, KEND)
  enddo
  close (iu)

 1000 format(1x,i0,1x,i4,1x,i4,1x,i4,1x,i5,1x,i5)
  end subroutine write_plates_format1
!=======================================================================



!=======================================================================
  subroutine write_plates_format3()

  integer, parameter :: iu = 78


  open  (iu, file=plates_file)
  close (iu)

  end subroutine write_plates_format3
!=======================================================================



!=======================================================================
  subroutine write_plates()
  use emc3_grid

  integer, parameter :: output_format = 1
  integer, parameter :: iu = 78

  integer, dimension(:), allocatable :: ntcell, ntcell_bundle
  integer :: iz, ir, ip, it, ig, np, np_bundle


  allocate (ntcell(maxval(ZON_TORO)), ntcell_bundle(maxval(ZON_TORO)))


  open  (iu, file='plates.dat')
  do iz=0,NZONET-1
     do ir=R_SURF_PL_TRANS_RANGE(1,iz),R_SURF_PL_TRANS_RANGE(2,iz)-1
     do ip=P_SURF_PL_TRANS_RANGE(1,iz),P_SURF_PL_TRANS_RANGE(2,iz)-1
        ! check number of plate cells in flux tube
        np     = 0
        ntcell = 0
        do it=0,ZON_TORO(iz)
           ig = ir+(ip+it*ZON_POLO(iz))*ZON_RADI(iz)+MESH_P_OS(iz)
           if (ID_TEM(ig) > 0) then
              np         = np + 1
              ntcell(np) = it
           endif
        enddo

        ! there is at least one plate cell
        if (np > 0) then
           select case (output_format)
           ! 1. bundle plate cells
           case(1)
              ! setup bundles of plate cells
              np_bundle     = 0
              it            = 0
              ntcell_bundle = 0
              flux_tube_loop: do 
                 ig = ir+(ip+it*ZON_POLO(iz))*ZON_RADI(iz)+MESH_P_OS(iz)

                 ! find 1st plate cell in bundle
                 if (ID_TEM(ig) > 0) then
                    np_bundle = np_bundle + 1
                    ntcell_bundle(np_bundle) = it

                    ! find last plate cell in bundle
                    bundle_loop: do
                       ! reached end of flux tube
                       if (it == ZON_TORO(iz)-1) exit

                       it = it + 1
                       ig = ir+(ip+it*ZON_POLO(iz))*ZON_RADI(iz)+MESH_P_OS(iz)
                       if (ID_TEM(ig) == 0) exit
                    enddo bundle_loop

                    np_bundle = np_bundle + 1
                    ntcell_bundle(np_bundle) = it
                 endif

                 it = it + 1
                 if (it >= ZON_TORO(iz)) exit
              enddo flux_tube_loop
              write (iu, 1001) iz, ir, ip, np_bundle, ntcell_bundle(1:np_bundle)

           ! 3. write each toroidal cell index for plate cells
           case(3)
              write (iu, *) iz, ir, ip, np, ntcell(1:np)

           case default
              write (6, 9001) output_format
           end select
        endif
     enddo
     enddo
  enddo
  close (iu)

  deallocate (ntcell, ntcell_bundle)
 1001 format(1x,i0,1x,i4,1x,i4,1x,i4,1x,i5,1x,i5)
 9001 format('error: invalid output format ',i0,' for plate cells!')
  end subroutine write_plates
!=======================================================================

end module plates
!=======================================================================





!=======================================================================
! GENERATE PLATE DEFINITIONS
!=======================================================================
  subroutine generate_plates_old()
  use iso_fortran_env
  use fieldline_grid, only: symmetry
  use emc3_grid
  use boundary, unused => outside_boundary
  use curve2D
  use math
  use run_control, only: Debug
  use string
  use dataset
  implicit none

  integer, parameter :: iu = 78
  real(real64), parameter :: l0 = 10.d-10

  type(t_curve), dimension(:,:), allocatable :: C

  real(real64), dimension(:), allocatable   :: RC_TEM, ZC_TEM
  !integer, dimension(:), allocatable   :: ID_TEM
  integer, dimension(:,:), allocatable :: iindex ! number of cells behind a plate
  integer, dimension(:), allocatable   :: knumb  ! cell index in flux tube for plate cells

  logical      :: plate_cell
  real(real64) :: x(3), phi
  integer      :: nr, np, nt, iz, i, j, k, l, l1, l2, irun, icut, ig(8), ic


  write (6, 9000)
 9000 format(3x,'- Set up plate surfaces')


  allocate (RC_TEM(0:MESH_P_OS(NZONET)-1))
  allocate (ZC_TEM(0:MESH_P_OS(NZONET)-1))
  ID_TEM = 0
  open  (iu, file='plates.dat')
  do iz=0,NZONET-1
     nr = ZON_RADI(iz)
     np = ZON_POLO(iz)
     nt = ZON_TORO(iz)
     allocate (iindex(0:nr-1, 0:np-1))
     iindex = 0

     write (6, 1000) iz, nr, np, nt
 1000 format (8x,'zone ',i0,' (',i0,' x ',i0,' x ',i0,')')

     ! setup slices for Q4-type surfaces
     if (n_quad > 0) then
     allocate (C(0:nt-1, n_quad))
     do k=0,nt-1
     do l=1,n_quad
         phi  = (PHI_PLANE(k+1+PHI_PL_OS(iz)) + PHI_PLANE(k+PHI_PL_OS(iz))) / 2.d0
         phi  = phi_sym(phi, symmetry)
         C(k,l) = S_quad(l)%slice(phi)

         if (Debug) then
            call C(k,l)%plot(filename='debug/Q4surf_'//trim(str(l))//'_zone'//trim(str(iz))// &
                                      '_t'//trim(str(k)))
         endif
     enddo
     enddo
     endif


     ! (0) count plate cells, (1) setup plate cells
     do irun=0,1
        icut = 0
        ! loop over all flux tubes
        do i=R_SURF_PL_TRANS_RANGE(1,iz),R_SURF_PL_TRANS_RANGE(2,iz)-1
        do j=P_SURF_PL_TRANS_RANGE(1,iz),P_SURF_PL_TRANS_RANGE(2,iz)-1
        do k=0,nt-1
           ig(1)   = i + (j + k*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
           ig(2)   = ig(1) + 1
           ig(3)   = ig(2) + SRF_RADI(iz)
           ig(4)   = ig(3) - 1
           ig(5:8) = ig(1:4) + SRF_POLO(iz)*SRF_RADI(iz)

           x(1)    = sum(RG(ig))/8.d0
           x(2)    = sum(ZG(ig))/8.d0
           x(3)    = (PHI_PLANE(k+1+PHI_PL_OS(iz)) + PHI_PLANE(k+PHI_PL_OS(iz))) / 2.d0
           x(3)    = phi_sym(x(3), symmetry)

           ic      = i + (j + k*ZON_POLO(iz))*ZON_RADI(iz) + MESH_P_OS(iz)
           RC_TEM(ic) = x(1)
           ZC_TEM(ic) = x(2)

           plate_cell = outside_boundary()

           if (plate_cell) then
              ID_TEM(ic) = 1
              icut = icut + 1
              if (irun == 1) then
                 iindex(i,j) = iindex(i,j) + 1
                 knumb(icut) = k
              endif
           endif
        enddo
        enddo
        enddo

        ! allocate knumb array after 1st run
        if (irun == 0) then
           allocate (knumb(icut))
           knumb = 0

        ! write plate cells after 2nd run
        else
           icut = 0
           do i=0,nr-1
           do j=0,np-1
              if (iindex(i,j) .ne. 0) then
                 l1 = icut + 1
                 l2 = icut + iindex(i,j)
                 write (iu, *) iz, i, j, iindex(i,j),(knumb(l),l=l1,l2)
                 icut = icut + iindex(i,j)
              endif
           enddo
           enddo

           write (6, 2000) icut
           2000 format (8x,i0,' plate cells')
           deallocate (knumb, iindex)
        endif
     enddo

     ! cleanup slices of Q4-type surfaces
     if (n_quad > 0) then
     do k=0,nt-1
     do l=1,n_quad
        call C(k,l)%destroy()
     enddo
     enddo
     deallocate (C)
     endif
  enddo
  close (iu)


  ! additional output for debugging
  if (Debug) then
     do iz=0,NZONET-1
     do k=0,ZON_TORO(iz)-1
        open  (99, file='debug/SlicePlasma_zone'//trim(str(iz))//'_t'//trim(str(k))//'.plt')
        open  (98, file='debug/SliceBoundary_zone'//trim(str(iz))//'_t'//trim(str(k))//'.plt')

        do i=R_SURF_PL_TRANS_RANGE(1,iz),R_SURF_PL_TRANS_RANGE(2,iz)-1
        do j=P_SURF_PL_TRANS_RANGE(1,iz),P_SURF_PL_TRANS_RANGE(2,iz)-1
           ic      = i + (j + k*ZON_POLO(iz))*ZON_RADI(iz) + MESH_P_OS(iz)
           if (ID_TEM(ic) == 1) then
              write (98, *) RC_TEM(ic), ZC_TEM(ic)
           else
              write (99, *) RC_TEM(ic), ZC_TEM(ic)
           endif
        enddo
        enddo
        close (99)
        close (98)
     enddo
     enddo

     call plate_check_all()
  endif
  deallocate (ID_TEM, RC_TEM, ZC_TEM)
  contains
  !-------------------------------------------------------------------
  function outside_boundary()
  ! input:
  !    k: toroidal index
  !    C: slices of Q4-type surfaces at phi(k)
  !    x: reference point
  logical :: outside_boundary, side


  ! set default
  outside_boundary = .false.


  ! check axisymmetric (L2-type) surfaces
  do l=1,n_axi
     side = boundary_side(l) == 1
     if (S_axi(l)%outside(x(1:2)) .eqv. side) then
        outside_boundary = .true.
        return
     endif
  enddo

  ! check Q4-type surfaces
  do l=1,n_quad
     if (C(k,l)%n_seg <= 0) cycle
     side = boundary_side(n_axi + n_block + l) == 1
     if (C(k,l)%outside(x(1:2)) .eqv. side) then
        outside_boundary = .true.
        return
     endif
  enddo

  ! check block limiters (CSG-type)
  do l=1,n_block
     if (bl_outside(l, x)) then
        outside_boundary = .true.
        return
     endif
  enddo

  end function outside_boundary
  !-------------------------------------------------------------------


  !-------------------------------------------------------------------
  subroutine outside_boundary_check

  integer, parameter :: iu = 99, iu2 = 98, n = 101

  character(len=72) :: s
  real(real64) :: x(2), Phi
  integer :: i

  allocate (C(0:n-1, n_quad))
  open  (iu, file='fl+.dat')
  open  (iu2, file='fl+B.dat')
  read  (iu, *) s
  do i=0,n-1
     read (iu, *) x, Phi

     phi  = phi_sym(phi, symmetry)
     C(i,1) = S_quad(1)%slice(phi)

     if (C(i,1)%outside(x)) then
        !write (iu2, *) Phi, 1
     else
        !write (iu2, *) Phi, 0
        write (iu2, *) x, Phi
     endif
  enddo
  close (iu)
  close (iu2)

  stop
  end subroutine outside_boundary_check
  !-------------------------------------------------------------------


  !-------------------------------------------------------------------
  subroutine plate_check_all
  integer :: iz, k, j, i, ic1, ic2, iplate


  write (6, *) 'running plate checks ...'

  ID_TEM = ID_TEM*2 - 1
  iplate = 0

  do iz=0,NZONET-1
  do i=0,ZON_RADI(iz)-1
  do j=0,ZON_POLO(iz)-1
     do k=0,ZON_TORO(iz)-2
        ic1 = i + (j +  k   *ZON_POLO(iz))*ZON_RADI(iz) + MESH_P_OS(iz)
        ic2 = i + (j + (k+1)*ZON_POLO(iz))*ZON_RADI(iz) + MESH_P_OS(iz)

        if (ID_TEM(ic1)*ID_TEM(ic2) < 0) then
           iplate = iplate + 1
           write (6, *) iplate
           call plate_check(iz,i,j,k+1)
        endif
     enddo
  enddo
  enddo
  enddo

  end subroutine plate_check_all
  !-------------------------------------------------------------------
  subroutine plate_check(iz, ir, ip, jt)
  use Q4
  use math
  use grid
  use fieldline
  use dataset
  integer, intent(in) :: iz, ir, ip, jt

  real(real64), parameter :: Limit = 360.d0

  type(t_Q4)   :: Q
  type(t_grid) :: G
  type(t_fieldline) :: F
  type(t_dataset)   :: D
  real(real64) :: x1(2), x2(2), x3(2), x4(2), phi, y(3), ts, Lc, Lcsav, dl
  integer      :: ig(4), i, idir, n, nsuccess


  ! calculate node indices
  ig(1)   = ir + (ip + jt*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
  ig(2)   = ig(1) + 1
  ig(3)   = ig(2) + SRF_RADI(iz)
  ig(4)   = ig(3) - 1

  ! get node coordinates
  x1(1)   = RG(ig(1)); x1(2)   = ZG(ig(1))
  x2(1)   = RG(ig(2)); x2(2)   = ZG(ig(2))
  x3(1)   = RG(ig(3)); x3(2)   = ZG(ig(3))
  x4(1)   = RG(ig(4)); x4(2)   = ZG(ig(4))
  phi     = PHI_PLANE(jt+PHI_PL_OS(iz))

  ! setup quadrilateral
  call Q%set_nodes(x1, x2, x3, x4)

  ! generate mesh
  n = 10
  G = Q%generate_mesh(n, n, phi)
  G%coordinates       = CYLINDRICAL
  !G%fixed_coord       = 3
  !G%fixed_coord_value = phi / 180.d0 * pi

  ! sample connection length on grid
  ts       = 1.d0
  Lcsav    = 0.d0
  nsuccess = 0
  call D%new(G%nodes(),5)
  do i=1,G%nodes()
     y = G%node(i)

     D%x(i,1) = y(1)
     D%x(i,2) = y(2)

     ! trace field line in both directions
     do idir=-1,1,2
        call F%init(y, idir*ts, NM_AdamsBashforth4, FL_ARC)

        trace_loop: do
           dl = F%trace_1step()
           Lc = F%phi_int * 180.d0 / pi

           if (abs(Lc) > Limit) exit trace_loop

           if (F%intersect_boundary()) exit trace_loop
        enddo trace_loop

        D%x(i,3 + (idir+1)/2) = Lc
     enddo
     D%x(i,5) = min(abs(D%x(i,3)), abs(D%x(i,4)))

     ! success frequency, average shortest toroidal distance
     if (D%x(i,5) < Limit) then
        nsuccess = nsuccess + 1
        Lcsav    = Lcsav    + D%x(i,5)
     endif
  enddo
  Lcsav = Lcsav / nsuccess
  write (6, *) 'success frequency = ', 100.d0 * nsuccess / G%nodes(), ' %'
  write (6, *) 'average toroidal distance to plates [deg] = ', Lcsav
  call G%store('plate1.grid')
  call D%plot(filename='plate1.dat')

  call D%destroy()
  call G%destroy()
  stop

  end subroutine plate_check
  !-------------------------------------------------------------------
  end subroutine generate_plates_old
!=======================================================================




!=======================================================================
  subroutine generate_plates()
  use fieldline_grid
  use emc3_grid
  use plates
  implicit none

  integer :: iz, kmin


  write (6, *)
  write (6, 1000)
  write (6, *)

  call initialize_plates()
  select case(plate_generator)
  case(PLATES_DEFAULT)
     call generate_plates_old()

  case(PLATES_RADIAL_INTERSECT)
     do iz=0,NZONET-1
        call minimum_radial_coordinate(iz)
     enddo
     kmin = 2
     call write_plates_format1(kmin)

  case default
     write (6, *) 'error: invalid plate generator ', trim(plate_generator)
     stop
  end select
!  call map_slice_to_target(0, ZON_TORO(0))
!  call write_plates()

 1000 format(3x,'- Generating plate surface approximation within flux-tube mesh')
  end subroutine generate_plates
!=======================================================================
