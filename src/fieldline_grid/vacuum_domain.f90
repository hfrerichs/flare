!===============================================================================
! Set up additional domain for EIRENE (vacuum region in far SOL and PFR)
!===============================================================================
subroutine vacuum_domain_for_EIRENE()
  use emc3_grid
  use fieldline_grid
  use divertor
  implicit none

  integer :: iz


  write (6, 1000)
 1000 format(3x,' - Setting up vacuum domain for EIRENE')
  do iz=0,NZONET-1
     ! set up far SOL
     if (Zone(iz)%isfr(1) == SF_VACUUM) then
        call setup_vacuum_domain(iz, nr_EIRENE_vac, 1)
     endif
     if (Zone(iz)%isfr(2) == SF_VACUUM) then
        call setup_vacuum_domain(iz, nr_EIRENE_vac, 2)
     endif


     ! adjust last cell of divertor legs to close simulatin domain
     if (Zone(iz)%isfp(1) == SF_VACUUM  .and.  Zone(iz)%isfp(2) == SF_VACUUM) then
        call close_grid_domain(iz)
     endif
  enddo

end subroutine vacuum_domain_for_EIRENE
!===============================================================================



!===============================================================================
subroutine setup_vacuum_domain(iz, nr_vac, boundary)
  use iso_fortran_env
  use fieldline_grid
  implicit none

  character(len=*), parameter :: s_boundary(2) = (/ 'lower', 'upper' /)
  integer, intent(in) :: iz, nr_vac, boundary

  integer, parameter :: &
     SCALE_BOUNDARY   = 1, &
     CELL_EXTEND      = 2, &
     RAY_ADJUST_SCALE = 3, &
     MANUAL_2D        = -2, &
     MANUAL_3D        = -3

  real(real64) :: dl
  integer      :: Method, ir0, idir, ir2


  select case(Zone(iz)%N0_method)
  case('scale_boundary','')
     Method = SCALE_BOUNDARY
  case('cell_extend')
     Method = CELL_EXTEND
  case('ray_adjust_scale')
     Method = RAY_ADJUST_SCALE
  case('manual', 'manual_2D')
     Method = MANUAL_2D
  case('manual_3D')
     Method = MANUAL_3D
  case default
     write (6, *) 'error: invalid method ', trim(Zone(iz)%N0_method), ' for N0 domain!'
     stop
  end select


  dl     = Zone(iz)%d_N0
!  if (Zone(iz)%N0_file .ne. '') then
!     Method = MANUAL_2D
!  endif


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
  case (SCALE_BOUNDARY)
      call scale_boundary_proc(iz, ir0, idir, ir2, dl)
  case (CELL_EXTEND)
      call cell_extend_proc(iz, ir0, idir, ir2, dl)
  case (RAY_ADJUST_SCALE)
      call ray_adjust_scale_proc(iz, ir0, idir, ir2, dl)
  case (MANUAL_2D)
      call vacuum_domain_manual(iz, ir0, idir, ir2, Zone(iz)%N0_file)
  case (MANUAL_3D)
      call vacuum_domain_manual_3D(iz, ir0, idir, ir2, Zone(iz)%N0_file, Zone(iz)%N0_filter)
  end select


end subroutine setup_vacuum_domain
!===============================================================================


!===============================================================================
! scale boundary of plasma transport domain by dl
! keep node spacings
!===============================================================================
subroutine scale_boundary_proc(iz, ir0, idir, ir2, dl)
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
     it = ZON_TORO(iz) / 2
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


     ! adjust ends
     call adjust_boundary(C)
     if (Debug) then
        call C%plot(filename='debug/VacuumBoundary_'//trim(str(iz))//'_'//trim(str(it))//'.adjust')
     endif
     call C%setup_length_sampling()


     ! poloidal loop (set new grid nodes)
  do it=0,SRF_TORO(iz)-1
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
  enddo

     ! cleanup
     call C%destroy()
  deallocate (en)

  contains
  !---------------------------------------------------------------------
  subroutine adjust_boundary(C)
  type(t_curve), intent(inout) :: C

  character(len=*), dimension(2), parameter :: Send = (/ 'right', 'left ' /)

  type(t_curve) :: C_split(3)
  real(real64)  :: x1(2), x2(2), tsplit(2), xh(2)
  integer       :: j, isplit(2)


  ! find position on C from where the ends will be split off
  ip = 0
  do j=1,2
     ig0   = ir0 + (ip +it*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
     x1(1) = RG(ig0)
     x1(2) = ZG(ig0)
     ig0   = ir0-idir + (ip +it*SRF_POLO(iz))*SRF_RADI(iz) +GRID_P_OS(iz)
     x2(1) = RG(ig0)
     x2(2) = ZG(ig0)

     if (intersect_curve(x2, x1, C, xh=xh, sh=tsplit(j), ish=isplit(j), intersect_mode=1)) then
        !write (6, 1000) Send(j)
        !write (6, 1001) xh
        isplit(j) = isplit(j) - 1
        !write (6, *) 'ish = ', isplit(j)
        !write (6, *) 'sh  = ', tsplit(j)
     else
        !write (6, 1002) Send(j)
        isplit(j) = (j-1) * (C%n_seg-1)
        tsplit(j) = (j-1) * 1.d0
     endif

     ip = ZON_POLO(iz)
  enddo

  ! split off the ends of C
  !write (6, *) 'isplit = ', isplit
  call C%splitnseg(3, isplit, tsplit, C_split)
  !call C_split(1)%plot(filename='split_right.plt')
  !call C_split(2)%plot(filename='split_center.plt')
  !call C_split(3)%plot(filename='split_left.plt')
  call C%copy(C_split(2))

 1000 format(8x,'adjusting ',a5,' end: ')
 1001 format(10x,2f12.6)
 1002 format(8x,'no adjustment on ',a5,' end')
  end subroutine adjust_boundary
  !---------------------------------------------------------------------
end subroutine scale_boundary_proc
!===============================================================================



!===============================================================================
! create boundary by "extending" cells by dl
!===============================================================================
subroutine cell_extend_proc(iz, ir0, idir, ir2, dl)
  use iso_fortran_env
  use emc3_grid
  use curve2D
  use string
  implicit none

  integer, intent(in)      :: iz, ir0, idir, ir2
  real(real64), intent(in) :: dl

  real(real64), dimension(:,:), allocatable :: v
  type(t_curve) :: C, Cout, directional_extend
  real(real64)  :: rho
  integer       :: it, ip, ir, ir1, ig, ig0
  logical       :: debug


  !if (iz .ne. 1) return

  call C%new(ZON_POLO(iz))
  allocate (v(0:SRF_POLO(iz)-1,2))
  ir1 = ir0 + idir


  do it=0,SRF_TORO(iz)-1
     ! set up boundary nodes and direction vectors
     do ip=0,SRF_POLO(iz)-1
        ig = ir0 + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        C%x(ip,1) = RG(ig)
        C%x(ip,2) = ZG(ig)
        v(ip,1)   = RG(ig) - RG(ig-idir)
        v(ip,2)   = ZG(ig) - ZG(ig-idir)
        v(ip,:)   = v(ip,:) / sqrt(sum(v(ip,:)**2)) * dl
     enddo

     debug = .false.
     !if (iz == 4  .and.  it == 0) debug = .true.
     Cout = directional_extend(C, v, debug)
     !call Cout%plot(filename='test.plt')


     ! set up new grid nodes
     do ip=0,SRF_POLO(iz)-1
        ig0 = ir0 + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        do ir=ir1,ir2,idir
           rho = 1.d0 * (ir-ir0) / (ir2-ir1+idir)
           ig  = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)

           RG(ig) = RG(ig0) + rho * (Cout%x(ip,1) - RG(ig0))
           ZG(ig) = ZG(ig0) + rho * (Cout%x(ip,2) - ZG(ig0))
        enddo
     enddo
  enddo
  deallocate (v)


end subroutine cell_extend_proc


  function directional_extend(C, v, debug) result(Cout)
  use iso_fortran_env
  use curve2D
  type(t_curve), intent(in) :: C
  real(real64),  intent(in) :: v(0:C%n_seg,2)
  logical,       intent(in) :: debug

  integer,      dimension(:),   allocatable :: icheck
  integer,      dimension(:,:), allocatable :: icheck2
  type(t_curve) :: Cout
  real(real64)  :: L1(2), L2(2), M1(2), M2(2), l, m
  integer       :: n, ip


  ! initialize output curve
  n = C%n_seg
  call Cout%new(n)
  do ip=0,n
     Cout%x(ip,:) = C%x(ip,:) + v(ip,:)
  enddo
  if (debug) then
     open  (70, file='debug1.plt')
     do ip=0,n
        write (70, *) C%x(ip,:)
        write (70, *) Cout%x(ip,:)
        write (70, *)
     enddo
     close (70)
  endif


  ! check misaligned cells
  allocate (icheck(0:n-1))
  icheck = 0
  do ip=0,n-1
     L1 = C%x(ip,:)
     L2 = L1 + v(ip,:)
     M1 = C%x(ip+1,:)
     M2 = M1 + v(ip+1,:)
     if (intersect_lines(L1, L2, M1, M2, l, m)) then
        if (l >= 0.d0  .and.  l <= 1.d0  .and.  m >= 0.d0  .and.  m <= 1.d0) then
           icheck(ip) = 1
        endif
     endif
  enddo


  ! mark unwanted segments
  allocate (icheck2(0:n,-1:1))
  icheck2 = 0
  do ip=1,n-1
     if (icheck(ip-1) > 0  .or.  icheck(ip) > 0) icheck2(ip,0) = 1
  enddo

  if (debug) then
     open  (70, file='debug2.plt')
     do ip=0,n
        if (icheck2(ip,0) == 0) then
           write (70, *) C%x(ip,:)
           write (70, *) Cout%x(ip,:)
           write (70, *)
        endif
     enddo
     close (70)
  endif
!  ! test output
!  write (99, *) C%x(0,:)
!  write (99, *) C%x(0,:) + en(0,:)
!  write (99, *)
!
!  do ip=1,SRF_POLO(iz)-2
!     !if (icheck(ip-1) == 0  .and.  icheck(ip) == 0) then
!     if (icheck2(ip,0) == 0) then
!        write (99, *) C%x(ip,:)
!        write (99, *) C%x(ip,:) + en(ip,:)
!        write (99, *)
!     endif
!  enddo
!  write (99, *) C%x(SRF_POLO(iz)-1,:)
!  write (99, *) C%x(SRF_POLO(iz)-1,:) + en(SRF_POLO(iz)-1,:)
!  write (99, *)




  ! now double check each zone of misaligned cells
  ip1 = -1
  ip2 = -1
  ip  = 0
  do
     ! stop at upper boundary
     if (ip >= n+1) exit

     ! find beginning of zone with misaligned segments
     if (icheck2(ip,0) > 0) then
        ! index of last good segment
        ip1 = ip - 1
        do
           ! reached upper boundary before zone of misaligned segments endes
           if (ip >= n+1) then
              write (6, *) 'error: boundary segment is not set up correctly!'
              stop
           endif

           ! find end of zone with misaligned segments
           if (icheck2(ip,0) == 0) then
              ip2 = ip
              exit
           endif

           ! continue search
           ip = ip + 1
        enddo
        if (debug) write (6, *) 'misaligned zone between segments ', ip1, ' and ', ip2


        ! update boundaries of misaligned zone
        ! adjust upper boundary
        ip2b = ip2
        do
           if (ip2b >= n) exit

           L1 = C%x(ip1,:)
           L2 = L1 + v(ip1,:)
           M1 = C%x(ip2b,:)
           M2 = M1 + v(ip2b,:)
           if (.not.intersect_lines(L1, L2, M1, M2, l, m)) exit

           ! found new index if segments are not intersecting
           if (l < 0.d0  .or.  m < 0.d0  .or.  m > 1.d0) exit

           ip2b = ip2b + 1
        enddo

        ! adjust lower boundary
        ip1b = ip1
        do
           if (ip1b <= 0) exit

           L1 = C%x(ip1b,:)
           L2 = L1 + v(ip1b,:)
           M1 = C%x(ip2,:)
           M2 = M1 + v(ip2,:)
           if (.not.intersect_lines(L1, L2, M1, M2, l, m)) exit

           ! found new index if segments are not intersecting
           if (l < 0.d0  .or.  l > 1.d0  .or.  m < 0.d0) exit

           ip1b = ip1b - 1
        enddo
        if (debug) write (6, *) 'adjusted region: ', ip1b, ' -> ', ip2b
        ! interpolate nodes between ip1b and ip2b
        do i=ip1b+1,ip2b-1
        !   write (97, *) C%x(i,:) + en(i,:)*dl

           l = 1.d0 * (i-ip1b) / (ip2b-ip1b)
           Cout%x(i,:) = Cout%x(ip1b,:) + l * (Cout%x(ip2b,:) - Cout%x(ip1b,:))
        enddo
     endif

     ip = ip + 1
  enddo

  deallocate (icheck, icheck2)
  end function directional_extend
!===============================================================================



!===============================================================================
! scale boundary of plasma transport domain by dl (same as scale_boundary_proc)
! but find node spacing in extension to last plasma cells
!===============================================================================
subroutine ray_adjust_scale_proc(iz, ir0, idir, ir2, dl)
  use iso_fortran_env
  use emc3_grid
  use curve2D
  use run_control, only: Debug
  implicit none

  integer, intent(in)      :: iz, ir0, idir, ir2
  real(real64), intent(in) :: dl

  real(real64), dimension(:,:), allocatable :: v
  type(t_curve) :: C1, C2, scale_with_ray_adjust_nodes
  real(real64)  :: rho
  integer       :: it, ip, ir, ir1, ig, ig0


  call C1%new(ZON_POLO(iz))
  !call C2%new(ZON_POLO(iz))

!  ! 1st toroidal slice
!  it = 0
!  do ip=0,SRF_POLO(iz)-1
!     ig = ir0 + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
!     C1%x(ip,1) = RG(ig)
!     C1%x(ip,2) = ZG(ig)
!  enddo
!  call C1%left_hand_shift(-idir*dl)
!  call C1%plot(filename='debug1.plt')
!
!  ! last toroidal slice
!  it = SRF_TORO(iz)-1
!  do ip=0,SRF_POLO(iz)-1
!     ig = ir0 + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
!     C2%x(ip,1) = RG(ig)
!     C2%x(ip,2) = ZG(ig)
!  enddo
!  call C2%left_hand_shift(-idir*dl)
!  call C2%plot(filename='debug2.plt')



  allocate (v(0:SRF_POLO(iz)-1,2))
  ir1 = ir0 + idir


  do it=0,SRF_TORO(iz)-1
     ! set up boundary nodes and orientation vectors
     do ip=0,SRF_POLO(iz)-1
        ig = ir0 + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        C1%x(ip,1) = RG(ig)
        C1%x(ip,2) = ZG(ig)
        v(ip,1)   = RG(ig) - RG(ig-idir)
        v(ip,2)   = ZG(ig) - ZG(ig-idir)
     enddo
     if (Debug) call C1%plot(filename=debug_str(iz,it,1))


     C2 = scale_with_ray_adjust_nodes(C1, v, -idir*dl, 0)
     if (Debug) call C2%plot(filename=debug_str(iz,it,2))


     ! set up new grid nodes
     do ip=0,SRF_POLO(iz)-1
        ig0 = ir0 + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)
        do ir=ir1,ir2,idir
           rho = 1.d0 * (ir-ir0) / (ir2-ir1+idir)
           ig  = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)

           RG(ig) = RG(ig0) + rho * (C2%x(ip,1) - RG(ig0))
           ZG(ig) = ZG(ig0) + rho * (C2%x(ip,2) - ZG(ig0))
        enddo
     enddo
  enddo


  ! cleanup
  deallocate (v)

  return
  contains
  !.....................................................................
  function debug_str(iz, it, isuffix)
  use string
  integer, intent(in) :: iz, it, isuffix
  character(len=80)   :: debug_str

  character(len=*), parameter :: SUFFIX(2) = (/ 'raw', 'plt' /)


  if (isuffix < 0  .or.  isuffix > 2) then
     write (6, *) 'error in debug_str: 1 <= isuffix <= 2 required!'
     stop
  endif

  debug_str = 'debug/VacuumBoundary_'//trim(str(iz))//'_'//trim(str(it))//'.'//SUFFIX(isuffix)
  end function debug_str
  !.....................................................................
end subroutine ray_adjust_scale_proc
!===============================================================================



!===============================================================================
  function scale_with_ray_adjust_nodes(C, v, dl, debug_mode) result(Cout)
  use iso_fortran_env
  use curve2D
  implicit none
  type(t_curve), intent(in) :: C
  real(real64),  intent(in) :: v(0:C%n_seg,2), dl
  integer,       intent(in) :: debug_mode
  type(t_curve)             :: Cout

  real(real64), dimension(:),   allocatable :: w, w_save
  integer,      dimension(:,:), allocatable :: ip_fix
  type(t_curve) :: Ctmp
  real(real64)  :: x(2), x1(2), x2(2), sh, r
  integer       :: ip, ip1, ip2, ish, i, n, n_fix

  logical       :: debug


  ! 0. debug mode
  debug = .false.
  if (debug_mode > 0) debug = .true.


  ! 1. initialize output
  n = C%n_seg
  allocate (w(0:n), w_save(0:n))
  call Ctmp%copy(C)
  call Ctmp%left_hand_shift(dl)
  call Ctmp%setup_segment_sampling()
  if (debug) call Ctmp%plot(filename='debug_1.plt')
  call Cout%new(n)


  ! 2. find initial position on Ctmp for each node on C
  w = -1.2345d0
  if (debug) open (99, file='weight.dbg')
  do ip=0,n
     x1 = C%x(ip,:)
     x2 = x1 + v(ip,:)
     if (intersect_curve(x1, x2, Ctmp, xh=x, sh=sh, ish=ish, intersect_mode=RAY)) then
        w(ip)        = 1.d0 * (ish-1) + sh
        Cout%x(ip,:) = x
        if (debug) write (96, *) x
     else
     endif

     if (debug) write (99, *) ip, w(ip)
  enddo
  if (debug) close (99)
  if (debug) call Cout%plot(filename='debug_1a.plt')
  w_save = w


  ! 3. correct positions (w -> strictly monotonically increasing)
  ! 3.1 fix lower boundary (if required)
  ip = 0
  do
     if (w(ip) >= 0.d0  .or.  ip == n) exit
     w(ip)        = 0.d0
     Cout%x(ip,:) = Ctmp%x(0,:)

     ip = ip+1
  enddo

  ! 3.2 fix upper boundary (if required)
  ip = n
  do
     if (w(ip) >= 0.d0  .or.  ip == 0) exit
     w(ip)        = Ctmp%n_seg
     Cout%x(ip,:) = Ctmp%x(Ctmp%n_seg,:)

     ip = ip-1
  enddo

  ! 3.3 fix inner nodes
  ! 3.3.1 find misaligned zones (i.e. w not monotonically increasing)
  ! ip_fix(:,1): first node of misaligned zone
  ! ip_fix(:,2): last node of misaligned zone
  allocate(ip_fix(n,2))
  n_fix = 0
  ip    = 1
  ip1   = -1
  ip2   = -1
  do
     ! stop at upper boundary
     if (ip > n) exit

     ! find beginning of non-monotonically increasing zone
     if (w(ip) <= w(ip-1)) then
        ! index of last good node
        ip1 = ip - 1
        do
           ! reached upper boundary at a non-monotonically increasing segment
           if (ip >= n) then
              ip2 = ip
              exit
           endif

           ! find end of non-monotonically increasing zone
           if (w(ip) > w(ip-1)) then
              ip2 = ip - 1
              exit
           endif

           ! continue search
           ip = ip + 1
        enddo
        if (debug) write (6, *) 'misaligned zone between nodes ', ip1, ' and ', ip2
        n_fix           = n_fix + 1
        ip_fix(n_fix,1) = ip1
        ip_fix(n_fix,2) = ip2
     endif

     ip = ip+1
  enddo

  ! 3.3.2 adjust boundary indices for misaligned zones
  ! adjust ip_fix(:,1) and ip_fix(:,2) so that
  ! w(ip_fix(:,1)) < w (ip_fix(:,2))
  do i=1,n_fix
     ! adjust upper boundary
     ip2 = ip_fix(i,2)
     do
        if (ip2 >= n) exit
        if (w(ip2) > w(ip_fix(i,1))) exit
        ip2 = ip2+1
     enddo

     ! adjust lower boundary
     ip1 = ip_fix(i,1)
     do
        if (ip1 <= 0) exit
        if (w(ip1) < w(ip_fix(i,2))) exit
        ip1 = ip1-1
     enddo

     ip_fix(i,1) = ip1
     ip_fix(i,2) = ip2
  enddo
  ! consistency check
  do i=2,n_fix
     if (ip_fix(i,1) < ip_fix(i-1,2)) then
        write (6, *) 'error: cannot adjust boundary!'
        write (6, *) 'check weights.plt, C.plt, Ctmp.plt and Cout.plt'
        do ip=1,n_fix
           write (6, *) ip_fix(i,:)
        enddo
        open  (99, file='weights.plt')
        do ip=0,n
           write (99, *) w(ip), w_save(ip)
        enddo
        close (99)
        call C%plot(filename='C.plt')
        call Ctmp%plot(filename='Ctmp.plt')
        call Cout%plot(filename='Cout.plt')
        stop
     endif
  enddo

  ! 3.3.3 fix misaligned zones
  !write (6, *) 'misaligned zones:'
  do i=1,n_fix
     !write (6, *) ip_fix(ip,:)
     if (debug) write (6, *) 'fix zone between nodes ', ip_fix(i,1), ' and ', ip_fix(i,2)
     ip1 = ip_fix(i,1)
     ip2 = ip_fix(i,2)
     do ip=ip1+1,ip2-1
        r     = 1.d0 * (ip-ip1) / (ip2-ip1)
        w(ip) = w(ip1) + r * (w(ip2)-w(ip1))
        !Cout%x(ip,:) = Cout%x(ip1,:) + r * (Cout%x(ip2,:)-Cout%x(ip1,:))
        call Ctmp%sample_at(w(ip)/Ctmp%n_seg, x)
        Cout%x(ip,:) = x
        if (debug) write (95, *) Cout%x(ip,:)
     enddo
  enddo


!  ! 4. set up corrected boundary nodes
!  call Ctmp%setup_segment_sampling()
!  do ip=0,n
!     r = w(ip) / n
!     if (debug) write (98, *) ip, r
!     call Ctmp%sample_at(r, x)
!     Cout%x(ip,:) = x
!     if (debug) write (97, *) x
!  enddo


  ! 99. cleanup
  deallocate (w, w_save, ip_fix)

  end function scale_with_ray_adjust_nodes
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
!===============================================================================



!===============================================================================
subroutine vacuum_domain_manual_3D(iz, ir0, idir, ir2, boundary_file, filter)
  use iso_fortran_env
  use math
  use emc3_grid
  use quad_ele
  use curve2D
  use magnetic_axis, only: get_magnetic_axis
  use run_control, only: Debug
  use string
  implicit none

  integer, intent(in)      :: iz, ir0, idir, ir2
  character(len=*), intent(in) :: boundary_file, filter

  type(t_quad_ele) :: S
  type(t_curve)    :: C
  character(len=256), dimension(:), allocatable :: apply_filter, filter_parameter
  character(len=72):: tmp
  real(real64)     :: phi, A(3), theta, xi, x1(2), x2(2), rho, dl
  integer          :: ir, ir1, ip, it, ig, ifilter, is, nfilter


  ! set up filter/processing routines for boundary surfaces ............
  ! 1. count filter
  nfilter = 0
  do
     tmp = parse_string(filter, nfilter+1)
     if (tmp == '') exit
     nfilter = nfilter + 1
  enddo
  allocate (apply_filter(nfilter), filter_parameter(nfilter))
  write (6, 1000) nfilter

  ! 2. set up filter
  do ifilter=1,nfilter
     tmp = parse_string(filter, ifilter)

     is = scan(tmp, '=')
     if (is == 0) then
        apply_filter(ifilter)     = tmp
        filter_parameter(ifilter) = ''
        write (6, 1001) ifilter, trim(apply_filter(ifilter))
     else
        apply_filter(ifilter)     = tmp(1:is-1)
        filter_parameter(ifilter) = tmp(is+1:len_trim(tmp))
        write (6, 1002) ifilter, trim(apply_filter(ifilter)), trim(filter_parameter(ifilter))
     endif
  enddo
  !.....................................................................


  call S%load(boundary_file)
  ir1 = ir0 + idir
  do it=0,SRF_TORO(iz)-1
     phi = PHI_PLANE(it+PHI_PL_OS(iz))

     ! initialize slice C from boundary surface
     C   = S%slice(phi)
     if (C%n_seg < 0) then
        write (6, *) 'error in S%slice'
        stop
     endif
     if (Debug) then
        write (tmp, 9001) iz, it
        call C%plot(filename=tmp)
     endif


     ! apply filter for slice C
     do ifilter=1,nfilter
        select case(apply_filter(ifilter))
        case('expand')
           read  (filter_parameter(ifilter), *) dl
           call C%left_hand_shift(dl)

        case default
           write (6, *) 'error: invalid filter type ', apply_filter(1:is-1)
           stop
        end select
        if (Debug) then
           write (tmp, 9002) iz, it, ifilter
           call C%plot(filename=tmp)
        endif
     enddo


     ! setup sampling on slice C
     A   = get_magnetic_axis(phi)
     call C%sort_loop(A)
     call C%setup_angular_sampling(A(1:2))
     if (Debug) then
        write (tmp, 9003) iz, it
        call C%plot(filename=tmp)
     endif


     ! poloidal loop (setup segment weights)
     do ip=0,SRF_POLO(iz)-1
        ig = ir0 + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)

        ! x1: reference point on last EMC3 surface -> theta
        x1(1) = RG(ig);  x1(2) = ZG(ig)
        theta = atan2(x1(2) - A(2), x1(1) - A(1));  if (theta < 0.d0) theta = theta + pi2
        xi    = theta / pi2

        ! x2: reference point on boundary (theta)
        call C%sample_at(xi, x2)

        ! generate vacuum domain
        do ir=ir1,ir2,idir
           rho = 1.d0 * (ir-ir0) / (ir2-ir1+idir)
           ig  = ir + (ip + it*SRF_POLO(iz))*SRF_RADI(iz) + GRID_P_OS(iz)

           RG(ig) = x1(1) + rho * (x2(1) - x1(1))
           ZG(ig) = x1(2) + rho * (x2(2) - x1(2))
        enddo
     enddo
  enddo


  ! cleanup
  deallocate (apply_filter, filter_parameter)

 1000 format(8x,'using ',i0,' filter to process boundary surface')
 1001 format(8x,i0,': ',a)
 1002 format(8x,i0,': ',a,' with parameter ',a)
 9001 format('debug/slice_Z',i0,'_T',i0,'.plt')
 9002 format('debug/slice_Z',i0,'_T',i0,'_filter',i0,'.plt')
 9003 format('debug/slice_Z',i0,'_T',i0,'_sample.plt')
end subroutine vacuum_domain_manual_3D
!===============================================================================
