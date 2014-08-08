!===============================================================================
! Calculate connection length
!
! Input (taken from run control file):
!    Grid_File          Provide initial positions for several field lines
!
!    Trace_Step         Size of trace steps (see Trace_Coords)
!    Limit              Maximum distance for field line tracing (in one direction)
!
!    Trace_Method	> 0: Integration (see module ODE_solver)
!                       = 0: Reconstruction from field aligned grid
!
!    Trace_Coords	= 1: solve d(x,y,z)/dl   = (Bx,By,Bx)/B
!                            x_start: Cartesian coordinates [cm]
!                            Trace_step: straight segment [cm]
!
!                       = 2: solve d(R,Z,phi)/ds = (BR,BZ,Bphi/R)/B
!                            x_start: Cylindrical coordinates [cm,deg]
!                            Trace_step: arc segment [cm]
!
!                       = 3: solve d(R,Z)/dphi   = (BR,BZ)*R/Bphi
!                            x_start: Cylindrical coordinates [cm,deg]
!                            Trace_step: arc segment [deg]
!
!    Output_File
!    Output_Format
!===============================================================================
subroutine connection_length
  use run_control, only: Grid_File, Output_File, Trace_Step, Trace_Method, Trace_Coords, &
                         Output_Format
  use grid
  use parallel
  use math
  use equilibrium
  use boundary
  use fieldline
  implicit none

  integer, parameter :: nout_max = 10
  integer, parameter :: iu = 42

  real*8, dimension(:,:), allocatable :: lc_data
  type(t_fieldline)  :: F

  character*12 :: fstr
  real*8  :: y(3), yc(3), yl(3), rc(3), rl(3), X(3), thetac, thetal, dtheta, maxis(3)
  real*8  :: Psi, Psi_min, Psi_av
  real*8  :: lc(-1:1), lpt(-1:1), dist_Psi0(-1:1)
  integer :: itrace, nout, iout(nout_max), i, i2, ig, iflag, idir
  real*8  :: xc(3), xl(3)


  if (firstP) then
     write (6, *) 'Calculate connection length, output in: ', adjustl(trim(Output_File))
     write (6, *)
  endif


! load grid for connection length calculation
  itrace = 1
  if (Trace_Coords > 1) itrace = 2
  call read_grid (Grid_File, log_progress=.false., use_coordinates=COORDINATES(itrace))


! prepare output data arrays (begin) ..................................
  nout = 0
  iout = 0
  do i=0,nout_max-1
     i2 = 2**i
     if (i2.eq.iand(i2, Output_Format)) then
        nout = nout + 1
        iout(nout) = i2
     endif
  enddo
  if (firstP) call additional_output_info()

  ! set format string
  write (fstr, '(i4)') 5+nout
  fstr = '('//trim(adjustl(fstr))//'e18.10)'

  allocate (lc_data(0:n_grid-1,5+nout))
  lc_data = 0.d0
! ... prepare output data arrays (end) .................................


! main loop (begin) ....................................................
  ig = mype
  grid_loop: do
     call get_next_grid_point (iflag, y)
     if (iflag.eq.-1) exit grid_loop

     ! initial values
     ! ... for connection length calculation
     lc        = 0.d0
     lpt       = 0.d0
     dist_Psi0 = 0.d0

     ! ... for average and maximum penetration of field lines
     call coord_trans (y, Trace_Coords, rc, CYLINDRICAL)
     Psi       = get_Psi(rc)
     Psi_min   = (Psi-Psi_axis)/(Psi_sepx-Psi_axis)
     Psi_av    =  0.d0

     ! trace each field line in positive and negative direction
     do idir=-1,1,2
        Trace_Step = - Trace_Step
        yl         = y
        call coord_trans (y, Trace_Coords, rl, CYLINDRICAL)
        call coord_trans (y, Trace_Coords, xl, CARTESIAN)
        maxis      = magnetic_axis(rl(3))
        thetal     = datan2(rl(2) - maxis(2), rl(1) - maxis(1))
        call F%init(y, Trace_Step, Trace_Method, Trace_Coords)

        ! start field line tracing
        trace_loop: do
           yc = F%next_step()
           call coord_trans (yc, Trace_Coords, rc, CYLINDRICAL)
           call coord_trans (yc, Trace_Coords, xc, CARTESIAN)
           !write (99, '(6e14.6)') yc, rc
           maxis      = magnetic_axis(rc(3))
           lc(idir)   = lc(idir) + Trace_Step
           if (abs(lc(idir)) .ge. Limit) exit trace_loop

           ! integrate poloidal angle
           thetac     = datan2(rc(2) - maxis(2), rc(1) - maxis(1))
           dtheta     = thetac - thetal
           if (abs(dtheta).gt.pi) dtheta = dtheta - dsign(pi2,dtheta)
           lpt(idir)  = lpt(idir)  + dtheta

           ! update field line penetration
           Psi        = get_Psi(rc)
           Psi        = (Psi-Psi_axis)/(Psi_sepx-Psi_axis)
           if (Psi.lt.Psi_min) then
              Psi_min = Psi
              !phi_min = phic
              !theta_min = thetac
           endif

           ! update average pol. flux
           Psi_av = Psi_av + Psi * abs(Trace_Step)
!               if (lRD2LCFS) then
!                  call calc_RD2LCFS(rc, xc(3), phic, d)
!                  if (d.lt.d_min) d_min = d
!               endif

           ! distance to Psi =1
           if (dist_Psi0(idir).le.0.d0 .and. Psi.le.1.d0) dist_Psi0(idir) = abs(lc(idir))


           ! check intersection with walls
           if (intersect_boundary(rl, rc, X)) then
                  !x_hit(idir,:) = xh
                  !s_hit(idir)   = s
                  !is_hit(idir)  = is
              exit trace_loop
           endif

           ! prepare trace next step
           yl     = yc
           xl     = xc
           rl     = rc
           thetal = thetac
        enddo trace_loop
     enddo

     ! update connection length data for selected field line
     lc(0)  = abs(lc(-1)) + abs(lc(1))
     lpt    = lpt / pi2

     Psi_av = Psi_av / lc(0)
     ! ...
     dist_Psi0(0) = dist_Psi0(-1)
     ! ....

     lc_data(ig,1) = lc(-1)
     lc_data(ig,2) = lc( 1)
     lc_data(ig,3) = lpt(-1)
     lc_data(ig,4) = lpt( 1)
     lc_data(ig,5) = Psi_min
     do i=1,nout
        select case (iout(i))
        case (1)
           lc_data(ig,5+i) = Psi_av
!        case (2)
!           lc_data(i_grid,5+i) = 1.d0 * is_hit(-1)
!        case (4)
!           lc_data(i_grid,5+i) = 1.d0 * is_hit(1)
!        case (8)
!           lc_data(ig,5+i) = phi_min
!        case (16)
!           lc_data(ig,5+i) = theta_min
        case (32)
           lc_data(ig,5+i) = dist_Psi0(0)
        end select
     enddo

     write (6, 4000) ig, lc(0)/1.d2
     ig = ig + nprs
  enddo grid_loop
! ... main loop (end) ..................................................


! finalize .............................................................
  call wait_pe()
  call sum_real_data (lc_data, (5+nout)*n_grid)

  if (firstP) then
     open  (iu, file=Output_File, err=5010)
     do ig=0,n_grid-1
        write (iu, fstr) lc_data(ig,:)
     enddo
     close (iu)
  endif

  deallocate (lc_data)
! ......................................................................


  return
 4000 format (5x,i8,',',8x,'L_c = ',f9.2,' m')
 5010 write (6,5011) Output_File
 5011 format ('error opening file for output: ', a120)
  stop
  contains
!.......................................................................
  subroutine additional_output_info

  character*120 :: info_file, text
  integer :: i, i2


  if (nout.gt.0) then
     write (6, 1000)
     info_file = trim(adjustl(Output_File))//'.info'
     open  (iu, file=info_file)
     write (iu, *) nout
  endif

  do i=1,nout
     i2 = iout(i)
     select case (iout(i))
     case (1)
        text = 'Field line averaged poloidal flux'
     case (2)
        text = 'Limiting surface # (neg. direction)'
     case (4)
        text = 'Limiting surface # (pos. direction)'
     case (8)
        text = 'Toroidal location of deepest penetration'
     case (16)
        text = 'Poloidal location of deepest penetration'
     case (32)
        text = 'Distance to Psi=1'
     end select
     write (iu, 2000) i2, text
  enddo
  if (nout.gt.0) close (iu)

 1000 format (3x,' - Additional output:')
 2000 format (i4,3x,a)
  end subroutine additional_output_info
end subroutine connection_length
