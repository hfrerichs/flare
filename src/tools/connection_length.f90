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
!    Output_Format      = 0: default operation mode, i.e. backward and forward
!                            connection length in cm and pol. turns and minimum
!                            pol. flux along field line
!                       > 0: provide additional information, the format value is
!                            binary coded, i.e. 6 = 2 + 4 will provide both
!                            additional data of type 2 and 4.
!                       = 1: Field line averaged poloidal flux
!
!                            id = 0 if connection length > Limit
!                       = 2: Limiting surface id in neg. direction
!                       = 4: Limiting surface id in pos. direction
!
!                            8-32 are set automatically if Psi(1) or Psi(2) > 0
!                       = 8: Backward distance along field line to Psi(1)
!                      = 16: Forward  distance along field line to Psi(1)
!                      = 32: Backward distance along field line to Psi(2)
!                      = 64: Forward  distance along field line to Psi(2)
!                     = 128: Radial distance to last closed flux surface (lcfs.dat)
!    Psi(1), Psi(2)    Reference radial coordinate for additional data type 8-32
!    Output_File
!===============================================================================
subroutine connection_length
  use run_control, only: Grid_File, Output_File, Trace_Step, Trace_Method, Trace_Coords, &
                         Output_Format, Limit, Psi
  use usr_grid
  use parallel
  use math
  use equilibrium
  use boundary
  use fieldline
  use flux_surface_3D
  implicit none

  integer, parameter :: nout_max = 10
  integer, parameter :: iu = 42

  real(real64), dimension(:,:), allocatable :: lc_data
  type(t_flux_surface_3D)                   :: LCFS
  type(t_fieldline)  :: F

  character(len=12)  :: fstr
  real(real64)       :: y(3), r(3), PsiN, Psi_min, Psi_av
  real(real64)       :: lc(-1:1), lpt(-1:1), dist2PsiN(-1:1,2), d, d_min
  logical :: distance_to_lcfs
  integer :: itrace, nout, iout(nout_max), i, i2, ig, iflag, idir, id, id_limit(-1:1)


  if (firstP) then
     write (6, *) 'Calculate connection length, output in: ', adjustl(trim(Output_File))
     write (6, *)
  endif


! load grid for connection length calculation
  itrace = 1
  if (Trace_Coords > 1) itrace = 2
  call read_grid (Grid_File, log_progress=.false., use_coordinates=COORDINATES(itrace))


! prepare output data arrays (begin) ..................................
  ! 1. distance along field line to Psi(1) and Psi(2)
  if (Psi(1) > 0.d0) then
     Output_Format = Output_Format - iand( 8,Output_Format) +  8
     Output_Format = Output_Format - iand(16,Output_Format) + 16
  endif
  if (Psi(2) > 0.d0) then
     Output_Format = Output_Format - iand(32,Output_Format) + 32
     Output_Format = Output_Format - iand(64,Output_Format) + 64
  endif

  ! 2. other additional output
  nout = 0
  iout = 0
  distance_to_lcfs = .false.
  do i=0,nout_max-1
     i2 = 2**i
     if (i2.eq.iand(i2, Output_Format)) then
        nout = nout + 1
        iout(nout) = i2

        ! distance to last closed flux surface
        if (i == 7) distance_to_lcfs = .true.
     endif
  enddo
  if (firstP) call additional_output_info()

  ! set format string
  write (fstr, '(i4)') 5+nout
  fstr = '('//trim(adjustl(fstr))//'e18.10)'

  allocate (lc_data(0:n_grid-1,5+nout))
  lc_data = 0.d0


  ! 2.7. distance to last closed flux surface
  if (distance_to_lcfs) then
     call LCFS%load('lcfs.dat')
     call LCFS%sample_distance_to(grid='distance.grid')
  endif
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

     ! ... for average and maximum penetration of field lines
     call coord_trans (y, Trace_Coords, r, CYLINDRICAL)
     Psi_min   = get_PsiN(r)
     Psi_av    =  0.d0

     ! ... for additional data (distance to flux surface, limiting surface id, ...)
     dist2PsiN = 0.d0
     id_limit  = 0

     ! ... for distance to last closed flux surface
     d_min = huge(d_min)

     ! trace each field line in positive and negative direction
     do idir=-1,1,2
        call F%init(y, idir*Trace_Step, Trace_Method, Trace_Coords)

        ! start field line tracing
        trace_loop: do
           call F%trace_1step()

           ! update connection length
           lc(idir)   = lc(idir) + Trace_Step
           if (abs(lc(idir)) .ge. Limit) exit trace_loop


           ! update field line penetration
           PsiN       = F%get_PsiN()
           if (PsiN.lt.Psi_min) Psi_min = PsiN


           ! update average pol. flux
           Psi_av = Psi_av + PsiN * abs(Trace_Step)


           ! distance (along field line) to PsiN(1:2)
           if (dist2PsiN(idir,1).le.0.d0 .and. F%cross_PsiN(Psi(1))) dist2PsiN(idir,1) = abs(lc(idir))
           if (dist2PsiN(idir,2).le.0.d0 .and. F%cross_PsiN(Psi(2))) dist2PsiN(idir,2) = abs(lc(idir))


           ! check intersection with walls
           if (F%intersect_boundary(id=id)) then
              id_limit(idir) = id
              exit trace_loop
           endif


           ! shortest distance (in RZ-plane) to last closed flux surface
           if (distance_to_lcfs) then
              d = abs(LCFS%get_distance_to(F%rc))
              if (d < d_min) d_min = d
           endif
        enddo trace_loop
        lpt(idir)  = F%theta_int
     enddo


     ! update connection length data for selected field line
     lc(0)  = abs(lc(-1)) + abs(lc(1))
     lpt    = lpt / pi2
     Psi_av = Psi_av / lc(0)
     lc_data(ig,1) = lc(-1)
     lc_data(ig,2) = lc( 1)
     lc_data(ig,3) = lpt(-1)
     lc_data(ig,4) = lpt( 1)
     lc_data(ig,5) = Psi_min
     do i=1,nout
        i2 = nint(log(1.d0*iout(i))/log(2.d0))
        select case (i2)
        case (0)
           lc_data(ig,5+i) = Psi_av
        case (1)
           lc_data(ig,5+i) = real(id_limit(-1))
        case (2)
           lc_data(ig,5+i) = real(id_limit( 1))
        case (3)
           lc_data(ig,5+i) = dist2PsiN(-1,1)
        case (4)
           lc_data(ig,5+i) = dist2PsiN( 1,1)
        case (5)
           lc_data(ig,5+i) = dist2PsiN(-1,2)
        case (6)
           lc_data(ig,5+i) = dist2PsiN( 1,2)
        case (7)
           lc_data(ig,5+i) = d_min
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
     i2 = nint(log(1.d0*iout(i))/log(2.d0))
     select case (i2)
     case (0)
        text = 'Field line averaged poloidal flux'
     case (1)
        text = 'Limiting surface id in neg. direction'
     case (2)
        text = 'Limiting surface id in pos. direction'
     case (3)
        write (text, '(f8.4)') Psi(1)
        text = 'Backward distance to Psi = '//trim(text)
     case (4)
        write (text, '(f8.4)') Psi(1)
        text = 'Forward  distance to Psi = '//trim(text)
     case (5)
        write (text, '(f8.4)') Psi(2)
        text = 'Backward distance to Psi = '//trim(text)
     case (6)
        write (text, '(f8.4)') Psi(2)
        text = 'Forward  distance to Psi = '//trim(text)
     case (7)
        text = 'Distance to last closed flux surface'
     case default
        write (6, *) 'error: ', 2**i2, ' is not a valid data id!'
        stop
     end select
     write (iu, 2000) i2, text
     write (6, 2001) trim(text)
  enddo
  if (nout.gt.0) close (iu)

 1000 format (3x,' - Additional output:')
 2000 format (i4,3x,a)
 2001 format (8x,a)
  end subroutine additional_output_info
!.......................................................................
end subroutine connection_length
