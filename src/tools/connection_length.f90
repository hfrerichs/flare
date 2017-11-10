!===============================================================================
! Calculate connection length
!
! Input (taken from run control file):
!    Grid_File          Provide initial positions for several field lines
!
!    Trace_Step         Size of trace steps (see Trace_Coords)
!    Limit              Maximum distance for field line tracing (in one direction)
!    max_pt             Maximum number of poloidal turns
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
!                     = 256: Final radial position (PsiN) in backward direction
!                     = 512: Final radial position (PsiN) in forward direction
!                    = 1024: Backward connection length [toroidal turns]
!                    = 2048: Forward connection length [toroidal turns]
!    Psi(1), Psi(2)    Reference radial coordinate for additional data type 8-32
!    Output_File
!===============================================================================
subroutine connection_length
  use run_control, only: Grid_File, Output_File, Trace_Step, Trace_Method, Trace_Coords, &
                         Output_Format, Limit, Psi, max_pt
  use grid
  use dataset
  use parallel
  use math
  use equilibrium
  use boundary
  use fieldline
  use flux_surface_3D
  implicit none

  integer, parameter :: nout_max = 12
  integer, parameter :: iu = 42

  type(t_flux_surface_3D)                   :: LCFS
  type(t_fieldline)  :: F
  type(t_grid)       :: G
  type(t_dataset)    :: D

  real(real64)       :: y(3), r(3), dl, PsiN, Psi_min, Psi_av, ntt(-1:1)
  real(real64)       :: lc(-1:1), npt(-1:1), dist2PsiN(-1:1,2), dist, d_min, PsiN_final(-1:1)
  logical :: distance_to_lcfs
  integer :: itrace, nout, iout(nout_max), i, i2, ig, ig1, iflag, idir, id, id_limit(-1:1)


  if (firstP) then
     write (6, *) 'Calculate connection length, output in: ', adjustl(trim(Output_File))
     write (6, *)
  endif


! load grid for connection length calculation
  itrace = 1
  if (Trace_Coords > 1) itrace = 2
  call G%load(Grid_File)


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


  call D%new(G%nodes(), 5+nout)
  call D%set_info(2, trim(Grid_File))
  call D%set_column_info(1, 'Lc_bwd', "Backward connection length [m]")
  call D%set_column_info(2, 'Lc_fwd', "Forward connection length [m]")
  call D%set_column_info(3, 'Lpt_bwd', "Backward connection length [poloidal turns]")
  call D%set_column_info(4, 'Lpt_fwd', "Forward connection length [poloidal turns]")
  call D%set_column_info(5, 'minPsiN', "Minimum(PsiN) along field line")
  if (firstP) call additional_output_info()


  ! 2.7. distance to last closed flux surface
  if (distance_to_lcfs) then
     call LCFS%load_distance_to()
  endif
! ... prepare output data arrays (end) .................................


! main loop (begin) ....................................................
  grid_loop: do ig=mype,G%nodes()-1,nprs
     ig1 = ig + 1
     y = G%node(ig1, coordinates=itrace)

     ! initial values
     ! ... for connection length calculation
     lc        = 0.d0
     npt       = 0.d0
     ntt       = 0.d0

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
           dl = F%trace_1step()

           ! update connection length
           lc(idir)   = lc(idir) + dl
           if (abs(lc(idir)) .ge. Limit) exit trace_loop
           if (max_pt > 0.d0  .and.  abs(F%theta_int) .ge. max_pt*pi2) exit trace_loop


           ! update field line penetration
           PsiN       = F%get_PsiN()
           if (PsiN.lt.Psi_min) Psi_min = PsiN


           ! update average pol. flux
           Psi_av = Psi_av + PsiN * abs(dl)


           ! distance (along field line) to PsiN(1:2)
           if (dist2PsiN(idir,1).le.0.d0 .and. F%cross_PsiN(Psi(1))) dist2PsiN(idir,1) = abs(lc(idir))
           if (dist2PsiN(idir,2).le.0.d0 .and. F%cross_PsiN(Psi(2))) dist2PsiN(idir,2) = abs(lc(idir))


           ! shortest distance (in RZ-plane) to last closed flux surface
           if (distance_to_lcfs) then
              dist = abs(LCFS%get_distance_to(F%rc))
              if (dist < d_min) d_min = dist
           endif


           ! check intersection with walls
           if (F%intersect_boundary(id=id)) then
              id_limit(idir) = id
              exit trace_loop
           endif
        enddo trace_loop
        npt(idir)        = F%theta_int / pi2
        ntt(idir)        = F%phi_int   / pi2
        PsiN_final(idir) = PsiN
     enddo


     ! update connection length data for selected field line
     lc(0)  = abs(lc(-1)) + abs(lc(1))
     Psi_av = Psi_av / lc(0)
     D%x(ig1, 1) = abs(lc(-1) / 100.d0)
     D%x(ig1, 2) = abs(lc( 1) / 100.d0)
     D%x(ig1, 3) = npt(-1)
     D%x(ig1, 4) = npt( 1)
     D%x(ig1, 5) = Psi_min
     do i=1,nout
        i2 = nint(log(1.d0*iout(i))/log(2.d0))
        select case (i2)
        case (0)
           D%x(ig1,5+i) = Psi_av
        case (1)
           D%x(ig1,5+i) = real(id_limit(-1))
        case (2)
           D%x(ig1,5+i) = real(id_limit( 1))
        case (3)
           D%x(ig1,5+i) = dist2PsiN(-1,1)
        case (4)
           D%x(ig1,5+i) = dist2PsiN( 1,1)
        case (5)
           D%x(ig1,5+i) = dist2PsiN(-1,2)
        case (6)
           D%x(ig1,5+i) = dist2PsiN( 1,2)
        case (7)
           D%x(ig1,5+i) = d_min
        case (8)
           D%x(ig1,5+i) = PsiN_final(-1)
        case (9)
           D%x(ig1,5+i) = PsiN_final(1)
        case (10)
           D%x(ig1,5+i) = ntt(-1)
        case (11)
           D%x(ig1,5+i) = ntt(1)
        end select
     enddo

     write (6, 4000) ig, lc(0)/1.d2
  enddo grid_loop
! ... main loop (end) ..................................................


! finalize .............................................................
  call D%mpi_allreduce()

  if (firstP) then
     call D%store(filename=Output_File)
  endif
! ......................................................................


  return
 4000 format (5x,i8,',',8x,'L_c = ',f9.2,' m')
 5010 write (6,5011) Output_File
 5011 format ('error opening file for output: ', a120)
  stop
  contains
!.......................................................................
  subroutine additional_output_info

  character(len=256) :: qkey, text
  integer :: i, i2


  if (nout.gt.0) then
     write (6, 1000)
  endif

  do i=1,nout
     i2 = nint(log(1.d0*iout(i))/log(2.d0))
     select case (i2)
     case (0)
        qkey = 'avPsiN'
        text = 'Field line averaged poloidal flux'
     case (1)
        qkey = 'surf_id_bwd'
        text = 'Limiting surface id in neg. direction'
     case (2)
        qkey = 'surf_id_fwd'
        text = 'Limiting surface id in pos. direction'
     case (3)
        qkey = 'dist_PsiN1_bwd'
        write (text, '(f8.4)') Psi(1)
        text = 'Backward distance to Psi = '//trim(text)
     case (4)
        qkey = 'dist_PsiN1_fwd'
        write (text, '(f8.4)') Psi(1)
        text = 'Forward  distance to Psi = '//trim(text)
     case (5)
        qkey = 'dist_PsiN2_bwd'
        write (text, '(f8.4)') Psi(2)
        text = 'Backward distance to Psi = '//trim(text)
     case (6)
        qkey = 'dist_PsiN2_fwd'
        write (text, '(f8.4)') Psi(2)
        text = 'Forward  distance to Psi = '//trim(text)
     case (7)
        qkey = 'dist_LCFS'
        text = 'Distance to last closed flux surface'
     case (8)
        qkey = 'PsiN_bwd'
        text = 'Final radial position (PsiN) in backward direction'
     case (9)
        qkey = 'PsiN_fwd'
        text = 'Final radial position (PsiN) in forward direction'
     case (10)
        qkey = 'Ltt_bwd'
        text = 'Backward connection length [toroidal turns]'
     case (11)
        qkey = 'Ltt_fwd'
        text = 'Forward connection length [toroidal turns]'
     case default
        write (6, *) 'error: ', 2**i2, ' is not a valid data id!'
        stop
     end select
     call D%set_column_info(5+i, trim(qkey), trim(text))
     write (6, 2001) trim(text)
  enddo

 1000 format (3x,'- Additional output:')
 2001 format (8x,a)
  end subroutine additional_output_info
!.......................................................................
end subroutine connection_length
