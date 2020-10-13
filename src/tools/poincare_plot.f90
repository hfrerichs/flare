!===============================================================================
! Generate Poincare plot
!
! Input (taken from run control file):
! start locations: 1.
!    R_start, R_end     Equidistant steps between R_start and R_end, Phi = 0
!    Z_start, Z_end     Equidistant steps between Z_start and Z_end (default Z = 0)
!    N_steps		Number of steps between start locations for field line tracing
! or. 2.
!    Grid_File          Provide file with start locations
! or. 3.
!    x_start            Single start location for field line tracing
!
!
!    N_points           Max. number of points for each Poincare plot
!    N_sym              Apply toroidal symmetry, i.e. Poincare sections at multiples of 2*pi/N_sym
!    N_mult		Generate multiple Poincare plots: the toroidal domain 2*pi/N_sym is split into N_mult sub-domains
!    Phi_Output         Toroidal position [deg] of (base) Poincare plot
!
!
!    Trace_Step         Size of trace steps (see Trace_Coords)
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
subroutine poincare_plot
  use run_control, only: R_start, R_end, N_steps, Grid_File, Output_File, Output_Format, &
                         Z_start, Z_end, &
                         Trace_Step, Trace_Method, Trace_Coords, &
                         N_points, N_sym, N_mult, Phi_output, x_start, stop_at_boundary
  use parallel
  use equilibrium
  use boundary
  use fieldline
  use grid
  use math
  implicit none

  integer, parameter :: iu = 40

  type t_poincare_data 
     ! number of sets (i.e. start points) and Poincare cuts (= N_mult)
     integer :: n_sets, n_cuts

     ! number of points for each set and each Poincare cut
     integer, dimension(:,:), allocatable :: n_points

     ! data points
     real*8, dimension(:,:,:,:), allocatable :: X
  end type t_poincare_data

  type(t_grid)          :: G
  type(t_fieldline)     :: F
  type(t_poincare_data) :: Pdata
  

  character*120 :: tmp_filename, suffix, smult
  real*8  :: lc, y(3), yc(3), X(3)
  real*8  :: dr, dz, Psi, theta
  integer :: imult, j, ig, icut, n_grid


  if (firstP) then
     write (6, *) 'Generate Poincare plot, output in: ', adjustl(trim(Output_File))
     write (6, 1001) N_sym
     if (N_mult > 1) write (6, 1002) N_mult
     write (6, 1003) Phi_output / pi * 180.d0
     write (6, 1004) COORDINATES(Trace_Coords)
     write (6, *)
  endif


! open output file(s) ..................................................
  if (firstP) then
  if (N_mult == 1) then
     open  (iu, file=Output_File)
     call write_header(iu)
  else
     do imult=0,N_mult-1
        tmp_filename = Output_File
        ! find and store suffix
        suffix = ''
        do j=len(tmp_filename),1,-1
           if (tmp_filename(j:j) == '.') then
              suffix = tmp_filename(j:len(tmp_filename))
              tmp_filename = tmp_filename(1:j-1)
              exit
           endif
        enddo

        ! insert multiplicity number
        write (smult, *) imult
        tmp_filename = trim(tmp_filename)//'_'//trim(adjustl(smult))//trim(suffix)
        open  (iu+imult, file=tmp_filename)
        call write_header(iu+imult)
     enddo
  endif
  endif
!.......................................................................


! prepare grid/radial range for poincare plot ..........................
  ! use radial range R_start -> R_end
  if (R_start > 0.d0) then
     dr   = 0.d0
     dz   = 0.d0
     call G%new(CYLINDRICAL, UNSTRUCTURED, FIXED_COORD3, N_steps+1)

     if (N_steps .ne. 0) then
        dr = (R_end - R_start) / N_steps
        dz = (Z_end - Z_start) / N_steps
        if (firstP) write (6,2001) R_start, R_end, N_steps
        do ig=0,N_steps
           G%x(ig+1,1) = R_start + ig*dr; G%x(ig+1,2) = Z_start + ig*dz
        enddo
     else
        if (firstP) write (6,2002) R_start
        G%x(1,1) = R_start; G%x(1,2) = Z_start
     endif
     G%x(:,3) = x_start(3)

  ! use x_start if N_steps is not specified
  elseif (N_steps == 0  .and.  x_start(1).ne.0.d0) then
     call G%new(CYLINDRICAL, UNSTRUCTURED, FIXED_COORD3, 1)
     G%x(1,:) = x_start
     if (firstP) write (6, 2003) x_start

  ! use points from grid file if R_start is missing
  elseif (Grid_File .ne. '') then
     call G%load(Grid_File)
     if (firstP) write (6,*)

  else
     if (firstP) write (6, *) 'error: starting points undefined!'
     stop
  endif
!.......................................................................


!.......................................................................
! prepare output data array
  n_grid       = G%nodes()
  Pdata%n_sets = n_grid
  Pdata%n_cuts = N_mult
  allocate (Pdata%n_points(0:n_grid-1, 0:N_mult-1))
  Pdata%n_points = 0
  allocate (Pdata%X       (0:n_grid-1, 0:N_mult-1, N_points, 4))
  Pdata%X        = 0.d0
!.......................................................................


! main loop (begin) ....................................................
  grid_loop: do ig=mype,n_grid-1,nprs
     yc = G%node(ig+1, coordinates=min(Trace_Coords,2))

     ! integrate connection length
     lc = 0.d0


     call F%init_toroidal_tracing(yc, Trace_Step, Trace_Method, Trace_Coords, N_sym*N_mult, Phi_output)
     icut = 0

     ! start field line tracing
     trace_loop: do
        lc = lc + F%trace_1step()
        if (F%ierr > 0) then
           write (6, 4002) ig
           exit trace_loop
        endif


        ! check intersection with Poincare plane
        if (F%intersect_sym_plane(icut, X, theta)) then
           imult = int(mod(icut,N_mult))
           if (imult < 0) imult = imult + N_mult
           j     = Pdata%n_points(ig,imult) + 1
           Pdata%n_points(ig,imult) = j

           Pdata%X(ig,imult,j,1:2)  = X(1:2)
           select case(Output_Format)
           case(1)
              theta           = get_poloidal_angle(X)
              if (theta.le.0.d0) theta = theta + pi2
           case(2)
              ! nothing to be done here, theta is already set in F%intersect_cym_plane
           case default
              write (6, *) 'error: invalid output format ', Output_Format
              stop
           end select
           Pdata%X(ig,imult,j,3  ) = theta / pi2 * 360.d0
           Psi             = get_Psi(X)
           Psi             = (Psi-Psi_axis)/(Psi_sepx-Psi_axis)
           Pdata%X(ig,imult,j,4  ) = Psi
        endif



        ! check intersection with boundaries
        if (stop_at_boundary  .and.  F%intersect_boundary(X)) then
           write (6,4000) ig, abs(lc/1.d2), abs(icut)/N_mult
           exit trace_loop
        endif


        ! check upper limit for field line tracing
        if (abs(icut) .ge. N_points*N_mult) then
           write (6,4001) ig, abs(lc/1.d2), abs(icut)/N_mult
           exit trace_loop
        endif
     enddo trace_loop
  enddo grid_loop
! ... main loop (end) ..................................................


! finalize .............................................................
  if (nprs > 1) then
     call wait_pe()
     call sum_inte_data (Pdata%n_points, n_grid*N_mult)
     call sum_real_data (Pdata%X       , n_grid*N_mult*N_points*4)
  endif

  if (firstP) then
  ! write data
  do ig=0,n_grid-1
     do imult=0,N_mult-1
        do j=1,Pdata%n_points(ig,imult)
           write (iu+imult, '(4(e25.18,1x))') Pdata%X(ig,imult,j,:)
        enddo
     enddo
  enddo

  ! close output file(s)
  do imult=0,N_mult-1
     close (iu+imult)
  enddo
  endif

  deallocate (Pdata%n_points, Pdata%X)
! ......................................................................


 1001 format (8x,'Toroidal symmetry number:         ',i4)
 1002 format (8x,'Number of slices to be generated: ',i4)
 1003 format (8x,'Reference location     ',f8.3,' deg')
 1004 format (8x,'Field line tracing coordinates: ',a)
 2001 format (8x,'radial domain:',5x,'R_start = ',f6.2,5x, &
              'R_end = ',f6.2,5x,'with ',i4,' steps'/)
 2002 format (8x,'position:',5x,'R = ',f6.2,5x/)
 2003 format (8x,'Initial coordinate: (',f7.3,', ',f7.3,', ',f7.3,')')

 4000 format (5x,i5,',',8x,'L_c = ',f9.2,' m,',8x,'n_points = ',i5)
 4001 format (5x,i5,',',8x,'L_c > ',f9.2,' m,',8x,'n_points = ',i5)
 4002 format (5x,i5,',',8x,'Field line leaves magnetic field domain')
  contains
  subroutine write_header(iu)
  integer, intent(in) :: iu

  write (iu, 5001)
  write (iu, 5002) 'R "Major radius [cm]"'
  write (iu, 5002) 'Z "Vertical coordinate [cm]"'
  write (iu, 5002) 'Theta "Poloidal angle [deg]"'
  write (iu, 5002) 'PsiN "Normalized poloidal flux"'
  write (iu, 5003)
  write (iu, 5004) "Theta-PsiN"
  write (iu, 5004) "R-Z"
 5001 format("# DATA DIMENSION 2D")
 5002 format("# DATA COLUMN ",a)
 5003 format("# TYPE POINT_DATA")
 5004 format("# GEOMETRY ",a)
  end subroutine write_header
end subroutine poincare_plot
