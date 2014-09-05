!===============================================================================
! Generate Poincare plot
!
! Input (taken from run control file):
! start locations: 1.
!    R_start, R_end     Equidistant steps between R_start and R_end at Z = 0, Phi = 0
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
                         Trace_Step, Trace_Method, Trace_Coords, &
                         N_points, N_sym, N_mult, Phi_output, x_start
  use parallel
  use equilibrium
  use boundary
  use fieldline
  use grid
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

  type(t_fieldline)     :: F
  type(t_poincare_data) :: Pdata
  real*8, dimension(:,:), pointer     :: my_grid
  

  character*120 :: tmp_filename, suffix, smult
  real*8  :: lc, y(3), yc(3), X(3)
  real*8  :: dr, Psi, theta
  integer :: imult, j, ig, iflag, icut


  if (firstP) then
     write (6, *) 'Generate Poincare plot, output in: ', adjustl(trim(Output_File))
     write (6, 1001) N_sym
     if (N_mult > 1) write (6, 1002) N_mult
     write (6, 1003) Phi_output
     write (6, *)
  endif


! open output file(s) ..................................................
  if (firstP) then
  if (N_mult == 1) then
     open  (iu, file=Output_File)
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
     enddo
  endif
  endif
!.......................................................................


! prepare grid/radial range for poincare plot ..........................
  ! use radial range R_start -> R_end
  if (R_start > 0.d0) then
     dr   = 0.d0
     my_grid => new_grid(N_steps+1, log_progress=.false.)

     if (N_steps .ne. 1) then
        dr = (R_end - R_start) / N_steps
        if (firstP) write (6,2001) R_start, R_end, N_steps
        do ig=0,N_steps
           my_grid(ig+1,1) = R_start + ig*dr
        enddo
     else
        if (firstP) write (6,2002) R_start
        my_grid(1,1) = R_start
     endif

  ! use x_start if N_steps is not specified
  elseif (N_steps == 0  .and.  x_start(1).ne.0.d0) then
     my_grid => new_grid(N_steps+1, log_progress=.false.)
     my_grid(1,:) = x_start

  ! use points from grid file if R_start is missing
  elseif (Grid_File .ne. '') then
     call read_grid (Grid_File, log_progress=.false., use_coordinates=COORDINATES(min(Trace_Coords,2)))
     if (firstP) write (6,*)

  else
     if (firstP) write (6, *) 'error: starting points undefined!'
     stop
  endif
!.......................................................................


!.......................................................................
! prepare output data array
  Pdata%n_sets = n_grid
  Pdata%n_cuts = N_mult
  allocate (Pdata%n_points(0:n_grid-1, 0:N_mult-1))
  Pdata%n_points = 0
  allocate (Pdata%X       (0:n_grid-1, 0:N_mult-1, N_points, 4))
  Pdata%X        = 0.d0
!.......................................................................


! main loop (begin) ....................................................
  ig = mype
  grid_loop: do
     call get_next_grid_point (iflag, yc)
     if (iflag.eq.-1) exit grid_loop

     ! integrate connection length
     lc = 0.d0


     call F%init_toroidal_tracing(yc, Trace_Step, Trace_Method, Trace_Coords, N_sym*N_mult, Phi_output)
     icut = 0

     ! start field line tracing
     trace_loop: do
        call F%trace_1step()
        lc = lc + Trace_Step


        ! check intersection with Poincare plane
        if (F%intersect_sym_plane(icut, X)) then
           imult = int(mod(icut,N_mult))
           if (imult < 0) imult = imult + N_mult
           j     = Pdata%n_points(ig,imult) + 1
           Pdata%n_points(ig,imult) = j

           Pdata%X(ig,imult,j,1:2)  = X(1:2)
           theta           = get_poloidal_angle(X)
           if (theta.le.0.d0) theta = theta + pi2
           Pdata%X(ig,imult,j,3  ) = theta / pi2 * 360.d0
           Psi             = get_Psi(X)
           Psi             = (Psi-Psi_axis)/(Psi_sepx-Psi_axis)
           Pdata%X(ig,imult,j,4  ) = Psi
        endif



        ! check intersection with boundaries
        if (F%intersect_boundary(X)) then
           write (6,4000) ig, abs(lc/1.d2), abs(icut)/N_mult
           exit trace_loop
        endif


        ! check upper limit for field line tracing
        if (abs(icut) .ge. N_points*N_mult) then
           write (6,4001) ig, abs(lc/1.d2), abs(icut)/N_mult
           exit trace_loop
        endif
     enddo trace_loop

     ig = ig + nprs
  enddo grid_loop
! ... main loop (end) ..................................................


! finalize .............................................................
  call wait_pe()
  call sum_inte_data (Pdata%n_points, n_grid*N_mult)
  call sum_real_data (Pdata%X       , n_grid*N_mult*N_points*4)

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
 1003 format (8x,'Reference location     ',f7.3,' deg')
 2001 format (8x,'radial domain:',5x,'R_start = ',f6.2,5x, &
              'R_end = ',f6.2,5x,'with ',i4,' steps'/)
 2002 format (8x,'position:',5x,'R = ',f6.2,5x/)

 4000 format (5x,i5,',',8x,'L_c = ',f9.2,' m,',8x,'n_points = ',i5)
 4001 format (5x,i5,',',8x,'L_c > ',f9.2,' m,',8x,'n_points = ',i5)
end subroutine poincare_plot
