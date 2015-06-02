!===============================================================================
! Generate homoclinic tangle for hyperbolic fixed point
!
! Input (taken from run control file):
! Grid_File             toroidal discretization of X-point (+ small offset for tracing)
! Phi_Output            toroidal reference position [deg]
! N_sym                 toroidal symmetry of the configuration
! N_mult                multiplicity of grid point (in terms of 2 pi / N_sym)
!
! Trace_Step, Trace_Method
! Output_File
!===============================================================================
subroutine homoclinic_tangle
  use iso_fortran_env
  use run_control, only: Grid_File, Output_File, Trace_Step, Trace_Method, Phi_Output, N_sym, N_mult
  use grid
  use fieldline
  use boundary
  use math
  use dataset
  use parallel
  implicit none

  integer, parameter :: iu = 54

  type(t_fieldline)  :: F
  type(t_grid)       :: G
  type(t_dataset)    :: D
  real(real64)       :: y(3), yh(3), Dphi
  integer            :: i, j, k, ind, ierr


  ! initialize
  if (firstP) then
     write (6, *) 'Generate homoclinic tangle for hyperbolic fixed point'
     write (6, *)
  endif
  call G%load(Grid_File)
  write (6, *)
  

  call D%new(G%nodes() * N_mult, 3)
  open  (iu, file='strike_point_'//trim(Output_File)//'.dat')
  grid_loop: do i=1,G%nodes()
     y = G%node(i)
     write (6, *) y(3) / pi * 180.d0

     ! set initial location
     call F%init(y, Trace_Step, Trace_Method, CYLINDRICAL)

     ! initial trace to reference plane
     Dphi = abs(Phi_Output / 180.d0 * pi - y(3))
     call F%trace_Dphi(Dphi, .false., yh, ierr)
     D%x(i,:) = yh


     ! trace to symmetry planes
     Dphi = pi2 / N_sym
     do j=1,N_mult-1
        call F%trace_Dphi(Dphi, .true., yh, ierr)
        if (ierr > 0) then
           write (iu, *) yh
           exit
        endif

        ind = i + j*G%nodes()
        D%x(ind,:) = yh
     enddo
  enddo grid_loop
  close (iu)


  ! write RZ-cut of manifold
  open  (iu, file='manifold_'//trim(Output_File)//'.dat')
  do j=0,N_mult-1
     do i=1,G%nodes()
        ind = i + j*G%nodes()
        if (D%x(ind,1) > 0.d0) write (iu, *) D%x(ind,1:2)
     enddo
  enddo
  close (iu)

end subroutine homoclinic_tangle
