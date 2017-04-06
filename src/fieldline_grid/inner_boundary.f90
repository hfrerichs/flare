module inner_boundary
  use iso_fortran_env
  use fieldline_grid, only: blocks, Block, x_in2, nr_perturbed
  use curve2D
  implicit none

  type(t_curve), dimension(:,:), allocatable :: C_in

  real(real64), dimension(:,:), allocatable :: DPsiN1
  real(real64), dimension(:),   allocatable :: PsiN1
  real(real64) :: PsiN_in


  contains
  !=====================================================================


  !=====================================================================
  subroutine load_inner_boundaries (theta0, sample_method)
  use magnetic_axis
  use equilibrium, only: get_PsiN
  use math
  use run_control, only: Debug
  real(real64), intent(in), optional :: theta0
  integer,      intent(in), optional :: sample_method

  character(len=72) :: filename
  real(real64)     :: phi, r3(3), Pmag(2), delta, PsiN_x
  integer          :: i, iblock, n, sort_method, apply_sample_method


  write (6, 1000)
  write (6, 1001)

  ! set reference direction
  sort_method = DISTANCE
  if (present(theta0)) then
     sort_method = ANGLE
  endif

  apply_sample_method = sort_method
  if (present(sample_method)) apply_sample_method = sample_method


  allocate (C_in(0:blocks-1,0:1))
  allocate (PsiN1(0:blocks-1), DPsiN1(0:blocks-1,-1:1))
  PsiN_in       = 0.d0
  PsiN1         = 0.d0
  DPsiN1        = 0.d0
  DPsiN1(:,-1)  = 1.d0
  do iblock=0,blocks-1
     phi  = Block(iblock)%phi_base / 180.d0 * pi
     r3   = get_magnetic_axis(phi)
     Pmag = r3(1:2)

     ! read boundary data from file and set up flux surface
     do i=0,nr_perturbed-1
        write (filename, 2000) i, iblock
        call C_in(iblock,i)%load(filename,output=SILENT)
        if (C_in(iblock,i)%n_seg <= 0) then
           write (6, 9000) iblock, i
           stop
        endif

        call C_in(iblock,i)%sort_loop(Pmag, theta0, sort_method)
        if (apply_sample_method .ne. sort_method) then
           call C_in(iblock,i)%setup_sampling_by_method(apply_sample_method, Pmag)
        endif
        if (Debug) then
           write (filename, 2001) i, iblock
           call C_in(iblock,i)%plot(filename=filename)
        endif
     enddo

     ! calculate mean and variance of PsiN on (perturbed) flux surface
     n = C_in(iblock,nr_perturbed-1)%n_seg
     do i=1,n
        PsiN_x            = get_PsiN(C_in(iblock,nr_perturbed-1)%x(i,:))
        delta             = PsiN_x - PsiN1(iblock)
        PsiN1(iblock)     = PsiN1(iblock) + delta/i
        DPsiN1(iblock, 0) = DPsiN1(iblock,0) + delta*(PsiN_x - PsiN1(iblock))
        DPsiN1(iblock,-1) = min(DPsiN1(iblock,-1), PsiN_x)
        DPsiN1(iblock, 1) = max(DPsiN1(iblock, 1), PsiN_x)

        PsiN_in           = PsiN_in + PsiN_x / n
     enddo
     DPsiN1(iblock,0) = sqrt(DPsiN1(iblock,0) / (n - 1.d0))
     write (6, 1002), iblock, PsiN1(iblock), DPsiN1(iblock,0), DPsiN1(iblock,-1), DPsiN1(iblock,1)
  enddo
  PsiN_in = PsiN_in / blocks

 1000 format (3x,'- Inner simulation boundaries:')
 1001 format (8x,'Block   PsiN_mean                Min(PsiN)   Max(PsiN)')
 1002 format (8x,i4,4x,f8.6,' +/- ',f8.6,2(4x,f8.6))
 2000 format ('fsin',i0,'_',i0,'.txt')
 2001 format ('fsin',i0,'_',i0,'.plt')
 9000 format ('error: boundary ',i0,' in block ',i0,' undefined!')
  end subroutine load_inner_boundaries
  !=====================================================================


  !=====================================================================
  ! return width of high pressure region (HPR) at poloidal angle of X-point
  !=====================================================================
  function get_d_HPR (Px, Pmag)
  use flux_surface_2D
  use equilibrium
  use math
  real(real64), intent(in) :: Px(2), Pmag(2)
  real(real64)             :: get_d_HPR(2)

  type(t_flux_surface_2D)  :: F
  real(real64), save       :: d_HPR(2) = 0.d0
  real(real64) :: x(2), d, dx(2), theta, theta0, r3(3), xi
  integer      :: i


  if (d_HPR(1) > 0) then
     get_d_HPR = d_HPR
     return
  endif

  r3      = x_in2
  x       = r3(1:2)
  theta0  = get_poloidal_angle(r3)
  call F%generate_closed(x, RIGHT_HANDED)
  call F%setup_angular_sampling(Pmag)

  r3(1:2) = Px
  r3(3)   = 0.d0
  theta   = get_poloidal_angle(r3)
  xi      = (-theta0 + theta)/pi2
  if (xi < 0) xi = xi + 1.d0

  call F%sample_at(xi, dx)
  d_HPR     = dx-Px
  get_d_HPR = d_HPR

  end function get_d_HPR
  !=====================================================================


  !=====================================================================
  ! Setup discretization of inner boundaries for block iblock
  ! Input:
  !    S    poloidal spacing function
  !    
  !=====================================================================
  subroutine setup_inner_boundaries(G, iblock, sampling_method, S)
  use grid
  use mesh_spacing
  type(t_grid),    intent(inout) :: G
  integer,         intent(in)    :: iblock, sampling_method
  type(t_spacing), intent(in)    :: S

  real(real64) :: xi, x(2)
  integer :: i, j, np


  np = G%n2-1
  do j=0,np
     xi = S%node(j,np)
     do i=0,nr_perturbed-1
        call C_in(iblock,i)%sample_at(xi, x)
        G%mesh(i, j, :) = x
     enddo
  enddo

  end subroutine setup_inner_boundaries
  !=====================================================================
  subroutine setup_inner_boundary(iblock, ir, S, G)
  use grid
  use mesh_spacing
  integer,         intent(in)    :: iblock, ir
  type(t_grid),    intent(inout) :: G
  type(t_spacing), intent(in)    :: S

  real(real64) :: xi, x(2)
  integer :: j, np


  if (ir < 0  .or.  ir > nr_perturbed-1) then
     write (6, *) 'error in setup_inner_boundary: ir out of range!'
     write (6, *) '0 <= ir <= ', nr_perturbed-1, ' required!'
     stop
  endif


  np = G%n2-1
  do j=0,np
     xi = S%node(j,np)
     call C_in(iblock,ir)%sample_at(xi, x)
     G%mesh(ir, j, :) = x
  enddo

  end subroutine setup_inner_boundary
  !=====================================================================
  subroutine setup_inner_boundary0(iblock, G)
  use grid
  use mesh_spacing
  integer,         intent(in)    :: iblock
  type(t_grid),    intent(inout) :: G

  real(real64) :: x1(2), x2(2), x(2)
  integer :: j, np


  np = G%n2-1
  do j=0,np
     x1 = G%mesh(1, j, :)
     x2 = G%mesh(2, j, :)
     if (intersect_curve(x1, x2, C_in(iblock,0), x, intersect_mode=-1)) then
        G%mesh(0, j, :) = x
     else
        write (6, *) 'error in setup_inner_boundary0 at grid node ', j
        write (6, *) 'x = ', x1
        stop
     endif
  enddo

  end subroutine setup_inner_boundary0
  !=====================================================================

end module inner_boundary







!=======================================================================
! Generate pair of innermost boundaries
!=======================================================================
  subroutine generate_innermost_boundaries
  use iso_fortran_env
  use fieldline_grid
  implicit none


  write (6, *)
  write (6, 1000)
  write (6, *)

  select case (Innermost_Flux_Surface)
  case (SF_EXACT)
     call exact_surface (x_in1, 'fsin0')
     call exact_surface (x_in2, 'fsin1')
  case (SF_QUASI)
     call quasi_surface (x_in1, 'fsin0')
     call quasi_surface (x_in2, 'fsin1')
  case default
     write (6, *) 'error: flux surface type ', trim(Innermost_Flux_Surface), ' not defined!'
     stop
  end select
  write (6, *) 'finished generating innermost boundaries'

  return
 1000 format(3x,'- Generate inner simulation boundaries (based on Poincare plots)')
 9000 write (6, *) 'error while reading input file ', trim(config_file), '!'
  stop
  contains
!-----------------------------------------------------------------------
  subroutine exact_surface (x, s5)
  use run_control, only: N_mult, N_sym, N_points, x_start, Phi_output, Output_File

  real(real64),     intent(in) :: x(3)
  character(len=5), intent(in) :: s5

  character(len=3) :: sblock
  integer :: iblock


  ! set default number of points for Poincare plot
  if (N_points == 0) N_points = 1000

  N_sym    = symmetry
  x_start  = x
  if (default_decomposition) then
     N_mult   = blocks
     Output_File = s5//'.txt'
     if (blocks == 1) Output_File = s5//'_0.txt'
     call poincare_plot()
  else
     N_mult   = 1
     do iblock = 0,blocks-1
        write (sblock, '(i3)') iblock
        Output_File = s5//'_'//trim(adjustl(sblock))//'.txt'
        Phi_output  = Block(iblock)%phi_base / 180.d0 * pi
        call poincare_plot()
     enddo
  endif


  end subroutine exact_surface
!-----------------------------------------------------------------------
  subroutine quasi_surface (x, s5)
  real(real64),     intent(in) :: x(3)
  character(len=5), intent(in) :: s5

  write (6, *) 'generation of quasi flux surfaces not yet implemented!'
  stop
  end subroutine quasi_surface
!-----------------------------------------------------------------------
  end subroutine generate_innermost_boundaries
!=======================================================================
