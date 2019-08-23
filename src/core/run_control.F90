#include "../../config.h"
module run_control
  use parallel
  use numerics
  use diffusion
  implicit none


  character(len=*), parameter :: data_dir = DATABASE_DIR

  character(len=*), parameter :: &
     AUTOMATIC = 'automatic', &
     MANUAL    = 'manual'


  integer, parameter :: &
     ZERO_TOLERANCE = 0, &
     INERVOUS   =  1, &
     IMODERATE  = 10, &
     IDONTPANIC = 42

  ! user defined variables
  character(len=256) :: &
     Machine        = ' ', &        ! select input directory (1st part)
     Configuration  = ' ', &        ! select input directory (2nd part)
     Boundary       = '', &
     Run_Type       = ' ', &        ! select sub-program to execute
     Run_Mode       = 'automatic', &
     Label          = '', &
     Output_File    = 'output.txt', &
     Grid_File      = 'grid.dat'

  real*8 :: &
     x_start(3)     = 0.d0, &       ! initial position for field line tracing
     Limit          = 2.d4, &       ! maximum distance for field line tracing (in one direction)
     R_start        = 0.d0, &       ! radial start- and
     R_end          = 0.d0, &       ! end position (e.g. for Poincare plots)
     Z_start        = 0.d0, &       ! vertical start- and
     Z_end          = 0.d0, &       ! end position (e.g. for Poincare plots)
     Phi_output     = 0.d0, &       ! Reference plane for Poincare plots
     Theta(2)       = 0.d0, &
     Psi(2)         = 0.d0, &
     offset         = 1.d-2, &
     tolerance      = 0.2d0, &
     max_pt         = 0.d0, &       ! max. number of poloidal turns
     Dfl            = 0.d0          ! field line diffusion coeff [cm**2 / cm]


  integer :: &
     N_steps        = 0, &          ! Number of discrete steps
     N_points       = 0, &          ! Max. number of points for Poincare plots
     N_sym          = 1, &          ! Toroidal symmetry factor (for Poincare plots)
     N_mult         = 1, &          !
     N_theta        = 0, &          ! Resolution in poloidal direction
     N_psi          = 1, &          ! Resolution in radial direction
     N_phi          = 1, &          ! Resolution in toroidal direction
     N_R            = 1, &          ! Resolution in R direction
     N_Z            = 1, &          ! Resolution in Z direction
     Run_Level(2)   = 0, &
     Side           = 1, &
     Input_Format   = 1, &
     Output_Format  = 1, &          ! See individual tools
     Panic_Level    = IMODERATE, &
     Surface_Id     = 1, &
     Surface_Type   = 1, &
     field_line_diffusion_model = DIFFUSION_RZ_CIRCLE


  logical :: &
     Debug          = .false., &
     use_boundary_from_equilibrium   = .true., &
     stop_at_boundary                = .true.



  ! internal variables
  character(len=256) :: Prefix = '', &
                        Boundary_Prefix, &
                        Bfield_input_file, &
                        run_control_file


  namelist /RunControl/ &
     Machine, Configuration, Boundary, &
     Run_Type, Output_File, Label, Grid_File, Input_Format, Output_Format, Panic_Level, &
     x_start, N_steps, Limit, max_pt, &
     R_start, R_end, Z_start, Z_end, Phi_output, N_points, N_sym, N_mult, Side, &
     Theta, Psi, N_theta, N_psi, N_phi, N_R, N_Z, offset, tolerance, &
     Run_Level, &
     Surface_Id, Surface_Type, &
     Debug, use_boundary_from_equilibrium, stop_at_boundary, Run_Mode, &
     Dfl, field_line_diffusion_model


  private :: Boundary, RunControl
  contains
!=======================================================================


!=======================================================================
  subroutine load_run_control(input_file)
  use math
  character(len=*), intent(in) :: input_file

  integer, parameter :: iu = 23


  run_control_file = input_file
  ! load run control on first processor
  if (firstP) then
     open  (iu, file=input_file, err=5000)
     read  (iu, RunControl, end=5000)
     close (iu)

     Boundary_Prefix = ''
     if (Machine .ne. ' ') then
        write (6, *)
        write (6, 1000)
        write (6, 1001) trim(Machine)
        write (6, 1002) trim(Configuration)
        Prefix = data_dir//'/'//trim(Machine)//'/'// &
                 trim(Configuration)//'/'
        if (Boundary .ne. '') then
           Boundary_Prefix = data_dir//'/'// &
                             trim(Machine)//'/'//trim(Boundary)//'/'
        endif
     else
        Prefix = './'
     endif

     Bfield_input_file = trim(Prefix)//'bfield.conf'


     ! post-processing of user-defined variables
     Phi_output = Phi_output / 180.d0 * pi
  endif


  call setup_diffusion(Dfl, field_line_diffusion_model)


  ! broadcast data to other processors
  if (nprs == 1) return
  call wait_pe()
  call broadcast_char   (Run_Type   , 120)
  call broadcast_char   (Run_Mode   , 120)
  call broadcast_char   (Grid_File  , 120)
  call broadcast_char   (Output_File, 120)
  call broadcast_real   (x_start    ,   3)
  call broadcast_real_s (Limit           )
  call broadcast_real_s (max_pt          )
  call broadcast_real_s (R_start         )
  call broadcast_real_s (R_end           )
  call broadcast_real_s (Z_start         )
  call broadcast_real_s (Z_end           )
  call broadcast_real_s (Phi_output      )
  call broadcast_real   (Theta      ,   2)
  call broadcast_real   (Psi        ,   2)
  call broadcast_real_s (offset          )
  call broadcast_real_s (tolerance       )
  call broadcast_inte_s (N_steps         )
  call broadcast_inte_s (N_points        )
  call broadcast_inte_s (N_sym           )
  call broadcast_inte_s (N_mult          )
  call broadcast_inte_s (N_theta         )
  call broadcast_inte_s (N_psi           )
  call broadcast_inte_s (N_phi           )
  call broadcast_inte_s (N_R             )
  call broadcast_inte_s (N_Z             )
  call broadcast_inte_s (Side            )
  call broadcast_inte_s (Input_Format    )
  call broadcast_inte_s (Output_Format   )
  call broadcast_inte_s (Panic_Level     )
  call broadcast_logi   (Debug           )
  call broadcast_logi   (stop_at_boundary)
  call broadcast_real_s (Dfl             )
  call broadcast_inte_s (field_line_diffusion_model)


  return
 1000 format(3x,'- Magnetic setup:')
 1001 format(8x,'Machine:                ',a)
 1002 format(8x,'Configuration:          ',a)
 5000 write  (6,5001)
 5001 format ('error reading control table from input file')
  stop
  end subroutine load_run_control
!=======================================================================



!=======================================================================
  subroutine run_control_main()
  !use grid
  !use dataset

  !type(t_grid)    :: G
  !type(t_dataset) :: D

  integer :: i, ind


  if (firstP) then
     write (6, 1000)
     write (6, 1000)
     write (6, *) 'Main program:'
     write (6, *)
  endif


  select case (Run_Type)
  case ('sample_bfield')
     call sample_bfield
  case ('trace_bline')
     call trace_bline
  case ('poincare_plot')
     call poincare_plot
  case ('connection_length')
     call connection_length
  case ('initialize_equilibrium')
     call initialize_equilibrium()
  case ('get_equi_info_2D')
     call get_equi_info_2D
  case ('generate_flux_surface_2D')
     call generate_flux_surface_2D
  case ('generate_flux_surface_3D')
     call generate_flux_surface_3D
  case ('plot_boundary')
     call plot_boundary
  case ('safety_factor')
     call safety_factor
  case ('transform_to_flux_coordinates')
     call transform_to_flux_coordinates
  case ('generate_mag_file')
     call generate_mag_file
  case ('generate_magnetic_axis')
     call generate_magnetic_axis()
  case ('flux_surface_grid')
     call flux_surface_grid
  case ('field_line_loss')
     call field_line_loss
  case ('generate_separatrix')
     call generate_separatrix
  case ('generate_rpath')
     call generate_rpath()
  case ('footprint_grid')
     call footprint_grid
  case ('setup_distance_to_surface')
     call setup_distance_to_surface
  case ('evaluate_distance_to_surface')
     call evaluate_distance_to_surface
  case ('separatrix_manifolds')
     call separatrix_manifolds()
  case ('generate_flux_tube')
     call generate_flux_tube
  case ('FLR_analysis')
     call FLR_analysis
  case ('melnikov_function')
     call melnikov_function()
  case ('generate_3D_fieldline_grid')
     call generate_3D_fieldline_grid(Run_Level(1), Run_Level(2))
  case ('critical_point_analysis')
     !call G%load(Grid_File)
     !call critical_point_analysis(G, D, ind)
     !call D%store(filename=Output_File)
  case ('export_gfile')
     call export_gfile()
  case ('quasi_surfaces')
     call quasi_surfaces()
  case ('perturbation_spectrum')
     call perturbation_spectrum()
  case ('equi2d_geometry')
     call equi2d_geometry(N_R)
  case default
!     if (Run_Type(1:27) == 'generate_field_aligned_grid') then
!        read (Run_Type(40:42), *) i
!        call generate_field_aligend_grid (i)
!     else
        call run_control_development(Run_Type)
!     endif
  end select

 1000 format (/ '========================================================================')
  end subroutine run_control_main
!=======================================================================

end module run_control
