module reconstruct
  implicit none
  private

  type t_fluxtube_coords
     real*8  :: tj, pj, rj
     integer :: ir, ip, it, iz

     contains
     procedure :: trace_1step
     !procedure :: to_real_coordinates
  end type t_fluxtube_coords


  character*120 :: &
     Grid_File   = 'grid.data', &
     Grid_Layout = 'grid.info', &
     Bfield_File = 'bfield.data'

  namelist /Reconstruct_Input/ &
     Grid_File, Grid_Layout, Bfield_File


  type (t_fluxtube_coords) :: yc
  real*8            :: sc, ds


  public :: &
     load_reconstruct_config, &
     broadcast_mod_reconstruct, &
     t_fluxtube_coords

  contains
!=======================================================================



!=======================================================================
  ! use Prefix as argument -> benchmark routines may load the grid in addition to the original field
  !subroutine load_reconstruct_config (iu, Prefix, iload)
  subroutine load_reconstruct_config (iu, iload)
  use run_control, only: Prefix, Trace_Coords
  use math
#if defined(EMC3)
  use GEOMETRY_PL
  use SURFACE_PL
  use MAPPING
#endif
  integer, intent(in)  :: iu
  integer, intent(out) :: iload

  integer, parameter :: iu_Grid_Layout = 25, &
                        iu_Grid_Data   = 26, &
                        iu_Bfield      = 27


  Trace_Coords = CYLINDRICAL

#if defined(EMC3)
  rewind (iu)
  read   (iu, Reconstruct_Input, end=1000)
  iload = 1
  write (6, *)
  write (6, 1001)

  open  (iu_Grid_Layout, file=trim(Prefix)//trim(Grid_Layout))
  open  (iu_Grid_Data,   file=trim(Prefix)//trim(Grid_File))
  open  (iu_Bfield,      file=trim(Prefix)//trim(Bfield_File))
  call READ_GEOMETRY_PL(iu_Grid_Layout, iu_Grid_Data, iu_Bfield)
  close (iu_Grid_Data)
  close (iu_Bfield)

  call SETUP_GEOMETRY_PL()
  call SETUP_SURFACE()

  call READ_NON_DEFAULT_SF(iu_Grid_Layout)
  call ACTIVIZE_MAPP_RAD()
  call ACTIVIZE_MAPP_POL()
  call ACTIVIZE_MAPP_TOR()

  close (iu_Grid_Layout)
  return
 1000 iload = 0
 1001 format ('   - Reconstruct from magnetic field aligned grid')
#endif
  end subroutine load_reconstruct_config
!=======================================================================



!=======================================================================
  subroutine init_fl_reconstruction (r, ds_)
  real*8, intent(in) :: r(3), ds_

  integer :: ierr

#if defined(EMC3)
  call real_to_ft (r(1), r(2), r(3), &
      yc%iz, yc%ir, yc%ip, yc%it, yc%rj, yc%pj, yc%tj, 0, ierr)
  ds = ds_

  if (ierr.ne.0) then
     write (6, *) 'could not find flux-tube coordinates for ', r
     write (6, *) 'ierr = ', ierr
     stop
  endif

#endif
  end subroutine init_fl_reconstruction
!=======================================================================



!=======================================================================
  function trace_fl()
  real*8 :: trace_fl(3)
  end function trace_fl
!=======================================================================



!=======================================================================
  subroutine trace_1step(this)
  class(t_fluxtube_coords) :: this
  end subroutine trace_1step
!=======================================================================



!=======================================================================
  subroutine broadcast_mod_reconstruct
  use parallel

  if (nprs > 1) then
     if (firstP) write (6, *) 'module reconstruct does not support parallel execution!'
     stop
  endif

  end subroutine broadcast_mod_reconstruct
!=======================================================================

end module reconstruct
