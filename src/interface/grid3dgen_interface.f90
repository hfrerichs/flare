module grid3dgen_interface
  use types
  implicit none

  !integer, parameter :: dp = selected_real_kind(15,307) ! double precision

  contains
  !---------------------------------------------------------------------


  ! initialize grid generation module ----------------------------------
  subroutine init()
  use fieldline_grid, only: setup_grid_configuration


  call setup_grid_configuration()

  end subroutine init
  !---------------------------------------------------------------------



  ! return number of toroidal blocks -----------------------------------
  function get_block_num()
  use fieldline_grid, only: blocks
  integer get_block_num

  get_block_num = blocks
  end function get_block_num
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  subroutine base_position(iblock, phi, ierr)
  use fieldline_grid, only: blocks, Block
  integer, intent(in)   :: iblock
  real(dp), intent(out) :: phi
  integer,  intent(out) :: ierr
!f2py intent(in)  iblock
!f2py intent(out) phi
!f2py intent(out) ierr


  ierr = 0
  if (iblock < 0  .or.  iblock >= blocks) then
     ierr = 1
     return
  endif
  phi = Block(iblock)%phi_base
  write (6, *) 'phi_base(',iblock,') = ', phi

  end subroutine base_position
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  subroutine level(ilevel, ierr)
  use base_mesh
  integer, intent(in)  :: ilevel
  integer, intent(out) :: ierr


  ierr = 0
  select case(ilevel)
  ! 1. generate innermost boundaries
  case(1)
     call generate_innermost_boundaries()

  ! 2. generate base mesh
  case(2)
     call setup_topology()
     call setup_geometry()
     call setup_interfaces()

  case default
     ierr = 1
  end select

  end subroutine level
  !---------------------------------------------------------------------

end module grid3dgen_interface
