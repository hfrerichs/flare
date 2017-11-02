module grid3dgen_interface
  use types
  implicit none

  integer :: nx, nrpath, ni

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



  ! return toroidal position of base mesh location ---------------------
  subroutine base_position(iblock, phi, ierr)
  use fieldline_grid, only: blocks, Block
  integer, intent(in)   :: iblock
  real(dp), intent(out) :: phi
  integer,  intent(out) :: ierr


  ierr = 0
  if (iblock < 0  .or.  iblock >= blocks) then
     ierr = 1
     return
  endif
  phi = Block(iblock)%phi_base

  end subroutine base_position
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  subroutine level(ilevel, isublevel, ierr)
  use fieldline_grid, only: blocks
  use base_mesh
  integer, intent(in)  :: ilevel, isublevel
  integer, intent(out) :: ierr

  integer :: iblock


  ierr = 0
  select case(ilevel)
  ! 1. generate innermost boundaries
  case(1)
     call generate_innermost_boundaries()

  ! 2. generate base mesh
  case(2)
     select case(isublevel)
     ! 2.1. set up geometry of computational domain (magnetic axis, separatrix, radial paths)
     case(1)
        call setup_topology(nx, nrpath)
        call setup_geometry()
     ! 2.2. set up interfaces between zones
     case(2)
        call setup_interfaces(ni)
     ! 2.3. generate mesh now
     case(3)
        do iblock=0,blocks-1
           call generate_base_mesh(iblock)
        enddo
     case default
        ierr = 2
     end select

  case default
     ierr = 1
  end select

  end subroutine level
  !---------------------------------------------------------------------

end module grid3dgen_interface
