module separatrix_interface
  use types
  implicit none

  ! store branches of separatrix
  ! B_lc:    clockwise direction          (left hand) along core
  ! B_rc:    counter-clockwise direction (right hand) along core
  ! B_ld:    clockwise direction          (left hand) to divertor
  ! B_rd:    counter-clockwise direction (right hand) to divertor
  real(dp), dimension(:,:), allocatable :: B_lc, B_rc, B_ld, B_rd

  contains
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  ! Generate separatrix for X-point ix
  !---------------------------------------------------------------------
  subroutine generate(ix, stop_at_boundary)
  use separatrix
  integer, intent(in) :: ix
  logical, intent(in) :: stop_at_boundary

  type(t_separatrix) ::S


  call S%generate_iX(ix, stop_at_boundary=stop_at_boundary)
  if (allocated(B_lc)) deallocate(B_lc, B_rc, B_ld, B_rd)
  allocate (B_lc(S%M1%n_seg+1, 2));   B_lc = S%M1%x
  allocate (B_rc(S%M2%n_seg+1, 2));   B_rc = S%M2%x
  allocate (B_ld(S%M3%n_seg+1, 2));   B_ld = S%M3%x
  allocate (B_rd(S%M4%n_seg+1, 2));   B_rd = S%M4%x

  end subroutine generate
  !---------------------------------------------------------------------

end module separatrix_interface
