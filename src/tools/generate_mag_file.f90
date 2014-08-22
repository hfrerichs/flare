!===============================================================================
! Generate data file for interpolation of the magnetic field
!
! Input (taken from run control file):
!    Grid_File          Input file (namelist Grid_Layout)
!    Output_File
!    Output_Format      = 1: ASCII data
!                       = 2: binary data
!
! Input (taken from Grid_File):
!    Title              Short (72 characters) description of the field configuration
!    N_sym              Toroidal symmetry number, output is generated for
!                       Phi = 0 -> 360/N_sym deg
!    N_phi, N_R, N_Z    Number of computational nodes in toroidal, radial and
!                       vertical direction
!    R_center           Center of computational box [cm] in radial direction
!    width, height      Width and height of computational box [cm], the lower
!                       left corner is at [R_center - width/2, -height/2], the
!                       upper right corner is at [R_center + width/2, height/2].
!===============================================================================
subroutine generate_mag_file
  use run_control, only: Grid_File, Output_File, Output_Format
  use interpolateB
  use parallel
  implicit none

  integer, parameter :: iu = 42

  integer, parameter :: &
     OUT_ASCII  = 1, &
     OUT_BINARY = 2

  character*72 :: &
     Title = 'pre-calculated magnetic field'

  real*8  :: &
     R_center = 1.d0, &        ! center of computational box [cm]
     width    = 1.d0, &        ! width and height of computational box [cm]
     height   = 1.d0
  integer :: N_sym, N_phi, N_R, N_Z

  namelist /Grid_Layout/ &
     Title, N_sym, N_phi, N_R, N_Z, R_center, width, height


  if (firstP) then
     write (6, *) 'Generate mag-file (pre-calculated magnetic field), output in: ', adjustl(trim(Output_File))
     write (6, *)
  endif


  ! read grid layout
  open  (iu, file=Grid_file)
  read  (iu, Grid_Layout)
  close (iu)
  write (6, 1000) N_sym, 360.d0 / N_sym
 1000 format(3x,'- Toroidal symmetry: ',i4,'    (Phi = 0 -> ',f5.1,' deg)')


  ! sample magnetic field
  call generate_field_from_coils (N_sym, N_R, N_Z, N_phi, R_center, width, height)


  ! select output format and write data
  if (Output_Format == OUT_ASCII) then
     call write_ascii_data (title, Output_File)
  elseif (Output_Format == OUT_BINARY) then
     call write_binary_data (title, Output_File)
  else
     write (6, *) 'error: undefined Output_Format ', Output_Format
     stop
  endif

end subroutine generate_mag_file
