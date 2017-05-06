! Magnetic field coils
module polygones
  use iso_fortran_env
  use mfc_polygon
  implicit none
  private


  integer, parameter :: max_coil_sets = 100


  character(len=120) :: &
     Poly_file(max_coil_sets)   = ''   ! input file(s) with data for polygones

  real(real64)       :: &
     I_scale(max_coil_sets)     = 1.d0 ! scale factor(s) for current through polygones

  integer            :: n_sets = 1


  namelist /Polygones_Input/ &
     n_sets, I_scale, Poly_file


  type(t_mfc_polygon), dimension(:), allocatable :: C
  integer :: n_coils


  public :: &
     read_polygones_config, &
     broadcast_mod_polygones, &
     get_Bcart_polygones, &
     get_Bcyl_polygones, &
     get_JBf_cyl_polygones

  contains
!=======================================================================



!=======================================================================
  subroutine read_polygones_config(iu, iconfig, Prefix)
  integer,          intent(in)  :: iu
  integer,          intent(out) :: iconfig
  character(len=*), intent(in)  :: Prefix


  rewind (iu)
  read   (iu, Polygones_Input, end=1000)
  iconfig = 1
  write (6, *)

  if (sum(abs(I_scale)).ne.0.d0) then
     call setup_polygones (Prefix)
  else
     iconfig = 0
  endif

  return
 1000 iconfig = 0
  end subroutine read_polygones_config
!=======================================================================



!=======================================================================
  subroutine setup_polygones(Prefix)
  character(len=*), intent(in) :: Prefix

  integer, parameter :: iu = 44

  character(len=256) :: input_file
  character(len=80)  :: vv80
  integer :: i, io, j, j0, n


  ! 1. scan through files to get total number of polygones
  n_coils = 0
  do i=1,n_sets
     input_file = trim(Prefix)//Poly_file(i)
     open  (iu, file=input_file, iostat=io)
     if (io.ne.0) then
        write (6, 2000) trim(input_file)
        stop
     endif

     read  (iu, 1001) vv80
     read  (iu, *) n
     n_coils = n_coils + n
     close (iu)
  enddo


  ! 2. set up polygones
  allocate (C(n_coils))
  j0 = 0
  do i=1,n_sets
     input_file = trim(Prefix)//Poly_file(i)
     open  (iu, file=input_file, iostat=io)
     read  (iu, 1001) vv80
     read  (iu, *) n
     write (6,  1002) vv80

     do j=1,n
        call C(j0+j)%load(iu, I_scale(i))
     enddo
     close (iu)

     j0 = j0 + n
  enddo

 2000 format ('error reading polygon file: ',a)
 1001 format (a80)
 1002 format (3x,'- ',a80)
  end subroutine setup_polygones
!=======================================================================



!=======================================================================
  subroutine broadcast_mod_polygones()
  use parallel

  integer :: i


  call broadcast_inte_s (n_coils)
  do i=1,n_coils
     call C(i)%broadcast()
  enddo

  end subroutine broadcast_mod_polygones
!=======================================================================



!=======================================================================
  function get_Bcart_polygones(x) result(Bf)
  real(real64), intent(in) :: x(3)
  real(real64)             :: Bf(3)

  integer :: i


  Bf = 0.d0
  do i=1,n_coils
     Bf = Bf + C(i)%get_Bf(x)
  enddo

  end function get_Bcart_polygones
!=======================================================================



!=======================================================================
  function get_Bcyl_polygones(r) result(Bf)
  real(real64), intent(in) :: r(3)
  real(real64)             :: Bf(3)

  real(real64) :: Bcart(3), x(3), sin_phi, cos_phi
  integer      :: i


  sin_phi = sin(r(3))
  cos_phi = cos(r(3))
  x(1)    = r(1) * cos_phi
  x(2)    = r(1) * sin_phi
  x(3)    = r(2)
  Bcart   = get_Bcart_polygones(x)

  Bf(1)   =  Bcart(1) * cos_phi  +  Bcart(2) * sin_phi
  Bf(2)   =  Bcart(3)
  Bf(3)   = -Bcart(1) * sin_phi  +  Bcart(2) * cos_phi

  end function get_Bcyl_polygones
!=======================================================================



!=======================================================================
  function get_JBf_cyl_polygones(r) result(JBf)
  real(real64), intent(in) :: r(3)
  real(real64)             :: JBf(3,3)

  integer :: i


  JBf = 0.d0
  do i=1,n_coils
     JBf = JBf + C(i)%get_JBf_cyl(r)
  enddo

  JBf = JBf * 1.d2 ! T/m -> Gauss/cm

  end function get_JBf_cyl_polygones
!=======================================================================

end module polygones
