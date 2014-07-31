!===============================================================================
! Magnetic field from coil sets (given as polygones)
!
! THIS MODULE NEEDS CLEANUP !!!!
!===============================================================================
      module polygones
      implicit none

      integer, parameter :: max_coil_sets = 100

      real*8, parameter :: min_dist_to_coil = 1.d-8

      integer :: n_sets = 1

      real*8 :: I_scale(max_coil_sets) = 1.d0	! scale factor for current through polygones
      character*120 ::
     .    Poly_file(max_coil_sets)   = ''	! input file with data for polygones

      namelist /Polygones_Input/
     .         n_sets, I_scale, Poly_file


      integer ::
     .    n_poly   =  0,	! number of polygons
     .    n_poly1  =  0,	! with current correction factor
     .    n_poly2  =  0		! without current correction factor

      type t_polygon
         integer  ::
     .      n_seg	! number of segments
         real*8 :: I_poly	! current
         real*8, dimension(:), allocatable :: X,Y,Z
      end type t_polygon

      type (t_polygon), dimension(:), allocatable :: polygon_data

      contains
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
      subroutine read_polygones_config (iun, iconfig, Prefix)
      integer, intent(in)  :: iun
      integer, intent(out) :: iconfig
      character*120, intent(in) :: Prefix

      rewind (iun)
      read   (iun, Polygones_Input, end=1000)
      iconfig = 1
      write (6, *)
      !write (6, 1001)
      !write (6, 1002) n_sets

      if (sum(abs(I_scale)).ne.0.d0) then
         call setup_polygones (Prefix)
      else
         iconfig = 0
      endif

      return
 1000 iconfig = 0
 1001 format ('   - Coils given by polygones:')
 1002 format ('        Number of coil sets: ',i4)

      end subroutine read_polygones_config
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
! scan each polygon set and collect total number of polygones
c-----------------------------------------------------------------------
      subroutine prepare_polygones (Prefix)
      character*120, intent(in) :: Prefix

      integer, parameter :: iun = 44

      character*120 :: Input_File
      character*80  :: vv80
      integer :: i, io, n1, n2


      n_poly1 = 0
      n_poly2 = 0
      do i=1,n_sets
         Input_File = Prefix(1:len_trim(Prefix))//Poly_file(i)
         open (iun, file=Input_File, iostat=io)
         if (io.ne.0) then
            write (6,2000) Input_file
            stop
         endif

         read (iun, 1001) vv80
         read (iun,*)  n1, n2
         n_poly1 = n_poly1 + n1
         n_poly2 = n_poly2 + n2
         close (iun)
      enddo
      n_poly = n_poly1 + n_poly2
      if (n_poly.gt.0) allocate (polygon_data(n_poly))

      return
 2000 format ('error reading polygon file: ',a120)
 1001 format (a80)
      end subroutine prepare_polygones
c-----------------------------------------------------------------------
      subroutine setup_polygones (Prefix)
      character*120, intent(in) :: Prefix

      integer, parameter :: iun = 44

      character*80 :: vv80
      integer      :: n1, n2, np, ns, nperio,
     .                ip, is, io, i
      real*8       :: I_poly
      character*120:: Input_File


      call prepare_polygones(Prefix)

      np = 0
      do i=1,n_sets
         Input_File = Prefix(1:len_trim(Prefix))//Poly_file(i)

         ! read polygon set i
         open (iun, file=Input_File, iostat=io)
         if (io.ne.0) then
            write (6,2000) Input_file
            stop
         endif
         read (iun, 1001) vv80
         write (6, 1002) vv80
         read (iun,*)  n1, n2
         write (6, 1003) n1+n2

         do ip=1,n1+n2
            np = np + 1
            read  (iun,*) ns, I_poly, nperio
            polygon_data(np)%n_seg  = ns
            polygon_data(np)%I_poly = I_poly * I_scale(i)
            allocate (polygon_data(np)%X(0:ns),
     .                polygon_data(np)%Y(0:ns),
     .                polygon_data(np)%Z(0:ns))
            write (6,1000) ip, ns, I_poly * I_scale(i) / 1.d3
            read  (iun,*)  (polygon_data(np)%X(is),
     .          polygon_data(np)%Y(is), polygon_data(np)%Z(is), is=0,ns)

         enddo
         close (iun)
      enddo

      return
! 1000 format (8x,i3,' segments,',5x,
!     1           'current = ',1pd14.6)
 1000 format (8x,i4,': ',i4,' segments,',5x,'Ic = ',f10.3,' kA')
 1001 format (a80)
 1002 format ('   - ',a80)
 1003 format (8x,'Number of coils: ',i4)
 2000 format ('error reading polygon file: ',a120)
      end subroutine setup_polygones
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
      subroutine broadcast_polygones
      use parallel
      integer :: ip, ns

      !call broadcast_real_s (I_scale)
      call broadcast_inte_s (n_poly)
      call broadcast_inte_s (n_poly1)
      call broadcast_inte_s (n_poly2)

      if (mype.gt.0) then
         allocate (polygon_data(n_poly))
      endif

      do ip=1,n_poly
         call broadcast_inte_s (polygon_data(ip)%n_seg)
         call broadcast_real_s (polygon_data(ip)%I_poly)

         ns = polygon_data(ip)%n_seg
         if (mype.gt.0) then
            allocate (polygon_data(ip)%X(0:ns),
     .                polygon_data(ip)%Y(0:ns),
     .                polygon_data(ip)%Z(0:ns))
         endif
         call broadcast_real   (polygon_data(ip)%X, ns+1)
         call broadcast_real   (polygon_data(ip)%Y, ns+1)
         call broadcast_real   (polygon_data(ip)%Z, ns+1)
      enddo

      return
      end subroutine broadcast_polygones
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
      function get_Bcart_polygones(x) result(Bf)
      !subroutine add_Bf_polygones (igrad, Bf_, dBfdx_)
      !use position

      !integer, intent(in)     :: igrad
      !real*8,  intent(inout)  :: Bf_(3), dBfdx_(3,3)
      real*8, intent(in) :: x(3)
      real*8             :: Bf(3)

      integer :: igrad = 0

      integer  :: ip, is, is1, ns, ip_start
      real*8   :: Bf1(4:8),dBfdx(3,3), scaleF(n_poly), rp12, ak, bk, ck,
     .            rp2, som, s1, rx, ry, rz, var,
     .            p(6), dkx, dky, dkz
      real*8, dimension(:), allocatable :: dx, dy, dz, rp, rp1

      ip_start = 1
      !if (I_scale.eq.0.d0) ip_start = n_poly1+1

      !scaleF(        1:n_poly1) = I_scale
      !scaleF(n_poly1+1:n_poly ) = 1.d0

      Bf  = 0.d0
      Bf1 = 0.d0
      dBfdx = 0.d0
      do ip=ip_start,n_poly
         ns  = polygon_data(ip)%n_seg
         var = polygon_data(ip)%I_poly * 0.1d0	! mu_0 / 4pi * 100 cm/m * 10000 Gauss/T
         allocate (dx(0:ns), dy(0:ns), dz(0:ns), rp(0:ns), rp1(0:ns))

         do is=0,ns
            dx(is) = x(1) - polygon_data(ip)%X(is)
            dy(is) = x(2) - polygon_data(ip)%Y(is)
            dz(is) = x(3) - polygon_data(ip)%Z(is)
            rp2    = dx(is)**2 + dy(is)**2 + dz(is)**2
            rp2    = max(rp2,min_dist_to_coil**2)
            rp(is) = dsqrt(rp2)
            rp1(is)= 1.d0/rp2
         enddo

         do is=1,ns
            is1    = is-1
            rp12   = rp(is1) * rp(is)
            som    = rp(is1) + rp(is)
            s1     = 1.d0 /  (rp12 * (rp12 + dx(is1)*dx(is)
     +                              + dy(is1)*dy(is)
     +                              + dz(is1)*dz(is)))
            ck     = var * s1
            ak     = som * ck
            rx     = dy(is1)*dz(is) - dz(is1)*dy(is)
            ry     = dz(is1)*dx(is) - dx(is1)*dz(is)
            rz     = dx(is1)*dy(is) - dy(is1)*dx(is)
            Bf(1)  = Bf(1)  +  ak * rx
            Bf(2)  = Bf(2)  +  ak * ry
            Bf(3)  = Bf(3)  +  ak * rz

            if (igrad.eq.1) then
            bk     = -ak *som * s1
            p(1)   = rp(is1) * dx(is)
            p(2)   = rp(is)  * dx(is1)
            dkx    = bk * (p(1)+p(2))  
     -               -  ck * (p(2)*rp1(is1) + p(1)*rp1(is))
            p(3)   = rp(is1) * dy(is)
            p(4)   = rp(is)  * dy(is1)
            dky    = bk * (p(3)+p(4))
     -               -  ck * (p(4)*rp1(is1) + p(3)*rp1(is))
            p(5)   = rp(is1) * dz(is)
            p(6)   = rp(is)  * dz(is1)
            dkz    = bk * (p(5)+p(6))
     -               -  ck * (p(6)*rp1(is1) + p(5)*rp1(is))
            Bf1(4)  = Bf1(4)  +  dkx*rx
            Bf1(5)  = Bf1(5)  +  dky*rx + ak*(dz(is)  - dz(is1))
            Bf1(6)  = Bf1(6)  +  dkz*rx + ak*(dy(is1) - dy(is))
            Bf1(7)  = Bf1(7)  +  dky*ry
            Bf1(8)  = Bf1(8)  +  dkz*ry + ak*(dx(is)  - dx(is1))
            endif	! igrad.eq.1
         enddo
         deallocate (dx, dy, dz, rp, rp1)

         !Bf_       =  Bf_     +  Bf(1:3)
         !Bf_       =  Bf_     +  Bf(1:3) * scaleF(ip)
         if (igrad.eq.1) then
            dBfdx(1,1)=  Bf1(4)
            dBfdx(1,2)=  Bf1(5)
            dBfdx(1,3)=  Bf1(6)
            dBfdx(2,1)=  Bf1(5)
            dBfdx(2,2)=  Bf1(7)
            dBfdx(2,3)=  Bf1(8)
            dBfdx(3,1)=  Bf1(6)
            dBfdx(3,2)=  Bf1(8)
            dBfdx(3,3)=  -Bf1(4)-Bf1(7)
            !dBfdx_    =  dBfdx_  +  dBfdx
            !dBfdx_    =  dBfdx_  +  dBfdx * scaleF(ip)
         endif
      enddo

      return
      end function get_Bcart_polygones
      !end subroutine add_Bf_polygones
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
!      subroutine get_Bf_polygones_at_RZphi (R, Z, phi, Bf)
!      use position, only: set_position_from_RZphi
!
!      real*8, intent(in)  :: R, Z, phi
!      real*8, intent(out) :: Bf(3)
!
!      real*8 :: dBfdx(3,3)
!
!      call set_position_from_RZphi (R, Z, phi)
!
!      Bf    = 0.d0
!      dBfdx = 0.d0
!      call add_Bf_polygones (0, Bf, dBfdx)
!
!      return
!      end subroutine get_Bf_polygones_at_RZphi
c-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! calculate R,phi,Z components of magnetic field vector (Gauss) at R,Z (cm) and P (rad)
!-----------------------------------------------------------------------
      function get_Bcyl_polygones(r) result(Bf)
      !subroutine Bcyl_polygones (R, P, Z, BR, BP, BZ)
      !real*8, intent(in)  :: R, P, Z
      !real*8, intent(out) :: BR, BP, BZ
      real*8, intent(in) :: r(3)
      real*8             :: Bf(3)

      real*8 :: Bcart(3), x(3), sin_phi, cos_phi

      !call get_Bf_polygones_at_RZphi (R, Z, P, Bcart)
      sin_phi = sin(r(3))
      cos_phi = cos(r(3))
      x(1) = r(1) * cos_phi
      x(2) = r(1) * sin_phi
      x(3) = r(2)
      Bcart = get_Bcart_polygones(x)

      Bf(1) =  Bcart(1) * cos_phi  +  Bcart(2) * sin_phi
      Bf(2) =  Bcart(3)
      Bf(3) = -Bcart(1) * sin_phi  +  Bcart(2) * cos_phi

      end function get_Bcyl_polygones
      !end subroutine Bcyl_polygones
c-----------------------------------------------------------------------

      end module polygones
