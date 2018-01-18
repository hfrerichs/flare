module marsf
  use iso_fortran_env
  use bfield_component
  use math
  implicit none
  private


  integer, parameter :: B_FIELD = 1, A_FIELD = 2, X_FIELD = 3


  ! for internal use
  integer::nr,nz,ns,mmaxe,nmax,mmaxpm,nsp
  integer,dimension(:),allocatable::nn,mmaxp,m1,m2
  real*8 ::raxis,zaxis
  real*8,dimension(:),  allocatable::cs,csm
  real*8,dimension(:,:),allocatable::r,z,s,chi
  complex*16,dimension(:,:),allocatable::rmi,zmi,rmm,zmm
  complex*16,dimension(:,:,:),allocatable::b1,b2,b3,x1
  real*8 :: I_scale = 1.d0


  public :: &
     marsf_load, &
     marsf_broadcast, &
     marsf_get_Bf, &
     marsf_close

  contains
  !=====================================================================


  !=====================================================================
  subroutine marsf_load(iu, iconfig)
  use run_control, only: Prefix
  integer, intent(in)  :: iu
  integer, intent(out) :: iconfig

  character(len=256) :: geo_file, data_file
  integer :: ierr

  namelist /MARSF_Input/ &
     geo_file, data_file, I_scale


  rewind (iu)
  read   (iu, MARSF_Input, end=1000)
  iconfig = 1
  write (6, *)
  write (6, 1001)
  geo_file  = trim(Prefix)//geo_file
  data_file = trim(Prefix)//data_file
  call load_module(geo_file, data_file, ierr)
  if (ierr /= 0) then
     iconfig = 0
     return
  endif

  return
 1000 iconfig = 0
 1001 format(3x,'- Loading MARS-F data')
  end subroutine marsf_load
  !=====================================================================



  !=====================================================================
  subroutine marsf_broadcast()
  use parallel


  call broadcast_inte_s(nr)
  call broadcast_inte_s(nz)
  call broadcast_inte_s(ns)
  call broadcast_inte_s(mmaxe)
  call broadcast_inte_s(nmax)
  call broadcast_inte_s(mmaxpm)
  call broadcast_inte_s(nsp)
  call broadcast_real_s(raxis)
  call broadcast_real_s(zaxis)
  call broadcast_real_s(I_scale)

  if (mype > 0) then
     allocate (r(nr,nz), z(nr,nz), s(nr,nz), chi(nr,nz))
     allocate (nn(nmax),mmaxp(nmax),m1(nmax),m2(nmax))
     allocate (cs(ns),csm(ns))
     allocate (rmi(ns,mmaxe),zmi(ns,mmaxe),rmm(ns,mmaxe),zmm(ns,mmaxe))
     allocate (b1(ns,mmaxpm,nmax),b2(ns,mmaxpm,nmax),b3(ns,mmaxpm,nmax))
     allocate (x1(nsp+1,mmaxpm,nmax))
  endif

  call broadcast_real(r, nr*nz)
  call broadcast_real(z, nr*nz)
  call broadcast_real(s, nr*nz)
  call broadcast_real(chi, nr*nz)
  call broadcast_inte(nn, nmax)
  call broadcast_inte(mmaxp, nmax)
  call broadcast_inte(m1, nmax)
  call broadcast_inte(m2, nmax)
  call broadcast_real(cs, ns)
  call broadcast_real(csm, ns)
  call broadcast_complex(rmi, ns*mmaxe)
  call broadcast_complex(zmi, ns*mmaxe)
  call broadcast_complex(rmm, ns*mmaxe)
  call broadcast_complex(zmm, ns*mmaxe)
  call broadcast_complex(b1, ns*mmaxpm*nmax)
  call broadcast_complex(b2, ns*mmaxpm*nmax)
  call broadcast_complex(b3, ns*mmaxpm*nmax)
  call broadcast_complex(x1, (nsp+1)*mmaxpm*nmax)

  end subroutine marsf_broadcast
  !=====================================================================



  !=====================================================================
  subroutine load_module(geo_file, data_file, ierr)
  character(len=*), intent(in)  :: geo_file, data_file
  integer,          intent(out) :: ierr

  integer, parameter :: iu = 99

  real(real64) :: rtmp(8)
  integer :: ir, iz, i, j, k, n, nsv


  ierr = 0
  ! read geometry
  open  (iu, file=geo_file, status='old', form='formatted')
  read  (iu, *) nr, nz, raxis, zaxis
  allocate (r(nr,nz), z(nr,nz), s(nr,nz), chi(nr,nz))
  do ir=1,nr
  do iz=1,nz
     read  (iu, *) r(ir,iz), z(ir,iz), s(ir,iz), chi(ir,iz)
  enddo
  enddo
  close (iu)


  nmax = 1
  allocate(nn(nmax),mmaxp(nmax),m1(nmax),m2(nmax))
  ! first read B-field data to get dimentsions 
  do n=1,nmax
  open  (iu, file=data_file, status='old', form='formatted')
  read  (iu, *) nn(n),mmaxe,m1(n),m2(n),nsp,nsv
  close (iu)

  mmaxp(n) = m2(n)-m1(n) + 1
  ns       = nsp + nsv
  enddo

  ! find max(mmaxp)
  mmaxpm = mmaxp(1)
  do n=1,nmax
     if (mmaxpm.lt.mmaxp(n)) mmaxpm = mmaxp(n)
  enddo

  ! allocate new arrays
  allocate(cs(ns),csm(ns))
  allocate(rmi(ns,mmaxe),zmi(ns,mmaxe),rmm(ns,mmaxe),zmm(ns,mmaxe))
  allocate(b1(ns,mmaxpm,nmax),b2(ns,mmaxpm,nmax),b3(ns,mmaxpm,nmax))
  allocate(x1(nsp+1,mmaxpm,nmax))

  ! read B-field data for all toroidal harmonics
  ! followed by X1 data for normal displacement of the plasma      
  do n=1,nmax
  open  (iu, file=data_file, status='old', form='formatted')
  read  (iu, *) nn(n),mmaxe,m1(n),m2(n),nsp,nsv
  do i=1,ns
     read  (iu, *) cs(i),csm(i),rtmp(1)
  enddo

  do j=1,mmaxe
  do i=1,ns
     read  (iu, *) (rtmp(k),k=1,8)
     rmi(i,j) = cmplx(rtmp(1),rtmp(2))
     zmi(i,j) = cmplx(rtmp(3),rtmp(4))
     rmm(i,j) = cmplx(rtmp(5),rtmp(6))
     zmm(i,j) = cmplx(rtmp(7),rtmp(8))
  enddo
  enddo

  do j=1,mmaxp(n)
  do i=1,ns
     read  (iu, *) (rtmp(k),k=1,6)
     b1(i,j,n) = cmplx(rtmp(1),rtmp(2))
     b2(i,j,n) = cmplx(rtmp(3),rtmp(4))
     b3(i,j,n) = cmplx(rtmp(5),rtmp(6))
  enddo
  enddo

  do j=1,mmaxp(n)
  do i=1,nsp+1
     read  (iu, *) (rtmp(k),k=1,2)
     x1(i,j,n) = cmplx(rtmp(1),rtmp(2))
  enddo
  enddo

  close (iu)
  enddo

  end subroutine load_module
  !=====================================================================



  !=====================================================================
      subroutine fio_eval_field_f(type_field,xrphz,brphz,ierr)
      implicit none
      integer::type_field,ierr
      integer::ir,iz,i,j,k,n
      real*8 ::r0,r1,r2,z0,z1,z2,s0,chi0,phi0
      real*8 ::drds,dzds,drdc,dzdc,jac,sg22,h1,h2,h3,h4,h5,d
      real*8 ::chi00,chi01,chi10,chi11
      real*8,dimension(*)::xrphz,brphz
      complex*16::br,bp,bz
      complex*16::y1,y2,y3,y0,yp,ctmp,b10,b20,b30,x10

      ierr = 0

!     for given (r0,z0) find corresponding (s0,chi0)
      r0   = xrphz(1)
      phi0 = xrphz(2)
      z0   = xrphz(3)
      

      k  = 0
      do ir=1,nr-1
         if ( (r(ir,1)-r0)*(r(ir+1,1)-r0).le.0. ) k=ir
      enddo
      
      if (k.eq.0) then
         ierr = 1
         !write(*,*) 'fio_eval_field_f: ierr=', ierr
         !stop
         return
      endif

      r1 = r(k,1)
      r2 = r(k+1,1)
      ir = k

      k  = 0
      do iz=1,nz-1
         if ( (z(1,iz)-z0)*(z(1,iz+1)-z0).le.0. ) k=iz
      enddo
      
      if (k.eq.0) then
         ierr = 2
         !write(*,*) 'fio_eval_field_f: ierr=', ierr
         !stop
         return
      endif

      z1 = z(1,k)
      z2 = z(1,k+1)
      iz = k

      s0 = ( s(ir,iz)*    (r2-r0)*(z2-z0) + &
             s(ir,iz+1)*  (r2-r0)*(z0-z1) + &
             s(ir+1,iz)*  (r0-r1)*(z2-z0) + &
             s(ir+1,iz+1)*(r0-r1)*(z0-z1) )/(r2-r1)/(z2-z1)

      chi00 = chi(ir,iz)
      chi01 = chi(ir,iz+1)
      chi10 = chi(ir+1,iz)
      chi11 = chi(ir+1,iz+1)

      if ( (max(chi00,chi01,chi10,chi11) - &
            min(chi00,chi01,chi10,chi11)).lt.pi ) then 
          
         chi0 = ( chi00*(r2-r0)*(z2-z0) + &
                  chi01*(r2-r0)*(z0-z1) + &
                  chi10*(r0-r1)*(z2-z0) + &
                  chi11*(r0-r1)*(z0-z1) )/(r2-r1)/(z2-z1)
      else
         ctmp = ( exp((0.,1.)*chi00)*(r2-r0)*(z2-z0) + &
                  exp((0.,1.)*chi01)*(r2-r0)*(z0-z1) + &
                  exp((0.,1.)*chi10)*(r0-r1)*(z2-z0) + &
                  exp((0.,1.)*chi11)*(r0-r1)*(z0-z1) )/(r2-r1)/(z2-z1)
         chi0 = datan2(imag(ctmp),real(ctmp))
         if (chi0.lt.0.) chi0 = chi0 + 2.*pi
      endif

      if (s0.lt.cs(2)) then
         chi0 = datan2(z0-zaxis,r0-raxis)
         if (chi0.lt.0.) chi0 = chi0 + 2.*pi
      endif

!     on the s-mesh, compute dR/ds,dZ/ds,dR/dchi,dZ/dchi, at (s0,chi0)
!     using 3 points: 2 integer  points and 1 middle point
      k = 0
      do i=1,ns-1
         if ( (cs(i)-s0)*(cs(i+1)-s0).le.0. ) k=i
      enddo
      if (k.eq.0) then
         ierr = 3
         !write(*,*) 'fio_eval_field_f: ierr=', ierr
         !write(*,*) 's0=',s0
         !stop
         return
      endif

      i    = k
      h1   = cs(i)   - s0
      h2   = csm(i)  - s0
      h3   = cs(i+1) - s0
      d    = (h1-h2)*(h2-h3)*(h3-h1)
      drds = 0.
      dzds = 0.
      drdc = 0.
      dzdc = 0.

      do j=1,mmaxe
         y1=rmi(i,j)
         y2=rmm(i,j)
         y3=rmi(i+1,j)
         y0=(h3*h2*(h3-h2)*y1+h1*h3*(h1-h3)*y2+h2*h1*(h2-h1)*y3)/d
         yp=((h2-h3)*(h2+h3)*y1+(h3-h1)*(h3+h1)*y2+(h1-h2)*(h1+h2)*y3)/d
         
         if (j.eq.1) then
            drds = drds + real(yp)
         else
            ctmp = exp((0.,1.)*(j-1)*chi0)
            drds = drds + 2.*real(yp*ctmp)
            drdc = drdc + 2.*real(y0*(0.,1.)*(j-1)*ctmp)
         endif

         y1=zmi(i,j)
         y2=zmm(i,j)
         y3=zmi(i+1,j)
         y0=(h3*h2*(h3-h2)*y1+h1*h3*(h1-h3)*y2+h2*h1*(h2-h1)*y3)/d
         yp=((h2-h3)*(h2+h3)*y1+(h3-h1)*(h3+h1)*y2+(h1-h2)*(h1+h2)*y3)/d
         
         if (j.eq.1) then
            dzds = dzds + real(yp)
         else
            dzds = dzds + 2.*real(yp*ctmp)
            dzdc = dzdc + 2.*real(y0*(0.,1.)*(j-1)*ctmp)
         endif
      enddo

!     compute jacobian
      jac = (drds*dzdc - drdc*dzds)*r0
      if (type_field.eq.3) sg22=sqrt(drdc**2+dzdc**2)

!     compute b10,b20,b30,x10 at (s0,chi0)
!     b30 is computed only if type_field=1
!     x10 is computed only if type_field=3
!     go through all toroidal harmonics
      if (s0.lt.1..and.i.gt.1) i=i-1 
      h1   = cs(i)    - s0
      h2   = cs(i+1)  - s0
      h3   = cs(i+2)  - s0
      h4   = csm(i)   - s0
      h5   = csm(i+1) - s0
      d    = (h1-h2)*(h2-h3)*(h3-h1)

      brphz(1:3) = 0.

      do n=1,nmax

      b10 = (0.,0.)
      b20 = (0.,0.)   
      b30 = (0.,0.)   
      x10 = (0.,0.)   
      do j=1,mmaxp(n)
         ctmp = exp((0.,1.)*(j-1+m1(n))*chi0)

         if (type_field.eq.1.or.type_field.eq.2) then
            y1 = b1(i,j,n)
            y2 = b1(i+1,j,n)
            y3 = b1(i+2,j,n)
            y0=(h3*h2*(h3-h2)*y1+h1*h3*(h1-h3)*y2+h2*h1*(h2-h1)*y3)/d
            b10 = b10 + y0*ctmp

            y1 = (h5*b2(i,j,n)-h4*b2(i+1,j,n))/(h5-h4)
            b20 = b20 + y1*ctmp
         endif   

         if (type_field.eq.1) then
            y1 = (h5*b3(i,j,n)-h4*b3(i+1,j,n))/(h5-h4)
            b30 = b30 + y1*ctmp
         endif

         if (type_field.eq.3.and.s0.le.1.) then
            y1 = x1(i,j,n)
            y2 = x1(i+1,j,n)
            y3 = x1(i+2,j,n)
            y0=(h3*h2*(h3-h2)*y1+h1*h3*(h1-h3)*y2+h2*h1*(h2-h1)*y3)/d
            x10 = x10 + y0*ctmp
         endif   
      enddo
      
!     compute toroidal harmonics of B- or A-field at (s0,chi0)
      if (type_field.eq.1) then     !B-field 
         br = (b10*drds+b20*drdc)/jac
         bp =-b30*r0/jac
         bz = (b10*dzds+b20*dzdc)/jac
      elseif (type_field.eq.2) then !A-field 
         br = (b10*dzds+b20*dzdc)*r0/jac/((0.,1.)*nn(n))
         bp = (0.,0.)
         bz = (b10*drds+b20*drdc)*r0/jac/(-(0.,1.)*nn(n))   
      elseif (type_field.eq.3) then !xi_n
         br = x10*jac*sg22/r0
         bp = (0.,0.)
         bz = (0.,0.)
      else
         br = (0.,0.)
         bp = (0.,0.)
         bz = (0.,0.)
      endif

!     compute brphz at real space point xrphz
      ctmp = exp(-(0.,1.)*nn(n)*phi0)


      if (type_field==1) then
         brphz(1) = brphz(1) + 2.*real(br*ctmp)
         brphz(2) = brphz(2) + 2.*real(bp*ctmp)
         brphz(3) = brphz(3) + 2.*real(bz*ctmp)
      elseif (type_field==2) then
         brphz(1) = brphz(1) + 2.*real(br*ctmp)
         brphz(3) = brphz(3) + 2.*real(bz*ctmp)
      elseif (type_field==3) then
         brphz(1) = brphz(1) + 2.*real(br*ctmp)
      endif

      enddo
   
      return
      end subroutine fio_eval_field_f
  !=====================================================================



  !=====================================================================
  function marsf_get_Bf(r) result(Bf)
  use numerics, only: cm_to_m
  real(real64), intent(in) :: r(3)
  real(real64)             :: Bf(3)

  real(real64) :: x(3), Bx(3)
  integer      :: ierr


  Bf   = 0.d0

  ! expected coordinate order: R[m], Phi[rad], Z[m]
  x(1) = r(1) * cm_to_m
  x(2) = r(3)
  x(3) = r(2) * cm_to_m
  call fio_eval_field_f(B_FIELD, x, Bx, ierr)
  if (ierr /= 0) return
  Bx    = Bx * I_scale * 1.d4! Tesla -> Gauss
  Bf(1) = Bx(1)
  Bf(2) = Bx(3)
  Bf(3) = Bx(2)

  end function marsf_get_Bf
  !=====================================================================



  !=====================================================================
  subroutine marsf_close()


  deallocate(r,z,s,chi)
  deallocate(cs,csm)
  deallocate(rmi,zmi,rmm,zmm)
  deallocate(nn,mmaxp,m1,m2)
  deallocate(b1,b2,b3,x1)

  end subroutine marsf_close
  !=====================================================================

end module marsf
