module quad_ele
  use iso_fortran_env
  implicit none
  private

  type, public :: t_quad_ele
     ! n_sym = toroidal symmetry, coordinates should be in the range 0 <= phi <= 2*pi / n_sym
     integer :: n_phi, n_RZ, n_sym
     real*8  :: dR, dZ
     real*8, dimension(:), allocatable   :: phi    ! dimension(0:n_phi)
     real*8, dimension(:,:), allocatable :: R, Z   ! dimension(0:n_phi, 0:n_RZ)

     ! working array
     real*8, dimension(:,:,:), allocatable :: cA, cB, cC, cD    ! dimension(nphi,nRZ,2)

     contains
     procedure :: new
     procedure :: load      => quad_ele_load
     procedure :: plot      => quad_ele_plot
     procedure :: plot_at   => quad_ele_plot_at
     procedure :: intersect => quad_ele_intersect
     procedure :: sample
     procedure :: sample_phi
     procedure :: slice
     procedure :: left_hand_shift
     procedure :: destroy
     procedure :: get_stellarator_symmetric_element
     procedure :: setup_coefficients
     procedure :: element
  end type t_quad_ele

  contains
!=======================================================================



!=======================================================================
  subroutine new(this, n_phi, n_RZ, n_sym)
  class(t_quad_ele)   :: this
  integer, intent(in) :: n_phi, n_RZ,  n_sym


  call this%destroy()
  this%n_phi = n_phi
  this%n_RZ  = n_RZ
  this%n_sym = n_sym
  allocate (this%phi(0:n_phi))
  allocate (this%R(0:n_phi, 0:n_RZ), this%Z(0:n_phi, 0:n_RZ))

  end subroutine new
!=======================================================================



!=======================================================================
  subroutine quad_ele_load(this, filename, title)
  use math
  class(t_quad_ele)         :: this
  character(len=*), intent(in) :: filename
  character*80, intent(out), optional :: title

  integer, parameter :: iu = 32

  real*8, dimension(:,:), allocatable :: Rtmp, Ztmp
  real*8, dimension(:), allocatable   :: phitmp
  integer, dimension(:), allocatable  :: Di
  character*80 :: s
  integer :: i, j, n, m


  ! read data from file
  open  (iu, file=filename)
  read  (iu, '(a)') s
  if (present(title)) title = s

  read  (iu, *) n, m, this%n_sym, this%dR, this%dZ
  n = n-1
  m = m-1
  this%n_phi = n
  this%n_RZ  = m
  allocate (this%phi(0:n))
  allocate (this%R(0:n, 0:m), this%Z(0:n, 0:m))
  do i=0,n
     read (iu, *) this%phi(i)
     do j=0,m
        read (iu, *) this%R(i, j), this%Z(i, j)
     enddo
  enddo
  close (iu)

  ! convert deg to rad
  this%phi = this%phi / 180.d0 * pi


  ! sanity check
  allocate (Di(n))
  do i=1,n
     Di = int(sign(1.d0, this%phi(i)-this%phi(i-1)))
  enddo
  if (abs(sum(Di)) .ne. n) then
     write (6, *) 'error: elements must be in increasing or decreasing order!'
     stop
  endif


  ! reverse order if necessary
  if (sum(di) < 0) then
     allocate (Rtmp(0:n, 0:m), Ztmp(0:n, 0:m))
     allocate (phitmp(0:n))
     do i=0,n
        Rtmp(n-i,:) = this%R(i,:)
        Ztmp(n-i,:) = this%Z(i,:)
        phitmp(n-i) = this%phi(i)
     enddo
     this%R   = Rtmp
     this%Z   = Ztmp
     this%phi = phitmp
     deallocate (Rtmp, Ztmp, phitmp)
  endif
  deallocate (Di)

  call this%setup_coefficients()

  end subroutine quad_ele_load
!=======================================================================



!=======================================================================
  subroutine setup_coefficients(this)
  class(t_quad_ele)         :: this

  integer :: i, j, n, m


  n = this%n_phi
  m = this%n_RZ
  allocate (this%cA(n,m,2), this%cB(n,m,2), this%cC(n,m,2), this%cD(n,m,2))

  ! setup coefficients
  do i=1,n
  do j=1,m
     this%cA(i,j,1) = 0.5d0 * (this%R(i-1,j-1) + this%R(i-1,j))
     this%cA(i,j,2) = 0.5d0 * (this%Z(i-1,j-1) + this%Z(i-1,j))

     this%cB(i,j,1) = 0.5d0 * (- this%R(i-1,j-1) + this%R(i-1,j))
     this%cB(i,j,2) = 0.5d0 * (- this%Z(i-1,j-1) + this%Z(i-1,j))

     this%cC(i,j,1) = 0.5d0 * (- this%R(i-1,j-1) + this%R(i,j-1) + this%R(i,j) - this%R(i-1,j))
     this%cC(i,j,2) = 0.5d0 * (- this%Z(i-1,j-1) + this%Z(i,j-1) + this%Z(i,j) - this%Z(i-1,j))

     this%cD(i,j,1) = 0.5d0 * (this%R(i-1,j-1) - this%R(i,j-1) + this%R(i,j) - this%R(i-1,j))
     this%cD(i,j,2) = 0.5d0 * (this%Z(i-1,j-1) - this%Z(i,j-1) + this%Z(i,j) - this%Z(i-1,j))
  enddo
  enddo

  end subroutine setup_coefficients
!=======================================================================



!=======================================================================
! Plot surface mesh
!=======================================================================
  subroutine quad_ele_plot(this, filename, Output_Format)
  class(t_quad_ele)         :: this
  character*120, intent(in) :: filename
  integer, intent(in)       :: Output_Format

  integer, parameter :: iu = 99

  integer :: i, j


  open  (iu, file=filename)
  if (Output_Format == 1) then
  ! plot slices at phi = const
  do i=0,this%n_phi
     do j=0,this%n_RZ
        write (iu, *) this%R(i, j), this%Z(i, j), this%phi(i)
     enddo
     write (iu, *)
  enddo

  elseif (Output_Format == 2) then
  ! plot toroidal profiles
  do j=0,this%n_RZ
     do i=0,this%n_phi
        write (iu, *) this%R(i, j), this%Z(i, j), this%phi(i)
     enddo
     write (iu, *)
  enddo
  endif
  close (iu)

  end subroutine quad_ele_plot
!=======================================================================



!=======================================================================
  function slice(this, phi0) result(C)
  use math
  use search
  use curve2D
  class(t_quad_ele)         :: this
  real(real64), intent(in)  :: phi0
  type(t_curve)             :: C

  real(real64) :: t, R, Z, dl, phi
  integer      :: i, ind, ierr


  ! check boundaries
  phi = phi0
  do
     if (phi >= this%phi(0)) exit

     phi = phi + pi2 / this%n_sym
  enddo

  do
     if (phi <= this%phi(this%n_phi)) exit

     phi = phi - pi2 / this%n_sym
  enddo

  if (phi.lt.this%phi(0) .or. phi.gt.this%phi(this%n_phi)) then
     return
  endif


  ! find index "ind" with phi(ind) <= phi0 <= phi(ind+1)
  ind = binary_interval_search (0, this%n_phi, this%phi, phi, ierr)
  if (ierr .ne. 0) then
     write (6, *) 'error: boundary check passed but could not find position!'
     write (6, *) ind, ierr
     stop
  endif
  t = (phi - this%phi(ind)) / (this%phi(ind+1) - this%phi(ind))


  ! sample slice at phi0 (ind, t)
  call C%new(this%n_RZ)


  do i=0,this%n_RZ
     R = this%R(ind,i) + t * (this%R(ind+1,i)-this%R(ind,i))
     Z = this%Z(ind,i) + t * (this%Z(ind+1,i)-this%Z(ind,i))
     C%x(i,1) = R
     C%x(i,2) = Z
  enddo

  ! check if slice is closed
  dl = sqrt((Z-C%x(0,2))**2 + (R-C%x(0,1))**2)
  if (dl < epsilon(real(1.0,real64))) C%closed = .true.

  end function slice
!=======================================================================



!=======================================================================
! Plot profile at toroidal location phi [rad]
!=======================================================================
  subroutine quad_ele_plot_at(this, phi, filename)
  use math
  use search
  use curve2D
  class(t_quad_ele)            :: this
  real(real64), intent(in)     :: phi
  character(len=*), intent(in) :: filename

  integer, parameter :: iu = 99

  type(t_curve) :: C


  ! get slice at phi0
  C = this%slice(phi)
  if (C%n_seg < 0) return

  ! write profile at phi0 (ind, t)
  open  (iu, file=filename)
  write (iu, 1000) phi / pi * 180.d0
  call C%plot(iu=iu)
  close (iu)

 1000 format ("# phi [deg] = ",e18.10)
  end subroutine quad_ele_plot_at
!=======================================================================



!=======================================================================
  function quad_ele_intersect(this, r1, r2, X, nelem, tau) result(l)
  use math
  use search
  class(t_quad_ele)   :: this
  real*8, intent(in)  :: r1(3), r2(3)
  real*8, intent(out) :: X(3)
  integer, intent(out), optional :: nelem
  real(real64), intent(out), optional :: tau
  logical             :: l

  logical :: ljump
  real*8  :: r1s(3), r2s(3), rAs(3), rBs(3), Dphi, xi
  integer :: k, n, Di, istat


  ! set position in 0 <= phi <= 2*pi/n_sym
  n      = this%n_sym
  r1s    = r1
  r1s(3) = phi_sym(r1s(3), n)
  r2s    = r2
  r2s(3) = phi_sym(r2s(3), n)

  if (present(nelem)) nelem = -1


  ! jump beyond symmetry domain? ASSUMPTION: |phi2 - phi1| < pi/n_sym
  ljump = .false.
  Dphi  = r2s(3) - r1s(3)
  if (abs(Dphi) > pi/n) then
     Dphi     = Dphi - sign(pi2/n, Dphi)
     Di       = int(sign(1.d0, Dphi))
     rBs(3)   = pi2/n * (Di+1)/2
     rAs(3)   = pi2/n * (1-Di)/2
     xi       = (rBs(3) - r1s(3)) / Dphi
     rAs(1:2) = r1s(1:2) + xi * (r2s(1:2)-r1s(1:2))
     rBs(1:2) = rAs(1:2)
     ljump    = .true.
#if defined(DEBUG)
     write (6, *) 'ljump = ', ljump
     write (6, *) 'Di    = ', Di
     write (6, *) 'rAs   = ', rAs
     write (6, *) 'rBs   = ', rBs
#endif
  endif


  ! now check for intersections
  if (ljump) then
     call check_intersection (r1s, rBs, xi, istat)
     if (istat.ne.1) call check_intersection (rAs, r2s, xi, istat)
     ! TODO: update xi?
  else
     call check_intersection (r1s, r2s, xi, istat)
  endif


  l = .false.
  if (present(tau)) tau = -1.d0
  if (istat.eq.1) then
     l = .true.
     X = r1 + (r2-r1) * xi
     if (present(tau)) tau = xi
     if (nelem < 0) then
        write (6, *) 'nelem = ', nelem
        write (6, *) 'ljump = ', ljump
        stop
     endif
  endif

  contains
  !---------------------------------------------------------------------
  subroutine check_intersection(r1, r2, t, istat)
  real*8, intent(in)   :: r1(3), r2(3)
  real*8, intent(out)  :: t
  integer, intent(out) :: istat

  real*8  :: Dphi, phil, phir, phi1, phi2
  real*8  :: xA(2), xB(2), xC(2), xD(2), x5(2), x6(2), y1(2), y2(2), tc(2), ts
  integer :: Di, i, iA, iB, j, n, is, js


  ! 0. get search direction
  phi1  = r1(3)
  phi2  = r2(3)
  Dphi  = phi2 - phi1
  istat = 0
  if (Dphi == 0.d0) return
  !if (abs(Dphi) > pi2/this%n_sym) Dphi = Dphi - sign(pi2/this%n_sym, Dphi)
  Di   = int(sign(1.d0,Dphi))
#if defined(DEBUG)
  if (Di > 0) then
     write (6, *) 'Forward search'
  else
     write (6, *) 'Backward search'
  endif
#endif


  ! 1. mark slices which need to be checked for intersections
  n          = this%n_phi
  phil       = this%phi(0)
  phir       = this%phi(n)

  ! 1.1. find start and end indices in mesh domain
  iA = binary_interval_search (0, this%n_phi, this%phi, phi1, istat)
  iB = binary_interval_search (0, this%n_phi, this%phi, phi2, istat)

  ! 1.2 check if mesh is somewhere in between phi1 and phi2
  ! left boundary in forward search domain (phi2 > phil > phi1)?
  if (phi2 > phil .and. phil > phi1) then
     iA = 0
     ! right boundary within search domain?
     if (phi2 > phir) then
        iB = n-1
     endif
  endif
  ! right boundary in forward search domain (phi2 > phir > phi1)?
  if (phi2 > phir .and. phir > phi1) then
     iB = n-1
     ! left boundary witin search domain
     if (phi1 < phil) then
        iA = 0
     endif
  endif

  ! left boundary in backward search domain (phi2 < phil < phi1)?
  if (phi2 < phil .and. phil < phi1) then
     iB = 0
     ! right boundary within search domain?
     if (phi1 > phir) then
        iA = n-1
     endif
  endif
  ! right boundary in backward search domain (phi2 < phir < phi1)?
  if (phi2 < phir .and. phir < phi1) then
     iA = n-1
     ! left boundary within search domain?
     if (phi2 < phil) then
        iB = 0
     endif
  endif

  ! return if no overlap is found
  if (iA.eq.-1 .or. iB.eq.-1) then
#if defined(DEBUG)
     ! sanity check
     if (iA*iB.eq.-1) then
        write (6, *) 'iA, iB = ', iA, iB
        write (6, *) 'this should not happen!'
        stop
     endif
#endif
     istat = 0
     return
  endif



  ! check slice [phi(i), phi(i+1)] for intersections
  istat = 0
  t     = 2.d0
  is    = -1
  js    = -1
  do i=iA+1,iB+1,Di
     do j=1,this%n_RZ
        x5 = r1(1:2) + (r2(1:2)-r1(1:2)) * (this%phi(i-1)-phi1) / Dphi
        x6 =           (r2(1:2)-r1(1:2)) * (this%phi(i)-this%phi(i-1)) / Dphi
        xA = this%cA(i,j,:) - x5
        xB = this%cB(i,j,:)
        xC = this%cC(i,j,:) - x6
        xD = this%cD(i,j,:)
        call solve_bilinear_system_bc (xA, xB, xC, xD, y1, y2, istat)

        ! for each intersection with the surface element check if ts in [0,1]
        tc(1) = y1(2)
        tc(2) = y2(2)
        do k=1,istat
           ts = (this%phi(i-1)-phi1) / Dphi + (this%phi(i)-this%phi(i-1)) / Dphi * tc(k)

           ! update intersection if new solution is closer to t=0
           if (ts < t .and. ts.ge.0.d0) then
              t = ts
              js = j
              is = i
           endif
        enddo
     enddo

     ! exit if an intersection has been found in this slice
     if (t < 2.d0) then
        istat = 1
        if (present(nelem)) nelem = (js-1) + (is-1)*this%n_RZ
        if (nelem < 0) then
           write (6, *) 'nelem1 = ', nelem, js, is, this%n_RZ
        endif
        return
     endif
  enddo


  ! no intersections found
  istat = 0
  if (istat > 0) then
     write (6, *) 'istat = ', istat
     write (6, *) 'nelem = ', nelem
     write (6, *) 't     = ', t
  endif


  end subroutine check_intersection
  !---------------------------------------------------------------------
  end function quad_ele_intersect
!=======================================================================



!=======================================================================
! Return stellarator symmetric element with coordinates:
! phi -> 2*pi/n_sym - phi
! Z   -> - Z
!=======================================================================
  function get_stellarator_symmetric_element(this) result(S)
  use math
  class(t_quad_ele)   :: this
  type(t_quad_ele)    :: S

  integer :: i, n, m

  S%n_phi = this%n_phi
  S%n_RZ  = this%n_RZ
  S%n_sym = this%n_sym
  S%dR    = this%dR
  S%dZ    = this%dZ

  n = S%n_phi
  m = S%n_RZ
  allocate (S%phi(0:n))
  allocate (S%R(0:n, 0:m), S%Z(0:n, 0:m))

  do i=0,n
     S%phi(i) = pi2/S%n_sym - this%phi(n-i)
     S%R(i,:) =   this%R(n-i,:)
     S%Z(i,:) = - this%Z(n-i,:)
  enddo

  call S%setup_coefficients()

  end function get_stellarator_symmetric_element
!=======================================================================



!=======================================================================
! Based on module curve2D
!=======================================================================
  subroutine left_hand_shift(this, dl)
  use curve2D
  class(t_quad_ele)        :: this
  real(real64), intent(in) :: dl

  type(t_curve) :: C
  integer :: i, n


  n = this%n_RZ + 1
  do i=0,this%n_phi
     call make_2D_curve (n, this%R(i,:), this%Z(i,:), C)
     call C%left_hand_shift(dl)

     this%R(i,:) = C%x(:,1)
     this%Z(i,:) = C%x(:,2)
  enddo

  end subroutine left_hand_shift
!=======================================================================



!=======================================================================
  function sample_phi(this, tau) result(phi)
  class(t_quad_ele)        :: this
  real(real64), intent(in) :: tau
  real(real64)             :: phi


  phi  = this%phi(0) + tau*(this%phi(this%n_phi) - this%phi(0))

  end function sample_phi
!=======================================================================
! sample surface position at relative coordinates (xi, tau)
!=======================================================================
  function sample(this, tau, xi) result(r)
  use search
  class(t_quad_ele)        :: this
  real(real64), intent(in) :: tau, xi
  real(real64)             :: r(3)

  real(real64) :: tau1, xi1
  integer      :: istat, i, j


  r(3) = this%phi(0) + tau*(this%phi(this%n_phi) - this%phi(0))
  i    = binary_interval_search (0, this%n_phi, this%phi, r(3), istat)

  tau1 = (r(3) - this%phi(i)) / (this%phi(i+1) - this%phi(i))
  xi1  = xi * this%n_RZ
  j    = int(xi1)
  xi1  = xi1 - j
  if (j == this%n_RZ) then
     j   = this%n_RZ - 1
     xi1 = 1.d0
  endif

  ! map [0,1] to [-1,1]
  xi1 = 2.d0 * xi1 - 1.d0
  j   = j + 1
  i   = i + 1


  r(1:2) = this%cA(i,j,:) + xi1*this%cB(i,j,:) + tau1*this%cC(i,j,:) + xi1*tau1*this%cD(i,j,:)


  end function sample
!=======================================================================



!=======================================================================
  subroutine destroy(this)
  class(t_quad_ele)         :: this


  if (allocated(this%cA)) deallocate (this%cA, this%cB, this%cC, this%cD)
  if (allocated(this%phi)) deallocate (this%phi, this%R, this%Z)

  end subroutine destroy
!=======================================================================



!=======================================================================
  subroutine element(this, ielem, i, j)
  class(t_quad_ele)    :: this
  integer, intent(in)  :: ielem
  integer, intent(out) :: i, j


  j = mod(ielem, this%n_RZ) + 1
  i = ielem / this%n_RZ + 1

  end subroutine element
!=======================================================================

end module quad_ele









!subroutine test_quad_ele
!  use quad_ele
!  use run_control, only: Output_File, Output_Format, Phi_Output
!
!  type(t_quad_ele) :: S
!  character*120 :: filename
!
!
!  filename = Output_File
!  call S%load(filename)
!  !call S%plot(Output_File, Output_Format)
!
!  filename = trim(Output_File)//'.plt'
!!  call S%plot(filename, 1)
!  call S%plot_at(Phi_output, filename)
!
!
!end subroutine test_quad_ele
