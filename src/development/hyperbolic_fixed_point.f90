! TODO
subroutine hyperbolic_fixed_point
  use iso_fortran_env
  use run_control, only: x_start, Trace_Step, Trace_Method
  use grid
  use fieldline
  use equilibrium
  use math
  implicit none

  integer, parameter :: iu = 72
  integer, parameter :: np = 10

  real(real64), dimension(np,np,3) :: x
  real(real64), dimension(np,np)   :: Delta

  type(t_fieldline)  :: F
  real(real64) :: Dphi, Xout(3), Dout, u0(2), v0(2), dl
  integer      :: n, m, ilevel

  n = 3
  m = 13
  Dphi       = pi2 * m / n
  Trace_Step = pi2 / 360.d0

  ilevel = 0
  dl     = 2.d0

  Xout  = x_start
  u0(1) = dl
  u0(2) = 0.d0
  v0(1) = 0.d0
  v0(2) = dl
  call evaluateX (Xout, u0, v0, Xout, Dout, ilevel)
  write (6, *) 'final approximation: ', Xout(1:2)
  write (6, *) 'quality parameter:   ', Dout


  contains

  recursive subroutine evaluateX (x0, u, v, xp, D, ilevel)
  real(real64), intent(in)  :: x0(3), u(2), v(2)
  real(real64), intent(out) :: xp(3), D
  integer, intent(inout)    :: ilevel

  character(len=80) :: filename
  real(real64) :: D2(2), X4, X3Y1, X2Y2, X1Y3, Y4, X3, X2Y1, X1Y2, Y3, X2, X1Y1, Y2, X1, Y1
  real(real64) :: DX2, DY2, DXY, DX, DY, D1
  real(real64) :: A(4,4), B(4), C(4), xa, ya, det
  real(real64) :: H(2,2), lambda(2), V1(2), V2(2)
  real(real64) :: xi, eta, xp1(3), theta
  integer      :: i, j, i0, j0, i1, j1


  write (6, *)
  write (6, *) 'ilevel = ', ilevel
  ! setup grid
  do j=1,np
  do i=1,np
     xi  = (-1.d0 + 2.d0*(i-1)/(np-1))
     eta = (-1.d0 + 2.d0*(j-1)/(np-1))

     x(i,j,1:2) = x0(1:2) + xi*u + eta*v
     x(i,j,3) = x0(3)
  enddo
  enddo

  open  (iu, file='box.dat', position='append')
  write (iu, *)
  write (iu, *) x( 1, 1,1:2)
  write (iu, *) x( 1,np,1:2)
  write (iu, *) x(np,np,1:2)
  write (iu, *) x(np, 1,1:2)
  write (iu, *) x( 1, 1,1:2)
  close (iu)

  D  = 1.d99
  i0 = 0
  j0 = 0
  i1 = 0
  j1 = 0
  xp1 = 0.d0
  ! trace field lines
  do j=1,np
  do i=1,np
     call F%init(x(i,j,:), Trace_Step, Trace_Method, FL_ANGLE)
     call F%trace(Dphi, .true.)
     D2         = F%rc(1:2) - x(i,j,1:2)
     Delta(i,j) = sqrt(sum(D2**2))

     if (Delta(i,j) < D) then
        i1  = i0
        j1  = j0
        xp1 = xp

        D  = Delta(i,j)
        i0 = i
        j0 = j
        xp = x(i,j,:)
     endif
  enddo
  enddo
  write (6, *) 'min. found at (i,j) = ', i0, j0
  write (6, *) 'min1. found at (i,j) = ', i1, j1
  write (6, *) '(R,     Z   ) = ', xp(1:2)
  theta = get_poloidal_angle(xp)/pi*180.d0
  if (theta < 0) theta = theta + 360.d0
  write (6, *) '(theta, PsiN) = ', theta, get_PsiN(xp)
  write (6, *) ' Delta        = ', D

  write (filename, '(i1)') ilevel
  filename = 'delta_'//trim(filename)//'.txt'
  open  (iu, file=filename)
  do j=1,np
     do i=1,np
        write (iu, *) x(i,j,1:2), Delta(i,j)
     enddo
     write (iu, *)
  enddo
  close (iu)





  xa = xp(1)
  ya = xp(2)
  x(:,:,1) = x(:,:,1) - xa
  x(:,:,2) = x(:,:,2) - ya

  ! evaluate
  X4   = sum(x(:,:,1)**4)
  X3Y1 = sum(x(:,:,1)**3  *  x(:,:,2)   )
  X2Y2 = sum(x(:,:,1)**2  *  x(:,:,2)**2)
  X1Y3 = sum(x(:,:,1)     *  x(:,:,2)**3)
    Y4 = sum(                x(:,:,2)**4)
  X3   = sum(x(:,:,1)**3)
  X2Y1 = sum(x(:,:,1)**2  *  x(:,:,2)   )
  X1Y2 = sum(x(:,:,1)     *  x(:,:,2)**2)
    Y3 = sum(                x(:,:,2)**3)
  X2   = sum(x(:,:,1)**2)
  X1Y1 = sum(x(:,:,1)     *  x(:,:,2)   )
    Y2 = sum(                x(:,:,2)**2)
  X1   = sum(x(:,:,1))
    Y1 = sum(                x(:,:,2)   )
  DX2  = sum(Delta(:,:)  *  x(:,:,1)**2)
  DY2  = sum(Delta(:,:)  *  x(:,:,2)**2)
  DXY  = sum(Delta(:,:)  *  x(:,:,1) * x(:,:,2))
  DX   = sum(Delta(:,:)  *  x(:,:,1))
  DY   = sum(Delta(:,:)  *  x(:,:,2))
  D1   = sum(Delta(:,:))

  A(1,1) = X4
  A(1,2) = X2Y2
  A(1,3) = X3Y1
  A(1,4) = X2

  A(2,1) = X2Y2
  A(2,2) =   Y4
  A(2,3) = X1Y3
  A(2,4) =   Y2

  A(3,1) = X3Y1
  A(3,2) = X1Y3
  A(3,3) = X2Y2
  A(3,4) = X1Y1

  A(4,1) = X2
  A(4,2) =   Y2
  A(4,3) = X1Y1
  A(4,4) = 1.d0

  B(1)   = DX2
  B(2)   = DY2
  B(3)   = DXY
  B(4)   = D1
  call elimi(4,C,A,B)
  !call linalglu(6,A,B,C)
  write (6, *) 'C = ', C

  ! Hessian
  H(1,1) = C(1)
  H(2,2) = C(2)
  H(1,2) = C(3)
  H(2,1) = C(3)
  lambda(1) = 0.5d0*(C(1)+C(2)) + sqrt(0.25d0*(C(1)+C(2))**2 - (C(1)*C(2) - C(3)**2))
  lambda(2) = 0.5d0*(C(1)+C(2)) - sqrt(0.25d0*(C(1)+C(2))**2 - (C(1)*C(2) - C(3)**2))
!  write (6, *) 'lambda = ', lambda

  ! eigenvectors
  V1(1) = 1.d0
  V1(2) = - (C(1) - lambda(1)) / C(2)
  !V1    = V1 / sqrt(sum(V1**2)) / sqrt(abs(lambda(1)))
  V1    = V1 / sqrt(sum(V1**2))

  V2(1) = 1.d0
  V2(2) = - (C(1) - lambda(2)) / C(2)
  !V2    = V2 / sqrt(sum(V2**2)) / sqrt(abs(lambda(2)))
  V2    = V2 / sqrt(sum(V2**2))

!  write (6, *) 'eigenvectors: '
!  write (6, *) 'v1 = ', v1
!  write (6, *) 'v2 = ', v2
!  open  (iu, file='eigenvectors.dat')
!  write (iu, *) xa, ya
!  write (iu, *) xa+V1(1), ya+V1(2)
!  write (iu, *)
!  write (iu, *) xa, ya
!  write (iu, *) xa+V2(1), ya+V2(2)
!  close (iu)


  ilevel = ilevel + 1
  !v1 = v1 * sqrt(sum(u**2)) * 4.d0/np * sqrt(abs(lambda(2)/lambda(1)))
  !v2 = v2 * sqrt(sum(v**2)) * 4.d0/np

  v1 = 0.5d0* (xp(1:2) - xp1(1:2))
  v2(1) =  0.5d0 * v1(2)
  v2(2) = -0.5d0 * v1(1)
  xp = 0.75d0*xp+ 0.25d0*xp1
  if (ilevel <= 4) call evaluateX (xp, v1, v2, xp, D, ilevel)

  end subroutine evaluateX

end subroutine hyperbolic_fixed_point






subroutine hyperbolic_fixed_point_v4
  use iso_fortran_env
  use run_control, only: x_start, Trace_Step, Trace_Method
  use grid
  use fieldline
  use equilibrium
  use math
  implicit none

  integer, parameter :: iu = 72
  integer, parameter :: np = 10

  real(real64), dimension(np,np,3) :: x
  real(real64), dimension(np,np)   :: Delta

  type(t_fieldline)  :: F
  real(real64) :: Dphi, Xout(3), Dout, u0(2), v0(2), dl
  integer      :: n, m, ilevel

  n = 3
  m = 13
  Dphi       = pi2 * m / n
  Trace_Step = pi2 / 360.d0

  ilevel = 0
  dl     = 2.d0
  !call evaluateX_grid (x_start, dl, 0.d0, Xout, Dout, ilevel)

  Xout  = x_start
  !u0(1) = 2.d0*dl/np
  u0(1) = dl
  u0(2) = 0.d0
  v0(1) = 0.d0
  !v0(2) = 2.d0*dl/np
  v0(2) = dl
  call evaluateX (Xout, u0, v0, 0.d0, Xout, Dout, ilevel)
  write (6, *) 'final approximation: ', Xout(1:2)
  write (6, *) 'quality parameter:   ', Dout


  contains

  recursive subroutine evaluateX_grid (x0, dl, alpha, xp, D, ilevel)
  real(real64), intent(in)  :: x0(3), dl, alpha
  real(real64), intent(out) :: xp(3), D
  integer, intent(inout)    :: ilevel

  real(real64) :: Delta(0:2)
  integer      :: i, j, i0, j0


  ! setup grid
  write (6, *)
  write (6, 1000) x0(1)-dl, x0(1)+dl, x0(2)-dl, x0(2)+dl
 1000 format ('box: ',f10.5,' -> ',f10.5,',  ',f10.5,' -> ',f10.5)
  do j=1,np
  do i=1,np
     x(i,j,1) = x0(1) + (-1.d0 + 2.d0*(i-1)/(np-1)) * dl
     x(i,j,2) = x0(2) + (-1.d0 + 2.d0*(j-1)/(np-1)) * dl
     x(i,j,3) = x0(3)
  enddo
  enddo


  D  = 1.d99
  i0 = 0
  j0 = 0
  ! trace field lines
  do j=1,np
  do i=1,np
     call F%init(x(i,j,:), Trace_Step, Trace_Method, FL_ANGLE)
     call F%trace(Dphi, .true.)
     Delta(1:2) = F%rc(1:2) - x(i,j,1:2)
     Delta(0)   = sqrt(sum(Delta(1:2)**2))

     !write (6, *) x(i,j,1:2), Delta(0)
     if (Delta(0) < D) then
        D  = Delta(0)
        i0 = i
        j0 = j
        xp = x(i,j,:)
     endif
  enddo
  enddo

  ! evaluate
  write (6, *) 'min. found at (i,j) = ', i0, j0
  write (6, *) '(R,     Z   ) = ', xp(1:2)
  write (6, *) '(theta, PsiN) = ', get_poloidal_angle(xp), get_PsiN(xp)
  write (6, *) ' Delta        = ', D

  !ilevel = ilevel + 1
  !if (ilevel < 8) call evaluateX_grid(xp, 2.d0*dl/np, alpha, xp, D, ilevel)

  end subroutine evaluateX_grid

  recursive subroutine evaluateX (x0, u, v, alpha, xp, D, ilevel)
  real(real64), intent(in)  :: x0(3), u(2), v(2), alpha
  real(real64), intent(out) :: xp(3), D
  integer, intent(inout)    :: ilevel

  real(real64) :: D2(2), X4, X3Y1, X2Y2, X1Y3, Y4, X3, X2Y1, X1Y2, Y3, X2, X1Y1, Y2, X1, Y1
  real(real64) :: DX2, DY2, DXY, DX, DY, D1
  real(real64) :: A(6,6), B(6), C(6), xa, ya, det
  real(real64) :: H(2,2), lambda(2), V1(2), V2(2)
  real(real64) :: xi, eta
  integer      :: i, j, i0, j0


  write (6, *)
  write (6, *) 'ilevel = ', ilevel
  ! setup grid
  do j=1,np
  do i=1,np
     xi  = (-1.d0 + 2.d0*(i-1)/(np-1))
     eta = (-1.d0 + 2.d0*(j-1)/(np-1))

     x(i,j,1:2) = x0(1:2) + xi*u + eta*v
     x(i,j,3) = x0(3)
  enddo
  enddo

  open  (iu, file='box.dat')
  write (iu, *) x( 1, 1,1:2)
  write (iu, *) x( 1,np,1:2)
  write (iu, *) x(np,np,1:2)
  write (iu, *) x(np, 1,1:2)
  write (iu, *) x( 1, 1,1:2)
  close (iu)

  D  = 1.d99
  i0 = 0
  j0 = 0
  ! trace field lines
  do j=1,np
  do i=1,np
     call F%init(x(i,j,:), Trace_Step, Trace_Method, FL_ANGLE)
     call F%trace(Dphi, .true.)
     D2         = F%rc(1:2) - x(i,j,1:2)
     Delta(i,j) = sqrt(sum(D2**2))

     if (Delta(i,j) < D) then
        D  = Delta(i,j)
        i0 = i
        j0 = j
        xp = x(i,j,:)
     endif
  enddo
  enddo

  ! evaluate
  X4   = sum(x(:,:,1)**4)
  X3Y1 = sum(x(:,:,1)**3  *  x(:,:,2)   )
  X2Y2 = sum(x(:,:,1)**2  *  x(:,:,2)**2)
  X1Y3 = sum(x(:,:,1)     *  x(:,:,2)**3)
    Y4 = sum(                x(:,:,2)**4)
  X3   = sum(x(:,:,1)**3)
  X2Y1 = sum(x(:,:,1)**2  *  x(:,:,2)   )
  X1Y2 = sum(x(:,:,1)     *  x(:,:,2)**2)
    Y3 = sum(                x(:,:,2)**3)
  X2   = sum(x(:,:,1)**2)
  X1Y1 = sum(x(:,:,1)     *  x(:,:,2)   )
    Y2 = sum(                x(:,:,2)**2)
  X1   = sum(x(:,:,1))
    Y1 = sum(                x(:,:,2)   )
  DX2  = sum(Delta(:,:)  *  x(:,:,1)**2)
  DY2  = sum(Delta(:,:)  *  x(:,:,2)**2)
  DXY  = sum(Delta(:,:)  *  x(:,:,1) * x(:,:,2))
  DX   = sum(Delta(:,:)  *  x(:,:,1))
  DY   = sum(Delta(:,:)  *  x(:,:,2))
  D1   = sum(Delta(:,:))

  A(1,1) = X4
  A(1,2) = X2Y2
  A(1,3) = X3Y1
  A(1,4) = X3
  A(1,5) = X2Y1
  A(1,6) = X2
  A(2,1) = X2Y2
  A(2,2) =   Y4
  A(2,3) = X1Y3
  A(2,4) = X1Y2
  A(2,5) =   Y3
  A(2,6) =   Y2
  A(3,1) = X3Y1
  A(3,2) = X1Y3
  A(3,3) = X2Y2
  A(3,4) = X2Y1
  A(3,5) = X1Y2
  A(3,6) = X1Y1
  A(4,1) = X3
  A(4,2) = X1Y2
  A(4,3) = X2Y1
  A(4,4) = X2
  A(4,5) = X1Y1
  A(4,6) = X1
  A(5,1) = X2Y1
  A(5,2) =   Y3
  A(5,3) = X1Y2
  A(5,4) = X1Y1
  A(5,5) =   Y2
  A(5,6) =   Y1
  A(6,1) = X2
  A(6,2) =   Y2
  A(6,3) = X1Y1
  A(6,4) = X1
  A(6,5) =   Y1
  A(6,6) = 1.d0
  B(1)   = DX2
  B(2)   = DY2
  B(3)   = DXY
  B(4)   = DX
  B(5)   = DY
  B(6)   = D1
  call elimi(6,C,A,B)
  !call linalglu(6,A,B,C)
  write (6, *) 'C = ', C
  det = 4.d0 * C(1)*C(2) - C(3)**2
  xa = (-2.d0 * C(2)*C(4) + C(3)*C(5)) / det
  ya = (C(3)*C(4) - 2.d0 * C(1)*C(5)) / det
  write (6, *) '(xa, ya) = ', xa, ya
  !xp(1) = xa
  !xp(2) = ya
  !xp(3) = 0.d0

  ! Hessian
  H(1,1) = C(1)
  H(2,2) = C(2)
  H(1,2) = C(3)
  H(2,1) = C(3)
  lambda(1) = 0.5d0*(C(1)+C(2)) + sqrt(0.25d0*(C(1)+C(2))**2 - (C(1)*C(2) - C(3)**2))
  lambda(2) = 0.5d0*(C(1)+C(2)) - sqrt(0.25d0*(C(1)+C(2))**2 - (C(1)*C(2) - C(3)**2))
  write (6, *) 'lambda = ', lambda

  ! eigenvectors
  V1(1) = 1.d0
  V1(2) = - (C(1) - lambda(1)) / C(2)
  V1    = V1 / sqrt(sum(V1**2)) / sqrt(abs(lambda(1)))

  V2(1) = 1.d0
  V2(2) = - (C(1) - lambda(2)) / C(2)
  V2    = V2 / sqrt(sum(V2**2)) / sqrt(abs(lambda(2)))

  write (6, *) 'eigenvectors: '
  write (6, *) 'v1 = ', v1
  write (6, *) 'v2 = ', v2
  open  (iu, file='eigenvectors.dat')
  write (iu, *) xa, ya
  write (iu, *) xa+V1(1), ya+V1(2)
  write (iu, *)
  write (iu, *) xa, ya
  write (iu, *) xa+V2(1), ya+V2(2)
  close (iu)




  do j=1,np
  do i=1,np
     X1 = x(i,j,1)
     Y1 = x(i,j,2)
     DX = C(1)*X1**2 + C(2)*Y1**2 + C(3)*X1*Y1 + C(4)*X1 + C(5)*Y1 + C(6)
     DY = C(1)*(X1-xa)**2 + C(2)*(Y1-ya)**2 + C(3)*(X1-xa)*(Y1-ya) + &
          C(6) - (xa**2*C(1) + ya**2*C(2) + xa*ya*C(3))
     write (96, '(5e18.10)') X1, Y1, DX, DY, Delta(i,j)
  enddo
  enddo

  Dout = C(6) - (xa**2*C(1) + ya**2*C(2) + xa*ya*C(3))

  ilevel = ilevel + 1
  !if (ilevel <= 4) call evaluateX (xp, v1, v2, 0.d0, xp, Dout, ilevel)

  end subroutine evaluateX

end subroutine hyperbolic_fixed_point_v4


!subroutine linalglu (n4, afi, bfi, xfo)
!  use fgsl
!  implicit none
!  !integer(fgsl_size_t), parameter :: n = 4_fgsl_size_t
!  integer              :: n4
!  real*8               :: afi(n4,n4), bfi(n4), xfo(n4)
!
!
!  integer(fgsl_size_t) :: n
!  integer(fgsl_int) :: status, signum
!  type(fgsl_matrix) :: a
!  type(fgsl_vector) :: b, x
!  real(fgsl_double), target :: af(n4, n4), bf(n4), xf(n4)
!  type(fgsl_permutation) :: p
!!
!  n = n4
!  a = fgsl_matrix_init(type=1.0_fgsl_double)
!  b = fgsl_vector_init(type=1.0_fgsl_double)
!  x = fgsl_vector_init(type=1.0_fgsl_double)
!  p = fgsl_permutation_alloc(n)
!
!!  af = reshape((/0.18d0, 0.60d0, 0.57d0, 0.96d0, &
!!                 0.41d0, 0.24d0, 0.99d0, 0.58d0, &
!!                 0.14d0, 0.30d0, 0.97d0, 0.66d0, &
!!                 0.51d0, 0.13d0, 0.19d0, 0.85d0/), (/ 4, 4 /))
!!  bf = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0 /) 
!  af = afi
!  bf = bfi
!  status = fgsl_matrix_align(af, n, n, n, a)
!  status = fgsl_vector_align(bf, n, b, n, 0_fgsl_size_t, &
!       1_fgsl_size_t)
!  status = fgsl_vector_align(xf, n, x, n, 0_fgsl_size_t, &
!       1_fgsl_size_t)
!
!  status = fgsl_linalg_LU_decomp (a, p, signum)
!  status = fgsl_linalg_LU_solve (a, p, b, x)
!
!  write(*, *) 'x = '
!  write(*, fmt='(F12.5)') xf
!  
!  call fgsl_matrix_free(a)
!  call fgsl_vector_free(b)
!  call fgsl_vector_free(x)
!  call fgsl_permutation_free(p)
!
!  xfo = xf
!end subroutine linalglu




subroutine hyperbolic_fixed_point_v3
  use iso_fortran_env
  use run_control, only: x_start, Trace_Step, Trace_Method
  use grid
  use fieldline
  use equilibrium
  use math
  implicit none

  integer, parameter :: iu = 72
  integer, parameter :: np = 10

  real(real64), dimension(np,np,3) :: x, Delta

  type(t_fieldline)  :: F
  real(real64) :: Dphi, Xout(3), Dout
  integer      :: n, m, ilevel

  n = 3
  m = 13
  Dphi       = pi2 * m / n
  Trace_Step = pi2 / 360.d0

  ilevel = 0
  call evaluateX (x_start, 2.d0, 0.d0, Xout, Dout, ilevel)
  write (6, *) 'final approximation: ', Xout(1:2)
  write (6, *) 'quality parameter:   ', Dout


  contains

  recursive subroutine evaluateX (x0, dl, alpha, xp, D, ilevel)
  real(real64), intent(in)  :: x0(3), dl, alpha
  real(real64), intent(out) :: xp(3), D
  integer, intent(inout)    :: ilevel

  real(real64) :: Delta(0:2)
  integer      :: i, j, i0, j0


  ! setup grid
  write (6, *)
  write (6, 1000) x0(1)-dl, x0(1)+dl, x0(2)-dl, x0(2)+dl
 1000 format ('box: ',f10.5,' -> ',f10.5,',  ',f10.5,' -> ',f10.5)
  do j=1,np
  do i=1,np
     x(i,j,1) = x0(1) + (-1.d0 + 2.d0*(i-1)/(np-1)) * dl
     x(i,j,2) = x0(2) + (-1.d0 + 2.d0*(j-1)/(np-1)) * dl
     x(i,j,3) = x0(3)
  enddo
  enddo


  D  = 1.d99
  i0 = 0
  j0 = 0
  ! trace field lines
  do j=1,np
  do i=1,np
     call F%init(x(i,j,:), Trace_Step, Trace_Method, FL_ANGLE)
     call F%trace(Dphi, .true.)
     Delta(1:2) = F%rc(1:2) - x(i,j,1:2)
     Delta(0)   = sqrt(sum(Delta(1:2)**2))

     !write (6, *) x(i,j,1:2), Delta(0)
     if (Delta(0) < D) then
        D  = Delta(0)
        i0 = i
        j0 = j
        xp = x(i,j,:)
     endif
  enddo
  enddo

  ! evaluate
  write (6, *) 'min. found at (i,j) = ', i0, j0
  write (6, *) '(R,     Z   ) = ', xp(1:2)
  write (6, *) '(theta, PsiN) = ', get_poloidal_angle(xp), get_PsiN(xp)
  write (6, *) ' Delta        = ', D

  ilevel = ilevel + 1
  if (ilevel < 8) call evaluateX(xp, 2.d0*dl/np, alpha, xp, D, ilevel)

  end subroutine evaluateX

end subroutine hyperbolic_fixed_point_v3






subroutine hyperbolic_fixed_point_v2
  use iso_fortran_env
  use run_control, only: Grid_File, Output_File, Trace_Step, Trace_Method
  use grid
  use fieldline
  use equilibrium
  implicit none

  integer, parameter :: iu = 72
  type(t_fieldline)  :: F
  real(real64) :: Dphi, x(3), theta, PsiN, Delta(2)
  integer      :: i, n, m, iflag

  n = 3
  m = 13
  Dphi       = pi2 * m / n
  Trace_Step = pi2 / 360.d0


  call read_grid(Grid_File, use_coordinates=COORDINATES(CYLINDRICAL))
  write (94, 1000)
  write (94, 1001) n_grid
  write (94, 1002) 0.d0
 1000 format('# grid_id = 1       (irregular RZ grid)')
 1001 format('# resolution: n_RZ           =  ',i8)
 1002 format('# angular position:  phi     =  ',f7.2)


  open  (iu, file=Output_File)
  i = 0
  grid_loop: do
     call get_next_grid_point(iflag, x)
     i = i + 1
     write (6, *) i
     if (iflag .ne. 0) exit

     theta = get_poloidal_angle (x)
     PsiN  = get_PsiN(x)
     write (94, *) theta, PsiN
     call F%init(x, Trace_Step, Trace_Method, FL_ANGLE)
     call F%trace(Dphi, .true.)
     Delta = F%rc(1:2) - x(1:2)
     write (iu, *) Delta, sum(Delta**2)
  enddo grid_loop
  close (iu)

end subroutine hyperbolic_fixed_point_v2




subroutine hyperbolic_fixed_point_v1
  use iso_fortran_env
  use run_control, only: x_start, Trace_Step, Trace_Method
  use grid
  use fieldline
  use equilibrium
  implicit none

  type(t_fieldline) :: F
  real(real64), dimension(:,:), pointer :: xi
  real(real64) :: Dphi, h, w, theta(5), PsiN(5), Delta(5,2)
  real(real64) :: X, Y, XY, X2, Y2, DX, DY, XDX, YDY, XDY, YDX, A(5,5), B(5), C(5), x0, y0, D, r(3)
  integer      :: i, n, m

  n = 3
  m = 13
  Dphi = pi2 * m / n

  xi => new_grid (5)

  h = 1.d0
  w = 1.d0
  xi(1,:) = x_start
  xi(2:5,3) = x_start(3)

  ! setup start points
  xi(2,1) = x_start(1) - w
  xi(2,2) = x_start(2)
  xi(3,1) = x_start(1) + w
  xi(3,2) = x_start(2)
  xi(4,1) = x_start(1)
  xi(4,2) = x_start(2) - h
  xi(5,1) = x_start(1)
  xi(5,2) = x_start(2) + h
  do i=1,5
     theta(i) = get_poloidal_angle(xi(i,:)) / pi * 180.d0
     if (theta(i) < 0) theta(i) = theta(i) + 360.d0
     PsiN(i)  = get_PsiN(xi(i,:))
  enddo
  do i=2,5
     write (98, *) theta(i), PsiN(i), xi(i,1:2)
     write (98, *) theta(1), PsiN(1), xi(1,1:2)
  enddo


  ! trace 1 period
  Trace_Step = pi2 / 360.d0
  do i=1,5
     write (6, *) i
     call F%init(xi(i,:), Trace_Step, Trace_Method, FL_ANGLE)
     call F%trace(Dphi, .true.)
     Delta(i,1:2) = F%rc(1:2) - xi(i,1:2)
     write (99, *) F%rc
     theta(i) = get_poloidal_angle(F%rc) / pi * 180.d0
     if (theta(i) < 0) theta(i) = theta(i) + 360.d0
     PsiN(i)  = get_PsiN(F%rc)
  enddo
  do i=2,5
     write (97, *) theta(i), PsiN(i), xi(i,1:2)+Delta(i,1:2)
     write (97, *) theta(1), PsiN(1), xi(1,1:2)+Delta(1,1:2)
  enddo


  ! calculate coefficients
  X   = sum(xi(:,1))
  Y   = sum(xi(:,2))
  XY  = sum(xi(:,1)*xi(:,2))
  X2  = sum(xi(:,1)**2)
  Y2  = sum(xi(:,2)**2)
  DX  = sum(Delta(:,1))
  DY  = sum(Delta(:,2))
  XDX = sum(xi(:,1)*Delta(:,1))
  YDY = sum(xi(:,2)*Delta(:,2))
  XDY = sum(xi(:,1)*Delta(:,2))
  YDX = sum(xi(:,2)*Delta(:,1))


  A      = 0.d0
  A(1,1) =  4.d0*X2
  A(3,1) = -4.d0*XY
  A(4,1) = -2.d0*X
  A(2,2) =  4.d0*Y2
  A(3,2) =  4.d0*XY
  A(5,2) =  2.d0*Y
  A(1,3) = -4.d0*XY
  A(2,3) =  4.d0*XY
  A(3,3) =  4.d0*(X2+Y2)
  A(4,3) =  2.d0*Y
  A(5,3) =  2.d0*X
  A(1,4) = -2.d0*X
  A(3,4) =  2.d0*Y
  A(4,4) =  1.d0
  A(2,5) =  2.d0*Y
  A(3,5) =  2.d0*X
  A(5,5) =  1.d0

  B(1)   = -2.d0*XDY
  B(2)   =  2.d0*YDX
  B(3)   =  2.d0*(XDX+YDY)
  B(4)   = DY
  B(5)   = DY
  call elimi(5,C,A,B)

  D  = 4.d0*C(1)*C(2) - C(3)**2
  x0 = (-2.d0*C(2)*C(4) +      C(3)*C(5)) / D
  y0 = (      C(3)*C(4) - 2.d0*C(1)*C(5)) / D

  write (6, *) 'test: ', x0, y0
  r(1) = x0
  r(2) = y0
  r(3) = 0.d0
  !write (6, *) 'test2: ', get_poloidal_angle(r), get_PsiN(r)

end subroutine hyperbolic_fixed_point_v1
