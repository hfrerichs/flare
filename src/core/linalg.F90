!===============================================================================
! Linear algebra module
!===============================================================================
module linalg
  use iso_fortran_env
  implicit none

  integer, parameter :: &
     GAUSS_ELIM    = 0, &
     LIB_FGSL      = 1

  contains
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  ! Solve system of linear equations (A*x = b)
  !---------------------------------------------------------------------
  subroutine solve(n, A, b, x, algorithm)
  integer,      intent(in)  :: n
  real(real64), intent(in)  :: A(n,n), b(n)
  real(real64), intent(out) :: x(n)
  integer,      intent(in), optional :: algorithm

  integer :: use_algorithm


  use_algorithm = GAUSS_ELIM
  if (present(algorithm)) use_algorithm = algorithm

  select case(use_algorithm)
  case(GAUSS_ELIM)
     call gauss_elim_solver(n, A, b, x)
  case(LIB_FGSL)
#if defined(FGSL)
     call fgsl_solver(n, A, b, x)
#else
     write (6, *) 'error: FLARE has been compiled without FGSL support!'
     stop
#endif
  case default
     write (6, *) 'error in subroutine solve (module linalg): invalid algorithm ', use_algorithm
     stop
  end select

  end subroutine solve
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  ! Solve system of linear equations by Gaussian elimination
  !---------------------------------------------------------------------
  subroutine gauss_elim_solver(n, b0, c0, a0)
  integer,intent (in)          :: n
  real(real64) ,dimension(n,n) :: b0
  real(real64) ,dimension(  n) :: a0,c0

  real(real64)  :: CVECT,FACT
  integer :: i,ie,j,M
  integer,dimension(2) :: piv
  intrinsic sum,maxloc,maxval

  integer,dimension(:  ),allocatable :: SEQ,NO_0  
  real(real64) ,dimension(:  ),allocatable :: a,c
  real(real64) ,dimension(:,:),allocatable :: b


  allocate(SEQ(N),NO_0(n),a(n),b(n,n),c(n))
  b = b0
  c = c0
  do i=1,n
     seq(i)=i
  enddo 


  do i=2,n
     ! find pivotting element
     piv= maxloc(abs(b(i-1:n,i-1:n)))
     piv = piv + i-2

     ! row interchange
     if(piv(1).ne.i-1) then
        a(       i-1:n) = b(   i-1,i-1:n)
        b(i-1   ,i-1:n) = b(piv(1),i-1:n)
        b(piv(1),i-1:n) = a(       i-1:n)

        cvect      = c(i-1)
        c(i-1)     = c(piv(1))
        c(piv(1))  = cvect
     endif

     ! column interchange
     if(piv(2).ne.i-1) then
        a(1:n       ) = b(1:n,i-1   )
        b(1:n,i-1   ) = b(1:n,piv(2))
        b(1:n,piv(2)) = a(1:n       )

        j          = seq(i-1)
        seq(i-1)   = seq(piv(2))
        seq(piv(2))=j
     endif

     ! Gaussian elimination
     M = 0
     do ie=i,n
        if(b(i-1,ie) /=0.) then
           M=M+1
           NO_0(M) = ie
        endif
     enddo

     do ie=i,n
        if(b(ie,i-1) /= 0.) then
           fact = -b(ie,i-1)/b(i-1,i-1)
           c(ie) = c(ie) + fact*c(i-1)
           b(ie,i-1) = 0.
           if(m>0) b(ie,No_0(1:M)) = b(ie,No_0(1:M)) +  fact*b(i-1,No_0(1:M))
        endif
     enddo
  enddo

  ! back substitution
  a(n) = C(n)/b(n,n)
  do i=n-1,1,-1
     a(i)=( c(i)- sum(b(i,i+1:n)*a(i+1:n)) )/b(i,i)
  enddo
     
  a0(seq(1:n)) = a(1:n)
  deallocate(seq,no_0,a,b,c)

  end subroutine gauss_elim_solver
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  ! Use FGSL to solve system of linear equations
  !---------------------------------------------------------------------
#if defined(FGSL)
  subroutine fgsl_solver(n4, afi, bfi, xfo)
  use fgsl
  implicit none
  integer              :: n4
  real*8               :: afi(n4,n4), bfi(n4), xfo(n4)

  integer(fgsl_size_t) :: n
  integer(fgsl_int) :: status, signum
  type(fgsl_matrix) :: a
  type(fgsl_vector) :: b, x
  real(fgsl_double), target :: af(n4, n4), bf(n4), xf(n4)
  type(fgsl_permutation) :: p


  n = n4
  a = fgsl_matrix_init(type=1.0_fgsl_double)
  b = fgsl_vector_init(type=1.0_fgsl_double)
  x = fgsl_vector_init(type=1.0_fgsl_double)
  p = fgsl_permutation_alloc(n)

  af = afi
  bf = bfi
  status = fgsl_matrix_align(af, n, n, n, a)
  status = fgsl_vector_align(bf, n, b, n, 0_fgsl_size_t, &
       1_fgsl_size_t)
  status = fgsl_vector_align(xf, n, x, n, 0_fgsl_size_t, &
       1_fgsl_size_t)

  status = fgsl_linalg_LU_decomp (a, p, signum)
  status = fgsl_linalg_LU_solve (a, p, b, x)

  call fgsl_matrix_free(a)
  call fgsl_vector_free(b)
  call fgsl_vector_free(x)
  call fgsl_permutation_free(p)

  xfo = xf

  end subroutine fgsl_solver
#endif
  !---------------------------------------------------------------------


end module linalg
