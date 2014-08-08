!-------------------------------------------------------------------------------
! module:	parallel
!
! description:	Supplemental routines for parallel execution
!-------------------------------------------------------------------------------
module parallel
implicit none

#if defined(parallelMPI)
    include 'mpif.h'
#endif

integer, save :: nprs = 1, mype = 0, er_par = 0
logical, save :: firstP = .true.

contains
!===============================================================================



!===============================================================================
subroutine initial_parallel

#if defined(parallelMPI)
    call MPI_INIT (er_par)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, nprs, er_par)
    call MPI_COMM_RANK (MPI_COMM_WORLD, mype, er_par)
#endif

if (mype.ne.0) firstP = .false.

return
end subroutine initial_parallel
!===============================================================================



!===============================================================================
subroutine finished_parallel

call wait_pe
#if defined(parallelMPI)
    call MPI_FINALIZE (er_par)
#endif

return
end subroutine finished_parallel
!===============================================================================



!===============================================================================
subroutine broadcast_real_s(A)
real*8, intent(inout) :: A

if (nprs.le.1) return
#if defined(parallelMPI)
    call MPI_BCAST (A, 1, MPI_REAL8, 0, MPI_COMM_WORLD, er_par)
#endif

return
end subroutine broadcast_real_s
!===============================================================================



!===============================================================================
subroutine broadcast_real(A,ndim)
integer, intent(in)   :: ndim
real*8, dimension(ndim), intent(inout) :: A

if (nprs.le.1) return
#if defined(parallelMPI)
    call MPI_BCAST (A, ndim, MPI_REAL8, 0, MPI_COMM_WORLD, er_par)
#endif

return
end subroutine broadcast_real
!===============================================================================



!===============================================================================
subroutine broadcast_inte_s(A)
integer, intent(inout) :: A

if (nprs.le.1) return
#if defined(parallelMPI)
    call MPI_BCAST (A, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, er_par)
#endif

return
end subroutine broadcast_inte_s
!===============================================================================



!===============================================================================
subroutine broadcast_inte(A,ndim)
integer, intent(in) :: ndim
integer, dimension(ndim), intent(inout) :: A

if (nprs.le.1) return
#if defined(parallelMPI)
    call MPI_BCAST (A, ndim, MPI_INTEGER, 0, MPI_COMM_WORLD, er_par)
#endif

return
end subroutine broadcast_inte
!===============================================================================



!===============================================================================
subroutine broadcast_logi(A)
logical, intent(inout) :: A

if (nprs.le.1) return
#if defined(parallelMPI)
    call MPI_BCAST (A, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, er_par)
#endif

return
end subroutine broadcast_logi
!===============================================================================



!===============================================================================
subroutine broadcast_char(A,ndim)
integer, intent(in) :: ndim
    character(ndim), intent(inout) :: A

if (nprs.le.1) return
#if defined(parallelMPI)
    call MPI_BCAST (A, ndim, MPI_CHARACTER, 0, MPI_COMM_WORLD, er_par)
#endif

return
end subroutine broadcast_char
!===============================================================================



!===============================================================================
subroutine sum_inte_data(A,ndim)
integer, intent(in)                    :: ndim
integer, dimension(ndim), intent(inout) :: A
integer, dimension(:), allocatable      :: B,C

if (nprs.le.1) return
#if defined(parallelMPI)
    call MPI_BARRIER(MPI_COMM_WORLD,ER_PAR)
    allocate (B(ndim))
    B = A
    call MPI_ALLREDUCE (B,A,NDIM,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ER_PAR)
    deallocate (B)
#endif

end subroutine sum_inte_data
!===============================================================================



!===============================================================================
subroutine sum_real_data(A,ndim)
integer, intent(in)                    :: ndim
real*8, dimension(ndim), intent(inout) :: A
real*8, dimension(:), allocatable      :: B,C

integer, parameter :: ndim_max = 100000
integer :: n1, n2

if (nprs.le.1) return
#if defined(parallelMPI)
    call MPI_BARRIER(MPI_COMM_WORLD,ER_PAR)
    if (ndim .le. ndim_max) then
        allocate (b(ndim))
        B = A
        call MPI_ALLREDUCE (B, A, ndim, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, er_par)
        deallocate (B)
    else
        allocate (b(ndim_max), c(ndim_max))
        n1 = 1
        do
            n2 = min(n1+ndim_max-1, ndim)
            B(1:n2-n1+1) = A(n1:n2)
            call MPI_ALLREDUCE (B, C, ndim_max, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, er_par)
            A(n1:n2) = C(1:n2-n1+1)
            if (n2.eq.ndim) exit
            n1 = n2 + 1
        enddo
        deallocate (B,C)
    endif
#endif

return
end subroutine sum_real_data
!===============================================================================



!===============================================================================
subroutine wait_pe()

#if defined(parallelMPI)
    call MPI_BARRIER (MPI_COMM_WORLD, er_par)
#endif

return
end subroutine wait_pe
!===============================================================================

end module parallel
