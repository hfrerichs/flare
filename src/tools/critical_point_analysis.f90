!=======================================================================
! Scan for critical points on grid G
!=======================================================================
  subroutine critical_point_analysis(G, D, ind)
  use equilibrium
  use grid
  use dataset
  implicit none

  type(t_grid),    intent(in)  :: G
  type(t_dataset), intent(out) :: D
  integer,         intent(out) :: ind

  ! minimum distance between critical points to be resolved
  ! this should be on the order of the grid resolution
  real(real64), parameter :: min_dist = 1.0d0

  type(t_Xpoint) :: C
  real(real64), dimension(:,:), allocatable :: xk
  real(real64)   :: x(2), H(2,2), r, lambda1, lambda2, v1(2), v2(2), r3(3)
  integer        :: i, j, k, ig, ierr


  write (6, 1000)
  call D%new(G%nodes(), 4)
  allocate (xk(G%nodes(), 2))

  ind = 0
  loop: do ig=1,G%nodes()
     r3 = G%node(ig)
     x  = r3(1:2)

     ! run Newton method to find critical point from x
     x = find_X(x, Hout=H)

     ! this is not a valid critical point
     if (x(1) < 0.d0) cycle

     ! check if present critical point is identical to previous ones
     do k=1,ind
        r = sqrt(sum((xk(k,:)-x)**2))
        if (r < min_dist) then
           D%x(ig,1) = k
           cycle loop
        endif
     enddo

     ! so this is a new critical point, run analysis
     C%X = x; C%H = H
     call C%analysis(lambda1, lambda2, v1, v2, ierr)

     ! add present point to list
     ind = ind + 1
     xk(ind,:) = x
     write (6, 1002) ind, x, lambda1, lambda2
     D%x(ig,1)    = ind
     D%x(ind,2:3) = x
     D%x(ind,4  ) = ierr
  enddo loop


  ! cleanup
  deallocate (xk)

 1000 format(3x,'- Analyzing critical (hyperbolic/extremal) points')
 1002 format(8x,i0,4x,'(',f10.4,', ',f10.4,')',4x,'l1 = ',e12.4,',',4x,'l2 = ',e12.4)
  end subroutine critical_point_analysis
