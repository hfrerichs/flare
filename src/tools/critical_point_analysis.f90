!=======================================================================
  subroutine critical_point_analysis(Grid_File, Output_File)
  use equilibrium
  use grid

  character(len=*), intent(in) :: Grid_File, Output_File


  integer, parameter :: iu = 54

  type(t_grid)   :: G
  type(t_Xpoint) :: C
  real(real64), dimension(:,:), allocatable :: xk
  real(real64)   :: x(2), H(2,2), r, lambda1, lambda2, v1(2), v2(2), r3(3)
  integer        :: i, j, k, ind, ierr


  write (6, 1000)
  call G%load(Grid_File)
  allocate (xk(G%nodes(),2))

  ind = 0
  open  (iu, file=Output_File)
  loop: do ig=1,G%nodes()
     r3 = G%node(ig)
     x  = r3(1:2)

     ! run Newton method to find critical point from x
     x = find_X(x, Hout=H)

     if (x(1) < 0.d0) then
        write (iu, *) 0
        cycle ! not a valid critical point
     endif

     ! check if present critical point is identical to previous ones
     do k=1,ind
        r = sqrt(sum((xk(k,:)-x)**2))
        if (r < 1.d-8) then
           write (iu, *) k
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
     write (iu     , *) ind
  enddo loop
  close (iu)


  ! cleanup
  deallocate (xk)

 1000 format(3x,'- Analyzing critical (hyperbolic/extremal) points')
 1002 format(8x,i0,4x,'(',f10.4,', ',f10.4,')',4x,'l1 = ',e12.4,',',4x,'l2 = ',e12.4)
  end subroutine critical_point_analysis
