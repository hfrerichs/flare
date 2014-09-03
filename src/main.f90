program main
  use iso_fortran_env
  use parallel
  use run_control
  use bfield
  use equilibrium
  use boundary
  implicit none

  real(real64) :: t1, t2, t3

  ! Initialize parallel execution (using MPI)
  call initial_parallel()
  if (firstP) then
     write (6, 1000)
     write (6, *) 'Running FLARE (tagname Monkey Island)'
     if (nprs.gt.1) then
        write (6, *) 'on ', nprs, ' processors'
     endif
     write (6, *)
  endif


  call cpu_time(t1)
  call load_run_control()

  call setup_bfield_configuration()

  call setup_boundary()

  call cpu_time(t2)
  call run_control_main()
  call cpu_time(t3)
  if (firstP) then
     write (6, *)
     write (6, *) 'Time taken for initialization: ', t2-t1, ' seconds'
     write (6, *) 'Time taken for computation:    ', t3-t2, ' seconds'
     write (6, *)
  endif

  call finished_parallel()
 1000 format ('========================================================================')
end program main
