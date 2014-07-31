program main
  use parallel
  use run_control
  use bfield
  use equilibrium
  use pfc
  implicit none


  ! Initialize parallel execution (using MPI)
  call initial_parallel()
  if (firstP) then
     write (6, 1000)
     write (6, *) 'Running FLARE (tagname Monkey Island)'
     write (6, *)
     if (nprs.gt.1) then
        write (6, *) 'on ', nprs, ' processors'
     endif
  endif


  call load_run_control()

  call setup_bfield_configuration()
  call setup_equilibrium()

  call setup_pfc()

  call run_control_main()

  call finished_parallel()
 1000 format ('========================================================================')
end program main
