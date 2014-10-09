!===============================================================================
! Two dimensional data sets (nrow x ncol)
!===============================================================================
module dataset
  use iso_fortran_env
  use parallel
  implicit none
  private


  integer, parameter, public :: &
     SILENT  = 0, &
     VERBOSE = 1


  type, public :: t_dataset
     integer :: nrow, ncol, nrow_offset

     real(real64), dimension(:,:), pointer :: x => null()

     contains
     procedure :: load, plot, new, destroy, mpi_allreduce
  end type t_dataset

  type(t_dataset), public, parameter :: Empty_dataset = t_dataset(0,0,0,null())

  contains
!=======================================================================



!=======================================================================
! load data from file "data_file"
! optional input/output:
!	columns		expected number of columns in data file
!	output=VERBOSE  write number of rows found in data file
!	header		return leading comment line "# ..." of data file
!       nrow_offset     the first element of x is designated 1+nrow_offset
!=======================================================================
  subroutine load(this, data_file, columns, output, header, nrow_offset)
  class (t_dataset), intent(inout)         :: this
  character(len=*),  intent(in)            :: data_file
  integer,           intent(in), optional  :: columns
  integer,           intent(in), optional  :: output
  character(len=*),  intent(out), optional :: header
  integer,           intent(in), optional  :: nrow_offset


  integer, parameter                       :: iu = 42

  real(real64), dimension(:), allocatable  :: tmp
  character(len=256) :: str
  integer            :: i, j, n0, ncount, ncol, icom
  logical            :: report


  ! display messages
  report = .true.
  if (present(output)) then
     if (output == SILENT) report = .false.
  endif


  ! number of data columns
  ncol = 2
  if (present(columns)) then
     ncol = columns
  endif


  ! header
  if (present(header)) then
     header = ''
  endif


  ! parse date file to get number of data lines
  open  (iu, file=data_file)
  ncount = 0
  icom   = 0
  parse_loop: do
     read (iu, 4000, end=2000) str
     do while (str(1:1) .eq. '#')
        if (icom == 0  .and. present(header)) header = str(3:82)
        read (iu, *, end=2000) str
        icom = icom + 1
     enddo
     ncount = ncount + 1
  enddo parse_loop
 2000 rewind(iu)
  if (report) write (6,1000) ncount, data_file(1:len_trim(data_file))


  ! allocate memory
  this%nrow = ncount
  n0        = 0
  if (present(nrow_offset)) n0 = nrow_offset
  this%nrow_offset = n0
  if (associated(this%x)) deallocate (this%x)
  allocate (this%x(1+n0:ncount+n0,ncol))
  allocate (tmp(ncol))


  ! read actual data
  j = 1+n0
  read_loop: do i=1,ncount+icom
     read (iu, 4000) str
     if (str(1:1) .ne. '#') then
        read (str, *) tmp
        this%x(j,:) = tmp
        j = j + 1
     endif
  enddo read_loop
  close (iu)
  this%ncol = ncol
  deallocate (tmp)


  return
 1000 format ('found ',i6,' data lines in file: ',a)
 4000 format (a)
  end subroutine load
!=======================================================================



!=======================================================================
! iu             output unit number
! filename       output filename
! nelem          number of elements to output
!=======================================================================
  subroutine plot(this, iu, filename, nelem, formatstr)
  class (t_dataset), intent(in)           :: this
  integer,           intent(in), optional :: iu, nelem
  character(len=*),  intent(in), optional :: filename, formatstr


  character(len=11) :: f
  integer           :: i, iu0, n0, n


  ! set default unit number for output
  iu0 = 99

  ! Unit number given for output?
  if (present(iu)) iu0 = iu

  ! Output_File given?
  if (present(filename)) then
     open  (iu0, file=filename)
  endif


  ! output format
  if (present(formatstr)) then
     write (f, 1001) this%ncol, formatstr
  else
     write (f, 1000) this%ncol
  endif
 1000 format ('(',i3,'e18.10)')
 1001 format ('(',i3,a,')')

  ! offset and number of elements
  n0 = this%nrow_offset
  n  = this%nrow
  if (present(nelem)) then
     if (nelem > this%nrow) then
        write (6, *) 'error in subroutine t_dataset%plot, nelem > nrow!'
        stop
     else
        n = nelem
     endif
  endif

  ! write data
  do i=1+n0,n+n0
     write (iu0, f) this%x(i,:)
  enddo


  ! Output_File given?
  if (present(filename)) then
     close (iu0)
  endif

  end subroutine plot
!=======================================================================



!=======================================================================
  subroutine new(this, nrow, ncol, nrow_offset)
  class(t_dataset)    :: this
  integer, intent(in) :: nrow, ncol
  integer, intent(in), optional :: nrow_offset

  integer :: n0


  n0 = 0
  if (present(nrow_offset)) n0 = nrow_offset

  call this%destroy()
  this%nrow        = nrow
  this%ncol        = ncol
  this%nrow_offset = n0
  allocate (this%x(1+n0:nrow+n0, ncol))
  this%x    = 0.d0

  end subroutine new
!=======================================================================



!=======================================================================
  subroutine destroy(this)
  class(t_dataset) :: this

  if (associated(this%x)) deallocate(this%x)
  end subroutine destroy
!=======================================================================



!=======================================================================
  subroutine mpi_allreduce(this)
  class(t_dataset) :: this

  call wait_pe()
  call sum_real_data (this%x, this%nrow*this%ncol)

  end subroutine mpi_allreduce
!=======================================================================

end module dataset
