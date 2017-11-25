!===============================================================================
! Two dimensional data sets (nrow x ncol)
!===============================================================================
module dataset
  use iso_fortran_env
  implicit none
  private


  integer, parameter, public :: &
     SILENT  = 0, &
     VERBOSE = 1


  type :: t_data_abstract
     character(len=256) :: key = '', label = ''
  end type t_data_abstract

  type :: t_derived_data_abstract
     character(len=256) :: key = '', recipe = '', label = ''
  end type t_derived_data_abstract


  type, public :: t_dataset
     integer :: nrow, ncol, nrow_offset

     integer :: ndim = -1
     character(len=256) :: geometry
     type(t_data_abstract), dimension(:), pointer :: col => null()
     type(t_derived_data_abstract), dimension(:), pointer :: der => null()

     real(real64), dimension(:,:), pointer :: x => null()

     character(len=256) :: title = ''

     contains
     procedure :: load
     procedure :: plot
     procedure :: store => plot
     procedure :: new
     procedure :: set_info
     procedure :: set_column_info
     procedure :: add_derived_data
     procedure :: extend
     procedure :: resize
     procedure :: destroy
     procedure :: mpi_allreduce
     procedure :: sort_rows
     procedure :: smooth
  end type t_dataset

  type(t_dataset), public, parameter :: Empty_dataset = t_dataset(0,0,0,0,'',null(),null(),null())

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
     elseif (i==1) then
        this%title = str
     endif
  enddo read_loop
  close (iu)
  this%ncol = ncol
  deallocate (tmp)
  if (this%title == '') this%title = data_file

  allocate (this%col(ncol))

  return
 1000 format ('found ',i6,' data lines in file: ',a)
 4000 format (a)
  end subroutine load
!=======================================================================



!=======================================================================
! iu             output unit number
! filename       output filename
! nelem          number of elements to output
! append         append to existing file (in conjunction with filename)
!=======================================================================
  subroutine plot(this, iu, filename, nelem, formatstr, append)
  class (t_dataset), intent(in)           :: this
  integer,           intent(in), optional :: iu, nelem
  character(len=*),  intent(in), optional :: filename, formatstr
  logical,           intent(in), optional :: append


  character(len=11) :: f
  logical           :: append_
  integer           :: i, iu0, n0, n


  ! set default unit number for output
  iu0 = 99

  ! Unit number given for output?
  if (present(iu)) iu0 = iu

  ! Output_File given?
  if (present(filename)) then
     append_ = .false.
     if (present(append)) append_ = append

     if (append_) then
        open  (iu0, file=filename, position='append')
        write (iu0, *)
     else
        open  (iu0, file=filename)
     endif
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


  ! writa data info
  if (this%ndim >= 0) write (iu0, 2001) this%ndim
  if (this%geometry /= '') write (iu0, 2002) trim(this%geometry)
 2001 format("# FLARE DATA DIMENSION ",i0,'D')
 2002 format("# FLARE GEOMETRY ",a)
  do i=1,this%ncol
     if (this%col(i)%key /= '') write (iu0, 2003) trim(this%col(i)%key), trim(this%col(i)%label)
  enddo
 2003 format("# FLARE DATA COLUMN ",a,' "',a,'"')
  if (associated(this%der)) then
  do i=1,size(this%der)
     write (iu0, 2004) trim(this%der(i)%key), trim(this%der(i)%recipe), &
                       trim(this%der(i)%label)
  enddo
  endif
 2004 format("# FLARE DERIVED DATA ",a,'="',a,'" "',a,'"')


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

  allocate (this%col(ncol))

  end subroutine new
!=======================================================================



!=======================================================================
  subroutine set_info(this, ndim, geometry)
  class(t_dataset)             :: this
  integer,          intent(in) :: ndim
  character(len=*), intent(in) :: geometry


  this%ndim     = ndim
  this%geometry = geometry

  end subroutine set_info
!=======================================================================



!=======================================================================
  subroutine set_column_info(this, col, key, label)
  class(t_dataset)             :: this
  integer,          intent(in) :: col
  character(len=*), intent(in) :: key, label


  if (col < 1  .or.  col > this%ncol) then
     write (6, 9000) col
 9000 format("error in t_dataset%set_column_info: col = ",i0," out of range!")
     stop
  endif
  this%col(col)%key   = key
  this%col(col)%label = label

  end subroutine set_column_info
!=======================================================================



!=======================================================================
  subroutine add_derived_data(this, key, recipe, label)
  class(t_dataset)             :: this
  character(len=*), intent(in) :: key, recipe, label

  type(t_derived_data_abstract), dimension(:), allocatable :: der
  character(len=256) :: q
  integer :: i, n


  if (.not.associated(this%der)) then
     n = 1
     allocate(this%der(n))
  else
     n = size(this%der)
     allocate(der(n))
     der = this%der
     deallocate(this%der)
     allocate(this%der(n+1))
     this%der(1:n) = der
     n = n + 1
  endif


  this%der(n)%key    = key
  this%der(n)%recipe = recipe
  this%der(n)%label  = label

  end subroutine add_derived_data
!=======================================================================



!=======================================================================
  subroutine extend(this, additional_rows)
  class(t_dataset)    :: this
  integer, intent(in) :: additional_rows

  real(real64), dimension(:,:), allocatable :: tmp
  integer :: rows, n, n0, m


  n    = this%nrow
  m    = this%ncol
  allocate (tmp(n, m))
  tmp  = this%x
  rows = n + additional_rows

  n0   = this%nrow_offset
  call this%new(rows, m, n0)
  this%x(1+n0:n+n0,:) = tmp

  deallocate (tmp)

  end subroutine extend
!=======================================================================



!=======================================================================
  subroutine resize(this, nrow)
  class(t_dataset)    :: this
  integer, intent(in) :: nrow

  real(real64), dimension(:,:), allocatable :: tmp
  integer :: n, n0, m


  n  = this%nrow
  n0 = this%nrow_offset
  m  = this%ncol

  if (nrow > n) then
     allocate (tmp(n,m))
     tmp = this%x
     call this%new(nrow, m, n0)
     this%x(1+n0:n+n0,:) = tmp
  elseif (nrow < n) then
     allocate (tmp(nrow,m))
     tmp = this%x(1+n0:nrow+n0,:)
     call this%new(nrow, m, n0)
     this%x = tmp
  else
     return
  endif
  deallocate (tmp)

  end subroutine resize
!=======================================================================



!=======================================================================
  subroutine destroy(this)
  class(t_dataset) :: this

  if (associated(this%x)) deallocate(this%x)
  if (associated(this%col)) deallocate(this%col)
  end subroutine destroy
!=======================================================================



!=======================================================================
  subroutine mpi_allreduce(this)
#if defined(MPI)
  use parallel
#endif
  class(t_dataset) :: this

#if defined(MPI)
  call wait_pe()
  call sum_real_data (this%x, this%nrow*this%ncol)
#endif

  end subroutine mpi_allreduce
!=======================================================================



!=======================================================================
  subroutine sort_rows(this, icol)
  class(t_dataset)    :: this
  integer, intent(in) :: icol

  integer :: n, n0, m


#ifdef FLARE
  m  = this%ncol
  if (icol < 1  .or.  icol > m) then
     write (6, *) 'error in t_dataset%sort_rows: invalid column number ', icol, '!'
     stop
  endif
  n  = this%nrow
  n0 = this%nrow_offset
  call quicksort(this%x, n, m, 1, n, icol)
#else
  write (6, *) 'error: support for sorting datasets is not included!'
  stop
#endif

  end subroutine sort_rows
!=======================================================================



!=======================================================================
! smooth data
! input:
!    column		select column for smoothing, <=0: smooth all columns
!    smooth_width
!=======================================================================
  subroutine smooth(this, smooth_width, column)
  class(t_dataset)    :: this
  integer, intent(in) :: smooth_width, column

  real(real64), dimension(:,:), allocatable :: tmp
  integer :: i, i1, i2, j, j1, j2, n, n0


  j2 = this%ncol
  if (column > j2) then
     write (6, *) 'error in t_dataset%smooth: column > ', j2, ' not allowed!'
     stop
  endif

  if (column < 0) then
     j1 = 1
  else
     j1 = column
     j2 = column
  endif


  n  = this%nrow
  n0 = this%nrow_offset
  allocate (tmp(1+n0:n+n0, j1:j2))
  do i=1+n0,n+n0
     i1 = max(i-smooth_width,1+n0)
     i2 = min(i+smooth_width,n+n0)

     do j=j1,j2
        tmp(i, j) = sum(this%x(i1:i2, j)) / (i2-i1+1)
     enddo
  enddo

  this%x(:,j1:j2) = tmp
  deallocate (tmp)

  end subroutine smooth
!=======================================================================

end module dataset
