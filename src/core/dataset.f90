!===============================================================================
! Two dimensional data sets (nrow x ncol)
!===============================================================================
module dataset
  use iso_fortran_env
  implicit none
  private

  type, public :: t_dataset
     integer :: nrow, ncol
     real(real64), dimension(:,:), allocatable :: x

     contains
     procedure :: load, store
  end type t_dataset

  contains
!=======================================================================



!=======================================================================
! load data from file "data_file"
! optional input/output:
!	columns		expected number of columns in data file
!	report		write number of rows found in data file
!	header		return leading comment line "# ..." of data file
!=======================================================================
  subroutine load(this, data_file, columns, report, header)
  class (t_dataset), intent(inout)         :: this
  character(len=*),  intent(in)            :: data_file
  integer,           intent(in), optional  :: columns
  logical,           intent(in), optional  :: report
  character(len=*),  intent(out), optional :: header


  integer, parameter                       :: iu = 42

  real(real64), dimension(:), allocatable  :: tmp
  character(len=256) :: str
  integer            :: i, j, ncount, ncol, icom
  logical            :: lreport


  ! display messages
  lreport = .true.
  if (present(report)) then
     if (report .eqv. .false.) lreport = .false.
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
  if (lreport) write (6,1000) ncount, data_file(1:len_trim(data_file))


  ! allocate memory
  this%nrow = ncount
  if (allocated(this%x)) deallocate (this%x)
  allocate (this%x(ncount,ncol))
  allocate (tmp(ncol))


  ! read actual data
  j = 1
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
  subroutine store(this, data_file)
  class (t_dataset), intent(in) :: this
  character(len=*),  intent(in) :: data_file

  integer, parameter            :: iu = 42


  character*11 :: f
  integer      :: i


  write (6, *) 'writing data: ', this%nrow, this%ncol

  write (f, 1000) this%ncol
 1000 format ('(',i3,'e18.10)')

  open  (iu, file=data_file)
  do i=1,this%nrow
     write (iu, f) this%x(i,:)
  enddo
  close (iu)
  end subroutine store
!=======================================================================

end module dataset
