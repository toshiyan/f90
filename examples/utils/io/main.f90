program testio
  use myconst, only: dlc
  use io, only: savefile, loadfile
  implicit none
  integer :: i, s, n=1000000
  double precision, allocatable :: dat(:,:)
  complex(dlc), allocatable :: dc(:,:)

  allocate(dat(2,n),dc(2,n))
  dat = 1d0
  dc = 1d0
  dc(1,2) = 0d0

  call savefile('1.dat',dat,ow=.true.)
  call loadfile('1.dat',dat)
  open(unit=20,file='test.dat',status='replace',form='unformatted')
  write(20) dat
  close(20)
  open(unit=20,file='test.dat',status='old',form='unformatted')
  read(20) dat
  close(20)

  write(*,*) sum(dat)
  deallocate(dat,dc)

end program testio

