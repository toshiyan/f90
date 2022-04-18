!///////////////////////////////////////////////////////////////////////!
! * Example code for save/load ascii files
!///////////////////////////////////////////////////////////////////////!

program testascii
  use myconst, only: dlc
  use myutils, only: savetxt, loadtxt, savearray
  implicit none
  integer :: i, s, n
  double precision :: dat(2,100)
  complex(dlc) :: dc(2,100)

  !data to be saved
  dat = 1d0
  dat(1,2) = 0d0

  !for speed check
  s   = 1 

  !save double
  write(*,*) '1'
  do i = 1,s
    call savetxt('1.dat',dat(1:1,:))
  end do
  write(*,*) '2'
  do i = 1,s
    open(unit=20,file='2.dat',status='replace')
    write(20,'(E14.7)') ((dat(1:1,n)),n=1,10000)
    close(20)
  end do
  write(*,*) '3'
  do i = 1,s
    open(unit=20,file='3.dat',status='replace')
    do n = 1, 10000
      write(20,'(E14.7)') dat(1:1,n)
    end do
    close(20)
  end do

  !load double
  write(*,*) '1.1'
  do i = 1,s
    call loadtxt('1.dat',dat(1:1,:))
  end do
  write(*,*) '1.2'
  do i = 1,s
    call loadtxt('1.dat',dat(1:1,:),fsize=[1,10000])
  end do
  write(*,*) '2'
  do i = 1,s
    open(unit=20,file='2.dat',status='old')
    read(20,*) ((dat(1:1,n)),n=1,10000)
    close(20)
  end do

  !save array with index
  write(*,*) '3'
  do i = 1,s
    call savearray('3.dat',dat)
  end do


end program testascii

