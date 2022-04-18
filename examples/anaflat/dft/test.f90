!///////////////////////////////////////////////////////////////////////!
! * DFT test
! * Toshiya Namikawa
! - Last Modified: Fri 30 Mar 2018 11:19:07 PM PDT
!///////////////////////////////////////////////////////////////////////!

program test_dft
  use myutils, only: savetxt, linspace
  use myconst, only: dlc, pi
  use anaflat, only: elarray_inv
  use myfftw
  implicit none
  integer :: i, el(2), nn(2), npix
  double precision :: D(2)
  complex(dlc), allocatable, dimension(:) :: alm

  !* variables
  nn   = [100,100]
  npix = nn(1)*nn(2)
  D    = [10d0,10d0]*pi/180d0

  allocate(alm(npix))
  write(*,*) 'TEST_dft: Delta(x=y=0)=', dble(nn(1)*nn(2))/D(1)/D(2)
  write(*,*) 'TEST_dft: Delta(lx=ly=0)=', D(1)*D(2)

!  1) alm = 1 -> map = Delta(x,y)
  !* inverse dft
  alm = 1d0
  call dft_1darray(alm,nn,D,-1)
  call savetxt('dft_alm.dat',alm)

!  2) map = 1 -> alm = Delta(lx,ly)
  !* dft
  alm = 1d0
  call dft_1darray(alm,nn,D,1)
  call savetxt('dft_map.dat',alm)

  call dft_1darray(alm,nn,D,-1)
  call savetxt('dft_unity.dat',alm)

! 3) red spectrum
  alm = cmplx(elarray_inv(nn,D))

  !* inverse dft
  call dft_1darray(alm,nn,D,-1)
  call savetxt('dft_redalm.dat',alm)

  deallocate(alm)

end program test_dft

