! * Testing elarray

program test_anaflat
  use myutils, only: savetxt, linspace
  use myconst, only: dlc, pi
  use anaflat
  implicit none
  integer :: i, el(2), nn(2), npix, simn
  double precision :: D(2)
  double precision, allocatable :: els(:), linv(:), e2d(:,:)

  !* variables
  simn = 10
  nn   = [236,100]
  npix = nn(1)*nn(2)
  D    = [59d0,25d0]*pi/180d0
  eL   = [1,500]

  allocate(els(npix),linv(npix),e2d(nn(1),nn(2)))
  ! check speed
  do i = 1,30000
    els = elarray(nn,D)
  end do
  deallocate(els,linv,e2d)

end program test_anaflat

