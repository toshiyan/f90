!////////////////////////////////////////////////////!
! * For galaxy survey
!////////////////////////////////////////////////////!

module utils_galaxy
  use funcs,   only: erfc, lnGamma
  use myutils, only: linspace
  implicit none

  !* Galaxy Survey Parameters
  type galaxy_survey
    double precision :: fsky, Ngal, zm
  end type galaxy_survey

  type HSC
    double precision :: fsky = 0.05d0, Ngal = 35.d0, zm = 1.d0
  end type HSC

  type DES
    double precision :: fsky = 0.125d0, Ngal = 12.d0, zm = 0.5d0
  end type DES

  type LSST
    double precision :: fsky = 0.5d0, Ngal = 100.d0, zm = 1.5d0
  end type LSST

  private erfc, lnGamma
  private linspace

contains


function z0_SF(a,b,zm)  result(f)
  implicit none
  double precision, intent(in) :: a,b,zm
  double precision :: f

  f = zm * dexp(lnGamma((a+1d0)/b)-lnGamma((a+2d0)/b))

end function z0_SF


function nz_SF_scal(z,a,b,z0)  result(f)
  implicit none
  double precision, intent(in) :: z,a,b,z0
  double precision :: f, N

  N = b / (z0*dexp(lnGamma((a+1d0)/b)))
  f = N*(z/z0)**a*dexp(-(z/z0)**b)

end function nz_SF_scal


function nz_SF_arr(z,a,b,z0) result(f)
  implicit none
  double precision, intent(in) :: z(:),a,b,z0
  integer :: i
  double precision :: N, f(size(z))

  N = b / (z0*dexp(lnGamma((a+1d0)/b)))
  do i = 1, size(z)
    f(i) = N*(z(i)/z0)**a*dexp(-(z(i)/z0)**b)
  end do

end function nz_SF_arr


function nz_21cm(z,nzp) result(f)
! - Oct 15, 2015
  implicit none
  double precision, intent(in) :: z, nzp(1:3)
  double precision :: f

  f = z**(nzp(1))*dexp(-nzp(2)*z**nzp(3))

end function nz_21cm


function nz_21cm_histgram(z,fluxcut) result(f)
! * older version
  implicit none
  character(*), intent(in) :: fluxcut
  double precision, intent(in) :: z
  double precision :: N(12), f

  select case(fluxcut)
  case("100nJy")
    N = [10815d0,25483d0,27587d0,23100d0,21093d0,16764d0,14190d0,12646d0,9927d0,8462d0,7405d0,6396d0]
  case("1uJy")
    N = [5445d0,10987d0,12171d0,10327d0,7655d0,5434d0,3870d0,2928d0,2216d0,1675d0,1344d0,1076d0]
  case("5uJy")
    N = [2642d0,4499d0,4405d0,3381d0,2155d0,1403d0,932d0,648d0,448d0,317d0,231d0,172d0]
  case("10uJy")
    N = [1701d0,2731d0,2548d0,1875d0,1126d0,690d0,431d0,282d0,186d0,125d0,88d0,65d0]
  end select

  f = N(int(2*z)+1)/sum(N)

end function nz_21cm_histgram


function nz_delta_scal(z,zmean,zwidth) result(f)
  implicit none
  double precision, intent(in) :: z, zmean, zwidth
  double precision :: f

  f = 0d0
  if(abs(z-zmean)<zwidth)  f = 1d0

end function nz_delta_scal


function pz_SF_scal(z,zi,sigma,zbias)  result(f)
  implicit none
  double precision, intent(in) :: z, zi(1:2), sigma, zbias
  double precision :: f, s, c

  f = 0d0
  c = zbias*(1d0+z)
  if(sigma==0d0) then
    if(z-c>=zi(1).and.z-c<zi(2)) f = 1d0
  else
    s = sigma*(1d0+z)*dsqrt(2d0)
    f = (erfc((zi(1)-z+c)/s)-erfc((zi(2)-z+c)/s))*0.5d0
  end if

end function pz_SF_scal


function pz_SF_arr(z,zi,sigma,zbias)  result(f)
  implicit none
  !I/O
  double precision, intent(in) :: z(:), zi(1:2), sigma, zbias
  !internal
  integer :: i
  double precision :: f(size(z)), s, c

  do i = 1, size(z)
    f(i) = pz_SF_scal(z(i),zi,sigma,zbias)
  end do

end function pz_SF_arr


subroutine zbin_SF(a,b,z0,zb)
!* divide total galaxy distribution (ns) into z-bins
  implicit none
  double precision, intent(in) :: a, b, z0
  double precision, intent(out) :: zb(:)
  integer :: i, j, n, nz
  integer, parameter :: jn = 5000
  double precision :: z(jn), g(jn), zz(size(zb)-1), gmax, zmax

  zmax = 20d0
  z = linspace(0d0,zmax,jn)
  nz = size(zb)-1

  gmax = sum(nz_SF_arr(z,a,b,z0))*(z(2)-z(1))
  do j = 1, jn
    z = linspace(0d0,j*zmax/dble(jn),jn)
    g(j) = sum(nz_SF_arr(z,a,b,z0))*(z(2)-z(1))
    do n = 1, nz
      if(g(j) < n*gmax/dble(nz)) zz(n) = j*zmax/dble(jn)
    end do
  end do

  zb(1) = 0d0
  zb(2:nz+1) = zz(1:nz)
  write(*,*) 'zbin =', zb

end subroutine zbin_SF


function integ_SF(zb,a,b,z0,sigma,zbias)  result(f)
  implicit none
  !I/O
  double precision, intent(in) :: zb(1:2), a, b, z0, sigma, zbias
  !internal
  integer :: i
  double precision :: f, z(1:100000)

  z = linspace([0d0,10d0],100000)
  f = 0d0
  do i = 1, 100000
    f = f + nz_SF_scal(z(i),a,b,z0) * pz_SF_scal(z(i),zb,sigma,zbias)
  end do
  f = f * (z(2)-z(1))  !multiply dz

end function integ_SF


subroutine ngal_SF(frac,zb,a,b,z0,sigma,zbias)
  implicit none
  !I/O
  double precision, intent(in) :: zb(:), a,b,z0,sigma,zbias
  double precision, intent(out) :: frac(:)
  !internal
  integer :: n, i, num

  num  = size(frac)
  frac = 0d0
  if (size(zb)/=num+1) stop 'error (ngal_SF): size of zb is strange'

  do n = 1, num
    frac(n) = integ_SF(zb(n:n+1),a,b,z0,sigma,zbias)
  end do
  frac = frac/sum(frac)
  write(*,*) 'frac ngal = ', frac

end subroutine ngal_SF


end module utils_galaxy

