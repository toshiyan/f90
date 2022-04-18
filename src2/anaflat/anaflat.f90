!////////////////////////////////////////////////////!
! * Utils for 2D Grid Analysis
!////////////////////////////////////////////////////!

module anaflat
  use myrandom, only: InitRandom, ranmar, Gaussian1
  use myconst,  only: dlc, iu, pi
  use myutils,  only: filelines, savetxt, loadtxt
  use io,       only: savefile
  use mycls,    only: cl2c2d
  implicit none

  interface map_smoothing
    module procedure map_smoothing_1D, map_smoothing_2D, map_smoothing_1D_cmplx, map_smoothing_2D_cmplx
  end interface map_smoothing

  interface spin_weight
    module procedure spin_weight_1darray, spin_weight_2darray
  end interface 

  interface Gaussian_alm
    module procedure gaussian_alm_1darray, gaussian_alm_2darray
  end interface Gaussian_alm

  private InitRandom, ranmar, gaussian1
  private dlc, iu, pi
  private filelines, savetxt, loadtxt
  private savefile
  private cl2c2d

  !local variable
  double precision, parameter :: twopi = 2d0*pi

contains 


subroutine array_fine(nn,s,arr1,arr2,f) 
!* make finer array with a same value
  implicit none
!
! [input]
! nn   --- x and y grids
! s    --- increasing factors 
! arr1 --- finer array with nn(1) x nn(2) grids
  integer, intent(in) :: nn(1:2), s(1:2)
  double precision, intent(in) :: arr1(:)
!
! [output]
! arr2 --- finer array with nn(1)*s(1) x nn(2)*s(2) grids
  double precision, intent(out) :: arr2(:)
!
! (optional)
! f    --- file names for input and output arrays
  character(*), intent(in), optional :: f(1:2)
!
! [internal]
  integer :: mm(2), i, j
  double precision, allocatable :: map1(:,:), map2(:,:)

  mm = nn*s
  allocate(map1(nn(1),nn(2)),map2(mm(1),mm(2)))
  if (present(f))  call savetxt(f(1),arr1)

  map1 = reshape(arr1,nn,order=[2,1])
  do i = 1, nn(1)
    do j = 1, nn(2)
      map2((i-1)*s(1)+1:(i-1)*s(1)+s(1),(j-1)*s(2)+1:(j-1)*s(2)+s(2)) = map1(i,j)
    end do
  end do
  arr2 = reshape(transpose(map2),[mm(1)*mm(2)])
  if (present(f))  call savetxt(f(2),arr2)

  deallocate(map1,map2)

end subroutine array_fine


subroutine map_write(nn,f,map1D,map2D,s) 
!* write (map1D or map2D) in a file (f) with a smoothing size (s)
  implicit none
  !I/O
  character(*), intent(in) :: f
  integer, intent(in) :: nn(2)
  integer, intent(in), optional :: s
  double precision, intent(in), optional :: map1D(:)
  double precision, intent(in), optional :: map2D(:,:)
  !internal
  integer :: rn(2), i, j, res
  double precision, allocatable :: smap(:,:), map(:,:)

  !* make 2D map
  allocate(map(nn(1),nn(2)))
  if(present(map1D))  map = reshape(map1D,nn,order=[2,1])
  if(present(map2D))  map = map2D

  !* nsides after smoothing
  if(present(s)) then 
    res = s
    rn = int(dble(nn)/dble(res))
    allocate(smap(rn(1),rn(2)))
    call map_smoothing(nn,rn,res,map,smap)
  else
    rn = nn
    allocate(smap(rn(1),rn(2)))
    smap = map
  end if
  deallocate(map)

  !* output
  call savefile(f,smap)
  deallocate(smap)

end subroutine map_write


subroutine map_smoothing_2D(nn,rn,s,map,smap)
!* input map (map,nn) is smoothed (smap,rn) with a smoothing size (s)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), rn(2), s
  double precision, intent(in) :: map(:,:)
  double precision, intent(out) :: smap(:,:)
  !internal
  integer :: i, j

  do i = 1, rn(1)
    do j = 1, rn(2)
      smap(i,j) = sum(map((i-1)*s+1:i*s,(j-1)*s+1:j*s))/dble(s**2)
    end do
  end do

end subroutine map_smoothing_2D


subroutine map_smoothing_2D_adv(nn,rn,s,map,smap)
!* An extention of map_smoothing_2D
!* set zero if one of pixels has zero value
  implicit none
  !I/O
  integer, intent(in) :: nn(2), rn(2), s
  double precision, intent(in) :: map(:,:)
  double precision, intent(out) :: smap(:,:)
  !internal
  integer :: i, j

  do i = 1, rn(1)
    do j = 1, rn(2)
      if(product(map((i-1)*s+1:i*s,(j-1)*s+1:j*s))==0) then
        smap(i,j) = 0d0
      else
        smap(i,j) = sum(map((i-1)*s+1:i*s,(j-1)*s+1:j*s))/dble(s**2)
      end if
    end do
  end do

end subroutine map_smoothing_2D_adv


subroutine map_smoothing_1D(nn,rn,s,map1D,smap1D)
!* input map (map,nn) is smoothed (smap,rn) with a smoothing size (s)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), rn(2), s
  double precision, intent(in) :: map1D(:)
  double precision, intent(out) :: smap1D(:)
  !internal
  double precision :: map(nn(1),nn(2)), smap(rn(1),rn(2))

  map = reshape(map1D,[nn(1),nn(2)],order=[2,1])
  call map_smoothing_2D(nn,rn,s,map,smap)
  smap1D = reshape(transpose(smap),[rn(1)*rn(2)])

end subroutine map_smoothing_1D


subroutine map_smoothing_1D_cmplx(nn,rn,s,map1D,smap1D)
!* input map (map,nn) is smoothed (smap,rn) with a smoothing size (s)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), rn(2), s
  complex(dlc), intent(in) :: map1D(:)
  complex(dlc), intent(out) :: smap1D(:)
  !internal
  complex(dlc) :: map(nn(1),nn(2)), smap(rn(1),rn(2))

  map = reshape(map1D,[nn(1),nn(2)],order=[2,1])
  call map_smoothing_2D_cmplx(nn,rn,s,map,smap)
  smap1D = reshape(transpose(smap),[rn(1)*rn(2)])

end subroutine map_smoothing_1D_cmplx


subroutine map_smoothing_2D_cmplx(nn,rn,s,map,smap)
!* input map (map,nn) is smoothed (smap,rn) with a smoothing size (s)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), rn(2), s
  complex(dlc), intent(in) :: map(:,:)
  complex(dlc), intent(out) :: smap(:,:)
  !internal
  integer :: i, j

  do i = 1, rn(1)
    do j = 1, rn(2)
      smap(i,j) = sum(map((i-1)*s+1:i*s,(j-1)*s+1:j*s))/dble(s**2)
    end do
  end do

end subroutine map_smoothing_2D_cmplx


!//// Arrays in Fourier plane ////!

function elarray(nn,D)  result(f)
! return absolute value of multipole in 2D
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  double precision, intent(in) :: D(1:2)
  !internal
  integer :: i, j, n
  double precision :: lx, ly, f(nn(1)*nn(2))

  f = 0d0
  n = 1
  do i = 1, nn(1)
    lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
    do j = 1, nn(2)
      ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
      if (lx/=0.or.ly/=0)  f(n) = dsqrt(lx**2+ly**2)
      n = n + 1
    end do
  end do

end function elarray


function elarray_inv(nn,D)  result(f)
! * return absolute value of inverse multipole in 2D
! * avoid pixel at ell=0
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  double precision, intent(in) :: D(1:2)
  !internal
  integer :: i, j, n
  double precision :: lx, ly, f(nn(1)*nn(2))

  f = 0d0
  n = 1
  do i = 1, nn(1)
    lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
    do j = 1, nn(2)
      ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
      if(lx/=0d0.or.ly/=0d0)  f(n) = 1d0/dsqrt(lx**2+ly**2)
      n = n + 1
    end do
  end do

end function elarray_inv


subroutine elarray_array(f,nn,D)
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  double precision, intent(in) :: D(1:2)
  double precision, intent(inout) :: f(:)
  !internal
  integer :: i, j, n
  double precision :: lx, ly

  f = 0d0
  n = 1
  do i = 1, nn(1)
    lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
    do j = 1, nn(2)
      ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
      if (lx/=0.or.ly/=0)  f(n) = dsqrt(lx**2+ly**2)
      n = n + 1
    end do
  end do

end subroutine elarray_array


function elarray_x(nn,D)  result(f)
! return ell_x in 2D
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  double precision, intent(in) :: D(1:2)
  !internal
  integer :: i, j, n
  double precision :: lx, ly, f(1:nn(1)*nn(2))

  n = 1
  do i = 1, nn(1)
    do j = 1, nn(2)
      f(n) = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
      n = n + 1
    end do
  end do

end function elarray_x


function elarray_y(nn,D)  result(f)
! return ell_y in 2D
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  double precision, intent(in) :: D(1:2)
  !internal
  integer :: i, j, n
  double precision :: lx, ly, f(nn(1)*nn(2))

  n = 1
  do i = 1, nn(1)
    do j = 1, nn(2)
      f(n) = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
      n = n + 1
    end do
  end do

end function elarray_y


subroutine elarray_wcut(nn,D,oL,els,lfac,elx,linv)
  implicit none
  integer, intent(in) :: nn(1:2), oL(1:2)
  double precision, intent(in) :: D(1:2)
  double precision, intent(out), optional :: els(:), lfac(:), elx(:), linv(:)
  integer :: n, npix
  double precision, allocatable :: els0(:), lfac0(:), linv0(:)

  npix = nn(1)*nn(2)
  allocate(els0(npix),lfac0(npix),linv0(npix)); lfac0=0d0; linv0=0d0
  els0 = elarray(nn,D)
  do n = 1, npix
    if (oL(1)<=els0(n).and.els0(n)<=oL(2)) lfac0(n) = 2d0/els0(n)**2
    if (oL(1)<=els0(n).and.els0(n)<=oL(2)) linv0(n) = 1d0/els0(n)
  end do
  if (present(els))  els  = els0
  if (present(lfac)) lfac = lfac0
  if (present(elx))  elx  = elarray_x(nn,D)
  if (present(linv)) linv = linv0
  deallocate(els0,lfac0,linv0)

end subroutine elarray_wcut


subroutine elarrays(nn,D,elx,ely,els,eli,ei2p)
!* Return elx, ely and absolute value of multipole in 2D
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  double precision, intent(in) :: D(1:2)
  double precision, intent(out), optional :: elx(:), ely(:), els(:), eli(:)
  complex(dlc), intent(out), optional :: ei2p(:)
  !internal
  integer :: i, j, n, npix
  double precision :: lx, ly
  double precision, dimension(:), allocatable :: elx0, ely0, els0, eli0
  complex(dlc), allocatable :: ei2p0(:)

  npix = nn(1)*nn(2)
  allocate(elx0(npix),ely0(npix),els0(npix),eli0(npix),ei2p0(npix))
  elx0=0d0;  ely0=0d0;  els0=0d0;  eli0=0d0;  ei2p0=0d0

  n = 1
  do i = 1, nn(1)
    lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
    do j = 1, nn(2)
      ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
      elx0(n) = lx
      ely0(n) = ly
      if (lx/=0d0.or.ly/=0d0) then
        els0(n) = dsqrt(lx**2+ly**2)
        eli0(n) = 1d0/els0(n)
        !compute exp(2i*phi) = cos(2phi) + i*sin(2phi)
        ei2p0(n) = ((lx**2-ly**2)+iu*2d0*lx*ly)/(lx**2+ly**2)
      end if
      n = n + 1
    end do
  end do

  if (present(elx)) elx=elx0
  if (present(ely)) ely=ely0
  if (present(els)) els=els0
  if (present(eli)) eli=eli0
  if (present(ei2p)) ei2p=ei2p0

  deallocate(elx0,ely0,els0,eli0,ei2p0)


end subroutine elarrays


subroutine gaussian_alm_2darray(nn,D,iL,alm,Cl,fix)
! Generate 1D-array random gaussian fields in 2D Fourier space for a given isotropic spectrum
! Note: satisfy a^*_l = a_{-l}
  implicit none

  ![input]
  ! nn(2) --- x and y grid number
  ! D(2)  --- x and y length
  ! iL(2) --- min/max multipoles of the random gaussian fields
  ! Cl(:) --- power spectrum
  integer, intent(in) :: iL(2), nn(2)
  double precision, intent(in) :: Cl(:), D(2)

  ![in-output]
  ! alm(:,:) --- random gaussian fields of a nn(1) x nn(2) array
  complex(dlc), intent(inout) :: alm(:,:)

  !(optional)
  ! fix   --- use sqrt{Cl} for alm (not random)
  logical, intent(in), optional :: fix

  !internal
  integer :: i, j, n, l0, l1
  integer :: ijmin(2), ijmax(2)
  double precision :: x, y, dx, dy, d0, l
  double precision, allocatable :: amp(:), amp2d(:,:)

  ! check
  if(size(iL)/=2) stop 'error (gaussian_alm) : size of iL is not 2'
  if(mod(nn(1),2)/=0.or.mod(nn(2),2)/=0) stop 'error (gaussian_alm) : nn(1) and/or nn(2) should be even integers'
 
  call InitRandom(-1)

  !* make cl on 2d grid
  allocate(amp(size(Cl)),amp2d(nn(1),nn(2)))
  amp = 0d0
  d0 = D(1)*D(2)
  do i = iL(1), iL(2)
    if (Cl(i)<0d0) stop 'error: cl is negative'
    amp(i) = dsqrt(d0*Cl(i)*0.5d0)  ! \sqrt(\delta(l=0)*Cl(l)/2)
  end do
  if(present(fix)) amp = dsqrt(d0*Cl)
  amp2d = 0d0
  do i = 1, nn(1)
    x = dx*dble(i-1-nn(1)*0.5d0)
    do j = 1, nn(2)
      y = dy*dble(j-1-nn(2)*0.5d0)
      l = dsqrt(x**2+y**2)
      if(iL(1)>l.or.l>iL(2)-1) cycle
      l0 = int(l)
      l1 = l0 + 1
      amp2d(i,j) = amp(l0) + (l-l0)*(amp(l1)-amp(l0))
    end do
  end do
  deallocate(amp)

  !* alm=0 if i=1 or j=1 for symmetry
  !* center: (ic,jc) = (nn(1)/2+1, nn(2)/2+1)
  alm = 0d0

  !* maximum nn
  !ijmin = max ( 2, int( nn(:)*0.5d0 + 1 - iL(2)*D(:)/twopi ) )
  !ijmax = min ( nn(:), int( nn(:)*0.5d0 + 1 + iL(2)*D(:)/twopi ) + 1 )
  ijmin = 2
  ijmax = nn

  !* dx, dy in l-space
  dx = twopi/D(1)
  dy = twopi/D(2)

  !* check
  if(dx*(nn(1)*0.5-1)<iL(2).or.dy*(nn(2)*0.5-1)<iL(2)) then
    write(*,*) 'error: inclusion of Fourier mode is incorrect'
    write(*,*) 'maximum ell should be lower than', dx*(nn(1)*0.5-1), 'or', dy*(nn(2)*0.5-1)
    stop
  end if

  ! half side (i < ic)
  do i = ijmin(1), nn(1)/2
    x = dx*dble(i-1-nn(1)*0.5d0)
    do j = ijmin(2), ijmax(2)
      y = dy*dble(j-1-nn(2)*0.5d0)
      l = dsqrt(x**2+y**2)
      if(l<iL(1).or.l>iL(2)) cycle
      if(present(fix)) then
        alm(i,j) = amp2d(i,j)
      else
        alm(i,j) = cmplx(Gaussian1(),Gaussian1())*amp2d(i,j)
      end if
    end do
  end do

  ! values on axis (i=ic) but avoid the origin (ell=0) 
  i = nn(1)/2+1
  ! x=0
  do j = ijmin(2), nn(2)/2
    y = dy*dble(j-1-nn(2)*0.5d0)
    l = abs(y)
    if(l<iL(1).or.l>iL(2)) cycle
    if(present(fix)) then
      alm(i,j) = amp2d(i,j)
    else
      alm(i,j) = cmplx(Gaussian1(),Gaussian1())*amp2d(i,j)
    end if
  end do
  do j = nn(2)/2+2, ijmax(2)
    alm(i,j) = conjg(alm(i,nn(2)-j+2))
  end do
  
  ! the other half side (i>ic)
  do i = nn(1)/2+2, ijmax(1)
    do j = ijmin(2), ijmax(2)
      alm(i,j) = conjg(alm(nn(1)-i+2,nn(2)-j+2))
    end do
  end do

  deallocate(amp2d)


end subroutine gaussian_alm_2darray


subroutine gaussian_alm_1darray(nn,D,iL,alm,Cl,fix)
! Generate 1D-array random gaussian fields in 2D Fourier space for a given isotropic spectrum
  implicit none

  ![input]
  ! nn(2) --- x and y grid number
  ! D(2)  --- x and y length
  ! iL(2) --- min/max multipoles of the random gaussian fields
  ! Cl(:) --- power spectrum
  integer, intent(in) :: iL(2), nn(2)
  double precision, intent(in) :: Cl(:), D(2)

  ![in-output]
  ! alm(:) --- nn(1) x nn(2) size random gaussian fileds
  complex(dlc), intent(inout) :: alm(:)

  !(optional)
  ! fix   --- use sqrt{Cl} for alm (not random)
  logical, intent(in), optional :: fix

  !internal
  integer :: i, j, n, npixc, npix, l0, l1
  double precision :: x, y, dd(1:2), d0, l
  double precision, allocatable :: els(:), amp(:), amp2d(:)

  ! check
  if(size(iL)/=2)   stop 'error (gaussian_alm) : size of iL is not 2'
  if(mod(nn(1),2)/=0.or.mod(nn(2),2)/=0) stop 'error (gaussian_alm) : nn(1) and/or nn(2) should be even integers'
 
  call InitRandom(-1)

  npix = nn(1)*nn(2)

  !* make cl on 2d grid
  allocate(amp(size(Cl)),els(npix),amp2d(npix))
  amp = 0d0
  d0 = D(1)*D(2)
  do i = iL(1), iL(2)
    if (Cl(i)<0d0) stop 'error: cl is negative'
    amp(i) = dsqrt(d0*Cl(i)*0.5d0)  ! \sqrt(\delta(l=0)*Cl(l)/2)
  end do
  if(present(fix)) amp = dsqrt(d0*Cl)
  els = elarray(nn,D)
  amp2d = 0d0
  do n = 1, npix
    if(iL(1)>els(n).or.els(n)>iL(2)-1) cycle
    l0 = int(els(n))
    l1 = l0 + 1
    amp2d(n) = amp(l0) + (els(n)-l0)*(amp(l1)-amp(l0))
  end do
  deallocate(amp,els)

  !* alm=0 if i=1 or j=1 for symmetry
  !* center: (ic,jc) = (nn(1)/2+1, nn(2)/2+1)
  alm = 0d0
  npixc = nn(2)*nn(1)/2 + nn(2)/2 + 1

  !* dx, dy in l-space
  dd = twopi/D

  !* check
  if(dd(1)*(nn(1)*0.5-1)<iL(2).or.dd(2)*(nn(2)*0.5-1)<iL(2)) then
    write(*,*) 'error: inclusion of Fourier mode is incorrect'
    write(*,*) 'maximum ell should be lower than', dd(1)*(nn(1)*0.5-1), 'or', dd(2)*(nn(2)*0.5-1)
    stop
  end if

  ! half side (i < ic)
  n = nn(2)+1 ! adding from j=1 to nn(2)
  do i = 2, nn(1)/2
    x = dd(1)*dble(i-1-nn(1)*0.5d0)
    n = n + 1 ! j=1
    do j = 2, nn(2)
      y = dd(2)*dble(j-1-nn(2)*0.5d0)
      l = dsqrt(x**2+y**2)
      if(l>=iL(1).and.l<=iL(2)) then
        if(present(fix)) then
          alm(n) = amp2d(n)
        else
          alm(n) = cmplx(Gaussian1(),Gaussian1())*amp2d(n)
        end if
      end if
      n = n + 1
    end do
  end do

  ! values on axis (i=ic) but avoid the origin (ell=0) 
  i = nn(1)/2+1
  ! x=0
  n = n + 1  ! j = 1
  do j = 2, nn(2)/2
    y = dd(2)*dble(j-1-nn(2)*0.5d0)
    l = abs(y)
    if(l>=iL(1).and.l<=iL(2)) then
      if(present(fix)) then
        alm(n) = amp2d(n)
      else
        alm(n) = cmplx(Gaussian1(),Gaussian1())*amp2d(n)
      end if
    end if
    n = n + 1
  end do

  !complex conjugate
  do n = npixc+1, nn(1)*nn(2)
    alm(n) = conjg(alm(2*npixc-n))
  end do

  deallocate(amp2d)
  

end subroutine gaussian_alm_1darray


subroutine gaussian_alm_te(nn,D,sL,TT,TE,EE,tlm,elm)
  implicit none
  integer, intent(in) :: nn(2), sL(2)
  double precision, intent(in) :: D(2), TT(:), TE(:), EE(:)
  complex(dlc), intent(out) :: tlm(:), elm(:)
  integer :: npix, n
  double precision, allocatable :: uni(:), els(:), TT2d(:), TE2d(:), EE2d(:)
  complex(dlc), allocatable :: ulm(:)

  npix = nn(1)*nn(2)

  allocate(uni(1:sL(2)),ulm(npix)); uni=1d0; ulm=0d0

  call gaussian_alm(nn,D,sL,tlm,TT)
  call gaussian_alm(nn,D,sL,ulm,uni)

  allocate(TT2d(npix),TE2d(npix),EE2d(npix),els(npix))
  els = elarray(nn,D)
  call cl2c2d(els,TT,sL,TT2d)
  call cl2c2d(els,TE,sL,TE2d)
  call cl2c2d(els,EE,sL,EE2d)

  do n = 1, npix 
    if (TT2d(n)==0d0) cycle
    elm(n) = (TE2d(n)/TT2d(n))*tlm(n) + ulm(n)*dsqrt(EE2d(n)-TE2d(n)**2/TT2d(n))
  end do

  deallocate(TT2d,TE2d,EE2d,els,ulm)

end subroutine gaussian_alm_te


subroutine spin_weight_1darray(sw,nn,D,trans,S)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), trans, S
  double precision, intent(in) :: D(2)
  complex(dlc), intent(out) :: sw(:)
  !internal
  integer :: i, j, n
  double precision :: ss, x, y

  n  = 1
  ss = -trans*S/abs(S)
  do i = 1, nn(1)
    x = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
    do j = 1, nn(2)
      y = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
      if(x/=0d0.or.y/=0d0) sw(n) = sw(n) * ((x**2-y**2)+iu*2d0*ss*x*y)/(x**2+y**2)
      n = n + 1
    end do
  end do

end subroutine spin_weight_1darray


subroutine spin_weight_2darray(sw,nn,D,trans,S)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), trans, S
  double precision, intent(in) :: D(2)
  complex(dlc), intent(out) :: sw(:,:)
  !internal
  integer :: i, j, n
  double precision :: ss, x, y

  ss = -trans*S/abs(S)
  do i = 1, nn(1)
    x = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
    do j = 1, nn(2)
      y = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
      if(x/=0d0.or.y/=0d0) sw(i,j) = sw(i,j) * ((x**2-y**2)+iu*2d0*ss*x*y)/(x**2+y**2)
    end do
  end do

end subroutine spin_weight_2darray


!//// Window function //////////////////////////////////////////////////////////////////////////////!

function WAP(mapsize,s,apfactor,skycut)
!* 1D window function
  implicit none
  !I/O
  double precision, intent(in) :: mapsize, s
  double precision, intent(in), optional :: skycut, apfactor
  !internal
  double precision :: s0, a, Wap, ss, x

  !map size
  a = mapsize*0.5d0
  if(present(skycut)) a = a*skycut

  !apodization length
  s0 = 1d0
  if(present(apfactor)) s0 = apfactor

  !normalized coordinate
  ss = abs(s)/a
  if(s0/=1d0) x = (1d0-ss)/(1d0-s0)
  
  !window function
  if (ss<s0) then
    Wap = 1d0
  else if (ss>=s0.and.ss<1d0) then
    Wap = x - dsin(2*pi*x)/(2*pi) 
  else 
    Wap = 0d0
  end if

end function WAP


subroutine window_sin(nn,D,W,ap,cut)
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  double precision, intent(in) :: D(1:2)
  double precision, intent(in), optional :: ap, cut
  double precision, intent(out) :: W(:)
  !internal
  integer :: i, j, n
  double precision :: a, c, xi, xj, sx ,sy

  sx = D(1)/dble(nn(1))
  sy = D(2)/dble(nn(2))

  a  = 1d0
  if (present(ap)) a = ap

  c  = 1d0
  if (present(cut)) c = cut

  n  = 1
  do i = 1, nn(1)
    do j = 1, nn(2)
      xi = dble(i-1-(nn(1)-1)*0.5)*sx
      xj = dble(j-1-(nn(2)-1)*0.5)*sy
      W(n) = Wap(D(1),abs(xi),a,c)*Wap(D(2),abs(xj),a,c)
      n = n + 1
    end do
  end do

end subroutine window_sin


subroutine maskpoint_generate(nn,D,Nm,px,py)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), Nm
  double precision, intent(in) :: D(2)
  double precision, intent(out) :: px(:), py(:)
  !internal
  integer :: i, j, n, m

  m = Nm
  write(*,*) 'generate masked points'
  call InitRandom(-1)
  n = 1
  do i = 1, nn(1)
    do j = 1, nn(2)
      if(m>0.or.dble(m)/dble(nn(1)*nn(2)-n)>ranmar()) then
        px(m) = dble(i-1-(nn(1)-1)*0.5)*D(1)/dble(nn(1))
        py(m) = dble(j-1-(nn(2)-1)*0.5)*D(2)/dble(nn(2))
        m = m - 1
      end if
      n = n + 1
    end do
  end do

end subroutine maskpoint_generate


function Mwin(mr,t,mfac) 
!* apodized mask window function
  implicit none
  !I/O
  double precision, intent(in) :: mr, t
  double precision, intent(in), optional :: mfac
  !internal
  double precision :: t0, Mwin, x, b

  b = mr
  if(present(mfac)) b = b*mfac

  if (t<mr) then
    Mwin = 1d0
  else if (t>=mr.and.t<b.and..not.b==mr) then
    x = (b-t)/(b-mr)
    Mwin = x - sin(2*pi*x)/(2*pi)
  else 
    Mwin = 0d0
  end if

end function Mwin


subroutine mask_generate(nn,D,px,py,mr,M)
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: D(2), mr, px(:), py(:)
  double precision, intent(inout) :: M(:)
  !internal
  integer :: i, j, p, n
  double precision :: xi, xj, dx, dy, sx ,sy

  sx = D(1)/dble(nn(1))
  sy = D(2)/dble(nn(2))

  n = 1
  do i = 1, nn(1)
    do j = 1, nn(2)
      do p = 1, size(px)
        xi = dble(i-1-(nn(1)-1)*0.5)*sx
        xj = dble(j-1-(nn(2)-1)*0.5)*sy
        dx = abs(xi-px(p))
        dy = abs(xj-py(p))
        M(n) = (1d0-Mwin(mr*sx,dx)*Mwin(mr*sy,dy))*M(n)
      end do
      n = n + 1
    end do
  end do

end subroutine mask_generate


subroutine window_generate(nn,D,W,mr,mn,f,ap,cut)
  implicit none
  !I/O
  character(*), intent(in), optional :: f
  integer, intent(in) :: nn(2)
  integer, intent(in), optional :: mn
  double precision, intent(in) :: D(2)
  double precision, intent(in), optional :: mr, ap, cut
  double precision, intent(inout) :: W(:)
  !internal
  integer :: i, j, n, maskn
  double precision :: a, c, xi, xj, dx, dy, sx ,sy
  double precision, allocatable :: px(:), py(:), M(:)

  sx = D(1)/dble(nn(1))
  sy = D(2)/dble(nn(2))

  a  = 1d0
  if (present(ap))  a = ap

  c  = 1d0
  if (present(cut)) c = cut

  call window_sin(nn,D,W,a,c)

  allocate(M(nn(1)*nn(2)));  M = 1d0
  if(present(mr).and..not.mr==0) then
    if(present(f)) then
      maskn = filelines(f)
      allocate(px(maskn),py(maskn))
      write(*,*) 'read masked points from a file, mask points =', maskn
      call loadtxt(f,px,py)
    else if(present(mn)) then
      maskn = mn
      allocate(px(maskn),py(maskn))
      call maskpoint_generate(nn,D,maskn,px,py)
      call savetxt('genpoints.dat',px,py)
    else
      stop 'error: require number of masks to be generated'
    end if
    call mask_generate(nn,D,px,py,mr,M)
    deallocate(px,py)
  end if

  W = W*M

  deallocate(M)

end subroutine window_generate


subroutine window_norm(W,Wn)
  implicit none
  double precision, intent(in) :: W(:)
  double precision, intent(out) :: Wn(:)
  integer :: n, npix

  npix = size(W)
  do n = 1, size(Wn)
    Wn(n) = sum(W**n)/dble(npix)
  end do

end subroutine window_norm


subroutine window_norm_multi(W1,W2,Wn)
  implicit none
  double precision, intent(in) :: W1(:), W2(:)
  double precision, intent(out) :: Wn(0:5,0:5)
  integer :: n, m, npix

  npix = size(W1)
  do n = 0, size(Wn,dim=1)-1
    do m = 0, size(Wn,dim=2)-1
      Wn(n,m) = sum(W1**n*W2**m)/dble(npix)
    end do
  end do

end subroutine window_norm_multi


subroutine gaussbeam_2d(t,els,beam,eL)
!* gaussain beam in 2D
  implicit none
  double precision, intent(in) :: els(:), t
  double precision, intent(out) :: beam(:)
  integer, intent(in), optional :: eL(1:2)
  !internal
  integer :: n, cL(2)

  cL = [0,int(maxval(abs(els)))+1]
  if (present(eL)) cL = eL

  do n = 1, size(els)
    if (cL(1)<=els(n).and.els(n)<=cL(2)) beam(n) = dexp(-els(n)**2*t**2/16d0/dlog(2d0))
  end do

end subroutine gaussbeam_2d


subroutine rotation(QU,rot,rtype)
  implicit none
  character(*), intent(in) :: rtype
  double precision, intent(in) :: rot(:)
  double precision, intent(inout) :: QU(:,:)
  integer :: npix
  double precision, allocatable :: tmp(:,:)  

  npix = size(rot)
  allocate(tmp(npix,2))
  select case(rtype)
  case('l')
    tmp(:,1) = QU(:,1) - QU(:,2)*2d0*rot
    tmp(:,2) = QU(:,2) + QU(:,1)*2d0*rot
  case('f')
    tmp(:,1) = QU(:,1)*dcos(2d0*dble(rot)) - QU(:,2)*dsin(2d0*dble(rot))
    tmp(:,2) = QU(:,2)*dcos(2d0*dble(rot)) + QU(:,1)*dsin(2d0*dble(rot))
  case default
    stop 'error (anaflat.f90/rotation): rotation type not specified'
  end select
  QU = tmp
  deallocate(tmp)

end subroutine rotation


end module anaflat

