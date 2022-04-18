!////////////////////////////////////////////////////!
! * bispectrum measurement from flatsky simulations
!////////////////////////////////////////////////////!

module bispflat
  use myconst, only: dlc
  use anaflat, only: elarray
  use myfftw,  only: dft

  private dlc
  private elarray
  private dft

contains


subroutine binfilter(nn,D,bmin,bmax,bf)
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: bmin(:), bmax(:), D(2)
  double precision, intent(out) :: bf(:,:)
  integer :: n, b, npix, bnum
  double precision, allocatable :: els(:)

  npix = nn(1)*nn(2)
  bnum = size(bmin)

  bf = 0d0

  allocate(els(npix))
  els = elarray(nn,D)

  do b = 1, bnum
    do n = 1, npix
      if (bmin(b)<=els(n).and.els(n)<bmax(b))  bf(b,n) = 1d0
    end do
  end do

  deallocate(els)

end subroutine binfilter


! for 3D bispectrum

subroutine bispec_norm_bin(nn,D,bp,dbmax,hl)
  implicit none
  integer, intent(in) :: nn(2), dbmax
  double precision, intent(in) :: D(2), bp(:)
  double precision, intent(out) :: hl(:)
  integer :: l, lmax, bmax, npix
  double complex, allocatable :: klm(:,:)
  double precision, allocatable :: fmap(:,:,:)

  npix = nn(1)*nn(2)
  bmax = size(bp) - 1

  allocate(klm(1,npix),fmap(1,bmax,npix)); klm=1d0; fmap=0d0

  !first prepare filtered map for each bin
  call prep_filtered_map(nn,D,bp,klm,fmap)

  !next compute the product of filtered maps for each combination of b1, b2, b3
  call bispec_bin(fmap,bp,dbmax,hl,1d0)

  hl = hl*D(1)*D(2)

  deallocate(klm,fmap)

end subroutine bispec_norm_bin


subroutine prep_filtered_map(nn,D,bp,alm,fmap,fl0)
  implicit none
  !I/O
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: bp(:), D(2)
  double precision, intent(in), optional :: fl0(:,:)
  double complex, intent(in) :: alm(:,:)
  double precision, intent(out) :: fmap(:,:,:)
  !internal
  integer :: npix, b, s, bmax, smax
  double precision, allocatable :: fl(:,:)
  double complex, allocatable :: falm(:)

  npix = nn(1)*nn(2)
  bmax = size(bp) - 1
  smax = size(alm,dim=1)

  ! create bin filter
  allocate(fl(bmax,npix)); fl=0d0
  if (present(fl0)) then
    fl = fl0
  else
    call binfilter(nn,D,bp(1:bmax),bp(2:bmax+1),fl)
  end if

  ! compute filtered map
  fmap = 0d0
  do b = 1, bmax
    do s = 1, smax
      allocate(falm(npix)); falm=0d0
      falm = fl(b,:)*alm(s,:)
      call dft(falm,nn,D,-1)
      fmap(s,b,:) = dble(falm)
      deallocate(falm)
    end do
  end do

  deallocate(fl)

end subroutine prep_filtered_map


subroutine bispec_bin(kmap,bp,dbmax,Bb,val)
  !unnormalized binned bispectrum for 3D shape
  implicit none
  !I/O
  integer, intent(in) :: dbmax
  double precision, intent(in) :: kmap(:,:,:)  !multipole-filtered map
  double precision, intent(in) :: bp(:)
  double precision, intent(in), optional :: val
  double precision, intent(out) :: Bb(:)
  !internal
  integer :: b, b0, b1, b2, kn, bmax
  double precision :: lbmin, lbmax
  double precision :: Bb111, Bb112, Bb121, Bb211

  kn   = size(kmap,dim=1)
  bmax = size(kmap,dim=2)

  if (present(val)) then 
    Bb = val
  else
    Bb = 0d0
  end if

  b = 0

  do b0 = 1, bmax
    do b1 = b0, bmax
      if (b1>b0+dbmax) cycle !discard far separate bins 
      !assuming bp(b0)<bp(b1)
      lbmin = bp(b1) - bp(b0+1)
      lbmax = bp(b1+1) + bp(b0+1)
      do b2 = b1, bmax
        if (b2>b1+dbmax) cycle !discard far separate bins
        b = b + 1
        if (bp(b2)<lbmin.or.lbmax<bp(b2+1)) cycle !triangle condition
        select case(kn)
        case(1)
          Bb111 = sum(kmap(1,b0,:)*kmap(1,b1,:)*kmap(1,b2,:))
          Bb(b) = Bb111
        case(2)
          Bb112 = sum(kmap(1,b0,:)*kmap(1,b1,:)*kmap(2,b2,:))
          Bb121 = sum(kmap(1,b0,:)*kmap(2,b1,:)*kmap(1,b2,:))
          Bb211 = sum(kmap(2,b0,:)*kmap(1,b1,:)*kmap(1,b2,:))
          Bb(b) = Bb112 + Bb121 + Bb211
        case default
          stop 'only kn=1 or 2 types are supported'
        end select
      end do
    end do
  end do

end subroutine bispec_bin



! for 1D or specific bispectrum

subroutine norm_fold(nn,D,bf,dlb)
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: bf(:,:,:), D(2)
  double precision, intent(out) :: dlb(:)
  integer :: b, npix
  complex(dlc), allocatable :: wlm(:,:)

  npix = nn(1)*nn(2)

  do b = 1, size(bf,dim=2)
    allocate(wlm(2,npix));  wlm=0d0
    !write(*,*) b
    wlm(1,:) = bf(1,b,:)
    wlm(2,:) = bf(2,b,:)
    call dft(wlm(1,:),nn,D,-1)
    call dft(wlm(2,:),nn,D,-1)
    dlb(b) = D(1)*D(2)*sum(wlm(1,:)*wlm(2,:)**2)
    deallocate(wlm)
  end do

end subroutine norm_fold


subroutine bisp_fold(nn,D,bf,dlb,klm,Bb)
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: bf(:,:,:), D(2), dlb(:)
  double precision, intent(out) :: Bb(:)
  complex(dlc), intent(in) :: klm(:)
  integer :: b, npix
  complex(dlc), allocatable :: wklm(:,:)

  npix = nn(1)*nn(2)
  
  do b = 1, size(bf,dim=2)
    allocate(wklm(2,npix));  wklm=0d0
    !write(*,*) 'binned kappa', b
    wklm(1,:) = bf(1,b,:)*klm
    wklm(2,:) = bf(2,b,:)*klm
    call dft(wklm(1,:),nn,D,-1)
    call dft(wklm(2,:),nn,D,-1)
    Bb(b) = sum(wklm(1,:)*wklm(2,:)**2)/dlb(b)
    deallocate(wklm)
  end do

end subroutine bisp_fold


subroutine norm_equi(nn,D,bf,dlb)
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: bf(:,:), D(2)
  double precision, intent(out) :: dlb(:)
  integer :: b, npix
  complex(dlc), allocatable :: wlm(:)

  npix = nn(1)*nn(2)

  do b = 1, size(bf,dim=1)
    allocate(wlm(npix));  wlm=0d0
    wlm = bf(b,:)
    call dft(wlm,nn,D,-1)
    dlb(b) = D(1)*D(2)*sum(wlm**3) !bispectrum 
    deallocate(wlm)
  end do

end subroutine norm_equi


subroutine bisp_equi(nn,D,bf,dlb,klm,Bb)
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: bf(:,:), D(2), dlb(:)
  double precision, intent(out) :: Bb(:)
  complex(dlc), intent(in) :: klm(:)
  integer :: b, npix
  complex(dlc), allocatable :: wklm(:)

  npix = nn(1)*nn(2)
  
  do b = 1, size(bf,dim=1)
    allocate(wklm(npix));  wklm=0d0
    wklm = bf(b,:)*klm
    call dft(wklm,nn,D,-1)
    Bb(b) = sum(wklm**3)/dlb(b)
    deallocate(wklm)
  end do

end subroutine bisp_equi


subroutine bisp_equi_x(nn,D,bf,dlb,klm1,klm2,klm3,Bb)
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: bf(:,:), D(2), dlb(:)
  double precision, intent(out) :: Bb(:)
  complex(dlc), intent(in) :: klm1(:), klm2(:), klm3(:)
  integer :: i, b, npix
  complex(dlc), allocatable :: wklm(:,:)

  npix = nn(1)*nn(2)
  
  do b = 1, size(bf,dim=1)
    allocate(wklm(3,npix));  wklm=0d0
    wklm(1,:) = bf(b,:)*klm1
    wklm(2,:) = bf(b,:)*klm2
    wklm(3,:) = bf(b,:)*klm3
    do i = 1, 3
      call dft(wklm(i,:),nn,D,-1)
    end do
    Bb(b) = sum(wklm(1,:)*wklm(2,:)*wklm(3,:))/dlb(b)
    deallocate(wklm)
  end do

end subroutine bisp_equi_x


end module bispflat

