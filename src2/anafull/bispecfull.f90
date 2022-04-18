!////////////////////////////////////////////////////!
! * Bispectrum in fullsky
!////////////////////////////////////////////////////!

module bispec_full
  use alm_tools, only: alm2map, map2alm
  use myconst,   only: dlc, pi

  private alm2map
  private dlc, pi

contains 


subroutine make_quad_gauss(nside,eL,alm)
  implicit none
  !I/O
  integer, intent(in) :: nside, eL(2)
  double complex, intent(inout), dimension(1,0:eL(2),0:eL(2)) :: alm
  !internal
  double precision, allocatable :: kmap(:)

  allocate(kmap(0:12*nside**2-1))
  call alm2map(nside,eL(2),eL(2),alm,kmap)
  kmap = kmap + kmap**2
  call map2alm(nside,eL(2),eL(2),kmap,alm)
 deallocate(kmap)

end subroutine make_quad_gauss


! for 3D bispectrum 

subroutine bispec_norm_bin(nside,bp,hl)
  implicit none
  integer, intent(in) :: nside
  double precision, intent(in) :: bp(:)
  double precision, intent(out) :: hl(:)
  integer :: l, lmax, bmax
  double complex, allocatable :: klm(:,:,:)
  double precision, allocatable :: fmap(:,:,:)

  lmax = int(maxval(bp))
  bmax = size(bp) - 1

  allocate(klm(1,0:lmax,0:lmax),fmap(1,bmax,0:12*nside**2-1)); klm=0d0; fmap=0d0

  do l = 1, lmax
    klm(1,l,0) = dsqrt(2d0*l+1d0)
  end do

  !first prepare filtered map for each bin
  call prep_filtered_map(nside,bp,klm,fmap)

  !next compute the product of filtered maps for each combination of b1, b2, b3
  call bispec_bin(fmap,hl)

  hl = hl/dsqrt(4d0*pi)

  deallocate(klm,fmap)

end subroutine bispec_norm_bin


subroutine create_bin_filter(bp,fl)
  implicit none
  double precision, intent(in) :: bp(:)
  double precision, intent(out) :: fl(:,0:,0:)
  integer :: l, b, bmax

  fl   = 0d0
  bmax = size(bp) - 1

  do b = 1, bmax
    do l = int(bp(b)), int(bp(b+1))
      fl(b,l,0:l) = 1d0
    end do
  end do

end subroutine create_bin_filter
 

subroutine prep_filtered_map(nside,bp,alm,fmap,fl0)
  implicit none
  !I/O
  integer, intent(in) :: nside
  double precision, intent(in) :: bp(:)
  double precision, intent(in), optional :: fl0(:,0:,0:)
  double complex, intent(in) :: alm(:,0:,0:)
  double precision, intent(out) :: fmap(:,:,0:)
  !internal
  integer :: b, s, lmax, ll, bmax, smax
  double precision, allocatable :: fl(:,:,:)
  double complex, allocatable :: falm(:,:,:)

  bmax = size(bp) - 1
  smax = size(alm,dim=1)
  lmax = int(maxval(bp))

  allocate(fl(bmax,0:lmax,0:lmax)); fl=0d0

  if (present(fl0)) then
    fl = fl0
  else
    call create_bin_filter(bp,fl)
  end if

  fmap = 0d0
  do b = 1, bmax
    do s = 1, smax
      ll = int(bp(b+1))
      allocate(falm(1,0:ll,0:ll)); falm=0d0
      falm(1,:,:) = fl(b,0:ll,0:ll)*alm(s,0:ll,0:ll)
      call alm2map(nside,ll,ll,falm(1:1,:,:),fmap(s,b,:))
      deallocate(falm)
    end do
  end do

  deallocate(fl)

end subroutine prep_filtered_map


subroutine bispec_bin(kmap,Bb)
  !unnormalized binned bispectrum for 3D shape
  implicit none
  !I/O
  double precision, intent(in) :: kmap(:,:,0:)  !multipole-filtered map
  double precision, intent(out) :: Bb(:)
  !internal
  integer :: b, b0, b1, b2, kn, bmax
  double precision :: Bb111, Bb112, Bb121, Bb211

  kn   = size(kmap,dim=1)
  bmax = size(kmap,dim=2)

  b = 0
  do b0 = 1, bmax
    do b1 = b0, bmax
      do b2 = b1, bmax
        b = b + 1
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
  !Bb = Bb*(4d0*pi)/(12d0*dble(nside)**2)

end subroutine bispec_bin




! for 1D bispectrum 

subroutine bispec_equi_norm(bp0,bp1,hl)
  implicit none
  double precision, intent(in) :: bp0(:), bp1(:)
  double precision, intent(out) :: hl(:)
  integer :: bmax, eL(2), lmax, l, b
  double complex, allocatable :: klm(:,:,:)

  bmax = size(hl)
  lmax = int(bp1(bmax))

  allocate(klm(1,0:lmax,0:lmax)); klm=0d0

  do l = 1, lmax
    klm(1,l,0) = dsqrt(2d0*l+1d0)
  end do

  do b = 1, bmax
    write(*,*) b
    eL = [int(bp0(b)),int(bp1(b))]
    call bispec_equi(eL,klm(1,0:eL(2),0:eL(2)),hl(b))
    hl(b) = hl(b)/dsqrt(4d0*pi)
  end do

  deallocate(klm)

end subroutine bispec_equi_norm


! this is not used for 3D bispectrum
subroutine bispec_norm(bp0,bp1,bp2,hl)
  implicit none
  double precision, intent(in) :: bp0(:), bp1(:), bp2(:)
  double precision, intent(out) :: hl
  integer :: lmax, l, b
  double complex, allocatable :: klm(:,:,:)

  lmax = int(maxval([bp0(2),bp1(2),bp2(2)]))

  allocate(klm(3,0:lmax,0:lmax)); klm=0d0

  do l = 1, lmax
    if (bp0(1)<=l.and.l<bp0(2))  klm(1,l,0) = dsqrt(2d0*l+1d0)
    if (bp1(1)<=l.and.l<bp1(2))  klm(2,l,0) = dsqrt(2d0*l+1d0)
    if (bp2(1)<=l.and.l<bp2(2))  klm(3,l,0) = dsqrt(2d0*l+1d0)
  end do

  call bispec_equi_x([1,lmax],klm(1,:,:),klm(2,:,:),klm(3,:,:),hl)
  hl = hl/dsqrt(4d0*pi)

  deallocate(klm)

end subroutine bispec_norm


subroutine bispec_equi_x(eL,alm1,alm2,alm3,bispec)
  implicit none
  !I/O
  integer, intent(in) :: eL(2)
  double complex, intent(in), dimension(0:eL(2),0:eL(2)) :: alm1, alm2, alm3
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside, lmax
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  lmax = eL(2)
  nside = 2**(int(dlog(dble(lmax))/dlog(2d0)))

  allocate(kmap(3,0:12*nside**2-1),klm(3,0:lmax,0:lmax)); klm=0d0

  klm(1,:,:) = alm1
  klm(2,:,:) = alm2
  klm(3,:,:) = alm3

  call alm2map(nside,eL(2),eL(2),klm(1:1,:,:),kmap(1,:))
  call alm2map(nside,eL(2),eL(2),klm(2:2,:,:),kmap(2,:))
  call alm2map(nside,eL(2),eL(2),klm(3:3,:,:),kmap(3,:))

  !unnormalized
  bispec = sum(kmap(1,:)*kmap(2,:)*kmap(3,:))*(4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine bispec_equi_x


subroutine bispec_equi(eL,alm,bispec)
  implicit none
  !I/O
  integer, intent(in) :: eL(2)
  double complex, intent(in), dimension(0:eL(2),0:eL(2)) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:)
  double complex, allocatable :: klm(:,:,:)

  nside = 2**(int(dlog(dble(eL(2)))/dlog(2d0)))

  allocate(kmap(0:12*nside**2-1),klm(1,0:eL(2),0:eL(2)))

  klm = 0d0
  do l = eL(1), eL(2) !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,eL(2),eL(2),klm,kmap)
  bispec = sum(kmap**3)*(4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine bispec_equi


subroutine bispec_fold(eL,bst,alm,bispec)
  implicit none
  !I/O
  integer, intent(in) :: eL(2), bst
  double complex, intent(in), dimension(0:eL(2),0:eL(2)) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  nside = eL(2)*bst

  allocate(kmap(0:12*nside**2-1,2),klm(2,0:eL(2),0:eL(2)))

  klm = 0d0
  do l = eL(1), eL(2) !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do
  do l = max(2,int(eL(1)/2d0)), int(eL(2)/2d0)
    klm(2,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,eL(2),eL(2),klm(1:1,:,:),kmap(:,1))
  call alm2map(nside,eL(2),eL(2),klm(2:2,:,:),kmap(:,2))
  bispec = sum(kmap(:,1)*kmap(:,2)**2) * (4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine bispec_fold


subroutine bispec_sque(eL,sL,bst,alm,bispec)
  implicit none
  !I/O
  integer, intent(in) :: eL(2), sL(2), bst
  double complex, intent(in), dimension(0:eL(2),0:eL(2)) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside, lmax
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  lmax = max(eL(2),sL(2))

  nside = lmax*bst

  allocate(kmap(0:12*nside**2-1,2),klm(2,0:lmax,0:lmax))

  klm = 0d0
  do l = sL(1), sL(2) !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do
  do l = eL(1), eL(2)
    klm(2,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,lmax,lmax,klm(1:1,:,:),kmap(:,1))
  call alm2map(nside,lmax,lmax,klm(2:2,:,:),kmap(:,2))
  bispec = sum(kmap(:,1)*kmap(:,2)**2) * (4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine bispec_sque


subroutine bispec_angl(eL,aL,l1,bst,alm,bispec)
  implicit none
  !I/O
  integer, intent(in) :: eL(2), aL(2), bst, l1
  double complex, intent(in), dimension(0:l1,0:l1) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  nside = l1*bst

  allocate(kmap(0:12*nside**2-1,2),klm(2,0:l1,0:l1))

  klm = 0d0
  do l = eL(1), eL(2) !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do
  do l = aL(1), aL(2)
    klm(2,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,l1,l1,klm(1:1,:,:),kmap(:,1))
  call alm2map(nside,l1,l1,klm(2:2,:,:),kmap(:,2))
  bispec = sum(kmap(:,1)*kmap(:,2)**2) * (4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine bispec_angl


end module bispec_full


