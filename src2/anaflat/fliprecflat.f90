!///////////////////////////////////////////////////////////!
! * Spin-Flip (EB mixing) Reconstruction Kernel
!///////////////////////////////////////////////////////////!

module fliprecflat
  use myconst, only: iu, i1, pi, dlc, twopi
  use anaflat, only: elarrays
  use myfftw,  only: dft
  implicit none

  private iu, i1, pi, dlc, twopi
  private elarrays
  private dft

contains 


subroutine quadeb_flip(nn,D,E,B,EE,rL,eL,fe,fb)
!
  implicit none
! [inputs]
!   nn(1:2) --- 2D grids
!   rL(1:2) --- multipole range of reconstruction
!   D(1:2)  --- map size (radian)
!   EE(:)   --- E-mode Cl in 2D grids
!   E(:)    --- filtered E-modes
!   B(:)    --- filtered B-modes
!   eL(1:2) --- multipole range of output
  integer, intent(in) :: nn(1:2), rL(1:2)
  double precision, intent(in) :: EE(:), D(1:2)
  complex(dlc), intent(in) :: E(:), B(:)
  integer, intent(in), optional :: eL(1:2)
!
! [output]
!   fe(:), fb(:) --- reconstructed fields
  complex(dlc), intent(out), optional :: fe(:), fb(:)
!
! [internal]
  integer :: i, n, npix
  double precision, allocatable :: l(:,:), els(:)
  complex(dlc), allocatable :: wB(:), alm(:,:), wE(:,:), ei2p(:)

  npix = nn(1)*nn(2)
  allocate(wB(npix),wE(2,npix),ei2p(npix),l(2,npix),els(npix));  wB=0d0;  wE=0d0
  call elarrays(nn,D,elx=l(1,:),ely=l(2,:),els=els,ei2p=ei2p)

  if (size(EE)/=npix)  stop 'error (quadeb): size of EE /= npix.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    wB(n)   = B(n)*conjg(ei2p(n))
    wE(:,n) = l(:,n)*EE(n)*E(n)*ei2p(n)
  end do
  deallocate(ei2p)

  !* convolution
  call dft(wB,nn,D,-1)
  call dft(wE(1,:),nn,D,-1)
  call dft(wE(2,:),nn,D,-1)
  allocate(alm(2,npix))
  alm(1,:) = dble(wB*wE(1,:))
  alm(2,:) = dble(wB*wE(2,:))
  deallocate(wB,wE)
  call dft(alm(1,:),nn,D,1)
  call dft(alm(2,:),nn,D,1)

  !* estimator 
  do n = 1, npix
    if(present(eL).and.(els(n)<eL(1).or.els(n)>eL(2))) alm(:,n)=0d0
    if(present(fe)) fe(n) = sum(l(:,n)*alm(:,n))/iu
    if(present(fb)) fb(n) = (l(2,n)*alm(1,n)-l(1,n)*alm(2,n))/iu
  end do

  deallocate(alm,l,els)

end subroutine quadeb_flip


subroutine alflat_eb_flip(nn,D,fE,fB,CE,rL,eL,Alg)
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2), el(1:2), rL(1:2)
  double precision, intent(in) :: CE(:), fE(:), fB(:), D(1:2)
  double precision, intent(out), optional :: Alg(:)
! [internal]
  integer :: n, npix
  double precision :: ll(2,nn(1)*nn(2)), els(nn(1)*nn(2)), iAlg, iAlc
  complex(dlc) :: vec(1:2)
  complex(dlc), allocatable :: Al(:,:), A1(:,:), A2(:,:), ei2p(:)

  write(*,*) 'EB flat pr'

  npix = nn(1)*nn(2)
  allocate(A1(3,npix),A2(3,npix),ei2p(npix));  A1=0d0;  A2=0d0
  call elarrays(nn,D,elx=ll(1,:),ely=ll(2,:),els=els,ei2p=ei2p)

  if (size(CE)/=npix)  stop 'error (alflat_eb_flip): CE is strange.'

  !* filtering
  do n = 1, npix
    if(rL(1)>els(n).or.els(n)>rL(2)) cycle
    A1(1,n)  = CE(n)**2*fE(n)*ei2p(n)**2
    A1(2,n)  = CE(n)**2*fE(n)
    A2(1,n)  = - fB(n)*conjg(ei2p(n)**2)
    A2(2,n)  = fB(n)
  end do
  deallocate(ei2p)

  !* convolution
  allocate(Al(2,npix))
  do n = 1, 2
    call dft(A1(n,:),nn,D,-1)
    call dft(A2(n,:),nn,D,-1)
    Al(n,:) = A1(n,:)*A2(n,:)
    if (n==1)  Al(n,:) = dble(Al(n,:))
    call dft(Al(n,:),nn,D,1)
  end do
  deallocate(A1,A2)

  !* normalization
  do n = 1, npix
    if(els(n)<eL(1).or.els(n)>eL(2))  cycle
    iAlg = sum(Al(1:2,n))*0.5d0
    if(iAlg<=0d0) cycle
    if (present(Alg))  Alg(n) = 1d0 / iAlg
  end do

  deallocate(Al)

end subroutine alflat_eb_flip


end module fliprecflat


