!////////////////////////////////////////////////////!
! * Tools for Fullsky Analysis
!////////////////////////////////////////////////////!

module anafull
  !* healpix modules
  use alm_tools, only: alm2map, alm2map_spin, map2alm, map2alm_spin, alm2map_der
  use fitstools, only: input_map, output_map
  use pix_tools, only: pix2ang_ring, ang2pix_ring
  !* my modules
  use myrandom, only: InitRandom, Gaussian1
  use myconst, only: pi, iu
  use myutils, only: str

  implicit none

  private initrandom, gaussian1
  private alm2map, alm2map_spin, map2alm, map2alm_spin, alm2map_der
  private input_map, output_map
  private pix2ang_ring, ang2pix_ring
  private pi, iu
  private str

contains


subroutine gaussianalm(alm,Cl,lmax)
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double precision, intent(in) :: Cl(:)
  double complex, intent(out) :: alm(0:lmax,0:lmax)
  !internal
  integer :: l,m
  
  call INITRANDOM(-1) 
  alm = 0
  do l = 1, lmax
    alm(l,0) = Gaussian1()* sqrt(Cl(l))
    do m = 1, l
      alm(l,m) = cmplx(Gaussian1(),Gaussian1())*sqrt(Cl(l)/2)
    end do 
  end do

end subroutine gaussianalm


subroutine gaussianTE(alm,TT,EE,TE,lmax)
  !* include correlation between T and E
  !* order of input Cl should be TT, EE, TE
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double precision, intent(in) :: TT(:), EE(:), TE(:)
  double complex, intent(out) :: alm(2,0:lmax,0:lmax)
  !internal
  integer :: l, m
  double precision :: tamp, corr, xamp

  alm = 0
  call gaussianalm(alm(1,:,:),TT,lmax)

  call initrandom(-1)
  do l = 2, lmax
    tamp = TT(l)
    if(TT(l)==0d0) tamp=1d0  !prevent divide by zero
    corr = TE(l)/tamp
    xamp = sqrt(EE(l) - corr*TE(l))
    alm(2,l,0) = corr*alm(1,l,0) + Gaussian1()*xamp
    xamp = xamp/sqrt(2d0)
    do m =1, l
      alm(2,l,m) = corr*alm(1,l,m) + cmplx(Gaussian1(),Gaussian1())*xamp
    end do
  end do

end subroutine gaussianTE


subroutine gaussianTEB(alm,TT,EE,BB,TE,lmax)
  !* include correlation between T and E
  !* order of input Cl should be TT, EE, BB, TE
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double precision, intent(in) :: TT(:), EE(:), BB(:), TE(:)
  double complex, intent(out) :: alm(3,0:lmax,0:lmax)
  !internal
  integer :: l, m
  double precision :: tamp, corr, xamp

  alm = 0
  call gaussianalm(alm(1,:,:),TT,lmax)
  call gaussianalm(alm(3,:,:),BB,lmax)

  call initrandom(-1)
  do l = 2, lmax
    tamp = TT(l)
    if(TT(l)==0d0) tamp=1d0  !prevent divide by zero
    corr = TE(l)/tamp
    xamp = sqrt(EE(l) - corr*TE(l))
    alm(2,l,0) = corr*alm(1,l,0) + Gaussian1()*xamp
    xamp = xamp/sqrt(2d0)
    do m =1, l
      alm(2,l,m) = corr*alm(1,l,m) + cmplx(Gaussian1(),Gaussian1())*xamp
    end do
  end do

end subroutine gaussianTEB


subroutine gaussianEB(alm,Cl,lmax)
  !* for any two independent cls
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double precision, intent(in) :: Cl(:,:)
  double complex, intent(out) :: alm(2,0:lmax,0:lmax)

  if(size(Cl,dim=1)<2) stop 'need more Cls'
  call gaussianalm(alm(1,:,:),Cl(1,:),lmax)
  call gaussianalm(alm(2,:,:),Cl(2,:),lmax)  

end subroutine gaussianEB


subroutine addnoiseEB(alm,Nl,lmax)
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double precision, intent(in) :: Nl(:)
  double complex, intent(inout), dimension(2,0:lmax,0:lmax) :: alm
  !internal
  double complex :: nlm(0:lmax,0:lmax)

  call gaussianalm(nlm,Nl,lmax)
  alm(1,:,:) = alm(1,:,:) + nlm(:,:)
  alm(2,:,:) = alm(2,:,:) + nlm(:,:)

end subroutine addnoiseEB


subroutine addnoisealm(alm,Nl,lmax)
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double precision, intent(in) :: Nl(:)
  double complex, intent(inout) :: alm(0:lmax,0:lmax)
  !internal
  double complex :: nlm(0:lmax,0:lmax)

  call gaussianalm(nlm,Nl,lmax)
  alm = alm + nlm

end subroutine addnoisealm


subroutine read_maps(map,f)
  implicit none
  !I/O
  character(*), intent(in) :: f(:)
  double precision, intent(out) :: map(:,:)
  !internal
  integer :: i, nmap, npix

  npix = size(map,dim=1)
  nmap = size(map,dim=2)

  if(.not.size(f)==nmap) stop 'error: incorrect num of filenames'

  do i = 1, nmap
    call read_map(npix,map(0:npix-1,i),f(i))
  end do

end subroutine read_maps


subroutine read_map(npix,map,f)
  implicit none
  !I/O
  integer, intent(in) :: npix
  character(*), intent(in) :: f
  double precision, intent(out) :: map(0:npix-1)
  !internal
  double precision, allocatable, dimension(:,:) :: imap

  allocate(imap(0:npix-1,1))
  call INPUT_MAP(f,imap,npix,1)
  map(0:npix-1) = imap(0:npix-1,1)
  deallocate(imap)

end subroutine read_map


subroutine read_qumap(npix,Pobs,root,simn,amp,add,pf)
  implicit none
  !I/O
  character(*), intent(in) :: root
  character(*), intent(in), optional :: pf
  logical, intent(in), optional :: add
  integer, intent(in), optional :: simn
  double precision, intent(in), optional :: amp
  integer, intent(in) :: npix
  double precision, intent(inout) :: Pobs(0:npix-1,2)
  !internal
  character(LEN=64) :: f, fs(2)
  double precision, allocatable, dimension(:,:) :: map

  if(.not.(present(amp).and.amp==0d0)) then
    allocate(map(0:npix-1,2))
    f = 'map'
    if(present(pf)) f = pf
    fs(1) = trim(root)//'Q'//trim(f)//'_r'//str(simn)//'.fits'
    fs(2) = trim(root)//'U'//trim(f)//'_r'//str(simn)//'.fits'
    call read_maps(map,fs)
    if(present(amp)) map = amp*map
    if(present(add)) then
      Pobs = Pobs + map
    else
      Pobs(0:npix-1,:) = map(0:npix-1,:)
    end if
    deallocate(map)
  end if

end subroutine read_qumap


subroutine READ_TMAP(Tobs,root,simn,amp,add,pf)
  implicit none
  !I/O
  character(*), intent(in) :: root
  character(*), intent(in), optional :: pf
  logical, intent(in), optional :: add
  integer, intent(in), optional :: simn
  double precision, intent(in), optional :: amp
  double precision, intent(inout) :: Tobs(:)
  !internal
  integer :: npix
  character(LEN=64) :: f
  double precision, allocatable, dimension(:) :: map

  npix = size(Tobs)
  if(.not.(present(amp).and.amp==0d0)) then
    allocate(map(0:size(Tobs)-1))
    f = 'map'
    if(present(pf)) f = pf
    call read_map(npix,map,trim(root)//'T'//trim(f)//'_r'//str(simn)//'.fits')
    if(present(amp)) map = amp*map
    if(present(add)) then
      Tobs = Tobs + map
    else
      Tobs = map
    end if
    deallocate(map)
  end if

end subroutine READ_TMAP


subroutine output_alm2map(f,ialm,lmax,nside)
  implicit none
  !I/O
  character(*), intent(in) :: f
  integer, intent(in) :: lmax, nside
  double complex, intent(in) :: ialm(0:lmax,0:lmax)
  !internal
  character(80) :: header(10)
  double precision :: map(0:12*nside**2-1,1)
  double complex :: alm(1,0:lmax,0:lmax)

  header = ' '
  alm(1,:,:) = ialm
  call alm2map(nside,lmax,lmax,alm,map(:,1))
  call output_map(map,header,f)

end subroutine output_alm2map


subroutine map_output(npix,f,imap)
  implicit none
  !I/O
  integer, intent(in) :: npix
  character(*), intent(in) :: f
  double precision, intent(in) :: imap(0:npix-1)
  !internal
  character(80) :: header(10)
  double precision :: map(0:npix-1,1)

  header = ' '
  map(:,1) = imap
  call output_map(map,header,f)

end subroutine map_output


subroutine map_output_p(npix,P,simn,pf)
  implicit none
  !I/O
  character(*), intent(in), optional :: pf
  integer, intent(in) :: npix
  integer, intent(in), optional :: simn
  double precision, intent(in) :: P(0:npix-1,2)
  !internal
  character(LEN=32) :: f

  f = 'map'
  if(present(pf)) f = pf
  call map_output(npix,'Q'//trim(f)//'_r'//str(simn)//'.fits',P(:,1))
  call map_output(npix,'U'//trim(f)//'_r'//str(simn)//'.fits',P(:,2))

end subroutine map_output_p


subroutine SIM_QUMAP(Pmap,Cl,lmax,nside)
  implicit none
  !order of cl shoud be TT, EE, BB, TE
  !I/O
  integer, intent(in) :: lmax,nside
  double precision, intent(in) :: Cl(:,:)
  double precision, intent(out) :: Pmap(:,:)
  !internal
  integer :: npix
  double precision, allocatable, dimension(:,:) :: map
  double complex, allocatable, dimension(:,:,:) :: alm

  npix = size(Pmap,dim=1)

  allocate(alm(3,0:lmax,0:lmax),map(0:npix-1,1:3))
  call gaussianTEB(alm,Cl(1,:),Cl(2,:),Cl(3,:),Cl(4,:),lmax)
  call alm2map(nside,lmax,lmax,alm,map)
  Pmap(:,1) = map(:,2)
  Pmap(:,2) = map(:,3)
  deallocate(alm,map)

end subroutine SIM_QUMAP


subroutine SIM_QUMAP_NOISE(Pmap,Nl,lmax,nside)
  implicit none
  !I/O
  integer, intent(in) :: lmax,nside
  double precision, intent(in) :: Nl(:)
  double precision, intent(inout) :: Pmap(:,:)
  !internal
  integer :: npix
  double complex, allocatable, dimension(:,:,:) :: alm

  npix = size(Pmap,dim=1)

  allocate(alm(2,0:lmax,0:lmax))
  call gaussianalm(alm(1,:,:),Nl,lmax)
  call gaussianalm(alm(2,:,:),Nl,lmax)
  call alm2map_spin(nside,lmax,lmax,2,alm,Pmap)
  deallocate(alm)

end subroutine SIM_QUMAP_NOISE


subroutine BASELINE_SUBTRACTION(BLmap,Pmap,nside,nside_subpatch)
  !//// Author: Ryo Nagata ////!
  use pix_tools
  implicit none
  !I/O
  integer, intent(in) :: nside,nside_subpatch
  double precision, intent(inout) :: Pmap(:,:)
  double precision, intent(out) :: BLmap(:,:)
  !internal
  integer :: i,j, npix_subpatch, npix_inpatch, npix
  integer, allocatable, dimension(:) :: patch_map
  double precision, allocatable, dimension(:,:) :: patch_pmean

  npix = size(Pmap,dim=1)
  npix_subpatch = nside2npix(nside_subpatch)
  npix_inpatch  = npix/nside2npix(nside_subpatch)

  allocate(patch_map(0:npix-1))
  allocate(patch_pmean(0:npix_subpatch-1,2))
  patch_pmean = 0d0

  !* compute baseline *!
  do i=0, npix-1
    call get_subpatch( nside_subpatch, nside, i, j)
    patch_map(i) = j
    patch_pmean(j,:) = patch_pmean(j,:) + Pmap(i,:)/dble(npix_inpatch)
  end do

  !* baseline subtraction *!
  do i=0, npix-1
    Pmap(i,:) = Pmap(i,:) - patch_pmean(patch_map(i),:)
    BLmap(i,:) = patch_pmean(patch_map(i),:)
  end do

  deallocate(patch_map,patch_pmean)

end subroutine BASELINE_SUBTRACTION


subroutine GET_SUBPATCH( nside_large, nside_small, ipix_in, ipix_out)
  !//// Written by Ryo Nagata ////!
  implicit none
  !I/O
  integer, intent(in) :: nside_large, nside_small, ipix_in
  integer, intent(out) :: ipix_out
  !internal
  double precision :: theta, phi

  call pix2ang_ring(nside_small, ipix_in, theta, phi)
  call ang2pix_ring(nside_large, theta, phi, ipix_out)

end subroutine GET_SUBPATCH


subroutine GET_WINMAP(nside_large,nside_small,ipix_pix,win_out)
  !//// Written by Ryo Nagata ////!
  implicit none
  !I/O
  integer, intent(in) :: nside_large, nside_small,ipix_pix
  double precision, intent(out) :: win_out
  !internal
  integer :: ipix_sub
  double precision :: dtheta, theta_pix, phi_pix, theta_sub, phi_sub

  call pix2ang_ring(nside_small, ipix_pix, theta_pix, phi_pix)
  call ang2pix_ring(nside_large, theta_pix, phi_pix, ipix_sub)
  call pix2ang_ring(nside_large, ipix_sub, theta_sub, phi_sub)

  dtheta = dacos(dcos(theta_pix)*dcos(theta_sub) &
      + dsin(theta_pix)*dcos(phi_pix)*dsin(theta_sub)*dcos(phi_sub)  &
      + dsin(theta_pix)*dsin(phi_pix)*dsin(theta_sub)*dsin(phi_sub) )
  call get_window( dtheta, win_out)

end subroutine GET_WINMAP


subroutine GET_WINDOW(s,w)
  !//// Written by Ryo Nagata ////!
  implicit none
  !I/O
  double precision, intent(in) ::s
  double precision, intent(out) :: w
  !internal
  double precision :: s_in,s_out,s0

  s_out = 1d0
  s0    = 1d0
  s_in  = s_out*s0
  if ( s < s_in  ) then 
    w=1d0
  else if ( s > s_out ) then
    w=0d0
  else
    w=(s_out-s)/(s_out-s_in)-dsin(2d0*pi*(s_out-s)/(s_out-s_in))/2d0/pi
  end if

end subroutine GET_WINDOW


subroutine cmplxmap2alm_spin(nside,lmax,mmax,spin,map,alm)
  implicit none
  !I/O
  integer, intent(in) :: nside, lmax, mmax, spin
  double complex, intent(in) :: map(0:12*nside**2-1)
  double complex, intent(out) :: alm(2,0:lmax,0:mmax)
  !internal
  double precision :: S(0:12*nside**2-1,2)

  S(:,1) = real(map)
  S(:,2) = aimag(map)
  call map2alm_spin(nside,lmax,lmax,spin,S,alm)

end subroutine cmplxmap2alm_spin


subroutine alm_to_map(nside,lmax,mmax,alm,map)
  implicit none
  !I/O
  integer, intent(in) :: nside, lmax, mmax
  double precision, intent(out) :: map(0:12*nside**2-1)
  double complex, intent(in) :: alm(0:lmax,0:lmax)
  !internal
  double complex :: tlm(1,0:lmax,0:lmax)

  tlm(1,:,:) = alm
  call alm2map(nside,lmax,mmax,tlm,map)

end subroutine alm_to_map


function cosin_healpix(nside,lmax,mmax) result(cosin)
  implicit none
  integer, intent(in) :: nside, lmax, mmax
  double complex :: alm(1,0:lmax,0:mmax)
  double precision :: cosin(0:12*nside**2-1)

  alm = 0d0
  alm(1,1,0) = 1d0
  cosin = 0d0
  !call alm_to_map(nside,lmax,mmax,alm,cosin)
  call alm2map(nside,lmax,mmax,alm,cosin)
  cosin = cosin*dsqrt(pi/0.75d0)

end function cosin_healpix


subroutine elm2map_spin(nside,lmax,mmax,spin,Elm,map)
! map = S1 + i*S2 = sum_{lm} E_{lm} Y_{lm} 
! (E_{lm}=a^+_{lm})
  implicit none
  !I/O
  integer, intent(in) :: nside, lmax, mmax, spin
  double complex, intent(in) :: Elm(0:lmax,0:mmax)
  double complex, intent(out) :: map(0:12*nside**2-1)
  !internal
  double precision :: S(0:12*nside**2-1,2)
  double complex :: alm(2,0:lmax,0:mmax)

  alm(1,:,:) = Elm
  alm(2,:,:) = 0d0 ! Blm = 0
  call alm2map_spin(nside,lmax,mmax,spin,alm,S)
  map = S(:,1) + iu*S(:,2)

end subroutine elm2map_spin


!* note: map^*(spin) = map(-spin)
subroutine blm2map_spin(nside,lmax,mmax,spin,Blm,map)
! map = S1 + i*S2 = sum_{lm} B_{lm} sY_{lm} 
! (E_{lm}=a^+_{lm})
  implicit none
  !I/O
  integer, intent(in) :: nside, lmax, mmax, spin
  double complex, intent(in) :: Blm(0:lmax,0:mmax)
  double complex, intent(out) :: map(0:12*nside**2-1)
  !internal
  double precision :: S(0:12*nside**2-1,2)
  double complex :: alm(2,0:lmax,0:mmax)

  alm(1,:,:) = 0d0
  alm(2,:,:) = Blm
  call alm2map_spin(nside,lmax,mmax,spin,alm,S)
  map = S(:,1) + iu*S(:,2)

end subroutine blm2map_spin


subroutine pureEB_full(nside,npix,lmax,Pobs,W0)
  implicit none
  !I/O
  integer, intent(in) :: nside, npix, lmax
  double precision, intent(inout) :: Pobs(0:npix-1,2),W0(0:npix-1,2)
  !internal
  integer :: i, l
  double precision :: al, n1, n2, pd
  double precision, allocatable, dimension(:,:) :: W1,W2,P1,P0
  double complex, allocatable, dimension(:,:,:) :: alm1, alm2, wlm, tlm

  !derivatives of window functions in harmonic space
  allocate(wlm(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,0,W0,wlm)
  wlm = -wlm
  allocate(W1(0:npix-1,2),W2(0:npix-1,2),tlm(2,0:lmax,0:lmax))
  tlm = 0d0
  do l = 2, lmax
    al = dble(l)
    pd = dsqrt((al+1d0)*al)
    tlm(1,l,:) = - wlm(1,l,:)*pd
    tlm(2,l,:) = 0d0
  end do
  call alm2map_spin(nside,lmax,lmax,1,tlm,W1)
  do l = 2, lmax
    al = dble(l)
    pd = dsqrt((al+2d0)*(al**2-1d0)*al)
    tlm(1,l,:) = - wlm(1,l,:)*pd
    tlm(2,l,:) = 0d0
  end do
  call alm2map_spin(nside,lmax,lmax,2,tlm,W2)
  deallocate(tlm,wlm)

  !compute P1=conjg(W1)*(Q+/-iU), P0=conjg(W2)*(Q+/-iU)
  allocate(P1(0:npix-1,2),P0(0:npix-1,2))
  P1(:,1) = W1(:,1)*Pobs(:,1)+W1(:,2)*Pobs(:,2)
  P1(:,2) = W1(:,1)*Pobs(:,2)-W1(:,2)*Pobs(:,1)
  P0(:,1) = W2(:,1)*Pobs(:,1)+W2(:,2)*Pobs(:,2)
  P0(:,2) = W2(:,1)*Pobs(:,2)-W2(:,2)*Pobs(:,1)
  deallocate(W1,W2)

  allocate(alm1(2,0:lmax,0:lmax),alm2(2,0:lmax,0:lmax))
  alm1 = 0d0
  call map2alm_spin(nside,lmax,lmax,1,P1,alm1)
  call map2alm_spin(nside,lmax,lmax,0,P0,alm2)
  deallocate(P1,P0)

  do l = 2, lmax
    al = dble(l)
    n1 = 2d0/dsqrt((al+2d0)*(al-1d0))
    n2 = 1d0/dsqrt((al+2d0)*(al**2-1d0)*al)
    alm1(:,l,:) = n1*alm1(:,l,:) + n2*alm2(:,l,:)
  end do
  deallocate(alm2)

  allocate(P1(0:npix-1,2))
  call alm2map_spin(nside,lmax,lmax,2,alm1,P1)
  deallocate(alm1)

  ! P2 = conjg(W0)*(Q+/-iU)
  Pobs(:,1) = Pobs(:,1)*W0(:,1)
  Pobs(:,2) = Pobs(:,2)*W0(:,1)

  ! total
  Pobs = Pobs + P1

end subroutine pureEB_full


subroutine CHI_FIELD(nside,npix,lmax,Pobs,W0)
  implicit none
  !I/O
  integer, intent(in) :: nside, npix, lmax
  double precision, intent(inout) :: Pobs(0:npix-1,2),W0(0:npix-1,2)
  !internal
  integer :: i, l
  double precision :: al, n1, n2, pd
  double precision, allocatable, dimension(:,:) :: P2,chi
  double complex, allocatable, dimension(:,:,:) :: alm,tlm

  allocate(P2(0:npix-1,2))
  P2 = Pobs

  !derivatives of window functions in harmonic space
  allocate(alm(2,0:lmax,0:lmax))
  call map2alm_spin(nside,lmax,lmax,2,P2,alm)
  deallocate(P2)
  allocate(tlm(2,0:lmax,0:lmax),chi(0:npix-1,2))
  do l = 2, lmax
    al = dble(l)
    pd = dsqrt((al+2d0)*(al**2-1d0)*al)
    tlm(1,l,:) = alm(1,l,:)*cmplx(pd)
    tlm(2,l,:) = cmplx(0d0)
  end do
  call alm2map_spin(nside,lmax,lmax,0,tlm,chi)
  Pobs(:,1) = chi(:,1)!*W0(:,1)
  do l = 2, lmax
    al = dble(l)
    pd = dsqrt((al+2d0)*(al**2-1d0)*al)
    tlm(1,l,:) = alm(2,l,:)*cmplx(pd)
    tlm(2,l,:) = cmplx(0d0)
  end do
  call alm2map_spin(nside,lmax,lmax,0,tlm,chi)
  Pobs(:,1) = chi(:,1)!*W0(:,1)

  deallocate(tlm,alm)
  deallocate(chi)

end subroutine CHI_FIELD


end module anafull

