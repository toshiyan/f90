!///////////////////////////////////////////////////////////////////////!
! * Reconstructing cosmic birefringence fluctuations from Q/U maps
!///////////////////////////////////////////////////////////////////////!

program main
  use readfile,  only: set_params_file, read_prm, read_str, read_log, read_int, read_dbl, read_val
  use myutils,   only: str, meanvar, savetxt, loadtxt, linspace
  use myconst,   only: dlc, pi, iu, twopi, Tcmb
  use mycls,     only: cl2c2d, alm2bcl_flat, binned_ells, calcbcl_flat, readcl_camb
  use anaflat,   only: elarray, window_generate, gaussian_alm, window_norm, elarrays
  use myfftw,    only: dft, dft_pol
  use rotrecflat, only: quadeb_cb, alflat_eb_cb
  use remapping
  use nldd_flat, only: alxy_flat
  implicit none
  character(1) :: rtype
  integer :: el(2), rL(2), oL(2), sL(2)
  integer :: i, j, n, nn(2), npix, simn, simn0, simn1, bmax, cln=5, lxcut
  double precision :: D(2), Wn(6)
  double precision, dimension(:), allocatable :: N2d, N0, N1, N0b, N1b, els, W, linv, bp, elx, Nb, alg, caa, Alc
  double precision, dimension(:,:), allocatable :: Fl, vCb, mCb, C1d, C2d, bc, FF, Il
  double precision, dimension(:,:,:), allocatable :: Cb
  complex(dlc), dimension(:), allocatable :: rot
  complex(dlc), dimension(:,:), allocatable :: alm, est

  call set_params_file

  !* variables
  simn  = read_int('simn')
  simn0 = read_int('simn0')
  simn1 = read_int('simn1')
  bmax  = read_int('bmax')

  !* multipoles
  call read_prm('rL',rL)
  call read_prm('oL',oL)
  call read_prm('sL',sL)
  eL(1) = min(rL(1),oL(1),sL(1))
  eL(2) = max(rL(2),oL(2),sL(2))
  lxcut = 0
  allocate(bp(bmax+1),bc(1,bmax))
  call binned_ells(oL,bp,bc(1,:))
  rtype = read_str('rtype')

  !nn = [236,100]
  !D  = [59d0,25d0]*pi/180d0
  call read_prm('nside',nn)
  call read_prm('pix',D)
  D    = D*pi/180d0
  npix = nn(1)*nn(2)
  write(*,*) 2d0*pi/D(1), D(1)*D(2)

  !* multipole factors
  allocate(els(npix),linv(npix),elx(npix),FF(2,npix));  els=0d0;  linv=0d0;  elx=0d0;  FF=0d0
  call elarrays(nn,D,elx=elx,els=els,eli=linv)

  !* read theoretical CMB power spectrum for filtering
  allocate(C1d(4,eL(2)),C2d(3,npix))
  call readcl_camb(C1d,read_str('cls'),eL,.true.,.true.)
  C1d(1,:) = C1d(2,:)
  C1d(2,:) = C1d(3,:)
  call cl2c2d(els,C1d(1,:),eL,C2d(1,:))
  call cl2c2d(els,C1d(2,:),eL,C2d(2,:))
  C1d(3,:) = 2d-4*pi/linspace(1,eL(2))**2

  !* make optimal diagonal filter including ell-cut for E/B-modes
  write(*,*) 'make filtering'
  allocate(Fl(2,npix),Il(2,eL(2)));  Fl=0d0; Il=0d0
  do n = 1, npix
    if(C2d(1,n)/=0d0) Fl(1,n) = 1d0/C2d(1,n)
    if(C2d(2,n)/=0d0) Fl(2,n) = 1d0/C2d(2,n)
  end do
  !do n = 1, eL(2)
  !  if(C1d(1,n)/=0d0) Il(1,n) = 1d0/C1d(1,n)
  !  if(C1d(2,n)/=0d0) Il(2,n) = 1d0/C1d(2,n)
  !end do
  deallocate(FF)

  write(*,*) 'normalization'
  allocate(N2d(npix),Nb(bmax));  N2d=0d0
  call alflat_eb_cb(nn,D,Fl(1,:),Fl(2,:),C2d(1,:),C2d(2,:),rL,N2d,oL)
  !call alflat_eb(nn,D,Fl(1,:),Fl(2,:),C2d(1,:),rL,eL,N2d)
  call calcbcl_flat(bmax,oL,els,N2d,Nb)
  call savetxt('n2d.dat',bc(1,:),Nb,ow=.true.)

  !write(*,*) 'analytic'
  !allocate(Alg(eL(2)),Alc(eL(2))); Alg=0d0
  !call ALEB_FLAT(rL,eL,Alg,Alc,C1d(1,:),Il(1,:),Il(2,:),weight='rotation')
  !call savetxt('al.dat',linspace(1,eL(2)),Alg)
  !deallocate(Alg,Alc)

  write(*,*) 'N0'
  allocate(N0(npix),N0b(bmax)); N0 = 0d0
  do i = 1, simn0
    allocate(rot(npix),est(2,npix),alm(4,npix))
    do j = 1, 2
      call gaussian_alm(nn,D,sL,rot,C1d(3,:)) !generate rotation
      call dft(rot,nn,D,-1)
      !* generate alm
      call gaussian_alm(nn,D,sL,alm(2*j-1,:),C1d(1,:))
      call gaussian_alm(nn,D,sL,alm(2*j,:),C1d(2,:))
      call rotation_alm(nn,D,alm(2*j-1:2*j,:),dble(rot),rtype) !rotate
    end do
    !* reconstruct
    call quadeb_cb(nn,D,alm(1,:)*Fl(1,:),alm(4,:)*Fl(2,:),C2d(1,:),C2d(2,:),rL,est(1,:))
    call quadeb_cb(nn,D,alm(3,:)*Fl(1,:),alm(2,:)*Fl(2,:),C2d(1,:),C2d(2,:),rL,est(2,:))
    N0 = N0 + abs(est(1,:)+est(2,:))**2*N2d**2/dble(simn0)*0.5d0/D(1)/D(2)
    deallocate(alm,rot,est)
  end do

  call calcbcl_flat(bmax,oL,els,N0,N0b)
  call savetxt('N0.dat',bc(1,:),N0b,ow=.true.)

  write(*,*) 'N0+N1'
  allocate(N1(npix),N1b(bmax)); N1 = 0d0
  do i = 1, simn1
    allocate(rot(npix),est(2,npix))
    call gaussian_alm(nn,D,sL,rot,C1d(3,:)) !generate rotation
    call dft(rot,nn,D,-1)
    allocate(alm(4,npix))
    do j = 1, 2
      !* generate alm
      call gaussian_alm(nn,D,sL,alm(2*j-1,:),C1d(1,:))
      call gaussian_alm(nn,D,sL,alm(2*j,:),C1d(2,:))
      call rotation_alm(nn,D,alm(2*j-1:2*j,:),dble(rot),rtype)
    end do
    !* reconstruct
    call quadeb_cb(nn,D,alm(1,:)*Fl(1,:),alm(4,:)*Fl(2,:),C2d(1,:),C2d(2,:),rL,est(1,:))
    call quadeb_cb(nn,D,alm(3,:)*Fl(1,:),alm(2,:)*Fl(2,:),C2d(1,:),C2d(2,:),rL,est(2,:))
    N1 = N1 + abs(est(1,:)+est(2,:))**2*N2d**2/(2d0*D(1)*D(2)*dble(simn1))
    deallocate(alm,rot,est)
  end do

  call calcbcl_flat(bmax,oL,els,N1,N1b)
  call savetxt('N1.dat',bc(1,:),N1b,ow=.true.)

  write(*,*) 'sim'
  allocate(Cb(simn,cln,bmax))
  do i = 1, simn
    !if (mod(i,100)==0) write(*,*) i
    write(*,*) i
    allocate(rot(npix),est(2,npix))
    call gaussian_alm(nn,D,sL,rot,C1d(3,:))
    call dft(rot,nn,D,-1)
    do j = 1, 1
      allocate(alm(2,npix))
      !* generate alm
      call gaussian_alm(nn,D,sL,alm(1,:),C1d(1,:))
      call gaussian_alm(nn,D,sL,alm(2,:),C1d(2,:))
      !* rotate
      call rotation_alm(nn,D,alm,dble(rot),rtype)
      !* reconstruct
      call quadeb_cb(nn,D,alm(1,:)*Fl(1,:),alm(2,:)*Fl(2,:),C2d(1,:),C2d(2,:),rL,est(j,:))
      est(j,:) = est(j,:)*N2d
      deallocate(alm)
    end do
    call dft(rot,nn,D,1) ! alm of rotation
    call alm2bcl_flat(bmax,oL,els,D,est(1,:),Cb=Cb(i,3,:))
    call alm2bcl_flat(bmax,oL,els,D,est(1,:),rot,Cb=Cb(i,4,:))
    call alm2bcl_flat(bmax,oL,els,D,rot,Cb=Cb(i,5,:))
    Cb(i,1,:) = Cb(i,3,:) - N0b
    Cb(i,2,:) = Cb(i,3,:) - N1b
    deallocate(est,rot)
  end do

  write(*,*) 'averaged cls'
  allocate(mCb(cln,bmax),vCb(cln,bmax))
  do i = 1, cln
    call meanvar(Cb(:,i,:),mCb(i,:),vCb(i,:))
  end do
  call savetxt('mCb.dat',bc,mCb,vCb,ow=.true.)

  deallocate(mCb,vCb)

end program main

