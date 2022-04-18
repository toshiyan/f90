! * Simulating Remapping

program main
  use readfile, only: set_params_file, read_prm, read_str, read_log, read_log, read_dbl, read_int, read_val
  use myconst,  only: TT, TE, EE, BB, dd, ac2rad, tcmb, pi, dlc, iu, twopi
  use myutils,  only: savetxt, linspace, meanvar
  use mycls,    only: readcl_camb, calccl_flat, calcbcl_flat, binned_ells, cl2c2d
  use mycmbexp, only: calcNl
  use anaflat,  only: window_generate, window_norm, gaussian_alm, elarrays
  use myfftw,   only: dft, dft_pol, pureEB
  use remapping, only: remap
  implicit none
  integer :: nn(1:2), mm(1:2), eL(1:2), npix, i, j, k, l, n, bmax, mapn
  double precision :: D(1:2), Wn(1:6)
  double precision, dimension(:), allocatable :: W, els, elx, ely, bc, bp, NP, bl, b2d
  double precision, dimension(:,:), allocatable :: UC, LC, QU, Cl, mCb, vCb, map
  double precision, dimension(:,:,:), allocatable :: Cb
  complex(dlc), dimension(:), allocatable :: T, phi
  complex(dlc), dimension(:,:), allocatable :: alm, nlm, def

  call set_params_file

  !map
  call read_prm('nside',nn)
  call read_prm('eL',eL)
  npix = nn(1)*nn(2)
  mm   = nn!*2
  mapn = read_int('mapn')

  !cl
  bmax = read_int('bmax')

  !* read theoretical CMB power spectrum
  allocate(UC(7,eL(2)),LC(7,eL(2)))
  call readcl_camb(UC,read_str('ucl'),eL,rowCl=.true.)
  call readcl_camb(LC,read_str('lcl'),eL,HasBB=.true.,rowCl=.true.)

  !* map characterization
  D(1) = read_dbl('psize')*pi/180d0
  D(2) = read_dbl('psize')*pi/180d0

  !* multipole
  allocate(els(npix),elx(npix),ely(npix))
  call elarrays(nn,D,elx,ely,els)

  !* noise model
  allocate(bl(eL(2)),NP(eL(2)),b2d(npix));  bl=0d0;  NP=0d0;  b2d=0d0
  !call calcnl(NP,eL,[read_dbl('sP')*ac2rad/Tcmb],[read_dbl('fwhm')*ac2rad])
  NP  = (read_dbl('sP')*ac2rad/Tcmb)**2
  bl  = dexp(-(linspace(1,eL(2))*read_dbl('fwhm')*ac2rad)**2/(16d0*log(2d0)))
  call cl2c2d(els,bl,eL,b2d)
  call savetxt('noise.dat',linspace(1,eL(2)),NP,bl**2*LC(BB,1:eL(2)))

  !* window function
  allocate(W(npix))
  call window_generate(nn,D,W,ap=read_dbl('afac'),cut=read_dbl('cutf'))
  call window_norm(W,Wn)
  write(*,*) Wn(2)

  allocate(Cb(mapn,6,bmax));  Cb=0d0

  do i = 1, mapn

    allocate(T(npix));  T=0d0
    !* generate temperature
    !write(*,*) 'generate T'
    !call gaussian_alm(nn,D,el,T,UC(TT,:))
    !call dft(T,nn,D,-1)

    !* generate QU maps
    allocate(QU(npix,2),alm(2,npix));  alm=0d0
    write(*,*) 'generate QU'
    call gaussian_alm(nn,D,el,alm(1,:),UC(EE,:))
    call dft_pol(QU,nn,D,alm,-1)
    deallocate(alm)

    !* generate deflection angle
    allocate(phi(npix),def(npix,2))
    write(*,*) 'generate deflection'
    call gaussian_alm(nn,D,el,phi,UC(dd,:))
    def(:,1) = -iu*elx*phi
    def(:,2) = -iu*ely*phi
    call dft(def(:,1),nn,D,-1)
    call dft(def(:,2),nn,D,-1)
    deallocate(phi)

    !* remapping by true phi
    if (read_log('lens')) then
      write(*,*) 'generate lensed map'
      write(*,*) dsqrt(dble(sum(def(:,1)**2+def(:,2)**2)/dble(npix))), D(1)/nn(1)
      allocate(map(npix,2))
      map = dble(def)
      !!$OMP PARALLEL DO
      do k=1, 2
        call remap(nn,D,map,QU(:,k),2)
      end do
      !!$OMP END PARALLEL DO
      deallocate(map)
    end if

    !* beam convolution + noise
    allocate(Cl(6,npix))
    if (read_log('beam')) then
      write(*,*) 'beam + noise'
      allocate(alm(2,npix),nlm(2,npix));  alm=0d0;  nlm=0d0
      call dft_pol(QU,nn,D,alm,1)
      call calccl_flat(D,alm(2,:),alm(2,:),Cl(TT,:),els)
      call gaussian_alm(nn,D,el,nlm(1,:),NP)
      call gaussian_alm(nn,D,el,nlm(2,:),NP)
      alm(1,:) = alm(1,:)*b2d + nlm(1,:)
      alm(2,:) = alm(2,:)*b2d + nlm(2,:)
      call calccl_flat(D,alm(2,:),alm(2,:),Cl(5,:),els)
      call dft_pol(QU,nn,D,alm,-1)
      deallocate(alm,nlm)
    end if

    !* inverse remapping by measured phi
    if (read_log('remap')) then
      write(*,*) 'inverse remapping'
      allocate(map(npix,2))
      map = -dble(def)
      call remap(nn,D,map,QU(:,1),2)
      call remap(nn,D,map,QU(:,2),2)
      deallocate(map)
    end if

    !* apodization
    write(*,*) 'apodization'
    QU(:,1) = QU(:,1)*W
    QU(:,2) = QU(:,2)*W

    !* T,Q,U -> T,E,B
    allocate(alm(2,npix))
    write(*,*) 'alms'
    !call dft(T,nn,D,1)
    !call dft(def(:,1),nn,D,1)
    !call dft(def(:,2),nn,D,1)
    !call pureEB(QU,nn,D,alm,W)
    call dft_pol(QU,nn,D,alm,1)
    deallocate(QU)

    !* calculate Cls at each pixel
    !call calccl_flat(D,T,T,Cl(TT,:),els)
    call calccl_flat(D,alm(1,:),alm(1,:),Cl(EE,:),els)
    call calccl_flat(D,alm(2,:),alm(2,:),Cl(BB,:),els)
    !call calccl_flat(D,def(:,1),def(:,1),Cl(5,:),els)
    !call calccl_flat(D,def(:,2),def(:,2),Cl(6,:),els)
    !Cl(5,:) = Cl(5,:) + Cl(6,:)
    Cl = Cl/Wn(2)
    deallocate(T,alm,def)

    !* calculate binned Cls
    do j = 1, 6
      call calcbcl_flat(bmax,el,els,Cl(j,:),Cb(i,j,:))
    end do
    deallocate(Cl)

  end do

  allocate(mCb(6,bmax),vCb(6,bmax),bc(bmax),bp(bmax+1))
  call binned_ells(eL,bp,bc)
  do j = 1, 6
    call meanvar(Cb(:,j,:),mCb(j,:),vCb(j,:)) 
  end do
  call savetxt('test.dat',bc,mCb(TT,:),mCb(EE,:),mCb(BB,:),mCb(5,:),vCb(TT,:),vCb(EE,:),vCb(BB,:),vCb(5,:))
  deallocate(W,els,NP,LC,UC,Cb,bc,bp,mCb,vCb,b2d)


end program main

