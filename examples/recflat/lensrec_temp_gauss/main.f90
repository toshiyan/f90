!///////////////////////////////////////////////////////////////////////!
! * Lensing reconstruction from Gaussian temperature map
!///////////////////////////////////////////////////////////////////////!

program main
  use readfile,  only: set_params_file, read_prm, read_str, read_log, read_int, read_dbl, read_val
  use myutils,   only: meanvar, savetxt, linspace, str
  use myconst,   only: dlc, pi, TT, ac2rad, Tcmb
  use mycls,     only: readcl_camb, cl2c2d, alm2bcl_flat, binned_ells, calcbcl_flat
  use anaflat,   only: elarray, gaussian_alm
  use recflat,   only: quadtt, alflat_tt
  implicit none
  integer :: eL(2), oL(2), rL(2), nn(2), npix, simn, i, j, l, n, bmax, cln=2
  double precision :: D(2)
  double precision, dimension(:), allocatable :: els, C2d, A2d, Ab, bp, nt, Il
  double precision, dimension(:,:), allocatable :: mCb, vCb, bc, Cl, OC
  double precision, dimension(:,:,:), allocatable :: Cb
  complex(dlc), dimension(:), allocatable :: nlm, T
  complex(dlc), dimension(:,:), allocatable :: est

  call set_params_file

  !* variables
  simn = read_int('simn')
  bmax = read_int('bmax')
  call read_prm('rL',rL)
  call read_prm('oL',oL)
  call read_prm('nside',nn)
  call read_prm('area',D)
  eL = [min(rL(1),oL(1)),max(rL(2),oL(2))]
  npix = nn(1)*nn(2)
  D = D*pi/180d0

  allocate(els(npix),bp(bmax+1),bc(1,bmax))
  call binned_ells(oL,bp,bc(1,:))
  els = elarray(nn,D)

  !* read cl
  allocate(Cl(7,eL(2)),OC(7,eL(2)),nt(eL(2)),C2d(npix))
  call readcl_camb(Cl,'../../dat/lensedfid_P15.dat',eL,.true.)
  call cl2c2d(els,Cl(TT,:),eL,C2d)

  !* noise
  nt = (6d0*ac2rad/Tcmb)**2*dexp((linspace(1,eL(2))*1.4d0*ac2rad)**2/(8d0*log(2d0)))
  OC(TT,:) = Cl(TT,:) + nt

  !* optimal diagonal filter
  allocate(Il(npix));  Il=0d0
  call cl2c2d(els,1d0/OC(TT,:),rL,Il)

  !* normalization
  allocate(A2d(npix),Ab(bmax))
  call alflat_tt(nn,D,Il,C2d,rL,oL,A2d)
  !do n = 1, npix
  !  if (A2d(n)/=0d0) A2d(n) = 1d0/A2d(n)
  !end do
  !call calcbcl_flat(bmax,oL,els,A2d,Ab)
  !call savetxt('ab_ana.dat',bc(1,:),Ab,ow=.true.)

  allocate(Cb(simn,cln,bmax));  Cb = 0d0;
  do i = 1, simn
    write(*,*) i
    !* generate unlensed CMB
    allocate(T(npix),nlm(npix),est(2,npix)); est=0d0
    call gaussian_alm(nn,D,eL,T,Cl(TT,:))
    !add noise
    call gaussian_alm(nn,D,rL,nlm,nt)
    T = T*Il+nlm
    call quadtt(nn,D,T,T,C2d,rL,est(1,:),est(2,:))
    call alm2bcl_flat(bmax,oL,els,D,est(1,:),Cb=Cb(i,1,:))
    deallocate(est,T,nlm)
  end do

  allocate(mCb(cln,bmax),vCb(cln,bmax))
  do j = 1, cln
    call meanvar(Cb(:,j,:),mCb(j,:),vCb(j,:))
  end do
  mCb(2,:) = Ab
  call savetxt('cl_l'//str(rL(1))//'-'//str(rL(2))//'_r'//str(simn)//'.dat',bc,mCb,vCb,ow=.true.)
  deallocate(mCb,vCb,Cb)

end program main

