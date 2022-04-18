!///////////////////////////////////////////////////////////////////////!
! * Lensing reconstruction with TT estimator
!///////////////////////////////////////////////////////////////////////!

program main
  use readfile,  only: set_params_file, read_prm, read_str, read_log, read_int, read_dbl, read_val
  use myutils,   only: meanvar, savetxt, linspace, str
  use myconst,   only: dlc, pi, TT, EE, BB, dd, ac2rad, Tcmb, iu, twopi
  use mycls,     only: readcl_camb, cl2c2d, alm2bcl_flat, binned_ells, calcbcl_flat
  use anaflat,   only: elarray, gaussian_alm
  use myfftw,    only: dft
  use recflat,   only: quadtt, alflat_tt
  use remapping, only: remap_lin
  implicit none
  integer :: eL(2), oL(2), rL(2), nn(2), npix, simn, a, b, i, j, l, n, bmax
  double precision :: D(2)
  double precision, dimension(:), allocatable :: lfac, lx, ly, Il, els, A2d, bp, nt, C2d
  double precision, dimension(:,:), allocatable :: mCb, vCb, bc, UC, Cl, OC
  double precision, dimension(:,:,:), allocatable :: Cb
  complex(dlc), dimension(:), allocatable :: nlm, T, klm
  complex(dlc), dimension(:,:), allocatable :: dT, def, est

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

  !* multipole arrays
  allocate(els(npix),lfac(npix),lx(npix),ly(npix),bp(bmax+1),bc(1,bmax)); els=0d0; lfac=0d0
  call binned_ells(oL,bp,bc(1,:))
  els = elarray(nn,D)
  n = 0
  do a = 1, nn(1)
    do b = 1, nn(2)
      n = n + 1
      lx(n) = twopi*dble(a-1-nn(1)*0.5d0)/D(1)
      ly(n) = twopi*dble(b-1-nn(2)*0.5d0)/D(2)
      if (lx(n)==0d0.and.ly(n)==0d0) cycle
      lfac(n) = 2d0/(lx(n)**2+ly(n)**2)
    end do
  end do

  allocate(UC(7,eL(2)),Cl(7,eL(2)),OC(7,eL(2)),nt(eL(2)),C2d(npix))
  call readcl_camb(Cl,'../../dat/lensedfid_P15.dat',eL,.true.)
  call cl2c2d(els,Cl(TT,:),eL,C2d)
  call readcl_camb(UC,'../../dat/fid_P15.dat',eL)
  UC(dd,:) = UC(dd,:)*linspace(1,eL(2))**2/(4d0*twopi) !to kappa

  !noise
  nt = (6d0*ac2rad/Tcmb)**2*dexp((linspace(1,eL(2))*1.4d0*ac2rad)**2/(8d0*log(2d0)))
  OC(TT,:) = Cl(TT,:) + nt

  !* optimal diagonal filter
  allocate(Il(npix));  Il=0d0
  call cl2c2d(els,1d0/OC(TT,:),rL,Il)

  !* normalization
  allocate(A2d(npix))
  call alflat_tt(nn,D,Il,C2d,rL,oL,A2d)
  A2d = A2d*els**4/4d0

  !* reconstructed kappa spectrum
  allocate(Cb(simn,2,bmax));  Cb = 0d0;
  do i = 1, simn
    write(*,*) i
    !* generate unlensed CMB, klm and noise
    allocate(T(npix),nlm(npix),klm(npix))
    call gaussian_alm(nn,D,eL,T,Cl(TT,:))
    call gaussian_alm(nn,D,eL,klm,UC(dd,:))
    call gaussian_alm(nn,D,rL,nlm,nt)

    !* linear phi remapping
    call remap_lin(nn,D,lx,ly,T,klm*lfac)

    !add noise and reconstruction
    allocate(est(2,npix));  est = 0d0
    T = T*Il+nlm
    call quadtt(nn,D,T,T,C2d,rL,est(1,:),est(2,:))
    call alm2bcl_flat(bmax,oL,els,D,A2d*est(1,:)*lfac,klm,Cb=Cb(i,1,:))
    call alm2bcl_flat(bmax,oL,els,D,klm,Cb=Cb(i,2,:))
    deallocate(est,T,nlm,klm)
  end do

  allocate(mCb(2,bmax),vCb(2,bmax))
  do j = 1, 2
    call meanvar(Cb(:,j,:),mCb(j,:),vCb(j,:))
  end do
  call savetxt('cl_l'//str(rL(1))//'-'//str(rL(2))//'_r'//str(simn)//'.dat',bc,mCb,vCb,ow=.true.)
  deallocate(mCb,vCb,Cb)

end program main

