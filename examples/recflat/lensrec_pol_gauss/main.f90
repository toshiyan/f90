!///////////////////////////////////////////////////////////////////////!
! * Example code for lensing reconstruction in flat sky
!///////////////////////////////////////////////////////////////////////!

program main
  use myutils,   only: meanvar, savetxt, linspace
  use myconst,   only: dlc, pi, TT, EE, BB
  use mycls,     only: readcl_camb, alm2bcl_flat, binned_ells, calcbcl_flat, cl2c2d
  use anaflat,   only: elarray, gaussian_alm
  use recflat,   only: quadtt, quadee, quadeb, alflat_ee, alflat_eb
  implicit none
  integer :: el(2), nn(2), npix, simn, i, j, l, n, bmax, cln=4
  integer, parameter :: gEE=1, gEB=2, cEE=3, cEB=4
  double precision :: D(2)
  double precision, dimension(:), allocatable :: els, bp
  double precision, dimension(:,:), allocatable :: mCb, c2d, Il, A2d, bc, Al, Cl
  double precision, dimension(:,:,:), allocatable :: Cb
  complex(dlc), dimension(:,:), allocatable :: alm, est

  !* variables
  simn  = 10
  bmax  = 100
  nn    = [708,300]
  npix  = nn(1)*nn(2)
  D     = [59d0,25d0]*pi/180d0
  eL    = [20,1500]

  allocate(bp(bmax+1),bc(1,bmax))
  call binned_ells(eL,bp,bc(1,:))
  allocate(els(npix));  els=0d0
  els = elarray(nn,D)

  allocate(Cl(7,eL(2)),c2d(2,npix))
  call readcl_camb(Cl,'../../dat/lensedfid_P15.dat',eL,.true.)
  call cl2c2d(els,Cl(EE,:),eL,c2d(1,:))

  !* optimal diagonal filter
  allocate(Il(2,npix));  Il=0d0
  do n = 1, npix
    if (eL(1)<=els(n).and.els(n)<=eL(2)) then
      Il(1,n) = 1d0/Cl(EE,int(els(n)))
      Il(2,n) = 1d0/Cl(BB,int(els(n)))
    end if
  end do

  allocate(A2d(4,npix),Al(4,bmax));  A2d=0d0;  Al=0d0
  call alflat_ee(nn,D,Il(1,:),c2d(1,:),eL,eL,A2d(gEE,:),A2d(cEE,:))
  call alflat_eb(nn,D,Il(1,:),Il(2,:),c2d(1,:),eL,eL,A2d(gEB,:),A2d(cEB,:))
  do i = 1, 4
    call calcbcl_flat(bmax,eL,els,A2d(i,:),Al(i,:))
  end do
  call savetxt('al_sim.dat',bc,Al,ow=.true.)
  deallocate(Al)

  allocate(Cb(simn,cln,bmax));  Cb=0d0;
  do i = 1, simn
    allocate(alm(2,npix),est(4,npix));  est=0d0
    call gaussian_alm(nn,D,eL,alm(1,:),Cl(EE,:))
    call gaussian_alm(nn,D,eL,alm(2,:),Cl(BB,:))
    alm = alm*Il
    call quadee(nn,D,alm(1,:),alm(1,:),c2d(1,:),eL,est(gEE,:),est(cEE,:))
    call quadeb(nn,D,alm(1,:),alm(2,:),c2d(1,:),eL,est(gEB,:),est(cEB,:))
    do j = 1, 4
      call alm2bcl_flat(bmax,eL,els,D,est(j,:),Cb=Cb(i,j,:))
      do n = 1, bmax
        if (Cb(i,j,n)/=0d0) Cb(i,j,n) = 1d0/Cb(i,j,n)
      end do
    end do
    deallocate(est,alm)
  end do

  allocate(mCb(cln,bmax))
  do j = 1, cln
    call meanvar(Cb(1:simn,j,:),mCb(j,:))
  end do
  call savetxt("cl.dat",bc,mCb,ow=.true.)
  deallocate(mCb,Cb)

end program main

