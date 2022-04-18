!////////////////////////////////////////////////////!
! Map Making with Inpainting
!////////////////////////////////////////////////////!

program main
  use readFile, only: set_params_file, read_prm, read_int, read_dbl, read_str
  use myconst,  only: pi, dlc, ac2rad, Tcmb, iu
  use myutils,  only: ones_d
  use mycls,    only: readcl_camb
  use anaflat,  only: window_generate, window_norm, gaussian_alm
  use myfftw,   only: dft, dft_pol, gaussian_map, gaussian_map_pol
  implicit none
  integer :: simn, eL(2), tL(2), mn, nu, n, npix, nn(2), mm(2), mpix, nc
  double precision :: s, mr, D(2), Wn(4), fsky
  double precision, allocatable :: Nl(:), Cl(:,:), W(:), QU(:,:)
  complex(dlc), allocatable :: T(:), EB(:,:), nois(:,:), vec(:,:,:)

  call set_params_file

  !parameters
  simn = read_int('snum')
  D    = read_dbl("pixsize")*pi/180d0 * nn
  fsky = read_dbl('fsky')
  nc   = read_int('nchan')
  call read_prm('eL',el)  
  call read_prm('tL',tL)
  call read_prm('nside',nn)
  call read_prm('mside',mm)

  mn = read_int('mask_num')
  mr = 0d0
  if (mn/=0)  mr = read_dbl('mask_size')

  npix = nn(1)*nn(2)
  mpix = mm(1)*mm(2)

  allocate(cl(4,tL(2)))
  call readcl_camb(cl,read_str('clfile'),tL,hasBB=.true.,rowCl=.false.)

  !Mask map
  allocate(W(npix))
  call window_generate(nn,D,W,mr,mn)
  call window_norm(W,Wn)

  do n = 1, simn

    !generate noise map
    allocate(nois(2,npix),Nl(tL(2)))
    do nu = 1, nc
      Nl = (read_dbl("sigmaT")*ac2rad/Tcmb)**2 * ones_d(npix)
      call gaussian_alm(nn,D,tL,nois(1,:),Nl)
      Nl = (read_dbl("sigmaP")*ac2rad/Tcmb)**2 * ones_d(npix)
      call gaussian_alm(nn,D,tL,nois(2,:),Nl)
    end do
    deallocate(Nl)

    !generate signal
    allocate(T(npix),QU(npix,2))
    call gaussian_map(nn,mm,D,tL,cl(1,:),T)
    call gaussian_map_pol(nn,mm,D,tL,cl(2,:),cl(3,:),QU)

    !add instrumental noise
    T       = (T+nois(1,:))*W
    QU(:,1) = (QU(:,1)+nois(2,:))*W
    QU(:,2) = (QU(:,2)+nois(2,:))*W

    !construct data vector
    allocate(vec(3,nc,npix),EB(2,npix))
    do nu = 1, nc 
      call dft(T,nn,D,1)
      call dft_pol(QU,nn,D,EB,1)
      vec(1,nu,:) = T
      vec(2,nu,:) = EB(1,:)
      vec(3,nu,:) = EB(2,:)
    end do
    deallocate(T,QU,EB)
    !! need to restore !!
    !call inpaint_interface(DATA,W,pmap,Read_Logical('Cinv'))
    deallocate(vec)
  end do

  deallocate(W,nois,cl)

end program main

