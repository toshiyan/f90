!////////////////////////////////////////////////////!
! * Quick estimate of reconstruction noise
!////////////////////////////////////////////////////!

program main
  use readfile,  only: SET_PARAMS_FILE, set_prm, read_prm, read_str, read_val
  use myconst,   only: TT, TE, EE, BB, dd
  use mylensrec, only: QTT, QTE, QTB, QEE, QEB, QBB, QMV
  use myutils,   only: linspace, savetxt
  use mycmbexp,  only: instnl
  use mycls,     only: readcl_camb
  use nldd_interface, only: Al_interface
  implicit none
  character(LEN=50) :: key(1:5)
  logical :: QDO(1:7)
  integer :: eL(2), rL(2), lmin, lmax
  double precision, allocatable :: l(:),UC(:,:),LC(:,:),OC(:,:),Nl(:,:),Al(:,:,:)

  call set_params_file
  call set_prm('lmin',lmin,10)
  call set_prm('lmax',lmax,100)

  call read_prm('rL',rL)
  QDO = .false.
  call read_prm('DO_QUAD',QDO)

  !* read theoretical CMB power spectrum
  eL = [2,3000]
  allocate(UC(7,eL(2)),LC(7,eL(2)))
  call readcl_camb(UC,read_str('ucl'),eL)
  call readcl_camb(LC,read_str('lcl'),eL,.true.)

  !* experiment 1
  allocate(Nl(2,eL(2)),OC(7,eL(2)),Al(QMV,2,eL(2)),l(eL(2)))
  l = linspace(1,eL(2))
  key(1) = 'nchan1'
  key(3:5) = ['fwhm','sT','sP']
  call instnl(eL,Nl(1,:),Nl(2,:),mykey=key)
  OC(TT,:) = LC(TT,:) + Nl(1,:)
  OC(EE,:) = LC(EE,:) + Nl(2,:)
  OC(BB,:) = LC(BB,:) + Nl(2,:)
  OC(TE,:) = LC(TE,:)
  call savetxt('nl.dat',l,Nl(2,:),ow=.true.)
  write(*,*) dsqrt(sum((l(lmin:lmax)+0.5)*LC(BB,lmin:lmax)**2/OC(BB,lmin:lmax)**2))
  deallocate(Nl)

  call al_interface(rL,eL,LC,OC,Al(:,1,:),Al(:,2,:),QDO)
  !call savetxt('al.dat',l,Al(QTT,1,:),Al(QTE,1,:),Al(QTB,1,:),Al(QEE,1,:),Al(QEB,1,:),Al(QMV,1,:),ow=.true.)
  call savetxt('al.dat',l,Al(QEB,1,:),ow=.true.)
  write(*,*) dsqrt(sum((l(lmin:lmax)+0.5)*UC(dd,lmin:lmax)**2/(UC(dd,lmin:lmax)+Al(QEB,1,lmin:lmax))**2))
  deallocate(LC,UC,OC,Al,l)

end program main

