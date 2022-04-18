!////////////////////////////////////////////////////!
! * Quick delensing efficiency calculation
!////////////////////////////////////////////////////!

program main
  use readfile,  only: set_params_file, read_prm, read_str, read_log
  use myconst,   only: TT, TE, EE, BB, dd, pi
  use mylensrec, only: QTT, QTE, QTB, QEE, QEB, QBB, QMV
  use myutils,   only: linspace, savetxt
  use mycmbexp,  only: instnl
  use mycls,     only: readcl_camb
  use nldd_interface, only: al_interface
  use nldd_delens, only: clbb_lin, res_clbb
  implicit none
  character(LEN=50) :: key(1:5)
  logical :: QDO(1:7)
  integer :: el(2), rL(2), dL(2)
  double precision, allocatable :: l(:),UC(:,:),LC(:,:),OC(:,:),Nl(:,:),Al(:,:,:), rBB(:,:)

  call set_params_file
  call read_prm('rL',rL)
  call read_prm('dL',dL)
  QDO = .false.
  call read_prm('DO_QUAD',QDO)

  !* read theoretical CMB power spectrum
  eL(1) = min(rL(1),dL(1))
  eL(2) = max(rL(2),rL(2))
  allocate(UC(7,eL(2)),LC(7,eL(2)),l(eL(2)))
  call readcl_camb(UC,read_str('ucl'),eL)
  call readcl_camb(LC,read_str('lcl'),eL,.true.)
  l = linspace(1,eL(2))
  !UC(dd,:) = UC(dd,:)*(l*(l+1d0)/(2d0*pi)/l**4) !for P13
  LC(dd,:) = UC(dd,:) !for P15

  !* experimental noise if needed
  allocate(OC(7,eL(2)),Al(QMV,2,dL(2)),Nl(2,eL(2)),rBB(2,eL(2)))
  OC = LC
  if(.not.read_log('CV')) then
    key(1) = 'nchan1'
    key(3:5) = ['fwhm','sT','sP']
    call instnl(eL,Nl(1,:),Nl(2,:),mykey=key)
    OC(TT,:) = LC(TT,:) + Nl(1,:)
    OC(EE,:) = LC(EE,:) + Nl(2,:)
    OC(BB,:) = LC(BB,:) + Nl(2,:)
  end if

  !* rec noise calculation
  !call al_interface(rL,dL,LC,OC,Al(:,1,:),Al(:,2,:),QDO,itern=10)
  call al_interface(rL,dL,LC,OC,Al(:,1,:),Al(:,2,:),QDO)

  !* computing residual Cl^BB
  call res_clbb(eL,dL,LC(EE,:),UC(dd,:),rBB(1,:),NE=Nl(2,:),Np=Al(QTT,1,:))
  call res_clbb(eL,dL,LC(EE,:),UC(dd,:),rBB(2,:),NE=Nl(2,:),Np=Al(QEB,1,:))
  !call clbb_lin(eL,dL,LC(EE,:),UC(dd,:),LC(BB,:))

  !* residual ClBB
  call savetxt('rClBB.dat',l,rBB(1,:),rBB(2,:),LC(BB,:),ow=.true.)
  deallocate(l,LC,UC,OC,Al,Nl,rBB)

end program main

