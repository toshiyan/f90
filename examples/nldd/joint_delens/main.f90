!////////////////////////////////////////////////////!
! * Quick estimate of joint delensing efficiency
!////////////////////////////////////////////////////!

program main
  use readfile,  only: SET_PARAMS_FILE, read_prm, read_str
  use myconst,   only: TT, TE, EE, BB, dd
  use mylensrec, only: QTT, QTE, QTB, QEE, QEB, QBB, QMV
  use myutils,   only: linspace, savetxt
  use mycmbexp,  only: instnl
  use mycls,     only: readcl_camb
  use nldd_interface, only: al_interface
  use nldd_lens, only: clbb_lin, res_clbb
  implicit none
  character(LEN=50) :: key(1:5)
  logical :: CV(2), QDO(1:7)
  integer :: l, el(2), rL(2), dL(2)
  double precision, dimension(:,:), allocatable :: UC,LC,OC
  double precision, dimension(:,:,:), allocatable :: Al,Nl,LCN

  call set_params_file

  call read_prm('rL',rL)
  call read_prm('dL',dL)
  call read_prm('CV',CV)
  QDO = .false.
  call read_prm('DO_QUAD',QDO(1:6))

  !* read theoretical CMB power spectrum
  eL(1) = min(rL(1),dL(1))
  eL(2) = max(rL(2),rL(2))
  allocate(UC(7,eL(2)),LC(7,eL(2)),LCN(2,7,eL(2)))
  call readcl_camb(UC,read_str('ucl'),eL)
  call readcl_camb(LC,read_str('lcl'),eL,.true.)
  LCN(1,:,:) = UC
  LCN(2,:,:) = LC

  !* experiment 1
  allocate(Nl(2,2,eL(2)))
  if(.not.CV(1)) then
    key(1) = 'nchan1'
    key(3:5) = ['fwhm_arcmin1','sigma_T1','sigma_P1']
    call instnl(eL,Nl(1,1,:),Nl(1,2,:),mykey=key)
    LCN(1,TT,:) = LC(TT,:) + Nl(1,1,:)
    LCN(1,EE,:) = LC(EE,:) + Nl(1,2,:)
    LCN(1,BB,:) = LC(BB,:) + Nl(1,2,:)
  end if

  !* experiment 2
  if(.not.CV(2)) then
    key(1) = 'nchan2'
    key(3:5) = ['fwhm_arcmin2','sigma_T2','sigma_P2']
    call instnl(eL,Nl(2,1,:),Nl(2,2,:),mykey=key)
    LCN(2,TT,:) = LC(TT,:) + Nl(2,1,:)
    LCN(2,EE,:) = LC(EE,:) + Nl(2,2,:) !using E-mode from ground exp.
    LCN(2,BB,:) = LC(BB,:) !B-mdoe noise 2 does not affect the results
  end if

  deallocate(Nl)

  allocate(OC(7,eL(2)),Al(QMV,2,dL(2)),Nl(QMV,2,dL(2)))
  OC       = LCN(1,:,:)
  OC(TE,:) = LC(TE,:)
  LC(dd,:) = UC(dd,:)

  !* phi noise calculation
  call al_interface(rL,dL,LC,OC,Al(:,1,:),Al(:,2,:),QDO)

  !* computing residual Cl^BB
  Al(QMV,1,:) = 1d0/(1d0/Al(QEB,1,:)+1d0/Al(QEE,1,:))
  call res_clbb(eL,dL,LC(EE,:),UC(dd,:),OC(BB,:),NE=LCN(2,EE,:)-LC(EE,:),Np=Al(QMV,1,:))
  call clbb_lin(eL,dL,LC(EE,:),UC(dd,:),LC(BB,:))

  !* residual ClBB
  call savetxt('rClBB.dat',linspace(1,eL(2)),Al(QMV,1,:),OC(BB,:),LC(BB,:),OC(BB,:)/LC(BB,:))
  deallocate(LC,UC,OC,Al,Nl,LCN)

end program main

