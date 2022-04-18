!////////////////////////////////////////////////////!
! * Quick delensing efficiency calculation
!////////////////////////////////////////////////////!

program main
  use readfile,  only: set_params_file, read_prm, read_dbl, read_str
  use myconst,   only: TT, TE, EE, BB, dd, pi, ac2rad, Tcmb
  use mylensrec, only: QTT, QTE, QTB, QEE, QEB, QBB, QMV
  use myutils,   only: linspace, savetxt
  use mycls,     only: readcl_camb
  use nldd_lens, only: clbb_lin, res_clbb, aleb
  implicit none
  character(LEN=50) :: key(1:5)
  logical :: QDO(1:7)
  integer :: el(2), rL(2), dL(2), nL(2), i
  double precision, allocatable :: l(:),UC(:,:),LC(:,:),OC(:,:),NT(:),Al(:,:),DB(:),fC(:)

  call set_params_file
  QDO = .false.
  QDO(5) = .true.
  
  !* read theoretical CMB power spectrum
  eL = [2,3000]
  dL = [2,3000]
  rL = [200,3000]
  allocate(UC(7,eL(2)),LC(7,eL(2)),l(eL(2)))
  call readcl_camb(UC,read_str('ucl'),eL)
  call readcl_camb(LC,read_str('lcl'),eL,.true.)
  l = linspace(1,eL(2))
  UC(dd,:) = UC(dd,:)*(l*(l+1d0)/(2d0*pi)/l**4) !for P13
  !LC(dd,:) = UC(dd,:) !for P15

  !* experimental noise if needed
  allocate(OC(7,eL(2)),Al(2,dL(2)),NT(eL(2)),DB(eL(2)),fC(eL(2))); DB=0d0
  OC = LC
  NT = (read_dbl('sT')*ac2rad/Tcmb)**2*dexp((linspace(1,eL(2))*read_dbl('fwhm')*ac2rad)**2/(8d0*log(2d0)))
  OC(EE,:) = LC(EE,:) + 2d0*NT
  OC(BB,:) = LC(BB,:) + 2d0*NT

  call aleb(rL,dL,Al(1,:),Al(2,:),LC(EE,:),OC(EE,:),OC(BB,:))
  call res_clbb(eL,dL,LC(EE,:),UC(dd,:),UC(BB,:),NE=2d0*NT,Np=Al(1,:))
  call savetxt('res_full.dat',l,UC(BB,:)/LC(BB,:),ow=.true.)
  do i = 1, 7
    write(*,*) i
    !* rec noise calculation
    nL = [200+(i-1)*400,200+i*400]
    fC = LC(EE,:)
    fC(nL(1)-50:nL(2)+50) = 0d0
    call aleb(rL,dL,Al(1,:),Al(2,:),fC,OC(EE,:),OC(BB,:))
    !* computing residual Cl^BB
    call res_clbb(eL,dL,LC(EE,:),UC(dd,:),UC(BB,:),NE=2d0*NT,Np=Al(1,:))
    DB(nL(1):nL(2)) = UC(BB,nL(1):nL(2))/LC(BB,nL(1):nL(2))
  end do

  !* residual ClBB
  call savetxt('res_lsplit.dat',l,DB,ow=.true.)
  !deallocate(l,LC,UC,OC,Al,DB,fC,NT)

end program main

