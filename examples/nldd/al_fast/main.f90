!////////////////////////////////////////////////////!
! * Rotation noise spectrum
!////////////////////////////////////////////////////!

program main
  use readfile,  only: SET_PARAMS_FILE, read_prm, read_str, read_dbl
  use myconst,   only: TT, TE, EE, BB, dd, Tcmb, ac2rad
  use myutils,   only: linspace, savetxt
  use mycls,     only: readcl_camb
  use nldd_rot,  only: AaaEB
  use nldd_tau,  only: AlTT_tau
  implicit none
  character(LEN=50) :: key(1:5)
  integer :: l, el(2), rL(2), oL(2)
  double precision, dimension(:,:), allocatable :: LC,OC
  double precision, dimension(:), allocatable :: Aaa,Nl,Att

  call set_params_file
  call read_prm('rL',rL)
  call read_prm('oL',oL)

  !* read theoretical CMB power spectrum
  eL(1) = min(rL(1),oL(1))
  eL(2) = max(rL(2),oL(2))
  allocate(LC(7,eL(2)))
  !call readcl_camb(UC,read_str('ucl'),eL)
  call readcl_camb(LC,read_str('lcl'),eL,.true.)
  LC = LC/2.72e6**2
  !call savetxt('test.dat',linspace(eL(1),eL(2),eL(2)-eL(1)+1),LC(TT,eL(1):eL(2)),LC(EE,eL(1):eL(2)),LC(BB,eL(1):eL(2)),LC(TE,eL(1):eL(2)),ow=.true.)

  allocate(Nl(eL(2)),OC(7,eL(2)),Aaa(oL(2)),Att(oL(2)))
  Nl = (read_dbl('sT')*ac2rad/Tcmb)**2*dexp(linspace(1,eL(2))*(linspace(1,eL(2))+1d0)*(read_dbl('fwhm')*ac2rad)**2/(8d0*dlog(2d0)))
  OC(TT,:) = LC(TT,:) + Nl
  OC(EE,:) = LC(EE,:) + 2d0*Nl
  OC(BB,:) = LC(BB,:)*(1d0-read_dbl('alpha')) + 2d0*Nl
  call AaaEB(rL,oL,Aaa,LC(EE,:),OC(EE,:),OC(BB,:))
  call AlTT_tau(rL,oL,Att,LC(TT,:),OC(TT,:))
  call savetxt('Al.dat',linspace(oL(1),oL(2),oL(2)-oL(1)+1),Aaa(oL(1):oL(2)),Att(oL(1):oL(2)),ow=.true.)
  deallocate(LC,OC,Aaa,Att,Nl)

end program main

