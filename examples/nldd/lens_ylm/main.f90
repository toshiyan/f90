!////////////////////////////////////////////////////!
! Phi noise with Hamonic-Based Calculation
!////////////////////////////////////////////////////!

program main
  use readFile, only: set_params_file, read_prm, read_str, read_log
  use myconst,  only: TT, EE, BB
  use myutils,  only: savetxt, linspace, str
  use mycls,    only: readcl_camb
  use mycmbexp, only: instnl
  use mylensrec, only: QTT, QTE, QTB, QEE, QEB, QBB, QMV
  !use nldd_ylm, only: al_harmonic
  implicit none
  integer :: l, i, rL(2), el(2)
  double precision, dimension(:), allocatable :: NLT, NLP
  double precision, dimension(:,:), allocatable :: UCl, LCl
  double precision, dimension(:,:,:), allocatable :: Nl, Al

  call set_params_file

  call read_prm('eL',eL)
  call read_prm('rL',rL)

  allocate(UCl(7,rL(2)),LCl(7,rL(2)))
  call readcl_camb(UCl,read_str('ucl'),rL)
  call readcl_camb(LCl,read_str('lcl'),rL)
  if(read_log('useLCl')) UCL = LCL

  if(.not.read_log('CV')) then
    allocate(NlT(rL(2)),NlP(rL(2)))
    call instnl(rL,NlT,NlP)
    LCl(TT,:) = LCl(TT,:) + NlT
    LCl(EE,:) = LCl(EE,:) + NlP
    LCl(BB,:) = LCl(BB,:) + NlP
    deallocate(NlT,NlP)
  end if

  allocate(Al(QMV,2,eL(2)),Nl(QMV-1,2,eL(2)))
  !call al_harmonic(eL,rL,UCl,LCl,Al,Nl)
  do i = 1, QMV
    call savetxt('Al'//str(i)//'.dat',linspace(1,eL(2)),Al(i,1,:),Al(i,2,:))
  end do
  do i = 1, QMV-1
    call savetxt('Nl'//str(i)//'.dat',linspace(1,eL(2)),Nl(i,1,:),Nl(i,2,:))
  end do
  deallocate(LCl,UCl,Al,Nl)

end program main

