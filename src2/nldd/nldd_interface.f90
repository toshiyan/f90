!////////////////////////////////////////////////////!
! * All and iterative estimator
!////////////////////////////////////////////////////!

module nldd_interface
  use myconst, only: pi, TT, TE, EE, BB, dd
  use mylensrec, only: QTT, QTE, QTB, QEE, QEB, QBB, QMV, QTTTE, QTTEE, QTEEE, QTBEB
  use mylapack95, only: inv_lapack
  use nldd_lens, only: AlTT, AlTE, AlTB, AlEE, AlEB, AlBB, ILTTTE, ILTTEE, ILTEEE, ILTBEB
  use nldd_delens, only: CLBB_EST
  implicit none

  private pi, TT, TE, EE, BB, dd
  private QTT, QTE, QTB, QEE, QEB, QBB, QMV, QTTTE, QTTEE, QTEEE, QTBEB
  private inv_lapack
  private AlTT, AlTE, AlTB, AlEE, AlEB, AlBB, CLBB_EST, ILTTTE, ILTTEE, ILTEEE, ILTBEB

contains


subroutine OPTCOMB(QDO,Al,Il,Nl)
  implicit none
  !I/O
  logical, intent(in) :: QDO(:)
  double precision, intent(in) :: Il(:)
  double precision, intent(inout) :: Al(:)
  double precision, intent(out), optional :: Nl(:)
  !internal
  integer :: X, Y, qmax, i, id(6)
  double precision, dimension(:,:), allocatable :: M

  id = 0

  !* set ids
  do X = 1, 6
    if(QDO(X)) id(X) = 1 + maxval(id)
  end do
  qmax = maxval(id)

  !* noise covariance
  allocate(M(qmax,qmax));  M = 0d0

  if(QDO(QTT).and.QDO(QTE)) M(id(QTT),id(QTE)) = Il(QTTTE)*Al(QTT)*Al(QTE)
  if(QDO(QTT).and.QDO(QEE)) M(id(QTT),id(QEE)) = Il(QTTEE)*Al(QTT)*Al(QEE)
  if(QDO(QTE).and.QDO(QEE)) M(id(QTE),id(QEE)) = Il(QTEEE)*Al(QTE)*Al(QEE)
  if(QDO(QTB).and.QDO(QEB)) M(id(QTB),id(QEB)) = Il(QTBEB)*Al(QTB)*Al(QEB)

  do X = 1, 6
    if (QDO(X)) M(id(X),id(X)) = Al(X)
    do Y = X + 1, 6
      if(QDO(X).and.QDO(Y)) M(id(Y),id(X)) = M(id(X),id(Y))
    end do
  end do 
  call INV_LAPACK(M)

  Al(QMV) = 1d0/sum(M)
  if (present(Nl)) Nl = sum(M,dim=2)
  deallocate(M)

end subroutine OPTCOMB


subroutine AL_INTERFACE(rL,dL,fC,OC,Alg,Alc,QDO,itern)
  implicit none
  !I/O
  logical, intent(in) :: QDO(1:7)
  integer, intent(in) :: rL(1:2), dL(1:2)
  integer, intent(in), optional :: itern
  double precision, intent(out) :: Alg(:,:), Alc(:,:)
  double precision, intent(in), dimension(:,:) :: fC, OC
  !internal
  integer :: i, n, l
  double precision :: conv, ratio
  double precision, dimension(:), allocatable :: AlgEB, rCBB
  double precision, dimension(:,:), allocatable :: Ilg, Ilc

  conv = 1d-6

  !//// interface ////!
  ! TT
  if (QDO(QTT))  call ALTT(rL,dL,Alg(QTT,:),Alc(QTT,:),fC(TT,:),OC(TT,:))
  ! TE
  if (QDO(QTE))  call ALTE(rL,dL,Alg(QTE,:),Alc(QTE,:),fC(TE,:),OC(TT,:),OC(EE,:))
  ! EE
  if (QDO(QEE))  call ALEE(rL,dL,Alg(QEE,:),Alc(QEE,:),fC(EE,:),OC(EE,:))
  ! TB
  if (QDO(QTB))  call ALTB(rL,dL,Alg(QTB,:),Alc(QTB,:),fC(TE,:),OC(TT,:),OC(BB,:))
  ! BB
  if (QDO(QBB))  call ALBB(rL,dL,Alg(QBB,:),Alc(QBB,:),fC(BB,:),OC(BB,:))
  ! EB
  if (QDO(QEB)) then
    if (present(itern).and.itern>0) then
      allocate(AlgEB(dL(2)),rCBB(rL(2)));  AlgEB = 0d0;  rCBB=0d0
      ratio = 1d0
      rCBB = OC(BB,1:rL(2))
      do n = 1, itern !loop for iteration 
        !* lensing reconstruction with B-mode
        call ALEB(rL,dL,AlgEB,Alc(QEB,:),fC(EE,:),OC(EE,:),rCBB)
        !* convergence check
        if (n>=2) then
          ratio = (sum(Alg(QEB,:))/sum(AlgEB)-1d0)/dble(dL(2))
          write(*,*) n, ratio
        end if
        Alg(QEB,1:dL(2)) = AlgEB
        if(abs(ratio) < conv) exit
        !* delensing with EB-estimator
        call CLBB_EST(rL,rL,fC(EE,:),fC(dd,:),OC(EE,:)-fC(EE,:),AlgEB(1:rL(2)),rCBB)
        rCBB = OC(BB,:) - rCBB !delensed B-mode
        if(n==itern) stop 'not converged'
      end do
      deallocate(AlgEB,rCBB)
    else
      call ALEB(rL,dL,Alg(QEB,:),Alc(QEB,:),fC(EE,:),OC(EE,:),OC(BB,:))
    end if
  end if

  !MV
  if(QDO(QMV)) then
    allocate(Ilg(4,dL(2)),Ilc(4,dL(2)))
    ! TT x TE
    if(QDO(QTT).and.QDO(QTE)) call ILTTTE(rL,dL,Ilg(QTTTE,:),Ilc(QTTTE,:),fC(TT,:),fC(TE,:),OC(TT,:),OC(EE,:),OC(TE,:))
    ! TT x EE
    if(QDO(QTT).and.QDO(QEE)) call ILTTEE(rL,dL,Ilg(QTTEE,:),Ilc(QTTEE,:),fC(TT,:),fC(EE,:),OC(TT,:),OC(EE,:),OC(TE,:))
    ! TE x EE
    if(QDO(QTE).and.QDO(QEE)) call ILTEEE(rL,dL,Ilg(QTEEE,:),Ilc(QTEEE,:),fC(EE,:),fC(TE,:),OC(TT,:),OC(EE,:),OC(TE,:))
    ! TB x EB
    if(QDO(QTB).and.QDO(QEB)) call ILTBEB(rL,dL,Ilg(QTBEB,:),Ilc(QTBEB,:),fC(EE,:),fC(BB,:),fC(TE,:),OC(TT,:),OC(EE,:),OC(BB,:),OC(TE,:))
    do l = dL(1), dL(2)
      call OPTCOMB(QDO,Alg(:,l),Ilg(:,l))
      call OPTCOMB(QDO,Alc(:,l),Ilc(:,l))
    end do
    deallocate(Ilg,Ilc)
  end if

end subroutine AL_INTERFACE


end module nldd_interface

