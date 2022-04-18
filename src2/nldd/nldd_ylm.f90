!////////////////////////////////////////////////////!
! * Nldd with Wigner 3j formalism
!////////////////////////////////////////////////////!

module nldd_ylm
  use myconst, only: pi
  use myutils, only: linspace
  use funcs
  use nldd_interface, only: optcomb

  private pi
  private linspace
  private optcomb

contains


subroutine al_harmonic(eL,tL,fC,OC,Al,Nl)
  implicit none
  !I/O
  integer, intent(in) :: eL(2), tL(2)
  double precision, intent(in) :: fC(:,:), OC(:,:)
  double precision, intent(out) :: Al(:,:,:), Nl(:,:,:)
  !internal
  integer :: jmax, i, j, ell, L1, L2max, n, L2, L, nmax
  integer, dimension(:), allocatable :: llens
  double precision :: Lcut, W3j(3,-1:3,-1:1)
  double precision, dimension(:,:), allocatable :: Aj, Ij

  !* multipole range
  jmax = eL(2)-eL(1) + 1
  allocate(llens(1:jmax))
  llens = linspace(eL,jmax)

  !//// calculate normalization and noise covariance ////!
  W3j = 0d0

  !* loop for L
  do j = 1, jmax
    allocate(Aj(Qn,3),Ij(4,3))
    Aj = 0d0
    Ij = 0d0
    l = llens(j)
    write(*,*) l

    Lcut = tL(1) + 1
    if (L*0.5d0 >= tL(1)+1) Lcut = L*0.5d0
    if (dint(Lcut).ne.Lcut) Lcut = Lcut+0.5d0
    if (Lcut.le.2d0) Lcut = 2d0

    do L1 = tL(1), tL(2)
      L2max = l + L1
      !* L2 = L2max
      L2 = L2max
      call ini_factors(l,L1,W3j(1,:,:))
      if (raneq(L2,tL)) call gweight(O,Q,l,L1,L2,W3j(1,:,:),fC,OC,Aj,Ij)
      !* L2 = L2max - 1
      L2 = L2max - 1
      call sec_factors(l,L1,W3j(1,:,:),W3j(2,:,:))
      if (raneq(L2,tL)) call gweight(O,Q,l,L1,L2,W3j(2,:,:),fC,OC,Aj,Ij)
      !* sum over other cases, down to max (L2 = L-L1, tlmin)
      if (L1 < int(Lcut))  nmax = min(L2max-L+L1-1, L2max-1-tL(1))
      if (L1 >= int(Lcut)) nmax = L-1
      do n = 1, nmax
        L2 = L2max-1-n
        call recursions(L,L1,L2,0,0,W3j)
        if (raneq(L2,tL)) call gweight(O,Q,L,L1,L2,W3j(3,:,:),fC,OC,Aj,Ij)
        call shift(W3j)
      end do 
    end do 
    Al(:,:,l) = dble(2*l+1)/Aj
    Ij(1:4,1) = Ij(1:4,1)/dble(2*l+1)
    Ij(1:4,3) = Ij(1:4,2)/dble(2*l+1)
    if(QMV>2) call OPTCOMB(QDO,Al(:,1,l),Ij(:,1),Nl(:,1,l))
    if(QMV>2) call OPTCOMB(QDO,Al(:,3,l),Ij(:,3),Nl(:,3,l))
    deallocate(Aj,Ij)
  end do
  deallocate(llens)

end subroutine al_harmonic


subroutine gweight(L0,L1,L2,W,fC,oC,Al,Il)
  implicit none
  !I/O
  integer, intent(in) :: L0, L1, L2
  double precision, intent(in) :: fC(:,:), oC(:,:), W(-1:3,-1:1)
  double precision, intent(inout) :: Al(:,:), Il(:,:)
  !internal
  double precision :: fac, D
  double precision :: Sg0_12, Sg2_12, Sc0_12, Sc2_12, Sg0_21, Sg2_21, Sc0_21, Sc2_21
  double precision :: fgTT, fgTE12,fgTE21, fgTB12,fgTB21, fgEE, fgEB12,fgEB21, fgBB
  double precision :: fcTT, fcTE12,fcTE21, fcTB12,fcTB21, fcEE, fcEB12,fcEB21, fcBB
  double precision :: ggTE12,ggTE21, gcTE12,gcTE21, ATE12, ATE21, BTE

  fac = dsqrt(dble(2*L0+1)*dble(2*L1+1)*dble(2*L2+1)/(16d0*pi))

  !//// Relations ////!
  ! grad: s=0
  !  * Wigner(L2,L0,L1,0,0,0) = W_00 
  !  * Wigner(L1,L0,L2,0,0,0) = W_00 
  ! curl: s=0 
  !  * Wigner(L2,L0,L1,0,-1,1) = Wigner(L1,L0,L2,-1,1,0) = W_m11
  !  * Wigner(L2,L0,L1,0,1,-1) = W_m11 ( with parity )
  !  * Wigner(L1,L0,L2,0,-1,1) = W_0m1
  !  * Wigner(L1,L0,L2,0,1,-1) = W_0m1 ( with parity ) 
  ! grad: s=2
  !  * Wigner(L2,L0,L1,2,0,-2) = W_20
  !  * Wigner(L1,L0,L2,2,0,-2) = W_20
  ! curl: s=2
  !  * Wigner(L2,L0,L1,2,-1,-1) = Wigner(L1,L0,L2,1,1,-2) = W_11 
  !  * Wigner(L2,L0,L1,2,1,-3)  = Wigner(L1,L0,L2,3,-1,-2) = W_3m1
  !  * Wigner(L1,L0,L2,2,-1,-1) = W_2m1
  !  * Wigner(L1,L0,L2,2,1,-3)  = W_21 

  Sg0_12 = fac*(-L2*(L2+1)+L0*(L0+1)+L1*(L1+1))*W(0,0)
  Sg2_12 = fac*(-L2*(L2+1)+L0*(L0+1)+L1*(L1+1))*W(2,0)
  Sc2_12 = fac*LLsq(L0,0)*(LLsq(L1,1)*W(1,1)-LLsq(L1,2)*W(3,-1))
  Sg0_21 = fac*(-L1*(L1+1)+L0*(L0+1)+L2*(L2+1))*W(0,0)
  Sg2_21 = fac*(-L1*(L1+1)+L0*(L0+1)+L2*(L2+1))*W(2,0)
  Sc2_21 = fac*LLsq(L0,0)*(LLsq(L2,1)*W(2,-1)-LLsq(L2,2)*W(2,1))
  if(L1==2) Sc2_12 = fac*LLsq(L0,0)*LLsq(L1,1)*W(1,1) !W_3m1 = 0
  if(L2==2) Sc2_21 = fac*LLsq(L0,0)*LLsq(L2,1)*W(2,-1) !W_21 = 0

  !//// weight function ////!
  !* vanish
  ! even ---> grad: TB, EB, 
  !           curl: TT, TE, EE, BB
  ! odd  ---> grad: TT, TE, EE, BB
  !           curl: TB, EB 
  !* symmetric
  ! fgTT, fgEE, fcTT, fcEE, fcBB
  ! ggTT, ggEE, gcTT, gcEE, gcBB

  if (mod(L0+L1+L2,2)==0) then !* L0 + L1 + L2 = even 
    fgTT   = Sg0_21*fC(TT,L2) + Sg0_12*fC(TT,L1)
    fgTE12 = Sg0_21*fC(TE,L2) + Sg2_12*fC(TE,L1)
    fgTE21 = Sg0_12*fC(TE,L1) + Sg2_21*fC(TE,L2)
    fgEE   = Sg2_21*fC(EE,L2) + Sg2_12*fC(EE,L1)
    fgBB   = Sg2_21*fC(BB,L2) + Sg2_12*fC(BB,L1)
    fcTB12 = Sc2_12*fC(TE,L1)
    fcTB21 = Sc2_21*fC(TE,L2)
    fcEB12 = Sc2_12*fC(EE,L1)
    fcEB21 = Sc2_21*fC(EE,L2)
    !TE (decompose into A and B)
    ATE12  = oC(TT,L1)*oC(EE,L2)
    ATE21  = oC(TT,L2)*oC(EE,L1)
    BTE    = oC(TE,L1)*oC(TE,L2)
    ggTE12 = (fgTE12*ATE21 - fgTE21*BTE)/(ATE12*ATE21 - BTE**2)
    ggTE21 = (fgTE21*ATE12 - fgTE12*BTE)/(ATE12*ATE21 - BTE**2)

    D = 1d0
    if(L1==L2) D = 0.5d0
    !* auto/cross normalization
    Al(QTT,1) = Al(QTT,1) + D*alTT_ylm(fgTT,oC(TT,L1),oC(TT,L2))
    Al(QTE,1) = Al(QTE,1) + D*(fgTE12*ggTE12 + fgTE21*ggTE21)
    Al(QTB,3) = Al(QTB,3) + D*(fcTB12**2/oC(TT,L1)/oC(BB,L2) + fcTB21**2/oC(TT,L2)/oC(BB,L1))
    Al(QEE,1) = Al(QEE,1) + D*fgEE**2/oC(EE,L1)/oC(EE,L2)
    Al(QEB,3) = Al(QEB,3) + D*(fcEB12**2/oC(EE,L1)/oC(BB,L2) + fcEB21**2/oC(EE,L2)/oC(BB,L1))
    Al(QBB,1) = Al(QBB,1) + D*fgBB**2/oC(BB,L1)/oC(BB,L2)
    Il(QTTTE,1) = Il(QTTTE,1) + D*(fgTT*(ggTE12*oC(TE,L2)/oC(TT,L2)+ggTE21*oC(TE,L1)/oC(TT,L1)))
    Il(QTTEE,1) = Il(QTTEE,1) + D*fgTT*fgEE*oC(TE,L1)*oC(TE,L2)/oC(TT,L1)/oC(TT,L2)/oC(EE,L1)/oC(EE,L2) 
    Il(QTEEE,1) = Il(QTEEE,1) + D*fgEE*(ggTE12*oC(TE,L1)/oC(EE,L1)+ggTE21*oC(TE,L2)/oC(EE,L2))
    Il(QTBEB,3) = Il(QTBEB,1) + D*(fcTB12*fcEB12*oC(TE,L1)/oC(TT,L1)/oC(BB,L2)/oC(EE,L1) &
       + fcTB12*fcEB21*oC(TE,L2)/oC(TT,L2)/oC(BB,L1)/oC(EE,L2))

  else  !* !L1 + L2 + L3 = odd
    Sc0_12 = fac*LLsq(L0,0)*LLsq(L1,0)*W(-1,1)*2d0
    Sc0_21 = fac*LLsq(L0,0)*LLsq(L2,0)*W(0,-1)*2d0
    fcTT   = Sc0_21*fC(TT,L2) - Sc0_12*fC(TT,L1)
    fgTB12 = Sg2_12*fC(TE,L1)
    fgTB21 = Sg2_21*fC(TE,L2)
    fgEB12 = Sg2_12*fC(EE,L1)
    fgEB21 = Sg2_21*fC(EE,L2)
    fcTE12 = Sc0_21*fC(TE,L2) - Sc2_12*fC(TE,L1)
    fcTE21 = Sc0_12*fC(TE,L1) - Sc2_21*fC(TE,L2)
    fcEE   = Sc2_21*fC(EE,L2) - Sc2_12*fC(EE,L1)
    fcBB   = Sc2_21*fC(BB,L2) - Sc2_12*fC(BB,L1)
    !TE (decompose into A and B)
    ATE12 = oC(TT,L2)*oC(EE,L1)
    ATE21 = oC(TT,L1)*oC(EE,L2)
    BTE   = oC(TE,L1)*oC(TE,L2)
    gcTE12 = (fcTE12*ATE12 + fcTE21*BTE)/(ATE12*ATE21 - BTE**2)
    gcTE21 = (fcTE21*ATE21 + fcTE12*BTE)/(ATE12*ATE21 - BTE**2)

    D = 1d0
    if(L1==L2) D = 0.5d0
    Al(QTT,3) = Al(QTT,3) + D*fcTT**2/oC(TT,L1)/oC(TT,L2)
    Al(QTE,3) = Al(QTE,3) + D*(fcTE12*gcTE12 + fcTE21*gcTE21)
    Al(QTB,1) = Al(QTB,1) + D*(fgTB12**2/oC(TT,L1)/oC(BB,L2) + fgTB21**2/oC(TT,L2)/oC(BB,L1))
    Al(QEE,3) = Al(QEE,3) + D*fcEE**2/oC(EE,L1)/oC(EE,L2)
    Al(QEB,1) = Al(QEB,1) + D*(fgEB12**2/oC(EE,L1)/oC(BB,L2) + fgEB21**2/oC(EE,L2)/oC(BB,L1))
    Al(QBB,3) = Al(QBB,3) + D*fcBB**2/oC(BB,L1)/oC(BB,L2)
    Il(QTTTE,3) = Il(QTTTE,3) + D*fcTT*(gcTE12*oC(TE,L2)/oC(TT,L2) - gcTE21*oC(TE,L1)/oC(TT,L1))
    Il(QTEEE,3) = Il(QTEEE,3) + D*fcEE*(gcTE12*oC(TE,L1)/oC(EE,L1) - gcTE21*oC(TE,L2)/oC(EE,L2))
    Il(QTBEB,1) = Il(QTBEB,1) + D*(fgTB12*fgEB12*oC(TE,L1)/oC(TT,L1)/oC(BB,L2)/oC(EE,L1) & 
          + fgTB21*fgEB21*oC(TE,L2)/oC(TT,L2)/oC(BB,L1)/oC(EE,L2))
  end if

end subroutine gweight


function alTT_ylm(fTT,oCTT_L1,oCTT_L2) result(al)
  implicit none
  double precision, intent(in) :: fTT, oCTT_L1, oCTT_L2
  double precision :: al

  al = fTT**2/(oCTT_L1*oCTT_L2)

end function alTT_ylm


subroutine ini_factors(L,L1,W3j)
  implicit none
  !I/O
  integer, intent(in) :: L, L1
  double precision, intent(out) :: W3j(-1:3,-1:1)
  !internal
  integer :: i
  double precision :: p, W00, L1L0

  !* other Wigner 3j
  W00 = w00_ini(L,L1)
  W3j(0,0) = W00
  W3j(-1,1) = W00*W3j_ini_ratio(L,L1,-1,1)
  W3j(0,-1) = W00*W3j_ini_ratio(L,L1,0,-1)
  W3j(2,0)  = W00*W3j_ini_ratio(L,L1,2,0)
  W3j(1,1)  = W00*W3j_ini_ratio(L,L1,1,1)
  W3j(3,-1) = W00*W3j_ini_ratio(L,L1,3,-1)
  W3j(2,-1) = W00*W3j_ini_ratio(L,L1,2,-1)
  W3j(2,1)  = W00*W3j_ini_ratio(L,L1,2,1)

end subroutine ini_factors


subroutine sec_factors(L0,L1,Wsup,Wmid)
  implicit none
  !I/O
  integer, intent(in) :: L0, L1
  double precision, intent(in) :: Wsup(-1:3,-1:1)
  double precision, intent(out) :: Wmid(-1:3,-1:1)

  Wmid(0,0)  = 0d0
  Wmid(0,-1) = Wsup(0,-1)*W3j_sec(L0,L1,0,-1) 
  Wmid(-1,1) = Wsup(-1,1)*W3j_sec(L0,L1,-1,1) 
  Wmid(2,0)  = Wsup(2,0)*W3j_sec(L0,L1,2,0) 
  Wmid(1,1)  = Wsup(1,1)*W3j_sec(L0,L1,1,1) 
  Wmid(3,-1) = Wsup(3,-1)*W3j_sec(L0,L1,3,-1)
  Wmid(2,-1) = Wsup(2,-1)*W3j_sec(L0,L1,2,-1)
  Wmid(2,1)  = Wsup(2,1)*W3j_sec(L0,L1,2,1) 

end subroutine sec_factors


subroutine recursions(iL0,iL1,iL2,iM1,iM2,W3j)
  implicit none
  !I/O
  integer, intent(in) :: iL0, iL1, iL2, iM1, iM2
  double precision, intent(inout) :: W3j(3,-1:3,-1:1)
  !internal

  call w3j_recursion(iL0,iL1,iL2,0,0,W3j(:,0,0))
  call w3j_recursion(iL0,iL1,iL2,-1,1,W3j(:,-1,1))
  call w3j_recursion(iL0,iL1,iL2,0,-1,W3j(:,0,-1))
  call w3j_recursion(iL0,iL1,iL2,2,0,W3j(:,2,0))
  call w3j_recursion(iL0,iL1,iL2,1,1,W3j(:,1,1))
  call w3j_recursion(iL0,iL1,iL2,3,-1,W3j(:,3,-1))
  call w3j_recursion(iL0,iL1,iL2,2,-1,W3j(:,2,-1))
  call w3j_recursion(iL0,iL1,iL2,2,1,W3j(:,2,1))

end subroutine recursions

subroutine shift(W3j)
  implicit none
  double precision, intent(inout) :: W3j(3,-1:3,-1:1)
  integer :: l1, l2

  do l1 = -1, 3
    do l2 = -1, 1
      W3j(1,l1,l2) = W3j(2,l1,l2)
      W3j(2,l1,l2) = W3j(3,l1,l2)
    end do
  end do

end subroutine shift


end module nldd_ylm

