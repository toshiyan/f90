!////////////////////////////////////////////////////!
! * Kernel Functions of Fast Al Algorithm
!////////////////////////////////////////////////////!

module nldd_lens
  use myconst, only: pi
  use funcs,   only: WIGd_INI!, WIGD_RECURSION
  use myutils, only: gauss_legendre_params, gl_initialize, gl_finalize, savetxt, linspace
  use nldd_kernel, only: zeta, zeta_l
  implicit none

  !local parameters
  double precision, parameter :: cx(1:2) = [1d0,-1d0], fourpi = 4d0*pi
  private cx, fourpi

  private pi
  private wigd_ini
  private gauss_legendre_params, gl_initialize, gl_finalize, savetxt, linspace
  private zeta, zeta_l

contains


subroutine Kernels(rL,eL,WA,WB,X,kernel)
!Kernels of lensing reconstruction noise
  implicit none
  !I/O
  character(*), intent(in) :: kernel
  integer, intent(in) :: rL(1:2), eL(1:2)
  double precision, intent(in), dimension(rL(1):rL(2)) :: WA, WB
  double precision, intent(out) :: X(:,:)
  !internal
  type(gauss_legendre_params) :: GL
  character(LEN=4) :: selec
  integer :: i, l
  double precision :: pm, mu, Ip, Im, al, c1_inv, c2p, c2m, c3
  double precision, dimension(rL(1):rL(2)) :: al0, alp2, alm2
  double precision, dimension(2) :: d11_sup, d11_mid, d11_inf, C
  double precision, dimension(2) :: d00_sup, d00_mid, d00_inf
  double precision, dimension(2) :: d10_sup, d10_mid, d10_inf
  double precision, dimension(2) :: ZA00, ZA10, ZA20, ZA21, ZA22, ZA30, ZA32
  double precision, dimension(2) :: ZB00, ZB10, ZB11, ZB31, ZB32, ZB33, ZB21

  !* initialize
  if (kernel=='Sp'.or.kernel=='Gp') pm = 1d0
  if (kernel=='Sm'.or.kernel=='Gm') pm = -1d0

  selec = kernel
  if (kernel=='Sp'.or.kernel=='Sm') selec = 'Spm'
  if (kernel=='Gp'.or.kernel=='Gm') selec = 'Gpm'

  X = 0d0

  do l = rL(1), rL(2)
    al = dble(l)
    al0(l) = dsqrt(al*(al+1))
    alp2(l) = dsqrt((al-2)*(al+3))
    alm2(l) = dsqrt((al+2)*(al-1))
  end do

  !* GL quadrature
  call gl_initialize(GL,int((3*max(rL(2),eL(2))+1)/2),1d-15)

  do i = 1, GL%n !Gauss-Legendre Integration
    mu = GL%z(i)
    select case(selec)
    case ('S0')
      call ZETA(0,0,rL,WA,mu,ZA00)
      call ZETA(1,1,rL,WB*al0**2,mu,ZB11)
      Ip = ZA00(1)*ZB11(1)
      Im = ZA00(1)*ZB11(2)
    case ('G0')
      call ZETA(1,0,rL,WA*al0,mu,ZA10)
      call ZETA(1,0,rL,WB*al0,mu,ZB10)
      Ip = -ZA10(1)*ZB10(1)
      Im = -Ip
    case ('Sc')
      call Zeta(2,0,rL,WA,mu,ZA20)
      call Zeta(1,1,rL,WB*al0*alm2,mu,ZB11)
      call Zeta(3,1,rL,WB*al0*alp2,mu,ZB31)
      Ip = ZA20(1)*(ZB11(2)+ZB31(1))*0.5d0
      Im = ZA20(1)*(ZB31(2)+ZB11(1))*0.5d0
    case ('Gc')
      call Zeta(3,0,rL,WA*alp2,mu,ZA30)
      call Zeta(1,0,rL,WA*alm2,mu,ZA10)
      call Zeta(2,1,rL,WB*al0,mu,ZB21)
      Ip = -(ZA30(1)*ZB21(2)+ZA10(1)*ZB21(1))*0.5d0
      Im = -(ZA30(1)*ZB21(1)+ZA10(1)*ZB21(2))*0.5d0
    case ('Spm')
      call Zeta(2,2,rL,WA,mu,ZA22)
      call Zeta(3,3,rL,WB*alp2**2,mu,ZB33)
      call Zeta(3,1,rL,WB*alp2*alm2,mu,ZB31)
      call Zeta(1,1,rL,WB*alm2**2,mu,ZB11)
      Ip = ((ZB33(1)+ZB11(1))*ZA22(1)+pm*2d0*ZB31(2)*ZA22(2))/4d0
      Im = (pm*(ZB33(2)+ZB11(2))*ZA22(2)+2d0*ZB31(1)*ZA22(1))/4d0
    case ('Gpm')
      call Zeta(3,2,rL,WA*alp2,mu,ZA32)
      call Zeta(2,1,rL,WA*alm2,mu,ZA21)
      call Zeta(3,2,rL,WB*alp2,mu,ZB32)
      call Zeta(2,1,rL,WB*alm2,mu,ZB21)
      Ip = -(ZA21(1)*ZB32(1)+ZA32(1)*ZB21(1)+pm*(ZA32(2)*ZB32(2)+ZA21(2)*ZB21(2)))/4d0
      Im = (ZA32(1)*ZB32(1)+ZA21(1)*ZB21(1)-pm*(ZA21(2)*ZB32(2)+ZA32(2)*ZB21(2)))/4d0
    case default
      stop 'error: no kernel'
    end select
    do l = 1, eL(2) ! loop for l
      al = dble(l)
      if(l==1) then 
        d00_sup = wigd_ini(0,0,mu)
        d00_inf = d00_mid
        d00_mid = d00_sup
        d00_sup(1:2) = mu*d00_mid(1:2)*(al*(2d0*al-1d0))/al**2
        d11_sup = wigd_ini(1,1,mu)
        d10_sup = wigd_ini(1,0,mu)
      else
        !* M1=M2=0
        c1_inv = (2d0*al-1d0)/al
        c3 = (al-1d0)**2/((al-1d0)*(2d0*al-1d0))
        d00_sup(1) = (mu*d00_mid(1) - c3*d00_inf(1))*c1_inv
        d00_sup(2) = (mu*d00_mid(2) - c3*d00_inf(2))*c1_inv
        !* M1=M2=1
        c1_inv = al*(2d0*al-1d0)/(al**2-1d0)
        c2p = 1d0/((aL-1d0)*aL) - mu
        c2m = -1d0/((aL-1d0)*aL) - mu
        c3 = ((al-1d0)**2-1d0)/((al-1d0)*(2d0*al-1d0))
        d11_sup(1) = -(c2p*d11_mid(1) + c3*d11_inf(1))*c1_inv
        d11_sup(2) = -(c2m*d11_mid(2) + c3*d11_inf(2))*c1_inv
        !* M1=1,M2=0
        c1_inv = al*(2d0*al-1d0)/(dsqrt((al**2-1d0))*al)
        c3 = dsqrt(((al-1d0)**2-1d0))*(al-1d0)/((al-1d0)*(2d0*al-1d0))
        d10_sup(1) = (mu*d10_mid(1) - c3*d10_inf(1))*c1_inv
        d10_sup(2) = (mu*d10_mid(2) - c3*d10_inf(2))*c1_inv
      end if
      C(:) = (Ip*d11_sup(1)+cx(:)*Im*d11_sup(2))*dble(l*(l+1))
      X(1:2,l) = X(1:2,l) + C(1:2)*GL%w(i)*pi
      d00_inf = d00_mid
      d00_mid = d00_sup
      d10_inf = d10_mid
      d10_mid = d10_sup
      d11_inf = d11_mid
      d11_mid = d11_sup
    end do
  end do

  call gl_finalize(GL)

end subroutine Kernels


subroutine AlTT(rL,eL,Agg,Acc,fC,OCTT)
  implicit none
  !I/O
  integer, intent(in) :: rL(2), eL(2)
  double precision, intent(in), dimension(:) :: fC, OCTT
  double precision, intent(out), dimension(:) :: Agg, Acc
  !internal
  double precision, dimension(rL(1):rL(2)) :: W1, W2
  double precision, dimension(2,eL(2)) :: S0, G0

  write(*,*) 'TT'
  if (eL(1)<1)  stop 'error (altt): lmin<1'

  W1 = 1d0 / OCTT(rL(1):rL(2))
  W2 = W1 * fC(rL(1):rL(2))**2
  S0 = 0d0
  call Kernels(rL,eL,W1,W2,S0,'S0')

  W2 = W1 * fC(rL(1):rL(2))
  G0 = 0d0
  call Kernels(rL,eL,W2,W2,G0,'G0')

  Agg(eL(1):eL(2)) = 1d0/(S0(1,eL(1):eL(2))+G0(1,eL(1):eL(2)))
  Acc(eL(1):eL(2)) = 1d0/(S0(2,eL(1):eL(2))+G0(2,eL(1):eL(2)))

end subroutine AlTT


subroutine AlTE(rL,eL,Agg,Acc,fC,OCTT,OCEE)
  implicit none
  !I/O
  integer, intent(in) :: rL(2), eL(2)
  double precision, intent(in), dimension(:) :: fC, OCTT, OCEE
  double precision, intent(out) :: Agg(:), Acc(:)
  !internal
  double precision, dimension(rL(1):rL(2)) :: W1, W2
  double precision, dimension(2,eL(2)) :: S0, Sp, Gc

  write(*,*) 'TE'
  if (eL(1)<1)  stop 'error (alte): lmin<1'

  W1 = 1d0/OCTT(rL(1):rL(2))
  W2 = fC(rL(1):rL(2))**2/OCEE(rL(1):rL(2))
  S0 = 0d0
  call kernels(rL,eL,W1,W2,S0,'S0')

  W1 = 1d0/OCEE(rL(1):rL(2))
  W2 = fC(rL(1):rL(2))**2/OCTT(rL(1):rL(2))
  Sp = 0d0
  call kernels(rL,eL,W1,W2,Sp,'Sp')

  W1 = fC(rL(1):rL(2))/OCTT(rL(1):rL(2))
  W2 = fC(rL(1):rL(2))/OCEE(rL(1):rL(2))
  Gc = 0d0
  call kernels(rL,eL,W1,W2,Gc,'Gc')

  Agg(eL(1):eL(2)) = 1d0/(S0(1,eL(1):eL(2))+Sp(1,eL(1):eL(2))+2d0*Gc(1,eL(1):eL(2)))
  Acc(eL(1):eL(2)) = 1d0/(S0(2,eL(1):eL(2))+Sp(2,eL(1):eL(2))+2d0*Gc(2,eL(1):eL(2)))

end subroutine AlTE


subroutine AlTB(rL,eL,Agg,Acc,fC,OCTT,OCBB)
  implicit none
  !I/O
  integer, intent(in) :: rL(2), eL(2)
  double precision, intent(in), dimension(:) :: fC, OCTT, OCBB
  double precision, intent(out) :: Agg(:), Acc(:)
  !internal
  double precision, dimension(rL(1):rL(2)) :: W1, W2
  double precision, dimension(2,eL(2)) :: Sm

  write(*,*) 'TB'
  if (eL(1)<1)  stop 'error (altb): lmin<1'

  W1 = 1d0/OCBB(rL(1):rL(2))
  W2 = fC(rL(1):rL(2))**2/OCTT(rL(1):rL(2))
  Sm = 0d0
  call kernels(rL,eL,W1,W2,Sm,'Sm')

  Agg(eL(1):eL(2)) = 1d0/Sm(1,eL(1):eL(2))
  Acc(eL(1):eL(2)) = 1d0/Sm(2,eL(1):eL(2))

end subroutine AlTB


subroutine AlEE(rL,eL,Agg,Acc,fC,OCEE)
! Al: 1 -> number of estimator power
!   : 2 -> number of multipoles
  implicit none
  !I/O
  integer, intent(in) :: rL(2), eL(2)
  double precision, intent(in), dimension(:) :: fC, OCEE
  double precision, intent(out) :: Agg(:), Acc(:)
  !internal
  double precision, dimension(rL(1):rL(2)) :: W1, W2
  double precision, dimension(2,eL(2)) :: Sp, Gp

  write(*,*) 'EE'
  if (eL(1)<1)  stop 'error (alee): lmin<1'

  W1 = 1d0/OCEE(rL(1):rL(2))
  W2 = W1 * fC(rL(1):rL(2))**2
  Sp = 0d0
  call kernels(rL,eL,W1,W2,Sp,'Sp')

  W2 = W1 * fC(rL(1):rL(2))
  Gp = 0d0
  call kernels(rL,eL,W2,W2,Gp,'Gp')

  Agg(eL(1):eL(2)) = 1d0/(Sp(1,eL(1):eL(2))+Gp(1,eL(1):eL(2)))
  Acc(eL(1):eL(2)) = 1d0/(Sp(2,eL(1):eL(2))+Gp(2,eL(1):eL(2)))

end subroutine AlEE


subroutine AlEB(rL,eL,Agg,Acc,fCEE,OCEE,OCBB)
  implicit none
  !I/O
  integer, intent(in) :: rL(2), eL(2)
  double precision, intent(in), dimension(:) :: fCEE, OCEE, OCBB
  double precision, intent(out) :: Agg(:), Acc(:)
  !internal
  double precision, dimension(rL(1):rL(2)) :: W1, W2
  double precision, dimension(2,eL(2)) :: Sm

  write(*,*) 'EB'
  if (eL(1)<1)  stop 'error (aleb): lmin<1'

  W1 = 1d0/OCBB(rL(1):rL(2))
  W2 = fCEE(rL(1):rL(2))**2 / OCEE(rL(1):rL(2))
  Sm = 0d0
  call kernels(rL,eL,W1,W2,Sm,'Sm')

  Agg(eL(1):eL(2)) = 1d0/Sm(1,eL(1):eL(2))
  Acc(eL(1):eL(2)) = 1d0/Sm(2,eL(1):eL(2))

end subroutine AlEB


subroutine AlBB(rL,eL,Agg,Acc,fC,OCBB)
  implicit none
  !I/O
  integer, intent(in) :: rL(2), eL(2)
  double precision, intent(in), dimension(:) :: fC, OCBB
  double precision, intent(out) :: Agg(:), ACc(:)
  !internal
  double precision, dimension(rL(1):rL(2)) :: W1, W2
  double precision, dimension(2,eL(2)) :: Sp, Gp

  write(*,*) 'BB'
  if (eL(1)<1)  stop 'error (albb): lmin<1'

  W1 = 1d0/OCBB(rL(1):rL(2))
  W2 = W1 * fC(rL(1):rL(2))**2
  call kernels(rL,eL,W1,W2,Sp,'Sp')

  W2 = W1 * fC(rL(1):rL(2))
  call kernels(rL,eL,W2,W2,Gp,'Gp')

  Agg(eL(1):eL(2)) = 1d0/(Sp(1,eL(1):eL(2))+Gp(1,eL(1):eL(2)))
  Acc(eL(1):eL(2)) = 1d0/(Sp(2,eL(1):eL(2))+Gp(2,eL(1):eL(2)))


end subroutine AlBB


subroutine IlTTTE(rL,eL,Ig,Ic,fCTT,fCTE,OCTT,OCEE,OCTE)
  implicit none
  !I/O
  integer, intent(in) :: rL(2), eL(2)
  double precision, intent(in), dimension(:) :: fCTT,fCTE,OCTT,OCEE,OCTE
  double precision, intent(out) :: Ig(:), Ic(:)
  !internal
  integer :: l
  double precision, dimension(rL(1):rL(2)) :: W1, W2
  double precision, dimension(2,eL(2)) :: S0, Gc, G0, Sc

  write(*,*) 'TTTE'
  if (eL(1)<1)  stop 'error (ilttte): lmin<1'
  Ig=0d0 ;  Ic=0d0

  do l = rL(1), rL(2)
    W1(l) = 1d0/OCTT(l)
    W2(l) = fCTT(l)*fCTE(l)*OCTE(l)/(OCTT(l)*OCEE(l))
  end do
  call kernels(rL,eL,W1,W2,S0,'S0')

  do l = rL(1), rL(2)
    W1(l) = fCTE(l)/OCTT(l)
    W2(l) = fCTT(l)*OCTE(l)/(OCTT(l)*OCEE(l))
  end do
  call kernels(rL,eL,W1,W2,Gc,'Gc')

  do l = rL(1), rL(2)
    W1(l) = fCTE(l)*OCTE(l)/(OCTT(l)*OCEE(l))
    W2(l) = fCTT(l)/OCTT(l)
  end do
  call kernels(rL,eL,W1,W2,G0,'G0')

  do l = rL(1), rL(2)
    W1(l) = OCTE(l)/(OCTT(l)*OCEE(l))
    W2(l) = fCTT(l)*fCTE(l)/OCTT(l)
  end do
  call kernels(rL,eL,W1,W2,Sc,'Sc')

  do l = eL(1), eL(2)
    Ig(l) = S0(1,l)+Gc(1,l)+G0(1,l)+Sc(1,l)
    Ic(l) = S0(2,l)+Gc(2,l)+G0(2,l)+Sc(2,l)
  end do

end subroutine IlTTTE


subroutine IlTTEE(rL,eL,Ig,Ic,fCTT,fCEE,OCTT,OCEE,OCTE)
  implicit none
  !I/O
  integer, intent(in) :: rL(2), eL(2)
  double precision, intent(in), dimension(:) :: fCTT,fCEE,OCTT,OCEE,OCTE
  double precision, intent(out) :: Ig(:), Ic(:)
  !internal
  integer :: l
  double precision, dimension(rL(1):rL(2)) :: W1, W2
  double precision, dimension(2,eL(2)) :: Sc, Gc

  write(*,*) 'TTEE'
  if (eL(1)<1)  stop 'error (ilttee): lmin<1'
  Ig=0d0 ;  Ic=0d0 ;  W1=0d0 ;  W2=0d0

  do l = rL(1), rL(2)
    W1(l) = OCTE(l)/(OCTT(l)*OCEE(l))
    W2(l) = fCTT(l)*fCEE(l)*OCTE(l)/(OCTT(l)*OCEE(l))
  end do
  call kernels(rL,eL,W1,W2,Sc,'Sc')

  do l = rL(1), rL(2)
    W1(l) = fCEE(l)*OCTE(l)/(OCTT(l)*OCEE(l))
    W2(l) = fCTT(l)*OCTE(l)/(OCTT(l)*OCEE(l))
  end do
  call kernels(rL,eL,W1,W2,Gc,'Gc')

  Ig(eL(1):eL(2)) = Sc(1,eL(1):eL(2)) + Gc(1,eL(1):eL(2))
  Ic = 0d0

end subroutine IlTTEE


subroutine IlTEEE(rL,eL,Ig,Ic,fCEE,fCTE,OCTT,OCEE,OCTE)
  implicit none
  !I/O
  integer, intent(in) :: rL(2), eL(2)
  double precision, intent(in), dimension(:) :: fCEE,fCTE,OCTT,OCEE,OCTE
  double precision, intent(out) :: Ig(:), Ic(:)
  !internal
  integer :: l
  double precision, dimension(rL(1):rL(2)) :: W1, W2
  double precision, dimension(2,eL(2)) :: Sc,Gp,Gc,Sp

  write(*,*) 'TEEE'
  if (eL(1)<1)  stop 'error (ilteee): lmin<1'
  Ig=0d0 ;  Ic=0d0

  do l = rL(1), rL(2)
    W1(l) = OCTE(l)/(OCTT(l)*OCEE(l))
    W2(l) = fCTE(l)*fCEE(l)/OCEE(l)
  end do
  call kernels(rL,eL,W1,W2,Sc,'Sc')

  do l = rL(1), rL(2)
    W1(l) = fCTE(l)*OCTE(l)/(OCTT(l)*OCEE(l))
    W2(l) = fCEE(l)/OCEE(l)
  end do
  call kernels(rL,eL,W1,W2,Gp,'Gp')

  do l = rL(1), rL(2)
    W1(l) = fCEE(l)*OCTE(l)/(OCTT(l)*OCEE(l))
    W2(l) = fCTE(l)/OCEE(l)
  end do
  call kernels(rL,eL,W1,W2,Gc,'Gc')

  do l = rL(1), rL(2)
    W1(l) = 1d0/OCEE(l)
    W2(l) = fCTE(l)*fCEE(l)*OCTE(l)/(OCTT(l)*OCEE(l))
  end do
  call kernels(rL,eL,W1,W2,Sp,'Sp')

  do l = eL(1), eL(2)
    Ig(l) = (Sc(1,l)+Gp(1,l)+Gc(1,l)+Sp(1,l))
    Ic(l) = (Sc(2,l)+Gp(2,l)+Gc(2,l)+Sp(2,l))
  end do

end subroutine IlTEEE


subroutine IlTBEB(rL,eL,Ig,Ic,fCEE,fCBB,fCTE,OCTT,OCEE,OCBB,OCTE)
  implicit none
  !I/O
  integer, intent(in) :: rL(2), eL(2)
  double precision, intent(in), dimension(:) :: fCEE,fCBB,fCTE,OCTT,OCEE,OCTE,OCBB
  double precision, intent(out) :: Ig(:), Ic(:)
  !internal
  integer :: l
  double precision, dimension(2,eL(2)) :: Gm, Sm
  double precision, dimension(rL(1):rL(2)) :: W1, W2

  write(*,*) 'TBEB'
  if (eL(1)<1)  stop 'error (iltbeb): lmin<1'
  Ig=0d0 ;  Ic=0d0

  do l = rL(1), rL(2)
    W1(l) = fCTE(l)*OCTE(l)/(OCTT(l)*OCEE(l))
    W2(l) = fCBB(l)/OCBB(l)
  end do
  call kernels(rL,eL,W1,W2,Gm,'Gm')

  do l = rL(1), rL(2)
    W1(l) = 1d0/OCBB(l)
    W2(l) = fCTE(l)*fCEE(l)*OCTE(l)/(OCTT(l)*OCEE(l))
  end do
  call kernels(rL,eL,W1,W2,Sm,'Sm')

  Ig(eL(1):eL(2)) = Sm(1,eL(1):eL(2)) + Gm(1,eL(1):eL(2))
  Ic(eL(1):eL(2)) = Sm(2,eL(1):eL(2)) + Gm(2,eL(1):eL(2))

end subroutine IlTBEB


end module nldd_lens

