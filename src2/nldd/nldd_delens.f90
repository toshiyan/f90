!////////////////////////////////////////////////////!
! * Delensing 
!////////////////////////////////////////////////////!

module nldd_delens
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


subroutine conv_egrad(oL,eL,WE,WP,X)
! * Kernel of lensing B-mode power spectrum
!
! The Kernel function is given by
!
!   (1/(2l+1)) \sum_{L1L2} (S^-_{lL1L2})^2 WE_L1 WP_L2
!
! see Eqs.(2.7) and (A.11) in Namikawa & Nagata 2014
!
! inputs  : oL --- multipole range of X
!         : eL --- multipole range of WE and WP
!         : WE, WP --- functions to be convolved
! outputs : X --- kernel function

  implicit none
  !I/O
  integer, intent(in) :: oL(1:2), eL(1:2)
  double precision, intent(in), dimension(eL(1):eL(2)) :: WE, WP
  double precision, intent(out) :: X(:)
  !internal
  type(gauss_legendre_params) :: GL
  integer :: i, l
  double precision :: mu, al, Ip, Im, c1_inv, c2p, c2m, c3
  double precision, dimension(2) :: ZE33, ZE31, ZE11, ZP, d22_sup, d22_mid, d22_inf
  double precision, dimension(eL(1):eL(2)) :: al0, alm, alp

  do l = eL(1), eL(2)
    al = dble(l)
    al0(l) = dsqrt(al*(al+1))
    alp(l) = dsqrt((al-2)*(al+3))
    alm(l) = dsqrt((al+2)*(al-1))
  end do

  X = 0d0

  !* GL quadrature
  call gl_initialize(GL,int((3*eL(2)+1)/2),1d-15)

  do i = 1, GL%n
    mu = GL%z(i)
    call ZETA(3,3,eL,WE*alp**2,mu,ZE33)
    call ZETA(3,1,eL,WE*alp*alm,mu,ZE31)
    call ZETA(1,1,eL,WE*alm**2,mu,ZE11)
    call ZETA(1,1,eL,WP*al0**2,mu,ZP)
    d22_sup = 0d0 ;  d22_mid = 0d0 ;  d22_inf = 0d0
    do l = 2, oL(2) !recursion of wigner d function
      if(l==2) then 
        d22_sup = wigd_ini(2,2,mu)
      else
        al = dble(l)
        c1_inv = (al*(2d0*al-1d0))/(al**2-4d0)
        c2p = 4d0/((al-1d0)*al) - mu
        c2m = -4d0/((al-1d0)*al) - mu
        c3 = ((al-1d0)**2-4d0)/((al-1d0)*(2d0*al-1d0))
        d22_sup(1) = -(c2p*d22_mid(1) + c3*d22_inf(1))*c1_inv
        d22_sup(2) = -(c2m*d22_mid(2) + c3*d22_inf(2))*c1_inv
      end if
      Ip = (ZE33(1)*ZP(1)+2*ZE31(1)*ZP(2)+ZE11(1)*ZP(1))
      Im = (ZE33(2)*ZP(2)+2*ZE31(2)*ZP(1)+ZE11(2)*ZP(2))
      X(l) = X(l) + (Ip*d22_sup(1)-Im*d22_sup(2))*GL%w(i)*pi*0.25d0
      d22_inf = d22_mid
      d22_mid = d22_sup
    end do
  end do

  call gl_finalize(GL)

end subroutine conv_egrad


subroutine RES_CLBB(oL,eL,CE,Cp,CB,WE,Wp,NE,Np,f)
!* residual ClBB = ClBB^lin - ClBB^est
!
  implicit none
! [inputs]  
!   oL --- multipole range of residual ClBB
!   eL --- multipole range of delensing
!   CE, Cp --- power spectrum of E-mode and lensing pontential
  integer, intent(in) :: oL(2), eL(2)
  double precision, intent(in), dimension(:) :: CE, Cp
!
! (optional)
!   WE, Wp --- Wiener filters of E-mode and lensing potential
!   NE, Np --- Noise power spectra of E-mode and lensing potential
!   f --- ClBB output filename
  double precision, intent(in), dimension(:), optional :: WE, Wp, NE, Np
  character(*), intent(in), optional :: f
!
! [outputs]
!   CB --- residual B-mode
  double precision, intent(out) :: CB(:)
!
! [internal]
  integer :: l
  double precision :: CBW(oL(2)), CBL(oL(2))
  double precision, dimension(eL(1):eL(2)) :: A, B

! Wiener filtered Lensing B-mode power spectrum
  CBW = 0d0

  if (present(WE).and.present(Wp)) then
    A = CE(eL(1):eL(2)) * WE(eL(1):eL(2))
    B = Cp(eL(1):eL(2)) * Wp(eL(1):eL(2))
    call conv_egrad(oL,eL,A,B,CBW)
  end if

  if (present(NE).and.present(Np)) then
    do l = eL(1), eL(2)
      A(l) = CE(l)**2/(CE(l)+NE(l))
      B(l) = Cp(l)**2/(Cp(l)+Np(l))
    end do
    call conv_egrad(oL,eL,A,B,CBW)
  end if

! Lensing B-mode power spectrum
  call conv_egrad(oL,eL,CE(eL(1):eL(2)),Cp(eL(1):eL(2)),CBL)

  CB = CBL
  CB(oL(1):oL(2)) = CBL(oL(1):oL(2)) - CBW(oL(1):oL(2))

  if (present(f))  call savetxt(trim(f),linspace(1,eL(2),eL(2)),CB)

end subroutine RES_CLBB


subroutine CLBB_LIN(oL,eL,CE,Cg,CB)
! * Lensing B-mode power spectrum as a convolution of ClEE and Clpp
  implicit none
!
! [input]
!   oL --- multipole range of residual ClBB
!   eL --- multipole range of delensing
!   CE, Cg --- power spectrum of E-mode and lensing pontential
  integer, intent(in) :: oL(2), eL(2)
  double precision, intent(in), dimension(:) :: CE, Cg
!
! [output]
!   CB --- residual B-mode
  double precision, intent(out) :: CB(:)

  call conv_egrad(oL,eL,CE(eL(1):eL(2)),Cg(eL(1):eL(2)),CB)

end subroutine CLBB_LIN


subroutine CLBB_EST_W(oL,eL,CE,Cp,WE,Wp,Cl)
! * Estimate of lensing B-mode power spectrum (Wiener filters as inputs)
  implicit none
!
! [input]
!   oL --- multipole range of residual ClBB
!   eL --- multipole range of delensing
!   CE --- E-mode power spectrum
!   Cp --- lensing potential power spectrum
!   WE --- E-mode Wiener filter
!   Wp --- lensing potential Wiener filter
  integer, intent(in) :: oL(2), eL(2)
  double precision, intent(in), dimension(:) :: CE, Cp, WE, Wp
!
! [output]
!   Cl --- estimated lensing B-mode power spectrum
  double precision, intent(out) :: Cl(:)
!
! [internal]
  double precision, dimension(eL(1):eL(2)) :: W1, W2

  W1 = CE(eL(1):eL(2))*WE(eL(1):eL(2))
  W2 = Cp(eL(1):eL(2))*Wp(eL(1):eL(2))
  call conv_egrad(oL,eL,W1,W2,Cl)

end subroutine CLBB_EST_W


subroutine CLBB_EST(oL,eL,CE,Cp,NE,Np,Cl)
! * Estimate of lensing B-mode power spectrum (Noise power spectra as inputs)
  implicit none
!
! [input]
!   oL --- multipole range of residual ClBB
!   eL --- multipole range of delensing
!   CE --- E-mode power spectrum
!   Cp --- lensing potential power spectrum
!   NE --- E-mode noise power spectrum
!   Np --- lensing potential noise power spectrum
  integer, intent(in) :: oL(2), eL(2)
  double precision, intent(in), dimension(:) :: CE, Cp, NE, Np
!
! [output]
!   Cl --- estimated lensing B-mode power spectrum
  double precision, intent(out) :: Cl(:)
!
! [internal]
  integer :: l
  double precision, dimension(eL(1):eL(2)) :: WE, Wp

  do l = eL(1), eL(2)
    WE(l) = CE(l)**2/(CE(l)+NE(l))
    Wp(l) = Cp(l)**2/(Cp(l)+Np(l))
  end do

  call conv_egrad(oL,eL,WE,Wp,Cl)

end subroutine CLBB_EST


subroutine Delensing_Bias(rL,eL,LE,LB,Cp,CE,DBL,NP1,Np,NP2,f)
!* Computing total corrections of the delensing bias
  implicit none
  integer, intent(in) :: rL(2), eL(2)
  character(*), intent(in), optional :: f
  double precision, intent(in), dimension(:) :: LE,LB,Cp,CE,Np,NP1
  double precision, intent(in), dimension(:), optional :: NP2
  double precision, intent(out) :: DBL(:)
  !internal
  integer :: l
  double precision, dimension(rL(2)) :: CBB_est, DLs

  !* delensing bias (Dl and sum of bias terms)
  DLs = 0d0
  DBl = 0d0
  call DLBB(rL,eL,CE,LE,LE+NP1,LB+NP1,Cp,Np,DLS,LE+NP2)
  do l = rL(1), rL(2)
    DBl(l) = DLS(l)*(-2d0*LB(l)+2d0*CBB_est(l)+DLS(l)*(LB(l)+NP1(l)))
  end do

  !* output result
  if(present(f))  call savetxt(f,linspace(1,rL(2),rL(2)),CBB_est,DLS,DBl)

end subroutine Delensing_Bias


subroutine DlBB(rL,eL,CE,LE,OE1,OB1,Cp,Np,Cl,OE2)
!* Computing delensing bias power spectrum (Dl)
  implicit none
!
! [input]
!   eL --- multipole range of residual ClBB
!   rL --- multipole range of delensing
!   CE --- E-mode power spectrum
!   Cp --- lensing potential power spectrum
!   NE --- E-mode noise power spectrum
!   Np --- lensing potential noise power spectrum
  integer, intent(in) :: rL(2), eL(2)
  double precision, intent(in), dimension(:) :: CE, LE, OE1, OB1, Cp, Np
!
! (optional)
  double precision, intent(in), dimension(:), optional :: OE2
!
! [output]
!   Cl --- delensing bias
  double precision, intent(out) :: Cl(:)
!
! [internal]
  integer :: l
  double precision, dimension(eL(1):eL(2)) :: WE, Wp

  do l = eL(1), eL(2)
    WE(l) = 0d0
    if(l>=rL(1)) WE(l) = (LE(l)/OE1(l))*CE(l)
    if(present(OE2)) WE(l) = WE(l)*LE(l)/OE2(l)
    Wp(l) = Cp(l)*Np(l)/(Cp(l)+Np(l))
  end do
  call conv_egrad(rL,eL,WE,Wp,Cl)

  Cl(rL(1):rL(2)) = Cl(rL(1):rL(2))/OB1(rL(1):rL(2))

end subroutine DlBB


!//// covariance of lensed/delensed ClBB ////!

subroutine dCBdCE(oL,eL,Cp,X,f)
! * Derivative of B-mode Cl in terms of the E-mode Cl
!
  implicit none
! [input]  
!   oL --- multipole range of residual ClBB
!   eL --- multipole range of delensing
!   Cp --- power spectrum of the lensing pontential
  integer, intent(in) :: oL(1:2), eL(1:2)
  double precision, intent(in), dimension(eL(1):eL(2)) :: Cp
!
! (optional)
!   f --- output filename
  character(*), intent(in), optional :: f
!
! [output]
!   X --- derivative of ClBB in terms of Clpp
  double precision, intent(out) :: X(eL(2),oL(2))
!
! [internal]
  type(gauss_legendre_params) :: GL
  integer :: i, l
  double precision :: mu, al, c1_inv, c2p, c2m, c3
  double precision, dimension(2) :: ZP, d22_sup, d22_mid, d22_inf
  double precision, dimension(eL(2)) :: Ip, Im
  double precision, dimension(eL(2),2) :: ZE33, ZE31, ZE11
  double precision, dimension(eL(1):eL(2)) :: al0, alm, alp

  do l = eL(1), eL(2)
    al = dble(l)
    al0(l) = dsqrt(al*(al+1))
    alp(l) = dsqrt((al-2)*(al+3))
    alm(l) = dsqrt((al+2)*(al-1))
  end do

  !* Gauss-Legendre Integration
  call gl_initialize(GL,int((3*eL(2)+1)/2),1d-15)
  X = 0d0
  do i = 1, GL%n
    mu = GL%z(i)
    call ZETA_L(3,3,eL,alp**2,mu,ZE33)
    call ZETA_L(3,1,eL,alp*alm,mu,ZE31)
    call ZETA_L(1,1,eL,alm**2,mu,ZE11)
    call ZETA(1,1,eL,Cp*al0**2,mu,ZP)
    Ip = (ZE33(:,1)*ZP(1)+2*ZE31(:,1)*ZP(2)+ZE11(:,1)*ZP(1))
    Im = (ZE33(:,2)*ZP(2)+2*ZE31(:,2)*ZP(1)+ZE11(:,2)*ZP(2))
    d22_sup = 0d0 ;  d22_mid = 0d0 ;  d22_inf = 0d0
    do l = 2, oL(2) !recursion of wigner d function
      if(l==2) then 
        d22_sup = wigd_ini(2,2,mu)
      else
        al = dble(l)
        c1_inv = (al*(2d0*al-1d0))/(al**2-4d0)
        c2p = 4d0/((al-1d0)*al) - mu
        c2m = -4d0/((al-1d0)*al) - mu
        c3 = ((al-1d0)**2-4d0)/((al-1d0)*(2d0*al-1d0))
        d22_sup(1) = -(c2p*d22_mid(1) + c3*d22_inf(1))*c1_inv
        d22_sup(2) = -(c2m*d22_mid(2) + c3*d22_inf(2))*c1_inv
      end if
      X(:,l) = X(:,l) + (Ip(:)*d22_sup(1)-Im(:)*d22_sup(2))*GL%w(i)*pi*0.25d0
      d22_inf = d22_mid
      d22_mid = d22_sup
    end do
  end do

  call gl_finalize(GL)

  if (present(f))  call savetxt(f,linspace(1,eL(2),eL(2)),X(:,oL(1)),X(:,oL(2)))

end subroutine dCBdCE


subroutine dCBdCp(oL,eL,CE,X,f)
! * derivative of B-mode Cl in terms of the lensing potential Cl
!
  implicit none
! [inputs]  
!   oL --- multipole range of residual ClBB
!   eL --- multipole range of delensing
!   CE --- power spectrum of the E-modes
  integer, intent(in) :: oL(1:2), eL(1:2)
  double precision, intent(in), dimension(eL(1):eL(2)) :: CE
!
! (optional)
!   f --- output filename
  character(*), intent(in), optional :: f
!
! [outputs]
!   X --- derivative of ClBB in terms of Clpp
  double precision, intent(out) :: X(eL(2),oL(2))
!
! [internal]
  type(gauss_legendre_params) :: GL
  integer :: i, l
  double precision :: mu, al, c1_inv, c2p, c2m, c3
  double precision, dimension(2) :: d22_sup, d22_mid, d22_inf, ZE33, ZE31, ZE11
  double precision, dimension(eL(2)) :: Ip, Im
  double precision, dimension(eL(2),2) :: ZP
  double precision, dimension(eL(1):eL(2)) :: al0, alm, alp

  do l = eL(1), eL(2)
    al = dble(l)
    al0(l) = dsqrt(al*(al+1))
    alp(l) = dsqrt((al-2)*(al+3))
    alm(l) = dsqrt((al+2)*(al-1))
  end do

  X = 0d0
  !* Gauss-Legendre Integration
  call gl_initialize(GL,int((3*eL(2)+1)/2),1d-15)

  do i = 1, GL%n
    mu = GL%z(i)
    call ZETA(3,3,eL,CE*alp**2,mu,ZE33)
    call ZETA(3,1,eL,CE*alp*alm,mu,ZE31)
    call ZETA(1,1,eL,CE*alm**2,mu,ZE11)
    call ZETA_L(1,1,eL,al0**2,mu,ZP)
    Ip = (ZE33(1)*ZP(:,1)+2*ZE31(1)*ZP(:,2)+ZE11(1)*ZP(:,1))
    Im = (ZE33(2)*ZP(:,2)+2*ZE31(2)*ZP(:,1)+ZE11(2)*ZP(:,2))
    d22_sup = 0d0 ;  d22_mid = 0d0 ;  d22_inf = 0d0
    do l = 2, oL(2) !recursion of wigner d function
      if(l==2) then 
        d22_sup = wigd_ini(2,2,mu)
      else
        al = dble(l)
        c1_inv = (al*(2d0*al-1d0))/(al**2-4d0)
        c2p = 4d0/((al-1d0)*al) - mu
        c2m = -4d0/((al-1d0)*al) - mu
        c3 = ((al-1d0)**2-4d0)/((al-1d0)*(2d0*al-1d0))
        d22_sup(1) = -(c2p*d22_mid(1) + c3*d22_inf(1))*c1_inv
        d22_sup(2) = -(c2m*d22_mid(2) + c3*d22_inf(2))*c1_inv
      end if
      X(:,l) = X(:,l) + (Ip(:)*d22_sup(1)-Im(:)*d22_sup(2))*GL%w(i)*pi*0.25d0
      d22_inf = d22_mid
      d22_mid = d22_sup
    end do
  end do

  call gl_finalize(GL)

  if (present(f))  call savetxt(f,linspace(1,eL(2),eL(2)),X(:,oL(1)),X(:,oL(2)))

end subroutine dCBdCp


end module nldd_delens

