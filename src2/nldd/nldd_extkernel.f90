!////////////////////////////////////////////////////!
! * Extention of Kernel Functions for Fast Al Algorithm
!////////////////////////////////////////////////////!

module nldd_extkernel
  use myconst, only: pi
  use funcs,   only: WIGd_INI
  use myutils, only: gauss_legendre_params, gl_initialize, gl_finalize
  use nldd_kernel, only: Zeta
  implicit none

  !local parameters
  double precision, parameter :: cx(1:2) = [1d0,-1d0], fourpi = 4d0*pi
  private cx, fourpi

  private dl, pi
  private wigd_ini
  private gauss_legendre_params, gl_initialize, gl_finalize
  private Zeta

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
  double precision :: pm, mu, Ip, Im, I0, I1, al, c1_inv, c2p, c2m, c3
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

  I0 = 0d0
  I1 = 0d0

  do i = 1, GL%n !Gauss-Legendre Integration
    mu = GL%z(i)
    select case(selec)
    case ('S0')
      call ZETA(1,0,rL,WB*al0,mu,ZB10)
      call ZETA(1,1,rL,WB*al0**2,mu,ZB11)
    case ('G0')
      call ZETA(1,0,rL,WA*al0,mu,ZA10)
      call ZETA(1,0,rL,WB*al0,mu,ZB10)
    case ('Sc')
    case ('Gc')
    case ('Spm')
    case ('Gpm')
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
      select case (selec)
      case ('S0')
        I0 = ZA00(1)*ZB00(1)
        I1 = ZA00(1)*ZB10(1)
      case ('G0')
        I0 = ZA00(1)*ZB00(1)
        I1 = ZA00(1)*ZB10(1)
      case ('Sc')
      case ('Gc')
      case ('Spm')
      case ('Gpm')
      case default
        stop 'error: no kernel 2'
      end select
      X(1,l) = X(1,l) + 2*(I0*d00_sup(1))*GL%w(i)*pi
      X(2,l) = X(2,l) - I1*d10_sup(1)*dsqrt(dble(l*(l+1)))*GL%w(i)*pi
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


function arr_cut(arr,cut)
  implicit none
  double precision, intent(in) :: arr(:)
  integer, intent(in) :: cut(2)
  double precision :: arr_cut(cut(1):cut(2))

  arr_cut(cut(1):cut(2)) = arr(cut(1):cut(2))

end function arr_cut


subroutine AlTT(rL,eL,Agg,Acc,fCTT,OCTT,Amm,Agm)
  implicit none
  !I/O
  integer, intent(in) :: rL(2), eL(2)
  double precision, intent(in), dimension(:) :: fCTT, OCTT
  double precision, intent(out), dimension(:) :: Agg, Acc
  double precision, intent(out), dimension(:), optional :: Amm, Agm
  !internal
  integer :: l
  double precision, dimension(rL(1):rL(2)) :: W1, W2
  double precision, dimension(4,eL(2)) :: S0, G0

  write(*,*) 'TT'

  do l = rL(1), rL(2)
    W1(l) = 1d0/OCTT(l)
    W2(l) = W1(l) * fCTT(l)**2
  end do
  S0 = 0d0
  call Kernels(rL,eL,W1,W2,S0,'S0')

  do l = rL(1), rL(2)
    W2(l) = W1(l) * fCTT(l)
  end do 
  G0 = 0d0
  call Kernels(rL,eL,W2,W2,G0,'G0')

  do l = eL(1), eL(2)
    Agg(l) = 1d0/(S0(1,l)+G0(1,l))
    Acc(l) = 1d0/(S0(2,l)+G0(2,l))
    if(present(Amm)) Amm(l) = 1d0/(S0(3,l)+G0(3,l))
    if(present(Agm)) Agm(l) = 2d0/(S0(4,l)+G0(4,l))
  end do

end subroutine AlTT


end module nldd_extkernel


