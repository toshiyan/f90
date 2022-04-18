!////////////////////////////////////////////////////!
! * Kernel Functions of Fast Al Algorithm
!////////////////////////////////////////////////////!

module nldd_rot
  use myconst, only: pi
  use funcs,   only: wigd_ini
  use myutils, only: gauss_legendre_params, gl_initialize, gl_finalize, linspace
  use nldd_kernel, only: zeta
  implicit none

  !local parameters
  double precision, parameter :: cx(1:2) = [1d0,-1d0], fourpi = 4d0*pi
  private cx, fourpi

  private pi
  private wigd_ini
  private gauss_legendre_params, gl_initialize, gl_finalize, linspace
  private zeta

contains


subroutine Kernels(rL,eL,WA,WB,X,kernel)
!Kernels of lensing reconstruction noise
  implicit none
  !I/O
  character(*), intent(in) :: kernel
  integer, intent(in) :: rL(1:2), eL(1:2)
  double precision, intent(in), dimension(rL(1):rL(2)) :: WA, WB
  double precision, intent(out) :: X(:)
  !internal
  type(gauss_legendre_params) :: GL
  character(LEN=4) :: selec
  integer :: i, l
  double precision :: mu, II, al, c1_inv, c2p, c2m, c3
  double precision, dimension(rL(1):rL(2)) :: al0, alp2, alm2
  double precision, dimension(2) :: d00_sup, d00_mid, d00_inf
  double precision, dimension(2) :: ZA22
  double precision, dimension(2) :: ZB22

  !* initialize
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
    select case(kernel)
    case ('Sm')
      call Zeta(2,2,rL,WA,mu,ZA22)
      call Zeta(2,2,rL,WB,mu,ZB22)
      II = ZA22(1)*ZB22(1)+ZA22(2)*ZB22(2)
    case ('Gm')
      call Zeta(2,2,rL,WA,mu,ZA22)
      call Zeta(2,2,rL,WB,mu,ZB22)
      II = ZA22(1)*ZB22(1)+ZA22(2)*ZB22(2)
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
      else
        c1_inv = (2d0*al-1d0)/al
        c3 = (al-1d0)**2/((al-1d0)*(2d0*al-1d0))
        d00_sup(1) = (mu*d00_mid(1) - c3*d00_inf(1))*c1_inv
        d00_sup(2) = (mu*d00_mid(2) - c3*d00_inf(2))*c1_inv
      end if
      X(l) = X(l) + 4d0*II*d00_sup(1)*GL%w(i)*pi
      d00_inf = d00_mid
      d00_mid = d00_sup
    end do
  end do

  call gl_finalize(GL)

end subroutine Kernels


subroutine AaaEB(rL,eL,Aaa,fCEE,OCEE,OCBB)
  implicit none
  !I/O
  integer, intent(in) :: rL(2), eL(2)
  double precision, intent(in), dimension(:) :: fCEE, OCEE, OCBB
  double precision, intent(out) :: Aaa(:)
  !internal
  double precision, dimension(rL(1):rL(2)) :: W1, W2
  double precision, dimension(eL(2)) :: Sm

  write(*,*) 'EB'

  W1 = 1d0/OCBB(rL(1):rL(2))
  W2 = fCEE(rL(1):rL(2))**2 / OCEE(rL(1):rL(2))
  Sm = 0d0
  call kernels(rL,eL,W1,W2,Sm,'Sm')

  Aaa(eL(1):eL(2)) = 1d0/Sm(eL(1):eL(2))

end subroutine AaaEB


end module nldd_rot

