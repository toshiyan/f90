!////////////////////////////////////////////////////!
! * Kernel Functions of Fast Al Algorithm
!////////////////////////////////////////////////////!

module nldd_src
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


subroutine Kernels_src(rL,eL,WA,WB,X,kernel)
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
  double precision :: pm, mu, Ip, Im, al, c1_inv, c3
  double precision, dimension(2) :: d00_sup, d00_mid, d00_inf
  double precision, dimension(2) :: ZA00
  double precision, dimension(2) :: ZB00

  !* initialize
  if (kernel=='Sp'.or.kernel=='Gp') pm = 1d0
  if (kernel=='Sm'.or.kernel=='Gm') pm = -1d0

  selec = kernel
  if (kernel=='Sp'.or.kernel=='Sm') selec = 'Spm'
  if (kernel=='Gp'.or.kernel=='Gm') selec = 'Gpm'

  X = 0d0

  !* GL quadrature
  call gl_initialize(GL,int((3*max(rL(2),eL(2))+1)/2),1d-15)

  do i = 1, GL%n !Gauss-Legendre Integration
    mu = GL%z(i)
    select case(selec)
    case ('S0')
      call ZETA(0,0,rL,WA,mu,ZA00)
      call ZETA(0,0,rL,WB,mu,ZB00)
      Ip = ZA00(1)*ZB00(1)
    case ('G0') !same as S0
      call ZETA(0,0,rL,WA,mu,ZA00)
      call ZETA(0,0,rL,WB,mu,ZB00)
      Ip = ZA00(1)*ZB00(1)
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
      else
        c1_inv = (2d0*al-1d0)/al
        c3 = (al-1d0)**2/((al-1d0)*(2d0*al-1d0))
        d00_sup(1) = (mu*d00_mid(1) - c3*d00_inf(1))*c1_inv
        d00_sup(2) = (mu*d00_mid(2) - c3*d00_inf(2))*c1_inv
      end if
      X(l) = X(l) + Ip*d00_sup(1)*GL%w(i)*2d0*pi
      d00_inf = d00_mid
      d00_mid = d00_sup
    end do
  end do

  call gl_finalize(GL)

end subroutine Kernels_src


subroutine AlTT_src(rL,eL,Att,OCTT,ITT)
  implicit none
  !I/O
  integer, intent(in) :: rL(2), eL(2)
  double precision, intent(in), dimension(:), optional :: OCTT, ITT
  double precision, intent(out), dimension(:) :: Att
  !internal
  double precision, dimension(rL(1):rL(2)) :: W1, W2
  double precision, dimension(eL(2)) :: S0, G0

  write(*,*) 'AlTT (src)'

  if (present(OCTT)) W1 = 1d0 / OCTT(rL(1):rL(2))
  if (present(ITT))  W1 = ITT(rL(1):rL(2))
  if (.not.present(OCTT).and..not.present(ITT)) stop 'need obs TT or inverse of obs TT'
  S0 = 0d0
  call Kernels_src(rL,eL,W1,W1,S0,'S0')

  G0 = 0d0
  call Kernels_src(rL,eL,W1,W1,G0,'G0')

  Att(eL(1):eL(2)) = 1d0/(S0(eL(1):eL(2))+G0(eL(1):eL(2)))

end subroutine AlTT_src


end module nldd_src

