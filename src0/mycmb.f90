!////////////////////////////////////////////////////!
! * Subroutines for CMB experimental specifications
!////////////////////////////////////////////////////!

module mycmbexp
  use readfile, only: read_int, read_prm, read_val
  use myutils, only: savetxt, linspace
  implicit none

  private read_int, read_prm, savetxt, linspace

  !* local variables
  integer, parameter :: dl = KIND(1d0)
  double precision, parameter :: pi = 3.1415926535897932384626433832795d0
  double precision, parameter :: ac2rad = pi/180d0/60d0
  double precision, parameter :: Tcmb = 2.726d6  ! micro K
  character(LEN=16) :: key(1:5)

  private dl, ac2rad, pi, Tcmb, key

  !* CMB experiment
  type CMBEXP
    integer :: nc
    double precision :: fsky
    double precision, dimension(8) :: fnu, fwhm, s_T, s_P
  end type CMBEXP

  !* Foreground parameters
  type FOREGROUND_PARAMETERS
    double precision :: nu0s=30d0,As=4.7d-5,el0s=350,alphas=-3,betas=-2.6
    double precision :: nu0d=94d0,Ad=1,el0d=10,alphad=2.2,betad=-2.5
    double precision :: p=0.05, fac=2.6662464d-12
    double precision :: fwhm_fid = 0.00872664625d0, s_P_fid = 0.06666666d0
  end type FOREGROUND_PARAMETERS


contains


subroutine set_expparams(expname,P)
! * default units:
!     fwhm    --- arcmin
!     s_T,s_P --- uK-arcmin
  implicit none
  character(*), intent(in) :: expname
  type(CMBEXP), intent(out) :: P

  select case (expname)

  case ("Planck")
    write(*,*) "set Planck params"
    P%nc = 7
    P%fsky = 0.65d0
    P%fwhm(1:7) = [33d0,23d0,14d0,9.5d0,7.1d0,5d0,5d0]
    P%s_T(1:7) = [4.4d0,6.5d0,9.8d0,6.8d0,6d0,13.1d0,40.1d0] !uk/pixel
    P%s_P(1:7) = [6.2d0,9.2d0,13.9d0,10.9d0,11.4d0,26.7d0,81.2d0] !uk/pixel

  case ("LiteBIRD")
    write(*,*) "set LiteBIRD params"
    P%nc = 1
    P%fsky = 1d0
    P%fwhm(1) = 30d0
    P%s_T(1) = 1d0
    P%s_P(1) = 2d0

  case ("ACTPol") !ACTPol Wide
    write(*,*) "set ACTPol params"
    P%nc = 1
    P%fsky = 0.1d0
    P%fwhm(1) = 1.4d0
    P%s_T(1) = 14d0
    P%s_P(1) = 20d0

  case ("SPT3G")
    P%nc = 3
    P%fwhm(1:3) = [1d0,1d0,1d0]
    P%s_T(1:3) = [4.243d0,2.475d0,4.243d0]
    P%s_P(1:3) = [6d0,3.5d0,6d0]

  case ("advACT")
    P%nc = 3
    P%fwhm(1:3) = [2.2d0,1.3d0,0.9d0]
    P%s_T(1:3) = [7.8d0,6.9d0,25d0]
    P%s_P(1:3) = [11.03d0,9.76d0,13.08d0]

  case ("CMBS4")
    P%nc = 1
    P%fwhm(1) = 3d0
    P%s_T(1) = 1d0
    P%s_P(1) = 1.4d0

  case DEFAULT
    stop "error: no experiment specified"

  end select

end subroutine set_expparams


subroutine read_expparams(key,P)
  implicit none
  !I/O
  character(*), intent(in) :: key(1:5)
  type(CMBEXP), intent(out) :: P

  P%nc = read_int(trim(key(1)))
  if(P%nc==0) stop "error: nchan is zero"
  if(read_val(trim(key(2)))) call read_prm(trim(key(2)),P%fnu(1:P%nc))
  call read_prm(trim(key(3)),P%fwhm(1:P%nc))
  if(read_val(trim(key(4)))) call read_prm(trim(key(4)),P%s_T(1:P%nc))
  if(read_val(trim(key(5)))) call read_prm(trim(key(5)),P%s_P(1:P%nc)) 

end subroutine read_expparams


subroutine INSTNL(el,NlT,NlP,expname,lknee,f,mykey)
  implicit none
  !I/O
  integer, intent(in) :: el(1:2)
  double precision, intent(out), optional :: NlT(:), NlP(:)
  character(*), intent(in), optional :: expname, f, mykey(:)
  integer, optional :: lknee
  !internal
  integer :: nc
  double precision :: Nls(1:2,1:el(2))
  double precision, dimension(:), allocatable :: theta, sigma
  type(CMBEXP) :: P

  if(present(expname)) then
    if(expname=="CV") then
      if(present(NlT)) NlT = 0d0
      if(present(NlP)) NlP = 0d0
      go to 11
    else
      call set_expparams(expname,P)
    end if
  else
    !read parameters from a file
    key = ["nchan","channel","fwhm_arcmin","sigma_T","sigma_P"]
    if(present(mykey))  key = mykey
    call read_expparams(key,P)
  end if

  nc = P%nc
  allocate(sigma(nc),theta(nc))
  theta = P%fwhm(1:nc) * ac2rad
  sigma = P%s_T(1:nc) * ac2rad / Tcmb
  if (present(NlT)) call CALCNL(NlT,el,sigma,theta)
  sigma = P%s_P(1:nc) * ac2rad / Tcmb
  if (present(NlP)) call CALCNL(NlP,el,sigma,theta)
  deallocate(sigma,theta)

11 Nls = 0d0
  if(present(NlT)) Nls(1,:) = NlT(:)
  if(present(NlP)) Nls(2,:) = NlP(:)
  if(present(f)) call savetxt(f,linspace(1,eL(2),eL(2)),Nls(1,:),Nls(2,:))

end subroutine INSTNl


subroutine Foreground(el,Fl,mykey)
  implicit none
  !I/O
  character(*), intent(in), optional :: mykey(:)
  integer, intent(in) :: el(2)
  double precision, intent(out) :: Fl(:)
  !internal
  character(LEN=16) :: key(1:5)
  integer :: l, n, m
  integer, parameter :: nc=6
  double precision :: al, dustfac, dnchan, ss, s_fg, nutps, nutpd
  double precision, allocatable, dimension(:) :: nu,Nlfid,Flnu,C,Nn
  double precision, allocatable, dimension(:,:) :: Nl
  double precision, dimension(:), allocatable :: freq, fwhm, s_T, s_P
  type(FOREGROUND_PARAMETERS) :: F
  type(CMBEXP) :: P

  key = ["nchan","channel","fwhm_arcmin","sigma_T","sigma_P"]
  if(present(mykey))  key = mykey
  allocate(fwhm(nc),s_P(nc),nu(nc))
  call read_expparams(key,P)

  allocate(C(2),Nn(2),Flnu(nc),Nl(el(2),nc),Nlfid(el(2)))

  s_fg = 0.1d0
  nutps = nu(1)
  nutpd = nu(nc)

  Fl = 0d0
  do l = el(1), el(2)
    al = dble(l)
    do n = 1, nc
      dnchan = dble(nc*(nc-1))/4d0
      !foreground spectrum
      C(1) = F%As*(nu(n)/F%nu0s)**(2*F%alphas)*(al/F%el0s)**(F%betas)
      C(2) = F%Ad*(nu(n)/F%nu0d)**(2*F%alphad)*(al/F%el0d)**(F%betad)
      dustfac = F%p*(dexp(F%nu0d*F%fac)-1d0)/(dexp(nu(n)*F%fac)-1d0)
      C(2) = C(2)*dustfac**2
      !noise spectrum
      ss = fwhm(n)**2/(8d0*dlog(2d0))
      Nl(l,n) = (s_P(n))**2*dexp(al*(al+1)*ss)
      Nn(1) = Nl(l,n)*(nu(n)/nutps)**(2*F%alphas)/dnchan
      Nn(2) = Nl(l,n)*(nu(n)/nutpd)**(2*F%alphad)/dnchan
      !total
      Flnu(n) = s_fg*sum(C)/Tcmb**2 + sum(Nn)
    end do
    ss = F%fwhm_fid**2/(8d0*dlog(2d0))
    Nlfid(l) = (F%fwhm_fid*F%s_P_fid/Tcmb)**2*dexp(al*(al+1)*ss)
    do n = 1, nc
      Fl(l) = Fl(l) + 1d0/(Nl(l,n)+Flnu(n))**2
      do m = n+1, nc
        Fl(l) = Fl(l) + 2d0/(Nl(l,n)+Flnu(n))/(Nl(l,m)+Flnu(m))
      end do
    end do
    Fl(l) = Fl(l)**(-0.5) - Nlfid(l)
  end do

  call savetxt("NlFl.dat",linspace(1,eL(2),eL(2)),Fl,Nlfid)

  deallocate(fwhm,s_P,Flnu,C,Nn,Nl,Nlfid)

end subroutine Foreground


subroutine CALCNL_TYPE(NlT,NlP,el,dtype)
  implicit none
  !I/O
  integer, intent(in) :: el(2)
  double precision, intent(out), optional :: NlT(:), NlP(:)
  type(CMBEXP), intent(in) :: dtype
  !internal
  integer :: l, n, nc
  double precision :: ss
  double precision, allocatable :: sigma(:), theta(:)

  nc = dtype%nc
  allocate(sigma(nc),theta(nc))
  theta = dtype%fwhm(1:nc) * ac2rad
  sigma = dtype%s_T(1:nc) * ac2rad / Tcmb
  if (present(NlT)) call CALCNL(NlT,el,sigma,theta)
  sigma = dtype%s_P(1:nc) * ac2rad / Tcmb
  if (present(NlP)) call CALCNL(NlP,el,sigma,theta)
  deallocate(sigma,theta)

end subroutine CALCNL_TYPE


subroutine CALCNL(Nl,el,sigma,theta,lknee)
!* computing noise power spectrum, assuming gaussian beam and white noise
  implicit none
  !I/O
  integer, intent(in) :: el(2)
  double precision, intent(in) :: sigma(:), theta(:)
  double precision, intent(out) :: Nl(:)
  integer, optional :: lknee
  !internal
  integer :: l, n, nc
  double precision :: ss

  nc = size(sigma)
  Nl = 0d0
  do l = el(1), el(2)
    do n = 1, nc
      ss = theta(n)**2/(8d0*dlog(2d0))
      Nl(l) = Nl(l) + dexp(dble(-l*(l+1))*ss)/sigma(n)**2
    end do
    Nl(l) = 1d0/Nl(l)
    if(present(lknee)) Nl(l) = Nl(l)*(1d0+(dble(lknee)/dble(l)))
  end do

end subroutine CALCNL


end module mycmbexp

