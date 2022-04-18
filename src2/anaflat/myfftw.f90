!////////////////////////////////////////////////////!
! * dft module
!////////////////////////////////////////////////////!

module myfftw
  use myconst, only: dlc, iu, pi, twopi
  !use array,   only: mapcut
  use anaflat, only: spin_weight, gaussian_alm
  implicit none

  INTEGER FFTW_ESTIMATE
  PARAMETER (FFTW_ESTIMATE=64)

  interface dft
    module procedure dft_1darray, dft_2darray
  end interface

  interface dft_pol
    module procedure dft_pol_1darray, dft_pol_2darray
  end interface

  interface dft_all
    module procedure dft_all_1darray, dft_all_2darray
  end interface

  interface derivemap
    module procedure derivemap_all, derivemap_nth_1d, derivemap_nth_2d, derivemap_alm_nth_1d
  end interface

  private FFTW_ESTIMATE
  private dlc, iu, pi, twopi
  !private mapcut
  private spin_weight, gaussian_alm

contains 

!////////////////////////////////////////////////////////////////////////!
!
! * FT convention of e.g., Hu & Okamoto (2002)
!
!   E(l) + iB(l) = - \int dx e^(-ilx) [ Q(x) + iU(x) ] e^(-2i\phi)
!   E(l) - iB(l) = - \int dx e^(-ilx) [ Q(x) - iU(x) ] e^(2i\phi)
! 
! Note: dx -> dx/2\pi in Lewis & Challinor (2006) 
!
! * Inversion
!
!   Q(x) + iU(x) = - \int (dl/2pi^2) e^(ilx) [ E(l) + iB(l) ] e^(2i\phi)
!   Q(x) - iU(x) = - \int (dl/2pi^2) e^(ilx) [ E(l) - iB(l) ] e^(-2i\phi)
!
! * Other quantities
!
!   Q(l) + iU(l) = - [ E(l) + iB(l) ] e^(2i\phi) 
!   Q(l) - iU(l) = - [ E(l) - iB(l) ] e^(-2i\phi) 
!
! where Q(l) and U(l) are the FT of Q(x) and U(x), respectively.
! This is equivalent to
! 
!   Q(l) = - E(l)*cos(2\phi) + B(l)*sin(2\phi)  
!   U(l) = - E(l)*sin(2\phi) - B(l)*cos(2\phi)  
!
! or
!
!   E(l) = - Q(l)*cos(2\phi) - U(l)*sin(2\phi)
!   B(l) = Q(l)*sin(2\phi) - U(l)*cos(2\phi)
!
! FT algorithm below assume 
!   
!   Q(x), U(x) <-> Q(l), U(l) <-> E(l), B(l)
!


!//////////////////////////////////////////////////////////////////////!
! * QU <-> EB transform in Fourier space

subroutine QU2EB(nn,D,QU,EB,trans)
! * trans = 1  :  QU -> EB
! * trans = -1 :  EB -> QU
  implicit none
  !I/O
  integer, intent(in) :: nn(2), trans
  double precision, intent(in) :: D(2)
  complex(dlc), intent(inout), dimension(:,:) :: QU,EB
  !internal
  integer :: i, j, n
  double precision :: x, y, sin2t, cos2t

  n = 1
  do i = 1, nn(1)
    x = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
    do j = 1, nn(2)
      y = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
      if (x==0d0.and.y==0d0) then
        cos2t = 0d0
        sin2t = 0d0
      else
        cos2t = (x**2-y**2)/(x**2+y**2)
        sin2t = 2d0*x*y/(x**2+y**2)
      end if
      if (trans==1) then
        EB(1,n) = QU(n,1)*cos2t + QU(n,2)*sin2t
        EB(2,n) = -QU(n,1)*sin2t + QU(n,2)*cos2t
      else if (trans==-1) then
        QU(n,1) = EB(1,n)*cos2t - EB(2,n)*sin2t
        QU(n,2) = EB(1,n)*sin2t + EB(2,n)*cos2t
      end if
      n = n + 1
    end do
  end do 

end subroutine QU2EB


!////////////////////////////////////////////////////////////////////////!
!
! DESCRIPTION OF dft ROUTINE
!
! [note]
!   trans = -1 -> alm to map
!   trans = +1 -> map to alm
!
! [from FFTW web (http://www.fftw.org/doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html) : Sec.4.8.1]
!   An FFTW_FORWARD transform corresponds to a sign of -1 in the exponent of the DFT. 
!   The standard "in-order" output ordering --- the k-th output corresponds to the frequency k/n (or k/T, where T is your total sampling period).
!
!
! * Definition of FT
!
!   FT:  F(lx,ly) = \int dxdy f(x,y)
!   IFT: f(x,y) = \int (dlx/2pi)(dly/2pi) F(lx,ly)
!
! * dx, dl ...
!
!   dx = D(1)/nn(1), dy = D(2)/nn(2)
!   dlx = 2pi/D(1), dly = 2pi/D(2)
!
! * Delta function 
!
!   Delta(lx,ly) = \int dxdy e^{-ixl} 
!       ->  dxdy*n(1)*n(2) = D(1)*D(2)        ( at lx=ly=0 )
!       ->  0                                 ( otherwise )
!
!   Delta(x,y)   = \int (dlx/2pi)*(dly/2pi) e^{ixl} 
!       ->  (dlx/2pi)(dly/2pi)*n(1)*n(2) 
!                   = n(1)*n(2)/D(1)/D(2)     ( at x=y=0 )
!       ->  0                                 ( otherwise )
!
! where the above Delta functions satisfy the usual condition
!
!   \int (dlx/2pi)(dly/2pi) Delta(lx,ly) = 1
!   \int dxdy Delta(x,y)     = 1
!
! The above should be satisfied because 
!
!   \int (dlx/2pi)(dly/2pi) Delta(lx,ly) 
!    = \int (dlx/2pi)(dly/2pi) \int dxdy e^{-ixl}
!    = \int dxdy Delta(x,y)
!

subroutine dft_2darray(map,nn,D,trans)
!* note (Jan 28, 2015)
! Even after FT and inverse FT, the resultant map has tiny error 
! but it biases the angular Cl significantly on large scale, 
! so these values are imposed to zero (add error control). 
  implicit none
  !I/O
  integer, intent(in) :: nn(2), trans
  double precision, intent(in) :: D(2)
  complex(dlc), intent(inout) :: map(:,:)
  !internal
  integer :: i, j, n, plan(8)
  double precision :: f, mean
  complex(dlc) :: amap(0:nn(1)-1,0:nn(2)-1)

  do i = 1, nn(1)
    do j = 1, nn(2)
      amap(i-1,j-1) = (-1)**(i+j)*map(i,j)
    end do
  end do

  call DFFTW_PLAN_dft_2D(plan,nn(1),nn(2),amap,amap,trans,FFTW_ESTIMATE)
  call DFFTW_EXEcutE_dft(plan,amap,amap)
  call DFFTW_DESTROY_PLAN(plan)

  if (trans==-1) f = 1d0/(D(1)*D(2))             ! multiply 1/dS * 1/(nn(1)*nn(2))
  if (trans==1 ) f = D(1)*D(2)/dble(nn(1)*nn(2)) ! area of one real-space pixel (dxdy)

  do i = 1, nn(1)
    do j = 1, nn(2)
      map(i,j) = f*amap(i-1,j-1)*(-1)**(i+j)
    end do
  end do

  ! error control
  mean = sum(abs(map))/dble(nn(1)*nn(2))
  do i = 1, nn(1)
    do j = 1, nn(2)
      if(abs(map(i,j))<mean*1d-15) map(i,j) = 0d0
    end do
  end do


end subroutine dft_2darray


subroutine dft_1darray(map,nn,D,trans)
!* revised (Jun 8, 2015)
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2), trans
  double precision, intent(in) :: D(1:2)
  complex(dlc), intent(inout) :: map(:)
  !internal
  integer :: i, j, n, plan(8)
  double precision :: f, mean
  complex(dlc), allocatable :: amap(:,:)

  allocate(amap(0:nn(1)-1,0:nn(2)-1))
  n = 1
  do i = 1, nn(1)
    do j = 1, nn(2)
      amap(i-1,j-1) = (-1)**(i+j)*map(n)
      n = n + 1
    end do
  end do

  call DFFTW_PLAN_dft_2D(plan,nn(1),nn(2),amap,amap,trans,FFTW_ESTIMATE)
  call DFFTW_EXEcutE_dft(plan,amap,amap)
  call DFFTW_DESTROY_PLAN(plan)

  if (trans==-1) f = 1d0/(D(1)*D(2))             ! multiply 1/dS * 1/(nn(1)*nn(2))
  if (trans==1 ) f = D(1)*D(2)/dble(nn(1)*nn(2)) ! area of one real-space pixel (dxdy)

  n = 1
  do i = 1, nn(1)
    do j = 1, nn(2)
      map(n) = f*amap(i-1,j-1)*(-1)**(i+j)
      n = n + 1
    end do
  end do
  deallocate(amap)

  ! error control (values smaller than double precision is assigned to zero)
  mean = sum(abs(map))/dble(nn(1)*nn(2))
  n = 1
  do i = 1, nn(1)
    do j = 1, nn(2)
      if(abs(map(n))<mean*1d-15) map(n) = 0d0
      n = n + 1
    end do
  end do

end subroutine dft_1darray


subroutine dft_pol_1darray(QU,nn,D,EB,trans)
  !trans=-1 : E,B -> Q+iU
  !trans=+1 : Q+iU -> E,B
  implicit none
  !I/O
  integer, intent(in) :: trans, nn(2)
  double precision, intent(in) :: D(2)
  double precision, intent(inout) :: QU(:,:)
  complex(dlc), intent(inout) :: EB(:,:)
  !internal
  complex(dlc), allocatable :: cP(:), P(:)

  select case(trans)
  case(-1)
    allocate(P(nn(1)*nn(2)))
    P = EB(1,:) + iu*EB(2,:)
    call spin_weight(P,nn,D,-1,2)
    call dft(P,nn,D,-1)
    QU(:,1) = dble(P)
    QU(:,2) = aimag(P)
    deallocate(P)
  case(1)
    allocate(P(nn(1)*nn(2)),cP(nn(1)*nn(2)))
    P = QU(:,1) + iu*QU(:,2)
    cP = conjg(P)
    call dft(P,nn,D,1)
    call spin_weight(P,nn,D,1,2)
    call dft(cP,nn,D,1)
    call spin_weight(cP,nn,D,1,-2)
    EB(1,:) = (P+cP)*0.5d0
    EB(2,:) = -iu*(P-cP)*0.5d0
    deallocate(P,cP)
  end select

end subroutine dft_pol_1darray


subroutine dft_pol_2darray(QU,nn,D,EB,trans)
  !trans=-1 : E,B -> Q+iU
  !trans=+1 : Q+iU -> E,B
  implicit none
  !I/O
  integer, intent(in) :: trans, nn(2)
  double precision, intent(in) :: D(2)
  double precision, intent(inout) :: QU(:,:,:)
  complex(dlc), intent(inout) :: EB(:,:,:)
  !internal
  complex(dlc), allocatable :: cP(:,:), P(:,:)

  select case(trans)
  case(-1)
    allocate(P(nn(1),nn(2)))
    P = EB(1,:,:) + iu*EB(2,:,:)
    call spin_weight(P,nn,D,-1,2)
    call dft(P,nn,D,-1)
    QU(:,:,1) = dble(P)
    QU(:,:,2) = aimag(P)
    deallocate(P)
  case(1)
    allocate(P(nn(1),nn(2)),cP(nn(1),nn(2)))
    P = QU(:,:,1) + iu*QU(:,:,2)
    cP = conjg(P)
    call dft(P,nn,D,1)
    call spin_weight(P,nn,D,1,2)
    call dft(cP,nn,D,1)
    call spin_weight(cP,nn,D,1,-2)
    EB(1,:,:) = (P+cP)*0.5d0
    EB(2,:,:) = -iu*(P-cP)*0.5d0
    deallocate(P,cP)
  end select

end subroutine dft_pol_2darray


subroutine dft_all_1darray(TQU,nn,D,TEB,trans)
  !trans=-1 : T,E,B -> T,Q+iU
  !trans=+1 : T,Q+iU -> T,E,B
  implicit none
  !I/O
  integer, intent(in) :: trans, nn(2)
  double precision, intent(in) :: D(2)
  double precision, intent(inout) :: TQU(:,:)
  complex(dlc), intent(inout) :: TEB(:,:)
  !internal
  integer :: npix
  complex(dlc), allocatable :: T(:), cP(:), P(:)

  npix = nn(1)*nn(2)

  select case(trans)
  case(-1)
    allocate(T(npix),P(npix))
    T = TEB(1,:)
    call dft(T,nn,D,-1)
    TQU(:,1) = dble(T)
    P = TEB(2,:) + iu*TEB(3,:)
    call spin_weight(P,nn,D,-1,2)
    call dft(P,nn,D,-1)
    TQU(:,2) = dble(P)
    TQU(:,3) = aimag(P)
    deallocate(T,P)
  case(1)
    allocate(T(npix),P(npix),cP(npix))
    T  = TQU(:,1)
    call dft(T,nn,D,-1)
    TQU(:,1) = dble(T)
    P  = TQU(:,2) + iu*TQU(:,3)
    cP = conjg(P)
    call dft(P,nn,D,1)
    call dft(cP,nn,D,1)
    call spin_weight(P,nn,D,1,2)
    call spin_weight(cP,nn,D,1,-2)
    TEB(1,:) = (P+cP)*0.5d0
    TEB(2,:) = -iu*(P-cP)*0.5d0
    deallocate(t,P,cP)
  end select

end subroutine dft_all_1darray


subroutine dft_all_2darray(TQU,nn,D,TEB,trans)
  !trans=-1 : T,E,B -> T,Q+iU
  !trans=+1 : T,Q+iU -> T,E,B
  implicit none
  !I/O
  integer, intent(in) :: trans, nn(2)
  double precision, intent(in) :: D(2)
  double precision, intent(inout) :: TQU(:,:,:)
  complex(dlc), intent(inout) :: TEB(:,:,:)
  !internal
  complex(dlc), allocatable :: T(:,:), cP(:,:), P(:,:)

  select case(trans)
  case(-1)
    allocate(T(nn(1),nn(2)),P(nn(1),nn(2)))
    T = TEB(1,:,:)
    call dft(T,nn,D,-1)
    TQU(:,:,1) = dble(T)
    P = TEB(2,:,:) + iu*TEB(3,:,:)
    call spin_weight(P,nn,D,-1,2)
    call dft(P,nn,D,-1)
    TQU(:,:,2) = dble(P)
    TQU(:,:,3) = aimag(P)
    deallocate(T,P)
  case(1)
    allocate(T(nn(1),nn(2)),P(nn(1),nn(2)),cP(nn(1),nn(2)))
    T  = TQU(:,:,1)
    call dft(T,nn,D,-1)
    TQU(:,:,1) = dble(T)
    P  = TQU(:,:,2) + iu*TQU(:,:,3)
    cP = conjg(P)
    call dft(P,nn,D,1)
    call dft(cP,nn,D,1)
    call spin_weight(P,nn,D,1,2)
    call spin_weight(cP,nn,D,1,-2)
    TEB(1,:,:) = (P+cP)*0.5d0
    TEB(2,:,:) = -iu*(P-cP)*0.5d0
    deallocate(t,P,cP)
  end select

end subroutine dft_all_2darray


subroutine pureEB(QU,nn,D,EB,W,Wd)
! * compute Smith's pure EB estimator
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  double precision, intent(in) :: D(1:2), W(:), QU(:,:)
  double precision, intent(in), optional :: Wd(:,:)
  complex(dlc), intent(out) :: EB(:,:)
  !internal
  integer :: i, j, n, npix
  double precision :: lx, ly
  double precision, allocatable :: WP(:,:)
  complex(dlc), allocatable, dimension(:) :: Wl,E1,E2x,E2y,B1,B2x,B2y
  complex(dlc), allocatable, dimension(:) :: Wx,Wy,Wxx,Wxy,Wyy
  complex(dlc), allocatable, dimension(:,:) :: WEB

  npix = nn(1)*nn(2)

  allocate(WP(npix,2),WEB(2,npix))
  WP(:,1) = QU(:,1)*W
  WP(:,2) = QU(:,2)*W
  call dft_pol(WP,nn,D,WEB,1)
  deallocate(WP)

  !//// derivatives of window functions ////!
  allocate(Wx(npix),Wy(npix),Wxx(npix),Wxy(npix),Wyy(npix))
  if (present(Wd)) then
    Wx  = Wd(1,:)
    Wy  = Wd(2,:)
    Wxx = Wd(3,:)
    Wyy = Wd(4,:)
    Wxy = Wd(5,:)
  else
    call array_deriv(nn,D,W,Wx,Wy,Wxx,Wxy,Wyy)
  end if 

  !//// correction terms ////!
  allocate(E1(npix),B1(npix),E2x(npix),E2y(npix),B2x(npix),B2y(npix))
  E1  = QU(:,1)*(Wyy-Wxx) - 2*QU(:,2)*Wxy
  B1  = QU(:,2)*(Wyy-Wxx) + 2*QU(:,1)*Wxy
  E2x =  2*iu*( QU(:,1)*Wx+QU(:,2)*Wy )
  E2y =  2*iu*( QU(:,2)*Wx-QU(:,1)*Wy )
  B2x =  2*iu*( QU(:,2)*Wx-QU(:,1)*Wy ) 
  B2y = -2*iu*( QU(:,1)*Wx+QU(:,2)*Wy ) 
  deallocate(Wx,Wy,Wxx,Wxy,Wyy)

  !* Transform to Fourier Space
  call dft(E1,nn,D,1)
  call dft(B1,nn,D,1)
  call dft(E2x,nn,D,1)
  call dft(E2y,nn,D,1)
  call dft(B2x,nn,D,1)
  call dft(B2y,nn,D,1)

  !* add corrections
  EB = 0d0
  n = 1
  do i = 1, nn(1)
    lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
    do j = 1, nn(2)
      ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
      if(lx==0.and.ly==0) cycle
      EB(1,n) = WEB(1,n) + (E1(n)+lx*E2x(n)+ly*E2y(n))/(lx**2+ly**2)
      EB(2,n) = WEB(2,n) + (B1(n)+lx*B2x(n)+ly*B2y(n))/(lx**2+ly**2)
      n = n + 1
    end do
  end do
  deallocate(E1,B1,E2x,E2y,B2x,B2y,WEB)

end subroutine pureEB


subroutine array_deriv(nn,D,W,Wx,Wy,Wxx,Wxy,Wyy)
  implicit none
  integer, intent(in) :: nn(2)
  double precision, intent(in) :: D(2), W(:)
  complex(dlc), intent(inout):: Wx(:), Wy(:), Wxx(:),Wxy(:), Wyy(:)
  integer :: i, j, n, npix
  double precision :: lx, ly
  complex(dlc), allocatable :: Wl(:)

  npix = nn(1)*nn(2)

  !* array in Fourier space
  allocate(Wl(npix))
  Wl = W
  call dft(Wl,nn,D,1)
  n = 1
  do i = 1, nn(1)
    lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
    do j = 1, nn(2)
      ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
      Wx(n) = iu*lx*Wl(n)
      Wy(n) = iu*ly*Wl(n)
      Wxx(n) = -lx**2*Wl(n)
      Wxy(n) = -lx*ly*Wl(n)
      Wyy(n) = -ly**2*Wl(n)
      n = n + 1
    end do
  end do 
  deallocate(Wl)

  !* Derivatives of windows functions
  call dft(Wx,nn,D,-1)
  call dft(Wy,nn,D,-1)
  call dft(Wxx,nn,D,-1)
  call dft(Wxy,nn,D,-1)
  call dft(Wyy,nn,D,-1)

end subroutine array_deriv


!///////////////////////////////////////////////////////////////////////!
!* generate Gaussian random fluctuations on 2D map

subroutine gaussian_map(nn,mm,D,eL,Cl,map,fix)
!simulate Gaussian alm with (mm), and mapping as (nn)-pixel map
  implicit none
  !I/O
  logical, intent(in), optional :: fix
  integer, intent(in) :: nn(1:2), mm(1:2), eL(1:2)
  double precision, intent(in) :: D(1:2), Cl(:)
  complex(dlc), intent(out) :: map(:)
  !internal
  integer :: iL(1:2)
  complex(dlc), dimension(:), allocatable :: alm

  allocate(alm(mm(1)*mm(2)))

  if(present(fix)) then
    call Gaussian_Alm(mm,D,eL,alm,Cl,fix)
  else
    write(*,*) 'generate gaussian alm'
    call Gaussian_Alm(mm,D,eL,alm,Cl)
  end if

  write(*,*) 'alm -> map'
  call dft(alm,mm,D,-1)
  if (mm(1)==nn(1).and.mm(2)==nn(2)) then
    map = alm
  else
    !write(*,*) 'cut map'
    stop 'mm and nn has to be equal'
    !call mapcut(alm,mm,map,nn)
  end if
  deallocate(alm)

end subroutine gaussian_map


subroutine gaussian_map_pol(nn,mm,D,eL,EE,BB,QU,P,fix)
!simulate Gaussian Elm and Blm with (mm), and mapping to (nn)-pixel QU map
  implicit none
  !I/O
  logical, intent(in), optional :: fix
  integer, intent(in) :: nn(1:2), mm(1:2), eL(1:2)
  double precision, intent(in) :: D(1:2)
  double precision, intent(in), optional :: EE(:), BB(:)
  double precision, intent(out), optional :: QU(:,:)
  complex(dlc), intent(out), optional :: P(:)
  !internal
  double precision, dimension(:,:), allocatable :: QUfull
  complex(dlc), dimension(:,:), allocatable :: alm

  allocate(alm(2,mm(1)*mm(2)))
  alm = 0d0

  if(present(fix)) then
    if(present(EE)) call gaussian_alm(mm,D,eL,alm(1,:),EE,fix)
    if(present(BB)) call gaussian_alm(mm,D,eL,alm(2,:),BB,fix)
  else
    if(present(EE)) call gaussian_alm(mm,D,eL,alm(1,:),EE)
    if(present(BB)) call gaussian_alm(mm,D,eL,alm(2,:),BB)
  end if

  allocate(QUfull(mm(1)*mm(2),2))
  call dft_pol_1darray(QUfull,mm,D,alm,-1)
  deallocate(alm)

  if (mm(1)==nn(1).and.mm(2)==nn(2)) then
    if(present(QU)) QU = QUfull
    if(present(P))  P  = QUfull(:,1) + iu*QUfull(:,2)
  else
    stop 'mm and nn has to be equal'
    !if(present(QU)) then
    !call mapcut(QUfull(:,1),mm,QU(:,1),nn)
    !  call mapcut(QUfull(:,2),mm,QU(:,2),nn)
    !else if(present(P)) then
    !  call mapcut(QUfull(:,1)+iu*QUfull(:,2),mm,P,nn)
    !end if
  end if
  deallocate(QUfull)

end subroutine gaussian_map_pol


!//// nth order derivatives ////!

subroutine derivemap_nth_1d(nn,D,map,dmap,nth)
! return nth order derivatives of input map
  implicit none
  !I/O
  integer, intent(in) :: nn(2), nth
  double precision, intent(in) :: D(2), map(:)
  double precision, intent(out) :: dmap(:,:)
  !internal
  integer :: i, j, k, n, npix
  double precision :: lx, ly
  complex(dlc), allocatable :: alm(:), dftmap(:)

  npix = nn(1)*nn(2)

  allocate(alm(npix))
  alm = cmplx(map)
  call dft(alm,nn,D,1)

  ! nth order (n>=1)
  do k = 0, nth
    n = 1
    allocate(dftmap(npix))
    do i = 1, nn(1)
      lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
      do j = 1, nn(2)
        ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
        dftmap(n) = (-iu)**nth*lx**(nth-k)*ly**k*alm(n) !need minus sign of iu for code consistency
        n = n + 1
      end do
    end do
    ! inverse-FT 
    call dft(dftmap,nn,D,-1)
    dmap(k+1,:) = dftmap
    deallocate(dftmap)
  end do

  deallocate(alm)

end subroutine derivemap_nth_1d


subroutine derivemap_nth_2d(D,map,derivmap,nth)
  implicit none
  !I/O
  integer, intent(in) :: nth
  double precision, intent(in) :: D(2), map(:,:)
  double precision, intent(out) :: derivmap(:,:,:)
  !internal
  integer :: i,j,k,nn(2)
  double precision :: lx,ly
  complex(dlc), allocatable :: ftmap(:,:), dftmap(:,:)

  ! get information
  nn(1) = size(map,dim=1)
  nn(2) = size(map,dim=2)

  ! prepare FT map
  allocate(ftmap(nn(1),nn(2))); ftmap=0d0
  ftmap = map
  call dft(ftmap,nn,D,1)

  ! nth order (n>=1)
  do k = 0, nth

    ! compute higher-order derivatives in F-space
    allocate(dftmap(nn(1),nn(2)))
    do i = 1, nn(1)
      lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
      do j = 1, nn(2)
        ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
        dftmap(i,j) = (-iu)**nth*lx**(nth-k)*ly**k*ftmap(i,j) !need minus sign of iu for code consistency
      end do
    end do

    ! inverse-FT 
    call dft(dftmap,nn,D,-1)
    derivmap(k+1,:,:) = dftmap
    deallocate(dftmap)

  end do

end subroutine derivemap_nth_2d


subroutine derivemap_alm_nth_1d(nn,D,alm,dmap,nth)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), nth
  double precision, intent(in) :: D(2)
  complex(dlc), intent(in) :: alm(:)
  double precision, intent(out) :: dmap(:,:)
  !internal
  integer :: i, j, k, n, npix
  double precision :: lx, ly
  complex(dlc), allocatable :: dftmap(:)

  npix = nn(1)*nn(2)

  ! nth order (n>=1)
  do k = 0, nth
    n = 1
    allocate(dftmap(npix))
    do i = 1, nn(1)
      lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
      do j = 1, nn(2)
        ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
        dftmap(n) = (-iu)**nth*lx**(nth-k)*ly**k*alm(n) !need minus sign of iu for code consistency
        n = n + 1
      end do
    end do
    ! inverse-FT 
    call dft(dftmap,nn,D,-1)
    dmap(k+1,:) = dftmap
    deallocate(dftmap)
  end do

end subroutine derivemap_alm_nth_1d


subroutine derivemap_all(D,map,derivmap)
! return derivatives up to nth order
  implicit none
  !I/O
  double precision, intent(in) :: D(2), map(:,:)
  double precision, intent(out) :: derivmap(:,:,:,:)
  !internal
  integer :: i,j,k,n,nth,nn(2)
  double precision :: lx,ly
  complex(dlc), allocatable :: ftmap(:,:), dftmap(:,:)

  ! get information
  nn(1) = size(map,dim=1)
  nn(2) = size(map,dim=2)
  nth   = size(derivmap,dim=1) - 1

  ! prepare FT map
  allocate(ftmap(nn(1),nn(2))); ftmap=0d0
  ftmap = map
  call dft(ftmap,nn,D,1)

  ! 0th order map
  derivmap(1,1,:,:) = map

  ! nth order (n>=1)
  do n = 1, nth
    do k = 0, n

      ! compute higher-order derivatives in F-space
      allocate(dftmap(nn(1),nn(2)))
      do i = 1, nn(1)
        lx = twopi*dble(i-1-nn(1)*0.5d0)/D(1)
        do j = 1, nn(2)
          ly = twopi*dble(j-1-nn(2)*0.5d0)/D(2)
          dftmap(i,j) = (-iu)**n*lx**(n-k)*ly**k*ftmap(i,j) !need minus sign of iu for code consistency
        end do
      end do

      ! inverse-FT 
      call dft(dftmap,nn,D,-1)
      derivmap(k+1,n+1,:,:) = dftmap
      deallocate(dftmap)

    end do
  end do

end subroutine derivemap_all


end module myfftw


