module forecast
  use readfile, only: read_prm, read_str, read_val
  use myutils,  only: str
  use array,    only: symmetric
  use mycls,    only: readcl_camb
  implicit none

  private str, read_prm, read_str, read_val, symmetric, readcl_camb

contains


subroutine MagCorrection(T,d,S,g,m,C)
  implicit none
  !I/O
  integer, intent(in) :: T, d, S(:), g(:)
  double precision, intent(in) :: m(:)
  double precision, intent(inout) :: C(:,:,:)
  !internal
  integer :: i, j, n

  n = size(m)
  do i = 1, n
    do j = 1, n
      ! ni*nj = (gi+mui)*(gj+muj) = gi*gj + gi*muj + mui*gj + mui*muj
      C(:,g(i),g(j)) = C(:,g(i),g(j)) + m(j)*C(:,g(i),S(j)) + m(i)*C(:,S(i),g(j)) + m(i)*m(j)*C(:,S(i),S(j))
      C(:,g(i),S(j)) = C(:,g(i),S(j)) + m(i)*C(:,S(i),S(j))
      C(:,S(i),g(j)) = C(:,S(i),g(j)) + m(j)*C(:,S(i),S(j))
    end do
    C(:,d,g(i)) = C(:,d,g(i)) + m(i)*C(:,d,S(i))
    C(:,g(i),d) = C(:,g(i),d) + m(i)*C(:,S(i),d)
    C(:,T,g(i)) = C(:,T,g(i)) + m(i)*C(:,T,S(i))
    C(:,g(i),T) = C(:,g(i),T) + m(i)*C(:,S(i),T)
  end do

end subroutine MagCorrection


subroutine cov_noise_cmb_l(Nl,Nov,T,E,B,d)
  implicit none
  !I/O
  integer, intent(in), optional :: T, E, B, d
  double precision, intent(in) :: Nl(:)
  double precision, intent(inout) :: Nov(:,:)
  !internal
  integer :: i

  if(present(T)) Nov(T,T) = Nl(1)
  if(present(E)) Nov(E,E) = Nl(2)
  if(present(B)) Nov(B,B) = Nl(2)
  if(present(d)) Nov(d,d) = Nl(3)

end subroutine cov_noise_cmb_l


subroutine cov_noise_gal_l(Ngal,eint,Nov,nz,S,g)
  implicit none
  !I/O
  integer, intent(in) :: nz
  integer, intent(in), optional :: S(:), g(:)
  double precision, intent(in) :: Ngal(:), eint
  double precision, intent(inout) :: Nov(:,:)
  !internal
  integer :: i, n

  n = 1
  do i = 1, nz
    if(present(S)) Nov(S(i),S(i)) = eint**2/Ngal(i)
    if(present(g)) Nov(g(i),g(i)) = 1d0/Ngal(i)
    n = n + nz - i + 1
  end do

end subroutine cov_noise_gal_l


!//// Fisher Matrix Subroutines ////!

function FISHER(iC,dCdp,l) result(FL)
!* compute Fisher matrix at one multipole using inverse covariance (iC) and derivative (dCdp)
!* Fl = (l+0.5)*Tr(C^-1 dC/dp C^-1 dCdq)
  implicit none
  !I/O
  integer, intent(in), optional :: l
  double precision, intent(in) :: iC(:,:),dCdp(:,:,:)
  !internal
  integer :: p, q, np, X, cn
  double precision, allocatable :: Fl(:,:), M(:,:)

  np = size(dCdp,dim=1)
  cn = size(iC,dim=2)

  allocate(M(cn,cn),Fl(np,np))

  !* Fisher matrix
  Fl = 0d0
  do p = 1, np
    do q = p, np
      ! matrix
      M = matmul(matmul(dCdp(p,:,:),iC),matmul(dCdp(q,:,:),iC))
      ! trace
      do X = 1, cn
        Fl(p,q) = Fl(p,q) + M(X,X)
      end do
    end do
  end do

  !* symmetric elements
  call symmetric(Fl)

  !* multiply the factor
  if(present(l)) FL = FL * (dble(l)+0.5d0)

  deallocate(M)

end function FISHER


function FISHER_SUM(el,iC,dCdp) result(F)
!* compute total Fisher matrix by summing up all multipoles
  implicit none
  !I/O
  integer, intent(in) :: el(2)
  double precision, intent(in) :: iC(:,:,:), dCdp(:,:,:,:)
  !internal
  integer :: l
  double precision :: F(size(dCdp,dim=1),size(dCdp,dim=1))

  F = 0d0
  do l = el(1), el(2)
    F = F + FISHER(iC(l,:,:),dCdp(:,l,:,:),l)
  end do

end function FISHER_SUM


function BIAS_PARAMS_L(l,iC,dCdp,dC) result(BP)
!* compute biased parameter at one multipole: 
!* Fl = (l+0.5)*Tr(C^-1 dC/dp C^-1 dC)
  implicit none
  !I/O
  integer, intent(in) :: l
  double precision, intent(in) :: iC(:,:),dCdp(:,:,:),dC(:,:)
  !internal
  integer :: p, np, X, cn
  double precision :: BP(size(dCdp,dim=1))
  double precision, allocatable :: M(:,:)

  np = size(dCdp,dim=1)
  cn = size(iC,dim=2)

  BP = 0d0
  !* Fisher matrix
  do p = 1, np
    allocate(M(cn,cn))
    M = matmul(matmul(dC,iC),matmul(dCdp(p,:,:),iC))
    do X = 1, cn
      BP(p) = BP(p) + M(X,X)
    end do
    deallocate(M)
  end do
  !* multiply the factor
  BP = BP * (dble(l)+0.5d0)

end function BIAS_PARAMS_L


!//// Derivatives ////!


!//// TO BE REMOVED ////!

subroutine DELTA_CL(dC,f1,f2,m)
!* difference of two Cls
  implicit none
  !I/O
  character(*), intent(in) :: f1, f2
  double precision, intent(out) :: dC(:,:,:)
  double precision, intent(in), optional :: m(:)
  !internal
  double precision, dimension(:,:,:), allocatable :: C1, C2
  integer :: l,n

  l = size(dC,dim=1)
  n = size(dC,dim=2)
  allocate(C1(l,n,n),C2(l,n,n))
  call CovCl(C1,f1)
  call CovCl(C2,f2)
  if(present(m)) call MagCorrection(O,m,C1)
  dC = C1 - C2
  deallocate(C1,C2)

end subroutine DELTA_CL


subroutine CALC_DERIVATIVES(q,deriv,fp,fm)
!* compute derivatives of the covariance (deriv) w.r.t. a parameter (q)
  implicit none
  !I/O
  character(*), intent(in) :: fp, fm
  integer, intent(in) :: q
  double precision, intent(out) :: deriv(:,:,:)
  !internal
  integer :: lmax, n
  double precision, allocatable :: m(:), dC(:,:,:)

  lmax = size(deriv,dim=1)
  n = size(deriv,dim=2)
  allocate(dC(lmax,n,n))
  if(p%LAB(q)==trim('s_z')) then
    !* with respect to slope
    allocate(m(p%sn))
    call DELTA_CL(dC,fp,fm,m)
    deallocate(m)
  else 
    !* other parameters
    call DELTA_CL(dC,fp,fm)
    step = p%der(q)
    if(q==p%O_L) step = (1d0-p%fid(p%O_L))*p%der(q)
    deriv = dC*0.5d0/step
  end if

  deallocate(dC)


end subroutine CALC_DERIVATIVES


subroutine addnoise(Cov,Nl,Ngal,eint)
  implicit none
  !I/O
  double precision, intent(in), optional :: Nl(:), Ngal(:), eint
  double precision, intent(inout) :: Cov(:,:)
  !internal
  integer :: n
  double precision, allocatable :: Nov(:,:)

  n = size(Cov,dim=1)
  allocate(Nov(n,n))
  Nov = 0d0
  if(present(Nl)) call cov_noise_cmb_l(O,Nl,Nov)
  if(present(Ngal)) call cov_noise_gal_l(O,Ngal,eint,Nov)
  Cov = Cov + Nov
  deallocate(Nov)

end subroutine addnoise


subroutine COVCL(Cov,f,nz,cln,T,E,B,d,g,S,TT,EE,TE,BB,dd,Td,TS,Tg,dg,dS,gg,SS,gS)
!* setting covariance matrix (Cov) using CAMB Cls
  implicit none
  !I/O
  integer, intent(in), optional :: TT,EE,TE,BB,dd,Td,TS(:),Tg(:),dg(:),dS(:),gg(:),SS(:),gS(:)
  integer, intent(in), optional :: T,E,B,d,g(:),S(:)
  integer, intent(in) :: nz, cln
  character(*), intent(in) :: f
  double precision, intent(out) :: Cov(:,:,:)
  !internal
  integer :: i, j, l, el(2), rn, n, m
  double precision, allocatable :: Cl(:,:), rCov(:,:,:)

  el(1) = lbound(Cov,dim=1)
  el(2) = ubound(Cov,dim=1)
  rn = size(Cov,dim=2)

  allocate(Cl(cln,el(2)),rCov(el(2),rn,rn))
  rCov = 0d0
  call READCL_CAMB(Cl,f,el)

  !//// Cl -> Covariance matrix ////!

  !* CMB diagonal
  if(present(TT)) rCov(:,T,T) = Cl(TT,:)
  if(present(EE)) rCov(:,E,E) = Cl(EE,:)
  if(present(BB)) rCov(:,B,B) = Cl(BB,:)
  if(present(dd)) rCov(:,d,d) = Cl(dd,:)

  !* CMB off-diagonal
  if(present(TE)) rCov(:,T,E) = Cl(TE,:)
  if(present(Td)) rCov(:,T,d) = Cl(Td,:)

  !* Galaxy
  n = 0
  m = 0
  do i = 1, nz
    if(present(TS)) rCov(:,T,S(i)) = Cl(TS(i),:)
    if(present(Tg)) rCov(:,T,g(i)) = Cl(Tg(i),:)
    if(present(dS)) rCov(:,d,S(i)) = Cl(dS(i),:)
    if(present(dg(i)) rCov(:,d,g(i)) = Cl(dg(i),:)
    do j = i, nz
      n = n + 1
      if(present(SS)) rCov(:,S(i),S(j)) = Cl(SS(n),:)
      if(present(gg)) rCov(:,g(i),g(j)) = Cl(gg(n),:)
    end do
    do j = 1, nz
      m = m + 1
      if(present(Sg)) rCov(:,S(i),g(j))=Cl(Sg(n),:)
    end do
  end do

  do i = el(1), el(2)
    call Symmetric(rCov(i,:,:))
  end do
  Cov = rCov
  
  deallocate(Cl,rCov)

end subroutine COVCL


end module forecast

