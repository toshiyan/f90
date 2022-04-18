!////////////////////////////////////////////////////!
! * Statistical Analysis Subroutines
!////////////////////////////////////////////////////!

module statistics
  use myconst, only: signif
  use myutils, only: savetxt, linspace
  use array, only: sort_1d, sort_2d
  implicit none

  private signif
  private savetxt, linspace
  private sort_1d, sort_2d

contains


subroutine cl_likelihood(el,oC,tC,LL,f)
! Calculate likelihood (LL) from observed (oC) and theoretical Cls (tCl)
! Log likelihood of Cls: 
!   log(LL) = -(1/2) * {(2l+1)*[oC/tC + log(tC)]-(2l-1)*log(oC)}
  implicit none
  !I/O
  character(*), intent(in), optional :: f
  integer, intent(in) :: el(1:2)
  double precision, intent(in) :: oC(:), tC(:,:)
  double precision, intent(out) :: LL(:)
  !internal
  integer :: i, l
  double precision :: o, t, logL, iniL, maxL

  !* calculate Likelihood
  do i = 1, size(LL)
    logL = 0d0
    do l = el(1), el(2)
      o = oC(l)
      t = tC(i,l)
      logL = logL - 0.5d0*(dble(2*l+1)*(o/t+dlog(t))-dble(2*l-1)*dlog(o))
    end do
    if (i==1)  iniL = logL
    LL(i) = dexp(logL - iniL)
    if (logL - iniL<-10d0)  LL(i) = 0d0
    if (i==1)  maxL = LL(i)
    if (LL(i) > maxL)  maxL = LL(i) 
  end do

  !* normalized by maxL
  LL = LL/maxL
  write(*,*) 'maximum of likelihood:', maxL

  !* output likelihood
  if(present(f)) call savetxt(f,linspace(1,size(LL)),LL)

end subroutine cl_likelihood


subroutine normed_likelihood_2D(L,p1,p2,f)
  implicit none
  !I/O
  character(*), intent(in), optional :: f
  double precision, intent(in) :: p1(:), p2(:)
  double precision, intent(inout) :: L(:,:)
  !internal
  integer :: i, j

  L = L/maxval(L)
  write(*,*) 'maximum of likelihood:', maxval(L)

  if(present(f)) then
    open(unit=20,file=trim(f),status='replace')
    do i = 1, size(p1)
      do j = 1, size(p2)
        write(20,'(3(E12.5,1X))') p1(i), p2(j), L(i,j)
      end do
    end do
    close(20)
  end if

end subroutine normed_likelihood_2D


subroutine save_likelihood(p,L,dS)
  implicit none
  !I/O
  double precision, intent(in) :: p(:), dS(:)
  double precision, intent(inout) :: L(:)
  !internal
  integer :: n, pn
  integer, allocatable :: ii(:)
  double precision :: LL
  double precision, allocatable :: ps(:,:)

  pn = size(p)
  allocate(ii(pn),ps(3,pn));  ps=0d0

  L = L/sum(L*dS)
  call sort_1d(L,ii)
  
  LL = 0d0
  do n = 1, pn
    LL = LL + L(n)*dS(ii(n))
    if(LL<signif(1))  ps(1,n) = p(ii(n))
    if(LL>=signif(1).and.LL<signif(2)) ps(2,n) = p(ii(n))
    if(LL>=signif(2).and.LL<signif(3)) ps(3,n) = p(ii(n))
  end do
  call savetxt('sigma.dat',ps(1,:),ps(2,:),ps(3,:),L)

  deallocate(ii)

end subroutine save_likelihood


subroutine save_likelihood_2D(p1,p2,L,dS)
  implicit none
  double precision, intent(inout) :: L(:,:), p1(:), p2(:), dS(:,:)
  !recover on 1-dimension
  integer :: n, pn1, pn2, i, j
  integer, allocatable :: ii(:), jj(:)
  double precision :: LL

  pn1 = size(p1)
  pn2 = size(p2)
  allocate(ii(pn1),jj(pn2))

  L = L/sum(L*dS)
  call sort_2d(L,ii,jj)

  open(unit=21,file='1sigma.dat',status='replace')
  open(unit=22,file='2sigma.dat',status='replace')
  open(unit=23,file='3sigma.dat',status='replace')
  LL = 0.d0
  do n = 1, pn1*pn2
    i = ii(n)
    j = jj(n)
    LL = LL + L(i,j)*dS(i,j)
    if(LL>=signif(1).and.LL<signif(2)) write(21,'(3(E12.5,1X))') p1(i), p2(j), L(i,j)
    if(LL>=signif(2).and.LL<signif(3)) write(22,'(3(E12.5,1X))') p1(i), p2(j), L(i,j)
    if(LL>=signif(3)) write(23,'(3(E12.5,1X))') p1(i), p2(j), L(i,j)
  end do
  close(21)
  close(22)
  close(23)

  deallocate(ii,jj)

end subroutine save_likelihood_2D


end module statistics

