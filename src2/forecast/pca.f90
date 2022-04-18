
module pca
  use readfile
  use myutils, only: str
  use mylapack95, only: inv_lapack, sygv_lapack

  implicit none

  private dl
  private str

contains


subroutine PCA_interface(el,wn,pF,Fl,pFl,n,nl,dz)
!* PCA analysis
  implicit none
  !I/O
  integer, intent(in) :: el(2), wn, n, nl
  double precision, intent(in) :: dz, pF(n,n), Fl(nl,n,n), pFl(nl,n,n)
  !internal
  integer :: i, j, l, r, i1, i2
  double precision, allocatable :: D(:), dw(:), Fpca(:,:), U(:,:), Fqr(:,:), iFrr(:,:), F(:,:)

  r = n - wn

  allocate(D(wn),Fpca(wn,wn),U(wn,wn),dw(wn),Fqr(wn,r),iFrr(r,r),F(n,n))
  F    = pF + sum(Fl(el(1):el(2),:,:)-pFl(el(1):el(2),:,:),dim=1)
  Fpca = F
  if (r>0) then
    Fqr  = F(r+1:n,1:r)
    iFrr = F(1:r,1:r)
    call inv_lapack(iFrr)
    Fpca = F(r+1:n,r+1:n) - matmul( matmul( Fqr, iFrr ), transpose(Fqr) )
  end if

  call sygv_lapack(Fpca,U,D)

  write(*,*) 'eigenvalues', D
  dw = dsqrt( matmul( U**2, D**(-1) ) )

  !check -------------------------------------------
  if (n == wn) then
    call inv_lapack(Fpca)
    do i = 1, wn
      if(abs(dw(i)/dsqrt(Fpca(i,i))-1.d0)>1e-3) then
        write(*,*) i, dw(i), dsqrt(Fpca(i,i))
        stop 'error in pca analysis'
      end if
    end do
  end if
  !-------------------------------------------------

  open(unit=20,file='pca.dat',status='replace')
  open(unit=21,file='risk_func.dat',status='replace')
  do i = 1, wn
    write(20,trim('(')//str(wn+3)//trim('(E12.5,1X))')) dz*dble(i-1), dw(i), D(i), U(i,:)
    write(21,*) i, risk_func(i,U)
  end do
  close(20)
  close(21)

  deallocate(F,Fpca,U,D,dw,Fqr,iFrr)

end subroutine PCA_interface


function risk_func(i,evec)
  implicit none
  !I/O
  integer, intent(in) :: i
  double precision, intent(in) :: evec(:,:)
  !internal
  integer :: wn
  double precision :: risk_func
  double precision, allocatable :: w_bias(:), w_fid(:), e(:,:)

  wn = size(evec,dim=1)
  allocate(w_bias(wn),w_fid(wn),e(i,wn))
  e = evec(1:i,1:wn)
  w_fid = -1d0
  w_bias = matmul( transpose(e), matmul(e, w_fid) )
  risk_func = sum( (w_bias-w_fid)**2 )
  deallocate(w_bias,w_fid,e)

end function risk_func


end module pca

