module forecast_results
  use myconst, only: dl, ac2rad
  use mylapack95, only: inv_lapack, sygv_lapack
  use myutils, only: arr, str
  implicit none

  private dl, ac2rad
  private inv_lapack, sygv_lapack
  private arr, str

contains


subroutine OUTPUT_1SIGMA(F,e,LAB,fid)
  implicit none
  !I/O
  character(*), intent(in) :: LAB(:)
  integer, intent(in) :: e(1:2)
  real(dl), intent(in) :: F(:,:), fid(:)
  !internal
  integer  :: i, np
  real(dl) :: phi, a, b
  real(dl), allocatable :: InvF(:,:)

  np = size(F,dim=1)
  allocate(InvF(np,np))

  open(90,file='1sigma.dat',status='replace')

  write(90,'(A35)') '#//// marginalized parameters ////#'
  do i = 1, np
    write(90,'(A3,1X,E11.4)') LAB(i), fid(i)
  end do

  InvF = F
  call inv_lapack(InvF)
  write(90,'(A34)') '#//// expected 1 sigma error ////#'
  do i = 1, np
    write(90,'(A3,1X,E11.4)') LAB(i), dsqrt(InvF(i,i))
  end do

  !* write parameters for gnuplot
  if(np>1) then
    write(90,'(A36)') '#//// error ellipse parameters ////#'
    call calc2dellipse(InvF,e,phi,a,b)
    write(90,'(A3,E11.4)') 'p1=', phi
    write(90,'(A3,E11.4)') 'a1=', a
    write(90,'(A3,E11.4)') 'b1=', b
  end if
  close(90)

  deallocate(InvF)

end subroutine OUTPUT_1SIGMA


subroutine ERROR_INTERFACE(el,Fl)
  implicit none
  !I/O
  integer, intent(in) :: el(2)
  real(dl), intent(in), dimension(:,:,:) :: Fl
  !internal
  integer i, j, l, m, pn
  real(dl), dimension(:), allocatable :: error, unmargin, CC
  real(dl), dimension(:,:), allocatable :: F, InvF, Cor_F

  pn = size(Fl,dim=2)

  write(*,*) 'calculate one sigma error vs lmax'
  allocate(error(pn),unmargin(pn),CC(pn*(pn+1)/2),F(pn,pn),InvF(pn,pn),Cor_F(pn,pn))
  F = 0d0

  open(20,file='err_margin.dat',status='replace')
  open(21,file='err_unmargin.dat',status='replace')
  open(22,file='degeneracy.dat',status='replace')
  do l = el(1), el(2)
    F = F + Fl(l,:,:)
    InvF = F
    if(l<10) goto 11
    call inv_lapack(InvF)
    do i = 1, pn
      error(i) = dsqrt(InvF(i,i))
      unmargin(i) = 1/dsqrt(F(i,i))
    end do
    write(20,trim('(I4,1X,')//str(pn)//trim('(E12.5,1X))')) l, error(:)
    write(21,trim('(I4,1X,')//str(pn)//trim('(E12.5,1X))')) l, unmargin(:)
    CC = 0.d0
    m=0
    do i = 1, pn
      do j = i, pn
        m = m + 1
        CC(m) = InvF(i,j)/(dsqrt(InvF(i,i)*InvF(j,j)))
      end do
    end do
    write(22,trim('(I4,1X,')//str(m)//trim('(E12.5,1X))')) l, CC(:)
11  end do
  close(20)
  close(21)
  close(22)

  deallocate(error,unmargin,CC,F,InvF,Cor_F)

end subroutine Error_interface


subroutine Calc2DEllipse(F,id,phi,a,b)
  implicit none
  !I/O
  integer, intent(in) :: id(2)
  real(dl), intent(out) :: phi, a, b
  real(dl), dimension(:,:), intent(in) :: F
  !internal
  real(dl) :: G(2,2), U(2,2), D(2)

  G(1,1) = F(id(1),id(1))
  G(1,2) = F(id(1),id(2))
  G(2,1) = F(id(2),id(1))
  G(2,2) = F(id(2),id(2))    

  call sygv_lapack(G,U,D)
  phi = atan(U(2,1)/U(1,1))

  a = dsqrt(D(1))*1.516575089
  b = dsqrt(D(2))*1.516575089

end subroutine Calc2DEllipse


subroutine FoM(F,id,np,FoM2D)
  implicit none
  !I/O
  integer, intent(in) :: id(2), np
  real(dl), intent(in) :: F(np,np)
  real(dl), intent(out), optional :: FoM2D
  !intenral
  real(dl) :: K(2,2), U(2,2), D(2), InvF(np,np)

  InvF = F
  call inv_lapack(InvF)
  K = SubFisher(InvF,id)
  call sygv_lapack(K,U,D)

  write(*,*) "FoM =", dsqrt(D(1)*D(2))
  if(present(FoM2D)) FoM2D = dsqrt(D(1)*D(2))

end subroutine FoM


function SubFisher(M,id)
  implicit none
  !I/O
  integer, intent(in) :: id(2)
  real(dl), dimension(:,:), intent(in) :: M
  real(dl), dimension(2,2) :: SubFisher
  !internal
  real(dl), dimension(2,2) :: G

  G(1,1) = M(id(1),id(1))
  G(1,2) = M(id(1),id(2))
  G(2,1) = M(id(2),id(1))
  G(2,2) = M(id(2),id(2))
  call inv_lapack(G)
  SubFisher = G

end function SubFisher


subroutine Biased_Parameter(el,bp,Fl,pF,pFl)
  implicit none
  !I/O
  integer, intent(in) :: el(2)
  real(dl), intent(in) :: bp(:,:), Fl(:,:,:), pF(:,:), pFl(:,:,:)
  !internal
  integer :: l, n
  real(dl), dimension(:), allocatable :: bpl
  real(dl), dimension(:,:), allocatable :: bF, ibF

  n = size(bp,dim=2)
  open(unit=30,file='bias.dat',status='replace')
  open(unit=31,file='delta_l.dat',status='replace')
  allocate(bF(n,n),bpl(n),ibF(n,n))
  bpl = 0.d0
  bF = pF
  do l = el(1), el(2)
    !Summation 1 --- Fisher part -----------------
    bF = bF + Fl(l,:,:) - pFl(l,:,:)
    !Summation 2 --- Delta part ------------------
    bpl = bpl + bp(l,:)
    !( Fisher part ) times ( Delta part ) --------
    ibF = bF
    call inv_lapack(ibF)
    write(30,'(I6,15(2X,E11.4))') l, matmul(ibF,bpl)
    write(31,'(I6,15(2X,E11.4))') l, bpl
  end do
  close(30)
  close(31)

  deallocate(bpl,bF,ibF)

end subroutine Biased_Parameter


end module forecast_results

