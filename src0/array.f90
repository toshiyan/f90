!////////////////////////////////////////////////////!
! * Statistical Analysis Subroutines
!////////////////////////////////////////////////////!

module array
  use myconst, only: dl, dlc, iu
  implicit none

  private dl, dlc, iu

  !* interfaces
  interface vec2mat
    module procedure vec2mat_dble, vec2mat_cmplx, vec2mat_diag
  end interface 

contains

!//// operator ////!

function vec2mat_dble(a,b) result(M)
  implicit none
  double precision, intent(in) :: a(:), b(:)
  double precision :: M(size(a),size(b))
  integer :: i,j

  do i = 1, size(a)
    do j = 1, size(b)
      M(i,j) = a(i)*b(j)
    end do
  end do

end function vec2mat_dble


function vec2mat_cmplx(a,b) result(M)
  implicit none
  complex(dlc), intent(in) :: a(:), b(:)
  complex(dlc) :: M(size(a),size(b))
  integer :: i,j

  do i = 1, size(a)
    do j =1, size(b)
      M(i,j) = a(i)*b(j)
    end do
  end do

end function vec2mat_cmplx


function vec2mat_diag(vec) result(M)
  implicit none
  double precision, intent(in) :: vec(:)
  integer :: i
  double precision :: M(size(vec),size(vec))

  M = 0d0
  do i = 1, size(vec)
    M(i,i) = vec(i)
  end do

end function vec2mat_diag


!//// 2D array as matrix ////!

!diagonal elements of n*n matrix
function mdiag(i,n)
  ! 1 + n + (n-1) + (n-2) + ...
  !   = n*(n+1)/2 - (n-i)(n-i+1)/2 - (n-i)
  implicit none
  !I/O
  integer, intent(in) :: i, n
  integer :: MDIAG

  MDIAG = n*(i-1) + i*(3-i)/2

end function mdiag


!elements of n*n matrix
!(1,1), (1,2), (1,3), ..., (2,2), (2,3), ... -> 1,2,3,...
function mlab(i,j,n)
  implicit none
  !I/O
  integer, intent(in) :: i,j,n
  integer :: MLAB

  if(j>=i) MLAB = mdiag(i,n) + (j-i)
  if(j<i)  MLAB = mdiag(j,n) + (i-j)

end function mlab


subroutine symmetric(M)
! symmetrize matrix (M_ji -> M_ij)
  implicit none
  !I/O
  double precision, intent(inout) :: M(:,:)
  !internal
  integer :: i, j, n

  n = size(M,dim=1)
  do i = 1, n
    do j = i + 1, n
      M(j,i) = M(i,j)
    end do
  end do

end subroutine symmetric


function mat_identity(N) result(M)
! get N*N indentity matrix
  implicit none
  integer, intent(in) :: N
  integer :: i
  double precision :: M(N,N)

  M = 0d0
  do i = 1, N
    M(i,i) = 1d0
  end do

end function mat_identity


subroutine Inv_Gauss_Jordan(K)
  implicit none
  double precision, intent(inout) :: K(:,:)
  integer :: IPIV(size(K,dim=1)), i, j, N
  double precision, allocatable, dimension(:) :: Lam 
  double precision, allocatable, dimension(:,:) ::L, M, U, MM

  N = size(K,dim=1)
  if(.not.size(K,dim=1)==size(K,dim=2)) stop 'size of matrix is incorrect'
  allocate(L(N,N),M(N,N),U(N,N),MM(N,N),LAM(N))

  do i = 1, N
    Lam(i) = K(i,i)
  end do

  do i = 1, N
    do j = 1, N
      K(i,j) = K(i,j)/(dsqrt(Lam(i)*Lam(j)))
    end do
  end do

  M = K

  call gaussj2(K,N,Mat_Identity(N),N) 

  do i = 1, N
    do j = 1, N
      K(i,j) = K(i,j)/(dsqrt(Lam(i)*Lam(j)))
    end do
  end do

  deallocate(L,M,U,MM,LAM)


end subroutine Inv_Gauss_Jordan


!+++++++++++++++++++++!
! Gauss-Jordan Method !
!+++++++++++++++++++++!

SUBROUTINE gaussj2(a,n,b,m)
  !implicit double precision (a-h,o-z)
  INTEGER::  m, mp, n, np
  REAL(dl) :: a(:,:), b(:,:) 
  ! Linear equation solution by Gauss-Jordan elimination, eq.(2.1.1) above. 
  ! a(1:n,1:n) is an input matrix stored in an array of physical dimensions
  ! np by np. 
  ! b(1:n,1:m) is an input matrix containing the m right-hand side 
  ! vectors, stored in an array of physical dimensions np by mp. 
  ! On output, a(1:n,1:n) is replaced by its matrix inverse, 
  ! and b(1:n,1:m) is replaced by the corresponding set of solution 
  ! vectors. Parameter: NMAX is the largest anticipated value of n. 
  integer, parameter :: NMAX = 50 
  INTEGER :: i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX) 
  ! The integer arrays ipiv, indxr,andindxc are used for bookkeeping 
  ! on the pivoting. 
  REAL(dl) :: big,dum,pivinv 

  np = size(a,dim=1)
  mp = size(b,dim=2)
  ipiv(1:n) = 0

  do i = 1, n  ! This is the main loop over the columns to be re-duced. 
    big = 0.d0 
    do j = 1, n  ! This is the outer loop of the search for a pivot element.
      if(ipiv(j).ne.1) then
        do k = 1, n 
          if (ipiv(k).eq.0) then 
            if (abs(a(j,k)).ge.big)then 
              big = abs(a(j,k)) 
              irow = j 
              icol = k 
            end if
          end if
        end do
      end if
    end do
    ipiv(icol) = ipiv(icol) + 1 
    if (irow.ne.icol) then 
      do l = 1, n 
        dum = a(irow,l) 
        a(irow,l) = a(icol,l) 
        a(icol,l) = dum 
      end do
      do l = 1, m 
        dum = b(irow,l) 
        b(irow,l) = b(icol,l) 
        b(icol,l) = dum 
      end do
    end if
    indxr(i) = irow 
    ! We are now ready to divide the pivot row by the pivot element, 
    ! located at irow and icol. 
    indxc(i) = icol 
    if (a(icol,icol).eq.0.) stop 'singular matrix in gaussj'  
    ! singular matrix in gaussj  
    pivinv = 1./a(icol,icol) 
    a(icol,icol) = 1. 
    a(icol,1:n) = a(icol,1:n)*pivinv 
    b(icol,1:m) = b(icol,1:m)*pivinv 
    do ll = 1, n  !Next, we reduce the rows... 
      if(ll.ne.icol)then  !...except for the pivot one, of course. 
        dum = a(ll,icol)
        a(ll,icol) = 0. 
        a(ll,1:n) = a(ll,1:n) - a(icol,1:n)*dum 
        b(ll,1:m) = b(ll,1:m) - b(icol,1:m)*dum 
      end if
    end do
  end do  !This is the end of the main loop over columns of the reduction. 
  do l = n, 1, -1 
  !It only remains to unscramble the solution in view of the 
  !column interchanges. We do this by in-terchanging pairs of 
  !columns in the reverse order that the permutation was built up. 
    if (indxr(l).ne.indxc(l))then 
      do k = 1, n 
        dum = a(k,indxr(l)) 
        a(k,indxr(l)) = a(k,indxc(l)) 
        a(k,indxc(l)) = dum 
      end do
    end if
  end do
  return 

END SUBROUTINE gaussj2


!//// Re-ordering array ////!

subroutine sort_1d(refval,ii)
!Sorting input array 
  implicit none
  !I/O
  integer, intent(out), optional :: ii(:) !label
  double precision, intent(inout) :: refval(:) !values
  !internal
  integer :: i, n, m, id, pn
  double precision :: dummy
  double precision, allocatable :: array(:)

  pn = size(refval)
  allocate(array(pn))

  !set array
  array = refval
  if (present(ii)) then
    do n = 1, pn
      ii(n) = n
    end do
  end if

  !sort
  do n = 1, pn
    do m = n + 1, pn
      if(array(n)>=array(m)) cycle
      dummy = array(n)
      array(n) = array(m)
      array(m) = dummy
      if (present(ii)) then
        id = ii(n)
        ii(n) = ii(m)
        ii(m) = id
      end if
    end do
  end do 

  refval = array

  deallocate(array)

end subroutine sort_1d


subroutine sort_2d(refval,ii,jj)
  implicit none
  !I/O
  integer, intent(inout) :: ii(:), jj(:)
  double precision, intent(inout) :: refval(:,:)
  !internal
  integer :: i, j, n, m, id, jd, pn1, pn2
  double precision :: dummy
  double precision, allocatable :: array(:)

  pn1 = size(refval,dim=1)
  pn2 = size(refval,dim=2)

  allocate(array(pn1*pn2))

  !set to array
  n = 1
  do i = 1, pn1
    do j = 1, pn2
      array(n) = refval(i,j) 
      ii(n) = i
      jj(n) = j
      n = n + 1
    end do
  end do

  !sort
  do n = 1, pn1*pn2
    do m = n + 1, pn1*pn2
      if(array(n)>=array(m)) cycle
      dummy = array(n)
      id = ii(n)
      jd = jj(n)
      array(n) = array(m)
      ii(n) = ii(m)
      jj(n) = jj(m)
      array(m) = dummy
      ii(m) = id
      jj(m) = jd
    end do
  end do

  do n=1, pn1*pn2
    refval(ii(n),jj(n)) = array(n)
  end do

  deallocate(array)

end subroutine sort_2d


end module array

