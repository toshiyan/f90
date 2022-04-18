!//////////////////////////////////////////////
! * Lensing reconstruction tools
!//////////////////////////////////////////////

module mylensrec
  implicit none

  !* public variables
  integer, parameter :: QTT = 1, QTE = 2, QTB = 3, QEE = 4, QEB = 5, QBB = 6, QMV = 7
  integer, parameter :: QTTTE = 1, QTTEE = 2, QTEEE = 3, QTBEB = 4

  !* local parameters
  integer, parameter :: dlc = KIND(0d0)
  private dlc

contains


subroutine READGN(f,Al)
! Read Gaussian Normalization
! file should start from ell=2
  implicit none
  !I/O
  character(*), intent(in) :: f
  double precision, intent(out) :: Al(:,:)
  !internal
  integer :: i, l, m, rmax
  double precision, allocatable :: rNl(:)

  rmax = size(Al,dim=1)
  allocate(rNl(rmax))
  open(unit=20,file=trim(f),status="old")
  m = size(Al,dim=2) - 1
  do i = 1, m
    read(20,*) l, rNl(1:rmax)
    Al(1:rmax,l) = rNl(1:rmax)
  end do
  close(20)
  deallocate(rNl)

end subroutine READGN


subroutine DIAG_FILTER(els,alm,Cl,el)
  implicit none
  !I/O
  integer, intent(in) :: el(1:2)
  double precision, intent(in) :: els(:), Cl(:)
  complex(dlc), intent(inout) :: alm(:)
  !internal
  integer :: n

  do n = 1, size(alm)
    if(els(n)>=el(1).and.els(n)<el(2)) then
      alm(n) = alm(n)/Cl(int(els(n)))
    else
      alm(n) = 0d0
    end if
  end do

end subroutine DIAG_FILTER


subroutine DIAG_FILTER_MAP(els,alm,Cl,el)
!* added Jun 12, 2015
  implicit none
  !I/O
  integer, intent(in) :: el(1:2)
  double precision, intent(in) :: els(:), Cl(:)
  complex(dlc), intent(inout) :: alm(:)
  !internal
  integer :: n

  do n = 1, size(alm)
    if(els(n)>=el(1).and.els(n)<el(2)) then
      alm(n) = alm(n)/Cl(n)
    else
      alm(n) = 0d0
    end if
  end do

end subroutine DIAG_FILTER_MAP



subroutine NORMALIZED_ESTIMATOR(alm,Al,Nl,lmax)
  implicit none
  !I/O
  integer, intent(in) :: lmax
  double precision, intent(in) :: Al(:,:), Nl(:,:)
  complex(dlc), intent(inout) :: alm(size(Al,dim=1),0:lmax,0:lmax)
  !intenral
  integer :: l, q, mv

  mv = size(alm,dim=1)
  do l = 2, lmax
    do q = 1, max(mv-1,1)
      alm(q,l,:) = Al(q,l)*alm(q,l,:)
      if(mv>=3) then 
        alm(mv,l,:) = alm(mv,l,:)+Al(mv,l)*Nl(q,l)*alm(q,l,:)
      else
        alm(mv,l,:) = alm(1,l,:)
      end if
    end do 
  end do

end subroutine NORMALIZED_ESTIMATOR


subroutine NEST_Q_FLAT(est,AL,Nl)
! normalize estimators for each quadratic combination
  implicit none
  !I/O 
  double precision, intent(in) :: AL(:,:), Nl(:,:)
  complex(dlc), intent(inout) :: est(:,:)
  !internal 
  integer :: mv, e, q, qmax, npix, emax

  qmax = max(size(AL,dim=1)-1,1)
  npix = size(est,dim=2)

  do q = 1, qmax
    est(q,:) = est(q,:)*AL(q,:)
  end do

  ! minimum variance
  if(qmax>1) then
    mv = qmax + 1
    est(mv,:) = 0.d0
    do q = 1, qmax
      est(mv,:) = est(mv,:) + AL(mv,:)*Nl(q,:)*est(q,:)
    end do
  end if

end subroutine NEST_Q_FLAT


subroutine NEST_E_FLAT(est,AL)
! normalize estimators for each estimator type
  implicit none
  !I/O 
  double precision, intent(in) :: AL(:,:)
  complex(dlc), intent(inout) :: est(:,:)
  !internal 
  integer :: e, m, emax

  emax = size(est,dim=1)
  do e = 1, emax
    m = emax*(e-1) + e*(3-e)/2  ! diagonal component index
    est(e,:) = est(e,:)*AL(m,:)
  end do

end subroutine NEST_E_FLAT


subroutine NEST_ALL_FLAT(est,AL,NL)
! normalize estimators for each estimator type and quadratic combination
  implicit none
  !I/O 
  double precision, intent(in) :: AL(:,:,:), NL(:,:,:)
  complex(dlc), intent(inout) :: est(:,:,:)
  !internal 
  integer :: e, m, emax

  emax = size(est,dim=1)
  do e = 1, emax
    m = emax*(e-1) + e*(3-e)/2  ! diagonal component index
    call NEST_Q_FLAT(est(e,:,:),AL(m,:,:),NL(m,:,:))
  end do

end subroutine NEST_ALL_FLAT


function estcl_thvar(M,CAL) result(tV)
  implicit none
  !I/O
  integer, intent(in) :: M
  double precision, intent(in) :: CAL(:,:)
  double precision :: tV(size(CAL,dim=1),size(CAL,dim=2))
  !internal
  integer :: i, j, ei, ej, n

  n = 1
  do i = 1, M
    ei = M*(i-1) + i*(3-i)/2  ! diagonal component index
    tV(n,:) = CAL(ei,:)
    n = n + 1
    do j = i + 1, M
      ej = M*(j-1) + j*(3-j)/2  ! diagonal component index
      tV(n,:) = dsqrt(CAL(ei,:)*CAL(ej,:))
      n = n + 1
    end do
  end do

end function estcl_thvar


end module mylensrec

