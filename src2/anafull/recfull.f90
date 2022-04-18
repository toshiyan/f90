!////////////////////////////////////////////////////!
! * Lensing Reconstruction in Fullsky
!////////////////////////////////////////////////////!

module recfull
  use alm_tools, only: alm2map, alm2map_spin, map2alm_spin
  use myconst, only: dl, dlc, iu

  private alm2map, alm2map_spin, map2alm_spin
  private dl, dlc, iu

contains 


!* Tlm,Elm,Blm: inverse-variance filtered alms (Cov^-1 (Tlm,Elm,Blm)^t)

subroutine quadtt(nsidet,Tlm1,Tlm2,fC,eL,rL,glm,clm)
!* fC = ClTT
  implicit none
  !I/O
  integer :: eL(2), rL(2), nsidet
  real(dl), intent(in) :: fC(:)
  complex(dlc), intent(in), dimension(0:rL(2),0:rL(2)) :: Tlm1, Tlm2
  complex(dlc), intent(inout), dimension(0:eL(2),0:eL(2)) :: glm, clm
  !internal
  integer :: l, npixt
  real(dl), allocatable :: at(:), map(:,:)
  complex(dlc), allocatable :: alm1(:,:,:), blm(:,:,:)

  write(*,*) 'calc TT-estimator'
  npixt = 12*nsidet**2

  !* convolution
  allocate(alm1(1,0:rL(2),0:rL(2)))
  alm1 = 0d0
  do l = rL(1), rL(2)
    alm1(1,l,0:l) = Tlm1(l,0:l)
  end do 
  allocate(at(0:npixt-1))
  call alm2map(nsidet,rL(2),rL(2),alm1,at)
  deallocate(alm1)

  allocate(alm1(2,0:rL(2),0:rL(2)))
  alm1 = 0d0
  do l = rL(1), rL(2)
    alm1(1,l,0:l) = fC(l)*Tlm2(l,0:l)*dsqrt(dble((l+1)*l))
  end do 
  allocate(map(0:npixt-1,2))
  call alm2map_spin(nsidet,rL(2),rL(2),1,alm1,map)
  map(:,1) = at*map(:,1)
  map(:,2) = at*map(:,2)
  deallocate(at,alm1)

  allocate(blm(2,0:eL(2),0:eL(2)))
  call map2alm_spin(nsidet,eL(2),eL(2),1,map,blm)
  deallocate(map)

  !* compute glm and clm
  write(*,*) 'compute estimator'; glm = 0d0; clm=0d0
  do l = 0, eL(2)
    glm(l,0:l) = dsqrt(dble(l*(l+1)))*blm(1,l,0:l)
    clm(l,0:l) = dsqrt(dble(l*(l+1)))*blm(2,l,0:l)
  end do

end subroutine quadtt


subroutine QUADTE(nsidet,Tlm,Elm,glm,clm,fC,eL,tL)
!* fC = ClTE
  implicit none 
  !I/O
  integer, intent(in) :: eL(2), tL(2), nsidet
  real(dl), intent(in) :: fC(:)
  complex(dlc), intent(in), dimension(0:tL(2),0:tL(2)) :: Tlm, Elm
  complex(dlc), intent(out), dimension(0:eL(2),0:eL(2)) :: glm, clm
  !internal
  integer :: l, npixt
  real(dl), dimension(:), allocatable :: AT
  real(dl), dimension(:,:), allocatable :: A, A1, A3, map
  complex(dlc), dimension(:,:,:), allocatable :: alm1, alm3, blm

  write(*,*) 'calc TE-estimator'
  npixt = 12*nsidet**2

  write(*,*) 'ET'

  !* convolution
  allocate(alm1(2,0:tL(2),0:tL(2)))
  do l = tL(2), tL(2)
    alm1(1,l,:) = Elm(l,:)
  end do 
  allocate(A(0:npixt-1,2))
  call alm2map_spin(nsidet,tL(2),tL(2),2,alm1,A)
  deallocate(alm1)

  allocate(alm1(2,0:tL(2),0:tL(2)),alm3(2,0:tL(2),0:tL(2)))
  alm1 = 0d0;  alm3 = 0d0
  do l = tL(2), tL(2)
    alm1(1,l,:) = fC(l)*Tlm(l,:)*dsqrt(dble((l+2)*(l-1)))
    alm3(1,l,:) = fC(l)*Tlm(l,:)*dsqrt(dble((l-2)*(l+3)))
  end do 
  allocate(A1(0:npixt-1,2))
  call alm2map_spin(nsidet,tL(2),tL(2),1,alm1,A1)
  deallocate(alm1)
  allocate(A3(0:npixt-1,2))
  call alm2map_spin(nsidet,tL(2),tL(2),3,alm3,A3)
  deallocate(alm3)

  allocate(map(0:npixt-1,2))
  !map = A*conjg(A1) - conjg(A)*A3
  map(:,1) = A(:,1)*(A1(:,1)-A3(:,1)) + A(:,2)*(A(:,2)-A3(:,2))
  map(:,2) = -A(:,1)*(A1(:,2)+A3(:,2)) + A(:,2)*(A(:,1)+A3(:,1))
  deallocate(A,A1,A3)

  allocate(blm(2,0:eL(2),0:eL(2)))
  call map2alm_spin(nsidet,eL(2),eL(2),1,map,blm)
  deallocate(map)

  !* compute glm and clm
  write(*,*) 'compute grad and curl'
  do l = el(1), eL(2)
    glm(l,:) = 0.5d0*dsqrt(dble(l*(l+1)))*blm(1,l,:)
    clm(l,:) = 0.5d0*dsqrt(dble(l*(l+1)))*blm(2,l,:)
  end do
  deallocate(blm)

  !//// second part ////!
  write(*,*) 'TE'
  !* convolution
  allocate(alm1(1,0:tL(2),0:tL(2)))
  alm1 = 0d0
  do l = tL(1), tL(2)
    alm1(1,l,:) = Tlm(l,:)
  end do
  allocate(AT(0:npixt-1))
  call alm2map(nsidet,tL(2),tL(2),alm1,AT)
  deallocate(alm1)

  allocate(alm1(2,0:tL(2),0:tL(2)))
  alm1 = 0d0
  do l = tL(1), tL(2)
    alm1(1,l,:) = Elm(l,:)*fC(l)*dsqrt(dble(l*(l+1)))
  end do
  allocate(map(0:npixt-1,2))
  call alm2map_spin(nsidet,tL(2),tL(2),1,alm1,map)
  map(:,1) = AT*map(:,1)
  map(:,2) = AT*map(:,2)
  deallocate(AT,alm1)

  allocate(blm(2,0:eL(2),0:eL(2)))
  call map2alm_spin(nsidet,eL(2),eL(2),1,map,blm)
  deallocate(map)

  !* compute glm and clm
  write(*,*) 'compute grad and curl'
  do l = el(1), eL(2)
    glm(l,:) = glm(l,:) + dsqrt(dble(l*(l+1)))*blm(1,l,:)
    clm(l,:) = clm(l,:) + dsqrt(dble(l*(l+1)))*blm(2,l,:)
  end do

end subroutine QUADTE


subroutine QuadTB(nsidet,Tlm,Blm,glm,clm,fC,eL,tL)
!* fC = ClTE
  implicit none 
  !I/O
  integer, intent(in) :: eL(2), tL(2), nsidet
  real(dl), intent(in) :: fC(:)
  complex(dlc), intent(in), dimension(0:tL(2),0:tL(2)) :: Tlm, Blm
  complex(dlc), intent(out), dimension(0:eL(2),0:eL(2)) :: glm, clm
  !internal
  integer :: l, npixt
  real(dl), dimension(:,:), allocatable :: A, A1, A3, map
  complex(dlc), dimension(:,:,:), allocatable :: alm1, alm3, zlm

  write(*,*) 'calc TB-estimator'
  npixt = 12*nsidet**2

  !* convolution
  allocate(alm1(2,0:tL(2),0:tL(2)))
  alm1 = 0d0
  do l = tL(1), tL(2)
    alm1(2,l,:) = Blm(l,:)
  end do 
  allocate(A(0:npixt-1,2))
  call alm2map_spin(nsidet,tL(2),tL(2),2,alm1,A)
  deallocate(alm1)

  allocate(alm1(2,0:tL(2),0:tL(2)),alm3(2,0:tL(2),0:tL(2)))
  alm1 = 0d0;  alm3 = 0d0
  do l = tL(1), tL(2)
    alm1(1,l,:) = fC(l)*Tlm(l,:)*dsqrt(dble((l+2)*(l-1)))
    alm3(1,l,:) = fC(l)*Tlm(l,:)*dsqrt(dble((l-2)*(l+3)))
  end do 
  allocate(A1(0:npixt-1,2))
  call alm2map_spin(nsidet,tL(2),tL(2),1,alm1,A1)
  deallocate(alm1)
  allocate(A3(0:npixt-1,2))
  call alm2map_spin(nsidet,tL(2),tL(2),3,alm3,A3)
  deallocate(alm3)

  allocate(map(0:npixt-1,2))
  !map = A*conjg(A1) - conjg(A)*A3
  map(:,1) = A(:,1)*(A1(:,1)-A3(:,1)) + A(:,2)*(A1(:,2)-A3(:,2))
  map(:,2) = -A(:,1)*(A1(:,2)+A3(:,2)) + A(:,2)*(A1(:,1)+A3(:,1))
  deallocate(A,A1,A3)

  allocate(zlm(2,0:eL(2),0:eL(2)))
  call map2alm_spin(nsidet,eL(2),eL(2),1,map,zlm)
  deallocate(map)

  !* compute glm and clm
  write(*,*) 'compute grad and curl'
  do l = el(1), eL(2)
    glm(l,:) = 0.5d0*dsqrt(dble(l*(l+1)))*zlm(1,l,:)
    clm(l,:) = 0.5d0*dsqrt(dble(l*(l+1)))*zlm(2,l,:)
  end do
  deallocate(zlm)

end subroutine QuadTB


subroutine QUADEE(nsidet,Elm1,Elm2,glm,clm,fC,eL,tL)
!* fC = ClEE
  implicit none 
  !I/O
  integer, intent(in) :: nsidet, eL(2), tL(2)
  real(dl), intent(in) :: fC(:)
  complex(dlc), intent(in), dimension(0:tL(2),0:tL(2)) :: Elm1, Elm2
  complex(dlc), intent(out), dimension(0:eL(2),0:eL(2)) :: glm, clm
  !internal
  integer :: l, npixt
  real(dl), dimension(:,:), allocatable :: A, A1, A3, map
  complex(dlc), dimension(:,:,:), allocatable :: alm, blm

  write(*,*) 'calc EE-estimator'
  npixt = 12*nsidet**2

  !* convolution
  write(*,*) 'convolution'
  allocate(A(0:npixt-1,2),alm(2,0:tL(2),0:tL(2)))
  alm = 0d0
  do l = tL(1), tL(2)
    alm(1,l,:) = Elm1(l,:)
  end do
  call alm2map_spin(nsidet,tL(2),tL(2),2,alm,A)
  deallocate(alm)

  allocate(alm(2,0:tL(2),0:tL(2)),blm(2,0:tL(2),0:tL(2)))
  alm = 0d0;  blm = 0d0
  do l = tL(1), tL(2)
    alm(1,l,:) = fC(l)*Elm2(l,:)*dsqrt(dble((l+2)*(l-1)))
    blm(1,l,:) = fC(l)*Elm2(l,:)*dsqrt(dble((l-1)*(l+3)))
  end do
  allocate(A1(0:npixt-1,2))
  call alm2map_spin(nsidet,tL(2),tL(2),1,alm,A1)
  deallocate(alm)
  allocate(A3(0:npixt-1,2))
  call alm2map_spin(nsidet,tL(2),tL(2),3,blm,A3)
  deallocate(blm)

  allocate(map(0:npixt-1,2))
  !map = A*conjg(A1) - conjg(A)*A3
  map(:,1) = A(:,1)*(A1(:,1)-A3(:,1)) + A(:,2)*(A1(:,2)-A3(:,2))
  map(:,2) = -A(:,1)*(A1(:,2)+A3(:,2)) + A(:,2)*(A1(:,1)+A3(:,1))
  deallocate(A,A1,A3)

  allocate(alm(2,0:eL(2),0:eL(2)))
  call map2alm_spin(nsidet,el(2),el(2),1,map,alm)
  deallocate(map)

  !* compute glm and clm 
  write(*,*) 'compute grad and curl'
  glm = 0d0
  clm = 0d0
  write(*,*) 'grad/curl EE'
  do l = eL(1), eL(2)
    glm(l,0:eL(2)) = 0.5d0*dsqrt(dble(l)*dble(l+1))*alm(1,l,:)
    clm(l,0:eL(2)) = 0.5d0*dsqrt(dble(l)*dble(l+1))*alm(2,l,:)
  end do
  deallocate(alm)
  write(*,*) 'finish EE'

end subroutine QuadEE


subroutine quadeb(nsidet,Elm,Blm,fC,eL,tL,glm,clm)
!* fC = ClEE
  implicit none
  !I/O
  integer, intent(in) :: nsidet, eL(2), tL(2)
  real(dl), intent(in) :: fC(:)
  complex(dlc), intent(in), dimension(0:tL(2),0:tL(2)) :: Elm, Blm
  complex(dlc), intent(out), dimension(0:eL(2),0:eL(2)) :: glm, clm
  !internal
  integer :: l, npixt
  real(dl), dimension(:,:), allocatable :: A,A1,A3,map
  complex(dlc), dimension(:,:,:), allocatable :: alm1,alm3,tlm

  write(*,*) 'calc EB-estimator'
  npixt = 12*nsidet**2

  !* convolution
  write(*,*) 'convolution'
  allocate(alm1(2,0:tL(2),0:tL(2)))
  alm1 = 0d0
  do l = tL(1), tL(2)
    alm1(2,l,:) = Blm(l,:)
  end do 
  allocate(A(0:npixt-1,2))
  call alm2map_spin(nsidet,tL(2),tL(2),2,alm1,A)
  deallocate(alm1)

  allocate(alm1(2,0:tL(2),0:tL(2)),alm3(2,0:tL(2),0:tL(2)))
  alm1 = 0d0;  alm3 = 0d0
  do l = tL(1), tL(2)
    alm1(1,l,:) = fC(l)*Elm(l,:)*dsqrt(dble((l+2)*(l-1)))
    alm3(1,l,:) = fC(l)*Elm(l,:)*dsqrt(dble((l-2)*(l+3)))
  end do 
  allocate(A1(0:npixt-1,2))
  call alm2map_spin(nsidet,tL(2),tL(2),1,alm1,A1)
  deallocate(alm1)

  allocate(A3(0:npixt-1,2))
  call alm2map_spin(nsidet,tL(2),tL(2),3,alm3,A3)
  deallocate(alm3)

  allocate(map(0:npixt-1,2))
  !map = A*conjg(A1) - conjg(A)*A3
  map(:,1) = A(:,1)*(A1(:,1)-A3(:,1)) + A(:,2)*(A1(:,2)-A3(:,2))
  map(:,2) = -A(:,1)*(A1(:,2)+A3(:,2)) + A(:,2)*(A1(:,1)+A3(:,1))
  deallocate(A,A1,A3)

  allocate(tlm(2,0:eL(2),0:eL(2)))
  call map2alm_spin(nsidet,eL(2),eL(2),1,map,tlm)
  deallocate(map)

  !* compute glm and clm 
  write(*,*) 'compute grad and curl'
  do l = eL(1), eL(2)
    glm(l,:) = 0.5d0*dsqrt(dble(l*(l+1)))*tlm(1,l,:)
    clm(l,:) = 0.5d0*dsqrt(dble(l*(l+1)))*tlm(2,l,:)
  end do
  deallocate(tlm)

!  allocate(alm1(0:tL(2),0:tL(2)),alm3(0:tL(2),0:tL(2)))
!  alm1 = 0d0;  alm3 = 0d0
!  do l = tL(1), tL(2)
!    alm1(l,:) = fCBB(l)*Blm(l,:)*llsq(l,-2)
!    alm3(l,:) = fCBB(l)*Blm(l,:)*llsq(l,2)
!  end do 

!  allocate(A(0:npixt-1))
!  call elm2map_spin(nsidet,tL(2),tL(2),2,Elm,A)

!  allocate(A1(0:npixt-1))
!  call blm2map_spin(nsidet,tL(2),tL(2),1,alm1,A1)
!  deallocate(alm1)

!  allocate(A3(0:npixt-1))
!  call blm2map_spin(nsidet,tL(2),tL(2),3,alm3,A3)
!  deallocate(alm3)

!  allocate(map(0:npixt-1))
!  map = A*conjg(A1) - conjg(AE)*A3
!  deallocate(A,A1,A3)

!  allocate(tlm(2,0:eL(2),0:eL(2)))
!  call cmplxmap2alm_spin(nsidet,eL(2),eL(2),1,map,tlm)
!  deallocate(map)

  !* compute glm and clm
!  write(*,*) 'compute grad and curl'
!  do l = eL(1), eL(2)
!    est(E%g,l,:) = est(E%g,l,:) + 0.5d0*llsq(l,0)*tlm(1,l,:)
!    est(E%c,l,:) = est(E%c,l,:) + 0.5d0*llsq(l,0)*tlm(2,l,:)
!  end do
!  deallocate(tlm)

end subroutine quadeb


subroutine QUADBB(nsidet,Blm1,Blm2,glm,clm,fC,eL,tL)
!* fC = ClBB
  implicit none
  !I/O
  integer, intent(in) :: nsidet, eL(2), tL(2)
  real(dl), intent(in) :: fC(:)
  complex(dlc), intent(in), dimension(0:tL(2),0:tL(2)) :: Blm1, Blm2
  complex(dlc), intent(out), dimension(0:eL(2),0:eL(2)) :: glm, clm
  !internal
  integer :: l, npixt
  real(dl), dimension(:,:), allocatable :: A, A1, A3, map
  complex(dlc), dimension(:,:,:), allocatable :: alm1, alm3, tlm

  write(*,*) 'calc BB-estimator'
  npixt = 12*nsidet**2

  !* convolution
  allocate(alm1(2,0:tL(2),0:tL(2)))
  alm1 = 0d0
  do l = tL(1), tL(2)
    alm1(2,l,:) = Blm1(l,:)
  end do
  allocate(A(0:npixt-1,2))
  call alm2map_spin(nsidet,tL(2),tL(2),2,alm1,A)
  deallocate(alm1)

  allocate(alm1(2,0:tL(2),0:tL(2)),alm3(2,0:tL(2),0:tL(2)))
  alm1 = 0d0;  alm3 = 0d0
  do l = tL(1), tL(2)
    alm1(2,l,:) = fC(l)*Blm2(l,:)*dsqrt(dble((l+2)*(l-1)))
    alm3(2,l,:) = fC(l)*Blm2(l,:)*dsqrt(dble((l-2)*(l+3)))
  end do
  allocate(A1(0:npixt-1,2))
  call alm2map_spin(nsidet,tL(2),tL(2),1,alm1,A1)
  allocate(A3(0:npixt-1,2))
  call alm2map_spin(nsidet,tL(2),tL(2),3,alm3,A3)
  deallocate(alm3)

  allocate(map(0:npixt-1,2))
  !map = A*conjg(A1) - conjg(A)*A3
  map(:,1) = A(:,1)*(A1(:,1)-A3(:,1)) + A(:,2)*(A1(:,2)-A3(:,2))
  map(:,2) = -A(:,1)*(A1(:,2)+A3(:,2)) + A(:,2)*(A1(:,1)+A3(:,1))
  deallocate(A,A1,A3)

  allocate(tlm(2,0:eL(2),0:eL(2)))
  call map2alm_spin(nsidet,eL(2),eL(2),1,map,tlm)
  deallocate(map)

  !* compute glm and clm
  write(*,*) 'compute grad and curl'
  do l = el(1), eL(2)
    glm(l,:) = 0.5d0*dsqrt(dble(l*(l+1)))*tlm(1,l,:)
    clm(l,:) = 0.5d0*dsqrt(dble(l*(l+1)))*tlm(2,l,:)
  end do
  deallocate(tlm)

end subroutine QuadBB


end module recfull


