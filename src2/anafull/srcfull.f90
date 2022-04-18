!////////////////////////////////////////////////////!
! * Point source reconstruction in Fullsky
!////////////////////////////////////////////////////!

module srcfull
  use alm_tools, only: alm2map, map2alm
  use myconst, only: dl, dlc, iu

  interface quadtt_src
    module procedure quadtt_src_sym, quadtt_src_asym
  end interface

  private alm2map, map2alm
  private dl, dlc, iu

contains 


!* Tlm,Elm,Blm: inverse-variance filtered alms (Cov^-1 (Tlm,Elm,Blm)^t)

subroutine quadtt_src_sym(nside,Tlm,eL,rL,xlm)
  implicit none
  !I/O
  integer :: eL(2), rL(2), nside
  complex(dlc), intent(in), dimension(0:eL(2),0:eL(2)) :: Tlm
  complex(dlc), intent(inout), dimension(0:eL(2),0:eL(2)) :: xlm
  !internal
  integer :: l, npix
  double precision, allocatable :: map(:)
  complex(dlc), allocatable :: alm(:,:,:)

  write(*,*) 'calc TT-estimator (src)'
  npix = 12*nside**2

  !* alm to map 
  allocate(alm(1,0:rL(2),0:rL(2))); alm = 0d0
  do l = rL(1), rL(2)
    alm(1,l,0:l) = Tlm(l,0:l)
  end do 
  allocate(map(0:npix-1))
  call alm2map(nside,rL(2),rL(2),alm(1:1,:,:),map)
  deallocate(alm)

  !* map to alm
  allocate(alm(1,0:rL(2),0:rL(2)))
  call map2alm(nside,rL(2),rL(2),map**2,alm)
  xlm = 0d0
  do l = 0, rL(2)
    xlm(l,0:l) = alm(1,l,0:l)
  end do
  deallocate(map,alm)

end subroutine quadtt_src_sym


subroutine quadtt_src_asym(nside,Tlm1,Tlm2,eL,rL,xlm)
  implicit none
  !I/O
  integer :: eL(2), rL(2), nside
  complex(dlc), intent(in), dimension(0:eL(2),0:eL(2)) :: Tlm1, Tlm2
  complex(dlc), intent(inout), dimension(0:eL(2),0:eL(2)) :: xlm
  !internal
  integer :: l, npix
  double precision, allocatable :: map(:,:)
  complex(dlc), allocatable :: alm(:,:,:)

  write(*,*) 'calc TT-estimator (src)'
  npix = 12*nside**2

  !* alm to map 
  allocate(alm(2,0:rL(2),0:rL(2))); alm = 0d0
  do l = rL(1), rL(2)
    alm(1,l,0:l) = Tlm1(l,0:l)
    alm(2,l,0:l) = Tlm2(l,0:l)
  end do 
  allocate(map(2,0:npix-1))
  call alm2map(nside,rL(2),rL(2),alm(1:1,:,:),map(1,:))
  call alm2map(nside,rL(2),rL(2),alm(2:2,:,:),map(2,:))
  deallocate(alm)

  !* map to alm
  allocate(alm(1,0:rL(2),0:rL(2)))
  call map2alm(nside,rL(2),rL(2),map(1,:)*map(2,:),alm)
  xlm = 0d0
  do l = 0, rL(2)
    xlm(l,0:l) = alm(1,l,0:l)
  end do
  deallocate(map,alm)

end subroutine quadtt_src_asym


end module srcfull


