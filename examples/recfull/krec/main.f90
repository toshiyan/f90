! example of fullsky lensing reconstruction

program main
  use readfile,  only: set_params_file, read_prm, read_str, read_int, read_dbl
  use myutils,   only: str, linspace, savetxt, loadtxt, save_average
  use io,        only: loadfile
  use myconst,   only: dlc, pi
  use mycls,     only: alm2bcl, binned_ells
  use nldd_lens, only: AlTT, AlEB
  use anafull,   only: gaussianTEB
  use recfull,   only: quadtt, quadeb
  implicit none
  integer :: i, e, l, bn, nside, rL(2), eL(2), sn(2)
  double precision, allocatable :: bc(:), ll(:), Fl(:,:,:), cl(:,:), cnl(:,:), cpp(:,:,:), Al(:,:,:)
  complex(dlc), allocatable :: alm(:,:,:), glm(:,:,:), clm(:,:,:)

  ! read parameters from params.ini
  call set_params_file
  nside = read_int('nside')
  bn    = read_int('bn')
  call read_prm('sn',sn)
  call read_prm('eL',eL)
  call read_prm('rL',rL)

  ! multipole
  allocate(ll(eL(2)),bc(bn))
  ll = linspace(1,eL(2))
  call binned_ells(eL,bc=bc) !binned multipole

  ! read CMB cl (change here to read your cl)
  allocate(cl(4,eL(2)),cnl(4,eL(2))); cl=0d0; cnl=0d0
  call loadtxt('../../dat/lensedfid_P15.dat',cl(1,2:),cl(2,2:),cl(3,2:),cl(4,2:),rows=[1,eL(2)-1],usecols=[2,3,4,5])
  do i = 1, 4
    cl(i,:) = cl(i,:) * 2d0*pi/(ll**2+ll) ! factor out
  end do

  ! read or define "measured" cl
  cnl = cl ! here no instrumental effects are included

  ! set filtering function
  allocate(Fl(3,0:eL(2),0:eL(2))); Fl=0d0
  do i = 1, 3
    do l = 2, eL(2)
      Fl(i,l,0:l) = 1d0/cnl(i,l)
    end do
  end do

  ! compute normalization
  allocate(Al(2,2,eL(2)))
  call AlTT(rL,eL,Al(1,1,:),Al(1,2,:),cl(1,:),cnl(1,:))
  call AlEB(rL,eL,Al(2,1,:),Al(2,2,:),cl(2,:),cnl(2,:),cnl(3,:))
  !call savetxt('Al.dat',ll,Al(1,1,:),Al(2,1,:),ow=.true.) !output normalization if needed

  ! calculate estimator
  allocate(alm(3,0:eL(2),0:eL(2)),glm(2,0:eL(2),0:eL(2)),clm(2,0:eL(2),0:eL(2)),cpp(sn(1):sn(2),2,bn))

  do i = sn(1), sn(2) ! iterate for realizations

    !//// you need some input alm and read it here ////!
    ! 1) read alm
    !call loadfile('alm_r'//str(i)//'.dat',alm) !change this line to read your data appropriately
    ! 2) or generate random gaussian alms
    call gaussianTEB(alm,cl(1,:),cl(2,:),cl(3,:),cl(4,:),eL(2))

    ! inverse variance filtering
    alm = alm*Fl

    ! compute unnormalized quadratic estimator
    call quadtt(nside,alm(1,0:rL(2),0:rL(2)),alm(1,0:rL(2),0:rL(2)),cl(1,:),eL,rL,glm(1,:,:),clm(1,:,:))
    call quadeb(nside,alm(2,0:rL(2),0:rL(2)),alm(3,0:rL(2),0:rL(2)),cl(2,:),eL,rL,glm(2,:,:),clm(2,:,:))

    ! correct normalization
    do e = 1, 2
      do l = 2, eL(2)
        glm(e,l,0:l) = glm(e,l,0:l)*Al(e,1,l)
        clm(e,l,0:l) = clm(e,l,0:l)*Al(e,2,l)
      end do
      ! compute binned power of the estimators
      call alm2bcl(bn,eL,glm(e,:,:),cb=cpp(i,e,:),oL=eL)
    end do

    ! save to file
    call savetxt('cpp_r'//str(i)//'.dat',bc,cpp(i,1,:),cpp(i,2,:),ow=.true.)

  end do

  ! output average of cl over realizations
  call save_average('cpp.dat',cpp,id=[2,3],bc=bc)

  deallocate(alm,glm,clm,Fl,cl,cnl,ll,cpp,Al,bc)

end program main

