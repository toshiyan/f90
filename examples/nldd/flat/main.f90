!////////////////////////////////////////////////////!
! * Normalization in Flat Sky
!////////////////////////////////////////////////////!

program main
  use readfile,  only: set_params_file, read_prm, read_str, read_log, read_dbl
  use myconst,   only: TT, TE, EE, BB, dd, Tcmb, ac2rad
  use myutils,   only: linspace, savetxt
  use mycls,     only: readcl_camb
  use nldd_flat, only: alxy_flat
  implicit none
  logical :: CV
  integer :: l, oL(2), eL(2), rL(2)
  double precision, dimension(:), allocatable :: Nl, Ag, Ac, At, Ar
  double precision, dimension(:,:), allocatable :: UC, LC, W

  call set_params_file
  call read_prm('rL',rL)
  call read_prm('oL',oL)
  CV  = read_log('CV')

  !* read theoretical CMB power spectrum
  eL(1) = min(rL(1),oL(1))
  eL(2) = max(rL(2),oL(2))
  allocate(UC(7,eL(2)),LC(7,eL(2)),W(3,eL(2)))
  call readcl_camb(UC,read_str('ucl'),eL)
  call readcl_camb(LC,read_str('lcl'),eL,.true.)

  !* experiment noise
  allocate(Nl(eL(2)))
  if(.not.CV) then
    Nl = (read_dbl('sT')*ac2rad/Tcmb)**2*dexp(linspace(1,eL(2))*(linspace(1,eL(2))+1d0)*(read_dbl('fwhm')*ac2rad)**2/(8d0*dlog(2d0)))
    do l = 2, eL(2)
      W(1,l) = 1d0/(LC(TT,l)+Nl(l))
      W(2,l) = 1d0/(LC(EE,l)+Nl(l)*2d0)
      W(3,l) = 1d0/(LC(BB,l)*read_dbl('alpha')+Nl(l)*2d0)
    end do
    !call savetxt('nlcl.dat',linspace(1,eL(2)),LC(BB,:),2d0*Nl,ow=.true.)
  end if
  deallocate(Nl)

  !Main Calculations
  allocate(Ag(oL(2)),Ac(oL(2)),At(oL(2)),Ar(oL(2))); Ag=0d0; Ac=0d0; At=0d0; Ar=0d0
  call alxy_flat('TB',rL,oL,Ag,Ac,LC(EE,:),W(1,:),W(3,:),weight='rotation')
  call alxy_flat('TT',rL,oL,At,Ar,LC(TT,:),W(1,:),W(1,:),weight='patchytau')
  call savetxt('al.dat',linspace(2,oL(2),oL(2)-1),Ag(2:),Ac(2:),At(2:),Ar(2:),ow=.true.)
  deallocate(Ag,Ac,At,Ar,W,LC,UC)

end program main

