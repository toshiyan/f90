module forecast_snr
  use myconst, only: dl
  use myutils, only: arr
  use forecast, only: covcl, cov_noise_cmb_l, cov_noise_gal_l, magcorrection
  implicit none

  private dl, arr
  private covcl, cov_noise_cmb_l, cov_noise_gal_l, magcorrection

contains


subroutine SNR_AUTO(Cl,Nl,el,SNRl,VARl,f)
  implicit none
  !I/O
  character(*), intent(in), optional :: f
  integer, intent(in) :: el(2)
  real(dl), intent(in) :: Cl(:), Nl(:)
  real(dl), intent(out), optional :: SNRl(:),VARl(:)
  !internal
  integer :: l
  real(dl) :: iVAR(el(2)), iSNR(el(2))
  
  iSNR = 0d0
  do l = el(1), el(2)
    iVAR(l) = (Cl(l)+Nl(l))/dsqrt(l+0.5d0)
    iSNR(l) = Cl(l)/iVAR(l)
    if(present(VARl)) VARl(l) = iVAR(l)
    if(present(SNRl)) SNRl(l) = iSNR(l)
  end do
  write(*,*) "S/N =", dsqrt(sum(iSNR**2))

  if(present(f)) CALL SNR_OUTPUT(f,el,iSNR,iVAR)

end subroutine SNR_AUTO


subroutine SNR_CROSS(Cl1,Cl2,Clc,Nl1,Nl2,el,SNRl,VARl,f)
  implicit none
  !I/O
  character(*), intent(in), optional :: f
  integer, intent(in) :: el(2)
  real(dl), intent(in), dimension(:) :: Cl1, Cl2, Clc, Nl1, Nl2
  real(dl), intent(out), optional :: SNRl(:),VARl(:)
  !internal
  integer :: l
  real(dl) :: iVAR(el(2)), iSNR(el(2))

  iSNR = 0d0
  do l = el(1), el(2)
    iVAR(l) = dsqrt((Clc(l)**2+(Cl1(l)+Nl1(l))*(Cl2(l)+Nl2(l)))/dble(2*l+1))
    iSNR(l) = Clc(l)/dsqrt(iVAR(l))
    if(present(VARl)) VARl(l) = dsqrt(iVAR(l))
    if(present(SNRl)) SNRl(l) = iSNR(l)
  end do
  write(*,*) "S/N =", dsqrt(sum(iSNR**2))

  if(present(f)) CALL SNR_OUTPUT(f,el,iSNR,iVAR)

end subroutine SNR_CROSS


subroutine SNR_OUTPUT(fname,el,SNRl,VARl)
  implicit none
  !I/O
  character(*), intent(in) :: fname
  integer, intent(in) :: el(2)
  real(dl), intent(in) :: SNRl(:)
  real(dl), intent(in), optional :: VARl(:)
  !internal
  integer :: l

  open(20,file=fname,status='replace')
  do l = el(1), el(2)
    if(present(VARl)) then 
      write(20,"(I6,2(1X,E12.5))") l, SNRl(l), VARl(l)
    else
      write(20,"(I6,1(1X,E12.5))") l, SNRl(l)
    end if
  end do
  close(20)

end subroutine SNR_OUTPUT


subroutine SNR_ALL_SCAL(el,Nls,Ngal,eint,fname,m)
  implicit none
  !I/O
  character(*), intent(in) :: fname
  integer, intent(in) :: el(2)
  real(dl), intent(in) :: Nls(:,:), Ngal(:), eint
  real(dl), intent(in), optional :: m(:)
  !internal
  integer :: l, X, Y, n
  real(dl), dimension(:,:,:), allocatable :: Cl, Nl

  allocate(O%S(O%nz),O%g(O%nz))
  O%T = 1
  O%E = 2
  O%d = 3
  O%S = arr((/4,O%nz/),O%nz)
  O%g = arr((/4+O%nz,2*O%nz/),O%nz)
  n = 3+2*O%nz
  O%DO = .true.

  allocate(Cl(el(2),n,n),Nl(el(2),n,n))
  call CovCl(O,Cl,fname)
  if(present(m)) call MagCorrection(O,m,Cl) !Magnification Correction
  do l = el(1), el(2)
    call cov_noise_cmb_l(O,Nls(:,l),Nl(l,:,:))
    call cov_noise_gal_l(O,Ngal,eint,Nl(l,:,:))
  end do

  do X = 1, n
    CALL SNR_AUTO(Cl(:,X,X),Nl(:,X,X),el)
    do Y = X+1, n
      CALL SNR_CROSS(Cl(:,X,Y),Cl(:,X,X),Cl(:,Y,Y),Nl(:,X,X),Nl(:,Y,Y),el)
    end do
  end do

  deallocate(O%S,O%g,Cl,Nl)

end subroutine SNR_ALL_SCAL


end module forecast_snr

