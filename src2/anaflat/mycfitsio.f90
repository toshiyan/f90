!////////////////////////////////////////////////////!
! * cfitsio tools
!////////////////////////////////////////////////////!

module mycfitsio
  use myconst, only: dl, dlc, pi, iu
  use anaflat, only: map_write
  implicit none

  interface READ_FITSDATA
    module procedure READFITS_2D_TO_1D_DBLE, READFITS_3D_TO_1D_CMPLX, READFITS_2D_TO_1D_CMPLX
  end interface READ_FITSDATA

  private dl, dlc, pi, iu
  private map_write

contains 


subroutine WRITEFITS_2D(nn,f,arr2d)
! taken from cfitsio example code "writeimage"
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  character(*), intent(in) :: f
  double precision, intent(in) :: arr2d(nn(1),nn(2))
  !internal
  logical :: simple = .true., extend = .true.
  integer :: u, stat, bitpix=-32

  stat = 0
  call deletefile(trim(f),stat)
  call FTGIOU(u,stat) !assign file unit 
  call FTINIT(u,trim(f),1,stat) !open new FITS file
  call FTPHPR(u,simple,bitpix,2,nn,0,1,extend,stat) !write header
  call FTPPRD(u,1,1,nn(1)*nn(2),arr2d,stat) !write values
  call printerror(stat)
  call FTCLOS(u,stat) !close FITS file
  call FTFIOU(u,stat) !free Memory 

end subroutine WRITEFITS_2D


subroutine READFITS_2D(nn,f,arr2d)
! read 2d fits data -> output as 2d array
! nn: size of read array
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  character(*), intent(in) :: f
  double precision, intent(inout) :: arr2d(nn(1),nn(2))
  !internal
  logical :: anynull
  integer :: u, rw=0, blocksize, s, rr(2), n

  call FTGIOU(u,s) !assign file unit 
  call FTOPEN(u,trim(f),rw,blocksize,s) !open FITS file with readmode
  call printerror(s)
  call FTGKNJ(u,'NAXIS',1,2,rr,n,s) !read naxis
  if ( n .ne. 2 )  stop 'error: NAXIS is not 2'
  if ( .not. (nn(1)==rr(1).and.nn(2)==rr(2)) ) then
    write(*,*) "expected =", nn
    write(*,*) "readed =", rr
    stop 'error: expected size is not equal to readed size'
  end if
  call FTG2DD(u,0,0.d0,nn(1),nn(1),nn(2),arr2d,anynull,s) !read values
  call FTCLOS(u,s) !close FITS file
  call FTFIOU(u,s) !free Memory 

end subroutine READFITS_2D


subroutine READFITS_2D_TO_1D_DBLE(nn,f,map,nxy,outf,nv,rr)
!read CMB map (MAP) from a fitsfile (f) whose size is (nn)
!(nxy) is the center of (MAP)
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2), nxy(1:2)
  integer, intent(in), optional :: rr(1:2)
  real(dl), intent(in), optional :: nv
  character(*), intent(in) :: f
  character(*), intent(in), optional :: outf
  real(dl), intent(inout) :: map(:)
  !internal
  logical :: anynull
  integer :: u, rw=0, naxes(1:2), blocksize, s, nfound, i, j, n
  real(dl) :: nullval
  real(dl), allocatable :: arr2d(:,:)

  naxes = nn
  if(present(rr)) naxes = rr
  allocate(arr2d(naxes(1),naxes(2)))
  CALL READFITS_2D(naxes,f,arr2d)
  !* cut MAP whose center is at nxy
  n = 1
  do i = 1, naxes(1)
    do j = 1, naxes(2)
      if(i<=nxy(1)+nn(1)/2.and.i>nxy(1)-nn(1)/2.and.j<=nxy(2)+nn(2)/2.and.j>nxy(2)-nn(2)/2) then
        map(n) = arr2d(i,j)
        n = n + 1
      end if
    end do
  end do
  deallocate(arr2d)

  !//// output map ////!
  if(present(outf)) call MAP_WRITE(nn,outf,map)

end subroutine READFITS_2D_TO_1D_DBLE


subroutine READFITS_2D_TO_1D_CMPLX(nn,f,map,nxy,outf,nv)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), nxy(2)
  real(dl), intent(in), optional :: nv
  character(*), intent(in) :: f
  character(*), intent(in), optional :: outf
  complex(dlc), intent(inout) :: map(:)
  !internal
  real(dl) :: arr2d(nn(1)*nn(2))

  if(present(outf)) then
    call READFITS_2D_TO_1D_DBLE(nn,f,arr2d,nxy,outf,nv)
  else
    call READFITS_2D_TO_1D_DBLE(nn,f,arr2d,nxy,nv=nv)
  end if
  map = cmplx(arr2d)

end subroutine READFITS_2D_TO_1D_CMPLX


subroutine READFITS_3D_TO_1D_CMPLX(nn,f,map,nxy,outf,nv)
  implicit none
  !I/O
  integer, intent(in) :: nn(2), nxy(2)
  real(dl), intent(in), optional :: nv
  character(*), intent(in) :: f
  character(*), intent(in), optional :: outf
  complex(dlc), intent(out) :: map(:,:)
  !internal
  logical :: anynull
  integer :: u, rw=0, blocksize, s, rr(3), nfound, i, j, n
  real(dl) :: nullval
  real(dl), allocatable :: arr3d(:,:,:)

  call FTGIOU(u,s)
  call FTOPEN(u,trim(f),rw,blocksize,s)
  call FTGKNJ(u,'NAXIS',1,3,rr,nfound,s)
  if(nfound.ne.3)  stop 'READIMAGE failed to read the NAXISn keywords.'
  nullval = 0.d0
  if(present(nv)) nullval = nv
  allocate(arr3d(rr(1),rr(2),rr(3)))
  call FTG3DD(u,0,nullval,rr(1),rr(2),rr(1),rr(2),rr(3),arr3d,anynull,s)
  call FTCLOS(u,s)
  call FTFIOU(u,s)

  !* cut MAP whose center is at nxy
  n = 1
  do i = 1, rr(1)
    do j = 1, rr(2)
      if(i<=nxy(1)+nn(1)/2.and.i>nxy(1)-nn(1)/2.and.j<=nxy(2)+nn(2)/2.and.j>nxy(2)-nn(2)/2) then
        map(n,1) = cmplx(arr3d(i,j,1))
        map(n,2) = cmplx(arr3d(i,j,2))
        map(n,3) = cmplx(arr3d(i,j,3))
        n = n + 1
      end if
    end do
  end do
  deallocate(arr3d)

  !//// output map ////!
  if(present(outf)) call MAP_WRITE(nn,outf,dble(map(:,1)))

end subroutine READFITS_3D_TO_1D_CMPLX


subroutine deletefile(f,stat)
  implicit none
  character(*), intent(in) :: f
  integer, intent(inout) :: stat
  integer :: u, blocksize

  if (stat>0) return
  CALL FTGIOU(u,stat)
  CALL FTOPEN(u,f,1,blocksize,stat)
  if (stat==0) then
    CALL FTDELT(u,stat)
  else if (stat==103) then
    stat = 0
    CALL FTCMSG
  else
    stat = 0
    CALL FTCMSG
    CALL FTDELT(u,stat)
  end if
  CALL FTFIOU(u,stat)

end subroutine deletefile


subroutine printerror(stat)
  implicit none
  integer, intent(in) :: stat
  character(LEN=30) :: err_text
  character(LEN=80) err_message

  if (stat<=0) return
  call FTGERR(stat,err_text) !get corresponding error text
  call FTGMSG(err_message)
  do while (err_message .ne. ' ')
    write(*,*) err_message
    call FTGMSG(err_message)
  end do
  stop "error in fitsio subroutine"

end subroutine printerror


end module mycfitsio

