module cosmoparams
  use readfile, only: read_prm, read_str, read_val
  use myutils,  only: str
  implicit none
  type cosmoparamsset
    integer :: maxn = 18
    integer :: n, e(2), sn, wn
    integer :: Obh,Omh,O_L,n_s,A_s,tau
    integer :: w_0,w_a,mnu,Nef
    integer :: fnl,b_0,b_z
    integer :: TSR,n_t
    integer :: A_L
    integer, allocatable :: s_z(:), w_p(:)
    logical, dimension(18) :: DO, ln
    character(LEN=3), dimension(18) :: lab
    double precision, dimension(18) :: fid, der
  end type cosmoparamsset

  private str, read_prm, read_str, read_val

contains


subroutine settype_cosmoparams(P)
!//// main routine for setting cosmological parameteres ////!
  implicit none
  type(cosmoparamsset), intent(out) :: p
  character(len=3) :: def(18) = ['Obh','Omh','O_L','n_s','A_s','tau','w_0','w_a','w_p','mnu','Nef','fnl','b_0','b_z','s_z','TSR','n_t','A_L']
  character(len=128) :: inline
  integer :: n, i

  !* read # of the lensing slope (sn) and DE-EoS parameteres (wn) 
  p%sn = 0
  p%wn = 0
  do i = 1, 100
    if(read_val('sz'//str(i))) p%sn = p%sn + 1
    if(read_val('wp'//str(i))) p%wn = p%wn + 1
  end do

  !* total number of marginalized parameters (n) and label (lab)
  n = 1
  do i = 1, p%maxn
    if(.not.read_val(def(i))) cycle
    p%lab(n) = def(i) 
    n = n + 1
  end do

  p%n = n
  do i = 1, p%maxn
    if(read_val(def(i))) cycle
    p%lab(n) = def(i) 
    n = n + 1
  end do

  !* read fiducial value (fid) and stepsize (der)
  do i = 1, p%n
    write(*,*) 'read ', p%lab(n)
    inline = read_str(p%lab(n))
    read(inline,*) p%fid(n), p%der(n)
  end do

  !* set parameter IDs
  !LCDM
  p%Obh = param_id('Obh',p%LAB)
  p%Omh = param_id('Omh',p%LAB)
  p%O_L = param_id('O_L',p%LAB)
  p%n_s = param_id('n_s',p%LAB)
  p%A_s = param_id('A_s',p%LAB)
  p%tau = param_id('tau',p%LAB)
  !DENU
  p%w_0 = param_id('w_0',p%LAB)
  p%w_a = param_id('w_a',p%LAB)
  p%mnu = param_id('mnu',p%LAB)
  p%Nef = param_id('Nef',p%LAB)
  !
  p%fnl = param_id('fnl',p%LAB)
  p%b_0 = param_id('b_0',p%LAB)
  p%b_z = param_id('b_z',p%LAB)
  !TENS
  p%TSR = param_id('TSR',p%LAB)
  p%n_t = param_id('n_t',p%LAB)
  !Lensing amplitude
  p%A_L = param_id('A_L',p%LAB)

  !* check
  write(*,*) ((p%lab(i),p%fid(i),p%der(i)),i=1,p%n)

  !* read two parameter ids of 2D contour ellipse 
  if(read_val('eID')) call read_prm('eID',p%e)


end subroutine settype_cosmoparams


function param_id(str,LAB)
  implicit none
  character(*), intent(in) :: str, LAB(:)
  integer :: param_id
  integer :: p

  do p = 1, size(LAB)
    if(trim(LAB(p))==str) param_id = p
  end do

end function param_id


end module cosmoparams

