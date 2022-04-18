!//////////////////////////////////////////////
! * To be removed
!//////////////////////////////////////////////

module mytypes
  use readFile
  use myconst, only: dl, dlc, pi, ac2rad
  use funcs, only: lnGamma
  use myutils, only: param_id, str, ran, FileColumns, FileLines
  use utils_galaxy, only: HSC, DES, LSST, zbin_SF
  implicit none

  private dl, dlc, pi, ac2rad
  private lnGamma
  private param_id, str, ran
  private filecolumns, filelines
  private HSC, DES, LSST, zbin_SF

  !* anisotropies
  type OBSERVABLES
    character(LEN=2) :: LAB(128)
    character(LEN=4) :: CLLAB(128)
    integer :: T, E, B, d, o
    integer :: TT, TE, TB, EE, EB, BB, Td, Ed, dd, oo, Bo, CLID(128)
    integer :: n, nz, maxn = 7
    integer :: cln, clmaxn
    integer, allocatable :: S(:), g(:)
    integer, allocatable, dimension(:) :: TS,SS,Tg,gg,dS,dg,Sg
    real(dl), allocatable, dimension(:,:) :: UC, LC, OC
    logical :: DO(128), CLDO(128)
    logical :: noTS, noSg, onlyTg, onlySg
  end type OBSERVABLES

  type QUADRATIC_COMBINATION
    character(LEN=2) :: LAB(7)
    integer :: TT, TE, TB, EE, EB, BB, MV, n
    logical :: DO(7), DO_TTTE, DO_TTEE, DO_TEEE, DO_TBEB
    integer :: TTTE=1, TTEE=2, TEEE=3, TBEB=4
  end type QUADRATIC_COMBINATION

  type LENSING_ESTIMATOR
    character(LEN=2) :: LAB(64)
    character(LEN=4) :: CLLAB(64)
    character(LEN=64), allocatable :: alfile(:), nlfile(:)
    integer :: g, c, m, s, t, p
    integer :: gg, gc, gm, gs, gt, gp, cc, cm, cs, mm, ms, ss, tt, pp
    integer :: n, maxn = 6
    integer :: cln, clmaxn
    integer :: dL(2), rL(2)
    integer :: arr(6,6)
    logical :: DO(64), CLDO(64), use_L
  end type LENSING_ESTIMATOR

  type GalaxySurvey
    integer :: nz
    real(dl) :: a, b, sigma, zbias, z0, zm, Ngal, fsky, eint
    real(dl), dimension(:), allocatable :: zb, Ng
  end type GalaxySurvey

contains

subroutine SET_PARAMS_OBS(n,DEF,LAB,OBSDO)
  implicit none
  integer, intent(inout) :: n
  character(*), intent(in) :: DEF
  character(*), intent(out) :: LAB
  logical, intent(out), optional :: OBSDO

  LAB = trim(DEF)
  if(present(OBSDO)) OBSDO = .true.
  n = n + 1

end subroutine SET_PARAMS_OBS


subroutine SET_PARAMS_CL(n,DEF,CLLAB,CLDO,CLID)
  implicit none
  integer, intent(inout) :: n
  character(*), intent(in) :: DEF
  character(*), intent(out) :: CLLAB
  integer, intent(out), optional :: CLID
  logical, intent(out), optional :: CLDO

  CLLAB = trim(DEF)
  if(present(CLDO)) CLDO = READ_VAL(DEF)
  if(present(CLID)) CLID = READ_INT(DEF)
  n = n + 1

end subroutine SET_PARAMS_CL



subroutine SET_TYPE_OBS_CMBDEF(O)
  implicit none
  type(OBSERVABLES), intent(out) :: O
  integer :: n, m, i, j
  character(2) :: OBS(11)
  character(4) :: DEF(11)

  OBS(1:3) = (/"T","E","B"/)
  n = 0
  do i = 1, 3
    do j = i, 3
      n = n + 1
      DEF(n) = trim(OBS(i))//trim(OBS(j))
    end do
  end do
  DEF(7) = "dd"
  DEF(8) = "Td"
  DEF(9) = "Ed"
  DEF(10) = "oo"
  DEF(11) = "Bo"
  O%cln = 11

  n = 1
  do i = 1, O%cln
    call set_params_cl(n,DEF(i),O%CLLAB(n))
  end do

  O%TT = PARAM_ID("TT",O%CLLAB)
  O%TE = PARAM_ID("TE",O%CLLAB)
  O%TB = PARAM_ID("TB",O%CLLAB)
  O%EE = PARAM_ID("EE",O%CLLAB)
  O%EB = PARAM_ID("EB",O%CLLAB)
  O%BB = PARAM_ID("BB",O%CLLAB)
  O%dd = PARAM_ID("dd",O%CLLAB)
  O%Td = PARAM_ID("Td",O%CLLAB)
  O%Ed = PARAM_ID("Ed",O%CLLAB)
  O%oo = PARAM_ID("oo",O%CLLAB)
  O%Bo = PARAM_ID("Bo",O%CLLAB)

end subroutine SET_TYPE_OBS_CMBDEF



subroutine SET_TYPE_OBS(O)
!* set label and number for each signals
  implicit none
  type(OBSERVABLES), intent(inout) :: O
  integer :: n, m, i, j, nobs
  integer, parameter :: ncmb = 5, ngal = 3
  logical :: OBSDO(128)
  character(2) :: OBS(128)
  character(4) :: DEF(128)

  OBS(1:5) = (/"T","E","B","d","o"/)
  nobs = ncmb + 2*ngal
  O%DO = .false.
  O%CLDO = .false.
  OBSDO = .false.

  !number of redshifts
  O%nz = 0
  do i = 1, ngal
    if(READ_VAL("S"//str(1)//"S"//str(i))) O%nz = O%nz + 1
    OBS(i+ncmb) = "S"//str(i) 
    OBS(i+ncmb+ngal) = "g"//str(i) 
  end do
  write(*,*) "number of redshifts =", O%nz

  !label
  allocate(O%S(O%nz),O%g(O%nz))
  n = 0
  do i = 1, nobs
    do j = 1, nobs
      n = n + 1
      DEF(n) = trim(OBS(i))//trim(OBS(j))
    end do
  end do

  !* Ordering Cl labels
  n = 1
  do i = 1, nobs**2
    if(READ_VAL(DEF(i))) call set_params_cl(n,DEF(i),O%CLLAB(n),O%CLDO(n),O%CLID(n))
  end do

  O%cln = n - 1

  do i = 1, nobs**2
    if(.not.READ_VAL(DEF(i))) call set_params_cl(n,DEF(i),O%CLLAB(n))
  end do


  !* ID's for Cls
  O%TT = PARAM_ID("TT",O%CLLAB)
  O%TE = PARAM_ID("TE",O%CLLAB)
  O%TB = PARAM_ID("TB",O%CLLAB)
  O%EE = PARAM_ID("EE",O%CLLAB)
  O%EB = PARAM_ID("EB",O%CLLAB)
  O%BB = PARAM_ID("BB",O%CLLAB)
  O%dd = PARAM_ID("dd",O%CLLAB)
  O%Td = PARAM_ID("Td",O%CLLAB)
  O%Ed = PARAM_ID("Ed",O%CLLAB)
  O%oo = PARAM_ID("oo",O%CLLAB)

  allocate(O%TS(O%nz),O%Tg(O%nz),O%dS(O%nz),O%dg(O%nz))
  allocate(O%SS(O%nz*(O%nz+1)/2),O%gg(O%nz*(O%nz+1)/2))
  allocate(O%Sg(O%nz**2))

  n = 0
  m = 0
  do i = 1, O%nz
    O%TS(i) = PARAM_ID("TS"//str(i),O%CLLAB)
    O%Tg(i) = PARAM_ID("Tg"//str(i),O%CLLAB)
    O%dS(i) = PARAM_ID("dS"//str(i),O%CLLAB)
    O%dg(i) = PARAM_ID("dg"//str(i),O%CLLAB)
    do j = i, O%nz
      n = n + 1
      O%SS(n) = PARAM_ID("S"//str(i)//"S"//str(j),O%CLLAB)
      O%gg(n) = PARAM_ID("g"//str(i)//"g"//str(j),O%CLLAB)
    end do
    do j = 1, O%nz
      m = m + 1
      O%Sg(m) = PARAM_ID("S"//str(i)//"g"//str(j),O%CLLAB)
    end do
  end do


  !* Type of signal used
  OBSDO = .false.
  do i = 1, nobs
    do j = 1, nobs
      if(READ_VAL(trim(OBS(i))//trim(OBS(j)))) then
        OBSDO(i) = .true.
        OBSDO(j) = .true.
      end if
    end do
  end do

  !* Ordering signal labels
  n = 1
  do i = 1, nobs
    if(OBSDO(i)) call set_params_obs(n,OBS(i),O%LAB(n),O%DO(n))
  end do
  O%n = n - 1
  do i = 1, nobs
    if(.not.OBSDO(i)) call set_params_obs(n,OBS(i),O%LAB(n))
  end do

  !* ID's for signals
  O%T = PARAM_ID("T",O%LAB)
  O%E = PARAM_ID("E",O%LAB)
  O%B = PARAM_ID("B",O%LAB)
  O%d = PARAM_ID("d",O%LAB)
  do i = 1, O%nz
    O%S(i) = PARAM_ID("S"//str(i),O%LAB)
    O%g(i) = PARAM_ID("g"//str(i),O%LAB)
  end do

  !check
  write(*,*) "read Cls"
  do i = 1, O%cln
    write(*,"(I2,1X,A4,1X,I2)") i, O%CLLAB(i), O%CLID(i)
  end do

  write(*,*) "elements of covariance"
  do i = 1, O%n
    write(*,"(I2,1X,A2)") i, O%LAB(i)
  end do

  if(O%cln<1) then
    stop "error: cannot read cls (specify which Cls to be used)"
  end if

end subroutine SET_TYPE_OBS


subroutine FREE_TYPE_OBS(O)
  implicit none
  type(OBSERVABLES), intent(inout) :: O

  deallocate(O%S,O%g)
  deallocate(O%TS,O%Tg,O%dS,O%dg,O%SS,O%gg,O%Sg)

end subroutine FREE_TYPE_OBS


subroutine SET_TYPE_QUAD(Q)
!* Label and ID for quadratic combination (Q%LAB)
  implicit none
  type(QUADRATIC_COMBINATION), intent(out) :: Q
  logical :: rDO(6)
  integer :: n, i, qc(6)
  character(LEN=2) :: DEF(7) = (/"TT","TE","TB","EE","EB","BB","MV"/)

  call READ_PARAMS('DO_QUAD',rDO)

  n = 0
  qc = 7
  Q%DO = .false.
  do i = 1, 6
    if(rDO(i)) then
      n = n + 1
      qc(i) = n
      Q%DO(qc(i)) = .true.
      Q%LAB(qc(i)) = trim(DEF(i))
    end if
  end do

  Q%TT = qc(1)
  Q%TE = qc(2)
  Q%TB = qc(3)
  Q%EE = qc(4)
  Q%EB = qc(5)
  Q%BB = qc(6)

  if(n>0) then
    Q%MV = n + 1
    Q%LAB(Q%MV) = trim(DEF(7))
  end if

  Q%n = n
  write(*,*) "qmax=", Q%n

  Q%DO_TTTE = Q%DO(Q%TT).and.Q%DO(Q%TE)
  Q%DO_TTEE = Q%DO(Q%TT).and.Q%DO(Q%EE)
  Q%DO_TEEE = Q%DO(Q%TE).and.Q%DO(Q%EE)
  Q%DO_TBEB = Q%DO(Q%TB).and.Q%DO(Q%EB)

end subroutine SET_TYPE_QUAD


subroutine SET_TYPE_EST(E)
!* set label and number for each estimators
  implicit none
  !I/O
  type(LENSING_ESTIMATOR), intent(inout) :: E
  !internal
  logical :: CLDO(64)
  integer :: n, m, i, j
  integer, parameter :: nobs = 6
  character(1), parameter :: OBS(1:6) = (/"g","c","m","s","t","p"/)
  character(2) :: DEF(64)

  E%DO = .false.
  E%CLDO = .false.

  !* Ordering signal labels
  n = 1
  do i = 1, nobs
    if(READ_VAL(OBS(i))) call set_params_obs(n,OBS(i),E%LAB(n),E%DO(n))
  end do
  E%n = n - 1
  do i = 1, nobs
    if(.not.READ_VAL(OBS(i))) call set_params_obs(n,OBS(i),E%LAB(n))
  end do

  !* ID's for signals
  E%g = PARAM_ID("g",E%LAB)
  E%c = PARAM_ID("c",E%LAB)
  E%m = PARAM_ID("m",E%LAB)
  E%s = PARAM_ID("s",E%LAB)
  E%t = PARAM_ID("t",E%LAB)
  E%p = PARAM_ID("p",E%LAB)

  !* Als to be used and Al-labels
  CLDO = .false.
  n = 1
  m = 1
  E%arr = 0
  do i = 1, nobs
    do j = i, nobs
      DEF(n) = trim(OBS(i))//trim(OBS(j))
      if(READ_VAL(OBS(i)).and.READ_VAL(OBS(j))) then
        CLDO(n) = .true.
        if(.not.DEF(n)=="gc") then
          E%arr(i,j) = m
          E%arr(j,i) = m
          m = m + 1
        end if
      end if
      n = n + 1
    end do
  end do

  !* Ordering Cl labels
  i = 1
  do n = 1, nobs**2
    if(CLDO(n)) call set_params_cl(i,DEF(n),E%CLLAB(i),E%CLDO(i))
  end do
  E%cln = i - 1
  do n = 1, nobs**2
    if(.not.CLDO(n)) call set_params_cl(i,DEF(n),E%CLLAB(i))
  end do

  !* ID's for Cls
  E%gg = PARAM_ID("gg",E%CLLAB)
  E%gc = PARAM_ID("gc",E%CLLAB)
  E%gm = PARAM_ID("gm",E%CLLAB)
  E%gs = PARAM_ID("gs",E%CLLAB)
  E%cc = PARAM_ID("cc",E%CLLAB)
  E%cm = PARAM_ID("cm",E%CLLAB)
  E%cs = PARAM_ID("cs",E%CLLAB)
  E%mm = PARAM_ID("mm",E%CLLAB)
  E%ms = PARAM_ID("ms",E%CLLAB)
  E%ss = PARAM_ID("ss",E%CLLAB)

  !check
  write(*,*) "read Cls"
  do i = 1, E%cln
    write(*,"(I2,1X,A4,1X,I2)") i, E%CLLAB(i)
  end do

  write(*,*) "elements of covariance"
  do i = 1, E%n
    write(*,"(I2,1X,A2)") i, E%LAB(i)
  end do

  if(E%n<1) stop "error: specify which estimator to be used"

  call READ_PARAMS('dL',E%dL)
  call READ_PARAMS('rL',E%rL)
  E%use_L = READ_LOGICAL("use_LCl")

end subroutine SET_TYPE_EST


subroutine SET_NORMALIZATION_FILES(Q,E)
  implicit none
  !I/O
  type(QUADRATIC_COMBINATION), intent(in) :: Q
  type(LENSING_ESTIMATOR), intent(inout) :: E
  !internal
  character(LEN=64) :: root
  integer :: i

  if(.not.allocated(E%alfile)) allocate(E%alfile(Q%MV))
  if(.not.allocated(E%nlfile)) allocate(E%nlfile(Q%n))

  root = read_str('Alroot')
  do i = 1, Q%n
    E%alfile(i) = trim(root)//"Al"//trim(Q%LAB(i))//".dat"
    E%nlfile(i) = trim(root)//"Nl"//trim(Q%LAB(i))//".dat"
  end do
  E%alfile(Q%MV) = trim(root)//"Al"//trim(Q%LAB(Q%MV))//".dat"

end subroutine SET_NORMALIZATION_FILES


subroutine READCL(O,Cl,f,el)
!//// read cls from a CAMB-formated file 
  implicit none 
  !I/O
  type(OBSERVABLES), intent(in) :: O
  character(*), intent(in) :: f
  integer, intent(in) :: el(:)
  real(dl), intent(out) :: Cl(:,:)
  !internal
  real(dl), allocatable, dimension(:) :: rCl
  integer :: l, ell, n, m, i, j, z, rn

  Cl = 0.d0
  write(*,*) "read Cls from", trim(adjustl(f))

  open(unit=20,file=trim(adjustl(f)),status="old")
  rn = FileColumns(20)

  !check
  do i = 1, O%cln
    if(O%CLID(i)>=rn) stop "error: input cl file has not enough columns"
  end do

  !read
  allocate(rCl(rn-1)) 
  m = FileLines(20)
  do l = 1, m
    read(20,*) ell, rCl(1:rn-1)
    rCl = rCl*2.d0*pi/dble(ell**2+ell)
    if(ran(ell,el).or.ell==el(2)) then 
      do i = 1, O%cln
        Cl(i,ell) = rCl(O%CLID(i))
      end do
    end if
  end do
  deallocate(rCl)

  close(20)


end subroutine READCL


subroutine READCL_CAMB(O,Cl,f,el,HasBB)
  implicit none 
  !I/O
  type(OBSERVABLES), intent(in) :: O
  character(*), intent(in) :: f
  logical, intent(in), optional :: HasBB
  integer, intent(in) :: eL(2)
  real(dl), intent(out) :: Cl(:,:)
  !internal
  real(dl) :: rC(10)
  integer :: i, l, n, m

  open(unit=20,file=trim(f),status="old")
  n = FileColumns(20)
  m = FileLines(20)
  do i = 1, m
    read(20,*) l, rC(1:n-1)
    rC = rC*2.d0*pi/dble(l**2+l)
    if(ran(l,el).or.l==el(2)) then 
      CL(O%TT,l) = rC(1)
      CL(O%EE,l) = rC(2)
      if(present(HasBB)) then
        CL(O%BB,l) = rC(3)
        CL(O%TE,l) = rC(4)
        if(n>=6) CL(O%oo,l) = rC(5)
      else
        CL(O%TE,l) = rC(3)
        CL(O%dd,l) = rC(4)
        CL(O%Td,l) = rC(5)
        if(n>=7) CL(O%Ed,l) = rC(6)
      end if
    end if
  end do
  close(20)

end subroutine READCL_CAMB



subroutine GalaxySurveyParamREAD(GS,survey)
  implicit none
  !I/O
  type(GalaxySurvey), intent(inout) :: GS
  character(*), intent(in), optional :: survey
  !internal
  type(HSC) :: tHSC
  type(DES) :: tDES
  type(LSST) :: tLSST

  !* shape of galaxy distribution function
  GS%a  = READ_DOUBLE('alpha')
  GS%b  = READ_DOUBLE('beta')

  !* photo-z error
  GS%sigma = READ_DOUBLE('sigma')

  !* redshift bias
  GS%zbias = Read_DOUBLE('zbias')

  if (present(survey)) then
    if (survey=="HSC") then
      GS%zm   = tHSC%zm
      GS%fsky = tHSC%fsky
      GS%Ngal = tHSC%Ngal
    else if (survey=="DES") then
      GS%zm   = tDES%zm
      GS%fsky = tDES%fsky
      GS%Ngal = tDES%Ngal
    else if (survey=="LSST") then
      GS%zm   = tLSST%zm
      GS%fsky = tLSST%fsky
      GS%Ngal = tLSST%Ngal
    end if
  else
    GS%zm = READ_DOUBLE('zm')  !* mean redshift
    GS%fsky = READ_DOUBLE("fsky")  !* sky coverage
    GS%Ngal = READ_DOUBLE("Ng")  !* projected numer density [/arcmin^2]
  end if

  GS%z0 = GS%zm*dexp(lnGamma((GS%a+1)/GS%b))/dexp(lnGamma((GS%a+2)/GS%b))
  GS%Ngal = GS%Ngal/ac2rad**2


  !* intrinsic ellipticity
  GS%eint = READ_DOUBLE('eint')


end subroutine GalaxySurveyParamREAD


subroutine GalaxySurveyParam(GS,survey)
  implicit none
  type(GalaxySurvey), intent(inout) :: GS
  character(*), intent(in), optional :: survey

  if(present(survey)) then
    CALL GalaxySurveyParamRead(GS,survey)
  else
    CALL GalaxySurveyParamRead(GS)
  end if

  !* z-binning and galaxy number density
  allocate(GS%Ng(GS%nz),GS%zb(GS%nz+1))
  call zbin_SF(GS%a,GS%b,GS%z0,GS%zb)
  GS%Ng = GS%Ngal/GS%nz
  write(*,*) GS%Ng
  deallocate(GS%zb)

  write(*,*) 'n(z) shape (alpha, beta, z0) =', GS%a,GS%b,GS%z0
  write(*,*) "fsky, Ngal =", GS%fsky,GS%Ngal*ac2rad**2
  write(*,*) "number of galaxies at each bin", GS%Ng*ac2rad**2

end subroutine GalaxySurveyParam

subroutine FREE_TYPE_GS(GS)
  implicit none
  type(GalaxySurvey), intent(inout) :: GS

  deallocate(GS%Ng)

end subroutine FREE_TYPE_GS


end module mytypes


