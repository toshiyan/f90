! * Module to read in name/value pairs from a file, with each line of the form line 'name = value'
! * Taken from CAMB originally written by Antony Lewis (ver. Apr 11)

module readfile
  implicit none
  public
  character(LEN=128) :: InLine !namikawa
  integer, parameter :: Ini_max_name_len=128, Ini_max_string_len=1024
  type TNameValue !no known way to make character string pointers..
    character(Ini_max_name_len)  :: Name
    character(Ini_max_string_len):: Value
  end type TNameValue
  type TNameValue_pointer
    Type(TNameValue), pointer :: P
  end type TNameValue_pointer
  Type TNameValueList
    integer Count
    integer Delta
    integer Capacity
    logical ignoreDuplicates
    type(TNameValue_pointer), dimension(:), pointer :: Items
  end Type TNameValueList
  Type TIniFile
    logical SlashComments
    Type (TNameValueList) :: L, ReadValues
  end Type TIniFile
  Type(TIniFile) :: DefIni

  interface read_prm
    module procedure read_prm_int, read_prm_dbl, read_prm_real, read_prm_log
  end interface read_prm

  interface set_prm
    module procedure set_prm_array_int, set_prm_array_dbl, set_prm_int, set_prm_dbl, set_prm_str, set_prm_log
  end interface set_prm

contains

!////////////////////////////////////////////////////////////////////////////////////////////////////!
! modified by namikawa

subroutine set_params_file
  implicit none
  logical :: bad
  character(LEN=128) :: numstr, InputFile, InLine

  InputFile = ''
  if (command_argument_count() /= 0)  call get_command_argument(1,InputFile)
  if (InputFile == '') stop 'No parameter input file'
  call Ini_Open(InputFile,1,bad,.false.)
  if (bad) stop 'Error opening parameter file'

end subroutine set_params_file

subroutine read_prm_int(string,params)
  implicit none
  character(*), intent(in) :: string
  integer, intent(out) :: params(:)

  inline = read_str(trim(adjustl(string)))
  read(inline,*) params(:)

end subroutine read_prm_int

subroutine read_prm_log(string,params)
  implicit none
  character(*), intent(in) :: string
  logical, intent(out) :: params(:)

  inline = read_str(trim(adjustl(string)))
  read(inline,*) params(:)

end subroutine read_prm_log

subroutine read_prm_dbl(string,params)
  implicit none
  character(*), intent(in) :: string
  double precision, intent(out) :: params(:)
  
  inline = read_str(trim(adjustl(string)))
  read(inline,*) params(:)

end subroutine read_prm_dbl

subroutine read_prm_real(string,params)
  implicit none
  character(*), intent(in) :: string
  real(4), intent(out) :: params(:)

  inline = read_str(trim(adjustl(string)))
  read(inline,*) params(:)

end subroutine read_prm_real

function read_str(Key)
  implicit none
  character(*), intent(IN) :: Key
  character(LEN=128) :: read_str

  read_str = Ini_Read_String_File(DefIni,Key)

end function read_str

function read_int(Key)
  implicit none
  character(*), intent(IN) :: Key
  integer :: read_int

  read_int = ini_read_int_file(DefIni,Key)

end function read_int

function read_real(Key)
  implicit none
  character (*), intent(IN) :: Key
  real read_real

  read_real = ini_read_real_file(DefIni,Key)

end function read_real

function read_log(Key)
  implicit none
  character(*), intent(IN) :: Key
  logical :: read_log

  read_log = Ini_read_log_File(DefIni,Key)

end function read_log

function read_dbl(Key)
  implicit none
  character(*), intent(IN) :: Key
  double precision :: read_dbl

  read_dbl = Ini_read_dbl_File(DefIni,Key)
  
end function read_dbl

function read_val(Key)
  implicit none
  character(*), intent(IN) :: Key
  logical :: read_val

  read_val = TNameValueList_HasKey(DefIni%L,Key)

end function read_val

!////////////////////////////////////////////////////////////////////////!
! added

subroutine set_prm_array_int(args,p,def)
  implicit none
  character(*), intent(in) :: args
  integer, intent(in) :: def(:)
  integer, intent(out) :: p(:)

  if (read_val(args)) then
    call read_prm(args,p)
    write(*,*) 'set '//trim(args), p
  else
    p = def
    write(*,*) 'use default for '//trim(args), p
  end if

end subroutine set_prm_array_int

subroutine set_prm_array_dbl(args,p,def)
  implicit none
  character(*), intent(in) :: args
  double precision, intent(in) :: def(:)
  double precision, intent(out) :: p(:)

  if (read_val(args)) then
    call read_prm(args,p)
    write(*,*) 'set '//trim(args), p
  else
    p = def
    write(*,*) 'use default for '//trim(args), p
  end if

end subroutine set_prm_array_dbl

subroutine set_prm_int(args,p,def)
  implicit none
  character(*), intent(in) :: args
  integer, intent(in) :: def
  integer, intent(out) :: p

  if (read_val(args)) then
    p = read_int(args)
    write(*,*) 'set '//trim(args), p
  else
    p = def
    write(*,*) 'use default for '//trim(args), p
  end if

end subroutine set_prm_int

subroutine set_prm_dbl(args,p,def)
  implicit none
  character(*), intent(in) :: args
  double precision, intent(in) :: def
  double precision, intent(out) :: p

  if (read_val(args)) then
    p = read_dbl(args)
    write(*,*) 'set '//trim(args), p
  else
    p = def
    write(*,*) 'use default for '//trim(args), p
  end if

end subroutine set_prm_dbl

subroutine set_prm_str(args,p,def)
  implicit none
  character(*), intent(in) :: args
  character(*), intent(in) :: def
  character(*), intent(out) :: p

  if (read_val(args)) then
    p = read_str(args)
    write(*,*) 'set '//trim(args), ', '//trim(p)
  else
    p = def
    write(*,*) 'use default for '//trim(args), ', '//trim(p)
  end if

end subroutine set_prm_str

subroutine set_prm_log(args,p,def)
  implicit none
  character(*), intent(in) :: args
  logical, intent(in) :: def
  logical, intent(out) :: p

  if (read_val(args)) then
    p = read_log(args)
    write(*,*) 'set '//trim(args), p
  else
    p = def
    write(*,*) 'use default for '//trim(args), p
  end if

end subroutine set_prm_log


!////////////////////////////////////////////////////////////////////////!

function Ini_Read_Int_File(Ini,Key)
  implicit none
  Type(TIniFile), intent(in) :: Ini
  character(*), intent(IN) :: Key
  character(LEN=Ini_max_string_len) :: S
  integer :: Ini_Read_Int_File
   
  S = Ini_Read_String_File(Ini,Key)
  if(S=='') then
    write(*,*) 'no value for key: '//Key
    stop
  else
    if (verify(trim(S),'-+0123456789') /= 0) goto 10
    read (S,*,err=10) Ini_Read_Int_File
  end if
  return

10 write (*,*) 'error reading integer for key: '//Key
  stop
  
end function Ini_Read_Int_File


function Ini_read_dbl_File(Ini,Key)
  implicit none
  Type(TIniFile), intent(in) :: Ini
  character(*), intent(IN) :: Key
  character(LEN=Ini_max_string_len) :: S
  double precision :: Ini_read_dbl_File 
   
  S = Ini_Read_String_File(Ini,Key)
  if(S=='') then
    write(*,*) 'no value for key: '//Key
    stop
  else
    read (S,*,err=10) Ini_read_dbl_File
  end if
  return

10 write (*,*) 'error reading double for key: '//Key
  stop

end function Ini_read_dbl_File


function Ini_Read_Real_File(Ini,Key)
  implicit none
  Type(TIniFile), intent(in) :: Ini
  real :: Ini_Read_Real_File 
  character(*), intent(IN) :: Key
  character(LEN=Ini_max_string_len) :: S

  S = Ini_Read_String_File(Ini,Key)
  if(S=='') then
    write(*,*) 'no value for key: '//Key
    stop
  else
    read (S,*,err=10) Ini_Read_Real_File
  end if
  return

10 write (*,*) 'error reading double for key: '//Key
  stop

end function Ini_Read_Real_File


function Ini_read_log_File(Ini,Key)
  implicit none
  Type(TIniFile), intent(in) :: Ini
  logical Ini_read_log_File
  character(*), intent(IN) :: Key
  character(LEN=Ini_max_string_len) :: S
   
  S = Ini_Read_String_File(Ini,Key)
  if(S=='') then
    write(*,*) 'no value for key: '//Key
    stop
  else
    if (verify(trim(S),'10TF') /= 0) goto 10  
    read (S,*,err=10) Ini_read_log_File
  end if
  return

10 write (*,*) 'error reading logical for key: '//Key
  stop

end function Ini_read_log_File


function Ini_Read_String_File(Ini,Key) result(f)
  implicit none
  Type(TIniFile), intent(in) :: Ini
  character(*), intent(IN) :: Key
  character(LEN=Ini_max_string_len) :: f

   call TNameValueList_ValueOf(Ini%L, Key, f)
    call  TNameValueList_Add(Ini%ReadValues, Key, f)

end function Ini_Read_String_File


subroutine TNameValueList_ValueOf(L,AName,AValue)
  implicit none
  Type(TNameValueList), intent(in) :: L
  character(*), intent(in) :: AName
  CHARACTER(*), intent(out) :: AValue
  integer :: i

  do i=1, L%Count
    if (L%Items(i)%P%Name == AName) then
      AValue = L%Items(i)%P%Value 
      return
    end if
  end do
  AValue = ''

end subroutine TNameValueList_ValueOf


function TNameValueList_HasKey(L,AName) result (AValue)
  implicit none
  Type (TNameValueList), intent(in) :: L
  character(*), intent(in) :: AName
  logical :: AValue
  integer i

  do i=1, L%Count
    if (L%Items(i)%P%Name==AName) then
      AValue = .true.
      return
    end if
  end do
  AValue = .false.
     
end function TNameValueList_HasKey
 

subroutine TNameValueList_Init(L,ignoreDuplicates)
  implicit none
  Type (TNameValueList) :: L
  logical, intent(in), optional :: ignoreDuplicates
    
  L%Count = 0
  L%Capacity = 0
  L%Delta = 128
  L%ignoreDuplicates = .false.
  if (present(ignoreDuplicates)) L%ignoreDuplicates=ignoreDuplicates
  nullify(L%Items)

end subroutine TNameValueList_Init


subroutine TNameValueList_Clear(L)
  implicit none
  Type (TNameValueList) :: L
  integer i, status
     
  do i = L%count, 1, -1
    deallocate(L%Items(i)%P, stat = status)
  end do
  deallocate(L%Items, stat = status)
  call TNameValueList_Init(L)

end subroutine TNameValueList_Clear

   
subroutine TNameValueList_Add(L,AName,AValue)
  implicit none
  Type(TNameValueList) :: L
  character(*), intent(in) :: AName, AValue

  if(TNameValueList_HasKey(L,AName)) then
    if (L%ignoreDuplicates) return
    write(*,*) 'IniFile,TNameValueList_Add: duplicate key name in file: '//trim(AName)
    stop
  end if
  if(L%Count==L%Capacity) call TNameValueList_SetCapacity(L,L%Capacity+L%Delta)
  L%Count = L%Count + 1
  allocate(L%Items(L%Count)%P)
  L%Items(L%Count)%P%Name = AName
  L%Items(L%Count)%P%Value = AValue

end subroutine TNameValueList_Add

subroutine TNameValueList_SetCapacity(L,C)
  implicit none
  Type (TNameValueList) :: L
  integer :: C
  type(TNameValue_pointer), dimension(:), pointer :: TmpItems
    
  if (L%Count > 0) then
    if (C < L%Count) stop 'TNameValueList_SetCapacity: smaller than Count'
    allocate(TmpItems(L%Count))
    TmpItems = L%Items(1:L%Count)
    deallocate(L%Items)
    allocate(L%Items(C))
    L%Items(1:L%Count) = TmpItems
    deallocate(TmpItems)
  else
    allocate(L%Items(C))
  end if  
  L%Capacity = C
  
end subroutine TNameValueList_SetCapacity

subroutine TNameValueList_Delete(L, i)
  implicit none
  Type (TNameValueList) :: L
  integer, intent(in) :: i
     
  deallocate(L%Items(i)%P)
  if (L%Count > 1) L%Items(i:L%Count-1) = L%Items(i+1:L%Count)
  L%Count = L%Count -1
     
end subroutine TNameValueList_Delete

subroutine Ini_NameValue_Add(Ini,AInLine)
  implicit none
  Type(TIniFile) :: Ini
  character(*), intent(IN) :: AInLine
  integer EqPos, slashpos, lastpos
  character(LEN=len(AInLine)) :: AName, S, InLine

  InLine = trim(adjustl(AInLine))
  EqPos = scan(InLine,'=')
  if(EqPos/=0.and.InLine(1:1)/='#') then
    AName = trim(InLine(1:EqPos-1))
    S = adjustl(InLine(EqPos+1:)) 
    if(Ini%SlashComments) then
      slashpos = scan(S,'/')
      if(slashpos/=0) S = S(1:slashpos-1)
    end if
    lastpos=len_trim(S)
    if(lastpos>1) then
      if(S(1:1)==''''.and.S(lastpos:lastpos)=='''') S = S(2:lastpos-1)
    end if
    call TNameValueList_Add(Ini%L,AName,S)
  end if

end subroutine Ini_NameValue_Add


subroutine Ini_Open(filename,unit_id,error,slash_comments)
  implicit none
  character(*), intent(IN) :: filename
  integer, intent(IN) :: unit_id
  logical, optional, intent(OUT) :: error
  logical, optional, intent(IN) :: slash_comments
  logical :: aerror

  call TNameValueList_Init(DefIni%L)
  call TNameValueList_Init(DefIni%ReadValues, .true.)
          
  if (present(slash_comments)) then
    call Ini_Open_File(DefIni,filename,unit_id,aerror,slash_comments)
  else
    call Ini_Open_File(DefIni,filename,unit_id,aerror)
  end if

  if (present(error)) then
    error = aerror
  else
    if (aerror) then
      write (*,*) 'Ini_Open: Error opening file ' // trim(filename)
      stop
    end if
  end if

end subroutine Ini_Open

function Ini_ExtractFilePath(aname)
  implicit none
  character(*), intent(IN) :: aname
  character(LEN=Ini_max_string_len) Ini_ExtractFilePath
  integer :: l, i

  l = len_trim(aname)
  do i = l, 1, -1
    if (aname(i:i)=='/') then
      Ini_ExtractFilePath = aname(1:i)
      return
    end if
  end do
  Ini_ExtractFilePath = ''

end function Ini_ExtractFilePath

recursive subroutine Ini_Open_File(Ini,filename,unit_id,error,slash_comments,append)
  implicit none
  Type(TIniFile) :: Ini
  character(*), intent(IN) :: filename
  integer, intent(IN) :: unit_id
  logical, intent(OUT) :: error
  logical, optional, intent(IN) :: slash_comments
  logical, optional, intent(in) :: append
  character (LEN=Ini_max_string_len) :: InLine, IncludeFile
  integer :: lastpos, i
  Type (TNameValueList) IncudeFiles
  logical :: doappend, FileExists
     
  if (present(append)) then
    doappend=append
  else
    doappend=.false.
  end if  
     
  if(.not.doappend) then
    call TNameValueList_Init(Ini%L)
    call TNameValueList_Init(Ini%ReadValues, .true.)
  end if
 
  call TNameValueList_Init(IncudeFiles) 
     
  if (present(slash_comments)) then
    Ini%SlashComments = slash_comments
  else
    Ini%SlashComments = .false.
  end if
         
  open(unit=unit_id,file=filename,form='formatted',status='old',err=500) 
  do 
    read (unit_id,'(a)',end=400) InLine
    if (InLine == 'END') exit;
    if (InLine(1:8) == 'INCLUDE(') then
      lastpos = scan(InLine,')')
      if (lastpos/=0) then
        call TNameValueList_Add(IncudeFiles,trim(adjustl(InLine(9:lastpos-1))),'')            
      else
        stop 'Ini_Open_File: error in INCLUDE line'
      end if 
    else if (InLine /= '') then
      call Ini_NameValue_Add(Ini,InLine) 
    end if
  end do

400 close(unit_id)
  error=.false.

  do i=1, IncudeFiles%Count
    if(error) exit
    IncludeFile=IncudeFiles%Items(i)%P%Name
    inquire(file=IncludeFile, exist = FileExists)
    if (.not. FileExists) then
      IncludeFile=trim(Ini_ExtractFilePath(filename))//trim(IncludeFile)
      inquire(file=IncludeFile, exist = FileExists)
      if (.not. FileExists) stop 'Ini_Open_File: INCLUDE file not found'
    end if
    call Ini_Open_File(Ini,IncludeFile,unit_id,error,slash_comments,append=.true.)      
  end do
  call TNameValueList_Clear(IncudeFiles)  
  return

500 error=.true.
  call TNameValueList_Clear(IncudeFiles)

end subroutine Ini_Open_File

subroutine Ini_Open_Fromlines(Ini, Lines, NumLines, slash_comments)
  implicit none
  Type(TIniFile) :: Ini
  integer, intent(IN) :: NumLines
  character(*), dimension(NumLines), intent(IN) :: Lines
  logical, intent(IN) :: slash_comments
  integer i

  call TNameValueList_Init(Ini%L)
  call TNameValueList_Init(Ini%ReadValues, .true.)

  Ini%SlashComments = slash_comments

  do i=1,NumLines
    call Ini_NameValue_Add(Ini,Lines(i))
  end do

end subroutine Ini_Open_Fromlines

subroutine Ini_Close
  implicit none
  Type(TIniFile) :: Ini

  call TNameValueList_Clear(Ini%L)
  call TNameValueList_Clear(Ini%ReadValues)

end subroutine Ini_Close

end module readfile

