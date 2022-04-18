!//////////////////////////////////////////////
! * SNeIa
!//////////////////////////////////////////////

module utils_SNeIa
  implicit none

contains


function nz_SNeIa(zs,zbin)  result(f)
  implicit none
  double precision, intent(in) :: zs, zbin(2)
  double precision :: f

  if(zbin(1)<zs.and.zbin(2)>zs) then 
    f = 1d0
  else 
    f = 0d0
  end if 

end function nz_SNeIa


end module utils_SNeIa

