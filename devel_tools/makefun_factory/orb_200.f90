!     THE  COSTANT

indorbp=indorb+1
indshellp=indshell+1

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=1.d0
end do
!           endif

if(typec.ne.1) then
  do i=1,3
    z(indorbp,indt+i)=0
  end do

  z(indorbp,indt+4)=0
  !endif for indt
end if

indshell=indshellp
indorb=indorbp


