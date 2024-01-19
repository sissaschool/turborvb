!     2s with cusp condition
!     ( r^3*exp(-dd2*r^2))  ! with no cusp condition

dd2=dd(indpar+1)
indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k)**2)*r(k)
end do
!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)*r(i)**2
end do
!           endif
if(typec.ne.1) then
  rp1=r(0)**2*dd2
  fun=(3.d0-2.d0*rp1)*distp(0,1)
  fun2=(6.d0-14.d0*rp1+4.d0*rp1**2)*distp(0,1)
  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do
  z(indorbp,indt+4)=2.d0*fun+fun2
  !endif for indt
end if
indpar=indpar+1
indshell=indshellp
indorb=indorbp
