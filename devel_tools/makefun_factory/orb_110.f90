!     2s without cusp condition
!     dd1*( dd3 +1/(1+dd2*r^3))


dd2=dd(indpar+1)
dd3=dd(indpar+2)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=1.d0/(1.d0+dd2*r(k)**3)
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)+dd3
end do
!           endif


if(typec.ne.1) then
  fun=-dd2*distp(0,1)**2*3.d0*r(0)
  fun2=fun*distp(0,1)*(2.d0-4.d0*dd2*r(0)**3)

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=fun2+2.d0*fun


  !endif for indt
end if

indpar=indpar+2
indshell=indshellp
indorb=indorbp

