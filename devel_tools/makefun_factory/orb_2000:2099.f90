!     s gaussian  -r**(2*(npower+1))*exp(-alpha*r**2) derivative of 1000

npower=iopt+1-2000

indorbp=indorb+1
indshellp=indshell+1


dd2=dd(indpar+1)
do k=indtmin,indtm
  distp(k,1)=-r(k)**(2*npower)*dexp(-dd2*r(k)**2)
end do

!        if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)
end do
!        endif


if(typec.ne.1) then


  rp1=r(0)**2
  fun0=distp(0,1)
  fun=(npower-dd2*rp1)*distp(0,1)*2.d0/rp1
  fun2=(npower*(2.d0*npower-1.d0)-                                   &
      (1.d0+4.d0*npower)*dd2*rp1+2.d0*(dd2*rp1)**2)*                 &
      distp(0,1)*2.d0/rp1

  !           if(iocc(indshellp).eq.1) then
  do i=1,3
    z(indorbp,indt+i)=rmu(i,0)*fun
  end do
  z(indorbp,indt+4)=2.d0*fun+fun2
  !           endif


  !endif for indt
end if

indpar=indpar+1
indshell=indshell+1
indorb=indorbp


