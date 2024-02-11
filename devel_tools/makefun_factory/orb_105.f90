!     2s double gaussian without  constant
!     (exp (-dd2 r^2)+dd4*exp(-dd5*r^2))



!        dd1=1.d0
dd2=dd(indpar+1)
!        dd3=dd(indpar+2)
!        dd4=dd(indpar+3)
!        dd5=dd(indpar+4)
dd4=dd(indpar+2)
dd5=dd(indpar+3)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k)*r(k))
  distp(k,2)=dexp(-dd5*r(k)*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)+dd4*distp(i,2)
!              write(6,*) ' function inside = ',z(indorbp,i)
end do
!           endif


if(typec.ne.1) then
  fun=-2.d0*(dd2*distp(0,1)+dd5*dd4*distp(0,2))
  fun2=r(0)**2

  !              write(6,*) ' fun inside = ',fun,fun2

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*(dd2*(-3.d0+2.d0*dd2*fun2)*                 &
      distp(0,1)+dd5*dd4*(-3.d0+2.d0*dd5*fun2)*distp(0,2))




  !endif for indt
end if

indpar=indpar+3
indshell=indshellp
indorb=indorbp

