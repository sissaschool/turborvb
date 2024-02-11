!     2s double lorentian with constant  parent of 102
!     (dd3+ L(dd2 r^2)+dd4*L(dd5*r^2)) ;  L(x)=1/1+x^2



dd2=dd(indpar+1)
dd3=dd(indpar+2)
dd4=dd(indpar+3)
dd5=dd(indpar+4)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=1.d0/(1.d0+dd2*r(k)*r(k))
  distp(k,2)=1.d0/(1.d0+dd5*r(k)*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=(distp(i,1)+dd3+dd4*distp(i,2))
!              write(6,*) ' function inside = ',z(indorbp,i)
end do
!           endif


if(typec.ne.1) then
  fun=-2.d0*(dd2*distp(0,1)**2+dd5*dd4*distp(0,2)**2)
  fun2=2.d0*dd2*distp(0,1)**3*(-1.d0+3.d0*dd2*r(0)**2)               &
      +2.d0*dd5*dd4*distp(0,2)**3*(-1.d0+3.d0*dd5*r(0)**2)

  !              write(6,*) ' fun inside = ',fun,fun2

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=fun2+2.d0*fun

  !           write(6,*) ' lap 106 =',z(indorbp,indt+4)

  !endif for indt
end if

indpar=indpar+4
indshell=indshellp
indorb=indorbp


