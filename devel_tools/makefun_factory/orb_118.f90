!     2s double lorentian with constant  parent of 102
!     (dd1+  1/ (1 + Exp[  dd2 (r^2 - r_0^2) ] )   | dd3=r_0
!      Fermi distribution with r^2


dd1=dd(indpar+1)
dd2=dd(indpar+2)
dd3=-dd2*dd(indpar+3)**2

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  arg=dd2*r(k)**2+dd3
  if(arg.gt.200) then
    distp(k,1)=dexp(200.d0)
  else
    distp(k,1)=dexp(arg)
  end if
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=dd1+1.d0/(1.d0+distp(i,1))
!              write(6,*) ' function inside = ',z(indorbp,i)
end do
!           endif


if(typec.ne.1) then

  fun= -2.d0*dd2*distp(0,1)/(1.d0+distp(0,1))**2
  fun2=-2.d0*dd2*(-distp(0,1)*(-1.d0-2.d0*dd2*r(0)**2)               &
      +distp(0,1)**2*(1.d0-2.d0*dd2*r(0)**2))/(1.d0+distp(0,1))**3


  !              write(6,*) ' fun inside = ',fun,fun2

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=fun2+2.d0*fun

  !           write(6,*) ' lap 106 =',z(indorbp,indt+4)

  !endif for indt
end if

indpar=indpar+3
indshell=indshellp
indorb=indorbp

