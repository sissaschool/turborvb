!     2s double lorentian with constant  parent of 102
!     (dd3+r^3/(1+dd5*r)^4;



dd3=dd(indpar+1)
dd5=dd(indpar+2)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=r(k)**3/(1.d0+dd5*r(k))**4
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=dd3+distp(i,1)
!              write(6,*) ' function inside = ',z(indorbp,i)
end do
!           endif


if(typec.ne.1) then

  fun=                                                               &
      -r(0)*(-3.d0+dd5*r(0))/(1.d0+dd5*r(0))**5
  fun2=                                                              &
      +2.d0*r(0)*(3.d0-6.d0*dd5*r(0)+(dd5*r(0))**2)                  &
      /(1.d0+dd5*r(0))**6


  !              write(6,*) ' fun inside = ',fun,fun2

  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=fun2+2.d0*fun

  !           write(6,*) ' lap 106 =',z(indorbp,indt+4)

  !endif for indt
end if

indpar=indpar+2
indshell=indshellp
indorb=indorbp

