  ! R(r)=c*exp(-z*r**2)*(3/4/z-r**2)

  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1

  dd1=dd(indpar+1)
  dd2=dsqrt(dd1)

  !         c=(2.d0*dd1/pi)**(3.d0/4.d0)*ratiocs
  c=dd1**0.75d0*ratiocs

  do k=indtmin,indtm
    cost=dd1*r(k)**2/(1.d0+dd2*r(k))
    distp(k,1)=c*dexp(-cost)
  end do

  do i=i0,indtm
    cost=(1.d0+0.5d0*dd2*r(i))/(1.d0+dd2*r(i))**2
    z(indorbp,i)=distp(i,1)*(3.d0/4.d0/dd1-r(i)**2*cost)
  end do

  if(typec.ne.1) then
    !              the first derivative /r

    rp1=dd1*r(0)**2
    rp2=dd2*r(0)
    rp3=(1.d0+rp2)**2

    fun=0.25d0*distp(0,1)*                                           &
        (-14.d0-29.d0*rp2-12.d0*rp1+3.d0*rp1*rp2+2.d0*rp1**2)/rp3**2
    !              the second derivative
    fun2=0.25d0*distp(0,1)*                                          &
        (-14.d0-30.d0*rp2+34.d0*rp1+118.d0*rp1*rp2+87.d0*rp1**2      &
        +18.d0*rp1**2*rp2-5.d0*rp1**3-2.d0*rp1**3*rp2)/rp3**3

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  indpar=indpar+1
  indshell=indshellp

