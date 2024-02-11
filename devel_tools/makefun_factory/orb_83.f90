  ! R(r)=x*exp(-z*r**2)*(5/4/z-r**2)

  dd1=dd(indpar+1)
  dd2=dsqrt(dd1)
  c=dd1**1.25d0*ratiocp

  do k=indtmin,indtm
    cost=dd1*r(k)**2/(1.d0+dd2*r(k))
    distp(k,1)=c*dexp(-cost)
  end do
  !        indorbp=indorb
  !
  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      cost=(1.d0+0.5d0*dd2*r(i))/(1.d0+dd2*r(i))**2
      z(indorbp,i)=rmu(ic,i)*distp(i,1)*(1.25d0/dd1-r(i)**2*cost)
    end do
    !           endif
  end do

  if(typec.ne.1) then
    rp1=dd1*r(0)**2
    rp2=dd2*r(0)
    rp3=(1.d0+rp2)**2
    cost=(1.d0+0.5d0*rp2)/rp3
    fun0=distp(0,1)*(1.25d0/dd1-r(0)**2*cost)
    fun=0.25d0*distp(0,1)*                                           &
        (-18.d0-39.d0*rp2-20.d0*rp1+rp1*rp2+2.d0*rp1**2)/rp3**2
    fun2=0.25d0*distp(0,1)*                                          &
        (-18.d0-42.d0*rp2+30.d0*rp1+138.d0*rp1*rp2+113.d0*rp1**2     &
        +30.d0*rp1**2*rp2-3.d0*rp1**3-2.d0*rp1**3*rp2)/rp3**3

    !              indorbp=indorb

    do ic=1,3
      !                if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
      !                 endif
    end do
    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp



