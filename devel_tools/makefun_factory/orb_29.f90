  ! derivative of (28)


  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1

  dd1=dd(indpar+1)

  !           if(iflagnorm.gt.2) then
  !             if(dd1.gt.0.) then
  c=cost1s*dd1**1.5d0
  !             else
  !               c=1.d0
  !             endif
  !           endif
  !             if(dd1.gt.0.) then
  c1=1.5d0/dd1
  !             else
  !               c1=0.d0
  !             endif

  do i=indtmin,indtm
    distp(i,1)=c*dexp(-dd1*r(i))
  end do

  do i=i0,indtm
    !              rp1=(b1s*r(i))**4*dd1**3
    !              rp4=rp1*dd1
    !              rp5=dd1*r(i)
    !              z(indorbp,i)=distp(i,1)*                          &
        !                                                            &(c1*rp4/(1+rp4)-rp1*(-4+rp5+rp4*rp5)/(1+rp4)**2)
        rp4=(b1s*dd1*r(i))**4
    rp5=dd1*r(i)
    z(indorbp,i)=distp(i,1)*rp4/(1+rp4)*                             &
        (c1 - (1.d0/dd1)*(-4+rp5+rp4*rp5)/(1+rp4))
  end do

  if(typec.ne.1) then
    rp1=dd1*b1s*r(0)
    rp2=rp1**2
    rp4=rp2**2
    rp5=rp4*rp1
    rp8=rp4*rp4

    fun=distp(0,1)* (dd1*rp2*(4*b1s**2*(11-5*rp4) +2*(rp1+rp5)**2    &
        -b1s*rp1*(21+26*rp4+5*rp8)))/(2.*(1+rp4)**3)

    fun2=distp(0,1)*(dd1*rp2*(b1s*(31 + 7*rp4)*(rp1 + rp5)**2        &
        - 2*(rp1 + rp5)**3 + 64*b1s**2*rp1*(-2 - rp4 + rp8) +        &
        4*b1s**3*(33 - 134*rp4 + 25*rp8)))/(2.*b1s*(1 + rp4)**4)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*fun+fun2

    !           endif

    indorb=indorbp

  end if

  indpar=indpar+1
  indshell=indshellp





