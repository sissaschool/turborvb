  ! R(r)=c (-r**5*exp(-z1*r**2)+c1 r**3 exp(-z1*r**2))


  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)

  !          if(iflagnorm.gt.2) then
  !     c=2.d0/pi**(3.d0/4.d0)*(2.d0*dd1)**(9.d0/4.d0)*dsqrt(2.d0/105.d0)
  c=dd1**2.25d0*.55642345640820284397d0
  !           endif

  c1=2.25d0/dd1

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)*r(k)
  end do

  do i=i0,indtm
    z(indorbp,i)=(-r(i)**4+c1*r(i)**2)*distp(i,1)
  end do


  if(typec.ne.1) then
    rp1=r(0)**2
    rp2=rp1*dd1

    fun=c1*distp(0,1)*(3.d0-2.d0*rp2)                                &
        +distp(0,1)*rp1*(-5.d0+2.d0*rp2)
    !              the second derivative
    fun2=c1*distp(0,1)*(6.d0-14.d0*rp2+4.d0*rp2**2)                  &
        +distp(0,1)*rp1*(-20.d0+22.d0*rp2-4.d0*rp2**2)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  !        endif
  indpar=indpar+1
  indshell=indshellp

  ! single gaussianx r  p orbitals
