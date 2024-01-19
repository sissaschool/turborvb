  ! R(r)=r**3*exp(-z*r**2) single zeta



  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)

  !           if(iflagnorm.gt.2) then
  !     c=2.d0/pi**(3.d0/4.d0)*(2.d0*dd1)**(9.d0/4.d0)*dsqrt(2.d0/105.d0)
  c=dd1**2.25d0*.55642345640820284397d0
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)*r(k)
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)*r(i)**2
  end do

  if(typec.ne.1) then
    rp1=r(0)**2*dd1
    !              the first derivative / r
    fun=distp(0,1)*(3.d0-2.d0*rp1)
    !              the second derivative
    fun2=distp(0,1)*(6.d0-14.d0*rp1+4.d0*rp1**2)
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  !        endif

  indpar=indpar+1
  indshell=indshellp

  ! 3s -derivative of 60 with respect to dd1
