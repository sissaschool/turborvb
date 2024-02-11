  ! R(r)=r**4*exp(-z*r**2) single zeta



  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)

  !           if(iflagnorm.gt.2) then
  !           c=(2.d0*dd1**11/pi)**(1.d0/4.d0)*(512.d0/945.d0/pi)
  c=dd1**2.75d0*0.1540487967684377d0
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  do i=i0,indtm
    z(indorbp,i)=r(i)**4*distp(i,1)
  end do

  if(typec.ne.1) then
    rp1=r(0)**2

    !              the first derivative
    fun=distp(0,1)*rp1*(4.d0-2.d0*dd1*rp1)

    !              the second derivative
    fun2=distp(0,1)*rp1*(12.d0-18.d0*dd1*rp1                         &
        +4.d0*dd1**2*rp1**2)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*fun+fun2

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+1
  indshell=indshellp

  ! derivative of 16 with respect to z
