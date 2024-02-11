  ! R(r)=c*r**2*exp(-z*r**2)*(7/4/z-r**2)


  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)

  !           if(iflagnorm.gt.2) then
  c=4.d0*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)/dsqrt(15.d0)
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)*(7.d0/4.d0/dd1*r(i)**2                   &
        -r(i)**4)
  end do

  if(typec.ne.1) then
    rp1=r(0)**2
    !              the first derivative / r
    fun=distp(0,1)*(7.d0-15.d0*dd1*rp1                               &
        +4.d0*(dd1*rp1)**2)/2.d0/dd1
    !              the second derivative
    fun2=distp(0,1)*(7.d0-59*dd1*rp1+50*(dd1*rp1)**2                 &
        -8*(dd1*rp1)**3)/2.d0/dd1
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  !        endif

  indpar=indpar+1
  indshell=indshellp

  ! 5s single zeta derivative of 12
