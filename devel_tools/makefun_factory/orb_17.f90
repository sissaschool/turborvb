  ! 2s gaussian for pseudo
  ! R(r)=r**2*exp(-z*r**2) single zeta

  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)

  !           if(iflagnorm.gt.2) then
  !        c=4.d0*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)/dsqrt(15.d0)
  c=.73607904464954686606d0*dd1**1.75d0
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)*r(i)**2
  end do

  if(typec.ne.1) then
    rp1=r(0)**2
    !              the first derivative / r
    fun=2.d0*distp(0,1)*(1.d0-dd1*rp1)
    !              the second derivative
    fun2=2.d0*distp(0,1)*(1.d0-5.d0*dd1*rp1+2.d0*dd1**2*rp1**2)
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  !        endif

  indpar=indpar+1
  indshell=indshellp

  ! 2s gaussian for pseudo
