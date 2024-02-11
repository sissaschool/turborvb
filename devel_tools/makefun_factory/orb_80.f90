  ! R(r)=exp(-z*r**2) single zeta

  indshellp=indshell+1
  indorbp=indorb+1

  dd1=dd(indpar+1)
  dd2=dsqrt(dd1)

  !           if(iflagnorm.gt.2) then
  !           c=(2.d0*dd1/pi)**(3.d0/4.d0)*ratiocs
  !           ratiocs--> ratiocs*(2/pi)**3/4
  c=dd1**0.75d0*ratiocs
  !           endif

  do k=indtmin,indtm
    cost=dd1*r(k)**2/(1.d0+dd2*r(k))
    distp(k,1)=c*dexp(-cost)
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)
  end do

  if(typec.ne.1) then
    ! the first derivative /r
    rp1=dd1*r(0)**2
    rp2=dd2*r(0)
    rp3=(1.d0+rp2)**2

    fun=-dd1*distp(0,1)*(2.d0+rp2)/rp3

    ! the second derivative
    fun2=dd1*distp(0,1)*(-2.d0-2.d0*rp2+4.d0*rp1+4.d0*rp1*rp2+rp1**2)/rp3**2

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp
  indpar=indpar+1
  indshell=indshellp

