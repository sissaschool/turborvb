  ! (r**2 + dd2*(1 + dd1*r))*exp(-dd1*r)  ! normalized


  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  dd2=dd(indpar+2)

  c=dsqrt(2.d0*dd1**7/pi/                                            &
      (45.d0+42.d0*dd1**2*dd2+14.d0*dd1**4*dd2**2))

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
  end do


  do i=i0,indtm
    z(indorbp,i)=(r(i)**2+dd2*(1.d0+dd1*r(i)))                       &
        *distp(i,1)
  end do

  if(typec.ne.1) then

    fun=distp(0,1)*r(0)*(2.d0-dd1**2*dd2-dd1*r(0))
    fun2=distp(0,1)*((1.d0-dd1*r(0))                                 &
        *(3.d0-dd1**2*dd2-dd1*r(0))-1.d0)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)/r(0)
    end do
    z(indorbp,indt+4)=2.d0*fun/r(0)+fun2

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+2
  indshell=indshellp

  ! 2s gaussian for pseudo
