  ! normalized

  indshellp=indshell+1


  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  dd2=dd(indpar+2)

  !           if(iflagnorm.gt.2) then
  c=dd1**2.5d0/dsqrt(3.d0*pi*(1.d0+dd2**2/3.d0))
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
  end do


  do i=i0,indtm
    z(indorbp,i)=(r(i)+dd2*rmu(3,i))*distp(i,1)
  end do

  if(typec.ne.1) then

    fun=distp(0,1)*(1.d0-dd1*r(0))
    funp=-dd2*dd1*distp(0,1)*rmu(3,0)
    fun2=distp(0,1)*(dd1**2*r(0)-2.d0*dd1)
    fun2p=dd1**2*dd2*distp(0,1)*rmu(3,0)

    do i=1,3
      z(indorbp,indt+i)=(fun+funp)*rmu(i,0)/r(0)
    end do
    z(indorbp,indt+3)=z(indorbp,indt+3)+dd2*distp(0,1)
    z(indorbp,indt+4)=(2.d0*fun+4.d0*funp)/r(0)                      &
        +(fun2+fun2p)

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+2
  indshell=indshellp


  ! 2s single Z NO CUSP
