  ! normalized IS WRONG!!!

  indshellp=indshell+1


  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd(indpar+3)

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
    distp(k,2)=dexp(-dd2*r(k))
  end do
  !           if(iflagnorm.gt.2) then
  c=                                                                 &
      1/dsqrt(1/(3.D0/4.D0/dd1**5+peff**2/dd2**3/4+12*peff/          &
      (dd1+dd2)**4))*1.d0/dsqrt(4.0*pi)
  !           endif

  do i=i0,indtm
    z(indorbp,i)=c*(distp(i,1)+r(i)*distp(i,2)*peff)
  end do

  if(typec.ne.1) then

    fun=-dd1*distp(0,1)+peff*distp(0,2)*(1.d0-dd2*r(0))
    fun2=distp(0,1)*dd1**2                                           &
        +peff*distp(0,2)*(dd2**2*r(0)-2.d0*dd2)

    do i=1,3
      z(indorbp,indt+i)=fun*c*rmu(i,0)/r(0)
    end do
    z(indorbp,indt+4)=c*(2.d0*fun/r(0)+fun2)

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+3
  indshell=indshellp



  ! 2s double Z WITH CUSP
