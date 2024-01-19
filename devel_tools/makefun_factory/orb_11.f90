  ! R(r)=r**2*(exp(-z1*r)+p*exp(-z2*r))



  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd(indpar+3)
  !           if(iflagnorm.gt.2) then
  c=1.d0/2.d0/dsqrt(pi*720.d0*(1.d0/(2.d0*dd1)**7+                   &
      2.d0*peff/(dd1+dd2)**7+peff**2/(2.d0*dd2)**7))
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=c*dexp(-dd2*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=(distp(i,1)+peff*distp(i,2))*r(i)**2
  end do

  if(typec.ne.1) then
    rp1=r(0)**2
    !              the first derivative
    fun=distp(0,1)*(2.d0*r(0)-dd1*rp1)                               &
        +peff*distp(0,2)*(2.d0*r(0)-dd2*rp1)
    !
    !              the second derivative
    fun2=distp(0,1)*(2.d0-4.d0*dd1*r(0)+dd1**2*rp1)                  &
        +peff*distp(0,2)*(2.d0-4.d0*dd2*r(0)+dd2**2*rp1)
    !
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)/r(0)
    end do

    z(indorbp,indt+4)=2.d0*fun/r(0)+fun2

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+3
  indshell=indshellp


  ! 4s single zeta
