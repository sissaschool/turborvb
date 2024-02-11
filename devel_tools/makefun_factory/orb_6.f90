  ! normalized

  indshellp=indshell+1


  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd(indpar+3)

  !           if(iflagnorm.gt.2) then
  !            c=                                              WRONG
  !                                                                  &0.5*dsqrt((1.d0/dd1**5/32.d0+2.d0*peff/(dd1+dd2)**5
      !                                                              &+peff**2/dd2**5/32.d0)/(24.d0*pi))

  c=1.d0/dsqrt((3.d0*pi)*                                            &
      (1.d0/dd1**5+ 64.d0*peff/(dd1+dd2)**5+peff**2/dd2**5))

  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=c*dexp(-dd2*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=r(i)*(distp(i,1)+distp(i,2)*peff)
  end do

  if(typec.ne.1) then

    fun=distp(0,1)*(1.d0-dd1*r(0))                                   &
        +peff*distp(0,2)*(1.d0-dd2*r(0))
    fun2=distp(0,1)*(dd1**2*r(0)-2.d0*dd1)+peff*distp(0,2)           &
        *(dd2**2*r(0)-2.d0*dd2)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)/r(0)
    end do
    z(indorbp,indt+4)=(2.d0*fun/r(0)+fun2)

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+3
  indshell=indshellp

  ! 2s double Z NO CUSP
