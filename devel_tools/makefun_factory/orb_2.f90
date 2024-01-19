  !



  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1

  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=(zeta(1)-dd1)/(dd2-zeta(1))

  !           if(iflagnorm.gt.2) then
  c=1.d0/2.d0/dsqrt(2.d0*pi*(1.d0/(2.d0*dd1)**3                      &
      +2.d0*peff/(dd1+dd2)**3+peff**2/(2.d0*dd2)**3))
  !           endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=c*dexp(-dd2*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)+peff*distp(i,2)
  end do

  if(typec.ne.1) then
    fun=(-dd1*distp(0,1)-dd2*distp(0,2)*peff)/r(0)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=(-2.d0*dd1/r(0)+dd1**2)                        &
        *distp(0,1)+peff*(-2.d0*dd2/r(0)+dd2**2)                     &
        *distp(0,2)


  end if

  indorb=indorbp

  !        endif

  indpar=indpar+2
  indshell=indshellp


  ! 1s double Z NO CUSP
