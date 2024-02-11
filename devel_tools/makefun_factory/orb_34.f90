  ! normalized
  ! exp(-dd1*r) + dd1*r*exp(-dd1*r)


  indshellp=indshell+1


  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  !           peff=dd1


  !           if(iflagnorm.gt.2) then
  !              c=1.d0/dsqrt(1.d0/4.d0/dd1**3+12.d0*dd1/(2.d0*dd1)**4+&
      !                                                              &3.d0*dd1**2/4.d0/dd1**5)/dsqrt(4.d0*pi)
      !             c=dd1*dsqrt(dd1)/dsqrt(7.d0*pi)
  c=dd1*dsqrt(dd1)*.2132436186229231d0
  !           endif

  do i=indtmin,indtm
    distp(i,1)=c*dexp(-dd1*r(i))
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)*(1.d0+r(i)*dd1)
  end do

  if(typec.ne.1) then

    fun=-dd1**2*distp(0,1)
    fun2=fun*(1.d0-dd1*r(0))

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=(2.d0*fun+fun2)

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+1
  indshell=indshellp




  ! 2s single  Z WITH CUSP
