  ! R(r)=(Z*b*r)^4/(1+(Z*b*r)^4)*exp(-z*r)     orbital 1s (no cusp)
  ! d -> b1s (defined in module constants)
  ! normadization: cost1s, depends on b1s


  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1

  dd1=dd(indpar+1)

  !           if(iflagnorm.gt.2) then
  !             if(dd1.gt.0.) then
  c=cost1s*dd1**1.5d0
  !             else
  !               c=1.d0
  !             endif
  !           endif

  do i=indtmin,indtm
    distp(i,1)=c*dexp(-dd1*r(i))
  end do

  do i=i0,indtm
    rp4=(dd1*b1s*r(i))**4
    z(indorbp,i)=distp(i,1)*rp4/(1.d0+rp4)
  end do

  if(typec.ne.1) then
    rp1=dd1*b1s*r(0)
    rp2=rp1**2
    rp4=rp2**2
    rp5=r(0)*dd1
    rp6=(b1s*dd1)**2*rp2
    !              the first derivative /r
    fun=-distp(0,1)*rp6*(-4.d0+rp5+rp4*rp5)/(1.d0+rp4)**2
    !              the second derivative derivative
    fun2=distp(0,1)*rp6*(12.d0-8*rp5+rp5**2-20*rp4-                  &
        8*rp4*rp5+2*rp4*rp5**2+(rp4*rp5)**2)/(1.d0+rp4)**3
    !    gradient: dR(r)/dr_i=r_i*fun
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    !    laplacian = 2*fun+fun2
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  !        endif

  indpar=indpar+1
  indshell=indshellp

