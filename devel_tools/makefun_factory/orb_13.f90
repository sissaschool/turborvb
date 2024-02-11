  ! R(r)=r**3*(exp(-z1*r)+z3*exp(-z2*r))
  !
  indshellp=indshell+1

  !
  !
  !        if(iocc(indshellp).eq.1) then
  !
  indorbp=indorb+1
  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  dd3=dd(indpar+3)
  !           if(iflagnorm.gt.2) then
  c=1.d0/2.d0/dsqrt(pi*40320.d0*(1.d0/(2.d0*dd1)**9+                 &
      2.d0*dd3/(dd1+dd2)**9+dd3**2/(2.d0*dd2)**9))
  !           endif

  !
  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=c*dexp(-dd2*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=(distp(i,1)+dd3*distp(i,2))*r(i)**3
  end do
  !
  if(typec.ne.1) then
    rp1=r(0)**3
    rp2=r(0)**2
    !
    !c              the first derivative
    fun=distp(0,1)*(3.d0*rp2-dd1*rp1)                                &
        +dd3*distp(0,2)*(3.d0*rp2-dd2*rp1)
    !c
    !              the second derivative
    fun2=distp(0,1)*(6.d0*r(0)-6.d0*dd1*rp2+dd1**2*rp1)              &
        +dd3*distp(0,2)*(6.d0*r(0)-6.d0*dd2*rp2+dd2**2*rp1)
    !c
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)/r(0)
    end do
    !
    z(indorbp,indt+4)=2.d0*fun/r(0)+fun2
    !
  end if

  indorb=indorbp

  !        endif
  indpar=indpar+3
  indshell=indshellp

  ! 1s single Z pseudo
