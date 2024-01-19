  ! R(r)=r**2*exp(-z1*r)



  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  !           if(iflagnorm.gt.2) then
  !              c=1.d0/dsqrt(1.d0/4.d0/dd1**3+12.d0*dd1/(2.d0*dd1)**4+&
      !                                                              &3.d0*dd1**2/4.d0/dd1**5)/dsqrt(4.d0*pi)
      !             c=dd1*dsqrt(dd1)/dsqrt(7.d0*pi)
  c=dd1*dsqrt(dd1)*0.21324361862292308211d0
  !           endif

  c0=-c*dd1

  c1=1.5d0*c/dd1



  do i=indtmin,indtm
    distp(i,1)=dexp(-dd1*r(i))
  end do

  do i=i0,indtm
    z(indorbp,i)=(c0*r(i)**2+c1*(1.d0+dd1*r(i)))                     &
        *distp(i,1)
  end do

  c1=c1*dd1**2

  if(typec.ne.1) then
    fun=(c0*(2.d0-dd1*r(0))-c1)*distp(0,1)
    fun2=(c0*(2.d0-4*dd1*r(0)+(dd1*r(0))**2)                         &
        +c1*(dd1*r(0)-1.d0))*distp(0,1)
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do
    z(indorbp,indt+4)=2.d0*fun+fun2
  end if

  indorb=indorbp

  !        endif
  indpar=indpar+1
  indshell=indshellp

  ! 4s single zeta derivative of 10
