  ! R(r)=r**3*exp(-z1*r)
  !
  indshellp=indshell+1



  !       if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)
  !           if(iflagnorm.gt.2) then
  !           c=dsqrt((2*dd1)**9/40320.d0/pi)/2.d0
  !           c=dsqrt((2*dd1)**7/720.d0/pi)/2.d0
  c=dd1**3.5d0*0.11894160774351807429d0
  !           c=-c
  !           endif

  c0=-c
  c1=3.5d0*c/dd1

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
  end do

  do i=i0,indtm
    z(indorbp,i)=(c0*r(i)**3+c1*r(i)**2)*distp(i,1)
  end do

  if(typec.ne.1) then
    rp1=r(0)**3
    rp2=r(0)**2

    !              fun=(2.d0-dd1*r(0))*distp(0,1)
    !              fun2=(2.d0-4*dd1*r(0)+(dd1*r(0))**2)*distp(0,1)
    !
    !c              the first derivative/r
    fun=distp(0,1)*(c0*(3.d0*r(0)-dd1*rp2)                           &
        +c1*(2.d0-dd1*r(0)))

    !c

    !c              the second derivative
    fun2=distp(0,1)*                                                 &
        (c0*(6.d0*r(0)-6.d0*dd1*rp2+dd1**2*rp1)                      &
        +c1*(2.d0-4*dd1*r(0)+(dd1*r(0))**2))
    !c
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*fun+fun2

  end if
  !
  indorb=indorbp
  !
  !        endif
  indpar=indpar+1
  indshell=indshellp
  !
  ! 3p single zeta
