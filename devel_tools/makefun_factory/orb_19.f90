  ! R(r)=c*exp(-z*r**2)*(3/4/z-r**2)




  indshellp=indshell+1

  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1

  dd1=dd(indpar+1)

  !         if(iflagnorm.gt.2) then
  !          if(dd1.ne.0.) then
  !            c=(2.d0*dd1/pi)**(3.d0/4.d0)
  c=0.71270547035499016d0*dd1**0.75d0
  !          else
  !          c=1.d0
  !          endif
  !         endif

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  do i=i0,indtm
    z(indorbp,i)=distp(i,1)*(3.d0/4.d0/dd1-r(i)**2)
  end do

  if(typec.ne.1) then
    !              the first derivative /r
    fun=distp(0,1)*(2.d0*dd1*r(0)**2-7.d0/2.d0)

    !              the second derivative
    fun2=distp(0,1)*(-4.d0*dd1**2*r(0)**4                            &
        +13.d0*dd1*r(0)**2-7.d0/2.d0)

    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*fun+fun2

  end if

  indorb=indorbp

  !        endif

  indpar=indpar+1
  indshell=indshellp



  ! 2p single zeta
