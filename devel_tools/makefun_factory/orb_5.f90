  ! normalized

  indshellp=indshell+1


  !        if(iocc(indshellp).eq.1) then

  indorbp=indorb+1
  dd1=dd(indpar+1)

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
  end do

  !           if(iflagnorm.gt.2) then
  !           c=dd1**2.5d0/dsqrt(3.d0*pi)
  c=dd1**2.5d0*0.32573500793527994772d0
  !           endif

  do i=i0,indtm
    z(indorbp,i)=c*r(i)*distp(i,1)
  end do

  if(typec.ne.1) then

    fun=distp(0,1)*(1.d0-dd1*r(0))
    fun2=distp(0,1)*(dd1**2*r(0)-2.d0*dd1)

    do i=1,3
      z(indorbp,indt+i)=c*fun*rmu(i,0)/r(0)
    end do
    z(indorbp,indt+4)=c*2.d0*fun/r(0)+c*fun2

  end if

  indorb=indorbp

  !        endif
  indpar=indpar+1
  indshell=indshellp


  ! 2s double  Z NO CUSP
