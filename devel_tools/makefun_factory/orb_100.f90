  !     2s single gaussian
  !     exp(-dd2*r^2)


  dd2=dd(indpar+1)


  indorbp=indorb+1
  indshellp=indshell+1
  do k=indtmin,indtm
    distp(k,1)=dexp(-dd2*r(k)*r(k))
  end do

  !           if(iocc(indshellp).eq.1) then
  do i=i0,indtm
    z(indorbp,i)=distp(i,1)
  end do
  !           endif


  if(typec.ne.1) then
    fun=-dd2*distp(0,1)*2.d0
    do i=1,3
      z(indorbp,indt+i)=fun*rmu(i,0)
    end do

    z(indorbp,indt+4)=2.d0*dd2*(-3.d0+2.d0*dd2*r(0)**2)*             &
        distp(0,1)
    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshellp
  indorb=indorbp



