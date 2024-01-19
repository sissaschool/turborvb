!     2p single exponential  r e^{-z r^2}



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k)**2)
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)*r(i)
  end do
  !           endif
end do


if(typec.ne.1) then

  fun0=distp(0,1)*r(0)
  cost=2.d0*dd2*r(0)**2
  fun=distp(0,1)*(1.d0-cost)/r(0)
  fun2=2.d0*dd2*fun0*(cost-3.d0)
  !              indorbp=indorb
  do ic=1,3
    !                if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                          &
          fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)                                      &
        *(4.d0*fun+fun2)
    !                 endif
  end do
  !endif for indt
end if
indpar=indpar+1
indshell=indshell+3
indorb=indorbp

