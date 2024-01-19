!     p gaussian  r**(2*npower)*exp(-alpha*r**2)

npower=iopt-1100

!        indorbp=indorb

dd2=dd(indpar+1)
do k=indtmin,indtm
  distp(k,1)=r(k)**(2*npower)*dexp(-dd2*r(k)**2)
end do

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)
  end do
  !           endif
end do


if(typec.ne.1) then


  rp1=r(0)**2
  fun0=distp(0,1)
  fun=(npower-dd2*rp1)*distp(0,1)*2.d0/rp1
  fun2=(npower*(2.d0*npower-1.d0)-                                   &
      (1.d0+4.d0*npower)*dd2*rp1+2.d0*(dd2*rp1)**2)*                 &
      distp(0,1)*2.d0/rp1

  !           indorbp=indorb
  do ic=1,3
    !              if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=1,3
      z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*fun
      if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
    end do
    z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
    !              endif
  end do



  !endif for indt
end if

indpar=indpar+1
indshell=indshell+3
indorb=indorbp


