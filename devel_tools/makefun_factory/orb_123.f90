!     2p  double  exp
!       dd1 * x_mu  (exp(-dd2 r)+dd3 * exp(-dd4*r))


dd2=dd(indpar+1)
dd3=dd(indpar+2)
dd4=dd(indpar+3)


do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))
  distp(k,2)=dexp(-dd4*r(k))
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*(distp(i,1)+dd3*distp(i,2))
  end do
  !           endif
end do


if(typec.ne.1) then

  fun0=distp(0,1)+dd3*distp(0,2)
  fun=-(dd2*distp(0,1)+dd3*dd4*distp(0,2))/r(0)
  fun2=dd2**2*distp(0,1)+dd3*dd4**2*distp(0,2)


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

indpar=indpar+3
indshell=indshell+3
indorb=indorbp

