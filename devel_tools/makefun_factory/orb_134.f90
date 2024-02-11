!     2p single exponential  r^3 e^{-z r}  !



dd2=dd(indpar+1)

do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))
end do

!        indorbp=indorb

do ic=1,3
  !           if(iocc(indshell+ic).eq.1) then
  indorbp=indorb+ic
  do i=i0,indtm
    z(indorbp,i)=rmu(ic,i)*distp(i,1)*r(i)**3
  end do
  !           endif
end do

if(typec.ne.1) then

  fun0=distp(0,1)*r(0)**3
  fun=distp(0,1)*(3.d0-dd2*r(0))*r(0)
  !           fun= derivative of fun0 respect to r divided dy r
  fun2=distp(0,1)*(dd2**2*r(0)**3-6*dd2*r(0)**2                      &
      +6*r(0))
  !           fun2= second derivative of fun0 respect to r
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


