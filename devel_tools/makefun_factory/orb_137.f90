!     2s with cusp condition
!     dd1*(exp(-dd2*r)*(1+dd2*r))


dd2=dd(indpar+1)

!           if(iflagnorm.gt.2) then
!           c=1.d0/dsqrt(1/4.d0/dd2**3+12*dd2/(2.d0*dd2)**4+
!                                                                    &3*dd2**2/4/dd2**5)/dsqrt(4.0*pi)
    !           endif

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,1)=dexp(-dd2*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)*(1.d0+dd2*r(i))
end do
!           endif


if(typec.ne.1) then
  fun=-dd2**2*distp(0,1)
  fun2=fun*(1.d0-dd2*r(0))


  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*fun+fun2


  !endif for indt
end if

indpar=indpar+1
indshell=indshellp
indorb=indorbp

