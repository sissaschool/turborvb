!     2s double exp  with constant and cusp cond.
!     (dd3+ exp (-dd2 r)*(1+dd2*r)+dd4*exp(-dd5*r)*(1+dd5*r))



dd2=dd(indpar+1)
dd3=dd(indpar+2)
dd4=dd(indpar+3)
dd5=dd(indpar+4)

indorbp=indorb+1
indshellp=indshell+1
do k=indtmin,indtm
  distp(k,3)=dexp(-dd2*r(k))
  distp(k,4)=dexp(-dd5*r(k))
  distp(k,1)=distp(k,3)*(1.d0+dd2*r(k))
  distp(k,2)=distp(k,4)*(1.d0+dd5*r(k))
end do

!           if(iocc(indshellp).eq.1) then
do i=i0,indtm
  z(indorbp,i)=distp(i,1)+dd3+dd4*distp(i,2)
!              write(6,*) ' function inside = ',z(indorbp,i)
end do
!           endif


if(typec.ne.1) then

  fun=-dd2**2*distp(0,3)-dd5**2*dd4*distp(0,4)
  fun2=-dd2**2*distp(0,3)*(1.d0-dd2*r(0))                            &
      -dd4*dd5**2*distp(0,4)*(1.d0-dd5*r(0))



  do i=1,3
    z(indorbp,indt+i)=fun*rmu(i,0)
  end do

  z(indorbp,indt+4)=2.d0*fun+fun2


  !endif for indt
end if

indpar=indpar+4
indshell=indshellp
indorb=indorbp

