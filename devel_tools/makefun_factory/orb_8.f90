  ! s orbital
  ! 
  ! - angmom = 0
  ! - type = Gaussian
  ! - normalized = yes
  ! - angtype = spherical
  ! - npar = 2
  ! - multiplicity = 1
  !
  ! = exp(-dd1*r) + (dd1-zeta) * r * exp(-dd2*r)
  !
  ! 2s double Z WITH CUSP
  !

  indshellp=indshell+1
  indorbp=indorb+1
  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd1-zeta(1)

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
    distp(k,2)=dexp(-dd2*r(k))
  end do

  c= 1.d0/dsqrt(1.d0/4.d0/dd1**3+12.d0*peff/(dd1+dd2)**4&
  &+ 3*peff**2/4/dd2**5)/dsqrt(4.0*pi)

  do i=i0,indtm
    z(indorbp,i)=c*(distp(i,1)+r(i)*distp(i,2)*peff)
  end do

  if(typec.ne.1) then
    fun=-dd1*distp(0,1)+peff*distp(0,2)*(1.d0-dd2*r(0))
    fun2=distp(0,1)*dd1**2&
      &+ peff*distp(0,2)*(dd2**2*r(0)-2.d0*dd2)
    do i=1,3
      z(indorbp,indt+i)=fun*c*rmu(i,0)/r(0)
    end do
    z(indorbp,indt+4)=c*(2.d0*fun/r(0)+fun2)
  end if

  indorb=indorbp
  indpar=indpar+2
  indshell=indshellp
