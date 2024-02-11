  ! s orbital
  ! 
  ! - angmom = 1
  ! - type = Slater
  ! - normalized = yes
  ! - angtype = spherical
  ! - npar = 5
  ! - multiplicity = 3
  !
  ! 2p with cusp conditions
  !

  dd1=dd(indpar+1)
  dd2=dd(indpar+2)
  peff=dd(indpar+3)
  dd3=dd(indpar+4)
  peff2=dd(indpar+5)

  c=1.d0/2.d0/dsqrt(8.d0*pi*(1.d0/(2.d0*dd1)**5                      &
      +2.d0*peff/(dd1+dd2)**5+peff**2/(2.d0*dd2)**5                  &
      +2.d0*peff2/(dd1+dd3)**5+peff2**2/(2.d0*dd3)**5                &
      +2.d0*peff2*peff/(dd2+dd3)**5))

  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k))
    distp(k,2)=c*dexp(-dd2*r(k))
    distp(k,3)=c*dexp(-dd3*r(k))
  end do

  do i=indtmin,indtm
    distp(i,4)=distp(i,1)+peff*distp(i,2)+peff2*distp(i,3)
  end do

  do ic=1,3
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,4)
    end do
  end do

  if(typec.ne.1) then
    fun=(-dd1*distp(0,1)-dd2*peff*distp(0,2)                         &
        -dd3*peff2*distp(0,3))/r(0)
    fun2=dd1**2*distp(0,1)+peff*dd2**2*distp(0,2)                    &
        +peff2*dd3**2*distp(0,3)

    do ic=1,3
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*fun
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+distp(0,4)
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
    end do
  end if

  indpar=indpar+5
  indshell=indshell+3
  indorb=indorbp
