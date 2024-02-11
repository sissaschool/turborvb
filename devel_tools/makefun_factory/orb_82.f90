
  dd1=dd(indpar+1)
  dd2=dsqrt(dd1)

  !        if(iflagnorm.gt.2) then
  !        ratiocp--> ratiocp*dsqrt(2.d0)*pi**(-0.75d0)*2**1.25
  !        c=dsqrt(2.d0)*pi**(-0.75d0)*(2.d0*dd1)**1.25d0*ratiocp
  c=dd1**1.25d0*ratiocp
  !        endif

  do k=indtmin,indtm
    cost=dd1*r(k)**2/(1.d0+dd2*r(k))
    distp(k,1)=c*dexp(-cost)
  end do
  ! indorbp=indorb
  !
  do ic=1,3
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*distp(i,1)
    end do
  end do

  if(typec.ne.1) then
    fun0=distp(0,1)
    ! fun=-2.d0*dd1*distp(0,1)
    ! fun2=fun*(1.d0-2.d0*dd1*r(0)**2)
    rp1=dd1*r(0)**2
    rp2=dd2*r(0)
    rp3=(1.d0+rp2)**2

    fun=-dd1*distp(0,1)*(2.d0+rp2)/rp3

    ! the second derivative
    fun2=dd1*distp(0,1)*(-2.d0-2.d0*rp2+4.d0*rp1+4.d0*rp1*rp2+rp1**2)/rp3**2
    ! indorbp=indorb
    do ic=1,3
      ! if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
    end do
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp

