


  dd1=dd(indpar+1)


  !        if(iflagnorm.gt.2) then
  !        c=2.d0/pi**0.75d0*(2.d0*dd1)**1.75d0/dsqrt(5.d0)
  c=dd1**1.75d0*1.2749263037197753d0
  !        endif


  do k=indtmin,indtm
    distp(k,1)=c*dexp(-dd1*r(k)**2)
  end do

  !        indorbp=indorb
  !
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
    cost=2.d0*dd1*r(0)**2
    fun=distp(0,1)*(1.d0-cost)/r(0)
    fun2=2.d0*dd1*fun0*(cost-3.d0)
    !              indorbp=indorb

    do ic=1,3
      !                if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2)
      !                 endif
    end do

    !endif for indt
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp

  ! derivative of 62 with respect zeta
