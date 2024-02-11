  !c     4p without cusp condition derivative of 22
  !c      r^2   e^{-z1 r }


  dd1=dd(indpar+1)
  !        if(iflagnorm.gt.2) then
  !        c=dsqrt((2.d0*dd1)**7/240.d0/pi)/2.d0
  c=dd1**3.5d0*0.2060129077457011d0
  !        endif
  c0=-c

  c1=3.5d0*c/dd1

  do k=indtmin,indtm
    distp(k,1)=dexp(-dd1*r(k))
  end do

  do i=indtmin,indtm
    distp(i,3)=r(i)**2*distp(i,1)
  end do

  !        indorbp=indorb

  do ic=1,3
    !           if(iocc(indshell+ic).eq.1) then
    indorbp=indorb+ic
    do i=i0,indtm
      z(indorbp,i)=rmu(ic,i)*(c0*distp(i,3)+c1*r(i)*distp(i,1))
    end do
    !           endif
  end do


  if(typec.ne.1) then
    !           fun=(1.d0-dd1*r(0))*distp(0,1)
    !           fun2=dd1*(dd1*r(0)-2.d0)*distp(0,1)

    fun0=c0*distp(0,3)+c1*r(0)*distp(0,1)
    fun=(c0*(2.d0-dd1*r(0))*r(0)                                     &
        +c1*(1.d0-dd1*r(0)))*distp(0,1)
    fun2=(c0*((dd1*r(0))**2+2.d0-4.d0*dd1*r(0))                      &
        +c1*dd1*(dd1*r(0)-2.d0))*distp(0,1)

    !              indorbp=indorb
    do ic=1,3
      !                if(iocc(indshell+ic).eq.1) then
      indorbp=indorb+ic
      do i=1,3
        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*                        &
            fun/r(0)
        if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0
      end do
      z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun/r(0)+fun2)
      !                 endif
    end do
  end if

  indpar=indpar+1
  indshell=indshell+3
  indorb=indorbp

